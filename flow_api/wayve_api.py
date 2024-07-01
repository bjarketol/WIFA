# General packages
import numpy as np
from scipy.interpolate import interp1d
from windIO.utils.yml_utils import load_yaml
import matplotlib.pyplot as plt
import mpmath
import warnings
import xarray as xr

# Atmospheric state setup
from wayve.abl.abl import ABL
from wayve.abl.abl_tools import Cg_cubic, alpha_cubic
from wayve.abl import ci_methods

# General APM setup
from wayve.apm import APM
from wayve.grid.grid import Stat2Dgrid
from wayve.forcing.wind_farms.wake_model_coupling.coupling_methods.pressure_based import PressureBased
from wayve.forcing.wind_farms.wake_model_coupling.wake_models.lanzilao_merging import UniDirectional
from wayve.forcing.wind_farms.wind_farm import WindFarm, Turbine
from wayve.forcing.apm_forcing import ForcingComposite
from wayve.momentum_flux_parametrizations import FrictionCoefficients
from wayve.pressure.gravity_waves.gravity_waves import Uniform
from wayve.solvers import FixedPointIteration


def run_wayve(yamlFile):

    ######################
    # construct APM grid
    ######################
    # Numerical parameters (values from Allaerts and Meyers, 2019)
    Nx = 2000   # grid points in x-direction
    Lx = 1.e6   # grid size in x-direction [m]
    Ny = 2000   # grid points in y-direction
    Ly = 1.e6   # grid size in y-direction [m]
    # Generate 2D grid object
    grid = Stat2Dgrid(Lx, Nx, Ly, Ny)

    #####################
    # Read out yaml file
    #####################
    # Yaml loading
    system_dat = load_yaml(yamlFile)
    # WindIO components
    farm_dat = system_dat['wind_farm']
    resource_dat = system_dat['site']['energy_resource']

    ####################
    # Set up WindFarm object
    ####################
    # Turbine geometry
    hh = farm_dat['turbines']['hub_height']
    rd = farm_dat['turbines']['rotor_diameter']
    # Ct curve data
    ct = farm_dat['turbines']['performance']['Ct_curve']['Ct_values']
    ct_ws = farm_dat['turbines']['performance']['Ct_curve']['Ct_wind_speeds']
    # Cp curve data
    air_density = 1.225     # Hard-coded for now
    if 'Cp_curve' in farm_dat['turbines']['performance']:
        # Read out Cp curve
        cp = farm_dat['turbines']['performance']['Cp_curve']['Cp_values']
        cp_ws = farm_dat['turbines']['performance']['Cp_curve']['Cp_wind_speeds']
        power_curve_type = 'cp'
    elif 'power_curve' in farm_dat['turbines']['performance']:
        # Convert power curve to Cp curve
        cp_ws = farm_dat['turbines']['performance']['power_curve']['power_wind_speeds']
        pows = farm_dat['turbines']['performance']['power_curve']['power_values']
        rotor_area = np.pi * (rd / 2) ** 2
        cp = np.divide(np.array(pows), 0.5 * air_density * np.array(cp_ws) ** 3 * rotor_area)
    else:
        raise Exception('Bad Power Curve')
    # Ct and Cp curves
    ct_curve = interp1d(ct_ws, ct, fill_value="extrapolate")
    cp_curve = interp1d(cp_ws, cp, fill_value="extrapolate")
    # Get x and y positions
    x = farm_dat['layouts']['initial_layout']['coordinates']['x']
    y = farm_dat['layouts']['initial_layout']['coordinates']['y']
    # Reposition to be at grid center
    x -= np.mean(x)
    y -= np.mean(y)
    # Number of turbines
    Nt = len(x)
    # Turbine setup
    turbines = []
    for t in range(Nt):
        turbine = Turbine(x[t], y[t], rd, hh, ct_curve, cp_curve)
        turbines.append(turbine)
    # Use wake merging method of Lanzilao and Meyers (2021)
    wakemodel = UniDirectional()
    # Set up coupling object
    # Use pressure-based method
    coupling = PressureBased(wakemodel)
    # Gaussian filter length (meaningless for uncoupled wake model run)
    Lfilter = 1000.
    # Generate wind farm object
    wind_farm = WindFarm(turbines, Lfilter, coupling)
    # Combined forcing object
    forcing = ForcingComposite([wind_farm])

    ##################
    # Read site data
    ##################
    # Note: we assume all variables are given, and do not set up any default behavior
    # Get raw data
    times = resource_dat['wind_resource']['time']

    #####################
    # Perform model runs
    #####################
    # Initialize crash counter
    crashes = 0
    # List of datasets
    ds_list = []
    # Loop over timeseries
    for time_index, time in enumerate(times):
        try:
            # Set up ABL
            abl = flow_io_abl(resource_dat['wind_resource'], time_index, hh)
            # Set up APM from components
            # Momentum flux parametrization
            mfp = FrictionCoefficients()
            # Pressure feedback parametrization
            pressure = Uniform(dynamic=True, rotating=False)
            # Create static 2D model
            model = APM(grid, forcing, abl, mfp, pressure)
            # Use a fixed-point iteration solver with a relaxation factor of 0.7
            solver = FixedPointIteration(tol=5.e-3, relax=0.7)
            # Solve model
            _ = model.solve(method=solver)  # APM run
            # wind_farm.preprocess(model)   # Wake model run
            # Evaluate wind-farm power output
            turbine_power = wind_farm.power_turbines(abl.rho)
            # NC setup
            ds = xr.Dataset(
                {
                    "power": ("turbine", turbine_power),
                    "rotor_effective_velocity": ("turbine", wind_farm.coupling.St)
                },
                coords={
                    "time": time,
                    "turbine": range(Nt)
                },
            )
            # Add to output list
            ds_list.append(ds)
        except Exception as exc:
            print(exc)
            # Update crash counter
            crashes += 1
            continue
    #     # Print timestep
    #     print(f"time {time_index}/{len(times)}")
    # print(f"crashes: {crashes}/{len(times)}")

    # Combine into total dataset
    ds_full = xr.concat(ds_list, dim="time")
    # ds_full = ds_full.fillna(0.)
    ds_full.to_netcdf("turbine_data.nc")

    return


def nieuwstadt83_profiles(zh, v, wd, z0=1.e-1, h=1.5e3, fc=1.e-4, ust=0.666):
    """Set up the cubic analytical profile from Nieuwstadt (1983), based on hub height velocity information"""
    # Constants #
    kappa = 0.41  # Von Karman constant
    # # We iterate until we find a profile that has the requested speed at hub height, by varying ust # #
    # Iteration settings
    ust_i = ust
    error = np.infty
    attempt = 0
    max_attempts = 30
    tolerance = 1.e-3
    # Iteration
    while error > tolerance and attempt < max_attempts:
        # # # Nieuwstadt solution # # #
        # Vertical grid
        Nz = 100
        zs = np.linspace(z0, h, Nz)
        # Dimensionless groups
        hstar = h*fc/ust_i
        z0_h = z0/h
        # Nieuwstadt relations
        Cg = Cg_cubic(hstar, z0_h, kappa)               # Geostrophic drag Cg = utau/G
        geo_angle = alpha_cubic(hstar, z0_h, kappa)     # Geostrophic wind angle
        # Nieuwstadt solution #
        C = h * fc / kappa / ust_i
        alpha = 0.5 + 0.5 * np.sqrt(1 + 4j * C)
        sigma_s = np.zeros(len(zs), dtype=np.complex128)
        wd_s = np.zeros(len(zs), dtype=np.complex128)
        with np.errstate(invalid='ignore'):  # z>=h will result in Nan. This is set to 0 below.
            for k in range(len(zs)):
                sigma_s[k] = alpha * (mpmath.gamma(alpha)) ** 2 / mpmath.gamma(2 * alpha) \
                             * np.power(1. - zs[k] / h, alpha) \
                             * mpmath.hyp2f1(alpha - 1, alpha, 2 * alpha, 1 - zs[k] / h)
                wd_s[k] = (1j * alpha ** 2 * (mpmath.gamma(alpha)) ** 2) / (kappa * C * mpmath.gamma(2 * alpha)) * (
                        1 - zs[k] / h) ** (alpha - 1) * mpmath.hyp2f1(alpha + 1, alpha - 1, 2 * alpha,
                                                                      1 - zs[k] / h)
        # Set Nan to 0
        sigma_s[np.isnan(sigma_s)] = np.complex128(0.)
        wd_s[np.isnan(wd_s)] = np.complex128(0.)
        # Velocity arrays
        us = ((Cg**-1)*np.cos(geo_angle) + np.real(wd_s)) * ust_i
        vs = ((Cg**-1)*np.sin(geo_angle) + np.imag(wd_s)) * ust_i
        # Error
        u_hh = np.interp(zh, zs, np.sqrt(np.square(us)+np.square(vs)))
        error = np.abs(u_hh - v) / v
        ust_i *= v / u_hh
        attempt += 1
    # Velocity arrays
    us = ((Cg**-1)*np.cos(geo_angle) + np.real(wd_s)) * ust_i
    vs = ((Cg**-1)*np.sin(geo_angle) + np.imag(wd_s)) * ust_i
    # Momentum flux arrays
    tauxs = np.real(sigma_s) * ust_i ** 2
    tauys = np.imag(sigma_s) * ust_i ** 2
    nus = kappa * ust * np.multiply(zs, (1 - zs / h) ** 2, out=np.zeros_like(zs), where=(zs <= h))
    # # Rotate to match wind direction at hub height # #
    # Current wind direction (angle w.r.t. x-axis)
    wd_hh_0 = np.arctan2(np.interp(zh, zs, vs), np.interp(zh, zs, us))
    # Rotation angle
    rotation_angle = -(wd_hh_0 + np.deg2rad(wd) + np.pi/2.)     # +pi/2 for wd convention
    # Velocity components
    us, vs = rotate_xy_arrays(us, vs, rotation_angle)
    tauxs, tauys = rotate_xy_arrays(tauxs, tauys, rotation_angle)
    # Upper atmosphere
    U3 = us[-1]
    V3 = vs[-1]
    return zs, us, vs, U3, V3, tauxs, tauys, nus


def rotate_xy_arrays(xs, ys, angle):
    """
    Rotate the given vectors around the given angle.

    Parameters
    ----------
    xs : array_like
        x-components of the vectors
    ys : array_like
        y-components of the vectors
    angle : float
        angle over which to rotate the vectors (in radians)
    """
    # Angle cosine and sine
    c, s = np.cos(angle), np.sin(angle)
    # Rotation matrix
    R = np.array(((c, -s),
                  (s, c)))
    # Output arrays
    xs_rot, ys_rot = 0.*xs, 0.*ys
    # Loop over vectors
    Ns = len(xs)
    for i in range(Ns):
        # Vector i
        vec = np.array([xs[i], ys[i]])
        # Multiply with rotation matrix
        rotated_vec = np.matmul(R, vec)
        # Store rotated vector
        xs_rot[i] = rotated_vec[0]
        ys_rot[i] = rotated_vec[1]
    return xs_rot, ys_rot


def ci_fitting(zs, ths, l_mo=5.e3, blh=1.e3, dh_max=300., serz=True, plot_fits=False):
    # Stable or unstable atmosphere
    stable = 0. < l_mo < 100
    # Estimate inversion parameters with RZ fit #
    # Relevant part of the vertical profiles
    max_z_fit = 5.e3
    z_ci = zs[zs <= max_z_fit]
    th_ci = ths[zs <= max_z_fit]
    # Surface-Extended RZ or regular RZ
    if serz:
        # Stable or unstable profile determines the initial guess for the CI height
        if stable:
            l_p0 = 1.e3
        else:
            l_p0 = blh
        # Initial estimate for MBL temperature in fit
        th_mbl = np.interp(l_p0, z_ci, th_ci)
        # Fitting procedure
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            ci_estimate = ci_methods.SERZ_fit(z_ci, th_ci,
                                              p0=[0.9, .1, 0., th_mbl, l_p0, 100, .05],
                                              initialGuess='RZ',
                                              dh_max=dh_max)
    else:
        # Stable or unstable profile
        if stable:
            # Ignore temperature decrease inside SBL
            # (we are trying to identify the mixing layer that preceded this SBL)
            # We want to capture the mixed layer (or residual layer since we are in
            # SBL) that preceded the SBL. Therefore, we take the potential temperature
            # at the top of the ABL where we have the mixed layer and we extrapolate
            # till the bottom. We use this constant value in the ABL.
            # p0 are the initial guess for [a,b,thm,l,dh] used in Ramp&Zar model
            th_ci[z_ci < blh] = np.interp(blh, z_ci, th_ci)
            l_p0 = 1.e3
        else:
            # Ignore temperature increase in CBL surface layer, therefore we take
            # the lowest value of potential temperature. We are able to capture the
            # mixed layer in this way, where the temperature is constant and equal
            # to the lowest theta. We use this constant value in the ABL.
            # p0 are the initial guess for [a,b,thm,l,dh] used in Ramp&Zar model
            th_ci[0:np.argmin(th_ci)] = np.min(th_ci)
            l_p0 = blh
        # Initial estimate for MBL temperature in fit
        th_mbl = np.interp(l_p0, z_ci, th_ci)
        # RZ fit
        ci_estimate = ci_methods.RZfit(z_ci, th_ci,
                                       p0=[0.9, 0.1, th_mbl, l_p0, 100.0],
                                       dh_max=dh_max)
    # Plot fitted potential temperature profile
    if plot_fits:
        fig, ax = plt.subplots()
        ax.plot(ths[zs <= max_z_fit], z_ci / 1.e3, 'b', label="Data")
        ax.plot(ci_estimate['thfit'], z_ci / 1.e3, '--k', label="RZ fit", zorder=-1)
        ax.set_xlim([285., 312.])
        ax.set_ylim([0., 4.])
        ax.set_ylabel('$z$ [km]')
        ax.set_xlabel('$\\theta$ [K]')
        plt.legend()
        plt.tight_layout()
        plt.show()
    # CI altitudes
    inv_bottom = ci_estimate['h0']
    H = ci_estimate['h1']
    inv_top = ci_estimate['h2']
    # Determine reference potential temperature
    th0 = np.interp(H, zs, ths)
    # Inversion strength
    if ci_estimate['a'] <= 0.2 or ci_estimate['a'] <= 2 * ci_estimate['b']:
        # No inversion strength in the following cases:
        # a<=0.2: encroachment (No inversion layer, so the entire profile is given by g and a=0
        #           (considered a,0.2 as in paper))
        # a<=2*b: inversion lapse rate is equal to or smaller than free lapse rate
        dth = 0.
    else:
        dth = ci_estimate['dth']
    # Lapse rate
    dthdz = ci_estimate['gamma']
    return inv_bottom, H, inv_top, th0, dth, dthdz


def flow_io_abl(wind_resource_dat, time_index, zh, dh_max=None, serz=True):
    '''
    Method to set up an ABL object based on FLOW IO

    Parameters
    ----------
    wind_resource_dat: dict
        Wind resource data
    time_index: int
        Index of the timestamp to set up ABL for
    zh: float
        Average turbine hub height
    dh_max (optional): float
        Maximum depth of the inversion layer used in the inversion curve fitting procedure (default: None)
    serz (optional): boolean
        Whether the surface-extended version of the RZ model is used (default: True)
    '''
    # Constants #
    gravity = 9.80665  # [m s-2]
    kappa = 0.41  # Von Karman constant
    omega = 7.2921159e-5  # angular speed of the Earth [rad/s]
    # Basic atmospheric scalars #
    air_density = 1.225     # Hard-coded for now
    # Surface roughness
    z0 = 1.e-1
    if 'z0' in wind_resource_dat.keys():
        z0 = wind_resource_dat['z0']["data"][time_index]
    # Monin-Obhukov length
    l_mo = 5.e3
    if 'LMO' in wind_resource_dat.keys():
        l_mo = wind_resource_dat['LMO']["data"][time_index]
    # Coriolis parameter #
    phi = 0.377     # Assume latitude location
    fc = 2 * omega * np.sin(phi)
    if 'fc' in wind_resource_dat.keys():
        fc = wind_resource_dat['fc']["data"][time_index]
    # Check if wind resource contains vertical profile
    profile_input = 'height' in wind_resource_dat.keys()
    if not profile_input:
        # Wind speed and direction
        v = wind_resource_dat['wind_speed']["data"][time_index]
        wd = wind_resource_dat['wind_direction']["data"][time_index]
        # Friction velocity
        ust = 0.666
        if 'friction_velocity' in wind_resource_dat.keys():
            ust = wind_resource_dat['friction_velocity']["data"][time_index]
        # Turbulence intensity
        TI = 0.04
        if 'z0' in wind_resource_dat.keys():
            TI = wind_resource_dat['turbulence_intensity']["data"][time_index] / 100.
        # Capping inversion information
        h = 1.5e3
        dh = 100.
        dth = 5.
        dthdz = 2.e-3
        th0 = 293.15
        if 'thermal_stratification' in wind_resource_dat.keys():
            thermal_data = wind_resource_dat['thermal_stratification']
            if 'capping_inversion' in thermal_data.keys():
                ci_data = thermal_data['capping_inversion']
                h = ci_data['ABL_height']["data"][time_index]
                dh = ci_data['dH']["data"][time_index]
                dth = ci_data['dtheta']["data"][time_index]
                dthdz = ci_data['lapse_rate']["data"][time_index]
        inv_bottom, inv_top = h - dh/2, h + dh/2
        # Nieuwstadt profiles for velocity and shear stress
        zs, us, vs, U3, V3, tauxs, tauys, nus = nieuwstadt83_profiles(zh, v, wd, z0=z0, h=h, ust=ust, fc=fc)
        # Potential temperature profile constant
        ths = th0 * np.ones_like(zs)
    else:
        # Read out vertical profile
        zs = np.array(wind_resource_dat['height'])
        vs = np.array(wind_resource_dat['wind_speed']["data"][time_index])
        wds = np.array(wind_resource_dat['wind_direction']["data"][time_index])
        ths = np.array(wind_resource_dat['potential_temperature']["data"][time_index])
        TIs = np.array(wind_resource_dat['turbulence_intensity']["data"][time_index])
        # Interpolate TI
        TI = np.interp(zh, zs, TIs)
        # Velocity components
        us = -vs * np.sin(np.deg2rad(wds))
        vs = -vs * np.cos(np.deg2rad(wds))
        # Check available inputs
        if 'k' in wind_resource_dat.keys():     # RANS-like inputs
            tkes = np.array(wind_resource_dat['k']["data"][time_index])
            eps = np.array(wind_resource_dat['epsilon']["data"][time_index])
            # Eddy viscosity
            C_mu = 0.09     # k-epsilon model value
            nus = C_mu * np.divide(np.square(tkes), eps, out=np.zeros_like(tkes), where=eps!=0)
            # Momentum fluxes
            dudz = np.gradient(us, zs, edge_order=2)
            dvdz = np.gradient(vs, zs, edge_order=2)
            tauxs = nus * dudz
            tauys = nus * dvdz
        else:   # Shear stress profile directly available
            tauxs = np.array(wind_resource_dat['tau_x']["data"][time_index])
            tauys = np.array(wind_resource_dat['tau_y']["data"][time_index])
            nus = None
        # Total momentum flux
        taus = np.sqrt(np.square(tauxs) + np.square(tauys))
        # Friction velocity
        ust = taus[0]   # Assume friction velocity is not given explicitly
        # Estimate boundary layer height based on momentum flux #
        f_tau = interp1d(taus, zs)
        blh = f_tau(0.01 * ust)
        # Capping inversion information
        if ('thermal_stratification' in wind_resource_dat.keys() and
                'capping_inversion' in wind_resource_dat['thermal_stratification'].keys()):
            thermal_data = wind_resource_dat['thermal_stratification']
            ci_data = thermal_data['capping_inversion']
            th0 = 293.15
            h = ci_data['ABL_height']["data"][time_index]
            dh = ci_data['dH']["data"][time_index]
            dth = ci_data['dtheta']["data"][time_index]
            dthdz = ci_data['lapse_rate']["data"][time_index]
            inv_bottom, inv_top = h - dh/2, h + dh/2
        else:
            inv_bottom, h, inv_top, th0, dth, dthdz = ci_fitting(zs, ths, l_mo, blh, dh_max=dh_max, serz=serz)
        # Geostrophic wind speed
        z = np.linspace(h, 15.e3, 1000)
        U3 = np.trapz(np.interp(z, zs, us), z) / (15.e3 - h)
        V3 = np.trapz(np.interp(z, zs, vs), z) / (15.e3 - h)
    # Lower layer thickness
    h1 = 2*zh
    # Upper layer thickness
    h2 = h - h1
    if inv_bottom <= h1 + 10.:  # H cannot be lower than H1 and the upper layer must be at least 10m
        raise RuntimeWarning(f'CI too low, CI bottom located at z={int(inv_bottom)}m')
    # CI check
    if dth == 0.:
        raise RuntimeWarning("No CI present!")
    # gprime and N
    gprime = gravity * dth / th0
    N = np.sqrt(gravity * dthdz / th0)
    # Set up ABL object
    return ABL(zs, us, vs, ths, tauxs, tauys,
               h1, h2,
               gprime, N, U3, V3,
               fc,
               nus=nus,
               rho=air_density, TI=TI, z0=z0, ust=ust,
               inv_bottom=inv_bottom, inv_top=inv_top)


if __name__ == '__main__':
    # Basic run using local file, for local testing #
    # Get given file location
    import argparse
    arg_parser = argparse.ArgumentParser()
    arg_parser.add_argument('-file', dest='file')
    args = arg_parser.parse_args()
    # Run wayve with given file
    run_wayve(args.file)

