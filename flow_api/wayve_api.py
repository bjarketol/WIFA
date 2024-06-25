# General packages
import numpy as np
import netCDF4
from scipy.interpolate import interp1d
from windIO.utils.yml_utils import load_yaml
import matplotlib.pyplot as plt
import pandas as pd

# Atmospheric state setup
from wayve.abl.abl import ABL
from wayve.abl import line_optimization, ci_methods
from wayve.abl.abl_tools import height_average

# General APM setup
from wayve.apm import APM
from wayve.grid.grid import Stat2Dgrid
from wayve.forcing.wind_farms.wake_model_coupling.coupling_methods.pressure_based import PressureBased
from wayve.forcing.wind_farms.wake_model_coupling.wake_models.lanzilao_merging import UniDirectional
from wayve.forcing.wind_farms.wind_farm import WindFarm, Turbine
from wayve.forcing.wind_farms.dispersive_stresses import DispersiveStresses
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
    zs = np.array(resource_dat['wind_resource']['height'])
    vs_list = resource_dat['wind_resource']['wind_speed']
    wds_list = resource_dat['wind_resource']['wind_direction']
    ks_list = resource_dat['wind_resource']['k']
    eps_list = resource_dat['wind_resource']['epsilon']
    ths_list = resource_dat['wind_resource']['potential_temperature']
    TIs_list = resource_dat['wind_resource']['turbulence_intensity']
    z0s_list = resource_dat['wind_resource']['z0']
    lmos_list = resource_dat['wind_resource']['LMO']

    #####################
    # Perform model runs
    #####################
    # Set up output list
    output_list = []
    # Initialize AEP
    aep = 0.
    # Initialize crash counter
    crashes = 0
    # Loop over timeseries
    for time_index, time in enumerate(times):
        # Select data
        vs = np.array(vs_list["data"][time_index])
        wds = np.array(wds_list["data"][time_index])
        ks = np.array(ks_list["data"][time_index])
        eps = np.array(eps_list["data"][time_index])
        ths = np.array(ths_list["data"][time_index])
        TIs = np.array(TIs_list["data"][time_index])
        z0 = z0s_list["data"][time_index]
        # Set up ABL
        H1 = 2 * hh
        try:
            abl = flow_io_abl(zs, vs, wds, ks, eps, ths, TIs, z0, H1=H1)
        except Exception as exc:
            print(exc)
            # Add to output list
            output_index = [time] + [None] * Nt
            output_list.append(output_index)
            # Update crash counter
            crashes += 1
            continue
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
        _ = model.solve(method=solver)
        # Evaluate wind-farm power output
        turbine_power = wind_farm.power_turbines(air_density)
        # Add to aep
        aep += turbine_power.sum()  # Not weighted by duration of timeseries
        # Add to output list
        output_index = [time] + turbine_power.tolist()
        output_list.append(output_index)
        # Print timestep
    #     print(f"time {time_index}/{len(times)}")
    # print(crashes)
    # print(aep)

    ##############
    # Save output
    ##############
    # Set up header
    header = ["# Time"]
    for t in range(Nt):
        header += [f"Power of turbine {t}"]
    # Set up dataframe
    output = pd.DataFrame(output_list, columns=header)
    # Save output as .csv
    output.to_csv("output_wayve", na_rep='NULL')

    return aep


def flow_io_abl(zs, Ms, wds, tke, eps, ths, TIs, z0,
                H1=240., wth=0., dh_max=None, serz=True, fit_trop=False, plot=False):
    '''
    Method to set up an ABL object based on FLOW IO

    Parameters
    ----------
    H1 (optional): float
        Height of the wind-farm layer (default: 240.)
    wth (optional): float
        Surface heat flux (default: 0.)
    dh_max (optional): float
        Maximum depth of the inversion layer used in the inversion curve fitting procedure (default: None)
    serz (optional): boolean
        Whether the surface-extended version of the RZ model is used (default: True)
    fit_trop (optional): boolean
        Whether or not the upper atmospheric theta profile is fitted (default: False)
    plot (optional): boolean
        Whether or not a plot of the upper atmosphere fitting procedure is shown (default: False)
    '''
    # Constants #
    gravity = 9.80665  # [m s-2]
    kappa = 0.41  # Von Karman constant
    omega = 7.2921159e-5  # angular speed of the Earth [rad/s]
    # Velocity components
    us = -Ms * np.sin(np.deg2rad(wds))
    vs = -Ms * np.cos(np.deg2rad(wds))
    # Eddy viscosity
    C_mu = 0.09     # k-epsilon model value
    nu = C_mu * np.divide(np.square(tke), eps, out=np.zeros_like(tke), where=eps!=0)
    # Momentum fluxes
    dudz = np.gradient(us, zs, edge_order=2)
    dvdz = np.gradient(vs, zs, edge_order=2)
    tauxs = nu * dudz
    tauys = nu * dvdz
    taus = np.sqrt(np.square(tauxs) + np.square(tauys))
    # Surface layer #
    # Friction velocity
    utau = taus[0]     # Assume friction velocity is not given explicitly
    # Monin-Obukhov theory (currently used for RZ fit, default behaviour / additional inputs should be discussed)
    T0 = ths[0]     # Surface temperature taken to be lowest theta datapoint
    zeta = -2.0 * kappa * gravity * wth / (T0 * utau ** 3)
    # Estimate boundary layer height based on momentum flux #
    f_tau = interp1d(taus, zs)
    blh = f_tau(0.01 * utau)
    fTi = interp1d(zs, TIs)
    TI_h = fTi(H1 / 2.0).item()
    # Coriolis parameter #
    phi = 0.377     # Assume latitude location
    fc = 2 * omega * np.sin(phi)
    # Estimate inversion parameters with RZ fit #
    # Relevant part of the vertical profiles
    zCI = zs[zs < 5000]
    thCI = ths[zs < 5000]
    # Surface-Extended RZ or regular RZ
    if serz:
        # Stable or unstable profile determines the initial guess for the CI height
        if zeta > 0.02:
            l_p0 = 1.e3
        else:
            l_p0 = blh
        CIestimate = ci_methods.SERZ_fit(zCI, thCI,
                                         p0=[0.9, .1, 0., T0, l_p0, 100, .05],
                                         initialGuess='RZ',
                                         dh_max=dh_max)
    else:
        # Stable or unstable profile
        if zeta > 0.02:
            # Ignore temperature decrease inside SBL
            # (we are trying to identify the mixing layer that preceded this SBL)
            # We want to capture the mixed layer (or residual layer since we are in
            # SBL) that preceded the SBL. Therefore, we take the potential temperature
            # at the top of the ABL where we have the mixed layer and we extrapolate
            # till the bottom. We use this constant value in the ABL.
            # p0 are the initial guess for [a,b,thm,l,dh] used in Ramp&Zar model
            thCI[zCI < blh] = interp1d(zCI, thCI)(blh)
            l_p0 = 1.e3
        else:
            # Ignore temperature increase in CBL surface layer, therefore we take
            # the lowest value of potential temperature. We are able to capture the
            # mixed layer in this way, where the temperature is constant and equal
            # to the lowest theta. We use this constant value in the ABL.
            # p0 are the initial guess for [a,b,thm,l,dh] used in Ramp&Zar model
            thCI[0:np.argmin(thCI)] = np.min(thCI)
            l_p0 = blh
        # RZ fit
        CIestimate = ci_methods.RZfit(zCI, thCI,
                                      p0=[0.9, 0.1, T0, l_p0, 100.0],
                                      dh_max=dh_max)
    # Inversion strength
    if CIestimate['a'] <= 0.2 or CIestimate['a'] <= 2 * CIestimate['b']:
        # No inversion strength in the following cases:
        # a<=0.2: encroachment (No inversion layer, so the entire profile is given by g and a=0
        #           (considered a,0.2 as in paper))
        # a<=2*b: inversion lapse rate is equal to or smaller than free lapse rate
        gprime = 0.
        raise Exception('No CI present')
    else:
        gprime = gravity * CIestimate['dth'] / T0
    # Brunt-Vaisala frequency
    N = np.sqrt(gravity * CIestimate['gamma'] / T0)
    # CI altitudes
    inv_bottom = CIestimate['h0']
    if inv_bottom <= H1 + 10.:  # H cannot be lower than H1 and the upper layer must be at least 10m
        raise Exception(f'CI too low, CI bottom located at z={int(inv_bottom)}m')
    H = CIestimate['h1']
    inv_top = CIestimate['h2']
    # Upper layer width
    H2 = H - H1
    # Compute height averaged quantities #
    # Layer 1
    U1 = height_average(us, zs, 0., H1)
    V1 = height_average(vs, zs, 0., H1)
    nu1 = height_average(nu, zs, 0., H1)
    # Layer 2
    U2 = height_average(us, zs, H1, H)
    V2 = height_average(vs, zs, H1, H)
    nu2 = height_average(nu, zs, H1, H)
    # Layer 3
    fu = interp1d(zs, us, fill_value='extrapolate')
    fv = interp1d(zs, vs, fill_value='extrapolate')
    z3 = np.linspace(H, 5000., 1000)
    U3 = np.trapz(fu(z3), z3) / (z3[-1] - z3[0])
    V3 = np.trapz(fv(z3), z3) / (z3[-1] - z3[0])
    # Upper atmosphere variable determination #
    u_strat = U3
    v_strat = V3
    h_strat = 10.e3
    N_strat = N
    # Fitting procedure
    if fit_trop:
        nLines = 2
        upper_atm = np.logical_and(inv_top < zs, zs < 15000.)
        altitudes = zs[upper_atm]
        theta_values = ths[upper_atm]
        u_values = us[upper_atm]
        v_values = vs[upper_atm]
        profile = np.stack((theta_values, altitudes), axis=1)
        # So have a np-array of dimensions n x 2 with profile1[:,0] the temperature and profile1[:,1] the heights
        profiles = [profile]
        results = line_optimization.smart_slsqp_optimization(nLines, profiles, verbose=0)
        theta_trop = results[0][0]
        a1 = results[0][1]
        dthetadz_trop = a1
        h_strat = results[0][2]
        a2 = results[0][3]
        theta_strat = theta_trop + a1 * h_strat
        u_strat = np.mean(u_values[altitudes > h_strat])
        v_strat = np.mean(v_values[altitudes > h_strat])
        dthetadz_strat = a1 + a2
        N_trop = np.sqrt(gravity * dthetadz_trop / theta_trop)
        N_strat = np.sqrt(gravity * dthetadz_strat / theta_strat)
        if plot:
            fig, (ax1, ax2, ax3) = plt.subplots(1, 3)
            ax1.plot(profile[:, 0], profile[:, 1] / 1.e3)
            ax1.plot(line_optimization.theta_approximation(nLines, profile[:, 1], results[0]), profile[:, 1] / 1.e3,
                     '--k')
            ax1.set_ylabel('$z$ [km]')
            ax1.set_xlabel('$\\theta$ [K]')
            ax2.plot(u_values, altitudes / 1.e3)
            ax2.plot([u_strat, u_strat], [h_strat / 1.e3, 15000. / 1.e3], '--k')
            ax2.set_xlabel('$u$ [m/s]')
            ax3.plot(v_values, altitudes / 1.e3)
            ax3.plot([v_strat, v_strat], [h_strat / 1.e3, 15000. / 1.e3], '--k')
            ax3.set_xlabel('$v$ [m/s]')
            plt.show()
        # Layer 3
        fu = interp1d(zs, us, fill_value='extrapolate')
        fv = interp1d(zs, vs, fill_value='extrapolate')
        z = np.linspace(H, h_strat, 1000)
        U3 = np.trapz(fu(z), z) / (h_strat - H)
        V3 = np.trapz(fv(z), z) / (h_strat - H)
        N = N_trop  # Use Brunt-Vaisala frequency as determined by the tropopauze fit
    # Set up ABL object
    return ABL(zs, us, vs, ths, tauxs, tauys,
                 H1, H2,
                 gprime, N, U3, V3,
                 fc,
                 nus=nu,
                 rho=1.225, TI=TI_h, z0=z0, ust=utau,
                 inv_bottom=inv_bottom, inv_top=inv_top,
                 h_strat=h_strat, Uinf=u_strat, Vinf=v_strat, Ninf=N_strat)


# NetCDF file location
file = r"C:\Users\u0138255\Documents\5 Other codes\FLOW API\FLOW_API\examples\cases\windio_4turbines_profiles_stable\wind_energy_system\FLOW_toy_study_wind_energy_system.yaml"

run_wayve(file)

