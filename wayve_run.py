# General packages
import numpy as np
import netCDF4
from scipy.interpolate import interp1d
from windIO.utils.yml_utils import load_yaml
import matplotlib.pyplot as plt

# Atmospheric state setup
from wayve.abl.abl import ABL
from wayve.abl import line_optimization, ci_methods
from wayve.abl.abl_tools import height_average

# General APM setup
from wayve.three_layer_model import TLM
from wayve.grid.grid import Stat2Dgrid
from wayve.forcing.wind_farms.wake_models import coupling_methods, wake_models
from wayve.forcing.wind_farms.wind_farm import WF
from wayve.momentum_flux_parametrizations import FrictionCoefficients
from wayve.pressure.gravity_waves.gravity_waves import Uniform
from wayve.solvers import FixedPointIteration


def NetCDF(dat, H1=240., dh_max=None, plot=False):
    '''
    Method to set up an ABL object based on a NetCDF file

    Parameters
    ----------
    nc_file: str
        Path to file containing the atmospheric data
    H1 (optional): float
        Height of the wind-farm layer (default: 240.)
    dh_max (optional): float
        Maximum depth of the inversion layer used in the inversion curve fitting procedure (default: None)
    plot (optional): boolean
        Whether or not a plot of the upper atmosphere fitting procedure is shown (default: False)
    '''
    # Constants #
    gravity = 9.80665  # [m s-2]
    kappa = 0.41  # Von Karman constant
    omega = 7.2921159e-5  # angular speed of the Earth [rad/s]
    # Read out data #
    # Load file
    # Get cell-centered arrays
    zs = np.array(dat['z'])
    ths = np.array(dat['ths'])
    us = np.array(dat['speed']) * np.cos(np.deg2rad(dat['veer']))
    vs = np.array(dat['speed']) * np.sin(np.deg2rad(dat['veer']))
    # Get staggered arrays
    zst = np.array(dat['z'])
    tauxst = np.array(dat['tauxst'])
    tauyst = np.array(dat['tauyst'])
    # Total velocity and momentum flux #
    Ms = np.sqrt(np.power(us, 2) + np.power(vs, 2))
    taust = np.sqrt(np.power(tauxst, 2) + np.power(tauyst, 2))
    # Surface layer #
    # Friction velocity
    utau = np.sqrt(tauxst[0]**2 + tauyst[0]**2)     # Assume friction velocity is not given explicitly
    if 'utau' in dat.keys():
        utau = dat['utau']
    # Monin-Obukhov theory (currently used for RZ fit, default behaviour / additional inputs should be discussed)
    wth = 0.        # Assume no surface heat flux
    if 'wth' in dat.keys():
        wth = dat['wth']
    T0 = 288.15     # Assume near-surface temperature of 15C
    if 'T0' in dat.keys():
        T0 = dat['T0']
    zeta = -2.0 * kappa * gravity * wth / (T0 * utau ** 3)
    # Estimate boundary layer height based on momentum flux #
    f_tau = interp1d(taust, zst)
    blh = f_tau(0.01 * utau)
    if 'blh' in dat.keys():
        blh = dat['blh']
    # Estimate eddy viscosity #
    f_s = interp1d(zs, Ms, bounds_error=False, fill_value=(0., Ms[-1]))
    ds = np.gradient(f_s(zst), zst, edge_order=2)
    nu = np.zeros(zst.shape)
    nu[zst<blh] = np.abs(taust[zst<blh] / ds[zst<blh])
    nu[zst<blh] = np.abs(taust[zst<blh] / (1e-9 + ds[zst<blh])) # TODO: CHECK
    # Turbulent intensity at hub height #
    if zeta > 0.0:  # stable
        # From Nieuwstadt (1984): q/sqrt(tau) = 3
        tke = 4.5 * taust
    else:  # unstable
        # From Stull (1988): q^2/tau = 8.5+2.5
        tke = 5.5 * taust
    TIs = np.sqrt(2. / 3. * tke) / Ms
    fTi = interp1d(zs, TIs)
    TI = fTi(H1 / 2.0).item()
    if 'TI' in dat.keys():
        TI = dat['TI']
    # Surface roughness (ignoring stability effects) #
    z0 = zs[0] / np.exp(kappa * Ms[0] / utau)
    if 'z0' in dat.keys():
        z0 = dat['z0']
    # Coriolis parameter #
    phi = 0.377     # Assume latitude location
    if 'phi' in dat.keys():
        phi = dat['phi']
    fc = 2 * omega * np.sin(phi)
    # Estimate inversion parameters with RZ fit #
    # Relevant part of the vertical profiles
    zCI = zs[zs < 5000]
    thCI = ths[zs < 5000]
    # Stable or unstable profile
    if zeta > 0.02:
        # Ignore temperature decrease inside SBL
        # (we are trying to indentify the mixing layer that preceded this SBL)
        # We want to capture the mixed layer (or residual layer since we are in
        # SBL) that preceded the SBL. Therefore, we take the potential temperature
        # at the top of the ABL where we have the mixed layer and we extrapolate
        # till the bottom. We use this constant value in the ABL.
        # p0 are the initial guess for [a,b,thm,l,dh] used in Ramp&Zar model
        thCI[zCI < blh] = interp1d(zCI, thCI)(blh)
        CIestimate = ci_methods.RZfit(zCI, thCI, p0=[0.9, 0.1, T0, 1000., 100.0],
                                      dh_max=dh_max)
    else:
        # Ignore temperature increase in CBL surface layer, therefore we take
        # the lowest value of potential temperature. We are able to capture the
        # mixed layer in this way, where the temperature is constant and equal
        # to the lowest theta. We use this constant value in the ABL.
        # p0 are the initial guess for [a,b,thm,l,dh] used in Ramp&Zar model
        thCI[0:np.argmin(thCI)] = np.min(thCI)
        CIestimate = ci_methods.RZfit(zCI, thCI, p0=[0.9, 0.1, T0, blh, 100.0],
                                      dh_max=dh_max)
    # Inversion strength
    if CIestimate['a'] <= 0.2 or CIestimate['a'] <= 2 * CIestimate['b']:
        # No inversion strength in the following cases:
        # a<=0.2: encroachment (No inversion layer, so the entire profile is given by g and a=0
        #           (considered a,0.2 as in paper))
        # a<=2*b: inversion lapse rate is equal to or smaller than free lapse rate
        gprime = 0.
    else:
        gprime = gravity * CIestimate['dth'] / T0
    # CI altitudes
    inv_bottom = CIestimate['h1']
    H = np.max([CIestimate['h1'], H1 + 10.])  # H cannot be lower than H1 and the upper layer must be at least 10m
    inv_top = CIestimate['h2']
    # Upper layer width
    H2 = H - H1
    # Upper atmosphere variable determination #
    nLines = 2
    altitudes_check = zs[zs < 15000.]
    altitudes = altitudes_check[altitudes_check > H]
    theta_values = ths[zs < 15000.]
    theta_values = theta_values[altitudes_check > H]
    u_values = us[zs < 15000.]
    u_values = u_values[altitudes_check > H]
    v_values = vs[zs < 15000.]
    v_values = v_values[altitudes_check > H]
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
    # Compute height averaged quantities
    # Layer 1
    U1 = height_average(us, zs, 0., H1)
    V1 = height_average(vs, zs, 0., H1)
    nu1 = height_average(nu, zst, 0., H1)
    # Layer 2
    U2 = height_average(us, zs, H1, H)
    V2 = height_average(vs, zs, H1, H)
    nu2 = height_average(nu, zst, H1, H)
    # Layer 3
    fu = interp1d(zs, us, fill_value='extrapolate')
    fv = interp1d(zs, vs, fill_value='extrapolate')
    z = np.linspace(H, h_strat, 1000)
    U3 = np.trapz(fu(z), z) / (h_strat - H)
    V3 = np.trapz(fv(z), z) / (h_strat - H)
    N = N_trop  # Use Brunt-Vaisala frequency as determined by the tropopauze fit
    # Set up ABL object
    return ABL(H1, H2, U1, V1, U2, V2, U3, V3, gprime, N, fc, TI, nu1, nu2,
               zs, zst, us, vs, ths, tauxst, tauyst,
               z0=z0, inv_bottom=inv_bottom, inv_top=inv_top, h_strat=h_strat, Uinf=u_strat, Vinf=v_strat, Ninf=N_strat)


def wf_setup(center_location):
    """
    Set up an example hard-coded wind farm setup.

    Long-term, this should be replaced by an interface to the Wind-IO framework.

    Parameters
    ----------
    center_location: array-like
        Grid location of the center of the farm

    Returns
    -------
    wind_farm:  WF object
        Wayve wind farm representation
    """
    # Input parameters
    D = 154.            # Turbine diameter [m]
    Ntx = 18            # Number of turbine rows
    Nty = 27            # Number of turbine columns
    Ct = 0.8            # Turbine thrust coefficient
    zh = 120.           # Turbine hub height [m]
    Lwfx = Ntx/9.0e-4   # Wind-farm length
    Lwfy = Nty/9.0e-4   # Wind-farm width
    # Turbine coordinates
    xs = np.linspace(0, Lwfx, Ntx, endpoint=False)
    ys, dy = np.linspace(0., Lwfy, Nty, endpoint=False, retstep=True)
    X, Y = np.meshgrid(xs, ys, indexing='ij')
    # Staggered layout
    Y[1::2, :] += dy/2.
    # Reposition to be at grid center
    X_grid = np.ravel(X) - (xs[-1]-xs[0]) / 2. + center_location[0]
    Y_grid = np.ravel(Y) - (ys[-1]-ys[0]+dy/2.) / 2. + center_location[1]
    # Use Gaussian wake model of Niayifar (2015)
    wakemodel = wake_models.Niayifar()
    # Set up coupling object
    # Use velocity at 10D upstream as inflow velocity
    coupling = coupling_methods.Upstream(wakemodel, 10. * D)
    # Gaussian filter length
    Lfilter = 1000.
    # Turbine parameter arrays
    Cts = Ct*np.ones(X.size)
    Ds = D*np.ones(X.size)
    zhs = zh*np.ones(X.size)
    # Generate wind farm object
    wind_farm = WF(X_grid, Y_grid, Ds, Cts, zhs, Lfilter, coupling)
    return wind_farm


def wayve_run(nc_file, verbose=False):
    """
    Run wayve for an example hard-coded wind farm setup, with atmospheric conditions defined by the given .nc file.

    Parameters
    ----------
    nc_file: str
        File location
    verbose (optional): bool
        flag for printing solver results (default: False)
    Returns
    -------
    turbine_power: numpy array
        List of turbine power outputs, scaled with air density
    """
    # ---------------------------------------------------- #
    # ------------- Step 1: define abl ------------------- #
    # ---------------------------------------------------- #

    # Set up abl
    abl = NetCDF(nc_file, plot=False)   # Set plot to True to show the free atmosphere fitting procedure.

    # ---------------------------------------------------- #
    # -------- Step 2: define numerical grid ------------- #
    # ---------------------------------------------------- #

    # Numerical parameters
    Nx = 2000           # grid points in x-direction
    Lx = 1.e6           # grid size in x-direction [m]
    Ny = 800            # grid points in y-direction
    Ly = 0.4e6          # grid size in y-direction [m]

    # Generate 2D grid object
    grid = Stat2Dgrid(Lx, Nx, Ly, Ny)

    # ---------------------------------------------------- #
    # ------------ Step 3: define wind farm -------------- #
    # ---------------------------------------------------- #

    # Place wind farm at center of grid
    center_location = [Lx/2., Ly/2.]

    # Wind farm setup
    wind_farm = wf_setup(center_location)

    # ------------------------------------------------------ #
    # - Step 4: create TLM model from components and solve - #
    # ------------------------------------------------------ #

    # Momentum flux parametrization
    mfp = FrictionCoefficients()

    # pressure feedback parametrization
    pressure = Uniform()

    # Create static 2D model
    model = TLM(grid, wind_farm, abl, mfp, pressure)

    # Use fixed-point iteration solver
    solver = FixedPointIteration()

    # Solve linear system
    _ = model.solve(method=solver, verbose=verbose)

    # Read out power outputs
    turbine_power = wind_farm.P_turbines()  # [W/(kg/m3)]

    return turbine_power


# NetCDF file location
file = "examples/flow_example_timeseries.yaml"
dat = load_yaml(file)
if 'timeseries' not in dat['site']['energy_resource']['wind_resource']:
   print('FLOW API only supports WAYVE timeseries')
   quit()

# Run wayve for each timeseries entry
powers = []
for tt in range(len(dat['site']['energy_resource']['wind_resource']['timeseries'])):
   powers.append(wayve_run(dat['site']['energy_resource']['wind_resource']['timeseries'][tt], verbose=True))
   print('powers: ',  powers)

# Print output
print('powers: ', np.array(powers))
#print(f"Turbine power [W/(kg/m3)]: {power}")
