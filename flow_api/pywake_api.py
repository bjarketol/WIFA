import numpy as np
from py_wake import HorizontalGrid
from collections import OrderedDict
import warnings
import pandas as pd
import matplotlib.pyplot as plt
import xarray as xr
import os
import yaml
from scipy.interpolate import interp1d
from py_wake.deficit_models.fuga import FugaDeficit
from py_wake.site import XRSite
from py_wake.wind_turbines import WindTurbine
from py_wake.wind_farm_models import PropagateDownwind, All2AllIterative
from py_wake.examples.data.hornsrev1 import Hornsrev1Site
from py_wake.wind_turbines.power_ct_functions import PowerCtTabular
from py_wake.deflection_models import JimenezWakeDeflection
from py_wake.deficit_models.noj import NOJLocalDeficit
from py_wake.deficit_models.gaussian import BastankhahGaussianDeficit
from py_wake.turbulence_models import STF2005TurbulenceModel, STF2017TurbulenceModel, CrespoHernandez
from py_wake import NOJ, BastankhahGaussian
from py_wake.deficit_models import SelfSimilarityDeficit2020
from py_wake.superposition_models import LinearSum, SquaredSum
from py_wake.rotor_avg_models import RotorCenter, GridRotorAvg, EqGridRotorAvg, GQGridRotorAvg, CGIRotorAvg, PolarGridRotorAvg, polar_gauss_quadrature, GaussianOverlapAvgModel
from windIO.utils.yml_utils import validate_yaml, Loader, load_yaml


# Define default values for wind_deficit_model parameters
DEFAULTS = {
    'wind_deficit_model': {
        'name': 'Jensen',
        'k': 0.04,  # Default wake expansion coefficient for Jensen
        'k2': 0.0,  # TI wake expansion modifier
    },
    'deflection_model': {
        'name': 'Jimenez',
        'beta': 0.1,  # Default Jimenez deflection coefficient
    },
    'turbulence_model': {
        'name': 'STF2005',
        'c1': 1.0,  # Default STF C1 value
        'c2': 1.0,  # Default STF C2 value
    },
    'superposition_model': {
        'name': 'Linear',
    },
    'rotor_averaging': {
        'name': 'Center',
    },
    'blockage': {
        'name': None,
        'ss_alpha': 0.8888888888888888
    }
}

def get_with_default(data, key, defaults):
    """
    Retrieve a value from a dictionary, using a default if the key is not present.
    If the value is a dictionary, apply the same process recursively.
    """
    if key not in data:
        return defaults[key]
    elif isinstance(data[key], dict):
        # For nested dictionaries, ensure all subkeys are checked for defaults
        return {sub_key: get_with_default(data[key], sub_key, defaults[key]) for sub_key in defaults[key]}
    else:
        return data[key]


def weighted_quantile(values, quantiles, sample_weight=None, 
                      values_sorted=False, old_style=False):
    """ Very close to numpy.percentile, but supports weights.
    NOTE: quantiles should be in [0, 1]!
    :param values: numpy.array with data
    :param quantiles: array-like with many quantiles needed
    :param sample_weight: array-like of the same length as `array`
    :param values_sorted: bool, if True, then will avoid sorting of
        initial array
    :param old_style: if True, will correct output to be consistent
        with numpy.percentile.
    :return: numpy.array with computed quantiles.
    """
    values = np.array(values)
    quantiles = np.array(quantiles)
    if sample_weight is None:
        sample_weight = np.ones(len(values))
    sample_weight = np.array(sample_weight)
    assert np.all(quantiles >= 0) and np.all(quantiles <= 1), \
        'quantiles should be in [0, 1]'

    if not values_sorted:
        sorter = np.argsort(values)
        values = values[sorter]
        sample_weight = sample_weight[sorter]

    weighted_quantiles = np.cumsum(sample_weight) - 0.5 * sample_weight
    if old_style:
        # To be convenient with numpy.percentile
        weighted_quantiles -= weighted_quantiles[0]
        weighted_quantiles /= weighted_quantiles[-1]
    else:
        weighted_quantiles /= np.sum(sample_weight)
    return np.interp(quantiles, weighted_quantiles, values)

def run_pywake(yamlFile, output_dir='output'):

    system_dat = load_yaml(yamlFile)

    # define turbine
    farm_dat = system_dat['wind_farm']
    hh = farm_dat['turbines']['hub_height']
    rd = farm_dat['turbines']['rotor_diameter']
    if 'Cp_curve' in farm_dat['turbines']['performance']:
       cp = farm_dat['turbines']['performance']['Cp_curve']['Cp_values']
       cp_ws = farm_dat['turbines']['performance']['Cp_curve']['Cp_wind_speeds']
       power_curve_type = 'cp'
    elif 'power_curve' in farm_dat['turbines']['performance']:
       cp_ws = farm_dat['turbines']['performance']['power_curve']['power_wind_speeds']
       pows = farm_dat['turbines']['performance']['power_curve']['power_values']
       power_curve_type = 'power'
    else: raise Exception('Bad Power Curve')
    ct = farm_dat['turbines']['performance']['Ct_curve']['Ct_values']
    ct_ws = farm_dat['turbines']['performance']['Ct_curve']['Ct_wind_speeds']
    speeds = np.linspace(np.min([cp_ws, ct_ws]), np.max([cp_ws, ct_ws]), 10000)
    cts_int = np.interp(speeds, ct_ws, ct)
    if power_curve_type == 'power':
       powers = np.interp(speeds, cp_ws, pows)
    else:
       cps_int = np.interp(speeds, cp_ws, cp)
       powers =  0.5 * cps_int * speeds ** 3 * 1.225 * (rd / 2) ** 2 * np.pi
    turbine = WindTurbine(name=farm_dat['turbines']['name'], diameter=rd, hub_height=hh, 
                          powerCtFunction=PowerCtTabular(speeds, powers, power_unit='W', ct=cts_int))


   #farm = 'examples/plant/plant_wind_farm/IEA37_case_study_3_wind_farm.yaml'
   #with open(farm, "r") as stream:
   #    try:
   #        farm_dat = yaml.safe_load(stream)
   #    except yaml.YAMLError as exc:
   #        print(exc)

   #resource = 'examples/plant/plant_energy_resource/IEA37_case_study_4_energy_resource.yaml'
   #with open(resource, "r") as stream:
    #    try:
    #        resource_dat = yaml.safe_load(stream)
    #    except yaml.YAMLError as exc:
    #        print(exc)
    resource_dat = system_dat['site']['energy_resource']
    WFXLB = np.min(system_dat['site']['boundaries']['polygons'][0]['x'])
    WFXUB = np.max(system_dat['site']['boundaries']['polygons'][0]['x'])
    WFYLB = np.min(system_dat['site']['boundaries']['polygons'][1]['y'])
    WFYUB = np.max(system_dat['site']['boundaries']['polygons'][1]['y'])
    
    ##################
    # construct site
    ##################
    if 'time' in resource_dat['wind_resource'].keys():
       timeseries = True
       wind_resource_timeseries = resource_dat['wind_resource']['time']
       times = wind_resource_timeseries
       if 'height' in resource_dat['wind_resource'].keys():
           heights = resource_dat['wind_resource']['height']
       else:
           heights = None
       ws = resource_dat['wind_resource']['wind_speed']['data']
       wd = resource_dat['wind_resource']['wind_direction']['data']
       print(np.array(ws).shape, np.array(heights).shape)
       if heights:
          ws = interp1d(heights, ws, axis=1)(hh)
          wd = interp1d(heights, wd, axis=1)(hh)
       assert(len(times) == len(ws))
       assert(len(wd) == len(ws))
       site = Hornsrev1Site()
       if 'turbulence_intensity' not in resource_dat['wind_resource']:
          TI = 0.02
       else:
          TI =  resource_dat['wind_resource']['turbulence_intensity']['data']
          if heights: TI = interp1d(heights, TI, axis=1)(hh)
       #ite = XRSite(xr.Dataset(
       # data_vars={'P': (('time'), np.ones(len(ws)) / len(speeds)), },
       # coords={'time': range(len(times)),
       #         'ws': speeds,
       #         'wd': wd}))
    
    elif 'weibull_k' in resource_dat['wind_resource'].keys():
       A = resource_dat['wind_resource']['weibull_a']
       k = resource_dat['wind_resource']['weibull_k']
       freq = resource_dat['wind_resource']['sector_probability']
       wd = resource_dat['wind_resource']['wind_direction']
       if 'wind_speed' in resource_dat['wind_resource']:
           ws = resource_dat['wind_resource']['wind_speed']
       else:
           ws = np.arange(2, 30, 1)
       site = XRSite(
              ds=xr.Dataset(data_vars=
                               {'Sector_frequency': ('wd', freq['data']), 
                                'Weibull_A': ('wd', A['data']), 
                                'Weibull_k': ('wd', k['data']), 
                                'TI': (resource_dat['wind_resource']['turbulence_intensity']['dims'], resource_dat['wind_resource']['turbulence_intensity']['data'])
                                },
                             coords={'wd': wd, 'ws': ws}))
       
       timeseries = False
       TI =  resource_dat['wind_resource']['turbulence_intensity']['data']
    else:
       timeseries = False
       ws = resource_dat['wind_resource']['wind_speed']
       wd = resource_dat['wind_resource']['wind_direction']
       P = np.array(resource_dat['wind_resource']['probability']['data'])
       site = XRSite(ds=xr.Dataset(data_vars={'P': (['wd', 'ws'], P)}, coords = {'ws': ws, 'wd': wd, 'TI': resource_dat['wind_resource']['turbulence_intensity']['data']}))
       TI = resource_dat['wind_resource']['turbulence_intensity']['data']

    #if 'name' in system_dat['attributes']['analysis']['outputs']:
    #   output_dir = system_dat['attributes']['analysis']['outputs']['name']
    #   if not os.path.exists(output_dir):
    #      os.makedirs(output_dir)




    # get x and y positions
    x = farm_dat['layouts']['initial_layout']['coordinates']['x']
    y = farm_dat['layouts']['initial_layout']['coordinates']['y']
    

    wind_deficit_model_data = get_with_default(system_dat['attributes']['analysis'], 'wind_deficit_model', DEFAULTS)

    deficit_args = {}
    deficit_param_mapping = {}
    if wind_deficit_model_data['name'] == 'Jensen':
       wakeModel = NOJLocalDeficit
       deficit_param_mapping = {'k': 'k', 'k2': 'k2'}
    elif wind_deficit_model_data['name'] == 'Bastankhah2014':
       wakeModel = BastankhahGaussianDeficit
       deficit_param_mapping = {'k': 'k', 'ceps': 'ceps'}
       #from py_wake.deficit_models.utils import ct2a_mom1d
       #deficit_args['ct2a'] = ct2a_mom1d
    elif wind_deficit_model_data['name'].upper() == 'FUGA':
       wakeModel = FugaDeficit
       from pyfuga import get_luts
       lut = get_luts(folder = 'luts', # Path where all files (intermediate and final) are stored
                                           zeta0=0, # Stability parameter
                                           nkz0=8, 
                                           nbeta=32, 
                                           diameter=rd, 
                                           zhub=hh, 
                                           z0=0.00001, 
                                           zi=500, 
                                           zlow=70, 
                                           zhigh=70, 
                                           lut_vars=['UL'], 
                                           nx=2048, 
                                           ny=512, n_cpu=1) 
       deficit_args['LUT_path'] = r'luts/LUTs_Zeta0=0.00e+00_8_32_D%.1f_zhub%.1f_zi500_z0=0.00001000_z69.2-72.8_UL_nx2048_ny512_dx44.575_dy11.14375.nc' % (rd, hh)
       #deficit_args['LUT_path'] = 'luts/LUTs_Zeta0=0.00e+00_8_32_D%i_zhub%i_zi500_z0=0.00001000_z70.0_UL_nx2048_ny512_dx20.0_dy5.0.nc' % (rd, hh)
    else:
       raise Exception('%s wake model not implemented in PyWake' % wind_deficit_model_data['name'])
    for key in wind_deficit_model_data.keys():
        if key == 'name': continue
        if key in deficit_param_mapping.keys():
            deficit_args[deficit_param_mapping[key]] = wind_deficit_model_data[key]

    if 'k2' in deficit_args:
        k = deficit_args.pop('k')
        k2 = deficit_args.pop('k2')
        deficit_args['a'] = [k, k2]

    # Continuing from the previous example...
    
    deflection_model_data = get_with_default(system_dat['attributes']['analysis'], 'deflection_model', DEFAULTS)
    turbulence_model_data = get_with_default(system_dat['attributes']['analysis'], 'turbulence_model', DEFAULTS)
    superposition_model_data = get_with_default(system_dat['attributes']['analysis'], 'superposition_model', DEFAULTS)
    rotor_averaging_data = get_with_default(system_dat['attributes']['analysis'], 'rotor_averaging', DEFAULTS)
    blockage_data = get_with_default(system_dat['attributes']['analysis'], 'blockage', DEFAULTS)
    
    # Map the deflection model
    if deflection_model_data['name'].lower() == 'none':
        deflectionModel = None
    elif deflection_model_data['name'].lower() == 'jimenez':
        deflectionModel = JimenezWakeDeflection(beta=deflection_model_data['beta'])  # Assuming Jimenez takes 'beta' as an argument
    elif deflection_model_data['name'] == "None": deflectionModel = None
    else:
        raise Exception('%s deflection model not implemented' % deflection_model_data['name'])
    deficit_args['use_effective_ws'] = True
    
    # Map the turbulence model
    if turbulence_model_data['name'].lower() == 'none':
        turbulenceModel = None
    elif turbulence_model_data['name'].upper() == 'STF2005':
        turbulenceModel = STF2005TurbulenceModel(c=[turbulence_model_data['c1'], turbulence_model_data['c2']])
    elif turbulence_model_data['name'].upper()  == 'STF2017':
        turbulenceModel = STF2017TurbulenceModel(c=[turbulence_model_data['c1'], turbulence_model_data['c2']])
        #turbulenceModel = STF(turbulence_model_data['c1'], turbulence_model_data['c2'])  # Assuming STF takes 'c1' and 'c2' as arguments
    elif turbulence_model_data['name'].upper()  == 'CRESPOHERNANDEZ':
        turbulenceModel = CrespoHernandez()
    else:
        raise Exception('%s turbulence model not implemented' % turbulence_model_data['name'])
    
    # Map the superposition model
    if superposition_model_data['name'].lower() == 'linear':
        superpositionModel = LinearSum()
    elif superposition_model_data['name'].lower() == 'squared':
        superpositionModel = SquaredSum()
    else:
        raise Exception('%s superposition model not implemented' % superposition_model_data['name'])
    
    # Map the rotor averaging model
    if rotor_averaging_data['name'].lower() == 'center':
        rotorAveraging = RotorCenter()
    elif rotor_averaging_data['name'].lower() == 'avg_deficit':
        rotorAveraging = GridRotorAvg()
    else:
        raise Exception('%s rotor averaging model not implemented' % rotor_averaging_data['name'])
    
    if blockage_data['name'] == 'None' or  blockage_data['name'] is None:
       blockage = None
    elif blockage_data['name'] == 'SelfSimilarityDeficit2020':
       blockage = SelfSimilarityDeficit2020(ss_alpha=blockage_data['ss_alpha'])
    elif blockage_data['name'].upper() == 'FUGA':
       blockage =  FugaDeficit(deficit_args['LUT_path'])
    else: raise Exception('Bad Blockage Specified')


    solver_args = {}
    if blockage is not None:
       solver = All2AllIterative
       solver_args['blockage_deficitModel'] = blockage
    else: 
       solver = PropagateDownwind

    windFarmModel = solver(
                           site,
                           turbine,
                           wake_deficitModel=wakeModel(
                              rotorAvgModel=rotorAveraging,
                              groundModel=None, **deficit_args),
                           superpositionModel=superpositionModel,
                           deflectionModel=deflectionModel,
                           turbulenceModel=turbulenceModel,
                           **solver_args
                           )
    #noj = NOJ(site, turbine, turbulenceModel=None)
#    sim_res = noj(x, y)
    sim_res = windFarmModel(x, y, time=timeseries, ws=ws, wd=wd, TI=TI, yaw=0, tilt=0)
    aep = sim_res.aep(normalize_probabilities=not timeseries).sum()
    print('aep is ', aep, 'GWh')
    #print('aep is ', sim_res.aep().sum(), 'GWh')
    print('(%.2f capcacity factor)' % ( aep / (len(x) * turbine.power(10000) * 8760 / 1e9)))


    ######################
    # Construct Outputs
    ####################
    # turbine specific AEP
    if timeseries:
       aep_per_turbine = sim_res.aep(normalize_probabilities=True).sum(['time']).to_numpy()
    else:
       aep_per_turbine = sim_res.aep(normalize_probabilities=True).sum(['ws', 'wd']).to_numpy()

    # names -- to be updated according to team consensus
    # (also, TODO: read these in from the input WES file!)
    name = "FLOW tool output"
    wind_deficit_model = "Bastankhah’s Gaussian wake model"
    statistical_description = 'PDF' # "time_series"
    statistical_dimensions = ['wind_direction','wind_velocity']
    
    # flow field NetCDF file
    # -----------------
    wind_output_file = "FarmFlow.nc"
    wind_output_variables = ['velocity_u', 'turbulence_intensity', 'turbulence_k', 'pressure', 'wind_loss_to_inflow']
    
    power_output_file = "FarmPower.nc"
    power_output_variables = ['power_per_turbine', 'load_per_turbine', 'Remaining_Useful_life_per_turbine', 'power_loss_to_inflow_per_turbine']

    # Create a dictionary to represent the YAML file structure
    # Create an OrderedDict to represent the YAML file structure
    data = OrderedDict([
      ('name', name),
      ('FLOW_simulation_config', 
        dict([
          ('tool', 'PyWake'),
          ('wind_deficit_model', wind_deficit_model),
          #('wind_energy_system', '!include recorded_inputs.yaml')
          ('wind_energy_system', 'INCLUDE_YAML_PLACEHOLDER')
          ])),
      ('FLOW_simulation_outputs', dict([
        ('statistical_description', statistical_description),
        #('statistical_dimensions', statistical_dimensions),
        #('power_percentiles', power_percentiles),
        #('AEP', aep),
        #('AEP_per_turbine', aep_per_turbine),
        #('lifetime_per_turbine', lifetime_per_turbine),
        #('wind_output_file', wind_output_file),
        #('wind_output_variables', wind_output_variables),
        #('power_output_file', power_output_file),
        #('power_output_variables', power_output_variables)
      ]))
      ])

    if 'AEP' in system_dat['attributes']['outputs']:
       data['FLOW_simulation_outputs']['AEP'] = float(aep)
    
    #if 'power_percentiles' in system_dat['attributes']['analysis']['outputs']:
    # if system_dat['attributes']['analysis']['outputs']['power_percentiles']['report']:
    #   # compute power percentiles
    #   percentiles = np.array(system_dat['attributes']['analysis']['outputs']['power_percentiles']['percentiles']) / 100
    #   if not timeseries:

    #      # Flatten the power and P arrays
    #      power_values = sim_res.Power.sum('wt').values.flatten()
    #      probabilities = sim_res.P.values.flatten()

    #      # Compute the weighted percentiles
    #      power_percentiles = weighted_quantile(power_values, percentiles, sample_weight=probabilities)

    #   else:
    #       total_power = sim_res.Power.sum(dim='wt')
    #       power_percentiles = total_power.quantile(percentiles, dim=['time']).to_numpy()
    #   data['FLOW_simulation_outputs']['computed_percentiles'] = system_dat['attributes']['analysis']['outputs']['power_percentiles']['percentiles']
    #   data['FLOW_simulation_outputs']['power_percentiles'] = power_percentiles
    
    if 'turbine_outputs' in system_dat['attributes']['outputs']:
       #print('aep per turbine', list(aep_per_turbine)); hey
       #data['FLOW_simulation_outputs']['AEP_per_turbine'] = [float(value) for value in aep_per_turbine]
       sim_res_formatted = sim_res[['Power', 'WS_eff']].rename({'Power': 'power', 'WS_eff': 'effective_wind_speed', 'wt': 'turbine'})
       sim_res_formatted['power'] /= 1e3 # Watts to kW
       sim_res_formatted.to_netcdf(output_dir + os.sep + 'PowerTable.nc')

    print(sim_res)

    # flow field handling
    flow_map = None
    if 'flow_field' in system_dat['attributes']['outputs'] and not timeseries:

       # compute flow map for specified directions (wd) and speeds (ws)
       flow_map = sim_res.flow_map(
                            x = np.linspace(WFXLB, WFXUB, 100),
                            y = np.linspace(WFYLB, WFYUB, 100),
                            wd=system_dat['attributes']['analysis']['outputs']['flow_field']['directions'],
                            ws=system_dat['attributes']['analysis']['outputs']['flow_field']['speeds'],)

       # remove unwanted data
       flow_map = flow_map.drop_vars(['WD', 'WS', 'TI', 'P'])

       # raise warning if user requests data we can not provide
       if any(element not in ['velocity_u', 'turbulence_intensity'] for element in system_dat['attributes']['analysis']['outputs']['flow_field']['output_variables']):
          warnings.warn('PyWake can only output velocity_u and turbulence_intensity')

       # remove TI or WS if they are not requested
       if 'turbulence_intensity' not in system_dat['attributes']['analysis']['outputs']['flow_field']['output_variables']:
          flow_map = flow_map.drop_vars(['TI_eff'], inplace=True)
       if 'velocity_u' not in system_dat['attributes']['analysis']['outputs']['flow_field']['output_variables']:
          flow_map = flow_map.drop_vars(['WS_eff'], inplace=True)

    elif 'flow_field' in system_dat['attributes']['outputs'] and timeseries:
       #time_to_plot = system_dat['attributes']['outputs']['flow_field']['time'][0]
       #print('TIMES ', time_to_plot)
       flow_map = sim_res.flow_map(HorizontalGrid(x = np.linspace(WFXLB, WFXUB, 100),
y = np.linspace(WFYLB, WFYUB, 100)),
                            wd=wd,
                            ws=ws)
     # power table    
    #if 'power_table' in system_dat['attributes']['analysis']['outputs']:
    #   # todo: more in depth stuff here, include loads
    #   #data['FLOW_simulation_outputs']['power_output_variables'] = tuple(system_dat['attributes']['analysis']['outputs']['flow_field']['output_variables'])
#
    if flow_map:
       # save data
       flow_map = flow_map[['WS_eff', 'TI_eff']].drop(['wd', 'ws']).rename({'h': 'z', 'WS_eff': 'effective_wind_speed', 'TI_eff': 'turbulence_intensity'})
       flow_map.to_netcdf(output_dir + os.sep + 'FarmFlow.nc')

       # record data
       data['FLOW_simulation_outputs']['wind_output_file'] = 'FarmFlow.nc'
       data['FLOW_simulation_outputs']['wind_output_variables'] = system_dat['attributes']['outputs']['flow_field']['output_variables']

    # Write out the YAML data
    output_yaml_nam = output_dir + os.sep + 'output.yaml'
#    with open(output_yaml_nam, 'w') as file:
#        yaml.dump(data, file, default_flow_style=False, allow_unicode=True)
#
#    with open(output_yaml_nam, 'r') as f:
#       yaml_content = f.read()
#
#    yaml_content = yaml_content.replace('INCLUDE_YAML_PLACEHOLDER', '!include recorded_inputs.yaml')
#
#    # Save the post-processed YAML
#    with open(output_yaml_nam, 'w') as f:
#       f.write(yaml_content)

######################
# Construct Outputs
####################
    # Create a dictionary to represent the simplified YAML file structure
    data = {
        'wind_energy_system': 'INCLUDE_YAML_PLACEHOLDER',
        'power_table': 'INCLUDE_POWER_TABLE_PLACEHOLDER',
        'flow_field': 'INCLUDE_FLOW_FIELD_PLACEHOLDER'
    }
    
    # Write out the YAML data to the specified output directory
    output_yaml_name = os.path.join(output_dir, 'output.yaml')
    with open(output_yaml_name, 'w') as file:
        yaml_content = yaml.dump(data, file, default_flow_style=False, allow_unicode=True)
    
    # Replace placeholders with actual include directives, ensuring no quotes are used
    with open(output_yaml_name, 'r') as file:
        yaml_content = file.read()

    #yaml_content = yaml_content.replace('INCLUDE_YAML_PLACEHOLDER', '!include recorded_inputs.yaml')
    #
    ## Save the post-processed YAML
    #with open(output_yaml_nam, 'w') as f:
    #   f.write(yaml_content)
    data = {
        'wind_energy_system': 'INCLUDE_YAML_PLACEHOLDER',
        'power_table': 'INCLUDE_POWER_TABLE_PLACEHOLDER',
        'flow_field': 'INCLUDE_FLOW_FIELD_PLACEHOLDER'
    }

    # Write out the YAML data to the specified output directory
    output_yaml_name = os.path.join(output_dir, 'output.yaml')
    with open(output_yaml_name, 'w') as file:
        yaml_content = yaml.dump(data, file, default_flow_style=False, allow_unicode=True)
    
    # Replace placeholders with actual include directives, ensuring no quotes are used
    with open(output_yaml_name, 'r') as file:
        yaml_content = file.read()
    
    yaml_content = yaml_content.replace('INCLUDE_YAML_PLACEHOLDER', '!include recorded_inputs.yaml')
    yaml_content = yaml_content.replace('INCLUDE_POWER_TABLE_PLACEHOLDER', '!include PowerTable.nc')
    yaml_content = yaml_content.replace('INCLUDE_FLOW_FIELD_PLACEHOLDER', '!include FarmFlow.nc')
    
    # Save the post-processed YAML
    with open(output_yaml_name, 'w') as file:
        file.write(yaml_content)

    return aep
