from foxes_api import runFoxes
from pywake_api import run_pywake
import os
import yaml
import sys

import sys
import windIO
from windIO.utils.yml_utils import validate_yaml, Loader, load_yaml
sys.path.append(windIO.__path__[0])



# parse input from command line
if len(sys.argv) != 2:
    raise Exception('Please supply the relevant yaml file as an input')

yaml_input = sys.argv[1]

# validate input
validate_yaml(yaml_input, windIO.__path__[0] + '/plant/wind_energy_system.yaml')


# get number of turbines
yaml_dat = load_yaml(yaml_input)
x = yaml_dat['wind_farm']['layouts']['initial_layout']['coordinates']['x']

yaml_input_no_ext = os.path.splitext(yaml_input)[0]  # Remove the file extension
output_dir_name = 'output_' + yaml_input_no_ext.replace(os.sep, '_')  # Replace directory separators
if not os.path.exists(output_dir_name):
    os.makedirs(output_dir_name)

# Specify input metadata file name
if 'name' in yaml_dat['attributes']['outputs']:
   output_dir_name = yaml_dat['attributes']['analyses']['outputs']['name']
   if not os.path.exists(output_dir_name):
      os.makedirs(output_dir_name)

output_file = output_dir_name + os.sep + 'recorded_inputs.yaml'

# Write yaml_dat to a YAML file
with open(output_file, 'w') as outfile:
    yaml.safe_dump(yaml_dat, outfile, default_flow_style=False)

# compute AEP (next step is to return a richer set of outputs)
pywake_aep = run_pywake(yaml_input, output_dir=output_dir_name)
#foxes_aep = runFoxes(yaml_input)

# report results
print('++++++')
print('PyWake AEP: ', pywake_aep, '(capacity factor %f)' % (pywake_aep / (len(x) * 15 * 8760 / 1e3)))
#print('FOXES AEP: ', foxes_aep, '(capacity factor %f)' % (foxes_aep / (len(x) * 15 * 8760 / 1e3)))
