from .foxes_api import run_foxes
from .pywake_api import run_pywake
from .cs_api.cs_modules.csLaunch.cs_run_function import run_code_saturne
from .wayve_api import run_wayve
import os
import yaml
import sys
import argparse

import windIO
from windIO.utils.yml_utils import validate_yaml, load_yaml
sys.path.append(windIO.__path__[0])

def run_api(yaml_input):

    # validate input
    validate_yaml(yaml_input, windIO.__path__[0] + '/plant/wind_energy_system.yaml')

    # get number of turbines
    yaml_dat = load_yaml(yaml_input)

    model_name = yaml_dat['attributes']['flow_model']['name']

    if model_name.lower() == 'pywake':
        pywake_aep = run_pywake(yaml_input)

    elif model_name.lower() == 'foxes':
        foxes_aep = run_foxes(yaml_input)

    elif model_name.lower() == 'wayve':
        run_wayve(yaml_input)

    elif model_name.lower() == 'code_saturne':
        run_code_saturne(yaml_input, test_mode=True)

    else:
        print('Invalid Model')

def run():

    parser = argparse.ArgumentParser()
    parser.add_argument("input_yaml", help="The input yaml file")
    args = parser.parse_args()

    run_api(args.input_yaml)

if __name__ == "__main__":
    run()

