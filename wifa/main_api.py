from .foxes_api import run_foxes
from .pywake_api import run_pywake
from .cs_api.cs_modules.csLaunch.cs_run_function import run_code_saturne
from .wayve_api import run_wayve
import os
import yaml
import sys
import argparse

import windIO
from windIO import validate as validate_yaml, load_yaml

sys.path.append(windIO.__path__[0])


def run_api(yaml_input):
    # validate input
    validate_yaml(yaml_input, windIO.__path__[0] + "/schemas/plant/wind_energy_system")

    # get number of turbines
    yaml_dat = load_yaml(yaml_input)

    model_name = yaml_dat["attributes"]["flow_model"]["name"]

    if model_name.lower() == "pywake":
        pywake_aep = run_pywake(yaml_input)

    elif model_name.lower() == "foxes":
        foxes_aep = run_foxes(yaml_input)

    elif model_name.lower() == "wayve":
        # Output directory
        # yaml_input_no_ext = os.path.splitext(yaml_input)[0]  # Remove the file extension
        # output_dir_name = 'output_wayve' + yaml_input_no_ext.replace(os.sep, '_')  # Replace directory separators
        output_dir_name = yaml_dat["attributes"]["model_outputs_specification"][
            "output_folder"
        ]
        if not os.path.exists(output_dir_name):
            os.makedirs(output_dir_name)

        run_wayve(yaml_input, output_dir_name)

    elif model_name.lower() == "codesaturne":
        run_code_saturne(yaml_input, test_mode=True)

    else:
        print("Invalid Model")


def run():
    parser = argparse.ArgumentParser()
    parser.add_argument("input_yaml", help="The input yaml file")
    args = parser.parse_args()

    run_api(args.input_yaml)


if __name__ == "__main__":
    run()
