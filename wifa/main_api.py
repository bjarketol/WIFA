import argparse
from pathlib import Path

import windIO
from windIO import load_yaml
from windIO import validate as validate_yaml

from .cs_api.cs_modules.csLaunch.cs_run_function import run_code_saturne
from .foxes_api import run_foxes
from .pywake_api import run_pywake
from .wayve_api import run_wayve

# Valid model names for error messages
VALID_MODELS = ["pywake", "foxes", "wayve", "codesaturne"]


def run_api(yaml_input):
    """Run wind farm simulation using the specified flow model.

    Args:
        yaml_input: Path to the input YAML file

    Raises:
        ValueError: If an invalid model is specified
    """
    # Validate input
    validate_yaml(yaml_input, windIO.__path__[0] + "/plant/wind_energy_system.yaml")

    # Load configuration
    yaml_dat = load_yaml(yaml_input)
    model_name = yaml_dat["attributes"]["flow_model"]["name"].lower()

    # Create output directory if specified
    output_spec = yaml_dat["attributes"].get("model_outputs_specification", {})
    output_dir = output_spec.get("output_folder")
    if output_dir:
        Path(output_dir).mkdir(parents=True, exist_ok=True)

    # Run the appropriate model
    if model_name == "pywake":
        run_pywake(yaml_input)

    elif model_name == "foxes":
        run_foxes(yaml_input)

    elif model_name == "wayve":
        run_wayve(yaml_input, output_dir or "output")

    elif model_name == "codesaturne":
        run_code_saturne(yaml_input, test_mode=True)

    else:
        raise ValueError(f"Invalid model '{model_name}'. Choose from: {VALID_MODELS}")


def run():
    parser = argparse.ArgumentParser()
    parser.add_argument("input_yaml", help="The input yaml file")
    args = parser.parse_args()

    run_api(args.input_yaml)


if __name__ == "__main__":
    run()
