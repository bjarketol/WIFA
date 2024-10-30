from flow_api.foxes_api import run_foxes
import windIO
from windIO.utils.yml_utils import validate_yaml
from pathlib import Path

from flow_api.tests import test_path, windIO_path

def test_foxes_KUL():

    yaml_input = test_path / '../../examples/cases/KUL_LES/wind_energy_system/FLOW_UQ_vnv_toy_study_wind_energy_system.yaml'

    # validate input
    validate_yaml(yaml_input, windIO_path / 'plant/wind_energy_system.yaml')

    # compute AEP (next step is to return a richer set of outputs)
    output_dir_name = Path('output_foxes_4wts')
    output_dir_name.mkdir(parents=True, exist_ok=True)
    foxes_aep = run_foxes(yaml_input, output_dir=output_dir_name)
    # print(foxes_aep)

    # Check result
    #aep_expected = 7431.00357285578 #5448.73538851
    #assert abs(foxes_aep-aep_expected) < 1e-6, f"AEP MISMATCH, EXPECTED {aep_expected}, got {foxes_aep}"
    
    
if __name__ == "__main__":
    test_foxes_KUL()