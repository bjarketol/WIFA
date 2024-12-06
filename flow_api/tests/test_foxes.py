from flow_api.foxes_api import run_foxes
import windIO
from windIO.utils.yml_utils import validate_yaml
from pathlib import Path

from flow_api.tests import test_path, windIO_path

def test_foxes_KUL():

    yaml_input = test_path / '../../examples/cases/KUL_LES/wind_energy_system/system.yaml'

    # validate input
    validate_yaml(yaml_input, windIO_path / 'plant/wind_energy_system.yaml')

    # compute AEP (next step is to return a richer set of outputs)
    output_dir_name = Path('output_test_foxes')
    output_dir_name.mkdir(parents=True, exist_ok=True)
    foxes_aep = run_foxes(yaml_input, output_dir=output_dir_name)
    # print(foxes_aep)

def test_foxes_4wts():

    yaml_input = test_path / '../../examples/cases/windio_4turbines/wind_energy_system/system_no_field.yaml'

    # validate input
    validate_yaml(yaml_input, windIO_path / 'plant/wind_energy_system.yaml')

    # compute AEP (next step is to return a richer set of outputs)
    output_dir_name = Path('output_test_foxes')
    output_dir_name.mkdir(parents=True, exist_ok=True)
    foxes_aep = run_foxes(yaml_input, output_dir=output_dir_name)
    # print(foxes_aep)

def test_foxes_4wts2():

    yaml_input = test_path / '../../examples/cases/windio_4turbines_2flowcases/wind_energy_system/system.yaml'

    # validate input
    validate_yaml(yaml_input, windIO_path / 'plant/wind_energy_system.yaml')

    # compute AEP (next step is to return a richer set of outputs)
    output_dir_name = Path('output_test_foxes')
    output_dir_name.mkdir(parents=True, exist_ok=True)
    foxes_aep = run_foxes(yaml_input, output_dir=output_dir_name)
    # print(foxes_aep)
    
def test_foxes_abl():

    yaml_input = test_path / '../../examples/cases/windio_4turbines_ABL/wind_energy_system/system.yaml'

    # validate input
    validate_yaml(yaml_input, windIO_path / 'plant/wind_energy_system.yaml')

    # compute AEP (next step is to return a richer set of outputs)
    output_dir_name = Path('output_test_foxes')
    output_dir_name.mkdir(parents=True, exist_ok=True)
    foxes_aep = run_foxes(yaml_input, output_dir=output_dir_name)
    # print(foxes_aep)
    
    
def test_foxes_abl_stable():

    yaml_input = test_path / '../../examples/cases/windio_4turbines_ABL_stable/wind_energy_system/system.yaml'
    
    # validate input
    validate_yaml(yaml_input, windIO_path / 'plant/wind_energy_system.yaml')

    # compute AEP (next step is to return a richer set of outputs)
    output_dir_name = Path('output_test_foxes')
    output_dir_name.mkdir(parents=True, exist_ok=True)
    foxes_aep = run_foxes(yaml_input, output_dir=output_dir_name)
    # print(foxes_aep)
  
def test_foxes_profiles():

    yaml_input = test_path / '../../examples/cases/windio_4turbines_profiles_stable/wind_energy_system/system.yaml'
    
    # validate input
    validate_yaml(yaml_input, windIO_path / 'plant/wind_energy_system.yaml')

    # compute AEP (next step is to return a richer set of outputs)
    output_dir_name = Path('output_test_foxes')
    output_dir_name.mkdir(parents=True, exist_ok=True)
    foxes_aep = run_foxes(yaml_input, output_dir=output_dir_name)
    # print(foxes_aep)
    
if __name__ == "__main__":
    test_foxes_KUL()
    test_foxes_4wts()
    test_foxes_4wts2()
    test_foxes_abl()
    test_foxes_abl_stable()
    test_foxes_profiles()
    
