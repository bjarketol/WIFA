from flow_api.pywake_api import run_pywake
from py_wake.tests import npt
from flow_api.tests import test_path
import os
import yaml
import sys

import windIO
from windIO.utils.yml_utils import validate_yaml, Loader, load_yaml
sys.path.append(windIO.__path__[0])


def test_pywake_KUL():

    yaml_input = test_path + '/../../examples/cases/KUL_LES/wind_energy_system/FLOW_UQ_vnv_toy_study_wind_energy_system.yaml'

    # validate input
    validate_yaml(yaml_input, windIO.__path__[0] + '/plant/wind_energy_system.yaml')

    # compute AEP (next step is to return a richer set of outputs)
    output_dir_name = 'output_pywake_4wts'
    if not os.path.exists(output_dir_name):
        os.makedirs(output_dir_name)
    pywake_aep = run_pywake(yaml_input, output_dir=output_dir_name)
    # print(pywake_aep)

    # Check result
    pywake_aep_expected = 5448.73538851
    npt.assert_array_almost_equal(pywake_aep, pywake_aep_expected, 6)


def test_pywake_4wts():

    yaml_input = test_path + '/../../examples/cases/windio_4turbines_2flowcases/wind_energy_system/FLOW_toy_study_wind_energy_system.yaml'

    # validate input
    validate_yaml(yaml_input, windIO.__path__[0] + '/plant/wind_energy_system.yaml')

    # compute AEP (next step is to return a richer set of outputs)
    output_dir_name = 'output_pywake_4wts'
    if not os.path.exists(output_dir_name):
        os.makedirs(output_dir_name)
    pywake_aep = run_pywake(yaml_input, output_dir=output_dir_name)
    # print(pywake_aep)

    # Check result
    pywake_aep_expected = 175.87182444
    npt.assert_array_almost_equal(pywake_aep, pywake_aep_expected, 6)


if __name__ == '__main__':
    test_pywake_KUL()
