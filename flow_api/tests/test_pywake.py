from flow_api.pywake_api import run_pywake
from py_wake.tests import npt
from flow_api.tests import test_path
import os
import yaml
import sys
from py_wake.examples.data.hornsrev1 import Hornsrev1Site
import windIO
from windIO.utils.yml_utils import validate_yaml, Loader, load_yaml
sys.path.append(windIO.__path__[0])
from py_wake.examples.data.dtu10mw._dtu10mw import DTU10MW
from py_wake.deficit_models.gaussian import BastankhahGaussian
from py_wake.superposition_models import LinearSum
from py_wake.rotor_avg_models import RotorCenter


# todo
# - set up KUL with constant thrust turbine / wake model
# - set up two turbines case
#    - test vertical profile options: new file(s)?

def two_turbine_site():
    x = [0, 1248.1, 2496.2, 3744.3]
    y = [0, 0, 0, 0]
    ws = [10.09, 10.23]
    wd = [271.8, 266.2]
    turbine = DTU10MW()
    site = Hornsrev1Site()
    #deficit = BastankhahGaussianDeficit()
    wfm = BastankhahGaussian(site, turbine, k=0.04, 
                             use_effective_ws=True, 
                             superpositionModel=LinearSum(),
                             rotorAvgModel=RotorCenter())
    return wfm(x, y, ws=ws, wd=wd, time=True)


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

    wfm = two_turbine_site()
 
    # Check result
    pywake_aep_expected = wfm.aep().sum()
    npt.assert_array_almost_equal(pywake_aep, pywake_aep_expected, 0)


if __name__ == '__main__':
    test_pywake_KUL()
    test_pywake_4wts()
