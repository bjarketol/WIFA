from flow_api.foxes_api import run_foxes
from windIO.utils.yml_utils import validate_yaml
from pathlib import Path
from windIO import __path__ as wiop
import pytest
import os

test_path = Path(os.path.dirname(__file__))
windIO_path = Path(wiop[0])


def _run_foxes(wes_dir):
    for yaml_input in wes_dir.glob("system*.yaml"):
        print("\nRUNNING FOXES ON", yaml_input, "\n")
        validate_yaml(yaml_input, windIO_path / "plant/wind_energy_system.yaml")
        output_dir_name = Path("output_test_foxes")
        output_dir_name.mkdir(parents=True, exist_ok=True)
        run_foxes(yaml_input, output_dir=output_dir_name)


def test_foxes_KUL():
    wes_dir = test_path / "../examples/cases/KUL_LES/wind_energy_system/"
    _run_foxes(wes_dir)


@pytest.mark.skip(reason="We need to skip the no xy grid case")
def test_foxes_4wts():
    wes_dir = test_path / "../examples/cases/windio_4turbines/wind_energy_system/"
    _run_foxes(wes_dir)


@pytest.mark.skip(reason="This case was deleted")
def test_foxes_4wts2():
    wes_dir = (
        test_path / "../examples/cases/windio_4turbines_2flowcases/wind_energy_system/"
    )
    _run_foxes(wes_dir)


def test_foxes_abl():
    wes_dir = test_path / "../examples/cases/windio_4turbines_ABL/wind_energy_system/"
    _run_foxes(wes_dir)


def test_foxes_abl_stable():
    wes_dir = (
        test_path / "../examples/cases/windio_4turbines_ABL_stable/wind_energy_system/"
    )
    _run_foxes(wes_dir)


def test_foxes_profiles():
    wes_dir = (
        test_path
        / "../examples/cases/windio_4turbines_profiles_stable/wind_energy_system/"
    )
    _run_foxes(wes_dir)


if __name__ == "__main__":
    test_foxes_KUL()
    test_foxes_4wts()
    test_foxes_4wts2()
    test_foxes_abl()
    test_foxes_abl_stable()
    test_foxes_profiles()
