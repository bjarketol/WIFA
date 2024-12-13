from wifa.foxes_api import run_foxes
from windIO.utils.yml_utils import validate_yaml
from pathlib import Path
from windIO import __path__ as wiop
from foxes import Engine

import os

test_path = Path(os.path.dirname(__file__))
windIO_path = Path(wiop[0])


class TestFoxes:

    def _run_foxes(self, wes_dir):
        if not hasattr(self, "_engine"):
            self._engine = Engine.new("threads", verbosity=0)
            self._engine.initialize()
            print("TestFoxes: Using engine", self._engine)

        for yaml_input in wes_dir.glob("system*.yaml"):
            if "_noXYgrid" not in str(yaml_input):
                print("\nRUNNING FOXES ON", yaml_input, "\n")
                validate_yaml(yaml_input, windIO_path / "plant/wind_energy_system.yaml")
                output_dir_name = Path("output_test_foxes")
                output_dir_name.mkdir(parents=True, exist_ok=True)
                try:
                    run_foxes(yaml_input, output_dir=output_dir_name, engine=None)
                except Exception as e:
                    self._engine.finalize()
                    raise e

    def test_foxes_KUL(self):
        wes_dir = test_path / "../examples/cases/KUL_LES/wind_energy_system/"
        self._run_foxes(wes_dir)

    def test_foxes_4wts(self):
        wes_dir = test_path / "../examples/cases/windio_4turbines/wind_energy_system/"
        self._run_foxes(wes_dir)

    def test_foxes_abl(self):
        wes_dir = test_path / "../examples/cases/windio_4turbines_ABL/wind_energy_system/"
        self._run_foxes(wes_dir)

    def test_foxes_abl_stable(self):
        wes_dir = (
            test_path / "../examples/cases/windio_4turbines_ABL_stable/wind_energy_system/"
        )
        self._run_foxes(wes_dir)

    def test_foxes_profiles(self):
        wes_dir = (
            test_path
            / "../examples/cases/windio_4turbines_profiles_stable/wind_energy_system/"
        )
        self._run_foxes(wes_dir)

    def __del__(self):
        if hasattr(self, "_engine"):
            self._engine.finalize()


if __name__ == "__main__":
    runner = TestFoxes()
    #test_foxes_KUL()
    runner.test_foxes_4wts()
    #test_foxes_4wts2()
    #test_foxes_abl()
    #test_foxes_abl_stable()
    #test_foxes_profiles()
