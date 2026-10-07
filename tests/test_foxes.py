import os
from pathlib import Path
from shutil import rmtree

import numpy as np
import pytest

pytest.importorskip(
    "foxes", reason="foxes not installed, install with: pip install wifa[foxes]"
)

from windIO import __path__ as wiop
from windIO import validate as validate_yaml

from wifa.foxes_api import _map_rotor_averaging, run_foxes

test_path = Path(os.path.dirname(__file__))
windIO_path = Path(wiop[0])


def _run_foxes(wes_dir):
    assert wes_dir.is_dir(), f"{wes_dir} is not a directory"

    for yaml_input in wes_dir.glob("system.yaml"):
        print("\nRUNNING FOXES ON", yaml_input, "\n")
        validate_yaml(yaml_input, Path("plant/wind_energy_system"))
        output_dir_name = Path("output_test_foxes")
        output_dir_name.mkdir(parents=True, exist_ok=True)
        run_foxes(yaml_input, output_dir=output_dir_name)
        rmtree(output_dir_name)


def test_foxes_KUL():
    wes_dir = test_path / "../examples/cases/KUL_LES/wind_energy_system/"
    _run_foxes(wes_dir)


def test_foxes_4wts():
    wes_dir = test_path / "../examples/cases/windio_4turbines/wind_energy_system/"
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


def test_foxes_heterogeneous_wind_rose_at_turbines():
    wes_dir = (
        test_path
        / "../examples/cases/heterogeneous_wind_rose_at_turbines/wind_energy_system/"
    )
    _run_foxes(wes_dir)


def test_foxes_heterogeneous_wind_rose_map():
    wes_dir = (
        test_path / "../examples/cases/heterogeneous_wind_rose_map/wind_energy_system/"
    )
    _run_foxes(wes_dir)


def test_foxes_simple_wind_rose():
    wes_dir = test_path / "../examples/cases/simple_wind_rose/wind_energy_system/"
    _run_foxes(wes_dir)


def test_foxes_timeseries_with_operating_flag():
    wes_dir = (
        test_path
        / "../examples/cases/timeseries_with_operating_flag/wind_energy_system/"
    )
    _run_foxes(wes_dir)


def test_timeseries_per_turbine_with_density(tmp_path=Path(".")):
    import foxes.variables as FV
    from conftest import make_timeseries_per_turbine_system_dict

    # Run with density
    system_dict = make_timeseries_per_turbine_system_dict("foxes")
    output_dir = tmp_path / "output_foxes_ts"
    farm_results = run_foxes(system_dict, verbosity=0, output_dir=str(output_dir))[0]
    farmP_with = farm_results[FV.P].sum()
    # print("Farm power with density:", farmP_with)
    assert np.isfinite(farmP_with) and farmP_with > 0

    # Run without density — same config but density removed
    system_dict_no = make_timeseries_per_turbine_system_dict("foxes")
    del system_dict_no["site"]["energy_resource"]["wind_resource"]["density"]
    output_dir_no = tmp_path / "output_foxes_ts_no_density"
    farm_results_no = run_foxes(
        system_dict_no, verbosity=0, output_dir=str(output_dir_no)
    )[0]
    farmP_without = farm_results_no[FV.P].sum()
    # print("Farm power without density:", farmP_without)

    rmtree(output_dir)
    rmtree(output_dir_no)

    # Density correction should change AEP (test data varies around 1.225)
    assert farmP_with != farmP_without


if __name__ == "__main__":
    test_foxes_KUL()
    test_foxes_4wts()
    test_foxes_abl()
    test_foxes_abl_stable()
    test_foxes_profiles()
    test_foxes_heterogeneous_wind_rose_at_turbines()
    test_foxes_heterogeneous_wind_rose_map()
    test_foxes_simple_wind_rose()
    test_timeseries_per_turbine_with_density()


def test_map_rotor_averaging_from_name():
    """foxes' partial wakes follow windIO's rotor_averaging.name (which pyWake
    reads) unless wake_averaging is given explicitly."""

    def system(**rotor_averaging):
        return {"attributes": {"analysis": {"rotor_averaging": rotor_averaging}}}

    def wake_averaging(wio):
        return _map_rotor_averaging(wio)["attributes"]["analysis"]["rotor_averaging"][
            "wake_averaging"
        ]

    assert wake_averaging(system(name="area_overlap")) == "top_hat"
    assert wake_averaging(system(name="gaussian_overlap")) == "gaussian"
    assert wake_averaging(system(name="center")) == "centre"
    assert wake_averaging(system(name="eq_grid")) == "auto"
    assert wake_averaging(system(name="area_overlap", wake_averaging="grid")) == "grid"

    wio = system(name="area_overlap")
    _map_rotor_averaging(wio)
    assert "wake_averaging" not in wio["attributes"]["analysis"]["rotor_averaging"]
    assert _map_rotor_averaging({"attributes": {"analysis": {}}}) == {
        "attributes": {"analysis": {}}
    }


def _two_turbine_blockage_system(listing, turbine_types, operating=None):
    """Two turbines in a row with Rathmann blockage, wind from 12 directions.

    A blockage model makes foxes' windIO reader pick its Iterative algorithm.
    ``listing`` orders the physical turbines (x = 0 m, x = 1000 m) in the
    layout; ``turbine_types`` and ``operating`` are per physical turbine.
    """
    from conftest import make_mixed_type_timeseries_system_dict

    system = make_mixed_type_timeseries_system_dict("foxes")
    wd = np.arange(0.0, 360.0, 30.0)
    resource = {
        "time": list(range(len(wd))),
        "wind_speed": {"data": [9.0] * len(wd), "dims": ["time"]},
        "wind_direction": {"data": wd.tolist(), "dims": ["time"]},
        "turbulence_intensity": {"data": [0.06] * len(wd), "dims": ["time"]},
    }
    if operating is not None:
        resource["wind_turbine"] = [0, 1]
        resource["operating"] = {
            "data": [[operating[i] for i in listing]] * len(wd),
            "dims": ["time", "wind_turbine"],
        }
    system["site"]["energy_resource"]["wind_resource"] = resource
    system["wind_farm"]["layouts"] = [
        {
            "coordinates": {"x": [[0.0, 1000.0][i] for i in listing], "y": [0.0, 0.0]},
            "turbine_types": [turbine_types[i] for i in listing],
        }
    ]
    system["attributes"]["analysis"]["blockage_model"] = {"name": "Rathmann"}
    return system


def _power_per_physical_turbine(listing, turbine_types, tmp_path, operating=None):
    import foxes.variables as FV

    system = _two_turbine_blockage_system(listing, turbine_types, operating)
    out = tmp_path / ("output_" + "".join(map(str, listing)))
    P = run_foxes(system, verbosity=0, output_dir=str(out))[0][FV.P].values
    return P[:, np.argsort(listing)]


def test_foxes_iterative_turbine_types_follow_turbines(tmp_path):
    """Listing the turbines in a different order must not change any
    turbine's power under blockage (foxes issue #65: Iterative swapped the
    power/Ct curves of turbines whose downwind and farm order differ)."""
    P_fwd = _power_per_physical_turbine([0, 1], [0, 1], tmp_path)
    P_rev = _power_per_physical_turbine([1, 0], [0, 1], tmp_path)
    np.testing.assert_allclose(P_rev, P_fwd, rtol=1e-10)


def test_foxes_iterative_operating_flags_follow_turbines(tmp_path):
    """Under blockage the switched-off turbine stays off in every wind
    direction, also when it is not first in the downwind order (foxes #65)."""
    P = _power_per_physical_turbine([0, 1], [0, 0], tmp_path, operating=[0, 1])
    assert np.all(P[:, 0] == 0.0)
    assert np.all(P[:, 1] > 0.0)
