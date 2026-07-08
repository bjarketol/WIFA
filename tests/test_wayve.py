import os
from pathlib import Path

import pytest

pytest.importorskip(
    "wayve", reason="wayve not installed, install with: pip install wifa[wayve]"
)

from windIO import __path__ as wiop
from windIO import validate as validate_yaml

from wifa.wayve_api import run_wayve

test_path = Path(os.path.dirname(__file__))
windIO_path = Path(wiop[0])


# @pytest.mark.skip()
def test_wayve_4wts():
    yaml_input = (
        test_path / "../examples/cases/windio_4turbines/wind_energy_system/system.yaml"
    )
    validate_yaml(yaml_input, Path("plant/wind_energy_system"))
    output_dir_name = Path("output_test_wayve")
    output_dir_name.mkdir(parents=True, exist_ok=True)
    run_wayve(yaml_input, output_dir=output_dir_name, debug_mode=True)


def test_wayve_k_convention_lanzilao():
    """windIO defines k = k_a + k_b*TI; wayve's Lanzilao computes
    kwake = ka*TI + kb, so the reader must map swapped: ka <- k_b, kb <- k_a.
    Regression for the unswapped pass-through that silently turned the
    constant expansion into a TI multiplier (and vice versa)."""
    from wifa.wayve_api import wake_model_setup

    wm = wake_model_setup(
        {
            "wind_deficit_model": {
                "wake_expansion_coefficient": {"k_a": 0.004, "k_b": 0.38},
                "ceps": 0.2,
            }
        }
    )
    assert wm.ka == 0.38  # Lanzilao TI multiplier holds windIO k_b
    assert wm.kb == 0.004  # Lanzilao constant holds windIO k_a


def test_wayve_scalar_k_is_constant_expansion():
    """A scalar k is a constant expansion: Lanzilao kb, with no TI term."""
    from wifa.wayve_api import wake_model_setup

    wm = wake_model_setup(
        {
            "wind_deficit_model": {
                "wake_expansion_coefficient": {"k": 0.05},
                "ceps": 0.2,
            }
        }
    )
    assert wm.ka == 0.0
    assert wm.kb == 0.05


if __name__ == "__main__":
    test_wayve_4wts()
