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


_LANZILAO_ANALYSIS = {
    "wind_deficit_model": {
        "wake_expansion_coefficient": {"k_a": 0.004, "k_b": 0.38},
        "ceps": 0.2,
    }
}


def test_wake_tool_read_from_wm_coupling():
    """The windIO schema nests `wake_tool` under `wm_coupling` (a top-level
    `analysis.wake_tool` fails validation). Reading only the analysis level made
    the foxes coupling unreachable from any schema-valid file."""
    from wifa.wayve_api import wake_model_setup

    analysis = {**_LANZILAO_ANALYSIS, "wm_coupling": {"wake_tool": "wayve"}}
    wm = wake_model_setup(analysis)
    assert type(wm).__name__ == "Lanzilao"

    analysis = {**_LANZILAO_ANALYSIS, "wm_coupling": {"wake_tool": "nonsense"}}
    with pytest.raises(NotImplementedError, match="nonsense"):
        wake_model_setup(analysis)


def test_wake_tool_analysis_level_fallback():
    """Legacy hand-written yamls put `wake_tool` at the analysis level."""
    from wifa.wayve_api import wake_model_setup

    with pytest.raises(NotImplementedError, match="nonsense"):
        wake_model_setup({**_LANZILAO_ANALYSIS, "wake_tool": "nonsense"})


def test_wake_tool_defaults_to_wayve():
    from wifa.wayve_api import wake_model_setup

    wm = wake_model_setup(dict(_LANZILAO_ANALYSIS))
    assert type(wm).__name__ == "Lanzilao"


def _scalar_resource(ti):
    """Minimal single-point (no `height`) wind resource for flow_io_abl."""
    return {
        "wind_speed": {"data": [9.0]},
        "wind_direction": {"data": [270.0]},
        "turbulence_intensity": {"data": [ti]},
        "z0": {"data": [0.03]},
    }


def test_scalar_resource_ti_is_a_fraction():
    """windIO carries TI as a fraction. The scalar branch used to divide it by
    100 (and read it only when `z0` happened to be present), collapsing the
    TI-proportional wake expansion by two orders of magnitude."""
    from wifa.wayve_api import flow_io_abl

    abl = flow_io_abl(_scalar_resource(0.08), 0, zh=78.0, h1=156.0)
    assert abl.TI == pytest.approx(0.08)


def test_scalar_resource_ti_read_without_z0():
    """TI must be read whenever it is present, not only when `z0` is."""
    from wifa.wayve_api import flow_io_abl

    resource = _scalar_resource(0.08)
    del resource["z0"]
    abl = flow_io_abl(resource, 0, zh=78.0, h1=156.0)
    assert abl.TI == pytest.approx(0.08)


def test_scalar_resource_ti_defaults_when_absent():
    from wifa.wayve_api import flow_io_abl

    resource = _scalar_resource(0.08)
    del resource["turbulence_intensity"]
    abl = flow_io_abl(resource, 0, zh=78.0, h1=156.0)
    assert abl.TI == pytest.approx(0.04)


def test_scalar_resource_air_density_is_per_state():
    """A per-state `air_density` in the wind resource sets abl.rho, which
    linearly scales the reported turbine power (the APM velocity solution
    does not depend on it)."""
    from wifa.wayve_api import flow_io_abl

    resource = _scalar_resource(0.08)
    resource["air_density"] = {"data": [1.29, 1.11]}
    resource["wind_speed"]["data"].append(9.0)
    resource["wind_direction"]["data"].append(270.0)
    resource["turbulence_intensity"]["data"].append(0.08)
    resource["z0"]["data"].append(0.03)
    assert flow_io_abl(resource, 0, zh=78.0, h1=156.0).rho == pytest.approx(1.29)
    assert flow_io_abl(resource, 1, zh=78.0, h1=156.0).rho == pytest.approx(1.11)


def test_scalar_resource_air_density_defaults_when_absent():
    from wifa.wayve_api import flow_io_abl

    abl = flow_io_abl(_scalar_resource(0.08), 0, zh=78.0, h1=156.0)
    assert abl.rho == pytest.approx(1.225)


def _timeseries_system(wind_speeds, subset):
    """Single-turbine system whose wind speed differs at every timestep."""
    n = len(wind_speeds)
    return {
        "site": {
            "energy_resource": {
                "wind_resource": {
                    "time": list(range(n)),
                    "wind_speed": {"data": list(wind_speeds)},
                    "wind_direction": {"data": [270.0] * n},
                    "turbulence_intensity": {"data": [0.08] * n},
                    "z0": {"data": [0.03] * n},
                    "fc": {"data": [1.0e-4] * n},
                }
            }
        },
        "wind_farm": {
            "layouts": [{"coordinates": {"x": [0.0], "y": [0.0]}}],
            "turbines": {
                "performance": {
                    "power_curve": {
                        "power_values": [0, 1e6, 3e6, 3e6],
                        "power_wind_speeds": [3, 8, 12, 25],
                    },
                    "Ct_curve": {
                        "Ct_values": [0.8, 0.8, 0.4, 0.2],
                        "Ct_wind_speeds": [3, 8, 12, 25],
                    },
                },
                "hub_height": 80.0,
                "rotor_diameter": 80.0,
            },
        },
        "attributes": {
            "flow_model": {"name": "wayve"},
            "analysis": {
                "wind_deficit_model": {
                    "wake_expansion_coefficient": {"k_a": 0.04, "k_b": 0.0},
                    "ceps": 0.2,
                },
                "superposition_model": {"ws_superposition": "Product"},
                "wm_coupling": {"method": "PB"},
            },
            "model_outputs_specification": {
                "output_folder": "output",
                "run_configuration": {
                    "times_run": {"all_occurences": False, "subset": subset}
                },
                "turbine_outputs": {
                    "turbine_nc_filename": "turbine_data.nc",
                    "output_variables": ["rotor_effective_velocity"],
                },
            },
        },
    }


def test_times_run_subset_selects_the_requested_rows(tmp_path):
    """Regression: the subsetted timestamps were re-enumerated from 0, so the
    ABL was built from rows 0..n-1 of the wind resource while being *labelled*
    with the requested timestamps — silently simulating the wrong states."""
    import xarray as xr

    from wifa.wayve_api import run_wayve

    wind_speeds = [6.0, 9.0, 12.0]
    subset = [2]
    # debug_mode skips the APM solve: the wake model still sees the full ABL.
    run_wayve(
        _timeseries_system(wind_speeds, subset),
        output_dir=tmp_path,
        debug_mode=True,
    )

    ds = xr.load_dataset(tmp_path / "turbine_data.nc")
    assert ds["states"].values.tolist() == subset
    # A lone turbine is unwaked, so its rotor-effective velocity is the
    # background speed of the *requested* state (12 m/s), not of row 0 (6 m/s).
    rews = float(ds["rotor_effective_velocity"].isel(states=0, turbine=0))
    assert rews == pytest.approx(wind_speeds[subset[0]], rel=0.02)


if __name__ == "__main__":
    test_wayve_4wts()
