"""Parametrized unit tests for PyWake submodel configuration functions.

Tests each _configure_*() function in wifa/pywake_api.py to verify:
- Correct PyWake class is returned for each model name
- Parameters are passed through correctly
- NotImplementedError raised for unsupported names
- Case-insensitive matching works
"""

import pytest
from py_wake.deficit_models import (
    HybridInduction,
    RankineHalfBody,
    SelfSimilarityDeficit,
    SelfSimilarityDeficit2020,
    VortexCylinder,
    VortexDipole,
)
from py_wake.deficit_models.gaussian import (
    BastankhahGaussianDeficit,
    BlondelSuperGaussianDeficit2020,
    BlondelSuperGaussianDeficit2023,
    CarbajofuertesGaussianDeficit,
    NiayifarGaussianDeficit,
    TurboGaussianDeficit,
    ZongGaussianDeficit,
)
from py_wake.deficit_models.gcl import GCLDeficit
from py_wake.deficit_models.noj import NOJDeficit, NOJLocalDeficit, TurboNOJDeficit
from py_wake.deficit_models.rathmann import Rathmann
from py_wake.deficit_models.utils import ct2a_madsen, ct2a_mom1d
from py_wake.deflection_models import JimenezWakeDeflection
from py_wake.deflection_models.gcl_hill_vortex import GCLHillDeflection
from py_wake.ground_models.ground_models import Mirror
from py_wake.rotor_avg_models import (
    AreaOverlapAvgModel,
    CGIRotorAvg,
    EqGridRotorAvg,
    GaussianOverlapAvgModel,
    GQGridRotorAvg,
    GridRotorAvg,
    PolarGridRotorAvg,
    RotorCenter,
)
from py_wake.superposition_models import (
    CumulativeWakeSum,
    LinearSum,
    MaxSum,
    SqrMaxSum,
    SquaredSum,
    WeightedSum,
)
from py_wake.turbulence_models import (
    CrespoHernandez,
    STF2005TurbulenceModel,
    STF2017TurbulenceModel,
)
from py_wake.turbulence_models.gcl_turb import GCLTurbulence

from wifa.pywake_api import (
    DEFAULTS,
    _configure_blockage_model,
    _configure_deficit_model,
    _configure_deflection_model,
    _configure_rotor_averaging,
    _configure_superposition_model,
    _configure_turbulence_model,
    configure_wake_model,
    get_with_default,
)

# Default rotor diameter and hub height for deficit model tests
_RD = 126.0
_HH = 90.0


def _call_deficit(name, analysis_extra=None, analysis_top=None):
    """Helper to call _configure_deficit_model with minimal boilerplate.

    Returns ``(wake_model_class, deficit_args)``; the third element
    (``deficit_post_attrs``) is dropped for the common case.  ``analysis_top``
    merges keys at the ``analysis`` level (e.g. ``axial_induction_model``).
    """
    wind_deficit_model = {"name": name, **(analysis_extra or {})}
    analysis = {"wind_deficit_model": wind_deficit_model, **(analysis_top or {})}
    cls, args, _post = _configure_deficit_model({"name": name}, analysis, _RD, _HH)
    return cls, args


def _call_deficit_full(name, analysis_extra=None, analysis_top=None):
    """Like ``_call_deficit`` but returns the full 3-tuple including post-attrs."""
    wind_deficit_model = {"name": name, **(analysis_extra or {})}
    analysis = {"wind_deficit_model": wind_deficit_model, **(analysis_top or {})}
    return _configure_deficit_model({"name": name}, analysis, _RD, _HH)


# ---------------------------------------------------------------------------
# Deficit model tests
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "name,expected_class",
    [
        ("Jensen", NOJLocalDeficit),
        ("jensen", NOJLocalDeficit),
        ("JENSEN", NOJLocalDeficit),
        ("Bastankhah2014", BastankhahGaussianDeficit),
        ("bastankhah2014", BastankhahGaussianDeficit),
        ("BASTANKHAH2014", BastankhahGaussianDeficit),
        ("SuperGaussian", BlondelSuperGaussianDeficit2020),
        ("supergaussian", BlondelSuperGaussianDeficit2020),
        ("SuperGaussian2023", BlondelSuperGaussianDeficit2023),
        ("TurboPark", TurboGaussianDeficit),
        ("TurbOPark", TurboGaussianDeficit),
        ("turbopark", TurboGaussianDeficit),
        ("Niayifar2016", NiayifarGaussianDeficit),
        ("niayifar2016", NiayifarGaussianDeficit),
        ("Zong2020", ZongGaussianDeficit),
        ("zong2020", ZongGaussianDeficit),
        ("Carbajofuertes2018", CarbajofuertesGaussianDeficit),
        ("carbajofuertes2018", CarbajofuertesGaussianDeficit),
        ("TurboNOJ", TurboNOJDeficit),
        ("turbonoj", TurboNOJDeficit),
        ("GCL", GCLDeficit),
        ("gcl", GCLDeficit),
        ("NOJLocalDeficit", NOJLocalDeficit),
        ("nojlocaldeficit", NOJLocalDeficit),
        ("NOJLOCALDEFICIT", NOJLocalDeficit),
        ("Jensen_1983", NOJDeficit),
        ("jensen_1983", NOJDeficit),
        ("JENSEN_1983", NOJDeficit),
        ("NOJDeficit", NOJDeficit),
        ("nojdeficit", NOJDeficit),
        ("NOJDEFICIT", NOJDeficit),
    ],
)
def test_configure_deficit_model(name, expected_class):
    cls, args = _call_deficit(name)
    assert cls is expected_class


@pytest.mark.parametrize(
    "name,expected_class",
    [
        ("Jensen", NOJLocalDeficit),
        ("Bastankhah2014", BastankhahGaussianDeficit),
        ("SuperGaussian", BlondelSuperGaussianDeficit2020),
        ("SuperGaussian2023", BlondelSuperGaussianDeficit2023),
        ("TurboPark", TurboGaussianDeficit),
        ("Niayifar2016", NiayifarGaussianDeficit),
        ("Zong2020", ZongGaussianDeficit),
        ("Carbajofuertes2018", CarbajofuertesGaussianDeficit),
        ("TurboNOJ", TurboNOJDeficit),
        ("GCL", GCLDeficit),
        ("NOJLocalDeficit", NOJLocalDeficit),
        ("Jensen_1983", NOJDeficit),
        ("NOJDeficit", NOJDeficit),
    ],
)
def test_configure_deficit_model_instantiation(name, expected_class):
    """Verify returned kwargs can actually instantiate the model without TypeError."""
    cls, args = _call_deficit(name)
    instance = cls(**args)
    assert isinstance(instance, expected_class)


def test_configure_deficit_model_bastankhah2014_params():
    """Verify wake expansion and ceps params are passed for Bastankhah2014."""
    cls, args = _call_deficit(
        "Bastankhah2014",
        {"wake_expansion_coefficient": {"k": 0.04}, "ceps": 0.2},
    )
    assert cls is BastankhahGaussianDeficit
    assert args["k"] == 0.04
    assert args["ceps"] == 0.2


def test_configure_deficit_model_bastankhah2014_k_b():
    """Verify k_b wake expansion is used when present."""
    _, args = _call_deficit(
        "Bastankhah2014",
        {"wake_expansion_coefficient": {"k_a": 0.38, "k_b": 0.004}},
    )
    assert args["k"] == 0.004


def test_configure_deficit_model_jensen_k_b():
    """Verify Jensen k_a/k_b expansion params."""
    _, args = _call_deficit(
        "Jensen",
        {"wake_expansion_coefficient": {"k_a": 0.38, "k_b": 0.004}},
    )
    assert args["a"] == [0.38, 0.004]


def test_configure_deficit_model_nojlocaldeficit_k_b():
    """Verify NOJLocalDeficit k_a/k_b expansion params."""
    _, args = _call_deficit(
        "NOJLocalDeficit",
        {"wake_expansion_coefficient": {"k_a": 0.38, "k_b": 0.004}},
    )
    assert args["a"] == [0.38, 0.004]


def test_configure_deficit_model_jensen_1983_k():
    """Verify Jensen_1983 passes scalar k and does not pass use_effective_ws."""
    cls, args = _call_deficit(
        "Jensen_1983",
        {"wake_expansion_coefficient": {"k": 0.04}},
    )
    assert cls is NOJDeficit
    assert args["k"] == 0.04
    assert "use_effective_ws" not in args


def test_configure_deficit_model_jensen_1983_k_b():
    """Jensen_1983 (NOJDeficit) accepts k via k_b, since windIO's
    wake_expansion_coefficient has no scalar k field."""
    cls, args = _call_deficit(
        "Jensen_1983",
        {"wake_expansion_coefficient": {"k_a": 0.0, "k_b": 0.1}},
    )
    assert cls is NOJDeficit
    assert args["k"] == 0.1


def test_configure_deficit_model_gaussian_params_niayifar():
    """Verify Gaussian params pass through for Niayifar2016."""
    cls, args = _call_deficit(
        "Niayifar2016",
        {
            "wake_expansion_coefficient": {
                "k_a": 0.38,
                "k_b": 0.004,
                "free_stream_ti": False,
            },
            "ceps": 0.3,
        },
    )
    assert cls is NiayifarGaussianDeficit
    assert args["a"] == [0.38, 0.004]
    assert args["ceps"] == 0.3
    assert args["use_effective_ti"] is True


def test_configure_deficit_model_zong_no_ceps():
    """Verify Zong2020 does not pass ceps (unsupported)."""
    _, args = _call_deficit(
        "Zong2020",
        {"ceps": 0.3, "wake_expansion_coefficient": {"free_stream_ti": True}},
    )
    assert "ceps" not in args
    assert args["use_effective_ti"] is False


def test_configure_deficit_model_bastankhah2014_no_effective_ti():
    """Verify Bastankhah2014 does not pass use_effective_ti (unsupported)."""
    _, args = _call_deficit(
        "Bastankhah2014",
        {"wake_expansion_coefficient": {"k": 0.04, "free_stream_ti": False}},
    )
    assert "use_effective_ti" not in args
    assert args["k"] == 0.04


def test_configure_deficit_model_a_param_warns_on_scalar_k():
    """Verify warning when scalar k is provided for a=[k_a, k_b] models."""
    with pytest.warns(UserWarning, match="uses a="):
        _, args = _call_deficit(
            "Niayifar2016", {"wake_expansion_coefficient": {"k": 0.05}}
        )
    assert "k" not in args
    assert "a" not in args


def test_configure_deficit_model_a_param_warns_on_missing_k_a():
    """Verify warning when k_b is provided without k_a."""
    with pytest.warns(UserWarning, match="k_a not specified"):
        _, args = _call_deficit(
            "Zong2020", {"wake_expansion_coefficient": {"k_b": 0.004}}
        )
    assert args["a"] == [0, 0.004]


@pytest.mark.parametrize(
    "name,extra,expected_class",
    [
        (
            "Bastankhah2014",
            {"wake_expansion_coefficient": {"k": 0.04}, "ceps": 0.2},
            BastankhahGaussianDeficit,
        ),
        (
            "Niayifar2016",
            {
                "wake_expansion_coefficient": {
                    "k_a": 0.38,
                    "k_b": 0.004,
                    "free_stream_ti": False,
                },
                "ceps": 0.3,
            },
            NiayifarGaussianDeficit,
        ),
        (
            "Zong2020",
            {
                "wake_expansion_coefficient": {
                    "k_a": 0.38,
                    "k_b": 0.004,
                    "free_stream_ti": True,
                }
            },
            ZongGaussianDeficit,
        ),
        (
            "Jensen",
            {"wake_expansion_coefficient": {"k_a": 0.38, "k_b": 0.004}},
            NOJLocalDeficit,
        ),
        (
            "NOJLocalDeficit",
            {"wake_expansion_coefficient": {"k_a": 0.38, "k_b": 0.004}},
            NOJLocalDeficit,
        ),
        (
            "Jensen_1983",
            {"wake_expansion_coefficient": {"k": 0.04}},
            NOJDeficit,
        ),
    ],
)
def test_configure_deficit_model_instantiation_with_params(name, extra, expected_class):
    """Verify models with user-specified params can be instantiated."""
    cls, args = _call_deficit(name, extra)
    assert isinstance(cls(**args), expected_class)


def test_configure_deficit_model_turbonoj_A_param():
    """Verify TurboNOJ passes through the A parameter and can be instantiated."""
    cls, args = _call_deficit("TurboNOJ", {"A": 0.6})
    assert cls is TurboNOJDeficit
    assert args["A"] == 0.6
    instance = cls(**args)
    assert isinstance(instance, TurboNOJDeficit)


@pytest.mark.parametrize("name", ["Bastankhah2016", "bastankhah2016"])
def test_configure_deficit_model_bastankhah2016_not_implemented(name):
    with pytest.raises(NotImplementedError, match="Bastankhah2016"):
        _call_deficit(name)


def test_configure_deficit_model_unknown():
    with pytest.raises(NotImplementedError, match="NonexistentModel"):
        _call_deficit("NonexistentModel")


# ---------------------------------------------------------------------------
# TI reference flag (windIO free_stream_ti -> PyWake use_effective_ti)
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "name",
    [
        "Jensen",
        "Niayifar2016",
        "Carbajofuertes2018",
        "Zong2020",
        "TurbOPark",
        "SuperGaussian",
        "SuperGaussian2023",
    ],
)
@pytest.mark.parametrize("free_stream_ti,expected", [(False, True), (True, False)])
def test_free_stream_ti_inverts_to_use_effective_ti(name, free_stream_ti, expected):
    """free_stream_ti maps to use_effective_ti with inverted polarity, for every
    TI-capable deficit (including SuperGaussian and TurbOPark, which the previous
    narrow handling missed)."""
    _, args = _call_deficit(
        name, {"wake_expansion_coefficient": {"free_stream_ti": free_stream_ti}}
    )
    assert args["use_effective_ti"] is expected
    # kwargs must actually instantiate the model
    cls, _ = _call_deficit(name)
    cls(**args)


@pytest.mark.parametrize("name", ["Bastankhah2014"])
def test_free_stream_ti_ignored_for_non_ti_capable(name):
    """Deficits without a use_effective_ti param must not receive it, even if
    free_stream_ti is present (would raise TypeError on instantiation).

    (GCL was moved to TI_CAPABLE — GCLDeficit accepts use_effective_ti, as
    GCLLocal demonstrates.)"""
    _, args = _call_deficit(
        name, {"wake_expansion_coefficient": {"k_b": 0.04, "free_stream_ti": True}}
    )
    assert "use_effective_ti" not in args


def test_use_effective_ti_key_is_ignored():
    """The deprecated top-level use_effective_ti key is no longer read."""
    _, args = _call_deficit(
        "Niayifar2016",
        {
            "wake_expansion_coefficient": {"k_a": 0.38, "k_b": 0.004},
            "use_effective_ti": False,
        },
    )
    assert "use_effective_ti" not in args


# ---------------------------------------------------------------------------
# Deflection model tests
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "name,expected_class",
    [
        ("Jimenez", JimenezWakeDeflection),
        ("jimenez", JimenezWakeDeflection),
        ("JIMENEZ", JimenezWakeDeflection),
        ("GCLHill", GCLHillDeflection),
        ("gclhill", GCLHillDeflection),
        ("GCLhill", GCLHillDeflection),
    ],
)
def test_configure_deflection_model(name, expected_class):
    model = _configure_deflection_model({"name": name, "beta": 0.1})
    assert isinstance(model, expected_class)


@pytest.mark.parametrize("name", [None, "None", "none", "NONE"])
def test_configure_deflection_model_none(name):
    assert _configure_deflection_model({"name": name, "beta": 0.1}) is None


def test_configure_deflection_model_bastankhah2016():
    with pytest.raises(NotImplementedError, match="Bastankhah2016"):
        _configure_deflection_model({"name": "Bastankhah2016", "beta": 0.1})


def test_configure_deflection_model_unknown():
    with pytest.raises(NotImplementedError, match="UnknownDeflection"):
        _configure_deflection_model({"name": "UnknownDeflection", "beta": 0.1})


# ---------------------------------------------------------------------------
# Turbulence model tests
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "name,expected_class",
    [
        ("STF2005", STF2005TurbulenceModel),
        ("stf2005", STF2005TurbulenceModel),
        ("STF2017", STF2017TurbulenceModel),
        ("stf2017", STF2017TurbulenceModel),
        ("CrespoHernandez", CrespoHernandez),
        ("crespohernandez", CrespoHernandez),
        ("CRESPOHERNANDEZ", CrespoHernandez),
        ("IEC-TI-2019", STF2017TurbulenceModel),
        ("iec-ti-2019", STF2017TurbulenceModel),
        ("GCL", GCLTurbulence),
        ("gcl", GCLTurbulence),
    ],
)
def test_configure_turbulence_model(name, expected_class):
    data = {"name": name, "c1": 1.0, "c2": 1.0}
    assert isinstance(_configure_turbulence_model(data), expected_class)


@pytest.mark.parametrize("name", [None, "None", "none", "NONE"])
def test_configure_turbulence_model_none(name):
    assert _configure_turbulence_model({"name": name, "c1": 1.0, "c2": 1.0}) is None


def test_configure_turbulence_model_unknown():
    with pytest.raises(NotImplementedError, match="UnknownTurb"):
        _configure_turbulence_model({"name": "UnknownTurb", "c1": 1.0, "c2": 1.0})


# ---------------------------------------------------------------------------
# Superposition model tests
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "name,expected_class",
    [
        ("Linear", LinearSum),
        ("linear", LinearSum),
        ("LINEAR", LinearSum),
        ("Squared", SquaredSum),
        ("squared", SquaredSum),
        ("Max", MaxSum),
        ("max", MaxSum),
        ("Weighted", WeightedSum),
        ("weighted", WeightedSum),
        ("Cumulative", CumulativeWakeSum),
        ("cumulative", CumulativeWakeSum),
    ],
)
def test_configure_superposition_model(name, expected_class):
    assert isinstance(
        _configure_superposition_model({"ws_superposition": name}), expected_class
    )


def test_configure_superposition_model_product_not_implemented():
    with pytest.raises(NotImplementedError, match="Product"):
        _configure_superposition_model({"ws_superposition": "Product"})


def test_configure_superposition_model_vector_not_implemented():
    """Vector superposition is foxes-only; the pyWake path rejects it."""
    with pytest.raises(NotImplementedError, match="Vector"):
        _configure_superposition_model({"ws_superposition": "Vector"})


def test_configure_superposition_model_unknown():
    with pytest.raises(NotImplementedError, match="UnknownSuper"):
        _configure_superposition_model({"ws_superposition": "UnknownSuper"})


# ---------------------------------------------------------------------------
# Rotor averaging tests
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "name,expected_class",
    [
        ("Center", RotorCenter),
        ("center", RotorCenter),
        ("CENTER", RotorCenter),
        ("grid", GridRotorAvg),
        ("avg_deficit", GridRotorAvg),
        ("Avg_Deficit", GridRotorAvg),
        ("gaussian_overlap", GaussianOverlapAvgModel),
        ("gaussianoverlap", GaussianOverlapAvgModel),
        ("area_overlap", AreaOverlapAvgModel),
        ("areaoverlap", AreaOverlapAvgModel),
        ("EqGrid", EqGridRotorAvg),
        ("eqgrid", EqGridRotorAvg),
        ("GQGrid", GQGridRotorAvg),
        ("gqgrid", GQGridRotorAvg),
        ("PolarGrid", PolarGridRotorAvg),
        ("polargrid", PolarGridRotorAvg),
        ("CGI", CGIRotorAvg),
        ("cgi", CGIRotorAvg),
    ],
)
def test_configure_rotor_averaging(name, expected_class):
    assert isinstance(_configure_rotor_averaging({"name": name}), expected_class)


def test_configure_rotor_averaging_eqgrid_n():
    assert isinstance(
        _configure_rotor_averaging({"name": "EqGrid", "n": 9}), EqGridRotorAvg
    )


def test_configure_rotor_averaging_gqgrid_params():
    data = {"name": "GQGrid", "n_x_grid_points": 3, "n_y_grid_points": 5}
    model = _configure_rotor_averaging(data)
    assert isinstance(model, GQGridRotorAvg)
    # Verify custom params produced different nodes than defaults
    default = GQGridRotorAvg()
    assert len(model.nodes_x) != len(default.nodes_x)


def test_configure_rotor_averaging_cgi_n():
    assert isinstance(_configure_rotor_averaging({"name": "CGI", "n": 7}), CGIRotorAvg)


def test_configure_rotor_averaging_unknown():
    with pytest.raises(NotImplementedError, match="UnknownRotor"):
        _configure_rotor_averaging({"name": "UnknownRotor"})


# ---------------------------------------------------------------------------
# Blockage model tests
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "name,expected_class",
    [
        ("SelfSimilarityDeficit2020", SelfSimilarityDeficit2020),
        ("selfsimilaritydeficit2020", SelfSimilarityDeficit2020),
        ("SelfSimilarityDeficit", SelfSimilarityDeficit),
        ("selfsimilaritydeficit", SelfSimilarityDeficit),
        ("RankineHalfBody", RankineHalfBody),
        ("rankinehalfbody", RankineHalfBody),
        ("Rathmann", Rathmann),
        ("rathmann", Rathmann),
        ("VortexCylinder", VortexCylinder),
        ("vortexcylinder", VortexCylinder),
        ("VortexDipole", VortexDipole),
        ("vortexdipole", VortexDipole),
        ("HybridInduction", HybridInduction),
        ("hybridinduction", HybridInduction),
    ],
)
def test_configure_blockage_model(name, expected_class):
    model = _configure_blockage_model({"name": name, "ss_alpha": 0.888}, {})
    assert isinstance(model, expected_class)


@pytest.mark.parametrize("name", [None, "None", "none", "NONE"])
def test_configure_blockage_model_none(name):
    assert _configure_blockage_model({"name": name}, {}) is None


def test_configure_blockage_model_unknown():
    with pytest.raises(NotImplementedError, match="UnknownBlockage"):
        _configure_blockage_model({"name": "UnknownBlockage"}, {})


# ---------------------------------------------------------------------------
# get_with_default preserves extra user keys
# ---------------------------------------------------------------------------


def test_get_with_default_preserves_extra_keys():
    """Verify that get_with_default merges defaults without dropping user keys."""
    analysis = {
        "rotor_averaging": {
            "name": "GQGrid",
            "n_x_grid_points": 3,
            "n_y_grid_points": 5,
        },
    }
    result = get_with_default(analysis, "rotor_averaging", DEFAULTS)
    assert result["name"] == "GQGrid"
    assert result["n_x_grid_points"] == 3
    assert result["n_y_grid_points"] == 5


def test_get_with_default_rotor_avg_eqgrid_n():
    """Verify EqGrid 'n' param survives through get_with_default."""
    analysis = {"rotor_averaging": {"name": "EqGrid", "n": 9}}
    result = get_with_default(analysis, "rotor_averaging", DEFAULTS)
    model = _configure_rotor_averaging(result)
    assert isinstance(model, EqGridRotorAvg)
    assert result["n"] == 9


def test_get_with_default_fills_missing_keys():
    """Verify that missing keys are filled from defaults."""
    # deflection_model defaults have beta=0.1; user only provides name
    analysis = {"deflection_model": {"name": "Jimenez"}}
    result = get_with_default(analysis, "deflection_model", DEFAULTS)
    assert result["name"] == "Jimenez"
    assert result["beta"] == 0.1


def test_get_with_default_recursive_nested_dicts():
    """Verify recursive merge fills deep missing keys while preserving user extras."""
    nested_defaults = {
        "model": {
            "params": {"a": 1, "b": 2},
            "name": "default",
        }
    }
    # User provides partial nested dict (missing key "b") plus an extra key "c"
    data = {
        "model": {
            "params": {"a": 10, "c": 99},
            "name": "custom",
        }
    }
    result = get_with_default(data, "model", nested_defaults)
    assert result["name"] == "custom"
    assert result["params"]["a"] == 10  # user value preserved
    assert result["params"]["b"] == 2  # missing key filled from defaults
    assert result["params"]["c"] == 99  # extra user key preserved


# ---------------------------------------------------------------------------
# configure_wake_model return contract
# ---------------------------------------------------------------------------


def test_configure_wake_model_returns_wake_deficit_key():
    """Verify configure_wake_model returns wake_deficit_key for API compat."""
    system_dat = {
        "attributes": {
            "analysis": {
                "wind_deficit_model": {"name": "Jensen"},
            }
        }
    }
    config = configure_wake_model(system_dat, rotor_diameter=126.0, hub_height=90.0)
    assert "wake_deficit_key" in config
    assert config["wake_deficit_key"] is None


def _weighted_system(rotor_name, ws_superposition="Weighted"):
    return {
        "attributes": {
            "analysis": {
                "wind_deficit_model": {
                    "name": "Zong2020",
                    "wake_expansion_coefficient": {"k_a": 0.38, "k_b": 0.004},
                },
                "superposition_model": {"ws_superposition": ws_superposition},
                "rotor_averaging": {"name": rotor_name},
                "deflection_model": {"name": "None"},
                "turbulence_model": {"name": "None"},
                "blockage_model": {"name": None},
            }
        }
    }


@pytest.mark.parametrize("rotor_name", ["gaussian_overlap", "area_overlap", "center"])
@pytest.mark.parametrize("ws_superposition", ["Weighted", "Cumulative"])
def test_weighted_superposition_requires_node_rotor_avg(rotor_name, ws_superposition):
    """WeightedSum/CumulativeWakeSum with a non-node rotor-averaging model raises
    a clear ValueError instead of PyWake's deep AssertionError."""
    with pytest.raises(ValueError, match="node"):
        configure_wake_model(
            _weighted_system(rotor_name, ws_superposition),
            rotor_diameter=126.0,
            hub_height=90.0,
        )


@pytest.mark.parametrize("rotor_name", ["grid", "eq_grid", "gq_grid", "cgi"])
def test_weighted_superposition_allows_node_rotor_avg(rotor_name):
    """A node rotor-averaging model is accepted with Weighted superposition."""
    config = configure_wake_model(
        _weighted_system(rotor_name), rotor_diameter=126.0, hub_height=90.0
    )
    assert isinstance(config["superposition_model"], WeightedSum)


def _weighted_deficit_system(deficit_name):
    """Weighted superposition + node rotor-avg, varying only the deficit."""
    return {
        "attributes": {
            "analysis": {
                "wind_deficit_model": {
                    "name": deficit_name,
                    "wake_expansion_coefficient": {"free_stream_ti": True},
                },
                "superposition_model": {"ws_superposition": "Weighted"},
                "rotor_averaging": {"name": "grid"},
                "deflection_model": {"name": "None"},
                "turbulence_model": {"name": "None"},
                "blockage_model": {"name": None},
            }
        }
    }


@pytest.mark.parametrize("deficit_name", ["SuperGaussian", "SuperGaussian2023", "GCL"])
def test_weighted_superposition_requires_convection_deficit(deficit_name):
    """WeightedSum/CumulativeWakeSum with a non-ConvectionDeficitModel deficit
    (super-Gaussian, GCL) raises a clear ValueError instead of PyWake's deep
    AssertionError, even when the rotor-averaging model is a node model."""
    with pytest.raises(ValueError, match="ConvectionDeficitModel"):
        configure_wake_model(
            _weighted_deficit_system(deficit_name),
            rotor_diameter=126.0,
            hub_height=90.0,
        )


@pytest.mark.parametrize("deficit_name", ["Zong2020", "Niayifar2016", "Bastankhah2014"])
def test_weighted_superposition_allows_convection_deficit(deficit_name):
    """Convection-based deficits are accepted with Weighted superposition."""
    config = configure_wake_model(
        _weighted_deficit_system(deficit_name), rotor_diameter=126.0, hub_height=90.0
    )
    assert isinstance(config["superposition_model"], WeightedSum)


# ---------------------------------------------------------------------------
# CrespoHernandez calibration coefficients (Phase 2 / Fix C)
# ---------------------------------------------------------------------------


def test_crespo_default_without_c():
    """No c -> the PyWake-default CrespoHernandez (ct2a_madsen)."""
    tm = _configure_turbulence_model({"name": "CrespoHernandez"})
    assert isinstance(tm, CrespoHernandez)
    assert tm.ct2a is ct2a_madsen


def test_crespo_with_c_uses_literature_recipe():
    """c -> CrespoHernandez with those coefficients, 1D induction and SqrMaxSum."""
    c = [0.73, 0.83, 0.03, -0.32]
    tm = _configure_turbulence_model({"name": "CrespoHernandez", "c": c})
    assert isinstance(tm, CrespoHernandez)
    assert list(tm.c) == c
    assert tm.ct2a is ct2a_mom1d
    assert isinstance(tm.addedTurbulenceSuperpositionModel, SqrMaxSum)


# ---------------------------------------------------------------------------
# 'none' rotor averaging + WeightedSum (Phase 2 / Fix D)
# ---------------------------------------------------------------------------


def test_rotor_averaging_none():
    assert _configure_rotor_averaging({"name": "none"}) is None


def test_weighted_superposition_allows_none_rotor():
    """WeightedSum accepts rotorAvgModel=None (rotor centre), as Zong (2020) uses."""
    sd = _weighted_deficit_system("Zong2020")
    sd["attributes"]["analysis"]["rotor_averaging"] = {"name": "none"}
    config = configure_wake_model(sd, rotor_diameter=126.0, hub_height=90.0)
    assert config["rotor_averaging"] is None
    assert isinstance(config["superposition_model"], WeightedSum)


def test_weighted_superposition_rejects_center_rotor():
    """A non-node, non-None rotor model is still rejected for WeightedSum."""
    sd = _weighted_deficit_system("Zong2020")
    sd["attributes"]["analysis"]["rotor_averaging"] = {"name": "center"}
    with pytest.raises(ValueError, match="node"):
        configure_wake_model(sd, rotor_diameter=126.0, hub_height=90.0)


# ---------------------------------------------------------------------------
# Zong ceps -> eps_coeff (Phase 2)
# ---------------------------------------------------------------------------


def test_zong_ceps_maps_to_eps_coeff():
    """Zong's near-wake epsilon is named eps_coeff (not ceps) in PyWake."""
    cls, args = _call_deficit(
        "Zong2020",
        {"wake_expansion_coefficient": {"k_a": 0.38, "k_b": 0.004}, "ceps": 0.35},
    )
    assert cls is ZongGaussianDeficit
    assert args["eps_coeff"] == 0.35
    assert "ceps" not in args


# ---------------------------------------------------------------------------
# Axial induction -> ct2a (honoring axial_induction_model)
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "axial, expected",
    [("1D", ct2a_mom1d), ("Madsen", ct2a_madsen), ("madsen", ct2a_madsen)],
)
def test_axial_induction_sets_ct2a(axial, expected):
    """axial_induction_model maps to the deficit's ct2a on a ct2a-capable model."""
    _, args = _call_deficit(
        "Bastankhah2014", analysis_top={"axial_induction_model": axial}
    )
    assert args["ct2a"] is expected


def test_axial_induction_default_absent_keeps_deficit_default():
    """No axial_induction_model -> no ct2a override (deficit keeps its default)."""
    _, args = _call_deficit("Bastankhah2014")
    assert "ct2a" not in args


def test_axial_induction_skipped_when_deficit_has_no_ct2a():
    """Blondel2020 has no ct2a parameter, so it is left untouched."""
    _, args = _call_deficit(
        "SuperGaussian", analysis_top={"axial_induction_model": "1D"}
    )
    assert "ct2a" not in args


def test_axial_induction_ct2a_instantiates():
    """The injected ct2a is a valid constructor argument."""
    cls, args = _call_deficit(
        "Niayifar2016", analysis_top={"axial_induction_model": "1D"}
    )
    assert isinstance(cls(**args), NiayifarGaussianDeficit)


# ---------------------------------------------------------------------------
# TurbOPark canonical recipe (Nygaard 2022)
# ---------------------------------------------------------------------------


def test_turbopark_recipe_ground_and_ctlim():
    """TurbOPark gets a Mirror ground model and ctlim=0.96 as constructor args."""
    cls, args, post = _call_deficit_full("TurbOPark")
    assert cls is TurboGaussianDeficit
    assert isinstance(args["groundModel"], Mirror)
    assert args["ctlim"] == 0.96


def test_turbopark_recipe_ws_key_post_attr():
    """TurbOPark scales by the downstream ambient WS via WS_key='WS_jlk'."""
    _, _, post = _call_deficit_full("TurbOPark")
    assert post == {"WS_key": "WS_jlk"}


def test_turbopark_recipe_instantiates_and_applies_ws_key():
    """The recipe builds a valid deficit and WS_key is set post-construction."""
    cls, args, post = _call_deficit_full("TurbOPark")
    deficit = cls(**args)
    for attr, value in post.items():
        setattr(deficit, attr, value)
    assert deficit.WS_key == "WS_jlk"


def test_non_turbopark_has_no_post_attrs():
    """Only TurbOPark carries post-construction attributes."""
    _, _, post = _call_deficit_full("Bastankhah2014")
    assert post == {}


# ---------------------------------------------------------------------------
# use_effective_ws / use_effective_ti for GCL (Phase 3)
# ---------------------------------------------------------------------------


def test_use_effective_ws_honored():
    """The windIO use_effective_ws flag is passed through (not hardcoded)."""
    _, args = _call_deficit("Bastankhah2014", {"use_effective_ws": False})
    assert args["use_effective_ws"] is False


def test_use_effective_ws_defaults_true():
    _, args = _call_deficit("Bastankhah2014")
    assert args["use_effective_ws"] is True


def test_gcl_honors_free_stream_ti():
    """GCLDeficit accepts use_effective_ti (GCLLocal); free_stream_ti is honored."""
    _, args = _call_deficit(
        "GCL", {"wake_expansion_coefficient": {"free_stream_ti": False}}
    )
    assert args["use_effective_ti"] is True
