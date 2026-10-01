import argparse
from pathlib import Path

from windIO import load_yaml

from wifa._optional import require

# windIO names the rotor averaging in rotor_averaging.name, which the pyWake
# adapter reads.  foxes' reader ignores it and takes its partial-wake model from
# wake_averaging alone, so the same file averaged wakes over the rotor area on
# pyWake and at the rotor centre on foxes.  Map the name to the foxes partial
# wakes that average the same way; anything else gets foxes' model default.
_FOXES_PARTIAL_WAKES = {
    "none": "centre",
    "center": "centre",
    "areaoverlap": "top_hat",
    "gaussianoverlap": "axiwake6",
}


def _map_rotor_averaging(wio):
    """Return *wio* with foxes' wake_averaging set from rotor_averaging.name.

    An explicit wake_averaging is left alone.  The input dict is not modified.
    """
    analysis = wio.get("attributes", {}).get("analysis", {})
    rotor_avg = analysis.get("rotor_averaging")
    if rotor_avg is None or "wake_averaging" in rotor_avg:
        return wio
    name = str(rotor_avg.get("name", "")).lower().replace("_", "").replace("-", "")
    rotor_avg = {**rotor_avg, "wake_averaging": _FOXES_PARTIAL_WAKES.get(name, "auto")}
    analysis = {**analysis, "rotor_averaging": rotor_avg}
    return {**wio, "attributes": {**wio["attributes"], "analysis": analysis}}


def run_foxes(
    input_yaml,
    input_dir=None,
    output_dir=None,
    engine="default",
    n_procs=None,
    chunksize_states=None,
    chunksize_points=None,
    verbosity=1,
    **kwargs,
):
    """
    Runs foxes based on windio yaml input

    Parameters
    ----------
    input_yaml: str or dict
        Path to the input data file, or the input data
    input_dir: str, optional
        The input base directory, for cases where
        input_yaml is a dict. In such cases it defaults to
        cwd, otherwise to the file containing directory
    output_dir: str, optional
        The output base directory, defaults to cwd
    engine: str
        The foxes engine choice
    n_procs; int, optional
        The number of processes to be used
    chunksize_states: int, optional
        The size of a states chunk
    chunksize_points: int, optional
        The size of a points chunk
    verbosity: int
        The verbosity level, 0 = silent
    kwargs: dict, optional
        Additional parameters for foxes.input.yaml.run_dict

    Returns
    -------
    farm_results: xarray.Dataset, optional
        The farm results
    point_results: xarray.Dataset, optional
        The point results, if requested by input_yaml
    outputs: list of tuple
        For each output enty, a tuple (dict, results),
        where results is a list that represents one
        entry per function call of the corresponding
        foxes output class

    """
    require("foxes")

    from foxes.input.yaml import run_dict
    from foxes.input.yaml.windio import read_windio_dict

    if isinstance(input_yaml, dict):
        wio = input_yaml
        idir = input_dir
    else:
        input_yaml = Path(input_yaml)
        # Keep an included wind_resource.nc as numpy arrays (foxes' reader uses
        # ndarrays directly), avoiding the dict-of-lists memory blow-up.
        wio = load_yaml(input_yaml, nc_data="array")
        idir = input_yaml.parent

    idict, algo, odir = read_windio_dict(_map_rotor_averaging(wio), verbosity=verbosity)

    if output_dir is not None:
        odir = output_dir

    if (
        engine is not None
        or n_procs is not None
        or chunksize_states is not None
        or chunksize_points is not None
    ):
        epars = dict(
            engine_type=engine,
            n_procs=n_procs,
            chunk_size_states=chunksize_states,
            chunk_size_points=chunksize_points,
            verbosity=verbosity,
        )
    else:
        epars = None

    return run_dict(
        idict,
        algo=algo,
        input_dir=idir,
        output_dir=odir,
        engine_pars=epars,
        verbosity=verbosity,
        **kwargs,
    )


def run():
    """
    Command line tool for running foxes from windio yaml file input.

    Examples
    --------
    >>> flow_api_foxes input.yaml

    """

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "input_yaml",
        help="The windio yaml file",
    )
    parser.add_argument(
        "-o",
        "--output_dir",
        help="The output directory",
        default=None,
    )
    parser.add_argument(
        "-e",
        "--engine",
        help="The engine",
        default=None,
    )
    parser.add_argument(
        "-n",
        "--n_procs",
        help="The number of processes",
        default=None,
        type=int,
    )
    parser.add_argument(
        "-c",
        "--chunksize_states",
        help="The chunk size for states",
        default=None,
        type=int,
    )
    parser.add_argument(
        "-C",
        "--chunksize_points",
        help="The chunk size for points",
        default=5000,
        type=int,
    )
    parser.add_argument(
        "-it",
        "--iterative",
        help="Use iterative algorithm",
        action="store_true",
    )
    parser.add_argument(
        "-v",
        "--verbosity",
        help="The verbosity level, 0 = silent",
        type=int,
        default=1,
    )
    args = parser.parse_args()

    run_foxes(
        input_yaml=args.input_yaml,
        output_dir=args.output_dir,
        engine=args.engine,
        n_procs=args.n_procs,
        chunksize_states=args.chunksize_states,
        chunksize_points=args.chunksize_points,
        verbosity=args.verbosity,
    )


if __name__ == "__main__":
    run()
