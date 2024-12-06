import argparse
from pathlib import Path
from windIO.utils.yml_utils import load_yaml

def run_foxes(
        input_yaml,
        input_dir=None, 
        output_dir=None, 
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
    verbosity: int
        The verbosity level, 0 = silent
    
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

    from foxes.input.yaml import read_windio, run_dict
    
    if isinstance(input_yaml, dict):
        wio = input_yaml
        idir = input_dir
    else:
        wio = load_yaml(input_yaml)
        idir = Path(input_yaml).parent

    idict, algo, odir = read_windio(wio, verbosity=verbosity)

    if output_dir is not None:
        odir = output_dir

    return run_dict(
        idict,
        algo=algo,
        input_dir=idir,
        output_dir=odir,
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

    if (
        args.engine is not None
        or args.n_procs is not None
        or args.chunksize_states is not None
        or args.chunksize_points is not None
    ):
        epars = dict(
            engine_type=args.engine,
            n_procs=args.n_procs,
            chunk_size_states=args.chunksize_states,
            chunk_size_points=args.chunksize_points,
            verbosity=args.verbosity,
        )
    else:
        epars = None

    run_foxes(
        input_yaml=args.input_yaml,
        output_dir=args.output_dir, 
        engine_pars=epars,
        verbosity=args.verbosity, 
    )

    
if __name__ == "__main__":
    run()
    