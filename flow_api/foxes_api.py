import argparse
from foxes.output import FarmResultsEval

def run_foxes(input_yaml, output_dir=None, ret_aep=True, verbosity=1):
    
    from foxes.input.windio import read_windio

    wio_runner = read_windio(input_yaml, output_dir=output_dir, verbosity=verbosity)

    with wio_runner as runner:
        runner.run()
    
        if ret_aep:
            o = FarmResultsEval(runner.farm_results)
            aep = o.calc_farm_yield(algo=runner.algo)
    
    return aep

def run():
    
    parser = argparse.ArgumentParser()
    parser.add_argument("input_yaml", help="The input yaml file")
    parser.add_argument(
        "-v", "--verbosity", 
        help="The verbosity level, 0 = silent", 
        type=int, default=1
    )
    args = parser.parse_args()
    
    aep = run_foxes(args.input_yaml, ret_aep=True, verbosity=args.verbosity)
    print(f"AEP = {aep:.2f} GWh")
    
if __name__ == "__main__":
    run()
    