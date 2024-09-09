import argparse

def run_foxes(input_yaml, verbosity=1):
    
    from foxes.input.windio import read_windio

    wio_runner = read_windio(input_yaml, verbosity=verbosity)

    with wio_runner as runner:
        runner.run()

def run():
    
    parser = argparse.ArgumentParser()
    parser.add_argument("input_yaml", help="The input yaml file")
    parser.add_argument(
        "-v", "--verbosity", 
        help="The verbosity level, 0 = silent", 
        type=int, default=1
    )
    args = parser.parse_args()
    
    run_foxes(args.input_yaml, verbosity=args.verbosity)
    
if __name__ == "__main__":
    run()
    