import argparse

def run_foxes(input_yaml):
    
    from foxes.input.windio import read_windio

    wio_runner = read_windio(input_yaml, verbosity=3)

    with wio_runner as runner:
        runner.run()

def run():
    
    parser = argparse.ArgumentParser()
    parser.add_argument("input_yaml", help="The input yaml file")
    args = parser.parse_args()
    
    run_foxes(args.input_yaml)
    
if __name__ == "__main__":
    run()
    