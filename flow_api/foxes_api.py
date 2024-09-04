import argparse
from foxes.input.windio import read_windio

def runFoxes(input_yaml):

    wio_runner = read_windio(input_yaml, verbosity=1)

    with wio_runner as runner:
        runner.run()

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "input_yaml",
        help="Path to the wind energy system yaml file",
    )
    args = parser.parse_args()
    
    runFoxes(args.input_yaml)
    