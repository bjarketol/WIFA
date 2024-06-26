from foxes.input.windio import read_windio

def runFoxes(input_yaml):

    wio_runner = read_windio(input_yaml, verbosity=1)

    with wio_runner as runner:
        runner.run()
