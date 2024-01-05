# Author: Jonas Schulte, Fraunhofer IWES
# Email: jonas.schulte@iwes.fraunhofer.de
# Date: 2023-12-20
#
# Please run this script on foxes branch eu_flow from github:
# https://github.com/FraunhoferIWES/foxes/tree/eu_flow

from pathlib import Path
import argparse
import foxes
import foxes.variables as FV

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--cases", help="The windio input folders", default=["windio_toy", "windio_toy_profiles"], nargs="+")
    parser.add_argument("-r", "--runs", help="The foxes parameter sets", default=["A", "B", "C"], nargs="+")
    args = parser.parse_args()

    with foxes.utils.runners.DaskRunner() as runner:
        for case in args.cases:
            for run in args.runs:

                ydir = Path(f"{case}/wind_energy_system")
                yfile = ydir/f"FLOW_UQ_vnv_toy_study_wind_energy_system_foxes_{run}.yaml"
                odir = Path(f"foxes_results/{case}")
                ofile1 = odir/f"power_{case}_foxes_{run}.csv"
                odir.mkdir(parents=True, exist_ok=True)

                mbook, farm, states, algo, outs = foxes.input.windio.read_case(yfile, runner=runner)

                farm_results = algo.calc_farm()
                sres = farm_results
                fres = sres.to_dataframe()[[FV.AMB_REWS, FV.REWS, FV.AMB_TI, FV.TI, FV.AMB_P, FV.P]]
                print(fres)
                print(fres.describe())
                print()

                for o in outs:
                    o.name = f"{o.name}_{case}_foxes_{run}"
                    print("Running output", o.name)
                    o.create(farm_results=sres, out_dir=odir, 
                            auto_fnames=lambda fname: f"{o.name}{fname.suffix}")
        