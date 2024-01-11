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
    parser.add_argument("-c", "--cases", help="The windio input folders", default=["windio_toy"], nargs="+")
    parser.add_argument("-r", "--runs", help="The foxes parameter sets", 
                        default=["A", "B", "C", "D", "E", "F", "G", "A1", "B1", "C1", "D1", "E1", "F1", "G1"], 
                        nargs="+")
    parser.add_argument("-sc", "--scheduler", help="The scheduler choice", default=None)
    parser.add_argument("-n", "--n_workers", help="The number of workers for distributed run", type=int, default=None)
    parser.add_argument("-tw", "--threads_per_worker", help="The number of threads per worker for distributed run", type=int,default=None)
    args = parser.parse_args()

    with foxes.utils.runners.DaskRunner(
        scheduler=args.scheduler,
        n_workers=args.n_workers,
        threads_per_worker=args.threads_per_worker
    ) as runner:
        for case in args.cases:
            for run in args.runs:

                ydir = Path(f"{case}/wind_energy_system/foxes")
                yfile = ydir/f"FLOW_UQ_vnv_toy_study_wind_energy_system_foxes_{run}.yaml"
                odir = Path(f"foxes_results/{case}")
                ofile1 = odir/f"power_{case}_foxes_{run}.csv"
                odir.mkdir(parents=True, exist_ok=True)

                mbook, farm, states, algo, outs = foxes.input.windio.read_case(yfile, runner=runner)

                farm_results = algo.calc_farm()
                fres = farm_results.to_dataframe()[[FV.AMB_REWS, FV.REWS, FV.AMB_TI, FV.TI, FV.AMB_P, FV.P]]
                print(fres)
                print(fres.describe())
                print()
                print(farm_results.to_dataframe()[[FV.WD, FV.YAW, FV.AMB_REWS, FV.REWS, FV.CT, FV.AMB_TI, FV.TI, FV.P]])
                
                for o in outs:
                    o.name = f"{o.name}_{case}_foxes_{run}"
                    print("Running output", o.name)
                    o.create(farm_results=farm_results, out_dir=odir, 
                            auto_fnames=lambda fname: f"{o.name}.{'.'.join(str(fname).split('.')[1:])}")
        