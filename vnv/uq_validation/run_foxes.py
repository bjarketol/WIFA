# Author: Jonas Schulte, Fraunhofer IWES
# Email: jonas.schulte@iwes.fraunhofer.de
# Date: 2023-12-20
#
# Please run this script on foxes branch eu_flow from github:
# https://github.com/FraunhoferIWES/foxes/tree/eu_flow

from pathlib import Path
from string import Template
import argparse
import foxes
import foxes.variables as FV

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--cases", help="The windio input folders", default=["windio_toy"], nargs="+")
    parser.add_argument("-r", "--runs", help="The foxes parameter sets", 
                        default=["A", "B", "C", "D", "E", "F", "G"], nargs="+")
    parser.add_argument("-single", help="Switch on single turbine farm", action="store_true")
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

                ydir = Path(f"{case}/wind_energy_system")

                with open(ydir/"FLOW_UQ_vnv_toy_study_wind_energy_system_foxes.yaml", "r") as f:
                    tmpl = Template(f.read())
                yfile = ydir/f"_temp.yaml"
                with open(yfile, "w") as f:
                    f.write(tmpl.substitute({
                        "FARM": "oneTurbine.yaml" if args.single else "FLOW_UQ_vnv_toy_study_wind_farm.yaml",
                        "CASE": case,
                        "RUN": run
                    }))
                
                odir = Path(f"foxes_results/{case}")
                ofile1 = odir/f"power_{case}_foxes_{run}.csv"
                odir.mkdir(parents=True, exist_ok=True)

                mbook, farm, states, algo, outs = foxes.input.windio.read_case(yfile, runner=runner)
                yfile.unlink()

                farm_results = algo.calc_farm()
                fres = farm_results.to_dataframe()[[FV.AMB_REWS, FV.REWS, FV.AMB_TI, FV.TI, FV.AMB_P, FV.P]]
                print(fres)
                print(fres.describe())
                print()
                print(farm_results.to_dataframe()[[FV.WD, FV.YAW, FV.AMB_REWS, FV.REWS, FV.CT, FV.AMB_TI, FV.TI, FV.P]])
                
                for o in outs:
                    o.name = f"{o.name}_{case}_foxes_{run}"
                    if args.single:
                        o.name += "1"
                    print("Running output", o.name)
                    o.create(farm_results=farm_results, out_dir=odir, 
                            auto_fnames=lambda fname: f"{o.name}.{'.'.join(str(fname).split('.')[1:])}")
        