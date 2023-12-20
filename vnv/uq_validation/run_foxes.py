from pathlib import Path
import argparse
import pandas as pd
import foxes
import foxes.variables as FV
import foxes.constants as FC

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--cases", help="The windio input folders", default=["windio_toy", "windio_toy_profiles"], nargs="+")
    parser.add_argument("-r", "--runs", help="The foxes parameter sets", default=["A", "B", "C"], nargs="+")
    args = parser.parse_args()

    for case in args.cases:
        for run in args.runs:

            ydir = Path(f"{case}/wind_energy_system")
            yfile = ydir/f"FLOW_UQ_vnv_toy_study_wind_energy_system_foxes_{run}.yaml"
            odir = Path(f"foxes_results/{case}")
            ofile1 = odir/f"power_{case}_foxes_{run}.csv"
            odir.mkdir(parents=True, exist_ok=True)

            mbook, farm, states, algo = foxes.input.windio.read_case(yfile)

            farm_results = algo.calc_farm()
            times = farm_results[FC.STATE].to_numpy()
            print(farm_results)
            fres = farm_results.sel(state=times[:100]).to_dataframe()
            fres = fres[[FV.AMB_REWS, FV.REWS, FV.AMB_TI, FV.TI, FV.AMB_P, FV.P]]
            print(fres)
            print(fres.describe())

            print("\nWriting file", ofile1)
            odata = None
            for ti, g in fres.reset_index().set_index("state").groupby(FC.TURBINE):
                if ti == 0:
                    odata = pd.DataFrame(index=g.index.to_numpy())
                    odata.index.name="# Time"
                odata[f" Power of turbine {ti+1}"] = g[FV.P].to_numpy()
            print(odata)
            odata.to_csv(ofile1, float_format="%.10e")

        