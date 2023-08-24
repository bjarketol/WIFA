import argparse
import foxes
import foxes.variables as FV

def runFoxes(input_yaml):

    site_pars = dict(
            fixed_vars={FV.RHO: 1.225}
    )

    keymap = {
            'Bastankhah’s Gaussian wake model (simplified version)': "Bastankhah_linear_k002"
    }
    ana_pars = dict(keymap=keymap)

    mbook, farm, states, algo = foxes.input.windio.read_case(
                                    input_yaml,
                                    site_pars=site_pars,
                                    ana_pars=ana_pars)

    farm_results = algo.calc_farm()

    # results by turbine
    o = foxes.output.FarmResultsEval(farm_results)
    o.add_capacity(algo)
    o.add_capacity(algo, ambient=True)
    o.add_efficiency()
    turbine_results = o.reduce_states(
        {
            FV.AMB_P: "mean",
            FV.P: "mean",
            FV.AMB_CAP: "mean",
            FV.CAP: "mean",
            FV.EFF: "mean",
        }
    )
    turbine_results[FV.AMB_YLD] = o.calc_turbine_yield(algo=algo, annual=True, ambient=True)
    turbine_results[FV.YLD] = o.calc_turbine_yield(algo=algo, annual=True)
    print("\nResults by turbine:\n")
    print(turbine_results)

    # power results
    P0 = o.calc_mean_farm_power(ambient=True)
    P = o.calc_mean_farm_power()
    print(f"\nFarm power        : {P/1000:.1f} MW")
    print(f"Farm ambient power: {P0/1000:.1f} MW")
    print(f"Farm efficiency   : {o.calc_farm_efficiency()*100:.2f} %")
    print(f"Annual farm yield : {turbine_results[FV.YLD].sum():.2f} GWh")
    return turbine_results[FV.YLD].sum()


