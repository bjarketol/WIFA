"""Work around FraunhoferIWES/foxes#65 in foxes' Iterative algorithm.

foxes' Iterative algorithm, which its windIO reader picks whenever a blockage
model is set, applies turbine types and operating flags to the wrong turbines
after the first iteration: FarmWakesCalculation gets the per-turbine model data
in farm order but the farm data in downwind order.

The fix is a subclass rather than a patched method so that it also reaches
foxes' worker processes when they are spawned rather than forked (Windows,
macOS, forkserver): they unpickle the model by its class, which imports this
module.  Importing foxes here is fine because only the foxes adapter imports
this module.
"""

import foxes.constants as FC
import foxes.variables as FV
import numpy as np
from foxes.algorithms.iterative import models as iterative_models
from foxes.algorithms.iterative.models.farm_wakes_calc import (
    FarmWakesCalculation as _FoxesFarmWakesCalculation,
)


class FarmWakesCalculation(_FoxesFarmWakesCalculation):
    """foxes' iterative FarmWakesCalculation with the model data in downwind
    order on every iteration.  Keeps foxes' class name, which foxes uses to
    look up calculation parameters."""

    def calculate(self, algo, mdata, fdata):
        # Iteration 0, and any iteration that re-runs the full model list,
        # went through InitFarmData, which already put this chunk's mdata in
        # downwind order; later iterations get a fresh chunk in farm order.
        if algo.iterations and not algo._reamb:
            order = fdata[FV.ORDER].astype(int)
            ssel = np.broadcast_to(np.arange(order.shape[0])[:, None], order.shape)
            for k in mdata.keys():
                if tuple(mdata.dims[k][:2]) == (FC.STATE, FC.TURBINE) and np.any(
                    mdata[k] != mdata[k][0, 0, None, None]
                ):
                    mdata[k][:] = mdata[k][ssel, order]
        return super().calculate(algo, mdata, fdata)


def install():
    """Make foxes' Iterative algorithm build the fixed FarmWakesCalculation.

    Iterative looks its helper models up by name in
    ``foxes.algorithms.iterative.models``.  Idempotent.
    """
    iterative_models.FarmWakesCalculation = FarmWakesCalculation
