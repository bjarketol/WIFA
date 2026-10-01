"""GCL wake deficit and added turbulence with the thrust coefficient capped.

GCL (Larsen 2009) takes ``m = 1 / sqrt(1 - Ct)`` from 1D momentum theory, so it
is undefined at Ct >= 1 and pyWake returns NaN there; its near-wake length
``x0`` already turns negative at Ct ~ 0.99 (low TI).  Manufacturer curves do
reach Ct > 1 near cut-in (1.11 at 4 m/s on the V82).  Ct is therefore capped
where 1D momentum theory stops holding (Glauert's turbulent-wake onset,
a = 0.4 -> Ct = 0.96), the limit pyWake's literature TurbOPark uses.  pyWake's
Gaussian deficits cap Ct the same way (``ctlim``); GCLDeficit has no such
parameter.
"""

from py_wake import np
from py_wake.deficit_models.gcl import GCLDeficit
from py_wake.turbulence_models.gcl_turb import GCLTurbulence

GCL_CT_LIMIT = 0.96


class CtLimitedGCLDeficit(GCLDeficit):
    def wake_radius(self, dw_ijlk, D_src_il, ct_ilk, **kwargs):
        ct_ilk = np.minimum(ct_ilk, GCL_CT_LIMIT)
        return super().wake_radius(dw_ijlk, D_src_il, ct_ilk, **kwargs)

    def calc_deficit(self, D_src_il, dw_ijlk, cw_ijlk, ct_ilk, **kwargs):
        ct_ilk = np.minimum(ct_ilk, GCL_CT_LIMIT)
        return super().calc_deficit(D_src_il, dw_ijlk, cw_ijlk, ct_ilk, **kwargs)


class CtLimitedGCLTurbulence(GCLTurbulence):
    def calc_added_turbulence(
        self, dw_ijlk, D_src_il, ct_ilk, wake_radius_ijlk, D_dst_ijl, cw_ijlk, **kwargs
    ):
        ct_ilk = np.minimum(ct_ilk, GCL_CT_LIMIT)
        return super().calc_added_turbulence(
            dw_ijlk, D_src_il, ct_ilk, wake_radius_ijlk, D_dst_ijl, cw_ijlk, **kwargs
        )
