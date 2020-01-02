from __future__ import annotations
from typing import Dict, Any, TYPE_CHECKING
import numpy as np
from pyNastran.bdf.cards.optimization import get_dvprel_key

if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf import BDF


def get_dvprel_ndarrays(model: BDF, nelements: int, pids: np.ndarray,
                        fdtype: str='float32', idtype: str='int32') -> Dict[str, Any]:
    """
    Creates arrays for dvprel results

    Parameters
    ----------
    nelements : int
        the number of elements
    pids : (nelements,) int ndarray
        properties array to map the results to
    fdtype : str; default='float32'
        the type of the init/min/max arrays
    idtype : str; default='int32'
        the type of the design_region

    Returns
    -------
    dvprel_dict[key] : (design_region, dvprel_init, dvprel_min, dvprel_max)
        key : str
            the optimization string
        design_region : (nelements,) int ndarray
            the DVPRELx id
        dvprel_init : (nelements,) float ndarray
            the initial values of the variable
        dvprel_min : (nelements,)float ndarray
            the min values of the variable
        dvprel_max : (nelements,)float ndarray
            the max values of the variable

    """
    dvprel_dict = {}
    def get_dvprel_data(key):
        if key in dvprel_dict:
            return dvprel_dict[key]

        dvprel_t_init = np.full(nelements, np.nan, dtype=fdtype)
        dvprel_t_min = np.full(nelements, np.nan, dtype=fdtype)
        dvprel_t_max = np.full(nelements, np.nan, dtype=fdtype)
        design_region = np.zeros(nelements, dtype=idtype)
        dvprel_dict[key] = (design_region, dvprel_t_init, dvprel_t_min, dvprel_t_max)
        return design_region, dvprel_t_init, dvprel_t_min, dvprel_t_max

    for dvprel_key, dvprel in model.dvprels.items():
        prop_type = dvprel.prop_type
        unused_desvars = dvprel.dvids
        if dvprel.pid_ref is not None:
            pid = dvprel.pid_ref.pid
        else:
            pid = dvprel.pid
        unused_var_to_change = dvprel.pname_fid

        prop = model.properties[pid]
        if not prop.type == prop_type:
            raise RuntimeError('Property type mismatch\n%s%s' % (str(dvprel), str(prop)))

        key, msg = get_dvprel_key(dvprel, prop)
        if dvprel.type == 'DVPREL1':
            if msg:
                model.log.warning(msg)
                continue

            i = np.where(pids == pid)[0]
            if len(i) == 0:
                continue
            assert len(i) > 0, i
            design_region, dvprel_init, dvprel_min, dvprel_max = get_dvprel_data(key)

            optimization_region = dvprel.oid
            assert optimization_region > 0, str(model)
            design_region[i] = optimization_region
            xinit, lower_bound, upper_bound = dvprel.get_xinit_lower_upper_bound(model)

            dvprel_init[i] = xinit
            dvprel_min[i] = lower_bound
            dvprel_max[i] = upper_bound
        #elif dvprel.type == 'DVPREL2':
            #print(dvprel.get_stats())
        else:
            msg = 'dvprel.type=%r; dvprel=\n%s' % (dvprel.type, str(dvprel))
            raise NotImplementedError(msg)

        # TODO: haven't quite decided what to do
        if dvprel.p_max != 1e20:
            dvprel.p_max

        # TODO: haven't quite decided what to do
        if dvprel.p_min is not None:
            dvprel.p_min

    #dvprel_dict['PSHELL']['T']  = dvprel_t_init, dvprel_t_min, dvprel_t_max
    return dvprel_dict
