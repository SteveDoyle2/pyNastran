from __future__ import annotations
from typing import TYPE_CHECKING
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.op2.op2 import OP2


def get_is_slot_saved(op2: OP2, result_name: str) -> tuple[bool, dict]:
    if op2._results.is_not_saved(result_name):
        # op2.log.info(f'skipping {op2.table_name!r} due to {result_name!r}')
        return False, {}
    op2._results._found_result(result_name)
    slot = op2.get_result(result_name)
    return True, slot

def get_eid_dt_from_eid_device(eid_device: int,
                               dt: int | float,
                               sort_method: int,
                               ) -> tuple[int | float, int]:
    """common unvectorized method for transforming SORT2 data into SORT1"""
    if sort_method == 1:
        eid = eid_device // 10
        #print("SORT1 dt=%s eid_device=%s eid=%s" % (dt, eid_device, eid))
    else:
        eid, dt = dt, eid_device
    return eid, dt
