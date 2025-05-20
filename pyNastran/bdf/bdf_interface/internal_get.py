from pyNastran.bdf.cards.bdf_sets import SET1
from pyNastran.bdf.cards.optimization_nx import GROUP


def set_group(set_group_ref: GROUP | SET1,
              set_group_id: int) -> int:
    if hasattr(set_group_ref, 'group_id'):
        return set_group_ref.group_id
    if hasattr(set_group_ref, 'sid'):
        return set_group_ref.sid
    return set_group_id
