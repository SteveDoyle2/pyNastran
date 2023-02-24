from __future__ import annotations
from typing import Callable, Union, Optional, Any, TYPE_CHECKING
import h5py
import pandas as pd
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf import BDF

def passer(*args):
    pass

def h5py_to_dataframe(group: h5py._hl.dataset.Dataset) -> pd.DataFrame:
    data = {}
    for key in group.dtype.names:
        data[key] = group[key]
    return pd.DataFrame(data)

def get_tree(h5_group: Union[h5py._hl.group.Group,
                             h5py._hl.dataset.Dataset],
             mydict: Optional[dict[str, str]]=None) -> dict[str, str]:
    """
    Dataset:
    attrs  : <Attributes of HDF5 object at 1238427909672>
    chunks : None
    compression : None
    compression_opts : None
    dims   : <Dimensions of HDF5 object at 1238427909672>
    dtype  : dtype([('DOMAIN_ID', '<i8'), ('POSITION', '<i8'), ('LENGTH', '<i8')])
    external : None
    file   : <HDF5 file "6+element-nastran-sol103.h5" (mode r)>
    fillvalue : (0, 0, 0)
    fletcher32 : False
    id     : <h5py.h5d.DatasetID object at 0x00000120580E3E28>
    is_virtual : False
    maxshape : (3,)
    name   : '/INDEX/NASTRAN/RESULT/ELEMENTAL/STRESS/QUAD4'
    nbytes : 72
    ndim   : 1
    parent : <HDF5 group "/INDEX/NASTRAN/RESULT/ELEMENTAL/STRESS" (1 members)>
    ref    : <HDF5 object reference>
    regionref : <h5py._hl.base._RegionProxy object at 0x0000012058134E88>
    scaleoffset : None
    shape  : (3,)
    shuffle : False
    size   : 3
    """
    if mydict is None:
        mydict = {}
    allowed = (
        h5py._hl.group.Group,
        #h5py._hl.dataset.Dataset,
    )
    assert isinstance(h5_group, allowed), h5_group
    headers = list(h5_group)
    for header in headers:
        mydict2 = {}
        header_group = h5_group.get(header)
        if isinstance(header_group, h5py._hl.dataset.Dataset):
            #print(f'found header={header}')
            mydict[header] = None
            continue
        get_tree(header_group, mydict2)
        mydict[header] = mydict2
    return mydict
#def show_tree(h5_group: h5py.File) -> dict[str, str]:
    #assert isinstance(h5_group, h5py.File), h5_group
    #out = {}
    #headers = list(h5_group)
    #for header in headers:
        #header_group = h5_group.get(header)
        #mydict = {}
        #get_tree(header_group, mydict)
        #out[header] = mydict
    #return out

BasicCallable = Callable[
    [int, int, list[int], str],
    Any]

def read_basic_element(group: h5py._hl.dataset.Dataset, geom_model: BDF,
                        add_card: BasicCallable):
    assert len(group.dtype.names) == 4, group.dtype.names
    EID = group['EID']
    PID = group['PID']
    NIDS = group['G']
    DOMAIN_ID = group['DOMAIN_ID']
    for eid, pid, nids in zip(EID, PID, NIDS):
        obj = add_card(eid, pid, nids, comment='')
        obj.validate()
