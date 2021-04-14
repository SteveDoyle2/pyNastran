from __future__ import annotations
from typing import TYPE_CHECKING
#import numpy as np
import h5py
from .h5_utils import read_basic_element, h5py_to_dataframe
if TYPE_CHECKING:
    from pyNastran.bdf.bdf import BDF

def read_cshear(name: str, group: h5py._hl.dataset.Dataset, geom_model: BDF) -> None:
    #('EID', 'PID', 'G', 'DOMAIN_ID')
    read_basic_element(group, geom_model, geom_model.add_cshear)
    #EID = group['EID']
    #PID = group['PID']
    #NIDS = group['G']
    #DOMAIN_ID = group['DOMAIN_ID']
    #for eid, pid, nids in zip(EID, PID, NIDS):
        #obj = geom_model.add_cshear(eid, pid, nids, comment='')
        #obj.validate()

def read_ctria3(name: str, group: h5py._hl.dataset.Dataset, geom_model: BDF) -> None:
    geom_model.card_count[name] = len(group['EID'])
    setattr(geom_model, name, group)
    return
    EID = group['EID']
    PID = group['PID']
    NIDS = group['G']
    THETA = group['THETA']
    ZOFFS = group['ZOFFS']
    TFLAG = group['TFLAG']
    T = group['T']
    MCID = group['MCID']
    DOMAIN_ID = group['DOMAIN_ID']
    for eid, pid, nids, theta, zoffs, tflag, t, mcid in zip(EID, PID, NIDS, THETA, ZOFFS, TFLAG, T, MCID):
        if mcid == -1:
            theta_mcid = theta
        else:
            asdf
        assert tflag == 0, tflag
        t1, t2, t3 = [ti if ti != -1.0 else None
                      for ti in t]
        obj = geom_model.add_ctria3(
            eid, pid, nids,
            theta_mcid=theta_mcid, zoffset=zoffs, tflag=tflag,
            T1=t1, T2=t2, T3=t3, comment='')
        obj.validate()

def read_cquad4(name: str, group: h5py._hl.dataset.Dataset, geom_model: BDF) -> None:
    """
    Dataset:
    attrs  : <Attributes of HDF5 object at 1739457403880>
    chunks : (270,)
    compression : 'gzip'
    compression_opts : 1
    dims   : <Dimensions of HDF5 object at 1739457403880>
    dtype  : dtype([('EID', '<i8'), ('PID', '<i8'), ('G', '<i8', (4,)), ('THETA', '<f8'), ('ZOFFS', '<f8'), ('TFLAG', '<i8'), ('T', '<f8', (4,)), ('MCID', '<i8'), ('DOMAIN_ID', '<i8')])
    external : None
    file   : <HDF5 file "6+element-nastran-sol103.h5" (mode r)>
    fillvalue : (0, 0, [0, 0, 0, 0], 0., 0., 0, [0., 0., 0., 0.], 0, 0)
    fletcher32 : False
    id     : <h5py.h5d.DatasetID object at 0x00000194FFBD9BE8>
    is_virtual : False
    maxshape : (None,)
    name   : '/NASTRAN/INPUT/ELEMENT/CQUAD4'
    nbytes : 4320
    ndim   : 1
    parent : <HDF5 group "/NASTRAN/INPUT/ELEMENT" (1 members)>
    ref    : <HDF5 object reference>
    regionref : <h5py._hl.base._RegionProxy object at 0x00000194FF9AA048>
    scaleoffset : None
    shape  : (36,)
    shuffle : True
    size   : 36
    """
    geom_model.card_count[name] = len(group['EID'])
    setattr(geom_model, name, group)
    return
    EID = group['EID']
    PID = group['PID']
    NIDS = group['G']
    THETA = group['THETA']
    ZOFFS = group['ZOFFS']
    TFLAG = group['TFLAG']
    T = group['T']
    MCID = group['MCID']
    DOMAIN_ID = group['DOMAIN_ID']
    for eid, pid, nids, theta, zoffs, tflag, t, mcid in zip(EID, PID, NIDS, THETA, ZOFFS, TFLAG, T, MCID):
        if mcid == -1:
            theta_mcid = theta
        else:
            asdf
        assert tflag == 0, tflag
        t1, t2, t3, t4 = [ti if ti != -1.0 else None
                          for ti in t]
        obj = geom_model.add_cquad4(
            eid, pid, nids,
            theta_mcid=theta_mcid, zoffset=zoffs, tflag=tflag,
            T1=t1, T2=t2, T3=t3, T4=t4, comment='')
        obj.validate()

def read_ctria6(name: str, group: h5py._hl.dataset.Dataset, geom_model: BDF) -> None:
    EID = group['EID']
    PID = group['PID']
    NIDS = group['G']
    THETA = group['THETA']
    ZOFFS = group['ZOFFS']
    TFLAG = group['TFLAG']
    T = group['T']
    MCID = group['MCID']
    DOMAIN_ID = group['DOMAIN_ID']
    for eid, pid, nids, theta, zoffs, tflag, t, mcid in zip(EID, PID, NIDS, THETA, ZOFFS, TFLAG, T, MCID):
        if mcid == -1:
            theta_mcid = theta
        else:
            asdf
        assert tflag == 0, tflag
        t1, t2, t3 = [ti if ti != -1.0 else None
                      for ti in t]
        obj = geom_model.add_ctria6(
            eid, pid, nids, theta_mcid=theta_mcid, zoffset=zoffs,
            tflag=tflag, T1=t1, T2=t2, T3=t3, comment='')
        obj.validate()

def read_cquad8(name: str, group: h5py._hl.dataset.Dataset, geom_model: BDF) -> None:
    EID = group['EID']
    PID = group['PID']
    NIDS = group['G']
    THETA = group['THETA']
    ZOFFS = group['ZOFFS']
    TFLAG = group['TFLAG']
    T = group['T']
    MCID = group['MCID']
    DOMAIN_ID = group['DOMAIN_ID']
    for eid, pid, nids, theta, zoffs, tflag, t, mcid in zip(EID, PID, NIDS, THETA, ZOFFS, TFLAG, T, MCID):
        if mcid == -1:
            theta_mcid = theta
        else:
            asdf
        assert tflag == 0, tflag
        t1, t2, t3, t4 = [ti if ti != -1.0 else None
                          for ti in t]
        obj = geom_model.add_cquad8(
            eid, pid, nids, theta_mcid=theta_mcid, zoffset=zoffs,
            tflag=tflag, T1=t1, T2=t2, T3=t3, T4=t4, comment='')
        obj.validate()

def read_ctriar(name: str, group: h5py._hl.dataset.Dataset, geom_model: BDF) -> None:
    EID = group['EID']
    PID = group['PID']
    NIDS = group['G']
    THETA = group['THETA']
    ZOFFS = group['ZOFFS']
    TFLAG = group['TFLAG']
    T = group['T']
    MCID = group['MCID']
    DOMAIN_ID = group['DOMAIN_ID']
    for eid, pid, nids, theta, zoffs, tflag, t, mcid in zip(EID, PID, NIDS, THETA, ZOFFS, TFLAG, T, MCID):
        if mcid == -1:
            theta_mcid = theta
        else:
            asdf
        assert tflag == 0, tflag
        t1, t2, t3= [ti if ti != -1.0 else None
                     for ti in t]
        obj = geom_model.add_ctriar(
            eid, pid, nids,
            theta_mcid=theta_mcid, zoffset=zoffs, tflag=tflag,
            T1=t1, T2=t2, T3=t3, comment='')
        obj.validate()

def read_cquadr(name: str, group: h5py._hl.dataset.Dataset, geom_model: BDF) -> None:
    #('EID', 'PID', 'G', 'THETA', 'ZOFFS', 'TFLAG', 'T', 'MCID', 'DOMAIN_ID')
    EID = group['EID']
    PID = group['PID']
    NIDS = group['G']
    THETA = group['THETA']
    ZOFFS = group['ZOFFS']
    TFLAG = group['TFLAG']
    T = group['T']
    MCID = group['MCID']
    DOMAIN_ID = group['DOMAIN_ID']
    for eid, pid, nids, theta, zoffs, tflag, t, mcid in zip(EID, PID, NIDS, THETA, ZOFFS, TFLAG, T, MCID):
        if mcid == -1:
            theta_mcid = theta
        else:
            asdf
        assert tflag == 0, tflag
        t1, t2, t3, t4 = [ti if ti != -1.0 else None
                          for ti in t]
        obj = geom_model.add_cquadr(
            eid, pid, nids,
            theta_mcid=theta_mcid, zoffset=zoffs, tflag=tflag,
            T1=t1, T2=t2, T3=t3, T4=t4, comment='')
        obj.validate()
