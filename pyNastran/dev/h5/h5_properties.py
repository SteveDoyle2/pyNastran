from __future__ import annotations
from typing import TYPE_CHECKING
import numpy as np
import h5py
from pyNastran.bdf.cards.properties.shell import map_failure_theory_int
from .h5_geometry import passer
if TYPE_CHECKING:
    from pyNastran.bdf.bdf import BDF

def read_pbush(name: str, group: h5py._hl.dataset.Dataset, geom_model: BDF):
    afsd

def read_pelas(name: str, group: h5py._hl.dataset.Dataset, geom_model: BDF):
    assert len(group.dtype.names) == 5, group.dtype.names
    PID = group['PID']
    K = group['K']
    GE = group['GE']
    S = group['S']
    DOMAIN_ID = group['DOMAIN_ID']
    for pid, k, ge, s in zip(PID, K, GE, S):
        obj = geom_model.add_pelas(pid, k, ge=ge, s=s, comment='')
        obj.validate()
        str(obj)

def read_pvisc(name: str, group: h5py._hl.dataset.Dataset, geom_model: BDF):
    assert len(group.dtype.names) == 4, group.dtype.names
    PID = group['PID']
    CE = group['CE']
    CR = group['CR']
    DOMAIN_ID = group['DOMAIN_ID']
    for pid, ce, cr in zip(PID, CE, CR):
        obj = geom_model.add_pvisc(pid, ce, cr, comment='')
        obj.validate()
        str(obj)

def read_pdamp(name: str, group: h5py._hl.dataset.Dataset, geom_model: BDF):
    #('PID', 'B', 'DOMAIN_ID')
    assert len(group.dtype.names) == 3, group.dtype.names
    PID = group['PID']
    B = group['B']
    DOMAIN_ID = group['DOMAIN_ID']
    for pid, b in zip(PID, B):
        obj = geom_model.add_pdamp(pid, b, comment='')
        obj.validate()

def read_prod(name: str, group: h5py._hl.dataset.Dataset, geom_model: BDF):
    assert len(group.dtype.names) == 7, group.dtype.names
    PID = group['PID']
    MID = group['MID']
    A = group['A']
    J = group['J']
    C = group['C']
    NSM = group['NSM']
    DOMAIN_ID = group['DOMAIN_ID']
    for pid, mid, a, j, c, nsm in zip(PID, MID, A, J, C, NSM):
        obj = geom_model.add_prod(pid, mid, a, j=j, c=c, nsm=nsm, comment='')
        obj.validate()

def read_ptube(name: str, group: h5py._hl.dataset.Dataset, geom_model: BDF):
    #('PID', 'MID', 'OD', 'T', 'NSM', 'DOMAIN_ID')
    assert len(group.dtype.names) == 6, group.dtype.names
    PID = group['PID']
    MID = group['MID']
    OD = group['OD']
    T = group['T']
    NSM = group['NSM']
    DOMAIN_ID = group['DOMAIN_ID']
    for pid, mid, OD1, t, nsm in zip(PID, MID, OD, T, NSM):
        obj = geom_model.add_ptube(pid, mid, OD1, t=t, nsm=nsm, OD2=None, comment='')
        obj.validate()
        str(obj)

def read_pshear(name: str, group: h5py._hl.dataset.Dataset, geom_model: BDF):
    #('PID', 'MID', 'T', 'NSM', 'F1', 'F2', 'DOMAIN_ID')
    assert len(group.dtype.names) == 7, group.dtype.names
    PID = group['PID']
    MID = group['MID']
    T = group['T']
    NSM = group['NSM']
    F1 = group['F1']
    F2 = group['F2']
    DOMAIN_ID = group['DOMAIN_ID']
    for pid, mid, t, nsm, f1, f2 in zip(PID, MID, T, NSM, F1, F2):
        obj = geom_model.add_pshear(pid, mid, t, nsm=nsm, f1=f1, f2=f2, comment='')
        obj.validate()
        str(obj)

def read_pshell(name: str, group: h5py._hl.dataset.Dataset, geom_model: BDF):
    """
    Dataset:
    attrs  : <Attributes of HDF5 object at 3135747352808>
    chunks : (340,)
    compression : 'gzip'
    compression_opts : 1
    dims   : <Dimensions of HDF5 object at 3135747352808>
    dtype  : dtype([('PID', '<i8'), ('MID1', '<i8'), ('T', '<f8'), ('MID2', '<i8'), ('BK', '<f8'), ('MID3', '<i8'), ('TS', '<f8'), ('NSM', '<f8'), ('Z1', '<f8'), ('Z2', '<f8'), ('MID4', '<i8'), ('DOMAIN_ID', '<i8')])
    external : None
    file   : <HDF5 file "6+element-nastran-sol103.h5" (mode r)>
    fillvalue : (0, 0, 0., 0, 0., 0, 0., 0., 0., 0., 0, 0)
    fletcher32 : False
    id     : <h5py.h5d.DatasetID object at 0x000002DA191B68E8>
    is_virtual : False
    maxshape : (None,)
    name   : '/NASTRAN/INPUT/PROPERTY/PSHELL'
    nbytes : 96
    ndim   : 1
    parent : <HDF5 group "/NASTRAN/INPUT/PROPERTY" (1 members)>
    ref    : <HDF5 object reference>
    regionref : <h5py._hl.base._RegionProxy object at 0x000002DA191E6D08>
    scaleoffset : None
    shape  : (1,)
    shuffle : True
    size   : 1
    """
    #'PID', 'MID1', 'T', 'MID2', 'BK', 'MID3', 'TS', 'NSM', 'Z1', 'Z2', 'MID4', 'DOMAIN_ID',
    assert len(group.dtype.names) == 12, group.dtype.names
    PID = group['PID']
    MID1 = group['MID1']
    MID2 = group['MID2']
    MID3 = group['MID3']
    MID4 = group['MID4']
    T = group['T']
    BK = group['BK']
    TS = group['TS']
    NSM = group['NSM']
    Z1 = group['Z1']
    Z2 = group['Z2']
    DOMAIN_ID = group['DOMAIN_ID']
    for pid, mid1, mid2, mid3, mid4, t, bk, ts, nsm, z1, z2 in zip(PID, MID1, MID2, MID3, MID4, T, BK, TS, NSM, Z1, Z2):
        if max([mid1, mid2, mid3, mid4]) > 100_000_000:
            continue
        obj = geom_model.add_pshell(
            pid, mid1=mid1, t=t, mid2=mid2, twelveIt3=bk, mid3=mid3,
            tst=ts, nsm=nsm, z1=z1, z2=z2, mid4=mid4, comment='')
        obj.validate()
        str(obj)

def read_pcomp(name: str, group: h5py._hl.dataset.Dataset, geom_model: BDF):
    identity = group.get('IDENTITY')
    ply = group.get('PLY')

    PID = identity['PID']
    NPLIES = identity['NPLIES']
    Z0 = identity['Z0']
    NSM = identity['NSM']
    SB = identity['SB']
    FT = identity['FT']
    TREF = identity['TREF']
    GE = identity['GE']
    PLY_POS = identity['PLY_POS']
    PLY_LEN = identity['PLY_LEN']
    DOMAIN_ID = identity['DOMAIN_ID']

    MID = ply['MID']
    T = ply['T']
    THETA = ply['THETA']
    SOUT = ply['SOUT']
    for (pid, nplies, z0, nsm, sb, ft, tref, ge, ply_pos, ply_len) in zip(
         PID, NPLIES, Z0, NSM, SB, FT, TREF, GE, PLY_POS, PLY_LEN):
        ft_str = map_failure_theory_int(ft)
        mids = MID[ply_pos:ply_pos+ply_len]
        thicknesses = T[ply_pos:ply_pos+ply_len]
        thetas = THETA[ply_pos:ply_pos+ply_len]
        souts = SOUT[ply_pos:ply_pos+ply_len]
        assert nplies > 0, nplies
        obj = geom_model.add_pcomp(
            pid, mids, thicknesses, thetas=thetas, souts=souts,
            nsm=nsm, sb=sb, ft=ft_str, tref=tref, ge=ge, lam=None, z0=z0, comment='')
        obj.validate()
        str(obj)

def read_pcompg(name: str, group: h5py._hl.dataset.Dataset, geom_model: BDF):
    asdf

def read_psolid(name: str, group: h5py._hl.dataset.Dataset, geom_model: BDF):
    """
    Dataset:
    attrs  : <Attributes of HDF5 object at 3135747352808>
    chunks : (340,)
    compression : 'gzip'
    compression_opts : 1
    dims   : <Dimensions of HDF5 object at 3135747352808>
    dtype  : dtype([('PID', '<i8'), ('MID1', '<i8'), ('T', '<f8'), ('MID2', '<i8'), ('BK', '<f8'), ('MID3', '<i8'), ('TS', '<f8'), ('NSM', '<f8'), ('Z1', '<f8'), ('Z2', '<f8'), ('MID4', '<i8'), ('DOMAIN_ID', '<i8')])
    external : None
    file   : <HDF5 file "6+element-nastran-sol103.h5" (mode r)>
    fillvalue : (0, 0, 0., 0, 0., 0, 0., 0., 0., 0., 0, 0)
    fletcher32 : False
    id     : <h5py.h5d.DatasetID object at 0x000002DA191B68E8>
    is_virtual : False
    maxshape : (None,)
    name   : '/NASTRAN/INPUT/PROPERTY/PSHELL'
    nbytes : 96
    ndim   : 1
    parent : <HDF5 group "/NASTRAN/INPUT/PROPERTY" (1 members)>
    ref    : <HDF5 object reference>
    regionref : <h5py._hl.base._RegionProxy object at 0x000002DA191E6D08>
    scaleoffset : None
    shape  : (1,)
    shuffle : True
    size   : 1
    """
    assert len(group.dtype.names) == 8, group.dtype.names
    PID = group['PID']
    MID = group['MID']
    CORDM = group['CORDM']
    IN = group['IN']
    STRESS = group['STRESS']
    ISOP = group['ISOP']
    FCTN = group['FCTN']
    DOMAIN_ID = group['DOMAIN_ID']
    for pid, mid, cordm, integ, stress, isop, fctn in zip(
        PID, MID, CORDM, IN, STRESS, ISOP, FCTN):

        integ = integ if integ != 0 else None
        stress = stress if stress != 0 else None
        isop = isop if isop != 0 else None
        obj = geom_model.add_psolid(
            pid, mid, cordm=cordm,
            integ=integ, stress=stress, isop=isop, fctn=fctn, comment='')
        obj.validate()
        str(obj)

def read_pbar(name: str, group: h5py._hl.dataset.Dataset, geom_model: BDF):
    assert len(group.dtype.names) == 20, group.dtype.names
    PID = group['PID']
    MID = group['MID']
    A = group['A']
    I1 = group['I1']
    I2 = group['I2']
    I12 = group['I12']
    J = group['J']
    NSM = group['NSM']

    FE = group['FE']
    C1 = group['C1']
    C2 = group['C2']
    D1 = group['D1']
    D2 = group['D2']
    E1 = group['E1']
    E2 = group['E2']
    F1 = group['F1']
    F2 = group['F2']

    K1 = group['K1']
    K2 = group['K2']

    DOMAIN_ID = group['DOMAIN_ID']
    for pid, mid, a, i1, i2, i12, j, nsm, fe, c1, c2, d1, d2, e1, e2, f1, f2, k1, k2 in zip(
        PID, MID, A, I1, I2, I12, J, NSM, FE, C1, C2, D1, D2, E1, E2, F1, F2, K1, K2):

        obj = geom_model.add_pbar(pid, mid, A=a, i1=i1, i2=i2, i12=i12, j=j, nsm=nsm,
                                  c1=c1, c2=c2, d1=d1, d2=d2, e1=e1, e2=e2, f1=f1, f2=f2,
                                  k1=k1, k2=k2, comment='')
        obj.validate()
        str(obj)

def read_pbeam(name: str, group: h5py._hl.dataset.Dataset, geom_model: BDF):
    assert len(group.dtype.names) == 38, group.dtype.names
    PID = group['PID']
    MID = group['MID']

    NSEGS = group['NSEGS']
    CCF = group['CCF']
    CWELD = group['CWELD']
    SO = group['SO']
    XXB = group['XXB']

    A = group['A']
    I1 = group['I1']
    I2 = group['I2']
    I12 = group['I12']
    J = group['J']
    NSM = group['NSM']

    C1 = group['C1']
    C2 = group['C2']
    D1 = group['D1']
    D2 = group['D2']
    E1 = group['E1']
    E2 = group['E2']
    F1 = group['F1']
    F2 = group['F2']

    K1 = group['K1']
    K2 = group['K2']
    S1 = group['S1']
    S2 = group['S2']
    NSIA = group['NSIA']
    NSIB = group['NSIB']
    CWA = group['CWA']
    CWB = group['CWB']

    M1A = group['M1A']
    M2A = group['M2A']
    M1B = group['M1B']
    M2B = group['M2B']

    N1A = group['N1A']
    N2A = group['N2A']
    N1B = group['N1B']
    N2B = group['N2B']

    DOMAIN_ID = group['DOMAIN_ID']
    for pid, mid, xxb, so, a, i1, i2, i12, j, nsm, c1, c2, d1, d2, e1, e2, f1, f2, k1, k2, s1, s2, nsia, nsib, cwa, cwb, m1a, m2a, m1b, m2b, n1a, n2a, n1b, n2b in zip(
        PID, MID, XXB, SO, A, I1, I2, I12, J, NSM, C1, C2, D1, D2, E1, E2, F1, F2, K1, K2, S1, S2, NSIA, NSIB, CWA, CWB, M1A, M2A, M1B, M2B, N1A, N2A, N1B, N2B):
        ixxb = (xxb > 0)
        ixxb[0] = True # save the first value (the 0 station)
        so_new = np.array(['NULL'] * len(xxb), dtype='|U4')
        so_new[np.where(so == 0)] = 'NO'
        so_new[np.where(so == 1)] = 'YES'
        so_list = so_new[ixxb].tolist()
        assert 'NULL' not in so_list, so_list
        #print(so_new)
        obj = geom_model.add_pbeam(
            pid, mid,
            xxb[ixxb].tolist(),
            so_list,
            a[ixxb].tolist(),
            i1[ixxb].tolist(),
            i2[ixxb].tolist(),
            i12[ixxb].tolist(),
            j[ixxb].tolist(),
            nsm[ixxb].tolist(),
            c1[ixxb].tolist(), c2[ixxb].tolist(),
            d1[ixxb].tolist(), d2[ixxb].tolist(),
            e1[ixxb].tolist(), e2[ixxb].tolist(),
            f1[ixxb].tolist(), f2[ixxb].tolist(),
            #xxb, so, a, i1, i2, i12, j, nsm=nsm,
            #c1=c1, c2=c2, d1=d1, d2=d2, e1=e1, e2=e2, f1=f1, f2=f2,
            k1=k1, k2=k2, s1=s1, s2=s2,
            nsia=nsia, nsib=nsib, cwa=cwa, cwb=cwb,
            m1a=m1a, m2a=m2a, m1b=m1b, m2b=m2b,
            n1a=n1a, n2a=n2a, n1b=n1b, n2b=n2b, comment='')
        #obj = geom_model.add_pbar(pid, mid, A=a, i1=i1, i2=i2, i12=i12, j=j, nsm=nsm,
                                  #c1=c1, c2=c2, d1=d1, d2=d2, e1=e1, e2=e2, f1=f1, f2=f2,
                                  #k1=k1, k2=k2, comment='')
        obj.validate()
        str(obj)

def read_pbarl(name: str, group: h5py._hl.group.Group, geom_model: BDF):
    IDENTITY = group['IDENTITY'] # ('PID', 'MID', 'GROUP', 'TYPE', 'INFO_POS', 'INFO_LEN', 'DOMAIN_ID')
    INFO = group['INFO'] # ('VALUE',)
    VALUE = INFO['VALUE']

    PID = IDENTITY['PID']
    MID = IDENTITY['MID']
    GROUP = IDENTITY['GROUP']
    TYPE = IDENTITY['TYPE']
    INFO_POS = IDENTITY['INFO_POS']
    INFO_LEN = IDENTITY['INFO_LEN']
    DOMAIN_ID = IDENTITY['DOMAIN_ID']

    properties = geom_model.properties
    for pid, mid, group, bar_type, ipos, ilen in zip(PID, MID, GROUP, TYPE, INFO_POS, INFO_LEN):
        if pid in properties and properties[pid].type == 'PBAR':
            del properties[pid]
        group_str = group.strip().decode('latin1')
        bar_type_str = bar_type.strip().decode('latin1')
        dim = VALUE[ipos:ipos+ilen-1].tolist()
        nsm = VALUE[ilen-1]
        obj = geom_model.add_pbarl(pid, mid, bar_type_str, dim, group=group_str, nsm=nsm, comment='')
        obj.validate()
        str(obj)

def read_pbeaml(name: str, group: h5py._hl.group.Group, geom_model: BDF):
    IDENTITY = group['IDENTITY'] # ('PID', 'MID', 'GROUP', 'TYPE', 'SECTION_POS', 'SECTION_LEN', 'DOMAIN_ID')
    DIMS = group['DIMS'] # ('DIM',)
    SECTION = group['SECTION'] # ('SO', 'RDIST', 'DIMS_POS', 'DIMS_LEN', 'NSM')

    DIM = DIMS['DIM']

    SO = SECTION['SO']
    RDIST = SECTION['RDIST']
    DIMS_POS = SECTION['DIMS_POS']
    DIMS_LEN = SECTION['DIMS_LEN']
    NSM = SECTION['NSM']

    PID = IDENTITY['PID']
    MID = IDENTITY['MID']
    GROUP = IDENTITY['GROUP']
    TYPE = IDENTITY['TYPE']
    SECTION_POS = IDENTITY['SECTION_POS']
    SECTION_LEN = IDENTITY['SECTION_LEN']
    DOMAIN_ID = IDENTITY['DOMAIN_ID']

    properties = geom_model.properties
    for pid, mid, group, beam_type, spos, slen in zip(PID, MID, GROUP, TYPE, SECTION_POS, SECTION_LEN):
        if pid in properties and properties[pid].type == 'PBEAM':
            del properties[pid]
        group_str = group.strip().decode('latin1')
        beam_type_str = beam_type.strip().decode('latin1')

        # section position/length
        si = slice(spos, spos+slen)
        xxb = RDIST[si]
        uxxb, ixxb = np.unique(xxb, return_index=True)
        xxb2 = xxb[ixxb]

        so = SO[si][ixxb]
        nsm = NSM[si][ixxb]
        dim_pos = DIMS_POS[si][ixxb]
        dim_len = DIMS_LEN[si][ixxb]

        # dim position/length
        dims = []
        for ipos, ilen in zip(dim_pos, dim_len):
            dimsi = DIM[ipos:ipos+ilen]
            dims.append(dimsi)

        #uxxb, ixxb = np.unique(xxb, return_index=True)
        #xxb2 = xxb[ixxb]
        #dims2 = [dims[i] for i in ixxb]
        #so2 = so[ixxb]
        #nsm2 = nsm[ixxb]
        obj = geom_model.add_pbeaml(pid, mid, beam_type_str,
                                    xxb2, dims, so=so, nsm=nsm,
                                    group=group_str, comment='')
        obj.validate()
        str(obj)

def read_paero1(name: str, group: h5py._hl.dataset.Dataset, geom_model: BDF):
    PID = group['PID']
    B1 = group['B1']
    B2 = group['B2']
    B3 = group['B3']
    B4 = group['B4']
    B5 = group['B5']
    B6 = group['B6']
    B = np.stack([B1, B2, B3, B4, B5, B6], axis=1)
    for pid, b in zip(PID, B):
        bi = b[np.where(b > 0)].tolist()
        if bi:
            obj = geom_model.add_paero1(pid, caero_body_ids=bi, comment='')
        else:
            obj = geom_model.add_paero1(pid, caero_body_ids=None, comment='')
        obj.validate()
        str(obj)

property_map = {
    #  mass
    'PMASS': passer,

    # 0d
    'PVISC': read_pvisc,
    'PBUSH': read_pbush,
    'PBUSH1D': passer,
    'PDAMP': read_pdamp,
    'PELAS': read_pelas,
    'PELAST': passer,

    #1d
    'PROD': read_prod,
    'PTUBE': read_ptube,
    'PBEND': passer,
    'PBAR': read_pbar,
    'PBEAM': read_pbeam,
    'PBARL': read_pbarl,
    'PBEAML': read_pbeaml,

    # 2D
    'PSHELL': read_pshell,
    'PSHEAR': read_pshear,
    'PCOMP': read_pcomp,
    'PCOMPG': read_pcompg,

    # 3d
    'PSOLID': read_psolid,

    # --------------------------
    # non-elements
    'PAERO1' : read_paero1,
    # --------------------------
    # thermal
    'PCONV' : passer,
}

