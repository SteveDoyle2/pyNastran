from __future__ import annotations
from typing import Union, Callable, TYPE_CHECKING
import numpy as np
import h5py
from ..h5_utils import get_tree, passer

if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf import BDF
    GeomCallable = Callable[[h5py._hl.dataset.Dataset, BDF],  # inputs
                            None]  # output

def read_grid(name: str, group: h5py._hl.dataset.Dataset, geom_model: BDF) -> None:
    """
    Dataset:
    attrs  : <Attributes of HDF5 object at 2221931292520>
    chunks : (450,)
    compression : 'gzip'
    compression_opts : 1
    dims   : <Dimensions of HDF5 object at 2221931292520>
    dtype  : dtype([('ID', '<i8'), ('CP', '<i8'), ('X', '<f8', (3,)), ('CD', '<i8'), ('PS', '<i8'), ('SEID', '<i8'), ('DOMAIN_ID', '<i8')])
    external : None
    file   : <HDF5 file "6+element-nastran-sol103.h5" (mode r)>
    fillvalue : (0, 0, [0., 0., 0.], 0, 0, 0, 0)
    fletcher32 : False
    id     : <h5py.h5d.DatasetID object at 0x00000205556CE768>
    is_virtual : False
    maxshape : (None,)
    name   : '/NASTRAN/INPUT/NODE/GRID'
    nbytes : 3528
    ndim   : 1
    parent : <HDF5 group "/NASTRAN/INPUT/NODE" (1 members)>
    ref    : <HDF5 object reference>
    regionref : <h5py._hl.base._RegionProxy object at 0x00000205569A4D48>
    scaleoffset : None
    shape  : (49,)
    shuffle : True
    size   : 49
    """
    geom_model.card_count[name] = len(group['ID'])
    setattr(geom_model, name, group)
    return
    ID = group['ID']
    CP = group['CP']
    X = group['X']
    CD = group['CD']
    PS = group['PS']
    SEID = group['SEID']
    DOMAIN_ID = group['DOMAIN_ID']

    u_seid = np.unique(SEID)
    u_domain = np.unique(DOMAIN_ID)
    assert len(u_seid) == 1, u_seid
    assert len(u_domain) == 1, u_domain
    assert u_seid[0] == 0, u_seid
    assert u_domain[0] == 1, u_domain
    for nid, cp, xyz, cd, ps, seid in zip(ID, CP, X, CD, PS, SEID):
        if ps == 0:
            ps = ''
        geom_model.add_grid(nid, xyz, cp=cp, cd=cd, ps=ps, seid=seid, comment='')

def read_spoint(name: str, group: h5py._hl.dataset.Dataset, geom_model: BDF) -> None:
    ID = group['ID']
    DOMAIN_ID = group['DOMAIN_ID']

    u_domain = np.unique(DOMAIN_ID)
    assert len(u_domain) == 1, u_domain
    assert u_domain[0] == 1, u_domain
    geom_model.add_spoint(ID)

def read_epoint(name: str, group: h5py._hl.dataset.Dataset, geom_model: BDF) -> None:
    ID = group['ID']
    DOMAIN_ID = group['DOMAIN_ID']

    u_domain = np.unique(DOMAIN_ID)
    assert len(u_domain) == 1, u_domain
    assert u_domain[0] == 1, u_domain
    geom_model.add_spoint(ID)

def read_desvar(name: str, group: h5py._hl.dataset.Dataset, geom_model: BDF) -> None:
    #('ID', 'LABEL', 'XINIT', 'XLB', 'XUB', 'DELX', 'DVID', 'DOMAIN_ID')
    ID = group['ID']
    LABEL = group['LABEL']
    XINIT = group['XINIT']
    XLB = group['XLB']
    XUB = group['XUB']
    DELX = group['DELX']
    DVID = group['DVID']
    DOMAIN_ID = group['DOMAIN_ID']
    for desvar_id, label, xinit, xlb, xub, delx, ddval in zip(ID, LABEL, XINIT, XLB, XUB, DELX, DVID):
        label_str = label.strip().decode('latin1')
        obj = geom_model.add_desvar(desvar_id, label_str, xinit, xlb=xlb, xub=xub,
                                    delx=delx, ddval=ddval, comment='')
        obj.validate()
        str(obj)

def read_dconstr(name: str, group: h5py._hl.dataset.Dataset, geom_model: BDF) -> None:
    #('DCID', 'RID', 'LALLOW', 'UALLOW', 'LOWFQ', 'HIGHFQ', 'DTYPE', 'DOMAIN_ID')
    DCID = group['DCID']
    RID = group['RID']
    LALLOW = group['LALLOW']
    UALLOW = group['UALLOW']
    LOWFQ = group['LOWFQ']
    HIGHFQ = group['HIGHFQ']
    DTYPE = group['DTYPE']
    DOMAIN_ID = group['DOMAIN_ID']
    for oid, dresp_id, lallow, uallow, lowfq, highfq, dtype in zip(
        DCID, RID, LALLOW, UALLOW, LOWFQ, HIGHFQ, DTYPE):
        obj = geom_model.add_dconstr(oid, dresp_id, lid=lallow, uid=uallow,
                                     lowfq=lowfq, highfq=highfq, comment='')
        obj.validate()
        str(obj)

def read_cord2c(name: str, group: h5py._hl.dataset.Dataset, geom_model: BDF) -> None:
    _read_cord2(name, group, geom_model.add_cord2c)

def read_cord2r(name: str, group: h5py._hl.dataset.Dataset, geom_model: BDF) -> None:
    _read_cord2(name, group, geom_model.add_cord2r)

def read_cord2s(name: str, group: h5py._hl.dataset.Dataset, geom_model: BDF) -> None:
    _read_cord2(name, group, geom_model.add_cord2s)

def _read_cord2(name: str, group: h5py._hl.dataset.Dataset,
                add_func: Callable) -> None:
    assert len(group.dtype.names) == 12, group.dtype.names
    CID = group['CID']
    RID = group['RID']
    A1 = group['A1']
    A2 = group['A2']
    A3 = group['A3']
    B1 = group['B1']
    B2 = group['B2']
    B3 = group['B3']
    C1 = group['C1']
    C2 = group['C2']
    C3 = group['C3']
    DOMAIN_ID = group['DOMAIN_ID']
    for cid, rid, a1, a2, a3, b1, b2, b3, c1, c2, c3 in zip(CID, RID, A1, A2, A3, B1, B2, B3, C1, C2, C3):
        origin = np.array([a1, a2, a3])
        zaxis = np.array([b1, b2, b3])
        xzplane = np.array([c1, c2, c3])
        obj = add_func(cid, origin, zaxis, xzplane, rid=rid, setup=True, comment='')
        obj.validate()
        str(obj)

def read_mat1(name: str, group: h5py._hl.dataset.Dataset, geom_model: BDF) -> None:
    """
    Dataset:
    attrs  : <Attributes of HDF5 object at 2553977821512>
    chunks : (310,)
    compression : 'gzip'
    compression_opts : 1
    dims   : <Dimensions of HDF5 object at 2553977821512>
    dtype  : dtype([('MID', '<i8'), ('E', '<f8'), ('G', '<f8'), ('NU', '<f8'), ('RHO', '<f8'), ('A', '<f8'), ('TREF', '<f8'), ('GE', '<f8'), ('ST', '<f8'), ('SC', '<f8'), ('SS', '<f8'), ('MCSID', '<i8'), ('DOMAIN_ID', '<i8')])
    external : None
    file   : <HDF5 file "6+element-nastran-sol103.h5" (mode r)>
    fillvalue : (0, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0, 0)
    fletcher32 : False
    id     : <h5py.h5d.DatasetID object at 0x00000252A4F0D948>
    is_virtual : False
    maxshape : (None,)
    name   : '/NASTRAN/INPUT/MATERIAL/MAT1'
    nbytes : 104
    ndim   : 1
    parent : <HDF5 group "/NASTRAN/INPUT/MATERIAL" (1 members)>
    ref    : <HDF5 object reference>
    regionref : <h5py._hl.base._RegionProxy object at 0x00000252A60614C8>
    scaleoffset : None
    shape  : (1,)
    shuffle : True
    size   : 1
    """
    #'MID', 'E', 'G', 'NU', 'RHO', 'A', 'TREF', 'GE', 'ST', 'SC', 'SS', 'MCSID', 'DOMAIN_ID',
    assert len(group.dtype.names) == 13, group.dtype.names
    MID = group['MID']
    E = group['E']
    G = group['G']
    NU = group['NU']
    RHO = group['RHO']
    A = group['A']
    TREF = group['TREF']
    GE = group['GE']
    ST = group['ST']
    SC = group['SC']
    SS = group['SS']
    MCSID = group['MCSID']
    DOMAIN_ID = group['DOMAIN_ID']
    for mid, e, g, nu, rho, a, tref, ge, st, sc, ss, mcsid in zip(MID, E, G, NU, RHO, A, TREF, GE, ST, SC, SS, MCSID):
        #if mcid == -1:
            #theta_mcid = theta
        #else:
            #asdf
        #assert tflag == 0, tflag
        #t1, t2, t3, t4 = [ti if ti != -1.0 else None
                          #for ti in t]
        assert mcsid == 0, mcsid
        obj = geom_model.add_mat1(mid, e, g, nu, rho=rho, a=a, tref=tref, ge=ge,
                                  St=st, Sc=sc, Ss=ss, mcsid=mcsid, comment='')
        obj.validate()
        str(obj)

def read_mat2(name: str, group: h5py._hl.dataset.Dataset, geom_model: BDF) -> None:
    assert len(group.dtype.names) == 18, group.dtype.names
    MID = group['MID']
    G11 = group['G11']
    G12 = group['G12']
    G13 = group['G13']
    G22 = group['G22']
    G23 = group['G23']
    G33 = group['G33']
    RHO = group['RHO']
    A1 = group['A1']
    A2 = group['A2']
    A3 = group['A12']
    TREF = group['TREF']
    GE = group['GE']
    ST = group['ST']
    SC = group['SC']
    SS = group['SS']
    MCSID = group['MCSID']
    DOMAIN_ID = group['DOMAIN_ID']
    for mid, g11, g12, g13, g22, g23, g33, rho, a1, a2, a3, tref, ge, st, sc, ss, mcsid in zip(
        MID, G11, G12, G13, G22, G23, G33, RHO, A1, A2, A3, TREF, GE, ST, SC, SS, MCSID):
        if mid > 100_000_000:
            continue
        assert mcsid == 0, mcsid
        obj = geom_model.add_mat2(mid, g11, g12, g13, g22, g23, g33,
                                  rho=rho, a1=a1, a2=a2, a3=a3,
                                  tref=tref, ge=ge,
                                  St=st, Sc=sc, Ss=ss, mcsid=None, comment='')
        obj.validate()
        str(obj)

def read_mat8(name: str, group: h5py._hl.dataset.Dataset, geom_model: BDF) -> None:
    #assert len(group.dtype.names) == 18, group.dtype.names
    #('MID', 'E1', 'E2', 'NU12', 'G12', 'G1Z', 'G2Z', 'RHO', 'A1', 'A2', 'TREF', 'XT', 'XC', 'YT', 'YC',
    # 'S', 'GE', 'F12', 'STRN', 'DOMAIN_ID')

    MID = group['MID']
    E1 = group['E1']
    E2 = group['E2']
    NU12 = group['NU12']
    G12 = group['G12']
    G1Z = group['G1Z']
    G2Z = group['G2Z']
    RHO = group['RHO']
    A1 = group['A1']
    A2 = group['A2']
    TREF = group['TREF']
    XT = group['XT']
    XC = group['XC']
    YT = group['YT']
    YC = group['YC']
    S = group['S']
    GE = group['GE']
    F12 = group['F12']
    STRN = group['STRN']
    DOMAIN_ID = group['DOMAIN_ID']
    for mid, e11, e22, nu12, g12, g1z, g2z, rho, a1, a2, tref, xt, xc, yt, yc, s, ge, f12, strn in zip(
        MID, E1,  E2,  NU12, G12, G1Z, G2Z, RHO, A1, A2, TREF, XT, XC, YT, YC, S, GE, F12, STRN):
        obj = geom_model.add_mat8(mid, e11, e22, nu12, g12=g12, g1z=g1z, g2z=g2z,
                                  rho=rho, a1=a1, a2=a2, tref=tref,
                                  Xt=xt, Xc=xc, Yt=yt, Yc=yc, S=s, ge=ge, F12=f12, strn=strn, comment='')
        obj.validate()
        str(obj)

def read_dvprel1(name: str, group: h5py._hl.group.Group, geom_model: BDF) -> None:
    # TODO: group, not dataset
    # {'ATTI': None, 'IDENTITY': None}

    ('ID', 'TYPE', 'PID', 'FID', 'PMIN', 'PMAX', 'C0', 'PNAME', 'START', 'LEN', 'DOMAIN_ID')
    IDENTITY = group.get('IDENTITY')
    RELATION = group.get('RELATION') # ('DVID', 'COEF')
    DVID = RELATION['DVID']
    COEF = RELATION['COEF']

    ID = IDENTITY['ID']
    TYPE = IDENTITY['TYPE']
    PID = IDENTITY['PID']
    FID = IDENTITY['FID']
    PMIN = IDENTITY['PMIN']
    PMAX = IDENTITY['PMAX']
    C0 = IDENTITY['C0']
    PNAME = IDENTITY['PNAME']

    START = IDENTITY['START']
    LEN = IDENTITY['LEN']
    DOMAIN_ID = IDENTITY['DOMAIN_ID']
    for oid, typei, pid, fid, pmin, pmax, c0, pname, ipos, ilen in zip(
        ID,  TYPE,  PID, FID, PMIN, PMAX, C0, PNAME, START, LEN):
        pname_str = pname.strip().decode('latin1')
        prop_type = typei.strip().decode('latin1')
        dvids = DVID[ipos:ipos+ilen]
        coeffs = COEF[ipos:ipos+ilen]
        if pname_str:
            pname_fid = pname_str
        elif fid != 0:
            pname_fid = fid
        else:
            out = (pname_str, fid)
            raise RuntimeError(out)

        obj = geom_model.add_dvprel1(oid, prop_type, pid, pname_fid, dvids, coeffs,
                                     p_min=pmin, p_max=pmax, c0=c0, validate=True, comment='')
        obj.validate()
        str(obj)

def read_dresp1(name: str, group: h5py._hl.group.Group, geom_model: BDF) -> None:
    # TODO: group, not dataset
    # {'ATTI': None, 'IDENTITY': None}
    #('ID', 'LABEL', 'RTYPE', 'PTYPE', 'RPSID', 'NTUSED', 'AFPMID', 'REGION', 'ATTA', 'ATTBI', 'ATTBR', 'ATTI_LEN', 'ATTI_POS', 'DOMAIN_ID')
    IDENTITY = group.get('IDENTITY')
    ID = IDENTITY['ID']
    LABEL = IDENTITY['LABEL']
    RTYPE = IDENTITY['RTYPE']
    PTYPE = IDENTITY['PTYPE']
    RPSID = IDENTITY['RPSID']
    NTUSED = IDENTITY['NTUSED']
    AFPMID = IDENTITY['AFPMID']
    REGION = IDENTITY['REGION']
    ATTA = IDENTITY['ATTA']
    ATTBI = IDENTITY['ATTBI']
    ATTBR = IDENTITY['ATTBR']
    ATTI_LEN = IDENTITY['ATTI_LEN']
    ATTI_POS = IDENTITY['ATTI_POS']
    DOMAIN_ID = IDENTITY['DOMAIN_ID']

    #('ATTI',)
    ATTI_group = group.get('ATTI')
    ATTI = ATTI_group['ATTI']

    for idi, label, rtype, ptype, rpsid, ntused, afpmid, region, atta, attbi, attbr, ilen, ipos in zip(
        ID,  LABEL, RTYPE, PTYPE, RPSID, NTUSED, AFPMID, REGION, ATTA, ATTBI, ATTBR, ATTI_LEN, ATTI_POS):
        # rpsid
        # ntused
        # attbi/r

        label_str = label.strip().decode('latin1')
        property_type = ptype.strip().decode('latin1')
        response_type = rtype.strip().decode('latin1')
        label_str = label.strip().decode('latin1')
        if property_type == '':
            property_type = None
        atta, attb = get_attb_from_atti(response_type, atta, attbi, attbr)
        atti = ATTI[ipos:ipos+ilen].tolist()
        obj = geom_model.add_dresp1(idi, label_str, response_type, property_type, region,
                                    atta, attb, atti, validate=True, comment='')
        assert afpmid == 0, afpmid
        obj.validate()
        str(obj)

def get_attb_from_atti(response_type: str, atta: int, attb_integer: int, attb_real: float) -> Union[int, float, None]:
    # 0,1
    if response_type in {'DISP', 'STRAIN', 'STRESS', 'ESE', 'FORCE', 'CSTRESS', 'CSTRAIN', 'CFAILURE', 'TOTSE',
                         'RMSDISP', 'RMSVELO', 'RMSACCL', 'STABDER'}:
        attb = attb_integer
        if attb == 0:
            attb = None
    elif response_type in {}:
        attb = attb_real
    elif response_type in ['VOLUME', 'WEIGHT']:
        if attb_integer == -9999:
            attb = None
        else:
            raise RuntimeError(f'response_type={response_type} attb_integer={attb_integer} attb_real={attb_real}; what type?')
    elif response_type in ['FLUTTER']:
        assert atta == 0, atta
        assert attb_integer == 0, attb_integer
        atta = None
        attb = None
    else:
        raise RuntimeError(f'response_type={response_type} attb_integer={attb_integer} attb_real={attb_real}; what type?')
    return atta, attb

def read_conm1(name: str, group: h5py._hl.dataset.Dataset, geom_model: BDF) -> None:
    #('EID', 'G', 'CID', 'M1', 'M2', 'M3', 'M4', 'M5', 'M6', 'DOMAIN_ID')
    EID = group['EID']
    G = group['G']
    CID = group['CID']
    M1 = group['M1']
    M2 = group['M2']
    M3 = group['M3']
    M4 = group['M4']
    M5 = group['M5']
    M6 = group['M6']

    # [M] = [M11 M21 M31 M41 M51 M61]
    #       [    M22 M32 M42 M52 M62]
    #       [        M33 M43 M53 M63]
    #       [            M44 M54 M64]
    #       [    Sym         M55 M65]
    #       [                    M66]

    DOMAIN_ID = group['DOMAIN_ID']
    for eid, nid, cid, m1, m2, m3, m4, m5, m6 in zip(EID, G, CID, M1, M2, M3, M4, M5, M6):
        mass_matrix = np.zeros((6, 6), dtype=m1.dtype)
        mass_matrix[0, 0] = m1
        mass_matrix[1, :2] = mass_matrix[:2, 1] = m2
        mass_matrix[2, :3] = mass_matrix[:3, 2] = m3
        mass_matrix[3, :4] = mass_matrix[:4, 3] = m4
        mass_matrix[4, :5] = mass_matrix[:5, 4] = m5
        mass_matrix[5, :]  = mass_matrix[:, 5]  = m6
        obj = geom_model.add_conm1(eid, nid, mass_matrix, cid=cid, comment='')
        obj.validate()
        str(obj)

def load_geometry_block(node: h5py._hl.group.Group,
                        name_map: dict[str, GeomCallable],
                        geom_model: BDF) -> None:
    assert node is not None, node
    assert isinstance(node, h5py._hl.group.Group)
    for name in list(node):
        group = node.get(name)
        assert group is not None, name
        if name in name_map:
            func = name_map[name]
            func(name, group, geom_model)
        else:
            raise RuntimeError(name)

coord_map = {
    'CORD2C': read_cord2c,
    'CORD2R': read_cord2r,
    'CORD2S': read_cord2s,
    'TRANSFORMATION': passer,
}
node_map = {
    'GRID': read_grid,
    'SPOINT': read_spoint,
    'EPOINT': read_epoint,
}

material_map = {
    'MAT1': read_mat1,
    'MAT2': read_mat2,
    'MAT4': passer,
    #'MAT5': passer,
    'MAT8': read_mat8,
    'NLMOPTS': passer,
}

design_map = {
    'DOPTPRM': passer,
    'DTABLE': passer,
    'DCONADD': passer,
    'DVPREL1': read_dvprel1,
    'DVPREL2': passer,
    'DRESP1': read_dresp1,
    'DRESP2': passer,
    #'DRESP3': passer,
    'DESVAR': read_desvar,
    'DCONSTR': read_dconstr,
}
dynamic_map = {
    'EIGB': passer,
    'EIGC': passer,
    'EIGR': passer,
    'EIGRL': passer,
    'FREQ': passer,
    'FREQ1': passer,
    'FREQ2': passer,
    'FREQ3': passer,
    'FREQ4': passer,
    'FREQ5': passer,
}

matrix_map = {
    'CONM1': read_conm1,
}
partition_map = {
    'SET1': passer,
    'AELIST': passer,
    'AESURF': passer,
}
