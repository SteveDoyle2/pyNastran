from __future__ import annotations
from typing import TYPE_CHECKING
import numpy as np
import h5py
from ..h5_utils import get_tree, passer
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf import BDF

def read_mpc(name: str, group: h5py._hl.group.Group, geom_model: BDF) -> None:
    #{'GCA': None, 'IDENTITY': None}
    GCA = group['GCA'] # ('G', 'C', 'A')
    GG = GCA['G']
    GC = GCA['C']
    GA = GCA['A']

    IDENTITY = group['IDENTITY']  # ('SID', 'G', 'C', 'A', 'GCA_POS', 'GCA_LEN', 'DOMAIN_ID')

    SID = IDENTITY['SID']
    IG = IDENTITY['G']
    IC = IDENTITY['C']
    IA = IDENTITY['A']
    GCA_POS = IDENTITY['GCA_POS']
    GCA_LEN = IDENTITY['GCA_LEN']
    DOMAIN_ID = IDENTITY['DOMAIN_ID']

    for sid, ig, ic, ia, gca_pos, gca_len in zip(SID, IG, IC, IA, GCA_POS, GCA_LEN):
        i_components = str(ic)
        g_nodes = GG[gca_pos:gca_pos+gca_len]
        g_components = GC[gca_pos:gca_pos+gca_len]
        g_coefficients = GA[gca_pos:gca_pos+gca_len]

        nodes = [ig] + g_nodes.tolist()
        components = [str(val) for val in [i_components] + g_components.tolist()]
        coefficients = [ia] + g_coefficients.tolist()
        obj = geom_model.add_mpc(sid, nodes, components, coefficients)
        obj.validate()
        str(obj)

def read_spc(name: str, group: h5py._hl.dataset.Dataset, geom_model: BDF) -> None:
    if geom_model.flags['constraint'] is False:
        return
    #('SID', 'G', 'C', 'D', 'DOMAIN_ID')
    SID = group['SID']
    G = group['G']
    C = group['C']
    D = group['D']
    DOMAIN_ID = group['DOMAIN_ID']
    for sid, g, c, d in zip(SID, G, C, D):
        nodes = [g]
        components = [str(c)]
        enforced = [d]
        obj = geom_model.add_spc(sid, nodes, components, enforced, comment='')
        obj.validate()
        str(obj)

def read_spc1(name: str, group: h5py._hl.group.Group, geom_model: BDF) -> None:
    if geom_model.flags['constraint'] is False:
        return
    #{'SPC1_G': {'G': None, 'IDENTITY': None}, 'SPC1_THRU': None}
    for name in list(group):
        groupi = group[name]
        if name == 'SPC1_G':
            G = groupi['G'] # ('ID',)
            IDENTITY = groupi['IDENTITY']  # ('SID', 'C', 'G_POS', 'G_LEN', 'DOMAIN_ID')

            ID = G['ID']

            SID = IDENTITY['SID']
            C = IDENTITY['C']
            G_POS = IDENTITY['G_POS']
            G_LEN = IDENTITY['G_LEN']
            DOMAIN_ID = IDENTITY['DOMAIN_ID']
            for sid, c, g_pos, g_len in zip(SID, C, G_POS, G_LEN):
                components = str(c)
                nodes = ID[g_pos:g_pos+g_len]
                obj = geom_model.add_spc1(sid, components, nodes, comment='')
                obj.validate()
        elif name == 'SPC1_THRU':
            # ('SID', 'C', 'FIRST', 'SECOND', 'DOMAIN_ID')
            SID = groupi['SID']
            C = groupi['C']
            FIRST = groupi['FIRST']
            SECOND = groupi['SECOND']
            for sid, c, first, second in zip(SID, C, FIRST, SECOND):
                components = str(c)
                nodes = [first, 'THRU', second]
                obj = geom_model.add_spc1(sid, components, nodes, comment='')
                obj.validate()
        else:  # pragma: no cover
            raise RuntimeError(name)

def read_mpcadd(name: str, group: h5py._hl.dataset.Dataset, geom_model: BDF) -> None:
    if geom_model.flags['constraint'] is False:
        return
    SETS = group['S'] # ('ID',)
    IDENTITY = group['IDENTITY']  # ('SID', 'C', 'S_POS', 'S_LEN', 'DOMAIN_ID')

    SET_ID = SETS['S']

    SID = IDENTITY['SID']
    S_POS = IDENTITY['S_POS']
    S_LEN = IDENTITY['S_LEN']
    DOMAIN_ID = IDENTITY['DOMAIN_ID']
    for sid, s_pos, s_len in zip(SID, S_POS, S_LEN):
        sets = SET_ID[s_pos:s_pos+s_len]
        obj = geom_model.add_spcadd(sid, sets, comment='')
        obj.validate()

def read_spcadd(name: str, group: h5py._hl.dataset.Dataset, geom_model: BDF) -> None:
    if geom_model.flags['constraint'] is False:
        return
    SETS = group['S'] # ('ID',)
    IDENTITY = group['IDENTITY']  # ('SID', 'C', 'S_POS', 'S_LEN', 'DOMAIN_ID')

    SET_ID = SETS['S']

    SID = IDENTITY['SID']
    S_POS = IDENTITY['S_POS']
    S_LEN = IDENTITY['S_LEN']
    DOMAIN_ID = IDENTITY['DOMAIN_ID']
    for sid, s_pos, s_len in zip(SID, S_POS, S_LEN):
        sets = SET_ID[s_pos:s_pos+s_len]
        obj = geom_model.add_spcadd(sid, sets, comment='')
        obj.validate()

def read_spline1(name: str, group: h5py._hl.dataset.Dataset, geom_model: BDF) -> None:
    if geom_model.flags['aero'] is False:
        return
    #('EID', 'CAERO', 'BOX1', 'BOX2', 'SETG', 'DZ', 'METHOD', 'USAGE', 'NELEM', 'MELEM', 'DOMAIN_ID')
    EID = group['EID']
    CAERO = group['CAERO']
    BOX1 = group['BOX1']
    BOX2 = group['BOX2']
    SETG = group['SETG']
    DZ = group['DZ']
    METHOD = group['METHOD']
    USAGE = group['USAGE']
    NELEM = group['NELEM']
    MELEM = group['MELEM']
    DOMAIN_ID = group['DOMAIN_ID']

    for eid, caero, box1, box2, setg, dz, method, usage, nelem, melem in zip(
        EID, CAERO, BOX1, BOX2, SETG, DZ, METHOD, USAGE, NELEM, MELEM):
                            method_str = method.strip().decode('latin1')
                            usage_str = usage.strip().decode('latin1')
                            assert isinstance(method_str, str), method_str
                            assert isinstance(usage_str, str), usage_str
                            obj = geom_model.add_spline1(
                                eid, caero, box1, box2, setg, dz=dz,
                                method=method_str, usage=usage_str,
                                nelements=nelem, melements=melem, comment='')
                            obj.validate()
                            str(obj)

def read_suport1(name: str, group: h5py._hl.group.Group, geom_model: BDF) -> None:
    if geom_model.flags['constraint'] is False:
        return
    COMPONENT = group['COMPONENT'] # ('ID', 'C')
    IDENTITY = group['IDENTITY'] # ('SID', 'COMPONENT_POS', 'COMPONENT_LEN', 'DOMAIN_ID')

    NODES = COMPONENT['ID']
    COMP = COMPONENT['C']

    SID = IDENTITY['SID']
    COMPONENT_POS = IDENTITY['COMPONENT_POS']
    COMPONENT_LEN = IDENTITY['COMPONENT_LEN']
    DOMAIN_ID = IDENTITY['DOMAIN_ID']
    for sid, ipos, ilen in zip(SID, COMPONENT_POS, COMPONENT_LEN):
        nodes = NODES[ipos:ipos+ilen]
        components = [str(val) for val in COMP[ipos:ipos+ilen]]
        obj = geom_model.add_suport1(sid, nodes, components, comment='')
        obj.validate()
        str(obj)

def read_trim(name: str, group: h5py._hl.group.Group, geom_model: BDF) -> None:
    TRIMS = group['TRIMS'] # ('LABEL', 'UX')
    IDENTITY = group['IDENTITY'] # ('ID', 'MACH', 'Q', 'AEQR', 'TRIMS_POS', 'TRIMS_LEN', 'DOMAIN_ID')

    LABEL = TRIMS['LABEL']
    UX = TRIMS['UX']

    SID = IDENTITY['ID']
    MACH = IDENTITY['MACH']
    Q = IDENTITY['Q']
    AEQR = IDENTITY['AEQR']
    TRIMS_POS = IDENTITY['TRIMS_POS']
    TRIMS_LEN = IDENTITY['TRIMS_LEN']
    DOMAIN_ID = IDENTITY['DOMAIN_ID']

    for sid, mach, q, aeqr, ipos, ilen in zip(SID, MACH, Q, AEQR, TRIMS_POS, TRIMS_LEN):
        labels = [val.strip().decode('latin1') for val in LABEL[ipos:ipos+ilen]]
        uxs = UX[ipos:ipos+ilen]
        obj = geom_model.add_trim(sid, mach, q, labels, uxs, aeqr=aeqr, trim_type=1, comment='')
        obj.validate()
        str(obj)

constraint_map = {
    'SPC': read_spc,
    'SPC1': read_spc1,
    'SPCD': passer,
    'SPCADD': read_spcadd,
    'MPC': read_mpc,
    'MPCADD': read_mpcadd,
    'SUPORT1': read_suport1,
    'TRIM': read_trim,
    'TEMP': passer,
    'TEMPD': passer,

    'SPLINE1': read_spline1,
    'SPLINE2': passer,
}
