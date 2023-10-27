from __future__ import annotations
from typing import TYPE_CHECKING
import numpy as np

from .utils import cast_encoding_strip, get_group_name
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.dev.bdf_vectorized3.bdf_interface.bdf_attributes import (
        PELAS, PROD, PTUBE, PVISC, PBUSH1D, PBAR, PSHEAR, PSHELL)
    from pyNastran.dev.bdf_vectorized3.cards.elements.beam import PBEAM
    from pyNastran.dev.bdf_vectorized3.cards.elements.solid import PSOLID
    from pyNastran.dev.bdf_vectorized3.bdf import BDF
    from tables import Group


def load_h5_property(model: BDF, input_group: Group):
    for h5_element in input_group._f_iter_nodes():
        if h5_element._c_classid == 'GROUP':
            class_name = get_group_name(h5_element)
            if class_name == 'PBARL':
                _load_h5_pbarl(model, h5_element)
            elif class_name == 'PBEAML':
                _load_h5_pbeaml(model, h5_element)
            elif class_name == 'PCOMP':
                _load_h5_pcomp(model, h5_element)
            else:
                print(f'skipping {class_name}')
                #print(class_name)
                continue
            #print(f'loading {class_name}')
            continue

        name = h5_element.name
        data = h5_element.read()
        property_id = data['PID']
        domain_id = data['DOMAIN_ID']
        nproperties = len(property_id)
        if name == 'PSHELL':
            prop = model.pshell
            _load_h5_pshell(prop, property_id, data, domain_id)
        elif name == 'PBAR':
            prop = model.pbar
            _load_h5_pbar(prop, property_id, data, domain_id)
        elif name == 'PBEAM':
            prop = model.pbeam
            _load_h5_pbeam(prop, property_id, data, domain_id)
        elif name == 'PSHEAR':
            prop = model.pshear
            _load_h5_pshear(prop, property_id, data, domain_id)

        elif name == 'PSOLID':
            prop = model.psolid
            _load_h5_psolid(prop, property_id, data, domain_id)

        elif name == 'PTUBE':
            prop = model.ptube
            _load_h5_ptube(prop, property_id, data, domain_id)
        elif name == 'PELAS':
            prop = model.pelas
            _load_h5_pelas(prop, property_id, data, domain_id)
        elif name == 'PDAMP':
            prop = model.pdamp
            b = data['B']
            prop._save(property_id, b)
            prop.domain_id = domain_id

        elif name == 'PROD':
            prop = model.prod
            _load_h5_prod(prop, property_id, data, domain_id)
        elif name == 'PVISC':
            prop = model.pvisc
            _load_h5_visc(prop, property_id, data, domain_id)
        elif name == 'PBUSH1D':
            prop = model.pbush1d
            _load_h5_pbush1d(prop, property_id, data, domain_id)
        elif name == 'PELAST':
            #('PID', 'TKID', 'TGEID', 'TKNID', 'DOMAIN_ID')
            prop = model.pelast
            table_k = data['TKID']
            table_ge = data['TGEID']
            table_k_nonlinear = data['TKNID']
            prop._save(property_id, table_k, table_ge, table_k_nonlinear)
        elif name in {'PAERO1', 'PBEND', 'PCONV'}:
            #print(f'skipping {name}')
            model.log.warning(f'skipping {name}')
            continue
        else:
            raise NotImplementedError(name)
        #print(f'loading {name}')
        assert prop is not None, name
        prop.property_id = property_id
        prop.n = nproperties
        #prop.write()
        x = 1
    model.pshell._filter_pcomp_pshells()
    x = 2

def _load_h5_pelas(prop: PELAS, property_id, data, domain_id):
    k = data['K']
    ge = data['GE']
    s = data['S']
    prop._save(property_id, k, ge, s)
    prop.domain_id = domain_id

def _load_h5_prod(prop: PROD, property_id, data, domain_id):
    material_id = data['MID']
    A = data['A']
    J = data['J']
    c = data['C']
    nsm = data['NSM']
    prop._save(property_id, material_id, A, J, c, nsm)
    prop.domain_id = domain_id

def _load_h5_visc(prop: PVISC, property_id, data, domain_id):
    #('PID', 'CE', 'CR', 'DOMAIN_ID')
    property_id = data['PID']
    ce = data['CE']
    cr = data['CR']
    prop._save(property_id, cr, ce)
    prop.domain_id = domain_id

def _load_h5_ptube(prop: PTUBE, property_id, data, domain_id) -> None:
    # dtype([('PID', '<i8'), ('MID', '<i8'), ('OD', '<f8'), ('T', '<f8'),
    #       ('NSM', '<f8'), ('DOMAIN_ID', '<i8')])
    #prop = model.ptube
    nsm = data['NSM']
    material_id = data['MID']
    t = data['T']
    diameter = data['OD']
    prop._save(property_id, material_id, diameter, t, nsm)
    prop.domain_id = domain_id

def _load_h5_pshear(prop: PSHEAR, property_id, data, domain_id) -> None:
    # dtype([('PID', '<i8'), ('MID', '<i8'), ('T', '<f8'), ('NSM', '<f8'),
    #       ('F1', '<f8'), ('F2', '<f8'), ('DOMAIN_ID', '<i8')])
    #prop = model.pshear
    material_id = data['MID']
    t = data['T']
    f1 = data['F1']
    f2 = data['F2']
    nsm = data['NSM']
    prop._save(property_id, material_id, t, nsm, f1, f2)
    prop.domain_id = domain_id
    #return prop

def _load_h5_pshell(prop: PSHELL,
                    property_id: np.ndarray,
                    data: np.ndarray,
                    domain_id: np.ndarray) -> None:
    #prop = model.pshell
    #dtype([('PID', '<i8'), ('MID1', '<i8'), ('T', '<f8'), ('MID2', '<i8'),
           #('BK', '<f8'), ('MID3', '<i8'), ('TS', '<f8'), ('NSM', '<f8'),
           #('Z1', '<f8'), ('Z2', '<f8'), ('MID4', '<i8'), ('DOMAIN_ID', '<i8')])
    material_id = np.stack([
        data['MID1'], data['MID2'],
        data['MID3'], data['MID4']], axis=1)
    z = np.stack([data['Z1'], data['Z2']], axis=1)
    nsm = data['NSM']
    twelveIt3 = data['BK']
    t = data['T']
    tst = data['TS']
    prop._save(property_id, material_id, t, twelveIt3, tst, nsm, z)
    prop.domain_id = domain_id
    #return prop

def _load_h5_psolid(prop: PSOLID, property_id, data, domain_id) -> None:
    # dtype({'names': ['PID', 'MID', 'CORDM', 'IN', 'STRESS',
    #                  'ISOP', 'FCTN', 'DOMAIN_ID'],
    #        'formats': ['<i8', '<i8', '<i8', '<i8', '<i8', '<i8', 'S4', '<i8'],
    #        'offsets': [0, 8, 16, 24, 32, 40, 48, 56], 'itemsize': 64})
    #prop = model.psolid
    material_id = data['MID']
    coord_id = data['CORDM']
    integ = data['IN']
    stress = data['STRESS']
    isop = data['ISOP']
    fctn = data['FCTN']
    prop._save(property_id, material_id, coord_id, integ, stress, isop, fctn)
    prop.domain_id = domain_id
    #return prop

def _load_h5_pbar(prop: PBAR,
                  property_id: np.ndarray,
                  data: np.ndarray,
                  domain_id: np.ndarray) -> None:
    #prop = model.pbar
    #nproperties = len(property_id)
    # dtype([('PID', '<i8'), ('MID', '<i8'), ('A', '<f8'),
    #        ('I1', '<f8'), ('I2', '<f8'), ('J', '<f8'), ('NSM', '<f8'),
    #        ('FE', '<f8'),
    #        ('C1', '<f8'), ('C2', '<f8'),
    #        ('D1', '<f8'), ('D2', '<f8'), ('E1', '<f8'), ('E2', '<f8'),
    #        ('F1', '<f8'), ('F2', '<f8'), ('K1', '<f8'), ('K2', '<f8'),
    #        ('I12', '<f8'), ('DOMAIN_ID', '<i8')])
    FE = data['FE']
    material_id = data['MID']
    c = np.stack([data['C1'], data['C2']], axis=1)
    d = np.stack([data['D1'], data['D2']], axis=1)
    e = np.stack([data['E1'], data['E2']], axis=1)
    f = np.stack([data['F1'], data['F2']], axis=1)
    k = np.stack([data['K1'], data['K2']], axis=1)
    I = np.stack([data['I1'], data['I2'], data['I12']], axis=1)
    A = data['A']
    J = data['J']
    nsm = data['NSM']
    prop._save(property_id, material_id, A, J, c, d, e, f, I, k, nsm)
    prop.domain_id = domain_id
    #return prop

def _load_h5_pbeam(prop: PBEAM,
                   property_id: np.ndarray,
                   data: np.ndarray,
                   domain_id: np.ndarray) -> PBEAM:
    #prop: PBEAM = model.pbeam
    nproperties = len(property_id)
    #('PID', 'MID', 'NSEGS', 'CCF', 'CWELD',
     #'SO', 'XXB', 'A', 'I1', 'I2', 'I12', 'J', 'NSM',
     #'C1', 'C2', 'D1', 'D2', 'E1', 'E2', 'F1', 'F2',
     #'K1', 'K2', 'S1', 'S2',
     #'NSIA', 'NSIB',
     #'CWA', 'CWB',
     #'M1A', 'M2A', 'M1B', 'M2B',
     #'N1A', 'N2A', 'N1B', 'N2B',
     #'DOMAIN_ID')

    #NSEGS I Number of segments (or intermediate stations)
    #CCF I Constant cross-section flag: 1=yes and 0=no
    nintermediate_segments = data['NSEGS']  # [1]
    #nstation = nintermediate_segments + 2
    nstation = np.zeros(nproperties, dtype=nintermediate_segments.dtype)

    # overly complicated way to filter cross-sections with no x/xb definition
    # also has to filter intermediate stations that have x/xb=1.0
    xxb_full = data['XXB']
    if nintermediate_segments.max() > 0:
        ikeep_ = []
        jkeep_ = []
        for i, nsegs in enumerate(nintermediate_segments):
            jkeepi_ = [0]
            val_10 = 10
            for iseg in range(1, nsegs+1):
                if xxb_full[i, iseg] == 1.0:
                    val_10 = iseg
                    continue
                jkeepi_.append(iseg)
            jkeepi_.append(val_10)

            nstationi = len(jkeepi_)
            ikeepi_ = [i] * nstationi
            nstation[i] = nstationi
            ikeep_.extend(ikeepi_)
            jkeep_.extend(jkeepi_)
        ikeep = np.array(ikeep_, dtype='int32')
        jkeep = np.array(jkeep_, dtype='int32')
        # filter xxb stations that have 1.0
        #jxxb = np.hstack([np.arange(1, nsegmenti+1) for nsegmenti in nintermediate_segments])
        #ixxb = np.repeat(np.arange(nproperties, dtype='int32'), nintermediate_segments)
        #assert len(ixxb) == len(jxxb)
        #xxb_intermediate = xxb_full[ixxb, jxxb]
        #ifilter = np.where(xxb_intermediate != 1.0)

        #jkeep = np.hstack([np.hstack([np.arange(nsegmenti+1), [10]]) for nsegmenti in nintermediate_segments])
        #ikeep = np.repeat(np.arange(nproperties, dtype='int32'), nstation)
        #jkeep = jkeep[ifilter]
        #ikeep = ikeep[ifilter]
    else:
        jkeep = np.array([0, 10] * nproperties, dtype='int32')
        ikeep = np.repeat(np.arange(nproperties, dtype='int32'), 2)
    assert len(jkeep) == len(ikeep)
    assert len(jkeep) == len(ikeep)
    section_flag = data['CCF']  # [2]
    cweld = data['CWELD']  # [0]
    #FE = data['FE']
    material_id = data['MID']
    #c = np.stack([data['C1'], data['C2']], axis=1)
    #d = np.stack([data['D1'], data['D2']], axis=1)
    #e = np.stack([data['E1'], data['E2']], axis=1)
    #f = np.stack([data['F1'], data['F2']], axis=1)
    #k = np.stack([data['K1'], data['K2']], axis=1)
    #s = np.stack([data['S1'], data['S2']], axis=1)
    #I = np.stack([data['I1'], data['I2'], data['I12']], axis=1)

    c1 = data['C1'][ikeep, jkeep]
    d1 = data['D1'][ikeep, jkeep]
    e1 = data['E1'][ikeep, jkeep]
    f1 = data['F1'][ikeep, jkeep]
    I1 = data['I1'][ikeep, jkeep]

    c2 = data['C2'][ikeep, jkeep]
    d2 = data['D2'][ikeep, jkeep]
    e2 = data['E2'][ikeep, jkeep]
    f2 = data['F2'][ikeep, jkeep]
    I2 = data['I2'][ikeep, jkeep]
    I12 = data['I12'][ikeep, jkeep]

    s1 = data['S1']
    s2 = data['S2']

    k1 = data['K1']
    k2 = data['K2']

    nsia = data['NSIA']
    nsib = data['NSIB']

    cwa = data['CWA']
    cwb = data['CWB']
    m1a = data['M1A']
    m1b = data['M1B']
    m2a = data['M2A']
    m2b = data['M2B']

    n1a = data['N1A']
    n1b = data['N1B']
    n2a = data['N2A']
    n2b = data['N2B']

    #nsi = np.stack([data['NSIA'], data['NSIB']], axis=1)
    #m = np.stack([data['M1A'], data['M1B']], axis=1)
    #cw = np.stack([data['CWA'], data['CWB']], axis=1)
    A = data['A'][ikeep, jkeep]
    J = data['J'][ikeep, jkeep]
    nsm = data['NSM'][ikeep, jkeep]
    so = data['SO'][ikeep, jkeep]
    xxb = xxb_full[ikeep, jkeep]
    nsm = data['NSM'][ikeep, jkeep]

    str_map = {
        1.0 : 'YES',
    }
    nstations_ = len(so)
    so_str = map_table_string_by_float_map(so, nstations_, str_map)

    prop._save(property_id, material_id,
              nstation, xxb, so_str,
              A, J, I1, I2, I12, nsm,
              c1, c2, d1, d2, e1, e2, f1, f2,
              s1, s2, k1, k2,
              nsia, nsib, cwa, cwb,
              m1a, m2a, m1b, m2b,
              n1a, n2a, n1b, n2b)
    ##property_id, material_id,
               ##nstation, xxb, so_str,
               ##A, J, I1, I2, I12, nsm,
               ##c1, c2, d1, d2, e1, e2, f1, f2,
               ##s1, s2, k1, k2)
    prop.write()
    prop.domain_id = domain_id
    #return prop

def map_table_string_by_float_map(float_array: np.ndarray,
                                  n: int,
                                  str_map: dict[float, str]) -> np.ndarray:
    str_array = np.zeros(n, dtype='|U4')
    for soi in np.unique(float_array):
        so_stri = str_map[soi]
        iso = np.where(float_array == soi)[0]
        str_array[iso] = so_stri
    return str_array

def _load_h5_pbush1d(prop: PBUSH1D,
                     property_id: np.ndarray,
                     data: np.ndarray,
                     domain_id: np.ndarray) -> None:
    #('PID', 'K', 'C', 'M', 'ALPHA', 'SA', 'EA', 'TYPEA', 'CVT', 'CVC',
    # 'EXPVT', 'EXPVC', 'IDTSU', 'IDTCU', 'IDTSUD', 'IDCSUD', 'TYPES', 'IDTS', 'IDCS',
    # 'IDTDU1', 'IDCDU1', 'TYPED', 'IDTD1', 'IDTD2', 'IDTDV1', 'IDCDV1', 'TYPEG', 'IDTG', 'IDCG',
    # 'IDTDU2', 'IDCDU2', 'IDTDV2', 'IDCDV2', 'TYPEF', 'IDTF', 'IDCF', 'UT', 'UC', 'DOMAIN_ID')
    nproperties = len(property_id)
    #prop = model.pbush1d

    #('PID', 'K', 'C', 'M', 'ALPHA', 'SA', 'EA',
    property_id = data['PID']
    k = data['K']
    c= data['C']
    mass = data['M']
    alpha = data['ALPHA']
    sa = data['SA']
    se = data['EA']

    #'TYPEA', 'CVT', 'CVC', 'EXPVT', 'EXPVC', 'IDTSU', 'IDTCU', 'IDTSUD', 'IDCSUD',
    # shock group
    shock_type = data['TYPEA']
    shock_cvt = data['CVT']
    shock_cvc = data['CVC']
    shock_exp_vt = data['EXPVT']
    shock_exp_vc = data['EXPVC']

    id_ecs = np.zeros(nproperties, dtype=property_id.dtype)
    id_tsu = data['IDTSU']
    id_tcu = data['IDTCU']
    id_tsud = data['IDTSUD']
    id_csud = data['IDCSUD']
    shock_table = [shock_cvt, shock_cvc, shock_exp_vt, shock_exp_vc,
                   id_tsu, id_tcu,
                   id_ecs, id_tsud, id_csud,
    ]

    # spring group
    #'TYPES', 'IDTS', 'IDCS', 'IDTDU1', 'IDCDU1',
    spring_type = data['TYPES']
    id_ts = data['IDTS']
    id_cs = data['IDCS']
    id_tdu1 = data['IDTDU1']
    id_cdu1 = data['IDCDU1']

    #SPRING TYPE IDT IDC IDTDU IDCDU
    spring_table = np.stack([id_ts, id_cs, id_tdu1, id_cdu1], axis=1)

    # damper group
    #'TYPED', 'IDTD1', 'IDTD2', 'IDTDV1', 'IDCDV1',
    damper_type = data['TYPED']
    id_td1 = data['IDTD1']
    id_td2 = data['IDTD2']

    # gener
    #'TYPEG', 'IDTG', 'IDCG', 'IDTDU2', 'IDCDU2', 'IDTDV2', 'IDCDV2',
    type_gener = data['TYPEG']
    id_tg = data['IDTG']
    id_cg = data['IDCG']
    id_tdu2 = data['IDTDU2']
    id_cdu2 = data['IDCDU2']

    id_tdv1 = data['IDTDV1']
    id_cdv1 = data['IDCDV1']
    id_tdv2 = data['IDTDV2']
    id_cdv2 = data['IDCDV2']

    #'TYPEF', 'IDTF', 'IDCF', 'UT', 'UC', 'DOMAIN_ID')
    type_f = data['TYPEF']
    id_tdtf = data['IDTF']
    id_tdcf = data['IDCF']
    id_uf = data['UT']
    id_uc = data['UC']

    #spring_idt = integer(card, istart + 2, 'springIDT')
    #spring_idc = integer_or_blank(card, istart + 3, 'springIDC', default=spring_idt)
    #spring_idtdu = integer_or_blank(card, istart + 4, 'springIDTDU', default=0)
    #spring_idcdu = integer_or_blank(card, istart + 5, 'springIDCDU', default=spring_idtdu)


    #DAMPER TYPE IDT IDC IDTDV IDCDV
    damper_table = np.stack([id_td1, id_td1, id_tdv1, id_tdv2, id_cdv1, id_cdv2], axis=1)
    #damper_idt = integer(card, istart + 2, 'damperIDT')
    #if damper_type == 'TABLE':
        #damper_idc = integer_or_blank(card, istart + 3, 'damperIDC', default=damper_idt)
        #damper_idtdv = integer_or_blank(card, istart + 4, 'damperIDTDV', default=0)
        #damper_idcdv = integer_or_blank(card, istart + 5, 'damperIDCDV', default=damper_idtdv)

    #shock_type = string_or_blank(card, istart + 1, 'shockType')
    #shock_cvt = double(card, istart + 2, 'shockCVT')
    #shock_cvc = double_or_blank(card, istart + 3, 'shockCVC', default=shock_cvt)
    #shock_exp_vt = double_or_blank(card, istart + 4, 'shockExpVT', default=1.0)
    #shock_exp_vc = double_or_blank(card, istart + 5,
                                   #'shockExpVC', default=shock_exp_vt)
    #shock_idets = 0
    #shock_idecs = 0
    #shock_idetsd = 0
    #shock_idecsd = 0
    #shock_idts = integer(card, istart + 6, 'shockIDTS')

    ncards = len(property_id)
    spring_equation = np.zeros((ncards, 4), dtype='int32')
    damper_equation = np.zeros((ncards, 4), dtype='int32')
    shock_equation = np.zeros((ncards, 4), dtype='int32')

    #GENER TYPE IDT IDC IDTDU IDCDU IDTDV IDCDV
    gener_equation = [id_tg, id_cg, id_tdu2, id_cdu2, id_tdv2, id_cdv2]
    prop._save(property_id, k, c, sa, se, mass,
               spring_type, spring_table, spring_equation,
               damper_type, damper_table, damper_equation,
               shock_type, shock_table, shock_equation,
               gener_equation)

    prop.domain_id = domain_id
    #return prop

def _load_h5_pbarl(prop: PBARL, group: Group):
    identity = group['IDENTITY'].read()
    #dtype([('PID', '<i8'), ('MID', '<i8'), ('GROUP', 'S8'), ('TYPE', 'S8'),
           #('INFO_POS', '<i8'), ('INFO_LEN', '<i8'), ('DOMAIN_ID', '<i8')])

    info = group['INFO'].read()
    #dtype([('VALUE', '<f8')])

    #prop = model.pbarl
    property_id = identity['PID']
    nproperties = len(property_id)

    bar_types = cast_encoding_strip(identity['TYPE'], 'latin1')
    bar_type_to_ndim = prop.valid_types
    nvalues = identity['INFO_LEN']
    is_nsm = [bar_type_to_ndim[bar_type] == n - 1
              for bar_type, n in zip(bar_types, nvalues)]
    assert np.all(is_nsm), is_nsm

    i0s = identity['INFO_POS']
    values = info['VALUE']
    insm = i0s + nvalues
    nsm = values[insm - 1]
    dims = np.hstack([values[i0:i0+n-1] for i0, n in zip(i0s, nvalues)])

    prop.property_id = property_id
    prop.material_id = identity['MID']
    prop.group = cast_encoding_strip(identity['GROUP'], 'latin1')
    prop.Type = bar_types
    prop.ndim = nvalues - 1
    prop.domain_id = identity['DOMAIN_ID']
    prop.dims = dims
    prop.nsm = nsm
    prop.n = nproperties
    prop.write()
    x = 1

def _load_h5_pbeaml(model: BDF, group: Group):
    #('PID', 'MID', 'GROUP', 'TYPE', 'SECTION_POS', 'SECTION_LEN', 'DOMAIN_ID')
    identity = group['IDENTITY'].read()
    property_id = identity['PID']
    material_id = identity['MID']
    Type = identity['TYPE']
    group_ = cast_encoding_strip(identity['GROUP'], 'latin1')

    section_pos = identity['SECTION_POS']
    section_len = identity['SECTION_LEN']
    nsection_ = len(section_len)

    domain_id = identity['DOMAIN_ID']
    #nproperties = len(property_id)

    dims_ = group['DIMS'].read()
    dims = dims_['DIM']

    #('SO', 'RDIST', 'DIMS_POS', 'DIMS_LEN', 'NSM')
    section = group['SECTION'].read()
    # so = section['SO']
    xxb = section['RDIST']
    dims_pos = section['DIMS_POS']
    dims_len = section['DIMS_LEN']
    ndim_ = len(dims_len)
    nsm = section['NSM']

    prop = model.pbeaml

    ## TODO: better cast int to string
    int_to_str_map_sout = {
        0: 'NO',
        1: 'YES',
    }
    so = np.array([int_to_str_map_sout[souti] for souti in section['SO']], dtype='|U8')

    model.log.warning(f'skipping PBEAML')

    #idim = dims_pos - 1
    idim = np.hstack([
        dims_pos.reshape(ndim_, 1),
        (dims_pos + dims_len).reshape(ndim_, 1)
    ])
    ndim = dims_len

    istation = np.hstack([
        section_pos.reshape(nsection_, 1),
        (section_pos + section_len).reshape(nsection_, 1)
    ])
    nstation = section_len
    return
    prop._save(property_id, material_id, idim, ndim, istation, nstation, Type, group_,
               xxb, dims, so, nsm)
    prop.write()
    x = 1

def _load_h5_pcomp(model: BDF, group: Group):
    identity = group['IDENTITY'].read()
    #dtype([('PID', '<i8'), ('NPLIES', '<i8'), ('Z0', '<f8'), ('NSM', '<f8'),
    #       ('SB', '<f8'), ('FT', '<i8'), ('TREF', '<f8'), ('GE', '<f8'),
    #       ('PLY_POS', '<i8'), ('PLY_LEN', '<i8'), ('DOMAIN_ID', '<i8')])

    ply = group['PLY'].read()
    #dtype([('MID', '<i8'), ('T', '<f8'), ('THETA', '<f8'), ('SOUT', '<i8')])

    prop = model.pcomp
    property_id = identity['PID']
    domain_id = identity['DOMAIN_ID']
    nproperties = len(property_id)

    property_id = property_id
    nlayer = identity['NPLIES']
    z0 = identity['Z0']
    nsm = identity['NSM']
    shear_bonding = identity['SB']

    ## TODO: better cast int to string
    int_to_str_map_ft = {
        0: '',
    }
    failure_theory = np.array([int_to_str_map_ft[souti] for souti in identity['FT']], dtype='|U8')

    tref = identity['TREF']
    ge = identity['GE']

    material_id = ply['MID']
    thickness = ply['T']
    theta = ply['THETA']

    ## TODO: better cast int to string
    int_to_str_map_sout = {
        0: 'NO',
        1: 'YES',
    }
    sout = np.array([int_to_str_map_sout[souti] for souti in ply['SOUT']], dtype='|U8')
    prop.sout = sout

    # identity
    #dtype=[('PID', '<i8'), ('NPLIES', '<i8'), ('Z0', '<f8'), ('NSM', '<f8'),
           #('SB', '<f8'), ('FT', '<i8'), ('TREF', '<f8'), ('GE', '<f8'),
           #('PLY_POS', '<i8'), ('PLY_LEN', '<i8'), ('DOMAIN_ID', '<i8')])

    nprop = len(property_id)
    lam = np.full(nprop, '', dtype='|U6')
    prop.n = nproperties
    prop._save(
        property_id, nlayer, material_id, thickness,
        sout, theta, z0, nsm, shear_bonding, failure_theory,
        tref, ge, lam)

    prop.domain_id = domain_id
    prop.write()
    x = 1

