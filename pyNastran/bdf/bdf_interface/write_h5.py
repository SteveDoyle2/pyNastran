from __future__ import annotations
from collections import defaultdict
from typing import Any, TYPE_CHECKING
from tables import File, Int64Col, Float64Col, StringCol
import numpy as np
from pyNastran.utils.numpy_utils import integer_types, float_types
if TYPE_CHECKING:
    from pyNastran.bdf.bdf import BDF


def write_h5_geometry(model: BDF, h5file, nastran_group, input_group=None):
    #node_group = nastran_group.create_group('NODE')
    #nastran_group = h5file.create_group('/', 'NASTRAN')
    if input_group is None:
        input_group = h5file.create_group(nastran_group, 'INPUT')

    write_h5_parameters(model.params, h5file, input_group)
    write_h5_nodes(model, h5file, input_group)
    write_h5_elements(model, h5file, input_group)
    write_h5_properties(model, h5file, input_group)
    write_h5_materials(model, h5file, input_group)

def write_h5_parameters(params: dict[str, Any],
                        h5file: File, input_group):
    if len(params) == 0:
        return
    parameter_group = h5file.create_group(input_group, 'PARAMETER')
    pvt_group = h5file.create_group(parameter_group, 'PVT')

    ints = []
    floats = []
    strings = []
    domain_id = 1
    for key, param in params.items():
        if len(param.values) == 1:
            value = param.values[0]
            if isinstance(value, str):
                strings.append([key, value, domain_id])
            elif isinstance(value, integer_types):
                ints.append([key, value, domain_id])
            elif isinstance(value, float_types):
                floats.append([key, value, domain_id])
            else:  # pragma: no cover
                raise NotImplementedError(param)
        else:  # pragma: no cover
            raise NotImplementedError(param)

    data_types = [
        ('CHAR', strings, StringCol(8, pos=1)),
        ('INT', ints, Int64Col(pos=1)),
        ('DOUBLE', floats, Float64Col(pos=1)), # TODO: verify name
    ]
    for name, data, h5_type in data_types:
        if len(data) == 0:
            continue
        nparam = len(data)
        h5_dict = {
            'NAME': StringCol(8, pos=0),
            'VALUE': h5_type,
            'DOMAIN_ID': Int64Col(pos=2),
        }
        table = h5file.create_table(pvt_group, name, h5_dict)
        arr = np.empty(nparam, dtype=table.dtype)
        arr["NAME"] = [item[0] for item in data]
        arr["VALUE"] = [item[1] for item in data]
        arr["DOMAIN_ID"] = np.ones(nparam, dtype='int64')
        table.append(arr)

def write_h5_nodes(model: BDF, h5file: File, input_group):
    node_group = h5file.create_group(input_group, 'NODE')
    write_grid(model, h5file, node_group)

    if len(model.coords) > 1:
        coord_group = h5file.create_group(input_group, 'COORDINATE_SYSTEM')
        write_coord(model, h5file, coord_group)


def write_h5_elements(model: BDF, h5file: File, input_group):
    log = model.log
    elements_by_type = defaultdict(list)
    #masses_by_type = defaultdict(list)
    for eid, elem in sorted(model.elements.items()):
        elements_by_type[elem.type].append(eid)
    for eid, elem in sorted(model.rigid_elements.items()):
        elements_by_type[elem.type].append(eid)
    for eid, elem in sorted(model.masses.items()):
        elements_by_type[elem.type].append(eid)
    if len(elements_by_type) == 0:
        return

    element_group = h5file.create_group(input_group, 'ELEMENT')
    cls_mapper = {
        'CROD': write_eid_pid_nodes,
        'CTUBE': write_eid_pid_nodes,
        'CSHEAR': write_eid_pid_nodes,
        'CONROD': write_conrod,
        'CBAR': write_cbar,
        'CTRIA3': write_ctria3_ctria6,
        'CTRIA6': write_ctria3_ctria6,
        'CQUAD4': write_cquad4_cquad8,
        'CQUAD8': write_cquad4_cquad8,
        'CTETRA': write_csolid,
        'CPYRAM': write_csolid,
        'CPENTA': write_csolid,
        'CHEXA': write_csolid,
        'CONM2': write_conm2,
    }
    for etype, eids in elements_by_type.items():
        if etype in cls_mapper:
            write_h5_func = cls_mapper[etype]
            write_h5_func(model, eids, h5file, element_group)
        else:
            log.warning(f'skipping h5 element {etype!r}')

def write_conrod(model: BDF, eids: list[int],
                 h5file: File, element_group):
    neid = len(eids)
    h5_dict = {
        'EID': Int64Col(pos=0),  # Signed 64-bit integer
        'PID': Int64Col(pos=1),
        'N': Int64Col(pos=2),
        'A': Float64Col(pos=3),
        'J': Float64Col(pos=4),
        'NSM': Float64Col(pos=5),
        'DOMAIN_ID': Int64Col(pos=5),
    }
    table = h5file.create_table(element_group, 'CONROD', h5_dict)
    arr = np.empty(neid, dtype=table.dtype)
    pids = np.array(neid, dtype='int32')
    nodes = np.array((neid, neid), dtype='int32')
    mids = np.array(neid, dtype='int32')
    a = np.full(neid, -1, dtype='float64')
    j = np.full(neid, np.nan, dtype='float64')
    nsm = np.full(neid, np.nan, dtype='float64')
    for i, eid in enumerate(eids):
        elem = model.elements[eid]
        pids[i] = elem.pid
        nodes[i, :] = elem.nodes
        mids[i] = elem.mid
        a[i] = elem.a
        j[i] = elem.j
        nsm[i] = elem.nsm
    arr["EID"] = eids
    arr["MID"] = mids
    arr["N"] = nodes
    arr["A"] = a
    arr["J"] = j
    arr["N"] = nsm
    arr["DOMAIN_ID"] = np.ones(neid, dtype='int64')
    table.append(arr)

def write_cbar(model: BDF, eids: list[int],
                h5file: File, element_group):
    neid = len(eids)
    if neid == 0:
        return
    h5_dict = {
        'EID': Int64Col(pos=0),  # Signed 64-bit integer
        'PID': Int64Col(pos=1),
        'GA': Int64Col(pos=2),
        'GB': Int64Col(pos=3),
        'FLAG': Int64Col(pos=4),
        'X1': Float64Col(pos=5),
        'X2': Float64Col(pos=6),
        'X3': Float64Col(pos=7),
        # dtype = [('EID', '<i8'), ('PID', '<i8'), ('GA', '<i8'), ('GB', '<i8'),
        #          ('FLAG', '<i8'), ('X1', '<f8'), ('X2', '<f8'), ('X3', '<f8'),
        #          ('GO', '<i8'), ('PA', '<i8'), ('PB', '<i8'),
        #          ('W1A', '<f8'), ('W2A', '<f8'), ('W3A', '<f8'),
        #          ('W1B', '<f8'), ('W2B', '<f8'), ('W3B', '<f8'),
        #          ('DOMAIN_ID', '<i8')])
        'GO': Int64Col(pos=8),
        'PA': Int64Col(pos=9),
        'PB': Int64Col(pos=10),
        'W1A': Int64Col(pos=11),
        'W2A': Float64Col(pos=12),
        'W3A': Float64Col(pos=13),
        'W1B': Int64Col(pos=14),
        'W2B': Float64Col(pos=15),
        'W3B': Float64Col(pos=16),
        'DOMAIN_ID': Int64Col(pos=17),
    }
    table = h5file.create_table(element_group, 'CBAR', h5_dict)
    arr = np.empty(neid, dtype=table.dtype)
    for i, eid in enumerate(eids):
        elem = model.elements[eid]
        nodes = elem.nodes
        arr["EID"][i] = eid
        arr["PID"][i] = elem.pid
        arr["GA"][i] = nodes[0]
        arr["GB"][i] = nodes[1]

        arr["FLAG"][i] = 1  # GGG
        if elem.g0 is None:
            arr["X1"][i] = elem.x[0]
            arr["X2"][i] = elem.x[1]
            arr["X3"][i] = elem.x[2]
            arr["GO"][i] = 0
        else:
            arr["X1"][i] = np.nan
            arr["X2"][i] = np.nan
            arr["X3"][i] = np.nan
            arr["GO"][i] = elem.g0

        arr["PA"][i] = elem.pa
        arr["PB"][i] = elem.pb

        arr["W1A"][i] = elem.wa[0]
        arr["W2A"][i] = elem.wa[1]
        arr["W3A"][i] = elem.wa[2]

        arr["W1B"][i] = elem.wb[0]
        arr["W2B"][i] = elem.wb[1]
        arr["W3B"][i] = elem.wb[2]

    arr["DOMAIN_ID"] = np.ones(neid, dtype='int64')
    table.append(arr)

def write_conm2(model: BDF, eids: list[int],
                h5file: File, element_group):
    neid = len(eids)
    if neid == 0:
        return
    h5_dict = {
        # dtype([('EID', '<i8'), ('G', '<i8'), ('CID', '<i8'), ('M', '<f8'),
        #        ('X1', '<f8'), ('X2', '<f8'), ('X3', '<f8'), ('I1', '<f8'), ('I2', '<f8', (2,)),
        #        ('I3', '<f8', (3,)), ('DOMAIN_ID', '<i8')])
        'EID': Int64Col(pos=0),  # Signed 64-bit integer
        'G': Int64Col(pos=1),
        'CID': Int64Col(pos=2),
        'M': Float64Col(pos=3),
        'X1': Float64Col(pos=4),
        'X2': Float64Col(pos=5),
        'X3': Float64Col(pos=6),
        # I11, [I21, I22], [I31, I32, I33]
        'I1': Float64Col(pos=7),
        'I2': Float64Col(shape=(2,), pos=8),
        'I3': Float64Col(shape=(3,), pos=9),
        'DOMAIN_ID': Int64Col(pos=10),
    }
    table = h5file.create_table(element_group, 'CONM2', h5_dict)
    arr = np.empty(neid, dtype=table.dtype)
    for i, eid in enumerate(eids):
        elem = model.masses[eid]
        arr["EID"][i] = eid
        arr["G"][i] = elem.nid
        arr["CID"][i] = elem.cid
        arr["X1"][i] = elem.X[0]
        arr["X2"][i] = elem.X[1]
        arr["X3"][i] = elem.X[2]

        arr["I1"][i] = elem.I[0]    # 1
        arr["I2"][i] = elem.I[1:3]  # 2
        arr["I3"][i] = elem.I[3:]   # 3

    arr["DOMAIN_ID"] = np.ones(neid, dtype='int64')
    table.append(arr)

def write_ctria3_ctria6(model: BDF, eids: list[int],
                        h5file: File, element_group):
    if len(eids) == 0:
        return
    neid = len(eids)
    eid0 = eids[0]
    elem0 = model.elements[eid0]
    # etype = elem0.type
    nnode_map = {
        'CTRIA3': 3,
        'CTRIA6': 6,
    }
    card_type = elem0.type
    nnode = nnode_map[card_type]
    # dtype([('EID', '<i8'), ('PID', '<i8'), ('G', '<i8', (4,)), ('THETA', '<f8'),
    #        ('ZOFFS', '<f8'), ('TFLAG', '<i8'), ('T', '<f8', (3,)),
    #        ('MCID', '<i8'), ('DOMAIN_ID', '<i8')])
    h5_dict = {
        'EID': Int64Col(pos=0),  # Signed 64-bit integer
        'PID': Int64Col(pos=1),
        'G': Int64Col(shape=(nnode,), pos=2),
        'THETA': Float64Col(pos=3),
        'ZOFFS': Float64Col(pos=4),
        'TFLAG': Int64Col(pos=5),
        'T': Float64Col(shape=(3,), pos=6),
        'MCID': Int64Col(pos=7),
        'DOMAIN_ID': Int64Col(pos=8),
    }
    table = h5file.create_table(element_group, card_type, h5_dict)
    arr = np.empty(neid, dtype=table.dtype)
    pids = np.zeros(neid, dtype='int32')
    nodes = np.zeros((neid, nnode), dtype='int32')
    tflag = np.full(neid, -1, dtype='int32')
    t = np.zeros((neid, nnode), dtype='float64')
    zoffset = np.full(neid, -1, dtype='float64')
    theta = np.full(neid, -1, dtype='float64')
    mcid = np.full(neid, np.nan, dtype='float64')
    is_theta = np.zeros(neid, dtype='bool')
    for i, eid in enumerate(eids):
        elem = model.elements[eid]
        pids[i] = elem.pid
        nodes[i, :] = elem.nodes
        zoffset[i] = elem.zoffset
        theta_mcid = elem.theta_mcid
        tflag[i] = elem.tflag
        ti = [elem.T1, elem.T2, elem.T3]
        t[i, :] = ti
        if isinstance(theta_mcid, integer_types):
            #is_theta[i] = False
            mcid[i] = theta_mcid
        else:
            is_theta[i] = True
            theta[i] = theta_mcid
    arr["EID"] = eids
    arr["PID"] = pids
    arr["G"] = nodes
    arr["THETA"] = theta
    arr["MCID"] = mcid
    arr['ZOFFS'] = zoffset
    arr['TFLAG'] = tflag
    arr['T'] = t
    # arr["IS_THETA"] = is_theta
    arr["DOMAIN_ID"] = np.ones(neid, dtype='int64')
    table.append(arr)

def write_cquad4_cquad8(model: BDF, eids: list[int],
                        h5file: File, element_group):
    neid = len(eids)
    eid0 = eids[0]
    elem0 = model.elements[eid0]
    etype = elem0.type
    nnode_map = {
        'CQUAD4': 4,
        'CQUAD8': 8,
    }
    card_type = elem0.type
    nnode = nnode_map[card_type]
    # dtype([('EID', '<i8'), ('PID', '<i8'), ('G', '<i8', (4,)), ('THETA', '<f8'),
    #        ('ZOFFS', '<f8'), ('TFLAG', '<i8'), ('T', '<f8', (4,)),
    #        ('MCID', '<i8'), ('DOMAIN_ID', '<i8')])

    h5_dict = {
        'EID': Int64Col(pos=0),  # Signed 64-bit integer
        'PID': Int64Col(pos=1),
        'G': Int64Col(shape=(nnode,), pos=2),
        'THETA': Float64Col(pos=3),
        'ZOFFS': Float64Col(pos=4),
        'TFLAG': Int64Col(pos=5),
        'T': Float64Col(shape=(4,), pos=6),
        'MCID': Int64Col(pos=7),
        'DOMAIN_ID': Int64Col(pos=8),
    }
    table = h5file.create_table(element_group, card_type, h5_dict)
    arr = np.empty(neid, dtype=table.dtype)
    pids = np.zeros(neid, dtype='int32')
    nodes = np.zeros((neid, nnode), dtype='int32')
    tflag = np.full(neid, -1, dtype='int32')
    t = np.zeros((neid, nnode), dtype='float64')
    zoffset = np.full(neid, -1, dtype='float64')
    theta = np.full(neid, -1, dtype='float64')
    mcid = np.full(neid, np.nan, dtype='float64')
    is_theta = np.zeros(neid, dtype='bool')
    for i, eid in enumerate(eids):
        elem = model.elements[eid]
        pids[i] = elem.pid
        nodes[i, :] = elem.nodes
        theta_mcid = elem.theta_mcid
        tflag[i] = elem.tflag
        zoffset[i] = elem.zoffset
        ti = [elem.T1, elem.T2, elem.T3, elem.T4]
        t[i, :] = ti
        if isinstance(theta_mcid, integer_types):
            #is_theta[i] = False
            mcid[i] = theta_mcid
        else:
            is_theta[i] = True
            theta[i] = theta_mcid
    arr["EID"] = eids
    arr["PID"] = pids
    arr["G"] = nodes
    arr["THETA"] = theta
    arr["MCID"] = mcid
    arr['ZOFFS'] = zoffset
    arr['TFLAG'] = tflag
    arr['T'] = t
    # arr["IS_THETA"] = is_theta
    arr["DOMAIN_ID"] = np.ones(neid, dtype='int64')
    table.append(arr)

def write_csolid(model: BDF, eids: list[int],
                 h5file: File, element_group):
    # {'CROD', 'PROD', 'CSHEAR'}:
    neid = len(eids)
    eid0 = eids[0]
    elem0 = model.elements[eid0]
    etype = elem0.type
    nnode_map = {
        'CTETRA': 10,
        'CPYRAM': 13,
        'CPENTA': 15,
        'CHEXA': 20,
    }
    card_type = elem0.type
    nnode = nnode_map[card_type]
    h5_dict = {
        'EID': Int64Col(pos=0),  # Signed 64-bit integer
        'PID': Int64Col(pos=1),
        'N': Int64Col(pos=nnode),
        'DOMAIN_ID': Int64Col(pos=4),
    }
    table = h5file.create_table(element_group, card_type, h5_dict)
    arr = np.empty(neid, dtype=table.dtype)
    pids = np.array(neid, dtype='int32')
    nodes = np.array((neid, nnode), dtype='int32')
    for i, eid in enumerate(eids):
        elem = model.elements[eid]
        pids[i] = elem.pid
        nodesi = elem.nodes
        nnodei = len(nodesi)
        nodes[i, :nnodei] = nodesi
    arr["EID"] = eids
    arr["PID"] = pids
    arr["N"] = nodes
    arr["DOMAIN_ID"] = np.ones(neid, dtype='int64')
    table.append(arr)

def write_eid_pid_nodes(model: BDF, eids: list[int],
                        h5file: File, element_group):
    # {'CROD', 'PROD', 'CSHEAR'}:
    neid = len(eids)
    eid0 = eids[0]
    elem0 = model.elements[eid0]
    nnode = len(elem0.nodes)
    h5_dict = {
        'EID': Int64Col(pos=0),  # Signed 64-bit integer
        'PID': Int64Col(pos=1),
        'N': Int64Col(shape=(nnode,), pos=2),
        'DOMAIN_ID': Int64Col(pos=3),
    }
    card_type = elem0.type
    table = h5file.create_table(element_group, card_type, h5_dict)
    arr = np.empty(neid, dtype=table.dtype)
    pids = []
    nodes = []
    for eid in eids:
        elem = model.elements[eid]
        pids.append(elem.pid)
        nodes.append(elem.nodes)
    arr["EID"] = eids
    arr["PID"] = pids
    arr["N"] = nodes
    arr["DOMAIN_ID"] = np.ones(neid, dtype='int64')
    table.append(arr)

def write_h5_properties(model: BDF, h5file: File, input_group):
    log = model.log
    properties_by_type = defaultdict(list)
    for pid, prop in model.properties.items():
        properties_by_type[prop.type].append(pid)
    for pid, prop in model.properties_mass.items():
        properties_by_type[prop.type].append(pid)

    if len(properties_by_type) == 0:
        return
    property_group = h5file.create_group(input_group, 'PROPERTY')
    cls_mapper = {
        'PROD': write_prod,
        # 'PBARL': write_pbarl,
        'PSHELL': write_pshell,
        'PCOMP': write_pcomp,
    }
    for prop_type, pids in properties_by_type.items():
        if prop_type in cls_mapper:
            write_h5_func = cls_mapper[prop_type]
            write_h5_func(model, pids, h5file, property_group)
        else:
            log.warning(f'skipping h5 property {prop_type!r}')


def write_h5_materials(model: BDF, h5file: File, input_group):
    log = model.log
    materials_by_type = defaultdict(list)
    for mid, mat in model.materials.items():
        materials_by_type[mat.type].append(mid)
    if len(materials_by_type) == 0:
        return

    material_group = h5file.create_group(input_group, 'MATERIAL')
    cls_mapper = {
        'MAT1': write_mat1,
        'MAT8': write_mat8,
    }
    for mat_type, mids in materials_by_type.items():
        if mat_type in cls_mapper:
            write_h5_func = cls_mapper[mat_type]
            write_h5_func(model, mids, h5file, material_group)
        else:
            log.warning(f'skipping h5 material {mat_type!r}')


def write_grid(model: BDF, h5file: File, node_group):
    h5_dict = {
        'ID': Int64Col(pos=0),  # Signed 64-bit integer
        'CP': Int64Col(pos=1),  # Signed 64-bit integer
        'X': Float64Col(shape=(3,), pos=2),  # double (double-precision)
        'CD': Int64Col(pos=3),  # Signed 64-bit integer
        'PS': Int64Col(pos=4),  # Signed 64-bit integer
        'SEID': Int64Col(pos=5),  # Signed 64-bit integer
        'DOMAIN_ID': Int64Col(pos=6),
    }
    node_id = []
    cp = []
    xyz = []
    cd = []
    ps = []
    seid = []
    for nid, node in sorted(model.nodes.items()):
        psi = 0 if node.ps == '' else int(node.ps)
        node_id.append(nid)
        cp.append(node.cp)
        xyz.append(node.xyz)
        cd.append(node.cd)
        seid.append(node.seid)
        ps.append(psi)

    table = h5file.create_table(node_group, 'GRID', h5_dict)
    ngrid = len(seid)
    arr = np.empty(ngrid, dtype=table.dtype)
    arr["ID"] = node_id
    arr["CP"] = cp
    arr["X"] = xyz
    arr["CD"] = cd
    arr["PS"] = ps
    arr["SEID"] = seid
    arr["DOMAIN_ID"] = np.ones(ngrid, dtype='int32')
    table.append(arr)


def write_coord(model: BDF, h5file: File, coord_group):
    ncoord = len(model.coords) - 1
    if ncoord == 0:
        return

    transform_group = h5file.create_group(coord_group, 'TRANSFORMATION')
    xform_data_dict = {'DATA': Float64Col()}
    h5_dict = {
        'CID': Int64Col(pos=0),  # Signed 64-bit integer
        'RID': Int64Col(pos=1),  # Signed 64-bit integer
        'A1': Float64Col(pos=2),  # double (double-precision)
        'A2': Float64Col(pos=3),
        'A3': Float64Col(pos=4),
        'B1': Float64Col(pos=5),
        'B2': Float64Col(pos=6),
        'B3': Float64Col(pos=7),
        'C1': Float64Col(pos=8),
        'C2': Float64Col(pos=9),
        'C3': Float64Col(pos=10),
        'DOMAIN_ID': Int64Col(pos=6),
    }
    i = -1
    rdata_list = []
    coord_id = np.zeros(ncoord, dtype='int32')
    coord_type_int = np.zeros(ncoord, dtype='int32')

    cids_by_type = defaultdict(list)
    icoords_by_type = defaultdict(list)
    assert 0 in model.coords, model.coords

    map_coord_type = {'R': 1, 'C': 2, 'S': 3}
    for cid, coord in sorted(model.coords.items()):
        if cid == 0:
            continue
        class_name = coord.type
        icoordi = int(class_name[-2])  # 2 in CORD2R
        coord_type_stri = class_name[-1]  # CORD2R
        coord_type_inti = map_coord_type[coord_type_stri]
        #class_name = f'CORD{icoord}{coord_type}'
        #coords_by_type[class_name].append(i)
        #xform = self.T[i, :, :].ravel()
        e1 = coord.e1
        xform = coord.beta().ravel()
        coord_id[i] = cid
        coord_type_int[i] = coord_type_inti
        rdata_list.append(e1)
        rdata_list.append(xform)
        cids_by_type[class_name].append(cid)
        icoords_by_type[class_name].append(i)
        i += 1
    rdata = np.hstack(rdata_list, dtype='float64')
    assert len(rdata) > 0
    table_rdata = h5file.create_table(transform_group, 'RDATA', xform_data_dict)
    arr_rdata = np.empty(len(rdata), dtype=table_rdata.dtype)
    arr_rdata['DATA'] = rdata
    table_rdata.append(arr_rdata)

    xform_identity_dict = {
        'CID': Int64Col(pos=0),  # Signed 64-bit integer
        'TYPE': Int64Col(pos=1),
        'INDEX': Int64Col(pos=2),
        'RINDEX': Int64Col(pos=3),
        'DOMAIN_ID': Int64Col(pos=4),
    }
    table_identity = h5file.create_table(transform_group, 'INDENTITY', xform_identity_dict)

    coord_id = np.array(coord_id, dtype='int32')
    arr_identity = np.empty(ncoord, dtype=table_identity.dtype)
    arr_identity['CID'] = coord_id
    arr_identity['TYPE'] = coord_type_int
    arr_identity['INDEX'] = coord_id * 0
    arr_identity['RINDEX'] = np.arange(ncoord) * 12 + 1
    arr_identity['DOMAIN_ID'] = np.ones(ncoord)
    table_identity.append(arr_identity)

    for coord_type, icoords in icoords_by_type.items():
        cids = cids_by_type[coord_type]
        table = h5file.create_table(coord_group, coord_type, h5_dict)
        ncoord = len(cids)
        arr = np.empty(ncoord, dtype=table.dtype)
        for icoord, cid in zip(icoords, cids):
            coord = model.coords[cid]
            arr["CID"][icoord] = cid  # self.coord_id[icoords]
            arr["RID"][icoord] = coord.rid  # self.ref_coord_id[icoords]
            arr["A1"][icoord] = coord.e1[0]  # self.e1[icoords, 0]
            arr["A2"][icoord] = coord.e1[1]  # self.e1[icoords, 1]
            arr["A3"][icoord] = coord.e1[2]  # self.e1[icoords, 2]
            arr["B1"][icoord] = coord.e2[0]  # self.e2[icoords, 0]
            arr["B2"][icoord] = coord.e2[1]  # self.e2[icoords, 1]
            arr["B3"][icoord] = coord.e2[2]  # self.e2[icoords, 2]
            arr["C1"][icoord] = coord.e3[0]  # self.e3[icoords, 0]
            arr["C2"][icoord] = coord.e3[1]  # self.e3[icoords, 1]
            arr["C3"][icoord] = coord.e3[2]  # self.e3[icoords, 2]
        arr["DOMAIN_ID"] = np.ones(ncoord)
        table.append(arr)


def write_prod(model: BDF, pids: list[int],
               h5file: File, property_group):
    nprop = len(pids)
    h5_dict = {
        'PID': Int64Col(pos=0),  # Signed 64-bit integer
        'MID': Int64Col(pos=1),
        'A': Float64Col(pos=2),
        'J': Float64Col(pos=3),
        'C': Float64Col(pos=4),
        'NSM': Float64Col(pos=5),
        'DOMAIN_ID': Int64Col(pos=6),
    }

    table = h5file.create_table(property_group, 'PROD', h5_dict)
    arr = np.empty(nprop, dtype=table.dtype)
    arr["PID"] = pids
    for i, pid in enumerate(pids):
        prop = model.properties[pid]
        arr["MID"][i] = prop.mid
        arr["A"][i] = prop.A
        arr["J"][i] = prop.J
        arr["C"][i] = prop.c
        arr["NSM"][i] = prop.nsm
    arr["DOMAIN_ID"] = np.ones(nprop, dtype='int64')
    table.append(arr)

def write_pbarl(model: BDF, pids: list[int],
                h5file: File, property_group):
    nprop = len(pids)
    h5_dict = {
        'PID': Int64Col(pos=0),  # Signed 64-bit integer
        'MID': Int64Col(pos=1),
        'A': Float64Col(pos=2),
        'J': Float64Col(pos=3),
        'NSM': Float64Col(pos=4),
        'DOMAIN_ID': Int64Col(pos=5),
    }

    table = h5file.create_table(property_group, 'PBARL', h5_dict)
    arr = np.empty(nprop, dtype=table.dtype)
    arr["PID"] = pids
    for i, pid in enumerate(pids):
        prop = model.properties[pid]
        arr["MID"][i] = prop.mid
        arr["A"][i] = prop.A
        # arr["J"][i] = prop.J
        arr["NSM"][i] = prop.nsm
    arr["DOMAIN_ID"] = np.ones(nprop, dtype='int64')
    table.append(arr)

def write_pshell(model: BDF, pids: list[int],
                 h5file: File, property_group):
    nprop = len(pids)
    h5_dict = {
        'PID': Int64Col(pos=0),
        'MID1': Int64Col(pos=1),
        'T': Float64Col(pos=2),
        'MID2': Int64Col(pos=3),
        'BK': Float64Col(pos=4),
        'MID3': Int64Col(pos=5),
        'TS': Float64Col(pos=6),
        'NSM': Float64Col(pos=7),
        'Z1': Float64Col(pos=8),
        'Z2': Float64Col(pos=9),
        'MID4': Int64Col(pos=10),
        'DOMAIN_ID': Int64Col(pos=11),
    }
    table = h5file.create_table(property_group, 'PSHELL', h5_dict)
    arr = np.empty(nprop, dtype=table.dtype)
    arr["PID"] = pids
    for i, pid in enumerate(pids):
        prop = model.properties[pid]
        mids = [prop.mid1, prop.mid2, prop.mid3, prop.mid4]
        mids = [mid if mid is not None else 0
                for mid in mids]
        arr["MID1"][i] = mids[0]  # self.material_id[:, 0]
        arr["MID2"][i] = mids[1]  # self.material_id[:, 1]
        arr["MID3"][i] = mids[2]  # self.material_id[:, 2]
        arr["MID4"][i] = mids[3]  # self.material_id[:, 3]
        arr["T"][i] = prop.t # self.t
        arr["BK"][i] = prop.twelveIt3 # self.twelveIt3
        arr["TS"][i] = prop.tst # self.tst
        arr["Z1"][i] = prop.z1 # self.z[:, 0]
        arr["Z2"][i] = prop.z2 # self.z[:, 1]
        arr["NSM"][i] = prop.nsm # self.nsm
    arr["DOMAIN_ID"] = np.ones(nprop, dtype='int64')
    table.append(arr)

def write_pcomp(model: BDF, pids: list[int],
                h5file: File, property_group):
    nprop = len(pids)
    if nprop == 0:
        return
    pcomp_group = h5file.create_group(property_group, 'PCOMP')

    h5_identity_dict = {
        'PID': Int64Col(pos=0),
        'NPLIES': Int64Col(pos=1),
        'Z0': Float64Col(pos=2),
        'NSM': Float64Col(pos=3),
        'SB': Float64Col(pos=4),
        'FT': Int64Col(pos=5),
        'TREF': Float64Col(pos=6),
        'GE': Float64Col(pos=7),
        'PLY_POS': Int64Col(pos=8),
        'PLY_LEN': Int64Col(pos=9),
        'DOMAIN_ID': Int64Col(pos=10),
    }
    #nlayer = len(model.pcomp.material_id)
    identity_table = h5file.create_table(pcomp_group, 'IDENTITY', h5_identity_dict)
    arr = np.empty(nprop, dtype=identity_table.dtype)

    arr["PID"] = pids
    nplies_total = 0
    material_id = []
    thickness = []
    theta = []
    sout_int = []
    for i, pid in enumerate(pids):
        prop = model.properties[pid]
        nplies = len(prop.thicknesses)
        arr["NPLIES"][i] = nplies
        arr["Z0"][i] = prop.z0
        arr["NSM"][i] = prop.nsm
        arr["SB"][i] = prop.sb
        arr["FT"][i] = prop.ft_int
        arr["TREF"][i] = prop.tref
        arr["GE"][i] = prop.ge
        arr["PLY_POS"][i] = nplies_total
        arr["PLY_LEN"][i] = nplies_total+nplies
        material_id.extend(prop.mids)
        thickness.extend(prop.thicknesses)
        theta.extend(prop.thetas)
        sout_int.extend(prop.sout_ints)
        # material_id.extend(prop.material_ids)
        nplies_total += nplies
    arr["DOMAIN_ID"] = np.ones(nprop, dtype='int64')
    identity_table.append(arr)

    #----------------------------------
    h5_dict = {
        'MID': Int64Col(pos=0),
        'T': Float64Col(pos=1),
        'THETA': Float64Col(pos=2),
        'SOUT': Int64Col(pos=3),
    }
    ply_table = h5file.create_table(pcomp_group, 'PLY', h5_dict)
    nlayer = len(material_id)
    arr = np.empty(nlayer, dtype=ply_table.dtype)
    assert len(material_id) == len(thickness), (len(material_id), len(thickness))
    arr["MID"] = material_id
    arr["T"] = thickness
    arr["THETA"] = theta
    arr["SOUT"] = sout_int
    ply_table.append(arr)


def write_mat1(model: BDF, mids: list[int],
               h5file: File, material_group):
    nmat = len(mids)
    h5_dict = {
        # dtype([('MID', '<i8'), ('E', '<f8'), ('G', '<f8'), ('NU', '<f8'),
        # ('RHO', '<f8'), ('A', '<f8'), ('TREF', '<f8'), ('GE', '<f8'),
        # ('ST', '<f8'), ('SC', '<f8'), ('SS', '<f8'), ('MCSID', '<i8'),
        # ('DOMAIN_ID', '<i8')])
        'MID': Int64Col(pos=0),
        'E': Float64Col(pos=1),
        'G': Float64Col(pos=2),
        'NU': Float64Col(pos=3),
        'RHO': Float64Col(pos=4),
        'A': Float64Col(pos=5),
        'TREF': Float64Col(pos=6),
        'GE': Float64Col(pos=7),
        'ST': Float64Col(pos=8),
        'SC': Float64Col(pos=9),
        'SS': Float64Col(pos=10),
        'MCSID': Int64Col(pos=11),
        'DOMAIN_ID': Int64Col(pos=12),
    }
    table = h5file.create_table(material_group, 'MAT1', h5_dict)
    arr = np.empty(nmat, dtype=table.dtype)
    arr["MID"] = mids
    for i, mid in enumerate(mids):
        mat = model.materials[mid]
        assert isinstance(mat.nu, float), mat.nu
        arr["E"][i] = mat.e
        arr["G"][i] = mat.g
        arr["NU"][i] = mat.nu
        arr["RHO"][i] = mat.rho
        arr["TREF"][i] = mat.tref
        arr["A"][i] = mat.a
        arr["GE"][i] = mat.ge
        arr["MCSID"][i] = mat.mcsid
    arr["DOMAIN_ID"] = np.ones(nmat, dtype='int64')
    table.append(arr)

def write_mat8(model: BDF, mids: list[int],
               h5file: File, material_group):
    nmat = len(mids)
    h5_dict = {
        'MID': Int64Col(pos=0),
        'E1': Int64Col(pos=1),
        'E2': Int64Col(pos=2),
        'NU12': Int64Col(pos=3),
        'G12': Float64Col(pos=4),
        'G1Z': Float64Col(pos=5),
        'G2Z': Float64Col(pos=6),
        'RHO': Float64Col(pos=7),
        # dtype([('MID', '<i8'), ('E1', '<f8'), ('E2', '<f8'), ('NU12', '<f8'),
        #        ('G12', '<f8'), ('G1Z', '<f8'), ('G2Z', '<f8'), ('RHO', '<f8'),
        #        ('A1', '<f8'), ('A2', '<f8'), ('TREF', '<f8'), ('XT', '<f8'),
        #        ('XC', '<f8'), ('YT', '<f8'), ('YC', '<f8'), ('S', '<f8'),
        #        ('GE', '<f8'), ('F12', '<f8'), ('STRN', '<f8'), ('DOMAIN_ID', '<i8')])
        'A1': Float64Col(pos=8),
        'A2': Float64Col(pos=9),
        'TREF': Float64Col(pos=10),
        'XT': Float64Col(pos=11),
        'XC': Float64Col(pos=12),
        'YT': Float64Col(pos=13),
        'YC': Float64Col(pos=14),
        'S': Float64Col(pos=15),
        'GE': Float64Col(pos=16),
        'F12': Float64Col(pos=17),
        'STRN': Float64Col(pos=18),
        'DOMAIN_ID': Int64Col(pos=19),
    }
    table = h5file.create_table(material_group, 'MAT8', h5_dict)
    arr = np.empty(nmat, dtype=table.dtype)
    arr["MID"] = mids
    for i, mid in enumerate(mids):
        mat = model.materials[mid]
        arr["E1"][i] = mat.e11
        arr["E2"][i] = mat.e22
        arr["NU12"][i] = mat.nu12
        arr["G12"][i] = mat.g12
        arr["G1Z"][i] = mat.g1z
        arr["G2Z"][i] = mat.g2z
        arr["RHO"][i] = mat.rho
        arr["TREF"][i] = mat.tref
        arr["A1"][i] = mat.a1
        arr["A2"][i] = mat.a2
        arr["GE"][i] = mat.ge
        arr["XT"][i] = mat.Xt
        arr["XC"][i] = mat.Xc
        arr["YC"][i] = mat.Yt
        arr["YC"][i] = mat.Yc
        arr["S"][i] = mat.S
        arr["STRN"][i] = mat.strn
        arr["F12"][i] = mat.F12
    arr["DOMAIN_ID"] = np.ones(nmat, dtype='int64')
    table.append(arr)

