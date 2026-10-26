from __future__ import annotations
from collections import defaultdict
from typing import Any, TYPE_CHECKING
from tables import File, Int64Col, Float64Col, StringCol
import numpy as np
from pyNastran.utils.numpy_utils import integer_types, float_types
if TYPE_CHECKING:
    from pyNastran.bdf.bdf import BDF


def write_h5_geometry(model: BDF, h5file, nastran_group):
    #node_group = nastran_group.create_group('NODE')
    #nastran_group = h5file.create_group('/', 'NASTRAN')
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
        ('FLOAT', floats, Float64Col(pos=1)), # TODO: verify name
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
    write_h5_grid(model, h5file, node_group)

    if len(model.coord) > 1:
        coord_group = h5file.create_group(input_group, 'COORDINATE_SYSTEM')
        model.write_h5_coord(h5file, coord_group)


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
        'CROD': write_h5_eid_pid_nodes,
        'CTUBE': write_h5_eid_pid_nodes,
        'CSHEAR': write_h5_eid_pid_nodes,
        'CONROD': write_h5_conrod,
        'CQUAD4': write_h5_cquad4,
        'CTETRA': write_h5_solid,
        'CPYRAM': write_h5_solid,
        'CPENTA': write_h5_solid,
        'CHEXA': write_h5_solid,
    }
    for etype, eids in elements_by_type.items():
        if etype in cls_mapper:
            write_h5_func = cls_mapper[etype]
            write_h5_func(model, eids, h5file, element_group)
        else:
            log.warning(f'skipping h5 element {etype!r}')

def write_h5_conrod(model: BDF, eids: list[int],
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

def write_h5_cquad4(model: BDF, eids: list[int],
                    h5file: File, element_group):
    neid = len(eids)
    eid0 = eids[0]
    elem0 = model.elements[eid0]
    etype = elem0.type
    nnode_map = {
        'CQUAD4': 4,
    #    'CPYRAM': 13,
    #    'CPENTA': 15,
    #    'CHEXA': 20,
    }
    card_type = elem0.type
    nnode = nnode_map[card_type]
    h5_dict = {
        'EID': Int64Col(pos=0),  # Signed 64-bit integer
        'PID': Int64Col(pos=1),
        'N': Int64Col(size=(nnode,), pos=2),
        'THETA': Float64Col(pos=3),
        'MCID': Int64Col(pos=4),
        'DOMAIN_ID': Int64Col(pos=5),
    }
    table = h5file.create_table(element_group, card_type, h5_dict)
    arr = np.empty(neid, dtype=table.dtype)
    pids = np.array(neid, dtype='int32')
    nodes = np.array((neid, nnode), dtype='int32')
    theta = np.full(neid, -1, dtype='int32')
    mcid = np.full(neid, np.nan, dtype='float64')
    is_theta = np.array(neid, dtype='bool')
    for i, eid in enumerate(eids):
        elem = model.elements[eid]
        pids[i] = elem.pid
        nodes[i, :] = elem.nodes
        theta_mcid = elem.theta_mcid
        if isinstance(theta_mcid, integer_types):
            #is_theta[i] = False
            mcid[i] = theta_mcid
        else:
            is_theta[i] = True
            theta[i] = theta_mcid
    arr["EID"] = eids
    arr["PID"] = pids
    arr["N"] = nodes
    arr["THETA"] = theta
    arr["MCID"] = mcid
    arr["IS_THETA"] = is_theta
    arr["DOMAIN_ID"] = np.ones(neid, dtype='int64')
    table.append(arr)

def write_h5_solid(model: BDF, eids: list[int],
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

def write_h5_eid_pid_nodes(model: BDF, eids: list[int],
                           h5file: File, element_group):
    # {'CROD', 'PROD', 'CSHEAR'}:
    neid = len(eids)
    eid0 = eids[0]
    elem0 = model.elements[eid0]
    nnode = len(elem0.nodes)
    h5_dict = {
        'EID': Int64Col(pos=0),  # Signed 64-bit integer
        'PID': Int64Col(pos=1),
        'N': Int64Col(size=(nnode,), pos=2),
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
        'PROD': write_h5_prod,
        'PSHELL': write_h5_pshell,
        'PCOMP': write_h5_pcomp,
    }
    for prop_type, pids in properties_by_type.items():
        if prop_type in cls_mapper:
            write_h5_func = cls_mapper[prop_type]
            write_h5_func(model, pids, h5file, property_group)
        else:
            log.warning(f'skipping h5 property {prop_type!r}')


def write_h5_materials(self, h5file: File, input_group):
    model = self.model
    log = model.log
    materials_by_type = defaultdict(list)
    for mid, mat in model.materials.items():
        materials_by_type[mat.type].append(mid)
    if len(materials_by_type) == 0:
        return

    material_group = h5file.create_group(input_group, 'MATERIAL')
    cls_mapper = {
        #'MAT1': write_h5_mat1,
    }
    for mat_type, mids in materials_by_type.items():
        if mat_type in cls_mapper:
            write_h5_func = cls_mapper[mat_type]
            write_h5_func(model, mids, h5file, material_group)
        else:
            log.warning(f'skipping h5 material {mat_type!r}')


def write_h5_grid(model: BDF, h5file: File, node_group):
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
    for nid, node in sorted(model.nodes):
        psi = 0 if ps == '' else int(ps)
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


def write_h5_coord(model: BDF, h5file: File, coord_group):
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
    coord_id = []
    coord_type_int = []

    cids_by_type = defaultdict(list)
    icoords_by_type = defaultdict(list)
    for cid, coord in sorted(model.coords.items()):
        if cid == 0:
            continue
        class_name = coord.type
        icoordi = int(class_name[-2])
        coord_type_inti = class_name[-1]
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


def write_h5_prod(model: BDF, pids: list[int],
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

    table = h5file.create_table(property_group, 'PROD', h5_dict)
    arr = np.empty(nprop, dtype=table.dtype)
    arr["PID"] = pids
    for i, pid in enumerate(pids):
        prop = model.properties[pid]
        arr["MID"][i] = prop.mid
        arr["A"][i] = prop.A
        arr["J"][i] = prop.J
        arr["NSM"][i] = prop.nsm
    arr["DOMAIN_ID"] = np.ones(nprop, dtype='int64')
    table.append(arr)

def write_h5_pshell(model: BDF, pids: list[int],
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
        mids = prop.mids
        arr["MID1"][i] = mids[0]  # self.material_id[:, 0]
        arr["MID2"][i] = mids[0]  # self.material_id[:, 1]
        arr["MID3"][i] = mids[0]  # self.material_id[:, 2]
        arr["MID4"][i] = mids[0]  # self.material_id[:, 3]
        arr["T"][i] = prop.t # self.t
        arr["BK"][i] = prop.twelveIt3 # self.twelveIt3
        arr["TS"][i] = prop.tst # self.tst
        arr["Z1"][i] = prop.z1 # self.z[:, 0]
        arr["Z2"][i] = prop.z2 # self.z[:, 1]
        arr["NSM"][i] = prop.nsm # self.nsm
    arr["DOMAIN_ID"] = np.ones(nprop, dtype='int64')
    table.append(arr)

def write_h5_pcomp(model: BDF, pids: list[int],
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
    nlayer = len(model.pcomp.material_id)
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
        sout_int.extend(prop.sout_int)
        material_id.extend(prop.material_id)
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
    arr = np.empty(nlayer, dtype=ply_table.dtype)
    arr["MID"] = material_id
    arr["T"] = thickness
    arr["THETA"] = theta
    arr["SOUT"] = sout_int
    ply_table.append(arr)
