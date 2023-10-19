from __future__ import annotations
from typing import TYPE_CHECKING
import numpy as np

from .utils import get_group_name, get_attributes
if TYPE_CHECKING:  # pragma: no cover
    #from pyNastran.dev.bdf_vectorized3.cards.aero.aero import CAERO1
    from pyNastran.dev.bdf_vectorized3.cards.elements.mass import CONM2
    #from pyNastran.dev.bdf_vectorized3.cards.elements.rod import CONROD
    #from pyNastran.dev.bdf_vectorized3.cards.elements.bar import CBAR
    #from pyNastran.dev.bdf_vectorized3.cards.elements.beam import CBEAM
    from pyNastran.dev.bdf_vectorized3.cards.elements.shell import CQUAD4
    from pyNastran.dev.bdf_vectorized3.bdf import BDF
    from tables import Group


def load_h5_element(model: BDF, input_group: Group):
    name = '???'
    for h5_element in input_group._f_iter_nodes():
        if h5_element._c_classid == 'GROUP':
            class_name = get_group_name(h5_element)
            print(f'skipping {class_name} in _load_h5_element')
            continue

        element_name = h5_element.name
        #print(f'checking {element_name}')
        data = h5_element.read()
        element_id = data['EID']
        nelements = len(element_id)

        #attributes = get_attributes(data)
        #assert len(attributes) == 1, attributes
        #version = attributes['version'][0]  # 20200

        skip_elements = {
            'CAERO2', 'CAERO3', 'CAERO4', 'CAERO5',
        }
        if element_name in skip_elements:
            if hasattr(model, element_name.lower()):
                raise NotImplementedError(element_name)
            model.log.warning(f'skipping {element_name} in _load_h5_element')
            continue

        elif element_name in {'CQUAD4', 'CQUAD8', 'CQUADR',
                              'CTRIA3', 'CTRIA6', 'CTRIAR'}:
            #model.log.info(f'loading {element_name}')
            #print(f'loading {element_name}')
            # CQUAD4 -> self.cquad4
            elem = getattr(model, element_name.lower())
            _load_h5_shell(elem, element_id, data)
        elif element_name == 'CONM2':
            elem = _load_h5_conm2(model, data, element_id)
            #elem._save()
        elif element_name == 'CROD':
            elem = model.crod
            _load_h5_crod(elem, data, element_id)
        elif element_name == 'CBAR':
            elem = model.cbar
            _load_h5_cbar(elem, data, element_id)
        elif element_name == 'CBEAM':
            elem = model.cbeam
            _load_h5_cbeam(elem, data, element_id)

        elif element_name in {'CHEXA', 'CPENTA', 'CTETRA'}:
            elem = getattr(model, element_name.lower())
            property_id = data['PID']
            nodes = data['G']
            elem._save(element_id, property_id, nodes)
            #elem.write()
        elif element_name == 'CELAS1':
            elem = model.celas1
            read_celas1(name, data, elem)
        elif element_name == 'CELAS2':
            elem = model.celas2
            read_celas2(name, data, elem)
        elif element_name == 'CELAS3':
            elem = model.celas3
            read_celas3(name, data, elem)
        elif element_name == 'CELAS4':
            elem = model.celas4
            read_celas4(name, data, elem)
        elif element_name == 'CDAMP1':
            elem = model.cdamp1
            read_cdamp1(name, data, elem)
        elif element_name == 'CDAMP2':
            elem = model.cdamp2
            read_cdamp2(name, data, elem)
        elif element_name == 'CDAMP3':
            elem = model.cdamp3
            read_cdamp3(name, data, elem)
        elif element_name == 'CDAMP4':
            elem = model.cdamp4
            read_cdamp4(name, data, elem)
        else:
            model.log.warning(f'skipping {element_name} in _load_h5_element')
            #print(f'skipping {element_name} in _load_h5_element')
            continue
        elem.domain_id = data['DOMAIN_ID']
        elem.n = nelements
        elem.write()
    return
    #for (element_name, h5_element) in input_group._v_children.items():
        #if element_name in {'RBE2', 'RBE3'}:
            #print(f'skipping {element_name}')
            #model.log.warning(f'skipping {element_name}')
            #continue

        #data = h5_element.read()

        #data.dtype
        #x = 1
    #y = 1

def _load_h5_conm2(model: BDF, data, element_id) -> CONM2:
    elem = model.conm2
    # dtype([('EID', '<i8'), ('G', '<i8'), ('CID', '<i8'), ('M', '<f8'),
    #        ('X1', '<f8'), ('X2', '<f8'), ('X3', '<f8'), ('I1', '<f8'), ('I2', '<f8', (2,)),
    #        ('I3', '<f8', (3,)), ('DOMAIN_ID', '<i8')])
    nelements = len(element_id)
    elem.element_id = element_id
    elem.node_id = data['G']
    elem.coord_id = data['CID']
    elem._mass = data['M']
    elem.xyz_offset = np.stack([data['X1'], data['X2'], data['X3']], axis=1)
    elem.inertia = np.hstack([
        data['I1'].reshape(nelements, 1),
        data['I2'].reshape(nelements, 2),
        data['I3'].reshape(nelements, 3),
    ])
    return elem

def _load_h5_crod(elem: CROD, data, element_id):
    """
    array([(14, 3, [20, 24], 1), (15, 3, [21, 25], 1)],
      dtype=[('EID', '<i8'), ('PID', '<i8'), ('G', '<i8', (2,)), ('DOMAIN_ID', '<i8')])
    """
    nelements = len(element_id)
    element_id = element_id
    property_id = data['PID']
    node_id = data['G']
    elem._save(element_id, property_id, node_id)
    return elem

def _load_h5_cbar(elem: CBAR, data, element_id):
    """
    array([(13, 1, 19, 23, 1, 0., 1., 0., 0, 0, 0, 0., 0., 0., 0., 0., 0., 1)],
      dtype=[('EID', '<i8'), ('PID', '<i8'), ('GA', '<i8'), ('GB', '<i8'),
      ('FLAG', '<i8'), ('X1', '<f8'), ('X2', '<f8'), ('X3', '<f8'),
      ('GO', '<i8'), ('PA', '<i8'), ('PB', '<i8'),
      ('W1A', '<f8'), ('W2A', '<f8'), ('W3A', '<f8'),
      ('W1B', '<f8'), ('W2B', '<f8'), ('W3B', '<f8'),
      ('DOMAIN_ID', '<i8')])
    """
    nelements = len(element_id)
    element_id = element_id
    property_id = data['PID']
    nodes = np.column_stack([data['GA'], data['GB'], ])
    wa = np.column_stack([data['W1A'], data['W2A'], data['W3A'], ])
    wb = np.column_stack([data['W1B'], data['W2B'], data['W3B'], ])
    g0 = data['GO']
    x = np.column_stack([data['X1'], data['X2'], data['X3'], ])

    pa = data['PA']
    pb = data['PB']
    flag = data['FLAG']

    offt = np.full(len(element_id), '', dtype='|U3')
    for uflag in np.unique(flag):
        iflag = (uflag == flag)
        if flag == 1:
            offt[iflag] = 'GGG'
        else:
            raise NotImplementedError(flag)

    #node_id = data['G']
    elem._save(element_id, property_id, nodes, g0, x,
               offt, pa, pb, wa, wb)
    #elem.coord_id = data['CID']
    #elem._mass = data['M']
    #elem.xyz_offset = np.stack([data['X1'], data['X2'], data['X3']], axis=1)
    #elem.inertia = np.hstack([
        #data['I1'].reshape(nelements, 1),
        #data['I2'].reshape(nelements, 2),
        #data['I3'].reshape(nelements, 3),
    #])
    return elem

def _load_h5_cbeam(elem: CBEAM, data, element_id):
    """
    array([(12, 5, 18, 22, 0, 0, [0., 1., 0.], 0, 1, 0, 0, [0., 0., 0.], [0., 0., 0.], 1)],
    dtype=[('EID', '<i8'), ('PID', '<i8'), ('GA', '<i8'), ('GB', '<i8'),
    ('SA', '<i8'), ('SB', '<i8'), ('X', '<f8', (3,)), ('G0', '<i8'),
    ('F', '<i8'), ('PA', '<i8'), ('PB', '<i8'),
    ('WA', '<f8', (3,)), ('WB', '<f8', (3,)), ('DOMAIN_ID', '<i8')])
    """
    nelements = len(element_id)
    element_id = element_id
    property_id = data['PID']
    #node_id = data['G']

    nodes = np.column_stack([data['GA'], data['GB'], ])
    wa = data['WA']
    wb = data['WB']
    g0 = data['G0']
    x = data['X']

    pa = data['PA']
    pb = data['PB']
    sa = data['SA']
    sb = data['SB']
    flag = data['F']

    offt = np.full(len(element_id), '', dtype='|U3')
    for uflag in np.unique(flag):
        iflag = (uflag == flag)
        if flag == 1:
            offt[iflag] = 'GGG'
        else:
            raise NotImplementedError(flag)

    bit = np.full(len(element_id), -1, dtype='int32')
    elem._save(element_id, property_id, nodes,
               offt, bit,
               g0, x,
               pa, pb, wa, wb, sa, sb)
    #elem.coord_id = data['CID']
    #elem._mass = data['M']
    #elem.xyz_offset = np.stack([data['X1'], data['X2'], data['X3']], axis=1)
    #elem.inertia = np.hstack([
        #data['I1'].reshape(nelements, 1),
        #data['I2'].reshape(nelements, 2),
        #data['I3'].reshape(nelements, 3),
    #])
    return elem

def _load_h5_shell(elem: CQUAD4, element_id, data):
    # dtype([('EID', '<i8'), ('PID', '<i8'), ('G', '<i8', (4,)), ('THETA', '<f8'),
    #        ('ZOFFS', '<f8'), ('TFLAG', '<i8'), ('T', '<f8', (4,)),
    #        ('MCID', '<i8'), ('DOMAIN_ID', '<i8')])
    property_id = data['PID']
    nodes = data['G']
    theta = data['THETA']
    element_id = data['EID']
    mcid = data['MCID']
    zoffset = data['ZOFFS']
    tflag = data['TFLAG']
    T = data['T']
    elem._save(element_id, property_id, nodes,
               zoffset=zoffset, theta=theta, mcid=mcid, tflag=tflag, T=T)


def read_celas1(name: str, group: h5py._hl.dataset.Dataset, elem: CELAS1) -> None:
    #('EID', 'PID', 'G1', 'G2', 'C1', 'C2', 'DOMAIN_ID')
    element_id = group['EID']
    property_id = group['PID']
    G1 = group['G1']
    G2 = group['G2']
    C1 = group['C1']
    C2 = group['C2']
    nodes = np.column_stack([G1, G2])
    components = np.column_stack([C1, C2])
    #assert NIDS.shape[1] == 2, NIDS.shape
    DOMAIN_ID = group['DOMAIN_ID']
    #for eid, pid, nids, c1, c2 in zip(EID, PID, NIDS, C1, C2):
        #obj = geom_model.add_celas1(eid, pid, nids, c1=c1, c2=c2, comment='')
        #obj.validate()
    elem._save(element_id, property_id, nodes, components)

def read_celas2(name: str, group: h5py._hl.dataset.Dataset, elem: CELAS2) -> None:
    #('EID', 'K', 'G1', 'G2', 'C1', 'C2', 'DOMAIN_ID')
    #dtype=[('EID', '<i8'), ('K', '<f8'), ('G1', '<i8'), ('G2', '<i8'), ('C1', '<i8'), ('C2', '<i8'),
           #('GE', '<f8'), ('S', '<f8'), ('DOMAIN_ID', '<i8')])
    element_id = group['EID']
    k = group['K']
    G1 = group['G1']
    G2 = group['G2']
    C1 = group['C1']
    C2 = group['C2']
    s = group['S']
    ge = group['GE']
    nodes = np.column_stack([G1, G2])
    components = np.column_stack([C1, C2])
    #assert NIDS.shape[1] == 2, NIDS.shape
    DOMAIN_ID = group['DOMAIN_ID']
    #for eid, k, nids, c1, c2 in zip(EID, K, NIDS, C1, C2):
        #obj = geom_model.add_celas2(eid, k, nids, c1=c1, c2=c2, comment='')
        #obj.validate()
    elem._save(element_id, nodes, components, k, ge, s)

def read_celas3(name: str, group: h5py._hl.dataset.Dataset, elem: CELAS3) -> None:
    element_id = group['EID']
    property_id = group['PID']
    G1 = group['S1']
    G2 = group['S2']
    spoints = np.stack([G1, G2], axis=1)
    #assert NIDS.shape[1] == 2, NIDS.shape
    DOMAIN_ID = group['DOMAIN_ID']
    #for eid, pid, nids in zip(EID, PID, NIDS):
        #obj = geom_model.add_celas3(eid, pid, nids, comment='')
        #obj.validate()
    elem._save(element_id, property_id, spoints)

def read_celas4(name: str, group: h5py._hl.dataset.Dataset, elem: CELAS4) -> None:
    element_id = group['EID']
    k = group['K']
    G1 = group['S1']
    G2 = group['S2']
    spoints = np.stack([G1, G2], axis=1)
    #assert spoints.shape[1] == 2, spoints.shape
    DOMAIN_ID = group['DOMAIN_ID']
    #for eid, k, nids in zip(EID, K, NIDS):
        #obj = geom_model.add_celas4(eid, k, nids)
        #obj.validate()
    elem._save(element_id, k, spoints)

def read_cdamp1(name: str, group: h5py._hl.dataset.Dataset, elem: CDAMP1) -> None:
    #('EID', 'PID', 'G1', 'G2', 'C1', 'C2', 'DOMAIN_ID')
    element_id = group['EID']
    property_id = group['PID']
    G1 = group['G1']
    G2 = group['G2']
    C1 = group['C1']
    C2 = group['C2']
    nodes = np.column_stack([G1, G2])
    components = np.column_stack([C1, C2])
    #assert NIDS.shape[1] == 2, NIDS.shape
    DOMAIN_ID = group['DOMAIN_ID']
    #for eid, pid, nids, c1, c2 in zip(EID, PID, NIDS, C1, C2):
        #obj = geom_model.add_cdamp1(eid, pid, nids, c1=c1, c2=c2, comment='')
        #obj.validate()
    elem._save(element_id, property_id, nodes, components)

def read_cdamp2(name: str, group: h5py._hl.dataset.Dataset, elem: CDAMP2) -> None:
    #('EID', 'B', 'G1', 'G2', 'C1', 'C2', 'DOMAIN_ID')
    element_id = group['EID']
    b = group['B']
    G1 = group['G1']
    G2 = group['G2']
    C1 = group['C1']
    C2 = group['C2']
    nodes = np.column_stack([G1, G2])
    components = np.column_stack([C1, C2])
    assert nodes.shape[1] == 2, nodes.shape
    DOMAIN_ID = group['DOMAIN_ID']
    #for eid, b, nids, c1, c2 in zip(EID, B, NIDS, C1, C2):
        #obj = geom_model.add_cdamp2(eid, b, nids, c1=c1, c2=c2, comment='')
        #obj.validate()
    elem._save(element_id, nodes, components, b)

def read_cdamp3(name: str, group: h5py._hl.dataset.Dataset, elem: CDAMP3) -> None:
    element_id = group['EID']
    property_id = group['PID']
    G1 = group['S1']
    G2 = group['S2']
    spoints = np.stack([G1, G2], axis=1)
    assert spoints.shape[1] == 2, spoints.shape
    DOMAIN_ID = group['DOMAIN_ID']
    elem._save(element_id, property_id, spoints)
    #for eid, pid, nids in zip(EID, PID, NIDS):
        #obj = geom_model.add_cdamp3(eid, pid, nids, comment='')
        #obj.validate()

def read_cdamp4(name: str, group: h5py._hl.dataset.Dataset, elem: CDAMP4) -> None:
    element_id = group['EID']
    b = group['B']
    G1 = group['S1']
    G2 = group['S2']
    spoints = np.column_stack([G1, G2])
    assert spoints.shape[1] == 2, spoints.shape
    DOMAIN_ID = group['DOMAIN_ID']
    #for eid, b, nids in zip(EID, B, NIDS):
        #obj = geom_model.add_cdamp4(eid, b, nids)
        #obj.validate()
    elem._save(element_id, b, spoints)
