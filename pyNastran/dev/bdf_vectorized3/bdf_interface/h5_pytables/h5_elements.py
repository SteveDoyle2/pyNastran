from __future__ import annotations
from typing import TYPE_CHECKING
import numpy as np

from .utils import get_group_name, get_attributes
if TYPE_CHECKING:
    #from pyNastran.dev.bdf_vectorized3.cards.aero.aero import CAERO1
    from pyNastran.dev.bdf_vectorized3.cards.elements.mass import CONM2
    #from pyNastran.dev.bdf_vectorized3.cards.elements.rod import CONROD
    #from pyNastran.dev.bdf_vectorized3.cards.elements.bar import CBAR
    #from pyNastran.dev.bdf_vectorized3.cards.elements.beam import CBEAM
    from pyNastran.dev.bdf_vectorized3.cards.elements.shell import CQUAD4
    from pyNastran.dev.bdf_vectorized3.bdf import BDF
    from tables import Group


def load_h5_element(model: BDF, input_group: Group):
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
        elif element_name in {'CHEXA', 'CPENTA', 'CTETRA'}:
            elem = getattr(model, element_name.lower())
            property_id = data['PID']
            nodes = data['G']
            elem._save(element_id, property_id, nodes)
            #elem.write()
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

def _load_h5_conm2(model, data, element_id) -> CONM2:
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

