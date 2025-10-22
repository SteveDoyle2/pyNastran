from collections import defaultdict
from typing import Optional

import numpy as np
from pyNastran.bdf.bdf import read_bdf, BDF, FORCE, PLOAD4
from pyNastran.utils import PathLike


def force_to_pressure(bdf_filename: PathLike | BDF,
                      bdf_filename_out: Optional[PathLike]=None,
                      clear_model: bool=False):
    """converts FORCE cards to PLOAD4s for a shell model"""
    if isinstance(bdf_filename, BDF):
        model = bdf_filename
    else:
        model = read_bdf(bdf_filename, validate=True, xref=False, punch=False,
                         encoding=None, log=None, debug=True,
                         mode='msc')

    #if 0:
        ##card_types = ['CQUAD4', 'CTRIA3']
        #card_ids_map = model.get_card_ids_by_card_types(card_types=None,
                                                        #reset_type_to_slot_map=False,
                                                        #stop_on_missing_card=False)

        #for eid in card_ids_map['CQUAD4']:
            #elem = model.elements[eid]
            ##for nid in elem.node_ids:
            ##raise NotImplementedError(elem)
        #for eid in card_ids_map['CTRIA3']:
            #elem = model.elements[eid]
            #raise NotImplementedError(elem)

    nid_elem_count = defaultdict(int)
    nid_elem_map = defaultdict(list)
    for eid, elem in model.elements.items():
        for nid in elem.nodes:
            nid_elem_count[nid] += 1
            nid_elem_map[nid].append(eid)
    #model.cross_reference()

    forces = {}
    for load_id, loads in model.loads.items():
        forces[load_id] = defaultdict(float)
        for load in loads:
            if load.type == 'FORCE':
                scale = load.mag
                # cid: 0
                # cid_ref: None
                # comment: ''
                # mag: 1.0
                # node: 1
                # node_id: 1
                # node_ref: None
                # scaled_vector: array([0., 0., 1.])
                # sid: 1
                # type_card: FORCE
                # xyz: array([0., 0., 1.])

                #loadi = FORCE(sid, node, cid, mag, xyz)

                #if load.node_id not in nids:
                    #continue
                node_id = load.node
                if load.Cid() != 0:  # pragma: no cover
                    cp = load.cid
                    forcei = load.mag * cp.transform_vector_to_global(load.xyz) * scale
                    raise NotImplementedError()
                else:
                    forcei = load.mag * load.xyz * scale

                elem_count = nid_elem_count[node_id]
                forcei /= elem_count
                for eid in nid_elem_map[node_id]:
                    forces[load_id][eid] += forcei
                #node = self.Node(load.node_id)
                #r = xyz[node.nid] - p
                #m = np.cross(r, f)
                #F += f
                #M += m
            else:  # pragma: no cover
                #print(load.get_stats())
                raise NotImplementedError(load)

    #pressures = {}
    model.cross_reference()
    model.loads = {}
    add_methods = model._add_methods
    # with open(pressure_filename, 'w') as pressure_file:
    for sid, forcesi in forces.items():
        for eid, press in forcesi.items():
            eids = [eid]
            forcei = forcesi[eid]
            elem = model.elements[eid]
            area = elem.Area()
            normal = elem.Normal()
            #print(f'normal = {normal}')
            #print(f'force  = {forcei}')
            force_normal = forcei * normal
            force_mag = np.linalg.norm(force_normal)
            if force_mag == 0.0:
                continue

            pressure = force_mag / area
            #print(f'force_normal = {force_normal}')
            #print(f'force_mag    = {force_mag}')
            #print(f'area         = {area}')
            #print(f'pressure     = {pressure}')
            assert isinstance(pressure, float), pressure
            pressures = [pressure, pressure, pressure, pressure]
            pload4 = PLOAD4(sid, eids, pressures,
                            g1=None, g34=None, cid=0, nvector=None,
                            surf_or_line='SURF', line_load_dir='NORM', comment='')
            #pressure_file.write(pload4.write_card(size=8, is_double=False))
            add_methods.add_load_object(pload4)

    if clear_model:
        model2 = BDF()
        model2.log = model.log
        model2.loads = model.loads
        model = model2

    if bdf_filename_out:
        model.write_bdf(bdf_filename_out)
    return model
