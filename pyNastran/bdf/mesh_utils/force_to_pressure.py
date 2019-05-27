from collections import defaultdict

from pyNastran.bdf.bdf import read_bdf, BDF, FORCE, PLOAD4


def force_to_pressure(bdf_filename, bdf_filename_out=None):
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

    forces = defaultdict(float)
    for load_id, load in model.loads.items():
        for loadi in loads:
            if loadi.type == 'FORCE':
                loadi = FORCE(sid, node, cid, mag, xyz)
                loadi.get

                #if load.node_id not in nids:
                    #continue
                if load.Cid() != 0:
                    cp = load.cid
                    forcei = load.mag * cp.transform_vector_to_global(load.xyz) * scale
                    raise NotImplementedError()
                else:
                    forcei = load.mag * load.xyz * scale

                elem_count = nid_elem_count[load.node]
                f /= elem_count
                for eid in nid_elem_map:
                    forces[eid] += forcei
                #node = self.Node(load.node_id)
                #r = xyz[node.nid] - p
                #m = cross(r, f)
                #F += f
                #M += m
            else:
                raise NotImplementedError(loadi)

    #pressures = {}
    model.cross_reference()
    model.loads = {}
    with open('pressures.out', 'w') as pressure_file:
        for eid, press in forces.items():
            eids = [eid]
            forcei = forces[eid]
            elem = model.elements[eid]
            area = elem.Area()
            pressures = [pressure, pressure, pressure, pressure]
            pload4 = PLOAD4(sid, eids, pressures,
                            g1=None, g34=None, cid=0, nvector=None,
                            surf_or_line='SURF', line_load_dir='NORM', comment='')
            #pressure_file.write(pload4.write_card(size=8, is_double=False))
            model._add_load_object(pload4)

    if bdf_filename_out:
        model.write_bdf(bdf_filename_out)
    return model
