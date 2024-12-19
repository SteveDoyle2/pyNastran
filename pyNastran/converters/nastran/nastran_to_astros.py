from pyNastran.converters.astros.astros import Astros
from pyNastran.converters.astros.astros_cards import Part

def nastran_to_astros(nastran_model, astros_filename_out):
    """
    Handles CTRIA3, CQUAD4
    """
    log = nastran_model.log
    model = Astros(log=log, debug=True)
    node_sets = None
    element_sets = None
    solid_sections = None

    ctria3s = []
    cquad4s = []
    for eid, elem in nastran_model.elements.items():
        if elem.type == 'CTRIA3':
            node_ids = elem.nodes
            ctria3 = [eid] + node_ids
            ctria3s.append(ctria3)
        elif elem.type == 'CQUAD4':
            node_ids = elem.nodes
            cquad4 = [eid] + node_ids
            cquad4s.append(cquad4)
        else:
            pass

    element_types = {
        'cpe3' : ctria3s,
        'cpe4' : cquad4s,
    }
    name = 'model'
    icd_transform, icp_transform, xyz_cp, nid_cp_cd = nastran_model.get_displacement_index_xyz_cp_cd(
        fdtype='float64', sort_ids=True)
    nids = nid_cp_cd[:, 0]
    xyz_cid0 = nastran_model.transform_xyzcp_to_xyz_cid(
        xyz_cp, nids, icp_transform, cid=0, in_place=False)
    nodes = xyz_cid0
    part = Part(name, nids, nodes, element_types, node_sets, element_sets,
                solid_sections, log)
    model.parts[name] = part
    model.write(astros_filename_out)
    return model