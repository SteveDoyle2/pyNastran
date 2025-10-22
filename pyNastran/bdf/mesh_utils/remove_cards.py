from pyNastran.bdf.bdf import BDF


def remap_cards(model: BDF,
                cross_reference: bool=True,
                remap_nodes: bool=True,
                remap_elements: bool=True,
                remap_properties: bool=True,
                remap_materials: bool=True):
    """
    Remap cards after a cross-reference; works on an uncross-referenced model

    Parameters
    ----------
    model : BDF()
        the BDF object
    cross_reference : bool; default=True
        cross_reference once the remapping is done
    remap_nodes : bool; default=True
        'x' may be nodes, elements, properties, materials

    """
    model.uncross_reference()

    if remap_nodes:
        nodes = {}
        for node in model.nodes.values():
            nid = node.nid
            nodes[nid] = node
        model.nodes = nodes

    if remap_elements:
        elements = {}
        for element in model.elements.values():
            eid = element.eid
            elements[eid] = element
        model.elements = elements

    if remap_properties:
        properties = {}
        for prop in model.properties.values():
            pid = prop.pid
            properties[pid] = prop
        model.properties = properties

    if remap_materials:
        materials = {}
        for material in model.materials.values():
            mid = material.mid
            materials[mid] = material
        model.materials = materials
    if cross_reference:
        model.cross_reference()


def delete_properties(bdf_model: BDF, property_types_to_save=None):
    """early version of way to delete specific property cards"""
    pids_to_delete = set()
    if property_types_to_save:
        for pid, prop in bdf_model.properties.items():
            ptype = prop.type
            if ptype not in property_types_to_save:
                pids_to_delete.add(pid)
                bdf_model._type_to_id_map[ptype].remove(pid)
    for pid in pids_to_delete:
        del bdf_model.properties[pid]


def delete_elements(bdf_model: BDF, element_types_to_save=None):
    """early version of way to delete specific element cards"""
    eids_to_delete = set()
    if element_types_to_save:
        for eid, element in bdf_model.elements.items():
            etype = element.type
            if etype not in element_types_to_save:
                eids_to_delete.add(eid)
                bdf_model._type_to_id_map[etype].remove(eid)
    for eid in eids_to_delete:
        del bdf_model.elements[eid]

#def delete_forces(bdf_model, eids_to_delete=None):
    #"""early version of way to delete specific force cards"""
    #eids_to_delete = set()
    #loads = {}
    #for load_id, load_set in bdf_model.loads.items():
        #load_set = []
        #for load in load_set:
            #print(load)
    #if element_types_to_save:
        #for eid, element in bdf_model.elements.items():
            #if element.type not in element_types_to_save:
                #eids_to_delete.add(eid)
    #for eid in eids_to_delete:
        #del bdf_model.elements[eid]
