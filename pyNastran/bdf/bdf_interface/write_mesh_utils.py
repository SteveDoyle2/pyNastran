# coding: utf-8
"""
This file defines:
  - write_dict(bdf_file, my_dict, size, is_double, is_long_ids)
  - write_aero_in_flutter, write_aero_in_gust = find_aero_location(model)
  - ptype_to_pid, property_type_to_property_class, ...
        properties_by_class = get_properties_by_element_type(model)

"""
from __future__ import annotations
from collections import defaultdict
from typing import List, Dict, Tuple, Any, TYPE_CHECKING

if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf import BDF


def write_dict(bdf_file, my_dict: Dict[int, Any], size: int,
               is_double: bool, is_long_ids: bool) -> None:
    """writes a dictionary that may require long format"""
    if is_long_ids:
        for (unused_nid, node) in sorted(my_dict.items()):
            bdf_file.write(node.write_card_16(is_double))
    else:
        for (unused_nid, node) in sorted(my_dict.items()):
            bdf_file.write(node.write_card(size, is_double))


def find_aero_location(model: BDF) -> Tuple[bool, bool]:
    """Determines where the AERO card should be written"""
    write_aero_in_flutter = False
    write_aero_in_gust = False
    if model.aero:
        if model.flfacts or model.flutters or model.mkaeros:
            write_aero_in_flutter = True
        elif model.gusts:
            write_aero_in_gust = True
        else:
            # an AERO card exists, but no FLUTTER, FLFACT, MKAEROx or GUST card
            write_aero_in_flutter = True
    return write_aero_in_flutter, write_aero_in_gust

def get_properties_by_element_type(model: BDF) -> Tuple[Dict[str, List[str]],
                                                        Dict[str, Any],
                                                        Dict[str, Any]]:
    """helper for ``_write_properties_by_element_type``"""
    propertys_class_to_property_types = {
        # prop_class -> property types
        'spring' : ['PELAS', 'PELAST'],
        'damper' : ['PDAMP', 'PDAMPT'],
        'rod' : ['PROD', 'PTUBE'],
        'bar' : ['PBAR', 'PBARL', 'PBRSECT'],
        'beam' : ['PBEAM', 'PBEAML', 'PBMSECT'],
        'bush' : ['PBUSH', 'PBUSH1D', 'PBUSH2D'],
        'shell' : ['PSHEAR', 'PSHELL', 'PCOMP', 'PCOMPG'],
        'solid' : ['PSOLID', 'PLSOLID'],
    }

    property_type_to_property_class = {
        #'other' : [],
    }
    # the inverse of propertys_class_to_property_types
    for prop_class, prop_types in propertys_class_to_property_types.items():
        for prop_type in prop_types:
            property_type_to_property_class[prop_type] = prop_class

    #if is_properties:

    # put each property object into a class (e.g., CQUAD4 -> PCOMP)
    properties_by_class = defaultdict(list)
    prop_groups = (model.properties, model.pelast, model.pdampt, model.pbusht)
    for properties in prop_groups:
        for unused_pid, prop in properties.items():
            prop_class = property_type_to_property_class[prop.type]
            properties_by_class[prop_class].append(prop)
    return propertys_class_to_property_types, property_type_to_property_class, properties_by_class
