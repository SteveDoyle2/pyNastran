# coding: utf-8
"""
This file defines:
  - write_dict(bdf_file, my_dict, size, is_double, is_long_ids)
  - write_list(bdf_file, my_list, size, is_double, is_long_ids)
  - write_aero_in_flutter, write_aero_in_gust = find_aero_location(model)
  - ptype_to_pid, property_type_to_property_class, ...
        properties_by_class = get_properties_by_element_type(model)

"""
from __future__ import annotations
from collections import defaultdict
from typing import TextIO, Any, TYPE_CHECKING

if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf import BDF
    from pyNastran.bdf.cards.base_card import BaseCard


def csv_card_8(bdf_file: TextIO, card: BaseCard) -> None:
    fields = card.raw_fields()
    fields_str = (f'{field},' for field in fields)
    out = ''
    for i, field in enumerate(fields_str):
        out += field
        if i > 0 and i % 8 == 0:
            out += '\n'
    lines = [line.strip(',') for line in out.strip(',\n').split('\n')]
    for line in lines:
        sline = line.split(',')
        assert len(sline) <= 9, lines
    bdf_file.write(card.comment + out)


def csv_card_16(bdf_file: TextIO, card: BaseCard) -> None:
    """not tested"""
    fields = card.raw_fields()
    fields[0] += '*'
    fields_str = (f'{field},' for field in fields)
    out = ''
    for i, field in enumerate(fields_str):
        out += field
        if i > 0 and i % 4 == 0:
            out += '\n,*,'
    n = len(fields) - 1


    if n % 8 == 0:
        # full card; no blank line
        out = out.strip('\n*,') + '\n'
    elif n % 8 < 3:
        # single line of card; add a blank line
        out += '\n,*,\n'
    else:
        # end of strip; no blank line
        out = out.rstrip('*\n') + '\n'

    lines = [line.strip(',*') for line in out.strip(',*\n').split('\n')]
    for line in lines:
        sline = line.split(',')
        assert len(sline) <= 5, lines
    bdf_file.write(card.comment + out)


def write_dict(bdf_file: TextIO, my_dict: dict[int, BaseCard],
               size: int, is_double: bool, is_csv: bool,
               is_long_ids: bool) -> None:
    """writes a dictionary that may require long format"""
    if is_csv:
        if size == 16 or is_long_ids:
            for nid, node in sorted(my_dict.items()):
                csv_card_16(bdf_file, node)
        else:
            for nid, node in sorted(my_dict.items()):
                csv_card_8(bdf_file, node)
    else:
        if is_long_ids:
            for (unused_nid, node) in sorted(my_dict.items()):
                bdf_file.write(node.write_card_16(is_double))
        else:
            for (unused_nid, node) in sorted(my_dict.items()):
                bdf_file.write(node.write_card(size, is_double))


def write_list(bdf_file: TextIO, my_list: list[BaseCard],
               size: int, is_double: bool, is_csv: bool, is_long_ids: bool):
    if is_csv:
        if size == 16 or is_long_ids:
            for card in my_list:  # list
                csv_card_16(bdf_file, card)
        else:
            for card in my_list:  # list
                csv_card_8(bdf_file, card)
    else:
        for card in my_list:  # list
            bdf_file.write(card.write_card(size, is_double))


def find_aero_location(model: BDF) -> tuple[bool, bool]:
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


def get_properties_by_element_type(model: BDF) -> tuple[dict[str, list[str]],
                                                        dict[str, Any],
                                                        dict[str, Any]]:
    """helper for ``_write_properties_by_element_type``"""
    propertys_class_to_property_types = {
        # prop_class -> property types
        'spring': ['PELAS', 'PELAST'],
        'damper': ['PDAMP', 'PDAMPT'],
        'rod': ['PROD', 'PTUBE'],
        'bar': ['PBAR', 'PBARL', 'PBRSECT'],
        'beam': ['PBEAM', 'PBEAML', 'PBMSECT'],
        'bush': ['PBUSH', 'PBUSH1D', 'PBUSH2D'],
        'shell': ['PSHEAR', 'PSHELL', 'PCOMP', 'PCOMPG'],
        'solid': ['PSOLID', 'PLSOLID'],
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
