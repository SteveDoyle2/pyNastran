from __future__ import annotations
from typing import Any, TYPE_CHECKING
import numpy as np

if TYPE_CHECKING:
    from tables import Group, Node


def cast_encoding_strip(myarray: np.ndarray, encoding: str) -> np.ndarray:
    """casts a byte array to a unicode string array"""
    out = np.array([val.strip().decode(encoding) for val in myarray])
    return out

def cast_encoding_strip0(myarray: np.ndarray, encoding: str) -> np.ndarray:
    """Casts a byte array to a unicode string array.  This one uses
    an (N,) array of tuples of byte strings..."""
    out = np.array([val[0].strip().decode(encoding) for val in myarray])
    return out

def cast_encoding(myarray: np.ndarray, encoding: str) -> np.ndarray:
    """casts a byte array to a unicode string array"""
    out = np.array([val.decode(encoding) for val in myarray])
    return out

def get_attributes(node: Node) -> dict[str, Any]:
    """
    Gets a dictionary of minor parameters that are attached to the pytables Node
    (e.g., version, date)
    """
    attributes = {}

    attrs = node._v_attrs
    for name in attrs._v_attrnames:
        attr = attrs[name]
        if isinstance(attr, bytes):
            attributes[name] = attr.decode('latin1')
        else:
            attributes[name] = attr
    return attributes

def get_group_name(h5_group: Group) -> str:
    """'/NASTRAN/INPUT//DYNAMIC' -> DYNAMIC"""
    class_name = h5_group._v_pathname.rsplit('/', 1)[1]  # PBARL
    return class_name

def print_group(group: Group, indent: str='') -> str:
    msg = ''
    for h5_element in group._f_iter_nodes():
        class_id = h5_element._c_classid
        if class_id == 'GROUP':
            group_name = get_group_name(group)
            msg += f'{indent}Group: {group_name}\n'
            indent += '    '
            print_group(group, indent=indent)
            continue
        elif class_id == 'TABLE':
            pass
        else:
            raise RuntimeError(class_id)
        name = h5_element.name
        msg += f'{indent}Table: {name}\n'
    print(msg)
    return msg
