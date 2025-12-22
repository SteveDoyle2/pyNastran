"""
defines:
    bdf_renumber(bdf_filename, bdf_filename_out, size=8, is_double=False,
                 starting_id_dict=None, round_ids=False, cards_to_skip=None,
                 log=None, debug=False)
    superelement_renumber(bdf_filename, bdf_filename_out=None, size=8, is_double=False,
                          starting_id_dict=None, cards_to_skip=None,
                          log=None, debug=False)

"""
from __future__ import annotations
from itertools import chain
from io import StringIO, IOBase
from typing import Optional, TYPE_CHECKING

import numpy as np

from pyNastran.bdf.bdf import BDF
from pyNastran.utils import PathLike
from pyNastran.utils.numpy_utils import integer_types
from pyNastran.bdf.mesh_utils.bdf_renumber import _write_bdf, _get_bdf_model
if TYPE_CHECKING:  # pragma: no cover
    from cpylog import SimpleLogger


SKIP_ATTRS = {
    'active_filenames', 'include_filenames',
    'allow_duplicate_element_rbe_mass',
    'allow_overwrites_set', 'card_count', 'has_enddata',
    'cards_to_read', 'case_control_deck', 'case_control_lines',
    'executive_control_lines', 'debug', 'is_bdf_vectorized',
    'is_msc', 'is_nx', 'is_optistruct',
    'is_zaero', 'is_mystran',
    'is_superelements',
    'is_strict_card_parser', 'nid_map',
    'dumplines', 'echo', 'force_echo_off',
    'rsolmap_to_str', 'special_cards', 'system_command_lines',
    'read_includes', 'save_file_structure', 'sol', 'sol_iline',
    'use_new_deck_parser', 'wtmass', 'xref_obj',
    'zona', 'punch',
}

def bdf_remove_comments(bdf_filename: PathLike | BDF | StringIO,
                        bdf_filename_out: str,
                        #xref: bool=True,
                        size: int=8, is_double: bool=False,
                        log: Optional[SimpleLogger]=None,
                        debug: bool=False) -> BDF:
    """
    Removes the comments from a BDF

    Parameters
    ----------
    bdf_filename : str / BDF
        str : a bdf_filename (string; supported)
        BDF : a BDF model that has been cross referenced and is
        fully valid (an equivalenced deck is not valid)
    bdf_filename_out : str / None
        str : a bdf_filename to write
        None : don't write the BDF
    size : int; {8, 16}; default=8
        the bdf write precision
    is_double : bool; default=False
        the field precision to write

    Returns
    -------
    model : BDF()
        a cleaned file

    """
    assert size in [8, 16], size
    assert isinstance(is_double, bool), is_double
    cards_to_skip = []
    model = _get_bdf_model(bdf_filename, punch=None,
                           xref=False,
                           cards_to_skip=cards_to_skip,
                           log=log, debug=debug)

    attrs = model.object_attributes()
    for attr in attrs:
        if attr in SKIP_ATTRS:
            continue
        dict_list_scalar = getattr(model, attr)
        if dict_list_scalar is None:
            continue

        if isinstance(dict_list_scalar, list):
            for value in dict_list_scalar:
                # list[card]
                assert hasattr(value, 'comment'), (attr, value)
                value.comment = ''
        elif isinstance(dict_list_scalar, dict):
            for key, value in dict_list_scalar.items():
                if isinstance(value, list):
                    # dict[list[card]]
                    for valuei in value:
                        assert hasattr(valuei, 'comment'), (attr, key, valuei)
                        valuei.comment = ''
                else:
                    # dict[card]
                    assert hasattr(value, 'comment'), (attr, value)
                    value.comment = ''

        elif isinstance(dict_list_scalar, PathLike):
            pass
        else:
            if not hasattr(dict_list_scalar, 'comment'):
                #print(attr)
                continue
            assert hasattr(dict_list_scalar, 'comment'), (attr, value, type(value))
            dict_list_scalar.comment = ''

    _write_bdf(model, bdf_filename_out, size=size, is_double=is_double)
    return model
