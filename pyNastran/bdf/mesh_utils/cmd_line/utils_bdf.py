from __future__ import annotations
from typing import TYPE_CHECKING
import os
import numpy as np

# from cpylog import SimpleLogger
if TYPE_CHECKING:
    from pyNastran.bdf.bdf import BDF


def read_lax_bdf(bdf_filename: str,
                 punch: bool=False,
                 validate: bool=True,
                 xref: bool=True,
                 is_strict_card_parser: bool=False,
                 cards_to_read: list[str]=None,
                 duplicate_cards: list[str]=None,
                 log=None) -> BDF:
    from pyNastran.bdf.bdf import BDF
    model = BDF(log=log)
    if not is_strict_card_parser:
        log.warning('using lax card parser')
        model.is_strict_card_parser = is_strict_card_parser
    if cards_to_read is not None:
        model.enable_cards(cards_to_read)
    if duplicate_cards is not None:
        #duplicate_cards = {'GRID', 'CONM2'}
        model.set_allow_duplicates(set(duplicate_cards))
    model.read_bdf(bdf_filename, punch=punch,
                   validate=validate, xref=xref)
    return model


def read_lax_obj(bdf_filename, obj_filename, is_obj: bool,
                 xref: bool=True,
                 is_strict_card_parser: bool=True,
                 punch: bool=False,
                 duplicate_cards=None,
                 log=None) -> BDF:
    if is_obj:
        from pyNastran.bdf.bdf import BDF
        if os.path.exists(obj_filename):
            model = BDF(log=log)
            model.load(obj_filename)
            model.cross_reference()
        else:
            model = read_lax_bdf(
                bdf_filename, punch=punch, xref=False,
                validate=False,
                is_strict_card_parser=is_strict_card_parser,
                duplicate_cards=duplicate_cards,
                log=log)
            model.save(obj_filename)
    else:
        model = read_lax_bdf(
            bdf_filename, punch=punch, xref=False,
            validate=False,
            is_strict_card_parser=is_strict_card_parser,
            duplicate_cards=duplicate_cards,
            log=log)
    return model

def get_ids(ids_str, default=None, dtype: str='int32'):
    """
    Supports FEMAP and CSV format
    FEMAP format:
        7100001,7101743,1;1,10,1
    CSV format:
        1,2,3,4,5,10
    Patran Format:
        7100001:7101743:1;1:10:1

    """
    ids_str = ids_str.strip(',;')
    if len(ids_str) == 0:
        return default

    ids_list = []
    if ';' in ids_str:
        ids_split = ids_str.split(';')
        for idi in ids_split:
            if idi.isdigit():
                ids_list.append(int(idi))
            elif ',' in idi:
                spliti = idi.split(',')
                if len(spliti) == 3:
                    start, stop, space = (int(val) for val in spliti)
                    ids_list.extend(list(range(start, stop + 1, space)))
    elif ':' in ids_str:
        raise NotImplementedError(': is not supported')
    elif ',' in ids_str:
        ids_list.extend(list(int(val) for val in ids_str.split(',')))
    else:
        raise NotImplementedError(f'ids_str={ids_str}')
    ids = np.array(ids_list, dtype=dtype)
    return ids
