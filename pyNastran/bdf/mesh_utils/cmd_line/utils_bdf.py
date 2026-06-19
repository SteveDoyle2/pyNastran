from __future__ import annotations
from typing import TYPE_CHECKING
import os
import numpy as np
from pyNastran.bdf.bdf import BDF

if TYPE_CHECKING:  # pragma: no cover
    from cpylog import SimpleLogger
    from pyNastran.utils import PathLike


def read_lax_bdf(bdf_filename: PathLike,
                 punch: bool=False,
                 validate: bool=True,
                 xref: bool=True,
                 read_includes: bool=True,
                 is_strict_card_parser: bool=False,
                 save_file_structure: bool=False,
                 cards_to_read: list[str]=None,
                 duplicate_cards: list[str]=None,
                 log: SimpleLogger | None=None) -> BDF:
    model = BDF(log=log)
    log = model.log
    model.read_includes = read_includes
    if not is_strict_card_parser:
        log.warning('using lax card parser')
        model.is_strict_card_parser = is_strict_card_parser
    if cards_to_read is not None:
        model.enable_cards(cards_to_read)
    if duplicate_cards is not None:
        log.warning(f'allowing duplicate cards={duplicate_cards}')
        #duplicate_cards = {'GRID', 'CONM2'}
        model.set_allow_duplicates(set(duplicate_cards))
    model.read_bdf(bdf_filename, punch=punch,
                   validate=validate, xref=xref, save_file_structure=save_file_structure)
    return model

def read_bdf_obj(bdf_filename: PathLike,
                 punch: bool = False,
                 xref: list[str] | bool=True,
                 load_obj: bool=True,
                 save_obj: bool=True,
                 is_strict_card_parser: bool=True,
                 duplicate_cards: list[str] | set[str] | None=None,
                 log: SimpleLogger | None =None) -> BDF:
    if isinstance(bdf_filename, PathLike):
        base = os.path.splitext(bdf_filename)[0]
        obj_filename = base + '.obj'
        model = read_lax_obj(
            bdf_filename, obj_filename, load_obj=load_obj, save_obj=save_obj,
            punch=punch, xref=xref,
            is_strict_card_parser=is_strict_card_parser,
            duplicate_cards=duplicate_cards,
            log=log,)
        return model
    assert isinstance(bdf_filename, BDF), bdf_filename
    return bdf_filename

def read_lax_obj(bdf_filename: PathLike,
                 obj_filename: PathLike,
                 xref: list[str] | bool=True,
                 read_includes: bool=True,
                 is_strict_card_parser: bool=True,
                 save_file_structure: bool=False,
                 punch: bool=False,
                 duplicate_cards: list[str] | set[str] | None=None,
                 load_obj: bool=True,
                 save_obj: bool=True,
                 log=None) -> BDF:
    """
    Parameters
    ----------
    bdf_filename : PathLike
        bdf filename
    obj_filename : PathLike
        object filename
    load_obj : bool; default=True
        consider the obj_filename to load from
    xref : bool; default=True
         cross-references the part
    is_strict_card_parser : bool; default=True
        ???
    punch : bool; default=False
        define a punch file
    duplicate_cards : set[str]
        set to {GRID, CONM2} ot allow replacing of duplicate cards
    log : SimpleLogger
        a logging object

    Returns
    -------
    model : BDF()
        the model object

    """
    if load_obj:
        if os.path.exists(obj_filename):
            model = BDF(log=log)
            model.load(obj_filename)
            # model.cross_reference(xref=xref)
        else:
            model = read_lax_bdf(
                bdf_filename, punch=punch, xref=False,
                validate=False,
                read_includes=read_includes,
                is_strict_card_parser=is_strict_card_parser,
                save_file_structure=save_file_structure,
                duplicate_cards=duplicate_cards,
                log=log)
            if save_obj:
                model.save(obj_filename)
            # model.cross_reference(xref=xref)
    else:
        model = read_lax_bdf(
            bdf_filename, punch=punch, xref=False,
            validate=False,
            read_includes=read_includes,
            is_strict_card_parser=is_strict_card_parser,
            save_file_structure=save_file_structure,
            duplicate_cards=duplicate_cards,
            log=log)
        if save_obj:
            model.save(obj_filename)
    model.apply_xref_list(xref)
    return model

def get_ids(ids_str: str, default=None, dtype: str='int32'):
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
