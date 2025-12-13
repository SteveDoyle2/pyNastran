from __future__ import annotations
from typing import TYPE_CHECKING

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
