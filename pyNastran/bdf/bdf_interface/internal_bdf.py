from __future__ import annotations
from typing import TYPE_CHECKING

from pyNastran.bdf.field_writer_16 import print_field_16
from pyNastran.bdf.errors import ReplicationError
from pyNastran.bdf.cards.utils import wipe_empty_fields

from .bdf_card import BDFCard
from .replication import (
    to_fields_replication, get_nrepeats, int_replication, float_replication,
    _field, repeat_cards)
from .utils import (
    to_fields, _parse_dynamic_syntax,)
if TYPE_CHECKING:
    from cpylog import SimpleLogger
    from pyNastran.bdf.bdf import BDF


def get_file_tag(model: BDF, ifile_iline: int) -> str:
    """get the filename

    ifile_iline = np.array([194, 1], dtype=int32)
    """
    # print(f"ifile_iline = {ifile_iline!r}")
    ifile, iline = ifile_iline
    filename = model.active_filenames[ifile-1]
    # print(f"model.save_file_structure = {model.save_file_structure}")
    # if not model.save_file_structure:
    #     return ''
    # tag = '\nfiles:\n'
    # for i, fname in enumerate(model.active_filenames):
    #     tag += f' - {i}: {fname}\n'
    # tag += f'\nline={iline+1} file={filename}\n'
    tag = f'\nline={iline+1} file={filename}\n'
    return tag


def expand_replication(card_name: str, icard: int,
                       cards_list: list,
                       card_lines_new: list,
                       dict_of_vars: dict[str, int | float | str],
                       log: SimpleLogger,
                       is_dynamic_syntax: bool=False,
                       dig: bool = True):
    """replication helper"""
    # dig_str = '  ' if dig is False else ''
    # print(dig_str, '-----------************---------')
    # print(dig_str, '--dig=%s--' % dig)
    # print(dig_str, 'card_lines_new=%s' % card_lines_new)
    card = []
    cards = []
    card_lines_old = cards_list[icard - 1][2]

    is_star_lines = any('*' in line for line in card_lines_old)
    if is_star_lines:
        # new_fields = to_fields_replication(card_lines_old)
        old_card = to_fields_replication(card_lines_old)
    else:
        # old_card, unused_card = self.create_card_object(
        #     card_lines_old, card_name,
        #     is_list=False, has_none=True)
        # print(card_lines_old)
        old_card = _old_card_fields(
            card_lines_old, card_name, dict_of_vars, log,
            is_list=False, has_none=True,
            is_dynamic_syntax=is_dynamic_syntax)
        # print(old_card)
        # assert '=' not in card_name

    nlines = len(card_lines_new)
    # print(dig_str, "card_lines_new =", card_lines_new)
    new_card = to_fields_replication(card_lines_new)
    assert len(card_lines_new) == nlines, card_lines_new

    # print(dig_str, 'old_card =', old_card)
    # print(dig_str, 'card_name = %r' % card_name)
    old_card_real = None
    if old_card[0] == '=':
        # print(dig_str, 'A!!!')
        if dig is False:
            raise ReplicationError(f'dig=False...old_card=\n{old_card}')

        cards2 = expand_replication(
            card_name, icard - 1, cards_list, card_lines_old,
            dict_of_vars, log,
            is_dynamic_syntax=is_dynamic_syntax,
            dig=False)
        assert len(cards2) == 1, f'cards2={cards2}; ncards={len(cards2)}'
        # print(dig_str, 'cards_equal =', cards2)
        old_card_fields = cards2[0]
        old_card_real = old_card
        # print(dig_str, 'old_card_fields =', old_card_fields)
        # print(dig_str, 'old_card_real =', old_card_real)
        # print(dig_str, 'card_lines_old =', card_lines_old)
        old_card = _old_card_fields(
            old_card_fields, card_name, dict_of_vars, log,
            is_list=True, has_none=True,
            is_dynamic_syntax=is_dynamic_syntax)
    elif '=' in card_name:
        # print(dig_str, 'B!!!')
        # print(dig_str, 'old_card =', old_card)
        # print(dig_str, 'card_lines_new =', card_lines_new)
        # print(dig_str, 'card_name = %r' % card_name)

        # good
        # new_card = [u'=3']
        # old_card = [u'CQUAD4', u'64', u'1', u'88', u'89', u'101', u'100']
        # old_card_real = [u'=', u'*1', u'=', u'*1', u'*1', u'*1', u'*1']

        # bad
        # card_name = u'=(7)'
        # new_card = [u'=(7)', u'*(10)', u'=', u'=', u'=', u'*(1.0)']
        # old_card = [u'grid', u'1001', None, u'0.', u'0.', u'0.']
        # old_card_real = [u'grid', u'1001', None, u'0.', u'0.', u'0.']
        old_card_real = new_card
    # else:
    #     print('old_card[0] %r' % old_card[0])

    # print(dig_str, "old_card =", old_card)
    # print(dig_str, "new_card =", new_card)
    for ifield, field in enumerate(new_card):
        if field is None:
            field2 = old_card.field(ifield)
            # print(' %i: %r -> %r' % (ifield, field, field2))
            # assert field2 is None, 'field=%s field2=%s' % (field, field2)
            card.append(field2)
            continue

        # if field == '':
        #     pass

        field = field.strip()
        if field == '=':
            field2 = old_card[ifield]
            # field2 = old_card.field(ifield)
        elif field == '==':
            # just append the remaining fields
            card.extend(old_card[ifield:])
            # print(dig_str, ' %i : extending %s' % (ifield, old_card[ifield:]))
            # print(dig_str, ' break _expand_replication...')
            break
        elif '=' in field:
            # =4
            assert ifield == 0, f'ifield={ifield} field={field!r} new_card={new_card}'
            nrepeats = get_nrepeats(field, old_card, new_card)
            if old_card_real is None:
                # old_card_real = old_card
                msg = (
                        'Invalid Replication Syntax (continuations arent supported)\n'
                        'old:\n%s\n'
                        'new:\n%s'
                        % (old_card, new_card))
                raise RuntimeError(msg)

            # new_card = [u'=3']
            # old_card = [u'CQUAD4', u'64', u'1', u'88', u'89', u'101', u'100']
            # old_card_real = [u'=', u'*1', u'=', u'*1', u'*1', u'*1', u'*1']

            # print('---')
            # print(dig_str, "nrepeats =", nrepeats)
            new_card[0] = '='
            # print(dig_str, "new_card =", new_card)
            # print(dig_str, "old_card =", old_card)
            # print(dig_str, "old_card_real =", old_card_real)
            for unused_irepeat in range(nrepeats):
                repeated_cards = repeat_cards(old_card, old_card_real)
                if len(repeated_cards) != 1:
                    for repeated_card in repeated_cards:
                        print("  repeated_card =", repeated_card)
                    raise RuntimeError('too many repeated cards')
                # for repeated_card in repeated_cards:
                #     print("  repeated_card =", repeated_card)
                repeated_card = repeated_cards[0]
                cards.append(repeated_card)
                old_card = repeated_card
            #     print(dig_str, "  repeated_card =", repeated_card)
            # print(dig_str, 'breaking...')
            return cards

        elif '*' in field:
            # this is an increment, not multiplication...
            old_field = _field(old_card, ifield)
            assert old_field is not None, f'old_card:{old_card}\nnew_card:\n{new_card}'
            try:
                if '.' in field:
                    field2 = float_replication(field, old_field)
                else:
                    field2 = int_replication(field, old_field)
            except Exception:
                log.error(f'old_card:{old_card}\nnew_card:\n{new_card}')
                raise
        else:
            assert '(' not in field, f'field={field!r}'
            assert '*' not in field, f'field={field!r}'
            assert '=' not in field, f'field={field!r}'
            field2 = field
        # print(dig_str, ' %i: %r -> %r' % (ifield, field, field2))
        card.append(field2)
    if card:
        cards.append(card)
        # print(dig_str, 'card_expanded = %s' % card)
    else:  # pragma: no cover
        raise RuntimeError(card)
    return cards


def check_replicated_cards(replicated_cards: list) -> None:
    """helper method for ``parse_cards_list``"""
    replicated_card_old = []
    try:
        for replicated_card in replicated_cards:
            assert replicated_card != replicated_card_old
            replicated_card_old = replicated_card
    except AssertionError:
        #print('card_list = %s' % card_list)
        #print('card_lines = %s' % card_lines)
        replicated_card_old = []
        for replicated_card in replicated_cards:
            #print('adding ', replicated_card)
            assert replicated_card != replicated_card_old
            replicated_card_old = replicated_card
        raise


def _old_card_fields(card_lines: list[str], card_name: str,
                     dict_of_vars: dict[str, int | float | str],
                     log: SimpleLogger,
                     is_list: bool=False, has_none: bool=True,
                     is_dynamic_syntax: bool=False) -> BDFCard:
    """replication helper"""
    if is_list:
        fields = card_lines
    else:
        fields = to_fields(card_lines, card_name)

    # apply OPENMDAO syntax
    if is_dynamic_syntax:
        fields = [print_field_16(_parse_dynamic_syntax(field, dict_of_vars, log))
                  if '%' in field.strip()[0:1] else print_field_16(field)
                  for field in fields]
        has_none = False

    if has_none:
        card = wipe_empty_fields([print_field_16(field) for field in fields])
    else:
        card = wipe_empty_fields(fields)
    card_obj = BDFCard(card, has_none=False)
    return card_obj
