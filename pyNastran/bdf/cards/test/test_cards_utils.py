"""Tests for pyNastran.bdf.cards.utils (wipe_empty_fields, build_table_lines)."""
import unittest
from pyNastran.bdf.cards.utils import wipe_empty_fields, wipe_empty_fields_str, build_table_lines


class TestWipeEmptyFields(unittest.TestCase):
    """Tests wipe_empty_fields: converts '' to None, strips strings, removes trailing Nones."""

    def test_trailing_nones_removed(self):
        """Trailing Nones are stripped; interior Nones preserved."""
        card = ['GRID', '1', None, '1.0', None, None]
        result = wipe_empty_fields(card)
        self.assertEqual(result, ['GRID', '1', None, '1.0'])

    def test_empty_strings_become_none(self):
        """Empty/whitespace strings become None, then trailing Nones stripped."""
        card = ['GRID', '1', '', '1.0', '  ', '']
        result = wipe_empty_fields(card)
        self.assertEqual(result, ['GRID', '1', None, '1.0'])

    def test_all_none_except_first(self):
        """Card with only the card name populated."""
        card = ['SET1', None, None, None]
        result = wipe_empty_fields(card)
        self.assertEqual(result, ['SET1'])

    def test_no_trailing_blanks(self):
        """Card with no trailing blanks unchanged (except string stripping)."""
        card = ['MAT1', '1', '3.0+7', None, '.3']
        result = wipe_empty_fields(card)
        self.assertEqual(result, ['MAT1', '1', '3.0+7', None, '.3'])

    def test_preserves_ints_and_floats(self):
        """Non-string fields (int, float) pass through unchanged."""
        card = ['GRID', 1, None, 1.0, 2.0, 3.0, None]
        result = wipe_empty_fields(card)
        self.assertEqual(result, ['GRID', 1, None, 1.0, 2.0, 3.0])

    def test_strings_stripped(self):
        """Strings with surrounding whitespace are stripped."""
        card = ['PSHELL', ' 1 ', '2', ' MAT ']
        result = wipe_empty_fields(card)
        self.assertEqual(result, ['PSHELL', '1', '2', 'MAT'])

    def test_single_field(self):
        """Single non-None element card."""
        self.assertEqual(wipe_empty_fields(['X']), ['X'])

    def test_interior_blanks_preserved(self):
        """Interior None/blank fields are kept as None."""
        card = ['FORCE', 1, None, None, 5.0]
        result = wipe_empty_fields(card)
        self.assertEqual(result, ['FORCE', 1, None, None, 5.0])

    def test_zero_int_at_end_preserved(self):
        """Trailing 0 (int) is not None, so it's preserved."""
        card = ['CARD', 1, 2, 0]
        result = wipe_empty_fields(card)
        self.assertEqual(result, ['CARD', 1, 2, 0])

    def test_zero_float_at_end_preserved(self):
        """Trailing 0.0 (float) is not None, so it's preserved."""
        card = ['GRID', 1, None, 0.0, 0.0, 0.0]
        result = wipe_empty_fields(card)
        self.assertEqual(result, ['GRID', 1, None, 0.0, 0.0, 0.0])

    def test_mixed_trailing_blanks(self):
        """Mixed trailing blank types (None, '', '  ') all removed."""
        card = ['MAT1', '1', '3.+7', None, '', '  ', None]
        result = wipe_empty_fields(card)
        self.assertEqual(result, ['MAT1', '1', '3.+7'])

    def test_all_fields_populated(self):
        """Card with all fields populated — nothing removed."""
        card = ['CQUAD4', '1', '2', '10', '20', '30', '40']
        result = wipe_empty_fields(card)
        self.assertEqual(result, ['CQUAD4', '1', '2', '10', '20', '30', '40'])

    def test_multiple_interior_blanks(self):
        """Multiple interior blank regions preserved."""
        # e.g., GRID: nid, cp(blank), x, y, z, cd(blank), ps(blank), seid
        card = ['GRID', '1', None, '0.0', '0.0', '0.0', None, None, '1']
        result = wipe_empty_fields(card)
        self.assertEqual(result, ['GRID', '1', None, '0.0', '0.0', '0.0', None, None, '1'])

    def test_only_blanks_after_first(self):
        """All fields blank except card name => single element."""
        card = ['DUMMY', '', '  ', None, None]
        result = wipe_empty_fields(card)
        self.assertEqual(result, ['DUMMY'])

    def test_last_field_is_string(self):
        """Last field is a non-empty string with spaces — stripped but kept."""
        card = ['CARD', None, ' VALUE ']
        result = wipe_empty_fields(card)
        self.assertEqual(result, ['CARD', None, 'VALUE'])

    def test_single_blank_string_returns_one_none(self):
        """A single blank string must return [None], not [].

        BDFCard relies on this to keep card.card non-empty so that
        print_card(card) doesn't crash with an IndexError.
        Do not change this behavior.
        """
        self.assertEqual(wipe_empty_fields(['']), [None])
        self.assertEqual(wipe_empty_fields(['  ']), [None])

    def test_single_none_returns_one_none(self):
        """A single None field must return [None], not [].

        Same rationale as test_single_blank_string_returns_one_none.
        Do not change this behavior.
        """
        self.assertEqual(wipe_empty_fields([None]), [None])


class TestWipeEmptyFieldsStr(unittest.TestCase):
    """Tests wipe_empty_fields_str: fast path for all-string inputs."""

    def test_single_blank_string_returns_one_none(self):
        """A single blank string must return [None], not [].

        Do not change this behavior.
        """
        self.assertEqual(wipe_empty_fields_str(['']), [None])
        self.assertEqual(wipe_empty_fields_str(['  ']), [None])

    def test_trailing_blanks_removed(self):
        result = wipe_empty_fields_str(['GRID', '1', '', '1.0', '', ''])
        self.assertEqual(result, ['GRID', '1', None, '1.0'])

    def test_all_blanks_after_first(self):
        result = wipe_empty_fields_str(['SET1', '', '  ', ''])
        self.assertEqual(result, ['SET1'])


class TestBuildTableLines(unittest.TestCase):
    """Tests build_table_lines: pads fields to 8-wide lines with nstart/nend offsets.

    build_table_lines inserts (nstart + nend) None pads every n = (8 - nstart - nend) fields,
    then calls wipe_empty_fields (strips trailing Nones), then pads to the next 8-boundary.
    """

    def test_exact_output_7_fields_nstart1(self):
        """7 fields with nstart=1, nend=0: n=7 fields per line.
        All 7 fit in one logical line, then pad to 8-boundary.
        Tolerances: exact output match.
        """
        # n = 8 - 1 - 0 = 7 fields per logical line
        # fields[0] appended, i=0 doesn't trigger pad
        # fields[1..6] appended, at i=7 triggers pad -> [None]*1 appended
        # but i only goes to 6 (7 fields, indices 0-6), so no pad triggered
        # wipe_empty_fields -> all non-None, nothing stripped
        # len=7, nspaces = 8 - 7%8 = 1, pad 1 None
        fields = [1, 2, 3, 4, 5, 6, 7]
        result = build_table_lines(fields, nstart=1, nend=0)
        self.assertEqual(result, [1, 2, 3, 4, 5, 6, 7, None])

    def test_exact_output_8_fields_nstart1(self):
        """8 fields with nstart=1: triggers one pad insertion at i=7.
        Tolerances: exact output match.
        """
        # n = 7, field indices 0-7
        # i=7: 7 % 7 == 0 and i > 0 -> insert [None]*1
        # fields_out = [1,2,3,4,5,6,7,8, None]
        # wipe: last non-None is index 7 (value 8) -> [1,2,3,4,5,6,7,8]
        # wait, None is appended AFTER 8, so fields_out = [1,2,3,4,5,6,7,None,8]
        # No -- let me re-read: append field, then check if i%n==0
        # i=0: append 1, 0%7=0 but i=0 so skip
        # i=1: append 2
        # ...
        # i=6: append 7
        # i=7: append 8, 7%7==0 and i>0 -> append [None]*1
        # fields_out = [1,2,3,4,5,6,7,8,None]
        # wipe_empty_fields: last non-None at index 7 -> [1,2,3,4,5,6,7,8]
        # nspaces = 8 - 8%8 = 8 -> not < 8, no pad
        fields = [1, 2, 3, 4, 5, 6, 7, 8]
        result = build_table_lines(fields, nstart=1, nend=0)
        self.assertEqual(len(result) % 8, 0)
        # 8 fields fills one full 8-wide line after wipe
        self.assertEqual(result, [1, 2, 3, 4, 5, 6, 7, 8])

    def test_exact_output_14_fields_nstart1(self):
        """14 fields with nstart=1: spans 2 logical lines.
        Interaction: pad inserted at line break, creating full 8-wide lines.
        Tolerances: exact output match.
        """
        # n = 7
        # i=0: append 1 (skip pad since i=0)
        # i=1..6: append 2-7
        # i=7: append 8, 7%7==0 -> append [None]
        # i=8..13: append 9-14
        # i=14: append... wait, only 14 fields so indices 0-13
        # i=13: append 14
        # fields_out = [1,2,3,4,5,6,7,8,None,9,10,11,12,13,14]
        # wipe: last non-None at index 14 -> same
        # nspaces = 8 - 15%8 = 8 - 7 = 1
        # pad 1 None -> len 16
        fields = list(range(1, 15))
        result = build_table_lines(fields, nstart=1, nend=0)
        expected = [1, 2, 3, 4, 5, 6, 7, 8, None, 9, 10, 11, 12, 13, 14, None]
        self.assertEqual(result, expected)
        self.assertEqual(len(result) % 8, 0)

    def test_exact_output_nstart2(self):
        """nstart=2: n=6 fields per logical line.
        Tolerances: exact output match.
        """
        # n = 8 - 2 - 0 = 6
        # fields = [1,2,3,4,5,6,7]
        # i=0: append 1 (skip)
        # i=1..5: append 2-6
        # i=6: append 7, 6%6==0 -> append [None]*2
        # fields_out = [1,2,3,4,5,6,7,None,None]
        # wipe: last non-None at index 6 -> [1,2,3,4,5,6,7]
        # nspaces = 8 - 7%8 = 1
        # pad -> [1,2,3,4,5,6,7,None]
        fields = [1, 2, 3, 4, 5, 6, 7]
        result = build_table_lines(fields, nstart=2, nend=0)
        expected = [1, 2, 3, 4, 5, 6, 7, None]
        self.assertEqual(result, expected)

    def test_exact_output_nend1(self):
        """nstart=1, nend=1: n=6 fields per logical line, pad=(nstart+nend)=2.
        Tolerances: exact output match.
        """
        # n = 8 - 1 - 1 = 6
        # fields = [1,2,3,4,5,6,7]
        # i=0: append 1 (skip)
        # i=1..5: append 2-6
        # i=6: append 7, 6%6==0 -> append [None]*2
        # fields_out = [1,2,3,4,5,6,7,None,None]
        # wipe: last non-None at index 6 -> [1,2,3,4,5,6,7]
        # nspaces = 8 - 7%8 = 1
        # pad -> [1,2,3,4,5,6,7,None]
        fields = [1, 2, 3, 4, 5, 6, 7]
        result = build_table_lines(fields, nstart=1, nend=1)
        expected = [1, 2, 3, 4, 5, 6, 7, None]
        self.assertEqual(result, expected)

    def test_exact_output_nstart1_nend1_multiline(self):
        """13 fields with nstart=1, nend=1: n=6, pads of 2 at each line break.
        Tolerances: exact output match.
        """
        # n = 6, pad_size = 2
        # i=0: append 1 (skip)
        # i=1-5: append 2-6
        # i=6: append 7, 6%6==0 -> append [None,None]
        # i=7-11: append 8-12
        # i=12: append 13, 12%6==0 -> append [None,None]
        # fields_out = [1,2,3,4,5,6,7,None,None,8,9,10,11,12,13,None,None]
        # wipe: last non-None at index 14 (value 13)
        # -> [1,2,3,4,5,6,7,None,None,8,9,10,11,12,13]
        # nspaces = 8 - 15%8 = 8 - 7 = 1
        # pad -> len 16
        fields = list(range(1, 14))
        result = build_table_lines(fields, nstart=1, nend=1)
        expected = [1, 2, 3, 4, 5, 6, 7, None, None, 8, 9, 10, 11, 12, 13, None]
        self.assertEqual(result, expected)
        self.assertEqual(len(result) % 8, 0)

    def test_single_field(self):
        """Single field pads to 8-boundary.
        Tolerances: exact output match.
        """
        fields = ['UM']
        result = build_table_lines(fields, nstart=1, nend=0)
        # fields_out = ['UM']
        # wipe -> ['UM']
        # nspaces = 8 - 1%8 = 7 -> pad 7 Nones
        expected = ['UM'] + [None] * 7
        self.assertEqual(result, expected)
        self.assertEqual(len(result), 8)

    def test_empty_fields(self):
        """Empty field list produces empty result."""
        result = build_table_lines([], nstart=1, nend=0)
        self.assertEqual(result, [])

    def test_exactly_n_fields(self):
        """Exactly n fields: no pad triggered (i%n==0 only at i=n which is field n+1).
        Tolerances: exact output match.
        """
        # n = 7 (nstart=1, nend=0), provide exactly 7 fields
        # i=0..6: append fields 0-6
        # i=6: 6%7 != 0, no pad
        # fields_out = [0,1,2,3,4,5,6]
        # wipe -> same (all non-None)
        # nspaces = 8 - 7%8 = 1 -> pad 1 None
        fields = list(range(7))
        result = build_table_lines(fields, nstart=1, nend=0)
        expected = [0, 1, 2, 3, 4, 5, 6, None]
        self.assertEqual(result, expected)

    def test_fields_already_8_boundary(self):
        """When fields naturally align to 8-boundary after padding, no extra pad added.
        Tolerances: exact output match, nspaces logic.
        """
        # 8 fields, nstart=1 -> n=7
        # After append+pad: [1..8, None]
        # wipe strips trailing None -> [1..8]
        # len=8, nspaces = 8 - 8%8 = 8, not < 8 so no pad
        fields = list(range(1, 9))
        result = build_table_lines(fields, nstart=1, nend=0)
        self.assertEqual(len(result) % 8, 0)
        self.assertEqual(len(result), 8)

    def test_all_fields_preserved_in_order(self):
        """Field ordering is preserved through padding.
        Tolerances: exact positional check.
        """
        fields = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I']
        result = build_table_lines(fields, nstart=1, nend=0)
        # Extract non-None fields and verify order
        non_none = [f for f in result if f is not None]
        self.assertEqual(non_none, fields)

    def test_none_pads_at_correct_positions(self):
        """None pads appear at correct positions for line breaks.
        Interaction: pad separates logical lines.
        Tolerances: exact positions of None insertions.
        """
        # n=7 (nstart=1, nend=0), 9 fields
        # After i=7: pad [None] inserted between field[7] and field[8]
        fields = list(range(1, 10))
        result = build_table_lines(fields, nstart=1, nend=0)
        # fields_out after loop: [1,2,3,4,5,6,7,8,None,9]
        # wipe: last non-None at 9 -> same
        # nspaces = 8 - 10%8 = 8-2 = 6
        # pad 6 Nones
        expected = [1, 2, 3, 4, 5, 6, 7, 8, None, 9] + [None] * 6
        self.assertEqual(result, expected)
        self.assertEqual(len(result), 16)

    def test_string_fields(self):
        """String fields work (used in DRESP2 for DESVAR/DTABLE labels).
        Tolerances: exact output match.
        """
        fields = ['DESVAR', 101, 102, 103, 104, 105, 106, 107, 108]
        result = build_table_lines(fields, nstart=1, nend=0)
        # n=7, i=7 triggers pad
        # fields_out = ['DESVAR',101,102,103,104,105,106,107,None,108]
        # wipe: last non-None at 9 -> same (but None at index 8 stays since it's interior)
        # wait: wipe removes trailing Nones only. None at index 8 is interior.
        # nspaces = 8 - 10%8 = 6
        expected = ['DESVAR', 101, 102, 103, 104, 105, 106, 107, None, 108] + [None] * 6
        self.assertEqual(result, expected)
        self.assertEqual(len(result), 16)


if __name__ == '__main__':  # pragma: no cover
    unittest.main()
