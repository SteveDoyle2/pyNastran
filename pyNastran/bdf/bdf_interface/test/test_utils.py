"""Tests for pyNastran.bdf.bdf_interface.utils (to_fields, expand_tabs)."""
import unittest
from io import StringIO
from pyNastran.bdf.bdf_interface.utils import to_fields, expand_tabs
from pyNastran.bdf.bdf import read_bdf
from pyNastran.bdf.errors import CardParseSyntaxError


class TestExpandTabs(unittest.TestCase):
    """Tests for expand_tabs: replaces tabs with spaces, rejects tab+comma mix."""

    def test_no_tabs(self):
        """Line without tabs passes through unchanged."""
        line = 'GRID           1       0     0.0     0.0     0.0'
        self.assertEqual(expand_tabs(line), line)

    def test_tabs_expanded(self):
        """Tabs are expanded to spaces."""
        line = 'GRID\t1\t\t1.0\t2.0\t3.0'
        result = expand_tabs(line)
        self.assertNotIn('\t', result)
        self.assertIn('GRID', result)
        self.assertIn('1.0', result)

    def test_tabs_and_commas_raises(self):
        """Mixing tabs and commas raises CardParseSyntaxError."""
        line = 'GRID,1,\t,1.0'
        with self.assertRaises(CardParseSyntaxError):
            expand_tabs(line)

    def test_single_tab(self):
        """Single tab in a line is expanded."""
        line = 'A\tB'
        result = expand_tabs(line)
        self.assertNotIn('\t', result)
        self.assertIn('A', result)
        self.assertIn('B', result)


class TestToFields(unittest.TestCase):
    """Tests for to_fields: converts card lines into field strings.
    Interactions tested: small/large field, CSV/fixed-width, continuation lines, tabs.
    """

    def test_small_field_fixed_width(self):
        """Standard 8-column small field format.
        Tolerances: exact field count (9 for first line).
        """
        # 80 chars: 8 + 8*8 = 72 (fields occupy cols 1-72)
        line = 'GRID           1       0     0.0     0.0     0.0       0               '
        fields = to_fields([line], 'GRID')
        self.assertEqual(fields[0].strip(), 'GRID')
        self.assertEqual(len(fields), 9)
        self.assertEqual(fields[1].strip(), '1')
        self.assertEqual(fields[2].strip(), '0')

    def test_small_field_csv(self):
        """CSV small field format (comma-separated).
        Tolerances: exact field values.
        """
        line = 'GRID,1,,1.0,2.0,3.0'
        fields = to_fields([line], 'GRID')
        self.assertEqual(fields[0], 'GRID')
        self.assertEqual(fields[1], '1')
        self.assertEqual(fields[2], '')
        self.assertEqual(fields[3], '1.0')
        self.assertEqual(fields[4], '2.0')
        self.assertEqual(fields[5], '3.0')
        # pads to 9 fields
        self.assertEqual(len(fields), 9)

    def test_large_field_fixed_width(self):
        """Large field format (16-char fields, * in card name).
        Tolerances: 5 fields for first line.
        """
        # Large field: 8 + 16 + 16 + 16 + 16 = 72
        line = 'GRID*                  1               0             0.0             0.0'
        fields = to_fields([line], 'GRID')
        self.assertEqual(len(fields), 5)
        self.assertEqual(fields[0].strip(), 'GRID*')
        self.assertEqual(fields[1].strip(), '1')

    def test_large_field_csv(self):
        """Large field CSV format.
        Tolerances: 5 fields for first line.
        """
        line = 'GRID*,1,,1.0,2.0'
        fields = to_fields([line], 'GRID')
        self.assertEqual(len(fields), 5)
        self.assertEqual(fields[0], 'GRID*')
        self.assertEqual(fields[1], '1')

    def test_continuation_small_field(self):
        """Small field with continuation line.
        Tolerances: field count = 9 + 8 = 17 for two lines.
        """
        line1 = 'PBAR           1       1    1.0    2.0    3.0    4.0                    '
        line2 = '             5.0    6.0    7.0    8.0    9.0   10.0   11.0   12.0'
        fields = to_fields([line1, line2], 'PBAR')
        self.assertEqual(len(fields), 9 + 8)
        self.assertEqual(fields[0].strip(), 'PBAR')
        self.assertEqual(fields[1].strip(), '1')

    def test_continuation_csv(self):
        """CSV with continuation line.
        Tolerances: correct field count.
        """
        line1 = 'GRID,1,,1.0,2.0,3.0,,,,'
        line2 = ',100'
        fields = to_fields([line1, line2], 'GRID')
        # first line: 9 fields, continuation: 8 fields
        self.assertEqual(len(fields), 9 + 8)
        self.assertEqual(fields[0], 'GRID')
        # field index 9 is the first continuation field
        self.assertEqual(fields[9], '100')

    def test_tabs_in_line(self):
        """Tabs in card lines are expanded (allow_tabs=True default).
        Tolerances: fields extracted correctly.
        """
        line = 'GRID\t       1       0     0.0     0.0     0.0'
        fields = to_fields([line], 'GRID')
        self.assertEqual(fields[0].strip(), 'GRID')

    def test_tabs_disallowed(self):
        """allow_tabs=False raises RuntimeError when tabs present."""
        line = 'GRID\t1'
        with self.assertRaises(RuntimeError):
            to_fields([line], 'GRID', allow_tabs=False)

    def test_equal_sign_raises(self):
        """Equal sign in card line raises CardParseSyntaxError (replication not supported here)."""
        line = 'GRID    =1'
        with self.assertRaises(CardParseSyntaxError):
            to_fields([line], 'GRID')

    def test_blank_fields_preserved(self):
        """Blank fields in fixed-width format are empty strings.
        Tolerances: blank fields are whitespace-only strings.
        """
        line = 'GRID           1                     0.0     0.0     0.0                '
        fields = to_fields([line], 'GRID')
        # field 2 (cp) should be blank (all spaces)
        self.assertEqual(fields[2].strip(), '')

    def test_csv_trailing_commas(self):
        """Trailing commas produce empty string fields.
        Tolerances: padded to 9 fields.
        """
        line = 'FORCE,1,2,0,1.0,,,'
        fields = to_fields([line], 'FORCE')
        self.assertEqual(len(fields), 9)
        self.assertEqual(fields[0], 'FORCE')
        self.assertEqual(fields[5], '')
        self.assertEqual(fields[6], '')

    def test_continuation_large_field(self):
        """Large field continuation line (16-char fields).
        Tolerances: 4 fields per continuation line.
        """
        line1 = 'GRID*                  1                             1.0             2.0'
        line2 = '*                    3.0                                                '
        fields = to_fields([line1, line2], 'GRID')
        # first line: 5 fields, continuation: 4 fields
        self.assertEqual(len(fields), 5 + 4)
        self.assertEqual(fields[5].strip(), '3.0')

    def test_continuation_marker_stripped_first_line(self):
        """A '+' continuation marker in field 9 (cols 65-72) of the first line
        is stripped so it doesn't appear as data.
        Tolerances: field 8 (0-indexed) is blank, continuation data appears in fields 9+.
        """
        # TRIM card with '+' continuation marker in cols 65-72
        line1 = 'TRIM    1       0.789   1.5     PITCH   0.0     URDD3   2.5     +'
        line2 = '+       URDD5   0.0'
        fields = to_fields([line1, line2], 'TRIM')
        self.assertEqual(fields[0].strip(), 'TRIM')
        self.assertEqual(fields[1].strip(), '1')
        self.assertEqual(fields[7].strip(), '2.5')
        # field 8 (cols 65-72) had '+' which should be stripped
        self.assertEqual(fields[8].strip(), '')
        # continuation line data
        self.assertEqual(fields[9].strip(), 'URDD5')
        self.assertEqual(fields[10].strip(), '0.0')

    def test_continuation_marker_stripped_continuation_line(self):
        """A '+' continuation marker in field 9 (cols 65-72) of a continuation
        line is stripped so it doesn't appear as data.
        Tolerances: the marker field is blank, subsequent continuation data is correct.
        """
        line1 = 'CQUAD4  1       1       1       2       3       4               +'
        line2 = '+       5.0     6.0     7.0     8.0     9.0     10.0    11.0    +'
        line3 = '+       12.0'
        fields = to_fields([line1, line2, line3], 'CQUAD4')
        # line1 field 8 (cols 65-72) = '+' -> stripped
        self.assertEqual(fields[8].strip(), '')
        # line2: 8 data fields in cols 9-72, last field (cols 65-72) = '+' -> stripped
        self.assertEqual(fields[16].strip(), '')
        # line3: first continuation field
        self.assertEqual(fields[17].strip(), '12.0')

    def test_named_continuation_marker_stripped(self):
        """A named continuation marker like '+ABC' in cols 65-72 is stripped.
        Tolerances: the marker field is blank.
        """
        line1 = 'GRID           1       0     0.0     0.0     0.0       0        +G1'
        line2 = '+G1      100'
        fields = to_fields([line1, line2], 'GRID')
        # field 8 had '+G1' -> stripped
        self.assertEqual(fields[8].strip(), '')
        # continuation data
        self.assertEqual(fields[9].strip(), '100')

    def test_positive_number_in_last_field_not_stripped(self):
        """A positive number like '+1.5' in cols 65-72 must NOT be stripped.
        Tolerances: the value is preserved exactly.
        """
        # 72 chars with '+1.5' in field 9 (cols 65-72)
        line = 'CARD    1       2       3.0     4.0     5.0     6.0     7.0     +1.5    '
        self.assertEqual(len(line), 72)
        fields = to_fields([line], 'CARD')
        self.assertEqual(fields[8].strip(), '+1.5')

    def test_positive_exponent_in_last_field_not_stripped(self):
        """A positive number with exponent like '+3E5' in cols 65-72 must NOT be stripped.
        Tolerances: the value is preserved exactly.
        """
        line = 'CARD    1       2       3.0     4.0     5.0     6.0     7.0     +3E5    '
        self.assertEqual(len(line), 72)
        fields = to_fields([line], 'CARD')
        self.assertEqual(fields[8].strip(), '+3E5')

    def test_positive_decimal_in_last_field_not_stripped(self):
        """A value like '+.5' in cols 65-72 must NOT be stripped.
        Tolerances: the value is preserved exactly.
        """
        line = 'CARD    1       2       3.0     4.0     5.0     6.0     7.0     +.5     '
        self.assertEqual(len(line), 72)
        fields = to_fields([line], 'CARD')
        self.assertEqual(fields[8].strip(), '+.5')

    def test_star_continuation_marker_stripped(self):
        """A '*' continuation marker in cols 65-72 is stripped (large-to-small transition).
        Tolerances: the marker field is blank.
        """
        line1 = 'CARD           1       2     3.0     4.0     5.0     6.0       *'
        line2 = '*        100'
        fields = to_fields([line1, line2], 'CARD')
        self.assertEqual(fields[8].strip(), '')

    def test_csv_line0_overflow_raises(self):
        """A free-field line-0 with more than 10 fields raises CardParseSyntaxError."""
        line = 'GRID,1,,1.0,2.0,3.0,,,1,999,888'
        with self.assertRaises(CardParseSyntaxError):
            to_fields([line], 'GRID')

    def test_csv_field10_data_raises(self):
        """Data in field 10 (the continuation slot) raises CardParseSyntaxError."""
        line = 'GRID,1,,1.0,2.0,3.0,,,1,999'
        with self.assertRaises(CardParseSyntaxError):
            to_fields([line], 'GRID')

    def test_csv_field10_marker_ok(self):
        """A continuation marker in field 10 is dropped, not an error."""
        line = 'GRID,1,,1.0,2.0,3.0,,,1,+'
        fields = to_fields([line], 'GRID')
        self.assertEqual(len(fields), 9)
        self.assertNotIn('+', fields)

    def test_csv_field10_numeric_marker_ok(self):
        """A numeric continuation marker like '+01' is allowed in field 10."""
        line = 'TABLED1,100,,,,,,,,+01'
        fields = to_fields([line], 'TABLED1')
        self.assertEqual(fields[0], 'TABLED1')
        self.assertEqual(fields[1], '100')

    def test_csv_trailing_comma_marker_ok(self):
        """A trailing comma after the field-10 marker is not counted as a field."""
        line = 'CORD2R,5,,0.,0.,0.,0.,0.,1.0,+C1,'
        fields = to_fields([line], 'CORD2R')
        self.assertEqual(len(fields), 9)
        self.assertEqual(fields[8], '1.0')

    def test_csv_continuation_overflow_raises(self):
        """A continuation line with 10 data tokens raises CardParseSyntaxError."""
        line1 = 'PBAR,1,1,1.0,2.0,3.0,4.0,5.0,6.0'
        line2 = ',8.0,9.0,10.0,11.0,12.0,13.0,14.0,15.0,16.0,17.0'
        with self.assertRaises(CardParseSyntaxError):
            to_fields([line1, line2], 'PBAR')

    def test_csv_continuation_9_data_raises(self):
        """9 data fields on a continuation line is overflow (field 10 = data)."""
        line1 = 'GRID,1,,1.0,2.0,3.0,,,1'
        line2 = ',100,101,102,103,104,105,106,107,108'
        with self.assertRaises(CardParseSyntaxError):
            to_fields([line1, line2], 'GRID')

    def test_csv_line0_marker_stripped(self):
        """A marker in the last data slot of a CSV line-0 is stripped."""
        line = 'GRID,1,,1.0,2.0,3.0,,,+'
        fields = to_fields([line], 'GRID')
        self.assertEqual(len(fields), 9)
        self.assertNotIn('+', fields)

    def test_csv_continuation_marker_stripped(self):
        """A marker in the last data slot of a CSV continuation is stripped."""
        line1 = 'GRID,1,,1.0,2.0,3.0,,,1'
        line2 = ',100,101,102,103,104,105,106,+B'
        fields = to_fields([line1, line2], 'GRID')
        self.assertEqual(len(fields), 9 + 8)
        self.assertNotIn('+B', fields)


class TestReadBdfFreeFieldOverflow(unittest.TestCase):
    """read_bdf rejects free-field lines that exceed the 10-field line grammar."""

    def test_read_bdf_grid_line0_overflow_raises(self):
        """An 11-token GRID line raises instead of silently dropping the overflow."""
        bdf = 'SOL 101\nCEND\nBEGIN BULK\nGRID,1,,1.0,2.0,3.0,,,1,999,888\nENDDATA\n'
        with self.assertRaises(CardParseSyntaxError):
            read_bdf(StringIO(bdf), validate=False, xref=False)

    def test_read_bdf_pbar_continuation_overflow_raises(self):
        """A PBAR continuation with 10 data tokens raises instead of using defaults."""
        bdf = (
            'SOL 101\nCEND\nBEGIN BULK\n'
            'PBAR,1,1,1.0,2.0,3.0,4.0,5.0,6.0,7.0\n'
            ',8.0,9.0,10.0,11.0,12.0,13.0,14.0,15.0,16.0,17.0\n'
            'ENDDATA\n')
        with self.assertRaises(CardParseSyntaxError):
            read_bdf(StringIO(bdf), validate=False, xref=False)

    def test_read_bdf_cord2c_one_line_rejected(self):
        """A one-line 12-token CORD2C is malformed and raises."""
        bdf = (
            'SOL 101\nCEND\nBEGIN BULK\n'
            'CORD2C,1,0,0.,0.,0.,0.,0.,1.,0.,1.,0.\n'
            'ENDDATA\n')
        with self.assertRaises(CardParseSyntaxError):
            read_bdf(StringIO(bdf), validate=False, xref=False)

    def test_read_bdf_cord2c_two_line_ok(self):
        """The same CORD2C in two lines parses and gives the written frame."""
        bdf = (
            'SOL 101\nCEND\nBEGIN BULK\n'
            'CORD2C,1,0,0.,0.,0.,0.,0.,1.\n'
            ',0.,1.,0.\n'
            'ENDDATA\n')
        model = read_bdf(StringIO(bdf), validate=False, xref=False)
        coord = model.coords[1]
        self.assertEqual(coord.e2.tolist(), [0., 0., 1.])
        self.assertEqual(coord.e3.tolist(), [0., 1., 0.])

    def test_read_bdf_dmi_row_overflow_ok(self):
        """DMI matrix data rows may exceed 10 fields (variable-length rows)."""
        bdf = (
            'DMI,ZRO,0,2,1,0,,11,1\n'
            'DMI,ZRO,1,1,1.0,1.0,1.0,1.0,1.0,2.0,3.0\n')
        model = read_bdf(StringIO(bdf), validate=False, xref=False, punch=True)
        self.assertIn('ZRO', model.dmi)


if __name__ == '__main__':  # pragma: no cover
    unittest.main()
