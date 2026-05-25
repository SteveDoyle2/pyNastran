import unittest
from io import StringIO
from pyNastran.bdf.bdf import BDF
from pyNastran.bdf.bdf_interface.dev.dmap import (
    make_matrix_export_alter, add_matrix_export_alter,
    make_matrix_export_malter, add_matrix_export_malter)

class TestExecutive(unittest.TestCase):
    def test_blank_id(self):
        lines = (
            'NASTRAN SYSTEM(442)=-1,SYSTEM(319)=1\n'
            'ID NASTRAN,Femap\n'
            'SOL SESTATIC\n'
            'GEOMCHECK, NONE\n'
            'CEND\n'
            'SUBCASE 1\n'
            '  LOAD = 1\n'
            'BEGIN BULK\n'
            '$*\n'
            'GRID,1\n'
        )
        bdf_file = StringIO()
        bdf_file.write(lines)
        bdf_file.seek(0)
        model = BDF()
        model.read_bdf(bdf_filename=bdf_file, validate=True, xref=True,
                       punch=False, read_includes=True, save_file_structure=False, encoding=None)
        assert len(model.nodes) == 1, model.nodes

    def test_dmap_alter_op4(self):
        """Test generating DMAP ALTER for OUTPUT4 export.

        Verifies ASSIGN, COMPILE, ALTER, and OUTPUT4 lines are generated.
        Default ALTER location is 'END'.
        """
        system_lines, exec_lines = make_matrix_export_alter(
            matrices=['KGG', 'MGG'],
            output_filename='stiffness.op4',
            output_format='op4',
            sol=101,
            unit=51,
            form='FORMATTED',
        )
        assert len(system_lines) == 1
        assert "ASSIGN OUTPUT4='stiffness.op4'" in system_lines[0]
        assert 'UNIT=51' in system_lines[0]
        assert 'FORM=FORMATTED' in system_lines[0]
        assert 'COMPILE SESTATIC' in exec_lines
        assert "ALTER 'END'" in exec_lines
        assert any('OUTPUT4' in line and 'KGG' in line and 'MGG' in line
                   for line in exec_lines)

    def test_dmap_alter_op2(self):
        """Test generating DMAP ALTER for OUTPUT2 export.

        Verifies ASSIGN OUTPUT2 and COMPILE for SOL 103.
        """
        system_lines, exec_lines = make_matrix_export_alter(
            matrices=['KGG', 'KDICT', 'GPDT'],
            output_filename='matrices.op2',
            output_format='op2',
            sol=103,
            unit=25,
        )
        assert len(system_lines) == 1
        assert "ASSIGN OUTPUT2='matrices.op2'" in system_lines[0]
        assert 'UNIT=25' in system_lines[0]
        assert 'COMPILE SEMODES' in exec_lines
        assert any('OUTPUT2' in line and 'KGG' in line for line in exec_lines)

    def test_dmap_alter_custom_compile_and_location(self):
        """Test overriding compile_module and alter_location.

        Verifies that custom parameters override SOL-based defaults.
        """
        system_lines, exec_lines = make_matrix_export_alter(
            matrices=['KELM', 'KDICT'],
            output_format='op4',
            sol=101,
            compile_module='SEMG',
            alter_location='345',
        )
        assert 'COMPILE SEMG' in exec_lines
        assert "ALTER '345'" in exec_lines

    def test_dmap_alter_add_to_model(self):
        """Test adding DMAP ALTER to a BDF model.

        Verifies:
          - ASSIGN line added to system_command_lines
          - COMPILE/ALTER/OUTPUT4 appended to executive_control_lines
          - Written BDF contains the ALTER
        """
        model = BDF(debug=False)
        model.sol = 101
        add_matrix_export_alter(
            model,
            matrices=['KGG', 'MGG', 'KDICT'],
            output_filename='kgg.op4',
            output_format='op4',
            unit=51,
        )
        # check system lines
        assert any("ASSIGN OUTPUT4='kgg.op4'" in line
                   for line in model.system_command_lines)

        # check executive control has ALTER
        exec_str = '\n'.join(model.executive_control_lines)
        assert 'COMPILE SESTATIC' in exec_str
        assert 'OUTPUT4' in exec_str

        # write and verify
        out = StringIO()
        model.write_bdf(out, close=False)
        written = out.getvalue()
        assert "ASSIGN OUTPUT4='kgg.op4'" in written
        assert 'COMPILE SESTATIC' in written
        assert 'OUTPUT4' in written
        assert 'KGG' in written

    def test_dmap_alter_many_matrices(self):
        """Test OUTPUT4 statement splitting for >5 matrices.

        OUTPUT4 supports at most 5 matrices per statement.
        7 matrices -> 2 OUTPUT4 lines (5 + 2).
        """
        system_lines, exec_lines = make_matrix_export_alter(
            matrices=['KGG', 'MGG', 'BGG', 'KDICT', 'MDICT', 'GPDT', 'EQEXIN'],
            output_filename='all.op4',
            output_format='op4',
            sol=101,
        )
        output4_lines = [l for l in exec_lines if l.startswith('OUTPUT4')]
        assert len(output4_lines) == 2
        assert 'KGG' in output4_lines[0]
        assert 'GPDT' in output4_lines[1]

    def test_dmap_alter_bad_sol(self):
        """Test that unsupported SOL raises ValueError."""
        with self.assertRaises(ValueError):
            make_matrix_export_alter(['KGG'], sol=999)

    def test_dmap_alter_bad_format(self):
        """Test that invalid output_format raises ValueError."""
        with self.assertRaises(ValueError):
            make_matrix_export_alter(['KGG'], output_format='hdf5')

    # --- NX MALTER tests ---

    def test_malter_kgg(self):
        """Test NX MALTER generation for KGG export point.

        Uses the predefined 'MALTER:(KGG, BGG, MGG, K4GG, PG)' label.
        """
        exec_lines = make_matrix_export_malter(
            matrices=['KGG', 'MGG'],
            malter_label='KGG',
        )
        assert exec_lines[0] == "MALTER 'MALTER:(KGG, BGG, MGG, K4GG, PG)' $"
        assert 'OUTPUT2 KGG//-3/-12/KGG// $' in exec_lines
        assert 'OUTPUT2 MGG//-3/-12/MGG// $' in exec_lines

    def test_malter_custom_label(self):
        """Test MALTER with a custom label string."""
        exec_lines = make_matrix_export_malter(
            matrices=['KAA'],
            malter_label='MALTER:AFTER PREFACE MODULES',
        )
        assert exec_lines[0] == "MALTER 'MALTER:AFTER PREFACE MODULES' $"
        assert 'OUTPUT2 KAA//-3/-12/KAA// $' in exec_lines

    def test_malter_custom_unit(self):
        """Test MALTER with custom output unit.

        Unit -45 means binary output to a user-assigned file.
        """
        exec_lines = make_matrix_export_malter(
            matrices=['KGG'],
            output_unit=-45,
        )
        assert 'OUTPUT2 KGG//-3/-45/KGG// $' in exec_lines

    def test_malter_add_to_model(self):
        """Test adding MALTER to BDF model and writing.

        Verifies MALTER and OUTPUT2 lines appear in written BDF.
        """
        model = BDF(debug=False)
        model.sol = 101
        add_matrix_export_malter(model, ['KGG', 'MGG'])

        out = StringIO()
        model.write_bdf(out, close=False)
        written = out.getvalue()
        assert "MALTER 'MALTER:(KGG, BGG, MGG, K4GG, PG)'" in written
        assert 'OUTPUT2 KGG//-3/-12/KGG//' in written
        assert 'OUTPUT2 MGG//-3/-12/MGG//' in written


if __name__ == '__main__':  # pragma: no cover
    unittest.main()
