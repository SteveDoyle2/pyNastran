import unittest
from io import StringIO
from pyNastran.bdf.bdf import BDF

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
        #print(model)

if __name__ == '__main__':  # pragma: no cover
    unittest.main()
