from __future__ import print_function
import os
import unittest
import pyNastran
from pyNastran.converters.nastran.dev_vectorized2.bdf_vectorized import read_bdf

pkg_path = pyNastran.__path__[0]
model_path = os.path.join(pkg_path, '../', 'models')


class TestVectorized(unittest.TestCase):
    def test_solid_bending(self):
        bdf_filename = os.path.join(model_path, 'solid_bending', 'solid_bending.bdf')
        model = read_bdf(bdf_filename, validate=True, xref=False, punch=False,
                         skip_cards=None, encoding=None, log=None, debug=True,
                         mode='msc')
        #print(model.get_bdf_stats())

        #model.grids[10] = GRID(10, [0., 0., 0.])
        print(model.grid)
        print(model.grid[10]) # nid or index?
        print(model.grid.get_grid_by_nid(10))
        out_filename = 'spike.bdf'
        model.write_bdf(out_filename, encoding=None, size=8, is_double=False,
                        interspersed=False, enddata=None,
                        close=True)

    def test_bwb(self):
        bdf_filename = os.path.join(model_path, 'bwb', 'BWB_saero.bdf')
        model = read_bdf(bdf_filename, validate=True, xref=False, punch=False,
                         skip_cards=None, encoding=None, log=None, debug=True,
                         mode='msc')
        #print(model.get_bdf_stats())

        out_filename = 'spike.bdf'
        model.write_bdf(out_filename, encoding=None, size=8, is_double=False,
                        interspersed=False, enddata=None,
                        close=True)

    def test_static_elements(self):
        bdf_filename = os.path.join(model_path, 'elements', 'static_elements.bdf')
        model = read_bdf(bdf_filename, validate=True, xref=False, punch=False,
                         skip_cards=None, encoding=None, log=None, debug=True,
                         mode='msc')

        print(model.cquad4)
        print(model.ctria3)
        print(model.shells)
        print(model.solids)
        print(model.elements)


        print(model.ctetra4)
        print(model.cpenta6)
        print(model.chexa8)
        print(model.cpyram5)
        print(model.ctetra10)
        print(model.cpenta15)
        print(model.chexa20)
        print(model.cpyram13)

        print(model.celas1)
        print(model.celas2)
        print(model.celas3)
        print(model.celas4)

        print(model.cdamp1)
        print(model.cdamp2)
        print(model.cdamp3)
        print(model.cdamp4)

        print(model.conrod)
        print(model.crod)
        print(model.ctube)

        print(model.elements2)
        print(len(model.elements2))
        print(model.get_bdf_stats())

        out_filename = 'spike.bdf'
        model.write_bdf(out_filename, encoding=None, size=8, is_double=False,
                        interspersed=False, enddata=None,
                        close=True)
        model.write_bdf(out_filename, encoding=None, size=16, is_double=False,
                        interspersed=False, enddata=None,
                        close=True)

if __name__ == '__main__':  # pragma: no cover
    unittest.main()
