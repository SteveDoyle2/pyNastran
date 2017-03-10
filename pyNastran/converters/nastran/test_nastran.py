from six.moves import range
import os
from numpy import array_equal, allclose
import unittest

import pyNastran
from pyNastran.bdf.bdf import read_bdf

from pyNastran.converters.nastran.nastran_to_tecplot import nastran_to_tecplot
from pyNastran.converters.nastran.nastran_to_cart3d import nastran_to_cart3d
from pyNastran.converters.nastran.nastran_to_stl import nastran_to_stl
from pyNastran.converters.nastran.nastran_to_ugrid import nastran_to_ugrid
from pyNastran.converters.aflr.ugrid.ugrid_reader import read_ugrid
from pyNastran.utils.log import get_logger

pkg_path = pyNastran.__path__[0]
model_path = os.path.join(pkg_path, '../', 'models')

class TestNastran(unittest.TestCase):

    def test_nastran_to_ugrid_01(self):
        bdf_filename = os.path.join(model_path, 'solid_bending', 'solid_bending.bdf')

        size = 8
        debug = False
        log = get_logger(log=None, level='warning', encoding='utf-8')
        model = read_bdf(bdf_filename, log=log, debug=debug)
        #log = model.log
        #model.get_element_faces()
        skin_bdf_filename = os.path.join(model_path, 'solid_bending', 'solid_skin.bdf')
        model.write_skin_solid_faces(skin_bdf_filename, write_solids=True,
                                    write_shells=True,
                                    size=size, is_double=False, encoding=None)

        bdf_model = read_bdf(skin_bdf_filename, log=log, debug=debug)
        ugrid_filename_out = os.path.join(model_path, 'solid_bending', 'solid_skin.b8.ugrid')
        nastran_to_ugrid(bdf_model, ugrid_filename_out, properties=None,
                         check_shells=True, check_solids=True)
        ugrid = read_ugrid(ugrid_filename_out, encoding=None, log=log,
                           debug=debug)

        skin_bdf_filename2 = os.path.join(model_path, 'solid_bending', 'solid_skin2.bdf')
        ugrid.write_bdf(skin_bdf_filename2, include_shells=True, include_solids=True,
                  convert_pyram_to_penta=True, encoding=None,
                  size=size, is_double=False)
        model2 = read_bdf(skin_bdf_filename2, log=log, debug=debug)

        os.remove(ugrid_filename_out)
        os.remove(skin_bdf_filename)
        os.remove(skin_bdf_filename2)

    def test_nastran_to_stl(self):
        """tests nastran_to_stl"""
        bdf_filename = os.path.join(model_path, 'plate', 'plate.bdf')
        stl_filename = os.path.join(model_path, 'plate', 'plate.stl')
        log = get_logger(log=None, level='warning', encoding='utf-8')
        nastran_to_stl(bdf_filename, stl_filename, is_binary=False, log=log)

if __name__ == '__main__':  # pragma: no cover
    import time
    t0 = time.time()
    unittest.main()
    #test_1()
    #test_2()
    #test_3()
    print("dt = %s" % (time.time() - t0))
