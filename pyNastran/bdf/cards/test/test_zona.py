import os
from io import StringIO
from pathlib import Path
# from collections import defaultdict
import unittest
# from typing import Optional, Any

# import numpy as np
from cpylog import SimpleLogger

import pyNastran
from pyNastran.utils import print_bad_path
from pyNastran.bdf.bdf import read_bdf
from pyNastran.bdf.cards.test.utils import save_load_deck


IS_MATPLOTLIB = False
if IS_MATPLOTLIB:
    import matplotlib.pyplot as plt

PKG_PATH = Path(pyNastran.__path__[0])
MODEL_PATH = PKG_PATH / '..' / 'models'
TEST_PATH = PKG_PATH / 'bdf' / 'cards' / 'test'
ZONA_PATH = MODEL_PATH / 'aero' / 'zona'

EXAMPLES_DIR = PKG_PATH / 'bdf' / 'cards' / 'aero' / 'examples'

ASE_DIR = EXAMPLES_DIR / 'ase'
FLUTTER_DIR = EXAMPLES_DIR / 'flutter'
GLOADS_DIR = EXAMPLES_DIR / 'gloads'
MLOADS_DIR = EXAMPLES_DIR / 'mloads'
assert FLUTTER_DIR.exists(), print_bad_path(FLUTTER_DIR)


class TestAeroZona(unittest.TestCase):
    def test_zona_1(self):
        """zona explicit test"""
        log = SimpleLogger(level='error', encoding='utf-8')  # lots of zona errors
        bdf_filename = ZONA_PATH / 'f16_ma41.bdf'
        model = read_bdf(bdf_filename, xref=False, debug=None, log=log)
        model.safe_cross_reference()
        save_load_deck(model, xref='safe',
                       run_renumber=False, run_convert=False, run_remove_unused=False,
                       run_save_load=False, run_save_load_hdf5=False, run_mass_properties=False,
                       run_test_bdf=False, run_op2_writer=False, run_export_caero=False,
                       stringify=True)
        with self.assertRaises(NotImplementedError):
            model.zona.convert_to_nastran()

    def test_zona_2(self):
        """zona explicit test"""
        log = SimpleLogger(level='error', encoding='utf-8')  # lots of zona errors
        bdf_filename = ZONA_PATH / 'ztran.bdf'
        model = read_bdf(bdf_filename, xref=False, debug=None, log=log)
        model.safe_cross_reference()
        save_load_deck(model, xref='safe',
                       run_renumber=False, run_convert=False, run_remove_unused=False,
                       run_save_load=False, run_save_load_hdf5=False, run_mass_properties=False,
                       run_export_caero=False, run_test_bdf=False, run_op2_writer=False)
        model.zona.convert_to_nastran()

    def test_zona_model_1(self):
        """totally fake zona model"""
        bdf_file = get_zona_model()

        model = read_bdf(bdf_filename=bdf_file, validate=True, xref=True, punch=False,
                         skip_cards=None, read_cards=None, encoding=None,
                         log=None, debug=None, mode='zona')
        #with self.assertRaises(AttributeError):

        model.uncross_reference()
        model.write_bdf('zona.bdf')
        model.safe_cross_reference()
        model.write_bdf('zona.bdf')

        bdf_file.seek(0)
        model.clear_attributes()
        model2 = read_bdf('zona.bdf', debug=None)
        os.remove('zona.bdf')
        model2.zona.convert_to_nastran()

    def test_zona_flutter_case1_in(self):
        zona_filename = FLUTTER_DIR / 'case1' / 'ha145e.inp'
        model = read_bdf(zona_filename, xref=True, debug=True)
        model.write_bdf(FLUTTER_DIR / 'zona.inp')

    def test_zona_flutter_case2_in(self):
        zona_filename = FLUTTER_DIR / 'case2' / 'crop.inp'
        model = read_bdf(zona_filename, xref=True, debug=True)
        model.write_bdf(FLUTTER_DIR / 'zona.inp')

    def test_zona_flutter_case3_in(self):
        zona_filename = FLUTTER_DIR / 'case3' / 'ha145fb.inp'
        model = read_bdf(zona_filename, xref=True, debug=True)
        model.write_bdf(FLUTTER_DIR / 'zona.inp')

    def test_zona_flutter_case4_in(self):
        zona_filename = FLUTTER_DIR / 'case4' / 'ha145g.inp'
        model = read_bdf(zona_filename, xref=True, debug=True)
        model.write_bdf(FLUTTER_DIR / 'zona.inp')

    def test_zona_flutter_case5_in(self):
        zona_filename = FLUTTER_DIR / 'case5' / 'f16ma41.inp'
        model = read_bdf(zona_filename, xref=False, debug=True)
        model.zona.cross_reference()
        model.write_bdf(FLUTTER_DIR / 'zona.inp')

    def test_zona_flutter_case6_in_trim(self):
        zona_filename = FLUTTER_DIR / 'case6' / 'agard_trim.inp'
        model = read_bdf(zona_filename, xref=False, debug=True)
        model.zona.cross_reference()
        model.write_bdf(FLUTTER_DIR / 'zona.inp')

    def test_zona_flutter_case6_in_tran(self):
        zona_filename = FLUTTER_DIR / 'case6' / 'agardztran.inp'
        model = read_bdf(zona_filename, xref=False, debug=True)
        model.zona.cross_reference()
        model.write_bdf(MLOADS_DIR / 'zona.inp')

    def test_zona_flutter_case7_in(self):
        zona_filename = FLUTTER_DIR / 'case7' / 'agardztaw.inp'
        model = read_bdf(zona_filename, xref=False, debug=True)
        model.zona.cross_reference()
        model.write_bdf(MLOADS_DIR / 'zona.inp')

    def test_zona_mloads_case1_in(self):
        zona_filename = MLOADS_DIR / 'case1' / 'm144open.inp'
        model = read_bdf(zona_filename, xref=True, debug=True)
        # model.zona.cross_reference()
        model.write_bdf(MLOADS_DIR / 'zona.inp')
        assert len(model.zona.mloads)
        plot = False
        # if plot:
        #     import matplotlib.pyplot as plt
        #     fig = plt.figure()
        #     for extid, mloads in model.zona.mloads.items():
        #         mloads.plot(fig)
        #         plt.show()
        #         break

    def test_zona_mloads_case2a_in(self):
        zona_filename = MLOADS_DIR / 'case2' / 'm144_trim.inp'
        model = read_bdf(zona_filename, xref=True, debug=True)
        model.write_bdf(MLOADS_DIR / 'zona.inp')

    def test_zona_mloads_case2b_in(self):
        zona_filename = MLOADS_DIR / 'case2' / 'm144clos.inp'
        model = read_bdf(zona_filename, xref=True, debug=True)
        model.write_bdf(MLOADS_DIR / 'zona.inp')

    def test_zona_ase_case1(self):
        zona_filename = ASE_DIR / 'case1' / 'cropase.inp'
        model = read_bdf(zona_filename, mode='zona', xref=True, debug=True)
        model.write_bdf(ASE_DIR / 'zona.inp')

    def test_zona_gloads_case1(self):
        zona_filename = GLOADS_DIR / 'case1' / 'kussner.inp'
        model = read_bdf(zona_filename, mode='zona', xref=False, debug=True)
        model.zona.cross_reference()
        model.write_bdf(GLOADS_DIR / 'zona.inp')

    def test_zona_gloads_case2(self):
        zona_filename = GLOADS_DIR / 'case2' / 'gbj_dgust.inp'
        model = read_bdf(zona_filename, mode='zona', xref=False, debug=True)
        model.zona.cross_reference()
        model.write_bdf(GLOADS_DIR / 'zona.inp')

    def test_zona_gloads_case3(self):
        zona_filename = GLOADS_DIR / 'case3' / 'gbj_cgust.inp'
        model = read_bdf(zona_filename, mode='zona', xref=False, debug=True)
        model.zona.cross_reference()
        model.write_bdf(GLOADS_DIR / 'zona.inp')

    # def test_zona_gloads_case4a(self):
    #     zona_filename = GLOADS_DIR / 'case4' / 'cgust_md.inp'
    #     model = read_bdf(zona_filename, mode='zona', xref=False, debug=True)
    #     model.zona.cross_reference()
    #     model.write_bdf(GLOADS_DIR / 'zona.inp')

    # def test_zona_gloads_case4b(self):
    #     zona_filename = GLOADS_DIR / 'case4' / 'cgust_sof.inp'
    #     model = read_bdf(zona_filename, mode='zona', xref=False, debug=True)
    #     model.zona.cross_reference()
    #     model.write_bdf(GLOADS_DIR / 'zona.inp')


def get_zona_model() -> StringIO:
    bdf_file = StringIO()
    bdf_file.write(
        '$ pyNastran: version=zona\n'
        'CEND\n'
        'BEGIN BULK\n'
        #'$       acsid, rcsid, cref, bref, sref, symxz, symxy\n'
        #'AEROZ, 10,     0,     1.,   10.,  100., YES\n'
        '$AEROZ  ACSID XZSYM FLIP FMMUNIT FMLUNIT REFC   REFB   REFS\n'
        '$       REFX  REFY  REFZ\n'
        'AEROZ,  0,    YES,  NO,  SLIN,   IN,      22.73,59.394,1175.8\n'
        ',       59.53,0.0,  0.0\n'

        '$       label, type, cid, PANLST, setg, actid\n'
        'AESURFZ,FLAP,  ASYM, 1,   10,       20,   0\n'
        #'AESURFZ,FLAP,  SYM,  1,  10,       20,   0\n'
        'CORD2R, 1,0, 0.,0.,0., 0.,0.,1.,\n'
        ',1.,0.,0.\n'
        '$BODY7,ID,LABEL,IPBODY7, ACOORD, NSEG, IDMESH1\n'
        'BODY7, 1, FUSE,        ,      2,     , 1\n'
        'PANLST3,10, FUSE, \n'
        '$       id,naxial,nradial, \n'
        'SEGMESH,1, 4,     3,       \n'

        # ITYPEi = 1 (Body of Revolution):
        #    Xi, CAMi, YRi
        # ITYPEi = 2 (Elliptical Body):
        #    Xi, YRi, ZRi
        # ITYPEi = 3 (Arbitrary Body):
        #    Xi, IDYi, IDZi
        '$       itype, x1, cam, yr1, zr1, idy1, idz1 \n'
        ',       1,        ,  1.,  1.,    ,\n'
        ',       2,      1.,    ,  1.,  2.,\n'
        ',       3,      2.,    ,    ,    , 13,   14   \n'
        ',       3,      3.,    ,    ,    , 13,   14   \n'

        # y
        'AEFACT,13, 1., 0.,  0.,-1.\n'
        'AEFACT,14, 0., 1., -1., 0.\n'
        '$ MKAEROZ, ID, MACH, METHOD, IDFLT\n'
        'MKAEROZ,   101, 0.8, -1,     -1,  \n'
        '$ TRIM, ID, MKAEROZ, Q,   LABEL1, UX1,    CGX, CGY,\n'
        'TRIM, 100,  101,     42., ALPHA,  5., 0., 0.,  0.,\n'
        '$CGZ, WEIGHT, Ixx, Ixy, Iyy, Ixz, Iyz, Izz\n'
        ',0.,  1e4,    1e3, 1e3, 1e5, 1e3, 1e3, 1e4\n'
        '$TRUE/G, NX,     NY,  NZ,  P,       Q,   R, \n'
        ', TRUE,  FREE, NONE,  32., FREE, NONE, 42., \n'
        '$var, value\n'
        ',17,  1.0,\n'
        '$\n'
        'TRIMVAR,17,VAR\n'
        '$\n'
        '$trimlnk,id,sym, ?,  ?\n'
        'TRIMLNK,10,SYM, -1, 17\n'
        'ACOORD, 2, 0.,0.,0., 1.0,0.\n'
        '$       ID,    MODEL, CP, PANLST, SETG, DZ, EPS\n'
        'SPLINE1,100,        ,   ,    422, 423,\n'
        '$,      NELEM, MELEM\n'
        '$,      10,    10\n'
        'PANLST3,422, FUSE, \n'
        '$       id,naxial,nradial, \n'

        #'$       ID,   MODEL, PANLST, SETG,\n'
        #'SPLINE2,1000,      ,    422,  423,\n'

        '$       ID,   MODEL, CP, PANLST, SETG,\n'
        'SPLINE3,1200,      ,   ,    422,  423,\n'

        'SET1,423,10\n'
        'GRID,10,,0.,0.,0.\n'
        'GRID,11,,1.,0.,0.\n'
        'CONROD,100, 10,11, 101,1.0\n'
        'MAT1,101,3.0e7,,0.3\n'
    )
    bdf_file.seek(0)
    return bdf_file


