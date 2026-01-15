import os
from io import StringIO
from pathlib import Path
import unittest

from cpylog import SimpleLogger

import pyNastran
from pyNastran.utils import print_bad_path
from pyNastran.bdf.bdf import read_bdf
from pyNastran.bdf.cards.test.utils import save_load_deck
from pyNastran.bdf.cards.aero.zaero import ZAERO, get_dicts


IS_MATPLOTLIB = False
if IS_MATPLOTLIB:
    import matplotlib.pyplot as plt

PKG_PATH = Path(pyNastran.__path__[0])
MODEL_PATH = PKG_PATH / '..' / 'models'
TEST_PATH = PKG_PATH / 'bdf' / 'cards' / 'test'
ZAERO_PATH = MODEL_PATH / 'aero' / 'zona'

EXAMPLES_DIR = PKG_PATH / 'bdf' / 'cards' / 'aero' / 'examples'

ASE_DIR = EXAMPLES_DIR / 'ase'
TRIM_DIR = EXAMPLES_DIR / 'trim'
FLUTTER_DIR = EXAMPLES_DIR / 'flutter'
GLOADS_DIR = EXAMPLES_DIR / 'gloads'
MLOADS_DIR = EXAMPLES_DIR / 'mloads'
assert FLUTTER_DIR.exists(), print_bad_path(FLUTTER_DIR)
dirname = Path(os.path.dirname(__file__))


class TestAeroZaero(unittest.TestCase):
    def test_zaero_1(self):
        """zaero explicit test"""
        log = SimpleLogger(level='error', encoding='utf-8')
        #log = SimpleLogger(level='debug', encoding='utf-8')
        bdf_filename = ZAERO_PATH / 'f16_ma41.bdf'
        model = read_bdf(bdf_filename, xref=False,
                         mode='zona', debug=None, log=log)
        model.zaero.uncross_reference()
        model.safe_cross_reference()
        # save_load_deck(
        #     model, xref='safe', nastran_format='zaero',
        #     run_renumber=False, run_convert=False, run_remove_unused=False,
        #     run_save_load=False, run_save_load_hdf5=False, run_mass_properties=False,
        #     run_test_bdf=False, run_op2_writer=False, run_export_caero=False,
        #     stringify=True)
        with self.assertRaises(NotImplementedError):
            model.zaero.convert_to_nastran()

    def _test_zaero_2(self):
        """zaero explicit test"""
        log = SimpleLogger(level='error', encoding='utf-8')  # lots of zaero errors
        bdf_filename = ZAERO_PATH / 'ztran.bdf'
        model = read_bdf(bdf_filename, xref=False, debug=None, log=log)
        model.safe_cross_reference()
        save_load_deck(model, xref='safe',
                       run_renumber=False, run_convert=False, run_remove_unused=False,
                       run_save_load=False, run_save_load_hdf5=False, run_mass_properties=False,
                       run_export_caero=False, run_test_bdf=False, run_op2_writer=False)
        model.zaero.convert_to_nastran()
        write_raw_fields(model.zaero)
        model.zaero.uncross_reference()

    def test_zaero_model_1(self):
        """totally fake zaero model"""
        bdf_file = get_zaero_model()

        model = read_bdf(
            bdf_filename=bdf_file, validate=True, xref=True, punch=False,
            skip_cards=None, read_cards=None, encoding=None,
            log=None, debug=False, mode='zaero')
        #with self.assertRaises(AttributeError):

        model.uncross_reference()
        model.write_bdf('zaero.bdf')
        model.safe_cross_reference()
        model.write_bdf('zaero.bdf')
        model.zaero.uncross_reference()
        # model.zaero.convert_to_nastran('zaero.bdf')

        bdf_file.seek(0)
        model.clear_attributes()
        model2 = read_bdf('zaero.bdf', debug=None)
        os.remove('zaero.bdf')
        write_raw_fields(model2.zaero)
        model2.zaero.convert_to_nastran()
        #model2.write_bdf('zaero2.bdf')

    def test_zaero_trim_case1_in(self):
        zaero_filename = TRIM_DIR / 'case1' / 'ha144d.inp'
        model = read_bdf(zaero_filename, xref=False, debug=False, mode='zaero')
        model.zaero.safe_cross_reference()
        model.zaero.uncross_reference()

        with self.assertRaises(AssertionError):
            model.cross_reference()
        write_raw_fields(model.zaero)
        model.write_bdf(TRIM_DIR / 'zaero.inp')

    def test_zaero_flutter_case1_in(self):
        zaero_filename = FLUTTER_DIR / 'case1' / 'ha145e.inp'
        model = read_bdf(zaero_filename, xref=True, debug=False)
        model.zaero.safe_cross_reference()
        write_raw_fields(model.zaero)
        model.write_bdf(FLUTTER_DIR / 'zaero.inp')
        model.zaero.uncross_reference()

    def test_zaero_flutter_case2_in(self):
        zaero_filename = FLUTTER_DIR / 'case2' / 'crop.inp'
        model = read_bdf(zaero_filename, xref=True, debug=False)
        model.zaero.safe_cross_reference()
        write_raw_fields(model.zaero)
        model.write_bdf(FLUTTER_DIR / 'zaero.inp')
        model.zaero.uncross_reference()

    def test_zaero_flutter_case3_in(self):
        zaero_filename = FLUTTER_DIR / 'case3' / 'ha145fb.inp'
        model = read_bdf(zaero_filename, xref=True, debug=False)
        model.zaero.safe_cross_reference()
        write_raw_fields(model.zaero)
        model.write_bdf(FLUTTER_DIR / 'zaero.inp')
        model.zaero.uncross_reference()

    def test_zaero_flutter_case4_in(self):
        zaero_filename = FLUTTER_DIR / 'case4' / 'ha145g.inp'
        model = read_bdf(zaero_filename, xref=True, debug=False)
        model.zaero.safe_cross_reference()
        write_raw_fields(model.zaero)
        model.write_bdf(FLUTTER_DIR / 'zaero.inp')
        model.zaero.uncross_reference()

    def test_zaero_flutter_case5_in(self):
        zaero_filename = FLUTTER_DIR / 'case5' / 'f16ma41.inp'
        model = read_bdf(zaero_filename, xref=False, debug=False)
        model.zaero.cross_reference()
        model.zaero.safe_cross_reference()
        write_raw_fields(model.zaero)
        model.write_bdf(FLUTTER_DIR / 'zaero.inp')
        model.zaero.uncross_reference()

    def test_zaero_flutter_case6_in_trim(self):
        zaero_filename = FLUTTER_DIR / 'case6' / 'agard_trim.inp'
        model = read_bdf(zaero_filename, xref=False, debug=False)
        model.zaero.cross_reference()
        model.zaero.safe_cross_reference()
        write_raw_fields(model.zaero)
        model.write_bdf(FLUTTER_DIR / 'zaero.inp')
        model.zaero.uncross_reference()

    def test_zaero_flutter_case6_in_tran(self):
        zaero_filename = FLUTTER_DIR / 'case6' / 'agardztran.inp'
        model = read_bdf(zaero_filename, xref=False, debug=False)
        model.zaero.cross_reference()
        model.zaero.safe_cross_reference()
        write_raw_fields(model.zaero)
        model.write_bdf(MLOADS_DIR / 'zaero.inp')
        model.zaero.uncross_reference()

    def test_zaero_flutter_case7_in(self):
        zaero_filename = FLUTTER_DIR / 'case7' / 'agardztaw.inp'
        model = read_bdf(zaero_filename, xref=False, debug=False)
        model.zaero.cross_reference()
        model.zaero.safe_cross_reference()
        write_raw_fields(model.zaero)
        model.write_bdf(MLOADS_DIR / 'zaero.inp')
        model.zaero.uncross_reference()

    def test_zaero_mloads_case1_in(self):
        zaero_filename = MLOADS_DIR / 'case1' / 'm144open.inp'
        model = read_bdf(zaero_filename, xref=True, debug=False)
        model.zaero.uncross_reference()
        # model.zaero.cross_reference()
        model.zaero.safe_cross_reference()
        write_raw_fields(model.zaero)
        model.write_bdf(MLOADS_DIR / 'zaero.inp')
        assert len(model.zaero.mloads)
        plot = False
        # if plot:
        #     import matplotlib.pyplot as plt
        #     fig = plt.figure()
        #     for extid, mloads in model.zaero.mloads.items():
        #         mloads.plot(fig)
        #         plt.show()
        #         break

    def test_zaero_mloads_case2a_in(self):
        zaero_filename = MLOADS_DIR / 'case2' / 'm144_trim.inp'
        model = read_bdf(zaero_filename, xref=True, debug=False)
        model.zaero.safe_cross_reference()
        write_raw_fields(model.zaero)
        model.write_bdf(MLOADS_DIR / 'zaero.inp')
        model.zaero.uncross_reference()

    def test_zaero_mloads_case2b_in(self):
        zaero_filename = MLOADS_DIR / 'case2' / 'm144clos.inp'
        model = read_bdf(zaero_filename, xref=True, debug=False)
        model.zaero.safe_cross_reference()
        write_raw_fields(model.zaero)
        model.write_bdf(MLOADS_DIR / 'zaero.inp')
        model.zaero.uncross_reference()

    def test_zaero_ase_case1(self):
        zaero_filename = ASE_DIR / 'case1' / 'cropase.inp'
        model = read_bdf(zaero_filename, mode='zaero', xref=True, debug=False)
        model.zaero.safe_cross_reference()
        write_raw_fields(model.zaero)
        model.write_bdf(ASE_DIR / 'zaero.inp')
        model.zaero.uncross_reference()

    def test_zaero_ase_case2(self):
        zaero_filename = ASE_DIR / 'case2' / 'gafa.inp'
        model = read_bdf(zaero_filename, mode='zaero', xref=True, debug=False)
        model.zaero.safe_cross_reference()
        write_raw_fields(model.zaero)
        model.write_bdf(ASE_DIR / 'zaero.inp')
        model.zaero.uncross_reference()

    def test_zaero_gloads_case1(self):
        zaero_filename = GLOADS_DIR / 'case1' / 'kussner.inp'
        model = read_bdf(zaero_filename, mode='zaero', xref=False, debug=False)
        model.zaero.cross_reference()
        model.zaero.safe_cross_reference()
        write_raw_fields(model.zaero)
        model.write_bdf(GLOADS_DIR / 'zaero.inp')
        model.zaero.uncross_reference()

    def test_zaero_gloads_case2(self):
        zaero_filename = GLOADS_DIR / 'case2' / 'gbj_dgust.inp'
        model = read_bdf(zaero_filename, mode='zaero', xref=False, debug=False)
        model.zaero.cross_reference()
        model.zaero.safe_cross_reference()
        write_raw_fields(model.zaero)
        model.write_bdf(GLOADS_DIR / 'zaero.inp')
        model.zaero.uncross_reference()

    def test_zaero_gloads_case3(self):
        zaero_filename = GLOADS_DIR / 'case3' / 'gbj_cgust.inp'
        model = read_bdf(zaero_filename, mode='zaero', xref=False, debug=False)
        model.zaero.cross_reference()
        model.zaero.safe_cross_reference()
        write_raw_fields(model.zaero)
        model.write_bdf(GLOADS_DIR / 'zaero.inp')
        model.zaero.uncross_reference()

    def test_zaero_gloads_case4a(self):
        zaero_filename = GLOADS_DIR / 'case4' / 'cgust_md.inp'
        model = read_bdf(zaero_filename, mode='zaero', xref=False, debug=False)
        model.zaero.cross_reference()
        model.write_bdf(GLOADS_DIR / 'zaero.inp')
        model.zaero.uncross_reference()

    def test_zaero_gloads_case4b(self):
        zaero_filename = GLOADS_DIR / 'case4' / 'cgust_sof.inp'
        model = read_bdf(zaero_filename, mode='zaero', xref=False, debug=False)
        model.zaero.cross_reference()
        model.write_bdf(GLOADS_DIR / 'zaero.inp')
        model.zaero.uncross_reference()


def get_zaero_model() -> StringIO:
    bdf_file = StringIO()
    bdf_file.write(
        '$ pyNastran: version=zona\n'
        'CEND\n'
        'BEGIN BULK\n'
        # ATMOS IDATM AMMUNIT AMLUNIT AMTUNIT
        #       ALT1  SOUND1  DEN1    TEMP1  ALT2 SOUND2 DEN2 TEMP2
        #       ALTi  SOUNDi  DENi    TEMPi -etc-
        'ATMOS, 12,   SLIN,   IN,     F\n'
        ',      0.0,  1000.,  1.0,    2.0,  10000.0, 900.,  1.0,    2.0\n'
        # FIXMATM SETID   IDMK  IDATM FTMUNIT FTLUNIT VREF FLUTTF PRINT CONT
        #         ALT1    ALT2  ALTi -etc-
        'FIXMATM, 100,    101,   12,  slug,   ft,     1.0, -1,    0\n'
        ',       -10000., 0., 10000., 20000., 30000.\n'
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
        '$ TRIM, ID, MKAEROZ, Q,   IDOBJ, IDCONS, CGX, CGY, CGZ\n'
        'TRIM, 100,  101,     42., 0,     0,      5., 0., 0.\n'
        '$WTMASS, WEIGHT, Ixx, Ixy, Iyy, Ixz, Iyz, Izz\n'
        ',1.,     1e4,    1e3, 1e3, 1e5, 1e3, 1e3, 1e4\n'
        '$TRUE/G, NX,     NY,  NZ,  PDOT, QDOT, RDOT,\n'
        ', TRUE,  FREE, NONE,  32., FREE, NONE, 42.,\n'
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


def write_raw_fields(zaero: ZAERO):
    singletons, dicts, dicts_list = get_dicts(zaero, 'write')
    for panlsts in zaero.panlsts.values():
        for panlst in panlsts:
            panlst.raw_fields()

    for dicti in dicts:
        for value in dicti.values():
            value.raw_fields()
