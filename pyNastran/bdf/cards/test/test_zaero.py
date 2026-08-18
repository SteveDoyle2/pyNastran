import os
from io import StringIO
from pathlib import Path
import unittest

import numpy as np

from cpylog import SimpleLogger
import pyNastran
from pyNastran.bdf.cards.aero.zaero_interface import graphviz_interface
from pyNastran.utils import print_bad_path
from pyNastran.bdf.bdf import read_bdf, BDF
from pyNastran.bdf.cards.test.utils import save_load_deck
from pyNastran.bdf.cards.aero.zaero import ZAERO, get_dicts
from pyNastran.bdf.cards.aero.zaero_cards.cards import AEROZ
from pyNastran.bdf.cards.aero.zaero_cards.geometry import CAERO7
from pyNastran.bdf.bdf_interface.bdf_card import BDFCard
from pyNastran.bdf.cards.aero.zaero_cards.ase import (
    MIMOTF,
    SENSR,
    GAIN,
    SUMBLK,
    DEADBN,
    DELAY_ZAERO,
    FILTFL,
    LIMTR,
)


PKG_PATH = Path(pyNastran.__path__[0])
MODEL_PATH = PKG_PATH / ".." / "models"
TEST_PATH = PKG_PATH / "bdf" / "cards" / "test"
ZAERO_PATH = MODEL_PATH / "aero" / "zona"

EXAMPLES_DIR = PKG_PATH / "bdf" / "cards" / "aero" / "examples"

ASE_DIR = EXAMPLES_DIR / "ase"
TRIM_DIR = EXAMPLES_DIR / "trim"
FLUTTER_DIR = EXAMPLES_DIR / "flutter"
GLOADS_DIR = EXAMPLES_DIR / "gloads"
MLOADS_DIR = EXAMPLES_DIR / "mloads"
assert FLUTTER_DIR.exists(), print_bad_path(FLUTTER_DIR)
dirname = Path(os.path.dirname(__file__))


class TestAeroZaero(unittest.TestCase):
    def setUp(self):
        self._original_digraph = graphviz_interface.Digraph
        self._original_enf = graphviz_interface.ExecutableNotFound
        graphviz_interface.Digraph = graphviz_interface.FakeDigraph
        graphviz_interface.ExecutableNotFound = graphviz_interface.FakeExecutableNotFound

    def tearDown(self):
        graphviz_interface.Digraph = self._original_digraph
        graphviz_interface.ExecutableNotFound = self._original_enf

    def test_zaero_mloads(self):
        bdf_file = get_mloads_file()
        log = SimpleLogger(level="warning", encoding="utf-8")
        model = read_bdf(bdf_file, xref=True, mode="zona", debug=None, log=log)
        model.zaero.view_block_diagram(subcase_id=-1)
        model.zaero.uncross_reference()
        model.safe_cross_reference()
        model.zaero.convert_to_nastran()

    def test_zaero_add_methods(self):
        """test all add_* convenience methods on the ZAERO class"""
        model = BDF(debug=False)
        model.sol = 145
        zaero = ZAERO(model)
        model.zaero = zaero

        # --- atmosphere ---
        atmos = zaero.add_atmos(1, "SLUG", "IN", "R", [0.0, 1000.0, 0.001, 500.0])
        assert atmos.type == "ATMOS"
        assert 1 in zaero.atmos

        fixhatm = zaero.add_fixhatm(2, 10000.0, 1, "SLUG", "IN", 1.0, 0, 0, [1])
        assert fixhatm.type == "FIXHATM"
        assert 2 in zaero.flutter_table

        fixmatm = zaero.add_fixmatm(3, 1, [5000.0, 10000.0], atm_id=1, mass_unit="SLUG", length_unit="IN", vref=1.0,
                                    fluttf_id=0, print_flag=0)
        assert fixmatm.type == "FIXMATM"
        assert 3 in zaero.flutter_table

        fixmach = zaero.add_fixmach(4, 1, "SLUG", "IN", 0, 0, [100.0, 200.0], [0.001, 0.002])
        assert fixmach.type == "FIXMACH"
        assert 4 in zaero.flutter_table

        fixmden = zaero.add_fixmden(5, 1, 0.001, "SLUG", "IN", 0, 0, [100.0, 200.0])
        assert fixmden.type == "FIXMDEN"
        assert 5 in zaero.flutter_table

        # --- geometry ---
        splinem = zaero.add_splinem("SAVE", "spline.dat")
        assert splinem.type == "SPLINEM"
        assert zaero.splinem is not None

        panlst1 = zaero.add_panlst1(10, 101, 101, 110)
        assert panlst1.type == "PANLST1"
        assert 10 in zaero.panlsts

        panlst2 = zaero.add_panlst2(11, 101, [101, 102, 103])
        assert panlst2.type == "PANLST2"
        assert 11 in zaero.panlsts

        panlst3 = zaero.add_panlst3(12, ["WING", "TAIL"])
        assert panlst3.type == "PANLST3"
        assert 12 in zaero.panlsts

        pafoil7 = zaero.add_pafoil7(20, 1, 2, 3, 0.01, 4, 5, 0.005)
        assert pafoil7.type == "PAFOIL7"
        assert 20 in zaero.pafoil

        pafoil8 = zaero.add_pafoil8(21, 1, 2, 3, 0.01, 4)
        assert pafoil8.type == "PAFOIL8"
        assert 21 in zaero.pafoil

        aesurfz = zaero.add_aesurfz("ELEV", "SYM", 0, 10, 0, 0)
        assert aesurfz.type == "AESURFZ"
        assert "ELEV" in model.aesurf

        aeslink = zaero.add_aeslink("LNK1", "SYM", 100, ["ELEV"], [1.0])
        assert aeslink.type == "AESLINK"
        assert "LNK1" in zaero.aeslink

        attach = zaero.add_attach(30, "WING", 10, 100)
        assert attach.type == "ATTACH"
        assert 30 in zaero.attach

        # --- plot ---
        pltmode = zaero.add_pltmode(40, "SYM", 1, "TECPLOT", "mode.plt", max_disp=10.0)
        assert pltmode.type == "PLTMODE"
        assert 40 in zaero.pltmode

        pltaero = zaero.add_pltaero(41, "YES", 0, "TECPLOT", "aero.plt")
        assert pltaero.type == "PLTAERO"
        assert 41 in zaero.pltaero

        pltvg = zaero.add_pltvg(42, 1, "V", "vg.plt")
        assert pltvg.type == "PLTVG"
        assert 42 in zaero.pltvg

        pltcp = zaero.add_pltcp(43, "SYM", 1, 1, 1, "TECPLOT", "cp.plt")
        assert pltcp.type == "PLTCP"
        assert 43 in zaero.pltcp

        pltmist = zaero.add_pltmist(44, 1, 1, 1, 1, "TECPLOT", "mist.plt")
        assert pltmist.type == "PLTMIST"
        assert 44 in zaero.pltmist

        pltsurf = zaero.add_pltsurf(45, "WING", "TECPLOT", "surf.plt")
        assert pltsurf.type == "PLTSURF"
        assert 45 in zaero.pltsurf

        pltflut = zaero.add_pltflut(46, "TECPLOT", "flut.plt")
        assert pltflut.type == "PLTFLUT"
        assert 46 in zaero.pltflut

        pltbode = zaero.add_pltbode(47, 1, 0.1, 100.0, 50, "bode.plt")
        assert pltbode.type == "PLTBODE"
        assert 47 in zaero.pltbode

        plttime = zaero.add_plttime(48, 1, 0.0, 1.0, 100, "FORCE", "time.plt", "aero_time.plt")
        assert plttime.type == "PLTTIME"
        assert 48 in zaero.plttime

        plttrim = zaero.add_plttrim(49, 1, "FORCE", "trim.plt", "aero_trim.plt")
        assert plttrim.type == "PLTTRIM"
        assert 49 in zaero.plttrim

        # --- flutter ---
        mkaeroz = zaero.add_mkaeroz(50, 0.8, 1, "aero.dat", 0, [0.01, 0.1, 1.0])
        assert mkaeroz.type == "MKAEROZ"
        assert 50 in zaero.mkaeroz

        # --- trim ---
        trimvar = zaero.add_trimvar(60, "ALPHA", -10.0, 10.0, 0, "NONE", 0)
        assert trimvar.type == "TRIMVAR"
        assert 60 in zaero.trimvar

        trimlnk = zaero.add_trimlnk(61, "SYM", [1.0, -1.0], [60, 60])
        assert trimlnk.type == "TRIMLNK"
        assert 61 in zaero.trimlnk

        trimfnc = zaero.add_trimfnc(62, "USER", "LIFT", "YES", 0, 0, "test")
        assert trimfnc.type == "TRIMFNC"
        assert 62 in zaero.trimfnc

        # --- maneuver loads ---
        mloads = zaero.add_mloads(70, 1, 1, 0, 0, 0, 0, 0, 100.0, "SAVE")
        assert mloads.type == "MLOADS"
        assert 70 in zaero.mloads

        extinp = zaero.add_extinp(71, "", 100, 1, "CMD1")
        assert extinp.type == "EXTINP"
        assert 71 in zaero.extinp

        extout = zaero.add_extout(72, "", 100, 1, "OUT1")
        assert extout.type == "EXTOUT"
        assert 72 in zaero.extout

        actu = zaero.add_actu(100, 0.8, 0.5, 0.2)
        assert actu.type == "ACTU"
        assert 100 in zaero.actu

        loadmod = zaero.add_loadmod(73, "WING", 0, 10, 20)
        assert loadmod.type == "LOADMOD"
        assert 73 in zaero.loadmod

        rbred = zaero.add_rbred(74, 1, "123", 100, "NO")
        assert rbred.type == "RBRED"
        assert 74 in zaero.rbred

        # --- gust ---
        gloads = zaero.add_gloads(
            80, 1, 1, 0, 0, 0, 0, 0, "SAVE", "UNFORM", "gloads.dat", "SAVE", "gloads_freq.dat"
        )
        assert gloads.type == "GLOADS"
        assert 80 in zaero.gloads

        dgust = zaero.add_dgust(81, "1MC", 100.0, 50.0, 0.0)
        assert dgust.type == "DGUST"
        assert 81 in zaero.dgust

        cgust = zaero.add_cgust(82, "VONK", 0.01, 0.0, 10.0)
        assert cgust.type == "CGUST"
        assert 82 in zaero.cgust

        # --- ASE ---
        ase = zaero.add_ase(90, 1, 1, 0, 0, 0)
        assert ase.type == "ASE"
        assert 90 in zaero.ase

        asecont = zaero.add_asecont(91, 0, 0, 0, 0, 0, 0, 0)
        assert asecont.type == "ASECONT"
        assert 91 in zaero.asecont

        asesnsr = zaero.add_asesnsr(92, 1, 100, 1, 1.0, "NO")
        assert asesnsr.type == "ASESNSR"
        assert 92 in zaero.asesnsr

        asesns1 = zaero.add_asesns1(93, "ACCEL", 1, 1.0)
        assert asesns1.type == "ASESNS1"
        assert 93 in zaero.asesns1

        asegain = zaero.add_asegain(94, 92, 1, 95, 1, 2.0, "Q")
        assert asegain.type == "ASEGAIN"
        assert 94 in zaero.asegain

        gainset = zaero.add_gainset(95, [94])
        assert gainset.type == "GAINSET"
        assert 95 in zaero.gainset

        cjunct = zaero.add_cjunct(96, 2, 1, [1.0, 1.0])
        assert cjunct.type == "CJUNCT"
        assert 96 in zaero.cjunct

        conct = zaero.add_conct(97, 92, 1, 95, 1)
        assert conct.type == "CONCT"
        assert 97 in zaero.conct

        tfset = zaero.add_tfset(98, [95])
        assert tfset.type == "TFSET"
        assert 98 in zaero.tfset

        mimoss = zaero.add_mimoss(99, 2, 1, 1, "DMI1", "S", 0, [1.0, 2.0], ["IN", "OUT"])
        assert mimoss.type == "MIMOSS"
        assert 99 in zaero.mimoss

        mimotf = zaero.add_mimotf(101, 1, 1, [95])
        assert mimotf.type == "MIMOTF"
        assert 101 in zaero.mimotf

        sisotf = zaero.add_sisotf(95, 1, 2, [1.0, 2.0], [1.0, 0.5])
        assert sisotf.type == "SISOTF"
        assert 95 in zaero.sisotf

        aerolag = zaero.add_aerolag(102, 2, [0.1, 0.2])
        assert aerolag.type == "AEROLAG"
        assert 102 in zaero.aerolag

        # --- classic transfer function blocks ---
        sensr = zaero.add_sensr(110, 1, 100, 1, 1.0, "NO")
        assert sensr.type == "SENSR"
        assert 110 in zaero.sensr

        gain = zaero.add_gain(111, 2.5)
        assert gain.type == "GAIN"
        assert 111 in zaero.gain

        sumblk = zaero.add_sumblk(112, 2, [1.0, -1.0])
        assert sumblk.type == "SUMBLK"
        assert 112 in zaero.sumblk

        deadbn = zaero.add_deadbn(113, 0.01)
        assert deadbn.type == "DEADBN"
        assert 113 in zaero.deadbn

        delay = zaero.add_delay_zaero(114, 0.05, 3)
        assert delay.type == "DELAY"
        assert 114 in zaero.delay_zaero

        filtfl = zaero.add_filtfl(115, 1, 10.0, 2, 0.7)
        assert filtfl.type == "FILTFL"
        assert 115 in zaero.filtfl

        limtr = zaero.add_limtr(116, -25.0, 25.0)
        assert limtr.type == "LIMTR"
        assert 116 in zaero.limtr

        # --- sets ---
        senset = zaero.add_senset(120, [92])
        assert senset.type == "SENSET"
        assert 120 in zaero.senset

        surfset = zaero.add_surfset(121, [10])
        assert surfset.type == "SURFSET"
        assert 121 in zaero.surfset

        cnctset = zaero.add_cnctset(122, [97])
        assert cnctset.type == "CNCTSET"
        assert 122 in zaero.cnctset

        setadd = zaero.add_setadd(123, [120, 121])
        assert setadd.type == "SETADD"
        assert 123 in zaero.setadd

        # --- MLD / transient ---
        mldprnt = zaero.add_mldprnt(130, "prnt.dat", "TABLE", 0.0, 1.0, ["DISP"], [1])
        assert 130 in zaero.mldprnt

        mldstat = zaero.add_mldstat(131, 0, "YES", "stat.dat", ["ALPHA"], [5.0])
        assert mldstat.type == "MLDSTAT"
        assert 131 in zaero.mldstat

        minstat = zaero.add_minstat(132, 102, 10, 0, 0, 0, 0, 0, 0, "SAVE", "min.dat", 0)
        assert minstat.type == "MINSTAT"
        assert 132 in zaero.minstat

        mldtrim = zaero.add_mldtrim(133, 386.4, 1.0, "NO", ["ALPHA"], [5.0])
        assert mldtrim.type == "MLDTRIM"
        assert 133 in zaero.mldtrim

        mldcomd = zaero.add_mldcomd(134, [71], [0])
        assert mldcomd.type == "MLDCOMD"
        assert 134 in zaero.mldcomd

        mldtime = zaero.add_mldtime(135, 0.0, 1.0, 0.01, 0, 0, "RK4")
        assert mldtime.type == "MLDTIME"
        assert 135 in zaero.mldtime

        # --- other ---
        extfile = zaero.add_extfile(140, "ext.dat")
        assert extfile.type == "EXTFILE"
        assert 140 in zaero.extfile

        dmil = zaero.add_dmil("DMI1", 1, 1, [1.0, 2.0])
        assert dmil.type == "DMIL"
        assert ("DMI1", 1, 1) in zaero.dmil

        conmlst = zaero.add_conmlst(141, [1.0, 2.0], [1, 2])
        assert conmlst.type == "CONMLST"
        assert 141 in zaero.conmlst

        cpfact = zaero.add_cpfact(142, 1, "SYM", "X", "MULT", "WING", 1.0, 0.0, [1, 2])
        assert cpfact.type == "CPFACT"
        assert 142 in zaero.cpfact

        apconst = zaero.add_apconst(143, 0, 0, 0, 1, 1, [1.0], [1.0])
        assert apconst.type == "APCONST"
        assert 143 in zaero.apconst

        spline0 = zaero.add_spline0(144, "WING", 0, 10)
        assert spline0.type == "SPLINE0"
        assert 144 in zaero.spline0

        pbody7 = zaero.add_pbody7(145, 1, 0, [], [],)
        assert pbody7.type == "PBODY7"
        assert 145 in zaero.pbody7

        trimflt = zaero.add_trimflt(146, "CASE1", 5.0, [1.0])
        assert trimflt.type == "TRIMFLT"
        assert 146 in zaero.trimflt

        cmargin = zaero.add_cmargin(147, 6.0, -20.0, 45.0, -45.0)
        assert cmargin.type == "CMARGIN"
        assert 147 in zaero.cmargin

    def test_zaero_add_methods_xref(self):
        """Test cross_reference/uncross_reference/safe_cross_reference
        on cards created via add_* convenience methods.
        """
        model = BDF(debug=False)
        model.sol = 145
        model.nastran_format = "zona"
        zaero = ZAERO(model)
        model.zaero = zaero

        # --- leaf cards (no xref deps) ---
        asesnsr_201 = zaero.add_asesnsr(201, 1, 1001, 1, 1.0, "NO")
        asesnsr_202 = zaero.add_asesnsr(202, 2, 1002, 3, 1.0, "NO")
        sisotf_301 = zaero.add_sisotf(301, 1, 2, [1.0, 2.0], [1.0, 0.5])
        cjunct_401 = zaero.add_cjunct(401, 2, 1, [1.0, 1.0])
        actu_801 = zaero.add_actu(801, 0.8, 0.5, 0.2)

        # --- ASEGAIN: output=ASESNSR(201), input=SISOTF(301) ---
        asegain_501 = zaero.add_asegain(501, 201, 1, 301, 1, 2.5, "Q")

        # --- CONCT: various topology connections ---
        conct_601 = zaero.add_conct(601, 202, 1, 301, 1)
        conct_602 = zaero.add_conct(602, 301, 1, 401, 1)
        conct_603 = zaero.add_conct(603, 401, 1, 801, 1)

        # --- set cards ---
        senset_20 = zaero.add_senset(20, [201, 202])
        tfset_30 = zaero.add_tfset(30, [301, 401])
        gainset_40 = zaero.add_gainset(40, [501])
        cnctset_50 = zaero.add_cnctset(50, [601, 602, 603])

        # --- aero geometry (needed by AESURFZ xref) ---
        aeroz = AEROZ("SLIN", "IN", 1.0, 1.0, 1.0)
        model._add_methods.add_aeros_object(aeroz)
        caero7 = CAERO7(
            101,
            "",
            np.array([0.0, 0.0, 0.0]),
            1.0,
            np.array([0.0, 1.0, 0.0]),
            1.0,
            nspan=2,
            nchord=2,
        )
        model._add_methods.add_caero_object(caero7)
        zaero.add_panlst1(209, 101, 101, 110)

        # --- AESURFZ + AESLINK ---
        zaero.add_aesurfz("ELEV", "ASYM", 0, 209, 0, 0)
        zaero.add_aesurfz("FLAP", "ASYM", 0, 209, 0, 0)
        aeslink = zaero.add_aeslink("CTRL", "SYM", 801, ["ELEV", "FLAP"], [1.0, 0.5])

        # --- MLDCOMD -> EXTINP -> ACTU (no ASECONT) ---
        extinp_901 = zaero.add_extinp(901, "", 801, 1, "CMD1")
        mldcomd_100 = zaero.add_mldcomd(100, [901], [1001])
        model.add_tabled1(1001, [0.0, 1.0], [0.0, 1.0])

        # --- cross_reference ---
        zaero.cross_reference()

        # verify ASEGAIN refs
        assert asegain_501.output_ref is asesnsr_201
        assert asegain_501.input_ref is sisotf_301

        # verify CONCT refs
        assert conct_601.output_ref is asesnsr_202
        assert conct_601.input_ref is sisotf_301
        assert conct_602.output_ref is sisotf_301
        assert conct_602.input_ref is cjunct_401
        assert conct_603.output_ref is cjunct_401
        assert conct_603.input_ref is actu_801

        # verify SET refs
        assert senset_20.ids_ref == [asesnsr_201, asesnsr_202]
        assert tfset_30.ids_ref == [sisotf_301, cjunct_401]
        assert gainset_40.ids_ref == [asegain_501]
        assert cnctset_50.ids_ref == [conct_601, conct_602, conct_603]

        # verify AESLINK refs
        assert aeslink.actu_ref is actu_801
        assert len(aeslink.independent_labels_ref) == 2

        # verify MLDCOMD -> EXTINP refs
        assert len(mldcomd_100.extinps_ref) == 1
        assert mldcomd_100.extinps_ref[0] is extinp_901
        assert len(mldcomd_100.tables_ref) == 1
        assert extinp_901.itf_ref is actu_801

        # --- uncross_reference ---
        zaero.uncross_reference()

        # AESLINK has proper uncross_reference
        assert aeslink.actu_ref is None
        assert aeslink.independent_labels_ref == []

        # MLDCOMD has proper uncross_reference
        assert mldcomd_100.extinps_ref == []
        assert mldcomd_100.tables_ref == []

        # EXTINP has proper uncross_reference
        assert extinp_901.itf_ref is None

        # --- safe_cross_reference ---
        zaero.safe_cross_reference()

        # verify refs are repopulated
        assert asegain_501.output_ref is asesnsr_201
        assert asegain_501.input_ref is sisotf_301
        assert conct_601.output_ref is asesnsr_202
        assert conct_602.input_ref is cjunct_401
        assert senset_20.ids_ref == [asesnsr_201, asesnsr_202]
        assert tfset_30.ids_ref == [sisotf_301, cjunct_401]
        assert gainset_40.ids_ref == [asegain_501]
        assert cnctset_50.ids_ref == [conct_601, conct_602, conct_603]
        assert aeslink.actu_ref is actu_801
        assert len(aeslink.independent_labels_ref) == 2
        assert mldcomd_100.extinps_ref[0] is extinp_901
        assert extinp_901.itf_ref is actu_801

        model.cross_reference()
        # save_load_deck(
        #     model, nastran_format='zaero',
        #     run_renumber=False, run_convert=False, run_remove_unused=False,
        #     run_save_load=False, run_save_load_hdf5=False, run_mass_properties=False,
        #     run_test_bdf=False, run_op2_writer=False, run_export_caero=False,
        #     stringify=True)

    def test_zaero_1(self):
        """zaero explicit test"""
        log = SimpleLogger(level="error", encoding="utf-8")
        # log = SimpleLogger(level='debug', encoding='utf-8')
        bdf_filename = ZAERO_PATH / "f16_ma41.bdf"
        model = read_bdf(bdf_filename, xref=False, mode="zona", debug=None, log=log)
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

    def test_zaero_2(self):
        """zaero explicit test"""
        log = SimpleLogger(level="error", encoding="utf-8")  # lots of zaero errors
        bdf_filename = ZAERO_PATH / "ztran.bdf"
        model = read_bdf(bdf_filename, xref=False, debug=None, log=log)
        model.safe_cross_reference()
        model.zaero.view_block_diagram()
        save_load_deck(
            model,
            xref="safe",
            run_renumber=False,
            run_convert=False,
            run_remove_unused=False,
            run_save_load=False,
            run_save_load_hdf5=False,
            run_mass_properties=False,
            run_export_caero=False,
            run_test_bdf=False,
            run_op2_writer=False,
        )
        model.zaero.convert_to_nastran()
        write_raw_fields(model.zaero)
        model.zaero.uncross_reference()

    def test_zaero_model_1(self):
        """totally fake zaero model"""
        bdf_file = get_zaero_model()

        model = read_bdf(
            bdf_filename=bdf_file,
            validate=True,
            xref=True,
            punch=False,
            skip_cards=None,
            read_cards=None,
            encoding=None,
            log=None,
            debug=False,
            mode="zaero",
        )
        # with self.assertRaises(AttributeError):

        model.zaero.view_block_diagram()
        model.uncross_reference()
        model.write_bdf("zaero.bdf")
        model.safe_cross_reference()
        model.write_bdf("zaero.bdf")
        model.zaero.uncross_reference()
        # model.zaero.convert_to_nastran('zaero.bdf')

        bdf_file.seek(0)
        model.clear_attributes()
        model2 = read_bdf("zaero.bdf", debug=None)
        os.remove("zaero.bdf")
        write_raw_fields(model2.zaero)
        model2.zaero.convert_to_nastran()
        # model2.write_bdf('zaero2.bdf')

    def test_zaero_trim_case1_in(self):
        zaero_filename = TRIM_DIR / "case1" / "ha144d.inp"
        model = read_bdf(zaero_filename, xref=False, debug=False, mode="zaero")
        model.zaero.safe_cross_reference()
        model.zaero.view_block_diagram()
        model.zaero.uncross_reference()

        with self.assertRaises((AssertionError, RuntimeError)):
            model.cross_reference()
        write_raw_fields(model.zaero)
        model.write_bdf(TRIM_DIR / "zaero.inp")

    def test_zaero_flutter_case1_in(self):
        zaero_filename = FLUTTER_DIR / "case1" / "ha145e.inp"
        model = read_bdf(zaero_filename, xref=True, debug=False)
        model.zaero.safe_cross_reference()
        model.zaero.view_block_diagram()
        write_raw_fields(model.zaero)
        model.write_bdf(FLUTTER_DIR / "zaero.inp")
        model.zaero.uncross_reference()

    def test_zaero_flutter_case2_in(self):
        zaero_filename = FLUTTER_DIR / "case2" / "crop.inp"
        model = read_bdf(zaero_filename, xref=True, debug=False)
        model.zaero.safe_cross_reference()
        model.zaero.view_block_diagram()
        write_raw_fields(model.zaero)
        model.write_bdf(FLUTTER_DIR / "zaero.inp")
        model.zaero.uncross_reference()

    def test_zaero_flutter_case3_in(self):
        zaero_filename = FLUTTER_DIR / "case3" / "ha145fb.inp"
        model = read_bdf(zaero_filename, xref=True, debug=False)
        model.zaero.safe_cross_reference()
        model.zaero.view_block_diagram()
        write_raw_fields(model.zaero)
        model.write_bdf(FLUTTER_DIR / "zaero.inp")
        model.zaero.uncross_reference()

    def test_zaero_flutter_case4_in(self):
        zaero_filename = FLUTTER_DIR / "case4" / "ha145g.inp"
        model = read_bdf(zaero_filename, xref=True, debug=False)
        model.zaero.safe_cross_reference()
        model.zaero.view_block_diagram()
        write_raw_fields(model.zaero)
        model.write_bdf(FLUTTER_DIR / "zaero.inp")
        model.zaero.uncross_reference()

    def test_zaero_flutter_case5_in(self):
        zaero_filename = FLUTTER_DIR / "case5" / "f16ma41.inp"
        model = read_bdf(zaero_filename, xref=False, debug=False)
        model.zaero.cross_reference()
        model.zaero.view_block_diagram()
        model.zaero.safe_cross_reference()
        write_raw_fields(model.zaero)
        model.write_bdf(FLUTTER_DIR / "zaero.inp")
        model.zaero.uncross_reference()

    def test_zaero_flutter_case6_in_trim(self):
        zaero_filename = FLUTTER_DIR / "case6" / "agard_trim.inp"
        model = read_bdf(zaero_filename, xref=False, debug=False)
        model.zaero.cross_reference()
        model.zaero.view_block_diagram()
        model.zaero.safe_cross_reference()
        write_raw_fields(model.zaero)
        model.write_bdf(FLUTTER_DIR / "zaero.inp")
        model.zaero.uncross_reference()

    def test_zaero_flutter_case6_in_tran(self):
        zaero_filename = FLUTTER_DIR / "case6" / "agardztran.inp"
        model = read_bdf(zaero_filename, xref=False, debug=False)
        model.zaero.cross_reference()
        model.zaero.view_block_diagram()
        model.zaero.safe_cross_reference()
        write_raw_fields(model.zaero)
        model.write_bdf(MLOADS_DIR / "zaero.inp")
        model.zaero.uncross_reference()

    def test_zaero_flutter_case7_in(self):
        zaero_filename = FLUTTER_DIR / "case7" / "agardztaw.inp"
        model = read_bdf(zaero_filename, xref=False, debug=False)
        model.zaero.cross_reference()
        model.zaero.view_block_diagram()
        model.zaero.safe_cross_reference()
        write_raw_fields(model.zaero)
        model.write_bdf(MLOADS_DIR / "zaero.inp")
        model.zaero.uncross_reference()

    def test_zaero_mloads_case1_in(self):
        zaero_filename = MLOADS_DIR / "case1" / "m144open.inp"
        model = read_bdf(zaero_filename, xref=True, debug=False)
        model.zaero.view_block_diagram()
        model.zaero.uncross_reference()
        # model.zaero.cross_reference()
        model.zaero.safe_cross_reference()
        write_raw_fields(model.zaero)
        model.write_bdf(MLOADS_DIR / "zaero.inp")
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
        zaero_filename = MLOADS_DIR / "case2" / "m144_trim.inp"
        model = read_bdf(zaero_filename, xref=True, debug=False)
        model.zaero.view_block_diagram()
        model.zaero.safe_cross_reference()
        write_raw_fields(model.zaero)
        model.write_bdf(MLOADS_DIR / "zaero.inp")
        model.zaero.uncross_reference()

    def test_zaero_mloads_case2b_in(self):
        zaero_filename = MLOADS_DIR / "case2" / "m144clos.inp"
        model = read_bdf(zaero_filename, xref=True, debug=False)
        model.zaero.view_block_diagram()
        model.zaero.safe_cross_reference()
        write_raw_fields(model.zaero)
        model.write_bdf(MLOADS_DIR / "zaero.inp")
        model.zaero.uncross_reference()

    def test_zaero_ase_case1(self):
        zaero_filename = ASE_DIR / "case1" / "cropase.inp"
        model = read_bdf(zaero_filename, mode="zaero", xref=True, debug=False)
        model.zaero.view_block_diagram()
        model.zaero.safe_cross_reference()
        write_raw_fields(model.zaero)
        model.write_bdf(ASE_DIR / "zaero.inp")
        model.zaero.uncross_reference()

    def test_zaero_ase_case2(self):
        zaero_filename = ASE_DIR / "case2" / "gafa.inp"
        model = read_bdf(zaero_filename, mode="zaero", xref=True, debug=False)
        model.zaero.view_block_diagram()
        model.zaero.safe_cross_reference()
        write_raw_fields(model.zaero)
        model.write_bdf(ASE_DIR / "zaero.inp")
        model.zaero.uncross_reference()

    def test_zaero_gloads_case1(self):
        zaero_filename = GLOADS_DIR / "case1" / "kussner.inp"
        model = read_bdf(zaero_filename, mode="zaero", xref=False, debug=False)
        model.zaero.cross_reference()
        model.zaero.view_block_diagram()
        model.zaero.safe_cross_reference()
        write_raw_fields(model.zaero)
        model.write_bdf(GLOADS_DIR / "zaero.inp")
        model.zaero.uncross_reference()

    def test_zaero_gloads_case2(self):
        zaero_filename = GLOADS_DIR / "case2" / "gbj_dgust.inp"
        model = read_bdf(zaero_filename, mode="zaero", xref=False, debug=False)
        model.zaero.cross_reference()
        model.zaero.view_block_diagram()
        model.zaero.safe_cross_reference()
        write_raw_fields(model.zaero)
        model.write_bdf(GLOADS_DIR / "zaero.inp")
        model.zaero.uncross_reference()

    def test_zaero_gloads_case3(self):
        zaero_filename = GLOADS_DIR / "case3" / "gbj_cgust.inp"
        model = read_bdf(zaero_filename, mode="zaero", xref=False, debug=False)
        model.zaero.cross_reference()
        model.zaero.view_block_diagram()
        model.zaero.safe_cross_reference()
        write_raw_fields(model.zaero)
        model.write_bdf(GLOADS_DIR / "zaero.inp")
        model.zaero.uncross_reference()

    def test_zaero_gloads_case4a(self):
        zaero_filename = GLOADS_DIR / "case4" / "cgust_md.inp"
        model = read_bdf(zaero_filename, mode="zaero", xref=False, debug=False)
        model.zaero.cross_reference()
        model.zaero.view_block_diagram()
        model.write_bdf(GLOADS_DIR / "zaero.inp")
        model.zaero.uncross_reference()

    def test_zaero_gloads_case4b(self):
        zaero_filename = GLOADS_DIR / "case4" / "cgust_sof.inp"
        model = read_bdf(zaero_filename, mode="zaero", xref=False, debug=False)
        model.zaero.cross_reference()
        model.zaero.view_block_diagram()
        model.write_bdf(GLOADS_DIR / "zaero.inp")
        model.zaero.uncross_reference()

    def test_zaero_new_ase_cards(self):
        """Test parsing and round-tripping MIMOTF, SENSR, GAIN, SUMBLK,
        DEADBN, DELAY, FILTFL, LIMTR cards.
        """
        # MIMOTF: 2x2 MIMO transfer function matrix referencing 4 SISOTFs
        card = BDFCard(["MIMOTF", "10", "2", "2", "101", "102", "103", "104"])
        obj = MIMOTF.add_card(card)
        assert obj.mimotf_id == 10
        assert obj.n_input == 2
        assert obj.n_output == 2
        assert obj.tf_ids == [101, 102, 103, 104]
        fields = obj.raw_fields()
        assert fields == ["MIMOTF", 10, 2, 2, 101, 102, 103, 104]
        obj.write_card()

        # SENSR: acceleration sensor at grid 100, DOF 3
        card = BDFCard(["SENSR", "20", "2", "100", "3", "0.5", "YES"])
        obj = SENSR.add_card(card)
        assert obj.sensr_id == 20
        assert obj.sensor_type == 2
        assert obj.sgid == 100
        assert obj.component == 3
        assert obj.factor == 0.5
        assert obj.sum_method == "YES"
        obj.write_card()

        # SENSR: defaults (factor=1.0, sum_method=NO)
        card = BDFCard(["SENSR", "21", "0", "200", "1"])
        obj = SENSR.add_card(card)
        assert obj.factor == 1.0
        assert obj.sum_method == "NO"

        # GAIN: scalar gain
        card = BDFCard(["GAIN", "30", "2.5"])
        obj = GAIN.add_card(card)
        assert obj.gain_id == 30
        assert obj.k == 2.5
        fields = obj.raw_fields()
        assert fields == ["GAIN", 30, 2.5]
        obj.write_card()

        # SUMBLK: 3-input summing junction
        card = BDFCard(["SUMBLK", "40", "3", "1.0", "-1.0", "1.0"])
        obj = SUMBLK.add_card(card)
        assert obj.sumblk_id == 40
        assert obj.nsignal == 3
        assert obj.signs == [1.0, -1.0, 1.0]
        obj.write_card()

        # DEADBN: deadband with threshold 0.25
        card = BDFCard(["DEADBN", "50", "0.25"])
        obj = DEADBN.add_card(card)
        assert obj.deadbn_id == 50
        assert obj.threshold == 0.25
        obj.write_card()

        # DELAY: 20ms delay, order 4 Pade approximation
        card = BDFCard(["DELAY", "60", "0.02", "4"])
        obj = DELAY_ZAERO.add_card(card)
        assert obj.delay_id == 60
        assert obj.tau == 0.02
        assert obj.order == 4
        obj.write_card()

        # DELAY: default order=3
        card = BDFCard(["DELAY", "61", "0.01"])
        obj = DELAY_ZAERO.add_card(card)
        assert obj.order == 3

        # FILTFL: low-pass filter, 15 Hz, order 3, zeta=0.5
        card = BDFCard(["FILTFL", "70", "1", "15.0", "3", "0.5"])
        obj = FILTFL.add_card(card)
        assert obj.filtfl_id == 70
        assert obj.filter_type == 1
        assert obj.freq == 15.0
        assert obj.order == 3
        assert obj.zeta == 0.5
        obj.write_card()

        # FILTFL: defaults (order=2, zeta=0.707)
        card = BDFCard(["FILTFL", "71", "2", "20.0"])
        obj = FILTFL.add_card(card)
        assert obj.order == 2
        assert abs(obj.zeta - 0.707) < 1e-10

        # LIMTR: saturation limits
        card = BDFCard(["LIMTR", "80", "-30.0", "30.0"])
        obj = LIMTR.add_card(card)
        assert obj.limtr_id == 80
        assert obj.lower == -30.0
        assert obj.upper == 30.0
        obj.write_card()

    def test_zaero_new_ase_cards_bdf_read(self):
        """Test reading new ASE cards from a BDF file and cross-referencing."""
        bdf_file = StringIO(
            "$ pyNastran: version=zona\n"
            "CEND\n"
            "BEGIN BULK\n"
            "SISOTF  101     2       1       1.0     0.5     1.0     0.0\n"
            "SISOTF  102     2       1       1.0     0.5     1.0     0.0\n"
            "SISOTF  103     2       1       1.0     0.5     1.0     0.0\n"
            "SISOTF  104     2       1       1.0     0.5     1.0     0.0\n"
            "MIMOTF  10      2       2       101     102     103     104\n"
            "SENSR   20      2       100     3       0.5     YES\n"
            "GAIN    30      2.5\n"
            "SUMBLK  40      3       1.0     -1.0    1.0\n"
            "DEADBN  50      0.25\n"
            "DELAY   60      0.02    4\n"
            "FILTFL  70      1       15.0    3       0.5\n"
            "LIMTR   80      -30.0   30.0\n"
            "ENDDATA\n"
        )
        model = read_bdf(bdf_file, debug=False, mode="zaero")

        # All cards parsed
        assert 10 in model.zaero.mimotf
        assert 20 in model.zaero.sensr
        assert 30 in model.zaero.gain
        assert 40 in model.zaero.sumblk
        assert 50 in model.zaero.deadbn
        assert 60 in model.zaero.delay_zaero
        assert 70 in model.zaero.filtfl
        assert 80 in model.zaero.limtr

        # MIMOTF cross-references 4 SISOTFs
        mimotf = model.zaero.mimotf[10]
        assert hasattr(mimotf, "tf_ids_ref")
        assert len(mimotf.tf_ids_ref) == 4
        assert all(ref is not None for ref in mimotf.tf_ids_ref)

        # Write and verify round-trip
        write_raw_fields(model.zaero)


def get_zaero_model() -> StringIO:
    bdf_file = StringIO()
    bdf_file.write(
        "$ pyNastran: version=zona\n"
        "CEND\n"
        "BEGIN BULK\n"
        # ATMOS IDATM AMMUNIT AMLUNIT AMTUNIT
        #       ALT1  SOUND1  DEN1    TEMP1  ALT2 SOUND2 DEN2 TEMP2
        #       ALTi  SOUNDi  DENi    TEMPi -etc-
        "ATMOS, 12,   SLIN,   IN,     F\n"
        ",      0.0,  1000.,  1.0,    2.0,  10000.0, 900.,  1.0,    2.0\n"
        # FIXMATM SETID   IDMK  IDATM FTMUNIT FTLUNIT VREF FLUTTF PRINT CONT
        #         ALT1    ALT2  ALTi -etc-
        "FIXMATM, 100,    101,   12,  slug,   ft,     1.0, -1,    0\n"
        ",       -10000., 0., 10000., 20000., 30000.\n"
        #'$       acsid, rcsid, cref, bref, sref, symxz, symxy\n'
        #'AEROZ, 10,     0,     1.,   10.,  100., YES\n'
        "$AEROZ  ACSID XZSYM FLIP FMMUNIT FMLUNIT REFC   REFB   REFS\n"
        "$       REFX  REFY  REFZ\n"
        "AEROZ,  0,    YES,  NO,  SLIN,   IN,      22.73,59.394,1175.8\n"
        ",       59.53,0.0,  0.0\n"
        "$       label, type, cid, PANLST, setg, actid\n"
        "AESURFZ,FLAP,  ASYM, 1,   10,       20,   0\n"
        #'AESURFZ,FLAP,  SYM,  1,  10,       20,   0\n'
        "CORD2R, 1,0, 0.,0.,0., 0.,0.,1.,\n"
        ",1.,0.,0.\n"
        "$BODY7,ID,LABEL,IPBODY7, ACOORD, NSEG, IDMESH1\n"
        "BODY7, 1, FUSE,        ,      2,     , 1\n"
        "PANLST3,10, FUSE, \n"
        "$       id,naxial,nradial, \n"
        "SEGMESH,1, 4,     3,       \n"
        # ITYPEi = 1 (Body of Revolution):
        #    Xi, CAMi, YRi
        # ITYPEi = 2 (Elliptical Body):
        #    Xi, YRi, ZRi
        # ITYPEi = 3 (Arbitrary Body):
        #    Xi, IDYi, IDZi
        "$       itype, x1, cam, yr1, zr1, idy1, idz1 \n"
        ",       1,        ,  1.,  1.,    ,\n"
        ",       2,      1.,    ,  1.,  2.,\n"
        ",       3,      2.,    ,    ,    , 13,   14   \n"
        ",       3,      3.,    ,    ,    , 13,   14   \n"
        # y
        "AEFACT,13, 1., 0.,  0.,-1.\n"
        "AEFACT,14, 0., 1., -1., 0.\n"
        "$ MKAEROZ, ID, MACH, METHOD, IDFLT\n"
        "MKAEROZ,   101, 0.8, -1,     -1,  \n"
        "$ TRIM, ID, MKAEROZ, Q,   IDOBJ, IDCONS, CGX, CGY, CGZ\n"
        "TRIM, 100,  101,     42., 0,     0,      5., 0., 0.\n"
        "$WTMASS, WEIGHT, Ixx, Ixy, Iyy, Ixz, Iyz, Izz\n"
        ",1.,     1e4,    1e3, 1e3, 1e5, 1e3, 1e3, 1e4\n"
        "$TRUE/G, NX,     NY,  NZ,  PDOT, QDOT, RDOT,\n"
        ", TRUE,  FREE, NONE,  32., FREE, NONE, 42.,\n"
        "$var, value\n"
        ",17,  1.0,\n"
        "$\n"
        "TRIMVAR,17,VAR\n"
        "$\n"
        "$trimlnk,id,sym, ?,  ?\n"
        "TRIMLNK,10,SYM, -1., 17\n"
        "ACOORD, 2, 0.,0.,0., 1.0,0.\n"
        "$       ID,    MODEL, CP, PANLST, SETG, DZ, EPS\n"
        "SPLINE1,100,        ,   ,    422, 423,\n"
        "$,      NELEM, MELEM\n"
        "$,      10,    10\n"
        "PANLST3,422, FUSE, \n"
        "$       id,naxial,nradial, \n"
        #'$       ID,   MODEL, PANLST, SETG,\n'
        #'SPLINE2,1000,      ,    422,  423,\n'
        "$       ID,   MODEL, CP, PANLST, SETG,\n"
        "SPLINE3,1200,      ,   ,    422,  423,\n"
        "SET1,423,10\n"
        "GRID,10,,0.,0.,0.\n"
        "GRID,11,,1.,0.,0.\n"
        "CONROD,100, 10,11, 101,1.0\n"
        "MAT1,101,3.0e7,,0.3\n"
    )
    bdf_file.seek(0)
    return bdf_file


def write_raw_fields(zaero: ZAERO):
    singletons, dicts, dicts_list = get_dicts(zaero, "write")
    for panlsts in zaero.panlsts.values():
        for panlst in panlsts:
            panlst.raw_fields()

    for dicti in dicts:
        for value in dicti.values():
            value.raw_fields()


class TestNewZaeroCards(unittest.TestCase):
    """Tests for newly added ZAERO cards: CONMLST, CPFACT, APCONST,
    SPLINE0, PBODY7, TRIMFLT, and generic cards (ASEOUT, OUTPUT4, etc.).
    """

    def test_conmlst(self):
        """CONMLST parses scale factor + CONM ID pairs."""
        bdf_file = StringIO(
            "$ pyNastran: version=zona\n"
            "CEND\n"
            "BEGIN BULK\n"
            "CONMLST 3                                                               +C1\n"
            "+C1     2.0     100\n"
            "ENDDATA\n"
        )
        model = read_bdf(bdf_file, debug=False, mode="zaero")
        assert 3 in model.zaero.conmlst
        card = model.zaero.conmlst[3]
        assert card.factors == [2.0]
        assert card.conm_ids == [100]
        card.raw_fields()

    def test_conmlst_multiple_pairs(self):
        """CONMLST with multiple (factor, CONM) pairs."""
        bdf_file = StringIO(
            "$ pyNastran: version=zona\n"
            "CEND\n"
            "BEGIN BULK\n"
            "CONMLST 5                                                               +C\n"
            "+C      3.5     200     1.2     300\n"
            "ENDDATA\n"
        )
        model = read_bdf(bdf_file, debug=False, mode="zaero")
        card = model.zaero.conmlst[5]
        assert card.factors == [3.5, 1.2]
        assert card.conm_ids == [200, 300]

    def test_cpfact(self):
        """CPFACT parses pressure scaling factor fields."""
        bdf_file = StringIO(
            "$ pyNastran: version=zona\n"
            "CEND\n"
            "BEGIN BULK\n"
            "CPFACT  90      90      SYM     ALL     FEM     ALL     0.0     0.0     +CP\n"
            "+CP     1\n"
            "ENDDATA\n"
        )
        model = read_bdf(bdf_file, debug=False, mode="zaero")
        assert 90 in model.zaero.cpfact
        card = model.zaero.cpfact[90]
        assert card.idmk == 90
        assert card.sym == "SYM"
        assert card.comp == "ALL"
        assert card.cptype == "FEM"
        assert card.factor1 == 0.0
        assert card.strips == [1]
        card.raw_fields()

    def test_apconst(self):
        """APCONST parses RFA constraint parameters."""
        bdf_file = StringIO(
            "$ pyNastran: version=zona\n"
            "CEND\n"
            "BEGIN BULK\n"
            "APCONST 30      1       -1      -1      1       1       1       1\n"
            "ENDDATA\n"
        )
        model = read_bdf(bdf_file, debug=False, mode="zaero")
        assert 30 in model.zaero.apconst
        card = model.zaero.apconst[30]
        assert card.sid == 30
        assert card.da0 == 1
        assert card.da1 == -1
        assert card.da2 == -1
        assert card.nrp == 1
        assert card.ncp == 1
        card.raw_fields()

    def test_spline0(self):
        """SPLINE0 parses rigid body spline fields."""
        bdf_file = StringIO(
            "$ pyNastran: version=zona\n"
            "CEND\n"
            "BEGIN BULK\n"
            "SPLINE0 102     FUSELAGE        112\n"
            "ENDDATA\n"
        )
        model = read_bdf(bdf_file, debug=False, mode="zaero")
        assert 102 in model.zaero.spline0
        card = model.zaero.spline0[102]
        assert card.eid == 102
        assert card.model_name == "FUSELAGE"
        assert card.setk == 112
        card.raw_fields()

    def test_pbody7(self):
        """PBODY7 parses body property fields."""
        bdf_file = StringIO(
            "$ pyNastran: version=zona\n"
            "CEND\n"
            "BEGIN BULK\n"
            "PBODY7  11      0                                               0\n"
            "ENDDATA\n"
        )
        model = read_bdf(bdf_file, debug=False, mode="zaero")
        assert 11 in model.zaero.pbody7
        card = model.zaero.pbody7[11]
        assert card.pid == 11
        assert card.wake == 0
        card.raw_fields()

    def test_trimflt(self):
        """TRIMFLT parses trim-flutter card."""
        bdf_file = StringIO(
            "$ pyNastran: version=zona\nCEND\nBEGIN BULK\nTRIMFLT 10              2.0\nENDDATA\n"
        )
        model = read_bdf(bdf_file, debug=False, mode="zaero")
        assert 10 in model.zaero.trimflt
        card = model.zaero.trimflt[10]
        assert card.trimflt_id == 10
        assert card.alpha == 2.0
        card.raw_fields()

    def test_trimflt_integer_field2(self):
        """TRIMFLT with integer in field 2 (reference ID format)."""
        bdf_file = StringIO(
            "$ pyNastran: version=zona\nCEND\nBEGIN BULK\nTRIMFLT 1       1\nENDDATA\n"
        )
        model = read_bdf(bdf_file, debug=False, mode="zaero")
        assert 1 in model.zaero.trimflt
        card = model.zaero.trimflt[1]
        assert card.trimflt_id == 1
        assert card.title == "1"

    def test_generic_cards_aseout(self):
        """ASEOUT generic card parses without error."""
        bdf_file = StringIO(
            "$ pyNastran: version=zona\n"
            "CEND\n"
            "BEGIN BULK\n"
            "ASEOUT  101     10      01      PLANT           PLANT01.DAT\n"
            "ENDDATA\n"
        )
        model = read_bdf(bdf_file, debug=False, mode="zaero")
        assert 101 in model.zaero.aseout
        card = model.zaero.aseout[101]
        assert card.card_id == 101
        card.raw_fields()

    def test_generic_cards_output4(self):
        """OUTPUT4 generic card parses without error."""
        bdf_file = StringIO(
            "$ pyNastran: version=zona\nCEND\nBEGIN BULK\nOUTPUT4 QHGS0101QHGS0101.DAT\nENDDATA\n"
        )
        model = read_bdf(bdf_file, debug=False, mode="zaero")
        assert len(model.zaero.output4) == 1

    def test_generic_cards_cellwng(self):
        """CELLWNG generic card parses without error."""
        bdf_file = StringIO(
            "$ pyNastran: version=zona\n"
            "CEND\n"
            "BEGIN BULK\n"
            "CELLWNG 10001   1001    1       3               20001\n"
            "ENDDATA\n"
        )
        model = read_bdf(bdf_file, debug=False, mode="zaero")
        assert 10001 in model.zaero.cellwng

    def test_all_examples_zero_rejects(self):
        """All ZAERO example .inp files parse with zero rejected cards."""
        all_rejects = {}
        for subdir in ["flutter", "ase", "gloads", "trim", "mloads"]:
            edir = EXAMPLES_DIR / ".." / subdir
            if not edir.exists():
                continue
            for case_dir in sorted(edir.iterdir()):
                if not case_dir.is_dir():
                    continue
                for inp in case_dir.glob("*.inp"):
                    if "test_bdf" in inp.name:
                        continue
                    try:
                        model = read_bdf(inp, mode="zaero", xref=False, debug=False)
                        for card, count in model.reject_count.items():
                            if card not in all_rejects:
                                all_rejects[card] = 0
                            all_rejects[card] += count
                    except Exception:
                        pass
        self.assertEqual(all_rejects, {}, f"Rejected cards: {all_rejects}")

    def test_rbred(self):
        """RBRED parses rigid body mode reduction card."""
        bdf_file = StringIO(
            "$ pyNastran: version=zona\n"
            "CEND\n"
            "BEGIN BULK\n"
            "RBRED   10      200     246     10\n"
            "ENDDATA\n"
        )
        model = read_bdf(bdf_file, debug=False, mode="zaero")
        assert 10 in model.zaero.rbred
        card = model.zaero.rbred[10]
        assert card.sid == 10
        assert card.id_ase == 200
        assert card.component == "246"
        assert card.node_id == 10
        card.raw_fields()

    def test_conct(self):
        """CONCT parses fixed connection between control elements."""
        bdf_file = StringIO(
            "$ pyNastran: version=zona\n"
            "CEND\n"
            "BEGIN BULK\n"
            "CONCT   425     209     1       210     1\n"
            "ENDDATA\n"
        )
        model = read_bdf(bdf_file, debug=False, mode="zaero")
        assert 425 in model.zaero.conct
        card = model.zaero.conct[425]
        assert card.conct_id == 425
        card.raw_fields()

    def test_cjunct(self):
        """CJUNCT parses junction element (gain splitter)."""
        bdf_file = StringIO(
            "$ pyNastran: version=zona\n"
            "CEND\n"
            "BEGIN BULK\n"
            "CJUNCT  214     1       2       0.033333-0.25\n"
            "ENDDATA\n"
        )
        model = read_bdf(bdf_file, debug=False, mode="zaero")
        assert 214 in model.zaero.cjunct
        card = model.zaero.cjunct[214]
        assert card.cjunct_id == 214
        card.raw_fields()

    def test_cnctset(self):
        """CNCTSET parses connection set card."""
        bdf_file = StringIO(
            "$ pyNastran: version=zona\n"
            "CEND\n"
            "BEGIN BULK\n"
            "CNCTSET 130     421     422     423     424     425\n"
            "ENDDATA\n"
        )
        model = read_bdf(bdf_file, debug=False, mode="zaero")
        assert 130 in model.zaero.cnctset
        card = model.zaero.cnctset[130]
        assert card.cnctset_id == 130
        card.raw_fields()

    def test_mldstat(self):
        """MLDSTAT parses airframe state card."""
        bdf_file = StringIO(
            "$ pyNastran: version=zona\n"
            "CEND\n"
            "BEGIN BULK\n"
            "MLDSTAT 100     100                                                     +M1\n"
            "+M1     ALPHA   -4.31-5 Q       0.      THETA   -4.31-5 H       0.0\n"
            "ENDDATA\n"
        )
        model = read_bdf(bdf_file, debug=False, mode="zaero", xref=False)
        assert 100 in model.zaero.mldstat
        card = model.zaero.mldstat[100]
        assert card.mldstat_id == 100
        assert card.mldtrim_id == 100
        assert "ALPHA" in card.states
        assert len(card.states) == 4
        card.raw_fields()

    def test_mldtrim(self):
        """MLDTRIM parses trim equilibrium card."""
        bdf_file = StringIO(
            "$ pyNastran: version=zona\n"
            "CEND\n"
            "BEGIN BULK\n"
            "MLDTRIM 100     32.2    1.0     YES     SMODAL                          +MA\n"
            "+MA     CANARD  0.13534\n"
            "ENDDATA\n"
        )
        model = read_bdf(bdf_file, debug=False, mode="zaero")
        assert 100 in model.zaero.mldtrim
        card = model.zaero.mldtrim[100]
        assert card.mldtrim_id == 100
        assert card.gravity == 32.2
        assert card.nz == 1.0
        card.raw_fields()

    def test_conmlst_from_example(self):
        """CONMLST from Kussner example parses mass scale correctly."""
        zaero_filename = GLOADS_DIR / "case1" / "kussner.inp"
        model = read_bdf(zaero_filename, mode="zaero", xref=False, debug=False)
        # CONMLST 1: factor=7956.982, CONM1=100
        assert 1 in model.zaero.conmlst
        card = model.zaero.conmlst[1]
        self.assertAlmostEqual(card.factors[0], 7956.982, places=2)
        assert card.conm_ids[0] == 100
        # CONMLST 3: factor=2.0, CONM1=100
        card3 = model.zaero.conmlst[3]
        assert card3.factors[0] == 2.0

    def test_dmil_from_example(self):
        """DMIL parses from GAFA ASE example (MIMO controller matrices)."""
        zaero_filename = ASE_DIR / "case2" / "gafa.inp"
        model = read_bdf(zaero_filename, mode="zaero", xref=False, debug=False)
        assert len(model.zaero.dmil) > 0
        # DMIL MIMO has multiple column entries
        for key, dmil in model.zaero.dmil.items():
            dmil.raw_fields()
            break  # just check one

    def test_gengust_from_example(self):
        """GENGUST parses from GBJ discrete gust example."""
        zaero_filename = GLOADS_DIR / "case2" / "gbj_dgust.inp"
        model = read_bdf(zaero_filename, mode="zaero", xref=False, debug=False)
        assert len(model.zaero.gengust_card) > 0

    def test_mftgust_from_example(self):
        """MFTGUST parses from GBJ discrete gust example."""
        zaero_filename = GLOADS_DIR / "case2" / "gbj_dgust.inp"
        model = read_bdf(zaero_filename, mode="zaero", xref=False, debug=False)
        assert len(model.zaero.mftgust) > 0

    def test_trimflt_from_flutter_example(self):
        """TRIMFLT parses from F-16 flutter example."""
        zaero_filename = FLUTTER_DIR / "case5" / "f16ma41.inp"
        model = read_bdf(zaero_filename, mode="zaero", xref=False, debug=False)
        assert len(model.zaero.trimflt) > 0
        card = next(iter(model.zaero.trimflt.values()))
        assert card.trimflt_id == 10
        assert card.alpha == 2.0

    def test_cmargin_cross_reference(self):
        """CMARGIN is properly cross-referenced from ASE card."""
        zaero_filename = ASE_DIR / "case2" / "gafa.inp"
        model = read_bdf(zaero_filename, mode="zaero", xref=False, debug=False)
        model.zaero.safe_cross_reference()
        # ASE 30 references CMARGIN 1
        ase30 = model.zaero.ase[30]
        assert ase30.cmargin_ref is not None
        assert ase30.cmargin_ref.cmargin_id == 1
        assert ase30.cmargin_ref.gm_high == 50.0
        assert ase30.cmargin_ref.pm_high == 60.0


def get_mloads_file():
    lines = [
        "$ pyNastran: version=zona",
        "CEND",
        "SUBCASE 1",
        "  MLOADS = 10",
        "BEGIN BULK",
        #
        # --- core cards ---
        # ASECONT SID SURFID SENSID TFID GAINID CONCTID EINPID EOUTID
        "ASECONT, 4,      0,    20,   30,    40,    50,      ,    60",
        #
        # MLOADS  ID ASECONT FLUTTER MINSTAT MLDSTAT MLDCOMD MLDTIME
        "MLOADS, 10,       4,      5,      ,       , 100,    101",
        #
        'MKAEROZ,7,0.8',
        ',1.0',
        #           idmk
        'FIXMATM,6, 7, 0, slug, ft,',
        ',0.0, 10000.',
        "FLUTTER, 5, SYM, 6",
        "MLDTIME,101, 0.0, 1.0, 0.1",
        #
        # --- sets ---
        "SENSET,  20, 201, 202",
        "TFSET,   30, 301, 401",
        "GAINSET, 40, 501",
        "CNCTSET, 50, 601, 602, 603",
        "SET1,    60, 701",
        #
        # --- sensors ---
        "$ ASESNSR ID TYPE SGID SC FACTOR SOF",
        "ASESNSR, 201, 1, 1001, 1",
        "ASESNSR, 202, 2, 1002, 3",
        #
        # --- transfer functions ---
        "$ SISOTF ID NDEN NNUM den1..denN num1..numN",
        "SISOTF,  301, 2, 1, 1.0, 0.5, 1.0, 2.0",
        #
        # --- junctions ---
        "$ CJUNCT ID NU NY values",
        "CJUNCT,  401, 2, 1, 1.0, 1.0",
        #
        # --- gains (ASESNSR output -> SISOTF input) ---
        "$ ASEGAIN ID OTFID CO ITFID CI GAIN TYPE",
        "ASEGAIN, 501, 201, 1, 301, 1, 2.5",
        #
        # --- connections ---
        "$ CONCT ID OTFID CO ITFID CI",
        "$ CONCT 601: ASESNSR(output) -> SISOTF(input)",
        "CONCT,   601, 202, 1, 301, 1",
        "$ CONCT 602: SISOTF(output) -> CJUNCT(input)",
        "CONCT,   602, 301, 1, 401, 1",
        "$ CONCT 603: CJUNCT(output) -> ACTU(input)",
        "CONCT,   603, 401, 1, 801, 1",
        #
        # --- actuator ---
        "ACTU, 801, 0.8, 0.5, 0.2",
        #
        # --- AESLINK (ACTU -> AESURFZ) ---
        "$ AESLINK LABEL TYPE ACTID",
        "$         COEFF1 AESURF1 COEFF2 AESURF2",
        "AESLINK, CTRL,  SYM, 801,,,,,",
        ",        1.0,  ELEV,  0.5, FLAP",
        #
        # --- AESURFZ control surfaces ---
        "AESURFZ, ELEV, ASYM, 0, 209",
        "AESURFZ, FLAP, ASYM, 0, 209",
        #
        # --- MLDCOMD -> EXTINP -> CJUNCT ---
        "$ MLDCOMD SETID CONT (then EXTINP/TABLE pairs on cont)",
        "MLDCOMD, 100,,,,,,,,",
        ",        901, 1001",
        #
        "$ EXTINP ID TYPE ITFID CI LABEL",
        "EXTINP,  901,   , 401, 1, CMD1",
        #
        "$ TABLED1 for MLDCOMD",
        "TABLED1, 1001",
        ",        0.0, 0.0, 1.0, 1.0, ENDT",
        #
        # --- EXTOUT (via SET1=60) -> CJUNCT ---
        "$ EXTOUT ID TYPE ITFID CI LABEL",
        "EXTOUT,  701,   , 401, 1, JOUT",
        #
        # --- aero geometry ---
        "PANLST1, 209, 101, 101, 110",
        "CAERO7       101                       2       2      0",
        "              0.      0.      0.      1.",
        "              0.      1.      0.      1.",
        "AEROZ,  , NO, NO, , IN",
        "ENDDATA",
    ]
    bdf_file = StringIO("\n".join(lines))
    bdf_file.bdf_filename = "cat.bdf"
    return bdf_file


if __name__ == "__main__":  # pragma: no cover
    unittest.main()
