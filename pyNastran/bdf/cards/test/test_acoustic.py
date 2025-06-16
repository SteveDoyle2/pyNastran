# coding: utf-8
"""tests acoustic cards"""
import os
from pathlib import Path
import unittest

from cpylog import get_logger
import pyNastran
from pyNastran.bdf.bdf import BDF, read_bdf
#from pyNastran.bdf.test.test_bdf import run_bdf
from pyNastran.bdf.cards.test.utils import save_load_deck

PKG_PATH = Path(pyNastran.__path__[0])
MODEL_PATH = PKG_PATH / '..' / 'models'
assert os.path.exists(MODEL_PATH), MODEL_PATH

class TestAcoustic(unittest.TestCase):
    """
    The cards are:
     * PMIC
     * MATPOR
     * ACPLNW
     * AMLREG

    """
    def test_panel(self):
        log = get_logger(level='warning')
        model = BDF(log=log)
        panel_id = 24
        names = ['CAT', 'FROG']
        set_ids = [37, 42]
        model.add_panel(panel_id, names, set_ids, comment='panel')

    def test_pacabs(self):
        log = get_logger(level='warning')
        model = BDF(log=log)
        pid = 10
        cutfr = 20.0
        m = 1.0
        b = 2.0
        k = 3.0
        model.add_pacabs(pid, cutfr,
                         b, k, m)
        save_load_deck(model)

    def test_acmodl_nx(self):
        log = get_logger(level='warning')

        model = BDF(log=log, mode='nx')
        infor = 'CAT'
        fset = 11
        sset = 12
        acmodl = model.add_acmodl(infor, fset, sset,
                                  nastran_version='nx')
        acmodl.raw_fields()
        model.cross_reference()
        save_load_deck(model)

    def test_acmodl_msc(self):
        log = get_logger(level='warning')
        model = BDF(log=log, mode='msc')
        infor = 'CAT'
        fset = 11
        sset = 12
        acmodl = model.add_acmodl(infor, fset, sset,
                                  nastran_version='msc')
        acmodl.raw_fields()
        model.cross_reference()
        save_load_deck(model)

    def test_acoustic1(self):
        log = get_logger(level='warning')
        bdf_filename = MODEL_PATH / 'nx' / 'test_vba' / 'ac108vatv5tc.bdf'
        model = read_bdf(bdf_filename, log=log)
        save_load_deck(model, run_renumber=False)

    def test_acoustic2(self):
        log = get_logger(level='warning')
        bdf_filename = MODEL_PATH / 'nx' / 'test_vba' / 'acssn108presvar.bdf'
        model = read_bdf(bdf_filename, log=log, mode='nx')
        save_load_deck(model, run_convert=False, run_renumber=False)

    def test_acplnw(self):
        log = get_logger(level='warning')
        model = BDF(debug=False, log=log, mode='msc')
        sid = 1
        form = 'REAL'
        scale = 2.0
        real = 3.0
        imag = 4.0
        cid1 = 5
        xyz = [6., 7., 8.]
        cid2 = 9
        nxyz = [10., 11., 12.]
        comment = 'acplnw'
        model.add_acplnw(sid, form, scale, real, imag, cid1, xyz, cid2, nxyz,
                         comment=comment)

    def test_amlreg(self):
        log = get_logger(level='warning')
        model = BDF(debug=False, log=log, mode='msc')
        rid = 1
        sid = 2
        name = 'AMLREG test'
        infid = [3, 4, 5]
        model.add_amlreg(rid, sid, name, infid, nlayers=5, radsurf='AML', comment='amlreg')
        save_load_deck(model)

    def test_micpnt(self):
        log = get_logger(level='warning')
        model = BDF(debug=False, log=log, mode='msc')
        eid = 1
        nid = 2
        name = 'micpnt test'
        model.add_micpnt(eid, nid, name, comment='micpnt')
        save_load_deck(model)

    def test_pmic(self):
        log = get_logger(level='warning')
        model = BDF(debug=False, log=log, mode='msc')
        pid = 42
        model.add_pmic(pid, comment='pmic')
        save_load_deck(model, run_remove_unused=False)

    def test_matpor(self):
        log = get_logger(level='warning')
        model = BDF(debug=False, log=log, mode='msc')
        mid = 1
        rho = 1.0
        c = 2.0
        resistivity = 3.0
        porosity = 4.0
        tortuosity = 5.0
        comment = 'matpor'
        model.add_matpor_craggs(mid, rho, c, resistivity,
                                porosity, tortuosity, comment=comment)
        mid += 1

        frame = 'FRAME?'
        mu = 0.1
        gamma = 1.4
        prandtl_number = 0.27
        L1 = 532.
        L2 = 234.
        model.add_matpor_delmiki(mid, rho, c, resistivity, porosity, frame, density=0.0, comment='matpor')
        mid += 1

        model.add_matpor_jca(mid, rho, c, resistivity, porosity, tortuosity,
                             frame, gamma, prandtl_number, mu, L1, L2, density=0.0, comment='')
        save_load_deck(model, run_remove_unused=False)

    def test_chacab(self):
        log = get_logger(level='warning')
        model = BDF(debug=False, log=log, mode='msc')
        eid = 10
        pid = 11
        nodes = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20]
        for nid in nodes:
            model.add_grid(nid+1, [float(nid), 0., 0.])
        elem = model.add_chacab(eid, pid, nodes, comment='chacbr')
        str(elem)
        model.pop_parse_errors()
        model.cross_reference()
        model.pop_xref_errors()
        save_load_deck(model, run_remove_unused=False, run_test_bdf=False, run_save_load_hdf5=False)

    def test_chacbr(self):
        log = get_logger(level='warning')
        model = BDF(debug=False, log=log, mode='msc')
        eid = 10
        pid = 11
        nodes = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20]
        for nid in nodes:
            model.add_grid(nid+1, [float(nid), 0., 0.])
        elem = model.add_chacbr(eid, pid, nodes, comment='chacbr')
        str(elem)
        model.pop_parse_errors()
        model.cross_reference()
        model.pop_xref_errors()
        save_load_deck(model, run_remove_unused=False, run_test_bdf=False, run_save_load_hdf5=False)

    def test_caabsf_paabsf(self):
        log = get_logger(level='warning')
        model = BDF(debug=False, log=log, mode='msc')
        eid = 10
        pid = 11
        nodes = [1, 2, 3, 4]
        for nid in nodes:
            model.add_grid(nid+1, [float(nid), 0., 0.])
        elem = model.add_caabsf(eid, pid, nodes, comment='caabsf')
        prop = model.add_paabsf(pid, tzreid=None, tzimid=None,
                                s=1.0, a=1.0, b=0.0, k=0.0, rhoc=1.0, comment='paabsf')
        str(elem)
        model.pop_parse_errors()
        model.cross_reference()
        model.pop_xref_errors()
        save_load_deck(model, run_remove_unused=False, run_test_bdf=False, run_save_load_hdf5=False, run_convert=False)

    def test_pacbar(self):
        log = get_logger(level='warning')
        model = BDF(debug=False, log=log, mode='msc')
        eid = 10
        pid = 11
        nodes = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20]
        for nid in nodes:
            model.add_grid(nid+1, [float(nid), 0., 0.])
        mback = 1.0
        mseptm = 2.0
        freson = 3.0
        kreson = 4.0
        elem = model.add_pacbar(pid, mback, mseptm, freson, kreson, comment='pacbar')
        str(elem)
        model.pop_parse_errors()
        model.cross_reference()
        model.pop_xref_errors()
        save_load_deck(model, run_remove_unused=False, run_test_bdf=False, run_save_load_hdf5=False, run_convert=False)

if __name__ == '__main__':  # pragma: no cover
    unittest.main()
