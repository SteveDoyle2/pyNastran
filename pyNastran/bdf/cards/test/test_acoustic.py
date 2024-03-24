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
    def test_acoustic(self):
        log = get_logger(level='warning')
        bdf_filename = MODEL_PATH / 'nx' / 'test_vba' / 'ac108vatv5tc.bdf'
        model = read_bdf(bdf_filename, log=log)
        save_load_deck(model, run_renumber=False)

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
        save_load_deck(model, run_remove_unused=False)


if __name__ == '__main__':  # pragma: no cover
    unittest.main()
