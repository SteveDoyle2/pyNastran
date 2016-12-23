# coding: utf-8
from __future__ import print_function
import os
import unittest

from six import iteritems, StringIO

import pyNastran
from pyNastran.bdf.bdf import BDF
from pyNastran.op2.op2 import OP2
from pyNastran.f06.test.f06_unit_tests import run_model

root_path = pyNastran.__path__[0]
#test_path = os.path.join(root_path, 'bdf', 'cards', 'test')
model_path = os.path.join(pyNastran.__path__[0], '..', 'models')

class TestOpt(unittest.TestCase):
    """
    The cards tested are:
     * DEQATN
    """
    @unittest.expectedFailure
    def test_opt_1(self):
        bdfname = os.path.join(model_path, 'sol200', 'model_200.bdf')
        op2name = os.path.join(model_path, 'sol200', 'model_200_nx.op2')
        bdf, op2 = run_model(bdfname, op2name,
                             f06_has_weight=False, vectorized=True,
                             encoding='utf-8')

        subcase_ids = op2.subcase_key.keys()
        for subcase_id in subcase_ids:
            assert isinstance(subcase_id, int), subcase_id
            for key, dresp in sorted(iteritems(bdf.dresps)):
                dresp.calculate(op2, subcase_id)

    def test_dvprel1(self):
        """tests a DESVAR, DVPREL1, DRESP1, DCONSTR"""
        model = BDF(debug=False)
        oid = 10
        desvar_id = 12
        desvar_ids = 12
        Type = 'PSHELL'
        pid = 20
        eid = 25
        mid = 30
        pname_fid = 'T'

        coeffs = 1.
        E = 30.e7
        G = None
        nu = 0.3
        nids = [1, 2, 3]

        label = 'T_SHELL'
        xinit = 0.1
        xlb = 0.01
        xub = 2.0

        model.add_grid(1, xyz=[0., 0., 0.])
        model.add_grid(2, xyz=[1., 0., 0.])
        model.add_grid(3, xyz=[1., 1., 0.])
        ctria3 = model.add_ctria3(eid, pid, nids, comment='ctria3')
        pshell = model.add_pshell(pid, mid1=30, t=0.1, comment='pshell')
        mat1 = model.add_mat1(mid, E, G, nu, rho=0.1, comment='mat1')
        desvar = model.add_desvar(desvar_id, label, xinit, xlb, xub, comment='desvar')
        dvprel1 = model.add_dvprel1(oid, Type, pid, pname_fid,
                                    desvar_ids, coeffs, p_min=None, p_max=1e20, c0=0.0,
                                    validate=True, comment='dvprel')

        deqation = 100
        dvids = desvar_id
        labels = None
        dvprel2 = model.add_dvprel2(oid+1, Type, pid, pname_fid, deqation,
                                   dvids, labels, p_min=None, p_max=1e20,
                                   validate=True, comment='')
        equation_id = 100
        name = 'fstress'
        eqs = ['fstress(x) = x + 10.']
        deqatn = model.add_deqatn(name, equation_id, eqs, comment='deqatn')
        #print(deqatn.object_attributes())
        #print(deqatn.func_str)
        #print(deqatn)

        dresp_id = 42
        label = 'STRESS1'
        response_type = 'STRESS'
        property_type = 'PSHELL'
        region = None
        atta = 9
        attb = None
        atti = pid
        dresp = model.add_dresp1(dresp_id, label, response_type,
                                 property_type, region,
                                 atta, attb, atti, validate=True, comment='dresp')
        dconstr = model.add_dconstr(oid, dresp_id, lid=-1.e20, uid=1.e20,
                                   lowfq=0., highfq=1.e20, comment='dconstr')
        desvar.write_card(size=8)
        desvar.write_card(size=16)
        dvprel1.write_card(size=8)
        dvprel1.write_card(size=16)
        dconstr.write_card(size=8)
        dconstr.write_card(size=16)
        dresp.write_card(size=8)
        dresp.write_card(size=16)

        model.validate()
        #model._verify_bdf(xref=False)
        model.cross_reference()

        desvar.write_card(size=8)
        desvar.write_card(size=16)
        desvar.raw_fields()
        dvprel1.write_card(size=8)
        dvprel1.write_card(size=16)
        dvprel1.raw_fields()
        dconstr.write_card(size=8)
        dconstr.write_card(size=16)
        dconstr.raw_fields()
        dresp.write_card(size=8)
        dresp.write_card(size=16)
        dresp.raw_fields()

        stringio = StringIO()
        model.write_bdf(stringio, close=False)
        stringio.getvalue()
        #model.uncross_reference()
        #model.cross_reference()
        #model._verify_bdf(xref=True)

    def test_dvmrel1(self):
        """tests a DVMREL1"""
        model = BDF()
        oid = 10
        Type = 'MAT1'
        mid = 4
        mp_name = 'E'
        mp_min = 1e6
        mp_max = 1e7
        dvids = 11
        coeffs = 1.0
        dvmrel = model.add_dvmrel1(oid, Type, mid, mp_name, dvids, coeffs,
                                   mp_min=mp_min, mp_max=mp_max, c0=0., validate=True,
                                   comment='dmvrel')

        E = 30.e7
        G = None
        nu = 0.3
        #nids = [1, 2, 3]
        mat1 = model.add_mat1(mid, E, G, nu, rho=0.1, comment='mat1')

        dvmrel.raw_fields()
        dvmrel.write_card(size=8)
        dvmrel.write_card(size=16)
        model.validate()
        model.cross_reference()

if __name__ == '__main__':  # pragma: no cover
    unittest.main()
