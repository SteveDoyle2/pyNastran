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

    def test_ddval(self):
        """tests a DDVAL"""
        model = BDF(debug=False)
        oid = 10
        ddvals = [0.1, 0.2, 0.3, 0.5]
        ddval = model.add_ddval(oid, ddvals, comment='ddval')
        ddval.write_card(size=8)
        ddval.write_card(size=16)
        ddval.raw_fields()
        model.validate()
        model.cross_reference()

    def test_dlink(self):
        """tests a DLINK"""
        model = BDF(debug=False)
        oid = 10 # optimization id
        ddvid = 11 # dlink id?
        IDv = [20, 21, 22]
        Ci = [1.0, 1.0, 1.0]
        dlink = model.add_dlink(oid, ddvid, IDv, Ci, c0=0., cmult=1.,
                                comment='dlink')
        dlink.raw_fields()
        dlink.comment = ''
        msg = dlink.write_card(size=8)
        dlink.write_card(size=16)
        assert '$' not in msg, msg
        lines = msg.split('\n')

        model2 = BDF()
        model2.add_card(lines, 'DLINK', is_list=False)
        model.dlinks[10]

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
        mid1 = 4
        mp_min = 1e6
        mp_max = 1e7
        dvids = 11
        coeffs = 1.0
        dvmrel_1 = model.add_dvmrel1(oid, 'MAT1', mid1, 'E', dvids, coeffs,
                                     mp_min=mp_min, mp_max=mp_max, c0=0., validate=True,
                                     comment='dmvrel')

        oid = 11
        mid8 = 8
        dvmrel_8 = model.add_dvmrel1(oid, 'MAT8', mid8, 'NU12', dvids, coeffs,
                                     mp_min=0.25, mp_max=0.3, c0=0., validate=True,
                                     comment='dmvrel')
        oid = 12
        mid10 = 10
        dvmrel_10 = model.add_dvmrel1(oid, 'MAT10', mid10, 'RHO', dvids, coeffs,
                                      mp_min=0.1, mp_max=0.2, c0=0., validate=True,
                                      comment='dmvrel')

        E = 30.e7
        G = None
        nu = 0.3
        #nids = [1, 2, 3]
        mat1 = model.add_mat1(mid1, E, G, nu, rho=0.1, comment='mat1')

        e11 = 3e7
        e22 = 0.2 * e11
        nu12 = 0.3
        mat8 = model.add_mat8(mid8, e11, e22, nu12, comment='mat8')
        #bulk = c ** 2. * rho
        bulk = None
        rho = 0.1
        c = 4000.
        mat10 = model.add_mat10(mid10, bulk, rho, c, ge=0.0, comment='mat10')

        dvmrel_1.raw_fields()
        dvmrel_1.write_card(size=8)
        dvmrel_1.write_card(size=16)

        dvmrel_8.raw_fields()
        dvmrel_8.write_card(size=8)
        dvmrel_8.write_card(size=16)

        dvmrel_10.raw_fields()
        dvmrel_10.write_card(size=8)
        dvmrel_10.write_card(size=16)

        mat8.raw_fields()
        mat8.write_card(size=8)
        mat8.write_card(size=16)

        mat10.raw_fields()
        mat10.write_card(size=8)
        mat10.write_card(size=16)

        model.validate()
        model.cross_reference()

        dvmrel_1.raw_fields()
        dvmrel_8.raw_fields()
        dvmrel_10.raw_fields()
        mat8.raw_fields()
        mat10.raw_fields()

    def test_dvcrel1(self):
        """tests a DVCREL"""
        model = BDF()
        oid = 10
        eid = 100
        cp_min = 0.01
        cp_max = 1.
        desvar_id = 11
        desvar_ids = 11
        coeffs = 1.0
        dvcrel1 = model.add_dvcrel1(oid, 'CONM2', eid, 'X2', desvar_ids, coeffs,
                                    cp_min, cp_max, c0=0., validate=True, comment='')

        label = 'X2_MASS'
        xinit = 0.1
        xlb = 0.01
        xub = 1.0
        desvar = model.add_desvar(desvar_id, label, xinit, xlb, xub, comment='desvar')

        mass = 1.
        nid1 = 100
        nid2 = 101
        conm2 = model.add_conm2(eid, nid1, mass, cid=0, X=None, I=None,
                                comment='conm2')
        model.add_grid(100, xyz=[1., 2., 3.])
        model.add_grid(101, xyz=[2., 2., 4.])

        eid = 101
        pid = 102
        x = [1., 0., 0.]
        g0 = None
        cbar = model.add_cbar(eid, pid, nid1, nid2, x, g0, offt='GGG', pa=0, pb=0,
                              wa=None, wb=None, comment='cbar')

        oid = 11
        equation_id = 100
        dvcrel2 = model.add_dvcrel2(oid, 'CBAR', eid, 'X3', equation_id, desvar_ids, labels=None,
                                    cp_min=2., cp_max=4.,
                                    validate=True, comment='dvcrel2')

        mid = 1000
        dim = [1., 2., 0.1, 0.2]
        model.add_pbarl(pid, mid, 'T', dim, group='MSCBMLO', nsm=0.1,
                        comment='pbarl')
        E = 30.e7
        G = None
        nu = 0.3
        mat1 = model.add_mat1(mid, E, G, nu, rho=0.1, comment='mat1')

        name = 'fx2'
        eqs = ['fx2(x) = x + 10.']
        deqatn = model.add_deqatn(name, equation_id, eqs, comment='deqatn')

        dvcrel1.raw_fields()
        dvcrel2.raw_fields()
        dvcrel1.write_card(size=16)
        dvcrel2.write_card(size=16)

        dvcrel1.comment = ''
        dvcrel2.comment = ''
        desvar.comment = ''
        dvcrel1_msg = dvcrel1.write_card(size=8)
        dvcrel2_msg = dvcrel2.write_card(size=8)
        desvar_msg = desvar.write_card(size=8)


        model.validate()
        model.cross_reference()
        dvcrel1.raw_fields()
        dvcrel1.write_card(size=16)
        dvcrel1.write_card(size=8)

        dvcrel2.raw_fields()
        dvcrel2.write_card(size=16)
        dvcrel2.write_card(size=8)

        assert cbar.Mass() > 0, cbar.Mass()

        #-------------------------------------------
        dvcrel1_lines = dvcrel1_msg.split('\n')
        dvcrel2_lines = dvcrel2_msg.split('\n')
        desvar_lines = desvar_msg.split('\n')
        model2 = BDF(debug=False)
        model2.add_card(dvcrel1_lines, 'DVCREL1', is_list=False)
        model2.add_card(dvcrel2_lines, 'DVCREL2', is_list=False)
        model2.add_card(desvar_lines, 'DESVAR', is_list=False)

if __name__ == '__main__':  # pragma: no cover
    unittest.main()
