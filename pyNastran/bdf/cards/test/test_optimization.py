"""
defines optimization tests
"""
# coding: utf-8
from __future__ import print_function
import os
import unittest

from six import integer_types
import numpy as np

import pyNastran
from pyNastran.bdf.bdf import BDF, read_bdf
from pyNastran.op2.op2 import read_op2
from pyNastran.bdf.cards.test.utils import save_load_deck
from pyNastran.utils.log import get_logger
from pyNastran.bdf.cards.optimization import break_word_by_trailing_integer

model_path = os.path.join(pyNastran.__path__[0], '..', 'models')


class TestOpt(unittest.TestCase):
    """
    The cards tested are:
     * DEQATN
    """
    def test_opt_1(self):
        """tests SOL 200"""
        log = get_logger(level='warning')
        bdf_filename = os.path.join(model_path, 'sol200', 'model_200.bdf')
        model = read_bdf(bdf_filename, xref=True, debug=False)
        op2_filename = os.path.join(model_path, 'sol200', 'model_200.op2')
        #bdf, op2 = run_model(bdf_filename, op2_filename,
                             #f06_has_weight=False, vectorized=True,
                             #encoding='utf-8')
        op2 = read_op2(op2_filename, log=log, debug=True, debug_file='temp.debug')
        subcase_ids = op2.subcase_key.keys()
        #for subcase_id in subcase_ids:
            #assert isinstance(subcase_id, integer_types), subcase_id
            #for key, dresp in sorted(model.dresps.items()):
                #print(dresp)
                #dresp.calculate(op2, subcase_id)
        os.remove('temp.debug')

    def test_opt_2(self):
        """tests updating model based on DESVARs"""
        model = BDF(debug=False)
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [1., 0., 0.])
        model.add_grid(3, [1., 1., 0.])
        model.add_grid(4, [0., 1., 0.])
        eid = 1
        pid_pshell = 1
        nids = [1, 2, 3, 4]
        model.add_cquad4(eid, pid_pshell, nids)

        delta = 0.1
        #---------------------------------------------------
        model.add_pshell(pid_pshell, mid1=1, t=0.1)
        c2 = 0.5
        t2 = 1.0
        c0p = 0.1
        #      0.1 + 0.5*1.0 = 0.6
        tnew = c0p + c2 * (t2 + delta)
        model.add_desvar(2, 'T2', t2)
        desvars = [2] # desvar
        coeffs = [c2]
        prop_type = 'PSHELL'
        pname_fid = 'T'
        model.add_dvprel1(1, prop_type, pid_pshell, pname_fid, desvars, coeffs,
                          p_min=None, p_max=1e20,
                          c0=c0p, validate=True,
                          comment='')
        #---------------------------------------------------

        mid = 1
        E = 1.0
        G = None
        nu = 0.3
        model.add_mat1(mid, E, G, nu)

        c1 = 0.5
        e1 = 110.
        c0m = 100.
        model.add_desvar(1, 'E1', e1)
        mat_type = 'MAT1'
        mp_name = 'E'
        coeffs = [c1]
        desvars = [1]
        # 100 + 0.5 * 110 = 100 + 55 = 155
        enew = c0m + c1 * (e1 + delta)
        model.add_dvmrel1(2, mat_type, mid, mp_name, desvars, coeffs,
                          mp_min=None, mp_max=1e20,
                          c0=c0m, validate=True,
                          comment='')
        #---------------------------------------------------

        eid_conm2 = 2
        nid = 1
        mass = 0.1
        model.add_conm2(eid_conm2, nid, mass, cid=0, X=None, I=None, comment='')
        c3 = 1.0
        mass3 = 12.
        model.add_desvar(3, 'MASS', mass3)
        element_type = 'CONM2'
        cp_name = 'M'
        coeffs = [c3]
        desvars = [3]
        mass0 = 0.
        # ???
        mass_new = mass0 + c3 * (mass3 + delta)
        model.add_dvcrel1(3, element_type, eid_conm2, cp_name, desvars, coeffs, cp_min=None,
                         cp_max=1e20, c0=mass0,
                         validate=True, comment='')
        #---------------------------------------------------

        x14 = 37.
        model.add_desvar(4, 'X1', x14)
        element_type = 'CONM2'
        cp_name = 'X1'
        c4 = 1.
        x10 = 0.
        coeffs = [c4]
        desvars = [4]
        x1_new = x10 + c4 * (x14 + delta)
        model.add_dvcrel1(4, element_type, eid_conm2, cp_name, desvars, coeffs, cp_min=None,
                          cp_max=1e20, c0=x10,
                          validate=True, comment='')
        #---------------------------------------------------

        pid_pcomp = 2
        mids = [1]
        thicknesses = [1.]
        model.add_pcomp(pid_pcomp, mids, thicknesses, thetas=None, souts=None,
                        nsm=0., sb=0., ft=None,
                        tref=0., ge=0., lam=None,
                        z0=None, comment='')
        t5 = 5.
        model.add_desvar(5, 'T5', t5)
        prop_type = 'PCOMP'
        pname_fid = 'T1'
        desvars = [5] # desvar
        coeffs = [1.0]
        tpcomp_new = t5 + delta
        model.add_dvprel1(5, prop_type, pid_pcomp, pname_fid, desvars, coeffs,
                          p_min=None, p_max=1e20,
                          c0=0., validate=True,
                          comment='')
        #---------------------------------------------------

        model.validate()
        model.cross_reference()
        model.update_model_by_desvars()

        assert model.properties[pid_pshell].t == tnew, 't=%s tnew=%s' % (model.properties[pid_pshell].t, tnew)
        assert model.materials[mid].e == enew, 'E=%s enew=%s' % (model.materials[mid].e, enew)
        assert model.Mass(eid_conm2).mass == mass_new, 'mass=%s mass_new=%s' % (model.Mass(eid_conm2).mass, mass_new)
        assert model.Mass(eid_conm2).X[0] == x1_new, 'X1=%s x1_new=%s' % (model.Mass(eid_conm2).mass, x1_new)
        assert model.properties[pid_pcomp].thicknesses[0] == tpcomp_new, 't=%s tnew=%s' % (model.properties[pid_pcomp].thicknesses[0], tpcomp_new)

    def test_ddval(self):
        """tests a DDVAL"""
        model = BDF(debug=False)
        oid = 10
        ddvals = [0.1, 0.2, 0.3, 0.5]
        ddval = model.add_ddval(oid, ddvals, comment='ddval')
        ddval.write_card(size=8)
        ddval.write_card(size=16)
        ddval.write_card(size=16, is_double=True)
        ddval.raw_fields()
        model.validate()
        model.cross_reference()
        save_load_deck(model)

    def test_doptprm(self):
        """tests a doptprm"""
        #DOPTPRM    CONV1  .00001  DELOBJ .000001  DESMAX     100      P1       1
        #              P2      13
        model = BDF(debug=False)

        params = {
            'CONV1' : 0.0001,
            'DELOBJ' : 0.000001,
            'DESMAX' : 100,
            'P1' : 1,
            'P2' : 13,
        }
        doptprm = model.add_doptprm(params, comment='doptprm')
        model.validate()
        model._verify_bdf(xref=False)
        model.cross_reference()
        model._verify_bdf(xref=True)

        doptprm.comment = ''
        doptprm.raw_fields()
        doptprm.write_card(size=8)
        doptprm.write_card(size=16)
        doptprm.write_card(size=16, is_double=True)
        save_load_deck(model)


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
        dlink.write_card(size=16, is_double=True)
        assert '$' not in msg, msg
        lines = msg.split('\n')

        model2 = BDF(debug=False)
        model2.add_card(lines, 'DLINK', is_list=False)
        dlink = model.dlinks[10]
        dlink.write_card()
        save_load_deck(model)

    def test_dvprel1(self):
        """tests a DESVAR, DVPREL1, DVPREL2, DRESP1, DRESP2, DRESP3, DCONSTR, DSCREEN, DCONADD"""
        model = BDF(debug=False)
        dvprel1_id = 10
        desvar_id = 12
        desvar_ids = 12
        prop_type = 'PSHELL'
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
        p_max = 1e20

        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [1., 0., 0.])
        model.add_grid(3, [1., 1., 0.])
        model.add_ctria3(eid, pid, nids, comment='ctria3')
        model.add_pshell(pid, mid1=30, t=0.1, comment='pshell')
        model.add_mat1(mid, E, G, nu, rho=0.1, comment='mat1')
        desvar = model.add_desvar(desvar_id, label, xinit, xlb, xub, comment='desvar')
        dvprel1 = model.add_dvprel1(dvprel1_id, prop_type, pid, pname_fid,
                                    desvar_ids, coeffs, p_min=None, p_max=p_max, c0=0.0,
                                    validate=True, comment='dvprel1')

        dvprel2_id = dvprel1_id + 1
        deqation = 100
        dvids = desvar_id
        labels = None
        dvprel2 = model.add_dvprel2(dvprel2_id, prop_type, pid, pname_fid, deqation,
                                    dvids, labels, p_min=None, p_max=p_max,
                                    validate=True, comment='dvprel2')

        equation_id = 100
        eqs = ['fstress(x) = x + 10.']
        model.add_deqatn(equation_id, eqs, comment='deqatn')
        #print(deqatn.object_attributes())
        #print(deqatn.func_str)
        #print(deqatn)

        dresp1_id = 42
        label = 'STRESS1'
        response_type = 'STRESS'
        property_type = 'PSHELL'
        region = None
        atta = 9
        attb = None
        atti = pid
        dresp1 = model.add_dresp1(dresp1_id, label, response_type,
                                  property_type, region,
                                  atta, attb, atti, validate=True, comment='dresp1')
        assert dresp1.rtype == dresp1.response_type
        assert dresp1.ptype == dresp1.property_type
        dresp1.rtype = response_type
        dresp1.ptype = property_type

        dconstr = model.add_dconstr(dresp1_id, dresp1_id, lid=-1.e20, uid=1.e20,
                                    lowfq=0., highfq=1.e20, comment='dconstr')

        params = {
            (0, 'DRESP1') : [42],
            (1, 'DESVAR') : [12],
            (3, 'DNODE') : [[100, 101], [1, 2]],
        }
        dresp2_id = 43
        dequation = equation_id
        label = 'dresp2'
        region = None
        dresp2 = model.add_dresp2(dresp2_id, label, dequation, region, params,
                                  method='MIN', c1=100., c2=0.005, c3=None,
                                  comment='dresp2')
        dresp2.raw_fields()

        dresp3_id = 44
        label = 'dresp3'
        group = 'cat'
        Type = 'dog'
        region = None
        params = {
            (0, 'DRESP1') : [42],
            (1, 'DESVAR') : [12],
            (2, 'DRESP2') : [dresp2_id],
            (2, 'DVPREL1') : [dvprel1_id],
            (3, 'DVPREL2') : [dvprel2_id],
            (3, 'DNODE') : [[100, 101], [1, 2]],
        }
        dresp3 = model.add_dresp3(dresp3_id, label, group, Type, region,
                                  params, comment='dresp3')
        dresp3.raw_fields()

        oid = 1001
        dconstr = model.add_dconstr(oid, dresp1_id, lid=-1.e20, uid=1.e20,
                                   lowfq=0., highfq=1.e20, comment='dconstr1')
        oid = 1002
        dconstr = model.add_dconstr(oid, dresp2_id, lid=-1.e20, uid=1.e20,
                                    lowfq=0., highfq=1.e20)
        oid = 1003
        dconstr = model.add_dconstr(oid, dresp3_id, lid=-1.e20, uid=1.e20,
                                    lowfq=0., highfq=1.e20)

        oid = 45
        dconstrs = [1001, 1002, 1003]
        dconadd = model.add_dconadd(oid, dconstrs, comment='dconadd')

        dscreen = model.add_dscreen('dunno', comment='dscreen')
        dscreen.raw_fields()

        #print(dresp3)
        grid = model.add_grid(100, [0., 0., 0.])
        model.add_grid(101, [0., 0., 0.])
        model.pop_parse_errors()

        desvar.write_card(size=8)
        desvar.write_card(size=16)
        dvprel1.write_card(size=8)
        dvprel1.write_card(size=16)
        dconstr.write_card(size=8)
        dconstr.write_card(size=16)
        dresp1.write_card(size=8)
        dresp1.write_card(size=16)
        dresp1.write_card(size=16, is_double=True)
        dresp2.write_card(size=8)
        dresp2.write_card(size=16)
        dresp2.write_card(size=16, is_double=True)
        dresp3.write_card(size=8)
        dresp3.write_card(size=16)
        dresp3.write_card(size=16, is_double=True)
        dvprel2.write_card(size=8)
        dvprel2.write_card(size=16)
        dvprel2.write_card(size=16, is_double=True)
        dconadd.write_card(size=8)
        dconadd.write_card(size=16)
        dconadd.write_card(size=16, is_double=True)
        dscreen.write_card(size=8)
        dscreen.write_card(size=16)
        dscreen.write_card(size=16, is_double=True)

        model.validate()
        model._verify_bdf(xref=False)
        model.cross_reference()
        model.pop_xref_errors()

        desvar.write_card(size=8)
        desvar.write_card(size=16)
        desvar.raw_fields()
        dvprel1.write_card(size=8)
        dvprel1.write_card(size=16)
        dvprel1.raw_fields()
        dconstr.write_card(size=8)
        dconstr.write_card(size=16)
        dconstr.raw_fields()
        dresp1.write_card(size=8)
        dresp1.write_card(size=16)
        dresp1.raw_fields()
        dresp2.write_card(size=8)
        dresp2.write_card(size=16)
        dresp2.write_card(size=16, is_double=True)
        dresp3.write_card(size=8)
        dresp3.write_card(size=16)
        dresp3.write_card(size=16, is_double=True)
        dvprel2.write_card(size=8)
        dvprel2.write_card(size=16)
        dvprel2.write_card(size=16, is_double=True)
        dconadd.write_card(size=8)
        dconadd.write_card(size=16)
        dconadd.write_card(size=16, is_double=True)

        grid.nid = 200
        assert '200' in str(dresp3), dresp3

        save_load_deck(model, run_convert=False)

    def test_dvprel1_02(self):
        model = BDF()
        oid = 1
        pid = 2
        prop_type = 'PCOMP'
        pname_fid = 'THETA11'
        dvids = [1, 2]
        coeffs = [1., 2.]
        dvprel1a = model.add_dvprel1(oid, prop_type, pid, pname_fid, dvids, coeffs,
                                     p_min=None, p_max=1e20,
                                     c0=0.0, validate=True,
                                     comment='')

        oid = 2
        pname_fid = 'T42'
        dvprel1b = model.add_dvprel1(oid, prop_type, pid, pname_fid, dvids, coeffs,
                                     p_min=None, p_max=1e20,
                                     c0=0.0, validate=True,
                                     comment='')
        assert 'THETA11' in dvprel1a.raw_fields(), dvprel1a
        assert 'T42' in dvprel1b.raw_fields(), dvprel1b

    def test_dvmrel1(self):
        """tests a DVMREL1"""
        model = BDF(debug=False)

        model.add_desvar(11, 'X11', 1.0)
        oid = 10
        mid1 = 4
        mp_min = 1e6
        mp_max = 1e7
        dvids = 11
        coeffs = 1.0
        dvmrel1_1 = model.add_dvmrel1(oid, 'MAT1', mid1, 'E', dvids, coeffs,
                                      mp_min=mp_min, mp_max=mp_max, c0=0., validate=True,
                                      comment='dmvrel')

        oid = 11
        mid8 = 8
        dvmrel1_8 = model.add_dvmrel1(oid, 'MAT8', mid8, 'NU12', dvids, coeffs,
                                      mp_min=0.25, mp_max=0.3, c0=0., validate=True,
                                      comment='dmvrel')
        oid = 12
        mid10 = 10
        dvmrel1_10 = model.add_dvmrel1(oid, 'MAT10', mid10, 'RHO', dvids, coeffs,
                                       mp_min=0.1, mp_max=0.2, c0=0., validate=True,
                                       comment='dmvrel')

        oid = 21
        deqation = 42
        mp_name = 'E'
        mat_type = 'MAT1'
        labels = []
        dvmrel2_1 = model.add_dvmrel2(oid, mat_type, mid1, mp_name, deqation,
                                      dvids, labels, mp_min=None, mp_max=1e20,
                                      validate=True,
                                      comment='dvmrel')
        E = 30.e7
        G = None
        nu = 0.3
        #nids = [1, 2, 3]
        model.add_mat1(mid1, E, G, nu, rho=0.1, comment='mat1')

        e11 = 3e7
        e22 = 0.2 * e11
        nu12 = 0.3
        mat8 = model.add_mat8(mid8, e11, e22, nu12, comment='mat8')
        #bulk = c ** 2. * rho
        bulk = None
        rho = 0.1
        c = 4000.
        mat10 = model.add_mat10(mid10, bulk, rho, c, ge=0.0, comment='mat10')

        equation_id = 42
        eqs = ['fstress(x) = x + 10.']
        model.add_deqatn(equation_id, eqs, comment='deqatn')


        dvmrel1_1.raw_fields()
        dvmrel1_1.write_card(size=8)
        dvmrel1_1.write_card(size=16)

        dvmrel1_8.raw_fields()
        dvmrel1_8.write_card(size=8)
        dvmrel1_8.write_card(size=16)

        dvmrel1_10.raw_fields()
        dvmrel1_10.write_card(size=8)
        dvmrel1_10.write_card(size=16)

        dvmrel2_1.raw_fields()
        dvmrel2_1.write_card(size=8)
        dvmrel2_1.write_card(size=16)

        mat8.raw_fields()
        mat8.write_card(size=8)
        mat8.write_card(size=16)

        mat10.raw_fields()
        mat10.write_card(size=8)
        mat10.write_card(size=16)

        model.validate()
        model.cross_reference()

        dvmrel1_1.raw_fields()
        dvmrel1_8.raw_fields()
        dvmrel1_10.raw_fields()
        dvmrel2_1.raw_fields()
        mat8.raw_fields()
        mat10.raw_fields()
        save_load_deck(model, run_convert=False)

    def test_dvcrel1(self):
        """tests a DVCREL1, DVCREL2, DVGRID"""
        model = BDF(debug=False)
        oid = 10
        conm2_eid = 100
        cp_min = 0.01
        cp_max = 1.
        desvar_id = 11
        desvar_ids = 11
        coeffs = 1.0
        dvcrel1 = model.add_dvcrel1(oid, 'CONM2', conm2_eid, 'X2', desvar_ids, coeffs,
                                    cp_min, cp_max, c0=0., validate=True, comment='dvcrel')
        dvcrel1.Type = 'CONM2'
        assert dvcrel1.Type == dvcrel1.element_type
        #dvcrel1.pMax = p_max
        #dvcrel1.pMin = None

        label = 'X2_MASS'
        xinit = 0.1
        xlb = 0.01
        xub = 1.0
        desvar = model.add_desvar(desvar_id, label, xinit, xlb, xub, comment='desvar')

        mass = 1.
        nid1 = 100
        nid2 = 101
        conm2 = model.add_conm2(conm2_eid, nid1, mass, cid=0, X=None, I=None,
                                comment='conm2')
        model.add_grid(100, [1., 2., 3.])
        model.add_grid(101, [2., 2., 4.])

        eid = 101
        pid = 102
        x = [1., 0., 0.]
        g0 = None
        cbar = model.add_cbar(eid, pid, [nid1, nid2], x, g0, offt='GGG', pa=0, pb=0,
                              wa=None, wb=None, comment='cbar')

        oid = 11
        equation_id = 100
        dvcrel2 = model.add_dvcrel2(oid, 'CBAR', eid, 'X3', equation_id, desvar_ids, labels=None,
                                    cp_min=2., cp_max=4.,
                                    validate=True, comment='dvcrel2')
        dvcrel2.Type = 'CBAR'
        assert dvcrel2.Type == dvcrel2.element_type

        mid = 1000
        dim = [1., 2., 0.1, 0.2]
        model.add_pbarl(pid, mid, 'T', dim, group='MSCBML0', nsm=0.1,
                        comment='pbarl')
        E = 30.e7
        G = None
        nu = 0.3
        model.add_mat1(mid, E, G, nu, rho=0.1, comment='mat1')

        eqs = ['fx2(x) = x + 10.']
        deqatn = model.add_deqatn(equation_id, eqs, comment='deqatn')

        nid = 100
        dvid = 10000
        dxyz = [1., 2., 3.]
        dvgrid1 = model.add_dvgrid(dvid, nid, dxyz, cid=0, coeff=1.0,
                                   comment='dvgrid')

        nid = 101
        dvid = 10001
        dxyz = np.array([1., 2., 3.])
        dvgrid2 = model.add_dvgrid(dvid, nid, dxyz, cid=0, coeff=1.0,
                                   comment='dvgrid')

        model.pop_parse_errors()

        dvcrel1.raw_fields()
        dvcrel1.write_card(size=16)
        dvcrel2.raw_fields()
        dvcrel2.write_card(size=16)
        dvgrid1.raw_fields()
        dvgrid1.write_card(size=16)

        dvcrel1.comment = ''
        dvcrel2.comment = ''
        desvar.comment = ''
        dvgrid1.comment = ''
        dvcrel1_msg = dvcrel1.write_card(size=8)
        dvcrel2_msg = dvcrel2.write_card(size=8)
        desvar_msg = desvar.write_card(size=8)
        dvgrid_msg = dvgrid1.write_card(size=8)

        model.validate()
        model.cross_reference()
        model.pop_xref_errors()

        dvcrel1.raw_fields()
        dvcrel1.write_card(size=16)
        dvcrel1.write_card(size=8)

        dvcrel2.raw_fields()
        dvcrel2.write_card(size=16)
        dvcrel2.write_card(size=8)

        dvgrid1.raw_fields()
        dvgrid1.write_card(size=16)
        dvgrid1.write_card(size=8)

        deqatn.write_card()
        assert cbar.Mass() > 0, cbar.Mass()
        model.uncross_reference()

        #-------------------------------------------
        dvcrel1_lines = dvcrel1_msg.split('\n')
        dvcrel2_lines = dvcrel2_msg.split('\n')
        desvar_lines = desvar_msg.split('\n')
        dvgrid_lines = dvgrid_msg.split('\n')
        model2 = BDF(debug=True)

        model2.add_card(dvcrel1_lines, 'DVCREL1', is_list=False)
        model2.add_card(dvcrel2_lines, 'DVCREL2', is_list=False)
        model2.add_card(desvar_lines, 'DESVAR', is_list=False)
        model2.add_card(dvgrid_lines, 'DVGRID', is_list=False)
        #model2.add_conm2(conm2_eid, nid1, mass, cid=0, X=None, I=None,
                         #comment='conm2')
        #model2.add_grid(100, [1., 2., 3.])
        #model2.add_grid(101, [2., 2., 4.])
        #save_load_deck(model2)

    def test_break_words(self):
        """tests break_word_by_trailing_integer"""
        assert break_word_by_trailing_integer('T11') == ('T', '11'), break_word_by_trailing_integer('T11')
        assert break_word_by_trailing_integer('THETA42') == ('THETA', '42'), break_word_by_trailing_integer('THETA42')
        assert break_word_by_trailing_integer('T3') == ('T', '3'), break_word_by_trailing_integer('T3')
        with self.assertRaises(SyntaxError):
            assert break_word_by_trailing_integer('THETA32X')

    def test_dvgrid(self):
        """tests DVGRID"""
        dvid = 1
        nid = 2
        dxyz = [1., 2., 3.]
        cid = 1
        model = BDF(debug=False)
        dvgrid = model.add_dvgrid(dvid, nid, dxyz, cid=cid, coeff=1.0, comment='')
        model.pop_parse_errors()
        model.cross_reference()
        save_load_deck(model)


if __name__ == '__main__':  # pragma: no cover
    unittest.main()

