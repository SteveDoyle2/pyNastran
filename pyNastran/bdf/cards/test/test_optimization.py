"""defines optimization tests"""
# coding: utf-8
import os
import unittest

import numpy as np
from cpylog import get_logger, SimpleLogger

import pyNastran
from pyNastran.bdf.bdf import BDF, read_bdf
from pyNastran.op2.op2 import OP2, read_op2
from pyNastran.bdf.cards.test.utils import save_load_deck
from pyNastran.bdf.cards.optimization import (
    break_word_by_trailing_integer,
    validate_dvcrel, validate_dvmrel, validate_dvprel)

MODEL_PATH = os.path.join(pyNastran.__path__[0], '..', 'models')


class TestOpt(unittest.TestCase):
    """
    The cards tested are:
     * DEQATN
    """
    def test_dvcrel2(self):
        model = BDF(debug=None)
        oid = 20

        equation_id = 100
        eqs = ['fstress(x,y) = x + 10.']
        model.add_deqatn(equation_id, eqs, comment='deqatn')

        model.add_desvar(10, 'VAR', 1.0)
        model.add_desvar(11, 'VAR2', 3.0)
        desvar_ids = [10, 11]
        eid = 100
        pid = 101
        mid = 102
        model.add_grid(20, [0., 0., 0.])
        model.add_grid(21, [1., 0., 0.])
        model.add_mat1(mid, 3.0e7, None, 0.3)
        model.add_pbarl(pid, mid, 'BAR', [1., 2.,])
        model.add_cbar(eid, pid, [20, 21], [0., 1., 0.], None)
        model.add_dvcrel2(oid, 'CBAR', eid, 'X3', equation_id,
                          desvar_ids, labels=None,
                          cp_min=2., cp_max=4.,
                          validate=True,)
        save_load_deck(model, run_convert=False)

    def test_validate_dvprel(self):
        prop_type = 'PSHELL'
        pname_fid = 'T'
        validate = False
        validate_dvprel(prop_type, pname_fid, validate)

        validate = True
        prop_type = 'PSHELL'
        pname_fid = 'T'
        validate_dvprel(prop_type, pname_fid, validate)
        pname_fid = 'Z1'
        validate_dvprel(prop_type, pname_fid, validate)
        pname_fid = 'Z2'
        validate_dvprel(prop_type, pname_fid, validate)

        prop_type = 'PCOMP'
        pname_fid = 'T1'
        validate_dvprel(prop_type, pname_fid, validate)
        pname_fid = 'THETA3'
        validate_dvprel(prop_type, pname_fid, validate)

        prop_type = 'PCOMPG'
        pname_fid = 'Z0'
        validate_dvprel(prop_type, pname_fid, validate)
        # pname_fid = 'THETA3'
        # validate_dvprel(prop_type, pname_fid, validate)

        prop_type = 'PSHEAR'
        pname_fid = 'T'
        validate_dvprel(prop_type, pname_fid, validate)
        pname_fid = 'NSM'
        validate_dvprel(prop_type, pname_fid, validate)

        prop_type = 'PROD'
        pname_fid = 'A'
        validate_dvprel(prop_type, pname_fid, validate)

        prop_type = 'PTUBE'
        pname_fid = 'T'
        validate_dvprel(prop_type, pname_fid, validate)

        prop_type = 'PBAR'
        pname_fid = 'A'
        validate_dvprel(prop_type, pname_fid, validate)

        prop_type = 'PBARL'
        pname_fid = 'DIM1'
        validate_dvprel(prop_type, pname_fid, validate)

        prop_type = 'PBEAM'
        pname_fid = 'A'
        validate_dvprel(prop_type, pname_fid, validate)

        prop_type = 'PBEAML'
        pname_fid = 'DIM1'
        validate_dvprel(prop_type, pname_fid, validate)

        prop_type = 'PMASS'
        pname_fid = 3
        validate_dvprel(prop_type, pname_fid, validate)

        prop_type = 'PVISC'
        pname_fid = 'CE1'
        validate_dvprel(prop_type, pname_fid, validate)

        prop_type = 'PBUSH'
        pname_fid = 'K1'
        validate_dvprel(prop_type, pname_fid, validate)

        prop_type = 'PBUSH1D'
        pname_fid = 'K'
        validate_dvprel(prop_type, pname_fid, validate)

        prop_type = 'PBUSHT'
        pname_fid = 'TBID1'
        validate_dvprel(prop_type, pname_fid, validate)

        prop_type = 'PGAP'
        pname_fid = 'KA'
        validate_dvprel(prop_type, pname_fid, validate)

        prop_type = 'PDAMP'
        pname_fid = 'B1'
        validate_dvprel(prop_type, pname_fid, validate)

        prop_type = 'PWELD'
        pname_fid = 'D'
        validate_dvprel(prop_type, pname_fid, validate)

        prop_type = 'PFAST'
        pname_fid = 'KT1'
        validate_dvprel(prop_type, pname_fid, validate)

    def test_validate_dvmrel(self):
        validate = False
        element_type = 'MAT1'
        cp_name = 'E'
        validate_dvmrel(validate, element_type, cp_name)

        validate = True
        validate_dvmrel(validate, element_type, cp_name)

        element_type = 'MAT2'
        cp_name = 'G11'
        validate_dvmrel(validate, element_type, cp_name)

        element_type = 'MAT3'
        cp_name = 'EX'
        validate_dvmrel(validate, element_type, cp_name)

        # element_type = 'MAT4'
        # cp_name = 'E'
        # validate_dvmrel(validate, element_type, cp_name)

        # element_type = 'MAT5'
        # cp_name = 'E'
        # validate_dvmrel(validate, element_type, cp_name)

        element_type = 'MAT8'
        cp_name = 'E1'
        validate_dvmrel(validate, element_type, cp_name)

        element_type = 'MAT9'
        cp_name = 'G11'
        validate_dvmrel(validate, element_type, cp_name)

        element_type = 'MAT10'
        cp_name = 'BULK'
        validate_dvmrel(validate, element_type, cp_name)

        element_type = 'MAT11'
        cp_name = 'E3'
        validate_dvmrel(validate, element_type, cp_name)

    def test_validate_dvcrel(self):
        validate = False
        element_type = 'CONM2'
        cp_name = 'X1'
        validate_dvcrel(validate, element_type, cp_name)

        validate = True
        validate_dvcrel(validate, element_type, cp_name)

        element_type = 'CQUAD4'
        cp_name = 'T1'
        validate_dvcrel(validate, element_type, cp_name)

        element_type = 'CTRIA3'
        cp_name = 'T1'
        validate_dvcrel(validate, element_type, cp_name)

        element_type = 'CMASS4'
        cp_name = 'M'
        validate_dvcrel(validate, element_type, cp_name)

        element_type = 'CELAS4'
        cp_name = 'K'
        validate_dvcrel(validate, element_type, cp_name)

        element_type = 'CBEAM'
        cp_name = 'X2'
        validate_dvcrel(validate, element_type, cp_name)

        element_type = 'CBUSH'
        cp_name = 'X2'
        validate_dvcrel(validate, element_type, cp_name)

        # element_type = 'CBEAM'
        # cp_name = 'X2'
        # validate_dvcrel(validate, element_type, cp_name)

    def test_dmncon_dvtrel1_topvar(self):
        log = SimpleLogger(level='warning')
        model = BDF(log=log)

        options = {
            'STRESS': {'allowable': 10.0},
        }
        options_stress = {
            'SYM': {
                'cid': 0,
                'mirror_symmetry1': 'ZX',
                'mirror_symmetry2': 'ZX',
                'mirror_symmetry3': 'ZX',
                'cyclic_symmetry': 'X',
                'num_cyclic_symmetries': 5,
            },
        }
        opt_id = 3
        xinit = 0.1
        pid = 42
        prop_type = 'FAKE'
        label = 'cat'
        topvar1 = model.add_topvar(
            opt_id, label, prop_type, xinit, pid,
            xlb=0.001, delxv=0.2, power=3.0,
            options=None, comment='topvar')
        model.add_pelas(pid, 3.14)
        topvar1.raw_fields()
        topvar1.write_card(size=8, is_double=False)
        topvar1.write_card(size=16, is_double=False)
        opt_id = 4
        topvar2 = model.add_topvar(
            opt_id, label, prop_type, xinit, pid,
            xlb=0.001, delxv=0.2, power=3.0,
            options=options, comment='topvar')
        topvar2.raw_fields()

        opt_id = 5
        topvar2 = model.add_topvar(
            opt_id, label, prop_type, xinit, pid,
            xlb=0.001, delxv=0.2, power=3.0,
            options=options_stress, comment='topvar')
        topvar2.raw_fields()

        constraint_id = 1
        constraint_type = 'SYMP'
        xyz = [1., 2., 3.]
        normal = [0., 0., 1.]
        dmncon = model.add_dmncon(
            constraint_id, constraint_type,
            xyz=xyz, normal=normal, comment='dmncon',
        )
        dmncon.raw_fields()
        dmncon.write_card(size=8, is_double=False)
        dmncon.write_card(size=16, is_double=False)

        constraint_id += 1
        constraint_type = 'ADDM'
        dmncon = model.add_dmncon(
            constraint_id, constraint_type,
            xyz=xyz, normal=normal,
            angle=223., mind=42.0,
            comment='dmncon',
        )
        dmncon.write_card(size=8, is_double=False)

        # constraint_id += 1
        # constraint_type = 'SYMC'
        # dmncon = model.add_dmncon(
        #     constraint_id, constraint_type,
        #     xyz=xyz, normal=normal,
        #     m=[1., 2., 3.], nsections=4, comment='dmncon',
        # )
        # dmncon.write_card(size=8, is_double=False)

        constraint_id += 1
        constraint_type = 'CDID'
        dmncon = model.add_dmncon(
            constraint_id, constraint_type,
            xyz=xyz, normal=normal, d=[1], comment='dmncon',
        )
        dmncon.write_card(size=8, is_double=False)

        constraint_id += 1
        constraint_type = 'EXTC'
        dmncon = model.add_dmncon(
            constraint_id, constraint_type,
            normal=normal, comment='dmncon',
        )
        dmncon.write_card(size=8, is_double=False)

        dvt_id = 3
        group_id = 1
        label = 'label'
        dvtrel1 = model.add_dvtrel1(
            dvt_id, label, group_id,
            state='ACTIVE', dsv_flag=0, dvid1=0,
            validate=True, comment='dvtrel1')
        dvtrel1.raw_fields()
        dvtrel1.write_card()
        save_load_deck(model, run_save_load_hdf5=False)

    def test_opt_1(self):
        """tests SOL 200"""
        log = get_logger(level='warning')
        bdf_filename = os.path.join(MODEL_PATH, 'sol200', 'model_200.bdf')
        unused_model = read_bdf(bdf_filename, xref=True, debug=False)
        op2_filename = os.path.join(MODEL_PATH, 'sol200', 'model_200.op2')
        #bdf, op2 = run_model(bdf_filename, op2_filename,
                             #f06_has_weight=False, vectorized=True,
                             #encoding='utf-8')

        op2 = OP2(log=log, debug=True, debug_file='temp.debug')
        op2.read_op2(op2_filename)
        unused_subcase_ids = op2.subcase_key.keys()
        #for subcase_id in subcase_ids:
            #assert isinstance(subcase_id, integer_types), subcase_id
            #for key, dresp in sorted(model.dresps.items()):
                #print(dresp)
                #dresp.calculate(op2, subcase_id)
        os.remove('temp.debug')

    def test_opt_2(self):
        """tests updating model based on DESVARs"""
        log = SimpleLogger(level='warning')
        model = BDF(log=log)

        modtrak_id = 100
        low_range = 3
        high_range = 6
        mt_filter = 0.9
        model.add_modtrak(modtrak_id, low_range, high_range, mt_filter, comment='modtrak')

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
        save_load_deck(model)

    def test_ddval(self):
        """tests a DDVAL"""
        log = SimpleLogger(level='warning')
        model = BDF(log=log)
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
        log = SimpleLogger(level='warning')
        model = BDF(log=log)

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
        log = SimpleLogger(level='warning')
        model = BDF(log=log)
        oid = 10 # optimization id
        dvid = 11 # dlink id?
        IDv = [20, 21, 22]
        Ci = [1.0, 1.0, 1.0]
        dlink = model.add_dlink(oid, dvid, IDv, Ci, c0=0., cmult=1.,
                                comment='dlink')

        xinit = 0.1

        desvar_id = 11
        model.add_desvar(desvar_id, 'DV1', xinit, xlb=-1e20, xub=1e20, delx=None, ddval=None, comment='')

        desvar_id = 20
        model.add_desvar(desvar_id, 'DV1', xinit, xlb=-1e20, xub=1e20, delx=None, ddval=None, comment='')
        desvar_id = 21
        model.add_desvar(desvar_id, 'DV1', xinit, xlb=-1e20, xub=1e20, delx=None, ddval=None, comment='')
        desvar_id = 22
        model.add_desvar(desvar_id, 'DV1', xinit, xlb=-1e20, xub=1e20, delx=None, ddval=None, comment='')

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

    def test_dresp1s(self):
        model = BDF(debug=False)
        model.set_as_nx()

        pid = 1
        mid = 1
        E = 3.0e7
        G = None
        nu = 0.3
        model.add_pshell(pid, mid1=mid, t=1.0)
        model.add_mat1(mid, E, G, nu)
        #----------------------------------------------
        #no args
        dresp_id = 1
        label = 'COMPLIA'
        #DRESP1, 12, CMPL1, CMPLNCE
        response_type = 'CMPLNCE' # 'Compliance'
        property_type = None
        region = None
        atta = None
        attb = None
        atti = []
        model.add_dresp1(dresp_id, label, response_type, property_type, region, atta, attb, atti)

        dresp_id += 1
        label = 'DWEIGHT'
        response_type = 'DWEIGHT'
        model.add_dresp1(dresp_id, label, response_type, property_type, region, atta, attb, atti)

        dresp_id += 1
        label = 'WEIGHT'
        response_type = 'WEIGHT'
        model.add_dresp1(dresp_id, label, response_type, property_type, region, atta, attb, atti)

        #----------------------------------------------
        dresp_id += 1
        label = 'FLUTTER'
        response_type = 'FLUTTER'
        atti = [1, 2, 3, 4]
        property_type = 'PKNL'
        model.add_dresp1(dresp_id, label, response_type, property_type, region, atta, attb, atti)
        model.add_set1(1, [10]) # nodes
        model.add_flfact(2, [0.1]) # rho
        model.add_flfact(3, [0.9]) # mach
        model.add_flfact(4, [1100.]) # vel
        property_type = None
        #----------------------------------------------
        # reset
        atti = []

        #----------------------------------------------
        # atta only
        dresp_id += 1
        atta = 2
        label = 'CEIG'
        response_type = 'CEIG'
        model.add_dresp1(dresp_id, label, response_type, property_type, region, atta, attb, atti)

        dresp_id += 1
        atta = 2
        label = 'LAMA'
        response_type = 'LAMA'
        model.add_dresp1(dresp_id, label, response_type, property_type, region, atta, attb, atti)

        dresp_id += 1
        atta = 2
        label = 'EIGN'
        response_type = 'EIGN'
        model.add_dresp1(dresp_id, label, response_type, property_type, region, atta, attb, atti)
        #----------------------------------------------
        # atti only
        dresp_id += 1
        #atta = 2
        label = 'STABDER'
        response_type = 'STABDER'
        atti = [1]
        model.add_dresp1(dresp_id, label, response_type, property_type, region, atta, attb, atti)

        # reset
        atti = []
        #----------------------------------------------
        # property_type, atta only
        property_type = None
        # ---------------------------------------------
        # property_type, atta, atti
        dresp_id += 1
        label = 'STRESS'
        response_type = 'STRESS'
        property_type = 'PSHELL'
        atta = 4
        atti = [1]
        model.add_dresp1(dresp_id, label, response_type, property_type, region, atta, attb, atti)

        dresp_id += 1
        label = 'STRAIN'
        response_type = 'STRAIN'
        property_type = 'PSHELL'
        atta = 4
        atti = [1]
        model.add_dresp1(dresp_id, label, response_type, property_type, region, atta, attb, atti)

        dresp_id += 1
        label = 'FORCE'
        response_type = 'FORCE'
        property_type = 'PROD'
        atta = 2
        atti = [1]
        model.add_dresp1(dresp_id, label, response_type, property_type, region, atta, attb, atti)

        # attb is the lamina number
        dresp_id += 1
        label = 'CSTRESS'
        response_type = 'CSTRESS'
        property_type = 'PCOMP'
        atta = 4
        atti = [1]
        model.add_dresp1(dresp_id, label, response_type, property_type, region, atta, attb, atti)

        dresp_id += 1
        label = 'CSTRAIN'
        response_type = 'CSTRAIN'
        property_type = 'PCOMP'
        atta = 4
        atti = [1]
        model.add_dresp1(dresp_id, label, response_type, property_type, region, atta, attb, atti)

        dresp_id += 1
        label = 'CFAILURE'
        response_type = 'CFAILURE'
        property_type = 'PCOMP'
        atta = 3
        attai = [1]
        model.add_dresp1(dresp_id, label, response_type, property_type, region, atta, attb, atti)
        # ---------------------------------------------

        #  reset
        property_type = None
        atta = None
        atti = []
        # ---------------------------------------------

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

        deqation2 = 101
        labels = ['CAT']
        dvprel2b = model.add_dvprel2(dvprel2_id+1, prop_type, pid, pname_fid, deqation2,
                                     dvids, labels, p_min=None, p_max=p_max,
                                     validate=True, comment='dvprel2')

        dvids = []
        labels = ['TVAL']
        dvprel2 = model.add_dvprel2(dvprel2_id+2, prop_type, pid, pname_fid, deqation,
                                    dvids, labels, p_min=None, p_max=p_max,
                                    validate=True, comment='dvprel2')

        equation_id = 101
        eqs = ['fstress(x,y) = x * y + 10.']
        model.add_deqatn(equation_id, eqs, comment='deqatn')
        default_values = {'CAT': 42.0,
                          'TVAL': 0.3,}
        model.add_dtable(default_values, comment='dtable')


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
        dresp1.response_type = response_type
        dresp1.property_type = property_type

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

        dscreen = model.add_dscreen('DISP', comment='dscreen')
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

        # DCONADD (45) is part of dconstrs
        dconstr_keys = list(model.dconstrs.keys())
        dconstr_keys.sort()
        assert dconstr_keys == [42, 45, 1001, 1002, 1003], f'actual={dconstr_keys}'

        model2 = save_load_deck(model, run_convert=False)

        # DCONADD (45) is not part of self.dconstrs
        dconstr_keys = list(model.dconstrs.keys())
        dconstr_keys.sort()
        assert dconstr_keys == [42, 45, 1001, 1002, 1003], f'actual={dconstr_keys}'

        # 42 is a free DCONSTR
        dconaddi = model2.dconstrs[45][0]
        dconstr_subkeys = dconaddi.dconstrs
        dconstr_subkeys.sort()
        assert dconaddi.dconstrs == dconstrs, f'model2.dconadds={model2.dconadds}'
        assert dconstr_subkeys == [1001, 1002, 1003], f'actual={dconstr_subkeys}'

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
        log = SimpleLogger(level='warning')
        model = BDF(log=log)

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
        dvids = []
        labels = ['NUVAL']
        mp_name = 'NU'
        model.add_dvmrel2(oid+1, mat_type, mid1, mp_name, deqation, dvids,
                          labels, mp_min=None, mp_max=1e20, validate=True, comment='')
        model.add_dtable({'NUVAL': 0.3,})

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
        log = SimpleLogger(level='warning')
        model = BDF(log=log)
        oid = 10
        conm2_eid = 100
        cp_min = 0.01
        cp_max = 1.
        desvar_id = 11
        desvar_ids = 11
        coeffs = 1.0

        mass = 2.0
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
        unused_conm2 = model.add_conm2(conm2_eid, nid1, mass, cid=0, X=None, I=None,
                                       comment='conm2')
        model.add_grid(nid1, [1., 2., 3.])
        model.add_grid(nid2, [2., 2., 4.])

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
        model.add_desvar(dvid, 'X10000', 0.127, xlb=-1e20, xub=1e20, delx=None, ddval=None, comment='desvar10000')

        nid = 101
        dvid = 10001
        dxyz = np.array([1., 2., 3.])
        unused_dvgrid2 = model.add_dvgrid(dvid, nid, dxyz, cid=0, coeff=1.0,
                                          comment='dvgrid')
        model.add_desvar(dvid, 'X10001', 0.127, xlb=-1e20, xub=1e20, delx=None, ddval=None, comment='desvar10001')

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
        save_load_deck(model, run_renumber=False, run_convert=False)

    def test_break_words(self):
        """tests break_word_by_trailing_integer"""
        assert break_word_by_trailing_integer('T11') == ('T', '11'), break_word_by_trailing_integer('T11')
        assert break_word_by_trailing_integer('THETA42') == ('THETA', '42'), break_word_by_trailing_integer('THETA42')
        assert break_word_by_trailing_integer('T3') == ('T', '3'), break_word_by_trailing_integer('T3')
        with self.assertRaises(SyntaxError):
            assert break_word_by_trailing_integer('THETA32X')

    def test_dvgrid(self):
        """tests DVGRID"""
        log = SimpleLogger(level='warning')
        model = BDF(log=log)
        desvar_id = 1
        nid = 2
        dxyz = [1., 2., 3.]
        cid = 1
        dvgrid = model.add_dvgrid(desvar_id, nid, dxyz, cid=cid, coeff=1.0, comment='')
        model.add_grid(nid, [0., 0., 0.])

        origin = [0., 0., 0.]
        zaxis = [0., 0., 1.]
        xzplane = [1., 0., 0.]
        model.add_cord2r(cid, origin, zaxis, xzplane)
        xinit = 0.1
        model.add_desvar(desvar_id, 'DV1', xinit, xlb=-1e20, xub=1e20, delx=None, ddval=None, comment='')
        dvgrid.raw_fields()
        model.pop_parse_errors()
        model.cross_reference()
        save_load_deck(model)

    def test_rod_dvprel(self):
        """tests CROD and DVPREL1"""
        log = get_logger(level='warning')
        model = BDF(log=log)
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [1., 0., 0.])

        mid = 10
        E = 3.0e7
        G = None
        nu = 0.3
        model.add_mat1(mid, E, G, nu)
        #---------------------------------------------
        eid = 2
        pid = 2
        nids = [1, 2]
        model.add_crod(eid, pid, nids, comment='')
        model.add_conrod(eid+1, mid, nids, A=0.0, j=0.0, c=0.0, nsm=0.0, comment='')

        A = 2.0
        model.add_prod(pid, mid, A, j=0., c=0., nsm=0., comment='')

        OD1 = 1.
        unused_t = 0.1
        model.add_ctube(eid+2, pid+1, nids, comment='')
        model.add_ptube(pid+1, mid, OD1, t=None, nsm=0., OD2=None, comment='')

        dvprels = [
            #oid, pid, prop_type, pname_fid
            (1, pid, 'PROD', 'A'),
            (2, pid, 'PROD', 'J'),
            #-----------
            (3, pid+1, 'PTUBE', 'OD'),
            (4, pid+1, 'PTUBE', 'T'),
        ]
        coeffs = [1.0]
        xinit = 1.0
        for dvprel_data in dvprels:
            oid, pid, prop_type, pname_fid = dvprel_data
            dvids = [oid]
            desvar_id = oid

            dim_label = pname_fid
            prop_label = prop_type
            #if pname_fid.startswith('DIM'):
                #dim_label = 'D%s' % pname_fid[-1]
            #if prop_type == 'PBEAML':
                #prop_label = 'PBML'
            #elif prop_type == 'PBEAM':
                #prop_label = 'PBM'
            label = '%s%s%s' % (prop_label, dim_label, pid)

            model.add_desvar(desvar_id, label, xinit,
                             xlb=-1e20, xub=1e20, delx=None, ddval=None, comment='')
            model.add_dvprel1(oid, prop_type, pid, pname_fid, dvids, coeffs,
                              p_min=None, p_max=1e20, c0=0.0, validate=True, comment='')

        model.cross_reference()
        save_load_deck(model, xref='standard', punch=True)

    def test_cbar_dvprel(self):
        """tests CBAR and DVPREL1"""
        log = get_logger(level='warning')
        model = BDF(log=log)
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [1., 0., 0.])
        model.add_grid(3, [0., 0., 1.])

        mid = 10
        E = 3.0e7
        G = None
        nu = 0.3
        model.add_mat1(mid, E, G, nu)
        #---------------------------------------------
        eid = 2
        pid = 2
        nids = [1, 2]
        x = None
        g0 = 3
        wa = [0., 0., 0.1]
        wb = [0., 0., 0.1]
        model.add_cbar(eid, pid, nids, x, g0, offt='GGG', pa=0, pb=0, wa=wa, wb=wb, comment='')
        model.add_cbar(eid+1, pid+1, nids, x, g0, offt='GGG', pa=0, pb=0, wa=wa, wb=wb, comment='')

        model.add_cbeam(eid+2, pid+2, nids, x, g0, offt='GGG', bit=None, pa=0, pb=0,
                        wa=wa, wb=wb, sa=0, sb=0, comment='')
        model.add_cbeam(eid+3, pid+3, nids, x, g0, offt='GGG', bit=None, pa=0, pb=0,
                        wa=wa, wb=wb, sa=0, sb=0, comment='')

        model.add_pbar(pid, mid, A=0., i1=0., i2=0., i12=0., j=0., nsm=0.,
                       c1=0., c2=0., d1=0., d2=0., e1=0., e2=0., f1=0., f2=0.,
                       k1=1.e8, k2=1.e8, comment='')
        Type = 'BAR'
        dim = [0.1, 0.2]
        model.add_pbarl(pid+1, mid, Type, dim, group='MSCBML0', nsm=0., comment='')

        #---------------------
        xxb = [0.]
        so = ['YES']
        area = [1.0]
        i1 = [1.]
        i2 = [2.]
        i12 = [0.4]
        j = [0.3]
        nsm = [0.]
        model.add_pbeam(pid+2, mid, xxb, so, area, i1, i2, i12, j, nsm,
                        c1=None, c2=None, d1=None, d2=None, e1=None, e2=None, f1=None, f2=None,
                        k1=1., k2=1., s1=0., s2=0.,
                        nsia=0., nsib=None, cwa=0., cwb=None,
                        m1a=0., m2a=0., m1b=None, m2b=None,
                        n1a=0., n2a=0., n1b=None, n2b=None,
                        comment='')
        beam_type = 'BAR'
        dims = [dim]
        model.add_pbeaml(pid+3, mid, beam_type, xxb, dims, so=None, nsm=None,
                         group='MSCBML0', comment='')

        dvprels = [
            #oid, pid, prop_type, pname_fid
            (1, pid, 'PBAR', 'A'),
            (2, pid, 'PBAR', 'J'),
            (3, pid, 'PBAR', 'I1'),
            (4, pid, 'PBAR', 'I2'),
            (5, pid, 'PBAR', 'I12'),
            (6, pid+1, 'PBARL', 'DIM1'),
            (7, pid+1, 'PBARL', 'DIM2'),
            #-----------
            (8, pid+2, 'PBEAM', 'A'),
            (9, pid+2, 'PBEAM', 'J'),
            (10, pid+2, 'PBEAM', 'I1'),
            (11, pid+2, 'PBEAM', 'I2'),
            (12, pid+2, 'PBEAM', 'I12'),
            (13, pid+2, 'PBEAM', 'J'),
            (14, pid+3, 'PBEAML', 'DIM1'),
            (15, pid+3, 'PBEAML', 'DIM2'),
        ]
        coeffs = [1.0]
        xinit = 1.0
        for dvprel_data in dvprels:
            oid, pid, prop_type, pname_fid = dvprel_data
            dvids = [oid]
            desvar_id = oid

            dim_label = pname_fid
            prop_label = prop_type
            if pname_fid.startswith('DIM'):
                dim_label = 'D%s' % pname_fid[-1]
            if prop_type == 'PBEAML':
                prop_label = 'PBML'
            elif prop_type == 'PBEAM':
                prop_label = 'PBM'
            label = '%s%s%s' % (prop_label, dim_label, pid)

            model.add_desvar(desvar_id, label, xinit,
                             xlb=-1e20, xub=1e20, delx=None, ddval=None, comment='')
            model.add_dvprel1(oid, prop_type, pid, pname_fid, dvids, coeffs,
                              p_min=None, p_max=1e20, c0=0.0, validate=True, comment='')

        model.cross_reference()
        save_load_deck(model, xref='standard', punch=True)

    def test_shell_dvprel(self):
        """tests CBAR and DVPREL1"""
        log = get_logger(level='warning')
        model = BDF(log=log)
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [1., 0., 0.])
        model.add_grid(3, [1., 1., 0.])
        model.add_grid(4, [0., 1., 0.])

        mid = 10
        E = 3.0e7
        G = None
        nu = 0.3
        model.add_mat1(mid, E, G, nu)
        #---------------------------------------------
        eid = 2
        pid = 2
        nids = [1, 2, 3, 4]
        model.add_cquad4(eid, pid, nids, theta_mcid=0.0, zoffset=0., tflag=0,
                         T1=None, T2=None, T3=None, T4=None, comment='')
        nids = [1, 2, 3]
        model.add_ctria3(eid+1, pid+1, nids, zoffset=0., theta_mcid=0.0, tflag=0,
                         T1=None, T2=None, T3=None, comment='')

        model.add_pshell(pid, mid1=mid, t=0.1, mid2=mid, twelveIt3=1.0, mid3=mid,
                         tst=0.833333, nsm=0.0, z1=None, z2=None, mid4=None, comment='')

        mids = [mid] * 3
        thicknesses = [0.1] * 3
        model.add_pcomp(pid+1, mids, thicknesses, thetas=None, souts=None,
                        nsm=0., sb=0., ft=None, tref=0., ge=0., lam=None, z0=None, comment='')
        #---------------------
        dvprels = [
            #oid, pid, prop_type, pname_fid
            (1, pid, 'PSHELL', 'T'),

            (2, pid+1, 'PCOMP', 'T1'),
            (3, pid+1, 'PCOMP', 'T2'),
            (4, pid+1, 'PCOMP', 'T3'),

            (5, pid+1, 'PCOMP', 'THETA1'),
            (6, pid+1, 'PCOMP', 'THETA2'),
            (7, pid+1, 'PCOMP', 'THETA3'),
            (8, pid+1, 'PCOMP', 'SB'),
            (9, pid+1, 'PCOMP', 'GE'),
            (10, pid+1, 'PCOMP', 'Z0'),
        ]
        coeffs = [1.0]
        xinit = 1.0
        for dvprel_data in dvprels:
            oid, pid, prop_type, pname_fid = dvprel_data
            dvids = [oid]
            desvar_id = oid

            dim_label = pname_fid
            prop_label = prop_type
            if pname_fid.startswith('THETA'):
                dim_label = 'TH%s' % pname_fid[-1]
            if prop_type == 'PCOMP':
                prop_label = 'COMP'
            #elif prop_type == 'PBEAM':
                #prop_label = 'PBM'
            label = '%s%s%s' % (prop_label, dim_label, pid)

            model.add_desvar(desvar_id, label, xinit,
                             xlb=-1e20, xub=1e20, delx=None, ddval=None, comment='')
            model.add_dvprel1(oid, prop_type, pid, pname_fid, dvids, coeffs,
                              p_min=None, p_max=1e20, c0=0.0, validate=True, comment='')

        model.cross_reference()
        save_load_deck(model, xref='standard', punch=True)

    def test_dtable(self):
        """
        tests:
         - DTABLE
         - DRESP2

        """
        log = get_logger(level='warning')
        model = BDF(log=log)
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [1., 0., 0.])
        model.add_grid(3, [1., 1., 0.])
        model.add_grid(4, [0., 1., 0.])
        default_values = {
            'A' : 1.0,
            'B1': 2.0,
            'C': -3.,
        }
        dtable = model.add_dtable(default_values, comment='table')
        str(dtable)
        #print(dtable)
        dresp_id = 42
        label = 'cat'
        dequation = 1000
        region = None
        params = {
            (0, 'DTABLE'): ['B1', 'C', 'A'],
        }
        dresp2 = model.add_dresp2(dresp_id, label, dequation, region, params,
                                  method='MIN', c1=1., c2=0.005, c3=10., validate=True,
                                  comment='')
        eqs = [
            'f(a,b,c)=a+b+c'
        ]
        model.add_deqatn(dequation, eqs)
        fields = dresp2.raw_fields()
        #print(fields)
        expected_fields = ['DRESP2', 42, 'cat', 1000, None, 'MIN', 1.0, 0.005, 10.0,
                           'DTABLE', 'B1', 'C', 'A', None, None, None, None]
        assert fields == expected_fields
        #DRESP2        42     cat    1000
        #          DTABLE      B1       C       A
        #---------------------------------
        dresp2.cross_reference(model)
        str(dresp2)
        fields = dresp2.raw_fields()
        assert fields == expected_fields, fields

if __name__ == '__main__':  # pragma: no cover
    unittest.main()
