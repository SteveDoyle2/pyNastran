"""tests static load cards"""
import os
import unittest
import numpy as np
from numpy import array, allclose, array_equal, set_printoptions
from cpylog import get_logger
set_printoptions(suppress=True, precision=3)

import pyNastran
from pyNastran.bdf.bdf import BDF, BDFCard, DAREA, PLOAD4, read_bdf, CaseControlDeck
from pyNastran.bdf.cards.base_card import expand_thru_by
from pyNastran.bdf.cards.collpase_card import collapse_thru_by
from pyNastran.bdf.cards.test.utils import save_load_deck
from pyNastran.bdf.mesh_utils.loads import sum_forces_moments, sum_forces_moments_elements
#from pyNastran.bdf.errors import DuplicateIDsError

from pyNastran.op2.op2 import read_op2 # OP2,

bdf = BDF(debug=False)
TEST_PATH = pyNastran.__path__[0]
MODEL_PATH = os.path.join(TEST_PATH, '..', 'models')
log = get_logger(level='warning')


class TestLoads(unittest.TestCase):
    def test_force(self):
        """tests CONROD, FORCE"""
        model = BDF(log=log)
        eid = 1
        mid = 100
        nids = [10, 11]
        A = 3.14
        model.add_conrod(eid, mid, nids, A, j=0.0, c=0.0, nsm=0.0,
                         comment='')
        model.add_grid(10, [10., 0., 0.])
        model.add_grid(11, [11., 0., 0.])
        E = 3.0e7
        G = None
        nu = 0.3
        model.add_mat1(mid, E, G, nu)

        sid = 10000
        node = 11
        mag = 42.
        xyz = [1., 1., 2.]
        force = model.add_force(sid, node, mag, xyz)
        force.raw_fields()
        model.validate()
        model.pop_parse_errors()
        assert np.array_equal(force.scaled_vector, np.array([42., 42., 84.])), force.scaled_vector
        model.cross_reference()
        force.raw_fields()
        save_load_deck(model)

    def test_moment(self):
        """tests CONROD, MOMENT"""
        model = BDF(log=log)
        eid = 1
        mid = 100
        nids = [10, 11]
        A = 3.14
        model.add_conrod(eid, mid, nids, A, j=0.0, c=0.0, nsm=0.0,
                         comment='')
        model.add_grid(10, [10., 0., 0.])
        model.add_grid(11, [11., 0., 0.])
        E = 3.0e7
        G = None
        nu = 0.3
        model.add_mat1(mid, E, G, nu)

        sid = 10000
        node = 11
        mag = 42.
        xyz = [1., 1., 2.]
        moment = model.add_moment(sid, node, mag, xyz)
        moment.raw_fields()
        model.validate()
        model.pop_parse_errors()
        assert np.array_equal(moment.scaled_vector, np.array([42., 42., 84.])), moment.scaled_vector
        model.cross_reference()
        moment.raw_fields()
        save_load_deck(model)

    def test_accel1(self):
        """tests ACCEL1"""
        model = BDF(log=log)
        sid = 42
        N = [0., 0., 1.]
        nodes = [10, 11]
        scale = 3.14
        accel1 = model.add_accel1(sid, scale, N, nodes, cid=0, comment='accel1')
        accel1.raw_fields()
        accel1.write_card(size=8)
        accel1.write_card(size=16)
        accel1.write_card(size=16, is_double=True)

        model.add_grid(10, [10., 0., 0.])
        model.add_grid(11, [11., 0., 0.])
        model.validate()
        model.pop_parse_errors()
        model.cross_reference()
        model.pop_xref_errors()

        accel1.raw_fields()
        accel1.write_card(size=8)
        accel1.write_card(size=16)
        accel1.write_card(size=16, is_double=True)
        save_load_deck(model)

    def test_accel1_2(self):
        """tests problematic ACCEL1 cards"""
        cards = [
            ['ACCEL1  1               -32.2   0.2672610.5345220.801784',
             '+       1       2       3       4       5'],
            ['ACCEL1  2               -64.4   0.2672610.5345220.801784',
             '+       6       THRU    9       10'],
            ['ACCEL1  3               -96.6   0.2672610.5345220.801784',
             '+       11      12      THRU    15'],
            ['ACCEL1  4               -128.8  0.2672610.5345220.801784',
             '+       1       THRU    10      BY      2       12      THRU    24',
             '+       BY      2'],
            ['ACCEL1  5               -161.0  0.2672610.5345220.801784',
             '+       14      THRU    24      BY      3'],
        ]
        fields = ['14', 'THRU', '24', 'BY', '2']
        fields = expand_thru_by(fields, set_fields=True, sort_fields=True)
        assert fields == [14, 16, 18, 20, 22, 24], 'fields=%s' % fields
        assert collapse_thru_by(fields) == [14, 'THRU', 24, 'BY', 2], collapse_thru_by(fields)

        fields = ['2', 'THRU', '5', 'BY', '1', '10']
        fields = expand_thru_by(fields, set_fields=True, sort_fields=True)
        assert fields == [2, 3, 4, 5, 10], 'fields=%s' % fields
        assert collapse_thru_by(fields) == [2, 'THRU', 5, 10], collapse_thru_by(fields)

        fields = ['14', 'THRU', '24', 'BY', '3']
        fields = expand_thru_by(fields, set_fields=True, sort_fields=True)
        assert fields == [14, 17, 20, 23, 24], 'fields=%s' % fields
        # [14, 'THRU', 24, 'BY', 3] - this is the ideal answer, but close enough...
        assert collapse_thru_by(fields) == [14, 17, 20, 23, 24], collapse_thru_by(fields)

        fields = ['14', 'THRU', '24', 'BY', '2']
        fields = expand_thru_by(fields, set_fields=True, sort_fields=True)
        assert collapse_thru_by(fields) == [14, 'THRU', 24, 'BY', 2], collapse_thru_by(fields)

        model = BDF(log=log)
        for nid in [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 20, 22, 23, 24]:
            model.add_grid(nid, [0., 0., 0.])

        for card_lines in cards:
            model.add_card(card_lines, 'ACCEL1', comment='',
                           is_list=False, has_none=True)

        for unused_key, loads in sorted(model.loads.items()):
            for load in loads:
                str(load)
        save_load_deck(model)

    def test_accel(self):
        """tests ACCEL"""
        model = BDF(log=log)
        sid = 42
        N = [0., 0., 1.]
        #nodes = [10, 11]
        #scale = 3.14
        direction = 'Z'
        locs = [11., 22., 33.]
        vals = [1., 2., 3.]
        accel = model.add_accel(sid, N, direction, locs, vals, cid=0,
                                comment='accel')
        accel.raw_fields()
        accel.write_card(size=8)
        accel.write_card(size=16)
        accel.write_card(size=16, is_double=True)

        model.add_grid(10, [10., 0., 0.])
        model.add_grid(11, [11., 0., 0.])
        model.validate()
        model.pop_parse_errors()
        model.cross_reference()
        model.pop_xref_errors()

        accel.raw_fields()
        accel.write_card(size=8)
        accel.write_card(size=16)
        accel.write_card(size=16, is_double=True)
        save_load_deck(model)

    def test_darea_01(self):
        """tests a DAREA"""
        #DAREA SID P1 C1 A1  P2 C2 A2
        #DAREA 3   6   2 8.2 15 1  10.1
        lines = ['DAREA,3,6,2,8.2,15,1,10.1']
        card = bdf._process_card(lines)
        cardi = BDFCard(card)

        size = 8
        card = DAREA.add_card(cardi)
        card.write_card(size, 'dummy')
        card.raw_fields()

    def test_grav(self):
        """tests a GRAV"""
        model = BDF(log=log)
        sid = 1
        scale = 1.0
        N = [1., 2., 3.]
        grav = model.add_grav(sid, scale, N, cid=0, mb=0, comment='grav')
        grav.raw_fields()
        grav.write_card(size=8)
        grav.write_card(size=16)
        grav.write_card(size=16, is_double=True)

        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [1., 0., 0.])
        model.add_grid(3, [1., 1., 0.])
        model.add_grid(4, [0., 1., 0.])

        eid = 4
        pid = 4
        nids = [1, 2, 3, 4]
        model.add_cquad4(eid, pid, nids, theta_mcid=0.0, zoffset=0.,
                         tflag=0, T1=1.0, T2=1.0, T3=1.0, T4=1.0, comment='')

        mid = 10
        # mass = (rho*t + nsm)*A = (0.2*0.5 + 0.3) * 0.5 = 0.4 * 0.5 = 0.2
        model.add_pshell(pid, mid1=None, t=0.5, mid2=mid, twelveIt3=1.0,
                         mid3=None, tst=0.833333,
                         nsm=0.3, z1=None, z2=None,
                         mid4=None, comment='')

        E = 3.0e7
        G = None
        nu = 0.3
        model.add_mat1(mid, E, G, nu, rho=0.2, a=0.0, tref=0.0, ge=0.0,
                       St=0.0, Sc=0.0, Ss=0.0, mcsid=0,
                       comment='')

        model.validate()
        model.cross_reference()

        eids = [4]
        p0 = [0., 0., 0.]
        loadcase_id = sid
        forces1, moments1 = sum_forces_moments(model, p0, loadcase_id,
                                               include_grav=False, xyz_cid0=None)
        forces2, moments2 = sum_forces_moments_elements(model, p0, loadcase_id, eids, nids,
                                                        include_grav=False, xyz_cid0=None)
        assert np.array_equal(forces1, forces2)
        assert np.array_equal(moments1, moments2)

        forces1, moments1 = model.sum_forces_moments(p0, loadcase_id,
                                                     include_grav=True, xyz_cid0=None)
        forces2, moments2 = sum_forces_moments_elements(model, p0, loadcase_id, eids, nids,
                                                        include_grav=True, xyz_cid0=None)
        assert np.array_equal(forces1, forces2)
        assert np.array_equal(moments1, moments2)
        save_load_deck(model)

    def test_gmload(self):
        """tests GMLOAD"""
        model = BDF(log=log)
        sid = 1
        normal = [1., 2., 3.]
        entity = 'cat'
        entity_id = 37
        method = 'frog'
        load_magnitudes = [31.]
        gmload = model.add_gmload(sid, normal, entity, entity_id, method,
                                  load_magnitudes, cid=0,
                                  comment='gmload')
        gmload.raw_fields()
        model.validate()
        model.cross_reference()
        save_load_deck(model, run_convert=False)

    def test_pload4_01(self):
        """tests a PLOAD4"""
        lines = ['PLOAD4  1000    1       -60.    -60.    60.             1']
        card = bdf._process_card(lines)
        cardi = BDFCard(card)

        size = 8
        card = PLOAD4.add_card(cardi)
        card.write_card(size, 'dummy')
        card.raw_fields()

    def test_pload4_02(self):
        """tests a PLOAD4"""
        lines = ['PLOAD4  1       101     1.                              10000   10011']
        card = bdf._process_card(lines)
        cardi = BDFCard(card)

        size = 8
        card = PLOAD4.add_card(cardi)
        card.write_card(size, 'dummy')
        card.raw_fields()

    def test_pload4_line(self):
        """tests a PLOAD4 LINE option"""
        #PLOAD4        10      10      0.819.2319
        #0      1.      0.      0.    LINE    NORM
        model = BDF(log=log, mode='msc')

        sid = 1
        eids = 1
        pressures = 1.
        dummy = model.add_pload4(sid, eids, pressures,
                                 g1=None, g34=None, cid=0, nvector=None,
                                 surf_or_line='SURFBAD', line_load_dir='NORMBAD', comment='')
        with self.assertRaises(RuntimeError):
            dummy.validate()
        dummy.surf_or_line = 'SURF'
        with self.assertRaises(RuntimeError):
            dummy.validate()
        dummy.line_load_dir = 'NORM'
        dummy.validate()
        model.clear_attributes()

        eid = 10
        pid = 20
        mid = 100
        nids = [1, 2, 3, 4]
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [1., 0., 0.])
        model.add_grid(3, [1., 1., 0.])
        model.add_grid(4, [0., 1., 0.])
        model.add_pshell(pid, mid1=mid, t=0.1, mid2=None, twelveIt3=1.0,
                         mid3=None, tst=0.833333, nsm=0.0, z1=None, z2=None, mid4=None, comment='')

        E = 3.0e7
        G = None
        nu = 0.3
        model.add_mat1(mid, E, G, nu,
                       rho=0.0, a=0.0, tref=0.0, ge=0.0,
                       St=0.0, Sc=0.0, Ss=0.0, mcsid=0, comment='')
        model.add_cquadr(eid, pid, nids,
                         theta_mcid=0.0, zoffset=0., tflag=0,
                         T1=None, T2=None, T3=None, T4=None, comment='')
        model.add_ctriar(eid+1, pid, nids[:-1],
                         theta_mcid=0.0, zoffset=0., tflag=0,
                         T1=None, T2=None, T3=None, comment='')

        # The SORL field is ignored by all elements except QUADR and TRIAR.
        #    For QUADR or TRIAR only, if SORL=LINE, the consistent edge loads
        #    are defined by the PLOAD4 entry. P1, P2, P3 and  P4 are load per
        #    unit length at the corner of the element.
        #
        # All four Ps are given:
        #    The line loads along all four edges of the element are defined.
        #
        # If any P is blank:
        #    The line loads for only two edges are defined. For example,
        #    if P1 is blank, the line loads of the two edges connecting to G1 are zero.
        #
        # If two Ps are given:
        #    The line load of the edge connecting to the two grid points is defined.
        #
        # If only one P is given:
        #    The second P value default to the first P value.  For example,
        #    P1 denotes that the line load along edge G1 and G2 has the
        #    constant value of P1.
        #
        sid = 10
        eids = [10, 11]
        pressures = [1., 0., 0., 0.]
        cid = 0

        # The direction of the line load (SORL=LINE) is defined by either (CID, N1, N2, N3) or
        # LDIR. Fatal error will be issued if both methods are given.  TANG denotes that the line
        # load is in tangential direction of the edge, pointing from G1 to G2 if the edge is
        # connecting G1 and G2.
        #
        # NORM denotes that the line load is in the mean plan, normal to the edge, and pointing
        # outward from the element.  X, Y, or Z denotes the line load is in the X, Y, or Z
        # direction of the element coordinate system.  If both (CID, N1, n2, N3) and LDIR are
        # blank, then the default is LDIR=NORM.
        nvector = [1., 0., 0.]
        pload4 = model.add_pload4(sid, eids, pressures, g1=None, g34=None,
                                  cid=cid, nvector=nvector,
                                  surf_or_line='LINE', line_load_dir='NORM', comment='pload4_line')
        assert pload4.raw_fields() == ['PLOAD4', 10, 10, 1.0, 0.0, 0.0, 0.0, 'THRU', 11,
                                       0, 1.0, 0.0, 0.0, 'LINE', 'NORM']
        str(pload4)

        unused_pload4_surf = model.add_pload4(sid, eids, pressures, g1=None, g34=None,
                                              cid=cid, nvector=nvector,
                                              surf_or_line='SURF', line_load_dir='NORM',
                                              comment='pload4_line')
        str(pload4.raw_fields())

        sid = 11
        eids = 10
        pressures = 1.0
        unused_pload4_surf = model.add_pload4(sid, eids, pressures, g1=None, g34=None,
                                              cid=cid, nvector=nvector,
                                              surf_or_line='SURF', line_load_dir='NORM',
                                              comment='pload4_line')
        model.validate()
        model.cross_reference()

        p0 = [0., 0., 0.]
        loadcase_id = sid
        sum_forces_moments(model, p0, loadcase_id, include_grav=False, xyz_cid0=None)

        eids = None
        nids = None
        sum_forces_moments_elements(model, p0, loadcase_id, eids, nids,
                                    include_grav=False, xyz_cid0=None)
        save_load_deck(model)

    def test_pload4_cpenta(self):
        """tests a PLOAD4 with a CPENTA"""
        bdf_filename = os.path.join(MODEL_PATH, 'unit', 'pload4', 'cpenta.bdf')
        op2_filename = os.path.join(MODEL_PATH, 'unit', 'pload4', 'cpenta.op2')
        op2 = read_op2(op2_filename, log=log)

        model = read_bdf(bdf_filename, log=log)
        # p0 = (model.nodes[21].xyz + model.nodes[22].xyz + model.nodes[23].xyz) / 3.
        p0 = model.nodes[21].xyz
        angles = [
            (23, 24), (24, 23),
            (21, 26), (26, 21),
        ]
        nx = [
            (23, 25), (25, 23),
            (22, 26), (26, 22),
        ]

        msg = ''
        for isubcase, subcase in sorted(model.subcases.items()):
            if isubcase == 0:
                continue
            #if isubcase != 17:
                #continue
            loadcase_id = subcase.get_parameter('LOAD')[0]
            load = model.Load(loadcase_id, consider_load_combinations=True)[0]
            elem = load.eids_ref[0]
            g1 = load.g1_ref.nid
            if load.g34_ref is None:
                g34 = None
                #print(load)
                face, area, centroid, normal = elem.get_face_area_centroid_normal(g1)
                assert area == 0.5, area
                if g1 in [21, 22, 23]:
                    assert face == (2, 1, 0), 'g1=%s face=%s' % (g1, face)
                    assert array_equal(centroid, array([2/3., 1/3., 0.])), 'fore g1=%s g34=%s face=%s centroid=%s\n%s' % (g1, g34, face, centroid, msg)
                    assert array_equal(normal, array([0., 0., 1.])), 'fore g1=%s g34=%s face=%s normal=%s\n%s' % (g1, g34, face, normal, msg)
                else:
                    assert face == (3, 4, 5), 'g1=%s face=%s' % (g1, face)
                    assert array_equal(centroid, array([2/3., 1/3., 2.])), 'aft g1=%s g34=%s face=%s centroid=%s\n%s' % (g1, g34, face, centroid, msg)
                    assert array_equal(normal, array([0., 0., -1.])), 'aft g1=%s g34=%s face=%s normal=%s\n%s' % (g1, g34, face, normal, msg)
            else:
                g34 = load.g34_ref.nid
                face, area, centroid, normal = elem.get_face_area_centroid_normal(g1, g34)
                if (g1, g34) in angles:
                    self.assertAlmostEqual(area, 2 * 2**0.5, msg='g1=%s g34=%s face=%s area=%s' % (g1, g34, face, area))
                elif (g1, g34) in nx:
                    self.assertEqual(area, 2.0, 'area=%s' % area)
                    msg = '%s%s%s%s\n' % (
                        elem.nodes[face[0]], elem.nodes[face[1]], elem.nodes[face[2]], elem.nodes[face[3]])
                    assert array_equal(centroid, array([1., .5, 1.])), 'Nx g1=%s g34=%s face=%s centroid=%g\n%s' % (g1, g34, face, centroid, msg)
                    assert array_equal(normal, array([-1., 0., 0.])), 'Nx g1=%s g34=%s face=%s normal=%g\n%s' % (g1, g34, face, normal, msg)
                else:
                    msg = '%s%s%s%s\n' % (
                        elem.nodes[face[0]], elem.nodes[face[1]], elem.nodes[face[2]], elem.nodes[face[3]])

                    assert array_equal(centroid, array([0.5, .0, 1.])), 'Ny g1=%s g34=%s face=%s centroid=%s\n%s' % (g1, g34, face, centroid, msg)
                    assert array_equal(normal, array([0., 1., 0.])), 'Ny g1=%s g34=%s face=%s normal=%s\n%s' % (g1, g34, face, normal, msg)
                    self.assertEqual(area, 2.0, 'area=%s' % area)

            forces1, moments1 = sum_forces_moments(model, p0, loadcase_id, include_grav=False)
            eids = None
            nids = None
            forces2, moments2 = sum_forces_moments_elements(
                model, p0, loadcase_id, eids, nids, include_grav=False)
            assert allclose(forces1, forces2), 'forces1=%s forces2=%s' % (forces1, forces2)
            assert allclose(moments1, moments2), 'moments1=%s moments2=%s' % (moments1, moments2)

            case = op2.spc_forces[isubcase]
            fm = -case.data[0, :3, :].sum(axis=0)
            assert len(fm) == 6, fm
            if not allclose(forces1[0], fm[0]):
                model.log.error('subcase=%-2i Fx f=%s fexpected=%s face=%s' % (
                    isubcase, forces1.tolist(), fm.tolist(), face))
            if not allclose(forces1[1], fm[1]):
                model.log.error('subcase=%-2i Fy f=%s fexpected=%s face=%s' % (
                    isubcase, forces1.tolist(), fm.tolist(), face))
            if not allclose(forces1[2], fm[2]):
                model.log.error('subcase=%-2i Fz f=%s fexpected=%s face=%s' % (
                    isubcase, forces1.tolist(), fm.tolist(), face))
            # if not allclose(moments1[0], fm[3]):
                # print('%i Mx m=%s fexpected=%s' % (isubcase, moments1, fm))
            # if not allclose(moments1[1], fm[4]):
                # print('%i My m=%s fexpected=%s' % (isubcase, moments1, fm))
            # if not allclose(moments1[2], fm[5]):
                # print('%i Mz m=%s fexpected=%s' % (isubcase, moments1, fm))

            #self.assertEqual(forces1[0], fm[0], 'f=%s fexpected=%s' % (forces1, fm[:3]))
            #self.assertEqual(forces1[1], fm[1], 'f=%s fexpected=%s' % (forces1, fm[:3]))
            #self.assertEqual(forces1[2], fm[2], 'f=%s fexpected=%s' % (forces1, fm[:3]))
            #self.assertEqual(moments1[0], fm[3], 'm=%s mexpected=%s' % (moments1, fm[3:]))
            #self.assertEqual(moments1[1], fm[4], 'm=%s mexpected=%s' % (moments1, fm[3:]))
            #self.assertEqual(moments1[2], fm[5], 'm=%s mexpected=%s' % (moments1, fm[3:]))
        save_load_deck(model, run_loads=False)

    def test_pload4_ctria3(self):
        """tests a PLOAD4 with a CTRIA3"""
        bdf_filename = os.path.join(MODEL_PATH, 'unit', 'pload4', 'ctria3.bdf')
        op2_filename = os.path.join(MODEL_PATH, 'unit', 'pload4', 'ctria3.op2')
        op2 = read_op2(op2_filename, log=log)

        model = read_bdf(bdf_filename, log=log)
        # p0 = (model.nodes[21].xyz + model.nodes[22].xyz + model.nodes[23].xyz) / 3.
        p0 = model.nodes[21].xyz

        subcase_ids = [1, 2, 4, 5, 6] #  no 3, 7, 8
        for isubcase in subcase_ids:
            subcase = model.subcases[isubcase]
            loadcase_id = subcase.get_parameter('LOAD')[0]
            loads = model.get_reduced_loads(loadcase_id, consider_load_combinations=True)[0]
            load = loads[0]
            assert len(loads) == 1, 'subcase:\n%s\nloads=\n%s' % (subcase, loads)
            elem = load.eids_ref[0]

            centroid = elem.Centroid()
            normal = elem.Normal()

            msg = '%s%s%s\n' % (
                elem.nodes[0], elem.nodes[1], elem.nodes[2])

            assert array_equal(centroid, array([2/3., 1/3., 0.])), 'centroid=%s\n%s' % (centroid, msg)
            assert array_equal(normal, array([0., 0., 1.])), 'normal=%s\n%s' % (normal, msg)

            forces1, moments1 = sum_forces_moments(model, p0, loadcase_id, include_grav=False)
            eids = None
            nids = None
            forces2, moments2 = sum_forces_moments_elements(
                model, p0, loadcase_id, eids, nids, include_grav=False)
            assert allclose(forces1, forces2), 'forces1=%s forces2=%s' % (forces1, forces2)
            assert allclose(moments1, moments2), 'moments1=%s moments2=%s' % (moments1, moments2)

            case = op2.spc_forces[isubcase]
            fm = -case.data[0, :, :].sum(axis=0)
            assert len(fm) == 6, fm
            if not allclose(forces1[0], fm[0]):
                model.log.error('subcase=%-2i Fx f=%s fm_expected=%s' % (
                    isubcase, forces1.tolist(), fm.tolist()))
            if not allclose(forces1[1], fm[1]):
                model.log.error('subcase=%-2i Fy f=%s fm_expected=%s' % (
                    isubcase, forces1.tolist(), fm.tolist()))
            if not allclose(forces1[2], fm[2]):
                model.log.error('subcase=%-2i Fz f=%s fm_expected=%s' % (
                    isubcase, forces1.tolist(), fm.tolist()))
        save_load_deck(model, punch=False)

    def test_pload4_cquad4(self):
        """tests a PLOAD4 with a CQUAD4"""
        bdf_filename = os.path.join(MODEL_PATH, 'unit', 'pload4', 'cquad4.bdf')
        op2_filename = os.path.join(MODEL_PATH, 'unit', 'pload4', 'cquad4.op2')
        op2 = read_op2(op2_filename, log=log)

        model = read_bdf(bdf_filename, log=log)
        # p0 = (model.nodes[21].xyz + model.nodes[22].xyz + model.nodes[23].xyz) / 3.
        p0 = model.nodes[21].xyz

        eids = None
        nids = None
        subcase_ids = [1, 2, 3, 4, 5, 6, 7, 8]
        for isubcase in subcase_ids:
            subcase = model.subcases[isubcase]

            loadcase_id = subcase.get_parameter('LOAD')[0]
            load = model.Load(loadcase_id)
            loadi = load[0]
            if loadi.type == 'PLOAD4':
                elem = loadi.eids_ref[0]
                #area = 1.0
                centroid = elem.Centroid()
                normal = elem.Normal()
                # centroid = [0.5, 0.5, 0.]
                # normal   = [0., 0., 1.]
                #print('centroid=%s normal=%s' % (centroid, normal))
                msg = '%s%s%s\n' % (elem.nodes[0], elem.nodes[1], elem.nodes[2])

                assert array_equal(centroid, array([0.5, 0.5, 0.])), 'centroid=%s\n%s' % (centroid, msg)
                assert array_equal(normal, array([0., 0., 1.])), 'normal=%s\n%s' % (normal, msg)

            f1, m1 = sum_forces_moments(model, p0, loadcase_id, include_grav=False)
            f2, m2 = sum_forces_moments_elements(model, p0, loadcase_id, eids, nids,
                                                 include_grav=False)
            assert allclose(f1, f2), 'f1=%s f2=%s' % (f1, f2)
            assert allclose(m1, m2), 'm1=%s m2=%s' % (m1, m2)

            case = op2.spc_forces[isubcase]
            fm = -case.data[0, :, :].sum(axis=0)
            assert len(fm) == 6, fm
            force = fm[:3]
            if not allclose(f1[0], force[0]):
                model.log.error('subcase=%-2i Fx f=%s force_expected=%s' % (
                    isubcase, f1.tolist(), force.tolist()))
            if not allclose(f1[1], force[1]):
                model.log.error('subcase=%-2i Fy f=%s force_expected=%s' % (
                    isubcase, f1.tolist(), force.tolist()))
            if not allclose(f1[2], force[2]):
                model.log.error('subcase=%-2i Fz f=%s force_expected=%s' % (
                    isubcase, f1.tolist(), force.tolist()))
        save_load_deck(model, punch=False)

    def test_pload4_ctetra(self):
        """tests a PLOAD4 with a CTETRA"""
        bdf_filename = os.path.join(MODEL_PATH, 'unit', 'pload4', 'ctetra.bdf')
        op2_filename = os.path.join(MODEL_PATH, 'unit', 'pload4', 'ctetra.op2')
        op2 = read_op2(op2_filename, log=log)

        model = read_bdf(bdf_filename, log=log)
        p0 = model.nodes[21].xyz

        nx_plus = [ # 21, 24, 23
            (21, 22), (24, 22), (23, 22), # (22, 23),
            #(21, 22), (23, 22),
        ]
        ny_plus = [ # 21, 22, 24
            #(21, 23), (22, 23), (24, 23),
            #(23, 21), (23, 22), (23, 24),
        ]
        nz_plus = [ # 21, 22, 23
            (21, 24), (22, 24), (23, 24),
            #(24, 21), (24, 22), (24, 23),
        ]

        for isubcase, subcase in sorted(model.subcases.items()):
            if isubcase == 0:
                continue
            loadcase_id = subcase.get_parameter('LOAD')[0]
            load = model.loads[loadcase_id][0]
            elem = load.eids_ref[0]
            g1 = load.g1_ref.nid

            # f, m = model.sum_forces_moments(p0, loadcase_id, include_grav=False)
            # case = op2.spc_forces[isubcase]
            # fm = case.data[0, 0, :]#.ravel()
            # if f[0] != fm[0]:
                # print('%i f=%s fexpected=%s' % (isubcase, f, fm))

            g34 = load.g34_ref.nid
            face, area, unused_centroid, normal = elem.get_face_area_centroid_normal(g1, g34)
            msg = '%s%s%s\n' % (
                elem.nodes[face[0]], elem.nodes[face[1]],
                elem.nodes[face[2]])

            if (g1, g34) in nx_plus:
                self.assertEqual(area, 0.5, '+Nx area=%s\n%s' % (area, msg))
                #assert array_equal(centroid, array([1., .5, 1.])), '-Nx g1=%s g34=%s face=%s centroid=%g\n%s' % (g1, g34, face, centroid, msg)
                assert array_equal(normal, array([1., 0., 0.])), '+Nx g1=%s g34=%s face=%s normal=%g\n%s' % (g1, g34, face, normal, msg)

            elif (g1, g34) in ny_plus:
                self.assertEqual(area, 0.5, '+Ny area=%s\n%s' % (area, msg))
                #assert array_equal(centroid, array([1., .5, 1.])), '-Nz g1=%s g34=%s face=%s centroid=%g\n%s' % (g1, g34, face, centroid, msg)
                assert array_equal(normal, array([0., 1., 0.])), '+Ny g1=%s g34=%s face=%s normal=%s\n%s' % (g1, g34, face, normal, msg)
            elif (g1, g34) in nz_plus:
                self.assertEqual(area, 0.5, '+Nz area=%s\n%s' % (area, msg))
                #assert array_equal(centroid, array([1., .5, 1.])), '-Nz g1=%s g34=%s face=%s centroid=%g\n%s' % (g1, g34, face, centroid, msg)
                assert array_equal(normal, array([0., 0., 1.])), '+Nz g1=%s g34=%s face=%s normal=%s\n%s' % (g1, g34, face, normal, msg)
                #assert array_equal(centroid, array([1., .5, 1.])),  'Nx g1=%s g34=%s face=%s centroid=%s\n%s' % (g1, g34, face, centroid, msg)
                #assert array_equal(normal, array([-1., 0., 0.])),  'Nx g1=%s g34=%s face=%s normal=%s\n%s' % (g1, g34, face, normal, msg)
            else:
                self.assertEqual(area, 0.75**0.5, 'slant g1=%s g34=%s face=%s area=%s\n%s' % (g1, g34, face, area, msg))
                #assert array_equal(centroid, array([1., .5, 1.])), '-Nz g1=%s g34=%s face=%s centroid=%g\n%s' % (g1, g34, face, centroid, msg)
                normal_expected = array([-0.57735027, -0.57735027, -0.57735027])
                diff = normal - normal_expected
                assert allclose(normal, normal_expected), 'slant g1=%s g34=%s face=%s normal=%s\ndiff=%s\n%s' % (g1, g34, face, normal, diff, msg)
                #raise RuntimeError('??? g1=%s g34=%s face=%s normal=%s\n%s' % (g1, g34, face, normal, msg))
            #self.assertEqual(f[0], fm[0], 'f=%s fexpected=%s' % (f, fm[:3]))
            #self.assertEqual(f[1], fm[1], 'f=%s fexpected=%s' % (f, fm[:3]))
            #self.assertEqual(f[2], fm[2], 'f=%s fexpected=%s' % (f, fm[:3]))
            #self.assertEqual(m[0], fm[3], 'm=%s mexpected=%s' % (m, fm[3:]))
            #self.assertEqual(m[1], fm[4], 'm=%s mexpected=%s' % (m, fm[3:]))
            #self.assertEqual(m[2], fm[5], 'm=%s mexpected=%s' % (m, fm[3:]))
            forces1, moments1 = sum_forces_moments(model, p0, loadcase_id, include_grav=False)
            eids = None
            nids = None
            forces2, moments2 = sum_forces_moments_elements(
                model, p0, loadcase_id, eids, nids, include_grav=False)
            assert allclose(forces1, forces2), 'forces1=%s forces2=%s' % (forces1, forces2)
            assert allclose(moments1, moments2), 'moments1=%s moments2=%s' % (moments1, moments2)

            case = op2.spc_forces[isubcase]
            fm = -case.data[0, :, :].sum(axis=0)
            assert len(fm) == 6, fm
            if not allclose(forces1[0], fm[0]):
                model.log.error('subcase=%-2i Fx g=(%s,%s) forces1=%s fexpected=%s '
                                'face=%s normal=%s' % (
                                    isubcase, g1, g34, forces1, fm, face, normal))
            if not allclose(forces1[1], fm[1]):
                model.log.error('subcase=%-2i Fy g=(%s,%s) forces1=%s fexpected=%s '
                                'face=%s normal=%s' % (
                                    isubcase, g1, g34, forces1, fm, face, normal))
            if not allclose(forces1[2], fm[2]):
                model.log.error('subcase=%-2i Fz g=(%s,%s) forces1=%s fexpected=%s '
                                'face=%s normal=%s' % (
                                    isubcase, g1, g34, forces1, fm, face, normal))
        save_load_deck(model, punch=False, run_loads=False)

    def test_pload4_chexa(self):
        """tests a PLOAD4 with a CHEXA"""
        bdf_filename = os.path.join(MODEL_PATH, 'unit', 'pload4', 'chexa.bdf')
        op2_filename = os.path.join(MODEL_PATH, 'unit', 'pload4', 'chexa.op2')
        op2 = read_op2(op2_filename, log=log)

        model = read_bdf(bdf_filename, log=log)
        p0 = model.nodes[21].xyz
        nx_minus = [
            (22, 27), (27, 22),
            (23, 26), (26, 23),
        ]
        nx_plus = [
            (24, 25), (25, 24),
            (21, 28), (28, 21),
            #(23, 25), (25, 23),
            #(22, 26), (26, 22),
        ]

        ny_minus = [
            (24, 27), (27, 24),
            (23, 28), (28, 23),
        ]
        ny_plus = [
            (21, 26), (26, 21),
            (22, 25), (25, 22),
        ]

        nz_minus = [
            (25, 27), (27, 25),
            (26, 28), (28, 26),
        ]
        nz_plus = [
            (21, 23), (23, 21),
            (24, 22), (22, 24),
        ]

        for isubcase, subcase in sorted(model.subcases.items()):
            if isubcase == 0:
                continue
            loadcase_id = subcase.get_parameter('LOAD')[0]
            load = model.loads[loadcase_id][0]
            elem = load.eids_ref[0]
            g1 = load.g1_ref.nid

            # f, m = model.sum_forces_moments(p0, loadcase_id, include_grav=False)
            # case = op2.spc_forces[isubcase]
            # fm = case.data[0, 0, :]#.ravel()
            # if f[0] != fm[0]:
                # print('%i f=%s fexpected=%s' % (isubcase, f, fm))

            g34 = load.g34_ref.nid
            face, area, unused_centroid, normal = elem.get_face_area_centroid_normal(g1, g34)
            msg = '%s%s%s%s\n' % (
                elem.nodes[face[0]], elem.nodes[face[1]],
                elem.nodes[face[2]], elem.nodes[face[3]])

            if (g1, g34) in nx_plus:
                self.assertEqual(area, 2.0, '+Nx area=%s' % area)
                #assert array_equal(centroid, array([1., .5, 1.])), '+Nx g1=%s g34=%s face=%s centroid=%g\n%s' % (g1, g34, face, centroid, msg)
                assert array_equal(normal, array([1., 0., 0.])), '+Nx g1=%s g34=%s face=%s normal=%g\n%s' % (g1, g34, face, normal, msg)
            elif (g1, g34) in nx_minus:
                self.assertEqual(area, 2.0, '-Nx area=%s' % area)
                #assert array_equal(centroid, array([1., .5, 1.])), '-Nx g1=%s g34=%s face=%s centroid=%g\n%s' % (g1, g34, face, centroid, msg)
                assert array_equal(normal, array([-1., 0., 0.])), '-Nx g1=%s g34=%s face=%s normal=%g\n%s' % (g1, g34, face, normal, msg)

            elif (g1, g34) in ny_plus:
                self.assertEqual(area, 2.0, '+Ny area=%s' % area)
                #assert array_equal(centroid, array([1., .5, 1.])), '+Nz g1=%s g34=%s face=%s centroid=%g\n%s' % (g1, g34, face, centroid, msg)
                assert array_equal(normal, array([0., 1., 0.])), '+Ny g1=%s g34=%s face=%s normal=%s\n%s' % (g1, g34, face, normal, msg)
            elif (g1, g34) in ny_minus:
                self.assertEqual(area, 2.0, '-Ny area=%s' % area)
                #assert array_equal(centroid, array([1., .5, 1.])), '-Nz g1=%s g34=%s face=%s centroid=%g\n%s' % (g1, g34, face, centroid, msg)
                assert array_equal(normal, array([0., -1., 0.])), '-Ny g1=%s g34=%s face=%s normal=%s\n%s' % (g1, g34, face, normal, msg)

            elif (g1, g34) in nz_plus:
                self.assertEqual(area, 1.0, '+Nz area=%s' % area)
                #assert array_equal(centroid, array([1., .5, 1.])), '+Nz g1=%s g34=%s face=%s centroid=%g\n%s' % (g1, g34, face, centroid, msg)
                assert array_equal(normal, array([0., 0., 1.])), '+Nz g1=%s g34=%s face=%s normal=%s\n%s' % (g1, g34, face, normal, msg)
            elif (g1, g34) in nz_minus:
                self.assertEqual(area, 1.0, '-Nz area=%s' % area)
                #assert array_equal(centroid, array([1., .5, 1.])), '-Nz g1=%s g34=%s face=%s centroid=%g\n%s' % (g1, g34, face, centroid, msg)
                assert array_equal(normal, array([0., 0., -1.])), '-Nz g1=%s g34=%s face=%s normal=%s\n%s' % (g1, g34, face, normal, msg)
                #assert array_equal(centroid, array([1., .5, 1.])),  'Nx g1=%s g34=%s face=%s centroid=%s\n%s' % (g1, g34, face, centroid, msg)
                #assert array_equal(normal, array([-1., 0., 0.])),  'Nx g1=%s g34=%s face=%s normal=%s\n%s' % (g1, g34, face, normal, msg)
            else:
                msg = '??? g1=%s g34=%s face=%s normal=%s\n%s' % (g1, g34, face, normal, msg)
                raise RuntimeError(msg)
            #self.assertEqual(forces1[0], fm[0], 'forces1=%s fexpected=%s' % (forces1, fm[:3]))
            #self.assertEqual(forces1[1], fm[1], 'forces1=%s fexpected=%s' % (forces1, fm[:3]))
            #self.assertEqual(forces1[2], fm[2], 'forces1=%s fexpected=%s' % (forces1, fm[:3]))
            #self.assertEqual(moments1[0], fm[3], 'moments1=%s mexpected=%s' % (moments1, fm[3:]))
            #self.assertEqual(moments1[1], fm[4], 'moments1=%s mexpected=%s' % (moments1, fm[3:]))
            #self.assertEqual(moments1[2], fm[5], 'moments1=%s mexpected=%s' % (moments1, fm[3:]))
            forces1, moments1 = sum_forces_moments(model, p0, loadcase_id, include_grav=False)
            eids = None
            nids = None
            forces2, moments2 = sum_forces_moments_elements(model, p0, loadcase_id, eids, nids,
                                                            include_grav=False)
            assert allclose(forces1, forces2), 'forces1=%s forces2=%s' % (forces1, forces2)
            assert allclose(moments1, moments2), 'moments1=%s moments2=%s' % (moments1, moments2)

            case = op2.spc_forces[isubcase]
            fm = -case.data[0, :4, :].sum(axis=0)
            assert len(fm) == 6, fm
            if not allclose(forces1[0], fm[0]):
                model.log.error('subcase=%-2i Fx forces1=%s fexpected=%s face=%s' % (
                    isubcase, forces1.tolist(), fm.tolist(), face))
            if not allclose(forces1[1], fm[1]):
                model.log.error('subcase=%-2i Fy forces1=%s fexpected=%s face=%s' % (
                    isubcase, forces1.tolist(), fm.tolist(), face))
            if not allclose(forces1[2], fm[2]):
                model.log.error('subcase=%-2i Fz forces1=%s fexpected=%s face=%s' % (
                    isubcase, forces1.tolist(), fm.tolist(), face))

    #@unittest.expectedFailure
    #def test_pload1_cbar(self):
        #bdf_filename = os.path.join(test_path, '..', 'models', 'pload4', 'pload1.bdf')
        #op2_filename = os.path.join(test_path, '..', 'models', 'pload4', 'pload1.op2')
        #op2 = OP2(debug=False)
        #op2.read_op2(op2_filename)

        #model = BDF(debug=False)
        #model.read_bdf(bdf_filename)
        ## p0 = (model.nodes[21].xyz + model.nodes[22].xyz + model.nodes[23].xyz) / 3.
        #p0 = model.nodes[1].xyz

        #fail = False
        #for isubcase, subcase in sorted(model.subcases.items()):
            #if isubcase == 0:
                #continue
            ##if isubcase != 17:
                ##continue
            #loadcase_id = subcase.get_parameter('LOAD')[0]
            #load = model.loads[loadcase_id][0]
            #elem = load.eid

            ##msg = '%s%s\n' % (elem.nodes[0], elem.nodes[1])

            #f, m = sum_forces_moments(model, p0, loadcase_id, include_grav=False)
            #eids = None
            #nids = None
            #f2, m2 = sum_forces_moments_elements(
                #model, p0, loadcase_id, eids, nids, include_grav=False)
            #assert allclose(f, f2), 'f=%s f2=%s' % (f, f2)
            #assert allclose(m, m2), 'm=%s m2=%s' % (m, m2)

            #case = op2.spc_forces[isubcase]
            #fm = -case.data[0, :, :].sum(axis=0)
            #assert len(fm) == 6, fm
            #if not allclose(f[0], fm[0]):
                #model.log.error('subcase=%-2i Fx f=%s fexpected=%s' % (
                    #isubcase, f.tolist(), fm.tolist()))
                #fail = True
            #if not allclose(f[1], fm[1]):
                #model.log.error('subcase=%-2i Fy f=%s fexpected=%s' % (
                    #isubcase, f.tolist(), fm.tolist()))
                #fail = True
            #if not allclose(f[2], fm[2]):
                #model.log.error('subcase=%-2i Fz f=%s fexpected=%s' % (
                    #isubcase, f.tolist(), fm.tolist()))
                #fail = True

            #if not allclose(m[0], fm[3]):
                #model.log.error('subcase=%-2i Mx m=%s fexpected=%s' % (
                    #isubcase, m.tolist(), fm.tolist()))
                #fail = True
            #if not allclose(m[1], fm[4]):
                #model.log.error('subcase=%-2i My m=%s fexpected=%s' % (
                    #isubcase, m.tolist(), fm.tolist()))
                #fail = True
            #if not allclose(m[2], fm[5]):
                #model.log.error('subcase=%-2i Mz m=%s fexpected=%s' % (
                    #isubcase, m.tolist(), fm.tolist()))
                #fail = True
        #if fail:
            #raise RuntimeError('incorrect loads')

    def test_ploadx1(self):
        """tests a CTRIAX, PLPLANE, MATHP, and PLOADX1"""
        model = BDF(debug=False)
        sid = 10
        eid1 = 11
        pa = 200.
        ga = 1
        gb = 2
        ploadx1 = model.add_ploadx1(sid, eid1, pa, [ga, gb], pb=None,
                                    theta=0., comment='ploadx1')
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [1., 0., 0.])
        model.add_grid(3, [1., 1., 0.])

        pid = 20
        nids = [1, 2, 3, None, None, None]
        ctriax = model.add_ctriax(eid1, pid, nids, theta_mcid=0., comment='ctriax')

        mid = 21
        plplane = model.add_plplane(pid, mid, cid=0,
                                    stress_strain_output_location='GRID',
                                    comment='plplane')

        #eid2 = 12
        #model.add_ctriax6(eid2, mid, nids, theta=0., comment='ctriax6')

        #E = 30.e7
        #G = None
        #nu = 0.3
        #mat1 = model.add_mat1(mid, E, G, nu, rho=0.1, comment='mat1')
        #mathe = model.add_mathe(mid, model, bulk, rho, texp, mus, alphas,
                                #betas, mooney, sussbat, comment='mathe')
        mathp = model.add_mathp(mid, comment='mathp')
        model.validate()

        ctriax.raw_fields()
        ctriax.write_card(size=8)
        ctriax.write_card(size=16)

        plplane.raw_fields()
        plplane.write_card(size=8)
        plplane.write_card(size=16)

        #mathe.raw_fields()
        #mathe.write_card(size=8)
        #mathe.write_card(size=16)

        mathp.raw_fields()
        mathp.write_card(size=8)
        mathp.write_card(size=16)

        ploadx1.raw_fields()
        ploadx1.write_card(size=8)
        ploadx1.write_card(size=16)
        ploadx1.write_card(size=16, is_double=True)

        model.validate()
        model._verify_bdf(xref=False)
        model.cross_reference()
        model._verify_bdf(xref=True)

        ctriax.write_card(size=8)
        plplane.write_card(size=8)
        #mathe.write_card(size=8)
        mathp.write_card(size=8)
        ploadx1.write_card(size=8)
        model.write_bdf('ploadx1.temp')

        model2 = read_bdf('ploadx1.temp', debug=None)
        model2._verify_bdf()
        os.remove('ploadx1.temp')
        save_load_deck(model2, run_convert=False, run_mass_properties=False)

    def test_loads_combo(self):
        r"""
        tests CONROD, CTRIA3-PSHELL, CQUAD4-PCOMP,
        CTETRA/CPENTA/CPYRAM/CHEXA-PSOLID
        FORCE, FORCE1,
        PLOAD4 (CHEXA)
        PLOAD2 (CTRIA3, CQUAD4)

        ^ y
        |
        4     3 12
        +-----+--+
        |     |     + 13
        |     |     |
        +-----+--+--+---S  -> x
        1     2  9  10  11
        """
        model = BDF(debug=False)
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [1., 0., 0.])
        model.add_grid(3, [1., 1., 0.])
        model.add_grid(4, [0., 1., 0.])

        model.add_grid(5, [0., 0., 1.])
        model.add_grid(6, [1., 0., 1.])
        model.add_grid(7, [1., 1., 1.])
        model.add_grid(8, [0., 1., 1.])

        model.add_grid(9, [5., 0., 0.])
        model.add_grid(10, [6., 0., 0.])
        model.add_grid(12, [2., 1., 0.])
        model.add_grid(13, [2., 0.5, 0.])

        eid = 1
        mid = 1
        A = 2.0
        nids = [1, 2]
        # L = 1; A=2
        # mass=(rho*A + nsm) * L = (0.2*2 + 1.0) * 1 = 1.4
        conrod = model.add_conrod(eid, mid, nids, A, j=0.0, c=0.0, nsm=1.0, comment='')
        model.add_conrod(eid, mid, nids, A, j=0.0, c=0.0, nsm=1.0, comment='')

        eid_tube = 2
        pid_tube = 2
        nids = [3, 12]
        unused_ctube = model.add_ctube(eid_tube, pid_tube, nids, comment='ctube')
        unused_ctube = model.add_ctube(eid_tube, pid_tube, nids, comment='ctube')
        OD1 = 0.1
        unused_ptube = model.add_ptube(pid_tube, mid, OD1, t=None, nsm=0., OD2=None,
                                       comment='ptube')
        model.add_ptube(pid_tube, mid, OD1, t=None, nsm=0., OD2=None,
                        comment='ptube')

        E = 3.0e7
        G = None
        nu = 0.3
        model.add_mat1(mid, E, G, nu, rho=0.2, a=0.0, tref=0.0, ge=0.0,
                       St=0.0, Sc=0.0, Ss=0.0, mcsid=0,
                       comment='')

        eid = 3
        pid = 3
        nids = [1, 2, 3]
        ctria3 = model.add_ctria3(eid, pid, nids, zoffset=0., theta_mcid=0.0, tflag=0,
                                  T1=None, T2=None, T3=None,
                                  comment='')

        # mass = (rho*t + nsm)*A = (0.2*0.5 + 0.3) * 0.5 = 0.4 * 0.5 = 0.2
        model.add_pshell(pid, mid1=None, t=0.5, mid2=mid, twelveIt3=1.0,
                         mid3=None, tst=0.833333,
                         nsm=0.3, z1=None, z2=None,
                         mid4=None, comment='')


        eid = 4
        pid = 4
        nids = [1, 2, 3, 4]
        unused_cquad4 = model.add_cquad4(eid, pid, nids, theta_mcid=0.0, zoffset=0.,
                                         tflag=0, T1=1.0, T2=1.0, T3=1.0, T4=1.0, comment='')
        mids = [mid, mid, mid]
        thicknesses = [0.1, 0.2, 0.3]
        model.add_pcomp(pid, mids, thicknesses, thetas=None, souts=None,
                        nsm=0., sb=0., ft=None,
                        tref=0., ge=0., lam=None,
                        z0=None, comment='pcomp')
        model.add_pcomp(pid, mids, thicknesses, thetas=None, souts=None,
                        nsm=0., sb=0., ft=None,
                        tref=0., ge=0., lam=None,
                        z0=None, comment='pcomp')

        pid = 5
        global_ply_ids = [5, 6, 7]
        mids = [mid, mid, mid]
        thicknesses = [0.1, 0.2, 0.3]
        unused_pcompg = model.add_pcompg(pid, global_ply_ids, mids, thicknesses, thetas=None,
                                         souts=None, nsm=0.0, sb=0.0, ft=None, tref=0.0, ge=0.0,
                                         lam=None, z0=None, comment='pcompg')
        model.add_pcompg(pid, global_ply_ids, mids, thicknesses, thetas=None,
                         souts=None, nsm=0.0, sb=0.0, ft=None, tref=0.0, ge=0.0,
                         lam=None, z0=None, comment='pcompg')


        pid = 40
        eid = 5
        nids = [1, 2, 3, 5]
        unused_ctetra = model.add_ctetra(eid, pid, nids, comment='ctetra')

        eid = 6
        nids = [1, 2, 3, 4, 5]
        unused_cpyram = model.add_cpyram(eid, pid, nids, comment='cpyram')

        eid = 7
        nids = [1, 2, 3, 5, 6, 7]
        unused_cpenta = model.add_cpenta(eid, pid, nids, comment='cpenta')

        eid = 8
        nids = [1, 2, 3, 4, 5, 6, 7, 8]
        chexa = model.add_chexa(eid, pid, nids, comment='chexa')
        # mass = rho*V = 0.2*1

        unused_psolid = model.add_psolid(pid, mid, cordm=0, integ=None, stress=None,
                                         isop=None, fctn='SMECH', comment='psolid')


        conid = 42
        nodes = [
            2, 2, 2, 2, 2, 2,
            9, 9, 9, 9, 9, 9,
        ]
        components = [
            '1', '2', '3', '4', '5', '6',
            '1', '2', '3', '4', '5', '6',
        ]
        coefficients = [
            1., 1., 1., 1., 1., 1.,
            1., 1., 1., 1., 1., 1.,
        ]
        unused_mpc = model.add_mpc(conid, nodes, components, coefficients, comment='mpc')

        conid = 43
        nodes = [10, 11]
        components = ['1', '0']
        coefficients = [1., 1.]
        unused_mpc = model.add_mpc(conid, nodes, components, coefficients)
        model.add_spoint(11, comment='spoint')
        conid = 44
        sets = [42, 43]
        unused_mpcadd = model.add_mpcadd(conid, sets, comment='mpcadd')
        #model.add_spoint([11, 'THRU', 42], comment='spoint3')
        str(model.spoints)

        sid = 12
        xyz = [2., 3., 4.]
        node = 7
        mag = 1.0
        unused_force = model.add_force(sid, node, mag, xyz, cid=0, comment='force')
        unused_moment = model.add_moment(sid, node, mag, xyz, comment='moment')

        node = 6
        mag = 1.0
        g1 = 2
        g2 = 3
        unused_force1 = model.add_force1(sid, node, mag, g1, g2, comment='force1')
        unused_moment1 = model.add_moment1(sid, node, mag, g1, g2, comment='moment1')

        g1 = 1
        g2 = 3
        g3 = 2
        g4 = 4
        unused_force2 = model.add_force2(sid, node, mag, g1, g2, g3, g4, comment='force2')
        unused_moment2 = model.add_moment2(sid, node, mag, g1, g2, g3, g4, comment='moment2')
        #g2, g3 = g3, g2
        #force2 = model.add_force2(sid, node, mag, g1, g2, g3, g4, comment='force2')

        load_id = 120
        scale = 1.
        scale_factors = [1.0, 2.0]
        load_ids = [12, 13] # force, pload4
        unused_load = model.add_load(load_id, scale, scale_factors, load_ids, comment='load')

        # loads
        sid = 13
        pressure = 2.0
        model.add_pload(sid, pressure, [1, 2, 3], comment='pload')
        model.add_pload2(sid, pressure, [3, 4], comment='pload2')  # ctria3, cquad4

        eids = [8] # hexa
        g1 = 6
        unused_g34 = 8
        pressures = [1., 1., 1., 1.]
        unused_pload4 = model.add_pload4(sid, eids, pressures, g1=1, g34=8,
                                         cid=None, nvector=None, surf_or_line='SURF',
                                         line_load_dir='NORM', comment='pload4')
        #print(model.loads)
        #print(model.load_combinations)

        #-----------------------------------------------------------------------
        # constraints
        conid = 42
        nodes = [1, 2]
        components = ['123', '123']
        enforced = [0., 0.]
        unused_spc = model.add_spc(conid, nodes, components, enforced, comment='spc')
        conid = 43
        nodes = [1, 2]
        components2 = '123456'
        unused_spc1 = model.add_spc1(conid, components2, nodes, comment='spc1')
        conid = 44
        sets = [42, 43]
        unused_spcadd = model.add_spcadd(conid, sets, comment='spcadd')

        #-----------------------------------------------------------------------
        sid = 12
        model.add_eigrl(sid, v1=None, v2=None, nd=None, msglvl=0,
                        maxset=None, shfscl=None, norm=None,
                        options=None, values=None, comment='eigrl')

        sid = 13
        model.add_eigr(sid, method='LAN', f1=None, f2=None, ne=None, nd=20,
                       norm='MASS', G=None, C=None,
                       comment='')
        #-----------------------------------------------------------------------
        model.validate()
        model._verify_bdf(xref=False)
        model.write_bdf('loads.temp')
        model.cross_reference()
        assert allclose(conrod.Mass(), 1.4)
        assert allclose(ctria3.Mass(), 0.2)
        assert allclose(chexa.Mass(), 0.2)


        model.write_bdf('loads.temp')
        model._verify_bdf(xref=True)
        model.write_bdf('loads.temp')
        model.uncross_reference()
        model.write_bdf('loads.temp')
        model.cross_reference()
        model.uncross_reference()
        model.safe_cross_reference()
        model.write_bdf('loads.temp', size=8, is_double=False)
        model.write_bdf('loads.temp', size=16, is_double=False)
        model.write_bdf('loads.temp', size=16, is_double=True)

        model2 = read_bdf('loads.temp', debug=None)
        os.remove('loads.temp')
        eids = list(model.elements.keys())
        nids = list(model.nodes.keys())
        p0 = [0., 0., 0.]
        loadcase_id = 120
        forces1, moments1 = sum_forces_moments_elements(model2, p0, loadcase_id, eids, nids,
                                                        include_grav=False, xyz_cid0=None)
        forces2, moments2 = sum_forces_moments(model2, p0, loadcase_id, include_grav=False,
                                               xyz_cid0=None)
        assert allclose(forces1, forces2), 'forces1=%s forces2=%s' % (forces1, forces2)
        assert allclose(moments1, moments2), 'moments1=%s moments2=%s' % (moments1, moments2)

        model2.get_area_breakdown()
        model2.get_volume_breakdown()
        model2.get_mass_breakdown()

        eids = list(model.elements.keys())
        p0 = [0., 0., 0.]
        forces1, moments1 = sum_forces_moments(model, p0, loadcase_id,
                                               include_grav=True, xyz_cid0=None)
        #forces2, moments2 = sum_forces_moments_elements(model, p0, loadcase_id, eids, nids,
                                                        #include_grav=True, xyz_cid0=None)
        #assert np.array_equal(forces1, forces2)
        #assert np.array_equal(moments1, moments2)


        model2.write_skin_solid_faces('skin.bdf', write_solids=False,
                                      write_shells=True)
        os.remove('skin.bdf')
        save_load_deck(model2)

    def test_load(self):
        """makes sure LOAD cards don't get sorted"""
        model = BDF(debug=False, log=log)

        load = model.add_load(sid=13, scale=1., scale_factors=[0.5, 0.1], load_ids=[11, 10])
        msg8 = load.write_card(size=8, is_double=False)
        load_expected = 'LOAD          13      1.      .5      11      .1      10'
        assert msg8.rstrip() == load_expected, '%r' % msg8

        load2_expected = 'LOAD          14      1.      .5      11      .1      10      .4      11'
        load2 = model.add_load(sid=14, scale=1.,
                               scale_factors=[0.5, 0.1, 0.4], load_ids=[11, 10, 11])
        msg8 = load2.write_card(size=8, is_double=False)
        assert msg8.rstrip() == load2_expected, '%r' % msg8
        model.validate()

    def test_load_sort(self):
        """makes sure LOAD cards don't get sorted"""
        model = BDF(debug=False, log=log)
        unused_load2 = model.add_load(sid=14, scale=1.,
                                      scale_factors=[0.5, 0.1, 0.4], load_ids=[11, 10])
        with self.assertRaises(IndexError):
            model.validate()

    def test_sload(self):
        """tests SLOAD"""
        log = get_logger(level='warning')
        model = BDF(log=log)
        model.add_spoint([11, 12])
        sid = 14
        nids = 11 # SPOINT
        mags = 20.
        sload = model.add_sload(sid, nids, mags, comment='an sload')
        sload.raw_fields()
        sload.repr_fields()

        sid = 14
        nids = [11, 12] # SPOINT, GRID
        mags = [20., 30.]
        unused_sload = model.add_sload(sid, nids, mags, comment='an sload')
        model.validate()

        model.cross_reference()
        model.uncross_reference()
        model.safe_cross_reference()

        p0 = [0., 0., 0.]
        eids = list(model.elements.keys())
        nids = list(model.nodes.keys())
        #print('nids =', nids) # empty
        loadcase_id = sid
        forces1, moments1 = sum_forces_moments(model, p0, loadcase_id,
                                               include_grav=True, xyz_cid0=None)
        forces2, moments2 = sum_forces_moments_elements(model, p0, loadcase_id, eids, nids,
                                                        include_grav=True, xyz_cid0=None)
        assert np.array_equal(forces1, forces2)
        assert np.array_equal(moments1, moments2)
        save_load_deck(model, run_convert=False)

    def test_loads_nonlinear(self):
        """tests a nonlinear variation on a FORCE and PLOAD4"""
        model = BDF()

        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [1., 0., 0.])
        model.add_grid(3, [1., 1., 0.])
        model.add_grid(4, [0., 1., 0.])

        eid = 10
        pid = 100
        mid = 1000
        nids = [1, 2, 3, 4]
        model.add_cquad4(eid, pid, nids)
        model.add_pshell(pid, mid1=mid, t=0.1, mid2=mid)
        E = 3.0e7
        G = None
        nu = 0.3
        model.add_mat1(mid, E, G, nu)
        spc_id = 10000
        model.add_spc1(spc_id, 123456, [1, 2])
        pload4_sid = 20000
        eids = 10
        pressures = 100.
        model.add_pload4(pload4_sid, eids, pressures)

        table_id = 30000
        theta = np.linspace(0., 1., num=50, endpoint=True)
        y = np.sin(theta)
        x = theta
        model.add_tabled1(table_id, x, y, xaxis='LINEAR', yaxis='LINEAR', extrap=0, comment='')


        dload_id = 40000
        scale = 1.
        scale_factors = 1.
        load_ids = [40001]
        model.add_dload(dload_id, scale, scale_factors, load_ids, comment='dload')
        #model.add_darea(sid, p, c, scale)

        sid = 40001
        excite_id = pload4_sid
        model.add_rload1(sid, excite_id, delay=0, dphase=0, tc=table_id, td=0, Type='LOAD', comment='rload')

        lines = [
            'SUBCASE 1',
            '  DLOAD = %s' % dload_id,
        ]
        cc = CaseControlDeck(lines, log=model.log)
        model.case_control_deck = cc
        model.validate()
        model.cross_reference()

    def test_loads_nonlinear2(self):
        """tests a nonlinear variation on a FORCE and PLOAD4"""
        model = BDF()

        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [1., 0., 0.])
        model.add_grid(3, [1., 1., 0.])
        model.add_grid(4, [0., 1., 0.])

        eid = 10
        pid = 100
        mid = 1000
        nids = [1, 2, 3, 4]
        model.add_cquad4(eid, pid, nids)
        model.add_pshell(pid, mid1=mid, t=0.1, mid2=mid)
        E = 3.0e7
        G = None
        nu = 0.3
        model.add_mat1(mid, E, G, nu)
        spc_id = 10000
        model.add_spc1(spc_id, 123456, [1, 2])
        #pload4_sid = 20000
        #eids = 10
        #pressures = 100.
        #model.add_pload4(pload4_sid, eids, pressures)

        table_id = 30000
        theta = np.linspace(0., 1., num=50, endpoint=True)
        y = np.sin(theta)
        x = theta
        model.add_tabled1(table_id, x, y, xaxis='LINEAR', yaxis='LINEAR', extrap=0, comment='')


        dload_id = 40000
        scale = 1.
        scale_factors = 1.
        load_ids = [40001]
        model.add_dload(dload_id, scale, scale_factors, load_ids, comment='dload')
        #model.add_darea(sid, p, c, scale)

        tload_id = 40001
        darea_id = 40002
        excite_id = darea_id
        model.add_tload1(tload_id, excite_id, table_id, delay=0, Type='LOAD', us0=0.0, vs0=0.0, comment='tload1')

        nid = 4
        component = 3
        scale = 2.0
        model.add_darea(darea_id, nid, component, scale, comment='darea')
        lines = [
            'SUBCASE 1',
            '  DLOAD = %s' % dload_id,
        ]
        cc = CaseControlDeck(lines, log=model.log)
        model.case_control_deck = cc
        model.validate()
        model.cross_reference()

    def test_loads_nonlinear3(self):
        """tests a nonlinear variation on a FORCE and PLOAD4"""
        model = BDF()

        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [1., 0., 0.])
        model.add_grid(3, [1., 1., 0.])
        model.add_grid(4, [0., 1., 0.])

        eid = 10
        pid = 100
        mid = 1000
        nids = [1, 2, 3, 4]
        model.add_cquad4(eid, pid, nids)
        model.add_pshell(pid, mid1=mid, t=0.1, mid2=mid)
        E = 3.0e7
        G = None
        nu = 0.3
        model.add_mat1(mid, E, G, nu)
        spc_id = 10000
        model.add_spc1(spc_id, 123456, [1, 2])
        #pload4_sid = 20000
        #eids = 10
        #pressures = 100.
        #model.add_pload4(pload4_sid, eids, pressures)

        table_id = 30000
        theta = np.linspace(0., 1., num=50, endpoint=True)
        y = np.sin(theta)
        x = theta
        model.add_tabled1(table_id, x, y, xaxis='LINEAR', yaxis='LINEAR', extrap=0, comment='')


        #dload_id = 40000
        #scale = 1.
        #scale_factors = 1.
        #load_ids = [40001]
        #model.add_dload(dload_id, scale, scale_factors, load_ids, comment='dload')
        #model.add_darea(sid, p, c, scale)

        tload_id = 40001
        darea_id = 40002
        excite_id = darea_id
        model.add_tload1(tload_id, excite_id, table_id, delay=0, Type='LOAD', us0=0.0, vs0=0.0, comment='tload1')

        nid = 4
        component = 3
        scale = 2.0
        model.add_darea(darea_id, nid, component, scale, comment='darea')
        lines = [
            'SUBCASE 1',
            '  DLOAD = %s' % tload_id,
        ]
        cc = CaseControlDeck(lines, log=model.log)
        model.case_control_deck = cc
        model.validate()
        model.cross_reference()

    def test_loads_nonlinear_thermal1(self):
        """tests a nonlinear variation on a FORCE and PLOAD4"""
        model = BDF()

        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [1., 0., 0.])
        model.add_grid(3, [1., 1., 0.])
        model.add_grid(4, [0., 1., 0.])

        eid = 10
        pid = 100
        mid = 1000
        nids = [1, 2, 3, 4]
        model.add_cquad4(eid, pid, nids)
        model.add_pshell(pid, mid1=mid, t=0.1, mid2=mid)
        E = 3.0e7
        G = None
        nu = 0.3
        model.add_mat1(mid, E, G, nu)
        spc_id = 10000
        model.add_spc1(spc_id, 123456, [1, 2])
        #pload4_sid = 20000
        #eids = 10
        #pressures = 100.
        #model.add_pload4(pload4_sid, eids, pressures)

        table_id = 30000
        theta = np.linspace(0., 1., num=50, endpoint=True)
        y = np.sin(theta)
        x = theta
        model.add_tabled1(table_id, x, y, xaxis='LINEAR', yaxis='LINEAR', extrap=0, comment='')


        #dload_id = 40000
        #scale = 1.
        #scale_factors = 1.
        #load_ids = [40001]
        #model.add_dload(dload_id, scale, scale_factors, load_ids, comment='dload')
        #model.add_darea(sid, p, c, scale)

        tload_id = 40001
        qvect_id = 40002
        excite_id = qvect_id
        model.add_tload1(tload_id, excite_id, table_id, delay=0, Type='LOAD', us0=0.0, vs0=0.0, comment='tload1')
        #nid = 4
        eids = [eid]
        q0 = 2.0
        model.add_qvect(qvect_id, q0, eids, t_source=None, ce=0, vector_tableds=None, control_id=0, comment='qvect')
        lines = [
            'SUBCASE 1',
            '  DLOAD = %s' % tload_id,
        ]
        cc = CaseControlDeck(lines, log=model.log)
        model.case_control_deck = cc
        model.validate()
        model.cross_reference()

if __name__ == '__main__':  # pragma: no cover
    unittest.main()
