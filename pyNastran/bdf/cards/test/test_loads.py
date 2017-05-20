from __future__ import print_function
import os
import unittest
from six import iteritems
import numpy as np
from numpy import array, allclose, array_equal, set_printoptions
set_printoptions(suppress=True, precision=3)

import pyNastran
from pyNastran.bdf.bdf import BDF, BDFCard, DAREA, PLOAD4, read_bdf, RROD
from pyNastran.bdf.errors import DuplicateIDsError
from pyNastran.op2.op2 import OP2

bdf = BDF(debug=False)
test_path = pyNastran.__path__[0]


class TestLoads(unittest.TestCase):
    def test_force(self):
        """CONROD, FORCE"""
        model = BDF(debug=False)
        eid = 1
        mid = 100
        nids = [10, 11]
        A = 3.14
        model.add_conrod(eid, mid, nids, A, j=0.0, c=0.0, nsm=0.0,
                         comment='')
        model.add_grid(10, xyz=[10., 0., 0.])
        model.add_grid(11, xyz=[11., 0., 0.])
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
        assert np.array_equal(force.F()[11], np.array([42., 42., 84.])), force.F()
        model.cross_reference()
        force.raw_fields()

    def test_moment(self):
        """CONROD, MOMENT"""
        model = BDF(debug=False)
        eid = 1
        mid = 100
        nids = [10, 11]
        A = 3.14
        model.add_conrod(eid, mid, nids, A, j=0.0, c=0.0, nsm=0.0,
                         comment='')
        model.add_grid(10, xyz=[10., 0., 0.])
        model.add_grid(11, xyz=[11., 0., 0.])
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
        assert np.array_equal(moment.M()[11], np.array([42., 42., 84.])), force.M()
        model.cross_reference()
        moment.raw_fields()

    def test_accel1(self):
        """tests ACCEL1"""
        model = BDF(debug=False)
        sid = 42
        N = [0., 0., 1.]
        nodes = [10, 11]
        scale = 3.14
        accel1 = model.add_accel1(sid, scale, N, nodes, cid=0, comment='accel1')
        accel1.raw_fields()
        accel1.write_card(size=8)
        accel1.write_card(size=16)
        accel1.write_card(size=16, is_double=True)

        model.add_grid(10, xyz=[10., 0., 0.])
        model.add_grid(11, xyz=[11., 0., 0.])
        model.validate()
        model.pop_parse_errors()
        model.cross_reference()
        model.pop_xref_errors()

        accel1.raw_fields()
        accel1.write_card(size=8)
        accel1.write_card(size=16)
        accel1.write_card(size=16, is_double=True)

    def test_accel(self):
        """tests ACCEL"""
        model = BDF(debug=False)
        sid = 42
        N = [0., 0., 1.]
        nodes = [10, 11]
        scale = 3.14
        direction = 'Z'
        locs = [11., 22., 33.]
        vals = [1., 2., 3.]
        accel = model.add_accel(sid, N, direction, locs, vals, cid=0,
                                comment='accel')
        accel.raw_fields()
        accel.write_card(size=8)
        accel.write_card(size=16)
        accel.write_card(size=16, is_double=True)

        model.add_grid(10, xyz=[10., 0., 0.])
        model.add_grid(11, xyz=[11., 0., 0.])
        model.validate()
        model.pop_parse_errors()
        model.cross_reference()
        model.pop_xref_errors()

        accel.raw_fields()
        accel.write_card(size=8)
        accel.write_card(size=16)
        accel.write_card(size=16, is_double=True)

    def test_darea_01(self):
        """tests a DAREA"""
        #DAREA SID P1 C1 A1  P2 C2 A2
        #DAREA 3   6   2 8.2 15 1  10.1
        lines = ['DAREA,3,6,2,8.2,15,1,10.1']
        card = bdf.process_card(lines)
        cardi = BDFCard(card)

        size = 8
        card = DAREA.add_card(cardi)
        card.write_card(size, 'dummy')
        card.raw_fields()

    def test_pload4_01(self):
        """tests a PLOAD4"""
        lines = ['PLOAD4  1000    1       -60.    -60.    60.             1']
        card = bdf.process_card(lines)
        cardi = BDFCard(card)

        size = 8
        card = PLOAD4.add_card(cardi)
        card.write_card(size, 'dummy')
        card.raw_fields()

    def test_pload4_02(self):
        """tests a PLOAD4"""
        lines = ['PLOAD4  1       101     1.                              10000   10011']
        card = bdf.process_card(lines)
        cardi = BDFCard(card)

        size = 8
        card = PLOAD4.add_card(cardi)
        card.write_card(size, 'dummy')
        card.raw_fields()

    def test_pload4_cpenta(self):
        """tests a PLOAD4 with a CPENTA"""
        bdf_filename = os.path.join(test_path, '..', 'models', 'pload4', 'cpenta.bdf')
        op2_filename = os.path.join(test_path, '..', 'models', 'pload4', 'cpenta.op2')
        op2 = OP2(debug=False)
        op2.read_op2(op2_filename)

        model_b = BDF(debug=False)
        model_b.read_bdf(bdf_filename)
        # p0 = (model_b.nodes[21].xyz + model_b.nodes[22].xyz + model_b.nodes[23].xyz) / 3.
        p0 = model_b.nodes[21].xyz
        angles = [
            (23, 24), (24, 23),
            (21, 26), (26, 21),
        ]
        nx = [
            (23, 25), (25, 23),
            (22, 26), (26, 22),
        ]

        msg = ''
        for isubcase, subcase in sorted(iteritems(model_b.subcases)):
            if isubcase == 0:
                continue
            #if isubcase != 17:
                #continue
            loadcase_id = subcase.get_parameter('LOAD')[0]
            load = model_b.loads[loadcase_id][0]
            elem = load.eids[0]
            g1 = load.g1.nid
            if load.g34 is None:
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
                g34 = load.g34.nid
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

            f, m = model_b.sum_forces_moments(p0, loadcase_id, include_grav=False)
            eids = None
            nids = None
            f2, m2 = model_b.sum_forces_moments_elements(p0, loadcase_id, eids, nids, include_grav=False)
            assert allclose(f, f2), 'f=%s f2=%s' % (f, f2)
            assert allclose(m, m2), 'm=%s m2=%s' % (m, m2)

            case = op2.spc_forces[isubcase]
            fm = -case.data[0, :3, :].sum(axis=0)
            assert len(fm) == 6, fm
            if not allclose(f[0], fm[0]):
                print('%-2i Fx f=%s fexpected=%s face=%s' % (isubcase, f, fm, face))
            if not allclose(f[1], fm[1]):
                print('%-2i Fy f=%s fexpected=%s face=%s' % (isubcase, f, fm, face))
            if not allclose(f[2], fm[2]):
                print('%-2i Fz f=%s fexpected=%s face=%s' % (isubcase, f, fm, face))
            # if not allclose(m[0], fm[3]):
                # print('%i Mx m=%s fexpected=%s' % (isubcase, m, fm))
            # if not allclose(m[1], fm[4]):
                # print('%i My m=%s fexpected=%s' % (isubcase, m, fm))
            # if not allclose(m[2], fm[5]):
                # print('%i Mz m=%s fexpected=%s' % (isubcase, m, fm))

            #self.assertEqual(f[0], fm[0], 'f=%s fexpected=%s' % (f, fm[:3]))
            #self.assertEqual(f[1], fm[1], 'f=%s fexpected=%s' % (f, fm[:3]))
            #self.assertEqual(f[2], fm[2], 'f=%s fexpected=%s' % (f, fm[:3]))
            #self.assertEqual(m[0], fm[3], 'm=%s mexpected=%s' % (m, fm[3:]))
            #self.assertEqual(m[1], fm[4], 'm=%s mexpected=%s' % (m, fm[3:]))
            #self.assertEqual(m[2], fm[5], 'm=%s mexpected=%s' % (m, fm[3:]))

    def test_pload4_ctria3(self):
        """tests a PLOAD4 with a CTRIA3"""
        bdf_filename = os.path.join(test_path, '..', 'models', 'pload4', 'ctria3.bdf')
        op2_filename = os.path.join(test_path, '..', 'models', 'pload4', 'ctria3.op2')
        op2 = OP2(debug=False)
        op2.read_op2(op2_filename)

        model_b = BDF(debug=False)
        model_b.read_bdf(bdf_filename)
        # p0 = (model_b.nodes[21].xyz + model_b.nodes[22].xyz + model_b.nodes[23].xyz) / 3.
        p0 = model_b.nodes[21].xyz

        for isubcase, subcase in sorted(iteritems(model_b.subcases)):
            if isubcase == 0:
                continue
            #if isubcase != 17:
                #continue
            loadcase_id = subcase.get_parameter('LOAD')[0]
            load = model_b.loads[loadcase_id][0]
            elem = load.eids[0]
            area = 0.5
            centroid = elem.Centroid()
            normal = elem.Normal()

            msg = '%s%s%s\n' % (
                elem.nodes[0], elem.nodes[1], elem.nodes[2])

            assert array_equal(centroid, array([2/3., 1/3., 0.])), 'centroid=%s\n%s' % (centroid, msg)
            assert array_equal(normal, array([0., 0., 1.])), 'normal=%s\n%s' % (normal, msg)

            f, m = model_b.sum_forces_moments(p0, loadcase_id, include_grav=False)
            eids = None
            nids = None
            f2, m2 = model_b.sum_forces_moments_elements(
                p0, loadcase_id, eids, nids, include_grav=False)
            assert allclose(f, f2), 'f=%s f2=%s' % (f, f2)
            assert allclose(m, m2), 'm=%s m2=%s' % (m, m2)

            case = op2.spc_forces[isubcase]
            fm = -case.data[0, :, :].sum(axis=0)
            assert len(fm) == 6, fm
            if not allclose(f[0], fm[0]):
                print('%-2i Fx f=%s fexpected=%s' % (isubcase, f, fm))
            if not allclose(f[1], fm[1]):
                print('%-2i Fy f=%s fexpected=%s' % (isubcase, f, fm))
            if not allclose(f[2], fm[2]):
                print('%-2i Fz f=%s fexpected=%s' % (isubcase, f, fm))

    def test_pload4_cquad4(self):
        """tests a PLOAD4 with a CQUAD4"""
        bdf_filename = os.path.join(test_path, '..', 'models', 'pload4', 'cquad4.bdf')
        op2_filename = os.path.join(test_path, '..', 'models', 'pload4', 'cquad4.op2')
        op2 = OP2(debug=False)
        op2.read_op2(op2_filename)

        model_b = BDF(debug=False)
        model_b.read_bdf(bdf_filename)
        # p0 = (model_b.nodes[21].xyz + model_b.nodes[22].xyz + model_b.nodes[23].xyz) / 3.
        p0 = model_b.nodes[21].xyz

        eids = None
        nids = None
        for isubcase, subcase in sorted(iteritems(model_b.subcases)):
            if isubcase == 0:
                continue
            #if isubcase != 17:
                #continue
            loadcase_id = subcase.get_parameter('LOAD')[0]
            load = model_b.loads[loadcase_id]
            loadi = load[0]
            if loadi.type == 'PLOAD4':
                elem = loadi.eids[0]
                area = 1.0
                centroid = elem.Centroid()
                normal = elem.Normal()
                msg = '%s%s%s\n' % (elem.nodes[0], elem.nodes[1], elem.nodes[2])

                assert array_equal(centroid, array([0.5, 0.5, 0.])), 'centroid=%s\n%s' % (centroid, msg)
                assert array_equal(normal, array([0., 0., 1.])), 'normal=%s\n%s' % (normal, msg)

            f, m = model_b.sum_forces_moments(p0, loadcase_id, include_grav=False)
            f2, m2 = model_b.sum_forces_moments_elements(p0, loadcase_id, eids, nids, include_grav=False)
            assert allclose(f, f2), 'f=%s f2=%s' % (f, f2)
            assert allclose(m, m2), 'm=%s m2=%s' % (m, m2)

            case = op2.spc_forces[isubcase]
            fm = -case.data[0, :, :].sum(axis=0)
            assert len(fm) == 6, fm
            if not allclose(f[0], fm[0]):
                print('%-2i Fx f=%s fexpected=%s' % (isubcase, f, fm))
            if not allclose(f[1], fm[1]):
                print('%-2i Fy f=%s fexpected=%s' % (isubcase, f, fm))
            if not allclose(f[2], fm[2]):
                print('%-2i Fz f=%s fexpected=%s' % (isubcase, f, fm))

    def test_pload4_ctetra(self):
        """tests a PLOAD4 with a CTETRA"""
        bdf_filename = os.path.join(test_path, '..', 'models', 'pload4', 'ctetra.bdf')
        op2_filename = os.path.join(test_path, '..', 'models', 'pload4', 'ctetra.op2')
        op2 = OP2(debug=False)
        op2.read_op2(op2_filename)

        model_b = BDF(debug=False)
        model_b.read_bdf(bdf_filename)
        p0 = model_b.nodes[21].xyz

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

        for isubcase, subcase in sorted(iteritems(model_b.subcases)):
            if isubcase == 0:
                continue
            loadcase_id = subcase.get_parameter('LOAD')[0]
            load = model_b.loads[loadcase_id][0]
            elem = load.eids[0]
            g1 = load.g1.nid

            # f, m = model_b.sum_forces_moments(p0, loadcase_id, include_grav=False)
            # case = op2.spc_forces[isubcase]
            # fm = case.data[0, 0, :]#.ravel()
            # if f[0] != fm[0]:
                # print('%i f=%s fexpected=%s' % (isubcase, f, fm))

            g34 = load.g34.nid
            face, area, centroid, normal = elem.get_face_area_centroid_normal(g1, g34)
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
            f, m = model_b.sum_forces_moments(p0, loadcase_id, include_grav=False)
            eids = None
            nids = None
            f2, m2 = model_b.sum_forces_moments_elements(p0, loadcase_id, eids, nids, include_grav=False)
            assert allclose(f, f2), 'f=%s f2=%s' % (f, f2)
            assert allclose(m, m2), 'm=%s m2=%s' % (m, m2)

            case = op2.spc_forces[isubcase]
            fm = -case.data[0, :, :].sum(axis=0)
            assert len(fm) == 6, fm
            if not allclose(f[0], fm[0]):
                print('%-2i Fx g=(%s,%s) f=%s fexpected=%s face=%s normal=%s' % (
                    isubcase, g1, g34, f, fm, face, normal))
            if not allclose(f[1], fm[1]):
                print('%-2i Fy g=(%s,%s) f=%s fexpected=%s face=%s normal=%s' % (
                    isubcase, g1, g34, f, fm, face, normal))
            if not allclose(f[2], fm[2]):
                print('%-2i Fz g=(%s,%s) f=%s fexpected=%s face=%s normal=%s' % (
                    isubcase, g1, g34, f, fm, face, normal))

    def test_pload4_chexa(self):
        """tests a PLOAD4 with a CHEXA"""
        bdf_filename = os.path.join(test_path, '..', 'models', 'pload4', 'chexa.bdf')
        op2_filename = os.path.join(test_path, '..', 'models', 'pload4', 'chexa.op2')
        op2 = OP2(debug=False)
        op2.read_op2(op2_filename)

        model_b = BDF(debug=False)
        model_b.read_bdf(bdf_filename)
        p0 = model_b.nodes[21].xyz
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

        for isubcase, subcase in sorted(iteritems(model_b.subcases)):
            if isubcase == 0:
                continue
            loadcase_id = subcase.get_parameter('LOAD')[0]
            load = model_b.loads[loadcase_id][0]
            elem = load.eids[0]
            g1 = load.g1.nid

            # f, m = model_b.sum_forces_moments(p0, loadcase_id, include_grav=False)
            # case = op2.spc_forces[isubcase]
            # fm = case.data[0, 0, :]#.ravel()
            # if f[0] != fm[0]:
                # print('%i f=%s fexpected=%s' % (isubcase, f, fm))

            g34 = load.g34.nid
            face, area, centroid, normal = elem.get_face_area_centroid_normal(g1, g34)
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
            #self.assertEqual(f[0], fm[0], 'f=%s fexpected=%s' % (f, fm[:3]))
            #self.assertEqual(f[1], fm[1], 'f=%s fexpected=%s' % (f, fm[:3]))
            #self.assertEqual(f[2], fm[2], 'f=%s fexpected=%s' % (f, fm[:3]))
            #self.assertEqual(m[0], fm[3], 'm=%s mexpected=%s' % (m, fm[3:]))
            #self.assertEqual(m[1], fm[4], 'm=%s mexpected=%s' % (m, fm[3:]))
            #self.assertEqual(m[2], fm[5], 'm=%s mexpected=%s' % (m, fm[3:]))
            f, m = model_b.sum_forces_moments(p0, loadcase_id, include_grav=False)
            eids = None
            nids = None
            f2, m2 = model_b.sum_forces_moments_elements(p0, loadcase_id, eids, nids, include_grav=False)
            assert allclose(f, f2), 'f=%s f2=%s' % (f, f2)
            assert allclose(m, m2), 'm=%s m2=%s' % (m, m2)

            case = op2.spc_forces[isubcase]
            fm = -case.data[0, :4, :].sum(axis=0)
            assert len(fm) == 6, fm
            if not allclose(f[0], fm[0]):
                print('%-2i Fx f=%s fexpected=%s face=%s' % (isubcase, f, fm, face))
            if not allclose(f[1], fm[1]):
                print('%-2i Fy f=%s fexpected=%s face=%s' % (isubcase, f, fm, face))
            if not allclose(f[2], fm[2]):
                print('%-2i Fz f=%s fexpected=%s face=%s' % (isubcase, f, fm, face))

    @unittest.expectedFailure
    def test_pload1_cbar(self):
        bdf_filename = os.path.join(test_path, '..', 'models', 'pload4', 'pload1.bdf')
        op2_filename = os.path.join(test_path, '..', 'models', 'pload4', 'pload1.op2')
        op2 = OP2(debug=False)
        op2.read_op2(op2_filename)

        model_b = BDF(debug=False)
        model_b.read_bdf(bdf_filename)
        # p0 = (model_b.nodes[21].xyz + model_b.nodes[22].xyz + model_b.nodes[23].xyz) / 3.
        p0 = model_b.nodes[1].xyz

        fail = False
        dashes = '-' * 80
        for isubcase, subcase in sorted(iteritems(model_b.subcases)):
            if isubcase == 0:
                continue
            #if isubcase != 17:
                #continue
            loadcase_id = subcase.get_parameter('LOAD')[0]
            load = model_b.loads[loadcase_id][0]
            elem = load.eid

            #msg = '%s%s\n' % (elem.nodes[0], elem.nodes[1])

            f, m = model_b.sum_forces_moments(p0, loadcase_id, include_grav=False)
            eids = None
            nids = None
            f2, m2 = model_b.sum_forces_moments_elements(p0, loadcase_id, eids, nids, include_grav=False)
            assert allclose(f, f2), 'f=%s f2=%s' % (f, f2)
            assert allclose(m, m2), 'm=%s m2=%s' % (m, m2)

            case = op2.spc_forces[isubcase]
            fm = -case.data[0, :, :].sum(axis=0)
            assert len(fm) == 6, fm
            if not allclose(f[0], fm[0]):
                print('%-2i Fx f=%s fexpected=%s' % (isubcase, f, fm))
                print(dashes)
                fail = True
            if not allclose(f[1], fm[1]):
                print('%-2i Fy f=%s fexpected=%s' % (isubcase, f, fm))
                print(dashes)
                fail = True
            if not allclose(f[2], fm[2]):
                print('%-2i Fz f=%s fexpected=%s' % (isubcase, f, fm))
                print(dashes)
                fail = True

            if not allclose(m[0], fm[3]):
                print('%-2i Mx m=%s fexpected=%s' % (isubcase, m, fm))
                print(dashes)
                fail = True
            if not allclose(m[1], fm[4]):
                print('%-2i My m=%s fexpected=%s' % (isubcase, m, fm))
                print(dashes)
                fail = True
            if not allclose(m[2], fm[5]):
                print('%-2i Mz m=%s fexpected=%s' % (isubcase, m, fm))
                print(dashes)
                fail = True
        if fail:
            raise RuntimeError('incorrect loads')

    def test_ploadx1(self):
        """tests a PLOADX1"""
        model = BDF(debug=False)
        sid = 10
        eid1 = 11
        pa = 200.
        ga = 1
        gb = 2
        ploadx1 = model.add_ploadx1(sid, eid1, pa, [ga, gb], pb=None,
                                    theta=0., comment='ploadx1')
        model.add_grid(1, xyz=[0., 0., 0.])
        model.add_grid(2, xyz=[1., 0., 0.])
        model.add_grid(3, xyz=[1., 1., 0.])

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

    def test_loads_combo(self):
        r"""
        tests CONROD, CTRIA3-PSHELL, CQUAD4-PCOMP,
        CTETRA/CPENTA/CPYRAM/CHEXA-PSOLID
        FORCE, FORCE1, PLOAD4-CHEXA

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
        model.add_grid(1, xyz=[0., 0., 0.])
        model.add_grid(2, xyz=[1., 0., 0.])
        model.add_grid(3, xyz=[1., 1., 0.])
        model.add_grid(4, xyz=[0., 1., 0.])

        model.add_grid(5, xyz=[0., 0., 1.])
        model.add_grid(6, xyz=[1., 0., 1.])
        model.add_grid(7, xyz=[1., 1., 1.])
        model.add_grid(8, xyz=[0., 1., 1.])

        model.add_grid(9, xyz=[5., 0., 0.])
        model.add_grid(10, xyz=[6., 0., 0.])
        model.add_grid(12, xyz=[2., 1., 0.])
        model.add_grid(13, xyz=[2., 0.5, 0.])

        eid = 1
        mid = 1
        A = 2.0
        nids = [1, 2]
        # L = 1; A=2
        # mass=(rho*A + nsm) * L = (0.2*2 + 1.0) * 1 = 1.4
        conrod = model.add_conrod(eid, mid, nids, A, j=0.0, c=0.0, nsm=1.0, comment='')
        model.add_conrod(eid, mid, nids, A, j=0.0, c=0.0, nsm=1.0, comment='')

        eid = 2
        pid = 2
        nids = [3, 12]
        ctube = model.add_ctube(eid, pid, nids, comment='ctube')
        ctube = model.add_ctube(eid, pid, nids, comment='ctube')
        OD1 = 0.1
        ptube = model.add_ptube(pid, mid, OD1, t=None, nsm=0., OD2=None,
                                comment='ptube')
        model.add_ptube(pid, mid, OD1, t=None, nsm=0., OD2=None,
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
                                  T1=1.0, T2=1.0, T3=1.0,
                                  comment='')

        # mass = (rho*t + nsm)*A = (0.2*0.5 + 0.3) * 0.5 = 0.4 * 0.5 = 0.2
        model.add_pshell(pid, mid1=None, t=0.5, mid2=mid, twelveIt3=1.0,
                         mid3=None, tst=0.833333,
                         nsm=0.3, z1=None, z2=None,
                         mid4=None, comment='')


        eid = 4
        pid = 4
        nids = [1, 2, 3, 4]
        cquad4 = model.add_cquad4(eid, pid, nids, theta_mcid=0.0, zoffset=0.,
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
        pcompg = model.add_pcompg(pid, global_ply_ids, mids, thicknesses, thetas=None,
                                  souts=None, nsm=0.0, sb=0.0, ft=None, tref=0.0, ge=0.0,
                                  lam=None, z0=None, comment='pcompg')
        model.add_pcompg(pid, global_ply_ids, mids, thicknesses, thetas=None,
                         souts=None, nsm=0.0, sb=0.0, ft=None, tref=0.0, ge=0.0,
                         lam=None, z0=None, comment='pcompg')


        pid = 40
        eid = 5
        nids = [1, 2, 3, 5]
        ctetra = model.add_ctetra(eid, pid, nids, comment='ctetra')

        eid = 6
        nids = [1, 2, 3, 4, 5]
        cpyram = model.add_cpyram(eid, pid, nids, comment='cpyram')

        eid = 7
        nids = [1, 2, 3, 5, 6, 7]
        cpenta = model.add_cpenta(eid, pid, nids, comment='cpenta')

        eid = 8
        nids = [1, 2, 3, 4, 5, 6, 7, 8]
        chexa = model.add_chexa(eid, pid, nids, comment='chexa')
        # mass = rho*V = 0.2*1

        psolid = model.add_psolid(pid, mid, cordm=0, integ=None, stress=None,
                                  isop=None, fctn='SMECH', comment='psolid')

        sid = 13
        eids = [eid]
        g1 = 6
        g34 = 8
        pressures = [1., 1., 1., 1.]
        pload4 = model.add_pload4(sid, eids, pressures, g1=1, g34=8,
                                  cid=0, nvector=None, surf_or_line='SURF',
                                  line_load_dir='NORM', comment='pload4')


        conid = 42
        gids = [
            2, 2, 2, 2, 2, 2,
            9, 9, 9, 9, 9, 9,
        ]
        components = [
            1, 2, 3, 4, 5, 6,
            1, 2, 3, 4, 5, 6,
        ]
        enforced = [
            1., 1., 1., 1., 1., 1.,
            1., 1., 1., 1., 1., 1.,
        ]
        mpc = model.add_mpc(conid, gids, components, enforced, comment='mpc')

        eid = 1
        ga = 9
        gb = 10
        cna = '123456'
        cnb = ''
        cma = ''
        cmb = ''
        rbar = model.add_rbar(eid, [ga, gb], cna, cnb, cma, cmb, alpha=0.,
                              comment='rbar')

        eid = 2
        ga = 10
        gb = 13
        rrod_a = RROD(eid, [ga, gb], cma='42', cmb='33')
        with self.assertRaises(RuntimeError):
            rrod_a.validate()
        rrod_b = model.add_rrod(eid, [ga, gb], cma='3', cmb=None, alpha=0.0, comment='')

        conid = 43
        gids = [10, 11]
        components = [1, 0]
        enforced = [1., 1.]
        mpc = model.add_mpc(conid, gids, components, enforced)
        model.add_spoint(11, comment='spoint')
        conid = 44
        sets = [42, 43]
        mpcadd = model.add_mpcadd(conid, sets, comment='mpcadd')
        #model.add_spoint([11, 'THRU', 42], comment='spoint3')
        str(model.spoints)

        sid = 14
        nids = 11
        mags = 20.
        sload = model.add_sload(sid, nids, mags, comment='an sload')

        sid = 12
        xyz = [2., 3., 4.]
        node = 7
        mag = 1.0
        force = model.add_force(sid, node, mag, xyz, cid=0, comment='force')
        moment = model.add_moment(sid, node, mag, xyz, comment='moment')

        node = 6
        mag = 1.0
        g1 = 2
        g2 = 3
        force1 = model.add_force1(sid, node, mag, g1, g2, comment='force1')
        moment1 = model.add_moment1(sid, node, mag, g1, g2, comment='moment1')

        g1 = 1
        g2 = 3
        g3 = 2
        g4 = 4
        force2 = model.add_force2(sid, node, mag, g1, g2, g3, g4, comment='force2')
        moment2 = model.add_moment2(sid, node, mag, g1, g2, g3, g4, comment='moment2')
        #g2, g3 = g3, g2
        #force2 = model.add_force2(sid, node, mag, g1, g2, g3, g4, comment='force2')

        load_id = 120
        scale = 1.
        scale_factors = [1.0, 2.0]
        load_ids = [12, 13]
        load = model.add_load(load_id, scale, scale_factors, load_ids, comment='load')

        #-----------------------------------------------------------------------
        # constraints
        conid = 42
        gids = [1, 2]
        components = ['123', '123']
        enforced = [0., 0.]
        spc = model.add_spc(conid, gids, components, enforced, comment='spc')
        conid = 43
        nodes = [1, 2]
        components = '123456'
        spc1 = model.add_spc1(conid, components, nodes, comment='spc1')
        conid = 44
        sets = [42, 43]
        spcadd = model.add_spcadd(conid, sets, comment='spcadd')
        #-----------------------------------------------------------------------
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
        F1, M1 = model2.sum_forces_moments_elements(p0, loadcase_id, eids, nids,
                                                    include_grav=False, xyz_cid0=None)
        F2, M2 = model2.sum_forces_moments(p0, loadcase_id, include_grav=False,
                                           xyz_cid0=None)
        assert allclose(F1, F2), 'F1=%s F2=%s' % (F1, F2)
        assert allclose(M1, M2), 'M1=%s M2=%s' % (M1, M2)

        model2.get_area_breakdown()
        model2.get_volume_breakdown()
        model2.get_mass_breakdown()

        model2.write_skin_solid_faces('skin.bdf', write_solids=False,
                                      write_shells=True)
        os.remove('skin.bdf')


if __name__ == '__main__':  # pragma: no cover
    unittest.main()

