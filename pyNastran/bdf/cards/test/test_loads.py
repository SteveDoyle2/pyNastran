from __future__ import print_function
import os
import unittest
from six import iteritems
from numpy import array, allclose, array_equal, set_printoptions
set_printoptions(suppress=True, precision=3)

import pyNastran
from pyNastran.bdf.bdf import BDF, BDFCard, DAREA, PLOAD4
from pyNastran.bdf.bdf import SET1, AESTAT, DMI, DMIG
from pyNastran.op2.op2 import OP2

bdf = BDF(debug=False)
test_path = pyNastran.__path__[0]


class TestLoads(unittest.TestCase):
    def test_darea_01(self):
        #
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
        lines = ['PLOAD4  1000    1       -60.    -60.    60.             1']
        card = bdf.process_card(lines)
        cardi = BDFCard(card)

        size = 8
        card = PLOAD4.add_card(cardi)
        card.write_card(size, 'dummy')
        card.raw_fields()

    def test_pload4_02(self):
        lines = ['PLOAD4  1       101     1.                              10000   10011']
        card = bdf.process_card(lines)
        cardi = BDFCard(card)

        size = 8
        card = PLOAD4.add_card(cardi)
        card.write_card(size, 'dummy')
        card.raw_fields()

    def test_pload4_cpenta(self):
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
                #print(load)
                face, area, centroid, normal = elem.getFaceAreaCentroidNormal(g1)
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
                face, area, centroid, normal = elem.getFaceAreaCentroidNormal(g1, g34)
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

    def test_pload4_cquad4(self):
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
            face, area, centroid, normal = elem.getFaceAreaCentroidNormal(g1, g34)
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
                print('%-2i Fx g=(%s,%s) f=%s fexpected=%s face=%s normal=%s' % (isubcase, g1, g34, f, fm, face, normal))
            if not allclose(f[1], fm[1]):
                print('%-2i Fy g=(%s,%s) f=%s fexpected=%s face=%s normal=%s' % (isubcase, g1, g34, f, fm, face, normal))
            if not allclose(f[2], fm[2]):
                print('%-2i Fz g=(%s,%s) f=%s fexpected=%s face=%s normal=%s' % (isubcase, g1, g34, f, fm, face, normal))

    def test_pload4_chexa(self):
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
            face, area, centroid, normal = elem.getFaceAreaCentroidNormal(g1, g34)
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
                raise RuntimeError('??? g1=%s g34=%s face=%s normal=%s\n%s' % (g1, g34, face, normal, msg))
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

            msg = '%s%s\n' % (elem.nodes[0], elem.nodes[1])

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

    def test_aestat_01(self):
        lines = ['AESTAT  502     PITCH']
        card = bdf.process_card(lines)
        card = BDFCard(card)

        size = 8
        card = AESTAT.add_card(card)
        card.write_card(size, 'dummy')
        card.raw_fields()

    def test_dmi_01(self):
        lines = ['DMI,Q,0,6,1,0,,4,4']
        card = bdf.process_card(lines)
        card = BDFCard(card)

        size = 8
        card = DMI(card)
        card.write_card(size, 'dummy')
        #card.rawFields()

    def test_set1_01(self):
        lines = ['SET1,    1100,    100,     101']
        card = bdf.process_card(lines)
        card = BDFCard(card)

        size = 8
        card = SET1.add_card(card)
        card.write_card(size, 'dummy')
        card.raw_fields()

        card2 = SET1(1100, [100, 101], is_skin=False, comment='')
        card2.write_card(size, 'dummy')

if __name__ == '__main__':  # pragma: no cover
    unittest.main()

