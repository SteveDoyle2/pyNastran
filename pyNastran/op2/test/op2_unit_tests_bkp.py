from __future__ import print_function
import os
import unittest
import numpy as np
from six import iteritems
from numpy import dot, array_equal

import pyNastran
test_path = pyNastran.__path__[0]

from pyNastran.bdf.bdf import BDF
from pyNastran.op2.op2 import OP2, FatalError, read_op2
from pyNastran.op2.op2_geom import read_op2_geom
from pyNastran.op2.test.test_op2 import run_op2

from pyNastran.bdf.test.bdf_unit_tests import Tester
from pyNastran.op2.tables.oef_forces.oef_force_objects import RealPlateBilinearForceArray, RealPlateForceArray
from pyNastran.op2.export_to_vtk import export_to_vtk_filename
from pyNastran.op2.vector_utils import filter1d


class TestOP2(Tester):
    #def _spike(self):
        #op2 = OP2()
        #op2.set_results('solidStress.oxx')
        #op2.read_op2(op2_filename, vectorized=False)

    def test_filter1d(self):
        a = np.array([1., 2., 0.1])
        i = filter1d(a, zero_tol=0.5)
        res = np.array([0, 1])
        self.assertTrue(np.array_equal(i, res), 'A i=%s res=%s' % (i, res))

        a = np.array([1., 2., 0.1])
        b = np.array([1., -0.1, 0.1])
        res = np.array([0, 1])
        i = filter1d(a, b, zero_tol=0.5)
        self.assertTrue(np.array_equal(i, res), 'B i=%s res=%s' % (i, res))

        a = np.array([1., 2., 0.1])
        b = np.array([1., -0.1, 0.1])
        i = filter1d(a, b, zero_tol=1.1)
        res = np.array([1])
        self.assertTrue(np.array_equal(i, res), 'C i=%s res=%s' % (i, res))

    def test_ibulk(self):
        """this test will fail if IBULK talble doesn't work"""
        bdf_filename = os.path.abspath(os.path.join(test_path, 'op2', 'test',
            'examples', 'ibulk', 'model1_sim1-solution_1.op2'))
        op2_filename = os.path.abspath(os.path.join(test_path, 'op2', 'test',
            'examples', 'ibulk', 'model1_sim1-solution_1.op2'))
        op2 = read_op2(op2_filename, debug=False)

    def test_set_results(self):
        folder = os.path.abspath(os.path.join(test_path, '..', 'models'))
        op2_filename = os.path.join(folder, 'solid_bending', 'solid_bending.op2')

        op2 = OP2(debug=False)
        op2.set_results('stress')
        op2.read_op2(op2_filename)
        self.assertEqual(len(op2.cpenta_stress), 0, len(op2.cpenta_stress))
        self.assertEqual(len(op2.chexa_stress), 0, len(op2.chexa_stress))
        self.assertEqual(len(op2.ctetra_stress), 1, len(op2.ctetra_stress))
        self.assertEqual(len(op2.displacements), 0, len(op2.displacements))

        op2 = OP2(debug=False)
        op2.set_results(['stress', 'displacements'])
        op2.read_op2(op2_filename)
        self.assertEqual(len(op2.cpenta_stress), 0, len(op2.cpenta_stress))
        self.assertEqual(len(op2.chexa_stress), 0, len(op2.chexa_stress))
        self.assertEqual(len(op2.ctetra_stress), 1, len(op2.ctetra_stress))
        self.assertEqual(len(op2.displacements), 1, len(op2.displacements))

    def test_op2_solid_bending_01(self):
        folder = os.path.abspath(os.path.join(test_path, '..', 'models', 'solid_bending'))
        op2_filename = os.path.join(folder, 'solid_bending.op2')
        make_geom = False
        write_bdf = False
        write_f06 = True
        debug = False
        #debug_file = 'solid_bending.debug.out'
        model, ext = os.path.splitext(op2_filename)
        debug_file = model + '.debug.out'

        if os.path.exists(debug_file):
            os.remove(debug_file)

        read_op2(op2_filename)
        run_op2(op2_filename, write_bdf=write_bdf, isubcases=[],
                write_f06=write_f06,
                debug=debug, stop_on_failure=True, binary_debug=True, quiet=True)
        assert os.path.exists(debug_file), os.listdir(folder)
        #os.remove(debug_file)

        make_geom = False
        write_bdf = False
        write_f06 = True
        run_op2(op2_filename, make_geom=make_geom, write_bdf=write_bdf, isubcases=[],
                write_f06=write_f06,
                debug=debug, stop_on_failure=True, binary_debug=True, quiet=True)
        assert os.path.exists(debug_file), os.listdir(folder)
        os.remove(debug_file)

    def test_op2_solid_bending_02(self):
        folder = os.path.abspath(os.path.join(test_path, '..', 'models', 'solid_bending'))
        op2_filename = os.path.join(folder, 'solid_bending.op2')
        op2 = read_op2(op2_filename, debug=False)

    def test_op2_solid_bending_02_geom(self):
        folder = os.path.abspath(os.path.join(test_path, '..', 'models', 'solid_bending'))
        op2_filename = os.path.join(folder, 'solid_bending.op2')
        op2 = read_op2_geom(op2_filename, debug=False)

    def test_op2_solid_shell_bar_01_geom(self):
        folder = os.path.abspath(os.path.join(test_path, '..', 'models', 'sol_101_elements'))
        op2_filename = os.path.join(folder, 'static_solid_shell_bar.op2')
        op2 = read_op2_geom(op2_filename, debug=False)


    def test_op2_solid_shell_bar_01_gpforce(self):
        folder = os.path.abspath(os.path.join(test_path, '..', 'models', 'sol_101_elements'))
        #bdf_filename = os.path.join(folder, 'static_solid_shell_bar.bdf')
        op2_filename = os.path.join(folder, 'static_solid_shell_bar.op2')
        op2 = read_op2_geom(op2_filename, debug=False)

        i_transform, beta_transforms = op2.get_displacement_index_transforms()
        op2.transform_displacements_to_global(i_transform, beta_transforms)

        gpforce = op2.grid_point_forces[1]

        #bdf_filename = os.path.join(folder, 'solid_shell_bar_xyz.bdf')
        #model = BDF(debug=False)
        #model.read_bdf(bdf_filename, xref=True)

        op2.cross_reference(xref_elements=False,
                            xref_nodes_with_elements=False,
                            xref_properties=False,
                            xref_masses=False,
                            xref_materials=False,
                            xref_loads=False,
                            xref_constraints=False,
                            xref_aero=False,
                            xref_sets=False,
                            xref_optimization=False)
        xyz_cid0 = op2.get_xyz_in_coord(cid=0)
        nid_cd = np.array([[nid, node.Cd()] for nid, node in sorted(iteritems(op2.nodes))])
        coords = op2.coords

        freebody_eids = [19]
        data = _get_gpforce_data()
        for datai in data:
            eids, nids, cid, summation_point, total_force_local_expected, total_moment_local_expected = datai
            if cid not in coords:
                continue
            op2.log.debug('*' * 30 + 'Next Test' + '*' * 30)
            coord_out = coords[cid]
            out = gpforce.extract_interface_loads(
                nids, eids,
                coord_out, coords,
                nid_cd, i_transform, beta_transforms,
                xyz_cid0, summation_point, itime=0, debug=False, logger=op2.log)
            total_force_global, total_moment_global, total_force_local, total_moment_local = out

            out = gpforce.extract_freebody_loads(
                freebody_eids,
                coord_out, coords, nid_cd, i_transform, beta_transforms,
                itime=0, debug=True, logger=op2.log)
            freebody_force, freebody_moment = out
            #print(freebody_force)
            #continue

            op2.log.debug('***********')
            op2.log.debug('force = %s; %s' % (total_force_global, np.linalg.norm(total_force_global)))
            op2.log.debug('moment = %s; %s' % (total_moment_global, np.linalg.norm(total_moment_global)))

            case = 'eids=%s nids=%s cid=%s summation_point=%s' % (
                eids, nids, cid, summation_point)
            msg = '%s\ntotal_force_local_expected=%s total_force_local=%s' % (
                case, total_force_local_expected, total_force_local)
            self.assertTrue(np.allclose(total_force_local_expected, total_force_local, atol=0.005), msg), msg

            msg = '%s\ntotal_moment_local_expected=%s total_moment_local=%s' % (
                case, total_moment_local_expected, total_moment_local)
            self.assertTrue(np.allclose(total_moment_local_expected, total_moment_local, atol=0.005), msg), msg

    def test_op2_solid_shell_bar_01_gpforce_xyz(self):
        folder = os.path.abspath(os.path.join(test_path, '..', 'models', 'sol_101_elements'))
        bdf_filename = os.path.join(folder, 'solid_shell_bar_xyz.bdf')
        op2_filename = os.path.join(folder, 'solid_shell_bar_xyz.op2')
        op2 = read_op2_geom(op2_filename)

        i_transform, beta_transforms = op2.get_displacement_index_transforms()
        op2.transform_displacements_to_global(i_transform, beta_transforms)

        gpforce = op2.grid_point_forces[1]
        op2.cross_reference(xref_elements=False,
                            xref_nodes_with_elements=False,
                            xref_properties=False,
                            xref_masses=False,
                            xref_materials=False,
                            xref_loads=False,
                            xref_constraints=False,
                            xref_aero=False,
                            xref_sets=False,
                            xref_optimization=False)
        xyz_cid0 = op2.get_xyz_in_coord(cid=0)
        nid_cd = np.array([[nid, node.Cd()] for nid, node in sorted(iteritems(op2.nodes))])

        model = BDF(debug=False)
        model.read_bdf(bdf_filename)

        data = _get_gpforce_data()
        coords = op2.coords
        for datai in data:
            eids, nids, cid, summation_point, total_force_local_expected, total_moment_local_expected = datai
            if cid not in coords:
                continue
            coord_out = coords[cid]
            op2.log.debug('*' * 30 + 'Next Test' + '*' * 30)
            out = gpforce.extract_interface_loads(
                nids, eids,
                coord_out, coords,
                nid_cd, i_transform, beta_transforms,
                xyz_cid0, summation_point, itime=0, debug=False, logger=op2.log)
            total_force_global, total_moment_global, total_force_local, total_moment_local = out

            op2.log.debug('***********')
            op2.log.debug('force = %s; %s' % (total_force_global, np.linalg.norm(total_force_global)))
            op2.log.debug('moment = %s; %s' % (total_moment_global, np.linalg.norm(total_moment_global)))

            case = 'eids=%s nids=%s cid=%s summation_point=%s' % (
                eids, nids, cid, summation_point)
            msg = '%s\ntotal_force_local_expected=%s total_force_local=%s delta=%s' % (
                case, total_force_local_expected, total_force_local,
                np.abs(total_force_local_expected - total_force_local))
            self.assertTrue(np.allclose(total_force_local_expected, total_force_local, atol=0.2), msg), msg

            msg = '%s\ntotal_moment_local_expected=%s total_moment_local=%s delta=%s' % (
                case, total_moment_local_expected, total_moment_local,
                np.abs(total_moment_local_expected - total_moment_local))
            self.assertTrue(np.allclose(total_moment_local_expected, total_moment_local, atol=0.005), msg), msg

    def test_op2_solid_shell_bar_01_gpforce_global_radial(self):
        folder = os.path.abspath(os.path.join(test_path, '..', 'models', 'sol_101_elements'))
        bdf_filename = os.path.join(folder, 'static_solid_shell_bar_global_radial.bdf')
        op2_filename = os.path.join(folder, 'static_solid_shell_bar_global_radial.op2')
        op2 = read_op2_geom(op2_filename)

        i_transform, beta_transforms = op2.get_displacement_index_transforms()
        op2.transform_displacements_to_global(i_transform, beta_transforms)

        gpforce = op2.grid_point_forces[1]
        op2.cross_reference(xref_elements=False,
                            xref_nodes_with_elements=False,
                            xref_properties=False,
                            xref_masses=False,
                            xref_materials=False,
                            xref_loads=False,
                            xref_constraints=False,
                            xref_aero=False,
                            xref_sets=False,
                            xref_optimization=False)
        xyz_cid0 = op2.get_xyz_in_coord(cid=0)
        nid_cd = np.array([[nid, node.Cd()] for nid, node in sorted(iteritems(op2.nodes))])

        # only uncommented entries are assumed to be correct
        data = [
            #eids, nids, cid, summation_point
            #[1], [1], 0, [0., 0., 0.],
            #[[1], [1, 2, 3, 4], 0, [0., 0., 0.], [0.0, 0.0, -10000.0], [-5000.0, 5000.0, 0.0],],

            # cid=0; eid=[1]; nid=[3]; sum=[0., 0., 0.] - done
            #               fmag     mmag       fx      fy       fz       mx       my       mz
            # F2      = [2589.95,     0.0, -12.59, -49.84, -2589.44,      0.0,     0.0,    0.0]  cid 0
            # F2Total = [2589.95, 3862.70,  26.34, -44.15, -2589.44, -2589.44, 2589.44, -70.49]  cid 0
            [[1], [3], 0, [0., 0., 0.], [26.34, -44.15, -2589.44], [-2589.44, 2589.44, -70.49],],

            # cid=0; eid=[1]; nid=[3]; sum=[0., 0., 0.] - done
            #               fmag     mmag       fx      fy       fz       mx       my       mz
            # F2      = [2589.95,     0.0, -12.59, -49.84, -2589.44,      0.0,     0.0,    0.0]  cid 0
            # F2Total = [2589.95, 3862.70,  26.34, -44.15, -2589.44, -2589.44, 2589.44, -70.49]  cid 1
            [[1], [3], 1, [0., 0., 0.], [44.15, -26.34, -2589.44], [2589.44, 2589.44, -70.49],],

            # cid=0; eid=[1]; nid=[3]; sum=[0., 0., 0.] - done
            #                            fx      fy       fz       mx       my       mz
            [[1], [1], 0, [0., 0., 0.], [-37.18, 32.00, -2589.44], [0.0, 0.0, 0.0],],

            # cid=0; eid=[1]; nid=[3]; sum=[0., 0., 0.] - done
            #                           fx      fy       fz       mx       my       mz
            [[1], [1], 1, [0., 0., 0.], [-32.00, 37.18, -2589.44], [0.0, 0.0, 0.0],],

            # conclusion from eid=1/nid=1/cid=0 and eid=1/nid=1/cid=1...forces are output in the XYZ frame
            # not RTZ as specified by freebody spreadsheet.  The RTZ table header indicates it's a
            # cylindrical frame, not that is a cylindrical theta.
        ]

        coords = op2.coords
        for datai in data:
            eids, nids, cid, summation_point, total_force_local_expected, total_moment_local_expected = datai
            if cid not in coords:
                continue
            coord_out = coords[cid]
            op2.log.debug('*' * 30 + 'Next Test' + '*' * 30)
            out = gpforce.extract_interface_loads(
                nids, eids,
                coord_out, coords,
                nid_cd, i_transform, beta_transforms,
                xyz_cid0, summation_point, itime=0, debug=True, logger=op2.log)
            total_force_global, total_moment_global, total_force_local, total_moment_local = out

            op2.log.debug('***********')
            op2.log.debug('force = %s; %s' % (total_force_global, np.linalg.norm(total_force_global)))
            op2.log.debug('moment = %s; %s' % (total_moment_global, np.linalg.norm(total_moment_global)))

            case = 'eids=%s nids=%s cid=%s summation_point=%s' % (
                eids, nids, cid, summation_point)
            msg = '%s\ntotal_force_local_expected=%s total_force_local=%s delta=%s' % (
                case, total_force_local_expected, total_force_local,
                np.abs(total_force_local_expected - total_force_local))
            self.assertTrue(np.allclose(total_force_local_expected, total_force_local, atol=0.2), msg), msg

            msg = '%s\ntotal_moment_local_expected=%s total_moment_local=%s delta=%s' % (
                case, total_moment_local_expected, total_moment_local,
                np.abs(total_moment_local_expected - total_moment_local))
            self.assertTrue(np.allclose(total_moment_local_expected, total_moment_local, atol=0.005), msg), msg

    def test_op2_solid_shell_bar_01_gpforce_global2_radial(self):
        folder = os.path.abspath(os.path.join(test_path, '..', 'models', 'sol_101_elements'))
        bdf_filename = os.path.join(folder, 'static_solid_shell_bar_global2_radial.bdf')
        op2_filename = os.path.join(folder, 'static_solid_shell_bar_global2_radial.op2')
        op2 = read_op2_geom(op2_filename)

        i_transform, beta_transforms = op2.get_displacement_index_transforms()
        op2.transform_displacements_to_global(i_transform, beta_transforms)

        gpforce = op2.grid_point_forces[1]
        op2.cross_reference(xref_elements=False,
                            xref_nodes_with_elements=False,
                            xref_properties=False,
                            xref_masses=False,
                            xref_materials=False,
                            xref_loads=False,
                            xref_constraints=False,
                            xref_aero=False,
                            xref_sets=False,
                            xref_optimization=False)
        xyz_cid0 = op2.get_xyz_in_coord(cid=0)
        nid_cd = np.array([[nid, node.Cd()] for nid, node in sorted(iteritems(op2.nodes))])

        # only uncommented entries are assumed to be correct
        data = [
            #eids, nids, cid, summation_point
            #[1], [1], 0, [0., 0., 0.],
            #[[1], [1, 2, 3, 4], 0, [0., 0., 0.], [0.0, 0.0, -10000.0], [-5000.0, 5000.0, 0.0],],

            # cid=0; eid=[1]; nid=[3]; sum=[0., 0., 0.] - done
            #               fmag     mmag       fx      fy       fz       mx       my       mz
            # F2      = [2589.95,     0.0, -12.59, -49.84, -2589.44,      0.0,     0.0,    0.0]  cid 0
            # F2Total = [2589.95, 3862.70,  26.34, -44.15, -2589.44, -2589.44, 2589.44, -70.49]  cid 0
            [[1], [3], 0, [0., 0., 0.], [26.34, -44.15, -2589.44], [-2589.44, 2589.44, -70.49],],

            # cid=0; eid=[1]; nid=[3]; sum=[0., 0., 0.] - done
            #               fmag     mmag       fx      fy       fz       mx       my       mz
            # F2      = [2589.95,     0.0, -12.59, -49.84, -2589.44,      0.0,     0.0,    0.0]  cid 0
            # F2Total = [2589.95, 3862.70,  26.34, -44.15, -2589.44, -2589.44, 2589.44, -70.49]  cid 1
            [[1], [3], 1, [0., 0., 0.], [-26.34, -44.15, 2589.44], [2589.44, 2589.44, -70.49],],
        ]

        coords = op2.coords
        for datai in data:
            eids, nids, cid, summation_point, total_force_local_expected, total_moment_local_expected = datai
            if cid not in coords:
                continue
            coord_out = coords[cid]
            op2.log.debug('*' * 30 + 'Next Test' + '*' * 30)
            out = gpforce.extract_interface_loads(
                nids, eids,
                coord_out, coords,
                nid_cd, i_transform, beta_transforms,
                xyz_cid0, summation_point, itime=0, debug=True, logger=op2.log)
            total_force_global, total_moment_global, total_force_local, total_moment_local = out

            op2.log.debug('***********')
            op2.log.debug('force = %s; %s' % (total_force_global, np.linalg.norm(total_force_global)))
            op2.log.debug('moment = %s; %s' % (total_moment_global, np.linalg.norm(total_moment_global)))

            case = 'eids=%s nids=%s cid=%s summation_point=%s' % (
                eids, nids, cid, summation_point)
            msg = '%s\ntotal_force_local_expected=%s total_force_local=%s delta=%s' % (
                case, total_force_local_expected, total_force_local,
                np.abs(total_force_local_expected - total_force_local))
            self.assertTrue(np.allclose(total_force_local_expected, total_force_local, atol=0.2), msg), msg

            msg = '%s\ntotal_moment_local_expected=%s total_moment_local=%s delta=%s' % (
                case, total_moment_local_expected, total_moment_local,
                np.abs(total_moment_local_expected - total_moment_local))
            self.assertTrue(np.allclose(total_moment_local_expected, total_moment_local, atol=0.005), msg), msg

    def test_op2_solid_shell_bar_01_gpforce_global_radial_cd(self):
        folder = os.path.abspath(os.path.join(test_path, '..', 'models', 'sol_101_elements'))
        bdf_filename = os.path.join(folder, 'static_solid_shell_bar_global_radial_cd.bdf')
        op2_filename = os.path.join(folder, 'static_solid_shell_bar_global_radial_cd.op2')
        op2 = read_op2_geom(op2_filename)

        i_transform, beta_transforms = op2.get_displacement_index_transforms()
        op2.transform_displacements_to_global(i_transform, beta_transforms)

        gpforce = op2.grid_point_forces[1]
        op2.cross_reference(xref_elements=False,
                            xref_nodes_with_elements=False,
                            xref_properties=False,
                            xref_masses=False,
                            xref_materials=False,
                            xref_loads=False,
                            xref_constraints=False,
                            xref_aero=False,
                            xref_sets=False,
                            xref_optimization=False)
        xyz_cid0 = op2.get_xyz_in_coord(cid=0)
        nid_cd = np.array([[nid, node.Cd()] for nid, node in sorted(iteritems(op2.nodes))])

        # only uncommented entries are assumed to be correct
        data = [
            #eids, nids, cid, summation_point
            #[1], [1], 0, [0., 0., 0.],
            #[[1], [1, 2, 3, 4], 0, [0., 0., 0.], [0.0, 0.0, -10000.0], [-5000.0, 5000.0, 0.0],],

            # cid=0; eid=[1]; nid=[3]; sum=[0., 0., 0.] - done
            #               fmag     mmag       fx      fy       fz       mx       my       mz
            # F2      = [2589.95,     0.0, -12.59, -49.84, -2589.44,      0.0,     0.0,    0.0]  cid 1
            # F2Total = [2589.95, 3862.70,  26.34, -44.15, -2589.44, -2589.44, 2589.44, -70.49]  cid 0
            [[1], [3], 0, [0., 0., 0.], [26.34, -44.15, -2589.44], [-2589.44, 2589.44, -70.49],],

            #[[3], [1], 0, [0., 0., 0.],],

            # cid=1; eid=[1]; nid=[1]; sum=[0., 0., 0.]
            #               fmag     mmag      fx      fy       fz       mx       my       mz
            # F1      = [2589.90,     0.0,1853.64, 567.74, 1717.35,     0.0,     0.0,     0.0]
            # F1Total = [2589.90,     0.0,1853.64, 567.74, 1717.35,     0.0,     0.0,     0.0]
            #[[1], [1], 1, [0., 0., 0.], [1853.64, 567.74, 1717.35], [0.0, 0.0, 0.0],],

            # cid=1; eid=[1]; nid=[2]; sum=[0., 0., 0.]
            #               fmag     mmag       fx      fy       fz       mx       my       mz
            # F2      = [2411.67,     0.0, 1710.67, 634.80, 1577.03,     0.0,     0.0,     0.0]
            # F2Total = [2411.67, 2410.58, 1710.67, 634.80, 1577.03, 1698.38, -570.22, -1612.84]
            #[[1], [2], 1, [0., 0., 0.], [1710.67, 634.80, 1577.03], [1698.38, -570.22, -1612.84],],

            # cid=1; eid=[1]; nid=[3]; sum=[0., 0., 0.]
            #           fmag          mmag     fx       fy       fz       mx        my       mz
            # F3      = [2589.95,     0.0, 1799.79, 645.58, 1746.94,     0.0,      0.0,     0.0]
            # F3Total = [2589.95, 3862.70, 1799.79, 645.58, 1746.94, 1880.85, -3035.07, -816.15]
            #[[1], [3], 1, [0., 0., 0.], [1799.79, 645.58, 1746.94], [1880.85, -3035.07, -816.15]],
        ]

        coords = op2.coords
        for datai in data:
            eids, nids, cid, summation_point, total_force_local_expected, total_moment_local_expected = datai
            if cid not in coords:
                continue
            coord_out = coords[cid]
            op2.log.debug('*' * 30 + 'Next Test' + '*' * 30)
            out = gpforce.extract_interface_loads(
                nids, eids,
                coord_out, coords,
                nid_cd, i_transform, beta_transforms,
                xyz_cid0, summation_point, itime=0, debug=True, logger=op2.log)
            total_force_global, total_moment_global, total_force_local, total_moment_local = out

            op2.log.debug('***********')
            op2.log.debug('force = %s; %s' % (total_force_global, np.linalg.norm(total_force_global)))
            op2.log.debug('moment = %s; %s' % (total_moment_global, np.linalg.norm(total_moment_global)))

            case = 'eids=%s nids=%s cid=%s summation_point=%s' % (
                eids, nids, cid, summation_point)
            msg = '%s\ntotal_force_local_expected=%s total_force_local=%s delta=%s' % (
                case, total_force_local_expected, total_force_local,
                np.abs(total_force_local_expected - total_force_local))
            self.assertTrue(np.allclose(total_force_local_expected, total_force_local, atol=0.2), msg), msg

            msg = '%s\ntotal_moment_local_expected=%s total_moment_local=%s delta=%s' % (
                case, total_moment_local_expected, total_moment_local,
                np.abs(total_moment_local_expected - total_moment_local))
            self.assertTrue(np.allclose(total_moment_local_expected, total_moment_local, atol=0.005), msg), msg

    @unittest.expectedFailure
    def test_op2_solid_shell_bar_01_gpforce_radial(self):
        folder = os.path.abspath(os.path.join(test_path, '..', 'models', 'sol_101_elements'))
        #bdf_filename = os.path.join(folder, 'static_solid_shell_bar_radial.bdf')
        op2_filename = os.path.join(folder, 'static_solid_shell_bar_radial.op2')
        op2 = read_op2_geom(op2_filename)

        i_transform, beta_transforms = op2.get_displacement_index_transforms()
        op2.transform_displacements_to_global(i_transform, beta_transforms)

        gpforce = op2.grid_point_forces[1]

        op2.cross_reference(xref_elements=False,
                            xref_nodes_with_elements=False,
                            xref_properties=False,
                            xref_masses=False,
                            xref_materials=False,
                            xref_loads=False,
                            xref_constraints=False,
                            xref_aero=False,
                            xref_sets=False,
                            xref_optimization=False)
        xyz_cid0 = op2.get_xyz_in_coord(cid=0)

        nid_cd = np.array([[nid, node.Cd()] for nid, node in sorted(iteritems(op2.nodes))])

        data = _get_gpforce_data()
        coords = op2.coords
        for i, datai in enumerate(data):
            eids, nids, cid, summation_point, total_force_local_expected, total_moment_local_expected = datai
            if cid not in coords:
                continue
            coord_out = coords[cid]
            op2.log.debug('*' * 30 + 'Next Test #%s' % i + '*' * 30)
            out = gpforce.extract_interface_loads(
                nids, eids,
                coord_out, coords,
                nid_cd, i_transform, beta_transforms,
                xyz_cid0, summation_point, itime=0, debug=False, logger=op2.log)
            total_force_global, total_moment_global, total_force_local, total_moment_local = out

            op2.log.debug('***********')
            op2.log.debug('force = %s; %s' % (total_force_global, np.linalg.norm(total_force_global)))
            op2.log.debug('moment = %s; %s' % (total_moment_global, np.linalg.norm(total_moment_global)))

            case = 'eids=%s nids=%s cid=%s summation_point=%s' % (
                eids, nids, cid, summation_point)
            msg = '%s\ntotal_force_local_expected=%s total_force_local=%s delta=%s' % (
                case, total_force_local_expected, total_force_local,
                np.abs(total_force_local_expected - total_force_local))
            self.assertTrue(np.allclose(total_force_local_expected, total_force_local, atol=0.2), msg), msg

            msg = '%s\ntotal_moment_local_expected=%s total_moment_local=%s delta=%s' % (
                case, total_moment_local_expected, total_moment_local,
                np.abs(total_moment_local_expected - total_moment_local))
            self.assertTrue(np.allclose(total_moment_local_expected, total_moment_local, atol=0.005), msg), msg

    def test_op2_solid_shell_bar_01(self):
        op2_filename = os.path.join('static_solid_shell_bar.op2')
        folder = os.path.abspath(os.path.join(test_path, '..', 'models', 'sol_101_elements'))
        op2_filename = os.path.join(folder, op2_filename)
        make_geom = False
        write_bdf = False
        write_f06 = True
        debug = False
        #debug_file = 'solid_bending.debug.out'
        model, ext = os.path.splitext(op2_filename)
        debug_file = model + '.debug.out'

        if os.path.exists(debug_file):
            os.remove(debug_file)
        read_op2(op2_filename)
        op2, is_passed = run_op2(op2_filename, make_geom=make_geom, write_bdf=write_bdf, isubcases=[],
                                 write_f06=write_f06,
                                 debug=debug, stop_on_failure=True, binary_debug=True, quiet=True)

        isubcase = 1
        rod_force = op2.crod_force[isubcase]
        rod_force.build_dataframe()
        assert rod_force.nelements == 2, rod_force.nelements
        assert rod_force.data.shape == (1, 2, 2), rod_force.data.shape

        rod_stress = op2.crod_stress[isubcase]
        rod_stress.build_dataframe()
        assert rod_stress.nelements == 2, rod_stress.nelements
        assert rod_stress.data.shape == (1, 2, 4), rod_stress.data.shape

        cbar_force = op2.cbar_force[isubcase]
        cbar_force.build_dataframe()
        assert cbar_force.nelements == 1, cbar_force.nelements
        assert cbar_force.data.shape == (1, 1, 8), cbar_force.data.shape

        cbar_stress = op2.cbar_stress[isubcase]
        cbar_stress.build_dataframe()
        assert cbar_stress.nelements == 1, cbar_stress.nelements
        assert cbar_stress.data.shape == (1, 1, 15), cbar_stress.data.shape

        cbeam_force = op2.cbeam_force[isubcase]
        cbeam_force.build_dataframe()
        assert cbeam_force.nelements == 1, cbeam_force.nelements
        assert cbeam_force.data.shape == (1, 2, 8), cbeam_force.data.shape

        cbeam_stress = op2.cbeam_stress[isubcase]
        cbeam_stress.build_dataframe()
        assert cbeam_stress.nelements == 11, cbeam_stress.nelements  # wrong
        assert cbeam_stress.data.shape == (1, 2, 8), cbeam_stress.data.shape

        cquad4_force = op2.cquad4_force[isubcase]
        cquad4_force.build_dataframe()
        assert cquad4_force.nelements == 4, cquad4_force.nelements
        assert cquad4_force.data.shape == (1, 20, 8), cquad4_force.data.shape

        cquad4_stress = op2.cquad4_stress[isubcase]
        cquad4_stress.build_dataframe()
        assert cquad4_stress.nelements == 20, cquad4_stress.nelements # TODO: should this be 4; yes by actual count...
        assert cquad4_stress.data.shape == (1, 20, 8), cquad4_stress.data.shape
        assert cquad4_stress.is_fiber_distance(), cquad4_stress
        assert cquad4_stress.is_von_mises(), cquad4_stress

        ctria3_force = op2.ctria3_force[isubcase]
        ctria3_force.build_dataframe()
        assert ctria3_force.nelements == 8, ctria3_force.nelements
        assert ctria3_force.data.shape == (1, 8, 8), ctria3_force.data.shape

        ctria3_stress = op2.ctria3_stress[isubcase]
        ctria3_stress.build_dataframe()
        assert ctria3_stress.nelements == 8, ctria3_stress.nelements
        assert ctria3_stress.data.shape == (1, 8, 8), ctria3_stress.data.shape
        assert ctria3_stress.is_fiber_distance(), ctria3_stress
        assert ctria3_stress.is_von_mises(), ctria3_stress

        ctetra_stress = op2.ctetra_stress[isubcase]
        ctetra_stress.build_dataframe()
        assert ctetra_stress.nelements == 2, ctetra_stress.nelements
        assert ctetra_stress.data.shape == (1, 10, 10), ctetra_stress.data.shape
        assert ctetra_stress.is_von_mises(), ctetra_stress

        cpenta_stress = op2.cpenta_stress[isubcase]
        cpenta_stress.build_dataframe()
        assert cpenta_stress.nelements == 2, cpenta_stress.nelements
        assert cpenta_stress.data.shape == (1, 14, 10), cpenta_stress.data.shape
        assert cpenta_stress.is_von_mises(), cpenta_stress

        chexa_stress = op2.chexa_stress[isubcase]
        chexa_stress.build_dataframe()
        assert chexa_stress.nelements == 1, chexa_stress.nelements
        assert chexa_stress.data.shape == (1, 9, 10), chexa_stress.data.shape
        assert chexa_stress.is_von_mises(), chexa_stress

        assert os.path.exists(debug_file), os.listdir(folder)
        os.remove(debug_file)

    def test_op2_solid_shell_bar_01_export(self):
        folder = os.path.abspath(os.path.join(test_path, '..', 'models', 'sol_101_elements'))
        bdf_filename = os.path.join(folder, 'static_solid_shell_bar.bdf')
        op2_filename = os.path.join(folder, 'static_solid_shell_bar.op2')
        vtk_filename = os.path.join(folder, 'static_solid_shell_bar.vtk')
        export_to_vtk_filename(bdf_filename, op2_filename, vtk_filename)

    def test_op2_solid_shell_bar_02(self):
        op2_filename = os.path.join('mode_solid_shell_bar.op2')
        folder = os.path.abspath(os.path.join(test_path, '..', 'models', 'sol_101_elements'))
        op2_filename = os.path.join(folder, op2_filename)
        make_geom = False
        write_bdf = False
        write_f06 = True
        debug = False
        #debug_file = 'solid_bending.debug.out'
        model, ext = os.path.splitext(op2_filename)
        debug_file = model + '.debug.out'

        if os.path.exists(debug_file):
            os.remove(debug_file)
        read_op2(op2_filename, debug=False)
        op2, is_passed = run_op2(op2_filename, make_geom=make_geom, write_bdf=write_bdf, isubcases=[],
                                 write_f06=write_f06,
                                 debug=debug, stop_on_failure=True, binary_debug=True, quiet=True)

        isubcase = 1
        rod_force = op2.crod_force[isubcase]
        assert rod_force.nelements == 2, rod_force.nelements
        assert rod_force.data.shape == (3, 2, 2), rod_force.data.shape

        rod_stress = op2.crod_stress[isubcase]
        assert rod_stress.nelements == 2, rod_stress.nelements
        assert rod_stress.data.shape == (3, 2, 4), rod_stress.data.shape

        cbar_force = op2.cbar_force[isubcase]
        cbar_force.build_dataframe()
        str(cbar_force.data_frame)
        assert cbar_force.nelements == 1, cbar_force.nelements
        assert cbar_force.data.shape == (3, 1, 8), cbar_force.data.shape

        cbar_stress = op2.cbar_stress[isubcase]
        assert cbar_stress.nelements == 3, cbar_stress.nelements  # TODO: wrong
        assert cbar_stress.data.shape == (3, 1, 15), cbar_stress.data.shape

        cbeam_stress = op2.cbeam_stress[isubcase]
        assert cbeam_stress.nelements == 11, cbeam_stress.nelements  # TODO: wrong
        assert cbeam_stress.data.shape == (3, 2, 8), cbeam_stress.data.shape

        cquad4_stress = op2.cquad4_stress[isubcase]
        assert cquad4_stress.nelements == 20, cquad4_stress.nelements
        assert cquad4_stress.data.shape == (3, 20, 8), cquad4_stress.data.shape

        ctria3_stress = op2.ctria3_stress[isubcase]
        assert ctria3_stress.nelements == 8, ctria3_stress.nelements
        assert ctria3_stress.data.shape == (3, 8, 8), ctria3_stress.data.shape

        ctetra_stress = op2.ctetra_stress[isubcase]
        assert ctetra_stress.nelements == 2, ctetra_stress.nelements
        assert ctetra_stress.data.shape == (3, 10, 10), ctetra_stress.data.shape

        cpenta_stress = op2.cpenta_stress[isubcase]
        assert cpenta_stress.nelements == 2, cpenta_stress.nelements
        assert cpenta_stress.data.shape == (3, 14, 10), cpenta_stress.data.shape

        chexa_stress = op2.chexa_stress[isubcase]
        assert chexa_stress.nelements == 1, chexa_stress.nelements
        assert chexa_stress.data.shape == (3, 9, 10), chexa_stress.data.shape

        assert os.path.exists(debug_file), os.listdir(folder)
        os.remove(debug_file)

    def test_op2_solid_shell_bar_02_export(self):
        folder = os.path.abspath(os.path.join(test_path, '..', 'models', 'sol_101_elements'))
        bdf_filename = os.path.join(folder, 'mode_solid_shell_bar.bdf')
        op2_filename = os.path.join(folder, 'mode_solid_shell_bar.op2')
        vtk_filename = os.path.join(folder, 'mode_solid_shell_bar.vtk')
        export_to_vtk_filename(bdf_filename, op2_filename, vtk_filename)

    def test_op2_solid_shell_bar_03(self):
        op2_filename = os.path.join('buckling_solid_shell_bar.op2')
        folder = os.path.abspath(os.path.join(test_path, '..', 'models', 'sol_101_elements'))
        op2_filename = os.path.join(folder, op2_filename)
        make_geom = False
        write_bdf = False
        write_f06 = True
        debug = False
        #debug_file = 'solid_bending.debug.out'
        model, ext = os.path.splitext(op2_filename)
        debug_file = model + '.debug.out'

        if os.path.exists(debug_file):
            os.remove(debug_file)
        read_op2(op2_filename, debug=False)
        op2, is_passed = run_op2(op2_filename, make_geom=make_geom, write_bdf=write_bdf, isubcases=[],
                                 write_f06=write_f06,
                                 debug=debug, stop_on_failure=True, binary_debug=True, quiet=True)

        isubcases = [(1, 1, 1, 0, 'DEFAULT1'), (1, 8, 1, 0, 'DEFAULT1')]
        isubcase = isubcases[1]

        try:
            rod_force = op2.crod_force[isubcase]
        except KeyError:
            msg = 'isubcase=%s was not found\nkeys=%s' % (isubcase, op2.crod_force.keys())
            raise KeyError(msg)
        assert rod_force.nelements == 2, rod_force.nelements
        assert rod_force.data.shape == (4, 2, 2), rod_force.data.shape

        rod_stress = op2.crod_stress[isubcase]
        assert rod_stress.nelements == 2, rod_stress.nelements
        assert rod_stress.data.shape == (4, 2, 4), rod_stress.data.shape

        cbar_force = op2.cbar_force[isubcase]
        cbar_force.build_dataframe()
        str(cbar_force.data_frame)
        assert cbar_force.nelements == 1, cbar_force.nelements
        assert cbar_force.data.shape == (4, 1, 8), cbar_force.data.shape

        cbar_stress = op2.cbar_stress[isubcase]
        assert cbar_stress.nelements == 4, cbar_stress.nelements  # TODO: wrong
        assert cbar_stress.data.shape == (4, 1, 15), cbar_stress.data.shape

        cbeam_stress = op2.cbeam_stress[isubcase]
        assert cbeam_stress.nelements == 11, cbeam_stress.nelements  # TODO: wrong
        assert cbeam_stress.data.shape == (4, 2, 8), cbeam_stress.data.shape

        cquad4_stress = op2.cquad4_stress[isubcase]
        assert cquad4_stress.nelements == 20, cquad4_stress.nelements
        assert cquad4_stress.data.shape == (4, 20, 8), cquad4_stress.data.shape

        ctria3_stress = op2.ctria3_stress[isubcase]
        assert ctria3_stress.nelements == 8, ctria3_stress.nelements
        assert ctria3_stress.data.shape == (4, 8, 8), ctria3_stress.data.shape

        ctetra_stress = op2.ctetra_stress[isubcase]
        assert ctetra_stress.nelements == 2, ctetra_stress.nelements
        assert ctetra_stress.data.shape == (4, 10, 10), ctetra_stress.data.shape

        cpenta_stress = op2.cpenta_stress[isubcase]
        assert cpenta_stress.nelements == 2, cpenta_stress.nelements
        assert cpenta_stress.data.shape == (4, 14, 10), cpenta_stress.data.shape

        chexa_stress = op2.chexa_stress[isubcase]
        assert chexa_stress.nelements == 1, chexa_stress.nelements
        assert chexa_stress.data.shape == (4, 9, 10), chexa_stress.data.shape

        assert os.path.exists(debug_file), os.listdir(folder)
        os.remove(debug_file)

    #@unittest.expectedFailure
    def test_op2_solid_shell_bar_04(self):
        op2_filename = os.path.join('freq_solid_shell_bar.op2')
        folder = os.path.abspath(os.path.join(test_path, '..', 'models', 'sol_101_elements'))
        op2_filename = os.path.join(folder, op2_filename)
        make_geom = False
        write_bdf = False
        write_f06 = True
        debug = False
        #debug_file = 'solid_bending.debug.out'
        model, ext = os.path.splitext(op2_filename)
        debug_file = model + '.debug.out'

        if os.path.exists(debug_file):
            os.remove(debug_file)
        read_op2(op2_filename, debug=False)
        op2, is_passed = run_op2(op2_filename, make_geom=make_geom, write_bdf=write_bdf, isubcases=[],
                                 write_f06=write_f06,
                                 debug=debug, stop_on_failure=True, binary_debug=True, quiet=True)
        isubcase = 1
        # rod_force = op2.crod_force[isubcase]
        # assert rod_force.nelements == 2, rod_force.nelements
        # assert rod_force.data.shape == (7, 2, 2), rod_force.data.shape


        # isubcases = [(1, 1, 1, 0, 'DEFAULT'), (1, 8, 1, 0, 'DEFAULT')]
        # isubcase = isubcases[1]

        rod_force = op2.crod_force[isubcase]
        rod_force.build_dataframe()
        assert rod_force.nelements == 2, rod_force.nelements
        assert rod_force.data.shape == (7, 2, 2), rod_force.data.shape

        rod_stress = op2.crod_stress[isubcase]
        rod_stress.build_dataframe()
        assert rod_stress.nelements == 2, rod_stress.nelements
        assert rod_stress.data.shape == (7, 2, 2), rod_stress.data.shape

        cbar_force = op2.cbar_force[isubcase]
        cbar_force.build_dataframe()
        assert cbar_force.nelements == 1, cbar_force.nelements
        assert cbar_force.data.shape == (7, 1, 8), cbar_force.data.shape

        cbar_stress = op2.cbar_stress[isubcase]
        cbar_stress.build_dataframe()
        assert cbar_stress.nelements == 1, cbar_stress.nelements
        assert cbar_stress.data.shape == (7, 1, 9), cbar_stress.data.shape

        #print(op2.cbeam_stress.keys())
        #cbeam_stress = op2.cbeam_stress[isubcase]
        #assert cbeam_stress.nelements == 2, cbeam_stress.nelements
        #assert cbeam_stress.data.shape == (7, 2, 8), cbeam_stress.data.shape

        cquad4_stress = op2.cquad4_stress[isubcase]
        #cquad4_stress.build_dataframe()
        assert cquad4_stress.nelements == 4, cquad4_stress.nelements # TODO: wrong
        assert cquad4_stress.data.shape == (7, 40, 3), cquad4_stress.data.shape

        #print(op2.ctria3_stress.keys())
        ctria3_stress = op2.ctria3_stress[isubcase]
        ctria3_stress.build_dataframe()
        assert ctria3_stress.nelements == 8, ctria3_stress.nelements # TODO: wrong
        assert ctria3_stress.data.shape == (7, 32, 3), ctria3_stress.data.shape

        ctetra_stress = op2.ctetra_stress[isubcase]
        ctetra_stress.build_dataframe()
        assert ctetra_stress.nelements == 2, ctetra_stress.nelements
        assert ctetra_stress.data.shape == (7, 10, 6), ctetra_stress.data.shape

        cpenta_stress = op2.cpenta_stress[isubcase]
        cpenta_stress.build_dataframe()
        assert cpenta_stress.nelements == 2, cpenta_stress.nelements
        assert cpenta_stress.data.shape == (7, 14, 6), cpenta_stress.data.shape

        chexa_stress = op2.chexa_stress[isubcase]
        chexa_stress.build_dataframe()
        assert chexa_stress.nelements == 1, chexa_stress.nelements
        assert chexa_stress.data.shape == (7, 9, 6), chexa_stress.data.shape

        grid_point_forces = op2.grid_point_forces[isubcase]
        grid_point_forces.build_dataframe()
        #print(grid_point_forces._ntotals)
        assert grid_point_forces.ntotal == 106, grid_point_forces.ntotal
        assert grid_point_forces.data.shape == (7, 106, 6), grid_point_forces.data.shape

        assert os.path.exists(debug_file), os.listdir(folder)
        os.remove(debug_file)

    def test_op2_solid_shell_bar_03_export(self):
        folder = os.path.abspath(os.path.join(test_path, '..', 'models', 'sol_101_elements'))
        bdf_filename = os.path.join(folder, 'freq_solid_shell_bar.bdf')
        op2_filename = os.path.join(folder, 'freq_solid_shell_bar.op2')
        vtk_filename = os.path.join(folder, 'freq_solid_shell_bar.vtk')
        export_to_vtk_filename(bdf_filename, op2_filename, vtk_filename)

    def test_op2_solid_shell_bar_05(self):
        """
        MSC 2005r2 Tables : GEOM1, GEOM2, GEOM3, GEOM4, EPT, MPTS, DYNAMICS, DIT
                            OQG1, OUGV1, OGPFB1, OEF1X, OES1X1, OSTR1X, OPG1
        NX 10 Tables : PVT0, CASECC, GEOM1S, GEOM2S, GEOM3S, GEOM4S, EPTS, MPTS,
                       DYNAMICS, BGPDTS, EQEXINS, DIT,
                       OQG1, OUGV1, OGPFB1, OEF1X, OES1X1, OSTR1X, OPG1
        """
        op2_filename = os.path.join('transient_solid_shell_bar.op2')
        folder = os.path.abspath(os.path.join(test_path, '..', 'models', 'sol_101_elements'))
        op2_filename = os.path.join(folder, op2_filename)
        make_geom = False
        write_bdf = False
        write_f06 = True
        debug = False
        #debug_file = 'solid_bending.debug.out'
        model, ext = os.path.splitext(op2_filename)
        debug_file = model + '.debug.out'

        if os.path.exists(debug_file):
            os.remove(debug_file)
        read_op2(op2_filename, debug=debug)
        op2, is_passed = run_op2(op2_filename, make_geom=make_geom, write_bdf=write_bdf, isubcases=[],
                                 write_f06=write_f06,
                                 debug=debug, stop_on_failure=True, binary_debug=True, quiet=True)
        isubcase = 1
        # rod_force = op2.crod_force[isubcase]
        # assert rod_force.nelements == 2, rod_force.nelements
        # assert rod_force.data.shape == (7, 2, 2), rod_force.data.shape


        # isubcases = [(1, 1, 1, 0, 'DEFAULT'), (1, 8, 1, 0, 'DEFAULT')]
        # isubcase = isubcases[1]

        rod_force = op2.crod_force[isubcase]
        rod_force.build_dataframe()
        assert rod_force.nelements == 2, rod_force.nelements
        assert rod_force.data.shape == (21, 2, 2), rod_force.data.shape

        rod_stress = op2.crod_stress[isubcase]
        rod_stress.build_dataframe()
        assert rod_stress.nelements == 2, rod_stress.nelements
        assert rod_stress.data.shape == (21, 2, 4), rod_stress.data.shape

        cbar_force = op2.cbar_force[isubcase]
        cbar_force.build_dataframe()
        assert cbar_force.nelements == 1, cbar_force.nelements
        assert cbar_force.data.shape == (21, 1, 8), cbar_force.data.shape

        cbar_stress = op2.cbar_stress[isubcase]
        cbar_stress.build_dataframe()
        assert cbar_stress.nelements == 21, cbar_stress.nelements # 1-wrong
        assert cbar_stress.data.shape == (21, 1, 15), cbar_stress.data.shape

        #print(op2.cbeam_stress.keys())
        # cbeam_stress = op2.cbeam_stress[isubcase]
        # assert cbeam_stress.nelements == 11, cbeam_stress.nelements  # TODO: wrong
        # assert cbeam_stress.data.shape == (7, 11, 8), cbeam_stress.data.shape

        cquad4_stress = op2.cquad4_stress[isubcase]
        #cquad4_stress.build_dataframe()
        assert cquad4_stress.nelements == 40, cquad4_stress.nelements # TODO: (840-wrong, 40-correct)
        assert cquad4_stress.data.shape == (21, 40, 8), cquad4_stress.data.shape

        #print(op2.ctria3_stress.keys())
        ctria3_stress = op2.ctria3_stress[isubcase]
        ctria3_stress.build_dataframe()
        assert ctria3_stress.nelements == 16, ctria3_stress.nelements # TODO: 8-wrong
        assert ctria3_stress.data.shape == (21, 16, 8), ctria3_stress.data.shape

        ctetra_stress = op2.ctetra_stress[isubcase]
        ctetra_stress.build_dataframe()
        assert ctetra_stress.nelements == 2, ctetra_stress.nelements
        assert ctetra_stress.data.shape == (21, 10, 10), ctetra_stress.data.shape

        cpenta_stress = op2.cpenta_stress[isubcase]
        cpenta_stress.build_dataframe()
        assert cpenta_stress.nelements == 2, cpenta_stress.nelements
        assert cpenta_stress.data.shape == (21, 14, 10), cpenta_stress.data.shape

        chexa_stress = op2.chexa_stress[isubcase]
        chexa_stress.build_dataframe()
        assert chexa_stress.nelements == 1, chexa_stress.nelements
        assert chexa_stress.data.shape == (21, 9, 10), chexa_stress.data.shape

        grid_point_forces = op2.grid_point_forces[isubcase]
        grid_point_forces.build_dataframe()
        #print(grid_point_forces._ntotals)
        assert grid_point_forces.ntotal == 130, grid_point_forces.ntotal
        assert grid_point_forces.data.shape == (21, 130, 6), grid_point_forces.data.shape

        assert os.path.exists(debug_file), os.listdir(folder)
        os.remove(debug_file)

    def test_op2_plate_py_01(self):
        op2_filename = os.path.join('plate_py', 'plate_py.op2')
        folder = os.path.abspath(os.path.join(test_path, '..', 'models'))
        make_geom = False
        write_bdf = False
        write_f06 = False
        debug = False
        op2file = os.path.join(folder, op2_filename)
        read_op2(op2file)
        run_op2(op2file, make_geom=make_geom, write_bdf=write_bdf, isubcases=[],
                write_f06=write_f06,
                debug=debug, stop_on_failure=True, quiet=True)

        make_geom = False
        write_bdf = False
        write_f06 = True
        run_op2(op2file, make_geom=make_geom, write_bdf=write_bdf, isubcases=[],
                write_f06=write_f06,
                debug=debug, stop_on_failure=True, quiet=True)

    def test_op2_good_sine_01(self):
        op2_filename = os.path.join('freq_sine', 'good_sine.op2')
        folder = os.path.abspath(os.path.join(test_path, '..', 'models'))
        make_geom = False
        write_bdf = False
        write_f06 = False
        debug = False
        op2file = os.path.join(folder, op2_filename)
        read_op2(op2file)
        op2i, is_passed = run_op2(op2file, make_geom=make_geom, write_bdf=write_bdf, isubcases=[],
                                  write_f06=write_f06,
                                  debug=debug, stop_on_failure=True,
                                  quiet=True)

        nids = [5]
        isubcase = 103
        try:
            acc = op2i.accelerations[isubcase]
        except KeyError:
            raise KeyError('getting accelerations; isubcase=%s; keys=%s' % (
                isubcase, op2i.accelerations.keys()))

        with self.assertRaises(AssertionError):
            # no index 0; fortran 1-based
            acc.extract_xyplot(nids, 0, 'real')

        acc.build_dataframe()
        accx = acc.extract_xyplot(nids, 1, 'real')
        accxi = acc.extract_xyplot(nids, 1, 'imag')
        #print(accx)
        #print(accxi)
        #make_geom = False
        #write_bdf = False
        #write_f06 = True
        #run_op2(op2file, make_geom=make_geom, write_bdf=write_bdf, iSubcases=[],
                #write_f06=write_f06,
                #debug=debug, stopOnFailure=True)

    def test_op2_good_sine_02(self):
        folder = os.path.abspath(os.path.join(test_path, '..', 'models'))
        bdf_filename = os.path.join(folder, 'freq_sine', 'good_sine.dat')
        op2_filename = os.path.join(folder, 'freq_sine', 'good_sine.op2')
        make_geom = False
        write_bdf = False
        write_f06 = True
        debug = False
        op2file = os.path.join(folder, op2_filename)
        bdf = BDF(debug=False)
        bdf.read_bdf(bdf_filename)

        debug = False
        debug_file = 'debug.out'

        read_op2(op2_filename)
        op2 = OP2(debug=debug, debug_file=debug_file)
        op2.read_op2(op2_filename)
        assert os.path.exists(debug_file), os.listdir('.')

        self._verify_ids(bdf, op2, isubcase=1)

    def test_op2_bcell_01(self):
        folder = os.path.abspath(os.path.join(test_path, '..', 'models'))
        bdf_filename = os.path.join(folder, 'other', 'bcell9p0.bdf')
        op2_filename = os.path.join(folder, 'other', 'bcell9p0.op2')
        make_geom = False
        write_bdf = False
        write_f06 = True
        debug = False
        op2file = os.path.join(folder, op2_filename)
        bdf = BDF(debug=False)
        bdf.read_bdf(bdf_filename, xref=False)

        debug = False
        debug_file = 'debug.out'
        read_op2(op2_filename)
        op2 = OP2(debug=debug, debug_file=debug_file)
        op2.read_op2(op2_filename)
        assert os.path.exists(debug_file), os.listdir('.')

        self._verify_ids(bdf, op2, isubcase=1)

        msg = ''
        for isubcase, keys in sorted(iteritems(op2.subcase_key)):
            if len(keys) != 1:
                msg += 'isubcase=%s keys=%s len(keys) != 1\n' % (isubcase, keys)
                if len(keys) == 0:
                    continue
            if isubcase != keys[0]:
                msg += 'isubcase=%s != key0=%s keys=%s\n' % (isubcase, keys[0], keys)
        if msg:
            assert msg == '', msg
        op2.write_f06('junk.f06', quiet=True)
        os.remove('junk.f06')

    def test_op2_cbush_01(self):
        op2_filename = os.path.join('cbush.op2')
        folder = os.path.abspath(os.path.join(test_path, '..', 'models', 'cbush'))
        op2_filename = os.path.join(folder, op2_filename)
        make_geom = False
        write_bdf = False
        write_f06 = True
        debug = False
        #debug_file = 'solid_bending.debug.out'
        model, ext = os.path.splitext(op2_filename)
        debug_file = model + '.debug.out'

        if os.path.exists(debug_file):
            os.remove(debug_file)
        read_op2(op2_filename)
        op2, is_passed = run_op2(op2_filename, make_geom=make_geom, write_bdf=write_bdf, isubcases=[],
                                 write_f06=write_f06,
                                 debug=debug, stop_on_failure=True, binary_debug=True, quiet=True)
        isubcase = 1

        cbush_stress = op2.cbush_stress[isubcase]
        cbush_stress.build_dataframe()
        assert cbush_stress.nelements == 1, cbush_stress.nelements
        assert cbush_stress.data.shape == (1, 1, 6), cbush_stress.data.shape

        cbush_strain = op2.cbush_strain[isubcase]
        cbush_strain.build_dataframe()
        assert cbush_strain.nelements == 1, cbush_strain.nelements
        assert cbush_strain.data.shape == (1, 1, 6), cbush_strain.data.shape

        cbush_force = op2.cbush_force[isubcase]
        cbush_force.build_dataframe()
        assert cbush_force.nelements == 1, cbush_force.nelements
        assert cbush_force.data.shape == (1, 1, 6), cbush_force.data.shape

        assert os.path.exists(debug_file), os.listdir(folder)
        os.remove(debug_file)

    def _verify_ids(self, bdf, op2, isubcase=1):
        types = ['CQUAD4', 'CTRIA3', 'CHEXA', 'CPENTA', 'CTETRA', 'CROD', 'CONROD', 'CTUBE']
        out = bdf.get_card_ids_by_card_types(types)

        card_type = 'CQUAD4'
        if op2.cquad4_stress:
            try:
                case = op2.cquad4_stress[isubcase]
            except KeyError:
                raise KeyError('getting cquad4_stress; isubcase=%s; keys=%s' % (
                    isubcase, op2.cquad4_stress.keys()))
            eids = np.unique(case.element_node[:, 0])
            for eid in eids:
                assert eid in out[card_type], 'eid=%s eids=%s card_type=%s'  % (eid, out[card_type], card_type)
        if op2.cquad4_strain:
            try:
                case = op2.cquad4_strain[isubcase]
            except KeyError:
                raise KeyError('getting cquad4_strain; isubcase=%s; keys=%s' % (
                    isubcase, op2.cquad4_strain.keys()))
            eids = np.unique(case.element_node[:, 0])
            for eid in eids:
                assert eid in out[card_type], 'eid=%s eids=%s card_type=%s'  % (eid, out[card_type], card_type)
        if op2.cquad4_composite_strain:
            try:
                case = op2.cquad4_composite_strain[isubcase]
            except KeyError:
                raise KeyError('getting cquad4_composite_strain; isubcase=%s; keys=%s' % (
                    isubcase, op2.cquad4_composite_strain.keys()))
            eids = np.unique(case.element_layer[:, 0])
            for eid in eids:
                assert eid in out[card_type], 'eid=%s eids=%s card_type=%s'  % (eid, out[card_type], card_type)
        if op2.cquad4_composite_stress:
            try:
                case = op2.cquad4_composite_stress[isubcase]
            except KeyError:
                raise KeyError('getting cquad4_composite_stress; isubcase=%s; keys=%s' % (
                    isubcase, op2.cquad4_composite_stress.keys()))

            eids = np.unique(case.element_layer[:, 0])
            for eid in eids:
                assert eid in out[card_type], 'eid=%s eids=%s card_type=%s'  % (eid, out[card_type], card_type)
        if op2.cquad4_force:
            try:
                case = op2.cquad4_force[isubcase]
            except KeyError:
                raise KeyError('getting cquad4_force; isubcase=%s; keys=%s' % (
                    isubcase, op2.cquad4_force.keys()))

            eids = np.unique(case.element)
            for eid in eids:
                assert eid in out[card_type], 'eid=%s eids=%s card_type=%s'  % (eid, out[card_type], card_type)


        card_type = 'CTRIA3'
        if op2.ctria3_stress:
            case = op2.ctria3_stress[isubcase]
            eids = np.unique(case.element_node[:, 0])
            for eid in eids:
                assert eid in out[card_type], 'eid=%s eids=%s card_type=%s'  % (eid, out[card_type], card_type)
        if op2.ctria3_strain:
            case = op2.ctria3_strain[isubcase]
            eids = np.unique(case.element_node[:, 0])
            for eid in eids:
                assert eid in out[card_type], 'eid=%s eids=%s card_type=%s'  % (eid, out[card_type], card_type)
        if op2.ctria3_composite_strain:
            case = op2.ctria3_composite_strain[isubcase]
            eids = np.unique(case.element_layer[:, 0])
            for eid in eids:
                assert eid in out[card_type], 'eid=%s eids=%s card_type=%s'  % (eid, out[card_type], card_type)
        if op2.ctria3_composite_stress:
            case = op2.ctria3_composite_stress[isubcase]
            eids = np.unique(case.element_layer[:, 0])
            for eid in eids:
                assert eid in out[card_type], 'eid=%s eids=%s card_type=%s'  % (eid, out[card_type], card_type)
        if op2.ctria3_force:
            case = op2.ctria3_force[isubcase]
            eids = np.unique(case.element)
            for eid in eids:
                assert eid in out[card_type], 'eid=%s eids=%s card_type=%s'  % (eid, out[card_type], card_type)

    def test_op2_dmi_01(self):
        folder = os.path.abspath(os.path.join(test_path, '..', 'models'))
        bdf_filename = os.path.join(folder, 'matrix', 'matrix.dat')
        op2_filename = os.path.join(folder, 'matrix', 'mymatrix.op2')
        matrices = {
            'A' : True,
            'B' : False,
            'ATB' : False,
            'BTA' : False,
            'MYDOF' : True,
        }
        model = BDF(debug=False)
        model.read_bdf(bdf_filename)

        dmi_a = model.dmis['A']
        assert dmi_a.shape == (4, 2), 'shape=%s' % (dmi_a.shape)
        #print('dmi_a\n', dmi_a)
        a, rows_reversed, cols_reversed = dmi_a.get_matrix(is_sparse=False, apply_symmetry=False)
        #print('model.dmi.A =\n%s' % dmi_a)
        #print('model.dmi.A =\n%s' % str(a))
        #return
        op2 = OP2(debug=False)
        op2.set_additional_matrices_to_read(matrices)
        try:
            op2.read_op2(op2_filename)
            raise RuntimeError('this is wrong...')
        except FatalError:
            # the OP2 doesn't have a trailing zero marker
            pass

        # M rows, Ncols
        A = np.array([
            [1., 0.],
            [3., 6.],
            [5., 0.],
            [0., 8.],
        ], dtype='float32')
        B = A
        mydof = np.array([
            -1.0, 1.0, 1.0, -1.0, 1.0,
            2.0, -1.0, 1.0, 3.0, -1.0, 1.0, 4.0, -1.0,
            1.0, 5.0, -1.0, 1.0, 6.0, -1.0, 2.0, 1.0,
            -1.0, 2.0, 2.0, -1.0, 2.0, 3.0, -1.0, 2.0,
            4.0, -1.0, 2.0, 5.0, -1.0, 2.0, 6.0,
        ])
        BTA = np.dot(B.T, A)
        ATB = np.dot(A.T, B)
        ATB_expected = np.array([
            [35., 18.],
            [18., 100.]
        ], dtype='float32')
        BTA_expected = ATB_expected

        expecteds = [A, ATB, B, BTA, mydof]
        matrix_names = sorted(matrices.keys())

        for table_name, expected in zip(matrix_names, expecteds):
            assert table_name in op2.matrices, table_name


            actual = op2.matrices[table_name].data
            if not (np.array_equal(expected, actual) or
                    np.array_equal(expected, np.squeeze(actual))):
                if table_name in model.dmis:
                    dmi = model.dmis[table_name]
                    table_array, rows_reversed, cols_reversed = dmi.get_matrix(is_sparse=False, apply_symmetry=False)
                    #stable_array, rows_reversed, cols_reversed = dmi.get_matrix(is_sparse=True, apply_symmetry=False)
                    print(table_array)
                #print(stable_array)
                msg = 'matrix %s was not read properly\n' % table_name
                msg += 'expected shape=%s\n%s\n' % (str(expected.shape), expected)
                msg += 'actual shape=%s\n%s' % (str(actual.shape), actual.ravel())
                #msg += '\n%s' % actual.ravel()
                print(msg)
                print('==========================')
                #raise RuntimeError(msg)

    def test_op2_dmi_02(self):
        folder = os.path.abspath(os.path.join(test_path, '..', 'models'))
        bdf_filename = os.path.join(folder, 'matrix', 'matrix.dat')
        op2_filename = os.path.join(folder, 'matrix', 'mymatrix.op2')
        matrices = {
            'A' : True,
            'B' : False,
            'ATB' : False,
            'BTA' : False,
            'MYDOF' : True,
        }
        model = BDF(debug=False)
        model.read_bdf(bdf_filename)

        dmi_a = model.dmis['A']
        a, rows_reversed, cols_reversed = dmi_a.get_matrix(is_sparse=False, apply_symmetry=False)
        #print('model.dmi.A =\n%s' % dmi_a)
        #print('model.dmi.A =\n%s' % str(a))
        #return
        op2 = OP2()
        try:
            op2.read_op2(op2_filename, skip_undefined_matrices=True)
            raise RuntimeError('this is wrong...')
        except FatalError:
            # the OP2 doesn't have a trailing zero marker
            pass

        # M rows, Ncols
        A = np.array([
            [1., 0.],
            [3., 6.],
            [5., 0.],
            [0., 8.],
        ], dtype='float32')
        B = A
        mydof = np.array([
            -1.0, 1.0, 1.0, -1.0, 1.0,
            2.0, -1.0, 1.0, 3.0, -1.0, 1.0, 4.0, -1.0,
            1.0, 5.0, -1.0, 1.0, 6.0, -1.0, 2.0, 1.0,
            -1.0, 2.0, 2.0, -1.0, 2.0, 3.0, -1.0, 2.0,
            4.0, -1.0, 2.0, 5.0, -1.0, 2.0, 6.0,
        ])
        BTA = dot(B.T, A)
        ATB = dot(A.T, B)
        ATB_expected = np.array([
            [35., 18.],
            [18., 100.]
        ], dtype='float32')
        BTA_expected = ATB_expected

        expecteds = [A, ATB, B, BTA, mydof]
        matrix_names = sorted(matrices.keys())

        for table_name, expected in zip(matrix_names, expecteds):
            assert table_name in op2.matrices, table_name

            actual = op2.matrices[table_name].data
            if not array_equal(expected, actual):
                if table_name in model.dmis:
                    dmi = model.dmis[table_name]
                    table_array, rows_reversed, cols_reversed = dmi.get_matrix(is_sparse=False, apply_symmetry=False)
                    #stable_array, rows_reversed, cols_reversed = dmi.get_matrix(is_sparse=True, apply_symmetry=False)
                    #print(table_array)
                #print(stable_array)
                msg = 'matrix %s was not read properly\n' % table_name
                msg += 'expected\n%s\n' % expected
                msg += 'actual\n%s' % actual
                print(msg)
                print('==========================')
                #raise RuntimeError(msg)

def _get_gpforce_data():
    data = [
        #eids, nids, cid, summation_point
        #[1], [1], 0, [0., 0., 0.],
        [[1], [1, 2, 3, 4], 0, [0., 0., 0.], [0.0, 0.0, -10000.0], [-5000.0, 5000.0, 0.0],],  # total

        # cid=0; eid=[1]; nid=[3]; sum=[0., 0., 0.] - done
        #               fmag     mmag       fx      fy       fz       mx       my       mz
        # F2      = [2589.95,     0.0,  26.34, -44.15, -2589.44,     0.0,     0.0,     0.0]  # ith
        # F2Total = [2589.95, 3862.70,  26.34, -44.15, -2589.44, -2589.44, 2589.44, -70.49]  # total
        [[1], [3], 0, [0., 0., 0.], [26.34, -44.15, -2589.44], [-2589.44, 2589.44, -70.49],],

        # cid=0; eid=[1]; nid=[1]; sum=[0., 0., 0.]
        #                            fx      fy       fz       mx       my       mz
        [[1], [1], 0, [0., 0., 0.], [-37.18, 32.00, -2589.44], [0.0, 0.0, 0.0],],  # only 1 line b/c no moment

        # cid=1; eid=[1]; nid=[1]; sum=[0., 0., 0.]
        #               fmag     mmag      fx      fy       fz       mx       my       mz
        # F1      = [2589.90,     0.0,1853.64, 567.74, 1717.35,     0.0,     0.0,     0.0]  # ith
        # F1Total = [2589.90,     0.0,1853.64, 567.74, 1717.35,     0.0,     0.0,     0.0]  # total
        [[1], [1], 1, [0., 0., 0.], [1853.64, 567.74, 1717.35], [0.0, 0.0, 0.0],],

        # cid=1; eid=[1]; nid=[2]; sum=[0., 0., 0.]
        #               fmag     mmag       fx      fy       fz       mx       my       mz
        # F2      = [2411.67,     0.0, 1710.67, 634.80, 1577.03,     0.0,     0.0,     0.0]  # ith
        # F2Total = [2411.67, 2410.58, 1710.67, 634.80, 1577.03, 1698.38, -570.22, -1612.84] # total
        [[1], [2], 1, [0., 0., 0.], [1710.67, 634.80, 1577.03], [1698.38, -570.22, -1612.84],],

        # cid=1; eid=[1]; nid=[3]; sum=[0., 0., 0.]
        #           fmag          mmag     fx       fy       fz       mx        my       mz
        # F3      = [2589.95,     0.0, 1799.79, 645.58, 1746.94,     0.0,      0.0,     0.0]  # ith
        # F3Total = [2589.95, 3862.70, 1799.79, 645.58, 1746.94, 1880.85, -3035.07, -816.15]  # total
        [[1], [3], 1, [0., 0., 0.], [1799.79, 645.58, 1746.94], [1880.85, -3035.07, -816.15]],
    ]
    return data


if __name__ == '__main__':  # pragma: no cover
    on_rtd = os.environ.get('READTHEDOCS', None) == 'True'
    if not on_rtd:
        unittest.main()
