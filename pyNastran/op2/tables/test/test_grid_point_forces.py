import os
import unittest
from io import StringIO

import numpy as np
from cpylog import SimpleLogger

import pyNastran

from pyNastran.op2.op2 import OP2
from pyNastran.op2.op2_geom import read_op2_geom

from pyNastran.op2.tables.ogf_gridPointForces.ogf_objects import RealGridPointForcesArray

test_path = pyNastran.__path__[0]
model_path = os.path.abspath(os.path.join(test_path, '..', 'models'))


class TestGridPointForces(unittest.TestCase):
    """various grid point force tests"""

    def test_gpforce_01(self):
        nids = np.array([1, 2, 3])
        xyz_cid0 = np.array([
            [1., 1., 1.],
            [4., 2., 5.],
            [3., 3., 3.],
        ])
        data_code = {
            'nonlinear_factor' : None,
            'sort_bits' : [0, 0, 0],
            'analysis_code' : 1,
            'is_msc' : True,
            'format_code' : 1,
            'table_code' : 1,
            'data_names' : 'cat',
            'device_code' : 1,
            'size' : 4,
            #'tcode' : 1,
        }
        is_sort1 = True
        isubcase = 1
        dt = 0.0
        gpforce = RealGridPointForcesArray(data_code, is_sort1, isubcase, dt)
        gpforce.ntimes = 1
        gpforce.ntotal = 3
        gpforce._ntotals = [3]

        gpforce.build()
        gpforce.data[0, :, :] = np.array([
            [3., 7., 11., 0., 0., 0.,], # fx, fy, fz, mx, my, mz
            [3., 7., 11., 0., 0., 0.,],
            [3., 7., 11., 0., 0., 0.,],
        ])
        gpforce.node_element[0, :, :] = np.array([
            [1, 1],
            [2, 1],
            [3, 1],
        ])
        op2 = OP2()
        summation_point = [0., 0., 0.]
        icd_transform = None
        nid_cd = np.array([
            [1, 0],
            [2, 0],
            [3, 0],
        ])
        from pyNastran.bdf.bdf import CORD2R
        coord_out = CORD2R(cid=0, origin=None, zaxis=None, xzplane=None)
        coords = {0 : coord_out}

        #eids = [1]
        #nids = [1]
        #gpforce.extract_interface_loads(
            #nids, eids, coord_out, coords, nid_cd,
            #icd_transform,
            #xyz_cid0,
            #summation_point,
            #itime=0,
            #debug=True,
            #log=op2.log)

        #print('------------')
        #eids = [1]
        #nids = [2]
        #gpforce.extract_interface_loads(
            #nids, eids, coord_out, coords, nid_cd,
            #icd_transform,
            #xyz_cid0,
            #summation_point,
            #itime=0,
            #debug=True,
            #log=op2.log)
        print('------------')

        eids = [1]
        nids = [1, 2]
        gpforce.extract_interface_loads(
            nids, eids, coord_out, coords, nid_cd,
            icd_transform,
            xyz_cid0,
            summation_point,
            itime=0,
            debug=True,
            log=op2.log)
        #print(gpforce)

    def test_op2_solid_shell_bar_01_gpforce(self):
        folder = os.path.join(model_path, 'sol_101_elements')
        #bdf_filename = os.path.join(folder, 'static_solid_shell_bar.bdf')
        op2_filename = os.path.join(folder, 'static_solid_shell_bar.op2')
        op2 = read_op2_geom(op2_filename, xref=False, debug=False)

        nids_all, nids_transform, icd_transform = op2.get_displacement_index()
        op2.transform_displacements_to_global(icd_transform, op2.coords)

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
        nid_cd = np.array([[nid, node.Cd()] for nid, node in sorted(op2.nodes.items())])
        coords = op2.coords

        data = _get_gpforce_data()
        for datai in data:
            eids, nids, cid, summation_point, total_force_local_expected, total_moment_local_expected = datai
            if cid not in coords:
                continue
            #op2.log.debug('*' * 30 + 'Next Test' + '*' * 30)
            coord_out = coords[cid]
            out = gpforce.extract_interface_loads(
                nids, eids,
                coord_out, coords,
                nid_cd, icd_transform,
                xyz_cid0, summation_point, itime=0, debug=False, log=op2.log)
            total_force_global, total_moment_global, total_force_local, total_moment_local = out

            #op2.log.debug('***********')
            #op2.log.debug('force = %s; %s' % (total_force_global, np.linalg.norm(total_force_global)))
            #op2.log.debug('moment = %s; %s' % (total_moment_global, np.linalg.norm(total_moment_global)))

            case = 'eids=%s nids=%s cid=%s summation_point=%s' % (
                eids, nids, cid, summation_point)
            msg = '%s\ntotal_force_local_expected=%s total_force_local=%s' % (
                case, total_force_local_expected, total_force_local)
            self.assertTrue(np.allclose(total_force_local_expected, total_force_local, atol=0.005), msg)

            msg = '%s\ntotal_moment_local_expected=%s total_moment_local=%s' % (
                case, total_moment_local_expected, total_moment_local)
            self.assertTrue(np.allclose(total_moment_local_expected, total_moment_local, atol=0.005), msg)

    def test_op2_solid_shell_bar_01_gpforce_xyz(self):
        folder = os.path.join(model_path, 'sol_101_elements')
        #bdf_filename1 = os.path.join(folder, 'static_solid_shell_bar_xyz.bdf')
        op2_filename1 = os.path.join(folder, 'static_solid_shell_bar_xyz.op2')
        op2_1 = read_op2_geom(op2_filename1, xref=False, debug=False)
        #print("disp_orig =\n", op2_1.displacements[1].data[0, :2, :])
        #print("spc_orig =\n", op2_1.spc_forces[1].data[0, -3:, :])
        print("gpf_orig =\n", op2_1.grid_point_forces[1].data[0, :2, :])

        nids_all, nids_transform_1, icd_transform_1 = op2_1.get_displacement_index()
        op2_1.transform_displacements_to_global(icd_transform_1, op2_1.coords)
        op2_1.transform_gpforce_to_global(nids_all, nids_transform_1, icd_transform_1, op2_1.coords)
        #print("disp_new =\n", op2_1.displacements[1].data[0, :2, :])
        #print("spc_new =\n", op2_1.spc_forces[1].data[0, -3:, :])
        print("gpf_new =\n", op2_1.grid_point_forces[1].data[0, :2, :])

        gpforce = op2_1.grid_point_forces[1]
        op2_1.cross_reference(xref_elements=False,
                              xref_nodes_with_elements=False,
                              xref_properties=False,
                              xref_masses=False,
                              xref_materials=False,
                              xref_loads=False,
                              xref_constraints=False,
                              xref_aero=False,
                              xref_sets=False,
                              xref_optimization=False)
        xyz_cid0 = op2_1.get_xyz_in_coord(cid=0)
        nid_cd = np.array([[nid, node.Cd()] for nid, node in sorted(op2_1.nodes.items())])

        #bdf_filename2 = os.path.join(folder, 'static_solid_shell_bar.bdf')
        op2_filename2 = os.path.join(folder, 'static_solid_shell_bar.op2')
        op2_2 = read_op2_geom(op2_filename2, debug=False)
        nids_all, nids_transform_2, icd_transform_2 = op2_2.get_displacement_index()
        op2_2.transform_displacements_to_global(icd_transform_2, op2_2.coords)
        op2_2.transform_gpforce_to_global(nids_all, nids_transform_2, icd_transform_2, op2_2.coords)

        #print("disp_goal =\n", op2_2.displacements[1].data[0, :2, :])
        #print("spc_goal =\n", op2_2.spc_forces[1].data[0, -3:, :])
        print("gpf_goal =\n", op2_2.grid_point_forces[1].data[0, :2, :])

        msg = 'displacements baseline=\n%s\ndisplacements xyz=\n%s' % (
            op2_1.displacements[1].data[0, :, :], op2_2.displacements[1].data[0, :, :])
        #print(msg)
        assert op2_1.displacements[1].assert_equal(op2_2.displacements[1])

        msg = 'grid_point_forces baseline=\n%s\ngrid_point_forces xyz=\n%s' % (
            op2_1.grid_point_forces[1].data[0, :, :], op2_2.grid_point_forces[1].data[0, :, :])
        #print(msg)

        assert op2_1.spc_forces[1].assert_equal(op2_2.spc_forces[1], atol=4.4341e-04), msg
        assert op2_1.mpc_forces[1].assert_equal(op2_2.mpc_forces[1]), msg
        assert op2_1.load_vectors[1].assert_equal(op2_2.load_vectors[1]), msg
        assert op2_1.grid_point_forces[1].assert_equal(op2_2.grid_point_forces[1], atol=0.000123), msg

        #-------------------------------------------------
        return
        data = _get_gpforce_data()
        coords = op2_1.coords
        #used_cds = np.unique(nid_cd[:, 1])
        #for cd in used_cds:
            #coord = op2_1.coords[cd]
            #print(coord)
            #print('origin = %s' % coord.origin)
            #print('beta =\n%s' % coord.beta())
            #print('-----------------------------')

        for datai in data:
            eids, nids, cid, summation_point, total_force_local_expected, total_moment_local_expected = datai
            coord_out = coords[cid]
            op2_1.log.debug('*' * 30 + 'Next Test' + '*' * 30)
            out = gpforce.extract_interface_loads(
                nids, eids,
                coord_out, coords,
                nid_cd, icd_transform_1,
                xyz_cid0, summation_point, itime=0, debug=False, log=op2_1.log)
            total_force_global, total_moment_global, total_force_local, total_moment_local = out

            op2_1.log.debug('***********')
            op2_1.log.debug('force = %s; %s' % (total_force_global, np.linalg.norm(total_force_global)))
            op2_1.log.debug('moment = %s; %s' % (total_moment_global, np.linalg.norm(total_moment_global)))

            case = 'eids=%s nids=%s cid=%s summation_point=%s' % (
                eids, nids, cid, summation_point)
            msg = '%s\ntotal_force_local_expected=%s total_force_local=%s delta=%s' % (
                case, total_force_local_expected, total_force_local,
                np.abs(total_force_local_expected - total_force_local))
            self.assertTrue(np.allclose(total_force_local_expected, total_force_local, atol=0.2), msg)

            msg = '%s\ntotal_moment_local_expected=%s total_moment_local=%s delta=%s' % (
                case, total_moment_local_expected, total_moment_local,
                np.abs(total_moment_local_expected - total_moment_local))
            self.assertTrue(np.allclose(total_moment_local_expected, total_moment_local, atol=0.005), msg)

    @unittest.expectedFailure
    def test_broken_op2_solid_shell_bar_01_gpforce_radial_global_cd(self):
        warning_log = SimpleLogger(level='warning')
        debug_log = SimpleLogger(level='debug')
        folder = os.path.join(model_path, 'sol_101_elements')
        op2_filename1 = os.path.join(folder, 'static_solid_shell_bar_global_radial_cd.op2')
        op2_1 = read_op2_geom(op2_filename1, xref=False, log=warning_log)
        op2_1.log = debug_log

        #print("disp_orig =\n", op2_1.displacements[1].data[0, :4, :])
        #print("spc_orig =\n", op2_1.spc_forces[1].data[0, -3:, :])
        print("gpf_orig =\n", op2_1.grid_point_forces[1].data[0, 2:8, :])

        op2_1.cross_reference(xref_elements=False,
                              xref_nodes_with_elements=False,
                              xref_properties=False,
                              xref_masses=False,
                              xref_materials=False,
                              xref_loads=False,
                              xref_constraints=False,
                              xref_aero=False,
                              xref_sets=False,
                              xref_optimization=False)
        xyz_cid0 = op2_1.get_xyz_in_coord(cid=0)

        nid_cd = np.array([[nid, node.Cd()] for nid, node in sorted(op2_1.nodes.items())])
        #-------------------------------------------------
        #coords = op2_1.coords
        #used_cds = np.unique(nid_cd[:, 1])
        #for cd in used_cds:
            #coord = op2_1.coords[cd]
            #print(coord)
            #print('origin = %s' % coord.origin)
            #print('beta =\n%s' % coord.beta())
            #print('-----------------------------')

        #disp = op2_1.displacements[1]
        #for line in list(disp.data[0, :, :3]):
            #print('%10.4e %10.4e %10.4e' % tuple(line))

        nids_all, nids_transform_1, icd_transform_1 = op2_1.get_displacement_index()
        op2_1.transform_displacements_to_global(icd_transform_1, op2_1.coords, xyz_cid0=xyz_cid0)
        op2_1.transform_gpforce_to_global(nids_all, nids_transform_1, icd_transform_1, op2_1.coords, xyz_cid0=xyz_cid0)
        #print('stuff...')
        #disp = op2_1.displacements[1]
        #for line in list(disp.data[0, :, :3]):
            #print('%10.4e %10.4e %10.4e' % tuple(line))
        #print(disp.data[0, :, :3])

        #print("disp_new =\n", op2_1.displacements[1].data[0, :4, :])
        #print("spc_new =\n", op2_1.spc_forces[1].data[0, -3:, :])
        print("gpf_new =\n", op2_1.grid_point_forces[1].data[0, 2:8, :])

        #-----------------------------------------------------------------------
        op2_filename2 = os.path.join(folder, 'static_solid_shell_bar.op2')
        op2_2 = read_op2_geom(op2_filename2, debug=False)
        nids_all, nids_transform_2, icd_transform_2 = op2_2.get_displacement_index()
        op2_2.transform_displacements_to_global(icd_transform_2, op2_2.coords)
        op2_2.transform_gpforce_to_global(nids_all, nids_transform_2, icd_transform_2, op2_2.coords, xyz_cid0=xyz_cid0)

        #print("disp_goal =\n", op2_2.displacements[1].data[0, :4, :])
        #print("spc_goal =\n", op2_2.spc_forces[1].data[0, -3:, :])
        print("gpf_goal =\n", op2_2.grid_point_forces[1].data[0, 2:8, :])

        #return
        #msg = 'displacements baseline=\n%s\ndisplacements xyz=\n%s' % (
            #op2_1.displacements[1].data[0, :, :], op2_2.displacements[1].data[0, :, :])
        #print(msg)
        #assert op2_1.displacements[1].assert_equal(op2_2.displacements[1])

        #msg = 'grid_point_forces baseline=\n%s\ngrid_point_forces xyz=\n%s' % (
            #op2_1.grid_point_forces[1].data[0, :, :], op2_2.grid_point_forces[1].data[0, :, :])
        #print(msg)

        csv_file = StringIO()
        op2_1.spc_forces[1].write_csv(csv_file)
        op2_1.grid_point_forces[1].write_csv(csv_file)
        #print(csv_file.getvalue())

        assert op2_1.spc_forces[1].assert_equal(op2_2.spc_forces[1], atol=4.4341e-04)
        assert op2_1.mpc_forces[1].assert_equal(op2_2.mpc_forces[1])
        assert op2_1.load_vectors[1].assert_equal(op2_2.load_vectors[1])
        #print('op2_2.grid_point_forces[1]\n', op2_2.grid_point_forces[1].data)
        assert op2_1.grid_point_forces[1].assert_equal(op2_2.grid_point_forces[1], atol=0.000123)
        return
        #-----------------------------------------------------------------------
        gpforce = op2_1.grid_point_forces[1]
        data = _get_gpforce_data()
        for i, datai in enumerate(data):
            eids, nids, cid, summation_point, total_force_local_expected, total_moment_local_expected = datai
            coord_out = op2_1.coords[cid]
            op2_1.log.debug('*' * 30 + 'Next Test #%s' % i + '*' * 30)
            out = gpforce.extract_interface_loads(
                nids, eids,
                coord_out, op2_1.coords,
                nid_cd, icd_transform_1,
                xyz_cid0, summation_point, itime=0, debug=False, log=op2_1.log)
            total_force_global, total_moment_global, total_force_local, total_moment_local = out

            op2_1.log.debug('***********')
            op2_1.log.debug('force = %s; %s' % (total_force_global, np.linalg.norm(total_force_global)))
            op2_1.log.debug('moment = %s; %s' % (total_moment_global, np.linalg.norm(total_moment_global)))

            case = 'eids=%s nids=%s cid=%s summation_point=%s' % (
                eids, nids, cid, summation_point)
            msg = '%s\ntotal_force_local_expected=%s total_force_local=%s delta=%s' % (
                case, total_force_local_expected, total_force_local,
                np.abs(total_force_local_expected - total_force_local))
            self.assertTrue(np.allclose(total_force_local_expected, total_force_local, atol=0.2), msg)

            msg = '%s\ntotal_moment_local_expected=%s total_moment_local=%s delta=%s' % (
                case, total_moment_local_expected, total_moment_local,
                np.abs(total_moment_local_expected - total_moment_local))
            self.assertTrue(np.allclose(total_moment_local_expected, total_moment_local, atol=0.005), msg)

    #@unittest.expectedFailure
    def test_op2_solid_shell_bar_01_gpforce_radial(self):
        warning_log = SimpleLogger(level='warning')
        debug_log = SimpleLogger(level='debug')
        folder = os.path.join(model_path, 'sol_101_elements')
        op2_filename = os.path.join(folder, 'static_solid_shell_bar_radial.op2')
        op2_1 = read_op2_geom(op2_filename, xref=False, log=warning_log)
        op2_1.log = debug_log

        print("disp_orig =\n", op2_1.displacements[1].data[0, :2, :])
        #print("spc_orig =\n", op2_1.spc_forces[1].data[0, -3:, :])
        #print("gpf_orig =\n", op2_1.grid_point_forces[1].data[0, :2, :])

        op2_1.cross_reference(xref_elements=False,
                              xref_nodes_with_elements=False,
                              xref_properties=False,
                              xref_masses=False,
                              xref_materials=False,
                              xref_loads=False,
                              xref_constraints=False,
                              xref_aero=False,
                              xref_sets=False,
                              xref_optimization=False)
        xyz_cid0 = op2_1.get_xyz_in_coord(cid=0)

        nid_cd = np.array([[nid, node.Cd()] for nid, node in sorted(op2_1.nodes.items())])
        #-------------------------------------------------
        #coords = op2_1.coords
        #used_cds = np.unique(nid_cd[:, 1])
        #for cd in used_cds:
            #coord = op2_1.coords[cd]
            #print(coord)
            #print('origin = %s' % coord.origin)
            #print('beta =\n%s' % coord.beta())
            #print('-----------------------------')

        #disp = op2_1.displacements[1]
        #for line in list(disp.data[0, :, :3]):
            #print('%10.4e %10.4e %10.4e' % tuple(line))

        nids_all, nids_transform_1, icd_transform_1 = op2_1.get_displacement_index()
        op2_1.transform_displacements_to_global(icd_transform_1, op2_1.coords, xyz_cid0=xyz_cid0)
        op2_1.transform_gpforce_to_global(
            nids_all, nids_transform_1, icd_transform_1, op2_1.coords, xyz_cid0=xyz_cid0)
        #print('stuff...')
        #disp = op2_1.displacements[1]
        #for line in list(disp.data[0, :, :3]):
            #print('%10.4e %10.4e %10.4e' % tuple(line))
        #print(disp.data[0, :, :3])

        print("disp_new =\n", op2_1.displacements[1].data[0, :2, :])
        #print("spc_new =\n", op2_1.spc_forces[1].data[0, -3:, :])
        #print("gpf_new =\n", op2_1.grid_point_forces[1].data[0, :2, :])

        #-----------------------------------------------------------------------
        op2_filename2 = os.path.join(folder, 'static_solid_shell_bar.op2')
        op2_2 = read_op2_geom(op2_filename2, debug=False)
        nids_all, nids_transform_2, icd_transform_2 = op2_2.get_displacement_index()
        op2_2.transform_displacements_to_global(icd_transform_2, op2_2.coords)
        op2_2.transform_gpforce_to_global(
            nids_all, nids_transform_2, icd_transform_2, op2_2.coords)

        print("disp_goal =\n", op2_2.displacements[1].data[0, :2, :])
        #print("spc_goal =\n", op2_2.spc_forces[1].data[0, -3:, :])
        #print("gpf_goal =\n", op2_2.grid_point_forces[1].data[0, :2, :])

        return
        msg = 'displacements baseline=\n%s\ndisplacements xyz=\n%s' % (
            op2_1.displacements[1].data[0, :, :], op2_2.displacements[1].data[0, :, :])
        #print(msg)
        assert op2_1.displacements[1].assert_equal(op2_2.displacements[1])

        msg = 'grid_point_forces baseline=\n%s\ngrid_point_forces xyz=\n%s' % (
            op2_1.grid_point_forces[1].data[0, :, :], op2_2.grid_point_forces[1].data[0, :, :])
        #print(msg)

        assert op2_1.spc_forces[1].assert_equal(op2_2.spc_forces[1], atol=4.4341e-04), msg
        assert op2_1.mpc_forces[1].assert_equal(op2_2.mpc_forces[1]), msg
        assert op2_1.load_vectors[1].assert_equal(op2_2.load_vectors[1]), msg
        assert op2_1.grid_point_forces[1].assert_equal(op2_2.grid_point_forces[1], atol=0.000123), msg
        #-----------------------------------------------------------------------
        gpforce = op2_1.grid_point_forces[1]
        data = _get_gpforce_data()
        for i, datai in enumerate(data):
            eids, nids, cid, summation_point, total_force_local_expected, total_moment_local_expected = datai
            coord_out = op2_1.coords[cid]
            op2_1.log.debug('*' * 30 + 'Next Test #%s' % i + '*' * 30)
            out = gpforce.extract_interface_loads(
                nids, eids,
                coord_out, op2_1.coords,
                nid_cd, icd_transform_1,
                xyz_cid0, summation_point, itime=0, debug=False, log=op2_1.log)
            total_force_global, total_moment_global, total_force_local, total_moment_local = out

            op2_1.log.debug('***********')
            op2_1.log.debug('force = %s; %s' % (total_force_global, np.linalg.norm(total_force_global)))
            op2_1.log.debug('moment = %s; %s' % (total_moment_global, np.linalg.norm(total_moment_global)))

            case = 'eids=%s nids=%s cid=%s summation_point=%s' % (
                eids, nids, cid, summation_point)
            msg = '%s\ntotal_force_local_expected=%s total_force_local=%s delta=%s' % (
                case, total_force_local_expected, total_force_local,
                np.abs(total_force_local_expected - total_force_local))
            self.assertTrue(np.allclose(total_force_local_expected, total_force_local, atol=0.2), msg)

            msg = '%s\ntotal_moment_local_expected=%s total_moment_local=%s delta=%s' % (
                case, total_moment_local_expected, total_moment_local,
                np.abs(total_moment_local_expected - total_moment_local))
            self.assertTrue(np.allclose(total_moment_local_expected, total_moment_local, atol=0.005), msg)

def _get_gpforce_data():
    data = [
        #eids, nids, cid, summation_point
        #[1], [1], 0, [0., 0., 0.],
        #[[1], [1, 2, 3, 4], 0, [0., 0., 0.], [0.0, 0.0, -10000.0], [-5000.0, 5000.0, 0.0],],  # total; good for gpforce

        # cid=0; eid=[1]; nid=[3]; sum=[0., 0., 0.] - done
        #               fmag     mmag       fx      fy       fz       mx       my       mz
        # F2      = [2589.95,     0.0,  26.34, -44.15, -2589.44,     0.0,     0.0,     0.0]  # ith
        # F2Total = [2589.95, 3862.70,  26.34, -44.15, -2589.44, -2589.44, 2589.44, -70.49]  # total
        #[[1], [3], 0, [0., 0., 0.], [26.34, -44.15, -2589.44], [-2589.44, 2589.44, -70.49],], # good for gpforce; failing for xyz (cid 11)

        # cid=0; eid=[1]; nid=[1]; sum=[0., 0., 0.]
        #                            fx      fy       fz       mx       my       mz
        [[1], [1], 0, [0., 0., 0.], [-37.18, 32.00, -2589.44], [0.0, 0.0, 0.0],],  # only 1 line b/c no moment; good for gpforce; failing for xyz (cid 11)

        # cid=0/1/2/3; eid=[1]; nid=[1]; sum=[0., 0., 0.]
        [[1], [1], 0, [0., 0., 0.], [-37.18, 32.00, -2589.44], [0.0, 0.0, 0.0],],  # only 1 line b/c no moment
        [[1], [1], 1, [0., 0., 0.], [-37.18, 32.00, -2589.44], [0.0, 0.0, 0.0],],  # only 1 line b/c no moment
        [[1], [1], 2, [0., 0., 0.], [-37.18, 32.00, -2589.44], [0.0, 0.0, 0.0],],  # only 1 line b/c no moment
        [[1], [1], 3, [0., 0., 0.], [-37.18, 32.00, -2589.44], [0.0, 0.0, 0.0],],  # only 1 line b/c no moment

        # cid=1; eid=[1]; nid=[1]; sum=[0., 0., 0.]
        #               fmag     mmag      fx      fy       fz       mx       my       mz
        # F1      = [2589.90,     0.0,1853.64, 567.74, 1717.35,     0.0,     0.0,     0.0]  # ith
        # F1Total = [2589.90,     0.0,1853.64, 567.74, 1717.35,     0.0,     0.0,     0.0]  # total
        [[1], [1], 11, [0., 0., 0.], [1853.64, 567.74, 1717.35], [0.0, 0.0, 0.0],], # good; failing for gpforce
        [[1], [1], 12, [0., 0., 0.], [1853.64, 567.74, 1717.35], [0.0, 0.0, 0.0],], # good; failing for gpforce
        [[1], [1], 13, [0., 0., 0.], [1853.64, 567.74, 1717.35], [0.0, 0.0, 0.0],], # good; failing for gpforce

        # cid=1; eid=[1]; nid=[2]; sum=[0., 0., 0.]
        #               fmag     mmag       fx      fy       fz       mx       my       mz
        # F2      = [2411.67,     0.0, 1710.67, 634.80, 1577.03,     0.0,     0.0,     0.0]  # ith
        # F2Total = [2411.67, 2410.58, 1710.67, 634.80, 1577.03, 1698.38, -570.22, -1612.84] # total
        [[1], [2], 11, [0., 0., 0.], [1710.67, 634.60, 1577.03], [1698.38, -570.22, -1612.84],],

        # cid=1; eid=[1]; nid=[3]; sum=[0., 0., 0.]
        #           fmag          mmag     fx       fy       fz       mx        my       mz
        # F3      = [2589.95,     0.0, 1799.79, 645.58, 1746.94,     0.0,      0.0,     0.0]  # ith
        # F3Total = [2589.95, 3862.70, 1799.79, 645.58, 1746.94, 1880.85, -3035.07, -816.15]  # total
        [[1], [3], 11, [0., 0., 0.], [1799.79, 645.58, 1746.94], [1880.85, -3035.07, -816.15]],

        #[[1], [1], 11, [0., 0., 0.], [1853.64, 567.74, 1717.35], [0., 0., 0.],],
        #[[1], [1], 12, [0., 0., 0.], [1938.05, -47.57, 1717.35], [0., 0., 0.],],
        #[[1], [1], 13, [0., 0., 0.], [2069.00, 1557.11, -47.57], [0., 0., 0.],],
    ]
    return data


if __name__ == "__main__":  # pragma: no cover
    unittest.main()

