from __future__ import print_function
import os
import unittest
from six import iteritems
import numpy as np

from pyNastran.bdf.bdf import BDF, CORD2R
from pyNastran.op2.op2 import OP2, FatalError
from pyNastran.op2.op2_geom import read_op2_geom
from pyNastran.op2.tables.ogf_gridPointForces.ogf_Objects import RealGridPointForcesArray
from pyNastran.bdf.test.bdf_unit_tests import Tester
#from pyNastran.utils.log import SimpleLogger

import pyNastran
PKG_PATH = pyNastran.__path__[0]
MODEL_PATH = os.path.abspath(os.path.join(PKG_PATH, '..', 'models'))


class TestOP2(Tester):
    """various OP2 tests"""

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
        i_transform = None
        nid_cd = np.array([
            [1, 0],
            [2, 0],
            [3, 0],
        ])
        coord_out = CORD2R(cid=0)
        coords = {0 : coord_out}

        #eids = [1]
        #nids = [1]
        #gpforce.extract_interface_loads(
            #nids, eids, coord_out, coords, nid_cd,
            #i_transform,
            #xyz_cid0,
            #summation_point,
            #itime=0,
            #debug=True,
            #logger=op2.log)

        #print('------------')
        #eids = [1]
        #nids = [2]
        #gpforce.extract_interface_loads(
            #nids, eids, coord_out, coords, nid_cd,
            #i_transform,
            #xyz_cid0,
            #summation_point,
            #itime=0,
            #debug=True,
            #logger=op2.log)
        #print('------------')

        eids = [1]
        nids = [1, 2]
        gpforce.extract_interface_loads(
            nids, eids, coord_out, coords, nid_cd,
            i_transform,
            xyz_cid0,
            summation_point,
            itime=0,
            debug=True,
            logger=op2.log)
        #print(gpforce)

    def test_op2_solid_shell_bar_01_gpforce(self):
        #bdf_filename = os.path.join(MODEL_PATH, 'sol_101_elements', 'static_solid_shell_bar.bdf')
        op2_filename = os.path.join(MODEL_PATH, 'sol_101_elements', 'static_solid_shell_bar.op2')
        op2 = read_op2_geom(op2_filename, debug=False)

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
        nid_cd = np.array([[nid, node.Cd()] for nid, node in sorted(iteritems(op2.nodes))])
        coords = op2.coords

        data = [
            #eids, nids, cid, summation_point
            #[1], [1], 0, [0., 0., 0.],
            [[1], [1, 2, 3, 4], 0, [0., 0., 0.], [0.0, 0.0, -10000.0], [-5000.0, 5000.0, 0.0],],

            # cid=0; eid=[1]; nid=[3]; sum=[0., 0., 0.] - done
            #               fmag     mmag       fx      fy       fz       mx       my       mz
            # F2      = [2589.95,     0.0,  26.34, -44.15, -2589.44,     0.0,     0.0,     0.0]
            # F2Total = [2589.95, 3862.70,  26.34, -44.15, -2589.44, -2589.44, 2589.44, -70.49]
            [[1], [3], 0, [0., 0., 0.], [26.34, -44.15, -2589.44], [-2589.44, 2589.44, -70.49],],

            # cid=1; eid=[1]; nid=[1]; sum=[0., 0., 0.]
            #               fmag     mmag      fx      fy       fz       mx       my       mz
            # F1      = [2589.90,     0.0,1853.64, 567.74, 1717.35,     0.0,     0.0,     0.0]
            # F1Total = [2589.90,     0.0,1853.64, 567.74, 1717.35,     0.0,     0.0,     0.0]
            [[1], [1], 1, [0., 0., 0.], [1853.64, 567.74, 1717.35], [0.0, 0.0, 0.0],],

            # cid=1; eid=[1]; nid=[2]; sum=[0., 0., 0.]
            #               fmag     mmag       fx      fy       fz       mx       my       mz
            # F2      = [2411.67,     0.0, 1710.67, 634.80, 1577.03,     0.0,     0.0,     0.0]
            # F2Total = [2411.67, 2410.58, 1710.67, 634.80, 1577.03, 1698.38, -570.22, -1612.84]
            [[1], [2], 1, [0., 0., 0.], [1710.67, 634.80, 1577.03], [1698.38, -570.22, -1612.84],],

            # cid=1; eid=[1]; nid=[3]; sum=[0., 0., 0.]
            #           fmag          mmag     fx       fy       fz       mx        my       mz
            # F3      = [2589.95,     0.0, 1799.79, 645.58, 1746.94,     0.0,      0.0,     0.0]
            # F3Total = [2589.95, 3862.70, 1799.79, 645.58, 1746.94, 1880.85, -3035.07, -816.15]
            [[1], [3], 1, [0., 0., 0.], [1799.79, 645.58, 1746.94], [1880.85, -3035.07, -816.15]],
        ]
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
                xyz_cid0, summation_point, itime=0, debug=False, logger=op2.log)
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
        folder = os.path.join(MODEL_PATH, 'sol_101_elements')
        bdf_filename1 = os.path.join(folder, 'static_solid_shell_bar_xyz.bdf')
        op2_filename1 = os.path.join(folder, 'static_solid_shell_bar_xyz.op2')
        op2_1 = read_op2_geom(op2_filename1, debug=False)

        nids_all, nids_transform, icd_transform_1 = op2_1.get_displacement_index()
        op2_1.transform_displacements_to_global(icd_transform_1, op2_1.coords)

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
        nid_cd = np.array([[nid, node.Cd()] for nid, node in sorted(iteritems(op2_1.nodes))])

        bdf_filename2 = os.path.join(folder, 'static_solid_shell_bar.bdf')
        model = BDF(debug=False)
        model.read_bdf(bdf_filename)

        data = _get_gpforce_data()
        coords = op2_1.coords
        for datai in data:
            eids, nids, cid, summation_point, total_force_local_expected, total_moment_local_expected = datai
            if cid not in coords:
                continue
            coord_out = coords[cid]
            op2_1.log.debug('*' * 30 + 'Next Test' + '*' * 30)
            out = gpforce.extract_interface_loads(
                nids, eids,
                coord_out, coords,
                nid_cd, icd_transform_1,
                xyz_cid0, summation_point, itime=0, debug=False, logger=op2_1.log)
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
        folder = os.path.abspath(os.path.join(PKG_PATH, '..', 'models', 'sol_101_elements'))
        bdf_filename = os.path.join(folder, 'static_solid_shell_bar_radial.bdf')
        op2_filename = os.path.join(folder, 'static_solid_shell_bar_radial.op2')
        op2 = read_op2_geom(op2_filename)

        op2.transform_displacements_to_global(icd_transform, op2.coords)


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
        ]

        gpforce = op2.grid_point_forces[1]
        coords = op2.coords
        for datai in data:
            eids, nids, cid, summation_point, total_force_local_expected, total_moment_local_expected = datai
            if cid not in coords:
                continue
            coord_out = coords[cid]
            op2.log.debug('*' * 30 + 'Next Test #%s' % i + '*' * 30)
            out = gpforce.extract_interface_loads(
                nids, eids,
                coord_out, coords,
                nid_cd, icd_transform_1,
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
            self.assertTrue(np.allclose(total_force_local_expected, total_force_local, atol=0.2), msg)

            msg = '%s\ntotal_moment_local_expected=%s total_moment_local=%s delta=%s' % (
                case, total_moment_local_expected, total_moment_local,
                np.abs(total_moment_local_expected - total_moment_local))
            self.assertTrue(np.allclose(total_moment_local_expected, total_moment_local, atol=0.005), msg)

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

        #[[3], [1], 0, [0., 0., 0.],],

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
    ON_RTD = os.environ.get('READTHEDOCS', None) == 'True'
    if not ON_RTD:
        unittest.main()
