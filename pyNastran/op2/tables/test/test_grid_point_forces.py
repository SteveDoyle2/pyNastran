import os
import getpass
from pathlib import Path
from io import StringIO
import unittest

import numpy as np
from cpylog import SimpleLogger, WarningRedirector

import pyNastran
from pyNastran.bdf.bdf import BDF, CORD2R, read_bdf
from pyNastran.op2.op2 import OP2, read_op2
from pyNastran.op2.op2_geom import OP2Geom, read_op2_geom

from pyNastran.op2.tables.ogf_gridPointForces.smt import smt_setup, plot_smt, create_shear_moment_torque
from pyNastran.op2.tables.ogf_gridPointForces.ogf_objects import RealGridPointForcesArray

from pyNastran.bdf.mesh_utils.cut.utils import p1_p2_zaxis_to_cord2r, get_stations
from pyNastran.bdf.mesh_utils.cut.cut_model_by_plane import (
    get_element_centroids)

PKG_PATH = Path(pyNastran.__path__[0])
MODEL_PATH = (PKG_PATH / '..' / 'models').resolve()

try:
    import matplotlib.pyplot as plt
    IS_MATPLOTLIB = True
except ModuleNotFoundError:
    IS_MATPLOTLIB = False

class TestGridPointForcesSMT(unittest.TestCase):
    def test_p1_p2_zaxis_to_cord2r_zaxis(self):
        p1 = np.array([0., 0., 0.])
        p2 = np.array([1., 0., 0.])
        zaxis = np.array([0., 0., 1.])
        model = BDF(debug=False)
        model.cross_reference()
        xyz1, xyz2, z_global, i, k, origin, zaxis2, xzplane = p1_p2_zaxis_to_cord2r(
            model, p1, p2, zaxis, method='Z-Axis Projection',
            cid_p1=0, cid_p2=0, cid_zaxis=0)
        assert np.array_equal(xyz1, p1), xyz1
        assert np.array_equal(xyz2, p2), xyz2
        assert np.array_equal(z_global, zaxis), z_global

        assert np.array_equal(i, [1., 0., 0.]), i
        assert np.array_equal(k, [0., 0., 1.]), k

        assert np.array_equal(origin, [0., 0., 0.]), origin
        assert np.array_equal(zaxis2, p1+zaxis), zaxis2
        assert np.array_equal(xzplane, [1., 0., 0.]), xzplane

    def test_p1_p2_zaxis_to_cord2r_cord2r(self):
        p1 = np.array([1., 1., 1.])
        p2 = np.array([2., 1., 1.5])  # xzplane
        zaxis = np.array([1., 1., 2.])
        model = BDF(debug=False)
        model.cross_reference()
        xyz1, xyz2, z_global, i, k, origin, zaxis2, xzplane = p1_p2_zaxis_to_cord2r(
            model, p1, p2, zaxis, method='CORD2R',
            cid_p1=0, cid_p2=0, cid_zaxis=0)
        assert np.array_equal(xyz1, p1), xyz1
        assert np.array_equal(xyz2, p2), xyz2
        assert np.array_equal(z_global, zaxis), z_global

        assert np.array_equal(i, [1., 0., 0.]), i
        assert np.array_equal(k, [0., 0., 1.]), k

        assert np.array_equal(origin, [1., 1., 1.]), origin
        assert np.array_equal(zaxis2, [1., 1., 2.]), zaxis2
        assert np.array_equal(xzplane, [2., 1., 1.5]), xzplane

    def test_wingbox(self):
        dirname = MODEL_PATH / 'wingbox'
        bdf_filename = dirname / 'wingbox_stitched_together-000.bdf'
        op2_filename = dirname / 'wingbox_stitched_together-000.op2'
        csv_filename = dirname / 'wingbox_stitched_together-000.csv'
        root_filename = dirname / 'wingbox'

        log = SimpleLogger(level='warning')
        include_results = 'grid_point_forces'
        bdf_model = read_bdf(bdf_filename, log=log)
        op2_model = read_op2(op2_filename, include_results=include_results, log=log)

        gpforce = op2_model.grid_point_forces[1]

        p1 = np.array([128., 6., 0.])  #  start; 127.95, 6.0, 19.559
        p3 = np.array([138., 96., 0.]) #  end; 119.95, 96., 21.832

        #method='Vector'
        method = 'vector'
        p2    = np.array([0., 1., 0.])  # xz-plane (normal direction)
        zaxis = np.array([0., 0., 1.])  # z-axis

        station_location = 'End-Origin'
        make_plot = IS_MATPLOTLIB
        forces_sum, moments_sum = create_shear_moment_torque(
            bdf_model, gpforce, p1, p2, p3, zaxis,
            method=method, station_location=station_location,
            cid_p1=0, cid_p2=0, cid_p3=0, cid_zaxis=0,
            nplanes=20, root_filename=None, csv_filename=None,
            length_scale=1.0, length_unit='in',
            force_scale=0.001, force_unit='kip',
            moment_scale=0.001, moment_unit='in-kip',
            make_plot=make_plot,
            show=False)

        if 0:
            station_location = 'y'
            forces_sum, moments_sum = create_shear_moment_torque(
                bdf_model, gpforce, p1, p2, p3, zaxis,
                method=method, station_location=station_location,
                cid_p1=0, cid_p2=0, cid_p3=0, cid_zaxis=0,
                nplanes=20,
                root_filename=root_filename,
                csv_filename=csv_filename,
                length_scale=1.0, length_unit='in',
                force_scale=0.001, force_unit='kip',
                moment_scale=0.001, moment_unit='in-kip',
                show=True)
            plt.close()

        #----------------------------------------------------------
        # spar 1
        p1 = np.array([110.95, 0.0, 22.8])
        p3 = np.array([110.95, 96., 22.8])
        station_location = 'y'

        method = 'vector'
        p2    = np.array([0., 1., 0.])  # xz-plane (normal direction)
        zaxis = np.array([0., 0., 1.])  # z-axis

        #1633:1856
        spar1_element_ids = np.arange(1633, 1856+1)
        forces_sum, moments_sum = create_shear_moment_torque(
            bdf_model, gpforce, p1, p2, p3, zaxis,
            method=method, station_location=station_location,
            cid_p1=0, cid_p2=0, cid_p3=0, cid_zaxis=0,
            nplanes=20,
            root_filename=root_filename,
            csv_filename=csv_filename,
            length_scale=1.0, length_unit='in',
            force_scale=0.001, force_unit='kip',
            moment_scale=0.001, moment_unit='in-kip',
            element_ids=spar1_element_ids,
            show=True)

    def test_cutting_plane_bwb(self):
        model = BDF(debug=False)
        model.cross_reference()
        p1 = np.array([1388.9, 1262.0, 86.3471])
        p2 = p1 + np.array([1., 0., 0.]) # xz-plane
        zaxis = np.array([0., 0., 1.])

        xyz1, xyz2, z_global, i, k, origin, zaxis2, xzplane = p1_p2_zaxis_to_cord2r(
            model, p1, p2, zaxis, method='Z-Axis Projection',
            cid_p1=0, cid_p2=0, cid_zaxis=0)
        assert np.array_equal(xyz1, p1), xyz1
        assert np.array_equal(xyz2, p2), xyz2
        assert np.array_equal(z_global, zaxis), z_global

        assert np.array_equal(i, [1., 0., 0.]), i
        assert np.array_equal(k, [0., 0., 1.]), k

        assert np.array_equal(origin, p1), origin
        assert np.array_equal(zaxis2, p1+zaxis), zaxis2
        assert np.array_equal(xzplane, p2), xzplane

    def test_get_stations_bwb(self):
        model = BDF(debug=False)
        model.cross_reference()
        p1 = np.array([1388.9, 1262.0, 86.3471])
        p3 = np.array([940.93, 1.7976e-07, 0.0])

        p2 = p1 + np.array([1., 0., 0.]) # xz-plane
        zaxis = np.array([0., 0., 1.])
        #xyz1, xyz2, z_global, i, k, origin, zaxis2, xzplane = p1_p2_zaxis_to_cord2r(
            #model, p1, p2, zaxis, method='Z-Axis Projection',
            #cid_p1=0, cid_p2=0, cid_zaxis=0)
        xyz1, xyz2, xyz3, i, k, coord_out, iaxis_march, stations = get_stations(
            model, p1, p2, p3, zaxis,
            method='Z-Axis Projection', cid_p1=0, cid_p2=0, cid_p3=0,
            cid_zaxis=0, nplanes=10)

        assert np.array_equal(xyz1, p1), xyz1
        assert np.array_equal(xyz2, p2), xyz2
        #assert np.array_equal(z_global, zaxis), z_global

        assert np.array_equal(i, [1., 0., 0.]), i
        assert np.array_equal(k, [0., 0., 1.]), k

        #assert np.array_equal(origin, p1), origin
        #assert np.array_equal(zaxis2, p1+zaxis), zaxis2
        #assert np.array_equal(xzplane, p2), xzplane
        #print(coord_out.e1, coord_out.e2, coord_out.e3)
        iaxis_march_expected = [-0.33382509, -0.94043632, -0.06434544]
        assert np.allclose(iaxis_march, iaxis_march_expected), f'iaxis_march={iaxis_march} expected={iaxis_march_expected}'
        #-----------------------------
        p1 = np.array([1388.9, 1262.0, 86.3471])
        p3 = np.array([940.93, 1.7976e-07, 0.0])

        p2 = np.array([1., 0., 0.]) # xz-plane
        zaxis = np.array([0., 0., 1.])
        xyz1, xyz2, xyz3, i, k, coord_out, iaxis_march, stations = get_stations(
            model, p1, p2, p3, zaxis,
            method='Vector', cid_p1=0, cid_p2=0, cid_p3=0,
            cid_zaxis=0, nplanes=10)
        assert np.allclose(iaxis_march, iaxis_march_expected), f'iaxis_march={iaxis_march} expected={iaxis_march_expected}'
        assert np.allclose(i, p2), f'i={i} expected={p2}'
        assert np.allclose(k, zaxis), f'k={k} expected={zaxis}'
        #-----------------------------
        p1 = np.array([1388.9, 1262.0, 86.3471])
        p3 = np.array([940.93, 1.7976e-07, 0.0])

        p2 = np.array([1., 0., 0.]) # xz-plane
        zaxis = np.array([0., 0., 0.])
        xyz1, xyz2, xyz3, i, k, coord_out, iaxis_march, stations = get_stations(
            model, p1, p2, p3, zaxis,
            method='Coord ID', cid_p1=0, cid_p2=0, cid_p3=0,
            cid_zaxis=0, nplanes=10)
        assert np.allclose(iaxis_march, iaxis_march_expected), f'iaxis_march={iaxis_march} expected={iaxis_march_expected}'
        assert np.allclose(i, p2), f'i={i} expected={p2}'
        assert np.allclose(k, [0., 0., 1.]), f'k={k} expected={[0., 0., 1.]}'

    def test_op2_bwb_smt_setup(self):  # pragma: no cover
        """how to plot an shear-moment-torque"""
        log = SimpleLogger(level='warning')
        bdf_filename = MODEL_PATH / 'bwb' / 'bwb_saero.bdf'
        model = BDF(debug=False, log=log)
        #model.load_as_h5 = True
        model.read_bdf(bdf_filename=bdf_filename)

        model.cross_reference()
        #gpforce = model.grid_point_forces[1]  # type: RealGridPointForcesArray

        #out = model.get_xyz_in_coord_array(cid=0)
        #nid_cp_cd, xyz_cid0, xyz_cp, icd_transform, icp_transform = out
        #nids = nid_cp_cd[:, 0]
        #nid_cd = nid_cp_cd[:, [0, 2]]
        #eids, element_centroids_cid0 = get_element_centroids(model)
        nids, nid_cd, xyz_cid0, icd_transform, eids, element_centroids_cid0 = smt_setup(model)
        nnodes = 10135
        nelements = 9424
        assert len(nids) == nnodes, len(nids)
        assert nid_cd.shape == (nnodes, 2), nid_cd.shape
        assert xyz_cid0.shape == (nnodes, 3), xyz_cid0.shape
        assert len(eids) == nelements, len(eids)
        assert element_centroids_cid0.shape == (nelements, 3), element_centroids_cid0.shape

    @unittest.skipIf(getpass.getuser() != 'sdoyle', 'local test')
    def test_op2_bwb_axial(self):  # pragma: no cover
        """how to plot an shear-moment-torque"""
        log = SimpleLogger(level='warning')
        BWB_PATH = MODEL_PATH / 'bwb'
        op2_filename = BWB_PATH / 'bwb_saero.op2'
        model = OP2Geom(debug=False, log=log, debug_file=None, mode=None)
        #model.load_as_h5 = True
        model.read_op2(op2_filename=op2_filename, combine=True,
                     build_dataframe=False, skip_undefined_matrices=False,
                     encoding=None)

        model.cross_reference()
        gpforce = model.grid_point_forces[1]  # type: RealGridPointForcesArray

        #out = model.get_xyz_in_coord_array(cid=0)
        #nid_cp_cd, xyz_cid0, xyz_cp, icd_transform, icp_transform = out
        #nids = nid_cp_cd[:, 0]
        #nid_cd = nid_cp_cd[:, [0, 2]]
        #eids, element_centroids_cid0 = get_element_centroids(model)
        nids, nid_cd, xyz_cid0, icd_transform, eids, element_centroids_cid0 = smt_setup(model)

        coord_out = model.coords[0]

        #cid_p1 = 0 # start
        #cid_p3 = 0 # end
        #cid_p2 = 0 # coord
        #p1-p2 defines the x-axis
        #k is defined by the z-axis
        #p1 = np.array([1354., 0., 0.]) # origin
        #p2 = np.array([1354., 1245., 0.]) # xaxis
        #p3 = np.array([1354., 1245., 0.]) # end
        #zaxis = np.array([0., 0., 1.])
        #method = 'Z-Axis Projection'
        #idir = 0

        #p1 = np.array([1354., 0., 0.]) # origin
        #p2 = np.array([1354., 0., 1.]) # xzplane
        #p3 = np.array([1354., 1245., 0.]) # end
        #zaxis = np.array([0., 0., 1.])
        #method = 'CORD2R'
        #idir = 1 # x-direction in this rotated system

        # axial
        p1 = np.array([0., 0., 0.]) # origin/start
        p2 = np.array([1600., 0., 0.]) # xaxis/end
        p3 = np.array([1600., 0., 0.]) # xzplane/end
        zaxis = np.array([0., 0., 1.])
        method = 'Z-Axis Projection'

        xyz1, xyz2, xyz3, i, k, coord_out, iaxis_march, stations = get_stations(
            model, p1, p2, p3, zaxis,
            method=method, cid_p1=0, cid_p2=0, cid_p3=0,
            cid_zaxis=0, nplanes=100)
        #print(f'stations = {stations}')

        # i/j/k vector is nan
        #print(f'origin: {coord_out.origin}')
        #print(f'zaxis: {coord_out.e2}')
        #print(f'xzplane: {coord_out.e3}')

        force_sum, moment_sum, new_coords, nelems, nnodes = gpforce.shear_moment_diagram(
            nids, xyz_cid0, nid_cd, icd_transform,
            eids, element_centroids_cid0,
            stations, model.coords, coord_out,
            iaxis_march=iaxis_march,
            itime=0, debug=True, log=model.log)
        plot_smt(stations, force_sum, moment_sum, nelems, nnodes, show=False)

    @unittest.skipIf(getpass.getuser() != 'sdoyle', 'local test')
    def test_op2_bwb_spanwise(self):  # pragma: no cover
        """how to plot an shear-moment-torque"""
        log_ = SimpleLogger(level='critical')
        BWB_PATH = MODEL_PATH / 'bwb'
        bdf_filename = BWB_PATH / 'bwb_saero_smt.bdf'
        op2_filename = BWB_PATH / 'bwb_saero.op2'
        with WarningRedirector(log_) as log:
            model = OP2Geom(debug=False, log=log, mode=None)
            model.include_exclude_results(include_results='grid_point_forces')
            model.read_op2(op2_filename=op2_filename)
            model.cross_reference()

        gpforce = model.grid_point_forces[1]  # type: RealGridPointForcesArray
        nids, nid_cd, xyz_cid0, icd_transform, eids, element_centroids_cid0 = smt_setup(model)

        #coord_out = model.coords[0]

        # spanwise
        # defines the starting point for the shear, moment, torque plot
        p1 = np.array([1388.9, 1262.0, 86.3471])  # 1/4c wing_tip

        # defines the end point for the shear, moment, torque plot
        p3 = np.array([910.93, 0.1, 0.0])  # 1/4c wing_root

        # y = 1262 (p3)
        # 0-0.1

        # defines the XZ plane for the shears/moments
        zaxis = p1 + np.array([0., 0., 1.])

        method = 'CORD2R'

        # coord_march:
        # - x: aft
        # - z: up

        # xz plane
        p2 = p1 + np.array([0., -1., 0.])

        xyz1, xyz2, xyz3, i, k, coord_out, iaxis_march, stations = get_stations(
            model, p1, p2, p3, zaxis,
            method=method, cid_p1=0, cid_p2=0, cid_p3=0,
            cid_zaxis=0, nplanes=50)

        nstations = len(stations)
        coord_origins = np.zeros((nstations, 3))
        for istation, station in enumerate(stations):
            coord_origin = coord_out.origin + iaxis_march * station
            coord_origins[istation, :] = coord_origin
        ylocations = coord_origins[:, 1]
        #coord_origins2 = coord_march.origin[:, np.newaxis] + coord_march.i[:, np.newaxis] * stations
        #assert np.array_equal(coord_origins, coord_origins2)

        #print(f'stations = {stations}')
        assert np.abs(stations).max() > 1350.
        assert np.array_equal(i, [0., -1., 0.])
        assert np.array_equal(k, [0., 0., 1.])

        # i/j/k vector is nan
        #print('coord_out:')
        #print(f'origin: {coord_out.origin}')
        #print(f'zaxis: {coord_out.e2}')
        #print(f'xzplane: {coord_out.e3}')
        #print(f'i: {coord_out.i}')
        #print(f'j: {coord_out.j}')
        #print(f'k: {coord_out.k}')

        #print('coord_march:')
        #print(f'i: {iaxis_march}')

        force_sum, moment_sum, new_coords, nelems, nnodes = gpforce.shear_moment_diagram(
            nids, xyz_cid0, nid_cd, icd_transform,
            eids, element_centroids_cid0,
            stations, model.coords, coord_out,
            iaxis_march=iaxis_march,
            itime=0, debug=True,
            nodes_tol=25., log=model.log)

        for cid, coord in new_coords.items():
            model.coords[cid] = coord
        model.sol = 144
        model.write_bdf(bdf_filename)
        plot_smt(ylocations, force_sum/1000, moment_sum/1000, nelems, nnodes,
                 plot_force_components=False,
                 plot_moment_components=False,
                 root_filename=os.path.join(BWB_PATH, 'bwb'),
                 show=False,
                 xtitle='Y', xlabel='Spanwise Location, Y (in)',
                 force_unit='kip', moment_unit='in-kip')
        y = 1


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

    def test_gpforce_02(self):
        IS_MATPLOTLIB = False
        log = SimpleLogger(level='warning')
        (model, coord_out, nid_cp_cd, icd_transform,
         all_nids, xyz_cid0,
         all_eids, element_centroids_cid0,
         gpforce, x1, x2, bending_moment2,
         stations) = _setup_bar_grid_point_forces(log)

        with self.assertRaises(AssertionError):
            log.error('problem with extract_freebody_loads...')
            fb_force, fb_moment = gpforce.extract_freebody_loads(
                all_eids,
                coord_out, model.coords,
                nid_cp_cd,
                icd_transform,
                itime=0, debug=True,
                log=log)

        if IS_MATPLOTLIB:  # pragma: no cover
            fig = plt.figure()
            ax = fig.gca()
            L = 10.0
            x = xyz_cid0[:, 0].copy()
            x.sort()
            M = x ** 2 / 2
            # F = wx
            # M = wx^2/2
            ax.plot(x1, bending_moment2[::2], 'o-', label='BM2', linewidth=3)
            ax.plot(x2, bending_moment2[1::2], 'o--', )
            ax.plot(L-x, M, 'o-', label='exact', linewidth=1)
            ax.grid(True)

        # trivial case
        #nids = [1]
        #eids = [1]
        force_out_sum, moment_out_sum = gpforce.extract_interface_loads(
            all_nids, all_eids,
            coord_out, model.coords,
            nid_cp_cd, icd_transform,
            xyz_cid0,
            #summation_point: Optional[NDArray3float]=None,
            consider_rxf=True, itime=0,
            stop_on_nan=True,
            debug=True, log=log)
        assert np.allclose(force_out_sum, [0., 0., 0.]), force_out_sum
        assert np.allclose(moment_out_sum, [0., 0., 0.]), moment_out_sum

        # some problematic case...
        eids = [1]
        nids = [2, 3, 4, 5, 6, 7, 8, 9, 10, 11]
        force_out_sum, moment_out_sum = gpforce.extract_interface_loads(
            nids, eids,
            coord_out, model.coords,
            nid_cp_cd, icd_transform,
            xyz_cid0,
            #summation_point: Optional[NDArray3float]=None,
            consider_rxf=True, itime=0,
            stop_on_nan=True,
            debug=True, log=log)
        assert np.allclose(force_out_sum, [0., 0., 9.5]), force_out_sum
        assert np.allclose(moment_out_sum, [0., -49.916668, 0.]), moment_out_sum

        # this one is empty...
        nids = [1]
        eids = [2]
        force_out_sum, moment_out_sum = gpforce.extract_interface_loads(
            nids, eids,
            coord_out, model.coords,
            nid_cp_cd, icd_transform,
            xyz_cid0,
            #summation_point: Optional[NDArray3float]=None,
            consider_rxf=False, itime=0,
            stop_on_nan=False,
            debug=False, log=log)

        # np.array([nan, nan, nan], dtype=float32)
        assert not np.any(np.isfinite(force_out_sum)), force_out_sum
        assert not np.any(np.isfinite(moment_out_sum)), moment_out_sum
        #coord0 = model.coords[0]
        #gpforce.extract_interface_loads(nids: np.ndarray,
                                        #eids: np.ndarray,
                                        #coord_out=coord0,
                                        #model.coords,
                                        #nid_cd,
                                        #icd_transform,
                                        #xyz_cid0,
                                        #summation_point=None,
                                        #consider_rxf=True,
                                        #itime=0,
                                        #debug=True,
                                        #log=model.log,
                                        #idtype='int32')
        # ----------------------------------------
        nodes_list = list(model.nodes.keys())
        nids = np.array(nodes_list, dtype='int32')
        nids.sort()
        #eids = np.ndarray(list(model.elements.keys()), dtype='int32')
        #eids.sort()
        # bar is [0,10] in x
        istations = np.where(stations > 0.5)[0]
        stations = stations[istations]
        force_sum, moment_sum, new_coords, nelems, nnodes = gpforce.shear_moment_diagram(
            nids, xyz_cid0, nid_cp_cd, icd_transform,
            all_eids, element_centroids_cid0,
            stations, model.coords, coord_out,
            #idir=0,
            itime=0, debug=True,
            nodes_tol=0.5,
            stop_on_nan=True, log=model.log)
        force_sum_expected = np.array([
            [0., 0., 9.5],
            [0., 0., 9.5],
            [0., 0., 9.5],
            [0., 0., 9.5],
            [0., 0., 9.5],
            [0., 0., 8.5],
            [0., 0., 8.5],
            [0., 0., 8.5],
            [0., 0., 8.5],
            [0., 0., 8.5],
            [0., 0., 7.5],
            [0., 0., 7.5],
            [0., 0., 7.5],
            [0., 0., 7.5],
            [0., 0., 7.5],
            [0., 0., 6.5],
            [0., 0., 6.5],
            [0., 0., 6.5],
            [0., 0., 6.5],
            [0., 0., 6.5],
            [0., 0., 5.5],
            [0., 0., 5.5],
            [0., 0., 5.5],
            [0., 0., 5.5],
            [0., 0., 5.5],
            [0., 0., 4.5],
            [0., 0., 4.5],
            [0., 0., 4.5],
            [0., 0., 4.5],
            [0., 0., 4.5],
            [0., 0., 3.5],
            [0., 0., 3.5],
            [0., 0., 3.5],
            [0., 0., 3.5],
            [0., 0., 3.5],
            [0., 0., 2.5],
            [0., 0., 2.5],
            [0., 0., 2.5],
            [0., 0., 2.5],
            [0., 0., 2.5],
            [0., 0., 1.5],
            [0., 0., 1.5],
            [0., 0., 1.5],
            [0., 0., 1.5],
            [0., 0., 1.5],
            [0., 0., 0.5],
            [0., 0., 0.5],
            [0., 0., 0.5],
        ])
        moment_sum_expected = np.array([
            [0., -4.42166672e+01, 0.],
            [0., -4.23166695e+01, 0.],
            [0., -4.04166679e+01, 0.],
            [0., -3.85166664e+01, 0.],
            [0., -3.66166687e+01, 0.],
            [0., -3.53166656e+01, 0.],
            [0., -3.36166649e+01, 0.],
            [0., -3.19166660e+01, 0.],
            [0., -3.02166653e+01, 0.],
            [0., -2.85166664e+01, 0.],
            [0., -2.74166660e+01, 0.],
            [0., -2.59166660e+01, 0.],
            [0., -2.44166660e+01, 0.],
            [0., -2.29166660e+01, 0.],
            [0., -2.14166660e+01, 0.],
            [0., -2.05166664e+01, 0.],
            [0., -1.92166653e+01, 0.],
            [0., -1.79166660e+01, 0.],
            [0., -1.66166668e+01, 0.],
            [0., -1.53166656e+01, 0.],
            [0., -1.46166668e+01, 0.],
            [0., -1.35166674e+01, 0.],
            [0., -1.24166670e+01, 0.],
            [0., -1.13166666e+01, 0.],
            [0., -1.02166672e+01, 0.],
            [0., -9.71666622e+00, 0.],
            [0., -8.81666660e+00, 0.],
            [0., -7.91666651e+00, 0.],
            [0., -7.01666641e+00, 0.],
            [0., -6.11666632e+00, 0.],
            [0., -5.81666660e+00, 0.],
            [0., -5.11666632e+00, 0.],
            [0., -4.41666651e+00, 0.],
            [0., -3.71666646e+00, 0.],
            [0., -3.01666641e+00, 0.],
            [0., -2.91666651e+00, 0.],
            [0., -2.41666651e+00, 0.],
            [0., -1.91666663e+00, 0.],
            [0., -1.41666663e+00, 0.],
            [0., -9.16666627e-01, 0.],
            [0., -1.01666665e+00, 0.],
            [0., -7.16666639e-01, 0.],
            [0., -4.16666657e-01, 0.],
            [0., -1.16666660e-01, 0.],
            [0.,  1.83333337e-01, 0.],
            [0., -1.16666667e-01, 0.],
            [0., -1.66666638e-02, 0.],
            [0.,  8.33333358e-02, 0.]])
        assert force_sum.shape == force_sum_expected.shape, force_out_sum
        assert moment_sum.shape == moment_sum_expected.shape, moment_sum.shape
        assert np.allclose(force_sum, force_sum_expected), force_out_sum
        assert np.allclose(moment_sum, moment_sum_expected), moment_sum
        if IS_MATPLOTLIB:  # pragma: no cover
            M2 = moment_sum[:, 1]
            ax.plot(stations, -M2, '*-', label='SMT')
            ax.legend()
            fig.show()
            x = 1

    def test_op2_solid_shell_bar_01_gpforce(self):
        folder = os.path.join(MODEL_PATH, 'sol_101_elements')
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
            out = gpforce._extract_interface_loads(
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
        folder = os.path.join(MODEL_PATH, 'sol_101_elements')
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
            out = gpforce._extract_interface_loads(
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

    def test_op2_solid_shell_bar_01_gpforce_radial_global_cd(self):
        warning_log = SimpleLogger(level='warning')
        debug_log = SimpleLogger(level='debug')
        folder = os.path.join(MODEL_PATH, 'sol_101_elements')
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
            out = gpforce._extract_interface_loads(
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
        folder = os.path.join(MODEL_PATH, 'sol_101_elements')
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

        assert op2_1.displacements[1].assert_equal(op2_2.displacements[1])
        assert op2_1.grid_point_forces[1].assert_equal(op2_2.grid_point_forces[1], atol=0.000123)

        # TODO: _extract_interface_loads has a separate bug with cylindrical coords
        return
        #-----------------------------------------------------------------------
        gpforce = op2_1.grid_point_forces[1]
        data = _get_gpforce_data()
        for i, datai in enumerate(data):
            eids, nids, cid, summation_point, total_force_local_expected, total_moment_local_expected = datai
            coord_out = op2_1.coords[cid]
            op2_1.log.debug('*' * 30 + 'Next Test #%s' % i + '*' * 30)
            out = gpforce._extract_interface_loads(
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

    def test_spherical_displacement_and_gpforce_transform(self):
        """Tests spherical coordinate transforms for displacements and GPF.

        Constructs a CORD2S at the global origin with identity beta, and one
        with a rotated beta. Places nodes at known spherical positions and
        verifies that vectors in spherical components (v_r, v_theta, v_phi)
        transform correctly to global rectangular (v_x, v_y, v_z).

        Tests:
        - Node on +z axis (theta=0): e_r=[0,0,1], e_theta=[1,0,0], e_phi=[0,1,0]
        - Node on +x axis (theta=90, phi=0): e_r=[1,0,0], e_theta=[0,0,-1], e_phi=[0,1,0]
        - Node on +y axis (theta=90, phi=90): e_r=[0,1,0], e_theta=[0,0,-1], e_phi=[-1,0,0]
        - Node at theta=90, phi=45: e_r=[c,c,0], e_theta=[0,0,-1], e_phi=[-c,c,0] (c=1/sqrt2)
        - Non-identity beta (90° rotation about z)
        """
        from pyNastran.op2.op2_interface.transforms import (
            _transform_spherical_displacement, _transform_spherical_gpforce)
        from pyNastran.bdf.cards.coordinate_systems import CORD2S

        # --- Case 1: identity beta (coord aligned with global) ---
        # CORD2S at origin with axes aligned to global
        cord2s = CORD2S(cid=10, rid=0,
                        origin=[0., 0., 0.],
                        zaxis=[0., 0., 1.],
                        xzplane=[1., 0., 0.])
        cord2s.setup()
        beta = cord2s.beta()
        assert np.allclose(beta, np.eye(3)), f'expected identity beta, got:\n{beta}'

        # 4 nodes at known spherical positions
        # Node 0: on +z axis (theta=0, phi=0) at r=1 -> xyz=(0, 0, 1)
        # Node 1: on +x axis (theta=90, phi=0) at r=1 -> xyz=(1, 0, 0)
        # Node 2: on +y axis (theta=90, phi=90) at r=1 -> xyz=(0, 1, 0)
        # Node 3: at theta=90, phi=45, r=1 -> xyz=(1/sqrt2, 1/sqrt2, 0)
        c = 1.0 / np.sqrt(2.0)
        xyz_cid0 = np.array([
            [0., 0., 1.],
            [1., 0., 0.],
            [0., 1., 0.],
            [c, c, 0.],
        ])
        inode = np.arange(4)

        # Displacement data: unit vector in e_r direction at each node
        # (1, 0, 0) in spherical = pure radial
        # ntimes=1, nnodes=4, 6 dof
        data = np.zeros((1, 4, 6), dtype='float64')
        data[0, :, 0] = 1.0  # v_r = 1 at all nodes

        _transform_spherical_displacement(
            inode, data, cord2s, xyz_cid0, beta, is_global_cid=True)

        # Expected: e_r in global for each node
        # Node 0 (theta=0, phi=0): e_r = [sin0*cos0, sin0*sin0, cos0] = [0, 0, 1]
        # Node 1 (theta=90, phi=0): e_r = [sin90*cos0, sin90*sin0, cos90] = [1, 0, 0]
        # Node 2 (theta=90, phi=90): e_r = [sin90*cos90, sin90*sin90, cos90] = [0, 1, 0]
        # Node 3 (theta=90, phi=45): e_r = [sin90*cos45, sin90*sin45, cos90] = [c, c, 0]
        expected_er = np.array([
            [0., 0., 1.],
            [1., 0., 0.],
            [0., 1., 0.],
            [c, c, 0.],
        ])
        assert np.allclose(data[0, :, :3], expected_er, atol=1e-14), (
            f'e_r transform failed:\n{data[0, :, :3]}\nexpected:\n{expected_er}')

        # Test e_theta: unit vector in theta direction
        data[:] = 0.
        data[0, :, 1] = 1.0  # v_theta = 1
        _transform_spherical_displacement(
            inode, data, cord2s, xyz_cid0, beta, is_global_cid=True)

        # Expected: e_theta in global
        # Node 0 (theta=0, phi=0): e_theta = [cos0*cos0, cos0*sin0, -sin0] = [1, 0, 0]
        # Node 1 (theta=90, phi=0): e_theta = [cos90*cos0, cos90*sin0, -sin90] = [0, 0, -1]
        # Node 2 (theta=90, phi=90): e_theta = [cos90*cos90, cos90*sin90, -sin90] = [0, 0, -1]
        # Node 3 (theta=90, phi=45): e_theta = [cos90*cos45, cos90*sin45, -sin90] = [0, 0, -1]
        expected_etheta = np.array([
            [1., 0., 0.],
            [0., 0., -1.],
            [0., 0., -1.],
            [0., 0., -1.],
        ])
        assert np.allclose(data[0, :, :3], expected_etheta, atol=1e-14), (
            f'e_theta transform failed:\n{data[0, :, :3]}\nexpected:\n{expected_etheta}')

        # Test e_phi: unit vector in phi direction
        data[:] = 0.
        data[0, :, 2] = 1.0  # v_phi = 1
        _transform_spherical_displacement(
            inode, data, cord2s, xyz_cid0, beta, is_global_cid=True)

        # Expected: e_phi in global
        # Node 0 (theta=0, phi=0): e_phi = [-sin0, cos0, 0] = [0, 1, 0]
        # Node 1 (theta=90, phi=0): e_phi = [-sin0, cos0, 0] = [0, 1, 0]
        # Node 2 (theta=90, phi=90): e_phi = [-sin90, cos90, 0] = [-1, 0, 0]
        # Node 3 (theta=90, phi=45): e_phi = [-sin45, cos45, 0] = [-c, c, 0]
        expected_ephi = np.array([
            [0., 1., 0.],
            [0., 1., 0.],
            [-1., 0., 0.],
            [-c, c, 0.],
        ])
        assert np.allclose(data[0, :, :3], expected_ephi, atol=1e-14), (
            f'e_phi transform failed:\n{data[0, :, :3]}\nexpected:\n{expected_ephi}')

        # --- Case 2: non-identity beta (coord rotated 90° about z) ---
        # x_coord = y_global, y_coord = -x_global, z_coord = z_global
        cord2s_rot = CORD2S(cid=11, rid=0,
                            origin=[0., 0., 0.],
                            zaxis=[0., 0., 1.],
                            xzplane=[0., 1., 0.])
        cord2s_rot.setup()
        beta_rot = cord2s_rot.beta()
        # beta should be [[0,1,0],[-1,0,0],[0,0,1]]
        expected_beta = np.array([[0., 1., 0.], [-1., 0., 0.], [0., 0., 1.]])
        assert np.allclose(beta_rot, expected_beta, atol=1e-14), (
            f'unexpected beta:\n{beta_rot}')

        # Node on +x axis in global: xyz=(1,0,0)
        # In coord-local rectangular: (1,0,0) @ beta.T = (0,-1,0)
        # In coord-local spherical: theta=90, phi=270 (or -90)
        # e_r at this point in coord-local rect = (0,-1,0) (unit vector toward the point)
        # e_r in global = beta.T @ (0,-1,0) = (0,-1,0) @ [[0,-1,0],[1,0,0],[0,0,1]] -> let's compute
        # Actually e_r in global: the point is at +x global, so radial = +x = (1,0,0)
        # Let's just verify the full transform numerically:
        xyz_single = np.array([[1., 0., 0.]])
        data_rot = np.zeros((1, 1, 6), dtype='float64')
        data_rot[0, 0, 0] = 1.0  # pure radial
        inode_single = np.array([0])
        _transform_spherical_displacement(
            inode_single, data_rot, cord2s_rot, xyz_single,
            beta_rot, is_global_cid=False)

        # The radial direction at (1,0,0) should point in +x global
        expected_rot = np.array([1., 0., 0.])
        assert np.allclose(data_rot[0, 0, :3], expected_rot, atol=1e-14), (
            f'rotated beta e_r failed: {data_rot[0, 0, :3]}, expected {expected_rot}')

        # --- Case 3: GPF transform (same math, different indexing) ---
        # Use identity beta case with inode_xyz mapping
        inode_xyz = np.arange(4)  # 1:1 mapping
        inode_gp = np.arange(4)
        data_gpf = np.zeros((1, 4, 6), dtype='float64')
        data_gpf[0, :, 0] = 1.0  # radial force = 1

        _transform_spherical_gpforce(
            inode_xyz, inode_gp, data_gpf, beta,
            cord2s, xyz_cid0, SimpleLogger(level='warning'))

        # Same expected as displacement e_r case
        assert np.allclose(data_gpf[0, :, :3], expected_er, atol=1e-14), (
            f'GPF spherical transform failed:\n{data_gpf[0, :, :3]}\nexpected:\n{expected_er}')

    def test_spherical_pressure_vessel(self):
        """Thin-walled pressure vessel: constant radial displacement.

        A constant internal pressure on a sphere produces uniform radial
        deformation. Every node has v_r=0.1, v_theta=0, v_phi=0 in spherical.
        After transform to global, each displacement should be 0.1 * r_hat,
        where r_hat = xyz / |xyz|.

        Uses 8 nodes distributed over the sphere (octants) with both
        identity beta and a rotated coord (offset origin + rotated axes).
        """
        from pyNastran.op2.op2_interface.transforms import _transform_spherical_displacement
        from pyNastran.bdf.cards.coordinate_systems import CORD2S

        # --- Case 1: identity beta, origin at global origin ---
        cord2s = CORD2S(cid=20, rid=0,
                        origin=[0., 0., 0.],
                        zaxis=[0., 0., 1.],
                        xzplane=[1., 0., 0.])
        beta = cord2s.beta()

        # 8 nodes at various positions on a unit sphere
        c = 1.0 / np.sqrt(3.0)
        xyz_cid0 = np.array([
            [1., 0., 0.],
            [0., 1., 0.],
            [0., 0., 1.],
            [-1., 0., 0.],
            [0., -1., 0.],
            [0., 0., -1.],
            [c, c, c],
            [-c, c, -c],
        ])
        nnodes = len(xyz_cid0)
        inode = np.arange(nnodes)

        # Constant radial displacement = 0.1
        radial_disp = 0.1
        data = np.zeros((1, nnodes, 6), dtype='float64')
        data[0, :, 0] = radial_disp  # v_r = 0.1

        _transform_spherical_displacement(
            inode, data, cord2s, xyz_cid0, beta, is_global_cid=True)

        # Expected: 0.1 * unit_radial_vector at each node
        norms = np.linalg.norm(xyz_cid0, axis=1, keepdims=True)
        expected = radial_disp * xyz_cid0 / norms
        assert np.allclose(data[0, :, :3], expected, atol=1e-14), (
            f'pressure vessel (identity) failed:\n{data[0, :, :3]}\nexpected:\n{expected}')
        # Rotational DOFs should be zero
        assert np.allclose(data[0, :, 3:], 0., atol=1e-14)

        # --- Case 2: offset origin + rotated axes ---
        # Coord centered at (10, 5, 3) with z along global y
        cord2s_off = CORD2S(cid=21, rid=0,
                            origin=[10., 5., 3.],
                            zaxis=[10., 6., 3.],   # z = +y global
                            xzplane=[11., 5., 3.]) # x on xz-plane -> +x global
        beta_off = cord2s_off.beta()

        # Nodes on a sphere of radius 2 centered at (10, 5, 3)
        r = 2.0
        xyz_off = np.array([
            [10. + r, 5., 3.],
            [10. - r, 5., 3.],
            [10., 5. + r, 3.],
            [10., 5. - r, 3.],
            [10., 5., 3. + r],
            [10., 5., 3. - r],
            [10. + r*c, 5. + r*c, 3. + r*c],
            [10. - r*c, 5. - r*c, 3. + r*c],
        ])
        nnodes_off = len(xyz_off)
        inode_off = np.arange(nnodes_off)

        data_off = np.zeros((1, nnodes_off, 6), dtype='float64')
        data_off[0, :, 0] = radial_disp

        _transform_spherical_displacement(
            inode_off, data_off, cord2s_off, xyz_off,
            beta_off, is_global_cid=False)

        # Expected: 0.1 * radial direction from the coord origin
        radial_vecs = xyz_off - cord2s_off.origin
        norms_off = np.linalg.norm(radial_vecs, axis=1, keepdims=True)
        expected_off = radial_disp * radial_vecs / norms_off
        assert np.allclose(data_off[0, :, :3], expected_off, atol=1e-14), (
            f'pressure vessel (offset) failed:\n{data_off[0, :, :3]}\nexpected:\n{expected_off}')
        assert np.allclose(data_off[0, :, 3:], 0., atol=1e-14)

    def test_spherical_gpforce_offset_rotated(self):
        """GPF transform for spherical coord with offset origin and rotated axes.

        Simulates uniform radial pressure on a sphere: each node has
        force = (F_r, 0, 0) in spherical. After transform to global,
        force at each node should be F_r * r_hat (pointing radially from
        the coord origin).

        Uses CORD2S with:
          - origin at (5, -3, 2)
          - z-axis along global (1, 1, 0)/sqrt(2) (45° tilt)
          - x on the xz-plane defined by a point

        This exercises origin subtraction, non-trivial beta, and the
        full GPF indexing path (inode_xyz -> inode_gp mapping).

        Tolerances: exact to machine precision (1e-14).
        """
        from pyNastran.op2.op2_interface.transforms import _transform_spherical_gpforce
        from pyNastran.bdf.cards.coordinate_systems import CORD2S

        # CORD2S with offset origin and tilted z-axis
        origin = np.array([5., -3., 2.])
        zaxis = origin + np.array([1., 1., 0.]) / np.sqrt(2.)
        xzplane = origin + np.array([1., 0., 0.])  # defines x in the xz-plane
        cord2s = CORD2S(cid=30, rid=0,
                        origin=origin.tolist(),
                        zaxis=zaxis.tolist(),
                        xzplane=xzplane.tolist())
        beta = cord2s.beta()

        # Verify beta is non-trivial
        assert not np.allclose(beta, np.eye(3)), 'beta should not be identity'

        # 6 nodes on a sphere of radius 3 centered at the coord origin
        r = 3.0
        c = 1.0 / np.sqrt(3.0)
        # Use directions that aren't aligned with coord axes
        directions = np.array([
            [1., 0., 0.],
            [0., 1., 0.],
            [0., 0., 1.],
            [-1., 0., 0.],
            [c, c, c],
            [-c, -c, c],
        ])
        xyz_cid0 = origin + r * directions
        nnodes = len(xyz_cid0)

        # GPF data: radial force = 100.0 at each node, 1 GPF entry per node
        # In a real model, multiple entries per node exist (one per element),
        # but this tests the core transform math.
        F_r = 100.0
        data = np.zeros((1, nnodes, 6), dtype='float64')
        data[0, :, 0] = F_r  # force in r-direction
        data[0, :, 4] = 50.0  # moment in theta-direction (arbitrary)

        # Index arrays: 1:1 mapping (each GPF row corresponds to one unique node)
        inode_xyz = np.arange(nnodes)
        inode_gp = np.arange(nnodes)

        _transform_spherical_gpforce(
            inode_xyz, inode_gp, data, beta,
            cord2s, xyz_cid0, SimpleLogger(level='warning'))

        # Expected forces: F_r * radial unit vector from coord origin
        radial_vecs = xyz_cid0 - origin
        norms = np.linalg.norm(radial_vecs, axis=1, keepdims=True)
        r_hat = radial_vecs / norms
        expected_force = F_r * r_hat
        assert np.allclose(data[0, :, :3], expected_force, atol=1e-13), (
            f'GPF radial force failed:\n{data[0, :, :3]}\nexpected:\n{expected_force}')

        # Expected moments: 50.0 in e_theta direction at each node
        # e_theta = [cos(theta)*cos(phi), cos(theta)*sin(phi), -sin(theta)]
        # in coord-local, then rotated to global via beta^T
        # Verify magnitude is preserved (50.0 at each node)
        moment_mags = np.linalg.norm(data[0, :, 3:], axis=1)
        assert np.allclose(moment_mags, 50.0, atol=1e-13), (
            f'moment magnitude not preserved: {moment_mags}')

        # Verify moments are perpendicular to radial direction
        # (e_theta is always perpendicular to e_r)
        dots = np.sum(data[0, :, 3:] * r_hat, axis=1)
        assert np.allclose(dots, 0., atol=1e-13), (
            f'moments not perpendicular to radial: dots={dots}')

    def test_spherical_gpforce_repeated_nodes(self):
        """GPF spherical transform with multiple entries per node.

        In a real model, each node has multiple GPF entries (one per
        connected element + SPC + total). All entries for the same node
        get the same rotation matrix. Verifies correct index handling.

        Uses CORD2S offset and rotated, 3 physical nodes, 2-3 GPF entries each.
        Tolerance: 1e-14 (exact math).
        """
        from pyNastran.op2.op2_interface.transforms import _transform_spherical_gpforce
        from pyNastran.bdf.cards.coordinate_systems import CORD2S

        origin = np.array([1., 2., 3.])
        cord2s = CORD2S(cid=31, rid=0,
                        origin=origin.tolist(),
                        zaxis=(origin + np.array([0., 0., 1.])).tolist(),
                        xzplane=(origin + np.array([1., 0., 0.])).tolist())
        beta = cord2s.beta()

        # 3 physical nodes
        r = 5.0
        xyz_cid0 = np.array([
            origin + [r, 0., 0.],
            origin + [0., r, 0.],
            origin + [0., 0., r],
        ])

        # GPF has 8 entries total: node0 has 3, node1 has 3, node2 has 2
        # Simulates: node0 (elem1, elem2, total), node1 (elem1, spc, total), node2 (elem1, total)
        ngpf = 8

        # inode_gp_xyz: maps each GPF row -> xyz_cid0 index (length = ngpf)
        # This is what searchsorted(nids_all, nids_all_gp) produces
        inode_gp_xyz = np.array([0, 0, 0, 1, 1, 1, 2, 2])

        # inode_gp: which GPF rows to transform (all of them in this case)
        inode_gp = np.arange(ngpf)

        # All entries get radial force = 1.0
        data = np.zeros((1, ngpf, 6), dtype='float64')
        data[0, :, 0] = 1.0  # F_r = 1.0

        _transform_spherical_gpforce(
            inode_gp_xyz, inode_gp, data, beta,
            cord2s, xyz_cid0, SimpleLogger(level='warning'))

        # Expected: all entries for node i get the same radial direction
        radial_vecs = xyz_cid0 - origin
        norms = np.linalg.norm(radial_vecs, axis=1, keepdims=True)
        r_hats = radial_vecs / norms

        # Check each GPF entry matches its node's radial direction
        node_map = [0, 0, 0, 1, 1, 1, 2, 2]
        for igpf, inode_phys in enumerate(node_map):
            assert np.allclose(data[0, igpf, :3], r_hats[inode_phys], atol=1e-14), (
                f'GPF entry {igpf} (node {inode_phys}): '
                f'{data[0, igpf, :3]} != {r_hats[inode_phys]}')

    def test_transform_displacement_full_pipeline(self):
        """Full pipeline test through transform_displacement_to_global.

        Exercises the full dispatch path:
          transform_displacement_to_global -> _transform_coord_nids
            -> _transform_cylindrical_displacement / _transform_spherical_displacement

        Tests:
        - Multiple time steps (ntimes=3)
        - Spherical coord with offset origin
        - Cylindrical coord with offset origin
        - Mixed: some nodes in cid=0 (no transform), others in cylindrical/spherical
        - Verifies each timestep independently

        Tolerance: 1e-13 (numerical precision).
        """
        from pyNastran.op2.op2_interface.transforms import transform_displacement_to_global
        from pyNastran.bdf.cards.coordinate_systems import CORD2C, CORD2S

        # --- Setup coords ---
        # cid=0: global (identity, no transform)
        coord0 = CORD2R(cid=0, origin=None, zaxis=None, xzplane=None)

        # cid=10: cylindrical, offset origin at (5, 0, 0), axis along global z
        origin_cyl = np.array([5., 0., 0.])
        cord2c = CORD2C(cid=10, rid=0,
                        origin=origin_cyl.tolist(),
                        zaxis=(origin_cyl + [0., 0., 1.]).tolist(),
                        xzplane=(origin_cyl + [1., 0., 0.]).tolist())

        # cid=20: spherical, offset origin at (-2, 3, 1), z along global (0,0,1)
        origin_sph = np.array([-2., 3., 1.])
        cord2s = CORD2S(cid=20, rid=0,
                        origin=origin_sph.tolist(),
                        zaxis=(origin_sph + [0., 0., 1.]).tolist(),
                        xzplane=(origin_sph + [1., 0., 0.]).tolist())

        coords = {0: coord0, 10: cord2c, 20: cord2s}

        # --- Setup nodes ---
        # 8 nodes total:
        #   nodes 0-1: cid=0 (no transform)
        #   nodes 2-4: cid=10 (cylindrical, around offset origin)
        #   nodes 5-7: cid=20 (spherical, around offset origin)
        r_cyl = 3.0
        angles_cyl = np.array([0., np.pi/3, 2*np.pi/3])
        xyz_cid0 = np.array([
            # cid=0 nodes
            [1., 2., 3.],
            [4., 5., 6.],
            # cid=10 nodes (r=3 around (5,0,0), z-axis = global z)
            [origin_cyl[0] + r_cyl, origin_cyl[1], origin_cyl[2]],  # theta=0
            [origin_cyl[0] + r_cyl*np.cos(angles_cyl[1]), origin_cyl[1] + r_cyl*np.sin(angles_cyl[1]), origin_cyl[2] + 1.],
            [origin_cyl[0] + r_cyl*np.cos(angles_cyl[2]), origin_cyl[1] + r_cyl*np.sin(angles_cyl[2]), origin_cyl[2] - 1.],
            # cid=20 nodes (on sphere of r=4 around (-2,3,1))
            [origin_sph[0] + 4., origin_sph[1], origin_sph[2]],  # theta=90, phi=0
            [origin_sph[0], origin_sph[1] + 4., origin_sph[2]],  # theta=90, phi=90
            [origin_sph[0], origin_sph[1], origin_sph[2] + 4.],  # theta=0
        ])
        nnodes = len(xyz_cid0)

        # icd_transform: maps cid -> node indices
        icd_transform = {
            0: np.array([0, 1]),
            10: np.array([2, 3, 4]),
            20: np.array([5, 6, 7]),
        }

        # --- Setup mock result data ---
        ntimes = 3

        # Create a mock result with .data attribute
        class MockResult:
            def __init__(self, data):
                self.data = data
                self.table_name = 'OUGV1'

        # All nodes get v_r (or v_x for rect) = 1.0 * (itime+1) as a scaling test
        data = np.zeros((ntimes, nnodes, 6), dtype='float64')
        for itime in range(ntimes):
            scale = float(itime + 1)
            # cid=0 nodes: displacement in x (stays as-is)
            data[itime, 0, 0] = 7.0 * scale
            data[itime, 1, 1] = -3.0 * scale
            # cid=10 nodes: radial displacement
            data[itime, 2, 0] = 0.1 * scale
            data[itime, 3, 0] = 0.1 * scale
            data[itime, 4, 0] = 0.1 * scale
            # cid=20 nodes: radial displacement
            data[itime, 5, 0] = 0.2 * scale
            data[itime, 6, 0] = 0.2 * scale
            data[itime, 7, 0] = 0.2 * scale

        result = MockResult(data)
        log = SimpleLogger(level='warning')

        # --- Run full transform ---
        transform_displacement_to_global(
            None, result, icd_transform, coords, xyz_cid0, log)

        # --- Verify ---
        for itime in range(ntimes):
            scale = float(itime + 1)

            # cid=0 nodes: unchanged
            assert np.isclose(data[itime, 0, 0], 7.0 * scale)
            assert np.isclose(data[itime, 1, 1], -3.0 * scale)

            # cid=10 nodes: radial displacement 0.1*scale in r -> global radial from cyl origin
            for i, a in zip([2, 3, 4], angles_cyl):
                r_hat = np.array([np.cos(a), np.sin(a), 0.])
                expected = 0.1 * scale * r_hat
                assert np.allclose(data[itime, i, :3], expected, atol=1e-13), (
                    f'itime={itime} node={i}: {data[itime, i, :3]} != {expected}')

            # cid=20 nodes: radial displacement 0.2*scale in r -> global radial from sph origin
            for i in [5, 6, 7]:
                r_vec = xyz_cid0[i] - origin_sph
                r_hat = r_vec / np.linalg.norm(r_vec)
                expected = 0.2 * scale * r_hat
                assert np.allclose(data[itime, i, :3], expected, atol=1e-13), (
                    f'itime={itime} node={i}: {data[itime, i, :3]} != {expected}')

    def test_cylindrical_pressure_vessel(self):
        """Cylindrical pressure vessel: constant radial displacement.

        Internal pressure on a cylinder produces uniform radial deformation
        in the r-direction. Every node has v_r=0.1, v_theta=0, v_z=0.
        After transform to global, each displacement should be
        0.1 * r_hat_xy, where r_hat_xy is the radial unit vector in the
        x-y plane (perpendicular to the cylinder axis).

        Tests identity beta, then offset + rotated coord.
        Tolerance: 1e-14 (exact math).
        """
        from pyNastran.op2.op2_interface.transforms import _transform_cylindrical_displacement
        from pyNastran.bdf.cards.coordinate_systems import CORD2C

        # --- Case 1: identity beta, origin at global origin ---
        cord2c = CORD2C(cid=40, rid=0,
                        origin=[0., 0., 0.],
                        zaxis=[0., 0., 1.],
                        xzplane=[1., 0., 0.])
        beta = cord2c.beta()
        assert np.allclose(beta, np.eye(3)), f'expected identity beta:\n{beta}'

        # 6 nodes around a cylinder of radius 2 at various z heights
        r = 2.0
        angles_deg = [0., 45., 90., 135., 225., 315.]
        angles_rad = np.radians(angles_deg)
        xyz_cid0 = np.array([
            [r * np.cos(a), r * np.sin(a), z]
            for a, z in zip(angles_rad, [0., 1., -1., 2., 0.5, -0.5])
        ])
        nnodes = len(xyz_cid0)
        inode = np.arange(nnodes)

        # Constant radial displacement = 0.1
        radial_disp = 0.1
        data = np.zeros((1, nnodes, 6), dtype='float64')
        data[0, :, 0] = radial_disp  # v_r = 0.1

        _transform_cylindrical_displacement(
            inode, data, cord2c, xyz_cid0, beta, is_global_cid=True)

        # Expected: 0.1 * radial unit vector in x-y plane
        xy = xyz_cid0[:, :2]
        r_xy = np.linalg.norm(xy, axis=1, keepdims=True)
        r_hat_xy = xy / r_xy
        expected = np.zeros((nnodes, 3))
        expected[:, :2] = radial_disp * r_hat_xy
        assert np.allclose(data[0, :, :3], expected, atol=1e-14), (
            f'cylinder (identity) failed:\n{data[0, :, :3]}\nexpected:\n{expected}')
        assert np.allclose(data[0, :, 3:], 0., atol=1e-14)

        # --- Case 2: tangential displacement (hoop stress) ---
        # v_theta = 0.05 at all nodes -> should be tangential (perpendicular to r_hat in x-y)
        data[:] = 0.
        data[0, :, 1] = 0.05  # v_theta

        _transform_cylindrical_displacement(
            inode, data, cord2c, xyz_cid0, beta, is_global_cid=True)

        # Expected: 0.05 * theta_hat, where theta_hat = [-sin(theta), cos(theta), 0]
        theta_hat = np.zeros((nnodes, 3))
        theta_hat[:, 0] = -np.sin(angles_rad)
        theta_hat[:, 1] = np.cos(angles_rad)
        expected_theta = 0.05 * theta_hat
        assert np.allclose(data[0, :, :3], expected_theta, atol=1e-14), (
            f'tangential (identity) failed:\n{data[0, :, :3]}\nexpected:\n{expected_theta}')

        # --- Case 3: axial displacement ---
        # v_z = 0.2 at all nodes -> should remain purely in z
        data[:] = 0.
        data[0, :, 2] = 0.2  # v_z

        _transform_cylindrical_displacement(
            inode, data, cord2c, xyz_cid0, beta, is_global_cid=True)

        expected_z = np.zeros((nnodes, 3))
        expected_z[:, 2] = 0.2
        assert np.allclose(data[0, :, :3], expected_z, atol=1e-14), (
            f'axial (identity) failed:\n{data[0, :, :3]}\nexpected:\n{expected_z}')

        # --- Case 4: offset origin + rotated axes ---
        # Cylinder axis along global x, centered at (3, -1, 5)
        origin = np.array([3., -1., 5.])
        cord2c_rot = CORD2C(cid=41, rid=0,
                            origin=origin.tolist(),
                            zaxis=(origin + np.array([1., 0., 0.])).tolist(),
                            xzplane=(origin + np.array([0., 0., 1.])).tolist())
        beta_rot = cord2c_rot.beta()
        assert not np.allclose(beta_rot, np.eye(3))

        # Nodes around the cylinder (axis = local z = global x)
        # In local frame: r=2, z varies along global x
        # local x -> global z, local y -> global y, local z -> global x
        r = 2.0
        local_angles_rad = np.array([0., np.pi/3, 2*np.pi/3, np.pi, 4*np.pi/3, 5*np.pi/3])
        # local coords: (r*cos(a), r*sin(a), z_local) -> transform to global via beta^T + origin
        z_locals = np.array([0., 1., 2., -1., 0.5, -0.5])
        xyz_local_rect = np.column_stack([
            r * np.cos(local_angles_rad),
            r * np.sin(local_angles_rad),
            z_locals
        ])
        # global = local_rect @ beta + origin
        xyz_rot = xyz_local_rect @ beta_rot + origin
        nnodes_rot = len(xyz_rot)
        inode_rot = np.arange(nnodes_rot)

        # Radial displacement = 0.1
        data_rot = np.zeros((1, nnodes_rot, 6), dtype='float64')
        data_rot[0, :, 0] = radial_disp

        _transform_cylindrical_displacement(
            inode_rot, data_rot, cord2c_rot, xyz_rot,
            beta_rot, is_global_cid=False)

        # Expected: 0.1 * radial direction in global
        # Radial direction in local rect = (cos(a), sin(a), 0), transform to global via beta^T
        r_hat_local = np.column_stack([
            np.cos(local_angles_rad),
            np.sin(local_angles_rad),
            np.zeros(nnodes_rot)
        ])
        r_hat_global = r_hat_local @ beta_rot  # local rect -> global
        expected_rot = radial_disp * r_hat_global
        assert np.allclose(data_rot[0, :, :3], expected_rot, atol=1e-14), (
            f'cylinder (offset/rotated) failed:\n{data_rot[0, :, :3]}\nexpected:\n{expected_rot}')

    def test_cylindrical_gpforce_offset_rotated(self):
        """GPF transform for cylindrical coord with offset origin and rotated axes.

        Simulates radial pressure on a cylinder: each node has force = (F_r, 0, 0)
        in cylindrical. After transform to global, force should point radially
        outward from the cylinder axis (perpendicular to the coord z-axis).

        Also tests tangential forces and the repeated-node indexing path.

        Tolerance: 1e-13.
        """
        from pyNastran.op2.op2_interface.transforms import _transform_cylindrical_gpforce
        from pyNastran.bdf.cards.coordinate_systems import CORD2C

        # Cylinder axis along global (1,1,0)/sqrt(2), centered at (2, 3, -1)
        origin = np.array([2., 3., -1.])
        z_dir = np.array([1., 1., 0.]) / np.sqrt(2.)
        cord2c = CORD2C(cid=42, rid=0,
                        origin=origin.tolist(),
                        zaxis=(origin + z_dir).tolist(),
                        xzplane=(origin + np.array([0., 0., 1.])).tolist())
        beta = cord2c.beta()
        assert not np.allclose(beta, np.eye(3))

        # 6 nodes around the cylinder
        r = 4.0
        local_angles_rad = np.linspace(0, 2*np.pi, 6, endpoint=False)
        z_locals = np.array([0., 1., -1., 2., -2., 0.5])
        xyz_local_rect = np.column_stack([
            r * np.cos(local_angles_rad),
            r * np.sin(local_angles_rad),
            z_locals
        ])
        xyz_cid0 = xyz_local_rect @ beta + origin
        nnodes = len(xyz_cid0)

        # GPF: 10 entries total (some nodes have multiple entries)
        # node 0: 3 entries, node 1: 2, node 2: 1, node 3: 2, node 4: 1, node 5: 1
        ngpf = 10
        inode_gp_xyz = np.array([0, 0, 0, 1, 1, 2, 3, 3, 4, 5])
        inode_gp = np.arange(ngpf)

        # Radial force = 100
        F_r = 100.0
        data = np.zeros((1, ngpf, 6), dtype='float64')
        data[0, :, 0] = F_r

        _transform_cylindrical_gpforce(
            inode_gp_xyz, inode_gp, data, beta,
            cord2c, xyz_cid0, SimpleLogger(level='warning'))

        # Expected: F_r * radial unit vector in global
        # Radial in local rect = (cos(a), sin(a), 0)
        r_hat_local = np.column_stack([
            np.cos(local_angles_rad),
            np.sin(local_angles_rad),
            np.zeros(nnodes)
        ])
        r_hat_global = r_hat_local @ beta

        # Map GPF entries to their node's expected direction
        node_map = [0, 0, 0, 1, 1, 2, 3, 3, 4, 5]
        for igpf, inode_phys in enumerate(node_map):
            expected_f = F_r * r_hat_global[inode_phys]
            assert np.allclose(data[0, igpf, :3], expected_f, atol=1e-13), (
                f'GPF entry {igpf} (node {inode_phys}): '
                f'{data[0, igpf, :3]} != {expected_f}')

        # --- Also test tangential force ---
        data[:] = 0.
        data[0, :, 1] = F_r  # F_theta

        _transform_cylindrical_gpforce(
            inode_gp_xyz, inode_gp, data, beta,
            cord2c, xyz_cid0, SimpleLogger(level='warning'))

        # theta_hat in local rect = (-sin(a), cos(a), 0)
        t_hat_local = np.column_stack([
            -np.sin(local_angles_rad),
            np.cos(local_angles_rad),
            np.zeros(nnodes)
        ])
        t_hat_global = t_hat_local @ beta

        for igpf, inode_phys in enumerate(node_map):
            expected_t = F_r * t_hat_global[inode_phys]
            assert np.allclose(data[0, igpf, :3], expected_t, atol=1e-13), (
                f'GPF tangential entry {igpf} (node {inode_phys}): '
                f'{data[0, igpf, :3]} != {expected_t}')

        # --- Axial force should align with coord z = global z_dir ---
        data[:] = 0.
        data[0, :, 2] = F_r  # F_z

        _transform_cylindrical_gpforce(
            inode_gp_xyz, inode_gp, data, beta,
            cord2c, xyz_cid0, SimpleLogger(level='warning'))

        # z_hat in local = (0, 0, 1), in global = (0,0,1) @ beta = beta[2,:]
        z_hat_global = beta[2, :]  # 3rd row of beta
        expected_axial = F_r * z_hat_global
        for igpf in range(ngpf):
            assert np.allclose(data[0, igpf, :3], expected_axial, atol=1e-13), (
                f'GPF axial entry {igpf}: {data[0, igpf, :3]} != {expected_axial}')

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

def _setup_bar_grid_point_forces(log):
    op2_filename = os.path.join(MODEL_PATH, 'grid_point_forces', 'bar_grid_point_forces.op2')
    #from pyNastran.bdf.bdf import read_bdf
    #bdf_model = read_bdf()
    model = read_op2(op2_filename, load_geometry=True, combine=True,
                     exclude_results=None, log=log)
    log = model.log
    gpforce = model.grid_point_forces[1]  # type: RealGridPointForcesArray
    forces = model.op2_results.force
    force = forces.cbar_force[1]
    #['station', 'bending_moment1', 'bending_moment2', 'shear1', 'shear2', 'axial', 'torque']
    headers = force.get_headers()
    #istation = headers.index('station')
    itime = 0
    #ibending_moment1 = headers.index('bending_moment1')
    ibending_moment2 = headers.index('bending_moment2')
    #station = force.data[itime, :, istation]
    #bending_moment1 = force.data[itime, :, ibending_moment1]
    bending_moment2 = force.data[itime, :, ibending_moment2]

    coord_out = model.coords[0]
    nid_cp_cd, xyz_cid0, xyz_cp, icd_transform, icp_transform = model.get_xyz_in_coord_array(
        cid=0, fdtype='float64', idtype='int32')
    all_nids = nid_cp_cd[:, 0]

    all_eids, element_centroids_cid0 = get_element_centroids(model, idtype='int32', fdtype='float64')
    #stations = element_centroids_cid0[:-1, 0]
    stations = np.linspace(0., 10., num=51)
    #model.log.level = 'warning'
    #print(stations)

    nids_bar = []
    nids_beam = []
    for eid, elem in sorted(model.elements.items()):
        if elem.type == 'CBAR':
            nids_bar.append(elem.nodes)
        elif elem.type == 'BEAM':
            nids_beam.append(elem.nodes)
    nids_bar = np.array(nids_bar, dtype='int32')
    nids_beam = np.array(nids_beam, dtype='int32')
    inid_bar = np.searchsorted(all_nids, nids_bar)
    x1 = xyz_cid0[inid_bar[:, 0], 0]
    x2 = xyz_cid0[inid_bar[:, 1], 0]
    out = (
        model, coord_out, nid_cp_cd, icd_transform,
        all_nids, xyz_cid0,
        all_eids, element_centroids_cid0,
        gpforce, x1, x2, bending_moment2,
        stations)
    return out

if __name__ == "__main__":  # pragma: no cover
    unittest.main()
