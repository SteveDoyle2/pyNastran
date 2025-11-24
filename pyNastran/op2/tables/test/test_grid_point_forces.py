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

from pyNastran.bdf.mesh_utils.cut_model_by_plane import (
    get_element_centroids, _p1_p2_zaxis_to_cord2r,
    get_stations,
)

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
        xyz1, xyz2, z_global, i, k, origin, zaxis2, xzplane = _p1_p2_zaxis_to_cord2r(
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
        xyz1, xyz2, z_global, i, k, origin, zaxis2, xzplane = _p1_p2_zaxis_to_cord2r(
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

        xyz1, xyz2, z_global, i, k, origin, zaxis2, xzplane = _p1_p2_zaxis_to_cord2r(
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
        #xyz1, xyz2, z_global, i, k, origin, zaxis2, xzplane = _p1_p2_zaxis_to_cord2r(
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

        # array([nan, nan, nan], dtype=float32)
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

    @unittest.expectedFailure
    def test_broken_op2_solid_shell_bar_01_gpforce_radial_global_cd(self):
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
