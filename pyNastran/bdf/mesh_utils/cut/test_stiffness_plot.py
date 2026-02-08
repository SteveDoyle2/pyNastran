"""defines cutting plane tests"""
import os
# import copy
import time
from pathlib import Path
import unittest
import numpy as np
from pyNastran.gui.matplotlib_backend import matplotlib_backend

try:
    import matplotlib  # pylint: disable=unused-import
    IS_MATPLOTLIB = True
except ModuleNotFoundError:  # pragma: no cover
    IS_MATPLOTLIB = False

if IS_MATPLOTLIB:
    matplotlib.use(matplotlib_backend)

import pyNastran
from pyNastran.bdf.bdf import read_bdf, BDF, CORD2R
from cpylog import SimpleLogger

from pyNastran.bdf.mesh_utils.cut.moi_plotter import (
    cut_and_plot_moi, plot_inertia)
from pyNastran.bdf.mesh_utils.cut.cut_model_by_plane import (
    _setup_faces)

PKG_PATH = pyNastran.__path__[0]
TEST_PATH = Path(__file__).parent
MODEL_PATH = Path(os.path.join(PKG_PATH, '..', 'models'))


class TestStiffnessPlot(unittest.TestCase):
    def test_cut_quad(self):
        """cut_and_plot_moi"""
        dirname = TEST_PATH
        log = SimpleLogger(level='debug', encoding='utf-8')
        model = BDF(log=log)
        model.add_grid(101, [0., 0., 0.])
        model.add_grid(102, [3., 0., 0.])
        model.add_grid(103, [3., 3., 0.])
        model.add_grid(104, [0., 3., 0.])
        model.add_cquad4(10, 11, [101, 102, 103, 104])
        t = 0.1
        E = 3.0e7
        model.add_pshell(11, 12, t=t)
        model.add_mat1(12, E=E, G=None, nu=0.3)
        dy = 0.5
        coord = model.add_cord2r(
            cid=1,
            origin=[0., dy, 0.],
            zaxis=[0., dy, 1.],
            xzplane=[1., dy, 0.])
        model.cross_reference()
        coords = [coord]
        ystations = [0.]
        normal_plane = np.array([0., 1., 0.])

        moi_data = cut_and_plot_moi(
            model, normal_plane, log,
            ystations, coords,
            ytol=2.0,
            dirname=dirname,
            plot=False, show=False, face_data=None,
            stop_on_failure=True,
            cut_data_span_filename='',
            beam_model_bdf_filename='',
            thetas_csv_filename='y_thetas.csv',
            # normalized_inertia_png_filename='y_normalized_inertia_vs_span.png',
            # area_span_png_filename='y_area_vs_span.png',
            # amoi_span_png_filename='y_amoi_vs_span.png',
            # e_amoi_span_png_filename='y_e_amoi_vs_span.png',
            # cg_span_png_filename='y_cg_vs_span.png',
            debug_vectorize=True,
        )
        (y, A, I, J,
         EI, GJ, avg_centroid,
         plane_bdf_filenames1, plane_bdf_filenames2) = moi_data
        # print(f'y = {y.tolist()}')
        # print(f'A = {A.tolist()}')
        # print(f'I = {I.tolist()}')
        # print(f'J = {J.tolist()}')
        # print(f'EI = {EI.tolist()}')
        # print(f'GJ = {GJ.tolist()}')
        print(f'avg_centroid = {avg_centroid.tolist()}')
        cut_length = 3
        # skipped this...1/12 * cut_length * t**3
        i_expected = 0.09375  # A*d^2 of 2 triangles
        y_expected = [0.0]
        A_expected = [t*cut_length]
        I_expected = [[i_expected, 0.0, 0.0, 0.0, 0.0, 0.0]]
        J_expected = [i_expected]
        EI_expected = [[E*i_expected, 0.0, 0.0, 0.0, 0.0, 0.0]]  # 2812500.0
        GJ_expected = [1.0]
        centroid_expected = [[1.5, 0.0, 0.0]]

        assert np.allclose(y, y_expected)
        assert np.allclose(A, A_expected)
        assert np.allclose(I, I_expected), (I, I_expected)
        assert np.allclose(J, J_expected)
        assert np.allclose(EI, EI_expected)
        assert np.allclose(GJ, GJ_expected)

    def test_cut_bwb(self):
        """cut_and_plot_moi"""
        # show = True
        t0 = time.time()
        show = False
        run_y_cuts = True
        run_x_cuts = False

        log = SimpleLogger(level='warning', encoding='utf-8')
        # log = SimpleLogger(level='debug', encoding='utf-8')
        dirname = MODEL_PATH / 'bwb'
        bdf_filename = dirname / 'bwb_saero.bdf'  # ymax~=1262.0
        model = read_bdf(bdf_filename, log=log)
        # model.log.level = 'debug'

        ymax = 1401.
        ncut = 10
        dy = (ymax-1.) / ncut
        i = np.arange(ncut, dtype='int32')
        # xstations = 50 * i + 1.
        ystations = dy * i + 1.

        # y is outboard
        origin = np.array([0., 0., 0.])
        zaxis = np.array([0., 0., 1.])   # z
        xzplane = np.array([1., 0., 0.]) # x

        coords = get_coords_bwb(
            ystations, cid=1,
            origin=origin,
            zaxis=zaxis,
            xzplane=xzplane)

        normal_plane = coords[0].j
        # normal_plane = np.array([0., 1., 0.])
        # assert np.allclose(normal_plane, normal_plane2)

        log, *face_data = _setup_faces(model)
        # nids, xyz_cid0, elements = face_data
        # y0, A0, I0, J0, EI0, J0, avg_centroid0, plane_bdf_filenames10, plane_bdf_filenames20 = cut_and_plot_moi(
        #     model, normal_plane, log,
        #     dys, coords,
        #     ytol=2.0,
        #     dirname=dirname,
        #     plot=False, show=False, face_data=face_data)

        if run_y_cuts:
            log.info('working on y-cuts')
            moi_data = cut_and_plot_moi(
                bdf_filename, normal_plane, log,
                ystations, coords,
                ytol=2.0,
                dirname=dirname,
                plot=True, show=False, face_data=face_data,
                cut_data_span_filename='y_cut_data_vs_span.csv',
                beam_model_bdf_filename='y_equivalent_beam_model.bdf',
                thetas_csv_filename='y_thetas.csv',
                normalized_inertia_png_filename='y_normalized_inertia_vs_span.png',
                area_span_png_filename='y_area_vs_span.png',
                amoi_span_png_filename='y_amoi_vs_span.png',
                e_amoi_span_png_filename='y_e_amoi_vs_span.png',
                cg_span_png_filename='y_cg_vs_span.png',
                debug_vectorize=True,
            )
            (y, A, I, J,
             EI, GJ, avg_centroid,
             plane_bdf_filenames1, plane_bdf_filenames2) = moi_data
            # assert np.allclose(avg_centroid, avg_centroid0)
            # print(f'y = {y.tolist()}')
            # print(f'A = {A.tolist()}')
            # print(f'I = {I.tolist()}')
            # print(f'J = {J.tolist()}')
            # print(f'EI = {EI.tolist()}')
            # print(f'GJ = {GJ.tolist()}')

            y_expected = [1.0, 141.0, 281.0, 421.0, 561.0, 701.0, 841.0, 981.0, 1121.0, 1261.0]
            A_expected = [27325.7035460869, 16791.75582131023, 3439.6593514378583, 1885.2005923508825, 1112.9466242866476, 902.0982620165437, 692.8361798071458, 528.0741288550782, 397.36522577804004, 286.59151276901025]
            I_expected = [
                [3403933746.398361, 8.911280281243883e-27, 87431348.83639666, -1.1524729838475616e-09, -8.115831147516347e-11, -113043322.10390155],
                [1333461947.3875391, 2.2380375734653757e-27, 41162522.72557803, -3.816419553672272e-11, 2.2490616251215358e-11, -16958977.338277496],
                [213026144.23155284, 2.421936632273299e-28, 4799720.5299459305, -2.789842444913242e-12, -3.8705946989251635e-14, -5341199.828407566],
                [41040943.73186442, 1.963256730304734e-28, 545179.3969024371, -1.2228253392113826e-11, -4.381822725941327e-13, 449897.1318028397],
                [10103921.858063448, 2.8359530263005473e-28, 145368.04628856102, 5.970341815480735e-12, 8.713118441540724e-13, 523741.2250703718],
                [5374113.109890341, 1.135309608747665e-28, 87417.81394571027, 1.069894873862433e-12, -1.6083252621051996e-13, 380076.317217249],
                [2845327.0112954406, 6.312458138991356e-29, 53886.33108173666, 2.2219427152776314e-12, 3.100104645122421e-13, 238803.7329846358],
                [1483124.33122053, 8.203848533496161e-29, 30610.86321437863, 1.167927493310214e-12, -8.845064798186607e-14, 142410.19408779056],
                [684110.0163668194, 8.231340803540231e-29, 16725.831058479114, 5.647390666346783e-12, 7.090225668624896e-13, 78866.11262783181],
                [250269.1763145865, 3.010055271280282e-29, 8678.51145346907, 1.2739523923193558e-12, 1.1959210151601072e-13, 38470.74439627658]]
            J_expected = [3491365095.234758, 1374624470.1131172, 217825864.76149878, 41586123.12876686, 10249289.90435201, 5461530.9238360515, 2899213.3423771774, 1513735.1944349087, 700835.8474252985, 258947.6877680556]
            EI_expected = [
                [3.295171322240443e+16, 1.0055942705118451e-19, 834725104132483.5, -0.016966902163427177, -0.0008386804640027322, -1305402543493497.0],
                [1.1424869043875884e+16, 2.436246336255071e-20, 346881503958272.25, -0.0014203565392383113, 0.00012928505137646288, -158076619489430.62],
                [1376220583763065.0, 1.4987514310082692e-21, 30847567187819.86, -6.739996447080615e-05, 2.686140388554683e-07, -37232864334729.03],
                [232849135638355.66, 1.163610376266595e-21, 3370933367866.2837, -6.256159445564982e-05, -2.130227964903126e-06, 2725803693395.542],
                [53434265948650.6, 1.5070800170642208e-21, 756376205291.4624, 2.8350427118543883e-05, 3.4861744237313086e-06, 2421036127784.6895],
                [28448320392572.47, 5.8563062886418695e-22, 452376143171.62006, 1.233289212310551e-05, -4.952364981663227e-07, 1838842284503.897],
                [11922707419038.45, 3.128149669510187e-22, 255514290711.5685, 1.4411389995587915e-05, 1.8297740043860598e-06, 1043116272629.821],
                [5081816000830.39, 3.04143608378661e-22, 110098452419.56664, 4.1518689145376174e-06, -2.915248889290899e-07, 483517837043.9404],
                [2127609200637.4373, 2.5769541338500123e-22, 52962312533.78477, 1.7402752844301667e-05, 2.1839472821258388e-06, 244954003117.0613],
                [778925018736.649, 9.542065160328333e-23, 27450908134.014584, 3.8466670979195605e-06, 3.395408960754233e-07, 119686895556.34503]]
            GJ_expected = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
            assert np.allclose(y, y_expected)
            assert np.allclose(A, A_expected)
            assert np.allclose(I, I_expected)
            assert np.allclose(J, J_expected)
            assert np.allclose(EI, EI_expected)
            assert np.allclose(GJ, GJ_expected)
            for plane_bdf_filename in plane_bdf_filenames1:
                os.remove(plane_bdf_filename)
            for plane_bdf_filename in plane_bdf_filenames2:
                os.remove(plane_bdf_filename)

            if IS_MATPLOTLIB:
                plot_inertia(y, A, I, J, EI, GJ, avg_centroid, show=show)

        if run_x_cuts:
            # xmax = 1001.
            ncut = 100
            # dx = (xmax-1.) / ncut
            i = np.arange(ncut, dtype='int32')
            xstations = 10 * i

            # y is outboard
            origin = np.array([0., 0., 0.])
            zaxis = np.array([0., 0., 1.])   # z
            xzplane = np.array([0., 1., 0.]) # y

            xcoords = get_coords_bwb(
                xstations, cid=1,
                origin=origin,
                zaxis=zaxis,
                xzplane=xzplane)

            log.info('working on x-cuts')
            normal_plane = coords[0].j
            log.debug(f'normal_plane = {normal_plane}')
            moi_data = cut_and_plot_moi(
                bdf_filename, normal_plane, log,
                xstations, xcoords,
                ytol=2.0,
                dirname=dirname, ifig=10,
                plot=True, show=False, face_data=face_data,
                cut_data_span_filename='x_cut_data_vs_span.csv',
                beam_model_bdf_filename='x_equivalent_beam_model.bdf',
                thetas_csv_filename='x_thetas.csv',
                normalized_inertia_png_filename='x_normalized_inertia_vs_span.png',
                area_span_png_filename='x_area_vs_span.png',
                amoi_span_png_filename='x_amoi_vs_span.png',
                e_amoi_span_png_filename='x_e_amoi_vs_span.png',
                cg_span_png_filename='x_cg_vs_span.png',
                debug_vectorize=True,
            )
            (x, A, I, J,
             EI, GJ, avg_centroid,
             plane_bdf_filenames1, plane_bdf_filenames2) = moi_data
            print(f'x = {x.tolist()}')
            print(f'A = {A.tolist()}')
            x_expected = []
            ax_expected = []
            assert np.allclose(x, x_expected)
            assert np.allclose(A, ax_expected)
            log.info('cleanup')
            for plane_bdf_filename in plane_bdf_filenames1:
                os.remove(plane_bdf_filename)
            for plane_bdf_filename in plane_bdf_filenames2:
                os.remove(plane_bdf_filename)
        if run_x_cuts:
            _cleanup_moi_files(dirname, 'x_')
        if run_y_cuts:
            _cleanup_moi_files(dirname, 'y_')
        print(f'dt = {time.time() - t0}')

def _cleanup_moi_files(dirname: Path,
                       tag: str):
    if IS_MATPLOTLIB:
        os.remove(dirname / f'{tag}normalized_inertia_vs_span.png')
        os.remove(dirname / f'{tag}area_vs_span.png')
        os.remove(dirname / f'{tag}amoi_vs_span.png')
        os.remove(dirname / f'{tag}e_amoi_vs_span.png')
        os.remove(dirname / f'{tag}cg_vs_span.png')

    # bdf_merge(plane_bdf_filenames, bdf_filename_out='merge.bdf', renumber=True,
    #           encoding=None, size=8, is_double=False, cards_to_skip=None,
    #           log=None, skip_case_control_deck=False)
    os.remove(dirname / f'{tag}thetas.csv')
    # os.remove(dirname / 'equivalent_beam_model.bdf')
    os.remove(dirname / f'{tag}cut_data_vs_span.csv')
    # os.remove('cut_face.csv')
    # if IS_MATPLOTLIB:
    #     os.remove('area_vs_span.png')
    #     os.remove('amoi_vs_span.png')
    #     os.remove('normalized_inertia_vs_span.png')
    #     os.remove('cg_vs_span.png')
    #     os.remove('e_amoi_vs_span.png')


def get_coords_bwb(ystations: np.ndarray,
                   cid: int=1,
                   base_coord: CORD2R=None,
                   origin: np.ndarray=None,
                   zaxis: np.ndarray=None,
                   xzplane: np.ndarray=None) -> list[CORD2R]:  # pragma: no cover
    """gets coords from y=0 to y=100*ncuts"""
    if base_coord:
        raise NotImplementedError('base_coord is not yet implemented')
    else:
        origin = np.asarray(origin, dtype='float64')
        zaxis = np.asarray(zaxis, dtype='float64')
        xzplane = np.asarray(xzplane, dtype='float64')

    coords = []
    nstation = len(ystations)
    dxyz = np.zeros((nstation, 3), dtype='float64')
    dxyz[:, 1] = ystations

    for dxyzi in dxyz:
        coord = CORD2R(1, rid=0,
                       origin=origin+dxyzi,
                       zaxis=zaxis+dxyzi,
                       xzplane=xzplane+dxyzi)
        coords.append(coord)
    return coords


if __name__ == '__main__':  # pragma: no cover
    unittest.main()
