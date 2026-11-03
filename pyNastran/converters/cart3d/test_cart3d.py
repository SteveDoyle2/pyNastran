"""tests non-gui related Cart3d class/interface"""
import os
from pathlib import Path
import unittest
import numpy as np
from cpylog import get_logger, SimpleLogger

import pyNastran
from pyNastran.converters.cart3d.cart3d import Cart3D, read_cart3d, comp2tri
from pyNastran.converters.cart3d.cart3d_to_nastran import cart3d_to_nastran_filename, cart3d_to_nastran_model
from pyNastran.converters.cart3d.cart3d_to_stl import cart3d_to_stl_filename
from pyNastran.converters.cart3d.cart3d_to_tecplot import cart3d_to_tecplot
from pyNastran.converters.cart3d.input_c3d_reader import read_input_c3d
from pyNastran.converters.cart3d.input_cntl_reader import read_input_cntl
from pyNastran.converters.format_converter import cmd_line_format_converter


PKG_PATH = Path(pyNastran.__path__[0])
MODEL_PATH = PKG_PATH / 'converters' / 'cart3d' / 'models'

class TestCart3d(unittest.TestCase):

    def test_cart3d_flat_geometry(self):
        """geometry"""
        lines = (
            "7 6\n"
            "0.000000 0.000000 0.000000\n"
            "1.000000 0.000000 0.000000\n"
            "2.000000 0.000000 0.000000\n"
            "1.000000 1.000000 0.000000\n"
            "2.000000 1.000000 0.000000\n"
            "1.000000 -1.000000 0.000000\n"
            "2.000000 -1.000000 0.000000\n"
            "1 4 2\n"
            "2 4 5\n"
            "2 5 3\n"
            "2 6 1\n"
            "5 6 2\n"
            "5 5 2\n"
            "1\n"
            "2\n"
            "3\n"
            "2\n"
            "4\n"
            "6\n"
        )
        log = get_logger(level='warning', encoding='utf-8')
        infile_name = os.path.join(MODEL_PATH, 'flat_full.tri')
        out_name = os.path.join(MODEL_PATH, 'combined.tri')
        with open(infile_name, 'w') as f:
            f.write(lines)

        model = comp2tri([infile_name, infile_name], out_name, is_binary=False, float_fmt='%6.7f',
                         log=log)
        read_cart3d(out_name, log=log, debug=False)

        cart3d = read_cart3d(infile_name, log=log, debug=False)
        assert len(cart3d.points) == 7, 'npoints=%s' % len(cart3d.points)
        assert len(cart3d.elements) == 6, 'nelements=%s' % len(cart3d.elements)
        assert len(cart3d.regions) == 6, 'nregions=%s' % len(cart3d.regions)
        assert len(cart3d.loads) == 0, 'nloads=%s' % len(cart3d.loads)

        response = np.sin(cart3d.points[:, 0])
        yslices = np.array([0.])
        cart3d.cut_model_node(response, yslices)

        edges = cart3d.get_free_edges(cart3d.elements)
        assert len(edges) == 8, len(edges)
        os.remove(out_name)

        argv = ['format_converter', 'cart3d', infile_name, 'cart3d', out_name, '-b']
        cmd_line_format_converter(argv=argv, quiet=True)

        os.remove(infile_name)
        os.remove(out_name)

    def test_cart3d_flat_results(self):
        """geometry + results"""
        lines = (
            "5 3 6\n"
            "0. 0. 0.\n"
            "1. 0. 0.\n"
            "2. 0. 0.\n"
            "1. 1. 0.\n"
            "2. 1. 0.\n"
            "1 4 2\n"
            "2 4 5\n"
            "2 5 3\n"
            "1\n"
            "2\n"
            "3\n"
            "1.\n"
            "1. 1. 1. 1. 1.\n"
            "2.\n"
            "2. 2. 2. 2. 2.\n"
            "3.\n"
            "3. 3. 3. 3. 3.\n"
            "4.\n"
            "4. 4. 4. 4. 4.\n"
            "5.\n"
            "5. 5. 5. 5. 5.\n"
        )
        cart3d_filename = os.path.join(MODEL_PATH, 'flat.tri')
        with open(cart3d_filename, 'w') as f:
            f.write(lines)

        log = get_logger(level='warning', encoding='utf-8')
        model = read_cart3d(cart3d_filename, log=log, debug=False,
                            result_names=None)

        assert len(model.points) == 5, 'npoints=%s' % len(model.points)
        assert len(model.elements) == 3, 'nelements=%s' % len(model.elements)
        assert len(model.regions) == 3, 'nregions=%s' % len(model.regions)

        assert len(model.loads) == 14, 'nloads=%s' % len(model.loads)  # was 10
        assert len(model.loads['Cp']) == 5, 'nCp=%s' % len(model.loads['Cp'])

        outfile_name1 = os.path.join(MODEL_PATH, 'flat.bin.tri')
        outfile_name2 = os.path.join(MODEL_PATH, 'flat2.tri')
        model.write_cart3d(outfile_name2, is_binary=False)
        with self.assertRaises(NotImplementedError):
            model.write_cart3d(outfile_name1, is_binary=True)

        # no results
        model.loads = {}
        model.write_cart3d(outfile_name1, is_binary=True)
        model.write_cart3d(outfile_name2, is_binary=False)

        area, centroid, normals = model.get_area_centroid_normals()
        cnormals = model.get_normals()
        nnormals = model.get_normals_at_nodes(cnormals)
        area, centroid, cnormals = model.get_area_centroid_normals()
        os.remove(cart3d_filename)
        os.remove(outfile_name1)
        os.remove(outfile_name2)

    def test_cart3d_io_plugs(self):
        """read/write geometry in ascii/binary"""
        log = get_logger(level='warning', encoding='utf-8')
        infile_name = os.path.join(MODEL_PATH, 'threePlugs.bin.tri')
        outfile_name = os.path.join(MODEL_PATH, 'threePlugs_out.tri')
        outfile_name_bin = os.path.join(MODEL_PATH, 'threePlugs_bin2.tri')
        outfile_name_bin_out = os.path.join(MODEL_PATH, 'threePlugs_bin_out.tri')

        cart3d = read_cart3d(infile_name, log=log, debug=False)
        cart3d.write_cart3d(outfile_name, is_binary=False)
        cart3d.write_cart3d(outfile_name_bin, is_binary=True)

        cart3d_ascii = read_cart3d(outfile_name, log=log, debug=False)
        check_array(cart3d.points, cart3d_ascii.points)
        check_array(cart3d.elements, cart3d_ascii.elements)

        cart3d_bin = read_cart3d(outfile_name_bin, log=log, debug=False)
        check_array(cart3d.points, cart3d_bin.points)
        check_array(cart3d.elements, cart3d_ascii.elements)
        #print(cart3d_bin.points)

        #print('---------------')
        #print(cart3d_bin.points)

        os.remove(outfile_name)
        os.remove(outfile_name_bin)

        cart3d.write_cart3d(outfile_name_bin_out, is_binary=False)
        os.remove(outfile_name_bin_out)

    def test_cart3d_to_stl(self):
        """convert to stl"""
        log = get_logger(level='warning', encoding='utf-8')
        cart3d_filename = os.path.join(MODEL_PATH, 'threePlugs.bin.tri')
        stl_filename = os.path.join(MODEL_PATH, 'threePlugs.stl')

        argv = ['format_converter', 'cart3d', cart3d_filename,
                'stl', stl_filename, '-b', '--scale', '2.0']
        cmd_line_format_converter(argv=argv, quiet=True)
        #cart3d_to_stl_filename(cart3d_filename, stl_filename, log=log)
        os.remove(stl_filename)

    def test_cart3d_to_tecplot(self):
        """convert to tecplot"""
        log = get_logger(level='warning', encoding='utf-8')
        cart3d_filename = os.path.join(MODEL_PATH, 'threePlugs.bin.tri')
        tecplot_filename = os.path.join(MODEL_PATH, 'threePlugs.plt')
        #cart3d_to_tecplot(cart3d_filename, tecplot_filename, log=log)

        argv = ['format_converter', 'cart3d', cart3d_filename,
                'tecplot', tecplot_filename, '--scale', '2.0']
        cmd_line_format_converter(argv=argv, quiet=True)
        os.remove(tecplot_filename)

    def test_cart3d_to_nastran_01(self):
        """convert to nastran small field"""
        log = get_logger(level='warning', encoding='utf-8')
        cart3d_filename = os.path.join(MODEL_PATH, 'threePlugs.bin.tri')
        bdf_filename = os.path.join(MODEL_PATH, 'threePlugs.bdf')
        #cart3d_to_nastran_filename(cart3d_filename, bdf_filename, log=log)

        argv = ['format_converter', 'cart3d', cart3d_filename,
                'nastran', bdf_filename, '-b', '--scale', '2.0']
        cmd_line_format_converter(argv=argv, quiet=True)
        os.remove(bdf_filename)

    def test_cart3d_to_nastran_02(self):
        """convert to nastran large field"""
        log = get_logger(level='warning', encoding='utf-8')
        cart3d_filename = os.path.join(MODEL_PATH, 'threePlugs.bin.tri')
        bdf_filename = os.path.join(MODEL_PATH, 'threePlugs2.bdf')
        model = cart3d_to_nastran_model(cart3d_filename, log=log)
        model.write_bdf(bdf_filename, size=16)
        self.assertAlmostEqual(model.nodes[1].xyz[0], 1.51971436,
                               msg='if this is 0.0, then the assign_type method had the float32 check removed')
        os.remove(bdf_filename)
        #model.write_bdf(out_filename=None, encoding=None, size=8,
                       #is_double=False,
                       #interspersed=True,
                       #enddata=None)

    #def test_cart3d_input_cntl(self):
        #"""tests the input.cntl reading"""
        #from pyNastran.converters.cart3d.input_cntl_reader import read_input_cntl
        #input_cntl_filename = os.path.join(MODEL_PATH, '')
        #read_input_cntl(input_cntl_filename, log=None, debug=False)

    def test_cart3d_input_c3d(self):
        """tests the input.c3d reading"""
        log = get_logger(level='warning', encoding='utf-8')
        input_c3d_filename = os.path.join(MODEL_PATH, 'input.c3d')
        nodes, elements = read_input_c3d(input_c3d_filename, log=log, debug=False, stack=True)
        assert isinstance(nodes, np.ndarray), nodes
        assert isinstance(elements, np.ndarray), elements

    def test_cart3d_cntl(self):
        input_cntl_filename = os.path.join(MODEL_PATH, 'business_jet', 'input.cntl')
        cntl = read_input_cntl(input_cntl_filename, log=None, debug=None)
        bcs = cntl.get_boundary_conditions()
        mach, alpha, beta, gamma = cntl.get_flow_conditions()
        post = cntl.get_post_processing()
        assert np.allclose(mach, 0.84), mach
        assert np.allclose(alpha, 2.81), alpha
        assert np.allclose(beta, 0.), beta
        assert gamma is None, gamma
        #assert np.allclose(gamma, 1.4), gamma
        assert post is None, post
        assert bcs == (0, 0, 0, 0, 0, 0, {}), bcs


    def test_bjet(self):
        log = get_logger(level='warning', encoding='utf-8')
        geometry_filename = os.path.join(MODEL_PATH, 'business_jet', 'bJet.a.tri')

        model = read_cart3d(geometry_filename, log=log)
        #print(model.points.min(axis=0))
        #print(model.points.max(axis=0))
        model.make_half_model(axis='z', remap_nodes=True)  # keep +z
        #print(model.points.min(axis=0))
        #print(model.points.max(axis=0))
        geometry_filename_half = os.path.join(MODEL_PATH, 'half_bjet.a.tri')
        model.write_cart3d(geometry_filename_half)
        read_cart3d(geometry_filename_half, log=log)

class TestCart3dRotateRegions(unittest.TestCase):
    def test_rotate_regions_split(self):
        """region 2 shares nodes 1/2 with region 1 -> split; loads copied"""
        model = _two_region_model()
        points0 = model.points.copy()
        model.rotate_regions([([2], [1., 0., 0.], np.eye(3), 'y', 90.)])

        # new nodes 6/7 are copies of 1/2 (on the hinge, so they don't move)
        assert model.points.shape == (8, 3), model.points.shape
        assert np.array_equal(model.points[:4], points0[:4])  # region 1 doesn't move
        assert np.allclose(model.points[6:], points0[[1, 2]]), model.points
        # +90 about y: x -> -z
        assert np.allclose(model.points[4:6], [[1., 0., -1.], [1., 1., -1.]]), model.points
        assert np.array_equal(model.elements, [[0, 1, 2], [0, 2, 3],
                                               [6, 4, 5], [6, 5, 7]]), model.elements
        assert np.array_equal(model.loads['Cp'], [0., 1., 2., 3., 4., 5., 1., 2.])
        # rigid: areas don't change
        assert np.allclose(model.get_area(), 0.5), model.get_area()

    def test_rotate_regions_local_axis(self):
        """local x = basic y, so 'x' in the hinge frame == 'y' in basic"""
        hinge_axes = np.array([[0., 1., 0.], [-1., 0., 0.], [0., 0., 1.]])
        model1 = _two_region_model()
        model1.rotate_regions([([2], [1., 0., 0.], np.eye(3), 'y', 30.)])
        model2 = _two_region_model()
        model2.rotate_regions([([2], [1., 0., 0.], hinge_axes, 'x', 30.)])
        assert np.allclose(model1.points, model2.points), model2.points

    def test_rotate_regions_errors(self):
        model = _two_region_model()
        points0 = model.points.copy()
        eye = np.eye(3)
        zero = np.zeros(3)
        with self.assertRaises(KeyError):
            model.rotate_regions([([99], zero, eye, 'x', 10.)])
        with self.assertRaises(ValueError):  # left-handed
            model.rotate_regions([([2], zero, np.diag([1., 1., -1.]), 'x', 10.)])
        with self.assertRaises(ValueError):  # region 2 is in 2 sets
            model.rotate_regions([([2], zero, eye, 'x', 10.), ([1, 2], zero, eye, 'y', 5.)])
        with self.assertRaises(ValueError):  # bad axis
            model.rotate_regions([([2], zero, eye, 'w', 10.)])
        # failed calls don't change the model
        assert np.array_equal(model.points, points0)

    def test_rotate_regions_velocity_warning(self):
        messages = []
        log = SimpleLogger(level='debug',
                           log_func=lambda typ, filename, lineno, msg: messages.append((typ, msg)))
        rotation = ([2], np.zeros(3), np.eye(3), 'z', 5.)
        model = _two_region_model(log=log)
        model.rotate_regions([rotation])
        assert not any(typ == 'WARNING' for typ, msg in messages), messages

        model = _two_region_model(log=log)
        model.loads['rhoU'] = np.zeros(6)
        model.rotate_regions([rotation])
        warnings_ = [msg for typ, msg in messages if typ == 'WARNING']
        assert len(warnings_) == 1 and 'rhoU' in warnings_[0], messages
        assert len(model.loads['rhoU']) == 8

    def test_rotate_regions_far_origin(self):
        """region 2 is 1.4 across; an origin ~1000 away warns, the hinge edge doesn't"""
        messages = []
        log = SimpleLogger(level='debug',
                           log_func=lambda typ, filename, lineno, msg: messages.append((typ, msg)))
        model = _two_region_model(log=log)
        model.rotate_regions([([2], [1., 0., 0.], np.eye(3), 'y', 10.)])  # on the hinge
        assert not any(typ == 'WARNING' for typ, msg in messages), messages

        model = _two_region_model(log=log)
        model.rotate_regions([([2], [916., 401., 0.], np.eye(3), 'z', 30.)])
        warnings_ = [msg for typ, msg in messages if typ == 'WARNING']
        assert len(warnings_) == 1 and 'far away' in warnings_[0], messages

    def test_business_jet_rotate(self):
        """verifies the origin isn't shifted by a rotation"""
        region = 3  # horizontal tail
        origin = [27.3, 10.0, 0.]  # ~midpoint of region 3 (bJet is ~68 long; y is up)
        # origin = [22.8, 10., 0.]

        cart3d_filename = MODEL_PATH / 'business_jet' / 'bJet.a.tri'
        cart3d_filename_out = MODEL_PATH / 'business_jet' / 'bJet2.tri'
        log = SimpleLogger(level='debug')
        model = read_cart3d(cart3d_filename, log=log)
        rotation = ([region], origin, np.eye(3), 'z', 30.)
        nodes = np.unique(model.elements[model.regions == region])
        center0 = model.points[nodes].mean(axis=0)
        model.rotate_regions([rotation])
        # rotating about its own midpoint keeps the tail in place
        nodes = np.unique(model.elements[model.regions == region])
        center = model.points[nodes].mean(axis=0)
        assert np.allclose(center, center0, atol=0.5), (center0, center)
        model.write_cart3d(cart3d_filename_out)
        cart3d_filename_out.unlink()


def _two_region_model(log=None) -> Cart3D:
    """region 1: square in x=[0, 1]; region 2: square in x=[1, 2]; share nodes 1/2"""
    if log is None:
        log = SimpleLogger(level='warning')
    model = Cart3D(log=log)
    model.points = np.array([
        [0., 0., 0.], [1., 0., 0.], [1., 1., 0.],
        [0., 1., 0.], [2., 0., 0.], [2., 1., 0.],
    ])
    model.elements = np.array([[0, 1, 2], [0, 2, 3], [1, 4, 5], [1, 5, 2]], dtype='int32')
    model.regions = np.array([1, 1, 2, 2], dtype='int32')
    model.loads = {'Cp': np.arange(6, dtype='float64')}
    return model


def check_array(points, points2):
    nnodes = points.shape[0]
    msg = ''
    nfailed = 0
    if not np.array_equal(points, points2):  # pragma: no cover
        for nid in range(nnodes):
            p1 = points[nid]
            p2 = points2[nid]
            abs_sum_delta = sum(abs(p1-p2))
            if not np.allclose(abs_sum_delta, 0.0, atol=1e-6):
                msg += 'n=%s p1=%s p2=%s diff=%s\nsum(abs(p1-p2))=%s\n' % (
                    nid, str(p1), str(p2), str(p1-p2), abs_sum_delta)
                nfailed += 1
                if nfailed == 10:
                    break
    if msg:  # pragma: no cover
        #print(msg)
        raise RuntimeError(msg)

if __name__ == '__main__':  # pragma: no cover
    import time
    time0 = time.time()
    unittest.main()
    print("dt = %s" % (time.time() - time0))
