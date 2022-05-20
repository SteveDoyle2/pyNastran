import os
import unittest
import numpy as np

import pyNastran
from pyNastran.bdf.bdf import BDF
from pyNastran.bdf.mesh_utils.refine import (
    refine_model, _quad_nids_to_node_ids, _hexa_nids_to_node_ids,
    _insert_quad_nodes)
pkg_path = pyNastran.__path__[0]
#model_path = os.path.join(pkg_path, '..', 'models')
#bwb_path = os.path.join(model_path, 'bwb')
DIRNAME = os.path.dirname(__file__)



class TestRefine(unittest.TestCase):

    def test_insert_tri_nodes(self):
        #[[  10196 1206947   10184]
            #[      0 1206949 1206937]
            #[      0       0   10195]]
        nodes = np.array([
            [1, 4, 2],
            [0, 6, 5],
            [0, 0, 3]], dtype='int32')
        n1, n4, n2 = nodes[0, :]
        n6, n5 = nodes[1, 1:]
        n3 = nodes[2, 2]
        unids = np.unique([n1, n2, n3, n4, n5, n6])
        assert len(unids) == 6

    def test_insert_quad_nodes(self):
        n1 = 1
        n2 = 2
        n3 = 3
        n4 = 4
        n5 = 5
        n6 = 6
        n7 = 7
        n8 = 8

        nodes = {
            1: np.array([0., 0., 0.]),
            2: np.array([1., 0., 0.]),
            3: np.array([1., 1., 0.]),
            4: np.array([0., 1., 0.]),
        }
        nid0 = max(nodes) + 1
        xyz1 = nodes[n1]
        xyz2 = nodes[n2]
        xyz3 = nodes[n3]
        xyz4 = nodes[n4]

        nnodes_to_add_with_ends = 3
        nids_array = np.zeros((nnodes_to_add_with_ends, nnodes_to_add_with_ends), dtype='int32')
        nids_array[0, 0] = n1
        nids_array[0, -1] = n2
        nids_array[-1, -1] = n3
        nids_array[-1, 0] = n4
        forward_edges = [(n1, n2), (n2, n3), (n3, n4), (n4, n1)]
        edges = forward_edges

        edges_to_center = {
            (n1, n2) : [n1, n5, n2],
            (n2, n3) : [n2, n6, n3],
            (n3, n4) : [n3, n7, n4],
            (n4, n1) : [n4, n8, n1],
        }
        n9 = 9
        faces_to_center = {
            (n1, n2, n3, n4): [n9],
        }
        face = (n1, n2, n3, n4)
        nid0 = _insert_quad_nodes(
            nodes, nids_array, nid0,
            edges, forward_edges,
            edges_to_center, faces_to_center,
            nnodes_to_add_with_ends,
            face,
            xyz1, xyz2, xyz3, xyz4,
        )

    def test_quad_nids_to_node_ids(self):
        nids = np.array([
            [1, 5, 2],
            [8, 9, 6],
            [4, 7, 3],
        ])
        n1, n2, n3, n4 = _quad_nids_to_node_ids(nids)
        assert np.array_equal(n1, [1, 5, 8, 9]), n1
        assert np.array_equal(n2, [5, 2, 9, 6]), n2
        assert np.array_equal(n3, [9, 6, 7, 3]), n3
        assert np.array_equal(n4, [8, 9, 4, 7]), n4

    def test_hexa_nids_to_node_ids(self):
        nids = np.array([
            [
                [ 1, 17,  5],
                [ 9, 21, 13],
                [ 2, 18,  6]],

            [
                [12, 22, 16],
                [23, 24, 25],
                [10, 26, 14]],

            [
                [ 4, 20,  8],
                [11, 27, 15],
                [ 3, 19,  7]]]
                        )
        nodes = _hexa_nids_to_node_ids(nids)
        n1 = nodes[:, 0]
        n2 = nodes[:, 1]
        n3 = nodes[:, 2]
        n4 = nodes[:, 3]
        n5 = nodes[:, 4]
        n6 = nodes[:, 5]
        n7 = nodes[:, 6]
        n8 = nodes[:, 7]
        un1 = np.unique(n1)
        un2 = np.unique(n2)
        un3 = np.unique(n3)
        un4 = np.unique(n4)
        un5 = np.unique(n5)
        un6 = np.unique(n6)
        un7 = np.unique(n7)
        un8 = np.unique(n8)
        #assert np.array_equal(un1, np.unique([1, 9, 23, 12,
                                              #17, 22, 24, 21])), un1
        #assert np.array_equal(un2, np.unique([9, 2, 10, 23,
                                              #22, 18, 27, 24])), un2
        #assert np.array_equal(un3, np.unique([23, 10, 3, 11,
                                              #24, 27, 19, 26])), un3
        #assert np.array_equal(un4, np.unique([12, 23, 11, 4, 21, 24, 26, 20])), un4
        #assert np.array_equal(n1, [1, 5, 8, 9]), n1
        #assert np.array_equal(n2, [5, 2, 9, 6]), n2
        #assert np.array_equal(n3, [9, 6, 7, 3]), n3
        #assert np.array_equal(n4, [8, 9, 4, 7]), n4
        x = 1

    def test_tri(self):
        bdf_filename_out = os.path.join(DIRNAME, 'tri.bdf')

        model = BDF()
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [1., 0., 0.])
        model.add_grid(3, [1., 1., 0.])
        model.add_grid(4, [0., 1., 0.])
        model.add_ctria3(1, 1, [1, 2, 3])
        model.add_ctria3(2, 1, [1, 3, 4])
        model.add_pshell(1, mid1=1, t=0.1)
        model.add_mat1(1, 3.0e7, None, 0.3)
        ntris = len(model.elements)

        model = refine_model(model, refinement_ratio=2)
        model.write_bdf(bdf_filename_out)
        model.validate()
        assert len(model.elements) == ntris * 4
        model.cross_reference()
        x = 1

    def test_quad_bar(self):
        bdf_filename_out = os.path.join(DIRNAME, 'quad_bar.bdf')

        model = BDF()
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [1., 0., 0.])
        model.add_grid(3, [1., 1., 0.])
        model.add_grid(4, [0., 1., 0.])
        model.add_cquad4(1, 1, [1, 2, 3, 4])
        x = [0., 0., 1.]
        g0 = None
        model.add_cbar(2, 2, [1, 2], x, g0,
                       offt='GGG', pa=0, pb=0, wa=None, wb=None, comment='', validate=False)
        model.add_pbarl(2, 1, 'BAR', [0., 1.])
        model.add_pshell(1, mid1=1, t=0.1)
        model.add_mat1(1, 3.0e7, None, 0.3)

        nquads = 1
        nbars = 1

        model = refine_model(model, refinement_ratio=2)
        model.write_bdf(bdf_filename_out)
        model.validate()
        assert len(model.elements) == nbars * 2 + nquads * 4
        assert len(model.nodes) == nquads * 9, len(model.nodes)
        model.cross_reference()
        x = 1

    def test_quad(self):
        bdf_filename_out = os.path.join(DIRNAME, 'quad.bdf')

        model = BDF()
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [1., 0., 0.])
        model.add_grid(3, [1., 1., 0.])
        model.add_grid(4, [0., 1., 0.])
        model.add_cquad4(1, 1, [1, 2, 3, 4])
        model.add_pshell(1, mid1=1, t=0.1)
        model.add_mat1(1, 3.0e7, None, 0.3)
        nquads = len(model.elements)

        model = refine_model(model, refinement_ratio=2)
        model.write_bdf(bdf_filename_out)
        model.validate()
        assert len(model.elements) == nquads * 4
        model.cross_reference()
        x = 1

    def test_hexa(self):
        bdf_filename_out = os.path.join(DIRNAME, 'hexa.bdf')

        model = BDF()
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [1., 0., 0.])
        model.add_grid(3, [1., 1., 0.])
        model.add_grid(4, [0., 1., 0.])

        model.add_grid(5, [0., 0., 1.])
        model.add_grid(6, [1., 0., 1.])
        model.add_grid(7, [1., 1., 1.])
        model.add_grid(8, [0., 1., 1.])

        model.add_chexa(1, 1, [1, 2, 3, 4, 5, 6, 7, 8])
        model.add_psolid(1, 1)
        model.add_mat1(1, 3.0e7, None, 0.3)
        nhexas = len(model.elements)

        model = refine_model(model, refinement_ratio=2)
        model.write_bdf(bdf_filename_out)
        model.validate()
        assert len(model.elements) == nhexas * 8
        model.cross_reference()
        x = 1

    def test_hexa_quad(self):
        bdf_filename_out = os.path.join(DIRNAME, 'hexa_quad.bdf')

        model = BDF()
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [1., 0., 0.])
        model.add_grid(3, [1., 1., 0.])
        model.add_grid(4, [0., 1., 0.])

        model.add_grid(5, [0., 0., 1.])
        model.add_grid(6, [1., 0., 1.])
        model.add_grid(7, [1., 1., 1.])
        model.add_grid(8, [0., 1., 1.])
        model.add_cquad4(2, 2, [4, 3, 2, 1])
        model.add_pshell(2, mid1=1, t=0.1)

        model.add_chexa(1, 1, [1, 2, 3, 4, 5, 6, 7, 8])
        model.add_psolid(1, 1)
        model.add_mat1(1, 3.0e7, None, 0.3)
        nquads = 1
        nhexas = 1
        nelements = 4 * nquads + 8 * nhexas

        nlayers = 3
        nnodes_nlayers = 9
        nnodes = nhexas * nlayers * nnodes_nlayers


        model = refine_model(model, refinement_ratio=2)
        model.write_bdf(bdf_filename_out)
        model.validate()
        assert len(model.elements) == nelements
        assert len(model.nodes) == nnodes
        model.cross_reference()

        x = 1

    def _test_tri_penta(self):
        bdf_filename_out = os.path.join(DIRNAME, 'penta_tri.bdf')

        model = BDF()
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [1., 0., 0.])
        model.add_grid(3, [1., 1., 0.])

        model.add_grid(4, [0., 0., 1.])
        model.add_grid(5, [1., 0., 1.])
        model.add_grid(6, [1., 1., 1.])
        model.add_ctria3(1, 1, [1, 2, 3])
        model.add_pshell(1, mid1=1, t=0.1)

        model.add_cpenta(2, 1, [1, 2, 3, 4, 5, 6])
        model.add_psolid(2, 1)
        model.add_mat1(1, 3.0e7, None, 0.3)
        ntri = 1
        npenta = 1
        nelements = 6 * npenta + ntri

        nlayers = 3
        nnodes_nlayers = 6
        nnodes = npenta * nlayers * nnodes_nlayers


        model = refine_model(model, refinement_ratio=2)
        model.write_bdf(bdf_filename_out)
        model.validate()
        assert len(model.elements) == npenta * 8
        model.cross_reference()

        assert len(model.elements) == nelements
        assert len(model.nodes) == nnodes
        x = 1

    def test_hexa2(self):
        bdf_filename_out = os.path.join(DIRNAME, 'hexa2.bdf')

        model = BDF()
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [1., 0., 0.])
        model.add_grid(3, [1., 1., 0.])
        model.add_grid(4, [0., 1., 0.])

        model.add_grid(5, [0., 0., 1.])
        model.add_grid(6, [1., 0., 1.])
        model.add_grid(7, [1., 1., 1.])
        model.add_grid(8, [0., 1., 1.])

        model.add_grid(15, [0., 0., 2.])
        model.add_grid(16, [1., 0., 2.])
        model.add_grid(17, [1., 1., 2.])
        model.add_grid(18, [0., 1., 2.])

        model.add_chexa(1, 1, [1, 2, 3, 4,
                               5, 6, 7, 8])
        model.add_chexa(2, 1, [5, 6, 7, 8,
                               15, 16, 17, 18])
        model.add_psolid(1, 1)
        model.add_mat1(1, 3.0e7, None, 0.3)
        nhexas = len(model.elements)

        model = refine_model(model, refinement_ratio=2)
        model.write_bdf(bdf_filename_out)
        model.validate()
        ncommon_faces = 1
        nnodes_layer = 9
        nlayers = 3
        nnodes_large_hexa = nnodes_layer * nlayers
        assert len(model.elements) == nhexas * 8
        assert len(model.nodes) == nhexas * nnodes_large_hexa - ncommon_faces * nnodes_layer
        model.cross_reference()
        x = 1

    def test_quad2(self):
        bdf_filename_out = os.path.join(DIRNAME, 'quad2.bdf')

        model = BDF()
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [1., 0., 0.])
        model.add_grid(3, [1., 1., 0.])
        model.add_grid(4, [0., 1., 0.])
        model.add_grid(5, [1., 2., 0.])
        model.add_grid(6, [0., 2., 0.])
        model.add_cquad4(1, 1, [1, 2, 3, 4])
        model.add_cquad4(2, 1, [4, 3, 5, 6])
        model.add_pshell(1, mid1=1, t=0.1)
        model.add_mat1(1, 3.0e7, None, 0.3)
        nquads = len(model.elements)

        model = refine_model(model, refinement_ratio=2)
        model.write_bdf(bdf_filename_out)
        model.validate()
        assert len(model.elements) == nquads * 4
        model.cross_reference()
        x = 1

    def test_split_quad2(self):
        nids_array = np.array([
            [1, 7, 2],
            [10, 14, 8],
            [4, 9, 3], ], dtype='int32')
        n1, n2, n3, n4 = _quad_nids_to_node_ids(nids_array)
        x = 1


    def _test_refine_bwb(self):
        model_path = os.path.join(pkg_path, '..', 'models')
        bwb_path = os.path.join(model_path, 'bwb')
        bdf_filename = os.path.join(bwb_path, 'bwb_saero.bdf')
        bdf_filename_out = os.path.join(bwb_path, 'bwb_saero_fine.bdf')

        #bdf_filename = os.path.join(model_path, 'plate', 'plate.bdf')
        model = refine_model(bdf_filename, refinement_ratio=2)
        model.write_bdf(bdf_filename_out)
        model.validate()
        model.cross_reference()
        x = 1

if __name__ == '__main__':
    unittest.main()
