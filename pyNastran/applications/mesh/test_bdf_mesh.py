from pyNastran.applications.mesh.bdf_mesh import MeshTools
import unittest

class TestMeshTools(unittest.TestCase):
    def test_get_free_nodes_01(self):
        """
        1      2      3      4
        *--11--*--22--*--33--*
               |
               44
               |
               *5
        """
        model = MeshTools()
        model.add_card(['GRID', 1, None, 0., 0., 0.], 'GRID')
        model.add_card(['GRID', 2, None, 1., 0., 0.], 'GRID')
        model.add_card(['GRID', 3, None, 2., 0., 0.], 'GRID')
        model.add_card(['GRID', 4, None, 3., 0., 0.], 'GRID')

        A = 1.0
        mid = 42
        E = 1e7
        nu = 0.3
        model.add_card(['CONROD', 11, 1, 2, mid, A], 'CONROD')
        model.add_card(['CONROD', 22, 2, 3, mid, A], 'CONROD')
        model.add_card(['CONROD', 33, 3, 4, mid, A], 'CONROD')
        model.add_card(['MAT1', mid, E, None, nu], 'MAT1')
        #model.cross_reference()

        eids = [11]
        """
        1      2      3      4
        *--11--*--22--*--33--*
        """
        used_free_eids, used_free_nids, free_eids, free_nids = model.get_free_nodes(eids)
        self.assertEqual(set([11, 22]), used_free_eids)
        self.assertEqual(set([22]), free_eids)
        self.assertEqual(set([1, 2, 3]), used_free_nids)
        self.assertEqual(set([3]), free_nids)

        model.add_card(['GRID', 5, None, 1., 1., 0.], 'GRID')
        model.add_card(['CONROD', 44, 2, 5, mid, A], 'CONROD')

        eids = [11]
        """
        1      2      3      4
        *--11--*--22--*--33--*
               |
               44
               |
               *5
        """
        used_free_eids, used_free_nids, free_eids, free_nids = model.get_free_nodes(eids)
        self.assertEqual(set([11, 22, 44]), used_free_eids)
        self.assertEqual(set([22, 44]), free_eids)
        self.assertEqual(set([1, 2, 3, 5]), used_free_nids)
        self.assertEqual(set([3, 5]), free_nids)

        eids = [22]
        """
        1      2      3      4
        *--11--*--22--*--33--*
               |
               44
               |
               *5
        """
        used_free_eids, used_free_nids, free_eids, free_nids = model.get_free_nodes(eids)
        self.assertEqual(set([11, 22, 33, 44]), used_free_eids)
        self.assertEqual(set([11, 33, 44]), free_eids)
        self.assertEqual(set([1, 2, 3, 4, 5]), used_free_nids)
        self.assertEqual(set([1, 4, 5]), free_nids)
        #model.read_bdf(bdf_filename)

    def test_get_free_nodes_02(self):
        """
        1      2      3      4
        *--11--*--22--*--33--*
               |
               44
               |
               *5
        """
        model = MeshTools()
        model.add_card(['GRID', 1, None, 0., 0., 0.], 'GRID')
        model.add_card(['GRID', 2, None, 1., 0., 0.], 'GRID')
        model.add_card(['GRID', 3, None, 2., 0., 0.], 'GRID')
        model.add_card(['GRID', 4, None, 3., 0., 0.], 'GRID')

        A = 1.0
        mid = 42
        E = 1e7
        nu = 0.3
        model.add_card(['CONROD', 11, 1, 2, mid, A], 'CONROD')
        model.add_card(['MAT1', mid, E, None, nu], 'MAT1')
        model.add_card(['MPC', 22, 2, 123456, -1., 3, 123456, 1., None,
                             None, 4, 123456, 2.,
                        ], 'MPC')
        eids = [11]
        used_free_eids, used_free_nids, free_eids, free_nids = model.get_free_nodes(eids)
        #print used_free_eids
        #print used_free_nids
        #print free_eids
        #print free_nids
        self.assertEqual(set([11, 22]), used_free_eids)
        self.assertEqual(set([22]), free_eids)
        self.assertEqual(set([1, 2, 3, 4]), used_free_nids)
        self.assertEqual(set([3, 4]), free_nids)

        eids = [22]
        used_free_eids, used_free_nids, free_eids, free_nids = model.get_free_nodes(eids)
        #print used_free_eids
        #print used_free_nids
        #print free_eids
        #print free_nids
        self.assertEqual(set([11, 22]), used_free_eids)
        self.assertEqual(set([11]), free_eids)
        self.assertEqual(set([1, 2, 3, 4]), used_free_nids)
        self.assertEqual(set([1]), free_nids)


if __name__ == '__main__':  # pragma: no cover
    unittest.main()
