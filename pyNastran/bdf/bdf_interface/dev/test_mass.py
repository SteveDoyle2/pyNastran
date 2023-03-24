import unittest

import numpy as np
from pyNastran.bdf.bdf import BDF
from pyNastran.bdf.bdf_interface.dev.mass import (
    make_gpwg, make_mass_matrix, make_reduced_mass_matrix)


class TestMassGeneration(unittest.TestCase):
    def test_make_gpwg_1(self):
        """
        per basic_dynamics.pdf

        O U T P U T F R O M G R I D P O I N T W E I G H T G E N E R A T O R
        REFERENCE POINT = 0
        M O
        *  1.300000E+01  0.000000E+00  0.000000E+00  0.000000E+00  5.000000E+00 -3.000000E+00 *
        *  0.000000E+00  1.300000E+01  0.000000E+00 -5.000000E+00  0.000000E+00 7.000000E+00 *
        *  0.000000E+00  0.000000E+00  1.300000E+01  3.000000E+00 -7.000000E+00 0.000000E+00 *
        *  0.000000E+00 -5.000000E+00  3.000000E+00  8.000000E+00 -1.500000E+00 -2.500000E+00 *
        *  5.000000E+00  0.000000E+00 -7.000000E+00 -1.500000E+00  1.000000E+01 0.000000E+00 *
        * -3.000000E+00  7.000000E+00  0.000000E+00 -2.500000E+00  0.000000E+00 8.000000E+00 *
        S
        * 1.000000E+00 0.000000E+00 0.000000E+00 *
        * 0.000000E+00 1.000000E+00 0.000000E+00 *
        * 0.000000E+00 0.000000E+00 1.000000E+00 *
        DIRECTION
        MASS AXIS SYSTEM (S) MASS X-C.G. Y-C.G. Z-C.G.
        X 1.300000E+01 0.000000E+00 2.307692E-01 3.846154E-01
        Y 1.300000E+01 5.384616E-01 0.000000E+00 3.846154E-01
        Z 1.300000E+01 5.384616E-01 2.307692E-01 0.000000E+00
        I(S)
        *  5.384615E+00 -1.153847E-01 -1.923079E-01 *
        * -1.153847E-01  4.307692E+00 -1.153846E+00 *
        * -1.923079E-01 -1.153846E+00  3.538461E+00 *
        I(Q)
        * 5.503882E+00 *
        * 5.023013E+00 *
        * 2.703873E+00 *
        Q
        * 8.702303E-01  4.915230E-01  3.323378E-02 *
        * 3.829170E-01 -7.173043E-01  5.821075E-01 *
        * 3.099580E-01 -4.938418E-01 -8.124324E-01 *
        """
        #make_gpwg(Mgg, reference_point, xyz_cid0, grid_cps, coords, log)

        model = BDF()

        eid = nid = 1
        mass = 2.
        model.add_grid(nid, [0., 0., 0.], cp=0, cd=0, ps='', seid=0, comment='')
        model.add_conm2(eid, nid, mass, cid=0, X=None, I=None, comment='')

        eid = nid = 2
        mass = 3.
        model.add_grid(nid, [1., 0., 0.], cp=0, cd=0, ps='', seid=0, comment='')
        model.add_conm2(eid, nid, mass, cid=0, X=None, I=None, comment='')

        eid = nid = 3
        mass = 3.
        model.add_grid(nid, [0.5, 1., 0.], cp=0, cd=0, ps='', seid=0, comment='')
        model.add_conm2(eid, nid, mass, cid=0, X=None, I=None, comment='')

        eid = nid = 4
        mass = 5.
        model.add_grid(nid, [0.5, 0., 1.], cp=0, cd=0, ps='', seid=0, comment='')
        model.add_conm2(eid, nid, mass, cid=0, X=None, I=None, comment='')

        reference_point = 0
        #print(model.elements)
        D, Mo, S, mass, cg, IS, IQ, Q = make_reduced_mass_matrix(model, reference_point, fdtype='float64', idtype='int32')
        Mo_expected = [
            [13.,   0.,  0.,  0. ,  5., -3. ],
            [ 0.,  13.,  0., -5. ,  0.,  7. ],
            [ 0.,   0., 13.,  3. , -7.,  0. ],
            [ 0.,  -5.,  3.,  8. , -1.5,-2.5],
            [ 5.,   0., -7., -1.5, 10.,  0. ],
            [-3.,   7.,  0., -2.5,  0.,  8. ],
        ]
        S_expected = np.eye(3)
        mass_expected = [13., 13., 13.]
        cg_expected = [
            [0.,         0.23076923, 0.3846154 ],
            [0.53846157, 0.,         0.3846154 ],
            [0.53846157, 0.23076923, 0.        ], ]
        IS_expected = [
            [ 5.3846154,  -0.11538471, -0.19230787],
            [-0.11538471,  4.307692,   -0.19230787],
            [-0.19230787, -1.1538461,   3.5384612 ], ]
        #IQ_expected = [
            #[ 0.99436104, -0.10022114,  0.07963255],
            #[-0.09721395, -0.19528615, -0.66057014],
            #[-0.04237364, -0.97561216,  0.7465291 ], ]
        IQ_expected = np.array([
            [ 5.404091  , -0.12931095,  0.507303  ],
            [-0.21254948,  3.2877433 , -2.7565107 ],
            [ 0.60399896, -1.996658  ,  4.538934  ]])
        #II_expected = np.array([
            #[ 5.3846154 , -0.11538471, -0.19230787],
            #[-0.11538471,  4.307692  , -0.19230787],
            #[-0.19230787, -1.1538461 ,  3.5384612 ]])
        Q_expected = np.array([
            [ 0.99436104, -0.10022114,  0.07963255],
            [-0.09721395, -0.19528615, -0.66057014],
            [-0.04237364, -0.97561216,  0.7465291 ]])
        assert np.allclose(D, D_expected)
        assert np.allclose(Mo, Mo_expected)
        assert np.allclose(S, S_expected)
        assert np.allclose(mass, mass_expected)
        assert np.allclose(cg, cg_expected)
        assert np.allclose(IS, IS_expected)
        assert np.allclose(IQ, IQ_expected)
        assert np.allclose(Q, Q_expected)

    def test_make_gpwg_2(self):
        """directional dependence"""
        model = BDF()

        masses = [2., 3., 5.]
        components = [1, 2, 3]
        xyz1 = [0., 0., 0.]
        xyz2 = [1., 0., 0.]
        xyz3 = [0.5, 1., 0.]
        xyz4 = [0.5, 0., 1.]

        cid = 1
        origin = [0., 0., 0.]
        zaxis = [0., 0., 1.]
        theta = np.radians(45.)
        xaxis = [np.cos(theta), np.sin(theta), 0.]
        xzplane = xaxis
        model.add_cord2r(cid, origin, zaxis, xzplane, rid=0, setup=True, comment='')

        cid = 3
        origin = [0., 0., 0.]
        zaxis = [0., 0., 1.]
        theta = np.radians(60.)
        xaxis = [np.cos(theta), np.sin(theta), 0.]
        xzplane = xaxis
        model.add_cord2r(cid, origin, zaxis, xzplane, rid=0, setup=True, comment='')

        nid = 1
        nids = [nid, nid]
        model.add_grid(nid, xyz1, cp=0, cd=1, ps='', seid=0, comment='')
        eid = 1
        for massi, ci in zip(masses, components):
            model.add_cmass2(eid, massi, nids, ci, ci, comment='')
            eid += 1
        eid += 7
        assert eid == 11, eid

        nid = 2
        nids = [nid, nid]
        model.add_grid(nid, xyz2, cp=0, cd=0, ps='', seid=0, comment='')
        for massi, ci in zip(masses, components):
            model.add_cmass2(eid, massi, nids, ci, ci, comment='')
            eid += 1
        eid += 7

        nid = 3
        nids = [nid, nid]
        model.add_grid(nid, xyz3, cp=0, cd=3, ps='', seid=0, comment='')
        for massi, ci in zip(masses, components):
            model.add_cmass2(eid, massi, nids, ci, ci, comment='')
            eid += 1
        eid += 7

        nid = 4
        nids = [nid, nid]
        model.add_grid(nid, xyz4, cp=0, cd=0, ps='', seid=0, comment='')
        for massi, ci in zip(masses, components):
            model.add_cmass2(eid, massi, nids, ci, ci, comment='')
            eid += 1
        eid += 7

        reference_point = 0
        #print(model.elements)
        Mgg, D, Mo, S, mass, cg, IS, IQ, Q = make_reduced_mass_matrix(
            model, reference_point, fdtype='float64', idtype='int32')
        x = 1
        Mo_expected = np.array([
            [ 9.24999974, -0.93301261,   0.,  0. , 2.  , -2.96650635],
            [-0.93301263, 10.74999982,   0., -3. , 0.  ,  6.05801277],
            [ 0.        ,  0.        ,  20.,  5. , -10.,  0.        ],
            [ 0.        , -3.        ,   5.,  8. , -2.5, -1.5       ],
            [ 2.        ,  0.        , -10., -2.5, 9.5 ,  0.        ],
            [-2.96650641,  6.05801277,   0., -1.5, 0.  ,  7.495513  ]])
        S_expected = np.array([ #  flipped n1 and n2...
            [-0.90180983, -0.43213312,  0.        ],
            [-0.43213312,  0.90180983,  0.        ],
            [-0.        ,  0.        ,  1.        ]])
        #assert np.allclose(D, D_expected)
        assert np.allclose(Mo, Mo_expected)

if __name__ == '__main__':
    unittest.main()
