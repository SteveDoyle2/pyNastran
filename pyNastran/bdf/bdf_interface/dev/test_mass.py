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
        I(S) - good
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
        Mgg, D, Mo, mass, cg, S, IS, II, IQ, Q = make_reduced_mass_matrix(
            model, reference_point, fdtype='float64', idtype='int32')

        Mo_expected = [
            [13.,   0.,  0.,  0. ,  5., -3. ],
            [ 0.,  13.,  0., -5. ,  0.,  7. ],
            [ 0.,   0., 13.,  3. , -7.,  0. ],
            [ 0.,  -5.,  3.,  8. , -1.5,-2.5],
            [ 5.,   0., -7., -1.5, 10.,  0. ],
            [-3.,   7.,  0., -2.5,  0.,  8. ],
        ]
        S_expected = np.eye(3)

        #MASS AXIS SYSTEM (S) MASS X-C.G. Y-C.G. Z-C.G.
        #X 1.300000E+01 0.000000E+00 2.307692E-01 3.846154E-01
        #Y 1.300000E+01 5.384616E-01 0.000000E+00 3.846154E-01
        #Z 1.300000E+01 5.384616E-01 2.307692E-01 0.000000E+00
        mass_expected = [13., 13., 13.]
        cg_expected = [
            [0.,         0.23076923, 0.3846154 ],
            [0.53846157, 0.,         0.3846154 ],
            [0.53846157, 0.23076923, 0.        ], ]

        #I(S) - good
        #*  5.384615E+00 -1.153847E-01 -1.923079E-01 *
        #* -1.153847E-01  4.307692E+00 -1.153846E+00 *
        #* -1.923079E-01 -1.153846E+00  3.538461E+00 *
        #
        # IS is II, but with flipped sign on diagonals
        #
        IS_expected = np.array([
            [ 5.3846154,  -0.11538471, -0.19230787],
            [-0.11538471,  4.307692,   -1.1538461 ],
            [-0.19230787, -1.1538461,   3.5384612 ], ])
        II_expected = np.array([
            [ 5.3846154,   0.11538471,  0.19230787],
            [ 0.11538471,  4.307692,    1.1538461 ],
            [ 0.19230787,  1.1538461,   3.5384612 ], ])

        #I(Q)
        #* 5.503882E+00 *
        #* 5.023013E+00 *
        #* 2.703873E+00 *
        IQ_expected = np.array([
            [ 5.5038824e+00,  3.6633622e-08,  8.1878504e-08],
            [ 8.8647887e-08,  5.0230141e+00, -3.5885233e-07],
            [ 5.3454769e-09,  1.8574326e-08,  2.7038724e+00]])

        # Q
        # * 8.702303E-01  4.915230E-01  3.323378E-02 *
        # * 3.829170E-01 -7.173043E-01  5.821075E-01 *
        # * 3.099580E-01 -4.938418E-01 -8.124324E-01 *
        Q_expected = np.array([
            [ 0.8702302 ,  0.4915236 , -0.03322647],
            [ 0.3829222 , -0.7173038 , -0.5821046 ],
            [ 0.30995163, -0.4938419 ,  0.81243473]])

        D_expected = np.array([
            [1., 0., 0., 0., 0. , 0. ],
            [0., 1., 0., 0., 0. , 0. ],
            [0., 0., 1., 0., 0. , 0. ],
            [0., 0., 0., 1., 0. , 0. ],
            [0., 0., 0., 0., 1. , 0. ],
            [0., 0., 0., 0., 0. , 1. ],
            [1., 0., 0., 0., 0. , 0. ],
            [0., 1., 0., 0., 0. , 1. ],
            [0., 0., 1., 0.,-1. , 0. ],
            [0., 0., 0., 1., 0. , 0. ],
            [0., 0., 0., 0., 1. , 0. ],
            [0., 0., 0., 0., 0. , 1. ],
            [1., 0., 0., 0., 0. ,-1. ],
            [0., 1., 0., 0., 0. , 0.5],
            [0., 0., 1., 1.,-0.5, 0. ],
            [0., 0., 0., 1., 0. , 0. ],
            [0., 0., 0., 0., 1. , 0. ],
            [0., 0., 0., 0., 0. , 1. ],
            [1., 0., 0., 0., 1. , 0. ],
            [0., 1., 0.,-1., 0. , 0.5],
            [0., 0., 1., 0.,-0.5, 0. ],
            [0., 0., 0., 1., 0. , 0. ],
            [0., 0., 0., 0., 1. , 0. ],
            [0., 0., 0., 0., 0. , 1. ],
        ])
        assert np.allclose(D, D_expected)
        assert np.allclose(Mo, Mo_expected)
        assert np.allclose(S, S_expected), f'S:\n{S}\nS_expected:\n{S_expected}'
        assert np.allclose(mass, mass_expected)
        assert np.allclose(cg, cg_expected)
        assert np.allclose(IS, IS_expected)
        assert np.allclose(II, II_expected)
        assert np.allclose(Q, Q_expected)
        assert np.allclose(IQ, IQ_expected)

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
        Mgg, D, Mo, mass, cg, S, IS, II, IQ, Q = make_reduced_mass_matrix(
            model, reference_point, fdtype='float64', idtype='int32')
        x = 1
        Mo_expected = np.array([
            [ 9.24999974, -0.93301261,   0.,  0. , 2.  , -2.96650635],
            [-0.93301263, 10.74999982,   0., -3. , 0.  ,  6.05801277],
            [ 0.        ,  0.        ,  20.,  5. , -10.,  0.        ],
            [ 0.        , -3.        ,   5.,  8. , -2.5, -1.5       ],
            [ 2.        ,  0.        , -10., -2.5, 9.5 ,  0.        ],
            [-2.96650641,  6.05801277,   0., -1.5, 0.  ,  7.495513  ]])
        D_expected = np.array([
            [ 0.70710677, 0.70710677, -0.,  0.        ,  0.        , 0.        ],
            [-0.70710677, 0.70710677,  0.,  0.        ,  0.        , 0.        ],
            [ 0.        , 0.        ,  1.,  0.        ,  0.        , 0.        ],
            [ 0.        , 0.        ,  0.,  0.70710677,  0.70710677, 0.        ],
            [ 0.        , 0.        ,  0., -0.70710677,  0.70710677, 0.        ],
            [ 0.        , 0.        ,  0.,  0.        ,  0.        , 1.        ],
            [ 1.        , 0.        ,  0.,  0.        ,  0.        , 0.        ],
            [ 0.        , 1.        ,  0.,  0.        ,  0.        , 1.        ],
            [ 0.        , 0.        ,  1.,  0.        , -1.        , 0.        ],
            [ 0.        , 0.        ,  0.,  1.        ,  0.        , 0.        ],
            [ 0.        , 0.        ,  0.,  0.        ,  1.        , 0.        ],
            [ 0.        , 0.        ,  0.,  0.        ,  0.        , 1.        ],
            [ 0.5       , 0.8660254 , -0.,  0.        ,  0.        , -0.0669873 ],
            [-0.8660254 , 0.5       ,  0.,  0.        ,  0.        , 1.1160254 ],
            [ 0.        , 0.        ,  1.,  1.        , -0.5       , 0.        ],
            [ 0.        , 0.        ,  0.,  0.5       ,  0.8660254 , -0.        ],
            [ 0.        , 0.        ,  0., -0.8660254 ,  0.5       , 0.        ],
            [ 0.        , 0.        ,  0.,  0.        ,  0.        , 1.        ],
            [ 1.        , 0.        ,  0.,  0.        ,  1.        , 0.        ],
            [ 0.        , 1.        ,  0., -1.        ,  0.        , 0.5       ],
            [ 0.        , 0.        ,  1.,  0.        , -0.5       , 0.        ],
            [ 0.        , 0.        ,  0.,  1.        ,  0.        , 0.        ],
            [ 0.        , 0.        ,  0.,  0.        ,  1.        , 0.        ],
            [ 0.        , 0.        ,  0.,  0.        ,  0.        , 1.        ]])
        mass_expected = np.array([8.80291473, 11.19708483, 20.])
        cg_expected = np.array([
            [-0.04426964, -0.00651564, -0.24841079],
            [ 0.602398  ,  0.03480387, -0.25124943],
            [ 0.5589382 , -0.0093859 ,  0.        ],])
        S_expected = np.array([ #  flipped n1 and n2...
            [-0.90180983, -0.43213312,  0.],
            [-0.43213312,  0.90180983,  0.],
            [-0.        ,  0.        ,  1.]])
        IS_expected = np.array([
            [ 5.623007  , -0.87682927,  0.34198812],
            [-0.87682927,  4.376953  , -0.6624477],
            [ 0.34198812, -0.6624477 ,  3.4319038 ]])
        II_expected = np.array([
            [ 5.623007  ,  0.87682927, -0.34198812],
            [ 0.87682927,  4.376953  ,  0.6624477 ],
            [-0.34198812,  0.6624477 ,  3.4319038 ]])
        Q_expected = np.array([
            [-8.8860166e-01, -2.6871899e-01, -3.7172186e-01],
            [-4.5867968e-01,  5.2056354e-01,  7.2015733e-01],
            [ 1.5101138e-05, -8.1043428e-01,  5.8582956e-01]])
        IQ_expected = np.array([
            [ 6.0756149e+00, -7.7424431e-08, -1.7855866e-07],
            [ 1.5370921e-08,  2.8930018e+00, -1.1259252e-07],
            [-3.5335850e-07,  1.3223152e-07,  4.4632468e+00]])
        assert np.allclose(D, D_expected)
        assert np.allclose(Mo, Mo_expected)
        assert np.allclose(S, S_expected), f'S:\n{S}\nS_expected:\n{S_expected}'
        assert np.allclose(mass, mass_expected)
        assert np.allclose(cg, cg_expected)
        assert np.allclose(IS, IS_expected)
        assert np.allclose(II, II_expected)
        assert np.allclose(IQ, IQ_expected)
        assert np.allclose(Q, Q_expected)

if __name__ == '__main__':
    unittest.main()
