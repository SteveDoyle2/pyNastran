from __future__ import print_function
from numpy import zeros

from pyNastran.gui.qt_files.result import ScalarResult, Results
from pyNastran.converters.cart3d.cart3dIO import Cart3dGeometry, Cart3dResult

class DispResult(object):
    def __init__(self, xyz_cid0, dxyz, uname='test', scale=40.,
                 location='node', title='Disp', fmt='%.3f'):
        self.xyz_cid0 = xyz_cid0
        self.dxyz = dxyz

        self.nnodes = dxyz.shape[0]
        self.uname = uname
        self.title = title
        self.location = location
        self.scale = scale
        self.fmt = fmt

    def get_methods(self, in_args):
        return [
            (self.location, 'x'),
            (self.location, 'y'),
            (self.location, 'z'),
            (self.location, 'xyz'),
        ]

    def get_data(self, in_args, method):
        print('method', method)
        location, var = method
        dxyz = zeros((self.nnodes, 3), dtype='float32')
        if var == 'x':
            dxyz[:, 0] = self.dxyz[:, 0]
        elif var == 'y':
            dxyz[:, 1] = self.dxyz[:, 1]
        elif var == 'z':
            dxyz[:, 2] = self.dxyz[:, 2]
        elif var == 'xyz':
            dxyz = self.dxyz
        else:
            raise NotImplementedError(method)

        if self.scale != 1.0:
            out = self.xyz_cid0 + dxyz[:, :3]
        else:
            out = self.xyz_cid0 + dxyz[:, :3] * self.scale
        return out

    def get_uname(self, i):
        return self.uname

    def __repr__(self):
        msg = 'DispResult\n'
        msg += '    uname=%r\n' % self.uname
        msg += '    locations=%r\n' % self.location
        msg += '    methods={x, y, z, xyz}\n'
        return msg


def spike():
    import numpy as np
    import numexpr as ne
    a = np.arange(10)
    b = np.arange(0, 20, 2)
    c = ne.evaluate("2*a+3*b")
    print(c)


def main():
    #spike()
    #return
    from numpy import zeros, ones
    #if new:
        #cases_new[0] = (ID, nids, 'NodeID', 'node', '%i')
        #cases_new[1] = (ID, eids, 'ElementID', 'centroid', '%i')
    #else:
        # current
        #cases[(ID, 0, 'NodeID', 1, 'node', '%i')] = nids

    nids = zeros(3, dtype='int32')
    eids = ones(3, dtype='int32') + 3

    nids2 = ScalarResult(nids, uname='a', title='NodeID',    dim=1, location='node', fmt='%i')
    eids2 = ScalarResult(eids, uname='b', title='ElementID', dim=1, location='centroid', fmt='%i')
    nids3 = ScalarResult(eids, uname='c', title='NodeID',    dim=1, location='node', fmt='%i')

    nodes = nids
    elements = eids
    regions = elements
    rho = nids
    rhoU = nids
    rhoV = nids
    rhoW = nids
    rhoE = nids
    cart3d_geo = Cart3dGeometry(nodes, elements, regions, labels=None)
    cart3d_res = Cart3dResult(rho, rhoU, rhoV, rhoW, rhoE, labels=None)
    #num_expr = ['a', ' + 2 * ', 'c']
    #nids2_nids3 = NumExprResult(num_expr, uname='comboID', title='combo', dim=1, location='node', fmt='%i')

    #def adder(results, method):
        #return results['a'].data + 2 * results['c'].data ** 2 / 3.

    #nids2_nids3 = FuncResult(adder, uname='comboID', title='combo', dim=1, location='node', fmt='%i')

    dxyz = zeros((5, 6), dtype='float32')
    xyz = zeros((5, 3), dtype='float32')
    disp = DispResult(xyz, dxyz, title='disp', fmt='%.3f')
    res = Results(
        #nids2, eids2, nids3, nids2_nids3,
        cart3d_geo,
        #cart3d_res,
        disp,
    )


    ID = 1
    #cases = [
        #(ID, 0, nids, 'NodeID', None, 'node', '%i'),
        #(ID, 1, eids, 'ElementID', None, 'centroid', '%i'),
        #(ID, 2, eids, 'ElementID', None, 'centroid', '%i'),
    #]

    cases2 = [
        # result, (result_index_in_results, name)
        (cart3d_geo, (0, 'nodeID')),
        (cart3d_geo, (0, 'elementID')),
        (cart3d_res, (1, 'rho')),

        (res, (0, 'NodeID')),
        (disp, (0, 'tx')),
        (disp, (0, 'ty')),
        (disp, (0, 'tz')),
        #(res, (1, 'rho')),
        #(res, (1, 'rhoU')),
        #(res, (1, 'rhoV')),
        #(res, (1, 'rhoW')),
        #(res, (1, 'rhoE')),
    ]
    form = [
        # just like the old form
        # title, case_id, next_level
        ('rho', 0, []),
        ('rhoU', 1, []),
        ('rhoV', 2, []),
        ('rhoW', 3, []),
        ('rhoE', 4, []),
    ]

    for case, in_args in cases2:
        print("case    =", case)
        print("in_args =", in_args)
        methods = case.get_methods(in_args)
        print("methods =", methods)
        for method in methods:
            print('****', in_args)
            val = case.get_data(in_args, method)
            print('%s = %s' % (case.get_uname(in_args), val))

if __name__ == '__main__':
    main()
