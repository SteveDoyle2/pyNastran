import os
import unittest
import numpy as np

import pyNastran
from pyNastran.bdf.bdf import BDF
from pyNastran.op2.op2 import OP2
from pyNastran.op2.data_in_material_coord import (data_in_material_coord,
        get_eids_from_op2_vector)
pkg_path = pyNastran.__path__[0]


prefix = os.path.join('test_dummy_wing_metallic', 'dummy_wing_metallic')

bdf = BDF(debug=False)
op2 = OP2(debug=False)
basepath = os.path.join(pkg_path, 'op2', 'test', prefix)
bdf.read_bdf(os.path.join(basepath + '.bdf'))
op2.read_op2(os.path.join(basepath + '.op2'))
op2_new = data_in_material_coord(bdf, op2)

vecname = 'cquad4_force'
#for vecname in stress_vectors:
subcase = 1
name = os.path.join('test_dummy_wing_metallic', '{0}_subcase_{1:02d}.txt'.format(vecname, subcase))
vector = getattr(op2_new, vecname)[subcase]
data = vector.data

eids = get_eids_from_op2_vector(vector)
check = eids != 0

ref_result = np.loadtxt(name)

print(ref_result)
print(data[:, check, :])

mcid = np.array([isinstance(bdf.elements[eid].thetaMcid, int) for eid in eids])

if not np.allclose(data[:, check, :][:, mcid, :], ref_result[mcid], rtol=0.02):
    print('MCID FAILED %r' % name)

if not np.allclose(data[:, check, :][:, ~mcid, :], ref_result[~mcid], rtol=0.02):
    print('THETA FAILED %r' % name)


