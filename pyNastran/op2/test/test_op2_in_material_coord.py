import numpy as np

from pyNastran.bdf.bdf import BDF
from pyNastran.op2.op2 import OP2
from pyNastran.op2.data_in_material_coord import (data_in_material_coord,
        get_eids_from_op2_vector, force_vectors, stress_vectors,
        strain_vectors)

if __name__ == '__main__':
    bdf = BDF()
    op2 = OP2()
    bdf.read_bdf('./test_flat_plate_metallic/flat_plate_metallic.bdf')
    op2.read_op2('./test_flat_plate_metallic/flat_plate_metallic.op2')
    op2_new = data_in_material_coord(bdf, op2)
    subcase = 3
    for vecname in force_vectors + stress_vectors + strain_vectors:
        name = ('./test_flat_plate_metallic/{0}_subcase_{1:02d}.txt'.
                format(vecname, subcase))
        result = np.loadtxt(name)
        vector = getattr(op2_new, vecname)[subcase]
        data = vector.data
        eids = get_eids_from_op2_vector(vector)
        check = eids != 0
        assert np.allclose(data[:, check], result, rtol=0.001)



