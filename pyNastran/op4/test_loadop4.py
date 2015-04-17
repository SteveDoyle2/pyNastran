#!/usr/bin/env python
from six import iteritems
import os
import sys

from numpy import ndarray
import pyNastran
from pyNastran.op4.op4 import OP4
from pyNastran.utils.mathematics import print_matrix

pkg_path = pyNastran.__path__[0]
test_path = os.path.join(pkg_path, 'op4', 'test')
assert os.path.exists(test_path)

for in_file in ['mat_b_dn.op4',
                'mat_b_s1.op4',
                'mat_b_s2.op4',
                'mat_t_dn.op4',
                'mat_t_s1.op4',
                'mat_t_s2.op4',
                #'b_sample.op4',
                ]:
    in_file = os.path.join(test_path, in_file)
    op4 = OP4()
    matrices = op4.read_op4(in_file)

    error = 0
    for name, matrix in sorted(iteritems(matrices)):
            Format, A = matrices[name]

            #print type(a),type(A)
            if isinstance(A, ndarray):
                pass
            else:  # sparse
                #a = a.todense()
                A = A.todense()

            if error > 0:
                print("Name = %r" % name)
                print("error[%s] = %s" % (name, error))

                print("OP4:")
                print(print_matrix(A))

                print('----------------------------')
                sys.exit('stopping')
print('done')