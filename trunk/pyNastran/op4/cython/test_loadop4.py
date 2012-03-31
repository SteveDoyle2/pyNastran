#!/usr/bin/env python

import numpy as np
import op4

"""
for in_file in [ '../test/mat_b_dn.op4' ,
                 '../test/mat_b_s2.op4' ,
                 '../test/mat_b_s1.op4' ,
                 '../test/mat_t_dn.op4' ,
                 '../test/mat_t_s1.op4' ,
                 '../test/mat_t_s2.op4' , ]:
"""
for in_file in [ '../test/mat_b_dn.op4' ,
                 '../test/mat_t_dn.op4' ]:
    try:
        fh = op4.File(in_file, 'r')
    except:
        print('Failed to get header of %s' % (in_file))
        raise SystemExit

    print('%s' % ('=' * 61))
    print('%s' % (in_file))
    op4.print_header(fh)

    for i in range(fh['nMat']):
        a = op4.Load(fh, nmat=1, skip=i)
        if a is None: 
            print('op4.Load returned None')
        else:
#           print('%2d. %-8s has type %s, contents of' % (
#                 i+1, fh['name'][i], type(a[0,0])))
            print(a)

# print 'It carries a reference to our deallocator: %s ' % a.base
# np.testing.assert_allclose(a, np.arange(10))
