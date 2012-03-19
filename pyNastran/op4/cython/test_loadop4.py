#!/usr/bin/env python

import numpy as np
import op4

type_map = { 1 : 'RS' ,   # real    single precision
             2 : 'RD' ,   # real    double precision
             3 : 'CS' ,   # complex single precision
             4 : 'CD' ,}  # complex double precision

"""
for in_file in [ '../test/mat_b_dn.op4' ,
                 '../test/mat_b_s2.op4' ,
                 '../test/mat_b_s1.op4' ,
                 '../test/mat_t_dn.op4' ,
                 '../test/mat_t_s1.op4' ,
                 '../test/mat_t_s2.op4' , ]:
"""
for in_file in [ '../test/mat_t_dn.op4' ]:
    try:
        fh = op4.File(in_file, 'r')
    except:
        print('Failed to get header of %s' % (in_file))
        raise SystemExit


    print('%s' % ('=' * 61))
    print('%s' % (in_file))
    print('    %-8s %5s %5s %8s %8s %2s %2s %2s %9s' % (
          'Name', 'nRow', 'nCol', 'nStr', 'nNnz', 'T', 'Fr', 'Dg', 'Offset'))
    for i in range(fh['nMat']):
        print('%2d. %-8s %5d %5d %8d %8d %2s %2d %2d %9d' % (
               i+1          ,
               fh['name'][i],
               fh['nRow'][i],
               fh['nCol'][i],
               fh['nStr'][i],
               fh['nNnz'][i],
               type_map[fh['type'][i]],
               fh['form'][i],
               fh['digits'][i],
               fh['offset'][i],))

    a = op4.Load(fh)
#   a = op4.Load(fh, n_mat=3, n_skip=5)

# a = op4.load2(10)

print 'The array created is %s' % a
# print 'It carries a reference to our deallocator: %s ' % a.base
# np.testing.assert_allclose(a, np.arange(10))
