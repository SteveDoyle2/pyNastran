## GNU Lesser General Public License
## 
## Program pyNastran - a python interface to NASTRAN files
## Copyright (C) 2011-2012  Steven Doyle, Al Danial
## 
## Authors and copyright holders of pyNastran
## Steven Doyle <mesheb82@gmail.com>
## Al Danial    <al.danial@gmail.com>
## 
## This file is part of pyNastran.
## 
## pyNastran is free software: you can redistribute it and/or modify
## it under the terms of the GNU Lesser General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
## 
## pyNastran is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU Lesser General Public License
## along with pyNastran.  If not, see <http://www.gnu.org/licenses/>.
## 
#!/usr/bin/env python

import sys
import numpy as np
from cop4 import OP4
from pyNastran.op4.op4 import OP4 as pyOP4
from pyNastran.general.generalMath import printMatrix

for in_file in [ 'test/mat_b_dn.op4' ,
                 'test/mat_b_s1.op4' ,
                 'test/mat_b_s2.op4' ,
                 'test/mat_t_dn.op4' ,
                 'test/mat_t_s1.op4' ,
                 'test/mat_t_s2.op4' , 
                 'test/b_sample.op4',
                 ]:
    try:
        op4fh = OP4(in_file, 'r',encoding='utf-8')
    except:
        print('Failed to get header of %s, ignoring.' % (in_file))
        continue
    
    op4 = pyOP4()
    matrices = op4.readOP4(in_file)
    print('%s\n%s' % ('=' * 61, in_file))
    op4fh.print_header()
    for i in range(op4fh.nmat):
        a = op4fh.Load(nmat=1, skip=i)
        if a is None: 
            print('Failed to get %d-th matrix' % i)
        else:
            Name = op4fh.name[i]
            
            Format,A = matrices[Name]
            error = 0.
            #print type(a),type(A)
            if isinstance(A,np.ndarray):
                pass
            else: # sparse
                a = a.todense()
                A = A.todense()

            error = abs(a-A).max()

            if error>0:
                print "Name = |%s|" %(Name)
                print('%s:' % op4fh.name[i])
                print "error[%s] = %s" %(Name,error)
                
                print "cOP4:"
                print printMatrix(a)

                print "pyOP4:"
                print printMatrix(A)
                
                print "diff:"
                print printMatrix(a-A)
                print '----------------------------'
                sys.exit('stopping')
