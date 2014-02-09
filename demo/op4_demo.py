# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

import os
from pyNastran.utils import print_bad_path
from pyNastran.op4.op4 import OP4
from numpy import float32, float64, int32, int64, product

# <codecell>

# define the OP4 object
op4 = OP4()

# print the docstring
#op4.read_op4?
help(op4.read_op4)

# <codecell>

# read the op4, will pop open a dialog box
#matrices = op4.read_op4()

#print matrices.keys()
#key = 'CMPLX'
#print matrices[key]

# <codecell>

op4_filename = r'C:\Users\Steve\Desktop\ISat_Launch_Sm_4pt.op4'
assert os.path.exists(op4_filename), print_bad_path(op4_filename)

#specify the file
matrices = op4.read_op4(op4_filename)

# <codecell>

# more ways to read a matrix

# only 1 matrix
matrices = op4.read_op4(op4_filename, matrix_names='FLAMA')

# 1 or more matrices
matrices = op4.read_op4(op4_filename, matrix_names=['FLAMA','UGEXT'])

# <codecell>

# extract a matrix
form, flama = matrices['FLAMA']
print "form =", form
print "type =", type(flama)

# <codecell>

from pyNastran.utils import object_methods
print "keys =", matrices.keys()

print "object_methods", object_methods(op4)
#op4.getTypeNWV?

# <codecell>

print matrices.keys()
form_flama, flama = matrices['FLAMA']
print "shape = ", flama.shape
print "flamat nvals =", product(flama.shape)

form_ugext, ugext = matrices['UGEXT']
print "form_ugext=%s type=%s" % (form_ugext, type(ugext[0,0]))
#print "ugext", ugext
print "ugext.shape =", ugext.shape
print "ugext nvals =", product(ugext.shape)

# <codecell>

print ugext[:,:]
#print flama

# <codecell>

#del matrices, op4

