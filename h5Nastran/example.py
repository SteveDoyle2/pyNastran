from __future__ import print_function

import numpy as np

from h5Nastran import H5Nastran

db = H5Nastran('./models/model_001.h5', 'w')  # , in_memory=True)
db.load_bdf('./models/model_001.bdf')
db.load_punch('./models/model_001.pch')

print(db.input.node.grid.identity)  # or db.input.node.grid.grid

domain_ids = range(10)
elements = [1, 5, 75]
f = db.result.elemental.element_force.quad4.search(elements, domain_ids)

print(f)

print(2 * f + 3 * f)

f = db.result.elemental.element_force.search(elements, domain_ids)

print(f.quad4)

# vec = db.input.coordinate_system.h5n_transformation.vector_to_basic([1, 1, 1], 1)

# print(vec)

# pynastran bdf
bdf = db.bdf

# currently, to modify the bdf and rewrite to h5,
# you'd need to modify the pynastran bdf, write it to a new file and create a new h5 database

# currently, the entire bdf is written to h5 as written by pynastran
# it can be loaded by doing db.load_bdf() without a filename
# the goal is to recreate the bdf file from only the h5 data

db.visualize()

db.close()
