from __future__ import print_function

from h5Nastran import H5Nastran

db = H5Nastran('./models/model_001.h5', 'w')
db.load_bdf('./models/model_001.bdf')
db.load_punch('./models/model_001.pch')

print(db.input.node.grid.identity)  # or db.input.node.grid.grid

domain_ids = [1, 2]
elements = [400002, 400111, 400198]
forces = db.result.elemental.element_force.quad4.search(domain_ids, elements)

print(forces)

# pynastran bdf
bdf = db.bdf

# currently, to modify the bdf and rewrite to h5,
# you'd need to modify the pynastran bdf, write it to a new file and create a new h5 database

# currently, the entire bdf is written to h5 as written by pynastran
# it can be loaded by doing db.load_bdf() without a filename
# the goal is to recreate the bdf file from only the h5 data

db.close()
