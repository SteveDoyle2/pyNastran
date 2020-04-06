from pyNastran.bdf.bdf import read_bdf
from pyNastran.op2.op2 import read_op2

bdf_filename = 'beam_modes.dat'
op2_filename = 'beam_modes_m2.op2'

#2) Open an op2
op2_model = read_op2(op2_filename)
print(op2_model.get_op2_stats())
subcase = 1

# 3a) Load grid point locations
model = read_bdf(bdf_filename)

# I'm assuming you don't have SPOINTs/EPOINTs
xyz_cid0 = {}
for nid, node in model.nodes.items():
    xyz_cid0[nid] = node.get_position()

# 3b) generalized mass, eigenvalues or frequencies, and mode shapes from a 103 solution
# it's rare you'd have more than 1 key, so we'll just grab the 0th key
eigenvalue_keys = list(op2_model.eigenvalues.keys())
title = eigenvalue_keys[0]

# grab the data
#self.mode = np.zeros(nmodes, dtype='int32')
#self.extraction_order = np.zeros(nmodes, dtype='int32')
#self.eigenvalues = np.zeros(nmodes, dtype='float32')
#self.radians = np.zeros(nmodes, dtype='float32')
#self.cycles = np.zeros(nmodes, dtype='float32')
#self.generalized_mass = np.zeros(nmodes, dtype='float32')
#self.generalized_stiffness = np.zeros(nmodes, dtype='float32')
eigenvalue_obj = op2_model.eigenvalues[title]
eigenvalues = eigenvalue_obj.eigenvalues
frequencies = eigenvalue_obj.cycles
generalized_mass = eigenvalue_obj.generalized_mass


#4) Put them into numpy arrays so I can do some math with them
eigenvector_obj = op2_model.eigenvectors[1]

# take your pick on which form, but probably eigenvectors2
eigenvectors1 = eigenvector_obj.data # (nmodes, nnodes, 6)
eigenvectors2 = eigenvector_obj.get_phi() # (ndof, nmodes)

#5) I will also need node ids
nids = eigenvector_obj.node_gridtype[:, 0]

print(eigenvalues, frequencies, generalized_mass)
