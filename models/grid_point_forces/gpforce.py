import os
import numpy as np
import matplotlib.pyplot as plt
import pyNastran
from pyNastran.op2.op2 import read_op2
from pyNastran.bdf.mesh_utils.cut_model_by_plane import get_element_centroids
PKG_PATH = pyNastran.__path__[0]
op2_filename = os.path.join(PKG_PATH, '..', 'models', 'grid_point_forces', 'bar_grid_point_forces.op2')
#from pyNastran.bdf.bdf import read_bdf
#bdf_model = read_bdf()
model = read_op2(op2_filename, load_geometry=True, combine=True,
                 exclude_results=None, log=None, debug=True)
gpforce = model.grid_point_forces[1]
force = model.cbar_force[1]
#['station', 'bending_moment1', 'bending_moment2', 'shear1', 'shear2', 'axial', 'torque']
headers = force.get_headers()
istation = headers.index('station')
itime = 0
#ibending_moment1 = headers.index('bending_moment1')
ibending_moment2 = headers.index('bending_moment2')
#station = force.data[itime, :, istation]
#bending_moment1 = force.data[itime, :, ibending_moment1]
bending_moment2 = force.data[itime, :, ibending_moment2]

coord_out = model.coords[0]
nid_cp_cd, xyz_cid0, xyz_cp, icd_transform, icp_transform = model.get_xyz_in_coord_array(
    cid=0, fdtype='float64', idtype='int32')
all_nids = nid_cp_cd[:, 0]

all_eids, element_centroids_cid0 = get_element_centroids(model, idtype='int32', fdtype='float64')
stations = element_centroids_cid0[:-1, 0]
stations = np.linspace(0., 10., num=51)
#model.log.level = 'warning'
#print(stations)

nids_bar = []
nids_beam = []
for eid, elem in sorted(model.elements.items()):
    if elem.type == 'CBAR':
        nids_bar.append(elem.nodes)
    elif elem.type == 'BEAM':
        nids_beam.append(elem.nodes)
nids_bar = np.array(nids_bar, dtype='int32')
nids_beam = np.array(nids_beam, dtype='int32')
inid_bar = np.searchsorted(all_nids, nids_bar)
x1 = xyz_cid0[inid_bar[:, 0], 0]
x2 = xyz_cid0[inid_bar[:, 1], 0]
fig = plt.figure()
ax = fig.gca()
L = 10.0
x = xyz_cid0[:, 0].copy()
x.sort()
M = x ** 2 / 2

# F = wx
# M = wx^2/2
ax.plot(x1, bending_moment2[::2], 'o-', label='BM2', linewidth=3)
ax.plot(x2, bending_moment2[1::2], 'o--', )
ax.plot(L-x, M, 'o-', label='exact', linewidth=1)
ax.grid(True)

#nids = [1]
#eids = [1]
#force_out, moment_out, force_out_sum, moment_out_sum = gpforce.extract_interface_loads(
    #all_nids, all_eids,
    #coord_out, model.coords,
    #nid_cp_cd, icd_transform,
    #xyz_cid0,
    ##summation_point: Optional[NDArray3float]=None,
    #consider_rxf=True, itime=0,
    #debug=True, log=model.log)

# this one is empty...
nids = [1]
eids = [2]
force_out, moment_out, force_out_sum, moment_out_sum = gpforce.extract_interface_loads(
    nids, eids,
    coord_out, model.coords,
    nid_cp_cd, icd_transform,
    xyz_cid0,
    #summation_point: Optional[NDArray3float]=None,
    consider_rxf=True, itime=0,
    debug=True, log=model.log)
#coord0 = model.coords[0]
#gpforce.extract_interface_loads(nids: np.ndarray,
                                #eids: np.ndarray,
                                #coord_out=coord0,
                                #model.coords,
                                #nid_cd,
                                #icd_transform,
                                #xyz_cid0,
                                #summation_point=None,
                                #consider_rxf=True,
                                #itime=0,
                                #debug=True,
                                #log=model.log,
                                #idtype='int32')
# ----------------------------------------
nodes_list = list(model.nodes.keys())
nids = np.array(nodes_list, dtype='int32')
nids.sort()
#eids = np.ndarray(list(model.elements.keys()), dtype='int32')
#eids.sort()
# bar is [0,10] in x
force_sum, moment_sum = gpforce.shear_moment_diagram(
                                xyz_cid0,
                                all_eids,
                                nids,
                                icd_transform,
                                element_centroids_cid0,
                                model.coords,
                                nid_cp_cd,
                                stations,
                                coord_out,
                                idir=0, itime=0,
                                debug=True, log=model.log)
M2 = moment_sum[:, 1]
ax.plot(stations, M2, '*-', label='SMT')
ax.legend()
fig.show()
x = 1
