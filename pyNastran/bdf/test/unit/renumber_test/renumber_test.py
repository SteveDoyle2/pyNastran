from pyNastran.bdf.bdf import BDF

# renumber model entities. cross-references should be maintained

startnum = 9000

mesh = BDF(debug=False)
mesh.read_bdf('before.bdf', punch=True)
mesh.write_bdf('before_i.bdf', interspersed=False)

# renumber elements

for k, n in enumerate(mesh.elements.keys()):
    mesh.elements[n].eid = k + startnum

for k, n in enumerate(mesh.rigid_elements.keys()):
    mesh.rigid_elements[n].eid = k + startnum

# renumber properties

for k, n in enumerate(mesh.properties.keys()):
    mesh.properties[n].pid = k + startnum

# renumber materials

for k, n in enumerate(mesh.materials.keys()):
    mesh.materials[n].mid = k + startnum

# renumber coordinate systems

for k, n in enumerate(mesh.coords.keys()):
    mesh.coords[n].cid = k + startnum

# renumber grids

for k, n in enumerate(mesh.nodes.keys()):
    mesh.nodes[n].nid = k + startnum

# write the mesh

mesh.write_bdf('after_i.bdf', interspersed=False)

meshC = BDF(debug=False)
meshC.read_bdf('after_i.bdf', punch=True)
meshC.write_bdf('after.bdf')
