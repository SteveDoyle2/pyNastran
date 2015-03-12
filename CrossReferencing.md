# Introduction #

Cross-Referencing the BDF is a way of linking nodes, elements, properties, materials, etc. together in a way that makes it easier for the code to access data members.  If you're not cross-referencing the BDF, accessing data takes a lot more work.

# Example With/Without Cross-Referencing #

## Initialize the Model ##
```
import numpy
from pyNastran.bdf.bdf import BDF

mesh = BDF() # instantiates the class with the default set of cards
```

## Getting the Position of a Grid Point ##
```
#---NODES---
# with cross-referencing
mesh.read_bdf(fname, xref=True) # reads the bdf
n1 = mesh.Node(1)   # gets element 1
print n1.type       # GRID
print n1.cp         # gets the Cp coordinate card
print n1.Cp()       # gets the Cp ID
print n1.cd         # gets the Cd coordinate card
print n1.Cd()       # gets the Cd ID
print n1.Position() # gets the position of the node in the global coordinate system (CORD2R.id=0)

# without cross-referencing
mesh.read_bdf(fname, xref=False) # reads the bdf
n1 = mesh.Node(1)   # gets element 1
print n1.type       # GRID
print n1.cp         # gets the Cp ID
print n1.Cp()       # gets the Cp ID
print n1.cd         # gets the Cd ID
print n1.Cd()       # gets the Cd ID
print n1.Position() # crashes


# how to get Position()
cp = n1.Cp()
coord = mesh.Coord(cp)
position = coord.transformToGlobal(n1.xyz) # equivalent to n1.Position()
```

## Getting the Mass of an Element ##
```
#----ELEMENTS---
# with cross-referencing
mesh.read_bdf(fname, xref=True)        # reads the bdf
element1 = mesh.Element(1)  # assume a CROD
print element1.Mass()       # gets the mass

# without cross-referencing
mesh.read_bdf(fname, xref=False) # reads the bdf
element1  = mesh.Element(1)  # assume a CROD with pid=1, nids=[1,2]
pid1 = element1.Pid()

nodes = element1.nodeIDs()
n1 = mesh.Node(nodes[0])
n2 = mesh.Node(nodes[1])
L = numpy.linalg.norm(n1.xyz-n2.xyz)  # assuming global XYZ coodinates

property1 = mesh.Element(pid1)  # assume a PROD with mid=1,nsm=1.0, A=2.0
mid1 = property1.Mid()
nsm  = property1.nsm  # non-structural mass
A    = property1.A    # area

material1 = mesh.Material(mid1) # assume n, with rho=3.0
rho = material1.rho  # density

mass = L*(A*rho+nsm)  # gets the mass; equivalent to element1.Mass()

```

So in short (assuming it works on your model), cross-referencing is really nice and lets you make use of the true functionality of the software.

## Why Cross Referencing Fails ##
  * You're missing cards in the BDF.
  * There is a reference to an unsupported card.  For example, the PBCOMP (CBEAM Property) is not supported.  It is referenced by a CBEAM, so the property cannot be found.  The PBEAM3 is also not supported, but the CBEAM3 is not supported either.  In this case you're alright.
  * All element midside nodes must exist.