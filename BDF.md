# Introduction #

This section discusses the BDF class and how to use it.  Make sure you read the CrossReferencing and the DeveloperDocs sections.

# Basic Usage #
```
from pyNastran.bdf.bdf import BDF

mesh = BDF(debug=True,log=None)
# instantiates the class with the default set of cards
# debug - prints messages
# log   - if you pass in a python logging object that will be used
# instead of the dummy logger that is created internally

mesh.cardsToRead = set(['GRID','CQUAD4','PSHELL','MAT1','CORD2R'])
# not required, but lets you limit the cards
# if there's a problem with one

node100 = ['GRID',id,cp,x,y,z,cd,ps,seid]
# Note that None may be used for default values

mesh.add_card(node100,'GRID')  # adds the new GRID card to the bdf

mesh.readBDF(bdfFileIn,includeDir=None,xref=True) # reads the bdf
#includeDir lets you append path information onto the include files
#xref=True links up nodes, elements, properties, materials to make it easier to reference data

mesh.writeBDF(bdfFileOut) # writes the bdf
mesh.writeBDFAsPatran(bdfFileOut2) # writes the bdf with the properties and elements interspersed


c1 = mesh.Coord(1) # gets coordinate system 1
print c1.type      # CORD2R
print c1 # prints the card in small field format version (8-character width)

n1 = mesh.Node(1)   # gets element 1
print n1.type       # GRID
print n1.cp         # gets the Cp coordinate card
print n1.Cp()       # gets the Cp ID
print n1.cd         # gets the Cd coordinate card
print n1.Cd()       # gets the Cd ID
print n1.Position() # gets the position of the node in the global coordinate system (CORD2R.id=0)


e1 = mesh.Element(1) # gets element 1
print e1.type      # CQUAD4
print e1.nodeIDs() # [1,2,3,4]
print e1.nodes     # [nodeObj1,nodeObj2,nodeObj3,nodeObj4]
print e1.zOffset   # 0.2 - only for shells
print e1.Length()  # length of the element (only for Bars - CONROD, CROD, CBEAM, CBAR)
print e1.Area()    # area of the element   (only for CTRIA, CQUAD)
print e1.Volume()  # volume of the element (only for CTETRA, CPENTA, CHEXA)
print e1.Nsm()     # gets the non-structural mass
print e1.Thickness() # gets the element thickness (only for CQUAD/CTRIA)
print e1.pid.nsm   # gets the non-structural mass
# the function version will work for any
# element (assuming it has the quantity)

print e1.Mass()  # returns the mass if a BAR or SHELL element


p1 = mesh.Property(1) # gets property 1
print p1.type # PSHELL
print p1.t    # thickness (per the nastran QRG)

m1 = mesh.Material(1) # gets material 1
print m1.type # MAT1
print m.E()   # predefined on card
print m.Nu()  # predefined on card
print m.G()   # calculated based on E, nu
```

# Specific Examples #
Check out the DeveloperDocs for more detail on the specific method that was implemented in the various examples.

## 1.  Print all the Grid details (Nodes, Coords) to an Output File ##
```
from pyNastran.bdf.bdf import BDF

model = BDF()
model.readBDF('model.bdf',xref=True)
f = open('junk.out','wb')

#Method 1 - using objects

# GRIDs
for nid,node in sorted(model.nodes.items()):
    f.write(node)

# GRIDSET
if model.gridSet:
    f.write(model.gridSet)

# SPOINTs
if model.spoints:
    f.write(model.spoints)

# CORDx
for cid,coord in sorted(model.coords.items()):
    if cid!=0:  # if CID=0 is the global frame, skip it
        f.write(coord)

#----------------
#Method 2 - using built-in methods
f.write(model.write_nodes())
f.write(model.write_coords())

```
## 2.  Print the Element ID and associated Node and Property to an Output File ##
```
# note this skips rigidElements

# Method 1 - using objects
for eid,element in sorted(model.elements.items()):
    f.write(element)

for pid,prop in sorted(model.properties.items()):
    f.write(prop)

# Method 2 - using built-in method
f.write(model.write_elements_properties())

# Method 3 - using built-in methods
f.write(model.write_elements())
f.write(model.write_properties())
```

## 3.  Print the Element ID and its Configuration (e.g. CQUAD4, CTRIA3, etc.) to an Output File ##
```
# note this skips rigidElements

# Method 1 - using objects
for eid,element in sorted(model.elements.items()):
    msg = 'eid=%s type=%s\n' %(eid,element.type)
    f.write(msg)
```

## 4.  Given a Node, get the Elements Attached to that Node ##
```
# assume node 8101
nidToEidsMap = model.getNodeIDToElementIDsMap() # doesnt support 0d/1d elements yet
eids = nidToEidsMap[8101]

# convert to elements instead of element IDs
elements = []
for eid in eids:
    elements.append(model.Element(eid))

print eids
print "elements =\n ",elements
```

## 5.  Creating a List of Elements based on a Property ID ##
```
# assume pid=5
pidToEidsMap = model.getPropertyIDToElementIDsMap()
eids5  = pidToEidsMap[5]
eids10 = pidToEidsMap[10]

# convert to elements instead of element IDs
elements5  = []
elements10 = []
for eid in eids5:
    elements5.append(model.Element(eid))
for eid in eids10:
    elements10.append(model.Element(eid))
```

## 6.  Create a List of Elements based on a Material ID ##
```
# assume you want the eids for material 10
pidToEidsMap = model.getPropertyIDToElementIDsMap()
midToPidsMap = model.getMaterialIDToPropertyIDsMap()

pids10 = midToPidsMap[10]

eids = []
for pid in pids10:
    eids += pidToEidsMap[pid]

# convert to elements instead of element IDs
elements = []
for eid in eids:
    elements.append(model.Element(eid))
```