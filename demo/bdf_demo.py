# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

import os
import pyNastran
print pyNastran.__file__
print pyNastran.__version__


from pyNastran.bdf.bdf import BDF
from pyNastran.utils import object_attributes, object_methods

# <codecell>

mode = 0

# create the BDF object
bdf = BDF()

# read the file from the GUI
# don't cross-reference
#bdf.read_bdf(xref=False)

bdf_filename = r'C:\Users\Steve\Dropbox\pyNastran_examples\iSat\ISat_Launch_Sm_Rgd.dat'
bdf.read_bdf(bdf_filename, xref=False)

# <codecell>

#bdf_filename = r'D:\work\pynastran_0.7.0_py27\models\iSat\ISat_Launch_Sm_Rgd.dat'
bdf_filename = r'C:\Users\Steve\Dropbox\pyNastran_examples\iSat\ISat_Launch_Sm_Rgd.dat'

# read the file as a path
bdf2 = BDF()
bdf2.read_bdf(bdf_filename, xref=True)

# <codecell>

print "attributes =", object_attributes(bdf)
print ""
print "methods =",object_methods(bdf)
print bdf.card_stats()
print bdf.card_count

# <codecell>

# explanation of cross-referencing

# no cross referencing (xref=False)
cquad = bdf.elements[1]
nid1 = cquad.nodes[0]
n1 = bdf.nodes[nid1]
cd4 = n1.cd
c4 = bdf.coords[cd4]
print "i xref=False", c4.i
#print object_attributes(c4)

# cross referenced (xref=True)
print "i xref=True", bdf2.elements[1].nodes[0].cd.i

# how is it done?
cquad.nodes[0] = n1
print cquad.nodes[0]

# <codecell>

# some Grid methods
n1 = bdf2.nodes[1]
print n1

# the comment
c1 = bdf2.nodes[1].comment()
c2 = bdf2.nodes[2].comment()
print "c1=%r" % c1
print "c2=%r" % c2


# get the position of a node
# in the local cooordinate system
print "xyz =", n1.xyz

# in the global frame
print "position =", n1.Position()

# in an arbitrary frame
print "wrt5 =", n1.PositionWRT(bdf2, 5)
print "wrt4 =", n1.PositionWRT(bdf2, 4)

# <codecell>

n1 = bdf2.nodes[1]
n1.xyz[1] = -7.5
print "repr  =", n1.reprFields()
print "raw   =", n1.rawFields()

#n1.xyz[1] = 100000000000.
print "repr2 =", n1.reprFields()
print n1

# <codecell>

print "mass = ", bdf2.Mass()

# <codecell>

eid100 = bdf2.elements[100]
print eid100
print "nodes =", eid100.nodes
print "--node0--\n", eid100.nodes[0]
print "--cd--\n", eid100.nodes[0].cd
print "cd.cid", eid100.nodes[0].cd.cid

print "area =", eid100.Area()
print "mass =", eid100.Mass()
print "--pid--\n", eid100.pid
print "pid.pid =", eid100.pid.pid
print "pid.Pid() =", eid100.Pid()

print eid100.pid.mid1
print "type =", eid100.pid.mid1.type
print "nu12 =", eid100.pid.mid1.nu12
print "mass =", eid100.Mass()

# <codecell>

import getpass
name = getpass.getuser()
os.chdir(os.path.join(r'C:\Users', name, 'Desktop'))

# <codecell>

pwd

# <codecell>

# write the bdf
bdf2.write_bdf('fem.bdf', interspersed=True, size=8)
!tail "fem.bdf"

