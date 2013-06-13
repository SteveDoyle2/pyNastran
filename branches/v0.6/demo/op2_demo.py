# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

import os
import copy
from numpy import abs

from pyNastran.utils import print_bad_path
from pyNastran.op2.op2 import OP2
from pyNastran.utils import object_methods

# <codecell>

# define the object
# define the input file using the GUI
#op2 = OP2()

# <codecell>

#op2_filename = r'D:\work\pynastran_0.7.0_py27\models\iSat\ISat_Launch_Sm_Rgd.op2'
#op2_filename = r'D:\work\pynastran_0.7.0_py27\models\iSat\ISat_Launch_Sm_4pt.op2'
op2_filename = r'C:\Users\Steve\Dropbox\pyNastran_examples\iSat\ISat_Launch_Sm_4pt.op2'
#assert os.path.exists(op2_filename), print_bad_path(op2_filename)

# define the input file with a file path
op2 = OP2(op2_filename)

# <codecell>

# read the op2
op2.read_op2()

# <codecell>

print op2.get_op2_stats()

# <codecell>

# what modes did we analyze:  1 to 167
print "loadcases =", op2.eigenvectors.keys()

# get subcase 1
eig1 = op2.eigenvectors[1]

modes = eig1.translations.keys()

print "modes =", modes, '\n'

mode2 = eig1.translations[2]
mode2_node10 = mode2[10]
print "translation mode2_node10 =", mode2_node10
print "rotations mode2_node10 =", eig1.rotations[2][10]

# <codecell>

plate_stress = op2.plateStress[1]
print "plate_stress_obj =", type(plate_stress), '\n'

# the set of variables in the PlateStressObject
print "plate_stress = ", plate_stress.__dict__.keys(), '\n'

# list of parameters that define the object (e.g. what is the nonlinear variable name
print "data_code_keys = ", plate_stress.data_code.keys(), '\n'

# nonlinear variable name
name = plate_stress.data_code['name']
print "name = %r" % plate_stress.data_code['name']

print "list-type variables =", plate_stress.data_code['dataNames']

# the special loop parameter
# for modal analysis, it's "modes"
# for transient, it's "dts"
print "modes = ", plate_stress.modes # name + 's'


# extra list-type parameter for modal analysis; see dataNames
#print "mode_cycles =", plate_stress.mode_cycles

# <codecell>

def abs_max_min(vals):
    absvals = list(abs(vals))
    maxval = max(absvals)
    i = absvals.index(maxval)
    return vals[i]

print "methods =", object_methods(plate_stress)
if plate_stress.isVonMises():  # True
    ovm = plate_stress.ovmShear
else:
    oMaxShear = plate_stress.ovmShear

#ovm.keys()  # modes 1-6
#ovm[6].keys()  # elements 1-3277
print ""
print "[layer1, layer2, ...] =", ovm[6][1000]['C'] 
ovm_mode6_eid1000 = abs_max_min(ovm[6][1000]['C'])
print "ovm_mode6_eid1000 =", ovm_mode6_eid1000

# <codecell>

# see the difference between "transient"/"modal"/"frequency"-style results
# and "nodal"/"elemental"-style results

imode = 6  # mode 6; could just as easily be dt
iele = 10  # element 10
inode = 'C'
ilayer = 1

# result[imode][ielement][ilayer]
print "ps.mode_cycle =", plate_stress.mode_cycles[imode]
print "oxx =", plate_stress.oxx[imode][iele][inode][ilayer]
print "oyy =", plate_stress.oyy[imode][iele][inode][ilayer]
print "txy =", plate_stress.txy[imode][iele][inode][ilayer]
print "omax =", plate_stress.majorP[imode][iele][inode][ilayer]
print "omin = ", plate_stress.minorP[imode][iele][inode][ilayer]
print "ovm/maxShear =", plate_stress.ovmShear[imode][iele][inode][ilayer]

# result[imode][ielement][inode][ilayer]
if plate_stress.isFiberDistance():
    print "fiberDistance =", plate_stress.fiberCurvature[iele][inode][ilayer]
else:
    print "curvature =", plate_stress.fiberCurvature[iele][inode][ilayer]

# <codecell>

# write the F06 with Real/Imaginary or Magnitude/Phase
# only matters for complex results
op2.write_f06(r'C:\Users\Steve\Desktop\isat.f06', is_mag_phase=False, make_file=True, delete_objects=True)

