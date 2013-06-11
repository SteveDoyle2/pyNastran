# -*- coding: utf-8 -*-
import os
import copy
from pyNastran.utils import print_bad_path
from pyNastran.op2.op2 import OP2

#op2_filename = r'D:\work\pynastran_0.7.0_py27\models\iSat\ISat_Launch_Sm_Rgd.op2'
op2_filename = r'D:\work\pynastran_0.7.0_py27\models\iSat\ISat_Launch_Sm_4pt.op2'
op2 = OP2(op2_filename)

#assert os.path.exists(op2_filename), print_bad_path(op2_filename)
op2.read_op2()

from pyNastran.utils import object_methods
from numpy import abs

#print op2.get_op2_stats()
print "loadcases =", op2.eigenvectors.keys()
eig1 = op2.eigenvectors[1]
modes = eig1.translations.keys()
#print "modes =", modes
mode2 = eig1.translations[2]
#print "mode2 =", mode2.keys()
mode2_node10 = mode2[10]
print "mode2_node10 =", mode2_node10
#m0 = modes[0]
#print "n1 = ", case3[m0]

plate_stress = op2.plateStress[1]
print "plate_stress_obj =", type(plate_stress)
print "plate_stress = ", plate_stress.__dict__.keys()
plate_stress.ovmShear.keys()  # 1-6
plate_stress.data_code

#type(plate_stress)

def abs_max_min(vals):
    maxval = max(abs(vals))
    i = vals.index(maxval)
    return vals[i]

object_methods(plate_stress)
plate_stress.isVonMises()  # True
ovm = plate_stress.ovmShear
#ovm.keys()  # modes 1-6
#ovm[6].keys()  # elements 1-3277
print "vals =", ovm[6][1000]['C']
ovm_mode6_eid1000 = abs_max_min(ovm[6][1000]['C'])
print "ovm_mode6_eid1000 =", ovm_mode6_eid1000

ps = copy.deepcopy(plate_stress)

ps = copy.deepcopy(plate_stress)
print "ps.mode_cycles", ps.mode_cycles
ps.mode_cycles = [0.0]

ps.oxx = {6: ps.oxx[6]}
ps.oyy = {6: ps.oyy[6]}
ps.txy = {6: ps.txy[6]}
ps.fiberCurvature = {1: ps.fiberCurvature[1]}
#print "fiberCurv", ps.fiberCurvature
ps.modes = [6]
ps.majorP = {6: ps.majorP[6]}
ps.minorP = {6: ps.minorP[6]}
ps.ovmShear = {6: ps.ovmShear[6]}

#object_methods(ps)
fname = os.path.join(r'C:\f.f06')
f = open(fname, 'wb')
ps.write_f06(['header1','header2'], 'stamp', f=f)
f.close()
#ps.write_f06?

