import os
configpath  = 'inputs'
workpath = 'outputsFinal'

aero_format = 'nastran'
structural_call = 'nastran scr=yes bat=no '
#structural_symmetry = None  # no symmetry
structural_symmetry = '+y'  # symmetry

aero_format = 'cart3d'
aero_input = 'Components.i.tri'
#aero_output = 'Components.i.triq'
aero_output = 'Cart3d_35000_0.825_10_0_0_0_0.i.triq'
aero_command = './COMMAND > command.out'
#aero_symmetry = '+y'
aero_symmetry = None  # full model

# unused
isubcase = 1


Mach = 0.825
pInf = 499.3        # psf, alt=35k (per Schaufele p. 11)
pInf = pInf / 144.  # convert to psi
qInf = 1.4 / 2. * pInf * Mach**2.
Sref = 1582876.  # inch^2
Lref = 623.  # inch
xref = 268.  # inch

