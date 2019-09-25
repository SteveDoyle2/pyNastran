from collections import Counter
import numpy as np

NX_RESULT_TABLES = [
    # GEOM table
    b'ICASE',

    #-----------------------
    # OESVM1  - OES Table of           element stresses
    # OESVM1C - OES Table of composite element stresses
    #           for frequency response analysis that includes von Mises stress
    #           output in SORT1 format.
    b'OESVM1', b'OESVM2', b'OSTRVM1', b'OSTRVM2',
    b'OESVM1C', b'OSTRVM1C',

    b'OES2C',
    b'OSTR2C',

    b'OESPSD2C', b'OSTPSD2C',
    b'OSTRRMS1', b'OSTRMS1C',
    b'OSTRNO1', b'OSTNO1C',

    # sol 401?
    b'OESNL2',
    b'OSTR1IN',  # OES output table of initial strains at corner grids in the basic coordinate system
    b'OSTR1G',   # Table of total strain at Gauss points in SORT1 format
    b'OSTR1PLG', # Table of plastic strain at Gauss points in SORT1 format
    b'OSTR1THG', # Table of thermal strain at Gauss points in SORT1 format
    b'OSTR1ELG', # Table of elastic strain at Gauss points in SORT1 format
    b'OSTR1TH',  # Table of thermal strain in SORT1 format
    b'OSTR1EL',  # Table of elastic strain in SORT1 format
    b'OSTR1ING', # OES output table of initial strains at corner Gauss points in the basic coordinate system

    b'OES1G',   # Grid point stress or strain table in SORT1 format and interpolated from the centroidal stress table, OES1M.

    #----------------------
    b'MDICT', b'BDICT', b'KDICTP', b'MDICTP',

    #----------------------
    # displacement/velocity/acceleration/eigenvector/temperature
    # OUGV1  - Displacements in the global coordinate system
    b'OUGV1',

    b'OPHIG', # Eigenvectors in the basic coordinate system.

    #----------------------
    # mpc forces - gset - sort 1
    b'OQMG1',

    #----------------------
    # forces
    # OEF1X - Element forces with intermediate (CBAR and CBEAM) station forces
    #         and forces on nonlinear elements
    b'OEF1X',

    # ---------------------------------------------
    # nx2019.2

    # geometry
    b'GPDTS',
    # results
    b'LAMAS', # Normal modes eigenvalue summary table for the structural portion of the model
    b'LAMAF', # Normal modes eigenvalue summary table for the fluid portion of the model
    b'GPLS',
    b'TRMBU',  # Transfomration matrices from undeformed to basic
    b'TRMBD',  # Transformation matrices from deformed to basic
    b'PVTS', # PVT0?
    b'OGPFB2',   # Table of grid point forces in SORT2 format
    b'OEFMXORD', # Data block that contains a list of element IDs with maximum frequency and element order diagnostics for FEMAO solution
    b'BOPHIGF',  # Eigenvectors in the basic coordinate system for the fluid portion of the model.
    b'BOPHIGS',  # Eigenvectors in the basic coordinate system for the structural portion of the model.
    b'OUGMC1',   # Table of modal contributions for displacements, velocities, accelerations.
    b'OUXY1',    # Table of displacements in SORT1 format for h-set or d-set.
    b'OUXY2',    # Table of displacements in SORT2 format for h-set or d-set.
    b'OSLIDE1',  # Incremental and total slide output for contact/glue.
    b'OSLIDEG1', # Glue slide distance output
    b'OBCKL',    # Table of load factor vs. cumulative arc-length in SORT2 format
    b'OJINT',    # Table of J-integral for a crack defined by CRAKTP.
    b'OPHSA',
    b'ONMD',     # Normalized material density for topology optimization output
    b'OACPERF',  # Performance data that indicates computation time in seconds and memory consumed in GB per frequency per subcase for FEMAO analysis.
    b'OGSTR1',   # Output table of grid point strains of superelement
    b'OELAR',    # Element status (active or inactive)
    b'OPNL2',   # Table of nonlinear loads in SORT2 format for the h-set or d-set.
    b'OSHT1',   # Output shell element thickness results.
    b'OTEMP1',  # Grid point temperature output
    b'OEFIIS',  # Data block for inter-laminar shear failure indices.
    b'OBOLT1',  # Bolt output data block
    b'PSDF',    # Power spectral density table.

    b'ODAMGPFE', # Table of damage energy for ply failure
    b'ODAMGPFD', # Table of damage values for ply failure
    b'ODAMGPFS', # Table of damage status for ply failure
    b'ODAMGPFR', # Table of crack density for ply failure EUD model from SOL 401. Crack density at corner grids on the middle of plies. The values are unitless

    b'OCONST1',  # Contact status in SORT1 format

    b'OCPSDF',  # Output table of cross-power-spectral-density functions.
    b'OCCORF',  # Output table of cross-correlation functions.
    b'OCPSDFC', # Table of cross-power spectral density functions for composites.
    b'OCCORFC', # Table of cross-correlation functions for composites.

    b'OESXRM1C', # Table of composite element RMS stresses in SORT1 format for random analysis that includes von Mises stress output.
    b'OESXNO1C',
    b'OESXNO1',

    b'OEDE1', # Elemental energy loss.


    # grid point pressure
    b'OPRNO1',   # SORT1 - NO
    b'OPRRMS1',  # SORT1 - RMS
    b'OPRPSD2',  # SORT2 - PSD
    b'OPRATO2',  # SORT2 - AUTO
    b'OPRCRM2',  # SORT2 - CRMS

    # modal contribution
    b'OUGMC2', # Table of modal contributions for displacements, velocities, accelerations.
    b'OQGMC1', # Table of modal contributions of single point constraint forces - SORT1
    b'OQGMC2', # Table of modal contributions of single point constraint forces - SORT2
    b'OESMC1', # Element stress modal contributions - SORT1
    b'OESMC2', # Element stress modal contributions - SORT2
    b'OSTRMC1', # Table of modal contributions of element strains - SORT1
    b'OSTRMC2', # Table of modal contributions of element strains - SORT2
    b'OEFMC1', # Table of modal contributions of element forces - SORT1
    b'OEFMC2', # Table of modal contributions of element forces - SORT2

    # modal strain energy
    b'OMSEC1',  # Constant modal strain energy - SORT1
    b'OMSEC2',  # Constant modal strain energy - SORT2
    b'OMECON1', # Constant total modal energies - SORT1
    b'OMECON2', # Constant total modal energies - SORT2
    b'OMEOSC1', # Oscillating total modal energies - SORT1
    b'OMEOSC2', # Oscillating total modal energies - SORT2
    b'OMKEC1',  # Constant modal kinetic energies - SORT1
    b'OMKEC2',  # Constant modal kinetic energies - SORT2
    b'OMKEO2',  # Oscillating modal kinetic energies - SORT2
    b'OMSEO1',  # Oscillating modal strain energies - SORT2
    b'OMSEO2',  # Oscillating modal strain energies - SORT2
    b'OMKEO1',  # Oscillating modal kinetic energies - SORT1

    # radiated power
    b'OERP',     # Equivalent radiated power output.
    b'OERPEL1',  # Element equivalent radiated power (element output)
    b'OERPEL2',  # Element equivalent radiated power output.

    # acoustic
    b'OUGPC1',  # Table of panel contributions - SORT1
    b'OUGPC2',  # Table of panel contributions - SORT2
    b'OUGF1',   # Acoustic pressures at microphone points in SORT1 format
    b'OUGF2',   # Acoustic pressures at microphone points in SORT2 format
    b'OUGGC1',  # Table of grid contributions - SORT1
    b'OUGGC2',  # Table of grid contributions - SORT2
    b'OUGRC1',  # Reciprocal panel contributions - SORT1
    b'OUGRC2',  # Reciprocal panel contributions - SORT2
    b'BOUGF1',  # Data block of acoustic pressures at microphone points in SORT1 format - basic frame

    # acoustic acceleration
    b'OACCQ',    # Acoustic coupling quality
    b'OACINT1',  # Acoustic intensities at microphone points - SORT1
    b'OACINT2',  # Acoustic intensities at microphone points - SORT2
    b'OACVELO1', # Acoustic velocities at microphone points - SORT1
    b'OACVELO2', # Acoustic velocities at microphone points - SORT2
    b'OACPWR2',  # Acoustic power for AML regions and GROUPs of 2D elements - SORT2
    b'OACPWRI2', # Acoustic incident power - SORT2
    b'OACPWRT2', # Transmitted acoustic power for AML regions and GROUPs of 2D elements - SORT2
    b'OACTRLS2', # Acoustic transmission loss - SORT2

    # random acoustic
    b'OAPPSD2', # Acoustic power for the PSD function - SORT2
]

if len(NX_RESULT_TABLES) != len(np.unique(NX_RESULT_TABLES)):
    counter = Counter(NX_RESULT_TABLES)
    _MSG = 'Invalid count:\n'
    for key, cvaluei in counter.items():
        if cvaluei != 1:
            _MSG += '%s = %s\n' % (key, cvaluei)
    raise RuntimeError(_MSG)

NX_TABLE_CONTENT = {
    # nx 8.5
    0: '',
    1: 'OUG - Displacement vector',
    2: 'OPG - Load vector',
    3: 'OQG - SPC/MPC Force vector',
    4: 'OEF - Element force/flux',
    5: 'OES - Element stress/strain',
    6: 'LAMA - Eigenvalue summary',
    7: 'OUG - Eigenvector',
    8: 'Grid Point Singularity Table (obsolete)',
    9: 'OEIGS - Eigenvalue analysis summary',
    10: 'OUG - Velocity vector',
    11: 'OUG - Acceleration vector',
    12: 'OPG - Nonlinear force vector',
    13: 'OGPWG - Grid point weight generator',
    14: 'OUG - Eigenvector (solution set)',
    15: 'OUG - Displacement vector (solution set)',
    16: 'OUG - Velocity vector (solution set)',
    17: 'OUG - Acceleration vector (solutin set)',
    18: 'OEE - Element strain energy',
    19: 'OGF - Grid point force balance',
    20: 'OES - Stresses at grid points',
    21: 'OES - Strain/curvature at grid points',
    22: 'OELOF1 - Element internal forces/moments',
    23: 'OELOP1 - Summation of element oriented forces on adjacent elements',
    24: 'OEP - Element pressures',
    25: 'OEF - Composite failure indices',
    26: 'OGS - Grid point stresses (surface)',
    27: 'OGS - Grid point stresses (volume - direct)',
    28: 'OGS - Grid point stresses (volume - princial)',
    29: 'OGS - Element stress discontinuities (surface)',
    30: 'OGS - Element stress discontinuities (volume - direct)',
    31: 'OGS - Element stress discontinuities (volume - princial)',
    32: 'OGS - Grid point stress discontinuities (surface)',
    33: 'OGS - Grid point stress discontinuities (volume - direct)',
    34: 'OGS - Grid point stress discontinuities (volume - princial)',
    35: 'OGS - Grid point stresses (plane strain)',
    36: 'OEE - Element kinetic energy',
    37: 'OEE - Element energy loss',

    38 : 'OMSEC - Constant modal strain energy',
    39 : 'OMSED - Oscillating modal strain energy',
    40 : 'OMKEC - Constant modal kinetic energy',
    41 : 'OMKED - Oscillating modal kinetic energy',
    42 : 'OMECON - Constant total modal energy',
    43 : 'OMEOSC - Oscillating total modal energy',
    44 : 'OUGMC - Displacement/velocity/acceleration modal contributions',
    45 : 'OEFMC - Element force modal contributions',
    46 : 'OESMC - Element stress modal contributions',
    47 : 'OSTRMC - Element strain modal contributions',
    48 : 'OQGMC - SPC force modal contributions',
    49 : 'OUGPC - Panel contributions',
    50 : 'OUGGC - Grid contributions',
    51 : 'OUGRC - Reciprocal panel contributions',
    61 : 'OGK - Gasket Element Results',
    62 : 'OBC - Contact Pressure and Traction',
    63 : 'OQG - Contact Force Results',
    64 : 'OSPDSI - Contact Separation Distance - Initial',
    65 : 'OSPDS - Contact Separation Distance',
    66 : 'OBG - Glue force results (normal and in-plane tractions)',
    67 : 'OQG - Glue force resutls',
    68 : '??? - Tosca normalized material properties',
    # nx 2019.2
    72 : 'OTEMP1 - grid point temperatures',
}
