import sys
import numpy as np
from pyNastran.converters.stl.stl import read_stl

stl_filename = sys.argv[1]
fld_filename = sys.argv[2]

model = read_stl(stl_filename)
xyz = model.nodes
pressure = np.arange(len(xyz))

with open(fld_filename, 'w') as fld_file:
    fld_file.write("""FIELD: [loadcase_name] : [TABLE]
FIELD LOCK STATE: [NO]
DUPLICATE_VALUE_OPTION: [0]
PARAMETERIZE INDEPENDENT DOMAIN: [NO]
PERSIST INTERPOL: [NO]
CREATE INTERPOLATION: [NO]
FALLBACK DEFAULT INTERPOLATOR: [YES]
INTERPOL [10]
VALUES OUTSIDE: [0]
REMOVE DELAUNAY SLIVERS: [NO]
MAP: [1]
INDEP VAR: [x] : [Length] : [in] : [0]
BOUNDS: [-100.000[ : [YES] : [100.000] : [YES] : [1000] : [207.000]
INDEP VAR: [y] : [Length] : [in] : [0]
BOUNDS: [-100.000[ : [YES] : [100.000] : [YES] : [1000] : [-10.000]
INDEP VAR: [z] : [Length] : [in] : [0]
BOUNDS: [-100.000[ : [YES] : [100.000] : [YES] : [1000] : [37.000]
DEP VAR: [pressure] : [Pressure] : [lbf/in^2(psi)] : [0]
START DATA
""")
    for xyzi, pressurei in zip(xyz, pressure):
        fld_file.write(f'{xyzi[0]}, {xyzi[1]}, {xyzi[2]}, {pressurei}\n')
    fld_file.write('END DATA\n')
