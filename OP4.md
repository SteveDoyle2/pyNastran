# Introduction #
The OP4 reader is a an ASCII or binary reader for Nastran OP4 files.  There are two readers.  One is pure Python and one is written in Cython, so it requires building.

See the InstallationGuide for how to install the Cython OP4 Reader.

# Matrix Forms #
The matrix form will be required when the OP4 writer is completed.

| **Form** | **Description** |
|:---------|:----------------|
| 1 | Square |
| 2 | Rectangular |
| 3 | Diagonal |
| 4 | Lower triangular factor (see note) |
| 5 | Upper triangular factor (see note) |
| 6 | Symmetric |
| 8 | Identity |
| 9 | Pseudoidentity |
| 10 | Cholesky factor (see note) |
| 11 | Trapezoidal factor |
| 13 | Sparse lower triangular factor (see note)|
| 15 | Sparse upper triangular factor (see note)|

Note:  not all matrix forms are supported by the Nastran OP4 Format.  This is not a limitation of the OP4 reader.



# Pure Python OP4 Reader (v0.5) #
The pure Python OP4 reader is slightly more robust and can read OP4s that are not formatted perfectly.  It can read ASCII or binary files of either big or little endian.  You may also specify the precision of the matrices that will be created.

Import the OP4 reader
```
from pyNastran.op4.op4 import OP4
```

Load the OP4 reader and print the matrix names
```
op4 = OP4()
matrices = op4.readOP4('sol103.op3')
matrices.keys()

>>> ['A','B','C','D']
```

Get matrix B.  A matrix will be a numpy.ndarray (dense matrix) or a scipy.coo\_matrix (sparse matrix).
```
(formatB,B) = matrices['B']
print B
```

You can also limit the matrices by specifying which ones you want
```
op4 = OP4()
matrices = op4.readOP4('sol103.op3',matrixNames=['A','B'])
matrices.keys()

>>> ['A','B']
```
Another way to limit the matrices
```
op4 = OP4()
matrices = op4.readOP4('sol103.op3',matrixNames='D')
matrices.keys()

>>> ['D']
```

You may also specify the precision of the matrix regardless of the data in the OP4.  The options are:  'single', 'double', and 'default'.  Double precision matrices take up twice as much memory, but are more accurate and therefore recommended.
```
op4 = OP4()
matrices = op4.readOP4('sol103.op3',matrixNames='D',precision='single')
matrices.keys()

>>> ['D']
```

# Cython OP4 Reader (v0.6) #
The Cython OP4 reader is faster and uses less memory than the pure Python version, but has a slightly different API.

Import the OP4 reader
```
import pyNastran.op4.cop4 as op4
```

## Examine OP4 File Contents ##
Print detailed information on the matrices stored in the file sol103.op4.
```
op4.OP4('sol103.op4').print_header()
```
Store the OP4 header information in the Python variable 'header'.
```
header = op4.OP4('sol103.op4')
```
Print the names of the matrices in the file 'sol103.op4' using the header variable created on the previous line:
```
for i in range(header.nmat): print('%2d. %s' % (i+1, header.name[i]))
```

## Load Matrices ##
Load the first matrix from the file kd.op4:
```
K = op4.Load('kd.op4')
```

Load the matrices named Kxx and Mxx from kd.op4:
```
K, M = op4.Load('kd.op4', 'Kxx', 'Mxx')
```

Load the first three matrices from file kd.op4:
```
K, D, M = op4.Load('kd.op4', nmat=3)
```

Skip the first matrix, then load the next two matrices from kd.op4:
```
D, M = op4.Load('kd.op4', nmat=2, skip=1)
```

Load all matrices from kd.op4.  Returns a dictionary of matrices keyed by matrix name:
```
Dict = op4.Load('kd.op4', '*')
```

Skip the first matrix, then load all remaining matrices from kd.op4.  Returns a dictionary of matrices keyed by matrix name:
```
Dict = op4.Load('kd.op4', '*', skip=1)
```

## Save Matrices ##
Write numpy arrays K and D to the file 'kd.op4' with OP4 matrix names "Kxx" and "M0000001" using native binary format (matrices that are not given explicit names, as with D in the line below, receive default names of "M0000001" through "M9999999"):
```
op4.Save('kd.op4', Kxx=K, D)
```

Write numpy arrays K and D to the file 'kd.op4' with names "Kxx" and "Dxx" as a text file with 5 mantissa digits:
```
op4.Save('kd.op4', Kxx=K, Dxx=D, digits=5)
```


# Cython OP4 Reader (v0.4, v0.5) #
The Cython OP4 reader is faster than the pure Python version, but has a slightly different API.  It requires that you know the order of your matrices in the OP4.

Import the OP4 reader
```
from pyNastran.op4.op4 import OP4  # v0.4
from pyNastran.op4.cop4 import OP4  # v0.5
```

Creates a read object associated with the file sol103.op4.
```
fh = OP4('sol103.op4', 'r')
```

Returns the number of matrices contained in the file.
```
fh.nmat
```

Prints detailed information on the matrices stored in the file.
```
fh.print_header()
```

Loads the second and third matrices from the file sol103.op4 into K and M.
```
K, M = fh.Load(skip=1, nmat=2)
```


## Additional COP4 Limitations ##
These will be resolved in future releases:
  1. Unable to byte-swap files created on a different endian machine