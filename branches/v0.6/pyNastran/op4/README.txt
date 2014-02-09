# Copyright (C) 1999-2012  Al Danial <al.danial@gmail.com>
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software Foundation,
# Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

Load and save matrices from/to NASTRAN Output 4 (.op4) files.

  >>> from pyNastran.op4 import cop4

    Examining OP4 File Contents
    ===========================

  >>> cop4.Scan('sol103.op4').print_header()

    Prints detailed information on the matrices stored in the file sol103.op4.

  >>> header = cop4.Scan('sol103.op4')

    Store the OP4 header information in the Python variable 'header'.

  >>> for i in range(header.nmat): print('%2d. %s' % (i+1, header.name[i]))

    Print the names of the matrices in the file 'sol103.op4'.


    Loading Matrices
    ================

  >>> K = cop4.Load('kd.op4')
        Load the first matrix from the file kd.op4.

  >>> K, D, M = cop4.Load('kd.op4', nmat=3)
        Load the first three matrices from file kd.op4.

  >>> D, M = cop4.Load('kd.op4', nmat=2, skip=1)
        Skip the first matrix, then load the next two matrices from kd.op4.

  >>> K, M = cop4.Load('kd.op4', 'Kxx', 'Mxx')
        Load the matrices named Kxx and Mxx from kd.op4.

  >>> Dict = cop4.Load('kd.op4', '*')
        Load all matrices from kd.op4.  Returns a dictionary of
        matrices keyed by matrix name.

  >>> Dict = cop4.Load('kd.op4', '*', skip=1)
        Skip the first matrix, then load all remaining matrices from 
        kd.op4.  Returns a dictionary of matrices keyed by matrix name.


    Saving Matrices
    ===============
  >>> cop4.Save('kd.op4', Kxx=K, D)
        Write numpy arrays K and D to the file 'kd.op4' with names 
        "Kxx" and "M0000001" using native binary format.
     
  >>> cop4.Save('kd.op4', Kxx=K, Dxx=D, digits=5)
        Write numpy arrays K and D to the file 'kd.op4' with names 
        "Kxx" and "Dxx" as a text file with 5 mantissa digits.


Limitation (to be resolved in a future release):
    - unable to byte-swap files created on a different-endian machine
