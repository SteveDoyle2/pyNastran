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

Load matrices from NASTRAN Output 4 (.op4) files.
Example:

  >>> from pyNastran.op4.op4 import OP4

    Load the OP4 reader
  
  >>> fh = OP4('sol103.op4', 'r')

    Creates a read object associated with the file sol103.op4.

  >>> fh.nmat

    Returns the number of matrices contained in the file.

  >>> fh.print_header()

    Prints detailed information on the matrices stored in the file.

  >>> K, M = fh.Load(skip=1, nmat=2)

    Loads the second and third matrices from the file sol103.op4 into K and M.

Limitations (these will be resolved in future releases):
    1. only handles dense matrices in text and binary op4 files
    2. unable to byte-swap files created on a different-endian machine
    3. unable to save Python arrays to new op4 files
