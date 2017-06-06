"""
defines:
 - nids_array, data_array = read_patran(patran_filename, fdtype='float64', idtype='int32')
"""
from __future__ import print_function
import os
from collections import defaultdict
import numpy as np
from pyNastran.utils import print_bad_path

def read_patran_format(patran_fmt_filename):
    """
    reads a file as shown below::

      /* mscnastran_op2_nod.res_tmpl */
      /* PATRAN 2.5 results file template for NASPAT OUTPUT2 .nod files */

      KEYLOC = 0

      TYPE = scalar
      COLUMN = 1
      PRI = Normals
      SEC = Normal (X)

      TYPE = scalar
      COLUMN = 2
      PRI = Normals
      SEC = Normal (Y)

      TYPE = scalar
      COLUMN = 3
      PRI = Normals
      SEC = Normal (Z)

      TYPE = END
    """
    assert os.path.exists(patran_fmt_filename), print_bad_path(patran_fmt_filename)
    with open(patran_fmt_filename, 'r') as patran_file:
        lines = patran_file.readlines()

    headers = defaultdict(list)
    for line in lines:
        if '=' in line:
            sline = line.strip().split('=')
            key = sline[0].strip()
            value = sline[1].strip()
            headers[key].append(value)
    return headers


def read_patran(patran_filename, fdtype='float64', idtype='int32'):
    """
    reads a patran .nod formatted file::

      KEY
         2214    0   0.000000E+00     0     3


         10.9126140E+00-.1825228E+00-.3658157E+00
         20.3790452E+00-.1844318E+00-.9068129E+00
         30.1772419E+00-.2389538E+00-.9547180E+00
         40.1056876E+00-.2347771E+00-.9662866E+00
         50.6002256E-01-.1678307E+00-.9839869E+00
    """
    base = os.path.splitext(patran_filename)[0]
    headers = read_patran_format(base + '.res_tmpl')

    assert os.path.exists(patran_filename), print_bad_path(patran_filename)
    with open(patran_filename, 'r') as patran_file:
        lines = patran_file.readlines()

    title = lines[0].strip()
    subtitle = (lines[2].strip() + ';' + lines[3].strip()).rstrip(';')

    sline = lines[1].strip().split()
    nnodes = int(sline[0])
    max_node = int(sline[1])
    fmt = sline[2]
    nvalues = int(sline[4])

    if 'e' in fmt or 'E' in fmt or '.' in fmt:
        assert '-' not in fmt, 'fmt=%r' % fmt
        dtype = fdtype
        width = len(fmt) + 1
    else:
        dtype = idtype
        width = len(fmt) + 1

    #print('fmt=%r; width=%s' % (fmt, width))
    nids = []

    #line0 = lines[0]
    #nid = line0[:8]
    #data

    data = []
    for line in lines[4:]:
        nid = line[:8]
        nids.append(nid)
        #print('nid = %r' % nid)
        i8 = 8
        i16 = i8 + width
        datai = []
        #print('i8=%s i16=%s' % (i8, i16))
        for ivalue in range(nvalues):
            value = line[i8:i16]
            i8 += width
            i16 += width
            datai.append(value)
        #print('datai = %r' % datai)
        #aaa
        data.append(datai)

    nids_array = np.array(nids, dtype=idtype)
    data_array = np.array(data, dtype=dtype)

    data_dict = {
        'title' : title,
        'subtitle' : subtitle,
        'nids' : nids_array,
        'data' : data_array,
        'headers' : headers,
    }
    return data_dict

if __name__ == '__main__':
    read_patran('normals.nod', fdtype='float64', idtype='int32')
