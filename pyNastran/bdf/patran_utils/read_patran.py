from __future__ import print_function
import numpy as np


def read_patran(patran_filename, fdtype='float64', idtype='int32'):
    """reads a patran .nod formatted file"""
    with open(patran_filename, 'r') as patran_file:
        lines = patran_file.readlines()
    lines = [line.rstrip() for line in lines if line.rstrip()]

    sline = lines[1].split()
    nnodes = int(sline[0])
    fmt = sline[2]
    nvalues = int(sline[4])

    if 'e' in fmt or 'E' in fmt or '.' in fmt:
        assert '-' not in fmt, 'fmt=%r' % fmt
        dtype = fdtype
        width = len(fmt) + 1
    else:
        dtype = idtype
        width = len(fmt) + 1

    print('fmt=%r; width=%s' % (fmt, width))
    nids = []

    #line0 = lines[0]
    #nid = line0[:8]
    #data

    data = []
    for line in lines[2:]:
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
    data_array = np.array(data, dtype=fdtype)
    return nids_array, data_array

if __name__ == '__main__':
    read_patran('normals.nod', fdtype='float64', idtype='int32')
