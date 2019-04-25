"""
defines:
 - version_latest, version_current, is_newer = check_for_newer_version(version_current=None)
 - res_obj, title = create_res_obj(
       islot, headers, header, A, fmt_dict, result_type,
       is_deflection=False, is_force=False,
       dim_max=None, xyz_cid0=None, colormap='jet')
 - B, fmt_dict_without_index, names_without_index = load_deflection_csv(
       out_filename, encoding='latin1')
 - A, fmt_dict, names = load_csv(out_filename, encoding='latin1')
 - grid_ids, xyz, bars, tris, quads = load_user_geom(fname, log=None, encoding='latin1')
"""
from __future__ import print_function
import os
#import re
import sys
import traceback
from codecs import open
from six.moves import urllib

import numpy as np
import pyNastran
from pyNastran.utils import _filename
from pyNastran.utils.numpy_utils import loadtxt_nice
from pyNastran.gui.gui_objects.gui_result import GuiResult
from pyNastran.converters.nastran.displacements import DisplacementResults, ForceTableResults
from pyNastran.converters.stl.stl import read_stl


def check_for_newer_version(version_current=None):
    """
    Checks to see if a newer version of pyNastran has been released.
    Only checks this for the GUI.

    Looks for:
        ## pyNastran v0.7.2 has been Released (4/25/2015)

    Specifically, it finds, 'has been released'
       then takes the the part that:
         - starts with 'v',
         - strips the 'v'
         - makes a version tuple:
           - (0,7,2)
       and compares that to the current version
    """
    is_newer = False
    if version_current is None:
        version_current = pyNastran.__version__
    target_url = 'https://raw.githubusercontent.com/SteveDoyle2/pyNastran/master/README.md'
    try:
        # it's a file like object and works just like a file
        data = urllib.request.urlopen(target_url)
    except (urllib.error.HTTPError, urllib.error.URLError):
    #except: #  urllib2.URLError
        #print(help(urllib))
        #raise
        return None, None, is_newer

    for btye_line in data: # files are iterable
        line_lower = btye_line.lower().decode('utf-8')
        if 'has been released' in line_lower:
            sline = line_lower.split()
            version_latest = [slot for slot in sline if slot.startswith('v')][0][1:]
            break

    is_dev = False
    if 'dev' in version_current:
        is_dev = True

    try:
        major, minor, rev = version_current.split('+')[0].split('.')
    except ValueError:
        print('sline = %s' % sline)
        print('version_current = %s' % version_current)
        raise
    major = int(major)
    minor = int(minor)
    rev = int(rev)
    tuple_current_version = (major, minor, rev)

    try:
        major, minor, rev = version_latest.split('_')[0].split('.')
    except ValueError:
        print('sline = %s' % sline)
        print('version_latest = %s' % version_latest)
        raise

    major = int(major)
    minor = int(minor)
    rev = int(rev)
    tuple_latest_version = (major, minor, rev)
    #print('tuple_latest_version = %s' % str(tuple_latest_version))  # (0,7,2)
    #print('tuple_current_version = %s' % str(tuple_current_version))  # (0,8,0)

    if (tuple_current_version < tuple_latest_version or
            (is_dev and tuple_current_version == tuple_latest_version)):
        print('pyNastran %s is now availible; current=%s' % (version_latest, version_current))
        is_newer = True
    return version_latest, version_current, is_newer

def create_res_obj(islot, headers, header, A, fmt_dict, result_type,
                   is_deflection=False, is_force=False,
                   dim_max=None, xyz_cid0=None, colormap='jet'):
    """
    Parameters
    ----------
    islot : int
        ???
    A : dict[key] = (n, m) array
        the numpy arrays
        key : str
            the name
        n : int
            number of nodes/elements
        m : int
            secondary dimension
            N/A : 1D array
            3 : deflection
    headers : List[str]???
        the titles???
    fmt_dict : dict[header] = fmt
        the format of the arrays
        header : str
            the name
        fmt : str
            '%i', '%f'
    result_type : str
        'node', 'centroid'
    dim_max : float
        required for forces/displacements
    xyz_cid0 : (nnodes, 3)
        the points
    """
    #print('create_res_object, header=%r' % header)
    datai = A[header]
    fmti = fmt_dict[header]
    title = header
    location = result_type

    dimension = len(datai.shape)
    if dimension == 1:
        vector_size = 1
    elif dimension == 2:
        vector_size = datai.shape[1]
    else:  # pramga: no cover
        raise RuntimeError('dimension=%s' % (dimension))

    if vector_size == 1:
        res_obj = GuiResult(
            islot, header, title, location, datai,
            nlabels=None, labelsize=None, ncolors=None,
            colormap=colormap, data_format=fmti,
        )
    elif vector_size == 3:
        # title is 3 values
        # header is 3 values
        # scale is 3 values
        titles = [header]
        headers = header

        norm_max = np.linalg.norm(datai, axis=1).max()
        scales = [dim_max / norm_max * 0.25]
        data_formats = [fmti] * 3
        scalar = None
        dxyz = datai
        if is_deflection:
            xyz = xyz_cid0
            res_obj = DisplacementResults(
                #subcase_id, titles, headers, xyz, dxyz, unused_scalar
                islot, titles, headers,
                xyz, dxyz, scalar, scales, data_formats=data_formats,
                nlabels=None, labelsize=None, ncolors=None,
                colormap=colormap,
                set_max_min=True)
        elif is_force:
            # scalar is unused
            res_obj = ForceTableResults(
                islot, titles, headers,
                dxyz, scalar, scales, data_formats=data_formats,
                nlabels=None, labelsize=None, ncolors=None,
                colormap=colormap,
                set_max_min=True)
        else:  # pramga: no cover
            raise RuntimeError('is_deflection=%s is_force=%s' % (is_deflection, is_force))
    else:  # pramga: no cover
        raise RuntimeError('vector_size=%s' % (vector_size))
    return res_obj, title

def load_deflection_csv(out_filename, encoding='latin1'):
    """
    The GUI deflection CSV loading function.

    Considers:
      - extension in determining how to load a file (e.g. commas or not)
      - header line of file for information regarding data types
    """
    ext = os.path.splitext(out_filename)[1].lower()
    if ext not in ['.csv', '.dat', '.txt']:
        raise NotImplementedError('extension=%r is not supported (use .dat, .txt, or .csv)' % ext)

    with open(_filename(out_filename), 'r', encoding=encoding) as file_obj:
        names, fmt_dict, dtype, delimiter = _load_format_header(file_obj, ext, force_float=False)

        try:
            #A = np.loadtxt(file_obj, dtype=dtype, delimiter=delimiter)
            A = loadtxt_nice(file_obj, comments='#', delimiter=delimiter)
        except:
            traceback.print_exc(file=sys.stdout)
            msg = 'extension=%r nheaders=%s delimiter=%r dtype=%s' % (
                ext, len(names), delimiter, dtype)
            raise RuntimeError(msg)

    names_without_index = names[1:]
    fmt_dict_without_index = {key:fmt_dict[key] for key in names_without_index}

    nnames_without_index = len(names_without_index)
    nexpected_results = 1 + 3 * nnames_without_index

    try:
        _nrows, ncols = A.shape
    except ValueError:
        msg = ('A should be (nnodes, 1+ndeflection_results); '
               'A.shape=%s nexpected_results=%s names=%s' % (
                   str(A.shape), nexpected_results, names))
        raise ValueError(msg)

    if ncols != nexpected_results:
        msg = 'A.shape=%s ncols=%s nexpected_results=%s names=%s nnames_without_index=%s' % (
            str(A.shape), ncols, nexpected_results, names, nnames_without_index)
        raise ValueError(msg)

    B = {}
    for i, name in enumerate(names_without_index):
        B[name] = A[:, 1+3*i:1+3*i+3]

    assert len(B) == len(fmt_dict_without_index), 'B.keys()=%s fmt_dict.keys()=%s' % (list(B.keys()), list(fmt_dict_without_index.keys()))
    assert len(B) == len(names_without_index), 'B.keys()=%s names.keys()=%s' % (list(B.keys()), names_without_index)
    return B, fmt_dict_without_index, names_without_index

def load_csv(out_filename, encoding='latin1'):
    """
    The GUI CSV loading function.

    Considers:
      - extension in determining how to load a file (e.g. commas or not)
      - header line of file for information regarding data types
    """
    ext = os.path.splitext(out_filename)[1].lower()
    if ext not in ['.csv', '.dat', '.txt']:
        raise NotImplementedError('extension=%r is not supported (use .dat, .txt, or .csv)' % ext)

    with open(_filename(out_filename), 'r', encoding=encoding) as file_obj:
        names, fmt_dict, dtype, delimiter = _load_format_header(file_obj, ext, force_float=False)

        try:
            #A = loadtxt(file_obj, dtype=dtype, delimiter=delimiter)
            A = loadtxt_nice(file_obj, dtype=dtype, comments='#', delimiter=delimiter)
        except:
            traceback.print_exc(file=sys.stdout)
            msg = 'extension=%r nheaders=%s delimiter=%r dtype=%s' % (
                ext, len(names), delimiter, dtype)
            raise RuntimeError(msg)
    return A, fmt_dict, names

def _load_format_header(file_obj, ext, force_float=False):
    header_line = file_obj.readline().strip()
    if not header_line.startswith('#'):
        msg = 'Expected file of the form:\n'
        if ext in ['.dat', '.txt']:
            msg += '# nodeid var1 var2\n'
            msg += '1 1 2\n'
            msg += '2 3 4\n'
            msg += '\nor:\n'
            msg += '# nodeid, var1(%i) var2(%f)\n'
            msg += '1 1 2.1\n'
            msg += '2 3 4.2\n'
        elif ext == '.csv':
            msg += '# nodeid, var1, var2\n'
            msg += '1, 1, 2\n'
            msg += '2, 3, 4\n'
            msg += '\nor:\n'
            msg += '# nodeid, var1(%i), var2(%f)\n'
            msg += '1, 1, 2.1\n'
            msg += '2, 3, 4.2\n'
        else:
            msg = 'extension=%r is not supported (use .dat, .txt, or .csv)' % ext
            raise NotImplementedError(msg)
        raise SyntaxError(msg)

    header_line = header_line.lstrip('# \t').strip()
    if ext in ['.dat', '.txt']:
        headers = header_line.split()
    elif ext == '.csv':
        headers = header_line.split(',')
    else:
        msg = 'extension=%r is not supported (use .dat, .txt, or .csv)' % ext
        raise NotImplementedError(msg)
    headers = [header.strip() for header in headers if header.strip()]

    fmt_dict = {}
    names = []
    dtype_fmts = []
    for iheader, header in enumerate(headers):
        # TODO: works for making a \n, but screws up the sidebar
        #       and scale
        header2 = header.strip()#.replace('\\n', '\n')
        dtype_fmt = 'float'

        str_fmt = '%.3f'
        header2_temp = None
        if iheader == 0:
            # the first column must be an integer
            if '(' in header2:
                header2_temp, fmt_temp = header2[:-1].rsplit('(', 1)
                header2_temp = header2_temp.strip()
                fmt = fmt_temp.strip()
                assert 'i' in fmt, 'header=%r must be an integer; fmt=%r' % (header2, fmt)
            header2_temp = 'index'
            dtype_fmt = 'int32'
            str_fmt = '%i'

        elif header2.endswith(')') and '%' in header2:
            header2_temp, fmt_temp = header2[:-1].rsplit('(', 1)
            header2_temp = header2_temp.strip()
            fmt = fmt_temp.strip()

            #('S1', 'i4', 'f4')
            if '%' in fmt:
                #fmt_temp = fmt_temp.replace('%', '%%')
                if force_float and 'i' in fmt:
                    fmt % 5
                    dtype_fmt = 'float32'
                    str_fmt = '%.0f'
                elif 'i' in fmt:
                    fmt % 5
                    dtype_fmt = 'int32'
                    str_fmt = '%i'

                elif 'g' in fmt or 'e' in fmt or 'f' in fmt or 's' in fmt:
                    dtype_fmt = 'float32'
                    str_fmt = fmt
                else:
                    raise TypeError('iheader=%s header=%r fmt=%r' % (iheader, header2, fmt))
            else:
                # default to float32
                dtype_fmt = 'float32'
        else:
            dtype_fmt = 'float32'
            header2_temp = header2

        assert header2_temp is not None
        names.append(header2_temp)
        dtype_fmts.append(dtype_fmt)
        fmt_dict[header2_temp] = str_fmt

    if ext in ['.dat', '.txt']:
        delimiter = None
    elif ext == '.csv':
        delimiter = ','
    else:
        raise NotImplementedError('extension=%r is not supported (use .dat, .txt, or .csv)' % ext)

    dtype = {
        'names': tuple(names),
        'formats': tuple(dtype_fmts),
    }
    return names, fmt_dict, dtype, delimiter

def load_user_geom(fname, log=None, encoding='latin1'):
    """
    Loads a file of the form:

    # all supported cards
    #  - GRID
    #  - BAR
    #  - TRI
    #  - QUAD
    #
    # doesn't support:
    #  - solid elements
    #  - element properties
    #  - custom colors
    #  - coordinate systems
    #  - materials
    #  - loads
    #  - results

    #    id  x    y    z
    GRID, 1, 0.2, 0.3, 0.3
    GRID, 2, 1.2, 0.3, 0.3
    GRID, 3, 2.2, 0.3, 0.3
    GRID, 4, 5.2, 0.3, 0.3
    grid, 5, 5.2, 1.3, 2.3  # case insensitive

    #    ID, nodes
    BAR,  1, 1, 2

    #   eid, n1,n2,n3,n4
    TRI,  2, 1, 2, 3
    # this is a comment

    #   eid, n1,n2,n3,n4
    QUAD, 3, 1, 5, 3, 4
    QUAD, 4, 1, 2, 3, 4  # this is after a blank line
    """
    if fname.lower().endswith('.stl'):
        stl_filename = fname
        stl = read_stl(stl_filename, log=log)
        nnodes = stl.nodes.shape[0]
        ntris = stl.elements.shape[0]
        grid_ids = np.arange(1, nnodes+1, dtype='int32')
        xyz = stl.nodes
        eids = np.arange(1, ntris+1, dtype='int32')
        tris = np.vstack([eids, stl.elements.T + 1]).T
        #tris = stl.elements + 1
        #print(tris)
        quads = np.array([], dtype='int32')
        bars = np.array([], dtype='int32')
        return grid_ids, xyz, bars, tris, quads

    with open(_filename(fname), 'r', encoding=encoding) as user_geom:
        lines = user_geom.readlines()

    grid_ids = []
    xyz = []
    bars = []
    tris = []
    quads = []
    #lines2 = []
    for line in lines:
        line2 = line.strip().split('#')[0].upper()
        if line2:
            sline = line2.split(',')
            if line2.startswith('GRID'):
                assert len(sline) == 5, sline
                grid_ids.append(sline[1])
                xyz.append(sline[2:])
            elif line2.startswith('BAR'):
                assert len(sline) == 4, sline
                bars.append(sline[1:])
            elif line2.startswith('TRI'):
                assert len(sline) == 5, sline
                tris.append(sline[1:])
            elif line2.startswith('QUAD'):
                assert len(sline) == 6, sline
                quads.append(sline[1:])
            else:
                log.warning(sline)

    grid_ids = np.array(grid_ids, dtype='int32')
    xyz = np.array(xyz, dtype='float32')
    tris = np.array(tris, dtype='int32')
    quads = np.array(quads, dtype='int32')
    bars = np.array(bars, dtype='int32')
    return grid_ids, xyz, bars, tris, quads

#def natural_sort(l):
    #convert = lambda text: int(text) if text.isdigit() else text.lower()
    #alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
    #return sorted(l, key=alphanum_key)

#def num_sort():
    #"""
    #in_values  = ['anti_main', 'main', 'Family 4', 'Family 3', 'Family 1',
                  #'Patch 119', 'Patch 118', 'Patch 19', 'Patch 18']
    #out_values = ['anti_main', 'main', 'Family 1', 'Family 3', 'Family 4',
                  #'Patch 118', 'Patch 119', 'Patch 18', 'Patch 19']

    #'Patch 19 cat 20' not handled
    #"""
    #convert = lambda text: int(text) if text.isdigit() else text.lower()
    #alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]

if __name__ == '__main__':  # pragma: no cover
    check_for_newer_version()
