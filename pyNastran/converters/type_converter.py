"""
Multi-input/output format converter
"""
from __future__ import print_function
import glob
from docopt import docopt

import pyNastran
from pyNastran.bdf.bdf import BDF
from pyNastran.converters.tecplot.tecplot import Tecplot
from pyNastran.converters.tecplot.utils import merge_tecplot_files
# from pyNastran.converters.stl.stl import STL
from pyNastran.converters.stl.utils import merge_stl_files
# from pyNastran.converters.cart3d.cart3d_reader import Cart3D
from pyNastran.converters.cart3d.cart3d import Cart3D

from pyNastran.converters.nastran.nastran_to_cart3d import nastran_to_cart3d
from pyNastran.converters.nastran.nastran_to_stl import nastran_to_stl
from pyNastran.converters.nastran.nastran_to_tecplot import nastran_to_tecplot
from pyNastran.converters.nastran.nastran_to_ugrid import nastran_to_ugrid

from pyNastran.converters.stl.stl_to_nastran import stl_to_nastran_filename
#from pyNastran.converters.stl.stl_to_cart3d import stl_to_cart3d, stl_to_cart3d_filename
from pyNastran.converters.cart3d.cart3d_to_nastran import cart3d_to_nastran_filename
from pyNastran.converters.cart3d.cart3d_to_stl import cart3d_to_stl_filename
from pyNastran.converters.ugrid.ugrid_reader import UGRID
from pyNastran.converters.ugrid.ugrid3d_to_tecplot import ugrid_to_tecplot
from pyNastran.converters.tecplot.tecplot_to_nastran import tecplot_to_nastran_filename


def process_nastran(bdf_filename, fmt2, fname2, data=None, debug=True):
    """
    Converts Nastran to STL/Cart3d/Tecplot
    """
    assert fmt2 in ['stl', 'cart3d', 'tecplot', 'ugrid'], 'format2=%s' % fmt2
    xref = True
    if fmt2 == 'ugrid':
        xref = False
    model = BDF(debug=debug)
    model.read_bdf(bdf_filename, xref=xref)

    if fmt2 == 'stl':
        nastran_to_stl(model, fname2, is_binary=data['--binary'])
    elif fmt2 == 'cart3d':
        cart3d = nastran_to_cart3d(model)
        cart3d.write_cart3d(fname2)
    elif fmt2 == 'tecplot':
        tecplot = nastran_to_tecplot(model)
        tecplot_filename = fname2
        tecplot.write_tecplot(tecplot_filename, adjust_nids=False)
    elif fmt2 == 'ugrid':
        ugrid = nastran_to_ugrid(model, fname2)
    else:
        raise NotImplementedError(fmt2)


def process_cart3d(cart3d_filename, fmt2, fname2, data):
    """
    Converts Cart3d to STL/Nastran
    """
    assert fmt2 in ['stl', 'nastran', 'tecplot'], 'format2=%s' % fmt2
    # model = Cart3D()
    # model.read_cart3d(cart3d_filename, fname2)
    if fmt2 == 'stl':
        cart3d_to_stl_filename(cart3d_filename, fname2, is_binary=data['--binary'])
    elif fmt2 == 'nastran':
        cart3d_to_nastran_filename(cart3d_filename, fname2)
    elif fmt2 == 'tecplot':
        cart3d_to_tecplot(cart3d_filename, fname2)
    # elif fmt2 == 'ugrid':
        # cart3d_to_ugrid(model, fname2)
    else:
        raise NotImplementedError(fmt2)

def cart3d_to_tecplot(cart3d_filename, tecplot_filename):
    """
    Converts Cart3d to Tecplot
    """
    model = Cart3D()
    model.read_cart3d(cart3d_filename)
    tecplot = Tecplot()
    tecplot.xyz = model.points
    tecplot.tri_elements = model.elements
    tecplot.write_tecplot(tecplot_filename, adjust_nids=False)
    return tecplot

def process_stl(stl_filename, fmt2, fname2, data=None):
    """
    Converts STL to Nastran/Cart3d
    """
    assert fmt2 in ['stl', 'nastran', 'cart3d'], 'format2=%s' % fmt2
    if '*' in stl_filename:
        stl_filenames = glob.glob(stl_filename)
    else:
        stl_filenames = [stl_filename]
    assert len(stl_filenames) > 0, stl_filenames
    model = merge_stl_files(stl_filenames, stl_out_filename=None)
    scale = data['--scale']
    if scale is not None:
        assert isinstance(scale, float), 'scale=%r type=%r' % (scale, type(scale))
        model.nodes *= scale

    # model = STL()
    # model.read_stl(stl_filename)
    if fmt2 == 'nastran':
        stl_temp_filename = '__temp__.stl'
        model.write_stl(stl_out_filename, is_binary=True, #float_fmt='%6.12f',
                        stop_on_failure=False)
        stl_to_nastran_filename(stl_temp_filename, fname2)
        os.remove(stl_temp_filename)
        # stl_to_cart3d(model, fname2 + '.cart')
        # process_cart3d(fname2 + '.cart', 'nastran', fname2)
        # stl_to_nastran(model, fname2)
    #elif fmt2 == 'cart3d':
        # we don't have an STL -> Cart3d, so we:
        #    - STL -> BDF
        #    - BDF -> Cart3D
        # stl_to_cart3d(model, fname2)
        #stl_to_nastran_filename(stl_filename, fname2 + '.bdf')
        #stl_to_cart3d_filename(fname2 + '.bdf', fname2)
    elif fmt2 == 'stl':
        is_binary = data['--binary']
        model.write_stl(fname2, is_binary=is_binary, float_fmt='%6.12f', stop_on_failure=False)
    # elif fmt2 == 'tecplot':
        # stl_to_tecplot(model, fname2)
    # elif fmt2 == 'ugrid':
        # stl_to_ugrid(model, fname2)
    else:
        raise NotImplementedError(fmt2)


def element_slice(tecplot, data):
    xslice = data['--xx']
    yslice = data['--yy']
    zslice = data['--zz']
    # if xslice is not None:
        # xslice = data['--xx']
        # tecplot.slice_x(xslice)
    # if yslice is not None:
        # yslice = data['--yy']
        # tecplot.slice_y(yslice)
    # if zslice is not None:
        # zslice = data['--zz']
        # tecplot.slice_z(zslice)
    tecplot.slice_xyz(xslice, yslice, zslice)

def process_tecplot(tecplot_filename, fmt2, fname2, data=None):
    """
    Converts Tecplot to Tecplot

    Globs all input tecplot files (e.g. tecplot*.plt)
    """
    assert fmt2 in ['stl', 'nastran', 'cart3d', 'tecplot'], 'format2=%s' % fmt2
    if '*' in tecplot_filename:
        tecplot_filenames = glob.glob(tecplot_filename)
    else:
        tecplot_filenames = [tecplot_filename]
    assert len(tecplot_filenames) > 0, tecplot_filename
    model = merge_tecplot_files(tecplot_filenames, tecplot_filename_out=None)
    #if fmt2 == 'cart3d':
        #tecplot_to_cart3d(model, fname2)
    #elif fmt2 == 'stl':
        #tecplot_to_stl(model, fname2)
    # elif fmt2 == 'ugrid':
        # tecplot_to_ugrid(model, fname2)
    res_types = data['RESTYPE']
    is_points = not data['--block']
    if fmt2 == 'tecplot':
        print(data)
        element_slice(model, data)

        # this is a good way to merge files
        model.write_tecplot(fname2, res_types=res_types, is_points=is_points)
    elif fmt2 == 'nastran':
        tecplot_to_nastran_filename(model, fname2)
    elif fmt2 == 'stl':
        tecplot_to_nastran_filename(model, fname2 + '.bdf')
        process_nastran(fname2 + '.bdf', fmt2, fname2, data=data)
    elif fmt2 == 'cart3d':
        tecplot_to_nastran_filename(model, fname2 + '.bdf')
        process_nastran(fname2 + '.bdf', fmt2, fname2, data=data)
    else:
        raise NotImplementedError(fmt2)

def process_ugrid(ugrid_filename, fmt2, fname2, data=None):
    """
    Converts UGRID to Nastran/Cart3d/STL/Tecplot
    """
    assert fmt2 in ['stl', 'nastran', 'cart3d', 'tecplot'], 'format2=%s' % fmt2
    model = UGRID()
    model.read_ugrid(ugrid_filename)
    if fmt2 == 'nastran':
        # ugrid_to_nastran(model, fname2
        include_shells = True
        include_solids = True
        bdf_filename = fname2
        model.write_bdf(bdf_filename, include_shells=include_shells, include_solids=include_solids)
    elif fmt2 == 'cart3d':
        include_shells = False
        include_solids = True
        bdf_filename = fname2 + '.bdf'
        model.write_bdf(bdf_filename, include_shells=include_shells, include_solids=include_solids)
        # ugrid_to_cart3d(model, fname2)
        process_nastran(bdf_filename, 'cart3d', fname2, data=None)
    elif fmt2 == 'stl':
        include_shells = False
        include_solids = True
        bdf_filename = fname2 + '.bdf'
        model.write_bdf(bdf_filename, include_shells=include_shells, include_solids=include_solids)
        process_nastran(bdf_filename, 'cart3d', fname2, data=None)
        # ugrid_to_stl(model, fname2)
    elif fmt2 == 'tecplot':
        # ugrid_to_tecplot(model, fname2)
        tecplot = ugrid_to_tecplot(model)
        element_slice(tecplot, data)
        tecplot_filename = fname2
        tecplot.write_tecplot(tecplot_filename)
    else:
        raise NotImplementedError(fmt2)

def run(fmt1, fname1, fmt2, fname2, data):
    """
    Runs the format converter
    """
    if fmt1 == 'nastran':
        process_nastran(fname1, fmt2, fname2, data)
    elif fmt1 == 'cart3d':
        process_cart3d(fname1, fmt2, fname2, data)
    elif fmt1 == 'stl':
        process_stl(fname1, fmt2, fname2, data)
    elif fmt1 == 'tecplot':
        process_tecplot(fname1, fmt2, fname2, data)
    elif fmt1 == 'ugrid':
        process_ugrid(fname1, fmt2, fname2, data)
    else:
        raise NotImplementedError(fmt1)


def main():
    """Interface for format_converter"""
    msg = "Usage:\n"
    msg += "  format_converter nastran <INPUT> <format2> <OUTPUT> [-o <OP2>]\n"
    msg += "  format_converter <format1> <INPUT> tecplot <OUTPUT> [-r RESTYPE...] [-b] [--block] [-x <X>] [-y <Y>] [-z <Z>] [--scale SCALE]\n"
    msg += "  format_converter <format1> <INPUT> stl     <OUTPUT> [-b]  [--scale SCALE]\n"
    msg += "  format_converter <format1> <INPUT> <format2> <OUTPUT> [--scale SCALE]\n"
    #msg += "  format_converter nastran  <INPUT> <format2> <OUTPUT>\n"
    #msg += "  format_converter cart3d   <INPUT> <format2> <OUTPUT>\n"
    msg += '  format_converter -h | --help\n'
    msg += '  format_converter -v | --version\n'
    msg += "\n"
    msg += "Options:\n"
    msg += "  format1        format type (nastran, cart3d, stl, ugrid, tecplot)\n"
    msg += "  format2        format type (nastran, cart3d, stl, ugrid, tecplot)\n"
    msg += "  INPUT          path to input file\n"
    msg += "  OUTPUT         path to output file\n"
    msg += "  -o OP2, --op2 OP2  path to results file (nastran-specific)\n"
    msg += "                 only used for Tecplot (not supported)\n"
    msg += "  -x X, --xx X   Creates a constant x slice; keeps points < X\n"
    msg += "  -y Y, --yy Y   Creates a constant y slice; keeps points < Y\n"
    msg += "  -z Z, --zz Z   Creates a constant z slice; keeps points < Z\n"
    msg += "  --scale SCALE  Apply a scale factor to the XYZ locations (default=1.0)\n"
    msg += "  --block        Writes the data in BLOCK (vs. POINT) format\n"
    msg += "  -r, --results  Specifies the results to write to limit output\n"
    msg += "  -b, --binary   writes the STL in binary (not supported for Tecplot)\n"
    msg += "  -h, --help     show this help message and exit\n"
    msg += "  -v, --version  show program's version number and exit\n"
    msg += '\n'
    msg += 'Notes:\n'
    msg += "  Nastran->Tecplot assumes sequential nodes and consistent types (shell/solid)\n"
    msg += "  STL/Tecplot supports globbing as the input filename\n"
    msg += "  Tecplot slicing doesn't support multiple slice values and will give bad results (not crash)\n"
    msg += "  UGRID outfiles must be of the form model.b8.ugrid, where b8, b4, lb8, lb4 are valid choices and periods are important\n"
    msg += "  Scale has only been tested on STL -> STL\n"


    ver = str(pyNastran.__version__)
    data = docopt(msg, version=ver)

    is_nastran = data['nastran']
    format1 = data['<format1>']
    if is_nastran:
        format1 = 'nastran'
        data['<format1>'] = format1


    format2 = data['<format2>']
    is_stl = data['stl']
    if is_stl:
        format2 = 'stl'
        data['<format2>'] = format2

    is_tecplot = data['tecplot']
    if is_tecplot:
        format2 = 'tecplot'
        data['<format2>'] = format2
    print(data)

    if data['--scale']:
        data['--scale'] = float(data['--scale'])

    input_filename = data['<INPUT>']
    output_filename = data['<OUTPUT>']
    run(format1, input_filename, format2, output_filename, data)

if __name__ == '__main__':
    main()
