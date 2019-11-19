"""Multi-input/output format converter"""
import os
import sys
import glob

# stl_to_plot3d ???


def process_nastran(bdf_filename, fmt2, fname2, log, data=None, debug=True, quiet=False):
    """
    Converts Nastran to STL/Cart3d/Tecplot/UGRID3d
    """
    assert fmt2 in ['stl', 'cart3d', 'tecplot', 'ugrid', 'nastran', 'abaqus'], 'format2=%s' % fmt2
    from pyNastran.bdf.bdf import BDF
    xref = True
    if fmt2 == 'ugrid':
        xref = False
    model = BDF(log=log, debug=debug)
    model.read_bdf(bdf_filename, validate=False, xref=xref)

    if data['--scale'] != 1.0:
        scale = data['--scale']
        data['--scale'] = 1.0
        for node in model.nodes.values():
            node.xyz = node.get_position() * scale
            node.cp = 0
            del node.cp_ref

    if fmt2 == 'stl':
        from pyNastran.converters.nastran.nastran_to_stl import nastran_to_stl
        nastran_to_stl(model, fname2, is_binary=data['--binary'])
    elif fmt2 == 'cart3d':
        from pyNastran.converters.nastran.nastran_to_cart3d import nastran_to_cart3d
        cart3d = nastran_to_cart3d(model)
        cart3d.write_cart3d(fname2)
    elif fmt2 == 'tecplot':
        from pyNastran.converters.nastran.nastran_to_tecplot import nastran_to_tecplot
        tecplot = nastran_to_tecplot(model)
        tecplot_filename = fname2
        tecplot.write_tecplot(tecplot_filename, adjust_nids=False)
    elif fmt2 == 'ugrid':
        from pyNastran.converters.nastran.nastran_to_ugrid import nastran_to_ugrid
        nastran_to_ugrid(model, fname2)
    elif fmt2 == 'abaqus':
        from pyNastran.converters.nastran.nastran_to_abaqus import nastran_to_abaqus
        nastran_to_abaqus(model, fname2)
    elif fmt2 == 'nastran':
        model.write_bdf(fname2, size=16)
    else:
        raise NotImplementedError('fmt2=%s is not supported by process_nastran' % fmt2)


def process_cart3d(cart3d_filename, fmt2, fname2, log, data, quiet=False):
    """
    Converts Cart3d to STL/Nastran/Tecplot/Cart3d
    """
    assert fmt2 in ['stl', 'nastran', 'tecplot', 'cart3d'], 'format2=%s' % fmt2
    from pyNastran.converters.cart3d.cart3d import read_cart3d

    model = read_cart3d(cart3d_filename, log=log)
    if data['--scale'] != 1.0:
        model.points *= data['--scale']
        data['--scale'] = 1.0

    if fmt2 == 'stl':
        from pyNastran.converters.cart3d.cart3d_to_stl import cart3d_to_stl_filename
        cart3d_to_stl_filename(model, fname2, log=log, is_binary=data['--binary'])
    elif fmt2 == 'nastran':
        from pyNastran.converters.cart3d.cart3d_to_nastran import cart3d_to_nastran_filename
        cart3d_to_nastran_filename(model, fname2, log=log)
    elif fmt2 == 'tecplot':
        from pyNastran.converters.cart3d.cart3d_to_tecplot import cart3d_to_tecplot
        cart3d_to_tecplot(model, fname2, log=log)
    elif fmt2 == 'cart3d':
        model.write_cart3d(fname2, is_binary=data['--binary'])
    # elif fmt2 == 'ugrid':
        # cart3d_to_ugrid(model, fname2)
    else:
        raise NotImplementedError('fmt2=%s is not supported by process_cart3d' % fmt2)


def process_stl(stl_filename, fmt2, fname2, log, data=None, quiet=False):
    """
    Converts STL to Nastran/Cart3d
    """
    assert fmt2 in ['stl', 'nastran', 'cart3d'], 'format2=%s' % fmt2
    if '*' in stl_filename:
        stl_filenames = glob.glob(stl_filename)
    else:
        stl_filenames = [stl_filename]
    assert len(stl_filenames) > 0, stl_filenames
    from pyNastran.converters.stl.utils import merge_stl_files

    model = merge_stl_files(stl_filenames, stl_out_filename=None, log=log)
    scale = data['--scale']
    if scale is not None:
        assert isinstance(scale, float), 'scale=%r type=%r' % (scale, type(scale))
        model.nodes *= scale

    # model = STL()
    # model.read_stl(stl_filename)
    if fmt2 == 'nastran':
        from pyNastran.converters.stl.stl_to_nastran import stl_to_nastran
        stl_to_nastran(model, fname2, log=log)
    elif fmt2 == 'cart3d':
        from pyNastran.converters.stl.stl_to_cart3d import stl_to_cart3d
        stl_to_cart3d(model, fname2, log=log)
    elif fmt2 == 'stl':
        is_binary = data['--binary']
        model.write_stl(fname2, is_binary=is_binary, float_fmt='%6.12f', stop_on_failure=False)
    # elif fmt2 == 'tecplot':
        # stl_to_tecplot(model, fname2)
    # elif fmt2 == 'ugrid':
        # stl_to_ugrid(model, fname2)
    else:
        raise NotImplementedError('fmt2=%s is not supported by process_stl' % fmt2)


def element_slice(tecplot, data):
    """removes solid elements from a tecplot model"""
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
    #print(tecplot)
    tecplot.slice_xyz(xslice, yslice, zslice)


def process_tecplot(tecplot_filename, fmt2, fname2, log, data=None, quiet=False):
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
    from pyNastran.converters.tecplot.utils import merge_tecplot_files
    from pyNastran.converters.tecplot.tecplot_to_nastran import tecplot_to_nastran_filename
    from pyNastran.converters.tecplot.tecplot_to_cart3d import tecplot_to_cart3d_filename

    model = merge_tecplot_files(tecplot_filenames, tecplot_filename_out=None, log=log)
    #if fmt2 == 'cart3d':
        #tecplot_to_cart3d(model, fname2)
    #elif fmt2 == 'stl':
        #tecplot_to_stl(model, fname2)
    # elif fmt2 == 'ugrid':
        # tecplot_to_ugrid(model, fname2)
    res_types = data['RESTYPE']
    unused_is_points = not data['--block']
    if fmt2 == 'tecplot':
        if not quiet:  # pragma: no cover
            print(data)
        element_slice(model, data)

        # this is a good way to merge files
        model.write_tecplot(fname2, res_types=res_types) # is_points=is_points
    elif fmt2 == 'nastran':
        tecplot_to_nastran_filename(model, fname2)
    elif fmt2 == 'stl':
        cart3d_filename = fname2 + '.tri'
        tecplot_to_cart3d_filename(model, cart3d_filename, log=log)
        process_cart3d(cart3d_filename, fmt2, fname2, log, data=data, quiet=quiet)
        os.remove(cart3d_filename)
        #tecplot_to_nastran_filename(model, fname2 + '.bdf')
        #process_nastran(fname2 + '.bdf', fmt2, fname2, log, data=data, quiet=quiet)
    elif fmt2 == 'cart3d':
        # supports tris/quads, not loads
        #tecplot_to_nastran_filename(model, fname2 + '.bdf')
        #process_nastran(fname2 + '.bdf', fmt2, fname2, log, data=data, quiet=quiet)

        # supports quads/loads, not tris
        tecplot_to_cart3d_filename(model, fname2, log=log)
    else:
        raise NotImplementedError('fmt2=%s is not supported by process_tecplot' % fmt2)


def process_ugrid(ugrid_filename, fmt2, fname2, log, data=None, quiet=False):
    """
    Converts UGRID to Nastran/Cart3d/STL/Tecplot
    """
    assert fmt2 in ['stl', 'nastran', 'cart3d', 'tecplot'], 'format2=%s' % fmt2
    read_shells = True
    read_solids = True
    if fmt2 in ['stl', 'cart3d']:
        read_shells = True
        read_solids = False

    from pyNastran.converters.aflr.ugrid.ugrid_reader import UGRID
    model = UGRID(read_shells=read_shells, read_solids=read_solids, log=log)
    model.read_ugrid(ugrid_filename)
    if fmt2 == 'nastran':
        # ugrid_to_nastran(model, fname2
        include_shells = True
        include_solids = True
        bdf_filename = fname2
        model.write_bdf(bdf_filename, include_shells=include_shells, include_solids=include_solids)
    elif fmt2 == 'cart3d':
        include_shells = True
        include_solids = False
        bdf_filename = fname2 + '.bdf'
        model.write_bdf(bdf_filename, include_shells=include_shells, include_solids=include_solids)
        # ugrid_to_cart3d(model, fname2)
        process_nastran(bdf_filename, 'cart3d', fname2, data=None)
    elif fmt2 == 'stl':
        include_shells = True
        include_solids = False
        bdf_filename = fname2 + '.bdf'
        model.write_bdf(bdf_filename, include_shells=include_shells, include_solids=include_solids)
        process_nastran(bdf_filename, 'cart3d', fname2, data=None)
        # ugrid_to_stl(model, fname2)
    elif fmt2 == 'tecplot':
        from pyNastran.converters.aflr.ugrid.ugrid3d_to_tecplot import ugrid_to_tecplot
        # ugrid_to_tecplot(model, fname2)
        tecplot, unused_zone = ugrid_to_tecplot(model)
        element_slice(tecplot, data)
        tecplot_filename = fname2
        tecplot.write_tecplot(tecplot_filename)
    else:
        raise NotImplementedError('fmt2=%s is not supported by process_ugrid' % fmt2)


def run_format_converter(fmt1, fname1, fmt2, fname2, data, log, quiet=False):
    """
    Runs the format converter
    """
    if fmt1 == 'nastran':
        process_nastran(fname1, fmt2, fname2, log, data=data, quiet=quiet)
    elif fmt1 == 'cart3d':
        process_cart3d(fname1, fmt2, fname2, log, data=data, quiet=quiet)
    elif fmt1 == 'stl':
        process_stl(fname1, fmt2, fname2, log, data=data, quiet=quiet)
    elif fmt1 == 'tecplot':
        process_tecplot(fname1, fmt2, fname2, log, data=data, quiet=quiet)
    elif fmt1 == 'ugrid':
        process_ugrid(fname1, fmt2, fname2, log, data=data, quiet=quiet)
    elif fmt1 == 'vrml':
        process_vrml(fname1, fmt2, fname2, log, data=data, quiet=quiet)
    else:
        format1s = ['nastran', 'cart3d', 'stl', 'tecplot', 'ugrid', 'vrml']
        #format2s = ['nastran', 'cart3d', 'stl', 'ugrid', 'tecplot']
        raise NotImplementedError(f'fmt1={fmt1} is not supported by run; '
                                  f'use {", ".join(format1s)}')


def cmd_line_format_converter(argv=None, quiet=False):
    """Interface for format_converter"""
    if argv is None:
        argv = sys.argv
    msg = "Usage:\n"
    #format1s = ['nastran', 'cart3d', 'stl', 'ugrid', 'tecplot', 'vrml']
    #format2s = ['nastran', 'cart3d', 'stl', 'ugrid', 'tecplot']
    msg += "  format_converter nastran   <INPUT> <format2> <OUTPUT> [-o <OP2>] --no_xref\n"
    msg += "  format_converter <format1> <INPUT> tecplot   <OUTPUT> [-r RESTYPE...] [-b] [--block] [-x <X>] [-y <Y>] [-z <Z>] [--scale SCALE]\n"
    msg += "  format_converter <format1> <INPUT> stl       <OUTPUT> [-b]  [--scale SCALE]\n"
    msg += "  format_converter cart3d    <INPUT> <format2> <OUTPUT> [-b]  [--scale SCALE]\n"
    msg += "  format_converter <format1> <INPUT> <format2> <OUTPUT> [--scale SCALE]\n"
    #msg += "  format_converter nastran  <INPUT> <format2> <OUTPUT>\n"
    #msg += "  format_converter cart3d   <INPUT> <format2> <OUTPUT>\n"
    msg += '  format_converter -h | --help\n'
    msg += '  format_converter -v | --version\n'
    msg += "\n"
    msg += "Required Arguments:\n"
    msg += "  format1        format type (nastran, cart3d, stl, ugrid, tecplot, vrml)\n"
    msg += "  format2        format type (nastran, cart3d, stl, ugrid, tecplot, abaqus)\n"
    msg += "  INPUT          path to input file\n"
    msg += "  OUTPUT         path to output file\n"

    msg += "\n"
    msg += "Nastran Options:\n"
    msg += "  -o OP2, --op2 OP2  path to results file (nastran-specific)\n"
    msg += "                 only used for Tecplot (not supported)\n"
    msg += "  --no_xref      Don't cross-reference (nastran-specific)\n"

    msg += "\n"
    msg += "Tecplot Options:\n"
    msg += "  -x X, --xx X   Creates a constant x slice; keeps points < X\n"
    msg += "  -y Y, --yy Y   Creates a constant y slice; keeps points < Y\n"
    msg += "  -z Z, --zz Z   Creates a constant z slice; keeps points < Z\n"
    msg += "  --block        Writes the data in BLOCK (vs. POINT) format\n"
    msg += "  -r, --results  Specifies the results to write to limit output\n"

    msg += "\n"
    msg += "Tecplot/Cart3d/STL Options:\n"
    msg += "  --scale SCALE  Apply a scale factor to the XYZ locations (default=1.0)\n"
    msg += "  -b, --binary   writes the STL in binary (not supported for Tecplot)\n"

    msg += "\n"
    msg += "Info:\n"
    msg += "  -h, --help     show this help message and exit\n"
    msg += "  -v, --version  show program's version number and exit\n"
    msg += '\n'
    msg += 'Notes:\n'
    msg += "  Nastran->Tecplot assumes sequential nodes and consistent types (shell/solid)\n"
    msg += "  STL/Tecplot supports globbing as the input filename\n"
    msg += "  Tecplot slicing doesn't support multiple slice values and will give bad results (not crash)\n"
    msg += "  UGRID outfiles must be of the form model.b8.ugrid, where\n"
    msg += "    b8, b4, lb8, lb4 are valid choices and periods are important\n"
    msg += "  Scale has only been tested on STL -> STL\n"

    from docopt import docopt
    import pyNastran
    ver = str(pyNastran.__version__)
    data = docopt(msg, version=ver, argv=argv[1:])

    # because we have special blocks for tecplot/stl/cart3d
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

    is_cart3d = data['cart3d']
    if is_cart3d:
        format1 = 'cart3d'
        data['<format1>'] = format1

    # common options
    if data['--scale']:
        data['--scale'] = eval(data['--scale'])
    else:
        data['--scale'] = 1.0

    if not quiet:  # pragma: no cover
        print(data)
    input_filename = data['<INPUT>']
    output_filename = data['<OUTPUT>']
    level = 'warning' if  quiet else 'debug'
    from cpylog import SimpleLogger
    log = SimpleLogger(level=level)

    run_format_converter(format1, input_filename, format2, output_filename, data, log=log, quiet=quiet)


def process_vrml(vrml_filename, fmt2, fname2, log, data, quiet=False):
    """
    Converts VRML to Nastran
    """
    assert fmt2 in ['nastran', 'stl'], 'format2=%s' % fmt2
    #if data['--scale'] != 1.0:
        #model.points *= data['--scale']
        #data['--scale'] = 1.0

    from pyNastran.converters.dev.vrml.vrml import vrml_to_nastran, vrml_to_stl
    if fmt2 == 'nastran':
        vrml_to_nastran(vrml_filename, fname2, log=log)
    if fmt2 == 'stl':
        vrml_to_stl(vrml_filename, fname2, log=log)
    else:
        raise NotImplementedError('fmt2=%s is not supported by process_vrml' % fmt2)

if __name__ == '__main__':  # pragma: no cover
    cmd_line_format_converter()
