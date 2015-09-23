from docopt import docopt

import pyNastran
from pyNastran.bdf.bdf import BDF
from pyNastran.converters.tecplot.tecplot_reader import Tecplot
from pyNastran.converters.stl.stl_reader import STL
# from pyNastran.converters.cart3d.cart3d_reader import Cart3DReader
from pyNastran.converters.cart3d.cart3d import Cart3D

from pyNastran.converters.nastran.nastran_to_cart3d import nastran_to_cart3d
from pyNastran.converters.nastran.nastran_to_stl import nastran_to_stl
# from pyNastran.converters.nastran.nastran_to_tecplot import nastran_to_tecplot

from pyNastran.converters.stl.stl_to_nastran import stl_to_nastran_filename
from pyNastran.converters.cart3d.cart3d_to_nastran import cart3d_to_nastran_filename
from pyNastran.converters.cart3d.cart3d_to_stl import cart3d_to_stl_filename
from pyNastran.converters.ugrid.ugrid_reader import UGRID


def process_nastran(bdf_filename, fmt2, fname2):
    """
    Converts Nastran to STL/Cart3d
    """
    model = BDF()
    model.read_bdf(bdf_filename)
    if fmt2 == 'stl':
        nastran_to_stl(model, fname2)
    elif fmt2 == 'cart3d':
        nastran_to_cart3d(model, fname2)
    #elif fmt2 == 'tecplot':
        #nastran_to_tecplot(model, fname2)
    # elif fmt2 == 'ugrid':
        # nastran_to_ugrid(model, fname2)
    else:
        raise NotImplementedError(fmt2)

def process_cart3d(cart3d_filename, fmt2, fname2):
    """
    Converts Cart3d to STL/Nastran
    """
    # model = Cart3DReader()
    # model.read_cart3d(cart3d_filename, fname2)
    if fmt2 == 'stl':
        cart3d_to_stl_filename(cart3d_filename, fname2)
    elif fmt2 == 'nastran':
        cart3d_to_nastran_filename(cart3d_filename, fname2)
    elif fmt2 == 'tecplot':
        cart3d_to_tecplot(cart3d_filename, fname2)
    # elif fmt2 == 'ugrid':
        # cart3d_to_ugrid(model, fname2)
    else:
        raise NotImplementedError(fmt2)

def cart3d_to_tecplot(cart3d_filename, tecplot_filename):
    model = Cart3D()
    model.read_cart3d(cart3d_filename)
    tecplot = TecplotReader()
    tecplot.xyz = model.points
    tecplot.tri_elements = model.elements
    tecplot.write_tecplot(tecplot_filename, adjust_nids=False)

def process_stl(stl_filename, fmt2, fname2):
    """
    Converts STL to Nastran/Cart3d
    """
    model = STLReader()
    model.read_stl(stl_filename)
    if fmt2 == 'nastran':
        stl_to_nastran_filename(stl_filename, fname2)
        # stl_to_cart3d(model, fname2 + '.cart')
        # process_cart3d(fname2 + '.cart', 'nastran', fname2)
        # stl_to_nastran(model, fname2)
    elif fmt2 == 'cart3d':
        # stl_to_cart3d(model, fname2)
        stl_to_nastran_filename(stl_filename, fname2 + '.bdf')
        stl_to_nastran_filename(fname2 + '.bdf', fname2)
    # elif fmt2 == 'tecplot':
        # stl_to_tecplot(model, fname2)
    # elif fmt2 == 'ugrid':
        # stl_to_ugrid(model, fname2)
    else:
        raise NotImplementedError(fmt2)

def process_tecplot(tecplot_filename, fmt2, fname2):
    """
    Converts Tecplot to ...
    """
    # model = Tecplot()
    # model.read_tecplot(tecplot_filename)
    # if fmt2 == 'nastran':
        # tecplot_to_nastran(model, fname2)
    #if fmt2 == 'cart3d':
        #tecplot_to_cart3d(model, fname2)
    #elif fmt2 == 'stl':
        #tecplot_to_stl(model, fname2)
    # elif fmt2 == 'ugrid':
        # tecplot_to_ugrid(model, fname2)
    raise NotImplementedError(fmt2)

def process_ugrid(ugrid_filename, fmt2, fname2):
    """
    Converts UGRID to Nastran/Cart3d/STL/Tecplot
    """
    model = UGRID()
    model.read_ugrid(ugrid_filename)
    if fmt2 == 'nastran':
        # ugrid_to_nastran(model, fname2
        include_shells = True
        include_solids = True
        bdf_filename = fname2
        model.export_bdf(bdf_filename, include_shells=include_shells, include_solids=include_solids)
    elif fmt2 == 'cart3d':
        include_shells = False
        include_solids = True
        bdf_filename = fname2 + '.bdf'
        model.export_bdf(bdf_filename, include_shells=include_shells, include_solids=include_solids)
        # ugrid_to_cart3d(model, fname2)
        process_nastran(bdf_filename, 'cart3d', fname2)
    elif fmt2 == 'stl':
        include_shells = False
        include_solids = True
        bdf_filename = fname2 + '.bdf'
        model.export_bdf(bdf_filename, include_shells=include_shells, include_solids=include_solids)
        process_nastran(bdf_filename, 'cart3d', fname2)
        # ugrid_to_stl(model, fname2)
    elif fmt2 == 'tecplot':
        # ugrid_to_tecplot(model, fname2)
        tecplot_filename = fname2
        model.write_tecplot(tecplot_filename)
    else:
        raise NotImplementedError(fmt2)

def run(fmt1, fname1, fmt2, fname2):
    if fmt1 == 'nastran':
        process_nastran(fname1, fmt2, fname2)
    elif fmt1 == 'cart3d':
        process_cart3d(fname1, fmt2, fname2)
    elif fmt1 == 'stl':
        process_stl(fname1, fmt2, fname2)
    elif fmt1 == 'tecplot':
        process_tecplot(fname1, fmt2, fname2)
    elif fmt1 == 'ugrid':
        process_ugrid(fname1, fmt2, fname2)
    else:
        raise NotImplementedError(fmt1)


def main():
    msg = "Usage:\n"
    msg += "  format_convert <format1> <INPUT> <format2> <OUTPUT>\n"
    #msg += "  format_convert nastran  <INPUT> <format2> <OUTPUT>\n"
    #msg += "  format_convert cart3d   <INPUT> <format2> <OUTPUT>\n"
    msg += '  format_convert -h | --help\n'
    msg += '  format_convert -v | --version\n'
    msg += "\n"
    msg += "Options:\n"
    msg += "  format1        format type (nastran, cart3d, stl, ugrid, tecplot)\n"
    msg += "  format2        format type (nastran, cart3d, stl, ugrid, tecplot)\n"
    msg += "  INPUT          path to input file\n"
    msg += "  OUTPUT         path to output file\n"
    msg += "  -h, --help     show this help message and exit\n"
    msg += "  -v, --version  show program's version number and exit\n"


    ver = str(pyNastran.__version__)
    data = docopt(msg, version=ver)
    #print(data)

    format1 = data['<format1>']
    format2 = data['<format2>']
    input_filename = data['<INPUT>']
    output_filename = data['<OUTPUT>']
    run(format1, input_filename, format2, output_filename)

if __name__ == '__main__':
    main()
