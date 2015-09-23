"""
Multi-input/output format converter
"""
from __future__ import print_function
from six import iteritems
import glob
from docopt import docopt

from numpy import zeros, array

import pyNastran
from pyNastran.bdf.bdf import BDF
from pyNastran.converters.tecplot.tecplot import Tecplot
from pyNastran.converters.tecplot.utils import merge_tecplot_files
from pyNastran.converters.stl.stl import STL
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
    elif fmt2 == 'tecplot':
        nastran_to_tecplot(model, fname2)
    # elif fmt2 == 'ugrid':
        # nastran_to_ugrid(model, fname2)
    else:
        raise NotImplementedError(fmt2)

def nastran_to_tecplot(model, tecplot_filename):
    """assumes sequential nodes"""
    tecplot = Tecplot()

    nnodes = model.nnodes
    xyz = zeros((nnodes, 3), dtype='float64')
    i = 0
    for nid, node in sorted(iteritems(model.nodes)):
        xyz[i, :] = node.get_position()
        i += 1
    tecplot.xyz = xyz

    nquads = model.card_count['CQUAD4'] if 'CQUAD4' in model.card_count else 0
    ntets = model.card_count['CTETRA'] if 'CTETRA' in model.card_count else 0
    #ntrias = model.card_count['CTRIA3'] if 'CTRIA3' in model.card_count else 0
    nhexas = model.card_count['CHEXA'] if 'CHEXA' in model.card_count else 0
    nelements = len(model.elements)
    tris = []
    quads = []
    tets = []
    hexas = []
    pentas = []
    #i = 0
    #pids = zeros(nelements, dtype='int32')
    #mids = zeros(nelements, dtype='int32')
    for eid, element in iteritems(model.elements):
        if element.type in ['CTRIA3']:
            tris.append(element.node_ids)
        elif element.type in ['CQUAD4']:
            quads.append(element.node_ids)
        elif element.type == 'CTETRA':
            tets.append(element.node_ids[:4])
        elif element.type == 'CPENTA':
            pentas.append(element.node_ids[:6])
        elif element.type == 'CHEXA':
            hexas.append(element.node_ids[:8])
        else:
            raise NotImplementedError(element.type)
        #pid = element.Pid()
        #mid = element.Mid()
        #pids[i] = pid
        #mids[i] = mid
        #i += 1

    # only supports nodal results
    #tecplot.results = vstack([pids, mids])#.T
    #print(tecplot.results.shape)
    #tecplot.result_names = ['PropertyID', 'MaterialID']

    ntris = len(tris)
    nquads = len(quads)
    nshells = ntris + nquads

    ntets = len(tets)
    npentas = len(pentas)
    nhexas = len(hexas)
    nsolids = ntets + npentas + nhexas
    nnot_tris = nquads
    nnot_quads = ntris
    nnot_tets = npentas + nhexas
    nnot_hexas = ntets + npentas
    if ntris and not nnot_tris and not nsolids:
        tecplot.tri_elements = array(tris, dtype='int32')
    elif nquads and not nnot_quads and not nsolids:
        tecplot.quad_elements = array(quads, dtype='int32')
    elif ntets and not nnot_tets and not nshells:
        tecplot.tet_elements = array(tets, dtype='int32')
    elif nhexas and not nnot_hexas and not nshells:
        tecplot.hexa_elements = array(hexas, dtype='int32')
    elif not nshells:
        elements = zeros((nelements, 8), dtype='int32')
        if ntets:
            tets = array(tets, dtype='int32')
            elements[:ntets, :4] = tets
            elements[:ntets, 4] = elements[:ntets, 3]
            elements[:ntets, 5] = elements[:ntets, 3]
            elements[:ntets, 6] = elements[:ntets, 3]
            elements[:ntets, 7] = elements[:ntets, 3]
        if npentas:
            # penta6
            pentas = array(pentas, dtype='int32')
            elements[ntets:ntets + npentas, :6] = pentas
            elements[ntets:ntets + npentas, 6] = elements[:ntets, 5]
            elements[ntets:ntets + npentas, 7] = elements[:ntets, 5]
        if nhexas:
            hexas = array(hexas, dtype='int32')
            elements[ntets + npentas:ntets + npentas + nhexas, :6] = pentas
            elements[ntets + npentas:ntets + npentas + nhexas, 6] = elements[:ntets, 5]
            elements[ntets + npentas:ntets + npentas + nhexas, 7] = elements[:ntets, 5]
        tecplot.hexa_elements = array(elements)
    elif not nsolids:
        elements = zeros((nelements, 4), dtype='int32')
        tris = array(tris, dtype='int32')
        elements[:ntris, :3] = tris
        elements[:ntris, 4] = elements[:ntets, 3]

        quads = array(quads, dtype='int32')
        elements[ntris:, :] = quads
    else:
        raise NotImplementedError()
    tecplot.write_tecplot(tecplot_filename, adjust_nids=False)


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
    """
    Converts Cart3d to Tecplot
    """
    model = Cart3D()
    model.read_cart3d(cart3d_filename)
    tecplot = Tecplot()
    tecplot.xyz = model.points
    tecplot.tri_elements = model.elements
    tecplot.write_tecplot(tecplot_filename, adjust_nids=False)

def process_stl(stl_filename, fmt2, fname2):
    """
    Converts STL to Nastran/Cart3d
    """
    model = STL()
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
    #model = Tecplot()
    if '*' in tecplot_filename:
        tecplot_filenames = glob.glob(tecplot_filename)
    else:
        tecplot_filenames = tecplot_filename
    tecplot = merge_tecplot_files(tecplot_filenames, tecplot_filename_out=None)
    # model.read_tecplot(tecplot_filename)
    # if fmt2 == 'nastran':
        # tecplot_to_nastran(model, fname2)
    #if fmt2 == 'cart3d':
        #tecplot_to_cart3d(model, fname2)
    #elif fmt2 == 'stl':
        #tecplot_to_stl(model, fname2)
    # elif fmt2 == 'ugrid':
        # tecplot_to_ugrid(model, fname2)
    if fmt2 == 'tecplot':
        # this is a good way to merge files
        tecplot.write_tecplot(fname2)
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
        model.write_bdf(bdf_filename, include_shells=include_shells, include_solids=include_solids)
    elif fmt2 == 'cart3d':
        include_shells = False
        include_solids = True
        bdf_filename = fname2 + '.bdf'
        model.write_bdf(bdf_filename, include_shells=include_shells, include_solids=include_solids)
        # ugrid_to_cart3d(model, fname2)
        process_nastran(bdf_filename, 'cart3d', fname2)
    elif fmt2 == 'stl':
        include_shells = False
        include_solids = True
        bdf_filename = fname2 + '.bdf'
        model.write_bdf(bdf_filename, include_shells=include_shells, include_solids=include_solids)
        process_nastran(bdf_filename, 'cart3d', fname2)
        # ugrid_to_stl(model, fname2)
    elif fmt2 == 'tecplot':
        # ugrid_to_tecplot(model, fname2)
        tecplot_filename = fname2
        model.write_tecplot(tecplot_filename)
    else:
        raise NotImplementedError(fmt2)

def run(fmt1, fname1, fmt2, fname2):
    """
    Runs the format converter
    """
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
    msg += "  format_convert nastran <INPUT> <format2> <OUTPUT> -o <OP2>\n"
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
    msg += "  OP2            path to results file (nastran-specific)\n"
    msg += "                 only used for Tecplot\n"
    msg += "  -h, --help     show this help message and exit\n"
    msg += "  -v, --version  show program's version number and exit\n\n"
    msg += 'Notes:\n'
    msg += "  Nastran->Tecplot assumes sequential nodes and consistent types (shell/solid)\n"
    msg += "  Tecplot supports globbing as the input filename\n"

    ver = str(pyNastran.__version__)
    data = docopt(msg, version=ver)
    #print(data)

    is_nastran = data['nastran']
    format1 = data['<format1>']
    if is_nastran:
        format1 = 'nastran'

    format2 = data['<format2>']
    input_filename = data['<INPUT>']
    output_filename = data['<OUTPUT>']
    op2_filename = data['<OP2>']
    run(format1, input_filename, format2, output_filename)

if __name__ == '__main__':
    main()
