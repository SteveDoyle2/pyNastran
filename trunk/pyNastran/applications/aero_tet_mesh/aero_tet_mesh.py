from pyNastran.utils.log import get_logger

log = None
debug = True
log = get_logger(log, 'debug' if debug else 'info')

def run_arg_parse():
    msg  = 'This program creates a CFD tet mesh from a nastran surface..\n'
    msg += 'Usage:\n'
    msg += '  stl_reader INPUT [-o OUTPUT]\n'
    msg += '             [-n] [-q]\n'
    msg += '  stl_reader -h | --help\n'
    msg += '  stl_reader -v | --version\n'
    msg += "  INPUT      path to input file\n"
    msg += "\n"
    msg += "Options:\n"
    msg += "  -h, --help                  show this help message and exit\n"
    msg += "  -f FORMAT, --format FORMAT  format type (usm3d, fun3d)\n"
    msg += "  -o OUTPUT, --output OUTPUT  path to output file\n"
    msg += "  -n, --normal                flip the element normals\n"
    #msg += "  -r XYZ, --rotation XYZ      [x, y, z, -x, -y, -z] default is ???\n"

    msg += "  -q, --quiet                 prints debug messages (default=True)\n"
    msg += "  -v, --version               show program's version number and exit\n"

    ver = str(pyNastran.__version__)
    data = docopt(msg, version=ver)
    #print data

    #format  = data['--format']
    input = data['INPUT']
    print "input =", input
    output = data['--output']

    reverse_normals = data['--normal']
    quiet = data['--quiet']
    return (input, output, reverse_normals, quiet)

from numpy import array, zeros, transpose
from pyNastran.bdf.bdf import BDF
#from pyNastran.converters.cart3d.cart3d_reader import Cart3DReader

from pyNastran.converters.nastran.nastran_to_cart3d import nastran_to_cart3d_filename
from pyNastran.converters.cart3d.cart3d_to_stl import cart3d_to_stl_filename
#from pyNastran.converters.cart3d.cart3d_to_nastran import cart3d_to_nastran_filename

from pyNastran.converters.nastran.nastran_to_stl import nastran_to_stl_filename
from pyNastran.converters.stl.stl_to_nastran import stl_to_nastran_filename
from pyNastran.converters.stl.stl_reader import STLReader



def main():
    bdf_filename = 'bay.bdf'
    bc_file = 'bc'

    bdf_filename2 = 'bay2.bdf'
    cart3d_filename = 'bay.i.tri'
    stl_filename = 'bay.stl'
    #in_format = 'nastran'
    flip_normals = True

    volume_bdfname = 'bay.vol.bdf'

    module = __import__(bc_file)
    pid_map = module.__dict__['pid_map']
    boundary_conditions = module.__dict__['boundary_conditions']

    #print dir(module)
    #e = execfile(bc_file, globals(), locals())
    #print globals()
    #print dir(e)
    
    pid_map
    boundary_conditions
    assert isinstance(pid_map, dict)
    assert isinstance(boundary_conditions, dict)
    
    get_bcs(bdf_filename, pid_map, boundary_conditions, log=None)
    #if flip_normals:
    #    bdf = BDF()
    #    bdf.read_bdf(bdf_filename, xref=False)
    #    bdf.flip_normals()
    #    bdf.write_bdf(bdf_filename2)
    #    del bdf
    #else:
    #bdf_filename2 = bdf_filename

    if 0:
        nastran_to_cart3d_filename(bdf_filename, cart3d_filename, log=log)
        cart3d_to_stl_filename(cart3d_filename, stl_filename, log=log)
    else:
        nastran_to_stl_filename(bdf_filename, stl_filename, log=log)
        stl_to_nastran_filename(stl_filename, bdf_filename2, log=log)
    print "----------"
    #nastran_to_cart3d(bdf_filename2, cart3d_filename)
    #cart3d = Cart3dReader()
    #cart3d.read_cart3d(cart3d_filename)

    stl = STLReader()
    stl.read_stl(stl_filename)

    if flip_normals:
        stl.flip_normals()
    stl.project_boundary_layer(stl.nodes, stl.elements, volume_bdfname)


def get_bcs(bdf_filename, pid_map, boundary_conditions, log=None):
    bdf = BDF(log=log)
    bdf.read_bdf(bdf_filename)
    
    ntris = len(bdf.elements)
    boundary_condition_ids = zeros(ntris, 'int32')
    i = 0
    
    for eid, element in sorted(bdf.elements.iteritems()):
        assert element.type == 'CTRIA3', element.type
        boundary_condition_ids[i] = pid_map[element.Pid()]
        i += 1
    bc_set = set(boundary_condition_ids)
    assert len(bc_set) == 1, 'len(bc_set) = %s' % len(bc_set)

if __name__ == '__main__':
    main()