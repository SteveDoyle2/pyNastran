import sys

import pyNastran
from docopt import docopt

def run_arg_parse(mode):
    msg  = "Usage:\n"
    msg += "  pyNastranGUI_%s [-f FORMAT] [-i INPUT] [-o OUTPUT]\n" % mode
    msg += '                  [-s SHOT] [-m MAGNIFY] [-p SCRIPT]\n'  #  [-r XYZ]
    msg += '                  [-q] [-e] [-n | -c]\n'
    msg += '  pyNastranGUI_%s -h | --help\n' % mode
    msg += '  pyNastranGUI_%s -v | --version\n' % mode
    msg += "\n"
    msg += "Options:\n"
    msg += "  -h, --help                  show this help message and exit\n"
    msg += "  -f FORMAT, --format FORMAT  format type (cart3d, lawgs, nastran, panair, \n"
    msg += "                                           plot3d, stl, tetgen, usm3d)\n"
    msg += "  -i INPUT, --input INPUT     path to input file\n"
    msg += "  -o OUTPUT, --output OUTPUT  path to output file\n"
    #msg += "  -r XYZ, --rotation XYZ      [x, y, z, -x, -y, -z] default is ???\n"

    #if mode != 'qt':
    msg += "  -p SCRIPT, --pyscript SCIPRT  path to script file\n"
    msg += "  -s SHOT, --shots SHOT       path to screenshot (only 1 for now)\n"
    msg += "  -m MAGNIFY, --magnify       MAGNIFY how much should the resolution on a picture be magnified (default=1)\n"

    msg += "  -q, --quiet                 prints debug messages (default=True)\n"
    #if mode != 'qt':
    msg += "  -e, --edges                 shows element edges as black lines (default=False)\n"
    msg += "  -n, --nodalResults          plots nodal results (default)\n"
    msg += "  -c, --centroidalResults     plots centroidal results\n"
    msg += "  -v, --version               show program's version number and exit\n"

    ver = str(pyNastran.__version__)
    data = docopt(msg, version=ver)
    #print data

    format  = data['--format']
    input   = data['--input']
    output  = data['--output']
    debug   = not(data['--quiet'])
    is_nodal = data['--nodalResults']
    is_centroidal = data['--centroidalResults']

    edges = False
    if '--edges' in data:
        edges = data['--edges']

    shots = []
    if '--shots' in data:
        shots = data['--shots']

    script = None
    if '--script' in data:
        script = data['--script']

    magnify = 1
    if '--magnify' in data and data['--magnify'] is not None:
        magnify = int(data['--magnify'])

    rotation = None
    if '--rotation' in data:
        rotation = data['--rotation']

    #print("is_nodal=%s is_centroidal=%s" % (is_nodal, is_centroidal))
    #print("shots", shots)
    #writeBDF = args.writeBDF
    if shots:
        #shots = shots[1]
        #print "shots2 = %r" % shots, type(shots)
        shots = shots.split(';')[0]
    return (edges, is_nodal, is_centroidal, format, input, output, shots,
            magnify, rotation, script, debug)


def get_inputs(mode='wx'):
    assert mode in ['wx', 'qt']
    is_edges = False
    format = None
    input = None
    output = None
    debug = True

    is_centroidal = True
    is_nodal = not(is_centroidal)
    magnify = 1
    rotation = None
    shots = None

    if sys.version_info < (2, 6):
        print("requires Python 2.6+ to use command line arguments...")
    else:
        if len(sys.argv) > 1:
            (is_edges, is_nodal, is_centroidal, format, input, output, shots, magnify,
             rotation, script, debug) = run_arg_parse(mode)

    if format is not None and is_centroidal == is_nodal:
        if format.lower() == 'nastran':
            is_nodal = False
            is_centroidal = True
        elif format.lower() == 'cart3d':
            is_nodal = False
            is_centroidal = True
        elif format.lower() == 'panair':
            is_nodal = False
            is_centroidal = True
    if is_centroidal == is_nodal:
        is_nodal = not(is_centroidal)

    print "is_centroidal =", is_centroidal
    inputs = {
        'is_edges' : is_edges,
        'is_nodal' : is_nodal,
        'is_centroidal' : is_centroidal,
        'format' : format,
        'input' : input,
        'output' : output,
        'shots' : shots,
        'magnify' : magnify,
        'rotation' : rotation,
        'debug' : debug,
        'script' : script,
        }
    return inputs