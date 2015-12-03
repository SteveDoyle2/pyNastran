import os
import sys
from six import iteritems

from docopt import docopt
import pyNastran
#from gui.formats import format_string

def determine_format(input):
    """
    Tries to map the input filename to an extension.

    .. note :: this function will not support generic extensions (e.g. .inp, .dat)
    """
    format_to_extension = {
        'nastran' : ['.bdf', '.ecd', '.nas'],
        'stl' : ['.stl'],
        'cart3d' : ['.tri', '.triq'],
        'tecplot' : ['.plt'],
        'ugrid' : ['.ugrid'],
        'plot3d' : ['.p3d', '.p3da'],
    }
    ext = os.path.splitext(input)[1].lower()
    extension_to_format = {val : key for key, value in iteritems(format_to_extension)
                           for val in value}
    try:
        format = extension_to_format[ext]
    except:
        msg = 'format=%r was not found; Specify it' % ext
        raise TypeError(msg)
    return format

def run_arg_parse():
    msg  = "Usage:\n"
    msg += "  pyNastranGUI [-f FORMAT] [-i INPUT] [-o OUTPUT...]\n"
    msg += '                  [-s SHOT] [-m MAGNIFY]\n'  #  [-r XYZ]
    msg += '                  [-g GSCRIPT] [-p PSCRIPT]\n'
    msg += '                  [-u POINTS_FNAME...]\n'
    msg += '                  [-q]\n'
    msg += '  pyNastranGUI -h | --help\n'
    msg += '  pyNastranGUI -v | --version\n'
    msg += "\n"
    msg += "Options:\n"
    msg += "  -h, --help                  show this help message and exit\n"
    msg += "  -f FORMAT, --format FORMAT  format type (cart3d, lawgs, nastran, panair, \n"
    msg += "                                           plot3d, stl, tetgen, usm3d)\n"
    msg += "  -i INPUT, --input INPUT     path to input file\n"
    msg += "  -o OUTPUT, --output OUTPUT  path to output file\n"
    #msg += "  -r XYZ, --rotation XYZ      [x, y, z, -x, -y, -z] default is ???\n"

    #if mode != 'qt':
    msg += "  -g GSCRIPT, --geomscript GSCRIPT  path to geometry script file (runs before load geometry)\n"
    msg += "  -p PSCRIPT, --postscript PSCRIPT  path to post script file (runs after load geometry)\n"
    msg += "  -s SHOT, --shots SHOT       path to screenshot (only 1 for now)\n"
    msg += "  -m MAGNIFY, --magnify       how much should the resolution on a picture be magnified [default: 5]\n"
    msg += "  -u POINTS_FNAME, --user_points POINTS_FNAME               add user specified points to an alternate grid (repeatable)\n"

    msg += "  -q, --quiet                 prints debug messages (default=True)\n"
    #if mode != 'qt':
    msg += "  -v, --version               show program's version number and exit\n"


    ver = str(pyNastran.__version__)
    data = docopt(msg, version=ver)
    #print(data)

    format  = data['--format']
    input = data['--input']
    output = data['--output']
    debug = not(data['--quiet'])

    if input and not format:
        format = determine_format(input)


    shots = []
    if '--shots' in data:
        shots = data['--shots']

    geom_script = None
    if '--geomscript' in data:
        geom_script = data['--geomscript']

    post_script = None
    if '--postscript' in data:
        post_script = data['--postscript']

    magnify = 1
    if '--magnify' in data and data['--magnify'] is not None:
        magnify = int(data['--magnify'])

    rotation = None
    if '--rotation' in data:
        rotation = data['--rotation']

    user_points = data['--user_points']

    for key, value in sorted(iteritems(data)):
        print(key, value)
    #print("shots", shots)
    if shots:
        #shots = shots[1]
        #print("shots2 = %r" % shots, type(shots))
        shots = shots.split(';')[0]
    return (format, input, output, shots,
            magnify, rotation, geom_script, post_script, debug, user_points)


def get_inputs():
    is_edges = False
    format = None
    input = None
    output = None
    debug = True

    magnify = 1
    rotation = None
    shots = None
    geom_script = None
    post_script = None
    user_points = None

    if sys.version_info < (2, 6):
        print("requires Python 2.6+ to use command line arguments...")
    else:
        if len(sys.argv) > 1:
            (format, input, output, shots, magnify,
             rotation, geom_script, post_script, debug, user_points) = run_arg_parse()

    inputs = {
        'format' : format,
        'input' : input,
        'output' : output,
        'shots' : shots,
        'magnify' : magnify,
        'rotation' : rotation,
        'debug' : debug,
        'geomscript' : geom_script,
        'postscript' : post_script,
        'user_points' : user_points
    }
    return inputs
