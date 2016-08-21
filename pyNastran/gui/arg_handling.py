from __future__ import print_function
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
        'nastran' : ['.bdf', '.ecd', '.nas', '.op2'],
        'stl' : ['.stl'],
        'cart3d' : ['.tri', '.triq'],
        'tecplot' : ['.plt'],
        'ugrid' : ['.ugrid'],
        'plot3d' : ['.p3d', '.p3da'],
        'surf' : ['.surf'],
        'lawgs' : ['.wgs'],
        'shabp' : ['.mk5'],
        'panair' : ['.inp'],
    }
    ext = os.path.splitext(input)[1].lower()
    extension_to_format = {val : key for key, value in iteritems(format_to_extension)
                           for val in value}
    try:
        format = extension_to_format[ext]
    except:
        msg = 'format=%r was not found; Specify the format' % ext
        raise TypeError(msg)
    return format

def run_docopt():
    msg  = "Usage:\n"
    msg += "  pyNastranGUI [-f FORMAT] INPUT [-o OUTPUT]\n"
    msg += '               [-s SHOT] [-m MAGNIFY]\n'  #  [-r XYZ]
    msg += '               [-g GSCRIPT] [-p PSCRIPT]\n'
    msg += '               [-u POINTS_FNAME...] [--user_geom GEOM_FNAME...]\n'
    msg += '               [-q] [--groups]\n'
    msg += "  pyNastranGUI [-f FORMAT] INPUT OUTPUT [-o OUTPUT]\n"
    msg += '               [-s SHOT] [-m MAGNIFY]\n'  #  [-r XYZ]
    msg += '               [-g GSCRIPT] [-p PSCRIPT]\n'
    msg += '               [-u POINTS_FNAME...] [--user_geom GEOM_FNAME...]\n'
    msg += '               [-q] [--groups]\n'
    msg += "  pyNastranGUI [-f FORMAT] [-i INPUT] [-o OUTPUT...]\n"
    msg += '               [-s SHOT] [-m MAGNIFY]\n'  #  [-r XYZ]
    msg += '               [-g GSCRIPT] [-p PSCRIPT]\n'
    msg += '               [-u POINTS_FNAME...] [--user_geom GEOM_FNAME...]\n'
    msg += '               [-q] [--groups]\n'
    msg += '  pyNastranGUI -h | --help\n'
    msg += '  pyNastranGUI -v | --version\n'
    msg += "\n"
    msg += "Primary Options:\n"
    msg += "  -f FORMAT, --format FORMAT  format type (cart3d, lawgs, nastran, panair, stl, surf, ugrid)\n"
    msg += "  -i INPUT, --input INPUT     path to input file\n"
    msg += "  -o OUTPUT, --output OUTPUT  path to output file\n"
    #msg += "  -r XYZ, --rotation XYZ      [x, y, z, -x, -y, -z] default is ???\n"
    msg += '\n'

    msg += "Secondary Options:\n"
    msg += "  -g GSCRIPT, --geomscript        path to geometry script file (runs before load geometry)\n"
    msg += "  -p PSCRIPT, --postscript        path to post script file (runs after load geometry)\n"
    msg += "  -s SHOT, --shots SHOT           path to screenshot (only 1 for now)\n"
    msg += "  -m MAGNIFY, --magnify           how much should the resolution on a picture be magnified [default: 5]\n"
    msg += "  --groups                        enables groups\n"
    msg += "  --user_geom GEOM_FNAME          add user specified points to an alternate grid (repeatable)\n"
    msg += "  -u POINTS_FNAME, --user_points  add user specified points to an alternate grid (repeatable)\n"
    msg += '\n'

    msg += "Info:\n"
    msg += "  -q, --quiet    prints debug messages (default=True)\n"
    msg += "  -h, --help     show this help message and exit\n"
    msg += "  -v, --version  show program's version number and exit\n"


    ver = str(pyNastran.__version__)
    data = docopt(msg, version=ver)
    #print(data)

    input_format = data['--format']

    input_filenames = []
    if data['INPUT']:
        input_filenames += [data['INPUT']]
    if data['--input']:
        input_filenames += [data['--input']]

    output_filenames = []
    if data['OUTPUT']:
        output_filenames += [data['OUTPUT']]
    if data['--output']:
        output_filenames += data['--output']
    debug = not(data['--quiet'])

    if input_filenames and not input_format:
        input_format = determine_format(input_filenames[0])

    # None is for custom geometry
    allowed_formats = [
        'nastran', 'stl', 'cart3d', 'tecplot', 'ugrid', 'panair', 'plot3d',
        'surf', 'lawgs', 'degen_geom', 'shabp', 'avus', 'fast', None]
    assert input_format in allowed_formats, 'format=%r is not supported' % input_format

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
    user_geom = data['--user_geom']

    for key, value in sorted(iteritems(data)):
        print(key, value)
    #print("shots", shots)
    if shots:
        #shots = shots[1]
        #print("shots2 = %r" % shots, type(shots))
        shots = shots.split(';')[0]

    is_groups = data['--groups']
    #assert data['--console'] == False, data['--console']
    return (input_format, input_filenames, output_filenames, shots,
            magnify, rotation, geom_script, post_script, debug, user_points, user_geom, is_groups)


def get_inputs():
    input_format = None
    input_filename = None
    output_filename = None
    debug = True

    magnify = 5
    rotation = None
    shots = None
    geom_script = None
    post_script = None
    user_points = None
    user_geom = None
    is_groups = False

    if sys.version_info < (2, 6):
        print("requires Python 2.6+ to use command line arguments...")
    else:
        if len(sys.argv) > 1:
            (input_format, input_filename, output_filename, shots, magnify,
             rotation, geom_script, post_script, debug, user_points, user_geom, is_groups) = run_docopt()

    inputs = {
        'format' : input_format,
        'input' : input_filename,
        'output' : output_filename,
        'shots' : shots,
        'magnify' : magnify,
        'rotation' : rotation,
        'debug' : debug,
        'geomscript' : geom_script,
        'postscript' : post_script,
        'user_points' : user_points,
        'user_geom' : user_geom,
        'is_groups' : is_groups,
    }
    return inputs
