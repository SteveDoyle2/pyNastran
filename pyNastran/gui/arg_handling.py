from __future__ import print_function
import os
import sys
from six import iteritems

from docopt import docopt
import pyNastran
from pyNastran.utils import print_bad_path
#from gui.formats import format_string

FORMAT_TO_EXTENSION = {
    'nastran' : ['.bdf', '.ecd', '.nas', '.op2', '.pch'],
    'stl' : ['.stl'],
    'cart3d' : ['.tri', '.triq'],
    'tecplot' : ['.plt'],
    'ugrid' : ['.ugrid'],
    'plot3d' : ['.p3d', '.p3da'],
    'surf' : ['.surf'],
    'lawgs' : ['.wgs'],
    'shabp' : ['.mk5'],
    'usm3d' : ['.front'],  # .cogsg
    'bedge' : ['.bedge'],
    'su2' : ['.su2'],
    'tetgen' : ['.smesh', '*.ele'],  # TODO: why does this have a *?

    # no duplicates are allowed
    #'panair' : ['.inp'],
    #'abaqus' : ['.inp'],
}

def determine_format(input_filename, allowed_formats=None):
    """
    Tries to map the input filename to an extension.

    .. note :: this function will not support generic extensions (e.g. .inp, .dat)
    """
    if allowed_formats is None:
        # used to include None...
        allowed_formats = [
            'nastran', 'stl', 'cart3d', 'tecplot', 'ugrid', 'panair',
            #'plot3d',
            'surf', 'lawgs', 'degen_geom', 'shabp', 'avus', 'fast', 'abaqus',
            'usm3d', 'bedge', 'su2', 'tetgen',
        ]

    ext = os.path.splitext(input_filename)[1].lower()
    extension_to_format = {val : key for key, value in iteritems(FORMAT_TO_EXTENSION)
                           for val in value}
    try:
        formati = extension_to_format[ext]
    except:
        print('allowed_formats =', allowed_formats)
        msg = 'format=%r was not found\nSpecify the format as [%s]' % (
            ext, ', '.join(allowed_formats))
        raise TypeError(msg)
    return formati

def run_docopt():
    msg  = "Usage:\n"

    # INPUT format may be explicitly or implicitly defined with or
    # without an output file
    msg += "  pyNastranGUI [-f FORMAT] INPUT [-o OUTPUT]\n"
    msg += '               [-s SHOT] [-m MAGNIFY]\n'  #  [-r XYZ]
    msg += '               [-g GSCRIPT] [-p PSCRIPT]\n'
    msg += '               [-u POINTS_FNAME...] [--user_geom GEOM_FNAME...]\n'
    msg += '               [-q] [--groups] [--noupdate] [--log LOG]\n'

    # You don't need to throw a -o flag
    msg += "  pyNastranGUI [-f FORMAT] INPUT OUTPUT [-o OUTPUT]\n"
    msg += '               [-s SHOT] [-m MAGNIFY]\n'  #  [-r XYZ]
    msg += '               [-g GSCRIPT] [-p PSCRIPT]\n'
    msg += '               [-u POINTS_FNAME...] [--user_geom GEOM_FNAME...]\n'
    msg += '               [-q] [--groups] [--noupdate] [--log LOG]\n'

    # no input/output files
    # can you ever have an OUTPUT, but no INPUT?
    msg += "  pyNastranGUI [-f FORMAT] [-i INPUT] [-o OUTPUT...]\n"
    msg += '               [-s SHOT] [-m MAGNIFY]\n'  #  [-r XYZ]
    msg += '               [-g GSCRIPT] [-p PSCRIPT]\n'
    msg += '               [-u POINTS_FNAME...] [--user_geom GEOM_FNAME...]\n'
    msg += '               [-q] [--groups] [--noupdate] [--log LOG]\n'
    msg += '  pyNastranGUI -h | --help\n'
    msg += '  pyNastranGUI -v | --version\n'
    msg += "\n"
    msg += "Primary Options:\n"
    msg += "  -f FORMAT, --format FORMAT  format type (avus, cart3d, lawgs, nastran, panair, \n"
    msg += "                                           stl, surf, tetgen, usm3d, ugrid)\n"
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
    msg += "  --noupdate                      disables the update check\n"
    msg += "  --user_geom GEOM_FNAME          add user specified points to an alternate grid (repeatable)\n"
    msg += "  -u POINTS_FNAME, --user_points  add user specified points to an alternate grid (repeatable)\n"
    msg += '\n'

    msg += "Info:\n"
    msg += "  -q, --quiet    prints debug messages (default=True)\n"
    msg += "  --log LOG      disables HTML logging; prints to the screen\n"
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
    for input_filename in input_filenames:
        assert os.path.exists(input_filename), print_bad_path(input_filename)

    output_filenames = []
    if data['OUTPUT']:
        output_filenames += [data['OUTPUT']]
    if data['--output']:
        output_filenames += data['--output']
    for output_filename in output_filenames:
        assert os.path.exists(output_filename), print_bad_path(output_filename)
    debug = not(data['--quiet'])

    if input_filenames and not input_format:
        input_format = determine_format(input_filenames[0])

    # None is for custom geometry
    allowed_formats = [
        'nastran', 'stl', 'cart3d', 'tecplot', 'ugrid', 'panair',
        'surf', 'lawgs', 'shabp', 'avus', 'fast',
        'usm3d', 'bedge', 'su2', 'tetgen',
        None,
    ]
    assert input_format in allowed_formats, 'format=%r is not supported' % input_format

    shots = []
    if '--shots' in data:
        shots = data['--shots']

    geom_script = data['--geomscript']
    if geom_script:
        assert os.path.exists(geom_script), print_bad_path(geom_script)

    post_script = data['--postscript']
    if post_script:
        assert os.path.exists(post_script), print_bad_path(post_script)

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
    no_update = data['--noupdate']
    #assert data['--console'] == False, data['--console']
    return (input_format, input_filenames, output_filenames, shots,
            magnify, rotation, geom_script, post_script, debug, user_points,
            user_geom, is_groups, no_update, data['--log'])


def get_inputs(argv=None):
    if argv is None:
        argv = sys.argv
        #print('get_inputs; argv was None -> %s' % argv)
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
    no_update = True
    log = None

    if sys.version_info < (2, 6):
        print("requires Python 2.6+ to use command line arguments...")
    else:
        if len(argv) > 1:
            (input_format, input_filename, output_filename, shots, magnify,
             rotation, geom_script, post_script, debug, user_points, user_geom,
             is_groups, no_update, log) = run_docopt()

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
        'no_update' : no_update,
        'log' : log,
    }
    return inputs
