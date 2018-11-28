from __future__ import print_function
import os
import sys

from docopt import docopt
import pyNastran
from pyNastran.utils import check_path
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
    'tetgen' : ['.smesh', '.ele'],
    'obj' : ['.obj'],
    'fast' : ['.fgrid'],
    'avl' : ['.avl'],
    #'abaqus' : []

    # no duplicate extensions are allowed; use the explicit --format option
    #'ugrid3d' : ['.ugrid'],
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
            'nastran', 'stl', 'cart3d', 'tecplot', 'ugrid', 'ugrid3d', 'panair',
            #'plot3d',
            'surf', 'lawgs', 'degen_geom', 'shabp', 'avus', 'fast', 'abaqus',
            'usm3d', 'bedge', 'su2', 'tetgen', 'obj',
            'openfoam_hex', 'openfoam_shell', 'openfoam_faces',
            'avl',
        ]

    ext = os.path.splitext(input_filename)[1].lower()
    extension_to_format = {val : key for key, value in FORMAT_TO_EXTENSION.items()
                           for val in value}
    try:
        formati = extension_to_format[ext]
    except KeyError:
        print('allowed_formats =', allowed_formats)
        msg = 'format=%r was not found\nSpecify the format as [%s]' % (
            ext, ', '.join(allowed_formats))
        raise TypeError(msg)
    return formati

def run_docopt():
    """helper for ``get_inputs_docopt``"""
    msg = "Usage:\n"

    # INPUT format may be explicitly or implicitly defined with or
    # without an output file
    test = ''
    qt = ''
    if not pyNastran.is_pynastrangui_exe:
        test = ' [--test]'
        qt = ' [--qt QT] [--plugin]'

    msg += "  pyNastranGUI INPUT [-f FORMAT] [-o OUTPUT]\n"
    msg += '               [-g GSCRIPT] [-p PSCRIPT]\n'
    msg += '               [-u POINTS_FNAME...] [--user_geom GEOM_FNAME...]\n'
    msg += '               [-q] [--groups] [--noupdate] [--log LOG]%s%s\n' % (test, qt)

    # You don't need to throw a -o flag
    msg += "  pyNastranGUI INPUT OUTPUT [-f FORMAT] [-o OUTPUT]\n"
    msg += '               [-g GSCRIPT] [-p PSCRIPT]\n'
    msg += '               [-u POINTS_FNAME...] [--user_geom GEOM_FNAME...]\n'
    msg += '               [-q] [--groups] [--noupdate] [--log LOG]%s%s\n' % (test, qt)

    # no input/output files
    # can you ever have an OUTPUT, but no INPUT?
    msg += "  pyNastranGUI [-f FORMAT] [-i INPUT] [-o OUTPUT...]\n"
    msg += '               [-g GSCRIPT] [-p PSCRIPT]\n'
    msg += '               [-u POINTS_FNAME...] [--user_geom GEOM_FNAME...]\n'
    msg += '               [-q] [--groups] [--noupdate] [--log LOG]%s%s\n' % (test, qt)
    msg += '  pyNastranGUI -h | --help\n'
    msg += '  pyNastranGUI -v | --version\n'
    msg += "\n"
    msg += "Primary Options:\n"
    msg += "  -f FORMAT, --format FORMAT  format type (avus, bedge, cart3d, lawgs, nastran,\n" # plot3d,
    msg += '                                  openfoam_hex, openfoam_shell, openfoam_faces,\n'
    msg += "                                  panair, stl, surf, tetgen, usm3d, ugrid, ugrid3d)\n"
    msg += "  -i INPUT, --input INPUT     path to input file\n"
    msg += "  -o OUTPUT, --output OUTPUT  path to output file\n"
    #msg += "  -r XYZ, --rotation XYZ      [x, y, z, -x, -y, -z] default is ???\n"
    msg += '\n'

    msg += "Secondary Options:\n"
    msg += "  --groups                        enables groups\n"
    msg += "  -g GSCRIPT, --geomscript        path to geometry script file (runs before load geometry)\n"
    msg += "  -p PSCRIPT, --postscript        path to post script file (runs after load geometry)\n"
    msg += "  --user_geom GEOM_FNAME          add user specified geometry (repeatable)\n"
    msg += "  -u POINTS_FNAME, --user_points  add user specified points (repeatable)\n"
    msg += '\n'

    msg += "Debug:\n"
    if not pyNastran.is_pynastrangui_exe:
        msg += "  --test         temporary dev mode (default=False)\n"
        msg += "  --qt QT        sets the qt version (default=QT_API)\n"
        msg += "  --plugin       disables the format check\n"
    msg += "  --noupdate     disables the update check\n"
    msg += "  --log LOG      disables HTML logging; prints to the screen\n"
    msg += '\n'

    msg += "Info:\n"
    msg += "  -q, --quiet    prints debug messages (default=True)\n"
    msg += "  -h, --help     show this help message and exits\n"
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
        check_path(input_filename, 'input file')

    output_filenames = []
    if data['OUTPUT']:
        output_filenames += [data['OUTPUT']]
    if data['--output']:
        output_filenames += data['--output']
    for output_filename in output_filenames:
        check_path(output_filename, 'output_filename')
    debug = not(data['--quiet'])

    if input_filenames and input_format is None:
        input_format = determine_format(input_filenames[0])

    plugin = False
    if '--plugin' in data:
        plugin = True

    if not plugin:
        # None is for custom geometry
        allowed_formats = [
            'nastran', 'stl', 'cart3d', 'tecplot', 'ugrid', 'ugrid3d', 'panair',
            #'plot3d',
            'surf', 'lawgs', 'degen_geom', 'shabp', 'avus', 'fast', 'abaqus',
            'usm3d', 'bedge', 'su2', 'tetgen',
            'openfoam_hex', 'openfoam_shell', 'openfoam_faces', 'obj', 'avl',
            None,
        ]
        assert input_format in allowed_formats, 'format=%r is not supported' % input_format

    geom_script = data['--geomscript']
    if geom_script:
        check_path(geom_script, name='geom_script')

    post_script = data['--postscript']
    if post_script:
        check_path(post_script, 'post_script')

    user_points = data['--user_points']
    user_geom = data['--user_geom']

    if data['--qt'] is not None:
        qt = data['--qt'].lower()
        assert qt in ['pyside', 'pyqt4', 'pyqt5', 'pyside2'], 'qt=%r' % qt
        os.environ.setdefault('QT_API', qt)

    for key, value in sorted(data.items()):
        print(key, value)
    #print("shots", shots)

    is_groups = data['--groups']
    no_update = data['--noupdate']
    test = ''
    if '--test' in data:
        test = data['--test']
    #assert data['--console'] == False, data['--console']
    return (input_format, input_filenames, output_filenames,
            geom_script, post_script, debug, user_points,
            user_geom, is_groups, no_update, data['--log'], test)


def get_inputs(print_inputs=False, argv=None):
    """Gets the inputs for pyNastranGUI using docopt."""
    if argv is None:
        argv = sys.argv
        #print('get_inputs; argv was None -> %s' % argv)
    input_format = None
    input_filename = None
    output_filename = None
    debug = True

    geom_script = None
    post_script = None
    user_points = None
    user_geom = None
    is_groups = False
    no_update = True
    log = None
    test = False

    if sys.version_info < (2, 7, 7):
        print("requires Python 2.7.7+...")
    else:
        if len(argv) > 1:
            if 0:
                (input_format, input_filename, output_filename,
                 geom_script, post_script, debug, user_points, user_geom,
                 is_groups, no_update, log, test) = run_docopt()
            else:
                argdict = run_argparse()
                if print_inputs:
                    for key, value in sorted(argdict.items()):
                        print(key, value)
                return argdict

    inputs = {
        'format' : input_format,
        'input' : input_filename,
        'output' : output_filename,
        'debug' : debug,
        'geomscript' : geom_script,
        'postscript' : post_script,
        'user_points' : user_points,
        'user_geom' : user_geom,
        'is_groups' : is_groups,
        'no_update' : no_update,
        'log' : log,
        'test' : test,
    }
    return inputs

def run_argparse():
    """Gets the inputs for pyNastranGUI using argparse."""
    import argparse
    #msg = "Usage:\n"

    # INPUT format may be explicitly or implicitly defined with or
    # without an output file
    #test = ' [--test]'
    #qt = ' [--qt QT] [--plugin]'

    #msg += "  pyNastranGUI INPUT [-f FORMAT] [-o OUTPUT]\n"
    #msg += '               [-g GSCRIPT] [-p PSCRIPT]\n'
    #msg += '               [-u POINTS_FNAME...] [--user_geom GEOM_FNAME...]\n'
    #msg += '               [-q] [--groups] [--noupdate] [--log LOG]%s%s\n' % (test, qt)

    # You don't need to throw a -o flag
    #msg += "  pyNastranGUI INPUT OUTPUT [-f FORMAT] [-o OUTPUT]\n"
    #msg += '               [-g GSCRIPT] [-p PSCRIPT]\n'
    #msg += '               [-u POINTS_FNAME...] [--user_geom GEOM_FNAME...]\n'
    #msg += '               [-q] [--groups] [--noupdate] [--log LOG]%s%s\n' % (test, qt)

    dev = ''
    dev_list = []
    if not pyNastran.is_pynastrangui_exe:
        #dev = ' [--noupdate] [--test] [--qt Qt] [--plugin]'
        dev_list = ['--noupdate', '--test', '--qt', '--plugin']
        dev = ''.join([' [%s]' % devi for devi in dev_list])

    # no input/output files
    # can you ever have an OUTPUT, but no INPUT?
    usage = "Usage:\n"
    usage += "  pyNastranGUI INPUT [-f FORMAT] [-o OUTPUT]\n"
    usage += '               [-g GSCRIPT] [-p PSCRIPT]\n'
    usage += '               [-u POINTS_FNAME...] [--user_geom GEOM_FNAME...]\n'
    usage += '               [-q] [--groups] [--noupdate] [--log LOG]%s\n' % (dev)

    #parent_parser.add_argument('-g', '--geomscript', type=str, help='path to geometry script file (runs before load geometry)', action='append')
    #parent_parser.add_argument('-p', '--postscript', type=str, help='path to post script file (runs after load geometry)', action='append')
    #parent_parser.add_argument('-u', '--points_fname', type=str, help='an (nrows, 3) comma/tab/space separated list of points')
    #parent_parser.add_argument('--user_geom', type=str, help='add user specified geometry (repeatable)')

    # You don't need to throw a -o flag
    usage += "  pyNastranGUI INPUT OUTPUT [-f FORMAT] [-o OUTPUT]\n"
    usage += '               [-g GSCRIPT] [-p PSCRIPT]\n'
    usage += '               [-u POINTS_FNAME...] [--user_geom GEOM_FNAME...]\n'
    usage += '               [-q] [--groups] [--log LOG]%s\n' % (dev)

    # no input/output files
    # can you ever have an OUTPUT, but no INPUT?
    usage += "  pyNastranGUI [-f FORMAT] [-i INPUT] [-o OUTPUT...]\n"
    usage += '               [-g GSCRIPT] [-p PSCRIPT]\n'
    usage += '               [-u POINTS_FNAME...] [--user_geom GEOM_FNAME...]\n'
    usage += '               [-q] [--groups] [--log LOG]%s\n' % (dev)
    #usage += '  pyNastranGUI -h | --help\n'
    usage += '  pyNastranGUI -v | --version\n'

    arg_msg = ''
    arg_msg += "\n"
    arg_msg += "Primary Options:\n"
    arg_msg += "  -f FORMAT, --format FORMAT  format type (avus, bedge, cart3d, lawgs, nastran,\n" # plot3d,
    arg_msg += '                                  openfoam_hex, openfoam_shell, openfoam_faces,\n'
    arg_msg += "                                  panair, stl, surf, tetgen, usm3d, ugrid, ugrid3d)\n"
    arg_msg += "  -i INPUT, --input INPUT     path to input file\n"
    arg_msg += "  -o OUTPUT, --output OUTPUT  path to output file\n"
    #arg_msg += "  -r XYZ, --rotation XYZ      [x, y, z, -x, -y, -z] default is ???\n"
    arg_msg += '\n'

    arg_msg += "Secondary Options:\n"
    arg_msg += "  --groups                        enables groups\n"
    arg_msg += "  -g GSCRIPT, --geomscript        path to geometry script file (runs before load geometry)\n"
    arg_msg += "  -p PSCRIPT, --postscript        path to post script file (runs after load geometry)\n"
    arg_msg += "  --user_geom GEOM_FNAME          add user specified geometry (repeatable)\n"
    arg_msg += "  -u POINTS_FNAME, --user_points  add user specified points (repeatable)\n"
    arg_msg += '\n'

    arg_msg += "Debug:\n"
    if not pyNastran.is_pynastrangui_exe:
        arg_msg += "  --test         temporary dev mode (default=False)\n"
        arg_msg += "  --qt QT        sets the qt version (default=QT_API)\n"
        arg_msg += "  --plugin       disables the format check\n"
    arg_msg += "  --noupdate     disables the update check\n"
    arg_msg += "  --log LOG      disables HTML logging; prints to the screen\n"
    arg_msg += '\n'

    arg_msg += "Info:\n"
    arg_msg += "  -q, --quiet    prints debug messages (default=True)\n"
    arg_msg += "  -h, --help     show this help message and exits\n"
    arg_msg += "  -v, --version  show program's version number and exit\n"
    arg_msg += '\n'

    #msg += "\n"
    #parser = argparse.ArgumentParser(
        #prog=None, usage=None, description=None, epilog=None,
        #version=None, parents=[], formatter_class=HelpFormatter,
        #prefix_chars='-', fromfile_prefix_chars=None, argument_default=None,
        #conflict_handler='error', add_help=True)

    #usage = '[options]'
    examples = (
        'Examples\n'
        '--------\n'
        '  pyNastranGUI\n'
        '  pyNastranGUI fem.bdf\n'
        '  pyNastranGUI fem.bdf fem.op2\n'
        '  pyNastranGUI fem.dat fem.op2 --format nastran -o fem2.op2\n\n'
    )
    import textwrap
    parent_parser = argparse.ArgumentParser(
        #prog = 'pyNastranGUI',
        #usage = usage,
        #description='A foo that bars',
        epilog="And that's how you'd foo a bar",
        #formatter_class=argparse.RawDescriptionHelpFormatter,
        #description=textwrap.dedent(text),
        #version=pyNastran.__version__,
        #add_help=False,
    )
    # positional arguments
    parent_parser.add_argument('INPUT', nargs='?', help='path to input file', type=str)
    parent_parser.add_argument('OUTPUT', nargs='?', help='path to output file', type=str)

    parent_parser.add_argument('-i', '--input', help='path to input file')
    parent_parser.add_argument('-o', '--output', help='path to output file')
    #parent_parser.add_argument('--user_geom', type=str, help='log msg')

    parent_parser.add_argument('-f', '--format', type=str,
                               help='format type (avus, bedge, cart3d, lawgs, nastran, '
                               'openfoam_hex, openfoam_shell, openfoam_faces, panair, '
                               'stl, surf, tetgen, usm3d, ugrid, ugrid3d, #plot3d)')

    # double args
    #parent_parser.add_argument('-f', '--format', type=str,
                               #help='format type (avus, bedge, cart3d, lawgs, nastran, '
                               #'openfoam_hex, openfoam_shell, openfoam_faces, panair, '
                               #'stl, surf, tetgen, usm3d, ugrid, ugrid3d, #plot3d)', action='append')
    parent_parser.add_argument('-g', '--geomscript', type=str, help='path to geometry script file (runs before load geometry)', action='append')
    parent_parser.add_argument('-p', '--postscript', type=str, help='path to post script file (runs after load geometry)', action='append')
    parent_parser.add_argument('-u', '--points_fname', type=str, help='an (nrows, 3) comma/tab/space separated list of points')
    parent_parser.add_argument('--user_geom', type=str, help='add user specified geometry (repeatable)')
    parent_parser.add_argument('--log', type=str, help='{debug, info, warning, error} msg')

    # no arguments
    if dev:
        parent_parser.add_argument('--qt', type=str, help='{pyqt4, pyqt5, pyside, pyside2} msg')
        parent_parser.add_argument('--test', help='test msg', action='store_true')
        parent_parser.add_argument('--noupdate', help='noupdate msg', action='store_true')
    parent_parser.add_argument('--groups', help='enables groups', action='store_true')
    parent_parser.add_argument('--plugin', help='disables the format check', action='store_true')

    parent_parser.add_argument('-q', '--quiet', help='prints debug messages (default=True)', action='store_true')
    #parent_parser.add_argument('-h', '--help', help='show this help message and exits', action='store_true')
    parent_parser.add_argument('-v', '--version', action='version',
                               version=pyNastran.__version__)

    #foo_parser = argparse.ArgumentParser(parents=[parent_parser])
    #foo_parser.parse_args(['INPUT', '--format', '--output',
                           #'--geomscript', '--postscript', '--points_fname', '--user_geom',
                           #'--quiet', '--groups', '--no_update', '--log' '--help'] + dev_list)

    #msg += "  pyNastranGUI INPUT [-f FORMAT] [-o OUTPUT]\n"
    #msg += '               [-g GSCRIPT] [-p PSCRIPT]\n'
    #msg += '               [-u POINTS_FNAME...] [--user_geom GEOM_FNAME...]\n'
    #msg += '               [-q] [--groups] [--noupdate] [--log LOG]%s%s\n' % (test, qt)
    #parser_no_output = p

    #parser = argparse.ArgumentParser(
        #description='A foo that bars',
        #epilog="And that's how you'd foo a bar",
        #version=pyNastran.__version__,
    #)
    #parser.add_argument("square", help="display a square of a given number",
                        #type=int)
    #parser.add_argument('-v', '--version', action='version',
                        #version='%%(prog)s %s' % pyNastran.__version__)
    #parser.add_argument("-w", "--verbosity", type=int, choices=[0, 1, 2],
                        #help="increase output verbosity")

    #mymsg = 'replacing argparse message'
    #def print_help(self, file=None):
        #if file is None:
            #file = _sys.stdout
        #self._print_message(self.format_help(), file)

    mymsg = usage + arg_msg + examples
    def _print_message(message, file=None):
        """overwrites the argparse print to get a better help message"""
        if message:
            if file is None:
                file = _sys.stderr
            file.write(mymsg)
    parent_parser._print_message = _print_message
    args = parent_parser.parse_args()
    #args.plugin = True
    argdict = argparse_to_dict(args)
    _update_argparse_argdict(argdict)
    return argdict

def _update_argparse_argdict(argdict):
    """converts to the pyNastranGUI argument format"""
    argdict['debug'] = not(argdict['quiet'])
    del argdict['quiet']

    swap_key(argdict, 'groups', 'is_groups')
    swap_key(argdict, 'noupdate', 'no_update')
    swap_key(argdict, 'points_fname', 'user_points')

    input_filenames = []
    if isinstance(argdict['input'], str):
        input_filenames += [argdict['input']]
    if isinstance(argdict['INPUT'], str):
        input_filenames += [argdict['INPUT']]
    del argdict['INPUT']
    argdict['input'] = input_filenames

    output_filenames = []
    if isinstance(argdict['output'], str):
        output_filenames += [argdict['output']]
    if isinstance(argdict['OUTPUT'], str):
        output_filenames += [argdict['OUTPUT']]
    del argdict['OUTPUT']
    argdict['output'] = output_filenames

    for output_filename in output_filenames:
        check_path(output_filename, name='output_filename')
    for input_filename in input_filenames:
        check_path(input_filename, name='input file')

    plugin = False
    if 'plugin' in argdict:
        plugin = True

    if not plugin:
        # None is for custom geometry
        allowed_formats = [
            'nastran', 'stl', 'cart3d', 'tecplot', 'ugrid', 'ugrid3d', 'panair',
            #'plot3d',
            'surf', 'lawgs', 'degen_geom', 'shabp', 'avus', 'fast', 'abaqus',
            'usm3d', 'bedge', 'su2', 'tetgen',
            'openfoam_hex', 'openfoam_shell', 'openfoam_faces', 'obj', 'avl',
            None,
        ]
        assert input_format in allowed_formats, 'format=%r is not supported' % input_format

    if input_filenames and argdict['format'] is None:
        input_format = determine_format(input_filenames[0])
        argdict['format'] = input_format

    if argdict['geomscript']:
        check_path(geom_script, name='geomscript')
    if argdict['postscript']:
        check_path(post_script, name='postscript')

    if argdict['qt'] is not None:
        qt = argdict['qt'].lower()
        assert qt in ['pyside', 'pyqt4', 'pyqt5', 'pyside2'], 'qt=%r' % qt
        os.environ.setdefault('QT_API', qt)

    #if argdict['input'] is None:
        #argdict['input'] = []

    #inputs = {
        #'format' : input_format,
        #'input' : input_filename,
        #'output' : output_filename,
        #'debug' : debug,
        #'geomscript' : geom_script,
        #'postscript' : post_script,
        #'user_points' : user_points,
        #'user_geom' : user_geom,
        #'is_groups' : is_groups,
        #'no_update' : no_update,
        #'log' : log,
        #'test' : test,
    #}
    #print(argdict)
    return argdict

def swap_key(mydict, key_orig, key_new):
    """replaces a key in a dictionary"""
    mydict[key_new] = mydict[key_orig]
    del mydict[key_orig]

def argparse_to_dict(args):
    """converts the argparse output into a dictionary"""
    argdict = {}
    for name, value in args._get_args():
        argdict[name] = value
    for name, value in args._get_kwargs():
        argdict[name] = value
    return argdict
