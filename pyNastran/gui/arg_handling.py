import os
import sys

import pyNastran
from pyNastran.utils import check_path
from pyNastran.utils.arg_handling import argparse_to_dict, swap_key, update_message

#from gui.formats import format_string

if sys.version_info < (3, 7):  # pragma: no cover
    sys.exit("requires Python 3.7+...")

SUPPORT_MULTIMODEL = False

# True: --groups flag
# False: --nogroups flag
GROUPS_DEFAULT = True
FORMAT_TO_EXTENSION = {
    # an extension should not be added to this list if it is shared with another type
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
    'vrml' : ['.wrl'],
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
            'avl', 'vrml',
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


def get_inputs(print_inputs=False, argv=None):
    """Gets the inputs for pyNastranGUI using docopt."""
    if argv is None:
        argv = sys.argv[1:]  # same as argparse
        #print('get_inputs; argv was None -> %s' % argv)
    else:
        # drop the pyNastranGUI; same as argparse
        argv = argv[1:]

    input_format = None
    input_filename = None
    output_filename = None
    debug = True

    geom_script = None
    post_script = None
    user_points = None
    user_geom = None
    is_groups = not GROUPS_DEFAULT
    log = None
    test = False

    if len(argv) >= 1:
        argdict = run_argparse(argv)
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
        'log' : log,
        'test' : test,
    }
    return inputs

def run_argparse(argv):
    """Gets the inputs for pyNastranGUI using argparse."""
    import argparse
    #msg = "Usage:\n"

    # INPUT format may be explicitly or implicitly defined with or
    # without an output file
    #test = ' [--test]'
    #qt = ' [--qt QT] [--plugin]'

    #msg += "  pyNastranGUI INPUT [-f FORMAT] [-o OUTPUT]\n"
    #msg += '               [-q] [--groups] [--noupdate] [--log LOG]%s%s\n' % (test, qt)

    # You don't need to throw a -o flag
    #msg += "  pyNastranGUI INPUT OUTPUT [-f FORMAT] [-o OUTPUT]\n"
    #msg += '               [-q] [--groups] [--noupdate] [--log LOG]%s%s\n' % (test, qt)

    dev = ''
    dev_list = []
    if not pyNastran.is_pynastrangui_exe:
        #dev = ' [--noupdate] [--test] [--qt Qt] [--plugin]'
        dev_list = ['--noupdate', '--test', '--qt', '--plugin']
        dev = ''.join([' [%s]' % devi for devi in dev_list])

    # no input/output files
    # can you ever have an OUTPUT, but no INPUT?
    usage = 'Usage:\n'
    usage += '  pyNastranGUI INPUT [-f FORMAT] [-o OUTPUT] [options]\n'

    # You don't need to throw a -o flag
    usage += '  pyNastranGUI INPUT OUTPUT [-f FORMAT] [-o OUTPUT] [options]\n'

    # no input/output files
    # can you ever have an OUTPUT, but no INPUT?
    usage += '  pyNastranGUI [-f FORMAT] [-i INPUT] [-o OUTPUT...] [options]\n'
    #usage += '  pyNastranGUI -h | --help\n'
    usage += '  pyNastranGUI -v | --version\n'

    usage += (
        '  [options] = [-g GSCRIPT] [-p PSCRIPT]\n'
        '              [-u POINTS_FNAME...] [--user_geom GEOM_FNAME...]\n'
    )
    if GROUPS_DEFAULT:
        usage += '              [-q] [--groups] [--log LOG]%s\n' % (dev)
    else:
        usage += '              [-q] [--nogroups] [--log LOG]%s\n' % (dev)

    arg_msg = ''
    arg_msg += '\n'
    arg_msg += 'Primary Options:\n'
     # plot3d,
    arg_msg += '  -f FORMAT, --format FORMAT  format type (avus, bedge, cart3d, lawgs, nastran,\n'
    arg_msg += '                                  openfoam_hex, openfoam_shell, openfoam_faces,\n'
    arg_msg += '                                  panair, stl, surf, tetgen, usm3d, ugrid, ugrid3d)\n'
    arg_msg += '  -i INPUT, --input INPUT     path to input file\n'
    arg_msg += '  -o OUTPUT, --output OUTPUT  path to output file\n'
    #arg_msg += "  -r XYZ, --rotation XYZ      [x, y, z, -x, -y, -z] default is ???\n"
    arg_msg += '\n'

    arg_msg += 'Secondary Options:\n'
    if GROUPS_DEFAULT:
        arg_msg += '  --groups                        enables groups\n'
    else:
        arg_msg += '  --nogroups                      disables groups\n'
    arg_msg += '  -g GSCRIPT, --geomscript        path to geometry script file (runs before load geometry)\n'
    arg_msg += '  -p PSCRIPT, --postscript        path to post script file (runs after load geometry)\n'
    arg_msg += '  --user_geom GEOM_FNAME          add user specified geometry (repeatable)\n'
    arg_msg += '  -u POINTS_FNAME, --user_points  add user specified points (repeatable)\n'
    arg_msg += '\n'

    arg_msg += "Debug:\n"
    if not pyNastran.is_pynastrangui_exe:
        arg_msg += '  --noupdate     disables the update check\n'
        arg_msg += '  --test         temporary dev mode (default=False)\n'
        arg_msg += '  --qt QT        sets the qt version (pyqt5, pyside2; default=QT_API)\n'
        arg_msg += '  --plugin       disables the format check\n'
    arg_msg += '  --log LOG      disables HTML logging; prints to the screen\n'
    arg_msg += '\n'

    arg_msg += 'Info:\n'
    arg_msg += '  -q, --quiet    prints debug messages (default=True)\n'
    arg_msg += '  -h, --help     show this help message and exits\n'
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
    #import textwrap
    parent_parser = argparse.ArgumentParser(
        #prog = 'pyNastranGUI',
        #usage = usage,
        #description='A foo that bars',
        #epilog="And that's how you'd foo a bar",
        #formatter_class=argparse.RawDescriptionHelpFormatter,
        #description=textwrap.dedent(text),
        #version=pyNastran.__version__,
        #add_help=False,
    )
    # pyNastranGUI INPUT OUTPUT [-f FORMAT] [-o OUTPUT]
    # positional arguments
    parent_parser.add_argument('INPUT', nargs='?', help='path to input file', type=str)
    parent_parser.add_argument('OUTPUT', nargs='?', help='path to output file', type=str)

    #nargs : str/int
    #   * : 0 or more
    #   + : one or more
    #   ? : optional
    #   int : int values
    #SUPPORT_MULTIMODEL = False
    #append_nargs = 1 if SUPPORT_MULTIMODEL else 1
    append_action = 'append' if SUPPORT_MULTIMODEL else None
    parent_parser.add_argument('-i', '--input', help='path to input file',
                               nargs=1, action=append_action)
    parent_parser.add_argument('-o', '--output', help='path to output file',
                               nargs=1, action=append_action)
    #parent_parser.add_argument('--user_geom', type=str, help='log msg')

    parent_parser.add_argument('-f', '--format', type=str, nargs=1, action=append_action,
                               help='format type (avus, bedge, cart3d, lawgs, nastran, '
                               'openfoam_hex, openfoam_shell, openfoam_faces, panair, '
                               'stl, surf, tetgen, usm3d, ugrid, ugrid3d, #plot3d)')

    # double args
    #parent_parser.add_argument('-f', '--format', type=str,
                               #help='format type (avus, bedge, cart3d, lawgs, nastran, '
                               #'openfoam_hex, openfoam_shell, openfoam_faces, panair, '
                               #'stl, surf, tetgen, usm3d, ugrid, ugrid3d, #plot3d)',
                               #action='append')
    parent_parser.add_argument('-g', '--geomscript', type=str,
                               help='path to geometry script file (runs before load geometry)')
    parent_parser.add_argument('-p', '--postscript', type=str,
                               help='path to post script file (runs after load geometry)')
    parent_parser.add_argument(
        '-u', '--points_fname', type=str, action='append',
        help='an (nrows, 3) comma/tab/space separated list of points (repeatable)')
    parent_parser.add_argument('--user_geom', type=str, action='append',
                               help='add user specified geometry (repeatable)')
    parent_parser.add_argument('--log', type=str, help='{debug, info, warning, error} msg')

    # no arguments
    if dev:
        parent_parser.add_argument('--qt', type=str, help='{pyqt5, pyside2} msg')
        parent_parser.add_argument('--test', help='test msg', action='store_true')
        parent_parser.add_argument('--noupdate', help='noupdate msg', action='store_true')
        parent_parser.add_argument('--plugin', help='disables the format check',
                                   action='store_true')

    if GROUPS_DEFAULT:
        parent_parser.add_argument('--groups', help='enables groups', action='store_true')
    else:
        parent_parser.add_argument('--nogroups', help='disables groups', action='store_false')

    parent_parser.add_argument('-q', '--quiet',
                               help='prints debug messages (default=True)', action='store_true')
    #parent_parser.add_argument('-h', '--help', help='show this help message and exits',
                               #action='store_true')
    parent_parser.add_argument('-v', '--version', action='version',
                               version=pyNastran.__version__)

    #foo_parser = argparse.ArgumentParser(parents=[parent_parser])
    #foo_parser.parse_args(['INPUT', '--format', '--output',
                           #'--geomscript', '--postscript', '--points_fname', '--user_geom',
                           #'--quiet', '--groups', '--log' '--help'] + dev_list)

    #msg += "  pyNastranGUI INPUT [-f FORMAT] [-o OUTPUT]\n"
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

    update_message(parent_parser, usage, arg_msg, examples)
    args = parent_parser.parse_args(args=argv)
    #args.plugin = True
    argdict = argparse_to_dict(args)
    _update_argparse_argdict(argdict)
    return argdict

def _add_inputs_outputs(positional_inputs, optional_inputs, word='input'):
    input_filenames = []
    if isinstance(optional_inputs, str):
        input_filenames.append(optional_inputs)
    elif isinstance(optional_inputs, list):
        for input_filenamesi in optional_inputs:
            if isinstance(input_filenamesi, str):
                input_filenames.append(input_filenamesi)
            elif isinstance(input_filenamesi, list):
                input_filenames.extend(input_filenamesi)
            else:
                raise TypeError('%s_filenamesi=%s type=%s' % (
                    word, input_filenames, type(input_filenamesi)))
    #print('input_filenames =', input_filenames)

    if isinstance(positional_inputs, str):
        input_filenames += [positional_inputs]
    return input_filenames

def _update_argparse_argdict(argdict):
    """converts to the pyNastranGUI argument format"""
    argdict['debug'] = not argdict['quiet']
    del argdict['quiet']

    if not GROUPS_DEFAULT:
        swap_key(argdict, 'nogroups', 'is_groups')
    else:
        argdict['is_groups'] = argdict['groups']

    swap_key(argdict, 'points_fname', 'user_points')

    input_filenames = _add_inputs_outputs(argdict['INPUT'], argdict['input'], word='input')
    del argdict['INPUT']
    argdict['input'] = input_filenames

    output_filenames = _add_inputs_outputs(argdict['OUTPUT'], argdict['output'], word='output')
    del argdict['OUTPUT']
    argdict['output'] = output_filenames

    for output_filename in output_filenames:
        check_path(output_filename, name='output_filename')
    for input_filename in input_filenames:
        check_path(input_filename, name='input file')

    plugin = False
    if 'plugin' in argdict:
        plugin = True

    formats = argdict['format']
    if input_filenames and formats is None:
        input_formats = []
        for input_filenamei in input_filenames:
            if isinstance(input_filenamei, str):
                formati = determine_format(input_filenamei)
            else:  # pragma: no cover
                raise TypeError('input_filenamei=%s type=%s' % (
                    input_filenamei, type(input_filenamei)))
            input_formats.append(formati)
        #input_formats = [determine_format(input_filenamei) for input_filenamei in input_filenames]
        argdict['format'] = input_formats
        del input_formats
    elif formats:
        input_formats = []
        for formati in formats:
            if isinstance(formati, str):
                input_formats.append(formati)
            else:
                input_formats.extend(formati)
        argdict['format'] = input_formats
        del input_formats
    del formats

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
        for input_format in input_formats:
            assert input_format in allowed_formats, 'format=%r is not supported' % input_format

    if argdict['geomscript']:
        geom_script = argdict['geomscript']
        check_path(geom_script, name='geomscript')
    if argdict['postscript']:
        post_script = argdict['postscript']
        check_path(post_script, name='postscript')

    if argdict['qt'] is not None:
        qt = argdict['qt'].lower()
        assert qt in ['pyqt5', 'pyside2'], 'qt=%r' % qt
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
        #'log' : log,
        #'test' : test,
    #}

    formats = argdict['format']
    ninput_files = len(input_filenames)
    if formats:
        if isinstance(formats, str):
            formats = [formats]
        nformats = len(formats)
        if nformats == 1 and ninput_files > 1:
            formats = formats * ninput_files
        argdict['format'] = formats
        if nformats != ninput_files:
            msg = (
                'nformats=%s formats=%s\n'
                'ninput_files=%s input_filenames=%s' % (
                    nformats, formats,
                    ninput_files, input_filenames))
            raise RuntimeError(msg)
    return argdict
