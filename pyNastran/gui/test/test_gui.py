# coding: latin1
"""defines the command line tool test_gui"""
import os
import sys
import time
import itertools
import traceback

from docopt import docopt

import pyNastran
from cpylog import get_logger

from pyNastran.gui.errors import NoGeometry
from pyNastran.bdf.errors import (CrossReferenceError, CardParseSyntaxError,
                                  DuplicateIDsError, MissingDeckSections)

from pyNastran.gui.testing_methods import FakeGUIMethods
from pyNastran.gui.formats import CLASS_MAP
from pyNastran.converters.nastran.gui.nastran_io import NastranIO

from pyNastran.gui.arg_handling import determine_format
from pyNastran.utils import print_bad_path, check_path
from pyNastran.utils.dev import get_files_of_type
from pyNastran.op2.test.op2_test import get_failed_files

FORMAT_TO_EXTENSION = {
    #'abaqus' : ['.inp'],
    'nastran' : ['.bdf', '.ecd', '.nas', '.op2', '.pch', '.dat'],
    'stl' : ['.stl'],
    'cart3d' : ['.tri', '.triq'],
    'tecplot' : ['.plt', '.dat'],
    'ugrid' : ['.ugrid'],
    #'plot3d' : ['.p3d', '.p3da'],
    'surf' : ['.surf'],
    'lawgs' : ['.wgs'],
    'shabp' : ['.mk5'],
    'fast' : ['.cogsg'],
    'usm3d' : ['.cogsg', '.front'],
    'bedge' : ['.bedge'],
    'su2' : ['.su2'],
    'tetgen' : ['.smesh', '.ele'],
    'avus' : ['.grd'],
    'avl' : ['.avl'],
    'vrml' : ['.wrl'],

    # no duplicates are allowed
    #'panair' : ['.inp'],
    #'abaqus' : ['.inp'],
}

EXTENSION_TO_OUPUT_FORMATS = {
    'nastran' : ['.op2'],
    'shabp' : ['.out'],
    #'usm3d' : ['.out'],
}


#pkg_path = pyNastran.__path__[0]
#model_path = os.path.join(pkg_path, '..', 'models')

class FakeGUI(FakeGUIMethods, NastranIO):
    """spoofs the gui for testing"""

    def __init__(self, formati, inputs=None):
        """
        Parameters
        ----------
        formati : str
            the file format (e.g., 'nastran', 'cart3d')
        inputs : None
            pass for the fake gui
        """
        self._formati = formati
        FakeGUIMethods.__init__(self, inputs=inputs)
        NastranIO.__init__(self)
        self.format_class_map = CLASS_MAP

    def load_geometry(self, input_filename):
        """loads a model"""
        load_geometry_name = 'load_%s_geometry' % self._formati
        if self._formati in self.format_class_map:
            cls = self.format_class_map[self._formati](self)
            getattr(cls, load_geometry_name)(input_filename)
        elif hasattr(self, load_geometry_name):
            # self.load_nastran_geometry(bdf_filename, None)
            getattr(self, load_geometry_name)(input_filename)
        else:
            msg = "load_geometry_name=%s doesn't exist" % load_geometry_name
            raise NotImplementedError(msg)

    def load_results(self, output_filename):
        """loads a model"""
        load_results_name = 'load_%s_results' % self._formati
        if self._formati in self.format_class_map:
            cls = self.format_class_map[self._formati](self)
            getattr(cls, load_results_name)(output_filename)
        elif hasattr(self, load_results_name):
            # self.load_nastran_ressults(op2_filename, None)
            getattr(self, load_results_name)(output_filename)
        else:
            msg = "load_results_name=%s doesn't exist" % load_results_name
            raise NotImplementedError(msg)


def run_docopt(argv=None):
    """
    The main function for the command line ``test_pynastran_gui`` script.
    """
    msg = "Usage:\n"
    # INPUT format may be explicitly or implicitly defined with or
    # without an output file
    msg += "  test_pynastrangui [-f FORMAT]           INPUT_FILENAME  OUTPUT_FILENAME [--log LOG] [--test]\n"
    msg += "  test_pynastrangui [-f FORMAT]           INPUT_FILENAME  [--log LOG] [--test]\n"
    msg += "  test_pynastrangui  -f FORMAT  [-r] [-d] INPUT_DIRECTORY [--log LOG] [--test]\n"
    msg += "  test_pynastrangui  -f FORMAT  [-r] [-d]                 [--log LOG] [--test]\n"

    msg += '  test_pynastrangui -h | --help\n'
    msg += '  test_pynastrangui -v | --version\n'
    msg += '\n'

    msg += 'Positional Arguments:\n'
    msg += '  INPUT_FILENAME   path to input file\n'
    msg += '  OUTPUT_FILENAME  path to output file\n'
    msg += '  INPUT_DIRECTORY  path to input directory\n'
    msg += '\n'

    msg += "Options:\n"
    msg += '  -f FORMAT, --format  format type (avus, cart3d, lawgs, nastran, panair,\n'
    msg += '                                    su2, stl, surf, tetgen, usm3d, ugrid)\n'
    msg += '  -d, --dir            directory to run tests on\n'
    msg += "  -r, --regenerate     Resets the tests\n"
    msg += '  --log LOG            debug, info, warning, error; default=debug\n'
    msg += '\n'

    msg += "Debug:\n"
    msg += "  --test    temporary dev mode (default=False)\n"

    msg += 'Info:\n'
    #msg += "  -q, --quiet    prints debug messages (default=True)\n"
    msg += '  -h, --help     show this help message and exit\n'
    msg += "  -v, --version  show program's version number and exit\n"

    if len(sys.argv) == 1:
        sys.exit(msg)
    ver = str(pyNastran.__version__)
    data = docopt(msg, argv=argv, version=ver, options_first=False) # help=True,

    isdir = data['INPUT_DIRECTORY'] or not data['INPUT_FILENAME'] or data['--dir']
    if isdir or data['--dir']:
        formati = data['--format'].lower()
        #print("sys.version_info[:2]) =", sys.version_info[:2])
        failed_cases_filename = 'failed_cases_%s_%s%s.in' % (
            formati, sys.version_info[0], sys.version_info[1])
        print(failed_cases_filename)

        if data['--regenerate'] or not os.path.exists(failed_cases_filename):
            dirname = data['INPUT_DIRECTORY']
            if not os.path.exists(dirname):
                msg = 'dirname=%r does not exist\n%s' % (
                    dirname, print_bad_path(dirname))
                raise RuntimeError(msg)
            if not os.path.isdir(dirname):
                msg = 'dirname=%r is not a directory' % dirname
                raise RuntimeError(msg)
            extensions = FORMAT_TO_EXTENSION[formati]
            input_filenames = [
                get_files_of_type(
                    dirname, extension=extension,
                    max_size=100., limit_file='no_dig.txt')
                for extension in extensions
            ]
            input_filenames = list(itertools.chain.from_iterable(input_filenames))
        else:
            input_filenames = get_failed_files(failed_cases_filename)
        output_filenames = [None] * len(input_filenames)
    else:
        failed_cases_filename = None
        input_filename = data['INPUT_FILENAME']
        output_filename = data['OUTPUT_FILENAME']
        check_path(input_filename, 'input_filename')
        if not os.path.isfile(input_filename):
            msg = 'input_filename=%r is not a file' % input_filename
            raise RuntimeError(msg)
        input_filenames = [input_filename]
        output_filenames = [output_filename]

        if data['--format']:
            formati = data['--format'].lower()
        else:
            formati = determine_format(input_filename)
            print('formati', formati)
    #assert formati == 'nastran', 'format=%r' % formati

    if data['--log']:
        log_method = data['--log'].lower()
        assert log_method in ['debug', 'info', 'warning', 'error'], 'log_method=%r' % log_method
    else:
        log_method = 'debug'
    return formati, input_filenames, output_filenames, failed_cases_filename, log_method, data['--test']


def main():
    """runs the gui"""
    (formati, input_filenames, output_filenames,
     failed_cases_filename, log_method, test) = run_docopt()
    log = get_logger(log=None, level=log_method, encoding='utf-8')
    npass = 0
    nfailed = 0
    failed_files = []
    ntotal = len(input_filenames)

    print('test =', test)
    test_gui = FakeGUI(formati)
    test_gui.log = log
    stop_on_failure = ntotal == 1
    time0 = time.time()
    for input_filename, output_filename in zip(input_filenames, output_filenames):
        input_filename = os.path.abspath(input_filename)
        #output_filename =
        print("filename = %s" % input_filename)
        is_passed = True
        try:
            test_gui.load_geometry(input_filename)
        except KeyboardInterrupt:
            sys.exit('KeyboardInterrupt...sys.exit()')
        except SystemExit:
            sys.exit('sys.exit...')

        except NoGeometry:
            traceback.print_exc(file=sys.stdout)
            print('failed test because NoGeometry...ignoring')
        except CrossReferenceError:
            traceback.print_exc(file=sys.stdout)
            print('failed test because CrossReferenceError...ignoring')
        except CardParseSyntaxError:
            traceback.print_exc(file=sys.stdout)
            print('failed test because CardParseSyntaxError...ignoring')
        except DuplicateIDsError:
            traceback.print_exc(file=sys.stdout)
            print('failed test because DuplicateIDsError...ignoring')
        except MissingDeckSections:
            traceback.print_exc(file=sys.stdout)
            print('failed test because MissingDeckSections...ignoring')
        except OverflowError:
            traceback.print_exc(file=sys.stdout)
            print('failed test because OverflowError...ignoring')

        except:
            is_passed = False
            traceback.print_exc(file=sys.stdout)
            if stop_on_failure:
                print('input_filenames =', input_filenames)
                raise

        if output_filename:
            test_gui.load_results(output_filename)

        if is_passed:
            sys.stderr.write('%i  %s' % (npass, input_filename))
            npass += 1
        else:
            sys.stderr.write('*%s ' % nfailed + input_filename)
            nfailed += 1
            failed_files.append(input_filename)
        sys.stderr.write('\n')

        if not stop_on_failure:
            print('-' * 80)
    dt = time.time() - time0
    nfailed = len(failed_files)
    npassed = ntotal - nfailed
    time_msg = 'dt = %.0f sec = %.1f min' % (dt, dt/60.)
    sys.stderr.write('%i/%i passed\n' % (npassed, ntotal))
    sys.stderr.write(time_msg)

    if ntotal > 1:
        with open(failed_cases_filename, 'w') as failed_cases_file:
            for fname in failed_files:
                failed_cases_file.write('%s\n' % fname)
        print(time_msg)
        sys.exit('finished...')

    #test_gui.load_nastran_geometry(bdf_filename, None)
    #test_gui.load_nastran_results(op2_filename, None)


if __name__ == '__main__':   # pragma: no cover
    main()
