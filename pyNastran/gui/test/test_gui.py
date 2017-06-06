"""
defines the command line tool test_gui
"""
from __future__ import print_function
import os
import sys
import time
import itertools
import traceback

from six import PY2
from docopt import docopt

import pyNastran
from pyNastran.utils.log import get_logger

from pyNastran.gui.errors import NoGeometry
from pyNastran.bdf.errors import (CrossReferenceError, CardParseSyntaxError,
                                  DuplicateIDsError, MissingDeckSections)

from pyNastran.gui.testing_methods import FakeGUIMethods
from pyNastran.converters.nastran.nastranIOv import NastranIO
from pyNastran.converters.cart3d.cart3d_io import Cart3dIO
from pyNastran.converters.panair.panair_io import PanairIO
from pyNastran.converters.fast.fast_io import FastIO
from pyNastran.converters.LaWGS.wgs_io import LaWGS_IO
from pyNastran.converters.shabp.shabp_io import ShabpIO
from pyNastran.converters.stl.stl_io import STL_IO
from pyNastran.converters.su2.su2_io import SU2_IO
from pyNastran.converters.tecplot.tecplot_io import TecplotIO
from pyNastran.converters.tetgen.tetgen_io import TetgenIO
from pyNastran.converters.usm3d.usm3d_io import Usm3dIO
from pyNastran.converters.abaqus.abaqus_io import AbaqusIO

from pyNastran.converters.aflr.aflr2.bedge_io import BEdge_IO
from pyNastran.converters.aflr.ugrid.surf_io import SurfIO
from pyNastran.converters.aflr.ugrid.ugrid_io import UGRID_IO

from pyNastran.gui.arg_handling import determine_format
from pyNastran.utils import print_bad_path
from pyNastran.utils.dev import get_files_of_type
from pyNastran.op2.test.test_op2 import get_failed_files

FORMAT_TO_EXTENSION = {
    #'abaqus' : ['.inp'],
    'nastran' : ['.bdf', '.ecd', '.nas', '.op2', '.pch', '.dat'],
    'stl' : ['.stl'],
    'cart3d' : ['.tri', '.triq'],
    'tecplot' : ['.plt', '.dat'],
    'ugrid' : ['.ugrid'],
    'plot3d' : ['.p3d', '.p3da'],
    'surf' : ['.surf'],
    'lawgs' : ['.wgs'],
    'shabp' : ['.mk5'],
    'fast' : ['.cogsg'],
    'usm3d' : ['.cogsg', '.front'],
    'bedge' : ['.bedge'],
    'su2' : ['.su2'],
    'tetgen' : ['.smesh', '.ele'],

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

class FakeGUI(FakeGUIMethods, NastranIO, AbaqusIO, Cart3dIO, ShabpIO,
              PanairIO, LaWGS_IO, STL_IO, TetgenIO, Usm3dIO,
              #Plot3d_io, ADB_IO, DegenGeomIO,
              # AbaqusIO, AvusIO,
              TecplotIO, FastIO, SurfIO, UGRID_IO, BEdge_IO, SU2_IO):
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
        #ADB_IO.__init__(self)
        #AvusIO.__init__(self)
        AbaqusIO.__init__(self)
        BEdge_IO.__init__(self)
        NastranIO.__init__(self)
        Cart3dIO.__init__(self)
        #DegenGeomIO.__init__(self)
        FastIO.__init__(self)
        LaWGS_IO.__init__(self)
        PanairIO.__init__(self)
        #Plot3d_io.__init__(self)
        STL_IO.__init__(self)
        ShabpIO.__init__(self)
        SurfIO.__init__(self)
        TetgenIO.__init__(self)
        TecplotIO.__init__(self)
        Usm3dIO.__init__(self)
        UGRID_IO.__init__(self)
        #AbaqusIO.__init__(self)
        SU2_IO.__init__(self)

    def load_geometry(self, input_filename):
        """loads a model"""
        load_geometry_name = 'load_%s_geometry' % self._formati
        if hasattr(self, load_geometry_name):
            # self.load_nastran_geometry(bdf_filename, None)
            dirname = None
            getattr(self, load_geometry_name)(input_filename, dirname)
        else:
            msg = "load_geometry_name=%s doesn't exist" % load_geometry_name
            raise NotImplementedError(msg)

    def load_results(self, output_filename):
        """loads a model"""
        load_results_name = 'load_%s_results' % self._formati
        if hasattr(self, load_results_name):
            # self.load_nastran_ressults(op2_filename, None)
            dirname = None
            getattr(self, load_results_name)(output_filename, dirname)
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
    msg += "  test_pynastrangui [-f FORMAT]           INPUT_FILENAME  OUTPUT_FILENAME [--log LOG]\n"
    msg += "  test_pynastrangui [-f FORMAT]           INPUT_FILENAME  [--log LOG]\n"
    msg += "  test_pynastrangui  -f FORMAT  [-r] [-d] INPUT_DIRECTORY [--log LOG]\n"
    msg += "  test_pynastrangui  -f FORMAT  [-r] [-d]                 [--log LOG]\n"

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

    msg += 'Info:\n'
    #msg += "  -q, --quiet    prints debug messages (default=True)\n"
    msg += '  -h, --help     show this help message and exit\n'
    msg += "  -v, --version  show program's version number and exit\n"

    if len(sys.argv) == 1:
        sys.exit(msg)
    ver = str(pyNastran.__version__)
    data = docopt(msg, argv=argv, help=True, version=ver, options_first=False)

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
        if not os.path.exists(input_filename):
            msg = 'input_filename=%r does not exist\n%s' % (
                input_filename, print_bad_path(input_filename))
            raise RuntimeError(msg)
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
    return formati, input_filenames, output_filenames, failed_cases_filename, log_method

def main():
    """runs the gui"""
    formati, input_filenames, output_filenames, failed_cases_filename, log_method = run_docopt()
    log = get_logger(log=None, level=log_method, encoding='utf-8')
    npass = 0
    nfailed = 0
    failed_files = []
    ntotal = len(input_filenames)

    test = FakeGUI(formati)
    test.log = log
    stop_on_failure = ntotal == 1
    t0 = time.time()
    for input_filename, output_filename in zip(input_filenames, output_filenames):
        input_filename = os.path.abspath(input_filename)
        #output_filename =
        print("filename = %s" % input_filename)
        is_passed = True
        try:
            test.load_geometry(input_filename)
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
            test.load_results(output_filename)

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
    dt = time.time() - t0
    nfailed = len(failed_files)
    npassed = ntotal - nfailed
    time_msg = 'dt = %.0f sec = %.1f min' % (dt, dt/60.)
    sys.stderr.write('%i/%i passed\n' % (npassed, ntotal))
    sys.stderr.write(time_msg)

    if PY2:
        write = 'wb'
    else:
        write = 'w'

    if ntotal > 1:
        with open(failed_cases_filename, write) as failed_cases_file:
            for fname in failed_files:
                failed_cases_file.write('%s\n' % fname)
        print(time_msg)
        sys.exit('finished...')

    #test.load_nastran_geometry(bdf_filename, None)
    #test.load_nastran_results(op2_filename, None)

if __name__ == '__main__':
    main()
