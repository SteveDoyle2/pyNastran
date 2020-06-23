"""Defines the command line tool `test_op2`"""
import os
import sys
import time
from traceback import print_exc
from typing import List, Tuple, Optional, Any

import numpy as np
np.set_printoptions(precision=3, threshold=20)

try:
    import pandas
    IS_PANDAS = True
except ImportError:
    IS_PANDAS = False

try:
    import h5py
    IS_HDF5 = True
except ImportError:
    IS_HDF5 = False

#import warnings
#warnings.filterwarnings('error')
#warnings.filterwarnings('error', category=UnicodeWarning)

import pyNastran
from pyNastran.op2.op2 import (
    OP2, FatalError, SixtyFourBitError, OverwriteTableError)
#SortCodeError, DeviceCodeError, FortranMarkerError

from pyNastran.op2.op2_geom import OP2Geom, DuplicateIDsError


def parse_table_names_from_f06(f06_filename):
    """gets the op2 names from the f06"""

    marker = 'NAME OF DATA BLOCK WRITTEN ON FORTRAN UNIT IS'
    names = []
    with open(f06_filename, 'r') as infile:
        for line in infile:
            if marker in line:
                word = line.replace(marker, '').strip().strip('.')
                names.append(word)
    return names


def run_lots_of_files(files, make_geom: bool=True, combine: bool=True,
                      write_bdf: bool=False, write_f06: bool=True,
                      delete_f06: bool=True, delete_op2: bool=True, delete_hdf5: bool=True,
                      delete_debug_out: bool=True, build_pandas: bool=True, write_op2: bool=False,
                      write_hdf5: bool=True, debug: bool=True, skip_files: Optional[List[str]]=None,
                      stop_on_failure: bool=False, nstart: int=0, nstop: int=1000000000,
                      short_stats: bool=False, binary_debug: bool=False,
                      compare: bool=True, quiet: bool=False, dev: bool=True, xref_safe: bool=False):
    """used by op2_test.py to run thousands of files"""
    if skip_files is None:
        skip_files = []
    #n = ''
    assert make_geom in [True, False]
    assert combine in [True, False]
    assert write_bdf in [True, False]
    assert write_f06 in [True, False]
    assert write_op2 in [True, False]
    assert write_hdf5 in [True, False]
    assert build_pandas in [True, False]
    if binary_debug in [True, False]:
        binary_debug = [binary_debug]

    subcases = []
    failed_cases = []
    nfailed = 0
    ntotal = 0
    npassed = 0
    #t0 = time.time()
    for i, op2file in enumerate(files[nstart:nstop], nstart):  # 149
        basename = os.path.basename(op2file)
        #if basename not in skip_files and not basename.startswith('acms') and i not in nskip:
        sys.stderr.write('%s file=%s\n' % (i, op2file))
        if basename not in skip_files and '#' not in op2file:
            print("%" * 80)
            print('file=%s\n' % op2file)
            #n = '%s ' % i
            ntotal += 1

            is_passed = True
            for binary_debugi in binary_debug:
                print('------running binary_debug=%s------' % binary_debugi)
                is_passedi = run_op2(op2file, make_geom=make_geom, combine=combine,
                                     write_bdf=write_bdf, write_f06=write_f06, write_op2=write_op2,
                                     is_mag_phase=False,
                                     delete_f06=delete_f06,
                                     delete_op2=delete_op2,
                                     delete_hdf5=delete_hdf5,
                                     delete_debug_out=delete_debug_out,
                                     build_pandas=build_pandas,
                                     write_hdf5=write_hdf5,
                                     short_stats=short_stats,
                                     subcases=subcases, debug=debug,
                                     stop_on_failure=stop_on_failure,
                                     binary_debug=binary_debug,
                                     compare=compare, dev=dev,
                                     xref_safe=xref_safe)[1]
                if not is_passedi:
                    is_passed = False
                    break

            if not is_passed:
                sys.stderr.write('**file=%s\n' % op2file)
                failed_cases.append(op2file)
                nfailed += 1
            else:
                npassed += 1
    return failed_cases


def run_op2(op2_filename: str, make_geom: bool=False, combine: bool=True,
            write_bdf: bool=False, read_bdf: Optional[bool]=None,
            write_f06: bool=True, write_op2: bool=False,
            write_hdf5: bool=True,
            is_mag_phase: bool=False, is_sort2: bool=False,
            is_nx: Optional[bool]=None, is_autodesk: Optional[bool]=None,
            is_nasa95: Optional[bool]=None,
            delete_f06: bool=False, delete_op2: bool=False, delete_hdf5: bool=False,
            delete_debug_out: bool=False,
            build_pandas: bool=True,
            subcases: Optional[str]=None, exclude: Optional[str]=None,
            short_stats: bool=False, compare: bool=True,
            debug: bool=False, log: Any=None,
            binary_debug: bool=False, quiet: bool=False,
            stop_on_failure: bool=True,
            dev: bool=False, xref_safe: bool=False,
            post: Any=None, load_as_h5: bool=False) -> Tuple[OP2, bool]:
    """
    Runs an OP2

    Parameters
    ----------
    op2_filename : str
        path of file to test
    make_geom : bool; default=False
        should the GEOMx, EPT, MPT, DYNAMIC, DIT, etc. tables be read
    combine : bool; default=True
        should the op2 tables be combined
    write_bdf : bool; default=False
        should a BDF be written based on the geometry tables
    write_f06 : bool; default=True
        should an F06 be written based on the results
    write_op2 : bool; default=False
        should an OP2 be written based on the results
    write_hdf5 : bool; default=True
        should an HDF5 be written based on the results
    is_mag_phase : bool; default=False
        False : write real/imag results
        True : write mag/phase results
        For static results, does nothing
    is_sort2 : bool; default=False
        False : writes "transient" data is SORT1
        True : writes "transient" data is SORT2
    is_nx : bool; default=None
        True : use NX Nastran
        False : use MSC Nastran
        None : guess
    is_autodesk : bool; default=None
        True : use Autodesk Nastran
        False : use MSC Nastran
        None : guess
    is_nasa95 : bool; default=None
        True : use NASA 95 Nastran
        False : use MSC Nastran
        None : guess
    delete_f06 : bool; default=False
        deletes the F06 (assumes write_f06 is True)
    delete_op2 : bool; default=False
        deletes the OP2 (assumes write_op2 is True)
    delete_hdf5 : bool; default=False
        deletes the HDF5 (assumes write_hdf5 is True)
    subcases : List[int, ...]; default=None
        limits subcases to specified values; default=None -> no limiting
    exclude : List[str, ...]; default=None
        limits result types; (remove what's listed)
    short_stats : bool; default=False
        print a short version of the op2 stats
    compare : bool
        True : compares vectorized result to slow vectorized result
        False : doesn't run slow vectorized result
    debug : bool; default=False
        debug flag for OP2
    log : logger; default=None
        a custom logger
        None : use debug
    binary_debug : bool; default=False
        creates a very cryptic developer debug file showing exactly what was parsed
    quiet : bool; default=False
        don't write debug messages
    stop_on_failure : bool; default=True
        is this used???
    dev : bool; default=False
        flag that is used by op2_test.py to ignore certain errors
        False : crash on errors
        True : don't crash

    Returns
    -------
    op2 : OP2()
        the op2 object
    is_passed : bool
        did the test pass

    """
    assert build_pandas in [True, False]

    if read_bdf is None:
        read_bdf = write_bdf
    op2 = None
    op2_nv = None
    if subcases is None:
        subcases = []
    if exclude is None:
        exclude = []
    if isinstance(is_sort2, bool):
        sort_methods = [is_sort2]
    else:
        sort_methods = is_sort2

    assert '.op2' in op2_filename.lower(), 'op2_filename=%s is not an OP2' % op2_filename
    is_passed = False

    fname_base = os.path.splitext(op2_filename)[0]
    bdf_filename = fname_base + '.test_op2.bdf'

    if isinstance(subcases, str):
        if '_' in subcases:
            subcases = [int(i) for i in subcases.split('_')]
        else:
            subcases = [int(subcases)]

    debug_file = None
    model = os.path.splitext(op2_filename)[0]
    if binary_debug or write_op2:
        debug_file = model + '.debug.out'
    #print('debug_file = %r' % debug_file, os.getcwd())

    if make_geom:
        op2 = OP2Geom(debug=debug, log=log)
        op2_nv = OP2Geom(debug=debug, log=log, debug_file=debug_file)
        op2_bdf = OP2Geom(debug=debug, log=log)
        set_versions([op2, op2_nv, op2_bdf], is_nx, is_autodesk, is_nasa95, post)

        if load_as_h5 and IS_HDF5:
            # you can't open the same h5 file twice
            op2.load_as_h5 = load_as_h5
            #op2_nv.load_as_h5 = load_as_h5
            #op2_bdf.load_as_h5 = load_as_h5

        op2_bdf.set_error_storage(nparse_errors=0, stop_on_parsing_error=True,
                                  nxref_errors=0, stop_on_xref_error=True)
    else:
        op2 = OP2(debug=debug, log=log)
        # have to double write this until ???
        op2_nv = OP2(debug=debug, log=log, debug_file=debug_file)

        set_versions([op2, op2_nv], is_nx, is_autodesk, is_nasa95, post)
        if load_as_h5 and IS_HDF5:
            # you can't open the same h5 file twice
            op2.load_as_h5 = load_as_h5
            #op2_nv.load_as_h5 = load_as_h5
        op2_bdf = None
    op2_nv.use_vector = False

    if not quiet:
        op2.log.debug('subcases = %s' % subcases)
    op2.set_subcases(subcases)
    op2_nv.set_subcases(subcases)
    op2.remove_results(exclude)
    op2_nv.remove_results(exclude)

    try:
        #op2.read_bdf(op2.bdf_filename, includeDir=None, xref=False)
        if compare:
            op2_nv.read_op2(op2_filename, combine=combine)
        op2.read_op2(op2_filename, combine=combine)
        #if not make_geom:  # TODO: enable this...
            #op2.save()

        #op2a.get_op2_stats()
        op2.get_op2_stats()
        op2.get_op2_stats(short=True)
        op2.object_attributes()
        op2.object_methods()
        if not quiet:
            print("---stats for %s---" % op2_filename)
            print(op2.get_op2_stats(short=short_stats))
            op2.print_subcase_key()

        write_op2_as_bdf(op2, op2_bdf, bdf_filename, write_bdf, make_geom, read_bdf, dev,
                         xref_safe=xref_safe)

        if compare:
            assert op2 == op2_nv

        if IS_HDF5 and write_hdf5:
            from pyNastran.op2.op2_interface.hdf5_interface import load_op2_from_hdf5_filename
            h5_filename = model + '.test_op2.h5'
            op2.export_hdf5_filename(h5_filename)
            load_op2_from_hdf5_filename(h5_filename, log=op2.log)
            if delete_hdf5:
                remove_file(h5_filename)
        if write_f06:
            for is_sort2i in sort_methods:
                f06_filename = model + '.test_op2.f06'
                op2.write_f06(f06_filename, is_mag_phase=is_mag_phase,
                              is_sort1=not is_sort2i, quiet=quiet, repr_check=True)
            if delete_f06:
                remove_file(f06_filename)

        # we put it down here so we don't blame the dataframe for real errors
        if IS_PANDAS and build_pandas:
            op2.build_dataframe()
        #if compare:
            #op2_nv.build_dataframe()

        if write_op2:
            model = os.path.splitext(op2_filename)[0]
            op2_filename2 = model + '.test_op2.op2'
            total_case_count = op2.write_op2(op2_filename2,
                                             #is_mag_phase=is_mag_phase,
                                             endian=b'<')
            if total_case_count > 0:
                #print('------------------------------')
                op2a = OP2(debug_file='debug.out', log=log)
                op2a.use_vector = False
                op2a.read_op2(op2_filename2)
                os.remove(op2_filename2)
            #read_op2(op2_filename2)
            if delete_op2:
                remove_file(op2_filename2)

        if debug_file is not None and delete_debug_out:
            os.remove(debug_file)

        #table_names_f06 = parse_table_names_from_F06(op2.f06FileName)
        #table_names_op2 = op2.getTableNamesFromOP2()
        #print("subcases = ", op2.subcases)

        #if table_names_f06 != table_names_op2:
            #msg = 'table_names_f06=%s table_names_op2=%s' % (table_names_f06, table_names_op2)
            #raise RuntimeError(msg)
        #op2.case_control_deck.sol = op2.sol
        #print(op2.case_control_deck.get_op2_data())
        #print(op2.case_control_deck.get_op2_data())
        is_passed = True
    except MemoryError:
        raise
    except KeyboardInterrupt:
        sys.stdout.flush()
        print_exc(file=sys.stdout)
        sys.stderr.write('**file=%s\n' % op2_filename)
        sys.exit('keyboard stop...')
    #except SortCodeError: # inherits from Runtime; comment this
        #is_passed = True

    #except RuntimeError: # the op2 is bad, not my fault; comment this
        #is_passed = True
        #if stop_on_failure:
            #raise
        #else:
            #is_passed = True
    #except RuntimeError:
        #pass
    #except ValueError:
        #pass
    #except IndexError:
        #pass
    #except FortranMarkerError:
        #pass
    except IOError: # missing file; this block should be uncommented
        #if stop_on_failure:
            #raise
        if not dev:
            raise
        is_passed = True
    #except UnicodeDecodeError:  # this block should be commented
        #is_passed = True
    #except NotImplementedError:  # this block should be commented
        #is_passed = True
    except SixtyFourBitError:
        #log.error('SixtyFourBitError')
        if not dev:
            raise
        is_passed = True
    except FatalError:  # this block should be commented
        #if stop_on_failure:
            #raise
        if not dev:
            raise
        is_passed = True
    #except KeyError:  # this block should be commented
        #is_passed = True
    #except DeviceCodeError:  # this block should be commented
        #is_passed = True

    #except AssertionError:  # this block should be commented
        #is_passed = True
    #except RuntimeError: #invalid analysis code; this block should be commented
        #is_passed = True
    #except ValueError:  # this block should be commented
        #is_passed = True
    #except NotImplementedError:  # this block should be commented
        #is_passed = True
    #except FortranMarkerError:  # this block should be commented
        #is_passed = True

    except DuplicateIDsError:
        if not dev:
            raise
        is_passed = True

    except SystemExit:
        #print_exc(file=sys.stdout)
        #sys.exit('stopping on sys.exit')
        raise
    #except NameError:  # variable isnt defined
    #    if stop_on_failure:
    #        raise
    #    else:
    #        is_passed = True
    #except IndexError: # this block should be commented
        #is_passed = True
    #except SyntaxError: #Param Parse; this block should be commented
        #if stop_on_failure:
            #raise
        #is_passed = True
    except:
        #print(e)
        if stop_on_failure:
            raise
        else:
            print_exc(file=sys.stdout)
            is_passed = False

    return op2, is_passed

def write_op2_as_bdf(op2, op2_bdf, bdf_filename, write_bdf, make_geom, read_bdf, dev,
                     xref_safe=False):
    if write_bdf:
        assert make_geom, f'write_bdf=False, but make_geom={make_geom!r}; expected make_geom=True'
        op2._nastran_format = 'msc'
        op2.executive_control_lines = ['CEND\n']
        op2.validate()
        op2.write_bdf(bdf_filename, size=8)
        op2.log.debug('bdf_filename = %s' % bdf_filename)
        xref = xref_safe is False
        if read_bdf:
            try:
                op2_bdf.read_bdf(bdf_filename, xref=xref)
                if xref_safe:
                    op2_bdf.safe_cross_reference()
            except:
                if dev and len(op2_bdf.card_count) == 0:
                    pass
                else:
                    raise
        #os.remove(bdf_filename)

def get_test_op2_data(argv):
    """defines the docopt interface"""
    from docopt import docopt
    ver = str(pyNastran.__version__)
    is_dev = 'dev' in ver

    msg = "Usage:  "
    nasa95 = '|--nasa95' if is_dev else ''
    version = f'[--nx|--autodesk{nasa95}]'
    options = f'[-p] [-d] [-z] [-w] [-t] [-s <sub>] [-x <arg>]... {version} [--safe] [--post POST] [--load_hdf5]'
    if is_dev:
        line1 = f"test_op2 [-q] [-b] [-c] [-g] [-n] [-f] [-o] [--profile] [--nocombine] {options} OP2_FILENAME\n"
    else:
        line1 = f"test_op2 [-q] [-b] [-c] [-g] [-n] [-f] [-o] {options} OP2_FILENAME\n"

    while '  ' in line1:
        line1 = line1.replace('  ', ' ')
    msg += line1
    msg += "        test_op2 -h | --help\n"
    msg += "        test_op2 -v | --version\n"
    msg += "\n"
    msg += "Tests to see if an OP2 will work with pyNastran %s.\n" % ver
    msg += "\n"
    msg += "Positional Arguments:\n"
    msg += "  OP2_FILENAME         Path to OP2 file\n"
    msg += "\n"
    msg += "Options:\n"
    msg += "  -b, --binarydebug      Dumps the OP2 as a readable text file\n"
    msg += "  -c, --disablecompare   Doesn't do a validation of the vectorized result\n"
    msg += "  -q, --quiet            Suppresses debug messages [default: False]\n"
    msg += "  -t, --short_stats      Short get_op2_stats printout\n"
    #if not is_release:
    msg += "  -g, --geometry         Reads the OP2 for geometry, which can be written out\n"
    # n is for NAS
    msg += "  -n, --write_bdf        Writes the bdf to fem.test_op2.bdf (default=False)\n"
    msg += "  -f, --write_f06        Writes the f06 to fem.test_op2.f06\n"
    msg += "  -d, --write_hdf5       Writes the h5 to fem.test_op2.h5\n"
    msg += "  -o, --write_op2        Writes the op2 to fem.test_op2.op2\n"
    msg += "  -z, --is_mag_phase     F06 Writer writes Magnitude/Phase instead of\n"
    msg += "                         Real/Imaginary (still stores Real/Imag); [default: False]\n"
    msg += "  --load_hdf5            Load as HDF5 (default=False)\n"
    msg += "  -p, --pandas           Enables pandas dataframe building; [default: False]\n"
    if is_dev:
        msg += "  --nocombine            Disables case combination\n"
    msg += "  -s <sub>, --subcase    Specify one or more subcases to parse; (e.g. 2_5)\n"
    msg += "  -w, --is_sort2         Sets the F06 transient to SORT2\n"
    msg += "  -x <arg>, --exclude    Exclude specific results\n"
    msg += "  --nx                   Assume NX Nastran\n"
    msg += "  --autodesk             Assume Autodesk Nastran\n"
    if is_dev:
        msg += "  --nasa95               Assume Nastran 95\n"
    msg += "  --post POST            Set the PARAM,POST flag\n"
    msg += "  --safe                 Safe cross-references BDF (default=False)\n"

    if is_dev:
        msg += "\n"
        msg += "Developer:\n"
        msg += '  --profile         Profiles the code (default=False)\n'

    msg += "\n"
    msg += "Info:\n"
    msg += "  -h, --help     Show this help message and exit\n"
    msg += "  -v, --version  Show program's version number and exit\n"

    if len(argv) == 1:
        sys.exit(msg)

    data = docopt(msg, version=ver, argv=argv[1:])
    if not is_dev:
        # just set the defaults for these so we don't need special code later
        data['--profile'] = False
        data['--write_xlsx'] = False
        data['--nocombine'] = False
        data['--nasa95'] = False

    if '--geometry' not in data:
        data['--geometry'] = False
    if '--write_bdf' not in data:
        data['--write_bdf'] = False
    data['--is_sort2'] = bool(data['--is_sort2'])
    #print("data", data)
    return data

def remove_file(filename):
    try:
        os.remove(filename)
    except:
        pass


def set_versions(op2s, is_nx, is_autodesk, is_nasa95, post):
    for op2 in op2s:
        op2.IS_TESTING = False

    if is_nx is None and is_autodesk is None and is_nasa95 is None:
        pass
    elif is_nx:
        for op2 in op2s:
            op2.set_as_nx()
    elif is_autodesk:
        for op2 in op2s:
            op2.set_as_autodesk()
    elif is_nasa95:
        for op2 in op2s:
            op2.set_as_nasa95()
    else:
        for op2 in op2s:
            op2.set_as_msc()

    if post is not None:
        for op2 in op2s:
            op2.post = -4

def main(argv=None, show_args=True):
    """the interface for test_op2"""
    if argv is None:
        argv = sys.argv
    data = get_test_op2_data(argv)

    if show_args:
        for key, value in sorted(data.items()):
            print("%-12s = %r" % (key.strip('--'), value))

    if os.path.exists('skippedCards.out'):
        os.remove('skippedCards.out')

    time0 = time.time()

    if data['--profile']:
        import pstats
        import cProfile
        prof = cProfile.Profile()
        prof.runcall(
            run_op2,
            data['OP2_FILENAME'],
            make_geom=data['--geometry'],
            combine=not data['--nocombine'],
            load_as_h5=data['--load_hdf5'],
            write_bdf=data['--write_bdf'],
            write_f06=data['--write_f06'],
            write_op2=data['--write_op2'],
            write_hdf5=data['--write_hdf5'],
            is_mag_phase=data['--is_mag_phase'],
            build_pandas=data['--pandas'],
            subcases=data['--subcase'],
            exclude=data['--exclude'],
            debug=not data['--quiet'],
            binary_debug=data['--binarydebug'],
            is_sort2=data['--is_sort2'],
            compare=not data['--disablecompare'],
            quiet=data['--quiet'],
            is_nx=data['--nx'],
            is_autodesk=data['--autodesk'],
            is_nasa95=data['--nasa95'],
            safe=data['--safe'],
            post=data['--post'],
        )
        prof.dump_stats('op2.profile')

        stats = pstats.Stats("op2.profile")
        stats.sort_stats('tottime')  # time in function
        #stats.sort_stats('cumtime')  # time in function & subfunctions
        stats.strip_dirs()
        stats.print_stats(40)
    else:
        run_op2(
            data['OP2_FILENAME'],
            make_geom=data['--geometry'],
            combine=not data['--nocombine'],
            load_as_h5=data['--load_hdf5'],
            write_bdf=data['--write_bdf'],
            write_f06=data['--write_f06'],
            write_op2=data['--write_op2'],
            write_hdf5=data['--write_hdf5'],
            is_mag_phase=data['--is_mag_phase'],
            build_pandas=data['--pandas'],
            subcases=data['--subcase'],
            exclude=data['--exclude'],
            short_stats=data['--short_stats'],
            debug=not data['--quiet'],
            binary_debug=data['--binarydebug'],
            is_sort2=data['--is_sort2'],
            compare=not data['--disablecompare'],
            quiet=data['--quiet'],
            is_nx=data['--nx'],
            is_autodesk=data['--autodesk'],
            is_nasa95=data['--nasa95'],
            xref_safe=data['--safe'],
            post=data['--post'],
        )
    print("dt = %f" % (time.time() - time0))

if __name__ == '__main__':  # pragma: no cover
    main(show_args=True)
