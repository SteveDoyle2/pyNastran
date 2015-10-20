"""
Defines the command line tool `test_op2`
"""
from __future__ import print_function
from six import string_types, iteritems, PY2
import os
import sys
import time
from traceback import print_exc

import pyNastran
from pyNastran import is_release
from pyNastran.f06.errors import FatalError
from pyNastran.op2.op2 import OP2, FatalError, SortCodeError, DeviceCodeError

try:
    #from pyNastran.op2.op2_geom import OP2Geom
    is_geom = True
except ImportError:
    is_geom = False
    #raise


# we need to check the memory usage
is_linux = None
is_memory = True
try:
    if os.name == 'nt':  # windows
        windows_flag = True
        is_linux = False
        import wmi
        comp = wmi.WMI()

        """Functions for getting memory usage of Windows processes."""

        __all__ = ['get_current_process', 'get_memory_info', 'get_memory_usage']

        import ctypes
        from ctypes import wintypes

        GetCurrentProcess = ctypes.windll.kernel32.GetCurrentProcess
        GetCurrentProcess.argtypes = []
        GetCurrentProcess.restype = wintypes.HANDLE

        SIZE_T = ctypes.c_size_t

        class PROCESS_MEMORY_COUNTERS_EX(ctypes.Structure):
            """
            Windows memory tool that has the same interface as
            `resource` on Linux/Mac.
            """
            _fields_ = [
                ('cb', wintypes.DWORD),
                ('PageFaultCount', wintypes.DWORD),
                ('PeakWorkingSetSize', SIZE_T),
                ('WorkingSetSize', SIZE_T),
                ('QuotaPeakPagedPoolUsage', SIZE_T),
                ('QuotaPagedPoolUsage', SIZE_T),
                ('QuotaPeakNonPagedPoolUsage', SIZE_T),
                ('QuotaNonPagedPoolUsage', SIZE_T),
                ('PagefileUsage', SIZE_T),
                ('PeakPagefileUsage', SIZE_T),
                ('PrivateUsage', SIZE_T),
            ]

        GetProcessMemoryInfo = ctypes.windll.psapi.GetProcessMemoryInfo
        GetProcessMemoryInfo.argtypes = [
            wintypes.HANDLE,
            ctypes.POINTER(PROCESS_MEMORY_COUNTERS_EX),
            wintypes.DWORD,
        ]
        GetProcessMemoryInfo.restype = wintypes.BOOL

        def get_current_process():
            """Return handle to current process."""
            return GetCurrentProcess()

        def get_memory_info(process=None):
            """Return Win32 process memory counters structure as a dict."""
            if process is None:
                process = get_current_process()
            counters = PROCESS_MEMORY_COUNTERS_EX()
            ret = GetProcessMemoryInfo(process, ctypes.byref(counters),
                                       ctypes.sizeof(counters))
            if not ret:
                raise ctypes.WinError()
            info = dict((name, getattr(counters, name))
                        for name, _ in counters._fields_)
            return info

        def get_memory_usage(process=None):
            """Return this process's memory usage in bytes."""
            info = get_memory_info(process=process)
            return info['PrivateUsage']

    elif os.name in ['posix', 'mac']:  # linux/mac
        import resource
        windows_flag = False
        is_linux = True
except:
    is_memory = False


def parse_table_names_from_F06(f06_filename):
    """gets the op2 names from the f06"""
    infile = open(f06_filename, 'r')
    marker = 'NAME OF DATA BLOCK WRITTEN ON FORTRAN UNIT IS'
    names = []
    for line in infile:
        if marker in line:
            word = line.replace(marker, '').strip().strip('.')
            names.append(word)
    infile.close()
    return names


def get_failed_files(filename):
    """Gets the list of failed files"""
    infile = open(filename, 'r')
    lines = infile.readlines()
    infile.close()

    files = []
    for line in lines:
        files.append(line.strip())
    return files


def run_lots_of_files(files, make_geom=True, write_bdf=False, write_f06=True,
                      delete_f06=True, write_op2=False,
                      is_vector=False, vector_stop=True,
                      debug=True, saveCases=True, skip_files=None,
                      stop_on_failure=False, nstart=0, nstop=1000000000, binary_debug=False,
                      compare=True):
    """used by op2_test.py to run thousands of files"""
    if skip_files is None:
        skip_files = []
    n = ''
    assert make_geom in [True, False]
    assert write_bdf in [True, False]
    assert write_f06 in [True, False]
    if is_vector in [True, False]:
        is_vector = [is_vector]
        vector_stop = [vector_stop]

    isubcases = []
    failed_cases = []
    nfailed = 0
    ntotal = 0
    npassed = 0
    t0 = time.time()
    for i, op2file in enumerate(files[nstart:nstop], nstart):  # 149
        basename = os.path.basename(op2file)
        #if baseName not in skip_files and not baseName.startswith('acms') and i not in nSkip:
        if basename not in skip_files and '#' not in op2file:
            print("%" * 80)
            print('file=%s\n' % op2file)
            n = '%s ' % i
            sys.stderr.write('%sfile=%s\n' % (n, op2file))
            ntotal += 1

            is_passed = True
            is_vector_failed = []
            for vectori, vector_stopi in zip(is_vector, vector_stop):
                print('------running is_vector=%s------' % vectori)
                is_passed_i = run_op2(op2file, make_geom=make_geom, write_bdf=write_bdf,
                                      write_f06=write_f06, write_op2=write_op2,
                                      is_mag_phase=False,
                                      delete_f06=delete_f06,
                                      isubcases=isubcases, debug=debug,
                                      stop_on_failure=stop_on_failure,
                                      binary_debug=binary_debug,
                                      compare=True)[1]
                if not is_passed_i and vector_stopi:
                    is_passed = False
                if not is_passed_i:
                    is_vector_failed.append(vectori)
            if not is_passed:
                sys.stderr.write('**file=%s vector_failed=%s\n' % (op2file, is_vector_failed))
                failed_cases.append(op2file)
                nfailed += 1
            else:
                npassed += 1
            #sys.exit('end of test...test_op2.py')

    if saveCases:
        if PY2:
            failed_cases_file = open('failedCases.in', 'wb')
        else:
            failed_cases_file = open('failedCases.in', 'w')
        for op2file in failed_cases:
            failed_cases_file.write('%s\n' % op2file)
        failed_cases_file.close()

    seconds = time.time() - t0
    minutes = seconds / 60.
    print("dt = %s seconds = %s minutes" % (seconds, minutes))

    msg = '-----done with all models %s/%s=%.2f%%  nFailed=%s-----' % (
        npassed, ntotal,
        100. * npassed / float(ntotal),
        ntotal - npassed)
    print(msg)
    sys.exit(msg)


def run_op2(op2_filename, make_geom=False, write_bdf=False,
            write_f06=True, write_op2=False, is_mag_phase=False, is_sort2=False,
            delete_f06=False,
            isubcases=None, exclude=None, compare=True, debug=False, binary_debug=False,
            quiet=False, stop_on_failure=True):
    """
    Runs an OP2

    Parameters
    ----------
    op2_filename : str
        path of file to test
    make_geom : bool; default=False
        should the GEOMx, EPT, MPT, DYNAMIC, DIT, etc. tables be read
    write_bdf : bool; default=False
        should a BDF be written based on the geometry tables
    write_f06 : bool; default=True
        should an F06 be written based on the results
    is_mag_phase : bool; default=False
        False : write real/imag results
        True : write mag/phase results
        For static results, does nothing
    is_sort2 : bool; default=False
        False : writes "transient" data is SORT1
        True : writes "transient" data is SORT2
    delete_f06 : bool; default=False
        deletes the F06 (assumes write_f06 is True)
    isubcases : List[int, ...]; default=None
        limits subcases to specified values; default=None -> no limiting
    exclude : List[str, ...]; default=None
        limits result types; (remove what's listed)
    compare : bool
        True : compares vectorized result to slow vectorized result
        False : doesn't run slow vectorized result
    debug : bool; default=False
        dunno???
    binary_debug : bool; default=False
        creates a very cryptic developer debug file showing exactly what was parsed
    quiet : bool; default=False
        dunno???
    stopOnFailure : bool; default=True
        is this used???
    """
    op2a = None
    op2b = None
    if isubcases is None:
        isubcases = []
    if exclude is None:
        exclude = []
    assert '.op2' in op2_filename.lower(), 'op2_filename=%s is not an OP2' % op2_filename
    is_passed = False

    fname_base = os.path.splitext(op2_filename)[0]
    bdf_filename = fname_base + '.test_op2.bdf'

    if isinstance(isubcases, string_types):
        if '_' in isubcases:
            isubcases = [int(i) for i in isubcases.split('_')]
        else:
            isubcases = [int(isubcases)]
    print('isubcases = %s' % isubcases)

    debug_file = None
    model = os.path.splitext(op2_filename)[0]
    if binary_debug or write_op2:
        debug_file = model + '.debug.out'
    #print('debug_file = %r' % debug_file, os.getcwd())

    if make_geom and not is_geom:
        raise RuntimeError('make_geom=%s is not supported' % make_geom)
    if make_geom:
        op2a = OP2Geom(debug=debug, debug_file=debug_file)
        op2b = OP2Geom()
    else:
        op2a = OP2(debug=debug, debug_file=debug_file)
        op2b = OP2()
    op2b.use_vector = False

    op2a.set_subcases(isubcases)
    op2b.set_subcases(isubcases)
    op2a.remove_results(exclude)
    op2b.remove_results(exclude)

    if is_memory:
        if is_linux: # linux
            kb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        else: # windows
            kb = get_memory_usage() / 1024
        mb = kb / 1024.
        #mbs.append(mb)
        print("Memory usage start: %s (KB); %.2f (MB)" % (kb, mb))

    try:
        #op2.read_bdf(op2.bdf_filename, includeDir=None, xref=False)
        op2a.read_op2(op2_filename)
        if compare:
            op2b.read_op2(op2_filename)

        #op2a.get_op2_stats()
        if quiet:
            op2a.get_op2_stats()
        else:
            print("---stats for %s---" % op2_filename)
            print(op2a.get_op2_stats())
            print(op2a.print_subcase_key())
        if write_bdf:
            op2a.write_bdf(bdf_filename)
            os.remove(bdf_filename)
        if compare:
            assert op2a == op2b

        if is_memory:
            if is_linux: # linux
                kb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
            else: # windows
                kb = get_memory_usage() / 1024
            mb = kb / 1024.
            #mbs.append(mb)
            print("Memory usage     end: %s (KB); %.2f (MB)" % (kb, mb))

        if write_f06:
            op2a.write_f06(model + '.test_op2.f06', is_mag_phase=is_mag_phase, is_sort1=not is_sort2)
            if delete_f06:
                try:
                    os.remove(model + '.test_op2.f06')
                except:
                    pass

        if write_op2:
            model = os.path.splitext(op2_filename)[0]
            op2a.write_op2(model + '.test_op2.op2', is_mag_phase=is_mag_phase)
            if delete_f06:
                try:
                    os.remove(model + '.test_op2.op2')
                except:
                    pass

        if is_memory:
            del op2a
            del op2b
            if is_linux: # linux
                kb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
            else: # windows
                kb = get_memory_usage() / 1024
            mb = kb / 1024.
            #mbs.append(mb)
            print("Memory usage cleanup: %s (KB); %.2f (MB)" % (kb, mb))


        #table_names_f06 = parse_table_names_from_F06(op2.f06FileName)
        #table_names_op2 = op2.getTableNamesFromOP2()
        #print("subcases = ", op2.subcases)

        #assert table_names_f06==table_names_op2, 'table_names_f06=%s table_names_op2=%s' % (table_names_f06, table_names_op2)
        #op2.case_control_deck.sol = op2.sol
        #print(op2.case_control_deck.get_op2_data())
        #print(op2.case_control_deck.get_op2_data())
        is_passed = True
    except KeyboardInterrupt:
        sys.stdout.flush()
        print_exc(file=sys.stdout)
        sys.stderr.write('**file=%s\n' % op2_filename)
        sys.exit('keyboard stop...')
    #except SortCodeError: # inherits from Runtime; comment this
        #is_passed = True

    #except RuntimeError: # the op2 is bad, not my fault; comment this
        #is_passed = True
        #if stopOnFailure:
            #raise
        #else:
            #is_passed = True

    except IOError: # missing file; this block should be commented
        #if stopOnFailure:
            #raise
        is_passed = True
    #except UnicodeDecodeError:  # this block should be commented
        #is_passed = True
    #except NotImplementedError:  # this block should be commented
        #is_passed = True
    except FatalError:  # this block should be commented
        #if stopOnFailure:
            #raise
        is_passed = True
    #except KeyError:  # this block should be commented
        #is_passed = True
    #except DeviceCodeError:  # this block should be commented
        #is_passed = True
    #except AssertionError:  # this block should be commented
        #is_passed = True
    #except RuntimeError: #invalid analysis code; this block should be commented
        #is_passed = True
    except SystemExit:
        #print_exc(file=sys.stdout)
        #sys.exit('stopping on sys.exit')
        raise
    #except NameError:  # variable isnt defined
    #    if stopOnFailure:
    #        raise
    #    else:
    #        is_passed = True
    #except IndexError: # this block should be commented
        #is_passed = True
    #except SyntaxError: #Param Parse; this block should be commented
        #if stopOnFailure:
            #raise
        #is_passed = True
    except:
        #print e
        if stop_on_failure:
            raise
        else:
            print_exc(file=sys.stdout)
            is_passed = False

    return op2a, is_passed


def main():
    """the interface for test_op2"""
    from docopt import docopt
    ver = str(pyNastran.__version__)

    msg = "Usage:\n"
    is_release = False
    if is_release:
        line1 = "test_op2 [-q] [-b] [-c]           [-f]      [-z] [-w] [-s <sub>] [-x <arg>]... OP2_FILENAME\n"
    else:
        line1 = "test_op2 [-q] [-b] [-c] [-g] [-w] [-f] [-o] [-z] [-w] [-s <sub>] [-x <arg>]... OP2_FILENAME\n"

    while '  ' in line1:
        line1 = line1.replace('  ', ' ')
    msg += line1
    msg += "  test_op2 -h | --help\n"
    msg += "  test_op2 -v | --version\n"
    msg += "\n"
    msg += "Tests to see if an OP2 will work with pyNastran %s.\n" % ver
    msg += "\n"
    msg += "Positional Arguments:\n"
    msg += "  OP2_FILENAME         Path to OP2 file\n"
    msg += "\n"
    msg += "Options:\n"
    msg += "  -b, --binarydebug     Dumps the OP2 as a readable text file\n"
    msg += "  -c, --disablecompare  Doesn't do a validation of the vectorized result\n"
    msg += "  -q, --quiet           Suppresses debug messages [default: False]\n"
    #if not is_release:
    msg += "  -g, --geometry        Reads the OP2 for geometry, which can be written out\n"
    #msg += "  -b, --write_bdf      Writes the bdf to fem.test_op2.bdf (default=False)\n"
    msg += "  -f, --write_f06       Writes the f06 to fem.test_op2.f06\n"
    if not is_release:
        msg += "  -o, --write_op2       Writes the op2 to fem.test_op2.op2\n"
    msg += "  -z, --is_mag_phase    F06 Writer writes Magnitude/Phase instead of\n"
    msg += "                        Real/Imaginary (still stores Real/Imag); [default: False]\n"
    msg += "  -s <sub>, --subcase   Specify one or more subcases to parse; (e.g. 2_5)\n"
    msg += "  -w, --is_sort2        Sets the F06 transient to SORT2\n"
    msg += "  -x <arg>, --exclude   Exclude specific results\n"
    msg += "  -h, --help            Show this help message and exit\n"
    msg += "  -v, --version         Show program's version number and exit\n"

    if len(sys.argv) == 1:
        sys.exit(msg)

    data = docopt(msg, version=ver)
    if '--geometry' not in data:
        data['--geometry'] = False
    if '--write_op2' not in data:
        data['--write_op2'] = False
    if '--write_bdf' not in data:
        data['--write_bdf'] = False
    data['--is_sort2'] = bool(data['--is_sort2'])
    #print("data", data)

    for key, value in sorted(iteritems(data)):
        print("%-12s = %r" % (key.strip('--'), value))

    if os.path.exists('skippedCards.out'):
        os.remove('skippedCards.out')

    t0 = time.time()
    run_op2(data['OP2_FILENAME'],
            make_geom=data['--geometry'],
            write_bdf=data['--write_bdf'],
            write_f06=data['--write_f06'],
            write_op2=data['--write_op2'],
            is_mag_phase=data['--is_mag_phase'],
            isubcases=data['--subcase'],
            exclude=data['--exclude'],
            debug=not data['--quiet'],
            binary_debug=data['--binarydebug'],
            is_sort2=data['--is_sort2'],
            compare=not data['--disablecompare'],
            quiet=data['--quiet'])
    print("dt = %f" % (time.time() - t0))

if __name__ == '__main__':  # op2
    main()
