from six import string_types
import os
import sys
import time
from traceback import print_exc

import pyNastran
from pyNastran.f06.errors import FatalError
from pyNastran.op2.op2 import OP2
from pyNastran.op2.op2_vectorized import OP2_Vectorized as OP2V

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


def parse_table_names_from_F06(f06Name):
    """gets the op2 names from the f06"""
    infile = open(f06Name,'r')
    marker = 'NAME OF DATA BLOCK WRITTEN ON FORTRAN UNIT IS'
    names = []
    for line in infile:
        if marker in line:
            word = line.replace(marker,'').strip().strip('.')
            names.append(word)

    infile.close()
    return names


def get_failed_files(filename):
    infile = open(filename, 'r')
    lines = infile.readlines()
    infile.close()

    files = []
    for line in lines:
        files.append(line.strip())
    return files


def run_lots_of_files(files ,make_geom=True, write_bdf=False, write_f06=True,
                   delete_f06=True, write_op2=True, is_vector=False,
                   debug=True, saveCases=True, skipFiles=[],
                   stopOnFailure=False, nStart=0, nStop=1000000000):
    n = ''
    assert make_geom in [True, False]
    assert write_bdf in [True, False]
    assert write_f06 in [True, False]
    iSubcases = []
    failedCases = []
    nFailed = 0
    nTotal  = 0
    nPassed = 0
    t0 = time.time()
    for (i, op2file) in enumerate(files[nStart:nStop], nStart):  # 149
        baseName = os.path.basename(op2file)
        #if baseName not in skipFiles and not baseName.startswith('acms') and i not in nSkip:
        if baseName not in skipFiles and '#' not in op2file:
            print("%"*80)
            print('file=%s\n' % op2file)
            n = '%s ' %(i)
            sys.stderr.write('%sfile=%s\n' %(n, op2file))
            nTotal += 1
            isPassed = run_op2(op2file, make_geom=make_geom, write_bdf=write_bdf,
                               write_f06=write_f06, write_op2=write_op2,
                               is_mag_phase=False,
                               is_vector=is_vector,
                               delete_f06=delete_f06,
                               iSubcases=iSubcases, debug=debug,
                               stopOnFailure=stopOnFailure) # True/False
            if not isPassed:
                sys.stderr.write('**file=%s\n' % op2file)
                failedCases.append(op2file)
                nFailed +=1
            else:
                nPassed +=1
            #sys.exit('end of test...test_op2.py')

    if saveCases:
        f = open('failedCases.in','wb')
        for op2file in failedCases:
            f.write('%s\n' % op2file)
        f.close()

    seconds = time.time()-t0
    minutes = seconds/60.
    print("dt = %s seconds = %s minutes" % (seconds, minutes))

    #op2 = OP2('test_tet10_subcase_1.op2')
    #op2.read_op2()

    msg = '-----done with all models %s/%s=%.2f%%  nFailed=%s-----' %(nPassed,nTotal,100.*nPassed/float(nTotal),nTotal-nPassed)
    print(msg)
    sys.exit(msg)


def run_op2(op2FileName, make_geom=False, write_bdf=False,
            write_f06=True, write_op2=True, is_mag_phase=False,
            is_vector=False, delete_f06=False,
            iSubcases=[], exclude=[], debug=False, stopOnFailure=True):
    assert '.op2' in op2FileName.lower(), 'op2FileName=%s is not an OP2' % op2FileName
    isPassed = False
    if isinstance(iSubcases, basestring):
        if '_' in iSubcases:
            iSubcases = [int(i) for i in iSubcases.split('_')]
        else:
            iSubcases = [int(iSubcases)]
    print('iSubcases =', iSubcases)

    #debug = True
    try:
        if is_vector:
            op2 = OP2V(make_geom=make_geom, debug=debug)
        else:
            op2 = OP2(make_geom=make_geom, debug=debug)

        op2.set_subcases(iSubcases)
        op2.remove_results(exclude)

        if is_memory:
            if is_linux: # linux
                kb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
            else: # windows
                kb = get_memory_usage() / 1024
            mb = kb / 1024.
            #mbs.append(mb)
            print("Memory usage start: %s (KB); %.2f (MB)" % (kb, mb))

        #op2.read_bdf(op2.bdfFileName,includeDir=None,xref=False)
        #op2.write_bdf_as_patran()
        op2.read_op2(op2FileName)
        print("---stats for %s---" % op2FileName)
        #op2.get_op2_stats()
        print(op2.get_op2_stats())
        if write_bdf and 0:
            op2.write_bdf('fem.bdf.out', interspersed=True)
        #tableNamesF06 = parse_table_names_from_F06(op2.f06FileName)
        #tableNamesOP2 = op2.getTableNamesFromOP2()

        if is_memory:
            if is_linux: # linux
                kb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
            else: # windows
                kb = get_memory_usage() / 1024
            mb = kb / 1024.
            #mbs.append(mb)
            print("Memory usage     end: %s (KB); %.2f (MB)" % (kb, mb))

        if write_f06:
            (model, ext) = os.path.splitext(op2FileName)
            op2.write_f06(model+'.test_op2.f06', is_mag_phase=is_mag_phase)
            if delete_f06:
                try:
                    os.remove(model+'.test_op2.f06')
                except:
                    pass

        if write_op2:
            (model, ext) = os.path.splitext(op2FileName)
            op2.write_op2(model+'.test_op2.op2', is_mag_phase=is_mag_phase)
            if delete_f06:
                try:
                    os.remove(model+'.test_op2.op2')
                except:
                    pass

        del op2
        if is_memory:
            if is_linux: # linux
                kb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
            else: # windows
                kb = get_memory_usage() / 1024
            mb = kb / 1024.
            #mbs.append(mb)
            print("Memory usage cleanup: %s (KB); %.2f (MB)" % (kb, mb))


        #print("subcases = ",op2.subcases)

        #assert tableNamesF06==tableNamesOP2,'tableNamesF06=%s tableNamesOP2=%s' %(tableNamesF06,tableNamesOP2)
        #op2.caseControlDeck.sol = op2.sol
        #print(op2.caseControlDeck.get_op2_data())
        #print(op2.caseControlDeck.get_op2_data())
        isPassed = True
    except KeyboardInterrupt:
        sys.stdout.flush()
        print_exc(file=sys.stdout)
        sys.stderr.write('**file=%s\n' % op2FileName)
        sys.exit('keyboard stop...')
    #except RuntimeError: # the op2 is bad, not my fault
    #    isPassed = True
    #    if stopOnFailure:
    #        raise
    #    else:
    #        isPassed = True

    except IOError: # missing file
        if stopOnFailure:
            raise
    except FatalError:
        if stopOnFailure:
            raise
        isPassed = True
    #except AssertionError:
    #    isPassed = True
    #except RuntimeError: #invalid analysis code
    #    isPassed = True
    except SystemExit:
        #print_exc(file=sys.stdout)
        #sys.exit('stopping on sys.exit')
        raise
    #except NameError:  # variable isnt defined
    #    if stopOnFailure:
    #        raise
    #    else:
    #        isPassed = True
    #except IndexError: # bad bdf
    #    isPassed = True
    except SyntaxError: #Param Parse
        if stopOnFailure:
            raise
        isPassed = True
    except:
        #print e
        if stopOnFailure:
            raise
        else:
            print_exc(file=sys.stdout)
            isPassed = False

    return isPassed


def main():
    from docopt import docopt
    ver = str(pyNastran.__version__)

    msg  = "Usage:\n"

    # all
    #msg += "test_op2 [-q] [-g] [-w] [-f] [-o] [-z] [-t] [-s <sub>] OP2_FILENAME\n"

    # release
    #msg += "test_op2 [-q] [-f] [-z] [-t] [-s <sub>] OP2_FILENAME\n"

    # current
    msg += "test_op2 [-q] [-f] [-o] [-z] [-t] [-s <sub>] [-x <arg>]... OP2_FILENAME\n"
    msg += "  test_op2 -h | --help\n"
    msg += "  test_op2 -v | --version\n"
    msg += "\n"
    msg += "Tests to see if an OP2 will work with pyNastran %s.\n" % ver
    msg += "\n"
    msg += "Positional Arguments:\n"
    msg += "  OP2_FILENAME         Path to OP2 file\n"
    msg += "\n"
    msg += "Options:\n"
    msg += "  -q, --quiet          Suppresses debug messages (default=False)\n"
    #msg += "  -g, --geometry       Reads the OP2 for geometry, which can be written out (default=False)\n"
    #msg += "  -w, --write_bdf      Writes the bdf to fem.test_op2.bdf (default=False)\n"
    msg += "  -f, --write_f06      Writes the f06 to fem.test_op2.f06 (default=True)\n"
    msg += "  -o, --write_op2      Writes the op2 to fem.test_op2.op2 (default=False)\n"
    msg += "  -z, --is_mag_phase   F06 Writer writes Magnitude/Phase instead of\n"
    msg += "                       Real/Imaginary (still stores Real/Imag); (default=False)\n"
    msg += "  -s <sub>, --subcase  Specify one or more subcases to parse; (e.g. 2_5)\n"
    msg += "  -t, --vector         Vectorizes the results (default=False)\n"
    msg += "  -x <arg>, --exclude  Exclude specific results\n"
    msg += "  -h, --help           Show this help message and exit\n"
    msg += "  -v, --version        Show program's version number and exit\n"

    if len(sys.argv) == 1:
        sys.exit(msg)

    data = docopt(msg, version=ver)
    #print("data", data)

    for key, value in sorted(data.iteritems()):
        print("%-12s = %r" % (key.strip('--'), value))

    if os.path.exists('skippedCards.out'):
        os.remove('skippedCards.out')

    import time
    t0 = time.time()
    run_op2(data['OP2_FILENAME'],
            #make_geom     = data['--geometry'],
            #write_bdf     = data['--write_bdf'],
            write_f06     = data['--write_f06'],
            write_op2     = data['--write_op2'],
            is_mag_phase  = data['--is_mag_phase'],
            is_vector     = data['--vector'],
            iSubcases     = data['--subcase'],
            exclude       = data['--exclude'],
            debug         = not(data['--quiet']))
    print("dt = %f" % (time.time() - t0))

if __name__=='__main__':  # op2
    main()
