from __future__ import print_function
from six import iteritems, PY2

import os
import sys
import time
from traceback import print_exc

import pyNastran
from pyNastran.f06.f06 import F06, FatalError
#from pyNastran.op2.test.test_op2 import parseTableNamesFromF06, getFailedFiles


def run_lots_of_files(files, debug=True, save_cases=True, skip_files=None,
                      stop_on_failure=False, nstart=0, nstop=1000000000):
    if skip_files is None:
        skip_files = []
    n = ''
    isubcases = []
    failed_cases = []
    nfailed = 0
    ntotal = 0
    npassed = 0
    t0 = time.time()
    for i, f06file in enumerate(files[nstart:nstop], nstart):  # 149
        base_name = os.path.basename(f06file)
        #if baseName not in skip_files and not base_name.startswith('acms') and i not in nSkip:
        if base_name not in skip_files:
            print("%" * 80)
            print('file=%s\n' % f06file)
            n = '%s ' % i
            sys.stderr.write('%sfile=%s\n' % (n, f06file))
            ntotal += 1
            is_passed = run_f06(f06file, isubcases=isubcases, debug=debug,
                                stop_on_failure=stop_on_failure)  # True/False
            if not is_passed:
                sys.stderr.write('**file=%s\n' % f06file)
                failed_cases.append(f06file)
                nfailed += 1
            else:
                npassed += 1
            #sys.exit('end of test...test_f06.py')

    if save_cases:
        with open('failed_cases.in', 'wb') as failed_cases_file:
            for f06file in failed_cases:
                failed_cases_file.write('%s\n' % (f06file))
    print("dt = %s seconds" % (time.time() - t0))

    sys.exit('-----done with all models %s/%s=%.2f%%  nfailed=%s-----' % (
        npassed, ntotal, 100. * npassed / float(ntotal), ntotal - npassed))


def run_f06(f06_filename, isubcases=None, write_f06=True, is_vector=False,
            debug=False, stop_on_failure=True):
    if isubcases is None:
        isubcases = []
    is_passed = False
    #stopOnFailure = False
    #debug = True
    try:
        f06 = F06(debug=debug)
        #f06.set_subcases(isubcases)  # TODO not supported

        #f06.read_bdf(f06.bdf_filename, includeDir=None, xref=False)
        f06.read_f06(f06_filename, vectorized=is_vector)
        #tableNamesF06 = parseTableNamesFromF06(f06.f06FileName)
        #tableNamesF06 = f06.getTableNamesFromF06()
        #assert write_f06 == True, write_f06
        if write_f06:
            model = os.path.splitext(f06_filename)[0]
            #print(f06.get_f06_stats())
            f06.write_f06(model + '.test_f06.f06')

        #print("subcases = ",f06.subcases)

        #assert tableNamesF06==tableNamesF06,'tableNamesF06=%s tableNamesF06=%s' %(tableNamesF06,tableNamesF06)
        #f06.case_control_deck.sol = f06.sol
        #print(f06.case_control_deck.getF06Data())
        #print(f06.print_results())
        #print(f06.case_control_deck.getF06Data())
        is_passed = True
    except KeyboardInterrupt:
        sys.stdout.flush()
        print_exc(file=sys.stdout)
        sys.stderr.write('**file=%r\n' % f06_filename)
        sys.exit('keyboard stop...')
    #except AddNewElementError:
    #    raise
    #except IOError: # missing file
        #pass
    #except AssertionError:
    #    is_passed = True
    except FatalError:  # remove this later...
        is_passed = True
    #except InvalidFormatCodeError:
    #    is_passed = True
    #except RuntimeError: #InvalidAnalysisCode
    #    is_passed = True
    #except SyntaxError: #Invalid Markers
    #    is_passed = True
    except SystemExit:
        #print_exc(file=sys.stdout)
        #sys.exit('stopping on sys.exit')
        raise
    #except NameError:  # variable isnt defined
    #    if stop_on_failure:
    #        raise
    #    else:
    #        is_passed = True
    #except AttributeError:  # missing function
    #    if stop_on_failure:
    #        raise
    #    else:
    #        is_passed = True
    #except KeyError:
    #    raise
    #except TypeError:  # numpy error
    #    is_passed = True
    #except IndexError: # bad bdf
    #    is_passed = True
    #except IOError:  # missing bdf file
        #is_passed = False
        #raise
    #except SyntaxError: #Invalid Subcase
    #    is_passed = True
    #except SyntaxError: # Param Parse:
    #    is_passed = True
    #except NotImplementedError:
        #is_passed = True
    #except InvalidFieldError: # bad bdf field
    #    is_passed = True
    except:
        #print(e)
        print_exc(file=sys.stdout)
        if stop_on_failure:
            raise
        else:
            is_passed = False
    #print("is_passed = %s" % is_passed)
    return is_passed


def main():
    from docopt import docopt

    msg = 'Tests to see if an F06 will work with pyNastran.\n'
    msg += 'Usage:\n'
    msg += '  f06.py [-f] [-q] [-t] F06_FILENAME' # [-p]
    msg += '  f06.py -h | --help\n'
    msg += '  f06.py -v | --version\n'
    msg += '\n'
    msg += 'Positional Arguments:\n'
    msg += '  F06_FILENAME         path to F06 file\n'
    msg += '\n'
    msg += 'Options:\n'
    msg += '  -q, --quiet      prints debug messages (default=False)\n'
    msg += '  -f, --write_f06  writes the f06 to fem.f06.out (default=True)\n'
    msg += "  -t, --vector     vectorizes the results (default=False)\n"
    msg += '  -h, --help       show this help message and exit\n'
    msg += "  -v, --version    show program's version number and exit\n"

    # disabled b/c the F06 doesn't support complex well
    #msg += '  -z, --is_mag_phase      F06 Writer writes Magnitude/Phase instead of\n'
    #msg += '                          Real/Imaginary (still stores Real/Imag)\n'

    if len(sys.argv) == 1:
        sys.exit(msg)

    ver = str(pyNastran.__version__)
    data = docopt(msg, version=ver)

    for key, value in sorted(iteritems(data)):
        print("%-12s = %r" % (key.strip('--'), value))

    if os.path.exists('skipped_cards.out'):
        os.remove('skipped_cards.out')
    run_f06(
        data['F06_FILENAME'],
        write_f06=data['--write_f06'],
        debug=not(data['--quiet']),
        is_vector=data['--vector'],
        stop_on_failure=True
    )

if __name__ == '__main__':  # pragma: no cover
    main()
