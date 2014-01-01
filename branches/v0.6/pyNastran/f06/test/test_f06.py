import os
import sys
import time
from traceback import print_exc

import pyNastran
from pyNastran.f06.f06 import F06
#from pyNastran.op2.test.test_op2 import parseTableNamesFromF06, getFailedFiles


def run_lots_of_files(files, debug=True, saveCases=True, skipFiles=[],
                      stopOnFailure=False, nStart=0, nStop=1000000000):
    n = ''
    iSubcases = []
    failedCases = []
    nFailed = 0
    nTotal = 0
    nPassed = 0
    t0 = time.time()
    for i, f06file in enumerate(files[nStart:nStop], nStart):  # 149
        baseName = os.path.basename(f06file)
        #if baseName not in skipFiles and not baseName.startswith('acms') and i not in nSkip:
        if baseName not in skipFiles:
            print("%" * 80)
            print('file=%s\n' % (f06file))
            n = '%s ' % (i)
            sys.stderr.write('%sfile=%s\n' % (n, f06file))
            nTotal += 1
            isPassed = run_f06(f06file, iSubcases=iSubcases, debug=debug,
                               stopOnFailure=stopOnFailure)  # True/False
            if not isPassed:
                sys.stderr.write('**file=%s\n' % (f06file))
                failedCases.append(f06file)
                nFailed += 1
            else:
                nPassed += 1
            #sys.exit('end of test...test_f06.py')

    if saveCases:
        f = open('failedCases.in', 'wb')
        for f06file in failedCases:
            f.write('%s\n' % (f06file))
        f.close()
    print("dt = %s seconds" % (time.time() - t0))

    #f06 = F06('test_tet10_subcase_1.f06')
    #f06.readF06()

    sys.exit('-----done with all models %s/%s=%.2f%%  nFailed=%s-----' % (nPassed, nTotal, 100. * nPassed / float(nTotal), nTotal - nPassed))


def run_f06(f06file, iSubcases=[], write_f06=True, print_f06=False, debug=False,
            stopOnFailure=True):
    isPassed = False
    stopOnFailure = False
    #debug = True
    try:
        f06 = F06(f06file, debug=debug)
        #f06.set_subcases(iSubcases)  # TODO not supported

        #f06.readBDF(f06.bdf_filename,includeDir=None,xref=False)
        f06.read_f06()
        #tableNamesF06 = parseTableNamesFromF06(f06.f06FileName)
        #tableNamesF06 = f06.getTableNamesFromF06()
        assert write_f06 == True, write_f06
        if write_f06:
            (model, ext) = os.path.splitext(f06file)
            f06.write_f06(model + '.test_f06.f06')

        if print_f06:
            print(f06.print_results())
        #print "subcases = ",f06.subcases

        #assert tableNamesF06==tableNamesF06,'tableNamesF06=%s tableNamesF06=%s' %(tableNamesF06,tableNamesF06)
        #f06.caseControlDeck.sol = f06.sol
        #print f06.caseControlDeck.getF06Data()
        #print f06.print_results()
        #print f06.caseControlDeck.getF06Data()
        isPassed = True
    except KeyboardInterrupt:
        sys.stdout.flush()
        print_exc(file=sys.stdout)
        sys.stderr.write('**file=%r\n' % f06file)
        sys.exit('keyboard stop...')
    #except AddNewElementError:
    #    raise
    #except IOError: # missing file
        #pass
    #except AssertionError:
    #    isPassed = True

    #except InvalidFormatCodeError:
    #    isPassed = True
    #except RuntimeError: #InvalidAnalysisCode
    #    isPassed = True
    #except SyntaxError: #Invalid Markers
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
    #except AttributeError:  # missing function
    #    if stopOnFailure:
    #        raise
    #    else:
    #        isPassed = True
    #except KeyError:
    #    raise
    #except TypeError:  # numpy error
    #    isPassed = True
    #except IndexError: # bad bdf
    #    isPassed = True
    except IOError:  # missing bdf file
        isPassed = False
        raise
    #except SyntaxError: #Invalid Subcase
    #    isPassed = True
    #except SyntaxError: # Param Parse:
    #    isPassed = True
    #except NotImplementedError:
        #isPassed = True
    #except InvalidFieldError: # bad bdf field
    #    isPassed = True
    except:
        #print e
        print_exc(file=sys.stdout)
        if stopOnFailure:
            raise
        else:
            isPassed = False
    print "isPassed =", isPassed
    return isPassed


def main():
    from docopt import docopt

    msg  = 'Tests to see if an F06 will work with pyNastran.\n'
    msg += 'Usage:\n'
    msg += '  f06.py [-f] [-p] [-q] F06_FILENAME'
    msg += '  f06.py -h | --help\n'
    msg += '  f06.py -v | --version\n'
    msg += '\n'
    msg += 'Positional Arguments:\n'
    msg += '  F06_FILENAME         path to F06 file\n'
    msg += '\n'
    msg += 'Options:\n'
    msg += '  -q, --quiet      prints debug messages (default=False)\n'
    msg += '  -f, --write_f06  writes the f06 to fem.f06.out (default=True)\n'
    msg += '  -p, --print_f06  prints objects to screen which can require lots of\n'
    msg += '                    memory (default=True)\n'
    msg += '  -h, --help       show this help message and exit\n'
    msg += "  -v, --version    show program's version number and exit\n"
   #msg += '  -z, --is_mag_phase      F06 Writer writes Magnitude/Phase instead of\n'
   #msg += '                          Real/Imaginary (still stores Real/Imag)\n'

    if len(sys.argv) == 1:
        sys.exit(msg)

    ver = str(pyNastran.__version__)
    data = docopt(msg, version=ver)
    
    for key, value in sorted(data.iteritems()):
        print("%-12s = %r" % (key.strip('--'), value))

    if os.path.exists('skippedCards.out'):
        os.remove('skippedCards.out')
    run_f06(data['F06_FILENAME'],
            write_f06 = data['--write_f06'],
            print_f06 = data['--print_f06'],
            debug     = not(data['--quiet'])
    )

if __name__ == '__main__':  # f06
    main()
