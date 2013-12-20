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


def run_f06(f06file, iSubcases=[], write_f06=True, printF06=False, debug=False,
            stopOnFailure=True):
    isPassed = False
    stopOnFailure = False
    #debug = True
    try:
        f06 = F06(f06file, debug=debug)
        #f06.set_subcases(iSubcases)  # TODO not supported

        #f06.readBDF(f06.bdf_filename,includeDir=None,xref=False)
        f06.readF06()
        #tableNamesF06 = parseTableNamesFromF06(f06.f06FileName)
        #tableNamesF06 = f06.getTableNamesFromF06()
        if write_f06:
            (model, ext) = os.path.splitext(f06file)
            f06.write_f06(model + '.f06.out')

        if printF06:
            f06.print_results()
        #print "subcases = ",f06.subcases

        #assert tableNamesF06==tableNamesF06,'tableNamesF06=%s tableNamesF06=%s' %(tableNamesF06,tableNamesF06)
        pass
        #f06.caseControlDeck.sol = f06.sol
        #print f06.caseControlDeck.getF06Data()
        #print f06.print_results()
        #print f06.caseControlDeck.getF06Data()
        isPassed = True
    except KeyboardInterrupt:
        sys.stdout.flush()
        print_exc(file=sys.stdout)
        sys.stderr.write('**file=%s\n' % (f06file))
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
    except NotImplementedError:
        isPassed = True
    #except InvalidFieldError: # bad bdf field
    #    isPassed = True
    except:
        #print e
        print_exc(file=sys.stdout)
        if stopOnFailure:
            raise
        else:
            isPassed = False
    return isPassed


def run_arg_parse():
    import argparse

    ver = str(pyNastran.__version__)
    parser = argparse.ArgumentParser(description='Tests to see if an F06 will work with pyNastran.', add_help=True)  # ,version=ver)
    parser.add_argument('f06FileName', metavar='f06FileName', type=str, nargs=1,
                        help='path to F06 file')

    group = parser.add_mutually_exclusive_group()
    group.add_argument('-q', '--quiet', dest='quiet', action='store_true',
                       help='Prints   debug messages (default=True)')
    parser.add_argument('-f', '--write_f06', dest='write_f06',
                        action='store_true', help='Writes the f06 to fem.f06.out')
    parser.add_argument('-p', '--printF06', dest='printF06',
                        action='store_true', help='Prints the F06 to the screen, slow & uses a lot of memory for large files')
    parser.add_argument('-v', '--version', action='version', version=ver)

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit()
    args = parser.parse_args()
    print("f06FileName = %s" % (args.f06FileName[0]))
    print("debug       = %s" % (not(args.quiet)))

    debug = not(args.quiet)
    write_f06 = args.write_f06
    printF06 = args.printF06
    f06FileName = args.f06FileName[0]

    return (f06FileName, write_f06, printF06, debug)


def main():
    (f06FileName, write_f06, printF06, debug) = run_arg_parse()
    if os.path.exists('skippedCards.out'):
        os.remove('skippedCards.out')
    run_f06(f06FileName, write_f06=write_f06, printF06=printF06, debug=debug)

if __name__ == '__main__':  # f06
    main()
