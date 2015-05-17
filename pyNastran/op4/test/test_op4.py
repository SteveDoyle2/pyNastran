from __future__ import print_function
from six import iteritems
import os
import sys
import time
from traceback import print_exc

import pyNastran
from pyNastran.op4.op4 import OP4


def run_lots_of_files(files, write_op4=True,
                   debug=True, saveCases=True, skipFiles=None,
                   stopOnFailure=False, nStart=0, nStop=1000000000):
    if skipFiles is None:
        skipFiles = []
    n = ''
    iSubcases = []
    failedCases = []
    nFailed = 0
    nTotal = 0
    nPassed = 0
    t0 = time.time()
    for (i, op4file) in enumerate(files[nStart:nStop], nStart):  # 149
        baseName = os.path.basename(op4file)
        #if baseName not in skipFiles and not baseName.startswith('acms') and i not in nSkip:
        if baseName not in skipFiles and '#' not in op4file:
            print("%"*80)
            print('file=%s\n' % op4file)
            n = '%s ' % i
            sys.stderr.write('%sfile=%s\n' %(n, op4file))
            nTotal += 1
            isPassed = run_op4(op4file,
                               debug=debug,
                               stopOnFailure=stopOnFailure) # True/False
            if not isPassed:
                sys.stderr.write('**file=%s\n' % op4file)
                failedCases.append(op4file)
                nFailed += 1
            else:
                nPassed += 1
            #sys.exit('end of test...test_op4.py')

    if saveCases:
        f = open('failedCases.in', 'wb')
        for op4file in failedCases:
            f.write('%s\n' % op4file)
        f.close()

    seconds = time.time()-t0
    minutes = seconds/60.
    print("dt = %s seconds = %s minutes" % (seconds, minutes))

    msg = '-----done with all models %s/%s=%.2f%%  nFailed=%s-----' %(nPassed,nTotal,100.*nPassed/float(nTotal),nTotal-nPassed)
    print(msg)
    sys.exit(msg)


def run_op4(op4_filename, write_op4=True, debug=True):
    assert '.op4' in op4_filename.lower(), 'op4_filename=%s is not an OP4' % op4_filename
    isPassed = False
    stopOnFailure = True
    delete_op4 = True

    #debug = True
    try:
        op4 = OP4()
        matrices = op4.read_op4(op4_filename)
        print(matrices)
        print('matrices =', matrices.keys())

        if write_op4:
            model, ext = os.path.splitext(op4_filename)
            op4.write_op4(model+'.test_op4_ascii.op4', matrices, is_binary=False)
            op4.write_op4(model+'.test_op4_binary.op4', matrices, is_binary=True)
            if delete_op4:
                try:
                    os.remove(model+'.test_op4_ascii.op4')
                    os.remove(model+'.test_op4_binary.op4')
                except:
                    pass

        del op4
        isPassed = True
    except KeyboardInterrupt:
        sys.stdout.flush()
        print_exc(file=sys.stdout)
        sys.stderr.write('**file=%s\n' % op4_filename)
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
    #except AssertionError:
    #    isPassed = True
    #except RuntimeError:
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
    #except IndexError:
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
    # release
    # current
    msg += "test_op4 [-q] [-o] OP4_FILENAME\n"
    msg += "  test_op4 -h | --help\n"
    msg += "  test_op4 -v | --version\n"
    msg += "\n"
    msg += "Tests to see if an OP4 will work with pyNastran %s.\n" % ver
    msg += "\n"
    msg += "Positional Arguments:\n"
    msg += "  OP4_FILENAME         Path to OP4 file\n"
    msg += "\n"
    msg += "Options:\n"
    msg += "  -q, --quiet          Suppresses debug messages (default=False)\n"
    msg += "  -o, --write_op4      Writes the op2 to fem.test_op4.op4 (default=True)\n"
    msg += "  -h, --help           Show this help message and exit\n"
    msg += "  -v, --version        Show program's version number and exit\n"

    if len(sys.argv) == 1:
        sys.exit(msg)

    data = docopt(msg, version=ver)
    #print("data", data)

    for key, value in sorted(iteritems(data)):
        print("%-12s = %r" % (key.strip('--'), value))

    import time
    t0 = time.time()
    run_op4(data['OP4_FILENAME'],
            write_op4 = data['--write_op4'],
            debug     = not(data['--quiet']),
    )
    print("dt = %f" % (time.time() - t0))


if __name__=='__main__':  # op4
    main()
