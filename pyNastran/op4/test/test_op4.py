"""defines the command line argument ``test_op4``"""
import os
import sys
import time
from traceback import print_exc

import pyNastran
from pyNastran.op4.op4 import read_op4


def run_lots_of_files(files, write_op4=True,
                      debug=True, save_cases=True, skip_files=None,
                      stop_on_failure=False, nstart=0, nstop=1000000000):
    """runs lots of op4 files"""
    if skip_files is None:
        skip_files = []
    n = ''
    failed_cases = []
    nfailed = 0
    ntotal = 0
    npassed = 0
    t0 = time.time()
    for (i, op4file) in enumerate(files[nstart:nstop], nstart):  # 149
        base_name = os.path.basename(op4file)
        #if baseName not in skipFiles and not base_name.startswith('acms') and i not in nSkip:
        if base_name not in skip_files and '#' not in op4file:
            print("%"*80)
            print('file=%s\n' % op4file)
            n = '%s ' % i
            sys.stderr.write('%sfile=%s\n' %(n, op4file))
            ntotal += 1
            is_passed = run_op4(op4file,
                                debug=debug,
                                stop_on_failure=stop_on_failure) # True/False
            if not is_passed:
                sys.stderr.write('**file=%s\n' % op4file)
                failed_cases.append(op4file)
                nfailed += 1
            else:
                npassed += 1
            #sys.exit('end of test...test_op4.py')

    if save_cases:
        with open('failed_cases.in', 'wb') as failed_file:
            for op4file in failed_cases:
                failed_file.write('%s\n' % op4file)

    seconds = time.time()-t0
    minutes = seconds/60.
    print("dt = %s seconds = %s minutes" % (seconds, minutes))

    msg = '-----done with all models %s/%s=%.2f%%  nFailed=%s-----' %(
        npassed, ntotal, 100.*npassed/float(ntotal), ntotal-npassed)
    print(msg)
    sys.exit(msg)


def run_op4(op4_filename, write_op4=True, debug=True,
            stop_on_failure=False):
    """run an op4"""
    #print('***debug=%s' % debug)
    assert '.op4' in op4_filename.lower(), 'op4_filename=%s is not an OP4' % op4_filename
    is_passed = False
    stop_on_failure = True
    delete_op4 = True

    #debug = True
    try:
        matrices = read_op4(op4_filename, debug=debug)
        keys = list(matrices.keys())
        keys.sort()
        print('matrices =', keys)

        #if 0:
            #matrices2 = op4.read_op4(op4_filename)

            #print(matrices)
            #print('matrices =', matrices.keys())

            #assert list(sorted(matrices.keys())) == list(sorted(matrices2.keys()))
            #for key, (form, matrix) in sorted(matrices.items()):
                #form2, matrix2 = matrices2[key]
                #assert form == form2
                #delta = matrix - matrix2
                #assert np.array_equal(matrix, matrix2), 'delta=\n%s' % delta

        if write_op4:
            model = os.path.splitext(op4_filename)[0]
            model.write_op4(model+'.test_op4_ascii.op4', matrices, is_binary=False)
            model.write_op4(model+'.test_op4_binary.op4', matrices, is_binary=True)
            if delete_op4:
                try:
                    os.remove(model+'.test_op4_ascii.op4')
                    os.remove(model+'.test_op4_binary.op4')
                except:
                    pass

        is_passed = True
    except KeyboardInterrupt:
        sys.stdout.flush()
        print_exc(file=sys.stdout)
        sys.stderr.write('**file=%s\n' % op4_filename)
        sys.exit('keyboard stop...')
    #except RuntimeError: # the op2 is bad, not my fault
    #    is_passed = True
    #    if stop_on_failure:
    #        raise
    #    else:
    #        is_passed = True

    except IOError: # missing file
        if stop_on_failure:
            raise
    #except AssertionError:
    #    is_passed = True
    #except RuntimeError:
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
    #except IndexError:
    #    is_passed = True
    except SyntaxError: #Param Parse
        if stop_on_failure:
            raise
        is_passed = True
    except:
        #print(e)
        if stop_on_failure:
            raise
        else:
            print_exc(file=sys.stdout)
            is_passed = False
    return is_passed


def main():
    """defines the command line argument ``test_op4``"""
    from docopt import docopt
    ver = str(pyNastran.__version__)

    msg = "Usage:\n"

    # all
    # release
    # current
    msg += "test_op4 [-o] [-d] OP4_FILENAME\n"
    msg += "  test_op4 -h | --help\n"
    msg += "  test_op4 -v | --version\n"
    msg += "\n"
    msg += "Tests to see if an OP4 will work with pyNastran %s.\n" % ver
    msg += "\n"
    msg += "Positional Arguments:\n"
    msg += "  OP4_FILENAME         Path to OP4 file\n"
    msg += "\n"
    msg += "Options:\n"
    msg += "  -d, --debug          Developer Debug (default=False)\n"
    msg += "  -o, --write_op4      Writes the op2 to fem.test_op4.op4 (default=True)\n"
    msg += "  -h, --help           Show this help message and exit\n"
    msg += "  -v, --version        Show program's version number and exit\n"

    if len(sys.argv) == 1:
        sys.exit(msg)

    data = docopt(msg, version=ver)
    #print("data", data)

    for key, value in sorted(data.items()):
        print("%-12s = %r" % (key.strip('--'), value))

    time0 = time.time()
    run_op4(
        data['OP4_FILENAME'],
        write_op4=data['--write_op4'],
        debug=data['--debug'],
    )
    print("dt = %f" % (time.time() - time0))


if __name__ == '__main__':   # pragma: no cover
    main()
