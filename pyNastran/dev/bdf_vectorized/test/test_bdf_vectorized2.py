# pylint: disable=W0612,C0103
"""
``test_bdf`` runs multiple checks on a BDF in order to make sure that:
  - no data is lost on IO
  - card field types are correct (e.g. node_ids are integers)
  - various card methods (e.g. Area) work correctly

As such, ``test_bdf`` is very useful for debugging models.

"""
import os
import sys
import numpy
import warnings
warnings.simplefilter('always')
numpy.seterr(all='raise')
#import traceback
#import resource

from pyNastran.utils import check_path
from pyNastran.dev.bdf_vectorized.bdf import BDF, read_bdf #, NastranMatrix

import pyNastran.dev.bdf_vectorized.test
test_path = pyNastran.dev.bdf_vectorized.test.__path__[0]
#print("test_path = ",test_path)


def run_lots_of_files(filenames, folder='', debug=False, xref=True, check=True,
                      punch=False, cid=None):
    """
    Runs multiple BDFs

    Parameters
    ----------
    folder : str
        the folder where the bdf_filename is
    filenames : list[str]
        the bdf files to analyze
    debug : bool, optional
        run with debug logging (default=False)
    xref : bool / str / list[bool/str], optional
        True : cross reference the model
        False  : don't cross reference the model
        'safe' : do safe cross referencing
    check : bool / list[bool], optional
        validate cards for things like mass, area, etc. (default=True)
    punch : bool / list[bool], optional
        this is a PUNCH file (no executive/case control decks; default=False)
    cid : int / None, optional
        convert the model grids to an alternate coordinate system (default=None; no conversion)
    size : int / list[int], optional
        The field width of the model (8/16)
    is_double : bool / list[bool], optional
        Is this a double precision model?
            True : size = 16
            False : size = {8, 16}
    nastran : str, optional
        the path to nastran (default=''; no analysis)
    post : int / list[int], optional
        the PARAM,POST,value to run
    sum_load : bool; default=True
        should the loads be summed
    dev : bool; default=True
        True : crashes if an Exception occurs
        False : doesn't crash; useful for running many tests
    crash_cards : list[str, str, ...]
        list of cards that are invalid and automatically crash the run

    Examples
    --------
    All control lists must be the same length.

    You can run xref=True and xref=False with::

        run_lots_of_files(filenames, xref=[True, False]) # valid
    """
    filenames = list(set(filenames))
    filenames.sort()

    if size is None:
        sizes = [8]
    elif isinstance(size, integer_types):
        sizes = [size]
    else:
        sizes = size

    if is_double is None:
        is_doubles = [8]
    elif isinstance(is_double, bool):
        is_doubles = [is_double]
    else:
        is_doubles = is_double

    #debug = True
    filenames2 = []
    diff_cards = []
    for filename in filenames:
        if(filename.endswith('.bdf') or filename.endswith('.dat') or
           filename.endswith('.nas') or filename.endswith('.nas')):
            filenames2.append(filename)

    failed_files = []
    n = 1
    for filename in filenames2:
        abs_filename = os.path.abspath(os.path.join(folder, filename))
        if folder != '':
            print("filename = %s" % abs_filename)
        is_passed = False
        #try:
        (fem1, fem2, diff_cards) = run_bdf(folder, filename, debug=debug,
                                           xref=xref, check=check, punch=punch,
                                           cid=cid, isFolder=True, dynamic_vars={},
                                           run_extract_bodies=False)
        del fem1
        del fem2
        diff_cards += diff_cards
        is_passed = True
        #except KeyboardInterrupt:
            #sys.exit('KeyboardInterrupt...sys.exit()')
        #except IOError:
            #pass
        #except RuntimeError:  # only temporarily uncomment this when running lots of tests
            #pass
        #except AttributeError:  # only temporarily uncomment this when running lots of tests
            #pass
        #except SyntaxError:  # only temporarily uncomment this when running lots of tests
            #pass
        #except SystemExit:
            #sys.exit('sys.exit...')
        #except Exception:
            #traceback.print_exc(file=sys.stdout)
            ##raise
        print('-' * 80)

        if is_passed:
            sys.stderr.write('%i %s' % (n, abs_filename))
            n += 1
        else:
            sys.stderr.write('*' + abs_filename)
            failed_files.append(abs_filename)
        sys.stderr.write('\n')

    print('*' * 80)
    try:
        print("diff_cards1 = %s" % list(set(diff_cards)))
    except TypeError:
        print("diff_cards2 = %s" % diff_cards)
    return failed_files


def run_bdf(folder, bdf_filename, debug=False, xref=True, check=True, punch=False,
            cid=None, mesh_form='combined', is_folder=False, print_stats=False,
            encoding=None, sum_load=False, size=8, is_double=False,
            quiet=False,
            reject=False, dynamic_vars=None, run_extract_bodies=True):
    """
    Runs a single BDF

    Parameters
    ----------
    folder : str
        the folder where the bdf_filename is
    bdf_filename : str
        the bdf file to analyze
    debug : bool, optional
        run with debug logging (default=False)
    xref : bool / str, optional
        True : cross reference the model
        False  : don't cross reference the model
        'safe' : do safe cross referencing
    check : bool, optional
        validate cards for things like mass, area, etc.
    punch : bool, optional
        this is a PUNCH file (no executive/case control decks)
    cid : int / None, optional
        convert the model grids to an alternate coordinate system (default=None; no conversion)
    mesh_form : str, optional, {'combined', 'separate'}
        'combined' : interspersed=True
        'separate' : interspersed=False
    is_folder : bool, optional
        attach the test path and the folder to the bdf_filename
    print_stats : bool, optional
        get a nicely formatted message of all the cards in the model
    sum_load : bool, optional
        Sum the static loads (doesn't work for frequency-based loads)
    size : int, optional, {8, 16}
        The field width of the model
    is_double : bool, optional
        Is this a double precision model?
            True : size = 16
            False : size = {8, 16}
    reject : bool, optional
        True : all the cards are rejected
        False : the model is read
    nastran : str, optional
        the path to nastran (default=''; no analysis)
    post : int, optional
        the PARAM,POST,value to run
    dynamic vars : dict[str]=int / float / str / None
        support OpenMDAO syntax  %myvar; max variable length=7
    quiet : bool; default=False
        suppresses print messages
    dumplines: bool; default=False
        writes pyNastran_dump.bdf
    dictsort : bool; default=False
        writes pyNastran_dict.bdf
    dev : bool; default=False
        True : crashes if an Exception occurs
        False : doesn't crash; useful for running many tests
    """
    if not quiet:
        print('debug = %s' % debug)
    if dynamic_vars is None:
        dynamic_vars = {}

    # TODO: why do we need this?
    bdf_model = str(bdf_filename)
    if not quiet:
        print("bdf_model = %s" % bdf_model)
    if is_folder:
        bdf_model = os.path.join(test_path, folder, bdf_filename)

    assert os.path.exists(bdf_model), '%r doesnt exist' % bdf_model

    #print("before read bdf, Memory usage: %s (Mb) " % memory_usage_psutil())
    #print('before read bdf, Memory usage: %s (kb)' % resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    fem1 = BDF(debug=debug, log=None)
    if dynamic_vars:
        fem1.set_dynamic_syntax(dynamic_vars)

    fem1.log.info('starting fem1')
    sys.stdout.flush()
    fem2 = None
    diff_cards = []
    #try:
    out_model = run_fem1(fem1, bdf_model, mesh_form, xref, punch, sum_load, size, is_double, cid)
    fem2 = run_fem2(bdf_model, out_model, xref, punch, sum_load, size, is_double, reject,
                    debug=debug, log=None)
    diff_cards = compare(fem1, fem2, xref=xref, check=check, print_stats=print_stats)

    #except KeyboardInterrupt:
        #sys.exit('KeyboardInterrupt...sys.exit()')
    #except IOError:
        #pass
    #except AttributeError:  # only temporarily uncomment this when running lots of tests
        #pass
    #except SyntaxError:  # only temporarily uncomment this when running lots of tests
        #pass
    #except AssertionError:  # only temporarily uncomment this when running lots of tests
        #pass
    #except SystemExit:
        #sys.exit('sys.exit...')
    #except Exception:
        #exc_type, exc_value, exc_traceback = sys.exc_info()
        #print("\n")
        #traceback.print_exc(file=sys.stdout)
        #print msg
        #print("-" * 80)
        #raise

    print("-" * 80)
    return (fem1, fem2, diff_cards)


def run_fem1(fem1, bdf_model, mesh_form, xref, punch, sum_load, size, is_double, cid):
    check_path(bdf_model, 'bdf_model')
    try:
        if '.pch' in bdf_model:
            fem1.read_bdf(bdf_model, xref=False, punch=True)
        else:
            fem1.read_bdf(bdf_model, xref=xref, punch=punch)
    except Exception:
        print("failed reading %r" % bdf_model)
        raise
    #fem1.sumForces()

    out_model = bdf_model + 'v_out'
    if cid is not None and xref:
        fem1.resolveGrids(cid=cid)
    if mesh_form == 'combined':
        fem1.write_bdf(out_model, interspersed=True, size=size, is_double=is_double)
    elif mesh_form == 'separate':
        fem1.write_bdf(out_model, interspersed=False, size=size, is_double=is_double)
    else:
        msg = "mesh_form=%r; allowedForms=['combined','separate']" % mesh_form
        raise NotImplementedError(msg)
    #fem1.writeAsCTRIA3(out_model)
    return out_model


def run_fem2(bdf_model, out_model, xref, punch,
             sum_load, size, is_double,
             reject, debug=False, log=None):
    """
    Reads/writes the BDF to verify nothing has been lost

    Parameters
    ----------
    bdf_model : str
        the filename to run
    xref : bool
       xrefs
    punch : bool
       punches
    """
    assert os.path.exists(bdf_model), bdf_model
    assert os.path.exists(out_model), out_model
    fem2 = BDF(debug=debug, log=log)
    fem2.log.info('starting fem2')
    sys.stdout.flush()
    try:
        fem2.read_bdf(out_model, xref=xref, punch=punch)
    except Exception:
        print("failed reading %r" % out_model)
        raise

    #fem2.sumForces()
    #fem2.sumMoments()
    out_model2 = bdf_model + '_out2'
    fem2.write_bdf(out_model2, interspersed=True)
    #fem2.writeAsCTRIA3(out_model_2)
    os.remove(out_model2)
    return fem2


def divide(value1, value2):
    if value1 == value2:  # good for 0/0
        return 1.0
    else:
        try:
            v = value1 / float(value2)
        except ZeroDivisionError:
            v = 0.
    return v


def compare_card_count(fem1, fem2, print_stats=False):
    cards1 = fem1.card_count
    cards2 = fem2.card_count
    if print_stats:
        print(fem1.get_bdf_stats())
    else:
        fem1.get_bdf_stats()
    return compute_ints(cards1, cards2, fem1)


def compute_ints(cards1, cards2, fem1):
    card_keys1 = set(cards1.keys())
    card_keys2 = set(cards2.keys())
    all_keys = card_keys1.union(card_keys2)
    diff_keys1 = list(all_keys.difference(card_keys1))
    diff_keys2 = list(all_keys.difference(card_keys2))

    list_keys1 = list(card_keys1)
    list_keys2 = list(card_keys2)
    if diff_keys1 or diff_keys2:
        print(' diff_keys1=%s diff_keys2=%s' % (diff_keys1, diff_keys2))

    for key in sorted(all_keys):
        msg = ''
        if key in list_keys1:
            value1 = cards1[key]
        else:
            value1 = 0

        if key in list_keys2:
            value2 = cards2[key]
        else:
            value2 = 0

        diff = abs(value1 - value2)
        star = ' '
        if diff and key not in ['INCLUDE']:
            star = '*'
        if key not in fem1.cards_to_read:
            star = '-'

        factor1 = divide(value1, value2)
        factor2 = divide(value2, value1)
        factor_msg = ''
        if factor1 != factor2:
            factor_msg = 'diff=%s factor1=%g factor2=%g' % (diff, factor1,
                                                            factor2)
        msg += '  %skey=%-7s value1=%-7s value2=%-7s' % (star, key, value1,
                                                         value2) + factor_msg
        msg = msg.rstrip()
        print(msg)
    return list_keys1 + list_keys2


def compute(cards1, cards2):
    card_keys1 = set(cards1.keys())
    card_keys2 = set(cards2.keys())
    all_keys = card_keys1.union(card_keys2)
    diff_keys1 = list(all_keys.difference(card_keys1))
    diff_keys2 = list(all_keys.difference(card_keys2))

    list_keys1 = list(card_keys1)
    list_keys2 = list(card_keys2)
    msg = ''
    if diff_keys1 or diff_keys2:
        msg = 'diff_keys1=%s diff_keys2=%s' % (diff_keys1, diff_keys2)

    for key in sorted(all_keys):
        msg = ''
        if key in list_keys1:
            value1 = cards1[key]
        else:
            value2 = 0

        if key in list_keys2:
            value2 = cards2[key]
        else:
            value2 = 0

        if key == 'INCLUDE':
            msg += '    key=%-7s value1=%-7s value2=%-7s' % (key,
                                                             value1, value2)
        else:
            msg += '   *key=%-7s value1=%-7s value2=%-7s' % (key,
                                                             value1, value2)
        msg = msg.rstrip()
        print(msg)


def get_element_stats(fem1, fem2):
    """verifies that the various element methods work"""
    if 0:
        for (key, loads) in sorted(fem1.loads.items()):
            for load in loads:
                try:
                    all_loads = load.get_loads()
                    if not isinstance(all_loads, list):
                        raise TypeError('allLoads should return a list...%s'
                                        % (type(all_loads)))
                except Exception:
                    print("load statistics not available - load.type=%s "
                          "load.sid=%s" % (load.type, load.sid))
                    raise

    fem1._verify_bdf()

   # for (key, e) in sorted(fem1.elements.items()):
   #     try:
   #         e._verify()
   #         #if isinstance(e, RigidElement):
   #             #pass
   #         #elif isinstance(e, DamperElement):
   #             #b = e.B()
   #         #elif isinstance(e, SpringElement):
   #             #L = e.Length()
   #             #K = e.K()
   #             #pid = e.Pid()
   #         #elif isinstance(e, PointElement):
   #             #m = e.Mass()
   #             #c = e.Centroid()
   #     except Exception as exp:
   #         #print("e=\n",str(e))
   #         print("*stats - e.type=%s eid=%s  element=\n%s"
   #             % (e.type, e.eid, str(exp.args)))
   #     except AssertionError as exp:
   #         print("e=\n",str(e))
   #         #print("*stats - e.type=%s eid=%s  element=\n%s"
   #             #% (e.type, e.eid, str(exp.args)))
   #
   #         #raise


def get_matrix_stats(fem1, fem2):
    for (key, dmig) in sorted(fem1.dmigs.items()):
        try:
            if isinstance(dmig, NastranMatrix):
                dmig.get_matrix()
            else:
                print("statistics not available - "
                      "matrix.type=%s matrix.name=%s" % (dmig.type, dmig.name))
        except Exception:
            print("*stats - matrix.type=%s name=%s  matrix=\n%s"
                  % (dmig.type, dmig.name, str(dmig)))
            raise


def compare(fem1, fem2, xref=True, check=True, print_stats=True):
    diff_cards = compare_card_count(fem1, fem2, print_stats=print_stats)
    return
    #if xref and check:
    if check:
        fem1.mass_properties()
        get_element_stats(fem1, fem2)
        get_matrix_stats(fem1, fem2)
    compare_card_content(fem1, fem2)
    #compare_params(fem1, fem2)
    #print_points(fem1, fem2)
    return diffCards


def compare_params(fem1, fem2):
    compute(fem1.params, fem2.params)


def print_points(fem1, fem2):
    for (nid, n1) in sorted(fem1.nodes.items()):
        print("%s   xyz=%s  n1=%s  n2=%s" % (nid, n1.xyz, n1.Position(True),
                                             fem2.Node(nid).get_position()))
        break
    coord = fem1.Coord(5)
    print(coord)
    #print coord.Stats()


def main():
    from docopt import docopt
    msg = "Usage:\n"
    msg += "  test_bdf [-q] [-x] [-p] [-c] [-L] BDF_FILENAME\n" #
    msg += "  test_bdf [-q] [-x] [-p] [-c] [-L] [-d] BDF_FILENAME\n" #
    msg += "  test_bdf [-q] [-x] [-p] [-c] [-L] [-l] BDF_FILENAME\n" #
    msg += "  test_bdf [-q] [-p] [-r] BDF_FILENAME\n" #
    #msg += "  test_bdf [-q] [-p] [-o [<VAR=VAL>]...] BDF_FILENAME\n" #
    msg += '  test_bdf -h | --help\n'
    msg += '  test_bdf -v | --version\n'
    msg += '\n'

    msg += "Positional Arguments:\n"
    msg += "  BDF_FILENAME   path to BDF/DAT/NAS file\n"
    msg += '\n'

    msg += 'Options:\n'
    msg += '  -q, --quiet    prints debug messages (default=False)\n'
    msg += '  -x, --xref     disables cross-referencing and checks of the BDF.\n'
    msg += '                  (default=False -> on)\n'
    msg += '  -p, --punch    disables reading the executive and case control decks in the BDF\n'
    msg += '                 (default=False -> reads entire deck)\n'
    msg += '  -c, --check    disables BDF checks.  Checks run the methods on \n'
    msg += '                 every element/property to test them.  May fails if a \n'
    msg += '                 card is fully not supported (default=False)\n'
    msg += '  -l, --large    writes the BDF in large field, single precision format (default=False)\n'
    msg += '  -d, --double   writes the BDF in large field, double precision format (default=False)\n'
    msg += '  -L, --loads    sums the forces/moments on the model for the different subcases (default=False)\n'
    msg += '  -r, --reject   rejects all cards with the appropriate values applied (default=False)\n'
    #msg += '  -o <VAR_VAL>, --openmdao <VAR_VAL>   rejects all cards with the appropriate values applied;\n'
    #msg += '                 Uses the OpenMDAO %var syntax to replace it with value.\n'
    #msg += '                 So test_bdf -r var1=val1 var2=val2\n'
    msg += "\n"
    msg += "Info:\n"

    msg += '  -h, --help     show this help message and exit\n'
    msg += "  -v, --version  show program's version number and exit\n"

    if len(sys.argv) == 1:
        sys.exit(msg)

    ver = str(pyNastran.__version__)
    data = docopt(msg, version=ver)

    for key, value in sorted(data.items()):
        print("%-12s = %r" % (key.strip('--'), value))

    import time
    time0 = time.time()

    is_double = False
    if data['--double']:
        size = 16
        is_double = True
    elif data['--large']:
        size = 16
    else:
        size = 8
    print('is_double =', is_double)

    debug = True
    if data['--quiet']:
        debug = None
    #model = read_bdf(data['BDF_FILENAME'])
    run_bdf('.',
            data['BDF_FILENAME'],
            debug=debug,
            xref=not(data['--xref']),
            check=not(data['--check']),
            punch=data['--punch'],
            reject=data['--reject'],
            size=size,
            is_double=is_double,
            run_extract_bodies=False,
            #sum_load=data['--loads'],
    )
    print("total time:  %.2f sec" % (time.time() - time0))


if __name__ == '__main__':  # pragma: no cover
    main()
