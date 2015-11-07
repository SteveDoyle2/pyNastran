# pylint: disable=W0612
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
from six import iteritems, integer_types
import os
import sys
import numpy
from numpy import array
import warnings
warnings.simplefilter('always')

numpy.seterr(all='raise')
import traceback

from pyNastran.op2.op2 import OP2
from pyNastran.utils import print_bad_path
#from pyNastran.bdf.errors import CrossReferenceError, CardParseSyntaxError, DuplicateIDsError
from pyNastran.bdf.bdf import BDF, DLOAD
from pyNastran.bdf.cards.dmig import NastranMatrix
from pyNastran.bdf.bdf_replacer import BDFReplacer
from pyNastran.bdf.test.compare_card_content import compare_card_content

import pyNastran.bdf.test
test_path = pyNastran.bdf.test.__path__[0]


def run_all_files_in_folder(folder, debug=False, xref=True, check=True,
                            punch=False, cid=None, nastran=''):
    print("folder = %s" % folder)
    filenames = os.listdir(folder)
    run_lots_of_files(filenames, debug=debug, xref=xref, check=check,
                      punch=punch, cid=cid, nastran=nastran)


def run_lots_of_files(filenames, folder='', debug=False, xref=True, check=True,
                      punch=False, cid=None, nastran='',
                      size=None, is_double=None, post=None, sum_load=True):
    """
    Runs multiple BDFs

    Parameters
    ----------
    folder : str
        the folder where the bdf_filename is
    filenames : List[str]
        the bdf files to analyze
    debug : bool, optional
        run with debug logging (default=False)
    xref : bool / str / List[bool/str], optional
        True : cross reference the model
        False  : don't cross reference the model
        'safe' : do safe cross referencing
    check : bool / List[bool], optional
        validate cards for things like mass, area, etc. (default=True)
    punch : bool / List[bool], optional
        this is a PUNCH file (no executive/case control decks; default=False)
    cid : int / None, optional
        convert the model grids to an alternate coordinate system (default=None; no conversion)
    size : int / List[int], optional
        The field width of the model (8/16)
    is_double : bool / List[bool], optional
        Is this a double precision model?
            True : size = 16
            False : six = {8, 16}
    nastran : str, optional
        the path to nastran (default=''; no analysis)
    post : int / List[int], optional
        the PARAM,POST,value to run
    sum_load : bool; default=True
        should the loads be summed

    Usage
    -----
    All control lists must be the same length.
    You can run xref=True and xref=False with:

    .. python ::

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

    if post is None:
        posts = [-1]
    elif isinstance(post, integer_types):
        posts = [post]
    else:
        posts = post

    size_doubles_post = []
    print('posts=%s' % posts)
    for size, is_double, post in zip(sizes, is_doubles, posts):
        size_doubles_post.append((size, is_double, post))

    #debug = True
    filenames2 = []
    diff_cards = []
    for filename in filenames:
        if filename.endswith(('.bdf', '.dat', '.nas')):
            filenames2.append(filename)

    failed_files = []
    n = 1
    for filename in filenames2:
        abs_filename = os.path.abspath(os.path.join(folder, filename))
        if folder != '':
            print("filename = %s" % abs_filename)
        is_passed = False
        try:
            for size, is_double, post in size_doubles_post:
                fem1, fem2, diff_cards2 = run_bdf(folder, filename, debug=debug,
                                                  xref=xref, check=check, punch=punch,
                                                  cid=cid, isFolder=True, dynamic_vars={},
                                                  nastran=nastran, size=size, is_double=is_double,
                                                  post=post, sum_load=sum_load)
                del fem1
                del fem2
            diff_cards += diff_cards
            is_passed = True
        except KeyboardInterrupt:
            sys.exit('KeyboardInterrupt...sys.exit()')
        #except IOError:
            #pass
        #except RuntimeError:  # only temporarily uncomment this when running lots of tests
            #pass
        #except AttributeError:  # only temporarily uncomment this when running lots of tests
            #pass
        #except SyntaxError:  # only temporarily uncomment this when running lots of tests
            #pass
        except SystemExit:
            sys.exit('sys.exit...')
        except:
            traceback.print_exc(file=sys.stdout)
            #raise
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


def memory_usage_psutil():
    # return the memory usage in MB
    try:
        import psutil
    except ImportError:
        return '???'
    process = psutil.Process(os.getpid())
    mem = process.get_memory_info()[0] / float(2 ** 20)
    return mem


def run_bdf(folder, bdf_filename, debug=False, xref=True, check=True, punch=False,
            cid=None, meshForm='combined', isFolder=False, print_stats=False,
            sum_load=False, size=8, is_double=False,
            reject=False, stop=False, nastran='', post=-1, dynamic_vars=None):
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
    meshForm : str, optional, {'combined', 'separate'}
        'combined' : interspersed=True
        'separate' : interspersed=False
    isFolder : bool, optional
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
            False : six = {8, 16}
    reject : bool, optional
        True : all the cards are rejected
        False : the model is read
    nastran : str, optional
        the path to nastran (default=''; no analysis)
    post : int, optional
        the PARAM,POST,value to run
    dynamic vars : dict[str]=int / float / str / None
        support OpenMDAO syntax  %myvar; max variable length=7
    """
    if dynamic_vars is None:
        dynamic_vars = {}

    # TODO: why do we need this?
    bdfModel = str(bdf_filename)
    print("bdfModel = %s" % bdfModel)
    if isFolder:
        bdfModel = os.path.join(test_path, folder, bdf_filename)

    assert os.path.exists(bdfModel), '%r doesnt exist' % bdfModel

    if reject:
        fem1 = BDFReplacer(bdfModel + '.rej', debug=debug, log=None)
    else:
        fem1 = BDF(debug=debug, log=None)
    fem1.set_error_storage(nparse_errors=100, stop_on_parsing_error=True,
                           nxref_errors=100, stop_on_xref_error=True)
    if dynamic_vars:
        fem1.set_dynamic_syntax(dynamic_vars)

    fem1.log.info('starting fem1')
    sys.stdout.flush()
    fem2 = None
    diffCards = []

    try:
        #nastran = 'nastran scr=yes bat=no old=no news=no '
        nastran = ''
        #try:
        outModel = run_fem1(fem1, bdfModel, meshForm, xref, punch, sum_load, size, is_double, cid)
        if stop:
            print('***stopping')
            return fem1, None, None
        fem2 = run_fem2(bdfModel, outModel, xref, punch, sum_load, size, is_double, reject, debug=debug, log=None)
        diffCards = compare(fem1, fem2, xref=xref, check=check, print_stats=print_stats)
        test_get_cards_by_card_types(fem2)
        #except:
            #return 1, 2, 3

        run_nastran(bdfModel, nastran, post, size, is_double)

    except KeyboardInterrupt:
        sys.exit('KeyboardInterrupt...sys.exit()')
    #except IOError:  # only temporarily uncomment this when running lots of tests
        #pass
    #except CardParseSyntaxError:  # only temporarily uncomment this when running lots of tests
        #print('failed test because CardParseSyntaxError...ignoring')
    #except DuplicateIDsError:  # only temporarily uncomment this when running lots of tests
        #print('failed test because DuplicateIDsError...ignoring')
    #except RuntimeError:  # only temporarily uncomment this when running lots of tests
        #if 'GRIDG' in fem1.card_count:
            #print('failed test because mesh adaption (GRIDG)...ignoring')
            #raise
    #except AttributeError:  # only temporarily uncomment this when running lots of tests
        #pass
    #except SyntaxError:  # only temporarily uncomment this when running lots of tests
        #pass
    #except AssertionError:  # only temporarily uncomment this when running lots of tests
        #pass
    except SystemExit:
        sys.exit('sys.exit...')
    except:
        #exc_type, exc_value, exc_traceback = sys.exc_info()
        #print "\n"
        traceback.print_exc(file=sys.stdout)
        #print msg
        print("-" * 80)
        raise

    print("-" * 80)
    return (fem1, fem2, diffCards)


def run_nastran(bdf_model, nastran, post=-1, size=8, is_double=False):
    """
    Verifies that a valid bdf was written by running nastran and parsing
    the OP2.  Many cards do not support double precision and since there
    is no list, a test is necessary.
    """
    if nastran:
        dirname = os.path.dirname(bdf_model)
        basename = os.path.basename(bdf_model).split('.')[0]

        f04_model = os.path.join(dirname, 'out_%s.f04' % basename)
        f06_model = os.path.join(dirname, 'out_%s.f06' % basename)
        op2_model = os.path.join(dirname, 'out_%s.f06' % basename)
        log_model = os.path.join(dirname, 'out_%s.log' % basename)
        xdb_model = os.path.join(dirname, 'out_%s.xdb' % basename)
        pch_model = os.path.join(dirname, 'out_%s.pch' % basename)
        asm_model = os.path.join(dirname, 'out_%s.asm' % basename)
        master_model = os.path.join(dirname, 'out_%s.master' % basename)
        #op2_model = os.path.join(dirname, 'out_%s.op2' % basename)

        #cwd = os.getcwd()
        cwd = dirname
        bdf_model2 = os.path.join(cwd, 'out_%s.bdf' % basename)
        op2_model2 = os.path.join(cwd, 'out_%s.op2' % basename)
        #f06_model2 = os.path.join(cwd, 'out_%s.f06' % basename)
        print(bdf_model2)
        #if os.path.exists(bdf_model2):
            #os.remove(bdf_model2)

        # make sure we're writing an OP2
        bdf = BDF()
        bdf.read_bdf(bdf_model)
        if 'POST' in bdf.params:
            param_post = bdf.params['POST']
            #print('post = %s' % post)
            param_post.update_values(value1=post)
            #print('post = %s' % post)
        else:
            card = ['PARAM', 'POST', post]
            bdf.add_card(card, 'PARAM', is_list=True)
        bdf.write_bdf(bdf_model2, size=size, is_double=is_double)

        #os.rename(outModel, outModel2)
        if not os.path.exists(f06_model):
            os.system(nastran + bdf_model2)
        for fnamei in [f04_model, log_model, xdb_model, pch_model, asm_model, master_model]:
            if os.path.exists(fnamei):
                os.remove(fnamei)
        op2 = OP2()
        if not os.path.exists(op2_model):
            raise RuntimeError('%s failed' % op2_model)
        op2.read_op2(op2_model2)
        print(op2.get_op2_stats())

def run_fem1(fem1, bdfModel, meshForm, xref, punch, sum_load, size, is_double, cid):
    """
    Reads/writes the BDF

    Parameters
    ----------
    fem1 : BDF()
        The BDF object
    bdfModel : str
        The root path of the bdf filename
    meshForm : str {combined, separate}
        'combined' : interspersed=True
        'separate' : interspersed=False
    xref : bool
        The xref mode
    punch : bool
        punch flag
    sum_load : bool
        static load sum flag
    size : int, {8, 16}
        size flag
    is_double : bool
        double flag
    cid : int / None
        cid flag
    """
    assert os.path.exists(bdfModel), print_bad_path(bdfModel)
    try:
        if '.pch' in bdfModel:
            fem1.read_bdf(bdfModel, xref=False, punch=True)
        else:
            fem1.read_bdf(bdfModel, xref=False, punch=punch)
            if xref:
                #fem1.uncross_reference()
                fem1.cross_reference()
                #fem1.uncross_reference()
                #fem1.cross_reference()
    except:
        print("failed reading %r" % bdfModel)
        raise
    #fem1.sumForces()

    if fem1._auto_reject:
        out_model = bdfModel + '.rej'
    else:
        out_model = bdfModel + '_out'
        if cid is not None and xref:
            fem1.resolve_grids(cid=cid)

        if meshForm == 'combined':
            fem1.write_bdf(out_model, interspersed=False, size=size, is_double=is_double)
        elif meshForm == 'separate':
            fem1.write_bdf(out_model, interspersed=False, size=size, is_double=is_double)
        else:
            msg = "meshForm=%r; allowedForms=['combined','separate']" % meshForm
            raise NotImplementedError(msg)
        #fem1.writeAsCTRIA3(out_model)

    fem1._get_maps()
    return out_model


def run_fem2(bdfModel, out_model, xref, punch,
             sum_load, size, is_double,
             reject, debug=False, log=None):
    """
    Reads/writes the BDF to verify nothing has been lost

    Parameters
    ----------
    bdfModel
    out_model
    xref : bool
       xrefs
    punch : bool
       punches
    sum_load : bool
       sums static load
    size : int
    is_double : bool
    reject : bool
        True : rejects the cards
    debug : bool
        debugs
    log : logger / None
        ignored
    """
    assert os.path.exists(bdfModel), bdfModel
    assert os.path.exists(out_model), out_model

    if reject:
        fem2 = BDFReplacer(bdfModel + '.rej', debug=debug, log=None)
    else:
        fem2 = BDF(debug=debug, log=None)
    fem2.log.info('starting fem2')
    sys.stdout.flush()
    try:
        fem2.read_bdf(out_model, xref=xref, punch=punch)
    except:
        print("failed reading %r" % out_model)
        raise

    outModel2 = bdfModel + '_out2'

    if sum_load:
        p0 = array([0., 0., 0.])

        subcase_keys = fem2.case_control_deck.get_subcase_list()
        subcases = fem2.subcases

        sol_200_map = fem2.case_control_deck.sol_200_map
        sol_base = fem2.sol
        for isubcase in subcase_keys[1:]:  # drop isubcase = 0
            subcase = subcases[isubcase]
            if sol_base == 200:
                analysis = subcase.get_parameter('ANALYSIS')[0]
                sol = sol_200_map[analysis]
            else:
                sol = sol_base

            if subcase.has_parameter('METHOD'):
                method_id = subcase.get_parameter('METHOD')[0]
                method = fem2.methods[method_id]
                assert sol in [5, 76, 101, 103, 105, 106, 107, 108, 110, 111, 112, 144, 145, 146, 187], 'sol=%s METHOD' % sol
            if subcase.has_parameter('CMETHOD'):
                method_id = subcase.get_parameter('CMETHOD')[0]
                method = fem2.cMethods[method_id]
                assert sol in [107, 110, 145], 'sol=%s CMETHOD' % sol

            if subcase.has_parameter('LOAD'):
                loadcase_id = fem2.case_control_deck.get_subcase_parameter(isubcase, 'LOAD')[0]
                F, M = fem2.sum_forces_moments(p0, loadcase_id, include_grav=False)
                print('  isubcase=%i F=%s M=%s' % (isubcase, F, M))
                assert sol in [1, 5, 24, 61, 64, 66, 101, 103, 105, 106, 107, 108, 110, 112,
                               144, 145, 153, 400, 601], 'sol=%s LOAD' % sol
            else:
                # print('is_load =', subcase.has_parameter('LOAD'))
                pass

            if subcase.has_parameter('FREQUENCY'):
                freq_id = subcase.get_parameter('FREQUENCY')[0]
                freq = fem2.frequencies[freq_id]
                assert sol in [26, 68, 76, 78, 88, 108, 101, 111, 112, 118, 146], 'sol=%s FREQUENCY' % sol
                # print(freq)

            # if subcase.has_parameter('LSEQ'):
                # lseq_id = subcase.get_parameter('LSEQ')[0]
                # lseq = fem2.loads[lseq_id]
                # assert sol in [], sol
                # print(lseq)
            if subcase.has_parameter('SPC'):
                spc_id = subcase.get_parameter('SPC')[0]
                fem2.get_spcs(spc_id)

            if subcase.has_parameter('DLOAD'):
                assert sol in [26, 68, 76, 78, 88, 99, 103, 108, 109, 111, 112, 118, 129, 146,
                               153, 159, 400, 601], 'sol=%s DLOAD' % sol
                if subcase.has_parameter('LOADSET'):
                    raise NotImplementedError('LOADSET & DLOAD -> LSEQ')
                if subcase.has_parameter('IC'):
                    raise NotImplementedError('IC & DLOAD -> TIC')

                # DLOAD (case)   -> dynamic loads -> DLOAD, RLOAD1, RLOAD2, TLOAD1, TLOAD2, ACSRCE
                # LOADSET (case) -> static load sequence - > LSEQ
                # LSEQ (bulk)    -> sequence of static load sets
                # IC (case)      -> points to TIC (initial conditions)
                #
                # TYPE 0 (LOAD)
                #  - no LOADSET -> DAREA, static, thermal load entry
                #  -    LOADSET -> static, thermal loads as specified by LSEQ
                # TYPE 1/2/3 (DISP, VELO, ACCE)
                #  - no LOADSET -> SPCD
                #  -    LOADSET -> SPCDs as specified by LSEQ
                dload_id = subcase.get_parameter('DLOAD')[0]
                if dload_id in fem2.dloads:
                    dload = fem2.dloads[dload_id]
                else:
                    dload = fem2.dload_entries[dload_id]
                # dload = DLOAD()
                # print(dload)
                # for
                # loads, sf = dload.get_loads()
                scale_factors2 = []
                loads2 = []
                for load in dload:
                    # print('DLOAD\n', load)
                    if isinstance(load, DLOAD):
                        scale = load.scale
                        scale_factors = []
                        loads = []
                        # scale_factors, loads = load.get_reduced_loads()
                        for load, sf in zip(load.loadIDs, load.scaleFactors):
                            if isinstance(load, list):
                                for loadi in load:
                                    assert not isinstance(loadi, list), loadi
                                    scale_factors.append(scale * sf)
                                    loads.append(loadi)
                            else:
                                scale_factors.append(scale * sf)
                                assert not isinstance(load, list), load
                                loads.append(load)
                        scale_factors2 += scale_factors
                        loads2 += loads
                    else:
                        scale_factors2.append(1.)
                        loads2.append(load)

                if sol in [108, 111]:  # direct frequency, modal frequency
                    for load2, scale_factor in zip(loads2, scale_factors2):
                        # for
                        #print(load2)
                        F = load2.get_load_at_freq(100.) * scale_factor
                elif sol in [109, 129]:  # direct transient (time linear), time nonlinear
                    for load2, scale_factor in zip(loads2, scale_factors2):
                        # for
                        #print(load2)
                        F = load2.get_load_at_time(0.) * scale_factor
                ### 111
                else:
                    fem2.log.debug('solution=%s; DLOAD is not supported' % sol)

                # print(loads)

    fem2.write_bdf(outModel2, interspersed=False, size=size, is_double=is_double)
    #fem2.writeAsCTRIA3(outModel2)
    os.remove(outModel2)
    return fem2


def divide(value1, value2):
    """
    Used to divide the number of cards to check that nothing was lost.
    Handles division by 0 by returning 0, which is the reciprocal.
    """
    if value1 == value2:  # good for 0/0
        return 1.0
    else:
        try:
            v = value1 / float(value2)
        except ZeroDivisionError:
            v = 0.
    return v


def test_get_cards_by_card_types(model):
    """
    Verifies the ``model.get_cards_by_card_types`` method works
    """
    # setup to remove hackish cards
    card_types = list(model.card_count.keys())
    removed_cards = []
    for card_type in ['ENDDATA', 'INCLUDE', 'JUNK']:
        if card_type in model.card_count:
            removed_cards.append(card_type)
    for removed_card in removed_cards:
        card_types.remove(removed_card)

    removed_cards = []
    for card_type in card_types:
        if card_type not in model.cards_to_read:
            try:
                removed_cards.append(card_type)
                #print('removed %s' % card_type)
            except ValueError:
                msg = 'card_type=%s cant be removed' % card_type
                raise ValueError(msg)
    for removed_card in removed_cards:
        card_types.remove(removed_card)

    # we now have a list of card types we would like to extract
    # we'll get the associated cards
    card_dict = model.get_cards_by_card_types(card_types, reset_type_to_slot_map=False)
    for card_type, cards in iteritems(card_dict):
        for card in cards:
            assert card_type == card.type, 'this should never crash here...card_type=%s card.type=%s' % (card_type, card.type)


def compare_card_count(fem1, fem2, print_stats=False):
    """
    Checks that no cards from fem1 are lost when we write fem2
    """
    cards1 = fem1.card_count
    cards2 = fem2.card_count
    for key in cards1:
        if key != key.upper():
            raise RuntimeError('Proper capitalization wasnt determined')
    if print_stats:
        print(fem1.get_bdf_stats())
    else:
        fem1.get_bdf_stats()
    return compute_ints(cards1, cards2, fem1)


def compute_ints(cards1, cards2, fem1):
    """
    computes the difference / ratio / inverse-ratio between
    fem1 and fem2 to verify the number of card are the same:

    Example
    -------

    name   fem1  fem2  diff  ratio  1/ratio
    ====   ====  ====  ==== ======  =======
    GRID      1     1     1     1.       1.
    *SPOINT  10     1     9    10.      0.1

    The * indicates a change, which may or may not be a problem.
    """
    cardKeys1 = set(cards1.keys())
    cardKeys2 = set(cards2.keys())
    allKeys = cardKeys1.union(cardKeys2)
    diffKeys1 = list(allKeys.difference(cardKeys1))
    diffKeys2 = list(allKeys.difference(cardKeys2))

    listKeys1 = list(cardKeys1)
    listKeys2 = list(cardKeys2)
    if diffKeys1 or diffKeys2:
        print(' diffKeys1=%s diffKeys2=%s' % (diffKeys1, diffKeys2))

    for key in sorted(allKeys):
        msg = ''
        if key in listKeys1:
            value1 = cards1[key]
        else:
            value1 = 0

        if key in listKeys2:
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
    #return listKeys1 + listKeys2
    return diffKeys1 + diffKeys2


def compute(cards1, cards2):
    """
    Computes the difference between two dictionaries to data is the same
    """
    cardKeys1 = set(cards1.keys())
    cardKeys2 = set(cards2.keys())
    allKeys = cardKeys1.union(cardKeys2)
    diffKeys1 = list(allKeys.difference(cardKeys1))
    diffKeys2 = list(allKeys.difference(cardKeys2))

    listKeys1 = list(cardKeys1)
    listKeys2 = list(cardKeys2)
    msg = ''
    if diffKeys1 or diffKeys2:
        msg = 'diffKeys1=%s diffKeys2=%s' % (diffKeys1, diffKeys2)

    for key in sorted(allKeys):
        msg = ''
        if key in listKeys1:
            value1 = cards1[key]
        else:
            value2 = 0

        if key in listKeys2:
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
    for (key, loads) in sorted(iteritems(fem1.loads)):
        for load in loads:
            try:
                allLoads = load.get_loads()
                if not isinstance(allLoads, list):
                    raise TypeError('allLoads should return a list...%s'
                                    % (type(allLoads)))
            except:
                print("load statistics not available - load.type=%s "
                      "load.sid=%s" % (load.type, load.sid))
                raise

    fem1._verify_bdf()

    mass, cg, I = fem1.mass_properties(reference_point=None, sym_axis=None)
    print("mass =", mass)
    print("cg   =", cg)
    #print("I    =", I)


def get_matrix_stats(fem1, fem2):
    """
    Verifies the dmig.get_matrix() method works.
    """
    for (key, dmig) in sorted(iteritems(fem1.dmigs)):
        try:
            if isinstance(dmig, NastranMatrix):
                dmig.get_matrix()
            else:
                print("statistics not available - "
                      "matrix.type=%s matrix.name=%s" % (dmig.type, dmig.name))
        except:
            print("*stats - matrix.type=%s name=%s  matrix=\n%s"
                  % (dmig.type, dmig.name, str(dmig)))
            raise


def compare(fem1, fem2, xref=True, check=True, print_stats=True):
    diffCards = compare_card_count(fem1, fem2, print_stats=print_stats)
    if xref and check:
        get_element_stats(fem1, fem2)
        get_matrix_stats(fem1, fem2)
    compare_card_content(fem1, fem2)
    #compare_params(fem1, fem2)
    #print_points(fem1, fem2)
    return diffCards


def compare_params(fem1, fem2):
    compute(fem1.params, fem2.params)


def print_points(fem1, fem2):
    for nid, node in sorted(iteritems(fem1.nodes)):
        print("%s   xyz=%s  n1=%s  n2=%s" % (nid, node.xyz, node.get_position(True),
                                             fem2.Node(nid).get_position()))
        break
    coord = fem1.Coord(5)
    print(coord)
    #print coord.Stats()


def main():
    """
    The main function for the command line ``test_bdf`` script.
    """
    from docopt import docopt
    msg = "Usage:\n"
    msg += "  test_bdf [-q] [-x] [-p] [-c] [-L] [-f] BDF_FILENAME\n" #
    msg += "  test_bdf [-q] [-x] [-p] [-c] [-L] [-d] [-f] BDF_FILENAME\n" #
    msg += "  test_bdf [-q] [-x] [-p] [-c] [-L] [-l] [-f] BDF_FILENAME\n" #
    msg += "  test_bdf [-q] [-p] [-r] [-f] BDF_FILENAME\n" #
    msg += "  test_bdf [-q] [-x] [-p] [-s] [-f] BDF_FILENAME\n" #

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
    msg += '  -L, --loads    Disables forces/moments summation for the different subcases (default=False)\n'
    msg += '  -r, --reject   rejects all cards with the appropriate values applied (default=False)\n'
    msg += '  -f, --profile  Profiles the code (default=False)\n'
    msg += '  -s, --stop     Stop after first read/write (default=False)\n'
    #msg += '  -o <VAR_VAL>, --openmdao <VAR_VAL>   rejects all cards with the appropriate values applied;\n'
    #msg += '                 Uses the OpenMDAO %var syntax to replace it with value.\n'
    #msg += '                 So test_bdf -r var1=val1 var2=val2\n'

    msg += '  -h, --help     show this help message and exit\n'
    msg += "  -v, --version  show program's version number and exit\n"

    if len(sys.argv) == 1:
        sys.exit(msg)

    ver = str(pyNastran.__version__)
    data = docopt(msg, version=ver)

    for key, value in sorted(iteritems(data)):
        print("%-12s = %r" % (key.strip('--'), value))

    import time
    t0 = time.time()

    is_double = False
    if data['--double']:
        size = 16
        is_double = True
    elif data['--large']:
        size = 16
    else:
        size = 8

    #print(data)
    if data['--profile']:
        #import cProfile
        import pstats

        import cProfile
        prof = cProfile.Profile()
        prof.runcall(
            run_bdf,
                '.',
                data['BDF_FILENAME'],
                debug=not(data['--quiet']),
                xref=not(data['--xref']),
                # xref_safe=data['--xref_safe'],
                check=not(data['--check']),
                punch=data['--punch'],
                reject=data['--reject'],
                size=size,
                is_double=is_double,
                sum_load=not data['--loads'],
                stop=data['--stop'],
            )
        prof.dump_stats('bdf.profile')

        stats = pstats.Stats("bdf.profile")
        stats.sort_stats('tottime')  # time in function
        #stats.sort_stats('cumtime')  # time in function & subfunctions
        stats.strip_dirs()
        stats.print_stats(40)

        #retval = prof.runcall(self.method_actual, *args, **kwargs)
        #print(prof.dump_stats(datafn))
        #cProfile.runctx(
            #code,
               #None, # globs
               #None,
               #'junk.stats',
               #1) # sort

        #p = pstats.Stats('restats')
        #p.strip_dirs().sort_stats(-1).print_stats()
    else:
        run_bdf(
            '.',
            data['BDF_FILENAME'],
            debug=not(data['--quiet']),
            xref=not(data['--xref']),
            # xref_safe=data['--xref_safe'],
            check=not(data['--check']),
            punch=data['--punch'],
            reject=data['--reject'],
            size=size,
            is_double=is_double,
            sum_load=not data['--loads'],
            stop=data['--stop'],
        )
    print("total time:  %.2f sec" % (time.time() - t0))


if __name__ == '__main__':  # pragma: no cover
    main()
