# pylint: disable=W0612
"""
``test_bdf`` runs multiple checks on a BDF in order to make sure that:
  - no data is lost on IO
  - card field types are correct (e.g. node_ids are integers)
  - various card methods (e.g. Area) work correctly

As such, ``test_bdf`` is very useful for debugging models.
"""
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
import os
import sys
import traceback
import warnings
from six import iteritems
import numpy as np
warnings.simplefilter('always')

np.seterr(all='raise')

from pyNastran.op2.op2 import OP2
from pyNastran.utils import print_bad_path, integer_types
from pyNastran.bdf.errors import (
    CrossReferenceError, CardParseSyntaxError, DuplicateIDsError, MissingDeckSections)
from pyNastran.bdf.bdf import BDF, DLOAD, read_bdf
from pyNastran.bdf.mesh_utils.extract_bodies import extract_bodies
from pyNastran.bdf.cards.dmig import NastranMatrix
from pyNastran.bdf.test.compare_card_content import compare_card_content
from pyNastran.bdf.mesh_utils.convert import convert
from pyNastran.bdf.mesh_utils.remove_unused import remove_unused

import pyNastran.bdf.test
test_path = pyNastran.bdf.test.__path__[0]

class DisabledCardError(RuntimeError):
    pass

def run_all_files_in_folder(folder, debug=False, xref=True, check=True,
                            punch=False, cid=None, nastran=''):
    """runs all the BDFs in a given folder"""
    print("folder = %s" % folder)
    filenames = os.listdir(folder)
    run_lots_of_files(filenames, debug=debug, xref=xref, check=check,
                      punch=punch, cid=cid, nastran=nastran)


def run_lots_of_files(filenames, folder='', debug=False, xref=True, check=True,
                      punch=False, cid=None, nastran='', encoding=None,
                      size=None, is_double=None, post=None, sum_load=True, dev=True,
                      crash_cards=None):
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
    dev : bool; default=True
        True : crashes if an Exception occurs
        False : doesn't crash; useful for running many tests
    crash_cards : List[str, str, ...]
        list of cards that are invalid and automatically crash the run
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
        if (filename.endswith(('.bdf', '.dat', '.nas')) and
                'pyNastran_crash' not in filename and
                'skin_file' not in filename):
            filenames2.append(filename)

    failed_files = []
    npass = 1
    nfailed = 1
    for filename in filenames2:
        abs_filename = os.path.abspath(os.path.join(folder, filename))
        if folder != '':
            print("filename = %s" % abs_filename)
        is_passed = False
        try:
            for size, is_double, post in size_doubles_post:
                fem1, fem2, diff_cards2 = run_bdf(folder, filename, debug=debug,
                                                  xref=xref, check=check, punch=punch,
                                                  cid=cid, encoding=encoding,
                                                  is_folder=True, dynamic_vars={},
                                                  nastran=nastran, size=size, is_double=is_double,
                                                  nerrors=0,
                                                  post=post, sum_load=sum_load, dev=dev,
                                                  crash_cards=crash_cards,
                                                  run_extract_bodies=False)
                del fem1
                del fem2
            diff_cards += diff_cards
            is_passed = True
        except KeyboardInterrupt:
            sys.exit('KeyboardInterrupt...sys.exit()')
        except DisabledCardError:
            #if dev:
                #pass
            raise
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
            sys.stderr.write('%i  %s' % (npass, abs_filename))
            npass += 1
        else:
            sys.stderr.write('*%s ' % nfailed + abs_filename)
            nfailed += 1
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
            cid=None, mesh_form='combined', is_folder=False, print_stats=False,
            encoding=None, sum_load=True, size=8, is_double=False,
            stop=False, nastran='', post=-1, dynamic_vars=None,
            quiet=False, dumplines=False, dictsort=False, run_extract_bodies=False,
            nerrors=0, dev=False, crash_cards=None):
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
    sum_load : bool; default=True
        Sum the static loads (doesn't work for frequency-based loads)
    size : int, optional, {8, 16}
        The field width of the model
    is_double : bool, optional
        Is this a double precision model?
            True : size = 16
            False : six = {8, 16}
    stop : bool; default=False
        stop reading the first BDF
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
    run_extract_bodies : bool; default=False
        isolate the fem bodies; typically 1 body; code is still buggy
    dev : bool; default=False
        True : crashes if an Exception occurs
        False : doesn't crash; useful for running many tests
    """
    if not quiet:
        print('debug = %s' % debug)
    if dynamic_vars is None:
        dynamic_vars = {}
    if crash_cards is None:
        crash_cards = []

    # TODO: why do we need this?
    bdf_model = str(bdf_filename)
    if not quiet:
        print("bdf_model = %s" % bdf_model)
    if is_folder:
        bdf_model = os.path.join(test_path, folder, bdf_filename)

    model, ext = os.path.splitext(bdf_model)
    out_model = '%s.test_bdf%s' % (model, ext)

    fem1, fem2, diff_cards = run_and_compare_fems(
        bdf_model, out_model, debug=debug, xref=xref, check=check,
        punch=punch, cid=cid, mesh_form=mesh_form,
        print_stats=print_stats, encoding=encoding,
        sum_load=sum_load, size=size, is_double=is_double,
        stop=stop, nastran=nastran, post=post,
        dynamic_vars=dynamic_vars,
        quiet=quiet, dumplines=dumplines, dictsort=dictsort,
        nerrors=nerrors, dev=dev, crash_cards=crash_cards,
        run_extract_bodies=run_extract_bodies,
    )
    return fem1, fem2, diff_cards

def run_and_compare_fems(
        bdf_model, out_model, debug=False, xref=True, check=True,
        punch=False, cid=None, mesh_form='combined',
        print_stats=False, encoding=None,
        sum_load=True, size=8, is_double=False,
        stop=False, nastran='', post=-1, dynamic_vars=None,
        quiet=False, dumplines=False, dictsort=False,
        nerrors=0, dev=False, crash_cards=None,
        run_extract_bodies=False,
    ):
    """runs two fem models and compares them"""
    assert os.path.exists(bdf_model), '%r doesnt exist' % bdf_model

    fem1 = BDF(debug=debug, log=None)

    fem1.set_error_storage(nparse_errors=nerrors, stop_on_parsing_error=True,
                           nxref_errors=nerrors, stop_on_xref_error=True)
    if dynamic_vars:
        fem1.set_dynamic_syntax(dynamic_vars)

    if not quiet:
        fem1.log.info('starting fem1')
    sys.stdout.flush()
    fem2 = None
    diff_cards = []

    try:
        #nastran_cmd = 'nastran scr=yes bat=no old=no news=no '
        nastran_cmd = ''
        #try:

        fem1 = run_fem1(fem1, bdf_model, out_model, mesh_form, xref, punch, sum_load,
                        size, is_double, cid,
                        run_extract_bodies=run_extract_bodies,
                        encoding=encoding, crash_cards=crash_cards)
        if stop:
            if not quiet:
                print('card_count:')
                print('-----------')
                for card_name, card_count in sorted(iteritems(fem1.card_count)):
                    print('key=%-8s value=%s' % (card_name, card_count))
            return fem1, None, None
        fem2 = run_fem2(bdf_model, out_model, xref, punch, sum_load, size, is_double, mesh_form,
                        encoding=encoding, debug=debug, quiet=quiet)

        diff_cards = compare(fem1, fem2, xref=xref, check=check,
                             print_stats=print_stats, quiet=quiet)
        test_get_cards_by_card_types(fem2)
        #except:
            #return 1, 2, 3

        run_nastran(bdf_model, nastran_cmd, post, size, is_double)

    except KeyboardInterrupt:
        sys.exit('KeyboardInterrupt...sys.exit()')
    except IOError:  # only temporarily uncomment this when running lots of tests
        if not dev:
            raise
    except CardParseSyntaxError:  # only temporarily uncomment this when running lots of tests
        if not dev:
            raise
        print('failed test because CardParseSyntaxError...ignoring')
    except MissingDeckSections:
        if not dev:
            raise
        print('failed test because MissingDeckSections...ignoring')
    except DuplicateIDsError as e:
        # only temporarily uncomment this when running lots of tests
        if 'GRIDG' in fem1.card_count or 'CGEN' in fem1.card_count or 'SPCG' in fem1.card_count:
            print('failed test because mesh adaption (GRIDG,CGEN,SPCG)...ignoring')
            print(e)
        elif not dev:
            raise
        else:
            print('failed test because DuplicateIDsError...ignoring')
    except DisabledCardError as e:
        if not dev:
            raise
    except RuntimeError as e:
        # only temporarily uncomment this when running lots of tests
        if not dev:
            raise
        if 'GRIDG' in fem1.card_count or 'CGEN' in fem1.card_count or 'SPCG' in fem1.card_count:
            print('failed test because mesh adaption (GRIDG,CGEN,SPCG)...ignoring')
            print(e)
        else:
            raise
    #except AttributeError:  # only temporarily uncomment this when running lots of tests
        #pass
    except SyntaxError as e:
        # only temporarily uncomment this when running lots of tests
        if not dev:
            raise
        if 'GRIDG' in fem1.card_count or 'CGEN' in fem1.card_count or 'SPCG' in fem1.card_count:
            print('failed test because mesh adaption (GRIDG,CGEN,SPCG)...ignoring')
            print(e)
        else:
            raise
    except KeyError as e:  # only temporarily uncomment this when running lots of tests
        if not dev:
            raise
        if 'GRIDG' in fem1.card_count or 'CGEN' in fem1.card_count or 'SPCG' in fem1.card_count:
            print('failed test because mesh adaption (GRIDG,CGEN,SPCG)...ignoring')
            print(e)
        else:
            raise
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

    if not quiet:
        print("-" * 80)
    return (fem1, fem2, diff_cards)


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
        bdf = BDF(debug=False)
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

def run_fem1(fem1, bdf_model, out_model, mesh_form, xref, punch, sum_load, size, is_double, cid,
             run_extract_bodies=False, encoding=None, crash_cards=None):
    """
    Reads/writes the BDF

    Parameters
    ----------
    fem1 : BDF()
        The BDF object
    bdf_model : str
        The root path of the bdf filename
    out_model : str
        The path to the output bdf
    mesh_form : str {combined, separate}
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
    run_extract_bodies : bool; default=False
        isolate the fem bodies; typically 1 body; code is still buggy
    encoding : str; default=None
        the file encoding
    crash_cards : ???
        ???
    """
    if crash_cards is None:
        crash_cards = []
    assert os.path.exists(bdf_model), print_bad_path(bdf_model)
    try:
        if '.pch' in bdf_model:
            fem1.read_bdf(bdf_model, xref=False, punch=True, encoding=encoding)
        else:
            fem1.read_bdf(bdf_model, xref=False, punch=punch, encoding=encoding)
            for card in crash_cards:
                if card in fem1.card_count:
                    raise DisabledCardError('card=%r has been disabled' % card)
            #fem1.geom_check(geom_check=True, xref=False)
            skin_filename = 'skin_file.bdf'
            fem1.write_skin_solid_faces(skin_filename, size=16, is_double=False)
            if os.path.exists(skin_filename):
                model = read_bdf(skin_filename, log=fem1.log)
                os.remove(skin_filename)
            if xref:
                if run_extract_bodies:
                    extract_bodies(fem1)

                #fem1.uncross_reference()
                fem1.cross_reference()
                #fem1.safe_cross_reference()
                fem1._xref = True
                spike_fem = read_bdf(fem1.bdf_filename, encoding=encoding,
                                     debug=fem1.debug, log=fem1.log)

                remake = False
                if remake:
                    #log = fem1.log
                    fem1.save('model.obj')
                    fem1.save('model.obj', unxref=False)
                    fem1.write_bdf('spike_out.bdf')
                    fem1.get_bdf_stats()

                    fem1 = BDF(debug=fem1.debug, log=fem1.log)
                    fem1.load('model.obj')
                    fem1.write_bdf('spike_in.bdf')
                    #fem1.log = log
                    fem1.get_bdf_stats()

                    fem1.cross_reference()
                    #fem1.get_bdf_stats()
                    fem1._xref = True

                #fem1.geom_check(geom_check=True, xref=True)
                #fem1.uncross_reference()
                #fem1.cross_reference()
    except:
        print("failed reading %r" % bdf_model)
        raise

    #out_model = bdf_model + '_out'
    #if cid is not None and xref:
        #fem1.resolve_grids(cid=cid)

    if mesh_form is None:
        pass
    elif mesh_form == 'combined':
        fem1.write_bdf(out_model, interspersed=False, size=size, is_double=is_double)
    elif mesh_form == 'separate':
        fem1.write_bdf(out_model, interspersed=False, size=size, is_double=is_double)
    else:
        msg = "mesh_form=%r; allowedForms=['combined','separate']" % mesh_form
        raise NotImplementedError(msg)
    #fem1.write_as_ctria3(out_model)

    fem1._get_maps()
    #remove_unused_materials(fem1)
    #remove_unused(fem1)
    units_to = ['m', 'kg', 's']
    units_from = ['m', 'kg', 's']
    #convert(fem1, units_to, units=units_from)
    if xref:
        try:
            if len(fem1.elements) + len(fem1.masses) > 0:
                fem1.get_area_breakdown()
                fem1.get_volume_breakdown()
                fem1.get_mass_breakdown()
            else:
                fem1.log.warning('no elements found')
        except RuntimeError:
            fem1.log.warning('there are elements, but none with mass?')
    return fem1


def run_fem2(bdf_model, out_model, xref, punch,
             sum_load, size, is_double, mesh_form,
             encoding=None, debug=False, quiet=False):
    """
    Reads/writes the BDF to verify nothing has been lost

    Parameters
    ----------
    bdf_model : str
        the filename to run
    out_model
    xref : bool
       xrefs
    punch : bool
       punches
    sum_load : bool
       sums static load
    size : int
        ???
    is_double : bool
        ???
    mesh_form : str {combined, separate}
        'combined' : interspersed=True
        'separate' : interspersed=False
    debug : bool
        debugs
    quiet : bool
        supress prints
    log : logger / None
        ignored
    """
    assert os.path.exists(bdf_model), bdf_model
    assert os.path.exists(out_model), out_model

    fem2 = BDF(debug=debug, log=None)
    if not quiet:
        fem2.log.info('starting fem2')
    sys.stdout.flush()
    try:
        fem2.read_bdf(out_model, xref=xref, punch=punch, encoding=encoding)
    except:
        print("failed reading %r" % out_model)
        raise

    out_model_2 = bdf_model + '_out2'

    if xref and sum_load:
        p0 = np.array([0., 0., 0.])

        subcase_keys = fem2.case_control_deck.get_subcase_list()
        subcases = fem2.subcases

        sol_200_map = fem2.case_control_deck.sol_200_map
        sol_base = fem2.sol
        is_restart = False
        for line in fem2.system_command_lines:
            if line.strip().upper().startswith('RESTART'):
                is_restart = True
        if not is_restart:
            validate_case_control(fem2, p0, sol_base, subcase_keys, subcases, sol_200_map)

    if mesh_form is not None:
        fem2.write_bdf(out_model_2, interspersed=False, size=size, is_double=is_double)
        os.remove(out_model_2)
    #fem2.write_as_ctria3(out_model_2)
    return fem2

def _assert_has_spc(subcase, fem):
    """
    SPCs may be defined on SPC/SPC1 cards or may be defined on
    the GRID PS field
    """
    has_ps = False
    for nid, node in iteritems(fem.nodes):
        if node.ps:
            has_ps = True
            break
    assert subcase.has_parameter('SPC', 'STATSUB') or has_ps, subcase

def validate_case_control(fem2, p0, sol_base, subcase_keys, subcases, sol_200_map):
    for isubcase in subcase_keys[1:]:  # drop isubcase = 0
        subcase = subcases[isubcase]
        str(subcase)
        assert sol_base is not None, sol_base
        #print('case\n%s' % subcase)
        #if sol_base == 200:
            #analysis = subcase.get_parameter('ANALYSIS')[0]
            #sol = sol_200_map[analysis]
            #if sol is None:
                #msg = 'sol=%s analysis=%r' % (sol, analysis)
                #raise NotImplementedError(msg)
        #else:
            #sol = sol_base
        check_case(sol_base, subcase, fem2, p0, isubcase, subcases)

def check_case(sol, subcase, fem2, p0, isubcase, subcases):
    """
    Checks to see if the case has all the required case control fields
    and that they are valid.

    For example, a SOL 145 deck could requires:

    SUBCASE 10
       METHOD = 42

       # one or both
       SPC = 20 (optional)
       GRID-PS (optional)

       MPC = 30 (optional)
       AEROS = 5

       # one or both
       SUPORT1 = 5 # implicit point to SUPORT1
       # implicit call to SUPORT
    """
    if sol == 24:
        if 'SPC' not in subcase:
            _assert_has_spc(subcase, fem2)
        assert True in subcase.has_parameter('LOAD'), subcase
    elif sol == 64:
        #assert 'NLPARM' in subcase, subcase
        #if 'SPC' not in subcase:
            #_assert_has_spc(subcase, fem2)
        assert True in subcase.has_parameter('LOAD'), subcase
    elif sol == 66:
        assert 'NLPARM' in subcase, subcase
        if 'SPC' not in subcase:
            _assert_has_spc(subcase, fem2)
        assert True in subcase.has_parameter('LOAD', 'TEMPERATURE'), subcase
    elif sol == 99:
        assert 'DLOAD' in subcase, subcase
        assert 'LOADSET' in subcase, subcase
        if 'SPC' not in subcase:
            _assert_has_spc(subcase, fem2)
        #assert True in subcase.has_parameter('LOAD', 'TEMPERATURE'), subcase
        assert True in subcase.has_parameter('TSTEP', 'TSTEPNL'), subcase
    elif sol == 101:
        if 'SPC' not in subcase:
            _assert_has_spc(subcase, fem2)
        assert True in subcase.has_parameter('LOAD', 'TEMPERATURE'), subcase
    elif sol == 103:
        assert True in subcase.has_parameter('METHOD', 'RSMETHOD'), subcase
    elif sol == 105: # buckling
        if 'SPC' not in subcase:
            _assert_has_spc(subcase, fem2)
        assert True in subcase.has_parameter('LOAD', 'METHOD'), subcase
        if 0:
            if 'METHOD' not in subcase:
                subcases = fem2.subcases
                subcase_ids = [isubcase for isubcase in subcases if isubcase > 0]
                assert len(subcases) == 2, 'METHOD not in subcase and not 2 subcases\n%s' % subcase
                subcase_id = subcase.subcase_id
                if subcase_id == 1 and 'METHOD' in subcases[2]:
                    pass
                else:
                    msg = 'METHOD not in subcase and not 2 subcases\n%s' % subcase
                    raise RuntimeError(msg)

        #assert True in subcase.has_parameter('LOAD', 'TEMPERATURE(LOAD)'), subcase
    elif sol == 106: # freq
        assert 'NLPARM' in subcase, subcase
        assert 'LOAD' in subcase, subcase
    elif sol == 107: # ???
        if 'SPC' not in subcase:
            _assert_has_spc(subcase, fem2)
        assert 'LOAD' in subcase, subcase
    elif sol == 108: # freq
        assert 'FREQUENCY' in subcase, subcase
    elif sol == 109:  # time
        if not any(subcase.has_parameter('TIME', 'TSTEP', 'TSTEPNL')):
            subcases = fem2.subcases
            subcase_ids = [isubcase for isubcase in subcases if isubcase > 0]
            has_flag = False
            for isubcase, subcasei in iteritems(fem2.subcases):
                if any(subcasei.has_parameter('TIME', 'TSTEP', 'TSTEPNL')):
                    has_flag = True
            if not has_flag:
                msg = 'sol=%r; [TIME, TSTEP, TSTEPNL] not in subcase\n' % fem2.sol
                for isubcase, subcasei in iteritems(fem2.subcases):
                    msg += str(subcasei)
                raise RuntimeError(msg)

    elif sol == 110:  # ???
        if 'SPC' not in subcase:
            _assert_has_spc(subcase, fem2)
        assert subcase.has_parameter('LOAD', 'STATSUB'), 'sol=%s\n%s' % (sol, subcase)
    elif sol == 111:  # modal frequency
        assert subcase.has_parameter('FREQUENCY'), 'sol=%s\n%s' % (sol, subcase)
        assert any(subcase.has_parameter('METHOD', 'RMETHOD')), 'sol=%s\n%s' % (sol, subcase)
    elif sol == 112:  # modal transient
        assert any(subcase.has_parameter('TIME', 'TSTEP', 'TSTEPNL')), 'sol=%s\n%s' % (sol, subcase)
    elif sol == 114:
        assert 'LOAD' in subcase, subcase
        assert 'HARMONICS' in subcase, subcase
        if 'SPC' not in subcase:
            _assert_has_spc(subcase, fem2)
    elif sol == 118:
        assert 'LOAD' in subcase, subcase
        assert 'HARMONICS' in subcase, subcase
        assert 'SDAMPING' in subcase, subcase
        assert 'FREQUENCY' in subcase, subcase
        assert 'DLOAD' in subcase, subcase
        if 'SPC' not in subcase:
            _assert_has_spc(subcase, fem2)

    elif sol == 129:  # nonlinear transient
        assert any(subcase.has_parameter('TIME', 'TSTEP', 'TSTEPNL')), 'sol=%s\n%s' % (sol, subcase)
    elif sol == 159:  # thermal transient
        assert any(subcase.has_parameter('TIME', 'TSTEP', 'TSTEPNL')), 'sol=%s\n%s' % (sol, subcase)

    elif sol == 144:
        assert any(subcase.has_parameter('TRIM', 'DIVERG')), subcase
        assert fem2.aeros is not None, 'An AEROS card is required for STATIC AERO - SOL %i' % sol
    elif sol == 145:
        if fem2.aero is None:
            msg = 'An AERO card is required for FLUTTER - SOL %i; %s' % (sol, fem2.aero)
            raise RuntimeError(msg)

        assert 'METHOD'in subcase, subcase  # EIGRL
        assert 'FMETHOD' in subcase, subcase  # FLUTTER
    elif sol == 146:
        assert 'METHOD'in subcase, subcase
        assert any(subcase.has_parameter('FREQUENCY', 'TIME', 'TSTEP', 'TSTEPNL')), 'sol=%s\n%s' % (sol, subcase)
        assert any(subcase.has_parameter('GUST', 'LOAD', 'DLOAD')), subcase
        assert fem2.aero is not None, 'An AERO card is required for GUST - SOL %i' % sol
    elif sol == 153: # heat?
        if 'SPC' not in subcase:
            _assert_has_spc(subcase, fem2)
        assert 'NLPARM' in subcase, subcase
        if 'ANALYSIS' in subcase and subcase.get_parameter('ANALYSIS')[0] == 'HEAT':
            assert 'TEMPERATURE' in subcase, subcase
        else:
            assert any(subcase.has_parameter('LOAD')), 'sol=%s\n%s' % (sol, subcase)

    elif sol == 159: #  nonlinear transient; heat?
        assert 'NLPARM' in subcase, 'sol=%s\n%s' % (sol, subcase)
        #assert any(subcase.has_parameter('TIME', 'TSTEP', 'TSTEPNL')), subcase
        #assert any(subcase.has_parameter('GUST', 'LOAD')), subcase
        if 'ANALYSIS' in subcase and subcase.get_parameter('ANALYSIS')[0] == 'HEAT':
            assert 'TEMPERATURE' in subcase, 'sol=%s\n%s' % (sol, subcase)

    elif sol == 200:
        # local level
        # DESSUB - Set constraints (DCONSTR, DCONADD) applied for subcase
        #          (e.g. STRESS, STRAIN, DAMP)
        #          optional locally
        # DESGLB - Set constraints (DCONSTR, DCONADD) applied globally
        #          (e.g. WEIGHT, VOLUME, WMPID, FRMASS)
        #          optional locally
        # DESOBJ - The objective function (DRESP1, DRESP2, DRESP3)
        #          required globally
        # 1 or more DESSUB/DESGLB are required globally
        # 1 DESOBJ is required
        assert 'ANALYSIS' in subcase, 'sol=%s\n%s' % (sol, subcase)

        analysis, options = subcase.get_parameter('ANALYSIS')
        if analysis != 'STATICS':
            # BUCKLING
            if 'DESOBJ' in subcase:
                value, options = subcase.get_parameter('DESOBJ')
                assert value in fem2.dresps, 'value=%s not in dresps' % value
            else:
                fem2.log.warning('no DESOBJ in this subcase; is this a buckling preload case?')
                fem2.log.warning('\n%s' % subcase)

            if 'DESSUB' not in subcase and 'DESGLB' not in subcase:
                fem2.log.warning('no DESSUB/DESGLB in this subcase;'
                                 ' is this a buckling preload case?')
                fem2.log.warning('\n%s' % subcase)

            #assert 'DESSUB' in subcase or 'DESGLB' in subcase, subcase
        if 'DESSUB' in subcase:
            value, options = subcase.get_parameter('DESSUB')
            if value not in fem2.dconstrs:
                msg = 'value=%s not in dconstrs; Allowed DCONSTRs=%s' % (
                    value, np.unique(list(fem2.dconstrs.keys())))
                raise RuntimeError(msg)

        if analysis == 'STATICS':
            sol = 101
            check_case(sol, subcase, fem2, p0, isubcase, subcases)
        elif analysis in ['MODE', 'MODES']:
            sol = 103
            check_case(sol, subcase, fem2, p0, isubcase, subcases)
        elif analysis in ['BUCK', 'BUCKLING']:
            sol = 105
            check_case(sol, subcase, fem2, p0, isubcase, subcases)
        elif analysis == 'DFREQ':
            sol = 108
            check_case(sol, subcase, fem2, p0, isubcase, subcases)
        elif analysis == 'MFREQ':
            sol = 111
            check_case(sol, subcase, fem2, p0, isubcase, subcases)
        elif analysis == 'MTRAN':
            sol = 112
            check_case(sol, subcase, fem2, p0, isubcase, subcases)
        elif analysis in ['SAERO', 'DIVERG', 'DIVERGE']:
            sol = 144
            check_case(sol, subcase, fem2, p0, isubcase, subcases)
        elif analysis == 'FLUTTER':
            sol = 145
            check_case(sol, subcase, fem2, p0, isubcase, subcases)
        elif analysis == 'DCEIG': # direct complex eigenvalues
            sol = 107
            check_case(sol, subcase, fem2, p0, isubcase, subcases)
        #elif analysis == 'MCEIG': # modal direct complex eigenvalues
        elif analysis == 'HEAT': # heat transfer analysis
            sol = 159
            check_case(sol, subcase, fem2, p0, isubcase, subcases)
        else:
            msg = 'analysis = %s\nsubcase =\n%s' % (analysis, subcase)
            raise NotImplementedError(msg)
    elif sol in [114, 115, 116, 118]:
        # cyclic statics, modes, buckling, frequency
        pass
    elif sol in [1, 5, 21, 61, 68, 76, 88, 100, 128, 187, 190, 400, 401, 601, 700, 701]:
        pass
    else:
        msg = 'SOL = %s\n' % (sol)
        msg += str(subcase)
        raise NotImplementedError(msg)

    if any(subcase.has_parameter('TIME', 'TSTEP')):
        if 'TIME' in subcase:
            value, options = subcase.get_parameter('TIME')
        elif 'TSTEP' in subcase:
            value, options = subcase.get_parameter('TSTEP')
        else:
            raise NotImplementedError(subcase)
        assert value in fem2.tsteps, fem2.tsteps

    if 'TSTEPNL' in subcase:
        value, options = subcase.get_parameter('TSTEPNL')
        assert value in fem2.tstepnls, fem2.tstepnls

    if 'SUPORT1' in subcase:
        value, options = subcase.get_parameter('SUPORT1')
        assert value in fem2.suport1, fem2.suport1

    if 'TRIM' in subcase:
        trim_id = subcase.get_parameter('TRIM')[0]
        assert trim_id in fem2.trims, fem2.trims
        trim = fem2.trims[trim_id]

        suport1 = None
        if 'SUPORT1' in subcase:
            suport_id = subcase.get_parameter('SUPORT1')[0]
            suport1 = fem2.suport1[suport_id]
        trim._verify(fem2.suport, suport1, fem2.aestats, fem2.aeparams,
                     fem2.aelinks, fem2.aesurf, xref=True)
        assert 'DIVERG' not in subcase, subcase

    if 'DIVERG' in subcase:
        value, options = subcase.get_parameter('DIVERG')
        assert value in fem2.divergs, fem2.divergs
        assert 'TRIM' not in subcase, subcase

    if 'METHOD' in subcase:
        method_id = subcase.get_parameter('METHOD')[0]
        if method_id in fem2.methods:
            method = fem2.methods[method_id]
        #elif method_id in fem2.cMethods:
            #method = fem2.cMethods[method_id]
        else:
            method_ids = list(fem2.methods.keys())
            raise RuntimeError('METHOD = %s not in method_ids=%s' % (method_id, method_ids))

        assert sol in [5, 76, 101, 103, 105, 106, 107, 108, 110, 111,
                       112, 144, 145, 146, 187], 'sol=%s METHOD\n%s' % (sol, subcase)

    if 'CMETHOD' in subcase:
        cmethod_id = subcase.get_parameter('CMETHOD')[0]
        if cmethod_id in fem2.cMethods:
            method = fem2.cMethods[cmethod_id]
        #elif method_id in fem2.cMethods:
            #method = fem2.cMethods[method_id]
        else:
            cmethod_ids = list(fem2.cMethods.keys())
            raise RuntimeError('CMETHOD = %s not in cmethod_ids=%s' % (cmethod_id, cmethod_ids))
        assert sol in [110, 111, 145], 'sol=%s CMETHOD\n%s' % (sol, subcase)

    if 'RMETHOD' in subcase:
        rmethod_id = subcase.get_parameter('RMETHOD')[0]
        #if method_id in fem2.methods:
            #method = fem2.methods[method_id]
        #elif method_id in fem2.cMethods:
            #method = fem2.cMethods[method_id]
        #else:
            #method_ids = list(fem2.methods.keys())
            #raise RuntimeError('METHOD = %s not in method_ids=%s' % (method_id, method_ids))

        assert sol in [110, 111], 'sol=%s RMETHOD\n%s' % (sol, subcase)

    if 'FMETHOD' in subcase:
        method_id = subcase.get_parameter('FMETHOD')[0]
        method = fem2.flutters[method_id]
        assert sol in [145], 'sol=%s FMETHOD\n%s' % (sol, subcase)
    if 'LOAD' in subcase:
        loadcase_id = subcase.get_parameter('LOAD')[0]
        force, moment = fem2.sum_forces_moments(p0, loadcase_id, include_grav=False)
        eids = None
        nids = None
        force2, moment2 = fem2.sum_forces_moments_elements(
            p0, loadcase_id, eids, nids, include_grav=False)
        assert np.allclose(force, force2), 'force=%s force2=%s' % (force, force2)
        assert np.allclose(moment, moment2), 'moment=%s moment2=%s' % (moment, moment2)
        print('  isubcase=%i F=%s M=%s' % (isubcase, force, moment))
        assert sol in [1, 5, 24, 61, 64, 66, 101, 103, 105, 106, 107,
                       108, 109, 110, 111, 112, 144, 145, 153, 400, 601,
                      ], 'sol=%s LOAD\n%s' % (sol, subcase)
    else:
        # print('is_load =', subcase.has_parameter('LOAD'))
        pass

    if 'FREQUENCY' in subcase:
        freq_id = subcase.get_parameter('FREQUENCY')[0]
        freq = fem2.frequencies[freq_id]
        assert sol in [26, 68, 76, 78, 88, 108, 101, 111, 112, 118, 146], 'sol=%s FREQUENCY' % sol

    # if 'LSEQ' in subcase:
        # lseq_id = subcase.get_parameter('LSEQ')[0]
        # lseq = fem2.loads[lseq_id]
        # assert sol in [], sol
        # print(lseq)
    if 'SPC' in subcase:
        spc_id = subcase.get_parameter('SPC')[0]
        fem2.get_spcs(spc_id)
    if 'MPC' in subcase:
        mpc_id = subcase.get_parameter('MPC')[0]
        fem2.get_mpcs(mpc_id)

    if 'SDAMPING' in subcase:
        sdamping_id = subcase.get_parameter('SDAMPING')[0]
        sdamping_table = fem2.tables_sdamping[sdamping_id]

    if 'LOADSET' in subcase:
        loadset_id = subcase.get_parameter('LOADSET')[0]
        lseq = fem2.loads[loadset_id]

    if 'DLOAD' in subcase:
        assert sol in [26, 68, 76, 78, 88, 99, 103, 108, 109, 111, 112, 118, 129, 146,
                       153, 159, 400, 401, 601], 'sol=%s DLOAD\n%s' % (sol, subcase)
        #if 'LOADSET' in subcase:
            #raise NotImplementedError('LOADSET & DLOAD -> LSEQ')
        if 'IC' in subcase:
            ic_val = subcase.get_parameter('IC')[0]
            if ic_val != 0:
                # TIC - SOL 109, 129, 601, 701 (structural/physical)
                #
                # TEMP/TEMPD - SOL 159
                #
                # TIC - SET for modal
                #
                # STATSUB - IC=subcase_id (SOL 109, 112)
                # For the PHYSICAL option, n is the set identification number
                # of TIC bulk entries for structural analysis (SOLs 109, 129,
                # 601, and 701) or TEMP and TEMPD bulk entries for heat
                # transfer analysis (SOL 159). For the MODAL option, n is
                # the set identification number of TIC bulk entries for modal
                # transient analysis (SOL 112). For the STATSUB option, n is
                # the ID of a static analysis subcase (SOL 109 and 112). For
                # the TZERO option, n is not used by the software, although a
                # dummy value must be defined. (Integer>0)
                if sol in [109, 129, 601, 701]:
                    raise NotImplementedError('IC & DLOAD; sol=%s -> TIC\n%s' % (sol, subcase))
                elif sol == 159:
                    fem2.loads[ic_val]
                else:
                    raise NotImplementedError('IC & DLOAD; sol=%s -> TIC\n%s' % (sol, subcase))

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

        scale_factors2 = []
        loads2 = []
        for load in dload:
            if isinstance(load, DLOAD):
                scale = load.scale
                scale_factors = []
                loads = []
                # scale_factors, loads = load.get_reduced_loads()
                for load, scale_factor in zip(load.load_ids, load.scale_factors):
                    if isinstance(load, list):
                        for loadi in load:
                            assert not isinstance(loadi, list), loadi
                            scale_factors.append(scale * scale_factor)
                            loads.append(loadi)
                    else:
                        scale_factors.append(scale * scale_factor)
                        assert not isinstance(load, list), load
                        loads.append(load)
                scale_factors2 += scale_factors
                loads2 += loads
            else:
                scale_factors2.append(1.)
                loads2.append(load)

        if sol in [108, 111]:  # direct frequency, modal frequency
            for load2, scale_factor in zip(loads2, scale_factors2):
                freq_id = subcase.get_parameter('FREQ')[0]
                freqs = fem2.frequencies[freq_id]
                for freq in freqs:
                    if freq.type in ['FREQ', 'FREQ1', 'FREQ2']:
                        fmax = freq.freqs[-1]
                        force = load2.get_load_at_freq(fmax) * scale_factor
        elif sol in [109, 129]:  # direct transient (time linear), time nonlinear
            for load2, scale_factor in zip(loads2, scale_factors2):
                force = load2.get_load_at_time(0.) * scale_factor
        else:
            fem2.log.debug('solution=%s; DLOAD is not supported' % sol)

        # print(loads)

def divide(value1, value2):
    """
    Used to divide the number of cards to check that nothing was lost.
    Handles division by 0 by returning 0, which is the reciprocal.
    """
    if value1 == value2:  # good for 0/0
        return 1.0
    else:
        try:
            div_value = value1 / float(value2)
        except ZeroDivisionError:
            div_value = 0.
    return div_value


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
    card_dict = model.get_cards_by_card_types(card_types,
                                              reset_type_to_slot_map=False)
    for card_type, cards in iteritems(card_dict):
        for card in cards:
            msg = 'this should never crash here...card_type=%s card.type=%s' % (
                card_type, card.type)
            if card_type != card.type:
                raise RuntimeError(msg)


def compare_card_count(fem1, fem2, print_stats=False, quiet=False):
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
    return compute_ints(cards1, cards2, fem1, quiet=quiet)


def compute_ints(cards1, cards2, fem1, quiet=True):
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
        if not quiet or not star or factor1 != factor2:
            if factor1 != factor2:
                factor_msg = 'diff=%s factor1=%g factor2=%g' % (
                    diff, factor1, factor2)
            msg += '  %skey=%-7s value1=%-7s value2=%-7s' % (
                star, key, value1, value2) + factor_msg
        if msg:
            msg = msg.rstrip()
            print(msg)
    #return list_keys1 + list_keys2
    return diff_keys1 + diff_keys2


def compute(cards1, cards2, quiet=False):
    """
    Computes the difference between two dictionaries to data is the same
    """
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
            if not quiet:
                msg += '    key=%-7s value1=%-7s value2=%-7s' % (
                    key, value1, value2)
        else:
            msg += '   *key=%-7s value1=%-7s value2=%-7s' % (
                key, value1, value2)
        msg = msg.rstrip()
        if msg:
            print(msg)


def get_element_stats(fem1, fem2, quiet=False):
    """verifies that the various element methods work"""
    for (key, loads) in sorted(iteritems(fem1.loads)):
        for load in loads:
            try:
                all_loads = load.get_loads()
                if not isinstance(all_loads, list):
                    raise TypeError('allLoads should return a list...%s'
                                    % (type(all_loads)))
            except:
                raise
                #print("load statistics not available - load.type=%s "
                      #"load.sid=%s" % (load.type, load.sid))

    fem1._verify_bdf()

    if fem1.elements:
        fem1.get_elements_nodes_by_property_type()
    mass, cg, I = fem1.mass_properties(reference_point=None, sym_axis=None)
    mass, cg, I = fem1._mass_properties_new(reference_point=None, sym_axis=None)
    if not quiet:
        print("mass = %s" % mass)
        print("cg   = %s" % cg)
        print("Ixx=%s, Iyy=%s, Izz=%s \nIxy=%s, Ixz=%s, Iyz=%s" % tuple(I))
        #mass, cg, I = fem1._mass_properties_new(reference_point=None, sym_axis=None)
        #print("mass_old =", mass)
        #print("cg_old   =", cg)
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
                      "dmig.type=%s matrix.name=%s" % (dmig.type, dmig.name))
        except:
            print("*stats - dmig.type=%s name=%s  matrix=\n%s"
                  % (dmig.type, dmig.name, str(dmig)))
            raise

    for (key, dmi) in sorted(iteritems(fem1.dmis)):
        try:
            if isinstance(dmi, NastranMatrix):
                dmi.get_matrix()
            else:
                print("statistics not available - "
                      "dmi.type=%s matrix.name=%s" % (dmi.type, dmi.name))
        except:
            print("*stats - dmi.type=%s name=%s  matrix=\n%s"
                  % (dmi.type, dmi.name, str(dmi)))
            raise

    for (key, dmij) in sorted(iteritems(fem1.dmijs)):
        try:
            if isinstance(dmij, NastranMatrix):
                dmij.get_matrix()
            else:
                print("statistics not available - "
                      "dmij.type=%s matrix.name=%s" % (dmij.type, dmij.name))
        except:
            print("*stats - dmij.type=%s name=%s  matrix=\n%s"
                  % (dmij.type, dmij.name, str(dmij)))
            raise

    for (key, dmiji) in sorted(iteritems(fem1.dmijis)):
        try:
            if isinstance(dmiji, NastranMatrix):
                dmiji.get_matrix()
            else:
                print("statistics not available - "
                      "dmiji.type=%s matrix.name=%s" % (dmiji.type, dmiji.name))
        except:
            print("*stats - dmiji.type=%s name=%s  matrix=\n%s"
                  % (dmiji.type, dmiji.name, str(dmiji)))
            raise

    for (key, dmik) in sorted(iteritems(fem1.dmiks)):
        try:
            if isinstance(dmik, NastranMatrix):
                dmik.get_matrix()
            else:
                print("statistics not available - "
                      "dmik.type=%s matrix.name=%s" % (dmik.type, dmik.name))
        except:
            print("*stats - dmik.type=%s name=%s  matrix=\n%s"
                  % (dmik.type, dmik.name, str(dmik)))
            raise

def compare(fem1, fem2, xref=True, check=True, print_stats=True, quiet=False):
    """compares two fem objects"""
    diff_cards = compare_card_count(fem1, fem2, print_stats=print_stats, quiet=quiet)
    if xref and check:
        get_element_stats(fem1, fem2, quiet=quiet)
        get_matrix_stats(fem1, fem2)
    compare_card_content(fem1, fem2)
    #compare_params(fem1, fem2)
    #print_points(fem1, fem2)
    return diff_cards


#def compare_params(fem1, fem2):
    #raise RuntimeError('is compare_parms used?')
    #compute(fem1.params, fem2.params)


#def print_points(fem1, fem2):
    #raise RuntimeError('is print_points used?')
    #for nid, node in sorted(iteritems(fem1.nodes)):
        #print("%s   xyz=%s  n1=%s  n2=%s" % (nid, node.xyz, node.get_position(True),
                                             #fem2.Node(nid).get_position()))
        #break
    #coord = fem1.Coord(5)
    #print(coord)
    #print coord.Stats()


def main():
    """
    The main function for the command line ``test_bdf`` script.
    """
    encoding = sys.getdefaultencoding()
    from pyNastran.utils.docopt_types import docopt_types
    msg = "Usage:\n"
    msg += "  test_bdf [-q] [-D] [-i] [-e E] [--crash C] [-x] [-p] [-c] [-L]      [-f] [--encoding ENCODE] BDF_FILENAME\n"
    msg += "  test_bdf [-q] [-D] [-i] [-e E] [--crash C] [-x] [-p] [-c] [-L] [-d] [-f] [--encoding ENCODE] BDF_FILENAME\n"
    msg += "  test_bdf [-q] [-D] [-i] [-e E] [--crash C] [-x] [-p] [-c] [-L] [-l] [-f] [--encoding ENCODE] BDF_FILENAME\n"
    msg += "  test_bdf [-q] [-D] [-i] [-e E] [--crash C]      [-p]                [-f] [--encoding ENCODE] BDF_FILENAME\n"
    msg += "  test_bdf [-q] [-D] [-i] [-e E] [--crash C] [-x] [-p] [-s]           [-f] [--encoding ENCODE] BDF_FILENAME\n"

    #msg += "  test_bdf [-q] [-p] [-o [<VAR=VAL>]...] BDF_FILENAME\n" #
    msg += '  test_bdf -h | --help\n'
    msg += '  test_bdf -v | --version\n'
    msg += '\n'

    msg += "Positional Arguments:\n"
    msg += "  BDF_FILENAME   path to BDF/DAT/NAS file\n"
    msg += '\n'

    msg += 'Options:\n'
    msg += '  -x, --xref     disables cross-referencing and checks of the BDF\n'
    msg += '                 (default=True -> on)\n'
    msg += '  -p, --punch    disables reading the executive and case control decks in the BDF\n'
    msg += '                 (default=False -> reads entire deck)\n'
    msg += '  -c, --check    disables BDF checks.  Checks run the methods on \n'
    msg += '                 every element/property to test them.  May fails if a \n'
    msg += '                 card is fully not supported (default=False)\n'
    msg += '  -l, --large    writes the BDF in large field, single precision format (default=False)\n'
    msg += '  -d, --double   writes the BDF in large field, double precision format (default=False)\n'
    msg += '  -L, --loads    Disables forces/moments summation for the different subcases (default=True)\n'
    msg += '  -e E, --nerrors E  Allow for cross-reference errors (default=100)\n'
    msg += '  --encoding ENCODE  the encoding method (default=None -> %r)\n' % encoding
    msg += '  -q, --quiet        prints debug messages (default=False)\n'

    msg += "\n"
    msg += "Developer:\n"
    msg += '  --crash C,       Crash on specific cards (e.g. CGEN,EGRID)\n'
    msg += '  -D, --dumplines  Writes the BDF exactly as read with the INCLUDES processed\n'
    msg += '                   (pyNastran_dump.bdf)\n'
    msg += '  -i, --dictsort   Writes the BDF with exactly as read with the INCLUDES processed\n'
    msg += '                   (pyNastran_dict.bdf)\n'
    msg += '  -f, --profile    Profiles the code (default=False)\n'
    msg += '  -s, --stop       Stop after first read/write (default=False)\n'
    msg += "\n"
    msg += "Info:\n"
    msg += '  -h, --help     show this help message and exit\n'
    msg += "  -v, --version  show program's version number and exit\n"

    if len(sys.argv) == 1:
        sys.exit(msg)

    ver = str(pyNastran.__version__)
    type_defaults = {
        '--nerrors' : [int, 100],
    }
    data = docopt_types(msg, version=ver, type_defaults=type_defaults)

    data['--xref'] = not data['--xref']
    data['--loads'] = not data['--loads']
    if not data['--encoding']:
        data['--encoding'] = None

    for key, value in sorted(iteritems(data)):
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

    crash_cards = []
    if data['--crash']:
        crash_cards = data['--crash'].split(',')


    #print(data)
    debug = True
    if data['--quiet']:
        debug = None
    if data['--profile']:
        #import cProfile
        import pstats

        import cProfile
        prof = cProfile.Profile()
        prof.runcall(
            run_bdf,
            '.',
            data['BDF_FILENAME'],
            debug=debug,
            xref=['--xref'],
            # xref_safe=data['--xref_safe'],
            check=not(data['--check']),
            punch=data['--punch'],
            size=size,
            is_double=is_double,
            sum_load=data['--loads'],
            stop=data['--stop'],
            quiet=data['--quiet'],
            dumplines=data['--dumplines'],
            dictsort=data['--dictsort'],
            nerrors=data['--nerrors'],
            encoding=data['--encoding'],
            crash_cards=crash_cards,
            run_extract_bodies=False,
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
            debug=debug,
            xref=data['--xref'],
            # xref_safe=data['--xref_safe'],
            check=not(data['--check']),
            punch=data['--punch'],
            size=size,
            is_double=is_double,
            sum_load=data['--loads'],
            stop=data['--stop'],
            quiet=data['--quiet'],
            dumplines=data['--dumplines'],
            dictsort=data['--dictsort'],
            nerrors=data['--nerrors'],
            encoding=data['--encoding'],
            crash_cards=crash_cards,
            run_extract_bodies=False,
        )
    print("total time:  %.2f sec" % (time.time() - time0))


if __name__ == '__main__':  # pragma: no cover
    main()
