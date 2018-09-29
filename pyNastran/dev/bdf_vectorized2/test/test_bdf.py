"""
``test_bdfv`` runs multiple checks on a BDF in order to make sure that:
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
from itertools import chain
from typing import List, Tuple, Optional
import numpy as np
warnings.simplefilter('always')

np.seterr(all='raise')

from pyNastran.utils.numpy_utils import integer_types
from pyNastran.utils import check_path
from pyNastran.bdf.errors import (
    #CrossReferenceError,
    CardParseSyntaxError, DuplicateIDsError, MissingDeckSections)
from pyNastran.dev.bdf_vectorized2.bdf_vectorized import BDF, read_bdf
from pyNastran.bdf.test.test_bdf import divide, get_matrix_stats, compare_card_content

import pyNastran.bdf.test
TEST_PATH = pyNastran.bdf.test.__path__[0]

class DisabledCardError(RuntimeError):
    """lets bdf_test.py flag cards as auto-crashing and then skipping the deck (e.g., CGEN)"""
    pass

def run_all_files_in_folder(folder, debug=False, xref=True, check=True,
                            punch=False, cid=None, nastran=''):
    # type: (str, bool, bool, bool, bool, Optional[int], str) -> None
    """runs all the BDFs in a given folder"""
    print("folder = %s" % folder)
    filenames = os.listdir(folder)
    run_lots_of_files(filenames, debug=debug, xref=xref, check=check,
                      punch=punch, cid=cid, nastran=nastran)


def run_lots_of_files(filenames, folder='', debug=False, xref=True, check=True,
                      punch=False, cid=None, nastran='', encoding=None,
                      size=None, is_double=None, post=None, sum_load=True, dev=True,
                      crash_cards=None, pickle_obj=True):
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
    pickle_obj : bool; default=True
        tests pickling

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
                                                  run_extract_bodies=False, pickle_obj=pickle_obj)
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
        import psutil  # type: ignore

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
            nerrors=0, dev=False, crash_cards=None, safe_xref=True, pickle_obj=False, safe=False,
            stop_on_failure=True):
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
    pickle_obj : bool; default=True
        tests pickling
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
        bdf_model = os.path.join(TEST_PATH, folder, bdf_filename)

    model, ext = os.path.splitext(bdf_model)
    out_model = '%s.test_bdfv%s' % (model, ext)

    fem1, fem2, diff_cards = run_and_compare_fems(
        bdf_model, out_model, debug=debug, xref=xref, check=check,
        punch=punch, cid=cid, mesh_form=mesh_form,
        print_stats=print_stats, encoding=encoding,
        sum_load=sum_load, size=size, is_double=is_double,
        stop=stop, nastran=nastran, post=post,
        dynamic_vars=dynamic_vars,
        quiet=quiet, dumplines=dumplines, dictsort=dictsort,
        nerrors=nerrors, dev=dev, crash_cards=crash_cards,
        safe_xref=safe_xref,
        run_extract_bodies=run_extract_bodies, pickle_obj=pickle_obj,
        stop_on_failure=stop_on_failure,
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
        safe_xref=True, run_extract_bodies=False, pickle_obj=False,
        stop_on_failure=True,
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
                        encoding=encoding, crash_cards=crash_cards, safe_xref=safe_xref,
                        pickle_obj=pickle_obj, stop=stop)
        if stop:
            if not quiet:
                print('card_count:')
                print('-----------')
                for card_name, card_count in sorted(fem1.card_count.items()):
                    print('key=%-8s value=%s' % (card_name, card_count))
            return fem1, None, None
        fem2 = run_fem2(bdf_model, out_model, xref, punch, sum_load, size, is_double, mesh_form,
                        encoding=encoding, debug=debug, quiet=quiet,
                        stop_on_failure=stop_on_failure)

        diff_cards = compare(fem1, fem2, xref=xref, check=check,
                             print_stats=print_stats, quiet=quiet)
        test_get_cards_by_card_types(fem2)

        #fem2.update_model_by_desvars(xref)
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
        from pyNastran.op2.op2 import read_op2
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
        bdf = read_bdf(bdf_model, debug=False)
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
        if not os.path.exists(op2_model):
            raise RuntimeError('%s failed' % op2_model)
        op2 = read_op2(op2_model2)
        print(op2.get_op2_stats())

def run_fem1(fem1, bdf_model, out_model, mesh_form, xref, punch, sum_load, size, is_double, cid,
             run_extract_bodies=False, encoding=None, crash_cards=None, safe_xref=True,
             pickle_obj=False, stop=False):
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
    safe_xref : bool; default=False
        ???
    run_extract_bodies : bool; default=False
        isolate the fem bodies; typically 1 body; code is still buggy
    encoding : str; default=None
        the file encoding
    crash_cards : ???
        ???
    """
    if crash_cards is None:
        crash_cards = []
    check_path(bdf_model, 'bdf_model')
    try:
        if '.pch' in bdf_model:
            fem1.read_bdf(bdf_model, xref=False, punch=True, encoding=encoding)
        else:
            fem1.read_bdf(bdf_model, xref=False, punch=punch, encoding=encoding)
            for card in crash_cards:
                if card in fem1.card_count:
                    raise DisabledCardError('card=%r has been disabled' % card)
            #fem1.geom_check(geom_check=True, xref=False)
            if not stop and not xref:
                skin_filename = 'skin_file.bdf'
                fem1.write_skin_solid_faces(skin_filename, size=16, is_double=False)
                if os.path.exists(skin_filename):
                    read_bdf(skin_filename, log=fem1.log)
                    os.remove(skin_filename)
            if xref:
                if run_extract_bodies:
                    extract_bodies(fem1)

                # 1. testing that these methods word without xref
                #fem1._get_rigid()
                #fem1.get_dependent_nid_to_components()
                #fem1._get_maps(eids=None, map_names=None,
                               #consider_0d=True, consider_0d_rigid=True,
                               #consider_1d=True, consider_2d=True, consider_3d=True)
                #fem1.get_dependent_nid_to_components()

                # 1. testing that these methods work with xref
                fem1._get_rigid()
                #common_node_ids = list(fem1.nodes.keys())
                #fem1.get_rigid_elements_with_node_ids(common_node_ids)

                #for spc_id in set(list(fem1.spcadds.keys()) + list(fem1.spcs.keys())):
                    #fem1.get_reduced_spcs(spc_id)
                #for mpc_id in set(list(fem1.mpcadds.keys()) + list(fem1.mpcs.keys())):
                    #fem1.get_reduced_mpcs(mpc_id)

                #fem1.get_dependent_nid_to_components()
                #fem1._get_maps(eids=None, map_names=None,
                               #consider_0d=True, consider_0d_rigid=True,
                               #consider_1d=True, consider_2d=True, consider_3d=True)
                #fem1.get_dependent_nid_to_components()
                #fem1.get_pid_to_node_ids_and_elements_array(pids=None, etypes=None, idtype='int32',
                                                            #msg=' which is required by test_bdf')
                #fem1.get_property_id_to_element_ids_map(msg=' which is required by test_bdf')
                #fem1.get_material_id_to_property_ids_map(msg=' which is required by test_bdf')
                #fem1.get_element_ids_list_with_pids(pids=None)
                #fem1.get_element_ids_dict_with_pids(pids=None, stop_if_no_eids=False,
                                                    #msg=' which is required by test_bdf')
                #fem1.get_node_id_to_element_ids_map()
                #fem1.get_node_id_to_elements_map()


                read_bdf(fem1.bdf_filename, encoding=encoding,
                         debug=fem1.debug, log=fem1.log)

                fem1 = remake_model(bdf_model, fem1, pickle_obj)
                #fem1.geom_check(geom_check=True, xref=True)
    except:
        print("failed reading %r" % bdf_model)
        raise

    #out_model = bdf_model + '_out'
    #if cid is not None and xref:
        #fem1.resolve_grids(cid=cid)

    if mesh_form is None:
        pass
    elif mesh_form == 'combined':
        fem1.write_bdf(out_model, interspersed=True, size=size, is_double=is_double)
    elif mesh_form == 'separate':
        fem1.write_bdf(out_model, interspersed=False, size=size, is_double=is_double)
    else:
        msg = "mesh_form=%r; allowed_mesh_forms=['combined','separate']" % mesh_form
        raise NotImplementedError(msg)
    #fem1.write_as_ctria3(out_model)

    fem1._get_maps()
    #remove_unused_materials(fem1)
    #remove_unused(fem1)
    #units_to = ['m', 'kg', 's']
    #units_from = ['m', 'kg', 's']
    #convert(fem1, units_to, units=units_from)
    if xref:
        check_for_cd_frame(fem1)

        #try:
            #fem1.get_area_breakdown()
            #fem1.get_volume_breakdown()
        #except:
            #if len(fem1.masses) > 0:
                #fem1.log.warning('no elements with area/volume found, but elements with mass were')
            #else:
                #fem1.log.warning('no elements found')

        #if len(fem1.elements) + len(fem1.masses) > 0:
            #try:
                #fem1.get_mass_breakdown()
            #except RuntimeError:
                #fem1.log.warning('no elements with mass found')
    return fem1


def remake_model(bdf_model, fem1, pickle_obj):
    """reloads the model if we're testing pickling"""
    remake = pickle_obj
    if remake:
        #log = fem1.log
        model_name = os.path.splitext(bdf_model)[0]
        obj_model = '%s.test_bdfv.obj' % (model_name)
        #out_model_8 = '%s.test_bdfv.bdf' % (model_name)
        #out_model_16 = '%s.test_bdfv.bdf' % (model_name)

        fem1.save(obj_model)
        fem1.save(obj_model, unxref=False)
        #fem1.write_bdf(out_model_8)
        fem1.get_bdf_stats()

        fem1 = BDF(debug=fem1.debug, log=fem1.log)
        fem1.load(obj_model)
        #fem1.write_bdf(out_model_8)
        #fem1.log = log
        os.remove(obj_model)
        fem1.get_bdf_stats()
    return fem1

def check_for_cd_frame(fem1):
    """
    A cylindrical/spherical CD frame will cause problems with the
    grid point force transformation
    """
    if any([card_name in fem1.card_count for card_name in ['GRID', 'SPOINT', 'EPOINT', 'RINGAX']]):
        icd_transform, icp_transform, xyz_cp, nid_cp_cd = fem1.get_displacement_index_xyz_cp_cd(
            fdtype='float64', idtype='int32', sort_ids=True)
        cds = np.unique(nid_cp_cd[:, 2])
        cd_coords = []
        for cd in cds:
            if cd == -1:
                continue
            coord = fem1.coords[cd]
            # coordRs work in op2 extraction
            if coord.type not in ['CORD2R', 'CORD1R']:
                cd_coords.append(cd)
        if cd_coords:
            msg = (
                'GRID-CD coords=%s can cause a problem in the OP2 results processing; '
                'be careful' % cd_coords
            )
            fem1.log.warning(msg)

def run_fem2(bdf_model, out_model, xref, punch,
             sum_load, size, is_double, mesh_form,
             encoding=None, debug=False, quiet=False,
             stop_on_failure=True):
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
    """
    assert os.path.exists(bdf_model), bdf_model
    assert os.path.exists(out_model), out_model

    fem2 = BDF(debug=debug, log=None)
    if not quiet:
        fem2.log.info('starting fem2')
    sys.stdout.flush()
    try:
        fem2.read_bdf(out_model, xref=False, punch=punch, encoding=encoding)
    except:
        print("failed reading %r" % out_model)
        raise

    out_model_2 = bdf_model + '_out2'

    if xref and sum_load:
        if 'POST' in fem2.params:
            value = fem2.params['POST'].values[0]
            if value >= 0:
                msg = 'PARAM,POST,%i is not supported by the OP2 reader' % value
                fem2.log.error(msg)
        else:
            msg = 'PARAM,POST,0 is not supported by the OP2 reader'
            fem2.log.error(msg)

        p0 = np.array([0., 0., 0.])

        subcase_keys = fem2.case_control_deck.get_subcase_list()
        subcases = fem2.subcases

        sol_200_map = fem2.case_control_deck.sol_200_map
        sol_base = fem2.sol
        is_restart = False
        for line in fem2.system_command_lines:
            if line.strip().upper().startswith('RESTART'):
                is_restart = True
        #if not is_restart:
            #validate_case_control(fem2, p0, sol_base, subcase_keys, subcases, sol_200_map,
                                  #stop_on_failure=stop_on_failure)

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
    if 'SPC' not in subcase:
        has_ps = False
        for nid, node in fem.nodes.items():
            if node.ps:
                has_ps = True
                break
        assert subcase.has_parameter('SPC', 'STATSUB') or has_ps, subcase

def require_cards(card_names, log, soltype, sol, subcase):
    nerrors = 0
    for card_name in card_names:
        if card_name not in subcase:
            log.error('A %s card is required for %s - SOL %i\n%s' % (
                card_name, soltype, sol, subcase))
            nerrors += 1
    return nerrors

def test_get_cards_by_card_types(model):
    # type: (BDF) -> None
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
    for card_type, cards in card_dict.items():
        for card in cards:
            msg = 'this should never crash here...card_type=%s card.type=%s' % (
                card_type, card.type)
            if card_type != card.type:
                raise RuntimeError(msg)


def compare_card_count(fem1, fem2, print_stats=False, quiet=False):
    # type: (BDF, BDF, bool, bool) -> List[str]
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
        print(fem1.loads)
    else:
        fem1.get_bdf_stats()
    return compute_ints(cards1, cards2, fem1, quiet=quiet)


def compute_ints(cards1, cards2, fem1, quiet=True):
    """
    computes the difference / ratio / inverse-ratio between
    fem1 and fem2 to verify the number of card are the same:

    Examples
    --------

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


def compare(fem1, fem2, xref=True, check=True, print_stats=True, quiet=False):
    """compares two fem objects"""
    diff_cards = compare_card_count(fem1, fem2, print_stats=print_stats, quiet=quiet)
    if xref and check:
        #get_element_stats(fem1, fem2, quiet=quiet)
        get_matrix_stats(fem1, fem2)
    #compare_card_content(fem1, fem2)
    return diff_cards


def get_test_bdf_data():
    """defines the docopt interface"""
    encoding = sys.getdefaultencoding()

    from pyNastran.utils.docopt_types import docopt_types
    options = '[-e E] [--encoding ENCODE] [-q] [-D] [-i] [--crash C] [-k] [-f] '
    msg = (
        "Usage:\n"
        '  test_bdfv [-x | --safe] [-p] [-c] [-L]      %sBDF_FILENAME\n' % options +
        '  test_bdfv [-x | --safe] [-p] [-c] [-L] [-d] %sBDF_FILENAME\n' % options +
        '  test_bdfv [-x | --safe] [-p] [-c] [-L] [-l] %sBDF_FILENAME\n' % options +
        '  test_bdfv               [-p]                %sBDF_FILENAME\n' % options +
        '  test_bdfv [-x | --safe] [-p] [-s]           %sBDF_FILENAME\n' % options +

        #"  test_bdf [-q] [-p] [-o [<VAR=VAL>]...] BDF_FILENAME\n"
        '  test_bdfv -h | --help\n'
        '  test_bdfv -v | --version\n'
        '\n'

        'Positional Arguments:\n'
        '  BDF_FILENAME   path to BDF/DAT/NAS file\n'
        '\n'

        'Options:\n'
        '  -x, --xref     disables cross-referencing and checks of the BDF\n'
        '                 (default=True -> on)\n'
        '  --safe         Use safe cross-reference (default=False)\n'
        '  -p, --punch    disables reading the executive and case control decks in the BDF\n'
        '                 (default=False -> reads entire deck)\n'
        '  -c, --check    disables BDF checks.  Checks run the methods on \n'
        '                 every element/property to test them.  May fails if a \n'
        '                 card is fully not supported (default=False)\n'
        '  -l, --large    writes the BDF in large field, single precision format (default=False)\n'
        '  -d, --double   writes the BDF in large field, double precision format (default=False)\n'
        '  -L, --loads    Disables forces/moments summation for the different subcases (default=True)\n'
        '  -e E, --nerrors E  Allow for cross-reference errors (default=100)\n'
        '  --encoding ENCODE  the encoding method (default=None -> %r)\n' % encoding +
        '  -q, --quiet        prints debug messages (default=False)\n'

        '\n'
        'Developer:\n'
        '  --crash C,       Crash on specific cards (e.g. CGEN,EGRID)\n'
        '  -D, --dumplines  Writes the BDF exactly as read with the INCLUDES processed\n'
        '                   (pyNastran_dump.bdf)\n'
        '  -i, --dictsort   Writes the BDF with exactly as read with the INCLUDES processed\n'
        '                   (pyNastran_dict.bdf)\n'
        '  -f, --profile    Profiles the code (default=False)\n'
        '  -s, --stop       Stop after first read/write (default=False)\n'
        '  -k, --pickle     Pickles the data objects (default=False)\n'
        '\n'
        'Info:\n'
        '  -h, --help     show this help message and exit\n'
        "  -v, --version  show program's version number and exit\n"
    )
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

    return data

def main():
    """
    The main function for the command line ``test_bdfv`` script.
    """
    data = get_test_bdf_data()
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

    crash_cards = []
    if data['--crash']:
        crash_cards = data['--crash'].split(',')

    #print(data)
    debug = True
    if data['--quiet']:
        debug = None
    if data['--profile']:
        import pstats

        import cProfile
        prof = cProfile.Profile()
        prof.runcall(
            run_bdf,
            '.',
            data['BDF_FILENAME'],
            debug=debug,
            xref=['--xref'],
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
            pickle_obj=data['--pickle'],
            safe_xref=data['--safe'],
            print_stats=True,
            stop_on_failure=False,
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
            pickle_obj=data['--pickle'],
            safe_xref=data['--safe'],
            print_stats=True,
            stop_on_failure=False,
        )
    print("total time:  %.2f sec" % (time.time() - time0))


if __name__ == '__main__':  # pragma: no cover
    main()
