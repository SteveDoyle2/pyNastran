"""
``test_bdf`` runs multiple checks on a BDF in order to make sure that:
  - no data is lost on IO
  - card field types are correct (e.g. node_ids are integers)
  - various card methods (e.g. Area) work correctly

As such, ``test_bdf`` is very useful for debugging models.

"""
import os
import sys
import traceback
import warnings
from itertools import chain
from typing import List, Any, Optional, Union
from io import StringIO

import numpy as np
from cpylog import get_logger2, WarningRedirector
#warnings.simplefilter('always')
warnings.simplefilter('default')

np.seterr(all='raise')

#from pyNastran.gui.qt_version import qt_version
#import PySide2
#import matplotlib
#matplotlib.use('Qt5Agg')

from pyNastran.utils import check_path
from pyNastran.utils.numpy_utils import integer_types
from pyNastran.bdf.errors import (
    #CrossReferenceError,
    CardParseSyntaxError, DuplicateIDsError, MissingDeckSections,
    UnsupportedCard,
    DisabledCardError,
    ReplicationError,
    EnvironmentVariableError,
)
from pyNastran.bdf.bdf import BDF, read_bdf
from pyNastran.bdf.subcase import Subcase
from pyNastran.bdf.mesh_utils.export_mcids import export_mcids, export_mcids_all
from pyNastran.bdf.mesh_utils.extract_bodies import extract_bodies
from pyNastran.bdf.mesh_utils.mass_properties import (
    mass_properties, mass_properties_nsm)  #, mass_properties_breakdown
from pyNastran.bdf.mesh_utils.forces_moments import get_temperatures_array
from pyNastran.bdf.mesh_utils.mpc_dependency import (
    get_mpc_node_ids, get_mpc_node_ids_c1,
    get_dependent_nid_to_components, get_mpcs)
from pyNastran.bdf.mesh_utils.loads import get_static_force_vector_from_subcase_id

from pyNastran.bdf.cards.dmig import NastranMatrix
from pyNastran.bdf.bdf_interface.compare_card_content import compare_card_content
#from pyNastran.bdf.mesh_utils.convert import convert
#from pyNastran.bdf.mesh_utils.remove_unused import remove_unused

import pyNastran.bdf.test
TEST_PATH = pyNastran.bdf.test.__path__[0]
#warnings.filterwarnings("ignore", category=DeprecationWarning)


def run_lots_of_files(filenames: List[str], folder: str='',
                      debug: bool=False,
                      xref: bool=True,
                      check: bool=True,
                      write_hdf5: bool=True,
                      punch: bool=False,
                      nastran: str='',
                      encoding: Optional[str]=None,
                      size: Union[int, List[int], None]=None,
                      post: Union[int, List[int], None]=None,
                      is_double: Union[bool, List[bool], None]=None,
                      sum_load: bool=True,
                      dev: bool=True,
                      crash_cards: Optional[List[str]]=None,
                      pickle_obj: bool=True, quiet: bool=False) -> List[str]:
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
    size : int / List[int], optional
        The field width of the model (8/16)
    is_double : bool / List[bool], optional
        Is this a double precision model?
            True : size = 16
            False : size = {8, 16}
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
        is_doubles = [False]
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
    #print('posts=%s' % posts)
    for sizei, is_doublei, posti in zip(sizes, is_doubles, posts):
        size_doubles_post.append((sizei, is_doublei, posti))

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
    log = get_logger2(log=None, debug=debug, encoding='utf-8')
    with WarningRedirector(log) as warn:
        for filename in filenames2:
            abs_filename = os.path.abspath(os.path.join(folder, filename))
            if folder != '':
                print("filename = %s" % abs_filename)
            is_passed = False
            try:
                for sizei, is_doublei, posti in size_doubles_post:
                    fem1, fem2, unused_diff_cards2 = run_bdf(
                        folder, filename, debug=debug,
                        xref=xref, check=check, punch=punch,
                        encoding=encoding,
                        is_folder=True, dynamic_vars={},
                        nastran=nastran, size=sizei, is_double=is_doublei,
                        nerrors=0,
                        post=posti, sum_load=sum_load, dev=dev,
                        crash_cards=crash_cards,
                        run_extract_bodies=False, pickle_obj=pickle_obj,
                        hdf5=write_hdf5, quiet=quiet, log=log)
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
            #except SyntaxError:  # only temporarily uncomment this when running lots of tests
                #pass

            #except AttributeError:  # only temporarily uncomment this when running lots of tests
                #pass
            #except KeyError:  # only temporarily uncomment this when running lots of tests
                #pass
            #except AssertionError:  # only temporarily uncomment this when running lots of tests
                #pass
            #except IndexError:  # only temporarily uncomment this when running lots of tests
                #pass
            #except ValueError:  # only temporarily uncomment this when running lots of tests
                #pass
            except ReplicationError:
                if not dev:
                    raise
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


def run_bdf(folder, bdf_filename, debug=False, xref=True, check=True, punch=False,
            mesh_form='separate', is_folder=False, print_stats=False,
            encoding=None, sum_load=True, size=8, is_double=False,
            hdf5=False,
            stop=False, nastran='', post=-1, dynamic_vars=None,
            quiet=False, dumplines=False, dictsort=False,
            run_extract_bodies=False, run_skin_solids=True,
            save_file_structure=False,
            nerrors=0, dev=False, crash_cards=None, safe_xref=False, pickle_obj=False,
            stop_on_failure=True, log=None):
    """
    Runs a single BDF

    Parameters
    ----------
    folder : str
        the folder where the bdf_filename is
    bdf_filename : str
        the bdf file to analyze
    debug : bool / None, default=False
        True : run with debug logging
        False : run with info logging
        None : run with warning logging
    xref : bool / str, optional
        True : cross reference the model
        False  : don't cross reference the model
        'safe' : do safe cross referencing
    check : bool, optional
        validate cards for things like mass, area, etc.
    punch : bool, optional
        this is a PUNCH file (no executive/case control decks)
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
            False : size = {8, 16}
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

    bdf_model = bdf_filename
    if not quiet:
        print("bdf_model = %s" % bdf_model)
    if is_folder:
        bdf_model = os.path.join(TEST_PATH, folder, bdf_filename)

    model, ext = os.path.splitext(bdf_model)
    out_model = '%s.test_bdf%s' % (model, ext)

    fem1, fem2, diff_cards = run_and_compare_fems(
        bdf_model, out_model, debug=debug, xref=xref, check=check,
        punch=punch, mesh_form=mesh_form,
        print_stats=print_stats, encoding=encoding,
        sum_load=sum_load, size=size, is_double=is_double,
        stop=stop, nastran=nastran, post=post, hdf5=hdf5,
        dynamic_vars=dynamic_vars,
        quiet=quiet, dumplines=dumplines, dictsort=dictsort,
        nerrors=nerrors, dev=dev, crash_cards=crash_cards,
        safe_xref=safe_xref,
        run_extract_bodies=run_extract_bodies,
        run_skin_solids=run_skin_solids,
        save_file_structure=save_file_structure,
        pickle_obj=pickle_obj,
        stop_on_failure=stop_on_failure,
        log=log,
    )
    return fem1, fem2, diff_cards

def run_and_compare_fems(
        bdf_model, out_model, debug=False, xref=True, check=True,
        punch=False, mesh_form='combined',
        print_stats=False, encoding=None,
        sum_load=True, size=8, is_double=False,
        save_file_structure=False,
        stop=False, nastran='', post=-1, hdf5=False,
        dynamic_vars=None,
        quiet=False, dumplines=False, dictsort=False,
        nerrors=0, dev=False, crash_cards=None,
        safe_xref=True,
        run_extract_bodies=False,
        run_skin_solids=True, pickle_obj=False,
        stop_on_failure=True, log=None,
    ):
    """runs two fem models and compares them"""
    assert os.path.exists(bdf_model), '%r doesnt exist' % bdf_model
    fem1 = BDF(debug=debug, log=log)
    fem1.dumplines = dumplines

    fem1.set_error_storage(nparse_errors=nerrors, stop_on_parsing_error=True,
                           nxref_errors=nerrors, stop_on_xref_error=True)
    if dynamic_vars:
        fem1.set_dynamic_syntax(dynamic_vars)

    if not quiet:
        fem1.log.info('starting fem1')
    sys.stdout.flush()
    fem2 = None
    diff_cards = []

    mesh_opt_cards = [
        #'GRIDG', 'CGEN', 'SPCG', 'EQUIV', 'FEEDGE', 'FEFACE', 'ADAPT',
        #'PVAL', 'GMCURV', 'GMSURF',
    ]
    #nastran_cmd = 'nastran scr=yes bat=no old=no news=no '
    nastran_cmd = ''
    is_mesh_opt = False
    try:
        #try:

        fem1 = run_fem1(fem1, bdf_model, out_model, mesh_form, xref, punch, sum_load,
                        size, is_double,
                        run_extract_bodies=run_extract_bodies,
                        run_skin_solids=run_skin_solids,
                        save_file_structure=save_file_structure,
                        hdf5=hdf5,
                        encoding=encoding, crash_cards=crash_cards, safe_xref=safe_xref,
                        pickle_obj=pickle_obj, stop=stop)
        is_mesh_opt = any([card_name in fem1.card_count for card_name in mesh_opt_cards])
        if dev and is_mesh_opt:
            return None, None, None

        if stop:
            if not quiet:
                print('card_count:')
                print('-----------')
                for card_name, card_count in sorted(fem1.card_count.items()):
                    print('key=%-8s value=%s' % (card_name, card_count))
            return fem1, None, None

        ierror = 0
        fem2 = run_fem2(bdf_model, out_model, xref, punch, sum_load, size, is_double, mesh_form,
                        safe_xref=safe_xref,
                        encoding=encoding, debug=debug, quiet=quiet,
                        ierror=ierror, nerrors=nerrors,
                        stop_on_failure=stop_on_failure, log=log)

        diff_cards = compare(fem1, fem2, xref=xref, check=check,
                             print_stats=print_stats, quiet=quiet)
        test_get_cards_by_card_types(fem2)

        fem2.update_model_by_desvars(xref)
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
    except UnsupportedCard:
        if not dev:
            raise
        print('failed test because UnsupportedCard /...ignoring')
    except MissingDeckSections:
        if not dev:
            raise
        print('failed test because MissingDeckSections...ignoring')
    except DuplicateIDsError as e:
        # only temporarily uncomment this when running lots of tests
        if not dev:
            raise
        elif is_mesh_opt:
            print('failed test because mesh adaption (GRIDG,CGEN,SPCG)...ignoring')
            print(e)
        else:
            print('failed test because DuplicateIDsError...ignoring')
    except DisabledCardError as e:
        if not dev:
            raise
    except EnvironmentVariableError:
        if not dev:
            raise
    except RuntimeError as e:
        # only temporarily uncomment this when running lots of tests
        if not dev:
            raise
        elif is_mesh_opt:
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
        elif is_mesh_opt:
            print('failed test because mesh adaption (GRIDG,CGEN,SPCG)...ignoring')
            print(e)
        else:
            raise
    except KeyError as e:  # only temporarily uncomment this when running lots of tests
        if not dev:
            raise
        elif is_mesh_opt:
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


def run_nastran(bdf_model: BDF, nastran: str, post: int=-1,
                size: int=8, is_double: bool=False) -> None:
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

def run_fem1(fem1: BDF, bdf_model: str, out_model: str, mesh_form: str,
             xref: bool, punch: bool, sum_load: bool,
             size: int, is_double: bool,
             run_extract_bodies: bool=False, run_skin_solids: bool=True,
             save_file_structure: bool=False, hdf5: bool=False,
             encoding: Optional[str]=None, crash_cards: Optional[List[str]]=None,
             safe_xref: bool=True, pickle_obj: bool=False, stop: bool=False) -> BDF:
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
            fem1.read_bdf(bdf_model, xref=False, punch=True, encoding=encoding,
                          save_file_structure=save_file_structure)
        else:
            fem1.read_bdf(bdf_model, xref=False, punch=punch, encoding=encoding,
                          save_file_structure=save_file_structure)
            for card in crash_cards:
                if card in fem1.card_count:
                    raise DisabledCardError('card=%r has been disabled' % card)
            #fem1.geom_check(geom_check=True, xref=False)
            if not stop and not xref and run_skin_solids:
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
                #get_dependent_nid_to_components(fem1)
                #fem1._get_maps(eids=None, map_names=None,
                               #consider_0d=True, consider_0d_rigid=True,
                               #consider_1d=True, consider_2d=True, consider_3d=True)
                #get_dependent_nid_to_components(fem1)

                #fem1.uncross_reference()
                if safe_xref:
                    fem1.safe_cross_reference()
                else:
                    fem1.cross_reference()

                # 1. testing that these methods work with xref
                fem1._get_rigid()
                common_node_ids = list(fem1.nodes.keys())
                fem1.get_rigid_elements_with_node_ids(common_node_ids)

                for spc_id in set(list(fem1.spcadds.keys()) + list(fem1.spcs.keys())):
                    fem1.get_reduced_spcs(spc_id, consider_spcadd=True)
                for mpc_id in set(list(fem1.mpcadds.keys()) + list(fem1.mpcs.keys())):
                    fem1.get_reduced_mpcs(mpc_id, consider_mpcadd=True)

                get_dependent_nid_to_components(fem1)
                fem1._get_maps(eids=None, map_names=None,
                               consider_0d=True, consider_0d_rigid=True,
                               consider_1d=True, consider_2d=True, consider_3d=True)
                get_dependent_nid_to_components(fem1)
                fem1.get_pid_to_node_ids_and_elements_array(pids=None, etypes=None, idtype='int32',
                                                            msg=' which is required by test_bdf')
                fem1.get_property_id_to_element_ids_map(msg=' which is required by test_bdf')
                fem1.get_material_id_to_property_ids_map(msg=' which is required by test_bdf')
                fem1.get_element_ids_list_with_pids(pids=None)
                fem1.get_element_ids_dict_with_pids(pids=None, stop_if_no_eids=False,
                                                    msg=' which is required by test_bdf')
                fem1.get_node_id_to_element_ids_map()
                fem1.get_node_id_to_elements_map()

                export_mcids(fem1, csv_filename=None, eids=None,
                             export_xaxis=True, export_yaxis=False, iply=0, log=None, debug=False)

                export_mcids_all(fem1)

                fem1._xref = True
                read_bdf(fem1.bdf_filename, encoding=encoding, xref=False,
                         debug=fem1.debug, log=fem1.log)
                if safe_xref:
                    fem1.safe_cross_reference()
                elif xref:
                    fem1.cross_reference()

                fem1 = remake_model(bdf_model, fem1, pickle_obj)
                #fem1.geom_check(geom_check=True, xref=True)
                #fem1.uncross_reference()
                #fem1.cross_reference()
    except:
        print("failed reading %r" % bdf_model)
        raise

    #out_model = bdf_model + '_out'
    #if cid is not None and xref:
        #fem1.resolve_grids(cid=cid)

    if hdf5:
        hdf5_filename = out_model + '.h5'
        fem1.export_hdf5_filename(hdf5_filename)
        fem1a = BDF(log=fem1.log)
        fem1a.load_hdf5_filename(hdf5_filename)
        fem1a.validate()
        bdf_stream = StringIO()
        fem1a.write_bdf(bdf_stream, encoding=None, size=8,
                        is_double=False, interspersed=False,
                        enddata=None, write_header=True, close=True) # hdf5
        for key, unused_value in fem1.card_count.items():
            if key in ['ECHOOFF', 'ECHOON']:
                continue
            #if key == 'ENDDATA':
            hdf5_msg = ''
            if key not in fem1a.card_count:
                hdf5_msg += 'key=%r was not loaded to hdf5\n' % key

            if hdf5_msg:
                hdf5_msg += 'expected=%s\nactual=%s' % (
                    fem1.card_count, fem1a.card_count)
                fem1a.log.error(hdf5_msg)
                raise RuntimeError(hdf5_msg)
        #sys.exit('hdf5')

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

        if len(fem1.elements) == 0 and len(fem1.masses) == 0:
            fem1.log.warning('no elements found')
        elif len(fem1.elements) == 0:
            fem1.get_mass_breakdown(stop_if_no_mass=False)
            fem1.log.warning('no elements with length/area/volume found, but elements with mass were')
        else:
            # len(elements) > 0 or len(masses) > 0
            fem1.get_length_breakdown(stop_if_no_length=False)
            fem1.get_area_breakdown(stop_if_no_area=False)
            fem1.get_volume_breakdown(stop_if_no_volume=False)
            fem1.get_mass_breakdown(stop_if_no_mass=False)
    return fem1


def remake_model(bdf_model: BDF, fem1: BDF, pickle_obj: bool) -> None:
    """reloads the model if we're testing pickling"""
    remake = pickle_obj
    if remake:
        #log = fem1.log
        model_name = os.path.splitext(bdf_model)[0]
        obj_model = '%s.test_bdf.obj' % (model_name)
        #out_model_8 = '%s.test_bdf.bdf' % (model_name)
        #out_model_16 = '%s.test_bdf.bdf' % (model_name)

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

        fem1.cross_reference()
        #fem1.get_bdf_stats()
        fem1._xref = True
    return fem1

def check_for_cd_frame(fem1: BDF) -> None:
    """
    A cylindrical/spherical CD frame will cause problems with the
    grid point force transformation

    """
    is_grid_points = any([card_name in fem1.card_count
                          for card_name in ['GRID', 'SPOINT', 'EPOINT', 'RINGAX']])
    if is_grid_points:
        out = fem1.get_displacement_index_xyz_cp_cd(
            fdtype='float64', idtype='int64', sort_ids=True)
        unused_icd_transform, unused_icp_transform, unused_xyz_cp, nid_cp_cd = out
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

def run_fem2(bdf_model: str, out_model: str, xref: bool, punch: bool,
             sum_load: bool, size: int, is_double: bool, mesh_form: str,
             safe_xref: bool=False,
             encoding: Optional[str]=None, debug: bool=False, quiet: bool=False,
             stop_on_failure: bool=True, ierror: int=0, nerrors: int=100, log=None) -> BDF:
    """
    Reads/writes the BDF to verify nothing has been lost

    Parameters
    ----------
    bdf_model : str
        the filename to run
    out_model : ???
        ???
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

    fem2 = BDF(debug=debug, log=log)
    if not quiet:
        fem2.log.info('starting fem2')
    sys.stdout.flush()
    try:
        fem2.read_bdf(out_model, xref=False, punch=punch, encoding=encoding)
        if safe_xref:
            fem2.safe_cross_reference()
        elif xref:
            fem2.cross_reference()
    except:
        print("failed reading %r" % out_model)
        raise

    out_model_2 = bdf_model + '_out2'

    if xref and sum_load:
        if 'POST' in fem2.params:
            value = fem2.params['POST'].values[0]
            if value >= 0:
                msg = 'PARAM,POST,%i is not supported by the OP2 reader' % value
                fem2.log.warning(msg)
        else:
            msg = 'PARAM,POST,0 is not supported by the OP2 reader'
            fem2.log.warning(msg)

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
            ierror = validate_case_control(
                fem2, p0, sol_base, subcase_keys, subcases, sol_200_map,
                ierror=ierror, nerrors=nerrors,
                stop_on_failure=stop_on_failure)

    if mesh_form is not None:
        fem2.write_bdf(out_model_2, interspersed=False, size=size, is_double=is_double)
        try:
            os.remove(out_model_2)
        except PermissionError:  # pragma: no cover
            fem2.log.warning('cannot remove %s due to a permissions error' % out_model_2)
    #fem2.write_as_ctria3(out_model_2)
    return fem2

def _assert_has_spc(subcase, fem):
    """
    SPCs may be defined on SPC/SPC1 cards or may be defined on
    the GRID PS field

    """
    if 'SPC' not in subcase:
        has_ps = False
        for unused_nid, node in fem.nodes.items():
            if node.ps:
                has_ps = True
                break
        assert subcase.has_parameter('SPC', 'STATSUB') or has_ps, subcase

def validate_case_control(fem2: BDF, p0: Any, sol_base: int, subcase_keys: List[int],
                          subcases: Any, unused_sol_200_map: Any,
                          stop_on_failure: bool=True,
                          ierror: int=0, nerrors: int=100) -> None:
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
        ierror = check_case(
            sol_base, subcase, fem2, p0, isubcase, subcases,
            ierror=ierror, nerrors=nerrors, stop_on_failure=stop_on_failure)
    return ierror

def check_for_flag_in_subcases(fem2: BDF, subcase: Any, parameters: List[str]) -> None:
    """
    For a multi-subcase deck, you can define specific required cards
    (e.g., TSTEP) in secondary cases, but not primary cases.  This
    results is a non-transient case being defined for the first (1 or more)
    subcase, but the final subcase must have the TSTEP card.

    You can also use this for things like post-buckling and dynamic time
    stepping.  In the dynamic time stepping case, for the first 3 seconds
    you can define one set of time stepping, but then you can switch to
    a much finer time step/increased number of convergence iterations.

    """
    if not any(subcase.has_parameter(*parameters)):
        unused_subcases = fem2.subcases
        #subcase_ids = [isubcase for isubcase in subcases if isubcase > 0]
        has_flag = False
        for unused_isubcase, subcasei in fem2.subcases.items():
            if any(subcasei.has_parameter(*parameters)):
                has_flag = True
        if not has_flag:
            msg = 'sol=%r; %s not in subcase\n' % (fem2.sol, str(parameters))
            for unused_isubcase, subcasei in fem2.subcases.items():
                msg += str(subcasei)
            raise RuntimeError(msg)

def stop_if_max_error(msg: str, error: Any, ierror: int, nerrors: int) -> int:
    """if the error count is greater than nerrors, stop"""
    if ierror == nerrors:
        raise error(msg)
    ierror += 1
    return ierror

def check_for_optional_param(keys: List[str], subcase: Any,
                             msg: str, error: Any, log: Any, ierror: int, nerrors: int) -> int:
    """one or more must be True"""
    if not any(subcase.has_parameter(*keys)):
        msg = 'Must have one of %s\n%s' % (str(keys), msg)
        log.error(msg)
        if ierror == nerrors:
            raise error(msg)
        ierror += 1
    return ierror

def check_sol(sol: int,
              subcase: Any,
              allowed_sols: List[int],
              case_control_key: str,
              log: Any, ierror: int, nerrors: int) -> int:
    """Checks that the solution is valid"""
    if sol not in allowed_sols:
        msg = '%s is not valid in sol=%s allowed_sols=%s\n%s' % (
            case_control_key, sol, allowed_sols, subcase)
        log.error(msg)
        if ierror == nerrors:
            raise RuntimeError(msg)
    if case_control_key not in subcase:
        msg = 'sol=%s is missing %r\n%s' % (sol, case_control_key, subcase)
        log.error(msg)
        if ierror == nerrors:
            raise RuntimeError(msg)
        ierror += 1
    return ierror

def check_case(sol, subcase, fem2, p0, isubcase, subcases,
               ierror=0, nerrors=100, stop_on_failure=True):
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
    log = fem2.log

    if sol == 24:
        _assert_has_spc(subcase, fem2)
        assert True in subcase.has_parameter('LOAD'), 'sol=%s\n%s' % (sol, subcase)
    elif sol == 64:
        #assert 'NLPARM' in subcase, subcase
        #_assert_has_spc(subcase, fem2)
        assert True in subcase.has_parameter('LOAD'), 'sol=%s\n%s' % (sol, subcase)
    elif sol == 66:
        assert 'NLPARM' in subcase, subcase
        _assert_has_spc(subcase, fem2)
        msg = 'sol=%s\n%s' % (sol, subcase)
        ierror = check_for_optional_param(('LOAD', 'TEMPERATURE(LOAD)'), subcase, msg,
                                          RuntimeError, log, ierror, nerrors)
    elif sol == 99:
        assert 'DLOAD' in subcase, subcase
        assert 'LOADSET' in subcase, subcase
        _assert_has_spc(subcase, fem2)
        #assert True in subcase.has_parameter('LOAD', 'TEMPERATURE'), 'sol=%s\n%s' % (sol, subcase)
        msg = 'sol=%s\n%s' % (sol, subcase)
        ierror = check_for_optional_param(('TSTEP', 'TSTEPNL'), subcase, msg,
                                          RuntimeError, log, ierror, nerrors)
    elif sol in [1, 101]:
        _assert_has_spc(subcase, fem2)
        msg = 'sol=%s\n%s' % (sol, subcase)
        ierror = check_for_optional_param(('LOAD', 'TEMPERATURE(LOAD)', 'P2G'), subcase, msg,
                                          RuntimeError, log, ierror, nerrors)
    elif sol in [3, 103]:
        msg = 'sol=%s\n%s' % (sol, subcase)
        ierror = check_for_optional_param(('METHOD', 'RSMETHOD', 'RIGID', 'BOLTID', 'BGSET'), subcase, msg,
                                          RuntimeError, log, ierror, nerrors)
    elif sol == 105: # buckling
        _assert_has_spc(subcase, fem2)
        msg = 'sol=%s\n%s' % (sol, subcase)
        ierror = check_for_optional_param(('LOAD', 'METHOD'), subcase, msg,
                                          RuntimeError, log, ierror, nerrors)
        #if 0:  # pragma: no cover
            #if 'METHOD' not in subcase:
                #subcases = fem2.subcases
                #subcase_ids = [isubcase for isubcase in subcases if isubcase > 0]
                #assert len(subcases) == 2, 'METHOD not in subcase and not 2 subcases\n%s' % subcase
                #subcase_id = subcase.subcase_id
                #if subcase_id == 1 and 'METHOD' in subcases[2]:
                    #pass
                #else:
                    #msg = 'METHOD not in subcase and not 2 subcases\n%s' % subcase
                    #raise RuntimeError(msg)

        #assert True in subcase.has_parameter('LOAD', 'TEMPERATURE(LOAD)'), 'sol=%s\n%s' % (sol, subcase)
    elif sol == 106: # freq
        assert 'NLPARM' in subcase, subcase
        msg = 'sol=%s\n%s' % (sol, subcase)
        ierror = check_for_optional_param(('LOAD', 'TEMPERATURE(LOAD)', 'CLOAD'), subcase, msg,
                                          RuntimeError, log, ierror, nerrors)
    elif sol == 107: # ???
        _assert_has_spc(subcase, fem2)
        msg = 'sol=%s\n%s' % (sol, subcase)
        ierror = check_for_optional_param(('LOAD', 'TEMPERATURE(LOAD)'), subcase, msg,
                                          RuntimeError, log, ierror, nerrors)
    elif sol in [8, 108]: # freq
        assert 'FREQUENCY' in subcase, subcase
    elif sol == 109:  # time
        check_for_flag_in_subcases(fem2, subcase, ('TIME', 'TSTEP', 'TSTEPNL'))

    elif sol == 110:  # ???
        _assert_has_spc(subcase, fem2)
        msg = 'sol=%s\n%s' % (sol, subcase)
        ierror = check_for_optional_param(('LOAD', 'STATSUB'), subcase, msg,
                                          RuntimeError, log, ierror, nerrors)
    elif sol == 111:  # modal frequency
        assert subcase.has_parameter('FREQUENCY'), 'sol=%s\n%s' % (sol, subcase)
        assert any(subcase.has_parameter('METHOD', 'RMETHOD')), 'sol=%s\n%s' % (sol, subcase)
    elif sol == 112:  # modal transient
        check_for_flag_in_subcases(fem2, subcase, ('TIME', 'TSTEP', 'TSTEPNL'))
        #assert any(subcase.has_parameter('TIME', 'TSTEP', 'TSTEPNL')), 'sol=%s\n%s' % (sol, subcase)
    elif sol == 114:
        soltype = 'CYCSTATX'
        ierror = require_cards(['LOAD', 'HARMONICS'], log, soltype, sol, subcase,
                               RuntimeError, ierror, nerrors)
        _assert_has_spc(subcase, fem2)
    elif sol == 118:
        soltype = 'CYCFREQ'
        ierror = require_cards(['LOAD', 'HARMONICS', 'SDAMPING', 'FREQUENCY'],
                               log, soltype, sol, subcase,
                               RuntimeError, ierror, nerrors)
        _assert_has_spc(subcase, fem2)

    elif sol == 129:  # nonlinear transient
        assert any(subcase.has_parameter('TIME', 'TSTEP', 'TSTEPNL')), 'sol=%s\n%s' % (sol, subcase)
    elif sol == 159:  # thermal transient
        assert any(subcase.has_parameter('TIME', 'TSTEP', 'TSTEPNL')), 'sol=%s\n%s' % (sol, subcase)

    elif sol == 144:
        ierror = _check_static_aero_case(fem2, log, sol, subcase, ierror, nerrors)
    elif sol == 145:
        ierror = _check_flutter_case(fem2, log, sol, subcase, ierror, nerrors)
    elif sol == 146:
        ierror = _check_gust_case(fem2, log, sol, subcase, ierror, nerrors)

    elif sol == 153: # heat?
        _assert_has_spc(subcase, fem2)
        if 'NLPARM' not in subcase:
            log.error('A NLPARM card is required for HEAT? - SOL %i\n%s' % (sol, subcase))

        if 'ANALYSIS' in subcase and subcase.get_parameter('ANALYSIS')[0] == 'HEAT':
            assert any(subcase.has_parameter('TEMPERATURE(LOAD)', 'TEMPERATURE(INITIAL)')), 'sol=%s\n%s' % (sol, subcase)
        else:
            assert any(subcase.has_parameter('LOAD', 'TEMPERATURE(LOAD)')), 'sol=%s\n%s' % (sol, subcase)

    elif sol == 159: #  nonlinear transient; heat?
        if 'NLPARM' not in subcase:
            msg = (
                'A NLPARM card is required for NONLINEAR_TRANSIENT? '
                '- SOL %i\n%s' % (sol, subcase))
            log.error(msg)
            ierror = stop_if_max_error(msg, RuntimeError, ierror, nerrors)

        #assert any(subcase.has_parameter('TIME', 'TSTEP', 'TSTEPNL')), subcase
        #assert any(subcase.has_parameter('GUST', 'LOAD')), subcase
        if 'ANALYSIS' in subcase and subcase.get_parameter('ANALYSIS')[0] == 'HEAT':
            assert any(subcase.has_parameter('TEMPERATURE(LOAD)', 'TEMPERATURE(INITIAL)')), 'sol=%s\n%s' % (sol, subcase)

    elif sol == 200:
        _check_case_sol_200(sol, subcase, fem2, p0, isubcase, subcases, log)
    elif sol in [114, 115, 116, 118]:
        # cyclic statics, modes, buckling, frequency
        pass
    elif sol in [5, 21, 26, 38, 47, 61, 63, 68, 76, 78, 81, 88,
                 100, 128, 187, 190,
                 400, 401, 402, 600, 601, 700, 701, 'AEDB2XDB', 'UPWARD']:
        pass
    else:
        msg = 'SOL = %s\n' % (sol)
        msg += str(subcase)
        raise NotImplementedError(msg)
    ierror = _check_case_parameters(
        subcase, fem2, p0, isubcase, sol,
        ierror=ierror, nerrors=nerrors,
        stop_on_failure=stop_on_failure)
    return ierror

def _check_static_aero_case(fem2: BDF, log: Any, sol: int,
                            subcase: Any, ierror: int, nerrors: int) -> int:
    """checks that TRIM/DIVERG is valid"""
    if not any(subcase.has_parameter('TRIM', 'DIVERG')):
        msg = 'A TRIM or DIVERG card is required for STATIC AERO - SOL %i\n%s' % (
            sol, subcase)
        log.error(msg)
        ierror = stop_if_max_error(msg, RuntimeError, ierror, nerrors)
    if fem2.aeros is None:
        msg = 'An AEROS card is required for STATIC AERO - SOL %i; AEROS=%s' % (sol, fem2.aeros)
        log.error(msg)
        ierror = stop_if_max_error(msg, RuntimeError, ierror, nerrors)
    if len(fem2.caeros) == 0:
        msg = 'An CAEROx card is required for STATIC AERO - SOL %i' % (sol)
        log.error(msg)
        ierror = stop_if_max_error(msg, RuntimeError, ierror, nerrors)
    if len(fem2.splines) == 0:
        msg = 'An SPLINEx card is required for STATIC AERO - SOL %i' % (sol)
        log.error(msg)
        ierror = stop_if_max_error(msg, RuntimeError, ierror, nerrors)
    return ierror

def _check_flutter_case(fem2, log, sol, subcase, ierror, nerrors):
    # type: (BDF, Any, int, Any, int, int) -> int
    """checks that FLUTTER is valid"""
    if fem2.aero is None:
        msg = 'An AERO card is required for FLUTTER - SOL %i; AERO=%s' % (sol, fem2.aero)
        log.error(msg)
        ierror = stop_if_max_error(msg, RuntimeError, ierror, nerrors)

    if len(fem2.caeros) == 0:
        msg = 'An CAEROx card is required for FLUTTER - SOL %i' % (sol)
        log.error(msg)
        ierror = stop_if_max_error(msg, RuntimeError, ierror, nerrors)
    if len(fem2.splines) == 0:
        msg = 'An SPLINEx card is required for FLUTTER - SOL %i' % (sol)
        log.error(msg)
        ierror = stop_if_max_error(msg, RuntimeError, ierror, nerrors)
    if len(fem2.mkaeros) == 0:
        msg = 'An MKAERO1/2 card is required for FLUTTER - SOL %i' % (sol)
        log.error(msg)
        ierror = stop_if_max_error(msg, RuntimeError, ierror, nerrors)
    unused_mklist = fem2.get_mklist()

    soltype = 'FLUTTER'
    # METHOD - EIGRL
    # CMETHOD - EIGC
    # FMETHOD - FLUTTER
    ierror = require_cards(['FMETHOD'], log, soltype, sol, subcase,
                           RuntimeError, ierror, nerrors)
    flutter_id = subcase.get_parameter('FMETHOD')[0]
    flutter = fem2.Flutter(flutter_id, msg=', which is required by test_bdf')

    #valid methods = [K, KE,
                     #PKS, PKNLS, PKNL, PKE]
    #if flutter.method in ['PK', 'PKNL']: # not supported in SOL 200
    if flutter.method == 'K':
        # EIGC
        ierror = require_cards(['CMETHOD'], log, soltype, sol, subcase,
                               RuntimeError, ierror, nerrors)
    else:
        # EIGRL
        ierror = require_cards(['METHOD'], log, soltype, sol, subcase,
                               RuntimeError, ierror, nerrors)
    return ierror

def _check_gust_case(fem2, log, sol, subcase, ierror, nerrors):
    # type: (BDF, Any, int, Any, int, int) -> int
    """checks that GUST is valid"""
    if 'METHOD' not in subcase:  # EIGRL
        msg = 'A METHOD card is required for FLUTTER - SOL %i\n%s' % (sol, subcase)
        log.error(msg)
        ierror = stop_if_max_error(msg, RuntimeError, ierror, nerrors)
    if not any(subcase.has_parameter('FREQUENCY', 'TIME', 'TSTEP', 'TSTEPNL')):
        msg = (
            'A FREQUENCY/TIME/TSTEP/TSTEPNL card is required for GUST'
            ' - SOL %i\n%s' % (sol, subcase))
        log.error(msg)
        ierror = stop_if_max_error(msg, RuntimeError, ierror, nerrors)
    if not any(subcase.has_parameter('GUST', 'LOAD', 'DLOAD')):
        msg = 'A GUST/LOAD/DLOAD card is required for GUST - SOL %i\n%s' % (sol, subcase)
        log.error(msg)
        ierror = stop_if_max_error(msg, RuntimeError, ierror, nerrors)
    if fem2.aero is None:
        msg = 'An AERO card is required for GUST - SOL %i' % sol
        log.error(msg)
        ierror = stop_if_max_error(msg, RuntimeError, ierror, nerrors)

    if len(fem2.caeros) == 0:
        msg = 'An CAEROx card is required for GUST - SOL %i' % (sol)
        log.error(msg)
        ierror = stop_if_max_error(msg, RuntimeError, ierror, nerrors)
    if len(fem2.splines) == 0:
        msg = 'An SPLINEx card is required for GUST - SOL %i' % (sol)
        log.error(msg)
        ierror = stop_if_max_error(msg, RuntimeError, ierror, nerrors)
    if len(fem2.mkaeros) == 0:
        msg = 'An MKAERO1/2 card is required for GUST - SOL %i' % (sol)
        log.error(msg)
        ierror = stop_if_max_error(msg, RuntimeError, ierror, nerrors)
    unused_mklist = fem2.get_mklist()
    return ierror

def _check_case_sol_200(sol: int,
                        subcase: Subcase,
                        fem2: BDF,
                        p0: Any,
                        isubcase: int, subcases: int,
                        log: Any):
    """
    helper method for ``check_case``

    local level
    DESSUB - Set constraints (DCONSTR, DCONADD) applied for subcase
             (e.g. STRESS, STRAIN, DAMP)
             optional locally
    DESGLB - Set constraints (DCONSTR, DCONADD) applied globally
             (e.g. WEIGHT, VOLUME, WMPID, FRMASS)
             optional locally
    DESOBJ - The objective function (DRESP1, DRESP2, DRESP3)
             required globally
    1 or more DESSUB/DESGLB are required globally
    1 DESOBJ is required

    """
    assert 'ANALYSIS' in subcase, 'sol=%s\n%s' % (sol, subcase)

    analysis = subcase.get_parameter('ANALYSIS')[0]
    # BUCKLING
    if 'DESOBJ' in subcase:
        value = subcase.get_parameter('DESOBJ')[0]
        assert value in fem2.dresps, 'value=%s not in dresps' % value
    else:
        fem2.log.warning('no DESOBJ (DRESPi) in this subcase; is this a buckling preload case?')
        fem2.log.warning('\n%s' % subcase)

    nopt = len(fem2.dvprels) + len(fem2.dvmrels) + len(fem2.dvcrels)
    if nopt == 0:
        fem2.log.error('no DVPRELs/DVMRELs/DVCRELs found')

    #--------------------------------------------------------------------------
    # DCONSTR
    if 'DESSUB' not in subcase and 'DESGLB' not in subcase:
        fem2.log.warning('no DESSUB/DESGLB (DCONSTR) in this subcase;'
                         ' is this a buckling preload case?')
        log.warning('\n%s' % subcase)

    if 'DESSUB' in subcase:
        value = subcase.get_parameter('DESSUB')[0]
        if value not in fem2.dconstrs:
            msg = 'value=%s not in dconstrs; Allowed DCONSTRs=%s' % (
                value, np.unique(list(fem2.dconstrs.keys())))
            raise RuntimeError(msg)
    if 'DESGLB' in subcase:
        value = subcase.get_parameter('DESGLB')[0]
        if value not in fem2.dconstrs:
            msg = 'value=%s not in dconstrs; Allowed DCONSTRs=%s' % (
                value, np.unique(list(fem2.dconstrs.keys())))
            raise RuntimeError(msg)
    #--------------------------------------------------------------------------

    if analysis in ['STATIC', 'STATICS']:
        solution = 101
        check_case(solution, subcase, fem2, p0, isubcase, subcases)
    elif analysis in ['MODE', 'MODES']:
        solution = 103
        check_case(solution, subcase, fem2, p0, isubcase, subcases)
    elif analysis in ['BUCK', 'BUCKLING']:
        solution = 105
        check_case(solution, subcase, fem2, p0, isubcase, subcases)
    elif analysis == 'DFREQ':
        solution = 108
        check_case(solution, subcase, fem2, p0, isubcase, subcases)
    elif analysis == 'MFREQ':
        if 'GUST' in subcase:
            solution = 146
        else:
            solution = 111
        check_case(solution, subcase, fem2, p0, isubcase, subcases)
    elif analysis == 'MTRAN':
        solution = 112
        check_case(solution, subcase, fem2, p0, isubcase, subcases)
    elif analysis in ['SAERO', 'DIVERG', 'DIVERGE']:
        solution = 144
        check_case(solution, subcase, fem2, p0, isubcase, subcases)
    elif analysis == 'FLUTTER':
        solution = 145
        check_case(solution, subcase, fem2, p0, isubcase, subcases)
    elif analysis == 'DCEIG': # direct complex eigenvalues
        solution = 107
        check_case(solution, subcase, fem2, p0, isubcase, subcases)
    #elif analysis == 'MCEIG': # modal direct complex eigenvalues
    elif analysis == 'HEAT': # heat transfer analysis
        solution = 159
        check_case(solution, subcase, fem2, p0, isubcase, subcases)
    elif analysis == 'MCEIG': # modal complex eigenvalues
        solution = 110
        check_case(solution, subcase, fem2, p0, isubcase, subcases)
    else:
        msg = 'analysis = %s\nsubcase =\n%s' % (analysis, subcase)
        raise NotImplementedError(msg)

def require_cards(card_names, log, soltype, sol, subcase,
                  error, ierror, nerrors):
    """all must be True"""
    for card_name in card_names:
        if card_name not in subcase:
            msg = 'A %s card is required for %s - SOL %i\n%s' % (
                card_name, soltype, sol, subcase)
            log.error(msg)
            if ierror == nerrors:
                raise error(msg)
            ierror += 1
    return ierror

def _tstep_msg(fem, subcase, tstep_id, tstep_type=''):
    """writes the TSTEP/TSTEPNL info"""
    msg = (
        'tstep%s_id=%s\n'
        'tsteps=%s\n'
        'tstepnls=%s\n'
        #'tstep1s=%s\n'
        'subcase:\n%s' % (
            tstep_type, tstep_id,
            str(fem.tsteps),
            str(fem.tstepnls),
            #str(fem.tstep1s),
            str(subcase))
    )
    return msg


def _check_case_parameters(subcase, fem2: BDF, p0, isubcase: int, sol: int,
                           ierror: int=0, nerrors: int=100, stop_on_failure: bool=True) -> int:
    """helper method for ``check_case``"""
    log = fem2.log
    if fem2.sol in [401, 402]:
        # TSTEP references a TSTEP1, but not a TSTEP
        # TSTEP1s are stored in tstepnls
        if any(subcase.has_parameter('TIME', 'TSTEP')):
            if 'TSTEP' in subcase:
                tstep_id = subcase.get_parameter('TSTEP')[0]
            else:  # pragma: no cover
                raise NotImplementedError(subcase)
            if tstep_id not in fem2.tstepnls:
                raise RuntimeError(_tstep_msg(fem2, subcase, tstep_id))
        else:
            raise RuntimeError(f'missing TSTEP in case control deck\n{subcase}')
    else:
        if any(subcase.has_parameter('TIME', 'TSTEP')):
            if 'TIME' in subcase:
                tstep_id = subcase.get_parameter('TIME')[0]
            elif 'TSTEP' in subcase:
                tstep_id = subcase.get_parameter('TSTEP')[0]
            else:  # pragma: no cover
                raise NotImplementedError(subcase)
            if tstep_id not in fem2.tsteps:
                raise RuntimeError(_tstep_msg(fem2, subcase, tstep_id))

    if 'TSTEPNL' in subcase:
        tstepnl_id = subcase.get_parameter('TSTEPNL')[0]
        assert tstepnl_id in fem2.tstepnls, _tstep_msg(fem2, subcase, tstepnl_id, tstep_type='nl')

    if 'SUPORT1' in subcase:
        suport1_id = subcase.get_parameter('SUPORT1')[0]
        assert suport1_id in fem2.suport1, 'suport1_id=%s\n suport1=%s\n subcase:\n%s' % (suport1_id, str(fem2.suport1), str(subcase))

    if 'TRIM' in subcase:
        trim_id = subcase.get_parameter('TRIM')[0]
        if trim_id not in fem2.trims:
            msg = (
                'TRIM = %s\n'
                'trims=%s\n'
                'subcase:\n%s' % (trim_id, str(fem2.trims), str(subcase)))
            raise RuntimeError(msg)
        trim = fem2.trims[trim_id]

        suport1 = None
        if 'SUPORT1' in subcase:
            suport_id = subcase.get_parameter('SUPORT1')[0]
            suport1 = fem2.suport1[suport_id]
        try:
            trim.verify_trim(
                fem2.suport, suport1, fem2.aestats, fem2.aeparams,
                fem2.aelinks, fem2.aesurf, xref=True)
        except RuntimeError:
            if stop_on_failure or ierror == nerrors:
                raise
            ierror += 1
            exc_info = sys.exc_info()
            traceback.print_exception(*exc_info)
            #traceback.print_stack()
            #fem2.log.error(e.msg)
            #raise
        assert 'DIVERG' not in subcase, subcase
        #allowed_sols = [144, 200]

    if 'DIVERG' in subcase:
        value = subcase.get_parameter('DIVERG')[0]
        assert value in fem2.divergs, 'value=%s\n divergs=%s\n subcase:\n%s' % (value, str(fem2.divergs), str(subcase))
        assert 'TRIM' not in subcase, subcase
        #allowed_sols = [144, 200]

    if 'METHOD' in subcase:
        method_id = subcase.get_parameter('METHOD')[0]
        if method_id in fem2.methods:
            unused_method = fem2.methods[method_id]
        #elif method_id in fem2.cMethods:
            #method = fem2.cMethods[method_id]
        else:
            method_ids = list(fem2.methods.keys())
            raise RuntimeError('METHOD = %s not in method_ids=%s' % (method_id, method_ids))
        allowed_sols = [3, 5, 76, 100, 101, 103, 105, 106, 107, 108, 110, 111,
                        112, 144, 145, 146, 187, 200]
        ierror = check_sol(sol, subcase, allowed_sols, 'METHOD', log, ierror, nerrors)

    if 'CMETHOD' in subcase:
        cmethod_id = subcase.get_parameter('CMETHOD')[0]
        if cmethod_id in fem2.cMethods:
            unused_method = fem2.cMethods[cmethod_id]
        #elif method_id in fem2.cMethods:
            #unused_method = fem2.cMethods[method_id]
        else:
            cmethod_ids = list(fem2.cMethods.keys())
            raise RuntimeError('CMETHOD = %s not in cmethod_ids=%s' % (cmethod_id, cmethod_ids))
        allowed_sols = [107, 110, 111, 144, 145, 200]
        ierror = check_sol(sol, subcase, allowed_sols, 'CMETHOD', log, ierror, nerrors)

    if 'RMETHOD' in subcase:
        unused_rmethod_id = subcase.get_parameter('RMETHOD')[0]
        #if method_id in fem2.methods:
            #method = fem2.methods[method_id]
        #elif method_id in fem2.cMethods:
            #method = fem2.cMethods[method_id]
        #else:
            #method_ids = list(fem2.methods.keys())
            #raise RuntimeError('METHOD = %s not in method_ids=%s' % (method_id, method_ids))

        allowed_sols = [101, 110, 111]
        ierror = check_sol(sol, subcase, allowed_sols, 'RMETHOD', log, ierror, nerrors)

    if 'FMETHOD' in subcase:
        # FLUTTER
        fmethod_id = subcase.get_parameter('FMETHOD')[0]
        unused_fmethod = fem2.flutters[fmethod_id]
        allowed_sols = [145, 200]
        ierror = check_sol(sol, subcase, allowed_sols, 'FMETHOD', log, ierror, nerrors)

    nid_map = fem2.nid_map
    if 'TEMPERATURE(LOAD)' in subcase:
        loadcase_id = subcase.get_parameter('TEMPERATURE(LOAD)')[0]
        get_temperatures_array(fem2, loadcase_id, nid_map=nid_map, dtype='float32')
    if 'TEMPERATURE(BOTH)' in subcase:
        loadcase_id = subcase.get_parameter('TEMPERATURE(BOTH)')[0]
        get_temperatures_array(fem2, loadcase_id, nid_map=nid_map, dtype='float32')
    if 'TEMPERATURE(INITIAL)' in subcase:
        loadcase_id = subcase.get_parameter('TEMPERATURE(INITIAL)')[0]
        get_temperatures_array(fem2, loadcase_id, nid_map=nid_map, dtype='float32')
    if 'TEMPERATURE(MATERIAL)' in subcase:
        loadcase_id = subcase.get_parameter('TEMPERATURE(MATERIAL)')[0]
        get_temperatures_array(fem2, loadcase_id, nid_map=nid_map, dtype='float32')

    if 'LOAD' in subcase:
        cid_new = 0
        cid_msg = '' if cid_new == 0 else f'(cid={cid_new:d})'
        loadcase_id = subcase.get_parameter('LOAD')[0]
        force, moment = fem2.sum_forces_moments(p0, loadcase_id, cid=cid_new, include_grav=False)
        unused_fvec = get_static_force_vector_from_subcase_id(fem2, isubcase)
        eids = None
        nids = None
        force2, moment2 = fem2.sum_forces_moments_elements(
            p0, loadcase_id, eids, nids, cid=cid_new, include_grav=False)
        assert np.allclose(force, force2), 'force=%s force2=%s' % (force, force2)
        assert np.allclose(moment, moment2), 'moment=%s moment2=%s' % (moment, moment2)
        print('  isubcase=%i F=%s M=%s%s' % (isubcase, force, moment, cid_msg))
        allowed_sols = [
            1, 5, 24, 61, 64, 66, 100, 101, 103, 105, 106, 107,
            108, 109, 110, 111, 112, 114, 144, 145, 153, 200, 400, 401, 600, 601,
        ]
        ierror = check_sol(sol, subcase, allowed_sols, 'LOAD', log, ierror, nerrors)
    else:
        # print('is_load =', subcase.has_parameter('LOAD'))
        pass

    if 'FREQUENCY' in subcase:
        freq_id = subcase.get_parameter('FREQUENCY')[0]
        freq = fem2.frequencies[freq_id]
        allowed_sols = [26, 68, 76, 78, 88, 108, 101, 111, 112, 118, 146, 200]
        ierror = check_sol(sol, subcase, allowed_sols, 'FREQUENCY', log, ierror, nerrors)

    #if 'LSEQ' in subcase:
        #lseq_id = subcase.get_parameter('LSEQ')[0]
        #lseq = fem2.loads[lseq_id]
        #assert sol in [], sol
        #print(lseq)
    if 'SPC' in subcase:
        spc_id = subcase.get_parameter('SPC')[0]
        fem2.get_spcs(spc_id, stop_on_failure=False)
        fem2.get_reduced_spcs(spc_id, stop_on_failure=False)
        fem2.get_SPCx_node_ids_c1(spc_id, stop_on_failure=False)
        fem2.get_SPCx_node_ids(spc_id, stop_on_failure=False)

    if 'MPC' in subcase:
        mpc_id = subcase.get_parameter('MPC')[0]
        get_mpcs(fem2, mpc_id, consider_mpcadd=True, stop_on_failure=False)
        fem2.get_reduced_mpcs(mpc_id, consider_mpcadd=True, stop_on_failure=False)
        get_mpc_node_ids_c1(fem2, mpc_id, consider_mpcadd=True, stop_on_failure=False)
        get_mpc_node_ids(fem2, mpc_id, consider_mpcadd=True, stop_on_failure=False)

    if 'NSM' in subcase:
        nsm_id = subcase.get_parameter('NSM')[0]
        fem2.get_reduced_nsms(nsm_id, stop_on_failure=False)

    if 'SDAMPING' in subcase:
        sdamping_id = subcase.get_parameter('SDAMPING')[0]
        sdamp_sols = [110, 111, 112, 145, 146, 200]
        if not sdamping_id in fem2.tables_sdamping and fem2.sol in sdamp_sols:
            msg = 'SDAMPING = %s; not in TABDMP1, but must since its SOL %i\n' % (
                sdamping_id, fem2.sol)
            msg += 'TABDMP1 = %s\n' % list(fem2.tables_sdamping.keys())
            raise RuntimeError(msg)
        if not(sdamping_id in fem2.tables_sdamping or fem2.tables_d):
            msg = 'SDAMPING = %s; not in TABDMP1/TABLEDi\n' % sdamping_id
            msg += 'TABDMP1 = %s\n' % list(fem2.tables_sdamping.keys())
            msg += 'TABLEDi = %s\n' % list(fem2.tables_d.keys())
            raise RuntimeError(msg)

    if 'LOADSET' in subcase:
        loadset_id = subcase.get_parameter('LOADSET')[0]
        unused_lseq = fem2.Load(loadset_id)
        fem2.get_reduced_loads(loadset_id)

    if 'DLOAD' in subcase:
        allowed_sols = [
            26, 68, 76, 78, 88, 99, 103, 108, 109, 111, 112, 118, 129, 146,
            153, 159, 200, 400, 401, 601, 700,
        ]
        ierror = check_sol(sol, subcase, allowed_sols, 'DLOAD', log, ierror, nerrors)
        dload_id = subcase.get_parameter('DLOAD')[0]
        fem2.get_reduced_dloads(dload_id)
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
                    fem2.Load(ic_val)
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
        loads, scale_factors = fem2.get_reduced_dloads(dload_id)
        #fem2.log.info('scale_factors=%s loads=%s' % (scale_factors, loads))

        if sol in [108, 111]:  # direct frequency, modal frequency
            for load2, scale_factor in zip(loads, scale_factors):
                freq_id = subcase.get_parameter('FREQ')[0]
                freqs = fem2.frequencies[freq_id]
                for freq in freqs:
                    if freq.type in ['FREQ', 'FREQ1', 'FREQ2']:
                        fmax = freq.freqs[-1]
                        force = load2.get_load_at_freq(fmax) * scale_factor
        elif sol in [109, 129, 401]:  # direct transient (time linear), time nonlinear, NX nonlinear???
            for load2, scale_factor in zip(loads, scale_factors):
                force = load2.get_load_at_time(0.) * scale_factor
        elif  sol == 200:
            pass
        else:
            fem2.log.debug('solution=%s; DLOAD is not supported' % sol)

        # print(loads)
    return ierror


def divide(value1: int, value2: int) -> float:
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


def test_get_cards_by_card_types(model: BDF) -> None:
    """Verifies the ``model.get_cards_by_card_types`` method works"""
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
            if card_type != card.type and card_type + '1' != card.type:
                raise RuntimeError(msg)


def compare_card_count(fem1: BDF, fem2: BDF,
                       print_stats: bool=False, quiet: bool=False) -> List[str]:
    """Checks that no cards from fem1 are lost when we write fem2"""
    cards1 = fem1.card_count
    cards2 = fem2.card_count
    for key in cards1:
        if key != key.upper():
            raise RuntimeError('Proper capitalization wasnt determined')
    if print_stats and not quiet:
        print(fem1.get_bdf_stats())
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
        value1 = 0
        if key in list_keys1:
            value1 = cards1[key]

        value2 = 0
        if key in list_keys2:
            value2 = cards2[key]

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
    """Computes the difference between two dictionaries to data is the same"""
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


def get_element_stats(fem1: BDF, unused_fem2: BDF, quiet: bool=False) -> None:
    """verifies that the various element methods work"""
    for (unused_key, loads) in sorted(fem1.loads.items()):
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
    check_mass(fem1, quiet=quiet)

def check_mass(fem1: BDF, quiet: bool=False):
    mass1, cg1, inertia1 = mass_properties(fem1, reference_point=None, sym_axis=None)
    mass2, cg2, inertia2 = mass_properties_nsm(fem1, reference_point=None, sym_axis=None)
    #mass3, cg3, inertia3 = mass_properties_breakdown(fem1)[:3]
    if not quiet:
        if fem1.wtmass != 1.0:
            print('weight = %s' % (mass1 / fem1.wtmass))
        print(f'mass = {mass1}')
        print(f'cg   = {cg1}')
        print('Ixx=%s, Iyy=%s, Izz=%s \nIxy=%s, Ixz=%s, Iyz=%s' % tuple(inertia1))
    assert np.allclose(mass1, mass2), f'mass1={mass1} mass2={mass2}'
    assert np.allclose(cg1, cg2), f'mass={mass1}\ncg1={cg1} cg2={cg2}'
    assert np.allclose(inertia1, inertia2, atol=1e-5), f'mass={mass1} cg={cg1}\ninertia1={inertia1}\ninertia2={inertia2}\ndinertia={inertia1-inertia2}'

    for nsm_id in chain(fem1.nsms, fem1.nsmadds):
        mass, unused_cg, unused_inertia = mass_properties_nsm(
            fem1, reference_point=None, sym_axis=None, nsm_id=nsm_id)
        print('nsm_id=%s' % nsm_id)
        print('  mass = %s' % mass)
        print('  cg = %s' % cg1)
        print('  Ixx=%s, Iyy=%s, Izz=%s \n  Ixy=%s, Ixz=%s, Iyz=%s' % tuple(inertia1))

    reference_point = [10., 10., 10.]
    mass1, cg1, inertia1 = mass_properties(fem1, reference_point=reference_point, sym_axis=None)
    mass2, cg2, inertia2 = mass_properties_nsm(fem1, reference_point=reference_point, sym_axis=None)
    assert np.allclose(mass1, mass2), f'reference_point=[10., 10., 10.]; mass1={mass1} mass2={mass2}'
    assert np.allclose(cg1, cg2), f'reference_point=[10., 10., 10.]; mass={mass1} cg1={cg1} cg2={cg2}'
    assert np.allclose(inertia1, inertia2, atol=1e-5), f'reference_point=[10., 10., 10.]; mass={mass1} cg={cg1} inertia1={inertia1} inertia2={inertia2}'


def get_matrix_stats(fem1: BDF, unused_fem2: BDF) -> None:
    """Verifies the dmig.get_matrix() method works."""
    for (unused_key, dmig) in sorted(fem1.dmigs.items()):
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

    for (unused_key, dmi) in sorted(fem1.dmis.items()):
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

    for (unused_key, dmij) in sorted(fem1.dmij.items()):
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

    for (unused_key, dmiji) in sorted(fem1.dmijis.items()):
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

    for (unused_key, dmik) in sorted(fem1.dmiks.items()):
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


def test_bdf_argparse(argv=None):
    """test_bdf argument parser"""
    if argv is None:
        argv = sys.argv[1:]  # same as argparse
        #print('get_inputs; argv was None -> %s' % argv)
    else:
        # drop the pyNastranGUI; same as argparse
        argv = argv[1:]

    encoding = sys.getdefaultencoding()
    import argparse
    parent_parser = argparse.ArgumentParser()
    parent_parser.add_argument('BDF_FILENAME', help='path to BDF/DAT/NAS file',
                               type=str)
    parent_parser.add_argument('-v', '--version', action='version',
                               version=pyNastran.__version__)

    #nargs : str/int
    #   * : 0 or more
    #   + : one or more
    #   ? : optional
    #   int : int values
    # --------------------------------------------------------------------------
    # Options
    xref_safe_group = parent_parser.add_mutually_exclusive_group()
    xref_safe_group.add_argument(
        '-x', '--xref', action='store_false',
        help='disables cross-referencing and checks of the BDF (default=True -> on)')
    xref_safe_group.add_argument(
        '--safe', action='store_true',
        help='Use safe cross-reference (default=False)')

    parent_parser.add_argument(
        '-p', '--punch', action='store_true',
        help='disables reading the executive and case control decks in the BDF\n'
        '(default=False -> reads entire deck)')

    stop_check_group = parent_parser.add_mutually_exclusive_group()
    stop_check_group.add_argument(
        '-c', '--check', action='store_true',
        help='disables BDF checks.  Checks run the methods on \n'
        '                 every element/property to test them.  May fails if a \n'
        '                 card is fully not supported (default=False)')
    stop_check_group.add_argument('--stop', action='store_true', # dev
                                  help='Stop after first read/write (default=False)\n')

    width_group = parent_parser.add_mutually_exclusive_group()
    width_group.add_argument(
        '-l', '--large', action='store_true',
        help='writes the BDF in large field, single precision format (default=False)')
    width_group.add_argument(
        '-d', '--double', action='store_true',
        help='writes the BDF in large field, double precision format (default=False)')

    parent_parser.add_argument(
        '-L', '--loads', action='store_false',
        help='Disables forces/moments summation for the different subcases (default=True)')

    parent_parser.add_argument('-e', '--nerrors', default=100, type=int,
                               help='Allow for cross-reference errors (default=100)')
    parent_parser.add_argument('--encoding', default=encoding, type=str,
                               help='the encoding method (default=%r)\n' % encoding)
    parent_parser.add_argument('-q', '--quiet', action='store_true',
                               help='prints debug messages (default=False)')
    # --------------------------------------------------------------------------
    #'Developer:\n'
    parent_parser.add_argument('--crash', nargs=1, type=str,
                               help='Crash on specific cards (e.g. CGEN,EGRID)')

    parent_parser.add_argument('--dumplines', action='store_true',
                               help='Writes the BDF exactly as read with the INCLUDEs processed\n'
                               '(pyNastran_dump.bdf)')
    parent_parser.add_argument('--dictsort', action='store_true',
                               help='Writes the BDF exactly as read with the INCLUDEs processed\n'
                               '(pyNastran_dict.bdf)')
    parent_parser.add_argument('--profile', action='store_true',
                               help='Profiles the code (default=False)\n')
    parent_parser.add_argument('--pickle', action='store_true',
                               help='Pickles the data objects (default=False)\n')
    parent_parser.add_argument('--hdf5', action='store_true',
                               help='Save/load the BDF in HDF5 format')

    usage, args, examples = get_test_bdf_usage_args_examples(encoding)

    # --------------------------------------------------------------------------

    #argv
    #print(argv)
    from pyNastran.utils.arg_handling import argparse_to_dict, update_message # swap_key
    update_message(parent_parser, usage, args, examples)

    #try:
        #args = parent_parser.parse_args(args=argv)
    #except SystemExit:
        #fobj = StringIO()
        ##args = parent_parser.format_usage()
        #parent_parser.print_usage(file=fobj)
        #args = fobj.getvalue()
        #raise
    args = parent_parser.parse_args(args=argv)

    args2 = argparse_to_dict(args)
    #optional_args = [
        #'double', 'large', 'crash', 'quiet', 'profile',
        #'xref', 'safe', 'check', 'punch', 'loads', 'stop', 'encoding',
        #'dumplines', 'dictsort', 'nerrors', 'pickle', 'hdf5',
    #]
    #for arg in optional_args:
        #swap_key(args2, arg, '--' + arg)
    return args2
# defaults
#check        = False
#crash        = None
#dictsort     = False
#double       = False
#dumplines    = False
#encoding     = None
#hdf5         = False
#help         = False
#large        = False
#loads        = True
#nerrors      = 100
#pickle       = False
#profile      = False
#punch        = False
#quiet        = False
#safe         = False
#stop         = False
#version      = False
#xref         = True

def get_test_bdf_usage_args_examples(encoding):
    """helper method"""
    options = (
        '\n  [options] = [-e E] [--encoding ENCODE] [-q] [--dumplines] [--dictsort]\n'
        '              [--crash C] [--pickle] [--profile] [--hdf5]\n')
    usage = (
        "Usage:\n"
        '  test_bdf [-x | --safe] [-p] [-c] [-L]      BDF_FILENAME [options]\n'
        '  test_bdf [-x | --safe] [-p] [-c] [-L] [-d] BDF_FILENAME [options]\n'
        '  test_bdf [-x | --safe] [-p] [-c] [-L] [-l] BDF_FILENAME [options]\n'
        '  test_bdf               [-p]                BDF_FILENAME [options]\n'
        '  test_bdf [-x | --safe] [-p] [--stop]       BDF_FILENAME [options]\n'
        '  test_bdf -h | --help\n'
        '  test_bdf -v | --version\n' +
        options

        #"  test_bdf [-q] [-p] [-o [<VAR=VAL>]...] BDF_FILENAME\n"
    )
    args = (
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
        '  --crash C,   Crash on specific cards (e.g. CGEN,EGRID)\n'
        '  --stop       Stop after first read/write (default=False)\n'
        '  --dumplines  Writes the BDF exactly as read with the INCLUDEs processed\n'
        '               (pyNastran_dump.bdf)\n'
        '  --dictsort   Writes the BDF exactly as read with the INCLUDEs processed\n'
        '               (pyNastran_dict.bdf)\n'
        '  --profile    Profiles the code (default=False)\n'
        '  --pickle     Pickles the data objects (default=False)\n'
        '  --hdf5       Save/load the BDF in HDF5 format\n'
        '\n'
        'Info:\n'
        '  -h, --help     show this help message and exit\n'
        "  -v, --version  show program's version number and exit\n"
        '\n'
    )
    examples = (
        'Examples\n'
        '--------\n'
        '  test_bdf fem.bdf\n'
        '  test_bdf --xref fem.bdf\n'
    )
    return usage, args, examples


def main(argv=None):
    """The main function for the command line ``test_bdf`` script."""
    if argv is None:
        argv = sys.argv

    data = test_bdf_argparse(argv)
    for key, value in sorted(data.items()):
        print("%-12s = %r" % (key.strip('--'), value))

    import time
    time0 = time.time()

    is_double = False
    if data['double']:
        size = 16
        is_double = True
    elif data['large']:
        size = 16
    else:
        size = 8

    crash_cards = []
    if data['crash']:
        crash_cards = data['crash'].split(',')

    #print(data)
    debug = True
    if data['quiet']:
        debug = None
    if data['profile']:
        import pstats

        import cProfile
        prof = cProfile.Profile()
        prof.runcall(
            run_bdf,
            '.',
            data['BDF_FILENAME'],
            debug=debug,
            xref=data['xref'],
            check=not(data['check']),
            punch=data['punch'],
            size=size,
            is_double=is_double,
            sum_load=data['loads'],
            stop=data['stop'],
            quiet=data['quiet'],
            dumplines=data['dumplines'],
            dictsort=data['dictsort'],
            nerrors=data['nerrors'],
            encoding=data['encoding'],
            crash_cards=crash_cards,
            run_extract_bodies=False,
            pickle_obj=data['pickle'],
            safe_xref=data['safe'],
            hdf5=data['hdf5'],
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
            xref=data['xref'],
            # xref_safe=data['xref_safe'],
            check=not(data['check']),
            punch=data['punch'],
            size=size,
            is_double=is_double,
            sum_load=data['loads'],
            stop=data['stop'],
            quiet=data['quiet'],
            dumplines=data['dumplines'],
            dictsort=data['dictsort'],
            nerrors=data['nerrors'],
            encoding=data['encoding'],
            crash_cards=crash_cards,
            run_extract_bodies=False,
            pickle_obj=data['pickle'],
            safe_xref=data['safe'],
            hdf5=data['hdf5'],
            print_stats=True,
            stop_on_failure=False,
        )
    print("total time:  %.2f sec" % (time.time() - time0))


if __name__ == '__main__':  # pragma: no cover
    main()
