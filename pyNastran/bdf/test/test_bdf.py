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
from typing import Optional, Union, Any, cast
from io import StringIO

import numpy as np
from cpylog import get_logger2, SimpleLogger, WarningRedirector
#warnings.simplefilter('always')
warnings.simplefilter('default')

np.seterr(all='raise')

#from pyNastran.gui.qt_version import qt_version
#import PySide2
#import matplotlib
#matplotlib.use('Qt5Agg')

from pyNastran.utils import check_path, print_bad_path
from pyNastran.utils.numpy_utils import integer_types
from pyNastran.bdf.errors import (
    #CrossReferenceError,
    CardParseSyntaxError, DuplicateIDsError, MissingDeckSections,
    UnsupportedCard,
    DisabledCardError,
    ReplicationError,
    EnvironmentVariableError,
)
from pyNastran.bdf.subcase import Subcase
from pyNastran.bdf.bdf import (
    BDF, read_bdf, map_version, CaseControlDeck,
    FREQ1, FREQ2, FREQ3, FREQ4, FREQ5)
from pyNastran.bdf.mesh_utils.export_mcids import export_mcids, export_mcids_all
from pyNastran.bdf.mesh_utils.extract_bodies import extract_bodies
from pyNastran.bdf.mesh_utils.forces_moments import get_temperatures_array
from pyNastran.bdf.mesh_utils.mpc_dependency import (
    get_mpc_node_ids, get_mpc_node_ids_c1,
    get_dependent_nid_to_components, get_mpcs)
from pyNastran.bdf.mesh_utils.loads import (
    sum_forces_moments, sum_forces_moments_elements,
    get_static_force_vector_from_subcase_id)
from pyNastran.bdf.mesh_utils.skin_solid_elements import write_skin_solid_faces
from pyNastran.bdf.mesh_utils.export_caero_mesh import export_caero_mesh

from pyNastran.bdf.test.compare import compare
#from pyNastran.bdf.mesh_utils.convert import convert
#from pyNastran.bdf.mesh_utils.remove_unused import remove_unused

import pyNastran.bdf.test
TEST_PATH = pyNastran.bdf.test.__path__[0]
#warnings.filterwarnings("ignore", category=DeprecationWarning)

MESH_OPT_CARDS = [
    'GRIDG', 'CGEN', 'SPCG', 'FEEDGE', 'FEFACE', 'ADAPT', # 'EQUIV',
    'PVAL', 'GMCURV', 'GMSURF',
]
FREQS = Union[FREQ1, FREQ2, FREQ3, FREQ4, FREQ5]
class MeshOptimizationError(RuntimeError):
    pass

def run_lots_of_files(filenames: list[str], folder: str='',
                      debug: bool=False,
                      xref: bool=True,
                      check: bool=True,
                      write_hdf5: bool=True,
                      punch: bool=False,
                      nastran: str='',
                      encoding: Optional[str]=None,
                      size: Union[int, list[int], None]=None,
                      post: Union[int, list[int], None]=None,
                      is_double: Union[bool, list[bool], None]=None,
                      sum_load: bool=True,
                      run_mass: bool=True,
                      dev: bool=True,
                      crash_cards: Optional[list[str]]=None,
                      pickle_obj: bool=True, quiet: bool=False) -> list[str]:
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
    run_mass : bool; default=True
        should the mass be summed
    dev : bool; default=True
        True : crashes if an Exception occurs
        False : doesn't crash; useful for running many tests
    crash_cards : list[str, str, ...]
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
    diff_cards: list[str] = []
    for filename in filenames:
        if (filename.endswith(('.bdf', '.dat', '.nas')) and
                'pyNastran_crash' not in filename and
                'skin_file' not in filename):
            filenames2.append(filename)
    nfiles = len(filenames2)

    failed_files = []
    npass = 1
    nfailed = 1
    log = get_logger2(log=None, debug=debug, encoding='utf-8')
    with WarningRedirector(log) as unused_warn:
        for ifile, filename in enumerate(filenames2):
            abs_filename = os.path.abspath(os.path.join(folder, filename))
            if folder != '':
                print(f'filename = {abs_filename}')
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
                        post=posti,
                        sum_load=sum_load, run_mass=run_mass,
                        run_extract_bodies=False,

                        dev=dev,
                        crash_cards=crash_cards,
                        limit_mesh_opt=True,
                        pickle_obj=pickle_obj,
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
            except (ReplicationError, MeshOptimizationError) as e:
                if not dev:
                    raise
            except SystemExit:
                sys.exit('sys.exit...')
            except Exception:
                traceback.print_exc(file=sys.stdout)
                #raise
            print('-' * 80)

            if is_passed:
                sys.stderr.write(f'{npass:d}  {abs_filename}')
                npass += 1
            else:
                sys.stderr.write(f'*{nfailed:d} {abs_filename}')
                nfailed += 1
                failed_files.append(abs_filename)
            if ifile % 500 == 0:
                sys.stderr.write(f'\n{npass:d}/{nfiles:d}; nfailed={nfailed:d}')
            sys.stderr.write('\n')

    print('*' * 80)
    try:
        print("diff_cards1 = %s" % list(set(diff_cards)))
    except TypeError:
        print(f'diff_cards2 = {diff_cards}')
    return failed_files


def run_bdf(folder: str, bdf_filename: str,
            debug: bool=False, xref: bool=True, check: bool=True,
            punch: bool=False,
            mesh_form='separate', is_folder=False, print_stats=False,
            encoding=None,
            size=8, is_double=False,
            hdf5=False,
            stop=False, nastran='', post=-1, dynamic_vars=None,
            quiet=False, dumplines=False, dictsort=False,
            limit_mesh_opt: bool=False,
            sum_load: bool=True,
            run_mass: bool=True,
            run_extract_bodies: bool=False,
            run_skin_solids: bool=True,
            run_export_caero: bool=True,
            save_file_structure: bool=False,
            nerrors=0, dev: bool=False, crash_cards=None,
            safe_xref: bool=False, pickle_obj: bool=False,
            version: Optional[str]=None,
            validate_case_control: bool=True,
            stop_on_failure: bool=True, log=None, name: str=''):
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
    run_mass : bool; default=True
        Sum the mass
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
        version=version,
        limit_mesh_opt=limit_mesh_opt,
        run_extract_bodies=run_extract_bodies,
        run_skin_solids=run_skin_solids,
        run_export_caero=run_export_caero,
        run_mass=run_mass,
        save_file_structure=save_file_structure,
        pickle_obj=pickle_obj,
        validate_case_control=validate_case_control,
        stop_on_failure=stop_on_failure,
        log=log,
        name=name,
    )
    return fem1, fem2, diff_cards

def run_and_compare_fems(
        bdf_model: str,
        out_model: str,
        debug: bool=False,
        xref: bool=True,
        check: bool=True,
        punch: bool=False,
        mesh_form: str='combined',
        print_stats: bool=False,
        encoding=None,
        sum_load: bool=True,
        size: int=8,
        is_double: bool=False,
        save_file_structure: bool=False,
        stop: bool=False,
        nastran: str='',
        post: int=-1,
        hdf5: bool=False,
        dynamic_vars=None,
        quiet: bool=False,
        dumplines: bool=False,
        dictsort: bool=False,
        nerrors: int=0,
        dev: bool=False,
        crash_cards=None,
        version: Optional[str]=None,
        limit_mesh_opt: bool=False,
        safe_xref: bool=True,
        run_extract_bodies: bool=False,
        run_skin_solids: bool=True,
        run_export_caero: bool=True,
        run_mass: bool=True,
        pickle_obj: bool=False,
        validate_case_control: bool=True,
        stop_on_failure: bool=True,
        log: Optional[SimpleLogger]=None,
        name: str=''):
    """runs two fem models and compares them"""
    assert os.path.exists(bdf_model), f'{bdf_model!r} doesnt exist\n%s' % print_bad_path(bdf_model)
    fem1 = BDF(debug=debug, log=log)
    fem1.use_new_deck_parser = True
    if version:
        map_version(fem1, version)
    fem1.dumplines = dumplines

    fem1.set_error_storage(nparse_errors=nerrors, stop_on_parsing_error=True,
                           nxref_errors=nerrors, stop_on_xref_error=True)
    if dynamic_vars:
        fem1.set_dynamic_syntax(dynamic_vars)

    if not quiet:
        fem1.log.info('starting fem1 (read/write)')
    sys.stdout.flush()
    fem2 = None
    diff_cards = []

    #nastran_cmd = 'nastran scr=yes bat=no old=no news=no '
    nastran_cmd = ''
    try:
        #try:
        fem1.log.info('running fem1 (read/write)')
        fem1 = run_fem1(
            fem1, bdf_model, out_model, mesh_form, xref, punch, sum_load,
            size, is_double,
            run_extract_bodies=run_extract_bodies,
            run_skin_solids=run_skin_solids,
            run_export_caero=run_export_caero,
            save_file_structure=save_file_structure,
            hdf5=hdf5,
            encoding=encoding, crash_cards=crash_cards, safe_xref=safe_xref,
            limit_mesh_opt=limit_mesh_opt,
            pickle_obj=pickle_obj, stop=stop, name=name)

        if stop:
            if not quiet:
                print('card_count:')
                print('-----------')
                for card_name, card_count in sorted(fem1.card_count.items()):
                    print('key=%-8s value=%s' % (card_name, card_count))
            return fem1, None, None

        ierror = 0
        fem1.log.info('running fem2')
        fem2 = run_fem2(
            bdf_model, out_model, xref, punch,
            sum_load, size, is_double, mesh_form,
            safe_xref=safe_xref,
            encoding=encoding, debug=debug, quiet=quiet,
            ierror=ierror, nerrors=nerrors,
            stop_on_failure=stop_on_failure,
            validate_case_control=validate_case_control, log=log)

        diff_cards = compare(fem1, fem2, xref=xref,
                             run_mass=run_mass, check=check,
                             print_stats=print_stats, quiet=quiet)
        test_get_cards_by_card_types(fem2)

        fem2.update_model_by_desvars(xref)
        #except Exception:
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
        #elif is_mesh_opt:
            #print('failed test because mesh adaption (GRIDG,CGEN,SPCG)...ignoring')
            #print(e)
        else:
            print('failed test because DuplicateIDsError...ignoring')
    except MeshOptimizationError:
        print('failed test because mesh adaption (GRIDG,CGEN,SPCG)...ignoring')
    except DisabledCardError as e:
        if not dev:
            raise
    except EnvironmentVariableError:
        if not dev:
            raise
    #except RuntimeError as e:
        # only temporarily uncomment this when running lots of tests
        #if not dev:
            #raise
        #elif is_mesh_opt:
            #print('failed test because mesh adaption (GRIDG,CGEN,SPCG)...ignoring')
            #print(e)
        #else:
            #raise
    #except AttributeError:  # only temporarily uncomment this when running lots of tests
        #pass
    #except SyntaxError as e:
        # only temporarily uncomment this when running lots of tests
        #if not dev:
            #raise
        #elif is_mesh_opt:
            #print('failed test because mesh adaption (GRIDG,CGEN,SPCG)...ignoring')
            #print(e)
        #else:
            #raise
    #except KeyError as e:  # only temporarily uncomment this when running lots of tests
        #if not dev:
            #raise
        #elif is_mesh_opt:
            #print('failed test because mesh adaption (GRIDG,CGEN,SPCG)...ignoring')
            #print(e)
        #else:
            #raise
    #except AssertionError:  # only temporarily uncomment this when running lots of tests
        #pass
    except SystemExit:
        sys.exit('sys.exit...')
    except Exception:
        #exc_type, exc_value, exc_traceback = sys.exc_info()
        #print "\n"
        traceback.print_exc(file=sys.stdout)
        #print msg
        print("-" * 80)
        raise

    if not quiet:
        print("-" * 80)
    return (fem1, fem2, diff_cards)


def run_nastran(bdf_model: str, nastran: str, post: int=-1,
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

        f04_model = os.path.join(dirname, f'out_{basename}.f04')
        f06_model = os.path.join(dirname, f'out_{basename}.f06')
        op2_model = os.path.join(dirname, f'out_{basename}.f06')
        log_model = os.path.join(dirname, f'out_{basename}.log')
        xdb_model = os.path.join(dirname, f'out_{basename}.xdb')
        pch_model = os.path.join(dirname, f'out_{basename}.pch')
        asm_model = os.path.join(dirname, f'out_{basename}.asm')
        master_model = os.path.join(dirname, f'out_{basename}.master')
        #op2_model = os.path.join(dirname, 'out_%s.op2' % basename)

        #cwd = os.getcwd()
        cwd = dirname
        bdf_model2 = os.path.join(cwd, f'out_{basename}.bdf')
        op2_model2 = os.path.join(cwd, f'out_{basename}.op2')
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

def run_fem1(fem1: BDF, bdf_filename: str, out_model: str, mesh_form: str,
             xref: bool, punch: bool, sum_load: bool,
             size: int, is_double: bool,
             run_extract_bodies: bool=False, run_skin_solids: bool=True,
             run_export_caero: bool=True,
             save_file_structure: bool=False, hdf5: bool=False,
             encoding: Optional[str]=None,
             crash_cards: Optional[list[str]]=None,
             limit_mesh_opt: bool=False,
             safe_xref: bool=True, pickle_obj: bool=False, stop: bool=False,
             name: str='') -> BDF:
    """
    Reads/writes the BDF

    Parameters
    ----------
    fem1 : BDF()
        The BDF object
    bdf_filename : str
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
    log = fem1.log
    if crash_cards is None:
        crash_cards = []
    check_path(bdf_filename, 'bdf_filename')
    base, ext = os.path.splitext(out_model)
    try:
        if '.pch' in bdf_filename:
            fem1.read_bdf(bdf_filename, xref=False, punch=True, encoding=encoding,
                          save_file_structure=save_file_structure)
        else:
            fem1.read_bdf(bdf_filename, xref=False, punch=punch, encoding=encoding,
                          save_file_structure=save_file_structure)
            for card in crash_cards:
                if card in fem1.card_count:
                    raise DisabledCardError(f'card={card!r} has been disabled')
            #fem1.geom_check(geom_check=True, xref=False)

            if not stop and run_export_caero and len(fem1.caeros):
                caero_bdf_filename = base + '.caero.bdf'
                export_caero_mesh(fem1, caero_bdf_filename,
                                  is_subpanel_model=True,
                                  pid_method='caero', write_panel_xyz=True)

            if not stop and not xref and run_skin_solids:
                log.info('fem1-write_skin_solid_faces')
                skin_filename = 'skin_file.bdf'
                write_skin_solid_faces(fem1, skin_filename, size=16, is_double=False)
                if os.path.exists(skin_filename):
                    read_bdf(skin_filename, log=fem1.log)
                    os.remove(skin_filename)

            if limit_mesh_opt:
                limit_mesh_optimization(fem1)
            if xref:
                if run_extract_bodies:
                    log.info('fem1-run_extract_bodies')
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
                    log.debug('fem1.safe_cross_reference()')
                    fem1.safe_cross_reference()
                else:
                    log.debug('fem1.cross_reference()')
                    fem1.cross_reference()

                _fem_xref_methods_check(fem1)

                fem1._xref = True
                if fem1._nastran_format not in ['optistruct', 'mystran']:
                    log.info(f'fem1.bdf_filename = {fem1.bdf_filename}')
                    log.info('trying read_bdf from the raw filename')
                    read_bdf(fem1.bdf_filename, encoding=encoding, xref=False,
                             debug=fem1.debug, log=fem1.log)
                if safe_xref:
                    fem1.safe_cross_reference()
                elif xref:
                    fem1.cross_reference()

                fem1 = remake_model(bdf_filename, fem1, pickle_obj)
                #fem1.geom_check(geom_check=True, xref=True)
                #fem1.uncross_reference()
                #fem1.cross_reference()
    except Exception:
        print(f'failed reading {bdf_filename!r}')
        raise

    #out_model = bdf_filename + '_out'
    #if cid is not None and xref:
        #fem1.resolve_grids(cid=cid)

    if hdf5:
        hdf5_filename = f'{out_model}{name}.h5'
        _test_hdf5(fem1, hdf5_filename)

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
        #for pid, prop in sorted(fem1.properties.items()):
        #    if prop.type == 'PBARL':
        #        #A, I1, I2, I12 = prop.A_I1_I2_I12()
        #        log.debug(f'  pid={pid:d} type={prop.type!r} A={prop.Area():g} '
        #                  f'I1={prop.I1():g} I2={prop.I2():g} I12={prop.I12():g}')

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

def limit_mesh_optimization(model: BDF):
    is_mesh_opt = [card_name in model.card_count for card_name in MESH_OPT_CARDS]
    mesh_opt_cards = [card_name for is_mesh_opti, card_name in zip(is_mesh_opt, MESH_OPT_CARDS)
                      if is_mesh_opti]
    _cards = ', '.join(mesh_opt_cards)
    if any(is_mesh_opt):
        raise MeshOptimizationError(f'model contains [{_cards}]; mesh optimization is not supported')

def _test_hdf5(fem1: BDF, hdf5_filename: str) -> None:
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
            hdf5_msg += f'key={key!r} was not loaded to hdf5\n'

        if hdf5_msg:
            hdf5_msg += 'expected=%s\nactual=%s' % (
                fem1.card_count, fem1a.card_count)
            fem1a.log.error(hdf5_msg)
            raise RuntimeError(hdf5_msg)
    #sys.exit('hdf5')

def _fem_xref_methods_check(fem1: BDF) -> None:
    """
    testing that these methods work with xref
    """
    log = fem1.log
    log.debug('_fem_xref_methods_check(fem1)')

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

def remake_model(bdf_model: str, fem1: BDF, pickle_obj: bool) -> BDF:
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
             stop_on_failure: bool=True, validate_case_control: bool=True,
             ierror: int=0, nerrors: int=100, log=None) -> BDF:
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
        suppress prints

    """
    assert os.path.exists(bdf_model), bdf_model
    assert os.path.exists(out_model), out_model
    assert isinstance(ierror, int), ierror

    fem2 = BDF(debug=debug, log=log)
    fem2.use_new_deck_parser = True
    if not quiet:
        fem2.log.info('starting fem2')
    sys.stdout.flush()
    try:
        fem2.read_bdf(out_model, xref=False, punch=punch, encoding=encoding)
        if safe_xref:
            fem2.safe_cross_reference()
        elif xref:
            fem2.cross_reference()
    except Exception:
        print(f'failed reading {out_model!r}')
        raise

    out_model_2 = f'{bdf_model}_out2'

    if xref:
        if 'POST' in fem2.params:
            value = fem2.params['POST'].values[0]
            if value >= 0:
                msg = f'PARAM,POST,{value:d} is not supported by the OP2 reader'
                fem2.log.warning(msg)
        else:
            msg = 'PARAM,POST,0 is not supported by the OP2 reader'
            fem2.log.warning(msg)

        p0 = np.array([0., 0., 0.])

        case_control: CaseControlDeck = fem2.case_control_deck
        subcase_keys = case_control.get_subcase_list()
        subcases = fem2.subcases

        sol_200_map = case_control.sol_200_map
        sol_base = fem2.sol
        is_restart = _has_restart(fem2)
        if validate_case_control and not is_restart:
            assert isinstance(ierror, int), ierror
            ierror = _validate_case_control(
                fem2, p0, sol_base, subcase_keys, subcases, sol_200_map,
                sum_load=sum_load,
                ierror=ierror, nerrors=nerrors,
                stop_on_failure=stop_on_failure)

    if mesh_form is not None:
        fem2.write_bdf(out_model_2, interspersed=False,
                       size=size, is_double=is_double)
        try:
            os.remove(out_model_2)
        except PermissionError:  # pragma: no cover
            fem2.log.warning(f'cannot remove {out_model_2} due to a permissions error')
    #fem2.write_as_ctria3(out_model_2)
    return fem2

def _has_restart(fem: BDF):
    is_restart = False
    for line in fem.system_command_lines:
        if line.strip().upper().startswith('RESTART'):
            is_restart = True
    return is_restart

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

def _validate_case_control(fem: BDF, p0: Any, sol_base: int,
                           subcase_keys: list[int],
                           subcases: Any, unused_sol_200_map: dict[str, str],
                           sum_load: bool=True,
                           stop_on_failure: bool=True,
                           ierror: int=0, nerrors: int=100) -> int:
    if len(subcase_keys) > 1:
        subcase_keys = subcase_keys[1:]  # drop isubcase = 0

    for isubcase in subcase_keys:
        subcase = subcases[isubcase]
        if len(subcase.params) == 0:
            fem.log.error(f'Subcase {isubcase:d} is empty')
            continue
        str(subcase)

        if fem.is_acoustic():
            pass
        elif sol_base is None:
            raise RuntimeError(f'subcase: {subcase}\n')
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
            sol_base, subcase, fem, p0, isubcase, subcases,
            sum_load=sum_load,
            ierror=ierror, nerrors=nerrors,
            stop_on_failure=stop_on_failure)
    return ierror

def check_for_flag_in_subcases(fem2: BDF, subcase: Any, parameters: list[str]) -> None:
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
            msg = f'sol={fem2.sol!r}; {str(parameters)} not in subcase\n'
            for unused_isubcase, subcasei in fem2.subcases.items():
                msg += str(subcasei)
            raise RuntimeError(msg)

def stop_if_max_error(msg: str, error: Any, ierror: int, nerrors: int) -> int:
    """if the error count is greater than nerrors, stop"""
    if ierror == nerrors:
        raise error(msg)
    ierror += 1
    return ierror

def check_for_optional_param(keys: list[str], subcase: Any,
                             msg: str, error: Any, log: Any,
                             ierror: int, nerrors: int) -> int:
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
              allowed_sols: list[int],
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
        msg = f'sol={sol} is missing {case_control_key!r}\n{subcase}'
        log.error(msg)
        if ierror == nerrors:
            raise RuntimeError(msg)
        ierror += 1
    return ierror

def check_subcase_dmig_matrices(fem: BDF, subcase: Subcase) -> None:
    """verifies that all the DMIG matrices exist"""
     # acoustic/fluid-structure coupling matrix
    check_subcase_dmig_matrix(fem, subcase, 'A2GG')

    # mass matrices
    check_subcase_dmig_matrix(fem, subcase, 'M2GG')
    # mass matrices, which are not included in normal modes.
    check_subcase_dmig_matrix(fem, subcase, 'M2PP', is_real=False)

    # stiffness matrices
    check_subcase_dmig_matrix(fem, subcase, 'K2GG')
    # stiffness matrices, which are not included in normal modes
    check_subcase_dmig_matrix(fem, subcase, 'K2PP', is_real=False)
    # structural damping matrices
    check_subcase_dmig_matrix(fem, subcase, 'K42GG')

    # damping matrices
    check_subcase_dmig_matrix(fem, subcase, 'B2GG')
    check_subcase_dmig_matrix(fem, subcase, 'B2PP')

    # load matrices
    check_subcase_dmig_matrix(fem, subcase, 'P2G')

def check_subcase_dmig_matrix(fem: BDF,
                              subcase: Subcase,
                              matrix_name: str,
                              is_real: bool=True) -> None:
    """
    K2GG=KDMIG1, KDMIG2, KDMIG3
    K2GG=1.25*KDMIG1, 1.0*KDMIG2, 0.75*KDMIG3
    """
    if not any(subcase.has_parameter(matrix_name)):
        return

    #([(1.0, 'MCB')], [(1.0, 'MCB')])
    #[(1.0, 1.0, 'UMASS')]
    scale_names: list[tuple[float, str]] = subcase.get_parameter(matrix_name)[0]
    print(f'{matrix_name} (scale,names) = {scale_names}')

    print(f'{matrix_name}_name')
    if is_real:
        for scale, name in scale_names:
            dmig = fem.dmigs[name]
    else:
        for scale1, scale2, name in scale_names:
            dmig = fem.dmigs[name]

    del dmig

def check_case(sol: int,
               subcase: Subcase,
               fem2: BDF,
               p0: np.ndarray,
               isubcase: int,
               subcases: dict[int, Subcase],
               sum_load:bool =True,
               ierror: int=0, nerrors: int=100,
               stop_on_failure: bool=True) -> int:
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
       AERO

       # one or both
       SUPORT1 = 5 # implicit point to SUPORT1
       # implicit call to SUPORT

    """
    log = fem2.log

    msg = f'sol={sol}\n{subcase}'

    if fem2.is_acoustic():
        pass
    elif sol == 24:
        _assert_has_spc(subcase, fem2)
        assert True in subcase.has_parameter('LOAD'), msg
    elif sol == 64:
        #assert 'NLPARM' in subcase, subcase
        #_assert_has_spc(subcase, fem2)
        assert True in subcase.has_parameter('LOAD'), msg
    elif sol == 66:
        assert 'NLPARM' in subcase, subcase
        _assert_has_spc(subcase, fem2)
        ierror = check_for_optional_param(('LOAD', 'TEMPERATURE(LOAD)'), subcase, msg,
                                          RuntimeError, log, ierror, nerrors)
    elif sol == 99:
        assert 'DLOAD' in subcase, subcase
        assert 'LOADSET' in subcase, subcase
        _assert_has_spc(subcase, fem2)
        #assert True in subcase.has_parameter('LOAD', 'TEMPERATURE'), 'sol=%s\n%s' % (sol, subcase)
        ierror = check_for_optional_param(('TSTEP', 'TSTEPNL'),
                                          subcase, msg,
                                          RuntimeError, log, ierror, nerrors)
    elif sol in {1, 101, 'SESTATIC', 'SESTATICS'}:
        _assert_has_spc(subcase, fem2)
        ierror = check_for_optional_param(
            ('LOAD', 'TEMPERATURE(LOAD)', 'P2G'),
            subcase, msg, RuntimeError, log, ierror, nerrors)
    elif sol in {3, 103, 'SEMODES'}:
        ierror = check_for_optional_param(
            ('METHOD', 'RSMETHOD', 'RIGID', 'BOLTID', 'BGSET'),
            subcase, msg, RuntimeError, log, ierror, nerrors)
    elif sol in {5, 105, 'BUCK', 'BUCKLING'}: # buckling
        _assert_has_spc(subcase, fem2)
        ierror = check_for_optional_param(
            ('LOAD', 'TEMPERATURE(LOAD)', 'METHOD'),
            subcase, msg, RuntimeError, log, ierror, nerrors)
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
    elif sol in {106, 'NLSTATIC', 'NLSTATICS'}: # freq
        assert 'NLPARM' in subcase, subcase
        ierror = check_for_optional_param(
            ('LOAD', 'TEMPERATURE(LOAD)', 'CLOAD'),
            subcase, msg, RuntimeError, log, ierror, nerrors)
    elif sol in {7, 107, 'SEDCEIG'}:
        # direct complex eigenvalue
        _assert_has_spc(subcase, fem2)
        ierror = check_for_optional_param(
            ('CMETHOD', ),
            subcase, msg, RuntimeError, log, ierror, nerrors)
        #ierror = check_for_optional_param(('LOAD', 'TEMPERATURE(LOAD)'), subcase, msg,
                                          #RuntimeError, log, ierror, nerrors)
    elif sol in {8, 108}: # freq
        assert 'FREQUENCY' in subcase, subcase
    elif sol in {109, 'SEDTRAN'}:  # time
        check_for_flag_in_subcases(fem2, subcase, ('TIME', 'TSTEP', 'TSTEPNL'))

    elif sol in {110, 'SEMCEIG'}:  # modal complex eigenvalues
        _assert_has_spc(subcase, fem2)
        ierror = check_for_optional_param(
            ('LOAD', 'STATSUB'),
            subcase, msg, RuntimeError, log, ierror, nerrors)
    elif sol in {11, 111, 'SEMFREQ'}:  # modal frequency
        assert subcase.has_parameter('FREQUENCY'), msg
        assert any(subcase.has_parameter('METHOD', 'RMETHOD')), msg
    elif sol in {112, 'SEMTRAN'}:  # modal transient
        check_for_flag_in_subcases(fem2, subcase, ('TIME', 'TSTEP', 'TSTEPNL'))
        #assert any(subcase.has_parameter('TIME', 'TSTEP', 'TSTEPNL')), 'sol=%s\n%s' % (sol, subcase)
    elif sol == 114:
        soltype = 'CYCSTATX'
        ierror = require_cards(['LOAD', 'HARMONICS'], log, soltype, sol, subcase,
                               RuntimeError, ierror, nerrors)
        _assert_has_spc(subcase, fem2)
    elif sol == 115:  # cyclic modes
        assert any(subcase.has_parameter('METHOD')), msg
    elif sol == 118:
        soltype = 'CYCFREQ'

        ierror = require_either_cards(['LOAD', 'DLOAD'],
                                      log, soltype, sol, subcase,
                                      RuntimeError, ierror, nerrors)
        ierror = require_cards(['HARMONICS', 'SDAMPING', 'FREQUENCY'],
                               log, soltype, sol, subcase,
                               RuntimeError, ierror, nerrors)
        _assert_has_spc(subcase, fem2)

    elif sol in {29, 129, 'NLTRAN'}:  # nonlinear transient
        assert any(subcase.has_parameter('TIME', 'TSTEP', 'TSTEPNL')), msg
    elif sol in {159, 'NLTCSH'}:  # thermal transient
        assert any(subcase.has_parameter('TIME', 'TSTEP', 'TSTEPNL')), msg

    elif sol == 144:  # aero - trim, diverg
        ierror = _check_static_aero_case(fem2, log, sol, subcase, ierror, nerrors)
    elif sol == 145:  # aero - flutter
        ierror = _check_flutter_case(fem2, log, sol, subcase, ierror, nerrors)
    elif sol == 146:  # aero - gust
        ierror = _check_gust_case(fem2, log, sol, subcase, ierror, nerrors)

    elif sol in {153, 'NLSCSH'}:
        # Static structural and/or steady state heat Transfer analysis with options:
        # Linear or nonlinear analysis.
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
                f'- SOL %{sol:d}\n{subcase}')
            log.error(msg)
            ierror = stop_if_max_error(msg, RuntimeError, ierror, nerrors)

        #assert any(subcase.has_parameter('TIME', 'TSTEP', 'TSTEPNL')), subcase
        #assert any(subcase.has_parameter('GUST', 'LOAD')), subcase
        if 'ANALYSIS' in subcase and subcase.get_parameter('ANALYSIS')[0] == 'HEAT':
            assert any(subcase.has_parameter('TEMPERATURE(LOAD)', 'TEMPERATURE(INITIAL)')), msg

    elif sol in {200, 'DESOPT'}:
        _check_case_sol_200(sol, subcase, fem2, p0, isubcase, subcases, log)
    elif sol in [114, 116, 118]:
        # cyclic statics, buckling, frequency
        pass
    elif sol in {401, 'NLSTEP'}:  # nx
        pass
    elif sol in {402, 'NLSTPKIN'}:  # nx
        pass
    elif sol in {187, 'RESDDAM'}:  ## DDAM
        pass
    elif sol in {21, 26, 27, 28, 30, 31, 38, 39,
                 47, 48, 61, 63, 67, 68,
                 75, 76, 78, 81, 83, 88, 91,
                 100, 128, 190,
                 400, 600, 601, 700, 701,
                 'AEDB2XDB', 'BUILDIT', 'CATALOG', 'CHKCPY', 'DBDIR',
                 'FTGRSTRT', 'FUNCTEST', 'INPUTT2',
                 'MAIN', 'MXDMAP50', 'MYDMAP', 'MYSOL',
                 'PARAMS', 'PRINT68', 'TABTSTB', 'TSTGINO', 'UPWARD', 'USERDMAP',
                 'U24'}:
        pass
    else:
        msg = f'SOL = {sol!r}\n'
        msg += str(subcase)
        raise NotImplementedError(msg)
    assert isinstance(ierror, int)
    ierror = _check_case_parameters(
        subcase, fem2, p0, isubcase, sol,
        sum_load=sum_load,
        ierror=ierror, nerrors=nerrors,
        stop_on_failure=stop_on_failure)
    assert isinstance(ierror, int)
    return ierror

def _check_static_aero_case(fem2: BDF, log: SimpleLogger, sol: int,
                            subcase: Subcase, ierror: int, nerrors: int) -> int:
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

def _check_flutter_case(fem2: BDF, log: SimpleLogger, sol: int, subcase: Subcase,
                        ierror: int, nerrors: int) -> int:
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

def _check_gust_case(fem2: BDF, log: SimpleLogger, sol: int, subcase: Subcase,
                     ierror: int, nerrors: int) -> int:
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
                        isubcase: int,
                        subcases: dict[int, Subcase],
                        log: SimpleLogger) -> None:
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
    log = fem2.log
    log.info('checking SOL 200...')
    assert 'ANALYSIS' in subcase, 'sol=%s\n%s' % (sol, subcase)

    analysis = subcase.get_str_parameter('ANALYSIS')[0]
    if analysis.startswith('SE'):
        analysis = analysis[2:]

    # BUCKLING
    if 'DESOBJ' in subcase:
        value = subcase.get_int_parameter('DESOBJ')[0]
        assert value in fem2.dresps, f'value={value} not in dresps'
    else:
        log.warning('no DESOBJ (DRESPi) in this subcase; '
                    'is this a buckling preload case?')
        log.warning('\n%s' % subcase)

    nopt = len(fem2.dvprels) + len(fem2.dvmrels) + len(fem2.dvcrels)
    if nopt == 0:
        log.error('no DVPRELs/DVMRELs/DVCRELs found')

    #--------------------------------------------------------------------------
    # DCONSTR
    if 'DESSUB' not in subcase and 'DESGLB' not in subcase:
        log.warning('no DESSUB/DESGLB (DCONSTR) in this subcase;'
                    ' is this a buckling preload case?')
        log.warning('\n%s' % subcase)

    if 'DESSUB' in subcase:
        value = subcase.get_int_parameter('DESSUB')[0]
        if value not in fem2.dconstrs:
            msg = 'value=%s not in dconstrs; Allowed DCONSTRs=%s' % (
                value, np.unique(list(fem2.dconstrs.keys())))
            raise RuntimeError(msg)
    if 'DESGLB' in subcase:
        value = subcase.get_int_parameter('DESGLB')[0]
        if value not in fem2.dconstrs:
            msg = 'value=%s not in dconstrs; Allowed DCONSTRs=%s' % (
                value, np.unique(list(fem2.dconstrs.keys())))
            raise RuntimeError(msg)
    #--------------------------------------------------------------------------

    if analysis in {'STATIC', 'STATICS'}:
        solution = 101
        check_case(solution, subcase, fem2, p0, isubcase, subcases)
    elif analysis in {'MODE', 'MODES'}:
        solution = 103
        check_case(solution, subcase, fem2, p0, isubcase, subcases)
    elif analysis in {'BUCK', 'BUCKLING'}:
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
    elif analysis in {'MTRAN', 'MTRANS'}:
        solution = 112
        check_case(solution, subcase, fem2, p0, isubcase, subcases)
    elif analysis in {'SAERO', 'DIVERG', 'DIVERGE'}:
        solution = 144
        check_case(solution, subcase, fem2, p0, isubcase, subcases)
    elif analysis in {'FLUT', 'FLUTTER', 'FLUTTR'}:
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

def require_cards(card_names: list[str], log: SimpleLogger,
                  soltype: str, sol: int, subcase: Subcase,
                  error, ierror: int, nerrors: int) -> int:
    """all must be True"""
    for card_name in card_names:
        if card_name not in subcase:
            msg = f'A {card_name} card is required for {soltype} - SOL {sol:d}\n{subcase}'
            log.error(msg)
            if ierror == nerrors:
                raise error(msg)
            ierror += 1
    return ierror

def require_either_cards(card_names: list[str], log: SimpleLogger,
                         soltype: str, sol: int, subcase: Subcase,
                         error, ierror: int, nerrors: int) -> int:
    """one or more must be True"""
    msg = ''
    nlocal_errors = 0
    for card_name in card_names:
        if card_name not in subcase:
            msg += f'A {card_name} card is required for {soltype} - SOL {sol:d}\n{subcase}'
            nlocal_errors += 1
    if nlocal_errors == len(card_names):
        ierror += 1
        log.error(msg)
        if ierror == nerrors:
            raise error(msg)
    return ierror

def _tstep_msg(fem: BDF,
               subcase: Subcase,
               tstep_id: int,
               tstep_type: str='') -> str:
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


def _check_case_parameters(subcase: Subcase,
                           fem: BDF,
                           p0: np.ndarray,
                           isubcase: int, sol: int,
                           sum_load: bool=True,
                           ierror: int=0, nerrors: int=100,
                           stop_on_failure: bool=True) -> int:
    """helper method for ``check_case``"""
    log = fem.log
    if fem.sol in {401, 402}:
        # TSTEP references a TSTEP1, but not a TSTEP
        # TSTEP1s are stored in tstepnls
        if any(subcase.has_parameter('TIME', 'TSTEP')):
            if 'TSTEP' in subcase:
                tstep_id = subcase.get_int_parameter('TSTEP')[0]
            else:  # pragma: no cover
                raise NotImplementedError(subcase)
            if tstep_id not in fem.tstepnls:
                raise RuntimeError(_tstep_msg(fem, subcase, tstep_id))
        else:
            raise RuntimeError(f'missing TSTEP in case control deck\n{subcase}')
    else:
        if any(subcase.has_parameter('TIME', 'TSTEP')):
            if 'TIME' in subcase:
                tstep_id = subcase.get_int_parameter('TIME')[0]
            elif 'TSTEP' in subcase:
                tstep_id = subcase.get_int_parameter('TSTEP')[0]
            else:  # pragma: no cover
                raise NotImplementedError(subcase)
            if tstep_id not in fem.tsteps:
                raise RuntimeError(_tstep_msg(fem, subcase, tstep_id))

    check_subcase_dmig_matrices(fem, subcase)
    if 'TSTEPNL' in subcase:
        tstepnl_id = subcase.get_int_parameter('TSTEPNL')[0]
        assert tstepnl_id in fem.tstepnls, _tstep_msg(fem, subcase, tstepnl_id, tstep_type='nl')

    if 'SUPORT1' in subcase:
        suport1_id = subcase.get_int_parameter('SUPORT1')[0]
        assert suport1_id in fem.suport1, f'suport1_id={suport1_id}\n suport1={fem.suport1}\n subcase:\n{subcase}'

    ierror = _check_case_parameters_aero(
        subcase, fem, sol,
        ierror=ierror, nerrors=nerrors, stop_on_failure=stop_on_failure)

    if 'METHOD' in subcase: # or 'CMETHOD' in subcase:
        method_id = subcase.get_int_parameter('METHOD')[0]
        if method_id in fem.methods:
            unused_method = fem.methods[method_id]
        #elif method_id in fem.cMethods:
            #method = fem.cMethods[method_id]
        else:
            method_ids = list(fem.methods.keys())
            raise RuntimeError('METHOD = %s not in method_ids=%s' % (method_id, method_ids))
        allowed_sols = [3, 5, 76, 100, 101, 103, 105, 106, 107, 108, 110, 111,
                        112, 115, 144, 145, 146, 187, 200, 400]
        ierror = check_sol(sol, subcase, allowed_sols, 'METHOD', log, ierror, nerrors)

    if 'CMETHOD' in subcase:
        cmethod_id = subcase.get_int_parameter('CMETHOD')[0]
        if cmethod_id in fem.cMethods:
            unused_method = fem.cMethods[cmethod_id]
        #elif method_id in fem.cMethods:
            #unused_method = fem.cMethods[method_id]
        else:
            cmethod_ids = list(fem.cMethods.keys())
            raise RuntimeError('CMETHOD = %s not in cmethod_ids=%s' % (cmethod_id, cmethod_ids))
        allowed_sols = [107, 110, 111, 144, 145, 200, 400]
        ierror = check_sol(sol, subcase, allowed_sols, 'CMETHOD', log, ierror, nerrors)

    if 'RMETHOD' in subcase:
        unused_rmethod_id = subcase.get_int_parameter('RMETHOD')[0]
        #if method_id in fem.methods:
            #method = fem.methods[method_id]
        #elif method_id in fem.cMethods:
            #method = fem.cMethods[method_id]
        #else:
            #method_ids = list(fem.methods.keys())
            #raise RuntimeError('METHOD = %s not in method_ids=%s' % (method_id, method_ids))

        allowed_sols = [101, 110, 111]
        ierror = check_sol(sol, subcase, allowed_sols, 'RMETHOD', log, ierror, nerrors)

    nid_map = fem.nid_map
    if 'TEMPERATURE(LOAD)' in subcase:
        loadcase_id = subcase.get_int_parameter('TEMPERATURE(LOAD)')[0]
        get_temperatures_array(fem, loadcase_id, nid_map=nid_map, fdtype='float32')
    if 'TEMPERATURE(BOTH)' in subcase:
        loadcase_id = subcase.get_int_parameter('TEMPERATURE(BOTH)')[0]
        get_temperatures_array(fem, loadcase_id, nid_map=nid_map, fdtype='float32')
    if 'TEMPERATURE(INITIAL)' in subcase:
        loadcase_id = subcase.get_int_parameter('TEMPERATURE(INITIAL)')[0]
        get_temperatures_array(fem, loadcase_id, nid_map=nid_map, fdtype='float32')
    if 'TEMPERATURE(MATERIAL)' in subcase:
        loadcase_id = subcase.get_int_parameter('TEMPERATURE(MATERIAL)')[0]
        get_temperatures_array(fem, loadcase_id, nid_map=nid_map, fdtype='float32')

    is_not_acoustic = not fem.is_acoustic()
    if 'LOAD' in subcase and is_not_acoustic:
        cid_new = 0
        cid_msg = '' if cid_new == 0 else f'(cid={cid_new:d})'
        loadcase_id = subcase.get_parameter('LOAD')[0]
        if sum_load:
            force, moment = sum_forces_moments(
                fem, p0, loadcase_id,
                cid=cid_new, include_grav=False)
            unused_fvec = get_static_force_vector_from_subcase_id(fem, isubcase)
            eids = None
            nids = None
            force2, moment2 = sum_forces_moments_elements(
                fem, p0, loadcase_id, eids, nids,
                cid=cid_new, include_grav=False)
            assert np.allclose(force, force2), 'force=%s force2=%s' % (force, force2)
            assert np.allclose(moment, moment2), 'moment=%s moment2=%s' % (moment, moment2)
            print('  isubcase=%i F=%s M=%s%s' % (isubcase, force, moment, cid_msg))
        allowed_sols = [
            1, 5, 24, 38, 61, 64, 66, 100, 101, 103, 105, 106, 107,
            108, 109, 110, 111, 112, 114, 144, 145, 153, 200, 400, 401, 600, 601,
            700,
        ]
        ierror = check_sol(sol, subcase, allowed_sols, 'LOAD', log, ierror, nerrors)
    else:
        # print('is_load =', subcase.has_parameter('LOAD'))
        pass

    if 'FREQUENCY' in subcase and is_not_acoustic:
        # not positive on 107 - complex eigenvalues
        freq_id = subcase.get_int_parameter('FREQUENCY')[0]
        freq = fem.frequencies[freq_id]
        allowed_sols = [26, 68, 76, 78, 88, 108, 101, 107, 111, 112, 118, 146, 200]
        ierror = check_sol(sol, subcase, allowed_sols, 'FREQUENCY', log, ierror, nerrors)

    #if 'LSEQ' in subcase:
        #lseq_id = subcase.get_int_parameter('LSEQ')[0]
        #lseq = fem.loads[lseq_id]
        #assert sol in [], sol
        #print(lseq)
    if 'SPC' in subcase:
        spc_id: int = subcase.get_int_parameter('SPC')[0]
        fem.get_spcs(spc_id, stop_on_failure=False)
        fem.get_reduced_spcs(spc_id, stop_on_failure=False)
        fem.get_SPCx_node_ids_c1(spc_id, stop_on_failure=False)
        fem.get_SPCx_node_ids(spc_id, stop_on_failure=False)

    if 'MPC' in subcase:
        mpc_id: int = subcase.get_int_parameter('MPC')[0]
        get_mpcs(fem, mpc_id, consider_mpcadd=True, stop_on_failure=False)
        fem.get_reduced_mpcs(mpc_id, consider_mpcadd=True, stop_on_failure=False)
        get_mpc_node_ids_c1(fem, mpc_id, consider_mpcadd=True, stop_on_failure=False)
        get_mpc_node_ids(fem, mpc_id, consider_mpcadd=True, stop_on_failure=False)

    if 'NSM' in subcase:
        nsm_id = subcase.get_int_parameter('NSM')[0]
        fem.get_reduced_nsms(nsm_id, stop_on_failure=False)

    if 'SDAMPING' in subcase:
        sdamping_id = subcase.get_int_parameter('SDAMPING')[0]
        sdamp_sols = [110, 111, 112, 145, 146, 200]
        if not sdamping_id in fem.tables_sdamping and fem.sol in sdamp_sols:
            msg = 'SDAMPING = %s; not in TABDMP1, but must since its SOL %i\n' % (
                sdamping_id, fem.sol)
            msg += 'TABDMP1 = %s\n' % list(fem.tables_sdamping.keys())
            raise RuntimeError(msg)
        if not(sdamping_id in fem.tables_sdamping or fem.tables_d):
            msg = 'SDAMPING = %s; not in TABDMP1/TABLEDi\n' % sdamping_id
            msg += 'TABDMP1 = %s\n' % list(fem.tables_sdamping.keys())
            msg += 'TABLEDi = %s\n' % list(fem.tables_d.keys())
            raise RuntimeError(msg)

    if 'LOADSET' in subcase:
        loadset_id = subcase.get_int_parameter('LOADSET')[0]
        unused_lseq = fem.Load(loadset_id)
        fem.get_reduced_loads(loadset_id)

    if 'DLOAD' in subcase and is_not_acoustic:
        allowed_sols = [
            26, 68, 76, 78, 88, 99, 103, 108, 109, 111, 112, 118, 129, 146,
            153, 159, 200, 400, 401, 601, 700,
        ]
        ierror = check_sol(sol, subcase, allowed_sols, 'DLOAD', log, ierror, nerrors)
        dload_id = subcase.get_int_parameter('DLOAD')[0]
        fem.get_reduced_dloads(dload_id)
        #if 'LOADSET' in subcase:
            #raise NotImplementedError('LOADSET & DLOAD -> LSEQ')
        if 'IC' in subcase:
            ic_val = subcase.get_int_parameter('IC')[0]
            if ic_val != 0:
                log.debug(f'IC={ic_val}')
                tic = fem.tics[ic_val]
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
                    fem.Load(ic_val)
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
        loads, scale_factors = fem.get_reduced_dloads(dload_id)
        #fem.log.info('scale_factors=%s loads=%s' % (scale_factors, loads))

        if sol in {108, 111}:  # direct frequency, modal frequency
            for load2, scale_factor in zip(loads, scale_factors):
                freq_id = subcase.get_int_parameter('FREQ')[0]
                freqs: list[FREQS] = fem.frequencies[freq_id]
                for freq in freqs:
                    if freq.type in {'FREQ', 'FREQ1', 'FREQ2'}:
                        fmax = freq.freqs[-1]
                        force = load2.get_load_at_freq(fmax) * scale_factor
        elif sol in {109, 112, 118, 129, 146, 159, 400, 401, 601, 700}:
            # 109: direct transient (time linear)
            # 112: modal transient
            # 118: Cyclic direct frequency response
            # 129: time nonlinear
            # 146: gust
            # 159: thermal nonlinear
            # 400:
            # 401: NX nonlinear???
            # 601:
            # 700:
            for load2, scale_factor in zip(loads, scale_factors):
                force = load2.get_load_at_time(0.) * scale_factor
        elif  sol == 200:
            pass
        else:
            # 112-
            fem.log.warning(f'solution={sol}; DLOAD is not supported')

        # print(loads)
    if 'NLPARM' in subcase:
        allowed_sols = [66, 101, 106, 153, 400, 600, 700]
        ierror = check_sol(sol, subcase, allowed_sols, 'NLPARM', log, ierror, nerrors)
        nlparm_id = subcase.get_int_parameter('NLPARM')[0]
        unused_nlparm = fem.NLParm(nlparm_id, f', which is required for {subcase}')
    return ierror

def _check_case_parameters_aero(subcase: Subcase, fem: BDF, sol: int,
                                ierror: int=0, nerrors: int=100,
                                stop_on_failure: bool=True) -> int:
    """helper method for ``_check_case_parameters``"""
    log = fem.log
    if 'TRIM' in subcase:
        trim_id = subcase.get_int_parameter('TRIM')[0]
        if trim_id not in fem.trims:
            msg = (
                f'SOL={sol}\n'
                f'TRIM = {trim_id}\n'
                f'trims={fem.trims}\n'
                f'subcase:\n{subcase}')
            log_error(sol, [144, 200], msg, log)
        else:
            trim = fem.trims[trim_id]
            suport1 = None
            if 'SUPORT1' in subcase:
                suport_id = subcase.get_int_parameter('SUPORT1')[0]
                suport1 = fem.suport1[suport_id]
            try:
                trim.verify_trim(
                    fem.suport, suport1, fem.aestats, fem.aeparams,
                    fem.aelinks, fem.aesurf, xref=True)
            except RuntimeError:
                if stop_on_failure or ierror == nerrors:
                    raise
                ierror += 1
                exc_info = sys.exc_info()
                traceback.print_exception(*exc_info)
                #traceback.print_stack()
                #fem.log.error(e.msg)
                #raise
            assert 'DIVERG' not in subcase, subcase
            #allowed_sols = [144, 200]

    if 'DIVERG' in subcase:
        value = subcase.get_int_parameter('DIVERG')[0]
        assert value in fem.divergs, 'value=%s\n divergs=%s\n subcase:\n%s' % (value, str(fem.divergs), str(subcase))
        assert 'TRIM' not in subcase, subcase
        #allowed_sols = [144, 200]

    if 'FMETHOD' in subcase:
        # FLUTTER
        fmethod_id = subcase.get_int_parameter('FMETHOD')[0]
        unused_fmethod = fem.flutters[fmethod_id]
        allowed_sols = [145, 200]
        ierror = check_sol(sol, subcase, allowed_sols, 'FMETHOD', log, ierror, nerrors)
    return ierror

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


def compute(cards1: dict[str, int],
            cards2: dict[str, int], quiet: bool=False):
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

    version_group = parent_parser.add_mutually_exclusive_group()

    version_group_map = {
        '--msc' : 'Assume MSC Nastran (default=True)',
        '--nx' : 'Assume NX Nastran/Simcenter (default=False)',
        '--optistruct' : 'Assume Altair OptiStruct (default=False)',
        '--nasa95' : 'Assume Nastran 95 (default=False)',
        '--mystran' : 'Assume Mystran (default=False)',
    }
    for version, help_msg in version_group_map.items():
        version_group.add_argument(
            version, action='store_true',
            help=help_msg)

    parent_parser.add_argument(
        '-L', '--loads', action='store_false',
        help='Disables forces/moments summation for the different subcases (default=True)')

    parent_parser.add_argument('-e', '--nerrors', default=100, type=int,
                               help='Allow for cross-reference errors (default=100)')
    parent_parser.add_argument('--encoding', default=encoding, type=str,
                               help='the encoding method (default=%r)\n' % encoding)
    #parent_parser.add_argument('--skip_nominal', action='store_true',
                               #help='skip the nominal model comparison (default=False)')
    parent_parser.add_argument('--skip_loads', action='store_true',
                               help='skip loads calcuations (default=False)')
    parent_parser.add_argument('--skip_mass', action='store_true',
                               help='skip mass calcuations (default=False)')
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
    _set_version(args2)
    #optional_args = [
        #'double', 'large', 'crash', 'quiet', 'profile',
        #'xref', 'safe', 'check', 'punch', 'loads', 'stop', 'encoding',
        #'dumplines', 'dictsort', 'nerrors', 'pickle', 'hdf5',
    #]
    #for arg in optional_args:
        #swap_key(args2, arg, '--' + arg)
    return args2


def log_error(sol: int, error_solutions, msg: str, log: SimpleLogger) -> None:
    if sol in error_solutions:
        raise RuntimeError(msg)
    else:
        log.warning(msg)

def _set_version(args: Any):
    """sets the version flag"""
    if args['msc']:
        version = 'msc'
    elif args['nx']:
        version = 'nx'
    elif args['optistruct']:
        version = 'optistruct'
    elif args['nasa95']:
        version = 'nasa95'
    elif args['mystran']:
        version = 'mystran'
    else:
        version = None
    args['version'] = version
    del args['msc'], args['nx'], args['nasa95'], args['mystran'], args['optistruct']

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
    formats = '--msc|--nx|--optistruct|--nasa95|--mystran'
    options = (
        '\n  [options] = [-e E] [--encoding ENCODE] [-q] [--dumplines] [--dictsort]\n'
        f'              [--crash C] [--pickle] [--profile] [--hdf5] [{formats}] [--filter]\n'
        '              [--skip_loads] [--skip_mass]\n'
    )
    usage = (
        "Usage:\n"
        '  test_bdf [-x | --safe] [-p] [-c]       BDF_FILENAME [options]\n'
        '  test_bdf [-x | --safe] [-p] [-c] [-d]  BDF_FILENAME [options]\n'
        '  test_bdf [-x | --safe] [-p] [-c] [-l]  BDF_FILENAME [options]\n'
        '  test_bdf               [-p]            BDF_FILENAME [options]\n'
        '  test_bdf [-x | --safe] [-p] [--stop]   BDF_FILENAME [options]\n'
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
        '  --loads        Disables forces/moments summation for the different subcases (default=True)\n'
        #'  --filter       Filters unused cards\n'

        '  -e E, --nerrors E  Allow for cross-reference errors (default=100)\n'
        f'  --encoding ENCODE  the encoding method (default=None -> {encoding!r})\n'  +
        '  -q, --quiet        prints debug messages (default=False)\n'

        '\n'
        'Developer:\n'
        '  --crash C     Crash on specific cards (e.g. CGEN,EGRID)\n'
        '  --stop        Stop after first read/write (default=False)\n'
        '  --dumplines   Writes the BDF exactly as read with the INCLUDEs processed\n'
        '                (pyNastran_dump.bdf)\n'
        '  --dictsort    Writes the BDF exactly as read with the INCLUDEs processed\n'
        '                (pyNastran_dict.bdf)\n'
        '  --profile     Profiles the code (default=False)\n'
        '  --pickle      Pickles the data objects (default=False)\n'
        '  --hdf5        Save/load the BDF in HDF5 format\n'
        '  --msc         Assume MSC Nastran\n'
        '  --nx          Assume NX Nastran\n'
        '  --optistruct  Assume OptiStruct\n'
        '  --nasa95      Assume Nastran 95\n'
        '  --mystran     Assume Mystran\n'
        '  --skip_loads   skip the loads summation calculations (default=False)\n'
        '  --skip_mass    skip the mass properties calculations (default=False)\n'
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

    #data['run_nominal'] = not data['skip_nominal']
    data['run_loads'] = not data['skip_loads']
    data['run_mass'] = not data['skip_mass']

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
            #filter_unused=data['filter'],
            check=not(data['check']),
            punch=data['punch'],
            size=size,
            is_double=is_double,
            sum_load=data['run_loads'],
            #run_nominal=data['run_nominal'],
            #run_loads=data['run_loads'],
            run_mass=data['run_mass'],
            run_extract_bodies=False,

            stop=data['stop'],
            quiet=data['quiet'],
            dumplines=data['dumplines'],
            dictsort=data['dictsort'],
            nerrors=data['nerrors'],
            encoding=data['encoding'],
            crash_cards=crash_cards,
            pickle_obj=data['pickle'],
            safe_xref=data['safe'],
            hdf5=data['hdf5'],
            version=data['version'],
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
            #filter_unused=data['filter'],
            xref=data['xref'],
            # xref_safe=data['xref_safe'],
            check=not(data['check']),
            punch=data['punch'],
            size=size,
            is_double=is_double,
            sum_load=data['run_loads'],
            run_mass=data['run_mass'],
            run_extract_bodies=False,

            stop=data['stop'],
            quiet=data['quiet'],
            dumplines=data['dumplines'],
            dictsort=data['dictsort'],
            nerrors=data['nerrors'],
            encoding=data['encoding'],
            crash_cards=crash_cards,
            pickle_obj=data['pickle'],
            safe_xref=data['safe'],
            hdf5=data['hdf5'],
            version=data['version'],
            print_stats=True,
            stop_on_failure=False,
        )
    print("total time:  %.2f sec" % (time.time() - time0))


if __name__ == '__main__':  # pragma: no cover
    main()
