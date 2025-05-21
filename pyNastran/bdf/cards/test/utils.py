"""defines testing utils"""
import os
from io import StringIO
from collections import ChainMap
import inspect

import numpy as np
from cpylog import SimpleLogger

from pyNastran.utils import PathLike
from pyNastran.bdf.bdf import BDF, read_bdf, DMIAX
from pyNastran.bdf.mesh_utils.delete_bad_elements import element_quality
from pyNastran.bdf.mesh_utils.remove_unused import remove_unused
from pyNastran.bdf.mesh_utils.convert import convert
from pyNastran.bdf.mesh_utils.bdf_renumber import bdf_renumber
from pyNastran.bdf.mesh_utils.mirror_mesh import bdf_mirror
from pyNastran.bdf.mesh_utils.mass_properties import (
    mass_properties, mass_properties_nsm, mass_properties_breakdown)
from pyNastran.bdf.mesh_utils.forces_moments import get_load_arrays, get_pressure_array
from pyNastran.bdf.mesh_utils.export_caero_mesh import export_caero_mesh

from pyNastran.bdf.test.test_bdf import run_bdf as test_bdf

from pyNastran.op2.tables.oug.oug_displacements import RealDisplacementArray
from pyNastran.op2.op2_geom import read_op2_geom, attach_op2_results_to_bdf, OP2Geom


try:
    import h5py
    IS_H5PY = True
except ModuleNotFoundError:  # pragma: no cover
    IS_H5PY = False

def save_load_deck(model: BDF, xref='standard', punch=True, run_remove_unused=True,
                   run_convert=True, run_renumber=True, run_mirror=True,
                   run_save_load=True, run_quality=True, write_saves=True,
                   run_save_load_hdf5=True, run_mass_properties=True, run_loads=True,
                   run_export_caero=True,
                   run_test_bdf=True, run_op2_writer=True, run_op2_reader=True,
                   remove_disabled_cards=True,
                   nastran_format: str='nx',
                   op2_log_level: str='warning') -> BDF:
    """writes, re-reads, saves an obj, loads an obj, and returns the deck"""
    if os.path.exists('junk.bdf'):
        os.remove('junk.bdf')
    if not remove_disabled_cards:
        model._add_disabled_cards()
    model.set_error_storage(nparse_errors=0, stop_on_parsing_error=True,
                            nxref_errors=0, stop_on_xref_error=True)
    model._remove_disabled_cards = remove_disabled_cards
    model.validate()
    model.pop_parse_errors()
    model.pop_xref_errors()
    bdf_file = StringIO()
    model.write_bdf(bdf_file, size=8, close=False)
    bdf_file.seek(0)
    model.write_bdf(bdf_file, size=16, close=False)
    bdf_file.seek(0)
    model.write_bdf(bdf_file, size=16, is_double=True, close=False)
    bdf_file.seek(0)
    get_matrices(model)

    if write_saves and model.save_file_structure:
        bdf_filenames = {0 : 'junk.bdf',}
        model.write_bdfs(bdf_filenames)
        os.remove('junk.bdf')

    if run_convert:
        units_to = ['m', 'kg', 's']
        units = ['ft', 'lbm', 's']
        convert(model, units_to, units)

    model2 = BDF(log=model.log)  # , mode=nastran_format
    #print(bdf_file.getvalue())
    if not remove_disabled_cards:
        model2._add_disabled_cards()
    model2.read_bdf(bdf_file, punch=punch, xref=False)
    _cross_reference(model2, xref)

    model2.pop_parse_errors()
    model2.get_bdf_stats()
    model2.write_bdf('model2.bdf')

    # pulls nastran_format from the header
    model_parse = BDF(log=model.log)
    model_parse._parse = False
    model_parse.read_bdf('model2.bdf')
    #head('model2.bdf')

    if run_test_bdf:
        folder = ''
        log_error = SimpleLogger(level='error', encoding='utf-8')
        test_bdf(folder, 'model2.bdf', stop_on_failure=True,
                 punch=punch,
                 quiet=True, log=log_error, is_lax_parser=True)
        os.remove('model2.test_bdf.bdf')

    nelements = len(model2.elements) + len(model2.masses)
    nnodes = len(model2.nodes) + len(model2.spoints) + len(model2.epoints)

    _run_mass_properties(model2, nnodes, nelements, run_mass_properties=run_mass_properties)
    _run_loads(model2, nelements, run_loads=run_loads)

    if run_save_load:
        model2.save(obj_filename='model.obj', unxref=True)
        model3 = BDF(debug=False, log=model.log, mode='msc')
        model3.load(obj_filename='model.obj')
        os.remove('model.obj')
    else:
        model2.uncross_reference()
        model3 = model2

    _run_hdf5(model2, model.log, run_save_load_hdf5=run_save_load_hdf5)

    cross_reference(model3, xref)
    if run_renumber:
        renumber('model2.bdf', model.log)
        if run_mirror:
            # we put this under renumber to prevent modifying an
            # existing model to prevent breaking tests
            #
            # shouldn't have any effect model2.bdf
            model_mirrored = bdf_mirror('model2.bdf', plane='xz', log=model.log)[0]
            model_mirrored.write_bdf('mirrored2.bdf')
            read_bdf('mirrored2.bdf', log=model.log)
            os.remove('mirrored2.bdf')
    os.remove('model2.bdf')

    if model.elements and run_quality:
        element_quality(model)

    if len(model.caeros) and run_export_caero:
        base = os.path.splitext(model.bdf_filename)[0]
        caero_bdf_filename = base + '.caero.bdf'
        export_caero_mesh(model, caero_bdf_filename,
                          is_subpanel_model=True,
                          pid_method='caero', write_panel_xyz=True)
        os.remove(caero_bdf_filename)

    read_write_op2_geom(
        model,
        run_op2_writer=run_op2_writer,
        run_op2_reader=run_op2_reader,
        nastran_format=nastran_format,
        op2_log_level=op2_log_level)

    if run_remove_unused:
        remove_unused(model)
    return model3


def head(filename: PathLike, nlines: int=20) -> None:  # pragma: no cover
    with open(filename, 'r') as file_obj:
        for i in range(nlines):
            line = file_obj.readline().rstrip()
            print(line)


def read_write_op2_geom(model: BDF,
                        run_op2_writer: bool=True,
                        run_op2_reader: bool=True,
                        nastran_format: str='nx',
                        op2_log_level: str='warning') -> None:
    """helper method for ``save_load_deck``"""
    if not run_op2_writer:
        return
    op2_geom_model = attach_displacement_to_bdf_model(model)

    op2_filename = 'spike.op2'
    bkp_log = op2_geom_model.log
    op2_geom_model.log = SimpleLogger(level=op2_log_level, encoding='utf-8')
    op2_geom_model.write_op2(op2_filename, post=-1, endian=b'<', skips=None,
                             nastran_format=nastran_format)
    if run_op2_reader:
        unused_op2_geom = read_op2_geom(op2_filename, log=op2_geom_model.log, xref=False)
    else:
        frame = inspect.currentframe()
        call_frame = inspect.getouterframes(frame, 2)
        op2_geom_model.log.warning('skipping op2 reader for %s' % call_frame[1][3])
    op2_geom_model.log = bkp_log
    os.remove(op2_filename)

def attach_displacement_to_bdf_model(model: BDF) -> OP2Geom:
    """helper method for ``save_load_deck``"""
    op2_geom_model = attach_op2_results_to_bdf(model, op2_model=None)

    table_name = 'OUGV1'
    node_gridtype = np.zeros((10, 2), dtype='int32')
    node_gridtype[:, 0] = np.arange(1, 11)
    data = np.zeros((1, 10, 6), dtype='float32')
    isubcase = 1
    disp = RealDisplacementArray.add_static_case(
        table_name, node_gridtype, data, isubcase, is_sort1=True)
    op2_geom_model.displacements[isubcase] = disp
    return op2_geom_model

def _run_hdf5(model2: BDF, log: SimpleLogger,
              run_save_load_hdf5: bool=True) -> None:
    """helper method"""
    if run_save_load_hdf5 and IS_H5PY:
        hdf5_filename = 'test.h5'
        model2.export_hdf5_filename(hdf5_filename)
        model4 = BDF(log=model2.log)
        model4.load_hdf5_filename(hdf5_filename)
        model4.validate()
        bdf_stream = StringIO()
        model4.write_bdf(bdf_stream, encoding=None, size=8, is_double=False,
                         interspersed=False, enddata=None, write_header=True, close=True)
        for key, unused_value in model2.card_count.items():
            if key == 'ENDDATA':
                continue
            if key not in model4.card_count:
                msg = 'key=%r was not loaded to hdf5\nexpected=%s\nactual=%s' % (
                    key, model2.card_count, model4.card_count)
                #raise RuntimeError(msg)
                log.error(msg)
        os.remove(hdf5_filename)

def _run_mass_properties(model2: BDF, nnodes: int, nelements: int,
                         run_mass_properties: bool=True) -> None:
    """helper method"""
    if not(run_mass_properties and nelements):
        return

    if nelements > 1 and nnodes == 0:  # pragma: no cover
        raise RuntimeError('no nodes exist')
    mass1, cg1, inertia1 = mass_properties(model2, reference_point=None, sym_axis=None)
    mass2, cg2, inertia2 = mass_properties_nsm(model2, reference_point=None, sym_axis=None)
    #if not quiet:
        #if model2.wtmass != 1.0:
            #print('weight = %s' % (mass1 / model2.wtmass))
        #print('mass = %s' % mass1)
        #print('cg   = %s' % cg1)
        #print('Ixx=%s, Iyy=%s, Izz=%s \nIxy=%s, Ixz=%s, Iyz=%s' % tuple(inertia1))
    assert np.allclose(mass1, mass2), 'mass1=%s mass2=%s' % (mass1, mass2)
    assert np.allclose(cg1, cg2), 'mass=%s\ncg1=%s cg2=%s' % (mass1, cg1, cg2)
    if not np.allclose(inertia1, inertia2, atol=1e-5):  # pragma: no cover
        raise ValueError('mass=%s cg=%s\ninertia1=%s\ninertia2=%s\ndinertia=%s' % (
            mass1, cg1, inertia1, inertia2, inertia1-inertia2))

    unused_mass3, unused_cg3, unused_inertia3 = mass_properties_breakdown(model2)[:3]
    #assert np.allclose(mass1, mass3), 'mass1=%s mass3=%s' % (mass1, mass3)
    #assert np.allclose(cg1, cg3), 'mass=%s\ncg1=%s cg3=%s' % (mass1, cg1, cg3)

def _run_loads(model: BDF, nelements: int, run_loads: bool=True) -> None:
    """helper method"""
    if not run_loads:
        return
    assert isinstance(nelements, int), nelements
    #node_ids = None

    #nnodes = model.nnodes
    #node_ids = [None] * nnodes
    try:
        out = model.get_displacement_index_xyz_cp_cd(
            fdtype='float64', idtype='int32', sort_ids=True)
    except ValueError:
        return
    unused_icd_transformi, unused_icp_transformi, unused_xyz_cpi, nid_cp_cd = out
    node_ids = nid_cp_cd[:, 0]
    normals, eids, eid_map = _normals_eid_map(model, nelements, fdtype='float64')

    nsubcases = len(model.subcases)
    if nsubcases == 0:
        load_ids = list(set(list(model.load_combinations) + list(model.loads)))
        for load_id in sorted(load_ids):
            unused_subcase = model.case_control_deck.add_subcase(load_id)
            get_pressure_array(model, load_id, eids, stop_on_failure=True,
                               fdtype='float32')

    load_ids = list(set(list(model.load_combinations) + list(model.loads)))
    for subcase_id in model.subcases:
        get_load_arrays(model, subcase_id, eid_map, node_ids, normals, nid_map=None,
                        fdtype='float32')

    if nsubcases == 0:
        del model.case_control_deck

def _normals_eid_map(model: BDF, nelements: int,
                     fdtype: str='float64') -> tuple[np.ndarray,
                                                     list[int],
                                                     dict[int, int]]:
    ieid = 0
    eids = []
    eid_map = {}
    normals = np.zeros((nelements, 3), dtype=fdtype)
    for eid, elem in model.elements.items():
        if hasattr(elem, 'Normal'):
            normals[ieid, :] = elem.Normal()
        eid_map[eid] = ieid
        eids.append(eid)
        ieid += 1
    return normals, eids, eid_map

def _cross_reference(model: BDF, xref: bool | str) -> None:
    """helper method for ``_cross_reference``"""
    if xref in [True, 'standard']:
        model.cross_reference()
        model.safe_cross_reference()
    elif xref in ['safe']:
        model.safe_cross_reference()

def cross_reference(model: BDF, xref: bool) -> None:
    """validate we're doing xref right"""
    _cross_reference(model, xref)
    model.pop_xref_errors()

    model.safe_cross_reference()
    model.pop_xref_errors()

    _cross_reference(model, xref)
    model.pop_xref_errors()


def renumber(bdf_filename: PathLike, log: SimpleLogger) -> None:
    bdf_filename_out = 'junk.bdf'
    #model3_copy = deepcopy(model3)
    #model3.cross_reference()
    bdf_renumber(bdf_filename, bdf_filename_out, size=8, is_double=False,
                 starting_id_dict=None, round_ids=False, cards_to_skip=None, log=None,
                 debug=False)
    model4 = BDF(debug=False, log=log)
    model4.read_bdf(bdf_filename_out)
    os.remove('junk.bdf')

def get_matrices(model: BDF) -> None:
    """tests the ``get_matrix`` method for every type of matrix"""
    dicts = ChainMap(model.dmig, model.dmij, model.dmiji, model.dmik,
                     model.dmi, model.dmiax)
    for key, matrix in dicts.items():
        if key == 'UACCEL' or isinstance(matrix, DMIAX):
            model.log.warning(f'skipping get_matrix() for name={key!r}; type={matrix.type}')
            continue
        sparse, *junk = matrix.get_matrix(is_sparse=True, apply_symmetry=False)
        dense, *junk = matrix.get_matrix(is_sparse=False, apply_symmetry=False)
        dense2 = sparse.toarray()
        assert np.array_equal(dense, dense2)
