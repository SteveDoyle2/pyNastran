"""defines testing utils"""
import os
from copy import deepcopy
from six import StringIO
from pyNastran.bdf.bdf import BDF
from pyNastran.bdf.mesh_utils.delete_bad_elements import element_quality
from pyNastran.bdf.mesh_utils.remove_unused import remove_unused
from pyNastran.bdf.mesh_utils.convert import convert
from pyNastran.bdf.mesh_utils.bdf_renumber import bdf_renumber
from pyNastran.bdf.mesh_utils.mirror_mesh import bdf_mirror

try:
    import h5py
    IS_H5PY = True
except ImportError:
    IS_H5PY = False

def save_load_deck(model, xref='standard', punch=True, run_remove_unused=True,
                   run_convert=True, run_renumber=True, run_mirror=True,
                   run_save_load=True, run_quality=True, write_saves=True,
                   run_save_load_hdf5=True):
    """writes, re-reads, saves an obj, loads an obj, and returns the deck"""
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

    if write_saves and model.save_file_structure:
        bdf_filenames = {0 : 'junk.bdf',}
        model.write_bdfs(bdf_filenames)
        os.remove('junk.bdf')

    if run_remove_unused:
        remove_unused(model)
    if run_convert:
        units_to = ['m', 'kg', 's']
        units = ['ft', 'lbm', 's']
        convert(model, units_to, units)

    model2 = BDF(log=model.log)
    #print(bdf_file.getvalue())
    model2.read_bdf(bdf_file, punch=punch, xref=False)
    _cross_reference(model2, xref)

    model2.pop_parse_errors()
    model2.get_bdf_stats()
    model2.write_bdf('model2.bdf')

    if run_save_load:
        model2.save(obj_filename='model.obj', unxref=True)
        model3 = BDF(debug=False, log=model.log, mode='msc')
        model3.load(obj_filename='model.obj')
        os.remove('model.obj')
    else:
        model2.uncross_reference()
        model3 = model2

    if run_save_load_hdf5 and IS_H5PY:
        model2.export_to_hdf5_filename('test.h5')
        model4 = BDF(log=model2.log)
        model4.load_hdf5_filename('test.h5')
        model4.validate()
        bdf_stream = StringIO()
        model4.write_bdf(bdf_stream, encoding=None, size=8, is_double=False,
                         interspersed=False, enddata=None, write_header=True, close=True)
        for key, value in model2.card_count.items():
            if key == 'ENDDATA':
                continue
            if key not in model4.card_count:
                msg = 'key=%r was not loaded to hdf5\nexpected=%s\nactual=%s' % (
                    key, model2.card_count, model4.card_count)
                #raise RuntimeError(msg)
                model.log.error(msg)

    cross_reference(model3, xref)
    if run_renumber:
        renumber('model2.bdf', model.log)
        if run_mirror:
            # we put embed this under renumber to prevent modifying an
            # existing model to prevent breaking tests
            #
            # shouldn't have any effect model2.bdf
            bdf_mirror('model2.bdf', plane='xz', log=model.log)
    os.remove('model2.bdf')

    if model.elements and run_quality:
        element_quality(model)
    return model3

def _cross_reference(model, xref):
    """helper method for ``_cross_reference``"""
    if xref in [True, 'standard']:
        model.cross_reference()
    elif xref in ['safe']:
        model.safe_cross_reference()

def cross_reference(model, xref):
    """validate we're doing xref right"""
    _cross_reference(model, xref)
    model.pop_xref_errors()

    model.safe_cross_reference()
    model.pop_xref_errors()

    _cross_reference(model, xref)
    model.pop_xref_errors()


def renumber(bdf_filename, log):
    bdf_filename_out = 'junk.bdf'
    #model3_copy = deepcopy(model3)
    #model3.cross_reference()
    bdf_renumber(bdf_filename, bdf_filename_out, size=8, is_double=False,
                 starting_id_dict=None, round_ids=False, cards_to_skip=None, log=None,
                 debug=False)
    model4 = BDF(debug=False, log=log)
    model4.read_bdf(bdf_filename_out)
