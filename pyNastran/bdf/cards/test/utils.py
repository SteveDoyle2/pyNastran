"""defines testing utils"""
import os
from copy import deepcopy
from six import StringIO
from pyNastran.bdf.bdf import BDF
from pyNastran.bdf.mesh_utils.remove_unused import remove_unused
from pyNastran.bdf.mesh_utils.convert import convert
from pyNastran.bdf.mesh_utils.bdf_renumber import bdf_renumber
from pyNastran.bdf.mesh_utils.mirror_mesh import bdf_mirror

def save_load_deck(model, punch=True, run_remove_unused=True,
                   run_convert=True, run_renumber=True, run_mirror=True):
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
    if run_remove_unused:
        remove_unused(model)
    if run_convert:
        units_to = ['m', 'kg', 's']
        units = ['ft', 'lbm', 's']
        convert(model, units_to, units)

    model2 = BDF(log=model.log)
    #print(bdf_file.getvalue())
    model2.read_bdf(bdf_file, punch=punch)
    model2.pop_parse_errors()
    model2.get_bdf_stats()
    model2.write_bdf('model2.bdf')

    model2.save(obj_filename='model.obj', unxref=True)
    model3 = BDF(debug=False, log=model.log, mode='msc')
    model3.load(obj_filename='model.obj')
    os.remove('model.obj')

    cross_reference(model3)
    if run_renumber:
        renumber('model2.bdf', model.log)
        if run_mirror:
            # we put embed this under renumber to prevent modifying an
            # existing model to prevent breaking tests
            #
            # shouldn't have any effect model2.bdf
            bdf_mirror('model2.bdf', plane='xz', log=model.log)
    return model3

def cross_reference(model):
    """validate we're doing xref right"""
    model.cross_reference()
    model.pop_xref_errors()

    model.safe_cross_reference()
    model.pop_xref_errors()

    model.cross_reference()
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
