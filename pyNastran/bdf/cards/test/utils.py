"""defines testing utils"""
import os
from six import StringIO
from pyNastran.bdf.bdf import BDF
from pyNastran.bdf.mesh_utils.remove_unused import remove_unused
from pyNastran.bdf.mesh_utils.convert import convert


def save_load_deck(model, punch=True, run_remove_unused=True, run_convert=True):
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

    model2.save(obj_filename='model.obj', unxref=True)
    model3 = BDF(debug=False, log=None, mode='msc')
    model3.load(obj_filename='model.obj')
    os.remove('model.obj')
    model3.cross_reference()
    model3.pop_xref_errors()
    return model3
