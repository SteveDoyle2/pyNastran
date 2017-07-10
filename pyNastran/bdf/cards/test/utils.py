"""defines testing utils"""
import os
from six import StringIO
from pyNastran.bdf.bdf import BDF


def save_load_deck(model):
    """writes and re-reads a deck"""
    bdf_file = StringIO()
    model.write_bdf(bdf_file, size=8, close=False)
    bdf_file.seek(0)
    model.write_bdf(bdf_file, size=16, close=False)
    bdf_file.seek(0)

    model2 = BDF(log=model.log)
    model2.read_bdf(bdf_file, punch=True)
    model2.pop_parse_errors()

    model2.save(obj_filename='model.obj', unxref=True)
    model3 = BDF(debug=False, log=None, mode='msc')
    model3.load(obj_filename='model.obj')
    os.remove('model.obj')
    model3.cross_reference()
    model3.pop_xref_errors()
    return model3
