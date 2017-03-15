"""defines testing utils"""
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
    return model2
