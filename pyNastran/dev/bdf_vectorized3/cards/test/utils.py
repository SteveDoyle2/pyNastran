#from pyNastran.bdf.cards.test.utils import save_load_deck
#from __future__ import annotations
import io
from typing import TYPE_CHECKING
import numpy as np
from pyNastran.dev.bdf_vectorized3.bdf import BDF
#from pyNastran.dev.bdf_vectorized3.mesh_utils.convert import convert
#if TYPE_CHECKING:  # pragma: no cover

def save_load_deck(model: BDF, xref='standard', punch=True, run_remove_unused=True,
                   run_convert=True, run_renumber=True, run_mirror=True,
                   run_save_load=True, run_quality=True, write_saves=True,
                   run_save_load_hdf5=True, run_mass_properties=True, run_loads=True,
                   run_test_bdf=True, run_op2_writer=True, run_op2_reader=True,
                   remove_disabled_cards=True,
                   run_read_write=True,
                   nastran_format: str='nx',
                   op2_log_level: str='warning') -> BDF:
    """writes, re-reads, saves an obj, loads an obj, and returns the deck"""
    model.setup(run_geom_check=True)
    if run_quality:
        model.quality()

    if run_mass_properties:
        eids_mass, mass1 = model.mass()
        eids_inertia, massi1, cgi1, inertiai1 = model.inertia()
        mass2 = model.mass_sum(element_id=None)
        massi2, cgi2, inertiai2 = model.inertia_sum(element_id=None)
        assert np.allclose(mass1, massi1)
        assert np.allclose(mass2, massi2)
        assert np.allclose(mass1.sum(), mass2.sum())

    if run_loads and 0:
        model.sum_forces_moments()

    if run_convert and 0:
        units_to = ['m', 'kg', 's']
        units = ['ft', 'lbm', 's']
        convert(model, units_to, units)

    if run_read_write:
        stringio = io.StringIO()
        model.write_bdf(stringio, close=False)
        stringio.seek(0)

        model2 = BDF(debug=False, log=model.log)
        model2.read_bdf(stringio, punch=model.punch)
    return model
