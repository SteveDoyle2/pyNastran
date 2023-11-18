#from pyNastran.bdf.cards.test.utils import save_load_deck
#from __future__ import annotations
import io
from typing import TYPE_CHECKING
import numpy as np
from pyNastran.dev.bdf_vectorized3.bdf import BDF
from pyNastran.dev.bdf_vectorized3.mesh_utils.convert import convert
from pyNastran.dev.bdf_vectorized3.mesh_utils.remove_unused import remove_unused
from pyNastran.dev.bdf_vectorized3.mesh_utils.bdf_equivalence import bdf_equivalence_nodes

#from pyNastran.dev.bdf_vectorized3.mesh_utils.convert import convert
#if TYPE_CHECKING:  # pragma: no cover

def save_load_deck(model: BDF,
                   xref: str='standard',
                   punch: bool=True,
                   run_remove_unused: bool=True,
                   run_convert: bool=True,
                   run_renumber: bool=True,
                   run_mirror: bool=True,
                   run_save_load: bool=True,
                   run_quality: bool=True,
                   write_saves: bool=True,
                   run_save_load_hdf5: bool=True,
                   run_mass_properties: bool=True,
                   run_loads: bool=True,
                   run_test_bdf: bool=True,
                   run_op2_writer: bool=True,
                   run_op2_reader: bool=True,
                   remove_disabled_cards: bool=True,
                   run_read_write: bool=True,
                   run_geom_check: bool=True,
                   run_equivalence: bool=True,
                   nastran_format: str='nx',
                   op2_log_level: str='warning') -> BDF:
    """writes, re-reads, saves an obj, loads an obj, and returns the deck"""
    model.setup(run_geom_check=run_geom_check)
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

    if run_convert:
        units_to = ['m', 'kg', 's']
        units = ['ft', 'lbm', 's']
        convert(model, units_to, units)

    if run_remove_unused:
        remove_unused(model, inplace=True)

    if run_read_write:
        stringio8 = io.StringIO()
        model.write_bdf(stringio8, close=False)
        stringio8.seek(0)

        stringio16 = io.StringIO()
        model.write_bdf(stringio16, size=16, close=False)
        stringio16.seek(0)

        stringio_double = io.StringIO()
        model.write_bdf(stringio_double, size=16, is_double=True, close=False)
        stringio_double.seek(0)

        model2 = BDF(debug=False, log=model.log)
        model2.read_bdf(stringio8, punch=model.punch)

        model3 = BDF(debug=False, log=model.log)
        model3.read_bdf(stringio16, punch=model.punch)

        model4 = BDF(debug=False, log=model.log)
        model4.read_bdf(stringio_double, punch=model.punch)

    if run_read_write and run_equivalence and len(model.grid):
        model2 = BDF(debug=False, log=model.log)
        model2.read_bdf(stringio8, punch=model.punch)
        bdf_filename_out = None

        tol = 0.0
        model_eq = bdf_equivalence_nodes(
            model2, bdf_filename_out, tol,
            renumber_nodes=False, neq_max=4, xref=True, node_set=None, size=8, is_double=False,
            remove_collapsed_elements=False, avoid_collapsed_elements=False,
            crash_on_collapse=False, log=None, debug=True, method='new')
    return model
