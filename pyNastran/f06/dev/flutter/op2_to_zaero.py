from __future__ import annotations
import os
from typing import TYPE_CHECKING

from pyNastran.utils import PathLike, print_bad_path
from pyNastran.op2.op2 import read_op2
if TYPE_CHECKING:
    from pyNastran.bdf.bdf import BDF


def nastran_to_zaero_f06(op2_filename: PathLike,
                         f06_filename: PathLike='',
                         bdf_filename: PathLike | BDF=''):
    assert os.path.exists(op2_filename), print_bad_path(op2_filename)

    if bdf_filename:
        assert os.path.exists(bdf_filename), print_bad_path(bdf_filename)
    if f06_filename == '':
        base = os.path.splitext(op2_filename)[0]
        f06_filename = base + '.f06'

    results_model = read_op2(op2_filename, load_geometry=True,
                             include_results=['eigenvectors', 'displacements'])
    # print(results_model.nodes)
    if bdf_filename and len(results_model.nodes) == 0:
        from pyNastran.bdf.bdf import read_bdf
        model = read_bdf(bdf_filename)
        results_model.replace_cards(model, write_log=False)
        asdf

    assert isinstance(f06_filename, str), (f06_filename, type(f06_filename))
    results_model.write_f06(f06_filename, write_cards=True)
