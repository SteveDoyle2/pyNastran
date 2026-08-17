from __future__ import annotations
import os
from typing import TYPE_CHECKING

from pyNastran.utils import PathLike, print_bad_path
from pyNastran.op2.op2 import read_op2
from pyNastran.op2.tables.lama_eigenvalues.lama_objects import RealEigenvalues
if TYPE_CHECKING:
    from pyNastran.bdf.bdf import BDF


def nastran_to_zaero_f06(op2_filename: PathLike,
                         f06_filename: PathLike='',
                         bdf_filename: PathLike | BDF='',
                         build: str | None='SIMCENTER NASTRAN 11/ 8/24'):
    assert os.path.exists(op2_filename), print_bad_path(op2_filename)

    if bdf_filename:
        assert os.path.exists(bdf_filename), print_bad_path(bdf_filename)
    if f06_filename == '':
        base = os.path.splitext(op2_filename)[0]
        f06_filename = base + '.f06'
    assert isinstance(f06_filename, PathLike), (f06_filename, type(f06_filename))

    results_model = read_op2(op2_filename, load_geometry=True,
                             include_results=['eigenvectors',])
    # print(results_model.nodes)
    if bdf_filename and len(results_model.nodes) == 0:
        from pyNastran.bdf.bdf import read_bdf
        model = read_bdf(bdf_filename)
        results_model.replace_cards(model, write_log=False)

    assert len(results_model.eigenvectors) == 1, len(results_model.eigenvectors)
    if len(results_model.eigenvalues) == 0:
        for i, eigenvector_obj in results_model.eigenvectors.items():
            eigenvalue_obj = RealEigenvalues.from_eigenvectors(eigenvector_obj)
            results_model.eigenvalues[''] = eigenvalue_obj
            break
    results_model.log.info(f'writing eigenvectors to {str(f06_filename)}')
    results_model.build = build
    results_model.write_f06(f06_filename, write_cards=True)
