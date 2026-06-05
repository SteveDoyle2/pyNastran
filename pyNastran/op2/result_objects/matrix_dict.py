"""Defines the MatrixDict class.

SIL-EQEXIN-KDICT-KELM Relationship
====================================

SIL (Scalar Index List)
    Nastran's internal DOF numbering. Each grid gets a contiguous block
    of scalar indices spaced by 6 (regardless of element DOF count):
      - SIL_base = (internal_sequence - 1) * 6 + 1
      - 6-DOF grids (shells/bars): uses SIL_base through SIL_base+5
      - 3-DOF grids (solids): uses SIL_base through SIL_base+2

EQEXIN (Equivalence External-Internal)
    Maps external grid IDs (BDF) to internal sequence numbers.
    Internal ordering may differ from external due to resequencing.
    Without EQEXIN, SIL values cannot be mapped back to external grid IDs.

KDICT (MatrixDict)
    Describes KELM structure per element group:
      - sils: (neids, max_nodes) — SIL base values for each element's grids.
        Zero entries are unused slots (e.g., CHEXA8 uses 8 of 20 slots).
      - dof_per_grid: 3 for solids, 6 for shells/bars/rods
      - xforms: 3x3 rotation matrices (element coord -> global coord).
        None for elements already in global coords (CBAR, CBEAM).
      - eids: element IDs per group
      - address: (neids, 2) — not needed for extraction (columns are sequential)

KELM
    Shape (max_tri_size, nelements). Each column stores one element's
    stiffness as a column-major packed lower triangle::

        ndof = count_nonzero(sil_row) * dof_per_grid
        tri_size = ndof * (ndof + 1) / 2

    Column order matches sequential element order through KDICT groups.

Assembly into KGG::

    global_dofs = [(sil - 1) ... (sil - 1 + dof_per_grid - 1)] per node
    KGG[global_dofs, global_dofs] += T^T @ Ke @ T

    where T = block_diag(xform, xform, ...) per node.

"""
from __future__ import annotations
from typing import TYPE_CHECKING

import numpy as np
from scipy import sparse

from pyNastran.op2.op2_interface.op2_codes import MSC_ELEMENTS

if TYPE_CHECKING:
    from pyNastran.op2.result_objects.matrix import Matrix


class MatrixDict:
    """storage object for KDICT, MDICT, BDICT, etc. is op2.matdicts"""
    def __init__(self, name):
        self.name = name
        self.element_types = []
        self.numwides = []
        self.numgrids = []
        self.dof_per_grids = []

        self.eids = []
        self.ge = []
        self.address = []
        self.forms = []
        self.sils = []
        self.xforms = []

    def add(self, eltype, numwids, numgrid, dof_per_grid, form,
            eids, ge, address, sil, xform=None):
        """Sets the next set of the KDICT"""
        self.element_types.append(eltype)
        self.numwides.append(numwids)
        self.numgrids.append(numgrid)
        self.dof_per_grids.append(dof_per_grid)
        self.forms.append(form)

        self.eids.append(eids)
        self.ge.append(ge)
        self.address.append(address)
        self.sils.append(sil)
        self.xforms.append(xform)

    #@property
    #def nodes(self):
        #return [sil // 10 for sil in self.sils]

    #@property
    #def dofs(self):
        #return [sil % 10 for sil in self.sils]

    @property
    def nelements(self):
        return sum([len(eids) for eids in self.eids])

    @property
    def element_names(self):
        return [MSC_ELEMENTS[etype] for etype in self.element_types]

    def get_element_matrices(self, kelm: Matrix) -> dict[int, np.ndarray]:
        """Extract individual element stiffness/mass matrices from KELM/MELM.

        Uses the SIL (Scalar Index List) to determine the actual number of
        DOFs per element. Each element's packed lower triangle is unpacked
        into a full symmetric matrix.

        Parameters
        ----------
        kelm : Matrix
            The KELM or MELM matrix object from op2.matrices

        Returns
        -------
        element_matrices : dict[int, np.ndarray]
            Mapping of element_id -> (ndof, ndof) symmetric numpy array.
            DOF ordering follows the SIL order (grid1_dofs, grid2_dofs, ...).

        Notes
        -----
        KELM storage: each column is one element (sequential across groups).
        Each column stores the lower triangle in column-major (Fortran) order.

        For shell/bar elements: 6 DOF per grid (T1,T2,T3,R1,R2,R3)
        For solid elements: 3 DOF per grid (T1,T2,T3)

        Matrix sizes reflect the full DOF space of the element's grids:
          - CHEXA:  24x24 (8 nodes * 3 DOF)
          - CPENTA: 18x18 (6 nodes * 3 DOF)
          - CTETRA: 12x12 (4 nodes * 3 DOF)
          - CQUAD4: 24x24 (4 nodes * 6 DOF)
          - CTRIA3: 18x18 (3 nodes * 6 DOF)
          - CBAR:   12x12 (2 nodes * 6 DOF)
          - CBEAM:  12x12 (2 nodes * 6 DOF)
          - CROD:   12x12 (2 nodes * 6 DOF)
          - CELAS:   2x2  (2 nodes * 1 DOF)

        Elements like CROD only use a subset of DOFs (axial + torsion),
        so most entries in the 12x12 are zero. The physical stiffness is
        a 4x4 submatrix (DOFs 1 and 4 at each grid). To get the reduced
        matrix, extract rows/columns corresponding to non-zero DOFs.

        """
        data = kelm.data.toarray() if sparse.issparse(kelm.data) else np.asarray(kelm.data)

        element_matrices = {}
        col_idx = 0
        for sil, eids, dof_per_grid in zip(self.sils, self.eids, self.dof_per_grids):
            for j, eid in enumerate(eids):
                nnodes = np.count_nonzero(sil[j])
                ndof = nnodes * dof_per_grid
                tri_size = ndof * (ndof + 1) // 2
                packed = data[:tri_size, col_idx]
                element_matrices[eid] = _unpack_lower_triangle(packed, ndof)
                col_idx += 1
        return element_matrices

    def get_element_dof_info(self) -> dict[int, tuple[np.ndarray, int]]:
        """Get DOF connectivity info for each element.

        Returns
        -------
        dof_info : dict[int, tuple[sil_values, dof_per_grid]]
            Mapping of element_id -> (non-zero SIL values, dof_per_grid).
            SIL values encode grid DOF positions in the global matrix.
            Grid IDs can be recovered as sil_value // 10, component as sil_value % 10.

        """
        dof_info = {}
        for sil, eids, dof_per_grid in zip(self.sils, self.eids, self.dof_per_grids):
            for j, eid in enumerate(eids):
                nonzero_mask = sil[j] != 0
                dof_info[eid] = (sil[j, nonzero_mask], dof_per_grid)
        return dof_info

    def __repr__(self):
        msg = 'MatrixDict(name=%r, nelements=%s element_types=%s, element_names=[%s])' % (
            self.name, self.nelements, self.element_types, ', '.join(self.element_names))
        return msg


def _unpack_lower_triangle(packed: np.ndarray, ndof: int) -> np.ndarray:
    """Unpack a column-major lower triangle into a full symmetric matrix.

    Parameters
    ----------
    packed : ndarray of shape (ndof*(ndof+1)/2,)
        Packed lower triangle in column-major (Fortran) order.
    ndof : int
        Number of degrees of freedom (matrix dimension).

    Returns
    -------
    mat : ndarray of shape (ndof, ndof)
        Full symmetric matrix.

    """
    mat = np.zeros((ndof, ndof), dtype=packed.dtype)
    idx = 0
    for col in range(ndof):
        for row in range(col, ndof):
            mat[row, col] = packed[idx]
            mat[col, row] = packed[idx]
            idx += 1
    return mat
