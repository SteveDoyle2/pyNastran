from __future__ import annotations
import getpass
import numpy as np
from typing import Union, TYPE_CHECKING

from pyNastran.femutils.utils import safe_norm
from .vector_results import DispForceVectorResults, _to_dense_vector

if TYPE_CHECKING:
    from pyNastran.op2.result_objects.table_object import (
        RealTableArray, ComplexTableArray)

translation = ['Magnitude', 'Tx', 'Ty', 'Tz']
rotation = ['Magnitude', 'Rx', 'Ry', 'Rz']
col_axis = 1


## TODO: better label to indicate Magnitude vs. Value
class DisplacementResults2(DispForceVectorResults):
    def __init__(self,
                 subcase_id: int,
                 node_id: np.ndarray,
                 xyz: np.ndarray,
                 dxyz: Union[RealTableArray, ComplexTableArray],
                 title: str,
                 t123_offset: int,
                 dim_max: float=1.0,
                 data_format: str='%g',
                 is_variable_data_format: bool=False,
                 nlabels=None, labelsize=None, ncolors=None,
                 colormap: str='',
                 set_max_min: bool=False,
                 uname: str='DisplacementResults2'):
        """
        Defines a Displacement/Eigenvector result

        Parameters
        ----------
        subcase_id : int
            the flag that points to self.subcases for a message
        headers : list[str]
            the sidebar word
        titles : list[str]
            the legend title
        xyz : (nnodes, 3)
            the nominal xyz locations
        case : RealTableArray, ComplexTableArray
            the delta xyz values
        scalars : (nnodes,n) float ndarray
            #the data to make a contour plot with
            does nothing
        scales : list[float]
            the deflection scale factors
            nominally, this starts as an empty list and is filled later
        data_formats : str
            the type of data result (e.g. '%i', '%.2f', '%.3f')
        ncolors : int; default=None
            sets the default for reverting the legend ncolors
        set_max_min : bool; default=False
            set default_mins and default_maxs

        Unused
        ------
        uname : str
            some unique name for ...
        """
        methods_txyz_rxyz = ['Tx', 'Ty', 'Tz', 'Rx', 'Ry', 'Rz']
        index_to_base_title_annotation = {
            0: {'title': 'T_', 'corner': 'T_'},
            3: {'title': 'R_', 'corner': 'R_'},
        }
        DispForceVectorResults.__init__(
            self,
            subcase_id,
            node_id,
            dxyz,
            t123_offset,
            methods_txyz_rxyz,
            index_to_base_title_annotation,
            dim_max,
            data_format=data_format,
            nlabels=nlabels, labelsize=labelsize, ncolors=ncolors,
            colormap=colormap,
            set_max_min=set_max_min,
            uname=uname)
        self.title = title

        self.is_variable_data_format = is_variable_data_format

        # setup the node mapping
        disp_nodes = dxyz.node_gridtype[:, 0]  #  local node id
        self.inode = np.searchsorted(node_id, disp_nodes)

        # dense -> no missing nodes in the results set
        self.is_dense = (len(node_id) == len(disp_nodes))

        self.xyz = xyz
        assert len(self.xyz.shape) == 2, self.xyz.shape
        assert isinstance(t123_offset, int), t123_offset
        self.location = 'node'
        str(self)

    def _calculate_scale(self, itime: int, res_name: str) -> float:
        fringe_data = self._get_fringe_data_sparse(itime, res_name)

        # Fringe is typically a 'Magnitude' and therefore positive.
        # if 'Value' is used, it can be negative.
        abs_maximax = np.abs(fringe_data).max()
        if abs_maximax > 0.0:
            scale = self.dim_max / abs_maximax * 0.10
        else:
            scale = 1.0
        return scale

    #-------------------------------------
    # unmodifyable getters

    def deflects(self, unused_i: int, unused_res_name: str) -> bool:
        """deflection is opt-in"""
        return True

    def has_output_checks(self, i: int, resname: str) -> tuple[bool, bool, bool,
                                                               bool, bool, bool]:
        is_enabled_fringe = True
        is_checked_fringe = True
        is_enabled_disp = True
        is_checked_disp = True
        is_enabled_vector = False
        is_checked_vector = False
        out = (
            is_enabled_fringe, is_checked_fringe,
            is_enabled_disp, is_checked_disp,
            is_enabled_vector, is_checked_vector)
        return out

    def get_force_vector_result(self, itime: int, res_name: str,
                                ) -> tuple[np.ndarray, np.ndarray]:
        dxyz, *unused_junk = self.get_vector_data_dense(itime, res_name)
        assert dxyz.ndim == 2, dxyz
        scale = 1.
        return self.xyz, dxyz * scale

    def get_vector_result(self, itime: int, res_name: str,
                          return_dense: bool=True) -> tuple[np.ndarray, np.ndarray]:
        dxyz, *unused_junk = self.get_vector_data_dense(itime, res_name)
        assert dxyz.ndim == 2, dxyz

        scale = self.get_scale(itime, res_name)
        deflected_xyz = self.xyz + scale * dxyz
        return self.xyz, deflected_xyz

    def get_vector_result_by_scale_phase(self, itime: int, res_name: str,
                                         scale: float,
                                         phase: float=0.) -> tuple[np.ndarray, np.ndarray]:
        """
        Gets the real/complex deflection result

        Parameters
        ----------
        i : int
            mode/time/loadstep number
        name : str
            unused; useful for debugging
        scale : float
            deflection scale factor; true scale
        phase : float; default=0.0
            phase angle (degrees); unused for real results

        Returns
        -------
        xyz : (nnodes, 3) float ndarray
            the nominal state
        deflected_xyz : (nnodes, 3) float ndarray
            the deflected state
        """
        assert self.dim == 3, self.dim
        assert len(self.xyz.shape) == 2, self.xyz.shape
        if self.is_real:
            dxyz, *unused_junk = self.get_vector_data_dense(itime, res_name)
            deflected_xyz = self.xyz + scale * dxyz[itime, :]
        else:
            assert isinstance(itime, int), (itime, phase)
            assert isinstance(phase, float), (itime, phase)
            dxyz = self._get_complex_displacements_by_phase(itime, phase)
            deflected_xyz = self.xyz + scale * dxyz
        assert len(deflected_xyz.shape) == 2, deflected_xyz.shape
        return self.xyz, deflected_xyz

    def __repr__(self) -> str:
        """defines str(self)"""
        msg = 'DisplacementResults2\n'
        #msg += f'    titles={self.titles!r}\n'
        msg += f'    subcase_id={self.subcase_id}\n'
        msg += f'    data_type={self.data_type!r}\n'
        msg += f'    is_real={self.is_real} is_complex={self.is_complex}\n'
        msg += f'    location={self.location!r}\n'
        msg += f'    header={self.headers!r}\n'
        msg += f'    data_format={self.data_formats!r}\n'
        msg += f'    uname={self.uname!r}\n'
        return msg
