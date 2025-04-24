from __future__ import annotations
from typing import TYPE_CHECKING
import numpy as np

#from pyNastran.femutils.utils import safe_norm
from pyNastran.gui.gui_objects.vector_results import (
    DispForceVectorResults)  # _to_dense_vector

if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.op2.result_objects.table_object import (
        RealTableArray, ComplexTableArray)

translation = ['Magnitude', 'Fx', 'Fy', 'Fz']
rotation = ['Magnitude', 'Mx', 'My', 'Mz']


class ForceResults2(DispForceVectorResults):
    def __init__(self,
                 subcase_id: int,
                 node_id: np.ndarray,
                 xyz: np.ndarray,
                 case: RealTableArray | ComplexTableArray,
                 title: str,
                 t123_offset: int,
                 methods_txyz_rxyz: list[str],
                 index_to_base_title_annotation:dict[int, tuple[str, str]],
                 dim_max: float=1.0,
                 data_format: str='%g',
                 is_variable_data_format: bool=False,
                 nlabels=None, labelsize=None, ncolors=None,
                 colormap: str='',
                 set_max_min: bool=False,
                 uname: str='ForceResults2'):
        """
        Defines a SPC Force/MPC Force/Applied Load result

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
        dxyz : (nnodes, 3)
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
        uname : str
            some unique name for ...
        """
        self.methods_txyz_rxyz = methods_txyz_rxyz
        DispForceVectorResults.__init__(
            self,
            subcase_id,
            title,
            node_id,
            case,
            t123_offset,
            methods_txyz_rxyz,
            index_to_base_title_annotation,
            dim_max,
            data_format=data_format,
            is_variable_data_format=is_variable_data_format,
            nlabels=nlabels, labelsize=labelsize, ncolors=ncolors,
            colormap=colormap,
            set_max_min=set_max_min,
            uname=uname)

        #linked_scale_factor = False
        #location = 'node'

        # setup the node mapping
        disp_nodes = case.node_gridtype[:, 0]  #  local node id
        self.common_nodes = np.intersect1d(node_id, disp_nodes)
        self.inode_common = np.searchsorted(node_id, self.common_nodes)
        self.inode_result = np.searchsorted(disp_nodes, self.common_nodes)
        assert disp_nodes.max() > 0, disp_nodes
        assert len(self.inode_result) > 0, self.inode_result

        # dense -> no missing nodes in the results set
        self.is_dense = (len(node_id) == len(disp_nodes))

        self.xyz = xyz
        assert len(self.xyz.shape) == 2, self.xyz.shape
        assert isinstance(t123_offset, int), t123_offset
        self.location = 'node'
        str(self)

    #-------------------------------------
    def _calculate_scale(self, i: int, resname: str) -> float:
        return 1.0

    def has_output_checks(self, i: int, resname: str) -> tuple[bool, bool, bool,
                                                               bool, bool, bool]:
        is_enabled_fringe = True
        is_checked_fringe = True
        is_enabled_disp = False
        is_checked_disp = False
        is_enabled_vector = True
        is_checked_vector = True
        out = (
            is_enabled_fringe, is_checked_fringe,
            is_enabled_disp, is_checked_disp,
            is_enabled_vector, is_checked_vector)
        return out

    def get_force_vector_result(self, itime: int, res_name: str,
                                ) -> tuple[np.ndarray, np.ndarray]:
        dxyz, *unused_junk = self.get_vector_data_dense(itime, res_name)
        #scale = 1.
        assert dxyz.ndim == 2, dxyz.shape
        return self.xyz, dxyz # * scale

    def get_vector_result(self, itime: int, res_name: str,
                          ) -> tuple[np.ndarray, np.ndarray]:
        scale = self.get_scale(itime, res_name)
        if self.is_real:
            dxyz, *unused_junk = self.get_vector_data_dense(itime, res_name)
            assert dxyz.ndim == 2, dxyz.shape
            xyz = self.xyz
            deflected_xyz = dxyz * scale
            #deflected_xyz = self.xyz + scale * dxyz
            #return self.xyz, deflected_xyz
        else:
            phase = self.get_phase(itime, res_name)
            xyz, deflected_xyz = self.get_vector_result_by_scale_phase(
                itime, res_name, scale, phase)
        return xyz, deflected_xyz

    def _get_complex_displacements_by_phase(self, itime: int, res_name: str,
                                            phase: float=0.) -> np.ndarray:
        """
        Get force for a complex result.

        e^(i*theta) = cos(theta) + 1j*sin(theta)
        """
        dxyz, *unused_junk = self.get_vector_data_dense(itime, res_name)
        assert dxyz.ndim == 2, dxyz

        theta = np.radians(phase)
        if self.is_real:
            dxyz_out = dxyz * np.cos(theta)
        else:
            dxyz_out = dxyz.real * np.cos(theta) + dxyz.imag * np.sin(theta)
        return dxyz_out

    def get_vector_result_by_scale_phase(self, itime: int, res_name: str,
                                         scale: float,
                                         phase: float=0.) -> tuple[np.ndarray, np.ndarray]:
        """
        Gets the real/complex deflection result

        Parameters
        ----------
        itime : int
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
        assert isinstance(itime, int), (itime, phase)
        assert self.dim == 3, self.dim
        assert len(self.xyz.shape) == 2, self.xyz.shape
        if self.is_real:
            dxyz, itime, case_flag = self.get_vector_data_dense(itime, res_name)
            deflected_xyz = self.xyz + scale * dxyz
        else:
            assert isinstance(phase, float), (itime, phase)
            dxyz = self._get_complex_displacements_by_phase(
                itime, res_name, phase)
            deflected_xyz = self.xyz + scale * dxyz.real
        assert len(deflected_xyz.shape) == 2, deflected_xyz.shape
        return self.xyz, deflected_xyz

    def __repr__(self) -> str:
        """defines str(self)"""
        msg = 'ForceResults2\n'
        #msg += f'    titles={self.titles!r}\n'
        msg += f'    subcase_id={self.subcase_id}\n'
        msg += f'    data_type={self.data_type!r}\n'
        msg += f'    is_real={self.is_real} is_complex={self.is_complex}\n'
        msg += f'    location={self.location!r}\n'
        msg += f'    header={self.headers!r}\n'
        msg += f'    data_format={self.data_formats!r}\n'
        msg += f'    uname={self.uname!r}\n'
        return msg
