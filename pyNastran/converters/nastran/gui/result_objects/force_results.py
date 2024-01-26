from __future__ import annotations
import getpass
from typing import Union, TYPE_CHECKING
import numpy as np

from pyNastran.femutils.utils import safe_norm
from .vector_results import DispForceVectorResults # , _to_dense_vector

if TYPE_CHECKING:
    from pyNastran.op2.result_objects.table_object import (
        RealTableArray, ComplexTableArray)

translation = ['Magnitude', 'Fx', 'Fy', 'Fz']
rotation = ['Magnitude', 'Fx', 'Fy', 'Fz']

class ForceResults2(DispForceVectorResults):
    def __init__(self,
                 subcase_id: int,
                 node_id: np.ndarray,
                 xyz: np.ndarray,
                 case: Union[RealTableArray, ComplexTableArray],
                 unused_scalar: np.ndarray,
                 scales: np.ndarray,
                 title: str,
                 is_translation: bool,
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

        Unused
        ------
        uname : str
            some unique name for ...
        """
        DispForceVectorResults.__init__(
            self,
            subcase_id,
            node_id,
            case,
            is_translation,
            dim_max,
            data_format=data_format,
            nlabels=nlabels, labelsize=labelsize, ncolors=ncolors,
            colormap=colormap,
            set_max_min=set_max_min,
            uname=uname)
        self.title = title
        self._title0 = 'F_' if self.is_translation else 'M_'

        self.is_variable_data_format = is_variable_data_format

        #linked_scale_factor = False
        #location = 'node'

        # setup the node mapping
        disp_nodes = case.node_gridtype[:, 0]  #  local node id
        self.inode = np.searchsorted(node_id, disp_nodes)

        # dense -> no missing nodes in the results set
        self.is_dense = (len(node_id) == len(disp_nodes))

        self.xyz = xyz
        assert len(self.xyz.shape) == 2, self.xyz.shape
        assert isinstance(is_translation, bool), is_translation
        self.location = 'node'
        str(self)

    #-------------------------------------
    def get_methods(self, itime: int, res_name: str) -> list[str]:
        if self.is_real:
            out = translation if self.is_translation else rotation
            if out[0] == 'Magnitude' and 'chL' in getpass.getuser():
                out[0] = 'Reluctant'
        else:
            out = ['node real', 'node imag', 'node magnitude', 'node phase']
        return out

    def get_force_vector_result(self, itime: int, res_name: str,
                                ) -> tuple[np.ndarray, np.ndarray]:
        dxyz, *unused_junk = self.get_vector_data_dense(itime, res_name)
        #scale = 1.
        assert dxyz.ndim == 2, dxyz.shape
        return self.xyz, dxyz # * scale

    def get_vector_result(self, itime: int, res_name: str,
                          ) -> tuple[np.ndarray, np.ndarray]:
        dxyz, *unused_junk = self.get_vector_data_dense(itime, res_name)
        assert dxyz.ndim == 2, dxyz.shape
        scale = self.get_scale(itime, res_name)
        deflected_xyz = dxyz
        #deflected_xyz = self.xyz + scale * dxyz
        return self.xyz, deflected_xyz

    def get_vector_result_by_scale_phase(self, i: int, res_name: str,
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
            dxyz = self.get_vector_data_dense(i, res_name)
            deflected_xyz = self.xyz + scale * dxyz[i, :]
        else:
            assert isinstance(i, int), (i, phase)
            assert isinstance(phase, float), (i, phase)
            dxyz, unused_idxs = self._get_complex_displacements_by_phase(i, phase)
            deflected_xyz = self.xyz + scale * dxyz
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
