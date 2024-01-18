import getpass
from collections import defaultdict
import numpy as np
from typing import Union, Optional # , TYPE_CHECKING

from pyNastran.femutils.utils import safe_norm
#if TYPE_CHECKING:
from pyNastran.op2.result_objects.table_object import (
    RealTableArray, ComplexTableArray)
from pyNastran.gui.gui_objects.gui_result import GuiResultCommon

#def safe_norm(data: np.ndarray, axis=None):
    #return np.linalg.norm(data, axis=axis)


translation = ['Magnitude', 'Tx', 'Ty', 'Tz']
rotation = ['Magnitude', 'Rx', 'Ry', 'Rz']
col_axis = 1

class VectorResults2(GuiResultCommon):
    def __init__(self,
                 subcase_id: int,
                 node_id: np.ndarray,
                 dxyz: Union[RealTableArray, ComplexTableArray],
                 is_translation: bool,
                 dim_max: float,
                 data_format: str='%g',
                 nlabels=None, labelsize=None, ncolors=None,
                 colormap: str='',
                 set_max_min: bool=False,
                 uname: str='VectorResults2'):
        GuiResultCommon.__init__(self)
        self.method_indices = (0, 1, 2)
        self.is_dense = False
        self.dim = dxyz.data.ndim
        assert self.dim == 3, dxyz.data.shape

        self.subcase_id = subcase_id
        self.is_translation = is_translation

        if dim_max == 0.0:
            dim_max = 1.0
        self.dim_max = dim_max
        self.linked_scale_factor = False

        self.data_format = data_format

        #  global node ids
        self.node_id = node_id

        # local case object
        self.dxyz = dxyz

        self.data_type = self.dxyz.data.dtype.str # '<c8', '<f4'
        self.is_real = True if self.data_type in ['<f4', '<f8'] else False
        #self.is_real = dxyz.data.dtype.name in {'float32', 'float64'}
        self.is_complex = not self.is_real

        ntimes = self.dxyz.data.shape[0]
        nscale = ntimes
        if self.linked_scale_factor:
            nscale = 1

        def fscales():
            return [None] * nscale
        def ftimes():
            return [None] * ntimes
        def fphases():
            return np.zeros(ntimes, dtype='float64')

        self.default_scales = defaultdict(fscales)
        self.scales = defaultdict(fscales)
        self.mins = defaultdict(ftimes)
        self.maxs = defaultdict(ftimes)
        self.phases = defaultdict(fphases)

        self.data_formats = [self.data_format]
        self.headers = ['headers'] * ntimes

        self.nlabels = None
        self.labelsize = None
        self.ncolors = None
        self.colormap = colormap

        self.uname = uname
        self.active_method = ''

    #def get_header(self, itime: int, name: str) -> int:
        #return 3

    def set_sidebar_args(self,
                         itime: str, name: str,
                         min_max_method: str='',
                         transform: str='',
                         methods_keys: Optional[list[int]]=None,
                         # unused
                         nodal_combine: str='',
                         **kwargs) -> None:
        assert len(kwargs) == 0, kwargs
        #sidebar_kwargs = {
            #'min_max_method': min_max_method,
            #'transform': coord,
            #'nodal_combine': nodal_combine,
            #'methods_keys': keys_b,
        #}
        assert transform in {'', 'Coord 0'}, transform

        # if Magnitude is selected, only use magnitude
        # methods = ['T_XYZ', 'TX', 'TY', 'TZ']
        default_indices = [1, 2, 3]
        if methods_keys is None or len(methods_keys) == 0:
            # default
            indices = default_indices
        elif 0 in methods_keys:
            # include all components b/c Magnitude is selected
            indices = default_indices
        else:
            # no 0 (magnitude) in methods_keys
            # update the indices to correspond to the array
            #methods_keys.sort()
            indices = methods_keys
        self.method_indices = tuple(np.array(indices, dtype='int32') - 1)

        # doesn't matter cause it's already nodal
        #nodal_combine
        #min_max_method
        self.transform = transform
        x = 1


    def has_methods_table(self, i: int, res_name: str) -> bool:
        return True
    def has_coord_transform(self, i: int, res_name: str) -> bool:
        return True
    def has_derivation_transform(self, i: int, res_name: str) -> bool:  # min/max/avg
        return False
    def has_nodal_combine_transform(self, i: int, res_name: str) -> bool:  # elemental -> nodal
        return False

    def get_vector_size(self, itime: int, res_name: str) -> int:
        return 3
    def is_normal_result(self, itime: int, res_name: str) -> int:
        return False

    def get_header(self, itime: int, res_name: str) -> str:
        """
        A header is the thingy that goes in the lower left corner
        header = 'Static'
        title = 'Displacement'
        method = 'T_XYZ'
        returns 'Displacement T_XYZ: Static'
        """
        #method = self.get_methods(itime, res_name)[0]
        method_ = 'T_' if self.is_translation else 'R_'
        index_map = {0: 'X', 1: 'Y', 2: 'Z',}
        self.method_indices
        method = method_ + ''.join(index_map[idx] for idx in self.method_indices)
        annotation_label = f'{self.title} {method}: {self.headers[itime]}'
        #return self.uname
        return annotation_label

    def get_default_scale(self, itime: int, res_name: str) -> float:
        if self.linked_scale_factor:
            itime = 0
        scales = self.default_scales[self.method_indices]
        if scales[itime] is not None:
            return scales[itime]

        unused_ntimes, nnodes = self.dxyz.data.shape[:2]
        datai = self._get_real_data(itime)
        tnorm_abs = safe_norm(datai, axis=col_axis)
        assert len(tnorm_abs) == nnodes
        tnorm_abs_max = np.nanmax(tnorm_abs)
        if tnorm_abs_max > 0.0:
            scale = self.dim_max / tnorm_abs_max * 0.10
        else:
            scale = 1.0
        scales[itime] = scale
        return scale
    def get_scale(self, itime: int, res_name: str) -> int:
        if self.linked_scale_factor:
            itime = 0

        scales = self.scales[self.method_indices]
        scale = scales[itime]
        if scale is None:
            scale = self.get_default_scale(itime, res_name)
            scales[itime] = scale
        return scale
    def set_scale(self, itime: int, res_name: str, scale: float) -> None:
        if self.linked_scale_factor:
            itime = 0
        scales = self.scales[self.method_indices]
        scales[itime] = scale

    def get_default_phase(self, itime: int, res_name: str) -> float:
        return 0.0
    def get_phase(self, itime: int, res_name: str) -> int:
        if self.is_real:
            return 0.0
        phases = self.phases[self.method_indices]
        return phases[itime]
    def set_phase(self, itime: int, res_name: str, phase: float) -> None:
        if self.is_real:
            return
        phases = self.phases[self.method_indices]
        phases[itime] = phase

    def get_data_format(self, itime: int, res_name: str) -> str:
        """TODO: currently single value for all results"""
        return self.data_format
    def set_data_format(self, itime: int, res_name: str, data_format: str) -> None:
        self.data_format = data_format
    def get_default_data_format(self, itime: int, res_name: str) -> str:
        return self.data_format

    def set_method(self, method: str):
        methods = self.get_methods()
        assert method in methods, (method, methods)
        self.active_method = method

    def get_default_min_max(self, itime: int, res_name: str,
                            method: str) -> tuple[float, float]:
        return None, None
    def get_min_max(self, itime, resname) -> tuple[float, float]:
        mins = self.mins[self.method_indices]
        maxs = self.maxs[self.method_indices]
        if mins[itime] is not None and maxs[itime] is not None:
            return mins[itime], maxs[itime]

        datai = self._get_real_data(itime)
        tnorm_abs = safe_norm(datai, axis=col_axis)
        if mins[itime] is None:
            mins[itime] = tnorm_abs.min()
        if maxs[itime] is None:
            maxs[itime] = tnorm_abs.max()
        return mins[itime], maxs[itime]

    def get_default_legend_title(self, itime: int, res_name: str) -> str:
        method_ = 'T_' if self.is_translation else 'R_'
        index_map = {0: 'X', 1: 'Y', 2: 'Z',}
        self.method_indices
        method = method_ + ''.join(index_map[idx] for idx in self.method_indices)
        title = f'{self.title} {method}'
        return title
    def set_legend_title(self, itime: int, res_name: str,
                         title: str) -> None:
        self.title = title
    def get_legend_title(self, itime: int, res_name: str):
        """displacement T_XYZ"""
        #method2 = self._update_method(method)
        #return f'{self.title} {method2}'
        method_ = 'T_' if self.is_translation else 'R_'
        index_map = {0: 'X', 1: 'Y', 2: 'Z',}
        self.method_indices
        method = method_ + ''.join(index_map[idx] for idx in self.method_indices)
        title = f'{self.title} {method}'
        return title


    def get_default_nlabels_labelsize_ncolors_colormap(self, itime: int, res_name: str,
                                                       ) -> tuple[None, None, None, str]:
        return None, None, None, ''
    def get_nlabels_labelsize_ncolors_colormap(self, i: int, res_name: str):
        return self.nlabels, self.labelsize, self.ncolors, self.colormap
    def set_nlabels_labelsize_ncolors_colormap(self, i: int, res_name: str,
                                               nlabels, labelsize,
                                               ncolors, colormap: str) -> None:
        self.nlabels = nlabels
        self.labelsize = labelsize
        self.ncolors = ncolors
        self.colormap = colormap
        pass

    def _get_real_data(self, itime: int) -> np.ndarray:
        if self.is_translation:
            datai = self.dxyz.data[itime, :, :3].copy()
            assert datai.shape[1] == 3, datai.shape
        else:
            datai = self.dxyz.data[itime, :, 3:].copy()
            assert datai.shape[1] == 3, datai.shape

        assert len(self.method_indices) > 0, self.method_indices
        for idx in [0, 1, 2]:
            if idx not in self.method_indices:
                datai[:, idx] = 0.
        if not self.is_real:
            asdf
        return datai

    def _get_complex_data(self, itime: int) -> np.ndarray:
        return self._get_real_data(itime)
        #if self.is_translation:
            #datai = self.dxyz.data[itime, :, :3]
            #assert datai.shape[1] == 3, datai.shape
        #else:
            #datai = self.dxyz.data[itime, :, 3:]
            #assert datai.shape[1] == 3, datai.shape
        #return datai

    def get_result(self, itime: int, res_name: str,
                   method: str='',
                   return_dense: bool=True) -> np.ndarray:
        """
        gets the 'typical' result which is a vector
         - GuiResult:           fringe; (n,)   array
         - DisplacementResults: vector; (n, 3) array

        Parameters
        ----------
        return_dense: bool
            Rreturns the data array in a way that the gui can use.
            Handles the null result case (e.g; SPC forces only
            at the SPC location).
        """
        method = self._update_method(itime, res_name, method)
        if self.is_real:
            # multiple results
            # .0006 -> 0.0
            # .057 -> 0.0123
            # min
            dxyz = self._get_real_data(itime)
            #dxyz = data[itime, :, :]
        else:
            dxyz = self._get_complex_data(itime)
        assert len(dxyz.shape) == 2, dxyz.shape

        return_sparse = not return_dense
        if return_sparse or self.is_dense:
            return dxyz

        nnodes = len(self.node_id)
        dxyz2 = np.full((nnodes, 3), np.nan, dtype=dxyz.dtype)
        dxyz2[self.inode, :] = dxyz
        return dxyz2

    def _update_method(self, i: int, res_name: str,
                       method: str) -> str:
        if method == '':
            method = 'Magnitude'
        self.active_method = method
        assert method in self.get_methods(i, res_name)
        return method
    #def _update_method(self, method: str) -> str:
        #assert method != '', method
        #self.active_method = method


class DisplacementResults2(VectorResults2):
    def __init__(self,
                 subcase_id: int,
                 node_id: np.ndarray,
                 xyz: np.ndarray,
                 dxyz: Union[RealTableArray, ComplexTableArray],
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
        VectorResults2.__init__(
            self,
            subcase_id,
            node_id,
            dxyz,
            is_translation,
            dim_max,
            data_format=data_format,
            nlabels=nlabels, labelsize=labelsize, ncolors=ncolors,
            colormap=colormap,
            set_max_min=set_max_min,
            uname=uname)
        self.title = title

        self.is_variable_data_format = is_variable_data_format

        #linked_scale_factor = False
        #location = 'node'

        # setup the node mapping
        disp_nodes = dxyz.node_gridtype[:, 0]  #  local node id
        self.inode = np.searchsorted(node_id, disp_nodes)

        # dense -> no missing nodes in the results set
        self.is_dense = (len(node_id) == len(disp_nodes))

        self.xyz = xyz
        assert len(self.xyz.shape) == 2, self.xyz.shape
        assert isinstance(is_translation, bool), is_translation
        self.location = 'node'
        str(self)

    #-------------------------------------
    # unmodifyable getters

    def deflects(self, unused_i: int, unused_res_name: str) -> bool:
        return True

    def get_location(self, unused_i: int, unused_res_name: str) -> str:
        """the result type"""
        return self.location

    def set_min_max(self, i: int, res_name: str,
                    min_value: float, max_value: float):
        raise TypeError((i, res_name))

    #-------------------------------------
    def get_methods(self, itime: int, res_name: str) -> list[str]:
        if self.is_real:
            out = translation if self.is_translation else rotation
            if out[0] == 'Magnitude' and 'chL' in getpass.getuser():
                out[0] = 'Reluctant'
        else:
            out = ['node real', 'node imag', 'node magnitude', 'node phase']
        return out

    def get_plot_value(self, itime: int, res_name: str, method: str) -> np.ndarray:
        """get_fringe_value"""
        dxyz = self.get_result(itime, res_name, method, return_dense=False)
        normi = safe_norm(dxyz, axis=col_axis)
        if self.is_dense:
            return normi

        #case.data.shape = (11, 43, 6)
        #nnodes = len(self.node_id) =  48
        #nnodesi = len(self.inode) = len(self.dxyz.node_gridtype) = 43
        normi2 = np.full(len(self.node_id), np.nan, dtype=normi.dtype)
        normi2[self.inode] = normi
        return normi2

    def get_force_vector_result(self, itime: int, res_name: str, method: str) -> np.ndarray:
        dxyz = self.get_result(itime, res_name, method, return_dense=True)
        scale = 1.
        return self.xyz, dxyz * scale

    def get_vector_result(self, itime: int, res_name: str, method: str) -> tuple[np.ndarray, np.ndarray]:
        dxyz = self.get_result(itime, res_name, method, return_dense=True)
        scale = self.get_scale(itime, res_name)
        deflected_xyz = self.xyz + scale * dxyz
        return self.xyz, deflected_xyz

    def get_vector_result_by_scale_phase(self, i: int, unused_name: str,
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
            deflected_xyz = self.xyz + scale * self.dxyz[i, :]
        else:
            assert isinstance(i, int), (i, phase)
            assert isinstance(phase, float), (i, phase)
            dxyz = self._get_complex_displacements_by_phase(i, phase)
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


class ForceResults2(VectorResults2):
    def __init__(self,
                 subcase_id: int,
                 node_id: np.ndarray,
                 xyz: np.ndarray,
                 dxyz: Union[RealTableArray, ComplexTableArray],
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
        VectorResults2.__init__(
            self,
            subcase_id,
            node_id,
            dxyz,
            is_translation,
            dim_max,
            data_format=data_format,
            nlabels=nlabels, labelsize=labelsize, ncolors=ncolors,
            colormap=colormap,
            set_max_min=set_max_min,
            uname=uname)
        self.title = title

        self.is_variable_data_format = is_variable_data_format

        #linked_scale_factor = False
        #location = 'node'

        # setup the node mapping
        disp_nodes = dxyz.node_gridtype[:, 0]  #  local node id
        self.inode = np.searchsorted(node_id, disp_nodes)

        # dense -> no missing nodes in the results set
        self.is_dense = (len(node_id) == len(disp_nodes))

        self.xyz = xyz
        assert len(self.xyz.shape) == 2, self.xyz.shape
        assert isinstance(is_translation, bool), is_translation
        self.location = 'node'
        str(self)

    #-------------------------------------
    # unmodifyable getters

    def deflects(self, unused_i: int, unused_res_name: str) -> bool:
        return False

    def get_location(self, unused_i: int, unused_res_name: str) -> str:
        """the result type"""
        return self.location

    def set_min_max(self, i: int, res_name: str,
                    min_value: float, max_value: float):
        raise TypeError((i, res_name))

    #-------------------------------------
    def get_methods(self, itime: int, res_name: str) -> list[str]:
        if self.is_real:
            out = translation if self.is_translation else rotation
            if out[0] == 'Magnitude' and 'chL' in getpass.getuser():
                out[0] = 'Reluctant'
        else:
            out = ['node real', 'node imag', 'node magnitude', 'node phase']
        return out

    def get_plot_value(self, itime: int, res_name: str, method: str) -> np.ndarray:
        """get_fringe_value"""
        dxyz = self.get_result(itime, res_name, method, return_dense=False)
        normi = safe_norm(dxyz, axis=col_axis)
        if self.is_dense:
            return normi

        #case.data.shape = (11, 43, 6)
        #nnodes = len(self.node_id) =  48
        #nnodesi = len(self.inode) = len(self.dxyz.node_gridtype) = 43
        normi2 = np.full(len(self.node_id), np.nan, dtype=normi.dtype)
        normi2[self.inode] = normi
        return normi2

    def get_force_vector_result(self, itime: int, res_name: str, method: str) -> np.ndarray:
        dxyz = self.get_result(itime, res_name, method, return_dense=True)
        scale = 1.
        return self.xyz, dxyz # * scale

    def get_vector_result(self, itime: int, res_name: str, method: str) -> tuple[np.ndarray, np.ndarray]:
        dxyz = self.get_result(itime, res_name, method, return_dense=True)
        scale = self.get_scale(itime, res_name)
        deflected_xyz = dxyz
        #deflected_xyz = self.xyz + scale * dxyz
        return self.xyz, deflected_xyz

    def get_vector_result_by_scale_phase(self, i: int, unused_name: str,
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
            deflected_xyz = self.xyz + scale * self.dxyz[i, :]
        else:
            assert isinstance(i, int), (i, phase)
            assert isinstance(phase, float), (i, phase)
            dxyz = self._get_complex_displacements_by_phase(i, phase)
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
