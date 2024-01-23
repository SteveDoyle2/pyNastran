import getpass
from collections import defaultdict
import numpy as np
from typing import Union, Optional, Any # , TYPE_CHECKING

from pyNastran.femutils.utils import safe_norm
#if TYPE_CHECKING:
from pyNastran.op2.result_objects.table_object import (
    RealTableArray, ComplexTableArray)
from pyNastran.gui.gui_objects.gui_result import GuiResultCommon


translation = ['Magnitude', 'Tx', 'Ty', 'Tz']
rotation = ['Magnitude', 'Rx', 'Ry', 'Rz']
col_axis = 1

class VectorResultsCommon(GuiResultCommon):
    def __init__(self,
                 subcase_id: int,
                 #node_id: np.ndarray,
                 case,
                 #dxyz: Union[RealTableArray, ComplexTableArray],
                 #is_translation: bool,
                 #dim_max: float,
                 data_format: str='%g',
                 nlabels=None, labelsize=None, ncolors=None,
                 colormap: str='',
                 #set_max_min: bool=False,
                 uname: str='VectorResults2'):
        GuiResultCommon.__init__(self)
        self.is_dense = False
        self.subcase_id = subcase_id

        # local case object
        self.case = case

        self.linked_scale_factor = False
        ntimes = case.data.shape[0]
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

        self.data_format = data_format
        self.data_formats = [self.data_format]
        self.headers = ['headers'] * ntimes

        self.nlabels = None
        self.labelsize = None
        self.ncolors = None
        self.colormap = colormap

        self.uname = uname
        self.location = ''

    # --------------------------------------------------------------------------
    # unmodifyable getters
    def set_min_max(self, i: int, res_name: str,
                    min_value: float, max_value: float):
        raise TypeError((i, res_name))

    def get_location(self, unused_i: int, unused_res_name: str) -> str:
        """the result location (node/centroid)"""
        assert self.location, self.location
        return self.location

    # --------------------------------------------------------------------------
    def get_default_scale(self, itime: int, res_name: str) -> float:
        if self.linked_scale_factor:
            itime = 0
        scales = self.default_scales[self.component_indices]
        if scales[itime] is not None:
            return scales[itime]

        unused_ntimes, nnodes = self.case.data.shape[:2]
        datai, idxs = self._get_real_data(itime)
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

        scales = self.scales[self.component_indices]
        scale = scales[itime]
        if scale is None:
            scale = self.get_default_scale(itime, res_name)
            scales[itime] = scale
        return scale
    def set_scale(self, itime: int, res_name: str, scale: float) -> None:
        if self.linked_scale_factor:
            itime = 0
        scales = self.scales[self.component_indices]
        scales[itime] = scale

    def get_default_phase(self, itime: int, res_name: str) -> float:
        return 0.0
    def get_phase(self, itime: int, res_name: str) -> int:
        if self.is_real:
            return 0.0
        phases = self.phases[self.component_indices]
        return phases[itime]
    def set_phase(self, itime: int, res_name: str, phase: float) -> None:
        if self.is_real:
            return
        phases = self.phases[self.component_indices]
        phases[itime] = phase

    def get_data_format(self, itime: int, res_name: str) -> str:
        """TODO: currently single value for all results"""
        return self.data_format
    def set_data_format(self, itime: int, res_name: str, data_format: str) -> None:
        self.data_format = data_format
    def get_default_data_format(self, itime: int, res_name: str) -> str:
        return self.data_format

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


class DispForceVectorResults(VectorResultsCommon):
    def __init__(self,
                 subcase_id: int,
                 node_id: np.ndarray,
                 case: Union[RealTableArray, ComplexTableArray],
                 is_translation: bool,
                 dim_max: float,
                 data_format: str='%g',
                 nlabels=None, labelsize=None, ncolors=None,
                 colormap: str='',
                 set_max_min: bool=False,
                 uname: str='VectorResults2'):
        VectorResultsCommon.__init__(
            self, subcase_id, case, data_format,
            nlabels, labelsize, ncolors,
            colormap, uname)

        self.dim = case.data.ndim
        assert self.dim == 3, case.data.shape
        assert case.data.shape[2] == 6, case.data.shape

        self.is_translation = is_translation
        self._title0 = 'T_' if self.is_translation else 'R_'

        if dim_max == 0.0:
            dim_max = 1.0
        self.dim_max = dim_max

        #  global node ids
        self.node_id = node_id

        self.data_type = self.case.data.dtype.str # '<c8', '<f4'
        self.is_real = True if self.data_type in ['<f4', '<f8'] else False
        #self.is_real = case.data.dtype.name in {'float32', 'float64'}
        self.is_complex = not self.is_real

        if self.is_real:
            self.component_indices = (0, 1, 2)
            self.index_map = {0: 'X', 1: 'Y', 2: 'Z',}
        else:
            self.component_indices = (7, 8, 9, 10, 11, 12)
            self.index_map = {
                0: 'Magnitude',
                1: 'X Mag', 2: 'Y Mag', 3: 'Z Mag',
                4: 'X Phase', 5: 'Y Phase', 6: 'Z Phase',
                7: 'X Real', 8: 'Y Real', 9: 'Z Real',
                10: 'X Imaginary', 11: 'Y Imaginary', 12: 'Z Imaginary',
            }
        self.min_max_method = 'Magnitude'

    def set_sidebar_args(self,
                         itime: str, res_name: str,
                         min_max_method: str='',
                         transform: str='',
                         methods_keys: Optional[list[int]]=None,
                         # unused
                         nodal_combine: str='',
                         **kwargs) -> None:
        assert len(kwargs) == 0, kwargs
        transforms = self.has_coord_transform(itime, res_name)[1]
        min_max_methods = self.has_derivation_transform(itime, res_name)[1]['derivation']
        transform = transform if transform else transforms[0]
        min_max_method = min_max_method if min_max_method else min_max_methods[0]
        #sidebar_kwargs = {
            #'min_max_method': min_max_method,
            #'transform': coord,
            #'nodal_combine': nodal_combine,
            #'methods_keys': keys_b,
        #}
        assert transform in transforms, transform

        if self.is_real:
            component_indices = _get_real_component_indices(methods_keys)
        else:
            component_indices = _get_complex_component_indices(methods_keys)
        self.component_indices = component_indices

        # handle confusing data (can't plot a scalar for 2 components)
        if len(self.component_indices) > 1 and min_max_method == 'Value':
            min_max_method = 'Magnitude'

        # doesn't matter cause it's already nodal
        #nodal_combine
        #min_max_method
        self.transform = transform
        self.min_max_method = min_max_method
        str(self)

    def has_methods_table(self, i: int, res_name: str) -> bool:
        return True
    def has_coord_transform(self, i: int, res_name: str) -> tuple[bool, list[str]]:
        return True, ['Global']
    def has_derivation_transform(self, i: int, res_name: str,
                                 ) -> tuple[bool, dict[str, Any]]:
        """min/max/avg"""
        out = {
            'label': 'Derivation',
            'derivation': ['Magnitude', 'Value'],
            'tooltip': 'Magnitude is automatically selected if multiple cmponents are selected',
        }
        return True, out
    def has_nodal_combine_transform(self, i: int, res_name: str) -> tuple[bool, list[str]]:
        """elemental -> nodal"""
        return False, []
    def has_output_checks(self, i: int, resname: str) -> tuple[bool, bool, bool]:
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

    def get_vector_size(self, itime: int, res_name: str) -> int:
        """vector_size=1 is the default and displacement has 3 components"""
        return 3

    def get_annotation(self, itime: int, res_name: str) -> str:
        """
        A header is the thingy that goes in the lower left corner
        header = 'Static'
        title = 'Displacement'
        method = 'T_XYZ'
        returns 'Displacement T_XYZ: Static'
        """
        #method = self.get_methods(itime, res_name)[0]
        self.component_indices
        method = self._title0 + ''.join(self.index_map[idx] for idx in self.component_indices)
        #annotation_label = f'{self.title} {method} ({self.min_max_method}): {self.headers[itime]}'
        annotation_label = f'{self.title} {method}: {self.headers[itime]}'
        #return self.uname
        return annotation_label

    def get_default_min_max(self, itime: int, res_name: str) -> tuple[float, float]:
        return None, None
    def get_case_flag(self, itime: int, res_name: str):
        case_flag = (self.component_indices, self.min_max_method)
        return itime, case_flag
    def get_min_max(self, itime, res_name) -> tuple[float, float]:
        itime, case_flag = self.get_case_flag(itime, res_name)
        mins = self.mins[case_flag]
        maxs = self.maxs[case_flag]
        if mins[itime] is not None and maxs[itime] is not None:
            return mins[itime], maxs[itime]

        data = self.get_result(itime, res_name, method='', return_dense=True)

        normi = self._get_normi(data)
        if mins[itime] is None:
            mins[itime] = np.nanmin(normi)
        if maxs[itime] is None:
            maxs[itime] = np.nanmax(normi)
        return mins[itime], maxs[itime]

    def _get_normi(self, data: np.ndarray) -> np.ndarray:
        if self.min_max_method == 'Value' and len(self.component_indices) == 1:
            tnorm_abs = data
        else:
            self.min_max_method = 'Magnitude'
            tnorm_abs = safe_norm(data, axis=col_axis)
        return tnorm_abs

    def get_default_legend_title(self, itime: int, res_name: str) -> str:
        self.component_indices
        method = self._title0 + ''.join(self.index_map[idx] for idx in self.component_indices)
        title = f'{self.title} {method}'
        return title
    def set_legend_title(self, itime: int, res_name: str,
                         title: str) -> None:
        self.title = title
    def get_legend_title(self, itime: int, res_name: str):
        """displacement T_XYZ"""
        #method2 = self._update_method(method)
        #return f'{self.title} {method2}'
        self.component_indices
        method = self._title0 + ''.join(self.index_map[idx] for idx in self.component_indices)
        title = f'{self.title} {method}'
        return title

    def _get_data(self, itime: int) -> np.ndarray:
        assert self.case.data.shape[2] == 6
        if self.is_translation:
            datai = self.case.data[itime, :, :3].copy()
            assert datai.shape[1] == 3, datai.shape
        else:
            datai = self.case.data[itime, :, 3:].copy()
            assert datai.shape[1] == 3, datai.shape
        assert len(self.component_indices) > 0, self.component_indices
        return datai


    def _get_real_data(self, itime: int) -> tuple[np.ndarray, list[int]]:
        """
        Gets the real data

        Returns
        -------
        data : (3, ) float np.ndarray
            real/complex ndarray
        idxs : list[int]
            the indicies that were used

        """
        data = self._get_data(itime)
        idxs = []
        if self.is_real:
            for idx in [0, 1, 2]:
                if idx not in self.component_indices:
                    data[:, idx] = 0.
                    idxs.append(idx)
        else:
            # real
            for idx in [7, 8, 9]:
                if idx not in self.component_indices:
                    data[:, idx].real = 0.
                    idxs.append(idx)
            # complex
            for idx in [10, 11, 12]:
                if idx not in self.component_indices:
                    data[:, idx].imag = 0.
                    idxs.append(idx)
        return data, idxs

    def _get_complex_data(self, itime: int) -> tuple[np.ndarray, list[int]]:
        return self._get_real_data(itime)
        #if self.is_translation:
            #datai = self.case.data[itime, :, :3]
            #assert datai.shape[1] == 3, datai.shape
        #else:
            #datai = self.case.data[itime, :, 3:]
            #assert datai.shape[1] == 3, datai.shape
        #return datai

    def _get_vector_result(self, itime: int, res_name: str,
                           return_dense: bool=True) -> tuple[np.ndarray, list[int]]:
        """
        gets the 'typical' result as a vector
         - DisplacementResults: vector; (n, 3) array

        Parameters
        ----------
        return_dense: bool
            Rreturns the data array in a way that the gui can use.
            Handles the null result case (e.g; SPC forces only
            at the SPC location).

        """
        if self.is_real:
            # multiple results
            # .0006 -> 0.0
            # .057 -> 0.0123
            # min
            dxyz, idxs = self._get_real_data(itime)
            #dxyz = data[itime, :, :]
        else:
            dxyz, idxs = self._get_complex_data(itime)
        assert len(dxyz.shape) == 2, dxyz.shape
        return dxyz, idxs

    def get_result(self, itime: int, res_name: str,
                   method: str='',
                   return_dense: bool=True) -> np.ndarray:
        """
        gets the 'typical' result which is a vector
         - GuiResult:           fringe; (n,)   array

        Parameters
        ----------
        return_dense: bool
            Rreturns the data array in a way that the gui can use.
            Handles the null result case (e.g; SPC forces only
            at the SPC location).

        """
        itime, case_flag = self.get_case_flag(itime, res_name)
        dxyz, idxs = self._get_vector_result(itime, res_name, return_dense=False)
        return_sparse = not return_dense
        if (return_sparse or self.is_dense) and len(self.component_indices) > 1:
            return dxyz

        # apply dense
        nnodes = len(self.node_id)
        if self.min_max_method == 'Value' and len(self.component_indices) == 1:
            dxyz2 = np.full(nnodes, np.nan, dtype=dxyz.dtype)
            dxyz2[self.inode] = dxyz[:, self.component_indices[0]]
        else:
            dxyz2 = np.full((nnodes, 3), np.nan, dtype=dxyz.dtype)
            dxyz2[self.inode, :] = dxyz

        # set the min/max if they're not set
        normi = self._get_normi(dxyz2)
        mins = self.mins[case_flag]
        maxs = self.maxs[case_flag]
        if mins[itime] is None:
            mins[itime] = np.nanmin(normi)
        if maxs[itime] is None:
            maxs[itime] = np.nanmax(normi)
        return normi


def _get_real_component_indices(methods_keys: list[int]) -> np.ndarray:
    """
    if Magnitude is selected, only use magnitude
    methods = ['T_XYZ', 'TX', 'TY', 'TZ']
    method_keys = [0,    1      2     3]
    """
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
    component_indices = tuple(np.array(indices, dtype='int32') - 1)
    return component_indices

def _get_complex_component_indices(methods_keys: list[str]) -> np.ndarray:
    """
    if Magnitude is selected, only use magnitude
    methods = [
        'Resultant',
        'X Magnitude', 'Y Magnitude', 'Z Magnitude',
        'X Phase', 'Y Phase', 'Z Phase',
        'X Real', 'Y Real', 'Z Real',
        'X Imaginary', 'Y Imaginary', 'Z Imaginary',
    ]
    """
    # the data is real/imag, so we use that...
    #default_indices = [1, 2, 3, 4, 5, 6]    # mag/phase
    default_indices = [7, 8, 9, 10, 11, 12]  # real/imag
    if methods_keys is None or len(methods_keys) == 0:
        # default
        indices = default_indices
    elif 0 in methods_keys:
        # include all components b/c Magnitude is selected
        indices = default_indices
    else:
        # no 0 (magnitude) in methods_keys
        # update the indices to correspond to the array
        #
        # first split names by type
        for key in methods_keys:
            assert isinstance(key, str), key
            assert isinstance(key, int), key
        indices = methods_keys
    component_indices = tuple(np.array(indices, dtype='int32') - 1)
    return component_indices


class DisplacementResults2(DispForceVectorResults):
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
        DispForceVectorResults.__init__(
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
        """deflection is opt-in"""
        return True

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
            out = [
                # if results of different type are selected (e.g., Real/Imag)
                # the "earliest" type is selected
                # Resultant -> Magnitude -> Phase -> Real -> Imaginary
                'Resultant',
                'X Magnitude', 'Y Magnitude', 'Z Magnitude',
                'X Phase', 'Y Phase', 'Z Phase',
                'X Real', 'Y Real', 'Z Real',
                'X Imaginary', 'Y Imaginary', 'Z Imaginary',
            ]
        return out

    def get_plot_value(self, itime: int, res_name: str, method: str) -> np.ndarray:
        """get_fringe_value"""
        dxyz = self.get_result(itime, res_name, method, return_dense=False)
        normi = safe_norm(dxyz, axis=col_axis)
        if self.is_dense:
            return normi
        #case.data.shape = (11, 43, 6)
        #nnodes = len(self.node_id) =  48
        #nnodesi = len(self.inode) = len(self.case.node_gridtype) = 43
        normi2 = np.full(len(self.node_id), np.nan, dtype=normi.dtype)
        normi2[self.inode] = normi
        return normi2

    def get_force_vector_result(self, itime: int, res_name: str, method: str) -> np.ndarray:
        dxyz = self._get_vector_result(itime, res_name, method, return_dense=True)
        scale = 1.
        return self.xyz, dxyz * scale

    def get_vector_result(self, itime: int, res_name: str, method: str,
                          return_dense: bool=True) -> tuple[np.ndarray, np.ndarray]:
        #dxyz = self.get_result(itime, res_name, method, return_dense=True)

        dxyzi, idxs = self._get_vector_result(itime, res_name, return_dense=False)

        # map vectored data to dense
        nnodes = len(self.node_id)
        return_sparse = not return_dense
        if (return_sparse or self.is_dense):
            dxyz = dxyzi
        else:
            dxyz = np.full((nnodes, 3), np.nan, dtype=dxyzi.dtype)
            dxyz[self.inode, :] = dxyzi

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


class ForceResults2(DispForceVectorResults):
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
        DispForceVectorResults.__init__(
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
        self._title0 = 'F_' if self.is_translation else 'M_'

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
        #nnodesi = len(self.inode) = len(self.case.node_gridtype) = 43
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
