import getpass
from abc import abstractmethod
from collections import defaultdict

import numpy as np
from typing import DefaultDict, Union, Optional, Any # , TYPE_CHECKING

from pyNastran.femutils.utils import safe_norm
#if TYPE_CHECKING:
from pyNastran.op2.result_objects.table_object import (
    RealTableArray, ComplexTableArray)
from pyNastran.gui.gui_objects.gui_result import GuiResultCommon


translation = ['Magnitude', 'Tx', 'Ty', 'Tz']
rotation = ['Magnitude', 'Rx', 'Ry', 'Rz']
col_axis = 1

#if TYPE_CHECKING:
    #from typing import TypedDict
    #class SidebarArgs(TypedDict):
        #label: str
        #tooltip: int
        #derivation: list[str]


class VectorResultsCommon(GuiResultCommon):
    def __init__(self,
                 subcase_id: int,
                 case,
                 #dxyz: Union[RealTableArray, ComplexTableArray],
                 data_format: str='%g',
                 nlabels=None, labelsize=None, ncolors=None,
                 colormap: str='',
                 #set_max_min: bool=False,
                 uname: str='VectorResults2'):
        GuiResultCommon.__init__(self)
        self.is_dense = False
        self.subcase_id = subcase_id

        # local case object
        if isinstance(case, list):  # TODO: might change this later...
            self.cases = case
            case = case[0]
        else:
            self.case = case

        self.linked_scale_factor = False
        ntimes = case.data.shape[0]
        nscale = ntimes
        if self.linked_scale_factor:
            nscale = 1

        def fscales() -> list[float]:
            return [np.nan] * nscale
        def ftimes() -> list[float]:
            return [np.nan] * ntimes
        def fphases() -> np.ndarray:
            return np.zeros(ntimes, dtype='float64')

        self.default_arrow_scales: DefaultDict[Any, list[float]] = defaultdict(fscales)
        self.arrow_scales:         DefaultDict[Any, list[float]] = defaultdict(fscales)

        self.default_scales: DefaultDict[Any, list[float]] = defaultdict(fscales)
        self.scales:         DefaultDict[Any, list[float]] = defaultdict(fscales)

        self.default_mins: DefaultDict[Any, list[float]] = defaultdict(ftimes)
        self.default_maxs: DefaultDict[Any, list[float]] = defaultdict(ftimes)
        self.mins:         DefaultDict[Any, list[float]] = defaultdict(ftimes)
        self.maxs:         DefaultDict[Any, list[float]] = defaultdict(ftimes)

        self.phases: DefaultDict[Any, np.ndarray] = defaultdict(fphases)

        self.data_format = data_format
        self.data_formats = [self.data_format]
        self.headers = ['VectorResultsCommon'] * ntimes

        self.nlabels = None
        self.labelsize = None
        self.ncolors = None
        self.colormap = colormap

        self.uname = uname
        self.location = ''
        self.min_max_method = ''
        self.nodal_combine = ''

    # --------------------------------------------------------------------------
    # abstractmethods - results access
    @abstractmethod
    def _get_fringe_data_sparse(self, itime: int, res_name) -> np.ndarray:
        raise NotImplementedError(f'{self.class_name}._get_fringe_data_sparse')
    @abstractmethod
    def _get_fringe_data_sparse(self, itime: int, res_name) -> np.ndarray:
        raise RuntimeError(f'{self.class_name}._get_fringe_data_sparse')

    # --------------------------------------------------------------------------
    # abstractmethods
    @abstractmethod
    def get_case_flag(self, itime: int, res_name):
        raise NotImplementedError(f'{self.class_name}.get_case_flag')

    @abstractmethod
    def set_sidebar_args(self,
                         itime: int, res_name: str, **kwargs) -> None:
        raise NotImplementedError(f'{self.class_name}.set_sidebar_args')

    # --------------------------------------------------------------------------
    # abstractmethods - sidebar stuff
    @abstractmethod
    def has_methods_table(self, i: int, res_name: str) -> bool:
        raise NotImplementedError(f'{self.class_name}.has_methods_table')
    @abstractmethod
    def has_coord_transform(self, i: int, res_name: str) -> tuple[bool, list[str]]:
        raise NotImplementedError(f'{self.class_name}.has_coord_transform')
    @abstractmethod
    def has_derivation_transform(self, i: int, res_name: str,
                                 ) -> tuple[bool, dict[str, Any]]:
        raise NotImplementedError(f'{self.class_name}.has_derivation_transform')
    @abstractmethod
    def has_nodal_combine_transform(self, i: int, res_name: str) -> tuple[bool, list[str]]:
        """elemental -> nodal"""
        raise NotImplementedError(f'{self.class_name}.has_output_checks')
    @abstractmethod
    def has_output_checks(self, i: int, resname: str) -> tuple[bool, bool, bool,
                                                               bool, bool, bool]:
        raise NotImplementedError(f'{self.class_name}.has_output_checks')

    # --------------------------------------------------------------------------
    # unmodifyable getters
    def set_min_max(self, itime: int, res_name,
                    min_value: float, max_value: float) -> None:
        itime, case_flag = self.get_case_flag(itime, res_name)

        mins = self.mins[case_flag]
        maxs = self.maxs[case_flag]
        mins[itime] = min_value
        maxs[itime] = max_value

    def get_default_min_max(self, itime: int,
                            res_name) -> tuple[float, float]:
        """Calculates the default min/max dynamically"""
        itime, case_flag = self.get_case_flag(itime, res_name)
        mins = self.default_mins[case_flag]
        maxs = self.default_maxs[case_flag]
        if is_value(mins[itime]) and is_value(maxs[itime]):
            return mins[itime], maxs[itime]

        fringe_data = self._get_fringe_data_sparse(itime, res_name)

        # save the defaults for the next time
        mins[itime] = np.nanmin(fringe_data)
        maxs[itime] = np.nanmax(fringe_data)
        return mins[itime], maxs[itime]

    def get_min_max(self, itime: int, res_name) -> tuple[float, float]:
        itime, case_flag = self.get_case_flag(itime, res_name)
        mins = self.mins[case_flag]
        maxs = self.maxs[case_flag]
        if is_value(mins[itime]) and is_value(maxs[itime]):
            return mins[itime], maxs[itime]

        # save the defaults if they're not None
        mini2, maxi2 = self.get_default_min_max(itime, res_name)
        if is_blank(mins[itime]):
            mins[itime] = mini2
        if is_blank(maxs[itime]):
            maxs[itime] = maxi2
        return mins[itime], maxs[itime]

    # --------------------------------------------------------------------------
    def get_default_arrow_scale(self, itime: int, res_name: str) -> float:
        if not hasattr(self, 'dim_max'):
            return 0.

        if self.linked_scale_factor:
            itime = 0
        itime, case_flag = self.get_case_flag(itime, res_name)
        scales = self.default_arrow_scales[case_flag]
        if is_value(scales[itime]):
            return scales[itime]

        fringe_data = self._get_fringe_data_sparse(itime, res_name)

        # Fringe is typically a 'Magnitude' and therefore positive.
        # if 'Value' is used, it can be negative.
        abs_maximax = np.abs(fringe_data).max()
        if abs_maximax > 0.0:
            scale = self.dim_max / abs_maximax * 0.10
        else:
            scale = 1.0
        scales[itime] = scale
        return scale

    def get_default_scale(self, itime: int, res_name: str) -> float:
        if not hasattr(self, 'dim_max'):
            return 0.
        if self.linked_scale_factor:
            itime = 0
        itime, case_flag = self.get_case_flag(itime, res_name)
        scales = self.default_scales[case_flag]
        if is_value(scales[itime]):
            return scales[itime]

        scale = self._calculate_scale(itime, res_name)
        scales[itime] = scale
        return scale

    def _calculate_scale(self, i: int, resname: str) -> float:
        raise NotImplementedError(self.class_name)

    def get_scale(self, itime: int, res_name: str) -> float:
        if not hasattr(self, 'dim_max'):
            return 0.
        itime, case_flag = self.get_case_flag(itime, res_name)
        if self.linked_scale_factor:
            itime = 0
        scales = self.scales[case_flag]
        scale = scales[itime]
        if is_blank(scale):
            scale = self.get_default_scale(itime, res_name)
            scales[itime] = scale
        return scale
    def get_arrow_scale(self, itime: int, res_name: str) -> float:
        if not hasattr(self, 'dim_max'):
            return 0.
        itime, case_flag = self.get_case_flag(itime, res_name)
        if self.linked_scale_factor:
            itime = 0
        scales = self.arrow_scales[case_flag]
        scale = scales[itime]
        if is_blank(scale):
            scale = self.get_default_arrow_scale(itime, res_name)
            scales[itime] = scale
        return scale

    def set_scale(self, itime: int, res_name: str, scale: float) -> None:
        if not hasattr(self, 'dim_max'):
            return 0.
        itime, case_flag = self.get_case_flag(itime, res_name)
        if self.linked_scale_factor:
            itime = 0
        scales = self.scales[case_flag]
        scales[itime] = scale
    def set_arrow_scale(self, itime: int, res_name: str, scale: float) -> None:
        if not hasattr(self, 'dim_max'):
            return 0.
        itime, case_flag = self.get_case_flag(itime, res_name)
        if self.linked_scale_factor:
            itime = 0
        scales = self.arrow_scales[case_flag]
        scales[itime] = scale

    def get_default_phase(self, itime: int, res_name: str) -> float:
        return 0.0
    def get_phase(self, itime: int, res_name: str) -> float:
        if self.is_real:
            return 0.0
        itime, case_flag = self.get_case_flag(itime, res_name)
        phases = self.phases[case_flag]
        return phases[itime]
    def set_phase(self, itime: int, res_name: str, phase: float) -> None:
        if self.is_real:
            return
        itime, case_flag = self.get_case_flag(itime, res_name)
        phases = self.phases[case_flag]
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


class DispForceVectorResults(VectorResultsCommon):
    def __init__(self,
                 subcase_id: int,
                 node_id: np.ndarray,
                 case: Union[RealTableArray, ComplexTableArray],
                 t123_offset: int,
                 methods_txyz_rxyz: list[str],
                 index_to_base_title_annotation: dict[int, dict[str, str]],
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
        assert isinstance(index_to_base_title_annotation, dict), index_to_base_title_annotation
        self.index_to_base_title_annotation = index_to_base_title_annotation
        self.component_indices: tuple[int, ...] = (0, )
        self.dim = case.data.ndim
        assert self.dim == 3, case.data.shape
        assert case.data.shape[2] in {3, 6}, case.data.shape

        #expected 0 / 3 indices
        #if len(index_to_base_title_annotation[i0]) != 2:
            #x = 1
        assert isinstance(index_to_base_title_annotation[t123_offset]['title'], str)
        assert isinstance(index_to_base_title_annotation[t123_offset]['corner'], str)
        self.t123_offset = t123_offset
        self.methods_txyz_rxyz = methods_txyz_rxyz

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
        self.inode = np.array([], dtype='int32')

    def set_sidebar_args(self,
                         itime: int, res_name: str,
                         min_max_method: str='',
                         transform: str='',
                         # Magnitude (0, )
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

    #--------------------------------------------------------------

    def get_methods(self, itime: int, res_name: str) -> list[str]:
        if self.is_real:
            i0 = self.t123_offset
            out = ['Magnitude'] + self.methods_txyz_rxyz[i0:i0+3]
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
        if 'chL' in getpass.getuser():
            out[0] = 'Reluctant'
        return out

    def get_location(self, unused_i: int, unused_res_name: str) -> str:
        """the result location (node/centroid)"""
        assert self.location, self.location
        return self.location

    def get_case_flag(self, itime: int, res_name) -> tuple[int, Any]:
        case_flag = (self.component_indices, self.min_max_method)
        return itime, case_flag

    #--------------------------------------------------------------
    def has_methods_table(self, i: int, res_name: str) -> bool:
        return True
    def has_coord_transform(self, i: int, res_name: str) -> tuple[bool, list[str]]:
        return True, ['Global']
    def has_derivation_transform(self, i: int, res_name: str,
                                 ) -> tuple[bool, dict[str, Any]]:
        """min/max/avg"""
        out = {
            #'label': 'Derivation Method:',
            'derivation': ['Magnitude', 'Value'],
            'tooltip': 'Magnitude is automatically selected if multiple components are selected',
        }
        return True, out
    def has_nodal_combine_transform(self, i: int, res_name: str) -> tuple[bool, list[str]]:
        """elemental -> nodal"""
        return False, []

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
        #title0 = self._title0
        title0 = self.index_to_base_title_annotation[self.t123_offset]['corner']
        method = title0 + ''.join(self.index_map[idx] for idx in self.component_indices)
        #annotation_label = f'{self.title} {method} ({self.min_max_method}): {self.headers[itime]}'
        annotation_label = f'{self.title} {method}: {self.headers[itime]}'
        #return self.uname
        return annotation_label

    def get_default_legend_title(self, itime: int, res_name: str) -> str:
        self.component_indices
        title0 = self.index_to_base_title_annotation[self.t123_offset]['title']
        method = title0 + ''.join(self.index_map[idx] for idx in self.component_indices)
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
        title0 = self.index_to_base_title_annotation[self.t123_offset]['title']
        method = title0 + ''.join(self.index_map[idx] for idx in self.component_indices)
        title = f'{self.title} {method}'
        return title

    def _get_data(self, itime: int) -> tuple[np.ndarray, list[int]]:
        """
        Gets the real data

        Returns
        -------
        data : (3, ) float np.ndarray
            real/complex ndarray
        idxs : list[int]
            the indicies that were used

        """
        # handles translation vs. rotation
        datai = self.case.data
        i0 = self.t123_offset
        if i0 == 1:
            assert datai.shape[2] == 6, datai.shape
        else:
            assert datai.shape[2] in {3, 6}, datai.shape

        data = datai[itime, :, i0:i0+3].copy()
        assert data.shape[1] == 3, data.shape
        assert len(self.component_indices) > 0, self.component_indices

        idxs = []
        if self.is_real:
            for idx in [0, 1, 2]:
                if idx not in self.component_indices:
                    data[:, idx] = 0.0
                    idxs.append(idx)
        else:
            real = data.real
            imag = data.imag
            # real
            for idx in [7, 8, 9]:
                if idx not in self.component_indices:
                    real[:, idx] = 0.0
                    idxs.append(idx)
            # complex
            for idx in [10, 11, 12]:
                if idx not in self.component_indices:
                    imag[:, idx] = 0.0
                    idxs.append(idx)
        return data, idxs

    def _get_vector_data_sparse(self, itime: int,
                                res_name: str) -> tuple[np.ndarray, list[int]]:
        """
        gets the 'typical' result as a vector
         - DisplacementResults: vector; (n, 3) array

        Parameters
        ----------
        return_dense: bool
            Returns the data array in a way that the gui can use.
            Handles the null result case (e.g; SPC forces only
            at the SPC location).

        """
        if self.min_max_method == 'Value':
            assert len(self.component_indices) == 1, self.component_indices

        # multiple results
        # .0006 -> 0.0
        # .057 -> 0.0123
        # min
        dxyz, idxs = self._get_data(itime)
        #dxyz = data[itime, :, :]
        assert len(dxyz.shape) == 2, dxyz.shape
        return dxyz, idxs

    def get_vector_data_dense(self, itime: int, res_name: str) -> tuple[np.ndarray, int, Any]:
        """calculates a dense vector"""
        itime, case_flag = self.get_case_flag(itime, res_name)
        dxyz_sparse, unused_idx = self._get_vector_data_sparse(itime, res_name)

        if self.is_dense:
            return dxyz_sparse, itime, case_flag

        # apply dense
        nnodes = len(self.node_id)
        dxyz = np.full((nnodes, 3), np.nan, dtype=dxyz_sparse.dtype)
        dxyz[self.inode, :] = dxyz_sparse
        return dxyz, itime, case_flag

    def _get_fringe_data_sparse(self, itime: int, res_name: str) -> np.ndarray:
        unused_ntimes, nnodes = self.case.data.shape[:2]
        vector_result_sparse, *unused_junk = self._get_vector_data_sparse(itime, res_name)
        if self.min_max_method == 'Value':
            assert len(self.component_indices) == 1, self.component_indices
            fringe_result_sparse = vector_result_sparse[:, self.component_indices[0]]
        else:
            fringe_result_sparse = safe_norm(vector_result_sparse, axis=1)
        assert len(fringe_result_sparse) == nnodes
        return fringe_result_sparse

    def get_fringe_result_dense(self, itime: int, res_name: str) -> np.ndarray:
        """get_fringe_value"""
        fringe_result_sparse = self._get_fringe_data_sparse(itime, res_name)

        if self.is_dense:
            return fringe_result_sparse
        #case.data.shape = (11, 43, 6)
        #nnodes = len(self.node_id) =  48
        #nnodesi = len(self.inode) = len(self.case.node_gridtype) = 43
        fringe_result_dense = np.full(len(self.node_id), np.nan, dtype=fringe_result_sparse.dtype)
        fringe_result_dense[self.inode] = fringe_result_sparse
        return fringe_result_dense

    def get_fringe_vector_result(self, itime: int,
                                 res_name: str) -> tuple[np.ndarray, np.ndarray]:
        """
        gets the 'typical' result which is a vector
         - GuiResult:           fringe; (n,)   array
         - Displacement:        vector; (n,3)  array

        Parameters
        ----------
        return_dense: bool
            Rreturns the data array in a way that the gui can use.
            Handles the null result case (e.g; SPC forces only
            at the SPC location).
        return_vector: bool
            Don't compress the result into a normalized vector (for min/max calculation)

        """
        vector_result, itime, case_flag = self.get_vector_data_dense(itime, res_name)

        if self.min_max_method == 'Value':
            assert len(self.component_indices) == 1, self.component_indices
            fringe_result = vector_result[:, self.component_indices[0]]
        else:
            fringe_result = safe_norm(vector_result, axis=1)

        # set the min/max if they're not set
        default_mins = self.default_mins[case_flag]
        default_maxs = self.default_maxs[case_flag]

        if is_blank(default_mins[itime]) and is_blank(default_maxs[itime]):
            mini = np.nanmin(fringe_result)
            maxi = np.nanmax(fringe_result)

            self.mins[case_flag][itime] = mini
            self.maxs[case_flag][itime] = maxi
            default_mins[itime] = mini
            default_maxs[itime] = maxi

        return fringe_result, vector_result

def _get_real_component_indices(methods_keys: Optional[list[int]]) -> tuple[int, ...]:
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

def _get_complex_component_indices(methods_keys: list[str]) -> tuple[int, ...]:
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


def _to_dense_vector(dxyz: np.ndarray,
                     inode: np.ndarray,
                     nnodes: int,
                     return_dense: bool,
                     is_dense: bool) -> np.ndarray:
    """map vectored data to dense"""
    return_sparse = not return_dense
    if (return_sparse or is_dense):
        dxyz2 = dxyz
    else:
        dxyz2 = np.full((nnodes, 3), np.nan, dtype=dxyz.dtype)
        dxyz2[inode, :] = dxyz
    return dxyz2

def is_value(value: Optional[float]) -> bool:
    return not is_blank(value)

def is_blank(value: Optional[float]) -> bool:
    if value is None or not np.isfinite(value):
        return True
    return False
