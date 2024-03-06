import getpass
from abc import abstractmethod
from collections import defaultdict

import numpy as np
from typing import DefaultDict, Union, Optional, Any # , TYPE_CHECKING

from pyNastran.utils.numpy_utils import integer_types
from pyNastran.femutils.utils import safe_norm
#if TYPE_CHECKING:  # pragma: no cover
from pyNastran.op2.result_objects.table_object import (
    RealTableArray, ComplexTableArray)
from pyNastran.gui.gui_objects.gui_result import GuiResultCommon
from pyNastran.gui.utils.utils import is_blank, is_value


col_axis = 1
#IDX_REAL = (0, 1, 2)
#IDX_IMAG = (3, 4, 5)

IDX_REAL = (0, 1, 2)
IDX_IMAG = (3, 4, 5)
#IDX_REAL = (7, 8, 9)
#IDX_IMAG = (10, 11, 12)
COMPLEX_DEFAULT_INDICES = tuple(list(IDX_REAL) + list(IDX_IMAG))

#if TYPE_CHECKING:  # pragma: no cover
    #from typing import TypedDict
    #class SidebarArgs(TypedDict):
        #label: str
        #tooltip: int
        #derivation: list[str]


class VectorResultsCommon(GuiResultCommon):
    def __init__(self,
                 subcase_id: int,
                 title: str,
                 case,
                 #dxyz: Union[RealTableArray, ComplexTableArray],
                 ntitles: int,
                 data_format: str='%g',
                 is_variable_data_format: bool=False,
                 nlabels=None, labelsize=None, ncolors=None,
                 colormap: str='',
                 #set_max_min: bool=False,
                 uname: str='VectorResults2'):
        GuiResultCommon.__init__(self)
        self.is_dense = False
        self.subcase_id = subcase_id
        self._title = title

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
        def findex() -> np.ndarray:
            # -1 -> nan
            # -2 -> default
            return np.ones(ntimes, dtype='int32') * -2

        # floats
        self.default_arrow_scales: DefaultDict[Any, list[float]] = defaultdict(fscales)
        self.arrow_scales:         DefaultDict[Any, list[float]] = defaultdict(fscales)

        self.default_scales: DefaultDict[Any, list[float]] = defaultdict(fscales)
        self.scales:         DefaultDict[Any, list[float]] = defaultdict(fscales)

        self.default_mins: DefaultDict[Any, list[float]] = defaultdict(ftimes)
        self.default_maxs: DefaultDict[Any, list[float]] = defaultdict(ftimes)
        self.mins:         DefaultDict[Any, list[float]] = defaultdict(ftimes)
        self.maxs:         DefaultDict[Any, list[float]] = defaultdict(ftimes)
        self.phases: DefaultDict[Any, np.ndarray] = defaultdict(fphases)

        # ints
        self.imins:        DefaultDict[Any, np.ndarray] = defaultdict(findex)
        self.imaxs:        DefaultDict[Any, np.ndarray] = defaultdict(findex)
        self.title = [''] * ntitles

        self.data_format = data_format
        self.data_formats = [self.data_format]
        self.is_variable_data_format = is_variable_data_format
        self.headers = ['VectorResultsCommon'] * ntimes

        self.nlabels = None
        self.labelsize = None
        self.ncolors = None
        self.colormap = colormap

        self.uname = uname
        self.location = ''
        self.min_max_method = ''
        self.nodal_combine = ''
        self.is_method_array = True

    # --------------------------------------------------------------------------
    # abstractmethods - results access
    @abstractmethod
    def _get_fringe_data_sparse(self, itime: int, res_name) -> np.ndarray:
        raise NotImplementedError(f'{self.class_name}._get_fringe_data_sparse')

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

    @abstractmethod
    def get_legend_tuple(self, i: int, resname: str) -> Any:
        raise NotImplementedError(f'{self.class_name}.get_legend_tuple')

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

        fringe_result = self._get_fringe_data_sparse(itime, res_name)
        mini, maxi = self._set_default_from_fringe(
            itime, case_flag, fringe_result, is_sparse=True)

        return mini, maxi

    def get_imin_imax(self, itime: int, res_name) -> tuple[int, int]:
        itime, case_flag = self.get_case_flag(itime, res_name)
        imins = self.imins[case_flag]
        imaxs = self.imaxs[case_flag]
        imin = imins[itime]
        imax = imaxs[itime]
        if imin == -2:
            self.get_default_min_max(itime, res_name)
        # -1 -> nan
        # -2 -> default
        assert imins[itime] != -2, imins
        assert isinstance(imin, integer_types), imin
        assert isinstance(imax, integer_types), imax
        return imins[itime], imaxs[itime]

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

    def get_location_arrays(self) -> tuple[np.ndarray, np.ndarray]:
        """used for _set_default_from_fringe"""
        raise NotImplementedError(f'{self.class_name}.get_location_arrays')

    def _set_default_from_fringe(self, itime: int, case_flag,
                                 fringe_result: np.ndarray,
                                 is_sparse: bool) -> tuple[float, float]:
        """
        set the min/max if they're not set
        also sets node_id imin/imax flags
        """
        default_mins = self.default_mins[case_flag]
        default_maxs = self.default_maxs[case_flag]

        if is_blank(default_mins[itime]) and is_blank(default_maxs[itime]):
            # save the defaults for the next time

            imins = self.imins[case_flag]
            imaxs = self.imaxs[case_flag]

            mins = self.mins[case_flag]
            maxs = self.maxs[case_flag]
            try:
                # gold standard
                mini1 = np.nanmin(fringe_result)
                maxi1 = np.nanmax(fringe_result)

                # solid_shell_bar -> node=13
                imin = np.nanargmin(fringe_result)
                imax = np.nanargmax(fringe_result)
                assert isinstance(imin, integer_types), imin
                assert isinstance(imax, integer_types), imax

                mini2 = fringe_result[imin]
                maxi2 = fringe_result[imax]
                assert np.allclose(mini1, mini2)
                assert np.allclose(maxi1, maxi2)

                default_mins[itime] = mini2
                default_maxs[itime] = maxi2
                mins[itime] = mini2
                maxs[itime] = maxi2

                if is_sparse:
                    # we found the id as sparse
                    # map it to dense
                    node_id, nids = self.get_location_arrays()
                    iimin, iimax = np.searchsorted(node_id, nids[[imin, imax]])
                    imins[itime] = iimin
                    imaxs[itime] = iimax
                else:
                    imins[itime] = imin
                    imaxs[itime] = imax
            except ValueError:
                # All NaN
                default_mins[itime] = np.nan
                default_maxs[itime] = np.nan
                mins[itime] = np.nan
                maxs[itime] = np.nan
                imins[itime] = 0
                imaxs[itime] = 0
        return default_mins[itime], default_maxs[itime]

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
            return
        itime, case_flag = self.get_case_flag(itime, res_name)
        if self.linked_scale_factor:
            itime = 0
        scales = self.scales[case_flag]
        scales[itime] = scale
    def set_arrow_scale(self, itime: int, res_name: str, scale: float) -> None:
        if not hasattr(self, 'dim_max'):
            return
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
    def set_legend_title(self, i: int, name: str, title: str) -> None:
        """if legend is default or blank, use the default legend"""
        legend_tuple = self.get_legend_tuple(i, name)
        legend_default = self.get_default_legend_title(i, name)
        if title == legend_default:
            title = ''
        self.title[legend_tuple] = title

    def get_legend_title(self, itime: int, res_name: str):
        """
        Uses self.title if it's set.
        Otherwise, use the default title.
        """
        legend_tuple = self.get_legend_tuple(itime, res_name)
        titlei = self.title[legend_tuple]
        #method2 = self._update_method(method)
        #return f'{self.title} {method2}'
        #self.component_indices (0, 1, 2)
        title = titlei if titlei else self.get_default_legend_title(itime, res_name)
        return title

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
        nlabels = -1 if nlabels is None else nlabels
        labelsize = -1 if labelsize is None else labelsize
        ncolors = -1 if ncolors is None else ncolors
        self.nlabels = nlabels
        self.labelsize = labelsize
        self.ncolors = ncolors
        self.colormap = colormap


class DispForceVectorResults(VectorResultsCommon):
    def __init__(self,
                 subcase_id: int,
                 title: str,
                 node_id: np.ndarray,
                 case: Union[RealTableArray, ComplexTableArray],
                 t123_offset: int,
                 methods_txyz_rxyz: list[str],
                 index_to_base_title_annotation: dict[int, dict[str, str]],
                 dim_max: float,
                 data_format: str='%g',
                 is_variable_data_format: bool=False,
                 nlabels=None, labelsize=None, ncolors=None,
                 colormap: str='',
                 set_max_min: bool=False,
                 uname: str='VectorResults2'):
        ntitles = 1
        VectorResultsCommon.__init__(
            self, subcase_id, title, case, ntitles,
            data_format=data_format,
            is_variable_data_format=is_variable_data_format,
            nlabels=nlabels, labelsize=labelsize, ncolors=ncolors,
            colormap=colormap, uname=uname)
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
            #self.component_indices = (7, 8, 9, 10, 11, 12)
            self.component_indices = COMPLEX_DEFAULT_INDICES # (1, 2, 3, 4, 5, 6)
            self.index_map = {
                #0: 'Magnitude',
                #1: 'X Mag', 2: 'Y Mag', 3: 'Z Mag',
                #4: 'X Phase', 5: 'Y Phase', 6: 'Z Phase',
                #7: 'X Real', 8: 'Y Real', 9: 'Z Real',
                #10: 'X Imaginary', 11: 'Y Imaginary', 12: 'Z Imaginary',
                0: 'X Real', 1: 'Y Real', 2: 'Z Real',
                3: 'X Imaginary', 4: 'Y Imaginary', 5: 'Z Imaginary',
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
                'Magnitude',
                #'X Magnitude', 'Y Magnitude', 'Z Magnitude',
                #'X Phase', 'Y Phase', 'Z Phase',
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
        header = 'mode = 10; freq = 2.85697 Hz'
        title = 'Displacement'
        method = 'T_XYZ'
        min_max_method = 'Value'
        returns 'Displacement (Static) T_XYZ'

        Magnitude vs. Value
        """
        #title0 = self.index_to_base_title_annotation[self.t123_offset]['corner']
        #method = title0 + ''.join(self.index_map[idx] for idx in self.component_indices)
        header = self.headers[itime]
        #annotation_label = f'{self.title} {method} ({self.min_max_method}): {self.headers[itime]}'


        # v1
        #title = f'{self._title} {method}'
        title = self.get_legend_title(itime, res_name)
        annotation_label = f'{title} ({self.min_max_method}, {header})'
        return annotation_label

    def get_legend_tuple(self, itime: int, res_name: str) -> int:
        return 0
    def get_default_legend_title(self, itime: int, res_name: str) -> str:
        #self.component_indices (0, 1, 2)
        # T_
        title0 = self.index_to_base_title_annotation[self.t123_offset]['title']
        if self.is_real:
            # T_XYZ
            method = title0 + ''.join(self.index_map[idx] for idx in self.component_indices)
        else:
            method = ''
            idxs = IDX_REAL
            if any(idi in self.component_indices for idi in idxs):
                # real
                method ='r'
                for idx in self.component_indices:
                    if idx in idxs:
                        mapped_index = self.index_map[idx].split()[0]
                        method += mapped_index
            idxs = IDX_IMAG
            if any(idi in self.component_indices for idi in idxs):
                # imag
                method += '_i'
                for idx in self.component_indices:
                    if idx in idxs:
                        mapped_index = self.index_map[idx].split()[0]
                        method += mapped_index
            method = method.replace('__', '_')
        title_out = f'{self._title} {method}'
        return title_out

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
        datai = self.case.data
        i0 = self.t123_offset
        if i0 == 1:
            # rotations
            assert datai.shape[2] == 6, datai.shape
        else:
            # translations; i0=0
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
            for idx in IDX_REAL:
                if idx not in self.component_indices:
                    idxi = idx - IDX_REAL[0]
                    real[:, idxi] = 0.0
                    #real[:, idx] = 0.0
                    idxs.append(idx)
            # complex
            for idx in IDX_IMAG:
                if idx not in self.component_indices:
                    idxi = idx - IDX_IMAG[0]
                    imag[:, idxi] = 0.0
                    #imag[:, idx] = 0.0
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
            fringe_result_sparse = get_fringe_value_array(
                self.component_indices, vector_result_sparse,
                self.is_real, self.is_complex)
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

    def get_location_arrays(self) -> tuple[np.ndarray, np.ndarray]:
        return self.node_id, self.case.node_gridtype[:, 0]

    def get_fringe_result(self, itime: int, res_name: str) -> np.ndarray:
        fringe, vector = self.get_fringe_vector_result(itime, res_name)
        return fringe

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
            fringe_result = get_fringe_value_array(
                self.component_indices, vector_result,
                self.is_real, self.is_complex)
        else:
            fringe_result = safe_norm(vector_result, axis=1)

        self._set_default_from_fringe(itime, case_flag, fringe_result,
                                      is_sparse=False)
        return fringe_result, vector_result


def get_fringe_value_array(component_indices: tuple[int, ...],
                           vector_result: np.ndarray,
                           is_real: bool, is_complex: bool, ) -> np.ndarray:
    """doesn't support magnitude/phase"""
    assert len(component_indices) == 1, component_indices
    component = component_indices[0]
    if is_real or is_complex and component <= 2:
        # real
        fringe_result = vector_result[:, component].real
    else:
        # imaginary
        fringe_result = vector_result[:, component-3].imag
    return fringe_result

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

def _get_complex_component_indices(methods_keys: list[int]) -> tuple[int, ...]:
    """
    if Magnitude is selected, only use magnitude
    methods = [
        'Resultant',
        #'X Magnitude', 'Y Magnitude', 'Z Magnitude',
        #'X Phase', 'Y Phase', 'Z Phase',
        'X Real', 'Y Real', 'Z Real',
        'X Imaginary', 'Y Imaginary', 'Z Imaginary',
    ]
    """
    # the data is real/imag, so we use that...
    #default_indices = [1, 2, 3, 4, 5, 6]    # mag/phase
    default_indices = COMPLEX_DEFAULT_INDICES  # real/imag
    if methods_keys is None or len(methods_keys) == 0:
        # default
        #indices = default_indices
        return default_indices
    elif 0 in methods_keys:
        # include all components b/c Magnitude is selected
        #indices = default_indices
        return default_indices
    else:
        # no 0 (magnitude) in methods_keys
        # update the indices to correspond to the array
        #
        # first split names by type
        for key in methods_keys:
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


def filter_ids(all_element_id: np.ndarray,
               eids: np.ndarray) -> tuple[np.ndarray, np.ndarray,
                                          int, bool]:
    """filters a set of elements"""
    neids = len(eids)
    is_filter = False
    intersect_eids = np.intersect1d(all_element_id, eids)
    if neids == len(intersect_eids):
        return intersect_eids, np.array([], dtype='int32'), neids, is_filter

    is_filter = True
    ieid_filter = np.searchsorted(eids, intersect_eids)
    nelement_filtered = len(intersect_eids)
    return intersect_eids, ieid_filter, nelement_filtered, is_filter
