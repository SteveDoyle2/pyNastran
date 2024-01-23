from __future__ import annotations
from collections import defaultdict
import numpy as np
from typing import Optional, TYPE_CHECKING

from pyNastran.gui.gui_objects.gui_result import GuiResultCommon
from pyNastran.femutils.utils import pivot_table, abs_nan_min_max # abs_min_max
from pyNastran.bdf.utils import write_patran_syntax_dict

from .displacement_results import VectorResultsCommon
if TYPE_CHECKING:
    from pyNastran.bdf.bdf import BDF
    from pyNastran.op2.tables.oes_stressStrain.real.oes_composite_plates import RealCompositePlateArray

pcomp_stress = ['o11', 'o22', 't12', 't1z', 't2z', 'oangle', 'Max Principal', 'minor', 'ovm', 'omax_shear']
pcomp_strain = ['e11', 'e22', 'e12', 'e1z', 'e2z', 'eangle', 'Max Principal', 'minor', 'evm', 'emax_shear']
col_axis = 1

class CompositeResults2(VectorResultsCommon):
    def __init__(self,
                 subcase_id: int,
                 model: BDF,
                 element_id: np.ndarray,
                 case: RealCompositePlateArray,
                 result: str,
                 #dim_max: float,
                 data_format: str='%g',
                 nlabels=None, labelsize=None, ncolors=None,
                 colormap: str='',
                 set_max_min: bool=False,
                 uname: str='CompositeResults2'):
        GuiResultCommon.__init__(self)
        self.layer_indices = (-1, )  # All
        i = -1
        name = None

        # slice off the methods (from the boolean) and then pull the 0th one
        self.min_max_method = self.has_derivation_transform(i, name)[1]['derivation'][0]
        self.transform = self.has_coord_transform(i, name)[1][0]
        self.nodal_combine = self.has_nodal_combine_transform(i, name)[1][0]

        self.is_dense = False
        self.dim = case.data.ndim
        assert self.dim == 3, case.data.shape

        self.subcase_id = subcase_id
        self.is_stress = case.is_stress
        if self.is_stress:
            self.iresult_map = {
                0: 'o11',
                1: 'o22',
                3: 't12',
                4: 't1z',
                5: 't2z',
                6: 'angle',
                7: 'major',
                8: 'minor',
                9: 'von_mises',
                10: 'max_shear',
                11: 'fiber_distance',
                12: 'fiber_curvature',
            }
        else:
            self.iresult_map = {
                0: 'e11',
                1: 'e22',
                3: 'e12',
                4: 'e1z',
                5: 'e2z',
                6: 'angle',
                7: 'major',
                8: 'minor',
                9: 'von_mises',
                10: 'max_shear',
                11: 'fiber_distance',
                12: 'fiber_curvature',
            }

        #if dim_max == 0.0:
            #dim_max = 1.0
        #self.dim_max = dim_max
        self.linked_scale_factor = False

        self.data_format = data_format

        #  global ids
        self.model = model
        self.element_id = element_id

        # local case object
        self.case = case
        self.result = result
        layers = case.element_layer[:, 1]
        ulayers = np.unique(layers)
        layers = {layer: f'Layer {layer}' for layer in ulayers}
        self.layer_map = {0: 'All Layers'}
        self.layer_map.update(layers)

        self.data_type = self.case.data.dtype.str # '<c8', '<f4'
        self.is_real = True if self.data_type in ['<f4', '<f8'] else False
        #self.is_real = dxyz.data.dtype.name in {'float32', 'float64'}
        self.is_complex = not self.is_real

        ntimes = case.data.shape[0]
        #nscale = ntimes
        #if self.linked_scale_factor:
            #nscale = 1

        #def fscales():
            #return [None] * nscale
        def ftimes():
            return [None] * ntimes
        def fphases():
            return np.zeros(ntimes, dtype='float64')

        #self.default_scales = defaultdict(fscales)
        #self.scales = defaultdict(fscales)
        self.default_mins = defaultdict(ftimes)
        self.default_maxs = defaultdict(ftimes)
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

    def _get_default_tuple_indices(self):
        out = tuple(np.array(self._get_default_layer_indicies()) - 1)
        return out

    def _get_default_layer_indicies(self):
        default_indices = list(self.layer_map.keys())
        default_indices.remove(0)
        return default_indices

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
        combine_methods = self.has_nodal_combine_transform(itime, res_name)[1]

        transform = transform if transform else transforms[0]
        min_max_method = min_max_method if min_max_method else min_max_methods[0]
        nodal_combine = nodal_combine if nodal_combine else combine_methods[0]

        assert transform in transforms, transform
        assert min_max_method in min_max_methods, min_max_method
        assert nodal_combine in combine_methods, nodal_combine

        #sidebar_kwargs = {
            #'min_max_method': min_max_method,
            #'transform': coord,
            #'nodal_combine': nodal_combine,
            #'methods_keys': keys_b,
        #}
        # if Magnitude is selected, only use magnitude
        # methods = ['T_XYZ', 'TX', 'TY', 'TZ']
        default_indices = self._get_default_layer_indicies()
        if methods_keys is None or len(methods_keys) == 0:
            # default; All
            indices = default_indices
        elif 0 in methods_keys:
            # include all components b/c All is selected
            indices = default_indices
        else:
            # no 0 (magnitude) in methods_keys
            # update the indices to correspond to the array
            #methods_keys.sort()
            indices = methods_keys
        self.layer_indices = tuple(np.array(indices, dtype='int32') - 1)
        #self.layer_indices = (1, )

        # doesn't matter cause it's already nodal
        assert min_max_method in min_max_methods, min_max_method
        assert nodal_combine in combine_methods, nodal_combine
        self.min_max_method = min_max_method
        self.nodal_combine = nodal_combine
        self.transform = transform

    def has_methods_table(self, i: int, res_name: str) -> bool:
        return True
    def has_coord_transform(self, i: int, res_name: str) -> tuple[bool, list[str]]:
        return True, ['Material']
    def has_derivation_transform(self, i: int, case_tuple: str) -> tuple[bool, list[str]]:
        """min/max/avg"""
        #(itime, iresult, header) = case_tuple
        out = {
            'tooltip': 'Method to reduce multiple layers into a single elemental value',
            'derivation': ['Absolute Max', 'Min', 'Max', 'Mean', 'Std. Dev.', 'Difference',
                           #'Derive/Average'
                           ],

        }
        return True, out
    def has_nodal_combine_transform(self, i: int, res_name: str) -> tuple[bool, list[str]]:
        """elemental -> nodal"""
        return True, ['Centroid']
        #return True, ['Absolute Max', 'Min', 'Max']

    def get_annotation(self, itime: int, case_tuple: str) -> str:
        """
        A header is the thingy that goes in the lower left corner
        title = 'Compostite Plate Stress'
        method = 'Absolute Max'
        header = 'Static'
        returns 'Compostite Plate Stress All Layers (Absolute Max; Static): sigma11'
        """
        # overwrite itime based on linked_scale factor
        (itime, iresult, header) = case_tuple
        itime, unused_case_flag = self.get_case_flag(case_tuple)

        default_indices = self._get_default_tuple_indices() # 0-based
        if self.layer_indices == default_indices:
            layer_str = 'All Layers'
        else:
            self.layer_indices
            dict_sets = {'': [(idx+1) for idx in self.layer_indices],}
            layer_str = 'Layers ' + write_patran_syntax_dict(dict_sets)

        results = list(self.result.keys())
        result = results[iresult]

        #'Compostite Plate Stress (Absolute Max; Static): sigma11'
        annotation_label = f'{self.title}; {layer_str} ({self.min_max_method}, {header}): {result}'
        #return self.uname
        return annotation_label

    def get_default_min_max(self, itime: int,
                            case_tuple: str) -> tuple[float, float]:
        #(itime, iresult, unused_header) = case_tuple
        itime, case_flag = self.get_case_flag(case_tuple)
        mins = self.default_mins[case_flag]
        maxs = self.default_maxs[case_flag]
        if mins[itime] is not None and maxs[itime] is not None:
            return mins[itime], maxs[itime]

        datai = self._get_real_data(case_tuple)
        mins[itime] = np.nanmin(datai)
        maxs[itime] = np.nanmax(datai)
        return mins[itime], maxs[itime]

    def get_min_max(self, itime, case_tuple) -> tuple[float, float]:
        #(itime, iresult, header) = case_tuple
        itime, case_flag = self.get_case_flag(case_tuple)
        mins = self.mins[case_flag]
        maxs = self.maxs[case_flag]
        if mins[itime] is not None and maxs[itime] is not None:
            return mins[itime], maxs[itime]

        # save the defaults if they're not None
        mini2, maxi2 = self.get_default_min_max(itime, case_tuple)
        if mini2 is not None:
            mins[itime] = mini2
        if maxi2 is not None:
            maxs[itime] = maxi2
        return mins[itime], maxs[itime]

    def set_min_max(self, itime, case_tuple, min_value, max_value) -> tuple[float, float]:
        #(itime, iresult, header) = case_tuple
        itime, case_flag = self.get_case_flag(case_tuple)

        mins = self.mins[case_flag]
        maxs = self.maxs[case_flag]
        mins[itime] = min_value
        maxs[itime] = max_value

    def get_case_flag(self, case_tuple: tuple[int, int, str]) -> tuple[int,
                                                                       tuple[int, int, tuple, str]]:
        """
        itime = 0
        iresult = 0 # o11
        layer_indices = (1, 2, 3, 4, 5)
        min_max_method = 'Absolute Max'
        """
        (itime, iresult, header) = case_tuple
        #if self.is_linked_scale_factor:
            #itime = 0

        return itime, (itime, iresult, self.layer_indices, self.min_max_method)

    def get_default_legend_title(self, itime: int, case_tuple: str) -> str:
        (itime, iresult, header) = case_tuple
        #method_ = 'Composite Stress Layers:' if self.is_stress else 'Composite Strain Layers:'
        #self.layer_indices
        results = list(self.result.values())
        #method = method_ + ', '.join(str(idx) for idx in (self.layer_indices+1))
        #method = method.strip()
        #title = f'{self.title} {method}'
        title = results[itime]
        return title
    def set_legend_title(self, itime: int, res_name: str,
                         title: str) -> None:
        self.title = title
    def get_legend_title(self, itime: int, case_tuple: str):
        """Composite Stress Layers: 1, 2, 3, 4"""
        (itime, iresult, header) = case_tuple
        #method_ = 'Composite Stress Layers:' if self.is_stress else 'Composite Strain Layers:'
        #self.layer_indices
        #self.result
        #method = method_ + ', '.join(str(idx) for idx in (self.layer_indices+1))
        #title = f'{self.title} {method}'
        results = list(self.result.values())
        return results[iresult]

    def _get_real_data(self, case_tuple: int) -> np.ndarray:
        (itime, iresult, header) = case_tuple

        #self.case.get_headers()
        #['o11', 'o22', 't12', 't1z', 't2z', 'angle', 'major', 'minor', 'max_shear']
        data = self.case.data[itime, :, iresult].copy()
        rows = self.case.element_layer[:, 0]
        cols = self.case.element_layer[:, 1]
        ulayer = np.unique(cols)

        if self.layer_indices == (-1, ):
            self.layer_indices = tuple(ulayer - 1)
        #if self.is_stress:
            #datai = self.dxyz.data[itime, :, :3].copy()
            #assert datai.shape[1] == 3, datai.shape
        #else:
            #datai = self.dxyz.data[itime, :, 3:].copy()
            #assert datai.shape[1] == 3, datai.shape
        #if self.result == 'o11':

            #case_res = self.case.data[0]

        #if imethod == 0:
            ## fiber distance
            #print(self.case.get_stats())
            #asdf
        #else:
        assert len(data.shape) == 1, data.shape
        data2, eids_new = pivot_table(data, rows, cols, shape=1)
        assert len(data2.shape) == 2, data2.shape

        # element_id is unique & sorted
        # eids_new   is unique & sorted
        #
        # bwb_saero.op2
        # data2.shape = (9236, 10)
        #if not hasattr(self, 'ielement'):
        ieids = np.searchsorted(self.element_id, eids_new)
        self.ielement = ieids
        neids = len(self.ielement)

        data3 = data2
        if data2.shape[1] != len(self.layer_indices):
            # slice out the correct layers
            data3 = data2[:, self.layer_indices]

        assert len(data3.shape) == 2, data3.shape
        if data3.shape[1] == 1 and 0:  # pragma: no cover
            # single ply
            data4 = data3[:, 0]
        else:
            # multiple plies
            # ['Absolute Max', 'Min', 'Max', 'Derive/Average']
            axis = 1
            if self.min_max_method == 'Absolute Max':
                data4 = abs_nan_min_max(data3, axis=axis)
            elif self.min_max_method == 'Min':
                data4 = np.nanmin(data3, axis=axis)
            elif self.min_max_method == 'Max':
                data4 = np.nanmax(data3, axis=axis)
            elif self.min_max_method == 'Mean':  #   (Derive/Average)???
                data4 = np.nanmean(data3, axis=axis)
            elif self.min_max_method == 'Std. Dev.':
                data4 = np.nanstd(data3, axis=axis)
            elif self.min_max_method == 'Difference':
                data4 = np.nanmax(data3, axis=axis) - np.nanmin(data2, axis=axis)
            #elif self.min_max_method == 'Max Over Time':
                #data4 = np.nanmax(data3, axis=axis) - np.nanmin(data2, axis=axis)
            #elif self.min_max_method == 'Derive/Average':
                #data3 = np.nanmax(data2, axis=1)
            else:  # pragma: no cover
                raise NotImplementedError(self.min_max_method)

        # TODO: hack to try and debug things...
        assert data4.shape == (neids, )
        #data4 = eids_new.astype('float32')
        return data4

    #def _get_complex_data(self, itime: int) -> np.ndarray:
        #return self._get_real_data(itime)
        #if self.is_translation:
            #datai = self.dxyz.data[itime, :, :3]
            #assert datai.shape[1] == 3, datai.shape
        #else:
            #datai = self.dxyz.data[itime, :, 3:]
            #assert datai.shape[1] == 3, datai.shape
        #return datai

    def get_result(self, itime: int, case_tuple: str,
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
        method = self._update_method(itime, case_tuple, method)
        assert self.is_real
        # multiple results
        # .0006 -> 0.0
        # .057 -> 0.0123
        # min
        data = self._get_real_data(case_tuple)
        #dxyz = data[itime, :, :]
        #else:
        #data = self._get_complex_data(case_tuple)
        assert len(data.shape) == 1, data.shape

        return_sparse = not return_dense
        if return_sparse or self.is_dense:
            return data

        nelements = len(self.element_id)
        result_out = np.full(nelements, np.nan, dtype=data.dtype)
        result_out[self.ielement] = data
        return result_out

    def _update_method(self, i: int, case_tuple: str,
                       method: str) -> str:
        if method == '':
            method = 'All Layers'
        self.active_method = method
        assert method in self.get_methods(i, case_tuple)
        return method
    #def _update_method(self, method: str) -> str:
        #assert method != '', method
        #self.active_method = method

    def get_default_scale(self, itime: int, res_name: str) -> float:
        return None
    def get_scale(self, itime: int, res_name: str) -> float:
        return 0.0
    def set_scale(self, itime: int, res_name: str) -> None:
        return

    def get_default_phase(self, itime: int, res_name: str) -> float:
        return 0.0
    def get_phase(self, itime: int, res_name: str) -> float:
        return 0.0
    def set_phase(self, itime: int, res_name: str) -> None:
        return
    #def get_phase(self, itime: int, case_tuple: str) -> int:
        #(itime, iresult, header) = case_tuple
        #if self.is_real:
            #return 0.0
        #phases = self.phases[self.layer_indices]
        #return phases[itime]
    #def set_phase(self, itime: int, case_tuple: str, phase: float) -> None:
        #(itime, iresult, header) = case_tuple
        #if self.is_real:
            #return
        #phases = self.phases[self.layer_indices]
        #phases[itime] = phase


class CompositeStrainStressResults2(CompositeResults2):
    def __init__(self,
                 subcase_id: int,
                 model: BDF,
                 element_id: np.ndarray,
                 case: RealCompositePlateArray,
                 result: str,
                 title: str,
                 #dim_max: float=1.0,
                 data_format: str='%g',
                 is_variable_data_format: bool=False,
                 nlabels=None, labelsize=None, ncolors=None,
                 colormap: str='',
                 set_max_min: bool=False,
                 uname: str='CompositeStressResults2'):
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
        CompositeResults2.__init__(
            self,
            subcase_id,
            model, element_id,
            case,
            result,
            #dim_max,
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
        elements = np.unique(case.element_layer[:, 0])  #  local element id
        #_ielement = np.searchsorted(element_id, elements)

        # dense -> no missing nodes in the results set
        self.is_dense = (len(element_id) == len(elements))

        #self.xyz = xyz
        #assert len(self.xyz.shape) == 2, self.xyz.shape
        self.location = 'centroid'
        if self.is_stress:
            self.headers = ['CompositeStress2']
        else:
            self.headers = ['CompositeStrain2']
        str(self)

    #-------------------------------------
    # unmodifyable getters

    def get_location(self, unused_i: int, unused_res_name: str) -> str:
        """the result type"""
        return self.location

    #-------------------------------------
    def get_methods(self, itime: int, res_name: str) -> list[str]:
        layers = list(self.layer_map.values())
        return layers

    def get_scalar(self, itime: int, res_name: str, method: str) -> np.ndarray:
        return self.get_plot_value(itime, res_name, method)

    def get_plot_value(self, itime: int, res_name: str, method: str) -> np.ndarray:
        """get_fringe_value"""
        normi = self.get_result(itime, res_name, method, return_dense=False)
        #normi = safe_norm(dxyz, axis=col_axis)
        if self.is_dense:
            return normi

        #case.data.shape = (11, 43, 6)
        #nnodes = len(self.node_id) =  48
        #nnodesi = len(self.inode) = len(self.dxyz.node_gridtype) = 43
        normi2 = np.full(len(self.element_id), np.nan, dtype=normi.dtype)
        normi2[self.ielement] = normi
        return normi2

    #def get_force_vector_result(self, itime: int, res_name: str, method: str) -> np.ndarray:
        #dxyz = self.get_result(itime, res_name, method, return_dense=True)
        #scale = 1.
        #return self.xyz, dxyz * scale

    #def get_vector_result(self, itime: int, res_name: str, method: str) -> tuple[np.ndarray, np.ndarray]:
        #dxyz = self.get_result(itime, res_name, method, return_dense=True)
        #scale = self.get_scale(itime, res_name)
        #deflected_xyz = self.xyz + scale * dxyz
        #return self.xyz, deflected_xyz

    #def get_vector_result_by_scale_phase(self, i: int, unused_name: str,
                                         #scale: float,
                                         #phase: float=0.) -> tuple[np.ndarray, np.ndarray]:
        #"""
        #Gets the real/complex deflection result

        #Parameters
        #----------
        #i : int
            #mode/time/loadstep number
        #name : str
            #unused; useful for debugging
        #scale : float
            #deflection scale factor; true scale
        #phase : float; default=0.0
            #phase angle (degrees); unused for real results

        #Returns
        #-------
        #xyz : (nnodes, 3) float ndarray
            #the nominal state
        #deflected_xyz : (nnodes, 3) float ndarray
            #the deflected state
        #"""
        #assert self.dim == 3, self.dim
        #assert len(self.xyz.shape) == 2, self.xyz.shape
        #if self.is_real:
            #deflected_xyz = self.xyz + scale * self.dxyz[i, :]
        #else:
            #assert isinstance(i, int), (i, phase)
            #assert isinstance(phase, float), (i, phase)
            #dxyz = self._get_complex_displacements_by_phase(i, phase)
            #deflected_xyz = self.xyz + scale * dxyz
        #assert len(deflected_xyz.shape) == 2, deflected_xyz.shape
        #return self.xyz, deflected_xyz

    def __repr__(self) -> str:
        """defines str(self)"""
        msg = 'CompositeStrainStressResults2\n'
        #msg += f'    titles={self.titles!r}\n'
        msg += f'    subcase_id={self.subcase_id}\n'
        msg += f'    data_type={self.data_type!r}\n'
        msg += f'    is_real={self.is_real} is_complex={self.is_complex}\n'
        msg += f'    location={self.location!r}\n'
        msg += f'    header={self.headers!r}\n'
        msg += f'    data_format={self.data_formats!r}\n'
        msg += f'    uname={self.uname!r}\n'
        return msg
