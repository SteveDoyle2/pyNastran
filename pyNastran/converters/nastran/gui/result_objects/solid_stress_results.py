from __future__ import annotations
from collections import defaultdict
import numpy as np
from typing import Optional, TYPE_CHECKING

from pyNastran.utils.mathematics import get_abs_max
from pyNastran.gui.gui_objects.gui_result import GuiResultCommon
from pyNastran.femutils.utils import abs_nan_min_max # , pivot_table,  # abs_min_max
#from pyNastran.bdf.utils import write_patran_syntax_dict

from .displacement_results import VectorResultsCommon
if TYPE_CHECKING:
    from pyNastran.bdf.bdf import BDF
    from pyNastran.op2.tables.oes_stressStrain.real.oes_plates import RealPlateArray


def abs_max_scalar(x: np.ndarray):
    mini = np.nanmin(x)
    maxi = np.nanmax(x)
    if np.abs(mini) > np.abs(maxi):
        return mini
    return maxi

def difference_scalar(x: np.ndarray):
    out = np.nanmax(x) - np.nanmin(x)
    return out

nodal_combine_map = {
    'Absolute Max': abs_max_scalar,
    'Mean': np.nanmean,
    'Max': np.nanmax,
    'Min': np.nanmin,
    'Difference': difference_scalar,
    'Std. Dev.': np.nanstd,
}
class SolidResults2(VectorResultsCommon):
    def __init__(self,
                 subcase_id: int,
                 model: BDF,
                 node_id: np.ndarray,
                 element_id: np.ndarray,
                 cases: list[RealPlateArray],
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
        self.min_max_method = '' # self.has_derivation_transform(i, name)[1]['derivation'][0]
        self.transform = self.has_coord_transform(i, name)[1][0]
        self.nodal_combine = self.has_nodal_combine_transform(i, name)[1][0]
        #assert len(element_id) >= self.case.

        self.is_dense = False
        self.dim = cases[0].data.ndim
        for case in cases:
            assert case.data.ndim == 3, case.data.shape

        self.subcase_id = subcase_id
        self.is_stress = case.is_stress

        self.layer_map = {
            0: 'Corner',  # default
            1: 'Centroid',
        }

        #if dim_max == 0.0:
            #dim_max = 1.0
        #self.dim_max = dim_max
        self.linked_scale_factor = False

        self.data_format = data_format

        #  global ids
        self.model = model
        self.node_id = node_id
        self.element_id = element_id

        # local case object
        self.cases = cases
        self.result = result

        self.data_type = case.data.dtype.str # '<c8', '<f4'
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
        self.headers = ['SolidResult2'] * ntimes

        self.nlabels = None
        self.labelsize = None
        self.ncolors = None
        self.colormap = colormap

        self.uname = uname
        #self.active_method = ''

    def _get_default_tuple_indices(self) -> tuple[int]:
        out = tuple(np.array(self._get_default_layer_indicies()) - 1)
        return out

    def _get_default_layer_indicies(self):
        return (0, )
        #default_indices = list(self.layer_map.keys())
        #default_indices.remove(0)
        #return default_indices

    def set_sidebar_args(self,
                         itime: str, res_name: str,
                         min_max_method: str='', # Absolute Max
                         transform: str='', # Material
                         methods_keys: Optional[list[int]]=None,
                         # unused
                         nodal_combine: str='', # Centroid
                         **kwargs) -> None:
        assert len(kwargs) == 0, kwargs
        transforms = self.has_coord_transform(itime, res_name)[1]
        #min_max_methods = self.has_derivation_transform(itime, res_name)[1]['derivation']
        combine_methods = self.has_nodal_combine_transform(itime, res_name)[1]

        transform = transform if transform else transforms[0]
        #min_max_method = min_max_method if min_max_method else min_max_methods[0]
        nodal_combine = nodal_combine if nodal_combine else combine_methods[0]

        assert transform in transforms, transform
        #assert min_max_method in min_max_methods, min_max_method
        assert nodal_combine in combine_methods, nodal_combine

        #sidebar_kwargs = {
            #'min_max_method': min_max_method,
            #'transform': coord,
            #'nodal_combine': nodal_combine,
            #'methods_keys': keys_b,
        #}
        # if Both is selected, only use Both
        # methods = ['Both', 'Top', 'Bottom']
        default_indices = self._get_default_layer_indicies()
        if methods_keys is None or len(methods_keys) == 0:
            # default; Corner
            indices = default_indices
        elif 0 in methods_keys: # Both
            # Corner takes precendence (it has more data)
            indices = (0, )
        else:
            # Centroid is secondary
            indices = (1, )
        self.layer_indices = indices
        #self.layer_indices = (1, )

        # doesn't matter cause it's already nodal
        #assert min_max_method in min_max_methods, min_max_method
        assert nodal_combine in combine_methods, nodal_combine
        #self.min_max_method = min_max_method
        self.nodal_combine = nodal_combine
        self.transform = transform

    def has_methods_table(self, i: int, res_name: str) -> bool:
        return True
    def has_coord_transform(self, i: int, res_name: str) -> tuple[bool, list[str]]:
        return True, ['Material']
    def has_derivation_transform(self, i: int, case_tuple: str) -> tuple[bool, dict[str, Any]]:
        """min/max/avg"""
        #(itime, iresult, header) = case_tuple
        #out = {
            #'tooltip': 'Method to reduce multiple nodes into a single nodal/elemental value',
            #'derivation': ['Centroid',
                           #'Absolute Max', 'Min', 'Max', 'Mean', 'Std. Dev.', 'Difference',
                           ##'Derive/Average'
                           #],
        #}
        #return True, out
        return False, {}
    def has_nodal_combine_transform(self, i: int, res_name: str) -> tuple[bool, list[str]]:
        """elemental -> nodal (only applies for Corner)"""
        return True, ['Absolute Max', 'Mean', 'Max', 'Min',
                      'Difference', 'Std. Dev.',]
    # 'Nodal Max'

    def get_annotation(self, itime: int, case_tuple: str) -> str:
        """
        A header is the thingy that goes in the lower left corner
        title = 'Solid Stress'
        method = 'Absolute Max'
        header = 'Static'
        nodal_combine = 'Nodal Max'
        returns 'Solid Stress Both (Absolute Max; Nodal Max, Static): sigma11'
        """
        # overwrite itime based on linked_scale factor
        (itime, iresult, header) = case_tuple
        itime, unused_case_flag = self.get_case_flag(case_tuple)

        default_indices = self._get_default_tuple_indices() # 0-based
        if self.layer_indices == default_indices:
            self.layer_indices = (0, )

        if self.layer_indices == (0, ):
            layer_str = self.layer_map[0]  # Corner
        elif self.layer_indices == (1, ):
            layer_str = self.layer_map[1]  # Centroid
        else:
            raise RuntimeError(self.layer_indices)
        self.layer_indices

        result = get_solid_result(self.result, iresult, index=1)

        if layer_str == 'Centroid':
            #'Solid Stress; Centroid (Static): von Mises'
            annotation_label = f'{self.title}; {layer_str} ({header}): {result}'
        else:
            #'Solid Stress; Corner (Mean; Static): von Mises'
            annotation_label = f'{self.title}; {layer_str} ({self.nodal_combine}, {header}): {result}'
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
                                                                       tuple[int, int, tuple, str, str]]:
        """
        itime = 0
        iresult = 0 # o11
        layer_indices = (1, 2)
        min_max_method = 'Absolute Max'
        nodal_combine = 'Centroid'
        """
        (itime, iresult, header) = case_tuple
        #if self.is_linked_scale_factor:
            #itime = 0

        return itime, (itime, iresult, self.layer_indices, self.min_max_method, self.nodal_combine)

    def get_default_legend_title(self, itime: int, case_tuple: str) -> str:
        (itime, iresult, header) = case_tuple
        #method_ = 'Composite Stress Layers:' if self.is_stress else 'Composite Strain Layers:'
        #self.layer_indices
        results = list(self.result.values())
        #method = method_ + ', '.join(str(idx) for idx in (self.layer_indices+1))
        #method = method.strip()
        #title = f'{self.title} {method}'
        title = results[iresult][0]  # sidebar label=legend
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
        result = get_solid_result(self.result, iresult, index=0)
        return result

    def _get_real_data(self, case_tuple: int) -> np.ndarray:
        (itime, iresult, header) = case_tuple

        # [itime, ielement, ilayer, iresult
        #self.centroid_eids = np.hstack(centroid_elements_list)
        #self.centroid_data = np.hstack(data_list)

        #ilayer = self.layer_indices
        if self.layer_indices == (-1, ):
            self.layer_indices = (0, )

        #self.case.get_headers()
        #[fiber_dist, 'oxx', 'oyy', 'txy', 'angle', 'omax', 'omin', ovm]
        #results = list(self.result.keys())
        neids = self.centroid_data.shape[1]

        #['oxx', 'oyy', 'ozz', 'txy', 'tyz', 'txz', 'omax', 'omid', 'omin', von_mises]
        if self.layer_indices == (0, ):
            data = self._get_nodal_result(itime, iresult)

        elif self.layer_indices == (1, ):
            data = self._get_centroid_result(itime, iresult)
        #elif self.nodal_combine == 'Nodal Max':
        #elif self.nodal_combine == 'Nodal Min':
        #elif self.nodal_combine == 'Nodal Mean':
        #elif self.nodal_combine == 'Nodal Abs Max':
        #elif self.nodal_combine == 'Nodal Std. Dev.':
        #elif self.nodal_combine == 'Nodal Difference':
        else:
            raise RuntimeError(self.nodal_combine)

        assert len(data.shape) == 1, data.shape
        return data

        # multiple plies
        # ['Absolute Max', 'Min', 'Max', 'Derive/Average']
        ## TODO: why is this shape backwards?!!!
        ## [ilayer, ielement] ???
        axis = 0
        if self.min_max_method == 'Absolute Max':
            data2 = abs_nan_min_max(data, axis=axis)
        elif self.min_max_method == 'Min':
            data2 = np.nanmin(data, axis=axis)
        elif self.min_max_method == 'Max':
            data2 = np.nanmax(data, axis=axis)
        elif self.min_max_method == 'Mean':  #   (Derive/Average)???
            data2 = np.nanmean(data, axis=axis)
        elif self.min_max_method == 'Std. Dev.':
            data2 = np.nanstd(data, axis=axis)
        elif self.min_max_method == 'Difference':
            data2 = np.nanmax(data, axis=axis) - np.nanmin(data, axis=axis)
        #elif self.min_max_method == 'Max Over Time':
            #data2 = np.nanmax(data, axis=axis) - np.nanmin(data2, axis=axis)
        #elif self.min_max_method == 'Derive/Average':
            #data2 = np.nanmax(data, axis=1)
        else:  # pragma: no cover
            raise NotImplementedError(self.min_max_method)

        # TODO: hack to try and debug things...
        assert data2.shape == (neids, )
        #data4 = eids_new.astype('float32')
        return data2

    def _get_centroid_result(self, itime: int,
                             iresult: Union[int, str]) -> np.ndarray:
        """
        Centroid
        Derive and no Averaging
        """
        ioxx = 0
        ioyy = 1
        iozz = 2
        itxy = 3
        ityz = 4
        itxz = 5
        imax = 6
        imin = 8
        ## nodal_combine == 'Centroid':

        if iresult == 'abs_principal': # abs max; should be good
            omax = self.centroid_data[itime, :, imax]
            omin = self.centroid_data[itime, :, imin]
            data = get_abs_max(omin, omax, dtype=omin.dtype)

        elif iresult == 'von_mises':
            # von mises
            #Derive
            oxx = self.centroid_data[itime, :, ioxx]
            oyy = self.centroid_data[itime, :, ioyy]
            ozz = self.centroid_data[itime, :, iozz]
            txy = self.centroid_data[itime, :, itxy]
            txz = self.centroid_data[itime, :, itxz]
            tyz = self.centroid_data[itime, :, ityz]
            data = von_mises_3d(oxx, oyy, ozz, txy, tyz, txz)
        elif iresult == 'max_shear': #  probably wrong
            # not checked for strain
            omax = self.centroid_data[itime, :, imax]
            omin = self.centroid_data[itime, :, imin]
            data = max_shear(omax, omin)
        else:
            data = self.centroid_data[itime, :, iresult].copy()
        return data

    def _get_nodal_result(self, itime: int,
                          iresult: Union[int, str]) -> np.ndarray:
        ioxx = 0
        ioyy = 1
        iozz = 2
        itxy = 3
        ityz = 4
        itxz = 5
        imax = 6
        imin = 8

        ## Corner
        ## TODO: consider implementing Average/Derive
        # ----------------------------------------------------------
        # setup
        element_node = self.element_node
        nids = np.unique(element_node[:, 1])

        nodal_combine_func = nodal_combine_map[self.nodal_combine]
        inid = np.searchsorted(self.node_id, nids)
        nid_to_inid_map = {nid: inidi for nid, inidi in zip(nids, inid)}

        ## Derive/Average
        if isinstance(iresult, int):
            datai = self.node_data[itime, :, iresult]
            data = nodal_average(
                nodal_combine_func,
                element_node, datai,
                nids, inid, nid_to_inid_map)

        elif iresult == 'von_mises':
            oxx = self.node_data[itime, :, ioxx]
            oyy = self.node_data[itime, :, ioyy]
            ozz = self.node_data[itime, :, iozz]
            txy = self.node_data[itime, :, itxy]
            txz = self.node_data[itime, :, itxz]
            tyz = self.node_data[itime, :, ityz]
            ovm_data = von_mises_3d(oxx, oyy, ozz, txy, tyz, txz)
            data = nodal_average(
                nodal_combine_func,
                element_node, ovm_data,
                nids, inid, nid_to_inid_map)

        elif iresult == 'max_shear':
            omax = self.node_data[itime, :, imax]
            omin = self.node_data[itime, :, imin]
            max_shear_data = max_shear(omax, omin)
            data = nodal_average(
                nodal_combine_func,
                element_node, max_shear_data,
                nids, inid, nid_to_inid_map)

        elif iresult == 'abs_principal':
            data_max = self.node_data[itime, :, imax]
            data_min = self.node_data[itime, :, imin]
            abs_data = get_abs_max(data_min, data_max, dtype=data_min.dtype)
            assert abs_data.shape == data_min.shape
            data = nodal_average(
                nodal_combine_func,
                element_node, abs_data,
                nids, inid, nid_to_inid_map)
        else:
            raise NotImplementedError(iresult)
        return data

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
        #method = self._update_method(itime, case_tuple, method)
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
        result_out[self.ielement_centroid] = data
        return result_out

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


class SolidStrainStressResults2(SolidResults2):
    def __init__(self,
                 subcase_id: int,
                 model: BDF,
                 node_id: np.ndarray,
                 element_id: np.ndarray,
                 cases: list[RealPlateArray],
                 result: str,
                 title: str,
                 #is_fiber_distance: bool,
                 #dim_max: float=1.0,
                 data_format: str='%g',
                 is_variable_data_format: bool=False,
                 nlabels=None, labelsize=None, ncolors=None,
                 colormap: str='',
                 set_max_min: bool=False,
                 uname: str='SolidStressStrainResults2'):
        """
        Defines a SolidStressResults/SolidStrainResults result

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
        SolidResults2.__init__(
            self,
            subcase_id,
            model, node_id, element_id,
            cases,
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
        centroid_data_list = []
        centroid_elements_list = []

        node_data_list = []
        element_node_list = []
        for case in cases:
            ntime, nelement_nnode, nresults = case.data.shape
            nelements_all = len(np.unique(case.element_node[:, 0]))
            nnode = case.nnodes_per_element
            if case.nnodes_per_element != 1:
                eids = case.element_node[0::nnode, 0]
                nelement = len(case.element_node) // nnode
                nnodal_node = nelement * (nnode - 1)

                element_node_3d = case.element_node.reshape(nelement, nnode, 2)
                element_node = element_node_3d[:, 1:, :].reshape(nnodal_node, 2)
                eids = element_node_3d[:, 0, 0]
                data = case.data.copy().reshape(ntime, nelement, nnode, nresults)
                centroid_elements_list.append(eids)
                #nodal_elements.append(case.element_node[1::nnodes_per_element, 0])
                assert nelement == nelements_all, (nelement, nelements_all)
                node_data = data[:, :, 1:, :]
                node_data2 = node_data.reshape(ntime, nnodal_node, nresults)
            else:
                eids = case.element_node[:, 0]
                centroid_elements_list.append(case.element_node)
                data = case.data.reshape(ntime, nelements_all, 1, nresults)
                node_data = data
                raise NotImplementedError('centroidal solid elements')
            # slice off the centroid
            centroid_datai = data[:, :, 0, :]
            centroid_data_list.append(centroid_datai)

            node_data_list.append(node_data2)
            element_node_list.append(element_node)

        self.centroid_eids = np.hstack(centroid_elements_list)
        # [ntimes, nelements, nresults]
        self.centroid_data = np.hstack(centroid_data_list)

        self.element_node = np.vstack(element_node_list)
        self.node_data = np.hstack(node_data_list)

        common_eids = np.intersect1d(self.centroid_eids, element_id)
        if len(common_eids) == 0:
            raise IndexError('no solid elements found...')
        elif len(common_eids) != len(self.centroid_eids):
            icommon = np.searchsorted(common_eids, self.centroid_eids)
            #self.centroid_data = self.centroid_data[:, icommon, :]
            raise RuntimeError('some common elements were found...but some are missing')


        self.ielement_centroid = np.searchsorted(element_id, self.centroid_eids)

        # dense -> no missing nodes in the results set
        self.is_dense = (len(element_id) == len(self.centroid_eids))
        #self.is_dense = False

        #self.xyz = xyz
        #assert len(self.xyz.shape) == 2, self.xyz.shape
        #self.location = 'centroid'
        if self.is_stress:
            self.headers = ['SolidStress2']
        else:
            self.headers = ['SolidStrain2']
        str(self)

    #-------------------------------------
    # unmodifyable getters

    def get_location(self, unused_i: int, unused_res_name: str) -> str:
        """the result type"""
        if self.layer_indices == (0, ):
            return 'node'
        return 'centroid'
        #return self.location

    #-------------------------------------
    def get_methods(self, itime: int, res_name: str) -> list[str]:
        layers = list(self.layer_map.values())
        return layers

    def get_scalar(self, itime: int, res_name: str, method: str) -> np.ndarray:
        return self.get_plot_value(itime, res_name, method)

    def get_plot_value(self, itime: int, res_name: str, method: str) -> np.ndarray:
        """get_fringe_value"""
        normi = self.get_result(itime, res_name, method, return_dense=False)
        if self.is_dense:
            return normi

        #case.data.shape = (11, 43, 6)
        #nnodes = len(self.node_id) =  48
        #nnodesi = len(self.inode) = len(self.dxyz.node_gridtype) = 43
        normi2 = np.full(len(self.element_id), np.nan, dtype=normi.dtype)
        normi2[self.ielement_centroid] = normi
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
        msg = 'SolidStrainStressResults2\n'
        #msg += f'    titles={self.titles!r}\n'
        msg += f'    subcase_id={self.subcase_id}\n'
        msg += f'    data_type={self.data_type!r}\n'
        msg += f'    is_real={self.is_real} is_complex={self.is_complex}\n'
        #msg += f'    location={self.location!r}\n'
        msg += f'    header={self.headers!r}\n'
        msg += f'    data_format={self.data_formats!r}\n'
        msg += f'    uname={self.uname!r}\n'
        return msg


def get_solid_result(result: dict[str, Any],
                     iresult: Union[int, str], index: int):
    """
    values
    0=title, 'annotation'
    ('sAbs Principal', 'Abs Principal')
    """
    assert index in (0, 1), index
    #if isinstance(iresult, int):
        #assert iresult >= 0, iresult
    results = result[iresult]
    return results[index]
    #elif iresult == 'von_mises':
        #word = 'von Mises'
    #elif iresult == 'max_shear':
        #word = 'Max Shear'
    #elif iresult == 'abs_principal':
        #word = 'Abs Principal'
    #else:
        #raise RuntimeError(iresult)
    #return word

def max_shear(omax, omin) -> np.ndarray:
    """
    not verified for stress/strain
    same as Tresca?
    """
    max_shear = (omax - omin) / 2.
    return max_shear

def von_mises_2d(oxx, oyy, txy) -> np.ndarray:
    """not verified for stress/strain"""
    ovm = np.sqrt(oxx**2 + oyy**2 - oxx*oyy +3*(txy**2) )
    return ovm

def von_mises_3d(oxx, oyy, ozz, txy, tyz, txz) -> np.ndarray:
    """not verified for stress/strain"""
    vm = np.sqrt(
        0.5 * ((oxx - oyy) ** 2 + (oyy - ozz) ** 2 +(oxx - ozz) ** 2) +
        3 * (txy**2 + tyz**2 + txz**2))
    return vm

def nodal_average(nodal_combine_func: Callable[np.ndarray],
                  element_node: np.ndarray,
                  data: np.ndarray,
                  nids: np.ndarray,
                  inid: np.ndarray,
                  nid_to_inid_map: dict[int, int]) -> np.ndarray:
    data_dict = defaultdict(list)
    nnode = len(nids)

    data2 = np.full(nnode, np.nan, dtype=data.dtype)
    for (eid, nid), datai in zip(element_node, data):
        data_dict[nid].append(datai)
    for nid, datasi in data_dict.items():
        collapsed_value = nodal_combine_func(datasi)
        inidi = nid_to_inid_map[nid]
        data2[inidi] = collapsed_value
    return data2
