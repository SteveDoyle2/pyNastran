from __future__ import annotations
#from collections import defaultdict
import numpy as np
from typing import Optional, Any, TYPE_CHECKING

from pyNastran.utils.mathematics import get_abs_max
from pyNastran.femutils.utils import abs_nan_min_max, safe_nanstd

from pyNastran.gui.gui_objects.vector_results import VectorResultsCommon, filter_ids
from pyNastran.op2.stress_reduction import von_mises_2d, max_shear
from .nodal_averaging import nodal_average, nodal_combine_map, derivation_map

if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf import BDF
    from pyNastran.op2.op2 import OP2
    from pyNastran.op2.tables.oes_stressStrain.real.oes_plates import RealPlateArray


col_axis = 1
DERIVATION_METHODS = list(derivation_map.keys())

NODAL_COMBINE_TRANSFORMS = [
    'Centroid', 'Absolute Max', 'Mean', 'Min', 'Max',
    'Std. Dev.', 'Difference']
NODAL_COMBINES = []
class PlateStrainStressResults2(VectorResultsCommon):
    def __init__(self,
                 subcase_id: int,
                 model: BDF,
                 node_id: np.ndarray,
                 element_id: np.ndarray,
                 cases: list[RealPlateArray],
                 result: str,
                 title: str,
                 is_fiber_distance: bool,
                 eid_to_nid_map: dict[int, list[int]],
                 #dim_max: float=1.0,
                 data_format: str='%g',
                 is_variable_data_format: bool=False,
                 nlabels=None, labelsize=None, ncolors=None,
                 colormap: str='',
                 set_max_min: bool=False,
                 uname: str='PlateStressStrainResults2'):
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
        assert isinstance(is_fiber_distance, bool), is_fiber_distance
        title = ''
        ntitles = 8
        VectorResultsCommon.__init__(
            self, subcase_id, title,
            cases,
            ntitles,
            data_format=data_format,
            is_variable_data_format=is_variable_data_format,
            nlabels=nlabels, labelsize=labelsize, ncolors=ncolors,
            colormap=colormap,
            #set_max_min: bool=False,
            uname=uname)
        self.layer_indices = (-1, )  # All
        self.is_fiber_distance = is_fiber_distance
        i = -1
        name = (0, 0, '')

        # slice off the methods (from the boolean) and then pull the 0th one
        self.min_max_method = self.has_derivation_transform(i, name)[1]['derivation'][0]
        self.transform = self.has_coord_transform(i, name)[1][0]
        self.nodal_combine = self.has_nodal_combine_transform(i, name)[1][0]

        self.dim = cases[0].data.ndim
        for case in cases:
            assert case.data.ndim == 3, case.data.shape

        self.is_stress = case.is_stress
        self.eid_to_nid_map = eid_to_nid_map
        if self.is_stress:
            self.iresult_map = {
                #0 : 'FiberCurvature',
                0 : 'FiberDistance',
                1 : 'XX Stress',
                2 : 'YY Stress',
                3 : 'XY Stress',
                4 : 'Theta',
                5 : 'MaxPrincipal Stress',
                6 : 'MinPrincipal Stress',
                7 : 'AbsPrincipal Stress',  # special
                8 : 'von Mises Stress',     # special
                9 : 'Max Shear Stress',     # special
            }
        else:
            self.iresult_map = {
                #0 : 'FiberCurvature',
                0 : 'FiberDistance',
                1 : 'XX Strain',
                2 : 'YY Strain',
                3 : 'XY Strain',
                4 : 'Theta',
                5 : 'MaxPrincipal Strain',
                6 : 'MinPrincipal Strain',
                7 : 'AbsPrincipal Strain', # special
                8 : 'von Mises Strain',    # special
                9 : 'Max Shear Strain',    # special
            }
        if self.is_fiber_distance:
            self.layer_map = {
                0: 'Both',
                1: 'Bottom',
                2: 'Top',
            }
        else:
            self.layer_map = {
                0: 'Both',  # wut?  This makes no sense...
                1: 'Mean',
                2: 'Slope',
            }
            self.iresult_map[0] = 'FiberCurvature'

        self.data_format = data_format

        #  global ids
        self.model = model
        self.node_id = node_id
        self.element_id = element_id

        # local case object
        #self.cases = cases
        assert isinstance(result, dict), result
        self.result = result

        self.data_type = case.data.dtype.str # '<c8', '<f4'
        self.is_real = True if self.data_type in ['<f4', '<f8'] else False
        #self.is_real = dxyz.data.dtype.name in {'float32', 'float64'}
        self.is_complex = not self.is_real

        #ntimes = case.data.shape[0]
        #self.headers = ['PlateResult2'] * ntimes

        self.location = 'centroid'

        #[ntime, nelement, nlayer, nresult]
        self.centroid_data = np.zeros((0, 0, 2, 8), dtype='float32')
        self.node_data = np.zeros((0, 0, 2, 8), dtype='float32')
        self.element_node = np.zeros((0, 2), dtype='int32')
        self.inode = np.zeros(0, dtype='int32')
        self.ielement_centroid = np.zeros(0, dtype='int32')
        #---------------------------------------------------------------------

        #linked_scale_factor = False

        out = setup_centroid_node_data(eid_to_nid_map, cases, element_id)
        centroid_eids, centroid_data, element_node, node_data = out

        # [ntime, nelement_nnode, nlayer, nresult]
        assert centroid_data.ndim == 4, centroid_data.shape
        assert node_data.ndim == 4, node_data.shape

        self.centroid_eids = centroid_eids
        # [ntime, nelement_nnode, nlayer, nresult]
        self.centroid_data = centroid_data

        # [ntime, nelement_nnode, nlayer, nresult]
        self.element_node = element_node
        self.node_data = node_data

        assert len(np.unique(self.centroid_eids)) == len(self.centroid_eids)

        common_eids = np.intersect1d(self.centroid_eids, element_id)
        if len(common_eids) == 0:
            raise IndexError('no plate elements found...')
        elif len(common_eids) != len(self.centroid_eids):
            icommon = np.searchsorted(common_eids, self.centroid_eids)
            #self.centroid_data = self.centroid_data[:, icommon, :]
            raise RuntimeError('some common elements were found...but some are missing')

        self.ielement_centroid = np.searchsorted(element_id, self.centroid_eids)

        nids = np.unique(self.element_node[:, 1])
        self.inode = np.searchsorted(node_id, nids)

        # dense -> no missing nodes in the results set
        self.is_dense = (
            (len(element_id) == len(self.centroid_eids)) and
            (len(node_id) == len(nids))
        )
        #self.is_dense = False

        #self.xyz = xyz
        #assert len(self.xyz.shape) == 2, self.xyz.shape
        ntimes = self.centroid_data.shape[0]
        if self.is_stress:
            self.headers = ['PlateStress2'] * ntimes
        else:
            self.headers = ['PlateStrain2'] * ntimes
        str(self)

    @classmethod
    def load_from_code(cls,
                       subcase_id: int,
                       model: BDF,
                       model_results: OP2,
                       element_id: np.ndarray,
                       #cases: list[RealPlateArray],
                       #subcase_key: int,
                       is_stress: bool,
                       eid_to_nid_map: dict[int, list[int]],
                       #is_variable_data_format: bool=False,
                       #prefix='stress',
                       require_results: bool=True):
        """the intention of this method is to not use the gui"""
        subcase_key0 = (subcase_id, -1)
        plate_cases, subcase_id = get_real_plate_cases(
            element_id,
            model_results,
            #subcase_key0,
            subcase_id,
            is_stress=True,
            prefix='',
            require_results=require_results)

        # get the xyzs
        out = model.get_xyz_in_coord_array(
            cid=0, fdtype='float64', idtype='int64')
        nid_cp_cd, xyz_cid, unused_xyz_cp, unused_icd_transform, unused_icp_transform = out
        node_id = nid_cp_cd[:, 0]
        #------------------------------------------
        out = plate_cases_to_iresult(plate_cases, is_stress)
        iresult_to_title_annotation_map, is_fiber_distance, word, max_sheari = out

        #------------------------------------------
        title = 'title'
        obj = PlateStrainStressResults2(
            subcase_id,
            model,
            node_id,
            element_id,
            plate_cases,
            iresult_to_title_annotation_map, title,
            is_fiber_distance,
            eid_to_nid_map,
            is_variable_data_format=False,  # ???
            set_max_min=False)
        return obj

    def set_centroid(self,
                     top_bottom_both: str='Both',
                     min_max_method: str='Absolute Max') -> None:
        """elemental"""
        self._set_top_bottom_both(top_bottom_both)
        self._set_derivation_method(min_max_method)
        self.transform = 'Material'
        self.nodal_combine = 'Centroid'

    def set_corner(self,
                   top_bottom_both: str='Both',
                   min_max_method: str='Absolute Max',
                   nodal_combine_method: str='Mean') -> None:
        """elemental -> nodal"""
        self._set_top_bottom_both(top_bottom_both)
        self._set_derivation_method(min_max_method)
        assert nodal_combine_method.title() != 'Centroid', f'nodal_combine_method={nodal_combine_method}'
        self.transform = 'Material'
        self._set_nodal_combine_method(nodal_combine_method)

    def _set_top_bottom_both(self, top_bottom_both: str) -> None:
        top_bottom_both = top_bottom_both.lower()
        assert top_bottom_both in {'top', 'bottom', 'both'}, top_bottom_both
        if top_bottom_both == 'both':
            layer_indices = (0, 1)
        elif top_bottom_both == 'bottom':
            layer_indices = (0,)
        elif top_bottom_both == 'top':
            layer_indices = (1,)
        self.layer_indices = layer_indices

    def _set_nodal_combine_method(self, nodal_combine_method: str) -> None:
       nodal_combine_method = nodal_combine_method.title()
       assert nodal_combine_method in NODAL_COMBINE_TRANSFORMS, f'nodal_combine_method={nodal_combine_method!r}'
       self.nodal_combine = nodal_combine_method

    def _set_derivation_method(self, min_max_method: str) -> None:
        """elemental -> nodal"""
        min_max_method = min_max_method.title()
        assert min_max_method in DERIVATION_METHODS, f'min_max_method={min_max_method!r}'
        self.min_max_method = min_max_method

    def _get_default_tuple_indices(self):
        out = tuple(np.array(self._get_default_layer_indicies()) - 1)
        return out

    def _get_default_layer_indicies(self):
        default_indices = list(self.layer_map.keys())
        default_indices.remove(0)
        return default_indices

    def set_sidebar_args(self,
                         itime: str, res_name: str,
                         min_max_method: str='', # Absolute Max
                         transform: str='', # Material
                         methods_keys: Optional[list[int]]=None,
                         # unused
                         nodal_combine: str='', # Centroid
                         **kwargs) -> None:
        assert len(kwargs) == 0, kwargs
        transform, min_max_method, nodal_combine = _check_sidebar_args(
            self, itime, res_name,
            transform, min_max_method, nodal_combine)

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
            # default; All
            indices = default_indices
        elif 0 in methods_keys: # Both
            # include all components b/c All is selected
            indices = default_indices
        else:
            # no 0 (Both) in methods_keys
            # update the indices to correspond to the array
            #methods_keys.sort()
            indices = methods_keys
        self.layer_indices = tuple(np.array(indices, dtype='int32') - 1)
        #self.layer_indices = (1, )

        # doesn't matter cause it's already nodal
        self.min_max_method = min_max_method
        self.nodal_combine = nodal_combine
        self.transform = transform

    def has_methods_table(self, i: int, res_name: str) -> bool:
        return True
    def has_coord_transform(self, i: int, res_name: str) -> tuple[bool, list[str]]:
        return True, ['Material']
    def has_derivation_transform(self, i: int, case_tuple: tuple) -> tuple[bool, list[str]]:
        """min/max/avg"""
        if not isinstance(case_tuple, tuple):
            return True, {}
        #(itime, iresult, header) = case_tuple
        imax = 5
        imin = 6
        iresult = case_tuple[1]
        mini = ['Max'] if iresult != imin else []
        maxi = ['Min'] if iresult != imax else []
        derivation = ['Absolute Max'] + mini + maxi + [
            'Mean', 'Std. Dev.', 'Difference']
        #'Derive/Average'

        out = {
            'tooltip': 'Method to reduce multiple layers (top/btm) into a single nodal/elemental value',
            'derivation': derivation,
        }
        return True, out
    def has_nodal_combine_transform(self, i: int, res_name: str) -> tuple[bool, list[str]]:
        """elemental -> nodal"""
        return True, NODAL_COMBINE_TRANSFORMS

    def has_output_checks(self, i: int, resname: str) -> tuple[bool, bool, bool,
                                                               bool, bool, bool]:
        is_enabled_fringe = is_checked_fringe = True
        is_enabled_disp = is_checked_disp = is_enabled_vector = is_checked_vector = False
        out = (is_enabled_fringe, is_checked_fringe,
               is_enabled_disp, is_checked_disp,
               is_enabled_vector, is_checked_vector)
        return out
    # --------------------------------------------------------------
    def get_location(self, itime: int, case_tuple: str) -> str:
        if self.nodal_combine == 'Centroid':
            return 'centroid'
        return 'node'

    def get_annotation(self, itime: int, case_tuple: str) -> str:
        """
        A header is the thingy that goes in the lower left corner
        title = 'Plate Stress'
        method = 'Absolute Max'
        header = 'Static'
        nodal_combine = 'Nodal Max'
        returns 'Plate Stress Both (Absolute Max; Nodal Max, Static): sigma11'
        """
        # overwrite itime based on linked_scale factor
        (itime, iresult, header) = case_tuple
        itime, unused_case_flag = self.get_case_flag(itime, case_tuple)

        default_indices = self._get_default_tuple_indices() # 0-based
        if self.layer_indices == default_indices:
            layer_str = 'Both'
        else:
            if self.layer_indices == (0, ):
                layer_str = self.layer_map[0]  # Bottom
            elif self.layer_indices == (1, ):
                layer_str = self.layer_map[1]  # Top
            else:
                raise RuntimeError(self.layer_indices)
            self.layer_indices

        result = get_plate_result(self.result, iresult, index=1)

        #'Compostite Plate Stress (Absolute Max; Static): sigma11'
        #annotation_label = f'{self.title}; {layer_str} ({self.min_max_method}, {self.nodal_combine}, {header}): {result}'

        title = self.get_legend_title(itime, case_tuple)
        annotation_label = f'Plate {title}; {layer_str} ({self.min_max_method}, {self.nodal_combine}, {header})'

        #return self.uname
        return annotation_label

    def get_case_flag(self, i,
                      case_tuple: tuple[int, int, str]) -> tuple[int,
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
        #method = method_ + ', '.join(str(idx) for idx in (self.layer_indices+1))
        #method = method.strip()
        #title = f'{self.title} {method}'

        # sidebar label=legend
        #('Normal XX', 'XX')
        title_label = self.result[iresult]
        title = title_label[0]
        return title

    def get_legend_tuple(self, itime: int, case_tuple: str) -> int:
        (itime, iresult, header) = case_tuple
        return iresult

    def get_legend_title(self, itime: int, case_tuple: str):
        """Composite Stress Layers: 1, 2, 3, 4"""
        (itime, iresult, header) = case_tuple
        #method_ = 'Composite Stress Layers:' if self.is_stress else 'Composite Strain Layers:'
        #self.layer_indices
        #self.result
        #method = method_ + ', '.join(str(idx) for idx in (self.layer_indices+1))
        #title = f'{self.title} {method}'
        result = get_plate_result(self.result, iresult, index=0)
        return result

    def _get_real_data(self, case_tuple: int) -> np.ndarray:
        (itime, iresult, header) = case_tuple

        # [itime, ielement, ilayer, iresult
        #self.centroid_eids = np.hstack(centroid_elements_list)
        #self.centroid_data = np.hstack(data_list)

        if self.layer_indices == (-1, ):
            self.layer_indices = (0, 1)
        ilayer = self.layer_indices

        #self.case.get_headers()
        #[fiber_dist, 'oxx', 'oyy', 'txy', 'angle', 'omax', 'omin', ovm]
        #results = list(self.result.keys())

        neids = self.centroid_data.shape[1]

        if self.nodal_combine == 'Centroid':
            data = self._get_centroid_result(itime, iresult, ilayer)
            assert len(data.shape) == 2, data.shape
            data2 = self._squash_layers(data, self.min_max_method)
            assert data2.shape == (neids, )
        else:
            data2 = self._get_nodal_result(itime, iresult, ilayer)
        #elif self.nodal_combine == 'Nodal Max':
        #elif self.nodal_combine == 'Nodal Min':
        #elif self.nodal_combine == 'Nodal Mean':
        #elif self.nodal_combine == 'Nodal Abs Max':
        #elif self.nodal_combine == 'Nodal Std. Dev.':
        #elif self.nodal_combine == 'Nodal Difference':
        #else:
            #raise RuntimeError(self.nodal_combine)

        # TODO: hack to try and debug things...
        #data4 = eids_new.astype('float32')
        return data2

    def _squash_layers(self, data: np.ndarray, min_max_method: str) -> np.ndarray:
        # multiple plies
        # ['Absolute Max', 'Min', 'Max', 'Derive/Average']
        ## TODO: why is this shape backwards?!!!
        ## [ilayer, ielement] ???
        axis = 0
        if min_max_method == 'Absolute Max':
            data2 = abs_nan_min_max(data, axis=axis)
        elif min_max_method == 'Min':
            data2 = np.nanmin(data, axis=axis)
        elif min_max_method == 'Max':
            data2 = np.nanmax(data, axis=axis)
        elif min_max_method == 'Mean':  #   (Derive/Average)???
            data2 = np.nanmean(data, axis=axis)
        elif min_max_method == 'Std. Dev.':
            data2 = safe_nanstd(data, axis=axis)
        elif min_max_method == 'Difference':
            data2 = np.nanmax(data, axis=axis) - np.nanmin(data, axis=axis)
        #elif min_max_method == 'Max Over Time':
            #data2 = np.nanmax(data, axis=axis) - np.nanmin(data2, axis=axis)
        #elif min_max_method == 'Derive/Average':
            #data2 = np.nanmax(data, axis=1)
        else:  # pragma: no cover
            raise NotImplementedError(min_max_method)
        return data2

    #def _get_complex_data(self, itime: int) -> np.ndarray:
        #return self._get_real_data(itime)

    def _get_nodal_result(self, itime: int,
                          iresult: int | str,
                          ilayer: np.ndarray) -> np.ndarray:
        #(ntime, neidsi_nnode, nlayer, nresult) = node_data.shape
        node_data = self.node_data
        element_node = self.element_node

        # make sure we never have an issue with these
        nodal_combine_func = nodal_combine_map[self.nodal_combine]
        derivation_func = derivation_map[self.min_max_method]

        nids = np.unique(element_node[:, 1])
        nnode = len(nids)
        nid_to_inid_map = {nid: i for i, nid in enumerate(nids)}

        #ioxx = 1
        #ioyy = 2
        #itxy = 3
        #imax = 5
        #imin = 6
        ## slice off the top/bottom layers
        if isinstance(iresult, int):
            data = node_data[itime, :, ilayer, iresult]
        elif iresult == 'abs_principal':
            omax = node_data[itime, :, ilayer, 5]
            omin = node_data[itime, :, ilayer, 6]
            abs_principal = get_abs_max(omin, omax, dtype=omin.dtype)
            data = abs_principal
        elif iresult == 'von_mises':
            oxx = node_data[itime, :, ilayer, 1]
            oyy = node_data[itime, :, ilayer, 2]
            txy = node_data[itime, :, ilayer, 3]
            data = von_mises_2d(oxx, oyy, txy, self.is_stress)
        elif iresult == 'max_shear':
            omax = node_data[itime, :, ilayer, 5]
            omin = node_data[itime, :, ilayer, 6]
            data = max_shear(omax, omin)
        else:  # pragma: no cover
            raise RuntimeError(iresult)
        assert data.ndim == 2, data.shape

        # ----------------------------------------------------------------------
        ## derivation step
        # squash the data from multiple layers into 1
        if len(ilayer) == 1:
            # we have one layer, so it's easy :)
            if self.min_max_method in {'Std. Dev.', 'Difference'}:
                # single layer -> by definition=0
                derived_data = np.zeros(data[0, :].shape, dtype=data.dtype)
            else:
                # Absolute Max, Min, Max, Mean
                derived_data = data[0, :]
        else:
            # 2 layers
            derived_data = derivation_func(data, axis=0)
        element_node2 = element_node[::2, :]

        # ----------------------------------------------------------------------
        ## nodal combine step
        # time to nodal average
        data2 = nodal_average(
            nodal_combine_func,
            element_node2, derived_data,
            nids, nid_to_inid_map)
        assert len(data2) == nnode, (len(data2), nnode)
        return data2

    def _get_centroid_result(self, itime: int,
                             iresult: int | str,
                             ilayer: np.ndarray) -> np.ndarray:
        centroid_data = self.centroid_data
        #(ntime, neidsi_nnode, nlayer, nresult) = centroid_data.shape
        #assert neids == neidsi_nnode

        # [itime, ielement, ilayer, iresult]
        if isinstance(iresult, int):
            data = centroid_data[itime, :, ilayer, iresult].copy()
        elif iresult == 'abs_principal': # abs max
            omax = centroid_data[itime, :, ilayer, 5]
            omin = centroid_data[itime, :, ilayer, 6]
            data = get_abs_max(omin, omax, dtype=omin.dtype)
            #'exx' : ('Strain XX', 1),
            #'eyy' : ('Strain YY', 2),
            #'exy' : ('Strain XY', 3),
        elif iresult == 'von_mises': # von mises
            oxx = centroid_data[itime, :, ilayer, 1]
            oyy = centroid_data[itime, :, ilayer, 2]
            txy = centroid_data[itime, :, ilayer, 3]
            data = von_mises_2d(oxx, oyy, txy, self.is_stress)
        elif iresult == 'max_shear':
            # not checked for strain
            omax = centroid_data[itime, :, ilayer, 5]
            omin = centroid_data[itime, :, ilayer, 6]
            data = max_shear(omax, omin)
        else:
            raise RuntimeError(iresult)
        return data

    def _get_fringe_data_sparse(self, itime: int, case_tuple: str) -> np.ndarray:
        """gets the sparse fringe"""
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
        return data

    def _get_fringe_data_dense(self, itime: int, case_tuple: str) -> np.ndarray:
        """gets the dense fringe"""
        data = self._get_fringe_data_sparse(itime, case_tuple)
        if self.is_dense:
            return data

        location = self.get_location(0, '')
        if location == 'node':
            nnode = len(self.node_id)
            result_out = np.full(nnode, np.nan, dtype=data.dtype)
            result_out[self.inode] = data
        else:
            nelement = len(self.element_id)
            result_out = np.full(nelement, np.nan, dtype=data.dtype)
            result_out[self.ielement_centroid] = data
        return result_out

    def get_location_arrays(self) -> tuple[np.ndarray, np.ndarray]:
        if self.nodal_combine == 'Centroid':
            all_ids = self.element_id
            ids = np.unique(self.element_node[:, 0])
        else:  # Corner
            all_ids = self.node_id
            ids = np.unique(self.element_node[:, 1])
        return all_ids, ids

    def get_fringe_result(self, itime: int, res_name: str) -> np.ndarray:
        """get_fringe_value"""
        fringe = self._get_fringe_data_dense(itime, res_name)
        return fringe
    def get_fringe_vector_result(self, itime: int, res_name: str) -> tuple[np.ndarray, None]:
        """get_fringe_value"""
        fringe = self.get_fringe_result(itime, res_name)
        return fringe, None

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

    #-------------------------------------
    # unmodifyable getters

    #def get_location(self, unused_i: int, unused_res_name: str) -> str:
        #"""the result type"""
        #return self.location

    #-------------------------------------
    def get_methods(self, itime: int, res_name: str) -> list[str]:
        layers = list(self.layer_map.values())
        return layers

    #def get_force_vector_result(self, itime: int, res_name: str) -> np.ndarray:
        #fringe, dxyz = self.get_fringe_vector_result(itime, res_name)
        #scale = 1.
        #return self.xyz, dxyz * scale

    #def get_vector_result(self, itime: int, res_name: str) -> tuple[np.ndarray, np.ndarray]:
        #frnige, dxyz = self.get_fringe_vector_result(itime, res_name)
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


def get_plate_result(result: dict[str, Any],
                     iresult: int | str, index: int):
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


def setup_centroid_node_data(eid_to_nid_map: dict[int, list[int]],
                             cases: list[RealPlateArray],
                             all_element_id: np.ndarray) -> tuple[np.ndarray, np.ndarray,
                                                                  np.ndarray, np.ndarray]:
    """
    Parameters
    ----------
    all_element_id : (nelement_all,) int array
        all the element ids

    Returns
    -------
    centroid_eids : (nelement, ) int array
        the element ids
    node_data : (ntime, nelement, nlayer, nresult) float array
        nlayers = 2
    element_node : (nelement, 2) int array
        centroids have been removed
    centroid_data : (ntime, nelement, nlayer, nresult) float array
        nlayers = 2

    """
    # setup the node mapping
    centroid_elements_list = []
    centroid_data_list = []

    element_node_list = []
    node_data_list = []
    nlayer = 2
    for case in cases:
        ntime, nelement_nnode_nlayer, nresult = case.data.shape
        nelement_nnode_og = nelement_nnode_nlayer // 2
        if case.is_bilinear():
            nnode = case.nnodes_per_element
            nelement = len(case.element_node) // (2 * nnode)

            # removed the centroid
            nplies = nelement * (nnode - 1) * nlayer

            centroid_eidsi_og = case.element_node[0::2*nnode, 0]
            centroid_eidsi, ieid_filter, nelement_filtered, is_filter = filter_ids(
                all_element_id, centroid_eidsi_og)
            if nelement_filtered == 0:
                continue

            element_node_4d = case.element_node.reshape(nelement, nnode, nlayer, 2)
            if is_filter:
                element_node_4d = element_node_4d[ieid_filter, :, :, :]
            element_node_3d = element_node_4d[:, 1:, :, :]
            element_nodei = element_node_3d.reshape(nplies, 2)

            centroid_datai_5d = case.data.copy().reshape(ntime, nelement, nnode, nlayer, nresult)  # reshape to be easier to slice

            # pull off the centroid
            #[ntime, nelement, nlayer, nresult]
            centroid_datai = centroid_datai_5d[:, :, 0, :, :]

            # pull off the nodes
            #[ntime, nelement, nnode, nlayer, nresult]
            node_datai     = centroid_datai_5d[:, :, 1:, :, :]
            if is_filter:
                node_datai = node_datai[:, ieid_filter, :, :, :]

            node_datai2 = node_datai.reshape(ntime, nelement_filtered*(nnode-1), nlayer, nresult)
            node_data_list.append(node_datai2)
            element_node_list.append(element_nodei)
            del node_datai, node_datai2, element_node_3d, element_node_4d, element_nodei, centroid_datai_5d
        else:
            # ctria3 - no nodal
            centroid_eidsi_og = case.element_node[::2, 0]
            #nelement = len(centroid_eidsi)

            # filter out unused ids
            centroid_eidsi, ieid_filter, nelement_filtered, is_filter = filter_ids(
                all_element_id, centroid_eidsi_og)
            if nelement_filtered == 0:
                continue

            #-------------------
            # this is long, but not that bad...

            # slice off for output
            centroid_datai = case.data.reshape(
                ntime, nelement_nnode_og, nlayer, nresult).copy()

            # setup for nodes
            node_data_5d = case.data.reshape(
                ntime, nelement_nnode_og, 1, nlayer, nresult).copy()

            if is_filter:
            # we reshape the data so we can slice off the desired rows
                centroid_datai = centroid_datai[:, ieid_filter, :, :]
                node_data_5d = node_data_5d[:, ieid_filter, :, :, :]

            eid0 = centroid_eidsi[0]
            nids0 = eid_to_nid_map[eid0]

            # ----------------------------------------------
            # Spoof the values for the tri
            # TODO: probably wrong for fancy CQUAD8/CTRIA6
            nnode = len(nids0)

            ##centroid_dataii (1, 278, 2, 8)

            # prove out stacking...
            #node_data_5d[0, 0, 0, 0, 0] = 1.
            #node_data_5d[0, 0, 0, 1, 0] = 2.
            node_data_3d = node_data_5d.reshape(ntime, nelement_filtered*nlayer, nresult)
            # ------------------------------------------------------------------------
            # good
            node_data_shape = (ntime, nelement_filtered*nnode*nlayer, nresult)
            node_datai = np.full(node_data_shape, np.nan, dtype=node_data_3d.dtype)

            ## gross...interleaving the centroid data to fake nodal data
            # alternate approaches are using np.stack, np.tile, ... I couldn
            # [0, 1, 2]   0
            # [3, 4, 5]   1
            #  +
            # [6, 7, 8]   2
            # [9, 10, 11] 3
            # -------------
            # [0, 1, 2]   0
            # [3, 4, 5]   1
            #
            # [0, 1, 2]   2
            # [3, 4, 5]   3
            #
            # [0, 1, 2]   4
            # [3, 4, 5]   5
            #  +
            # [6, 7, 8]   6
            # [9, 10, 11] 7
            #
            # [6, 7, 8]   8
            # [9, 10, 11] 9
            #
            # [6, 7, 8]   10
            # [9, 10, 11] 11

            # every 3rd layer; should be reasonbly fast
            ilayer0 =  np.arange(nelement_filtered*nlayer)
            ilayer0a = ilayer0[::nlayer]
            ilayer0b = ilayer0a + 1
            ilayer1 = nnode * ilayer0a
            ilayer2 = ilayer1 + 1
            node_datai2 = np.full(node_data_shape, np.nan, dtype=node_data_3d.dtype)
            for inode in range(nnode):
                node_datai2[:, nlayer*inode+ilayer1, :] = node_data_3d[:, ilayer0a, :]
                node_datai2[:, nlayer*inode+ilayer2, :] = node_data_3d[:, ilayer0b, :]
            node_datai = node_datai2.reshape(ntime, nelement_filtered*nnode, nlayer, nresult)

            # not nearly as critical as the previous, but could be sped up
            element_nodei = []
            for eid in centroid_eidsi:
                nids = eid_to_nid_map[eid]
                for nid in nids:
                    # two layers per node
                    element_nodei.append((eid, nid))
                    element_nodei.append((eid, nid))
            element_node_list.append(np.array(element_nodei))
            node_data_list.append(node_datai)
            del eid0, nids, nnode, element_nodei
            del ilayer0, ilayer0a, ilayer0b, ilayer1, ilayer2

        # slice off the centroid
        centroid_elements_list.append(centroid_eidsi)
        centroid_data_list.append(centroid_datai)
        del centroid_datai, centroid_eidsi

    centroid_eids = np.hstack(centroid_elements_list)
    centroid_data = np.hstack(centroid_data_list)

    element_node = np.vstack(element_node_list)
    node_data = np.hstack(node_data_list)
    return centroid_eids, centroid_data, element_node, node_data

def _check_sidebar_args(plate: PlateStrainStressResults2,
                        itime: int, res_name: tuple,
                        transform: str, # Material
                        min_max_method: str, # Absolute Max
                        nodal_combine: str): # Centroid
    transforms = plate.has_coord_transform(itime, res_name)[1]
    min_max_methods = plate.has_derivation_transform(itime, res_name)[1]['derivation']
    combine_methods = plate.has_nodal_combine_transform(itime, res_name)[1]

    transform = transform if transform else transforms[0]
    min_max_method = min_max_method if min_max_method else min_max_methods[0]
    nodal_combine = nodal_combine if nodal_combine else combine_methods[0]

    assert transform in transforms, transform
    assert min_max_method in min_max_methods, min_max_method
    assert nodal_combine in combine_methods, nodal_combine
    return transform, min_max_method, nodal_combine

def get_real_plate_cases(element_id: np.ndarray,
                         model: OP2,
                         key,
                         is_stress: bool,
                         prefix: str='',
                         require_results: bool=True) -> tuple[list, int]:
    """
    Parameters
    ----------
    key : subcase id or tuple[subcase_id, analysis_code, ...]
        subcase_id: source of subcase_id
        analysis_code: unused here, but used somewhere else...
    element_id : np.ndarray
        the elements in the model to consider
        - may be all the elements in the model
        - may be a subset of the plate elements
    prefix: str; default=''
        ''                   : get stress or strain
        'modal_contribution' : get model_contribution.stress or model_contribution.strain
        #'stress', 'model_contribution.stress', 'strain', 'model_contribution.strain'
    is_stress : bool
        flag for selecting stress/strain

    Returns
    -------
    plate_cases : list[PlateStressArray] | list[PlateStrainArray]
        the results to add
    subcase_id : int
        #an echo of key if key is subcase_id
        otherwise, pulls out the subcase id

    """
    plates, word, subcase_id, analysis_code = _get_plates(
        model, key, is_stress, prefix)

    plate_cases = []
    if len(plates) == 0:
        if require_results:
            raise RuntimeError('no results found')
        return plate_cases, subcase_id

    for iplate, result in enumerate(plates):
        # if result:
        # print(f'keys[{iplate}] = {result.keys()}')
        if key not in result:
            continue
        plate_case = result[key]
        if plate_case.is_complex:
            continue
        eids = np.unique(plate_case.element_node[:, 0])
        common_eids = np.intersect1d(element_id, eids)
        if len(common_eids) == 0:
            continue
        # ieids = np.unique(np.searchsorted(element_id, eids))
        # if element_id[ieids[0]] != eids[0]
        plate_cases.append(plate_case)

    if len(plate_cases) == 0  and require_results:
        raise RuntimeError('no results found')
    return plate_cases, subcase_id

def _get_plates(model: OP2,
                key: tuple[int, int],
                is_stress: bool,
                prefix: str) -> tuple[str, list, int, int]:
    if isinstance(key, int):
        subcase_id = key
        analysis_code = -1
    else:
        subcase_id = key[0]
        analysis_code = key[1]
        #print("***stress eids=", eids)

    assert isinstance(prefix, str), f'prefix={prefix}'
    if prefix == 'modal_contribution':
        results = model.op2_results.modal_contribution
        preword = 'Modal Contribution '
    elif prefix == '':
        if is_stress:
            results = model.op2_results.stress
        else:
            results = model.op2_results.strain
        preword = ''
    else:  # pragma: no cover
        raise NotImplementedError(f'prefix={prefix}')

    if is_stress:
        plates = [
            results.ctria3_stress, results.cquad4_stress,
            results.ctria6_stress, results.cquad8_stress,
            results.ctriar_stress, results.cquadr_stress, # results.cquad_stress,
        ]
        word = preword + 'Stress'
    else:
        plates = [
            results.ctria3_strain, results.cquad4_strain,
            results.ctria6_strain, results.cquad8_strain,
            results.ctriar_strain, results.cquadr_strain, # results.cquad_strain,
        ]
        word = preword + 'Strain'

    plates2 = [plate for plate in plates if len(plate)]
    return plates2, word, subcase_id, analysis_code

def plate_cases_to_iresult(plate_cases: list,
                           is_stress: bool) -> tuple[dict[int, tuple[str, str]],
                                                     bool, str, int | str]:
    plate_case = plate_cases[0]
    is_fiber_distance = plate_case.is_fiber_distance
    is_von_mises = plate_case.is_von_mises
    assert isinstance(is_von_mises, bool), is_von_mises
    von_misesi = 7 if is_von_mises else 'von_mises'
    max_sheari = 7 if not is_von_mises else 'max_shear'
    # print(case_headers)
    if is_stress:
        # [fiber_dist, 'oxx', 'oyy', 'txy', 'angle', 'omax', 'omin', ovm]
        iresult_to_title_annotation_map = {
            # iresult: (sidebar_label, annotation)
            0: ('FiberDistance', 'Fiber Distance'),
            1: ('Normal XX', 'XX'),
            2: ('Normal YY', 'YY'),
            3: ('Shear XY', 'XY'),
            4: ('Theta', 'Theta'),
            5: ('Max Principal', 'Max Principal'),
            6: ('Min Principal', 'Min Principal'),
            'abs_principal': ('Abs Principal', 'Abs Principal'),
            von_misesi: ('von Mises', 'von Mises'),  # the magnitude is large
            max_sheari: ('Max Shear', 'Max Shear'),  # the magnitude is large
        }
        word = 'Stress'
    else:
        iresult_to_title_annotation_map = {
            # iresult: (sidebar_label, annotation)
            # 'fiber_curvature' : 'FiberCurvature',
            0: ('FiberDistance', 'Fiber Distance'),
            1: ('Normal XX', 'XX'),
            2: ('Normal YY', 'YY'),
            3: ('Shear XY', 'XY'),
            4: ('Theta', 'Theta'),
            5: ('Max Principal', 'Max Principal'),
            6: ('Min Principal', 'Min Principal'),
            'abs_principal': ('Abs Principal', 'Abs Principal'),
            von_misesi: ('von Mises', 'von Mises'),  # the magnitude is small
            max_sheari: ('Max Shear', 'Max Shear'),  # the magnitude is small
        }
        word = 'Strain'
    return iresult_to_title_annotation_map, is_fiber_distance, word, max_sheari
