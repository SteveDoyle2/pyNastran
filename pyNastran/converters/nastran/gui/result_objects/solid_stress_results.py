from __future__ import annotations
#from collections import defaultdict
import numpy as np
from typing import Optional, Any, TYPE_CHECKING

from pyNastran.utils.mathematics import get_abs_max
from pyNastran.op2.stress_reduction import max_shear, von_mises_3d

#from pyNastran.femutils.utils import abs_nan_min_max # , pivot_table,  # abs_min_max
#from pyNastran.bdf.utils import write_patran_syntax_dict

from .nodal_averaging import nodal_average, nodal_combine_map
from pyNastran.gui.gui_objects.vector_results import VectorResultsCommon, filter_ids
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf import BDF
    from pyNastran.op2.op2 import OP2
    from pyNastran.op2.tables.oes_stressStrain.real.oes_solids import RealSolidArray
    CaseTuple = tuple[int, int, str]


class SolidStrainStressResults2(VectorResultsCommon):
    def __init__(self,
                 subcase_id: int,
                 model: BDF,
                 node_id: np.ndarray,
                 element_id: np.ndarray,
                 cases: list[RealSolidArray],
                 result_dict: dict[int | str, Any],
                 title: str,
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
        model : BDF
            the geometry model
        node_id: (nnode,) int np.ndarray
            all the nodes ids in the model
        element_id: (nelement,) int np.ndarray
            all the element ids in the model
        cases: list[RealSolidArray]
            all the solid stress/strain cases for a given subcase
        result_dict: dict[int | str, Any]
            mapping of the result types
        title: str
            the legend title
        data_format: str='%g'
            ???
        is_variable_data_format: bool=False
            ???
        headers : list[str]
            the sidebar word
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
        title = ''
        ntitles = 10
        VectorResultsCommon.__init__(
            self, subcase_id, title,
            cases, ntitles,
            data_format=data_format,
            is_variable_data_format=is_variable_data_format,
            nlabels=nlabels, labelsize=labelsize, ncolors=ncolors,
            colormap=colormap,
            #set_max_min: bool=False,
            uname=uname)

        self.centroid_data = np.zeros((0, 0, 0), dtype='float32')
        self.node_data = np.zeros((0, 0, 0), dtype='float32')
        self.element_node = np.zeros((0, 2), dtype='int32')
        self.inode = np.zeros(0, dtype='int32')
        self.ielement_centroid = np.zeros(0, dtype='int32')

        self.layer_indices = (-1, )  # All
        i = -1
        name = ''

        # slice off the methods (from the boolean) and then pull the 0th one
        self.min_max_method = '' # self.has_derivation_transform(i, name)[1]['derivation'][0]
        self.transform = self.has_coord_transform(i, name)[1][0]
        self.nodal_combine = self.has_nodal_combine_transform(i, name)[1][0]

        self.dim = cases[0].data.ndim
        for case in cases:
            assert case.data.ndim == 3, case.data.shape

        self.is_stress = case.is_stress

        self.layer_map = {
            0: 'Centroid',  # default
            1: 'Corner',
        }

        #  global ids
        self.model = model
        self.node_id = node_id
        self.element_id = element_id

        # local case object
        #self.cases = cases
        self.result_dict = result_dict

        self.data_type = case.data.dtype.str # '<c8', '<f4'
        self.is_real = True if self.data_type in ['<f4', '<f8'] else False
        #self.is_real = dxyz.data.dtype.name in {'float32', 'float64'}
        self.is_complex = not self.is_real

        #ntimes = case.data.shape[0]
        # ------------------------------------------------------------------
        #linked_scale_factor = False
        #location = 'node'

        out = setup_centroid_node_data(cases, element_id)
        centroid_eids, centroid_data, element_node, node_data = out

        self.centroid_eids = centroid_eids
        # [ntimes, nelements, nresults]
        self.centroid_data = centroid_data

        self.element_node = element_node
        self.node_data = node_data

        #common_eids = np.intersect1d(self.centroid_eids, element_id)
        #if len(common_eids) == 0:
            #raise IndexError('no solid elements found...')
        #elif len(common_eids) != len(self.centroid_eids):
            #icommon = np.searchsorted(common_eids, self.centroid_eids)
            ##self.centroid_data = self.centroid_data[:, icommon, :]
            #raise RuntimeError('some common elements were found...but some are missing')

        nids = np.unique(self.element_node[:, 1])
        self.inode = np.searchsorted(node_id, nids)
        self.ielement_centroid = np.searchsorted(element_id, self.centroid_eids)

        assert np.array_equal(element_id[self.ielement_centroid], self.centroid_eids)
        # dense -> no missing nodes in the results set
        self.is_dense = (
            (len(element_id) == len(self.centroid_eids)) and
            (len(node_id) == len(nids)))
        #self.is_dense = False

        #self.xyz = xyz
        #assert len(self.xyz.shape) == 2, self.xyz.shape
        #self.location = 'centroid'
        ntimes = self.centroid_data.shape[0]
        if self.is_stress:
            self.headers = ['SolidStress2'] * ntimes
        else:
            self.headers = ['SolidStrain2'] * ntimes
        str(self)

    @classmethod
    def add_stress(cls,
                   subcase_id: int,
                   model: BDF,
                   model_results: OP2,
                   element_id: np.ndarray,):
                   #is_variable_data_format: bool=False,
                   #prefix='stress'):
        out = cls.load_from_code(
            subcase_id, model, model_results, element_id,
            is_stress=True,
            require_results=True)
        return out

    @classmethod
    def add_strain(cls,
                   subcase_id: int,
                   model: BDF,
                   model_results: OP2,
                   element_id: np.ndarray,):
                   #is_variable_data_format: bool=False,
                   #prefix='stress'):
        out = cls.load_from_code(
            subcase_id, model, model_results, element_id,
            is_stress=False,
            require_results=True)
        return out

    @classmethod
    def load_from_code(cls,
                       subcase_id: int,
                       model: BDF,
                       model_results: OP2,
                       element_id: np.ndarray,
                       is_stress: bool,
                       #is_variable_data_format: bool=False,
                       #prefix='stress',
                       require_results: bool=True):
        """the intention of this method is to not use the gui"""
        assert is_stress, is_stress
        solid_cases, subcase_id = get_real_solid_cases(
            element_id,
            model_results,
            subcase_id,
            is_stress=is_stress,
            prefix='',
            require_results=require_results)

        # get the xyzs
        out = model.get_xyz_in_coord_array(
            cid=0, fdtype='float64', idtype='int64')
        nid_cp_cd, xyz_cid, unused_xyz_cp, unused_icd_transform, unused_icp_transform = out
        node_id = nid_cp_cd[:, 0]
        #------------------------------------------
        out = solid_cases_to_iresult(solid_cases, is_stress)
        iresult_to_title_annotation_map, word, max_sheari = out

        #------------------------------------------
        title = 'title'
        obj = SolidStrainStressResults2(
            subcase_id,
            model,
            node_id,
            element_id,
            solid_cases,
            iresult_to_title_annotation_map, title,
            is_variable_data_format=False,  # ???
            set_max_min=False)
        # data_format: str = '%g',
        # nlabels = None, labelsize = None, ncolors = None,
        # colormap: str = '',
        # uname: str = 'SolidStressStrainResults2'):

        return obj

    def get_methods(self, itime: int, res_name: str) -> list[str]:
        layers = list(self.layer_map.values())
        return layers

    def _get_default_tuple_indices(self) -> tuple[int, ...]:
        out = tuple(np.array(self._get_default_layer_indicies()) - 1)
        return out

    def _get_default_layer_indicies(self) -> tuple[int, ...]:
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
        self._method_keys_to_layer_indices(methods_keys)

        # doesn't matter cause it's already nodal
        #assert min_max_method in min_max_methods, min_max_method
        assert nodal_combine in combine_methods, nodal_combine
        #self.min_max_method = min_max_method
        self.nodal_combine = nodal_combine
        self.transform = transform

    def _method_keys_to_layer_indices(self,
                                      methods_keys: Optional[list[int]]=None,
                                      ) -> tuple[int, ...]:
        default_indices = self._get_default_layer_indicies()
        if methods_keys is None or len(methods_keys) == 0:
            # default; Centroid
            indices = default_indices
        elif 1 in methods_keys: # Centroid
            # Corner takes precendence (it has more data)
            indices = (1, )
        else:
            # Centroid is secondary (even though it's first)
            indices = (0, )
        self.layer_indices = indices
        #self.layer_indices = (1, )
    def set_to_centroid(self):
        self.layer_indices = (0, )
    def set_to_node(self, nodal_combine: str):
        self.layer_indices = (1, )
        self._set_nodal_combine_transform(nodal_combine)

    def has_methods_table(self, i: int, res_name: str) -> bool:
        return True
    def has_coord_transform(self, i: int, res_name: str) -> tuple[bool, list[str]]:
        return True, ['Material']
    def has_derivation_transform(self, i: int, case_tuple: CaseTuple,
                                 ) -> tuple[bool, dict[str, Any]]:
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
    def _set_nodal_combine_transform(self, nodal_combine: str) -> None:
        itime = -1
        res_name = ''
        junk, allowed = self.has_nodal_combine_transform(itime, res_name)
        assert nodal_combine in allowed, f'nodal_combine={nodal_combine!r} is not allowed; allowed={allowed}'
        self.nodal_combine = nodal_combine
    def has_nodal_combine_transform(self, i: int, res_name: str) -> tuple[bool, list[str]]:
        """elemental -> nodal (only applies for Corner)"""
        return True, ['Absolute Max', 'Mean', 'Max', 'Min',
                      'Difference', 'Std. Dev.',]
    # 'Nodal Max'
    def has_output_checks(self, i: int, resname: str) -> tuple[bool, bool, bool,
                                                               bool, bool, bool]:
        is_enabled_fringe = is_checked_fringe = True
        is_enabled_disp = is_checked_disp = is_enabled_vector = is_checked_vector = False
        out = (is_enabled_fringe, is_checked_fringe,
               is_enabled_disp, is_checked_disp,
               is_enabled_vector, is_checked_vector)
        return out

    def get_annotation(self, itime: int, case_tuple: CaseTuple) -> str:
        """
        A header is the thingy that goes in the lower left corner
        title = 'Solid Stress'
        method = 'Absolute Max'
        header = 'Static'
        nodal_combine = 'Nodal Max'
        returns 'Solid Stress (Center, Absolute Max; Max, Static): sigma11'

        returns 'Solid Stress (Centroid, Static): sigma11'
        """
        # overwrite itime based on linked_scale factor
        (itime, iresult, header) = case_tuple
        itime, unused_case_flag = self.get_case_flag(itime, case_tuple)

        default_indices = self._get_default_tuple_indices() # 0-based
        if self.layer_indices == default_indices:
            self.layer_indices = (0, )

        if self.layer_indices == (0, ):
            layer_str = self.layer_map[0]  # Centroid
        elif self.layer_indices == (1, ):
            layer_str = self.layer_map[1]  # Corner
        else:  # pragma: no cover
            raise RuntimeError(self.layer_indices)
        self.layer_indices

        result = get_solid_result(self.result_dict, iresult, index=1)

        title = self.get_legend_title(itime, case_tuple)
        if layer_str == 'Centroid':
            #'Solid Stress; Centroid (Static): von Mises'
            annotation_label = f'Solid {title} ({layer_str}, {header}): {result}'
        else:
            #'Solid Stress; Corner (Mean; Static): von Mises'
            annotation_label = f'Solid {title} ({layer_str}, {self.nodal_combine}, {header}): {result}'
        #return self.uname
        return annotation_label

    def get_case_flag(self, i: int,
                      case_tuple: CaseTuple) -> tuple[int,
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

    #def get_default_legend_title(self, itime: int, case_tuple: CaseTuple) -> str:
        #(itime, iresult, header) = case_tuple
        ##method_ = 'Composite Stress Layers:' if self.is_stress else 'Composite Strain Layers:'
        ##self.layer_indices
        #results = list(self.result.values())
        ##method = method_ + ', '.join(str(idx) for idx in (self.layer_indices+1))
        ##method = method.strip()
        ##title = f'{self.title} {method}'
        #title = results[iresult][0]  # sidebar label=legend
        #return title
    def get_legend_tuple(self, itime: int, case_tuple: CaseTuple) -> int:
        (itime, iresult, header) = case_tuple
        return iresult
    def get_default_legend_title(self, itime: int, case_tuple: CaseTuple) -> str:
        """Composite Stress Layers: 1, 2, 3, 4"""
        (itime, iresult, header) = case_tuple
        #method_ = 'Composite Stress Layers:' if self.is_stress else 'Composite Strain Layers:'
        #self.layer_indices
        #self.result
        #method = method_ + ', '.join(str(idx) for idx in (self.layer_indices+1))
        #title = f'{self.title} {method}'
        result = get_solid_result(self.result_dict, iresult, index=0)
        return result

    def _get_real_data(self, case_tuple: CaseTuple) -> np.ndarray:
        (itime, iresult, header) = case_tuple

        #ilayer = self.layer_indices
        if self.layer_indices == (-1, ):
            self.layer_indices = (0, )

        #self.case.get_headers()
        #['oxx', 'oyy', 'ozz', 'txy', 'tyz', 'txz', 'omax', 'omid', 'omin', von_mises]
        if self.layer_indices == (0, ):  # Centroid
            data = self._get_centroid_result(itime, iresult)
        elif self.layer_indices == (1, ):  # Corner
            data = self._get_nodal_result(itime, iresult)
        #elif self.nodal_combine == 'Nodal Max':
        #elif self.nodal_combine == 'Nodal Min':
        #elif self.nodal_combine == 'Nodal Mean':
        #elif self.nodal_combine == 'Nodal Abs Max':
        #elif self.nodal_combine == 'Nodal Std. Dev.':
        #elif self.nodal_combine == 'Nodal Difference':
        else:  # pragma: no cover
            raise RuntimeError(self.nodal_combine)

        assert len(data.shape) == 1, data.shape
        return data

    def get_case_tuple(self, itime: int, result_name: str,
                       ) -> tuple[int, int | str, str]:
        assert itime <= self.centroid_data.shape[0]
        result_name = result_name.lower()
        if result_name in self.result_dict:
            iresult = result_name

        # result_dict = {
        #     0: ('Normal XX', 'XX'),
        #     1: ('Normal YY', 'YY'),
        #     2: ('Normal ZZ', 'ZZ'),
        #     3: ('Shear XY', 'XY'),
        #     4: ('Shear YZ', 'YZ'),
        #     5: ('Shear XZ', 'XZ'),
        #     6: ('Max Principal', 'Max Principal'),
        #     8: ('Min Principal', 'Min Principal'),
        #     7: ('Mid Principal', 'Mid Principal'),
        #     9: ('von Mises', 'von Mises'),
        #     'max_shear': ('Max Shear', 'Max Shear'),
        # }
        elif result_name in (0, 'xx', 'normal xx'):
            iresult = 0
        elif result_name in (1, 'yy', 'normal yy'):
            iresult = 1
        elif result_name in (2, 'zz', 'normal zz'):
            iresult = 2
        elif result_name in (3, 'xy', 'shear xy'):
            iresult = 3
        elif result_name in (4, 'yz', 'shear xy'):
            iresult = 4
        elif result_name in (5, 'xz', 'shear xz'):
            iresult = 5
            #     6: ('Max Principal', 'Max Principal'),
            #     8: ('Min Principal', 'Min Principal'),
            #     7: ('Mid Principal', 'Mid Principal'),
            #     9: ('von Mises', 'von Mises'),
            #     'max_shear': ('Max Shear', 'Max Shear'),
        elif result_name in {'von_mises', 'max_shear'}:
            iresult = result_name if result_name in self.result_dict else 9
            # is_von_mises = solid_case.is_von_mises
            # assert isinstance(is_von_mises, bool), is_von_mises
            # von_misesi = 9 if is_von_mises else 'von_mises'
            # max_sheari = 9 if not is_von_mises else 'max_shear'
            # iresult = result_name
        else:  # pragma: no cover
            raise RuntimeError(result_name)
        return (itime, iresult, '')

    def _get_centroid_result(self, itime: int,
                             iresult: int | str) -> np.ndarray:
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
        centroid_data = self.centroid_data
        if iresult == 'abs_principal': # abs max; should be good
            omax = centroid_data[itime, :, imax]
            omin = centroid_data[itime, :, imin]
            data = get_abs_max(omin, omax, dtype=omin.dtype)

        elif iresult == 'von_mises':
            # von mises
            #Derive
            oxx = centroid_data[itime, :, ioxx]
            oyy = centroid_data[itime, :, ioyy]
            ozz = centroid_data[itime, :, iozz]
            txy = centroid_data[itime, :, itxy]
            txz = centroid_data[itime, :, itxz]
            tyz = centroid_data[itime, :, ityz]
            data = von_mises_3d(oxx, oyy, ozz, txy, tyz, txz, self.is_stress)
        elif iresult == 'max_shear': #  probably wrong
            # not checked for strain
            omax = centroid_data[itime, :, imax]
            omin = centroid_data[itime, :, imin]
            data = max_shear(omax, omin)
        else:
            data = centroid_data[itime, :, iresult].copy()
        return data

    def _get_nodal_result(self, itime: int,
                          iresult: int | str) -> np.ndarray:
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
        nodal_combine_func = nodal_combine_map[self.nodal_combine]

        element_node = self.element_node
        nids = np.unique(element_node[:, 1])
        #ids = np.arange(len(unode_id))

        #inid = self.inode   # this is the global mapper
        #inid = np.searchsorted(unode_id, nids)
        #nid_to_inid_map = {nid: inidi for nid, inidi in zip(nids, inid)}
        nid_to_inid_map = {nid: i for i, nid in enumerate(nids)}

        ## Derive/Average
        node_data = self.node_data
        if isinstance(iresult, int):
            datai = node_data[itime, :, iresult]
        elif iresult == 'von_mises':
            oxx = node_data[itime, :, ioxx]
            oyy = node_data[itime, :, ioyy]
            ozz = node_data[itime, :, iozz]
            txy = node_data[itime, :, itxy]
            txz = node_data[itime, :, itxz]
            tyz = node_data[itime, :, ityz]
            datai = von_mises_3d(oxx, oyy, ozz, txy, tyz, txz, self.is_stress)
        elif iresult == 'max_shear':
            omax = node_data[itime, :, imax]
            omin = node_data[itime, :, imin]
            datai = max_shear(omax, omin)
        elif iresult == 'abs_principal':
            data_max = node_data[itime, :, imax]
            data_min = node_data[itime, :, imin]
            datai = get_abs_max(data_min, data_max, dtype=data_min.dtype)
            assert datai.shape == data_min.shape
        else:  # pragma: no cover
            raise NotImplementedError(iresult)

        data = nodal_average(
            nodal_combine_func,
            element_node, datai,
            nids, nid_to_inid_map)
        return data

    #def _get_complex_data(self, itime: int) -> np.ndarray:
        #return self._get_real_data(itime)

    def _get_fringe_data_sparse(self, itime: int,
                                case_tuple: CaseTuple) -> np.ndarray:
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

    def _get_fringe_data_dense(self, itime: int,
                               case_tuple: CaseTuple) -> np.ndarray:
        data = self._get_fringe_data_sparse(itime, case_tuple)
        if self.is_dense:
            return data

        if self.get_location(0, 0) == 'node':
            nnode = len(self.node_id)
            result_out = np.full(nnode, np.nan, dtype=data.dtype)
            result_out[self.inode] = data
        else:
            nelements = len(self.element_id)
            result_out = np.full(nelements, np.nan, dtype=data.dtype)
            result_out[self.ielement_centroid] = data
        return result_out

    def get_location_arrays(self) -> tuple[np.ndarray, np.ndarray]:
        if self.layer_indices == (0, ):  # Centroid
            all_ids = self.element_id
            ids = np.unique(self.element_node[:, 0])
        elif self.layer_indices == (1, ):  # Corner
            all_ids = self.node_id
            ids = np.unique(self.element_node[:, 1])
        return all_ids, ids

    def get_fringe_result(self, itime: int,
                          case_tuple: CaseTuple) -> np.ndarray:
        """
        gets the 'typical' result which is a vector
         - GuiResult:           fringe; (n,)   array
         - DisplacementResults: vector; (n, 3) array
        """
        fringe = self._get_fringe_data_dense(itime, case_tuple)
        return fringe

    def get_fringe_vector_result(self, itime: int,
                                 case_tuple: CaseTuple) -> tuple[np.ndarray, None]:
        """
        gets the 'typical' result which is a vector
         - GuiResult:           fringe; (n,)   array
         - DisplacementResults: vector; (n, 3) array
        """
        fringe = self.get_fringe_result(itime, case_tuple)
        return fringe, None

    def get_default_scale(self, itime: int, res_name: str) -> float:
        return 0.0
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

    def get_location(self, unused_i: int, unused_res_name: str) -> str:
        """the result type"""
        if self.layer_indices == (0, ):
            return 'centroid'
        return 'node'
        #return self.location

    #-------------------------------------

    #def get_force_vector_result(self, itime: int, res_name: str) -> np.ndarray:
        #dxyz = self._get_fringe_data_dense(itime, res_name)
        #scale = 1.
        #return self.xyz, dxyz * scale

    #def get_vector_result(self, itime: int, res_name: str) -> tuple[np.ndarray, np.ndarray]:
        #dxyz = self._get_fringe_data_dense(itime, res_name)
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


def get_solid_result(result_dict: dict[int | str, Any],
                     iresult: int | str, index: int) -> str:
    """
    values
    0=title, 'annotation'
    ('sAbs Principal', 'Abs Principal')
    """
    assert index in (0, 1), index
    results = result_dict[iresult]
    return results[index]

def setup_centroid_node_data(cases: list[RealSolidArray],
                             element_id: np.ndarray) -> tuple[np.ndarray, np.ndarray,
                                                              np.ndarray, np.ndarray]:
    # setup the node mapping
    centroid_data_list = []
    centroid_elements_list = []

    node_data_list = []
    element_node_list = []
    for case in cases:
        ntime, nelement_nnode, nresults = case.data.shape
        ueids = np.unique(case.element_node[:, 0])
        nelementi = len(ueids)
        nnode = case.nnodes_per_element
        if case.nnodes_per_element > 1:
            eids_og = case.element_node[0::nnode, 0]
            nelement = len(case.element_node) // nnode
            nnodal_node = nelement * (nnode - 1)

            eids, ifilter, nelement_filtered, is_filter = filter_ids(element_id, eids_og)
            if nelement_filtered == 0:
                continue
            nnodal_node_filtered = nelement_filtered * (nnode - 1)
            element_node_3d = case.element_node.reshape(nelement, nnode, 2)
            data = case.data.copy().reshape(ntime, nelement, nnode, nresults)
            if is_filter:
                element_node_3d = element_node_3d[ifilter, :, :]
                data = data[:, ifilter, :, :]
            element_node = element_node_3d[:, 1:, :].reshape(nnodal_node_filtered, 2)
            eids = element_node_3d[:, 0, 0]

            centroid_elements_list.append(eids)
            #nodal_elements.append(case.element_node[1::nnodes_per_element, 0])
            assert nelement == nelementi, (nelement, nelementi)
            node_data = data[:, :, 1:, :]
            node_data2 = node_data.reshape(ntime, nnodal_node_filtered, nresults)
        else:
            eids_og = case.element_node[:, 0]
            eids, ifilter, nelement_filtered, is_filter = filter_ids(element_id, eids_og)
            if nelement_filtered == 0:
                continue
            centroid_elements_list.append(case.element_node)
            data = case.data.reshape(ntime, nelementi, 1, nresults)
            node_data = data
            raise NotImplementedError('centroidal solid elements')
        # slice off the centroid
        centroid_datai = data[:, :, 0, :]
        centroid_data_list.append(centroid_datai)

        node_data_list.append(node_data2)
        element_node_list.append(element_node)

    centroid_eids = np.hstack(centroid_elements_list)
    # [ntimes, nelements, nresults]
    centroid_data = np.hstack(centroid_data_list)

    element_node = np.vstack(element_node_list)
    node_data = np.hstack(node_data_list)

    #intersect_eids, ieid_filter, nelement_filtered, is_filter = filter_ids(element_id, centroid_eids)
    return centroid_eids, centroid_data, element_node, node_data

def solid_cases_to_iresult(solid_cases: list,
                           is_stress: bool) -> tuple[dict[int, tuple[str, str]],
                                                     str, int | str]:
    solid_case = solid_cases[0]
    #solid_case_headers = solid_case.get_headers()
    is_von_mises = solid_case.is_von_mises
    assert isinstance(is_von_mises, bool), is_von_mises
    von_misesi = 9 if is_von_mises else 'von_mises'
    max_sheari = 9 if not is_von_mises else 'max_shear'

    if is_stress:
        #['oxx', 'oyy', 'ozz', 'txy', 'tyz', 'txz', 'omax', 'omid', 'omin', von_mises]
        iresult_to_title_annotation_map = {
            # iresult: (sidebar_label, annotation)
            0: ('Normal XX', 'XX'),
            1: ('Normal YY', 'YY'),
            2: ('Normal ZZ', 'ZZ'),
            3: ('Shear XY', 'XY'),
            4: ('Shear YZ', 'YZ'),
            5: ('Shear XZ', 'XZ'),

            6: ('Max Principal', 'Max Principal'),
            8: ('Min Principal', 'Min Principal'),
            7: ('Mid Principal', 'Mid Principal'),
            #'abs_principal': ('sAbs Principal', 'Abs Principal'),
            von_misesi: ('von Mises', 'von Mises'), # the magnitude is large
            max_sheari: ('Max Shear', 'Max Shear'), # the magnitude is large
        }
        word = 'Stress'
    else:
        iresult_to_title_annotation_map = {
            # iresult: (sidebar_label, annotation)
            0: ('Normal XX', 'XX'),
            1: ('Normal YY', 'YY'),
            2: ('Normal ZZ', 'ZZ'),
            3: ('Shear XY', 'XY'),
            4: ('Shear YZ', 'YZ'),
            5: ('Shear XZ', 'XZ'),

            6: ('Max Principal', 'Max Principal'),
            8: ('Min Principal', 'Min Principal'),
            7: ('Mid Principal', 'Mid Principal'),
            von_misesi: ('von Mises', 'von Mises'), # the magnitude is small
            max_sheari: ('Max Shear', 'Max Shear'), # the magnitude is small
        }
        word = 'Strain'
    return iresult_to_title_annotation_map, word, max_sheari


def get_real_solid_cases(element_id: np.ndarray,
                         model: OP2,
                         key: int,
                         is_stress: bool=True,
                         prefix: str='',
                         require_results: bool=True) -> tuple[list, int]:
    assert isinstance(require_results, bool), require_results
    solids, word, subcase_id, analysis_code = _get_solids(
        model, key, is_stress, prefix)
    #assert len(solids) > 0, solids

    solid_cases = []
    for solid_case in solids:
        if solid_case.is_complex:
            continue
        solid_eids = np.unique(solid_case.element_node[:, 0])
        common_eids = np.intersect1d(element_id, solid_eids)
        if len(common_eids) == 0:
            continue
        solid_cases.append(solid_case)

    #assert len(solid_cases) > 0, solid_cases
    return solid_cases, subcase_id

def _get_solids(results: OP2,
                key,
                is_stress: bool,
                prefix: str) -> tuple[str, list, int, int]:
    if isinstance(key, tuple):
        analysis_code = key[1]
        #print("***stress eids=", eids)
        subcase_id = key[0]
    else:
        subcase_id = key
        analysis_code = -1
    #if prefix == 'modal_contribution':
        #results = model.op2_results.modal_contribution
        #preword = 'Modal Contribution '
    #elif prefix == '':
        #results = model
        #preword = ''
    #else:  # pragma: no cover
        #raise NotImplementedError(prefix)

    if is_stress:
        stress = results.op2_results.stress
        cards = [
            stress.ctetra_stress, stress.cpenta_stress, stress.chexa_stress, # stress.cpyram_stress,
        ]
        word = 'Stress'
    else:
        strain = results.op2_results.strain
        cards = [
            strain.ctetra_strain, strain.cpenta_strain, strain.chexa_strain, # strain.cpyram_strain,
        ]
        word = 'Strain'
    cards = [card for card in cards if card]
    #assert len(cards) > 0, cards

    cards2 = []
    for result in cards:
        if key not in result:
            continue
        cards2.append(result[key])
    #assert len(cards2) > 0, cards2
    return cards2, word, subcase_id, analysis_code
