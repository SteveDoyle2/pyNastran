# pylint: disable=C0301,C0103,R0913,R0914,R0904,C0111,R0201,R0902
from itertools import count
from struct import pack
import warnings
from typing import TextIO, Optional

import numpy as np
from numpy.linalg import eigh  # type: ignore

from pyNastran.utils.numpy_utils import float_types, integer_float_types, type_integer_float
from pyNastran.utils.mathematics import get_abs_max

from pyNastran.f06.f06_formatting import (
    write_floats_13e, write_floats_13e_long, _eigenvalue_header)
from pyNastran.op2.result_objects.op2_objects import (
    get_times_dtype, combination_inplace)
from pyNastran.op2.result_objects.utils_pandas import (
    build_dataframe_transient_header, build_pandas_transient_element_node)
from pyNastran.op2.tables.oes_stressStrain.real.oes_objects import (
    StressObject, StrainObject, OES_Object,
    oes_real_data_code, set_static_case, set_modal_case,
    set_transient_case, set_post_buckling_case)
from pyNastran.op2.op2_interface.write_utils import to_column_bytes
from pyNastran.op2.stress_reduction import (
    ovm_shear_3d, principal_components_3d,
    von_mises_3d, max_shear_3d)


ELEMENT_NAME_TO_ELEMENT_TYPE = {
    'CTETRA' : 39,
    'CHEXA' : 67,
    'CPENTA' : 68,
    'CPYRAM' : 255,
}
# SOLID_PRINCIPAL_METHOD = '132'
SOLID_PRINCIPAL_METHOD = 'og'


class RealSolidArray(OES_Object):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        OES_Object.__init__(self, data_code, isubcase, apply_data_code=False)
        #self.code = [self.format_code, self.sort_code, self.s_code]

        #self.ntimes = 0  # or frequency/mode
        #self.ntotal = 0
        self.nelements = 0  # result specific
        self.element_cid = np.zeros((0, 2), dtype='int32')
        self.element_node = np.zeros((0, 2), dtype='int32')

        # if is_sort1:
        #     #sort1
        #     self.add_node = self.add_node_sort1
        #     self.add_eid = self.add_eid_sort1
        # else:
        #     raise NotImplementedError('SORT2')

    @property
    def is_real(self) -> bool:
        return True

    @property
    def is_complex(self) -> bool:
        return False

    def _reset_indices(self) -> None:
        self.itotal = 0
        self.ielement = 0

    def update_data_components(self):
        ntimes, nelements_nnodes = self.data.shape[:2]
        # vm
        oxx = self.data[:, :, 0].reshape(ntimes * nelements_nnodes)
        oyy = self.data[:, :, 1].reshape(ntimes * nelements_nnodes)
        ozz = self.data[:, :, 2].reshape(ntimes * nelements_nnodes)
        txy = self.data[:, :, 3].reshape(ntimes * nelements_nnodes)
        tyz = self.data[:, :, 4].reshape(ntimes * nelements_nnodes)
        txz = self.data[:, :, 5].reshape(ntimes * nelements_nnodes)

        #I1 = oxx + oyy + ozz
        #txyz = txy**2 + tyz**2 + txz ** 2
        #I2 = oxx * oyy + oyy * ozz + ozz * oxx - txyz
        #I3 = oxx * oyy * ozz + 2 * txy * tyz * txz + oxx * tyz**2 - oyy * txz**2 - ozz * txy

        # (n_subarrays, nrows, ncols)
        omax, omid, omin = principal_components_3d(
            ntimes, nelements_nnodes,
            oxx, oyy, ozz, txy, tyz, txz,
            self.is_stress)

        if self.is_strain:
            import warnings
            warnings.warn('verify solid principal stress/strains; solid von-mises stress/strain')

        ovm_sheari = ovm_shear_3d(oxx, oyy, ozz, txy, tyz, txz, omax, omin,
                                  self.is_von_mises, self.is_stress)
        ovm_sheari2 = ovm_sheari.reshape(ntimes, nelements_nnodes)

        self.data[:, :, 6] = omax.reshape(ntimes, nelements_nnodes)
        self.data[:, :, 7] = omid.reshape(ntimes, nelements_nnodes)
        self.data[:, :, 8] = omin.reshape(ntimes, nelements_nnodes)
        self.data[:, :, 9] = ovm_sheari2

        #A = [[doxx, dtxy, dtxz],
             #[dtxy, doyy, dtyz],
             #[dtxz, dtyz, dozz]]
        #(_lambda, v) = eigh(A)  # a hermitian matrix is a symmetric-real matrix

    def __iadd__(self, factor):
        """[A] += b"""
        #[oxx, oyy, ozz, txy, tyz, txz, o1, o2, o3, ovmShear]
        if isinstance(factor, float_types):
            self.data[:, :, :6] += factor
        else:
            # TODO: should support arrays
            raise TypeError(f'factor={factor} and must be a float')
        self.update_data_components()

    def __isub__(self, factor):
        """[A] -= b"""
        if isinstance(factor, float_types):
            self.data[:, :, :6] -= factor
        else:
            # TODO: should support arrays
            raise TypeError(f'factor={factor} and must be a float')
        self.update_data_components()

    def __imul__(self, factor: type_integer_float):
        """[A] *= b"""
        assert isinstance(factor, integer_float_types), f'factor={factor} and must be a float'
        self.data[:, :, :6] *= factor
        self.update_data_components()

    def linear_combination(self, factor: type_integer_float,
                           data: Optional[np.ndarray]=None,
                           update: bool=True):
        """[A] * b"""
        combination_inplace(self.data, data, factor, ires=slice(None, 6))
        # if data is None:
        #     self.data[:, :, :6] *= factor
        # else:
        #     self.data[:, :, :6] += data[:, :, :6] * factor
        if update:
            self.update_data_components()

    # def __rmul__(self, factor: integer_float_types):
    #     """[A] * b"""
    #     assert isinstance(factor, integer_float_types), f'factor={factor} and must be a float'
    #     self.data[:, :, :6] *= factor
    #     self.update_data_components()

    def __idiv__(self, factor: type_integer_float):
        """[A] *= b"""
        assert isinstance(factor, float_types), f'factor={factor} and must be a float'
        self.data[:, :, :6] *= 1. / factor
        self.update_data_components()

    def slice_data(self, slice_elements: np.ndarray) -> int:
        assert slice_elements is not None, slice_elements
        eids = self.element_cid[:, 0]
        common_eids = np.intersect1d(eids, slice_elements)
        icommon = np.searchsorted(eids, common_eids)
        ntotal = len(icommon)
        ncommon = ntotal

        neids = self.element_cid.shape[0]
        nnodes_per_element = self.nnodes_per_element
        self.element_cid = self.element_cid[icommon, :]
        element_node = self.element_node.reshape(neids, nnodes_per_element, 2)
        self.element_node = element_node[icommon, :].reshape(ncommon*nnodes_per_element, 2)

        ntime, neids_nnodes_per_element, nresults = self.data.shape
        data = self.data.reshape(ntime, neids, nnodes_per_element, nresults)
        data_slice = data[:, icommon, :, :]
        self.data = data_slice.reshape(ntime, ncommon * nnodes_per_element, nresults)
        self.ntotal = ntotal
        return ntotal


    def build(self):
        """sizes the vectorized attributes of the RealSolidArray"""
        #print('ntimes=%s nelements=%s ntotal=%s' % (self.ntimes, self.nelements, self.ntotal))
        assert self.ntimes > 0, 'ntimes=%s' % self.ntimes
        assert self.nelements > 0, 'nelements=%s' % self.nelements
        assert self.ntotal > 0, 'ntotal=%s' % self.ntotal
        #self.names = []
        self.nelements //= self.ntimes
        self.itime = 0
        self.ielement = 0
        self.itotal = 0
        #self.ntimes = 0
        #self.nelements = 0

        dtype, idtype, fdtype = get_times_dtype(self.nonlinear_factor, self.size, self.analysis_fmt)

        if self.is_sort1:
            ntimes = self.ntimes
            ntotal = self.ntotal
            nelements = self.nelements
        else:
            #print(f'ntimes={self.ntimes} nelements={self.nelements} ntotal={self.ntotal}')
            ntimes = self.nelements
            ntotal = self.ntotal
            nelements = self.ntimes
            #print(f'ntimes={ntimes} nelements={nelements} ntotal={ntotal}')
        #self.ntimes = ntimes
        #self.ntotal = ntotal
        #self.nelements = nelements

        _times = np.zeros(ntimes, dtype=self.analysis_fmt)

        # TODO: could be more efficient by using nelements for cid
        element_node = np.zeros((ntotal, 2), dtype=idtype)
        element_cid = np.zeros((nelements, 2), dtype=idtype)
        #if nelements > 5000:
            #raise RuntimeError(nelements)

        #if self.element_name == 'CTETRA':
            #nnodes = 4
        #elif self.element_name == 'CPENTA':
            #nnodes = 6
        #elif self.element_name == 'CHEXA':
            #nnodes = 8
        #self.element_node = np.zeros((self.ntotal, nnodes, 2), 'int32')

        #[oxx, oyy, ozz, txy, tyz, txz, o1, o2, o3, ovmShear]
        data = np.zeros((ntimes, ntotal, 10), fdtype)
        self.nnodes = element_node.shape[0] // self.nelements
        #self.data = np.zeros((self.ntimes, self.nelements, nnodes+1, 10), 'float32')

        if self.load_as_h5:
            #for key, value in sorted(self.data_code.items()):
                #print(key, value)
            group = self._get_result_group()
            self._times = group.create_dataset('_times', data=_times)
            self.element_node = group.create_dataset('element_node', data=element_node)
            self.element_cid = group.create_dataset('element_cid', data=element_cid)
            self.data = group.create_dataset('data', data=data)
        else:
            self._times = _times
            self.element_node = element_node
            self.element_cid = element_cid
            self.data = data

    def build_dataframe(self):
        """creates a pandas dataframe"""
        import pandas as pd

        headers = self.get_headers()
        # TODO: cid?
        #element_node = [self.element_node[:, 0], self.element_node[:, 1]]
        if self.nonlinear_factor not in (None, np.nan):
            column_names, column_values = build_dataframe_transient_header(self)
            data_frame = build_pandas_transient_element_node(
                self, column_values, column_names,
                headers, self.element_node, self.data)
        else:
            # Static            sxc  sxd  sxe  sxf  smax  smin    MS_tension  MS_compression
            # ElementID NodeID
            # 12        22      0.0  0.0  0.0  0.0   0.0   0.0  1.401298e-45    1.401298e-45
            #           26      0.0  0.0  0.0  0.0   0.0   0.0  1.401298e-45    1.401298e-45
            index = pd.MultiIndex.from_arrays(self.element_node.T, names=['ElementID', 'NodeID'])
            data_frame = pd.DataFrame(self.data[0], columns=headers, index=index)
            data_frame.columns.names = ['Static']
        self.data_frame = data_frame

    def envelope(self,
                 eids: np.ndarray,
                 result_name: str,
                 consider_corner_nodes: bool=False) -> np.ndarray:
        is_min = (result_name == 'min')
        func_name = 'min' if is_min else 'max'
        max_min_func = getattr(np, func_name)

        # print(''.join(self.get_stats()))
        element = self.element_cid[:, 0]
        ielement = np.searchsorted(element, eids)
        assert np.array_equal(element[ielement], eids)

        #[oxx, oyy, ozz, txy, tyz, txz, o1, o2, o3, ovmShear]
        solid_stress_strain_dict_map = {
            'max': 6,
            'min': 8,
        }
        ntime, nelement_node, _ = self.data.shape
        if result_name == 'von_mises':
            data = self.von_mises()
        elif result_name == 'max_shear':
            data = self.max_shear()
        elif result_name == 'abs_max':
            # TODO: wrong, but ok...
            data = np.abs(self.abs_principal())
            assert data.shape == (ntime, nelement_node), data.shape
        elif result_name in solid_stress_strain_dict_map:
            iresult = solid_stress_strain_dict_map[result_name]
            data = self.data[:, :, iresult]
        else:  # pragma: No cover
            raise NotImplementedError(result_name)

        is_corner_nodes = len(self.element_node) > 1 and self.element_node[1, 1] != 0
        neid = len(self.element_cid)
        assert data.ndim == 2, data.shape
        ntime, nelement_nnode = data.shape

        if is_corner_nodes:
            data = data.reshape(ntime, neid, self.nnodes_per_element)
            if consider_corner_nodes:
                # arr = np.random.rand(3, 4, 5)  # 3 matrices, 4 rows, 5 columns
                # col_max = arr.max(axis=1)  # Result shape: (3, 5)
                # eid_data = data[:, :, :].max(axis=(0, 2))
                eid_data = max_min_func(data, axis=(0, 2))
                assert eid_data.shape == (neid,), eid_data.shape
            else:
                # eid_data = data[:, :, 0].max(axis=0)
                eid_data = max_min_func(data[:, :, 0], axis=0)
                assert eid_data.shape == (neid,), eid_data.shape
        else:
            # eid_data = data.max(axis=0)
            eid_data = max_min_func(data, axis=0)
            assert eid_data.shape == (neid, ), eid_data.shape
        # TODO: move ielement to earlier
        return eid_data[ielement]

    def von_mises(self) -> np.ndarray:
        if self.is_von_mises:
            von_mises = self.data[:, :, 9]
            return von_mises
        oxx = self.data[:, :, 0]
        oyy = self.data[:, :, 1]
        ozz = self.data[:, :, 2]
        txy = self.data[:, :, 3]
        tyz = self.data[:, :, 4]
        txz = self.data[:, :, 5]
        von_mises = von_mises_3d(
            oxx, oyy, ozz, txy, tyz, txz,
            self.is_stress)
        return von_mises

    def max_shear(self) -> np.ndarray:
        if self.is_max_shear:
            max_shear = self.data[:, :, 9]
            return max_shear

        omax = self.data[:, :, 6]
        omin = self.data[:, :, 8]
        max_shear = max_shear_3d(
            omax, omin, self.is_stress)
        return max_shear

    def abs_principal(self) -> np.ndarray:
        omax = self.data[:, :, 6]
        omin = self.data[:, :, 8]
        abs_max = get_abs_max(omin, omax, dtype=omax.dtype)
        return abs_max

    @classmethod
    def _add_case(cls,
                  table_name, element_name, isubcase,
                  is_sort1, is_random, is_msc,
                  random_code, title, subtitle, label):
        num_wide = 3
        is_strain = 'Strain' in cls.__name__
        data_code = oes_real_data_code(table_name,
                                       element_name, num_wide,
                                       is_sort1=is_sort1, is_random=is_random,
                                       random_code=random_code,
                                       title=title, subtitle=subtitle, label=label,
                                       is_msc=is_msc)

        # I'm only sure about the 1s in the strains and the
        # corresponding 0s in the stresses.
        #stress / strain -> 1, 3
        # stress_bits[4] = 1 # von mises (vs. max shear)
        if is_strain:
            stress_bits = [0, 1, 0, 1, 1]
            s_code = 1
        else:
            stress_bits = [0, 0, 0, 0, 1]
            s_code = 0
            assert stress_bits[1] == 0, 'stress_bits=%s' % str(stress_bits)

        # stress
        assert stress_bits[1] == stress_bits[3], 'stress_bits=%s' % str(stress_bits)
        data_code['stress_bits'] = stress_bits
        data_code['s_code'] = s_code

        element_type = ELEMENT_NAME_TO_ELEMENT_TYPE[element_name]
        data_code['element_name'] = element_name
        data_code['element_type'] = element_type
        return data_code

    @classmethod
    def add_static_case(cls, table_name, element_name, element_node, element_cid, data, isubcase,
                        is_sort1=True, is_random=False, is_msc=True,
                        random_code=0, title='', subtitle='', label=''):
        data_code = cls._add_case(
            table_name, element_name,
            isubcase, is_sort1, is_random, is_msc,
            random_code, title, subtitle, label)
        obj = set_static_case(cls, is_sort1, isubcase, data_code,
                              set_element_cid_case, (element_node, element_cid, data))
        return obj

    @classmethod
    def add_modal_case(cls, table_name, element_name: str, element_node, element_cid, data,
                       isubcase: int,
                       modes: np.ndarray, eigns: np.ndarray, cycles: np.ndarray,
                       is_sort1: bool=True, is_random: bool=False, is_msc: bool=True,
                       random_code: int=0, title: str='', subtitle: str='', label: str=''):
        data_code = cls._add_case(
            table_name, element_name,
            isubcase, is_sort1, is_random, is_msc,
            random_code, title, subtitle, label)
        obj = set_modal_case(cls, is_sort1, isubcase, data_code,
                             set_element_cid_case, (element_node, element_cid, data),
                             modes, eigns, cycles)
        return obj

    @classmethod
    def add_transient_case(cls, table_name, element_name, element_node, element_cid, data, isubcase,
                           times,
                           is_sort1=True, is_random=False, is_msc=True,
                           random_code=0, title='', subtitle='', label=''):
        data_code = cls._add_case(
            table_name, element_name,
            isubcase, is_sort1, is_random, is_msc,
            random_code, title, subtitle, label)
        obj = set_transient_case(cls, is_sort1, isubcase, data_code,
                                 set_element_cid_case, (element_node, element_cid, data),
                                 times)
        return obj

    @classmethod
    def add_post_buckling_case(cls, table_name, element_name, element_node, element_cid, data,
                               isubcase: int,
                               modes: np.ndarray, eigrs: np.ndarray, eigis: np.ndarray,
                               is_sort1: bool=True, is_random: bool=False, is_msc: bool=True,
                               random_code: int=0, title: str='', subtitle: str='', label: str=''):
        data_code = cls._add_case(
            table_name, element_name,
            isubcase, is_sort1, is_random, is_msc,
            random_code, title, subtitle, label)
        obj = set_post_buckling_case(cls, is_sort1, isubcase, data_code,
                                     set_element_cid_case, (element_node, element_cid, data),
                                     modes, eigrs, eigis)
        return obj

    def add_eid_sort1(self, unused_etype, cid, dt, eid, unused_node_id,
                      oxx, oyy, ozz, txy, tyz, txz, o1, o2, o3,
                      unused_acos, unused_bcos, unused_ccos, unused_pressure, ovm):
        # See the CHEXA, CPENTA, or CTETRA entry for the definition of
        #    the element coordinate systems
        #   The material coordinate system (CORDM) may be:
        #    - basic system (0 or blank)
        #    - any defined system
        # (Integer > 0), or the standard internal coordinate system of the element:
        #   -1: element coordinate system (-1)
        #   -2: element system based on eigenvalue techniques to insure non bias in the
        #       element formulation
        #     C:\MSC.Software\msc_nastran_runs\ecs-2-rg.op2
        assert cid >= -2, cid
        assert eid >= 0, eid

        #print(f'dt={dt} eid={eid}')
        self._times[self.itime] = dt
        self.element_node[self.itotal, :] = [eid, 0]  # 0 is center
        # o1 += 1000
        omax, omid, omin = o123_to_max_mid_min(o1, o2, o3)
        # print(f'eid={eid} nid=0   s1={o1  :+6e}   s2={o2  :+6e}   s3={o3:+6e}')
        # print(f'eid={eid} nid=0 omax={omax:+6e} omid={omid:+6e} omin={omin:+6e}')
        # passed = (omax >= omid >= omin)
        # assert passed
        self.data[self.itime, self.itotal, :] = [oxx, oyy, ozz, txy, tyz, txz, omax, omid, omin, ovm]

        #print('element_cid[%i, :] = [%s, %s]' % (self.ielement, eid, cid))
        if self.ielement == self.nelements:
            self.ielement = 0
        self.element_cid[self.ielement, :] = [eid, cid]
        self.itotal += 1
        self.ielement += 1

    def add_node_sort1(self, dt, eid, unused_inode, node_id,
                       oxx, oyy, ozz, txy, tyz, txz, o1, o2, o3,
                       unused_acos, unused_bcos, unused_ccos, unused_pressure, ovm):
        # skipping aCos, bCos, cCos, pressure
        omax, omid, omin = o123_to_max_mid_min(o1, o2, o3)
        # print(f'eid={eid} nid={node_id}   s1={o1  :+6e}   s2={o2  :+6e}   s3={o3  :+6e}')
        # print(f'eid={eid} nid={node_id} omax={omax:+6e} omid={omid:+6e} omin={omin:+6e}')
        # passed = (omax >= omid >= omin)
        # assert passed
        self.data[self.itime, self.itotal, :] = [oxx, oyy, ozz, txy, tyz, txz, omax, omid, omin, ovm]
        #print('data[%s, %s, :] = %s' % (self.itime, self.itotal, str(self.data[self.itime, self.itotal, :])))

        #print('eid=%i node_id=%i exx=%s' % (eid, node_id, str(oxx)))
        self.element_node[self.itotal, :] = [eid, node_id]
        #self.element_node[self.ielement-1, self.inode-1, :] = [eid, node_id]
        self.itotal += 1

    def add_eid_sort2(self, unused_etype, cid, dt, eid, unused_node_id,
                      oxx, oyy, ozz, txy, tyz, txz, o1, o2, o3,
                      unused_acos, unused_bcos, unused_ccos, unused_pressure, ovm):
        #itime = self.ielement
        #ielement = self.itotal
        #itotal = self.itime
        #print(self.ntimes, self.nelements, self.ntotal, self.nnodes)
        itime = self.itotal // self.nnodes
        ielement = self.itime
        itotal = self.itotal
        assert cid >= -2, cid
        assert eid >= 0, eid

        #try:
        self._times[itime] = dt
            #print(f'dt={dt} eid={eid} ielement={ielement} -> itime={itime} itotal={itotal}')
        #except IndexError:
            #print(f'*dt={dt} eid={eid} ielement={ielement} -> itime={itime} itotal={itotal}')
            #self.itime += 1
            #self.ielement += 1
            #return
        self.element_node[itotal, :] = [eid, 0]  # 0 is center

        omax, omid, omin = o123_to_max_mid_min(o1, o2, o3)
        self.data[itime, itotal, :] = [oxx, oyy, ozz, txy, tyz, txz, omax, omid, omin, ovm]

        #print('element_cid[%i, :] = [%s, %s]' % (self.ielement, eid, cid))
        #if self.ielement == self.nelements:
            #self.ielement = 0
        self.element_cid[ielement, :] = [eid, cid]
        #self.itime += 1
        self.itotal += 1
        self.ielement += 1
        #print('self._times', self._times)

    def add_node_sort2(self, dt, eid, unused_inode, node_id,
                       oxx, oyy, ozz, txy, tyz, txz, o1, o2, o3,
                       unused_acos, unused_bcos, unused_ccos, unused_pressure, ovm):
        #ielement = self.ielement
        #itotal = self.itotal
        #itime = self.itime
        #itime=0 ielement=1 itotal=1
        #itime=0 ielement=1 itotal=2
        #itime=0 ielement=1 itotal=3
        #itime=0 ielement=1 itotal=4

        #ielement = self.ielement
        #itime = (self.itime - 1) % self.nelements
        #itime = self.itime - 1
        nnodes = self.nnodes
        itime = self.itotal // nnodes
        itotal = self.itotal
        #ielement = self.ielement - 1
        #ielement = self.itime
        #inode = self.itotal % nnodes
        #itotal2 = (self.ielement - 1) * nnodes + inode
        #print(f'  itime={itime} itotal={itotal}; nid={node_id}; '
              #f'ielement={ielement} inode={inode} -> itotal2={itotal2}')

        # skipping aCos, bCos, cCos, pressure
        omax, omid, omin = o123_to_max_mid_min(o1, o2, o3)
        self.data[itime, itotal, :] = [oxx, oyy, ozz, txy, tyz, txz, omax, omid, omin, ovm]
        #print('data[%s, %s, :] = %s' % (self.itime, self.itotal, str(self.data[self.itime, self.itotal, :])))

        #print('eid=%i node_id=%i exx=%s' % (eid, node_id, str(oxx)))
        self.element_node[itotal, :] = [eid, node_id]
        #self.element_node[ielement-1, inode-1, :] = [eid, node_id]
        self.itotal += 1

    def __eq__(self, table):  # pragma: no cover
        assert self.is_sort1 == table.is_sort1
        self._eq_header(table)
        if not np.array_equal(self.data, table.data):
            msg = 'table_name=%r class_name=%s\n' % (self.table_name, self.__class__.__name__)
            msg += '%s\n' % str(self.code_information())
            i = 0
            for itime in range(self.ntimes):
                for ieid, eid_nid in enumerate(self.element_node):
                    (eid, nid) = eid_nid

                    s1 = self.data[itime, ieid, :6]  # component stresses
                    s2 = table.data[itime, ieid, :6]

                    t1 = self.data[itime, ieid, :]
                    t2 = table.data[itime, ieid, :]
                    (oxx1, oyy1, ozz1, txy1, tyz1, txz1, o11, o21, o31, ovm1) = t1
                    (oxx2, oyy2, ozz2, txy2, tyz2, txz2, o12, o22, o32, ovm2) = t2
                    if not np.array_equal(s1, s2):
                        msg += (
                            '(%s, %s)    (%s, %s, %s, %s, %s, %s)\n'
                                '%s      (%s, %s, %s, %s, %s, %s)\n' % (
                                eid, nid,
                                oxx1, oyy1, ozz1, txy1, tyz1, txz1,
                                ' ' * (len(str(eid)) + len(str(nid)) + 2),
                                oxx2, oyy2, ozz2, txy2, tyz2, txz2))
                        i += 1
                        if i > 10:
                            print(msg)
                            raise ValueError(msg)
                        continue

                    if not np.allclose(t1, t2):
                        msg += (
                            '(%s, %s)    (%s, %s, %s, %s)\n'
                            '%s      (%s, %s, %s, %s)\n' % (
                                eid, nid,
                                o11, o21, o31, ovm1,
                                ' ' * (len(str(eid)) + len(str(nid)) + 2),
                                o12, o22, o32, ovm2))
                        i += 1
                        if i > 10:
                            print(msg)
                            raise ValueError(msg)
                #print(msg)
                if i > 0:
                    raise ValueError(msg)
        return True

    @property
    def nnodes_per_element(self) -> int:
        return self.nnodes_per_element_no_centroid + 1

    @property
    def nnodes_per_element_no_centroid(self) -> int:
        if self.element_type == 39: # CTETRA
            nnodes = 4
        elif self.element_type == 67: # CHEXA
            nnodes = 8
        elif self.element_type == 68: # CPENTA
            nnodes = 6
        elif self.element_type == 255: # CPYRAM
            nnodes = 5
        else:  # pragma: no cover
            raise NotImplementedError(f'element_name={self.element_name} self.element_type={self.element_type}')
        return nnodes

    def get_stats(self, short: bool=False) -> list[str]:
        if not self.is_built:
            return [
                f'<{self.__class__.__name__}>; table_name={self.table_name!r}\n',
                f'  ntimes: {self.ntimes:d}\n',
                f'  ntotal: {self.ntotal:d}\n',
            ]

        nelements = self.nelements
        ntimes = self.ntimes
        #ntotal = self.ntotal
        try:
            nnodes_per_element = self.element_node.shape[0] // nelements
        except ZeroDivisionError:
            nnodes_per_element = '???'
        nnodes = self.element_node.shape[0]

        msg = []

        if self.nonlinear_factor not in (None, np.nan):  # transient
            msg.append('  type=%s ntimes=%i nelements=%i nnodes=%i\n  nnodes_per_element=%s (including centroid)\n'
                       % (self.__class__.__name__, ntimes, nelements, nnodes, nnodes_per_element))
            ntimes_word = 'ntimes'
        else:
            msg.append('  type=%s nelements=%d nnodes=%d\n  nodes_per_element=%d (including centroid)\n'
                       % (self.__class__.__name__, nelements, nnodes, nnodes_per_element))
            ntimes_word = '1'
        msg.append('  eType, cid\n')
        headers = self.get_headers()
        n = len(headers)
        msg.append('  data: [%s, nnodes, %d] where %d=[%s]\n' % (
            ntimes_word, n, n, str(', '.join(headers))))
        msg.append(f'  element_node.shape = {self.element_node.shape}\n')
        msg.append(f'  element_cid.shape = {self.element_cid.shape}\n')
        msg.append(f'  data.shape = {self.data.shape}\n')
        msg.append('  element name: %s\n' % self.element_name)
        msg += self.get_data_code()
        #print(''.join(msg))
        return msg


    def get_element_index(self, eids):
        # elements are always sorted; nodes are not
        itot = np.searchsorted(eids, self.element_node[:, 0])  #[0]
        return itot

    def eid_to_element_node_index(self, eids):
        #ind = np.ravel([np.searchsorted(self.element_node[:, 0] == eid) for eid in eids])
        ind = np.searchsorted(eids, self.element_node[:, 0])
        #ind = ind.reshape(ind.size)
        #ind.sort()
        return ind

    def write_csv(self, csv_file: TextIO,
                  is_exponent_format: bool=False,
                  is_mag_phase: bool=False, is_sort1: bool=True,
                  write_header: bool=True):
        """
        Stress Table - Solid CHEXA
        --------------------------
        Flag,  SubcaseID, iTime, EID, NID,     CID,     Sxx,    Syy,     Szz,    Sxy,     Syz,    Szx
        10,    1,             0, 309,   0,       0, 142.500, 31.219, 589.597, 676.54, 1138.79, 213.68
        10,    1,             0, 309, 101,       0, 54.3342, 21.410, 553.325, 354.90, 87.0544, 20.192
        10,    1,             0, 309, 102,       0, 113.506, 80.846, 53.6931, 29.766, 1033.05, 19.109
        10,    1,             0, 309, 103,       0, 176.472, 721.43, 17.1733, 301.05, 374.726, 372.35
        10,    1,             0, 309, 104,       0, 21.1607, 81.748, 66.6382, 21.331, 783.494, 796.67
        10,    1,             0, 309, 105,       0, 114.81,  833.38, 271.391, 65.490, 1773.04, 74.355
        10,    1,             0, 309, 106,       0, 90.7118, 84.456, 783.545, 573.50, 1623.11, 08.347
        10,    1,             0, 309, 107,       0, 46.5565, 97.237, 540.913, 777.87, 1824.50, 13.681
        10,    1,             0, 309, 108,       0, 87.779,  923.94, 211.281, 4.9608, 1387.49, 945.86
        """
        name = str(self.__class__.__name__)
        if write_header:
            csv_file.write('# %s\n' % name)
            headers = ['Flag', 'SubcaseID', 'iTime', 'Eid', 'Nid', 'CID',
                       'Sxx', 'Syy', 'Szz', 'Sxy', 'Syz', 'Sxz']
            csv_file.write('# ' + ','.join(headers) + '\n')

        # stress vs. strain
        flag = 10 if 'Stress' in name else 11

        #nnodes, msg_temp = _get_f06_header_nnodes(self, is_mag_phase)
        nnodes = self.nnodes
        isubcase = self.isubcase

        # write the f06
        ntimes = self.data.shape[0]
        eids2 = self.element_node[:, 0]
        nodes = self.element_node[:, 1]

        eids3 = self.element_cid[:, 0]
        cids3 = self.element_cid[:, 1]

        #fdtype = self.data.dtype
        # oxx = self.data[:, :, 0]
        # oyy = self.data[:, :, 1]
        # ozz = self.data[:, :, 2]
        # txy = self.data[:, :, 3]
        # tyz = self.data[:, :, 4]
        # txz = self.data[:, :, 5]
        #o1 = self.data[:, :, 6]
        #o2 = self.data[:, :, 7]
        #o3 = self.data[:, :, 8]
        #ovm = self.data[:, :, 9]
        #p = (o1 + o2 + o3) / -3.
        eid_len = '%d' % len(str(eids2.max()))
        nid_len = '%d' % len(str(nodes.max()))
        #cid_len = '%d' % len(str(cids3.max()))

        #zero = ' 0.000000E+00'
        for itime in range(ntimes):
            #dt = self._times[itime]
            #header = _eigenvalue_header(self, header, itime, ntimes, dt)
            #f06_file.write(''.join(header + msg_temp))

            #print("self.data.shape=%s itime=%s ieids=%s" % (str(self.data.shape), itime, str(ieids)))
            oxx = self.data[itime, :, 0]
            oyy = self.data[itime, :, 1]
            ozz = self.data[itime, :, 2]
            txy = self.data[itime, :, 3]
            tyz = self.data[itime, :, 4]
            txz = self.data[itime, :, 5]
            #o1 = self.data[itime, :, 6]
            #o2 = self.data[itime, :, 7]
            #o3 = self.data[itime, :, 8]
            #ovm = self.data[itime, :, 9]
            #vi = v[itime, :, :, :]
            #pi = p[itime, :]

            cnnodes = nnodes + 1
            for i, deid, node_id, oxxi, oyyi, ozzi, txyi, tyzi, txzi in zip(
                    count(), eids2, nodes, oxx, oyy, ozz, txy, tyz, txz):

                if is_exponent_format:
                    [oxxi, oyyi, ozzi, txyi, tyzi, txzi] = write_floats_13e_long(
                        [oxxi, oyyi, ozzi, txyi, tyzi, txzi])

                if i % cnnodes == 0:
                    j = np.where(eids3 == deid)[0][0]
                    cid = cids3[j]
                csv_file.write(f'{flag}, {isubcase}, {itime}, {deid:{eid_len}d}, '
                               f'{node_id:{nid_len}d}, {cid}, '
                               f'{oxxi}, {oyyi}, {ozzi}, {txyi}, {tyzi}, {txzi}\n')
        return

    def check_update(self, check_code: bool=False):  # pragma: no cover
        if not check_code:
            return
        i1 = 6
        i2 = 7
        i3 = 8
        iovm = 9
        ovm1 = self.data[:, :, iovm].copy()
        o11 = self.data[:, :, i1].copy()
        o21 = self.data[:, :, i2].copy()
        o31 = self.data[:, :, i3].copy()

        self.update_data_components()
        ovm2 = self.data[:, :, iovm]
        o12 = self.data[:, :, i1]
        o22 = self.data[:, :, i2]
        o32 = self.data[:, :, i3]
        stress = 'stress' if self.is_stress else 'strain'
        von_mises = f'von_mises_{stress}' if self.is_von_mises else f'max_shear_{stress}'
        assert np.allclose(ovm1, ovm2), (von_mises, ovm1.ravel(), ovm2.ravel())

        for i, o11i, o21i, o31i, o12i, o22i, o32i in zip(count(), o11, o21, o31, o12, o22, o32):
            assert np.allclose(o11i, o12i), str((i, (o11i, o21i, o31i), (o12i, o22i, o32i)))
            assert np.allclose(o21i, o22i), str((i, (o11i, o21i, o31i), (o12i, o22i, o32i)))
            assert np.allclose(o31i, o32i), str((i, (o11i, o21i, o31i), (o12i, o22i, o32i)))

        assert np.allclose(o11, o12), (f'o1_{stress}', o11.ravel(), o12.ravel(), np.abs(o11-o12).max())
        assert np.allclose(o21, o22), (f'o2_{stress}', o21.ravel(), o22.ravel(), np.abs(o21-o22).max())
        assert np.allclose(o31, o32), (f'o3_{stress}', o31.ravel(), o32.ravel(), np.abs(o31-o32).max())

    def write_f06(self, f06_file: TextIO, header=None, page_stamp: str='PAGE %s',
                  page_num: int=1, is_mag_phase: bool=False, is_sort1: bool=True):
        self.check_update()

        calculate_directional_vectors = True
        if header is None:
            header = []
        nnodes, msg_temp = _get_f06_header_nnodes(self, is_mag_phase)

        # write the f06
        ntimes = self.data.shape[0]

        eids2 = self.element_node[:, 0]
        nodes = self.element_node[:, 1]

        eids3 = self.element_cid[:, 0]
        cids3 = self.element_cid[:, 1]

        fdtype = self.data.dtype
        oxx = self.data[:, :, 0]
        oyy = self.data[:, :, 1]
        ozz = self.data[:, :, 2]
        txy = self.data[:, :, 3]
        tyz = self.data[:, :, 4]
        txz = self.data[:, :, 5]
        # omax = self.data[:, :, 6]
        # omid = self.data[:, :, 7]
        # omin = self.data[:, :, 8]
        # ovm = self.data[:, :, 9]

        nnodes_total = self.data.shape[1]
        if calculate_directional_vectors:

            # (1, 930, 3)(1, 930, 3, 3)
            eigenvalues, v = calculate_principal_eigenvectors4(
                ntimes, nnodes_total,
                oxx.ravel(), oyy.ravel(), ozz.ravel(), txy.ravel(), txz.ravel(), tyz.ravel(),
                self.is_stress, fdtype)
            assert len(eigenvalues.shape) == 3
            # print(eigenvalues[0, :, :])
            omax = eigenvalues[:, :, 2]  # max
            omid = eigenvalues[:, :, 1]  # mid
            omin = eigenvalues[:, :, 0]  # min
        else:  # pragma: no cover
            v = np.zeros((ntimes, nnodes, 3, 3), dtype=fdtype)
            omax = self.data[:, :, 6]
            omid = self.data[:, :, 7]
            omin = self.data[:, :, 8]
        p = (omax + omid + omin) / -3.
        assert p.shape == omax.shape

        for itime in range(ntimes):
            dt = self._times[itime]
            header = _eigenvalue_header(self, header, itime, ntimes, dt)
            f06_file.write(''.join(header + msg_temp))

            #print("self.data.shape=%s itime=%s ieids=%s" % (str(self.data.shape), itime, str(ieids)))
            oxx = self.data[itime, :, 0]
            oyy = self.data[itime, :, 1]
            ozz = self.data[itime, :, 2]
            txy = self.data[itime, :, 3]
            tyz = self.data[itime, :, 4]
            txz = self.data[itime, :, 5]

            omaxi1 = omax[itime, :]
            omidi1 = omid[itime, :]
            omini1 = omin[itime, :]

            # omaxi1 = self.data[itime, :, 6]
            # omidi1 = self.data[itime, :, 7]
            # omini1 = self.data[itime, :, 8]
            #assert not np.allclose(omaxi1, omidi1)
            #assert not np.allclose(omidi1, omini1)
            assert omaxi1.shape == oxx.shape, (omaxi1.shape, oxx.shape)

            # assert np.all(omax >= omid), 'omax < omid'
            # assert np.all(omid >= omin), 'omid < omin'
            ovm = self.data[itime, :, 9]
            vi = v[itime, :, :, :]
            # vmaxi = v[itime, :, :, :]
            # vmidi = v[itime, :, :, :]
            # vmini = v[itime, :, :, :]
            # pi = (omax + omid + omin) / -3.
            pi = p[itime, :]

            # %(node_id, oxxi, txyi, omaxi, dv[0, 1], dv[0, 2], dv[0, 0], pii, ovmi,
            #        '', oyyi, tyzi, omini, dv[2, 1], dv[2, 2], dv[2, 0],
            #        '', ozzi, txzi, omidi, dv[1, 1], dv[1, 2], dv[1, 0]))

            cnnodes = nnodes + 1
            for i, deid, node_id, doxx, doyy, dozz, dtxy, dtyz, dtxz, domax, domid, domin, dp, dv, dovm in zip(
                    count(), eids2, nodes, oxx, oyy, ozz, txy, tyz, txz, omaxi1, omidi1, omini1, pi, vi, ovm):
                # assert domax >= domid, (domax, domid, domin)
                # assert domid >= domin, (domax, domid, domin)

                # assert domax >= domid >= domin, 'o1 >= o2 >= o3; eid=%s o1=%e o2=%e o3=%e' % (deid, domax, domid, domin)
                [oxxi, oyyi, ozzi, txyi, tyzi, txzi, omaxi, omidi, omini, pii, ovmi] = write_floats_13e(
                    [doxx, doyy, dozz, dtxy, dtyz, dtxz, domax, domid, domin, dp, dovm])
                i0 = 0
                i1 = 1
                i2 = 2
                if i % cnnodes == 0:
                    j = np.where(eids3 == deid)[0][0]
                    cid = cids3[j]
                    f06_file.write('0  %8s    %8iGRID CS  %i GP\n' % (deid, cid, nnodes))
                    f06_file.write(
                        '0              %8s  X  %-13s  XY  %-13s   A  %-13s  LX%5.2f%5.2f%5.2f  %-13s   %s\n'
                        '               %8s  Y  %-13s  YZ  %-13s   B  %-13s  LY%5.2f%5.2f%5.2f\n'
                        '               %8s  Z  %-13s  ZX  %-13s   C  %-13s  LZ%5.2f%5.2f%5.2f\n'
                        % ('CENTER', oxxi, txyi, omaxi, dv[0, i0], dv[0, i2], dv[0, i1], pii, ovmi,
                                 '', oyyi, tyzi, omini, dv[1, i0], dv[1, i2], dv[1, i1],
                                 '', ozzi, txzi, omidi, dv[2, i0], dv[2, i2], dv[2, i1]))
                else:
                    f06_file.write(
                        '0              %8s  X  %-13s  XY  %-13s   A  %-13s  LX%5.2f%5.2f%5.2f  %-13s   %s\n'
                        '               %8s  Y  %-13s  YZ  %-13s   B  %-13s  LY%5.2f%5.2f%5.2f\n'
                        '               %8s  Z  %-13s  ZX  %-13s   C  %-13s  LZ%5.2f%5.2f%5.2f\n'
                        % (node_id, oxxi, txyi, omaxi, dv[0, i0], dv[0, i2], dv[0, i1], pii, ovmi,
                                '', oyyi, tyzi, omini, dv[1, i0], dv[1, i2], dv[1, i1],
                                '', ozzi, txzi, omidi, dv[2, i0], dv[2, i2], dv[2, i1]))
                i += 1
            f06_file.write(page_stamp % page_num)
            page_num += 1
        return page_num - 1

    def check_stress(self, itime: int=0):
        msg = ''
        # eids = self.element_node[:, 0]
        # nodes = self.element_node[:, 1]
        nelement = len(np.unique(self.element_node[:, 0]))
        nnodes_centroid = self.nnodes_per_element
        # print(len(self.element_node[:, 0]), nelement, nnodes_centroid, self.element_node.shape)
        # print(self.element_node.tolist())
        nelement_node = len(self.element_node[:, 0])
        if nelement_node != (nelement * nnodes_centroid):
            warnings.warn('incorrect data in RealSolidArray.element_node')
            return
        eids2 = self.element_node[:, 0].reshape(nelement, nnodes_centroid)
        nids2 = self.element_node[:, 1].reshape(nelement, nnodes_centroid)
        eid_max = eids2.max()
        nid_max = nids2.max()
        len_eid = max(3, len(str(eid_max)))
        len_nid = max(3, len(str(nid_max)))
        eid_nid1_fmt = f'%{len_eid}d %{len_nid}d'
        eid_nid2_fmt = f'%{len_eid}s %{len_nid}s'
        omax = self.data[:, :, 6]
        omid = self.data[:, :, 7]
        omin = self.data[:, :, 8]
        nfailed = 0
        for eid, nid, omaxi, omidi, omini in zip(
                eids2.ravel(), nids2.ravel(),
                omax[itime, :].ravel(),
                omid[itime, :].ravel(),
                omin[itime, :].ravel(),):
            # passed = (omaxi >= omidi >= omini)
            passed = (
                (omaxi >= omidi or np.isclose(omaxi, omidi)) and
                (omidi >= omini or np.isclose(omidi, omini))
            )
            # passed = (omaxi+1e-10 >= omidi >= omini-1e-10)
            if not passed:
                # passed1 = (omaxi+1e-10 >= omidi)
                # passed2 = (omidi+1e-10 >= omini)
                passed1 = (omaxi >= omidi)
                passed2 = (omidi >= omini)
                eid_nid1 = eid_nid1_fmt % (eid, nid)
                msg += f'({eid_nid1}) {omaxi:+15.8e} {omidi:+15.8e} {omini:+15.8e} {passed1} {passed2}\n'
                nfailed += 1
        if msg:
            eid_nid2 = eid_nid2_fmt % ("eid", "nid")
            titles = f'({eid_nid2}) {"omax":15s} {"omid":15s} {"omin":15s}\n'
            nfailed_max = 20
            if nfailed > nfailed_max:
                msg = '\n'.join(msg.split('\n')[:nfailed_max])
                msg += f'\n limiting to {nfailed_max}/{nfailed}'
            raise RuntimeError(f'Incorrect stress order for {SOLID_PRINCIPAL_METHOD}\n' + titles + msg)

    def write_op2(self, op2_file, op2_ascii, itable, new_result,
                  date, is_mag_phase=False, endian='>'):
        """writes an OP2"""
        self.check_update()
        import inspect
        calculate_directional_vectors = True
        frame = inspect.currentframe()
        call_frame = inspect.getouterframes(frame, 2)
        op2_ascii.write(f'{self.__class__.__name__}.write_op2: {call_frame[1][3]}\n')

        if itable == -1:
            #print('***************', itable)
            self._write_table_header(op2_file, op2_ascii, date)
            itable = -3

        # if isinstance(self.nonlinear_factor, float):
        #     op2_format = '%sif' % (7 * self.ntimes)
        #     raise NotImplementedError()
        # else:
        #     op2_format = 'i21f'
        # s = Struct(op2_format)
        # nnodes_expected = self.nnodes

        eids2 = self.element_node[:, 0]
        nodes = self.element_node[:, 1]
        # nelements_nodes = len(nodes)

        eids3 = self.element_cid[:, 0]
        cids3 = self.element_cid[:, 1]
        element_device = eids3 * 10 + self.device_code

        # table 4 info
        #ntimes = self.data.shape[0]
        nnodes = self.data.shape[1]
        ueids2 = np.unique(eids2)
        nelements = len(ueids2)
        if len(ueids2) == 1 and len(eids3) != 1:
            raise RuntimeError(f'SORT2: ueids2={ueids2}')

        # 21 = 1 node, 3 principal, 6 components, 9 vectors, 2 p/ovm
        #ntotal = ((nnodes * 21) + 1) + (nelements * 4)
        nnodes_centroid = self.nnodes_per_element
        nnodes_no_centroid = self.nnodes_per_element_no_centroid
        ntotali = 4 + 21 * nnodes_no_centroid
        ntotali = self.num_wide
        ntotal = ntotali * nelements

        #print('shape = %s' % str(self.data.shape))
        assert nnodes > 1, nnodes
        #assert self.ntimes == 1, self.ntimes

        op2_ascii.write(f'  ntimes = {self.ntimes}\n')
        ntimes = self.ntimes

        #print('ntotal=%s' % (ntotal))
        if not self.is_sort1:
            raise NotImplementedError('SORT2')
        #op2_format = endian + b'2i6f'

        idtype = self.element_cid.dtype
        fdtype = self.data.dtype
        if self.size == fdtype.itemsize:
            grid_bytes = b'GRID'
        else:
            # warnings.warn(f'downcasting {self.class_name}...')
            idtype = np.int32(1)
            fdtype = np.float32(1.0)
            grid_bytes = b'GRID'

        cen_array = np.full(nelements, grid_bytes, dtype='|S4')
        nnodes_no_centroid_array = np.full(nelements, nnodes_no_centroid, dtype=idtype)

        assert len(element_device) == nelements, f'element_device.shape={element_device.shape} nelements={nelements}'
        assert len(cen_array) == nelements, f'cen_array.shape={cen_array.shape} cen_array={nelements}'
        element_wise_data = to_column_bytes([
            element_device, # ints
            cids3, # ints
            cen_array, # bytes
            nnodes_no_centroid_array, # ints
        ], fdtype, debug=False)

        oxx = self.data[:, :, 0]
        oyy = self.data[:, :, 1]
        ozz = self.data[:, :, 2]
        txy = self.data[:, :, 3]
        tyz = self.data[:, :, 4]
        txz = self.data[:, :, 5]

        ovm = self.data[:, :, 9]

        # speed up transient cases, but slightly slows down static cases
        data_out = np.empty((nelements, 4+21*nnodes_centroid), dtype=fdtype)

        # setting:
        #  - CTETRA: [element_device, cid, 'CEN/', 4]
        #  - CPYRAM: [element_device, cid, 'CEN/', 5]
        #  - CPENTA: [element_device, cid, 'CEN/', 6]
        #  - CHEXA:  [element_device, cid, 'CEN/', 8]
        data_out[:, :4] = element_wise_data

        # we could tack the nodes on, so we don't have to keep stacking it
        # but we run into issues with datai
        #
        # total=nelements_nodes
        #nodes_view = nodes.view(fdtype).reshape(nelements, nnodes_centroid)
        #inode = np.arange(nnodes_centroid)
        #data_out[:, 4+inode*21] = nodes_view[:, inode]

        # v is the (3, 3) eigenvector for every time and every element
        if calculate_directional_vectors:
            eigenvalues, v = calculate_principal_eigenvectors4(
                ntimes, nnodes,
                oxx.ravel(), oyy.ravel(), ozz.ravel(), txy.ravel(), txz.ravel(), tyz.ravel(),
                self.is_stress, fdtype)
            omax = eigenvalues[:, :, 2]  # max
            omid = eigenvalues[:, :, 1]  # mid
            omin = eigenvalues[:, :, 0]  # min
        else:
            v = np.zeros((ntimes, nnodes, 3, 3), dtype=fdtype)
            omax = self.data[:, :, 6]
            omid = self.data[:, :, 7]
            omin = self.data[:, :, 8]
        p = (omax + omin + omid) / -3.

        op2_ascii.write(f'nelements={nelements:d}\n')
        for itime in range(self.ntimes):
            vi = v[itime, :, :, :]
            self._write_table_3(op2_file, op2_ascii, new_result, itable, itime)

            # record 4
            #print('stress itable = %s' % itable)
            itable -= 1
            header = [4, itable, 4,
                      4, 1, 4,
                      4, 0, 4,
                      4, ntotal, 4,
                      4 * ntotal]
            op2_file.write(pack('%ii' % len(header), *header))
            op2_ascii.write('r4 [4, 0, 4]\n')
            op2_ascii.write(f'r4 [4, {itable:d}, 4]\n')
            op2_ascii.write(f'r4 [4, {4 * ntotal:d}, 4]\n')

            col_inputs = [
                nodes,
                oxx[itime, :], txy[itime, :], omax[itime, :], vi[:, 0, 1], vi[:, 0, 2], vi[:, 0, 0], p[itime, :], ovm[itime, :],
                oyy[itime, :], tyz[itime, :], omin[itime, :], vi[:, 1, 1], vi[:, 1, 2], vi[:, 1, 0],
                ozz[itime, :], txz[itime, :], omid[itime, :], vi[:, 2, 1], vi[:, 2, 2], vi[:, 2, 0],
            ]

            # stack each output by columns and fix any dtypes
            datai = to_column_bytes(col_inputs, fdtype)
            #datai2 = datai.reshape(nelements, 21*nnodes_centroid)
            #data_out = np.hstack([element_wise_data, datai2])
            #data_out[:, 4:] = datai2

            # switch datai to element format and put it in the output buffer
            data_out[:, 4:] = datai.reshape(nelements, 21*nnodes_centroid)
            op2_file.write(data_out)

            itable -= 1
            header = [4 * ntotal,]
            op2_file.write(pack('i', *header))
            op2_ascii.write('footer = %s\n' % header)
            new_result = False
        return itable


class RealSolidStressArray(RealSolidArray, StressObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        RealSolidArray.__init__(self, data_code, is_sort1, isubcase, dt)
        StressObject.__init__(self, data_code, isubcase)

    @property
    def headers(self) -> list[str]:
        if self.is_von_mises:
            von_mises = 'von_mises'
        else:
            von_mises = 'max_shear'
        headers = ['oxx', 'oyy', 'ozz', 'txy', 'tyz', 'txz', 'omax', 'omid', 'omin', von_mises]
        return headers


class RealSolidStrainArray(RealSolidArray, StrainObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        RealSolidArray.__init__(self, data_code, is_sort1, isubcase, dt)
        StrainObject.__init__(self, data_code, isubcase)

    @property
    def headers(self) -> list[str]:
        if self.is_von_mises:
            von_mises = 'von_mises'
        else:
            von_mises = 'max_shear'
        headers = ['exx', 'eyy', 'ezz', 'exy', 'eyz', 'exz', 'emax', 'emid', 'emin', von_mises]
        return headers

def _get_solid_msgs(self: RealSolidArray):
    if self.is_von_mises:
        von_mises = 'VON MISES'
    else:
        von_mises = 'MAX SHEAR'

    if self.is_stress:
        base_msg = [
            '0                CORNER        ------CENTER AND CORNER POINT STRESSES---------       DIR.  COSINES       MEAN                   \n',
            '  ELEMENT-ID    GRID-ID        NORMAL              SHEAR             PRINCIPAL       -A-  -B-  -C-     PRESSURE       %s \n' % von_mises]
        tetra_msg = ['                   S T R E S S E S   I N    T E T R A H E D R O N   S O L I D   E L E M E N T S   ( C T E T R A )\n', ]
        penta_msg = ['                    S T R E S S E S   I N   P E N T A H E D R O N   S O L I D   E L E M E N T S   ( P E N T A )\n', ]
        hexa_msg = ['                      S T R E S S E S   I N   H E X A H E D R O N   S O L I D   E L E M E N T S   ( H E X A )\n', ]
        pyram_msg = ['                      S T R E S S E S   I N   P Y R A M I D   S O L I D   E L E M E N T S   ( P Y R A M )\n', ]
    else:
        base_msg = [
            '0                CORNER        ------CENTER AND CORNER POINT  STRAINS---------       DIR.  COSINES       MEAN                   \n',
            '  ELEMENT-ID    GRID-ID        NORMAL              SHEAR             PRINCIPAL       -A-  -B-  -C-     PRESSURE       %s \n' % von_mises]
        tetra_msg = ['                     S T R A I N S   I N    T E T R A H E D R O N   S O L I D   E L E M E N T S   ( C T E T R A )\n', ]
        penta_msg = ['                      S T R A I N S   I N   P E N T A H E D R O N   S O L I D   E L E M E N T S   ( P E N T A )\n', ]
        hexa_msg = ['                        S T R A I N S   I N   H E X A H E D R O N   S O L I D   E L E M E N T S   ( H E X A )\n', ]
        pyram_msg = ['                        S T R A I N S   I N   P Y R A M I D   S O L I D   E L E M E N T S   ( P Y R A M )\n', ]

    tetra_msg += base_msg
    penta_msg += base_msg
    hexa_msg += base_msg
    pyram_msg += base_msg
    return tetra_msg, penta_msg, hexa_msg, pyram_msg

def _get_f06_header_nnodes(self: RealSolidArray, is_mag_phase=True):
    tetra_msg, penta_msg, hexa_msg, pyram_msg = _get_solid_msgs(self)
    if self.element_type == 39: # CTETRA
        msg = tetra_msg
        nnodes = 4
    elif self.element_type == 67: # CHEXA
        msg = hexa_msg
        nnodes = 8
    elif self.element_type == 68: # CPENTA
        msg = penta_msg
        nnodes = 6
    elif self.element_type == 255: # CPYRAM
        msg = pyram_msg
        nnodes = 5
    else:  # pragma: no cover
        msg = f'element_name={self.element_name} self.element_type={self.element_type}'
        raise NotImplementedError(msg)
    return nnodes, msg


def calculate_principal_eigenvectors5(self,
                                      ntimes: int, nelements: int, nnodes: int,
                                      oxx: np.ndarray, oyy: np.ndarray, ozz: np.ndarray,
                                      txy: np.ndarray, txz: np.ndarray, tyz: np.ndarray,
                                      dtype):
    """
    For 10 CTETRA elements (5 nodes) with 2 times, the shape would be:
    >>> (ntimes, nelements, nnodes, 3, 3)
    >>> (2,      10,        5,      3, 3)

    TODO: scale by 2 for strain
    """
    a_matrix = np.empty((ntimes, nelements, nnodes, 3, 3), dtype=self.analysis_fmt)

    # we're only filling the lower part of the A matrix
    a_matrix[:, :, :, 0, 0] = oxx
    a_matrix[:, :, :, 1, 1] = oyy
    a_matrix[:, :, :, 2, 2] = ozz
    a_matrix[:, :, :, 1, 0] = txy
    a_matrix[:, :, :, 2, 0] = txz
    a_matrix[:, :, :, 2, 1] = tyz

    # _lambda: ntimes, nelements, nnodes, (3)
    # v:       ntimes, nelements, nnodes, (3, 3)
    (_lambda, v) = eigh(a_matrix)  # a hermitian matrix is a symmetric-real matrix
    return _lambda, v


def calculate_principal_eigenvectors4(ntimes: int, nelements: int,
                                      oxx: np.ndarray, oyy: np.ndarray, ozz: np.ndarray,
                                      txy: np.ndarray, txz: np.ndarray, tyz: np.ndarray,
                                      is_stress: bool, dtype):
    """
    For 10 CTETRA elements (5 nodes) with 2 times, the shape would be:
    >>> (ntimes, nelements*nnodes, 3, 3)
    >>> (2,      10*5,             3, 3)

    TODO: scale by 2 for strain

    Parameters
    ----------
    oxx : (ntimes, nnodes) np.ndarray

    Returns
    -------
    eigenvalues : (ntimes*nnodes, 3)
        the eigenvalues
    eigenvectors : (ntimes*nnodes, 3, 3)
        the eigenvectors

    """
    a_matrix = np.zeros((ntimes*nelements, 3, 3), dtype=dtype)

    # we're only filling the lower part of the A matrix
    a_matrix[:, 0, 0] = oxx
    a_matrix[:, 1, 1] = oyy
    a_matrix[:, 2, 2] = ozz
    if is_stress:
        a_matrix[:, 1, 0] = txy
        a_matrix[:, 2, 0] = txz
        a_matrix[:, 2, 1] = tyz
    else:
        a_matrix[:, 1, 0] = txy / 2.
        a_matrix[:, 2, 0] = txz / 2.
        a_matrix[:, 2, 1] = tyz / 2.

    # eigenvalues:  ntimes, nnodes, (3)
    # eigenvectors: ntimes, nnodes, (3, 3)
    eigenvalues, eigenvectors = eigh(a_matrix)  # a hermitian matrix is a symmetric-real matrix
    eigenvalues = eigenvalues.reshape(ntimes, nelements, 3)
    eigenvectors = eigenvectors.reshape(ntimes, nelements, 3, 3)
    # eigenvectors[np.where(eigenvectors < 0.1)] = 0.
    return eigenvalues, eigenvectors


def set_element_cid_case(cls: RealSolidArray, data_code, is_sort1, isubcase,
                         element_node, element_cid, data, times):
    assert element_node.ndim == 2, element_node.shape
    assert element_cid.ndim == 2, element_cid.shape
    assert data.ndim == 3, data.shape
    ntimes = data.shape[0]
    nnodes = data.shape[1]
    dt = times[0]
    obj = cls(data_code, is_sort1, isubcase, dt)

    nresults = data.shape[1]
    nelements = element_cid.shape[0]

    # including centroid
    nnodes = nresults // nelements
    assert nresults % nelements == 0

    obj.element_node = element_node
    obj.element_cid = element_cid
    obj.nnodes = nnodes
    obj.data = data

    obj.ntimes = ntimes
    obj.ntotal = nnodes
    obj.nelements = nelements
    obj._times = times
    obj.update_data_components()
    return obj

def o123_to_max_mid_min(o1, o2, o3):
    if SOLID_PRINCIPAL_METHOD == 'og':
        omax_mid_min = [o1, o2, o3]
        omin = min(omax_mid_min)
        omax_mid_min.remove(omin)

        omax = max(omax_mid_min)
        omax_mid_min.remove(omax)

        omid = omax_mid_min[0]
    # elif SOLID_PRINCIPAL_METHOD == 'raw':
    #     omax = o1
    #     omid = o2
    #     omin = o3
    elif SOLID_PRINCIPAL_METHOD == '132':
        omax = o1
        omid = o3
        omin = o2
    else:  # pragma: no cover
        raise RuntimeError(SOLID_PRINCIPAL_METHOD)
    return omax, omid, omin
