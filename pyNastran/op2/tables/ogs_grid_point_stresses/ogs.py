from __future__ import annotations
from typing import TYPE_CHECKING
from struct import Struct
from numpy import frombuffer
from pyNastran.op2.op2_interface.op2_reader import mapfmt
from pyNastran.op2.tables.ogs_grid_point_stresses.ogs_surface_stresses import (
    GridPointSurfaceStressesArray,
    GridPointStressesVolumeDirectArray, GridPointStressesVolumePrincipalArray,
    GridPointStressesSurfaceDiscontinutiesArray,
    GridPointStressesVolumeDiscontinutiesArray,
    # strains
    GridPointSurfaceStrainsArray, GridPointStrainsVolumeDirectArray,
    GridPointStrainsVolumePrincipalArray,
    GridPointStrainsSurfaceDiscontinutiesArray
)
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.op2.op2 import OP2

class OGS:
    def __init__(self, op2: OP2):
        self.op2 = op2

    @property
    def size(self) -> int:
        return self.op2.size
    @property
    def factor(self) -> int:
        return self.op2.factor

    def _read_ogstr1_3(self, data: bytes, ndata: int):
        """OGSTR1 - grid point strains"""
        self._read_ogs1_3(data, ndata)

    def _read_ogs1_3(self, data: bytes, ndata: int):
        """OGS1 - grid point stresses"""
        op2 = self.op2
        op2.parse_approach_code(data)  # field 3
        op2.words = [
            'aCode', 'tCode', '???', 'isubcase',
            '???', '???', '???', 'dLoadID',
            'format_code', 'num_wide', 'o_code', '???',
            'acoustic_flag', '???', '???', '???',
            '???', '???', '???', '???',
            '???', '???', 'thermal', '???',
            '???', 'Title', 'subtitle', 'label']

        op2.parse_approach_code(data)
        #isubcase = self.get_values(data, b'i', 4)

        ## surface/volumeID
        op2.ogs = op2.add_data_parameter(data, 'ogs_id', b'i', 3, False)

        #: Reference coordinate system ID
        op2.refid = op2.add_data_parameter(data, 'refid', b'i', 8, False)

        ## format code
        op2.format_code = op2.add_data_parameter(data, 'format_code', b'i', 9, False)

        ## number of words per entry in record
        op2.num_wide = op2.add_data_parameter(data, 'num_wide', b'i', 10, False)

        ## Stress/Strain code
        op2.sCode = op2.add_data_parameter(data, 'sCode', b'i', 11, False)

        ## Output Coordinate System
        op2.oCoord = op2.add_data_parameter(data, 'oCoord', b'i', 12, False)

        ## Axis Specification code
        op2.axis = op2.add_data_parameter(data, 'axis', b'i', 13, False)

        #: Normal Specification Code
        op2.normal = op2.add_data_parameter(data, 'normal', b'i', 14, False)

        op2.fix_format_code()
        if not op2.is_sort1:
            raise NotImplementedError('OGS sort2...')

        ## assuming tCode=1
        if op2.analysis_code == 1:   # statics
            ## load set number
            op2.lsdvmn = op2.add_data_parameter(data, 'lsdvmn', b'i', 5, False)
            op2.data_names = op2.apply_data_code_value('data_names', ['lsdvmn'])
            op2.setNullNonlinearFactor()
        elif op2.analysis_code == 2:  # normal modes/buckling (real eigenvalues)
            ## mode number
            op2.mode = op2.add_data_parameter(data, 'mode', b'i', 5)
            ## real eigenvalue
            op2.eign = op2.add_data_parameter(data, 'eign', b'f', 6, False)
            op2.mode_cycle = 0.0
            op2._op2_readers.reader_oug.update_mode_cycle('mode_cycle')
            op2.data_names = op2.apply_data_code_value('data_names', ['mode', 'eign', 'mode_cycle'])
        #elif op2.analysis_code == 3: # differential stiffness
        #elif op2.analysis_code == 4: # differential stiffness
        #elif op2.analysis_code == 5: # frequency
        elif op2.analysis_code == 6:  # transient
            ## time step
            op2.time = op2.add_data_parameter(data, 'time', b'f', 5)
            op2.data_names = op2.apply_data_code_value('data_names', ['time'])
        #elif op2.analysis_code == 7:  # pre-buckling
        #elif op2.analysis_code == 8:  # post-buckling
        #elif op2.analysis_code == 9:  # complex eigenvalues
        elif op2.analysis_code == 10:  # nonlinear statics
            ## load step
            op2.lftsfq = op2.add_data_parameter(data, 'lftsfq', b'f', 5)
            op2.data_names = op2.apply_data_code_value('data_names', ['lftsfq'])
        #elif op2.analysis_code == 11:  # old geometric nonlinear statics
        #elif op2.analysis_code == 12:  # contran ? (may appear as aCode=6)  --> straight from DMAP...grrr...
        else:
            raise RuntimeError('invalid analysis_code...analysis_code=%s' % op2.analysis_code)

        #print "*isubcase=%s" % (op2.isubcase)
        #print "analysis_code=%s table_code=%s thermal=%s" %(op2.analysis_code,op2.table_code,self.thermal)

        #print op2.code_information()
        if op2.is_debug_file:
            op2.binary_debug.write('  approach_code  = %r\n' % op2.approach_code)
            op2.binary_debug.write('  tCode          = %r\n' % op2.tCode)
            op2.binary_debug.write('  isubcase       = %r\n' % op2.isubcase)
        op2._read_title(data)
        op2._write_debug_bits()

    def _read_ogstr1_4(self, data: bytes, ndata: int) -> int:
        """OGSTR1 - grid point strains"""
        return self._read_ogs1_4(data, ndata, restype='strains')

    def _read_ogs1_4(self, data: bytes, ndata: int, restype: str='stresses') -> int:
        """OGS1 - grid point stresses"""
        op2 = self.op2
        if op2.table_code == 26:
            # OGS1 - grid point stresses - surface
            assert op2.table_name in [b'OGS1', b'OGSTR1'], f'table_name={op2.table_name} table_code={op2.table_code}'
            n = self._read_ogs1_table26(data, ndata, restype)
        elif op2.table_code == 27:
            #OGS1 - grid point stresses - volume direct
            assert op2.table_name in [b'OGS1', b'OGSTR1', b'OGS1X'], f'table_name={op2.table_name} table_code={op2.table_code}'
            n = self._read_ogs1_table27(data, ndata, restype)
        elif op2.table_code == 28:
            #OGS1- grid point stresses - principal
            assert op2.table_name in [b'OGS1', b'OGSTR1'], f'table_name={op2.table_name} table_code={op2.table_code}'
            n = self._read_ogs1_table28(data, ndata, restype)
        elif op2.table_code == 35:
            # OGS - Grid point stress discontinuities (plane strain)
            assert op2.table_name in [b'OGS1', b'OGSTR1'], f'table_name={op2.table_name} table_code={op2.table_code}'
            n = self._read_ogs1_table35(data, ndata, restype)
        else:
            #msg = op2.code_information()
            raise RuntimeError(op2.code_information())
            #n = self._not_implemented_or_skip(data, ndata, msg)
        del op2.ogs
        return n

    def _read_ogs1_table28(self, data, ndata, restype: str):
        op2 = self.op2
        if op2.num_wide == 15:
            n = self._read_ogs1_table28_numwide15(data, ndata, restype)
        else:
            raise RuntimeError(op2.code_information())
        return n

    def _read_ogs1_table28_numwide15(self, data, ndata, restype: str):
        """
        TCODE =28 Volume with principal
        1 EKEY I 10*grid point identification number + device code
        2 LXA RS Direction cosine from x to a
        3 LXB RS Direction cosine from x to b
        4 LXC RS Direction cosine from x to c

        5 LYA RS Direction cosine from y to a
        6 LYB RS Direction cosine from y to b
        7 LYC RS Direction cosine from y to c

        8 LZA RS Direction cosine from z to a
        9 LZB RS Direction cosine from z to b
        10 LZC RS Direction cosine from z to c

        11 SA RS Principal in a
        12 SB RS Principal in b
        13 SC RS Principal in c
        14 EPR RS Mean pressure
        15 EHVM RS Hencky-von Mises or octahedral
        """
        op2 = self.op2
        result_name = f'grid_point_{restype}_volume_principal'
        if 'strain' in restype:
            obj_vector_real = GridPointStrainsVolumePrincipalArray
        else:
            obj_vector_real = GridPointStressesVolumePrincipalArray

        if op2._results.is_not_saved(result_name):
            op2.log.warning(f'skipping {result_name}')
            return ndata
        op2._results._found_result(result_name)
        slot = getattr(op2, result_name)
        n = 0

        #result_name, is_random = self._apply_oes_ato_crm_psd_rms_no(result_name)
        ntotal = 60 * self.factor # 15 * 4
        nelements = ndata // ntotal
        assert ndata % ntotal == 0
        auto_return, is_vectorized = op2._create_oes_object4(
            nelements, result_name, slot, obj_vector_real)
        if auto_return:
            return nelements * ntotal

        obj = op2.obj
        dt = op2.nonlinear_factor
        if op2.use_vector and is_vectorized and 0:
            n = nelements * ntotal
            #itotal = obj.ielement
            #ielement2 = obj.itotal + nelements
            #itotal2 = ielement2

            #floats = frombuffer(data, dtype=op2.fdtype).reshape(nelements, 11).copy()
            #obj._times[obj.itime] = dt
            #if obj.itime == 0:
                #ints = frombuffer(data, dtype=op2.idtype).reshape(nelements, 11).copy()
                #nids = ints[:, 0] // 10
                #eids = ints[:, 1]
                #assert nids.min() > 0, nids.min()
                #obj.node_element[itotal:itotal2, 0] = nids
                #obj.node_element[itotal:itotal2, 1] = eids

            ##[lxa, lxb, lxc, lya, lyb, lyc, lza, lzb, lzc, sa, sb, sc, epr, ovm]
            #strings = frombuffer(data, dtype=op2._uendian + 'S4').reshape(nelements, 11)[:, 2].copy()
            #obj.location[itotal:itotal2] = strings
            #obj.data[obj.itime, itotal:itotal2, :] = floats[:, 3:]#.copy()
            #obj.itotal = itotal2
            #obj.ielement = ielement2
            #n = ndata
        else:
            s = Struct(mapfmt(op2._endian + b'i14f', self.size))
            #nelements = ndata // 60  # 15*4
            for unused_i in range(nelements):
                edata = data[n:n+ntotal]
                out = s.unpack(edata)
                (eid_device, lxa, lxb, lxc, lya, lyb, lyc, lza, lzb, lzc, sa, sb, sc, epr, ovm) = out
                eid = eid_device // 10
                assert eid > 0, eid
                #op2.obj.add_sort1(dt, eid, lxa, lxb, lxc, lya, lyb, lyc, lza, lzb, lzc,
                                   #sa, sb, sc, epr, ovm)
                n += ntotal

        assert ndata > 0, ndata
        assert nelements > 0, f'nelements={nelements} element_type={op2.element_type} element_name={op2.element_name!r}'
        #assert ndata % ntotal == 0, '%s n=%s nwide=%s len=%s ntotal=%s' % (op2.element_name, ndata % ntotal, ndata % op2.num_wide, ndata, ntotal)
        assert op2.num_wide * 4 * self.factor == ntotal, 'numwide*4=%s ntotal=%s' % (op2.num_wide * 4, ntotal)
        assert n > 0, f'n = {n} result_name={result_name}'
        return n

    #-----------------------------------------------------------------------------------
    def _read_ogs1_table26(self, data: bytes, ndata: int, restype: str) -> int:
        """reads grid point stresses"""
        op2 = self.op2
        if op2.num_wide == 11:  # real/random
            n = self._read_ogs1_table26_numwide11(data, ndata, restype)
        else:
            msg = f'only num_wide=11 is allowed  num_wide={op2.num_wide}'
            raise NotImplementedError(msg)
        return n

    def _read_ogs1_table26_numwide11(self, data: bytes, ndata: int, restype: str) -> int:
        """surface stresses"""
        op2 = self.op2
        result_name = f'grid_point_surface_{restype}'
        if 'strain' in restype:
            obj_vector_real = GridPointSurfaceStrainsArray
        else:
            obj_vector_real = GridPointSurfaceStressesArray
        if op2._results.is_not_saved(result_name):
            op2.log.warning(f'skipping {result_name}')
            return ndata
        op2._results._found_result(result_name)
        slot = getattr(op2, result_name)
        n = 0

        #result_name, is_random = self._apply_oes_ato_crm_psd_rms_no(result_name)
        ntotal = 44 * self.factor # 4*11
        nelements = ndata // ntotal
        auto_return, is_vectorized = op2._create_oes_object4(
            nelements, result_name, slot, obj_vector_real)
        if auto_return:
            return nelements * ntotal

        obj = op2.obj
        dt = op2.nonlinear_factor

        if op2.use_vector and is_vectorized:
            n = nelements * ntotal
            itotal = obj.ielement
            ielement2 = obj.itotal + nelements
            itotal2 = ielement2

            floats = frombuffer(data, dtype=op2.fdtype8).reshape(nelements, 11).copy()
            obj._times[obj.itime] = dt
            if obj.itime == 0:
                ints = frombuffer(data, dtype=op2.idtype8).reshape(nelements, 11).copy()
                nids = ints[:, 0] // 10
                eids = ints[:, 1]
                assert nids.min() > 0, nids.min()
                obj.node_element[itotal:itotal2, 0] = nids
                obj.node_element[itotal:itotal2, 1] = eids

            #[fiber, nx, ny, txy, angle, major, minor, tmax, ovm]
            s4 = 'S%i' % self.size
            strings = frombuffer(data, dtype=op2._uendian + s4).reshape(nelements, 11)[:, 2].copy()
            obj.location[itotal:itotal2] = strings
            obj.data[obj.itime, itotal:itotal2, :] = floats[:, 3:]#.copy()
            obj.itotal = itotal2
            obj.ielement = ielement2
            n = ndata
        else:
            fmt = op2._endian + (b'2i4s8f' if self.size == 4 else b'2q8s8d')
            s = Struct(fmt)
            nelements = ndata // ntotal  # 11*4
            for unused_i in range(nelements):
                edata = data[n:n+ntotal]
                out = s.unpack(edata)
                (nid_device, eid, fiber, nx, ny, txy, angle, major, minor, tmax, ovm) = out
                nid = nid_device // 10
                fiber = fiber.decode('utf-8').strip()
                assert nid > 0, nid
                op2.obj.add_sort1(dt, nid, eid, fiber, nx, ny, txy,
                                   angle, major, minor, tmax, ovm)
                n += ntotal

        assert ndata > 0, ndata
        assert nelements > 0, 'nelements=%r element_type=%s element_name=%r' % (nelements, op2.element_type, op2.element_name)
        #assert ndata % ntotal == 0, '%s n=%s nwide=%s len=%s ntotal=%s' % (op2.element_name, ndata % ntotal, ndata % op2.num_wide, ndata, ntotal)
        #assert op2.num_wide * 4 * self.factor == ntotal, 'numwide*4=%s ntotal=%s' % (op2.num_wide * 4, ntotal)
        assert n > 0, f'n = {n} result_name={result_name}'
        return n

    def _read_ogs1_table27(self, data: bytes, ndata: int, restype: str) -> int:
        """OGS1 - grid point stresses - volume direct"""
        op2 = self.op2
        #is_sort1 = op2.is_sort1
        if op2.num_wide == 9:  # real/random
            #result_name = 'grid_point_stresses_volume_direct'
            n = self._read_ogs1_table27_numwide9(data, ndata, restype)
        else:
            msg = op2.code_information()
            #msg = 'only num_wide=9 is allowed  num_wide=%s' % op2.num_wide
            raise RuntimeError(msg)
        return n

    def _read_ogs1_table27_numwide9(self, data: bytes, ndata: int, restype: str) -> int:
        """
        TCODE =27 Volume with direct
        1 EKEY I 10*grid point identification number + Device Code
        2 NX RS Normal in x
        3 NY RS Normal in y
        4 NZ RS Normal in z
        5 TXY RS Shear in xy
        6 TYZ RS Shear in yz
        7 TZX RS Shear in zx
        8 PR RS Mean pressure
        9 HVM RS Hencky-von Mises or Octahedral
        """
        op2 = self.op2
        result_name = f'grid_point_{restype}_volume_direct'
        if op2._results.is_not_saved(result_name):
            op2.log.warning(f'skipping {result_name}')
            return ndata

        if 'strain' in restype:
            obj_vector_real = GridPointStrainsVolumeDirectArray
        else:
            obj_vector_real = GridPointStressesVolumeDirectArray
        op2._results._found_result(result_name)
        slot = getattr(op2, result_name)
        n = 0

        #result_name, is_random = self._apply_oes_ato_crm_psd_rms_no(result_name)
        ntotal = 36 * self.factor  # 9 * 4
        nelements = ndata // ntotal
        assert ndata % (nelements * ntotal) == 0, ndata % (nelements * ntotal)
        auto_return, is_vectorized = op2._create_oes_object4(
            nelements, result_name, slot, obj_vector_real)
        if auto_return:
            return nelements * ntotal

        obj = op2.obj
        dt = op2.nonlinear_factor

        if op2.use_vector and is_vectorized:
            n = nelements * ntotal
            itotal = obj.ielement
            ielement2 = obj.itotal + nelements
            itotal2 = ielement2

            floats = frombuffer(data, dtype=op2.fdtype8).reshape(nelements, 9)#.copy()
            obj._times[obj.itime] = dt
            if obj.itime == 0:
                ints = frombuffer(data, dtype=op2.idtype8).reshape(nelements, 9)
                nids = ints[:, 0] // 10
                assert nids.min() > 0, nids.min()
                obj.node[itotal:itotal2] = nids

            #[nid, nx, ny, nz, txy, tyz, txz, pressure, ovm]
            #strings = frombuffer(data, dtype=op2._uendian + 'S4').reshape(nelements, 11)[:, 2].copy()
            #obj.location[itotal:itotal2] = strings
            obj.data[obj.itime, itotal:itotal2, :] = floats[:, 1:]#.copy()
            obj.itotal = itotal2
            obj.ielement = ielement2
            n = ndata
        else:
            fmt = mapfmt(op2._endian + b'i8f', self.size)
            s = Struct(fmt)
            for unused_i in range(nelements):
                edata = data[n:n+ntotal]
                out = s.unpack(edata)
                (nid_device, nx, ny, nz, txy, tyz, txz, pressure, ovm) = out
                nid = nid_device // 10
                assert nid > 0, nid
                op2.obj.add_sort1(dt, nid, nx, ny, nz, txy, tyz, txz, pressure, ovm)
                n += ntotal
        return n


    def _read_ogs1_table35(self, data: bytes, ndata: int, restype: str) -> int:
        """
        grid point stress discontinuities (plane stress/strain)

        TCODE =35 Grid point stresses for surfaces with plane strain
        1 EKEY I 10*grid point identification number and grid code
        2 NX RS Normal in x
        3 NY RS Normal in y
        4 NZ RS Normal in z (always -1)
        5 TXY RS Shear in xy
        6 PR RS Mean pressure (always -1)
        """
        op2 = self.op2
        if restype in 'strains':
            result_name = 'grid_point_strain_discontinuities'
        else:
            result_name = 'grid_point_stress_discontinuities'
        if op2._results.is_not_saved(result_name):
            op2.log.warning(f'skipping {result_name}')
            return ndata
        op2._results._found_result(result_name)
        slot = getattr(op2, result_name)
        n = 0

        if op2.num_wide == 6:
            if 'strain' in restype:
                obj_vector_real = GridPointStrainsSurfaceDiscontinutiesArray
            else:
                obj_vector_real = GridPointStressesSurfaceDiscontinutiesArray

            #result_name, is_random = self._apply_oes_ato_crm_psd_rms_no(result_name)
            ntotal = 6 * 4 * self.factor
            nelements = ndata // ntotal
            assert ndata % (nelements * ntotal) == 0, ndata % (nelements * ntotal)
            auto_return, is_vectorized = op2._create_oes_object4(
                nelements, result_name, slot, obj_vector_real)
            if auto_return:
                return nelements * ntotal

            obj = op2.obj
            dt = op2.nonlinear_factor

            if op2.use_vector and is_vectorized:
                n = nelements * ntotal
                itotal = obj.ielement
                ielement2 = obj.itotal + nelements
                itotal2 = ielement2

                floats = frombuffer(data, dtype=op2.fdtype).reshape(nelements, 6)#.copy()
                obj._times[obj.itime] = dt
                if obj.itime == 0:
                    ints = frombuffer(data, dtype=op2.idtype).reshape(nelements, 6)
                    nids = ints[:, 0] // 10
                    assert nids.min() > 0, nids.min()
                    obj.node[itotal:itotal2] = nids

                #[nid, nx, ny, nz, txy, pressure]
                obj.data[obj.itime, itotal:itotal2, :] = floats[:, 1:]#.copy()
                obj.itotal = itotal2
                obj.ielement = ielement2
                n = ndata
            else:
                s = Struct(mapfmt(op2._endian + b'i5f', self.size))
                nelements = ndata // ntotal  # 6*4
                for unused_i in range(nelements):
                    out = s.unpack(data[n:n+ntotal])
                    (nid_device, nx, ny, nz, txy, pressure) = out
                    nid = nid_device // 10
                    assert nid > 0, nid
                    op2.obj.add_sort1(dt, nid, nx, ny, nz, txy, pressure)
                    n += ntotal
        else:
            msg = 'only num_wide=11 is allowed  num_wide=%s' % op2.num_wide
            raise RuntimeError(msg)
        return n
