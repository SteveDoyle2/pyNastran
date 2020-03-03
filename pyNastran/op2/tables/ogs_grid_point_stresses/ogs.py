from struct import Struct
from numpy import frombuffer
from pyNastran.op2.op2_interface.op2_common import OP2Common
from pyNastran.op2.op2_interface.op2_reader import mapfmt
from pyNastran.op2.tables.ogs_grid_point_stresses.ogs_surface_stresses import (
    GridPointSurfaceStressesArray,
    GridPointStressesVolumeDirectArray, GridPointStressesVolumePrincipalArray,
    GridPointStressesSurfaceDiscontinutiesArray,
    GridPointStressesVolumeDiscontinutiesArray)


class OGS(OP2Common):
    def __init__(self):
        OP2Common.__init__(self)

    def _read_ogs1_3(self, data, ndata):
        unused_three = self.parse_approach_code(data)
        self.words = [
            'aCode', 'tCode', '???', 'isubcase',
            '???', '???', '???', 'dLoadID',
            'format_code', 'num_wide', 'o_code', '???',
            'acoustic_flag', '???', '???', '???',
            '???', '???', '???', '???',
            '???', '???', 'thermal', '???',
            '???', 'Title', 'subtitle', 'label']

        self.parse_approach_code(data)
        #isubcase = self.get_values(data, b'i', 4)

        ## surface/volumeID
        self.ogs = self.add_data_parameter(data, 'ogs_id', b'i', 3, False)

        #: Reference coordinate system ID
        self.refid = self.add_data_parameter(data, 'refid', b'i', 8, False)

        ## format code
        self.format_code = self.add_data_parameter(data, 'format_code', b'i', 9, False)

        ## number of words per entry in record
        self.num_wide = self.add_data_parameter(data, 'num_wide', b'i', 10, False)

        ## Stress/Strain code
        self.sCode = self.add_data_parameter(data, 'sCode', b'i', 11, False)

        ## Output Coordinate System
        self.oCoord = self.add_data_parameter(data, 'oCoord', b'i', 12, False)

        ## Axis Specification code
        self.axis = self.add_data_parameter(data, 'axis', b'i', 13, False)

        #: Normal Specification Code
        self.normal = self.add_data_parameter(data, 'normal', b'i', 14, False)

        self.fix_format_code()
        if not self.is_sort1:
            raise NotImplementedError('OGS sort2...')

        ## assuming tCode=1
        if self.analysis_code == 1:   # statics
            ## load set number
            self.lsdvmn = self.add_data_parameter(data, 'lsdvmn', b'i', 5, False)
            self.data_names = self.apply_data_code_value('data_names', ['lsdvmn'])
            self.setNullNonlinearFactor()
        elif self.analysis_code == 2:  # normal modes/buckling (real eigenvalues)
            ## mode number
            self.mode = self.add_data_parameter(data, 'mode', b'i', 5)
            ## real eigenvalue
            self.eign = self.add_data_parameter(data, 'eign', b'f', 6, False)
            self.mode_cycle = 0.0
            self.update_mode_cycle('mode_cycle')
            self.data_names = self.apply_data_code_value('data_names', ['mode', 'eign', 'mode_cycle'])
        #elif self.analysis_code == 3: # differential stiffness
        #elif self.analysis_code == 4: # differential stiffness
        #elif self.analysis_code == 5: # frequency
        elif self.analysis_code == 6:  # transient
            ## time step
            self.time = self.add_data_parameter(data, 'time', b'f', 5)
            self.data_names = self.apply_data_code_value('data_names', ['time'])
        #elif self.analysis_code == 7:  # pre-buckling
        #elif self.analysis_code == 8:  # post-buckling
        #elif self.analysis_code == 9:  # complex eigenvalues
        elif self.analysis_code == 10:  # nonlinear statics
            ## load step
            self.lftsfq = self.add_data_parameter(data, 'lftsfq', b'f', 5)
            self.data_names = self.apply_data_code_value('data_names', ['lftsfq'])
        #elif self.analysis_code == 11:  # old geometric nonlinear statics
        #elif self.analysis_code == 12:  # contran ? (may appear as aCode=6)  --> straight from DMAP...grrr...
        else:
            raise RuntimeError('invalid analysis_code...analysis_code=%s' % self.analysis_code)

        #print "*isubcase=%s" % (self.isubcase)
        #print "analysis_code=%s table_code=%s thermal=%s" %(self.analysis_code,self.table_code,self.thermal)

        #print self.code_information()
        if self.is_debug_file:
            self.binary_debug.write('  approach_code  = %r\n' % self.approach_code)
            self.binary_debug.write('  tCode          = %r\n' % self.tCode)
            self.binary_debug.write('  isubcase       = %r\n' % self.isubcase)
        self._read_title(data)
        self._write_debug_bits()

    def _read_ogs1_4(self, data, ndata):
        if self.table_code == 26:
            # OGS1 - grid point stresses - surface
            assert self.table_name in [b'OGS1'], 'table_name=%s table_code=%s' % (self.table_name, self.table_code)
            n = self._read_ogs1_table26(data, ndata)
        elif self.table_code == 27:
            #OGS1 - grid point stresses - volume direct
            assert self.table_name in [b'OGS1'], 'table_name=%s table_code=%s' % (self.table_name, self.table_code)
            n = self._read_ogs1_table27(data, ndata)
        elif self.table_code == 28:
            #OGS1- grid point stresses - principal
            assert self.table_name in [b'OGS1'], 'table_name=%s table_code=%s' % (self.table_name, self.table_code)
            n = self._read_ogs1_table28(data, ndata)
        elif self.table_code == 35:
            # OGS - Grid point stress discontinuities (plane strain)
            assert self.table_name in [b'OGS1'], 'table_name=%s table_code=%s' % (self.table_name, self.table_code)
            n = self._read_ogs1_table35(data, ndata)
        else:
            msg = self.code_information()
            raise RuntimeError(self.code_information())
            #n = self._not_implemented_or_skip(data, ndata, msg)
        del self.ogs
        return n

    def _read_ogs1_table28(self, data, ndata):
        if self.num_wide == 15:
            n = self._read_ogs1_table28_numwide15(data, ndata)
        else:
            raise RuntimeError(self.code_information())
        return n

    def _read_ogs1_table28_numwide15(self, data, ndata):
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
        result_name = 'grid_point_stresses_volume_principal'
        obj_vector_real = GridPointStressesVolumePrincipalArray
        if self._results.is_not_saved(result_name):
            return ndata
        self._results._found_result(result_name)
        slot = getattr(self, result_name)
        n = 0

        #result_name, is_random = self._apply_oes_ato_crm_psd_rms_no(result_name)
        ntotal = 60 * self.factor # 15 * 4
        nelements = ndata // ntotal
        assert ndata % ntotal == 0
        auto_return, is_vectorized = self._create_oes_object4(
            nelements, result_name, slot, obj_vector_real)
        if auto_return:
            return nelements * ntotal

        obj = self.obj
        dt = self.nonlinear_factor
        if self.use_vector and is_vectorized and 0:
            n = nelements * ntotal
            #itotal = obj.ielement
            #ielement2 = obj.itotal + nelements
            #itotal2 = ielement2

            #floats = frombuffer(data, dtype=self.fdtype).reshape(nelements, 11).copy()
            #obj._times[obj.itime] = dt
            #if obj.itime == 0:
                #ints = frombuffer(data, dtype=self.idtype).reshape(nelements, 11).copy()
                #nids = ints[:, 0] // 10
                #eids = ints[:, 1]
                #assert nids.min() > 0, nids.min()
                #obj.node_element[itotal:itotal2, 0] = nids
                #obj.node_element[itotal:itotal2, 1] = eids

            ##[lxa, lxb, lxc, lya, lyb, lyc, lza, lzb, lzc, sa, sb, sc, epr, ovm]
            #strings = frombuffer(data, dtype=self._uendian + 'S4').reshape(nelements, 11)[:, 2].copy()
            #obj.location[itotal:itotal2] = strings
            #obj.data[obj.itime, itotal:itotal2, :] = floats[:, 3:]#.copy()
            #obj.itotal = itotal2
            #obj.ielement = ielement2
            #n = ndata
        else:
            s = Struct(mapfmt(self._endian + b'i14f', self.size))
            #nelements = ndata // 60  # 15*4
            for unused_i in range(nelements):
                edata = data[n:n+ntotal]
                out = s.unpack(edata)
                (eid_device, lxa, lxb, lxc, lya, lyb, lyc, lza, lzb, lzc, sa, sb, sc, epr, ovm) = out
                eid = eid_device // 10
                assert eid > 0, eid
                #self.obj.add_sort1(dt, eid, lxa, lxb, lxc, lya, lyb, lyc, lza, lzb, lzc,
                                   #sa, sb, sc, epr, ovm)
                n += ntotal

        assert ndata > 0, ndata
        assert nelements > 0, 'nelements=%r element_type=%s element_name=%r' % (nelements, self.element_type, self.element_name)
        #assert ndata % ntotal == 0, '%s n=%s nwide=%s len=%s ntotal=%s' % (self.element_name, ndata % ntotal, ndata % self.num_wide, ndata, ntotal)
        assert self.num_wide * 4 * self.factor == ntotal, 'numwide*4=%s ntotal=%s' % (self.num_wide * 4, ntotal)
        assert n > 0, "n = %s result_name=%s" % (n, result_name)
        return n
    #-----------------------------------------------------------------------------------
    def _read_ogs1_table26(self, data, ndata):
        """reads grid point stresses"""
        if self.num_wide == 11:  # real/random
            n = self._read_ogs1_table26_numwide11(data, ndata)
        else:
            msg = 'only num_wide=11 is allowed  num_wide=%s' % self.num_wide
            raise NotImplementedError(msg)
        return n

    def _read_ogs1_table26_numwide11(self, data, ndata):
        """surface stresses"""
        result_name = 'grid_point_surface_stresses'
        obj_vector_real = GridPointSurfaceStressesArray
        if self._results.is_not_saved(result_name):
            return ndata
        self._results._found_result(result_name)
        slot = getattr(self, result_name)
        n = 0

        #result_name, is_random = self._apply_oes_ato_crm_psd_rms_no(result_name)
        ntotal = 44 * self.factor # 4*11
        nelements = ndata // ntotal
        auto_return, is_vectorized = self._create_oes_object4(
            nelements, result_name, slot, obj_vector_real)
        if auto_return:
            return nelements * ntotal

        obj = self.obj
        dt = self.nonlinear_factor

        if self.use_vector and is_vectorized:
            n = nelements * ntotal
            itotal = obj.ielement
            ielement2 = obj.itotal + nelements
            itotal2 = ielement2

            floats = frombuffer(data, dtype=self.fdtype8).reshape(nelements, 11).copy()
            obj._times[obj.itime] = dt
            if obj.itime == 0:
                ints = frombuffer(data, dtype=self.idtype8).reshape(nelements, 11).copy()
                nids = ints[:, 0] // 10
                eids = ints[:, 1]
                assert nids.min() > 0, nids.min()
                obj.node_element[itotal:itotal2, 0] = nids
                obj.node_element[itotal:itotal2, 1] = eids

            #[fiber, nx, ny, txy, angle, major, minor, tmax, ovm]
            s4 = 'S%i' % self.size
            strings = frombuffer(data, dtype=self._uendian + s4).reshape(nelements, 11)[:, 2].copy()
            obj.location[itotal:itotal2] = strings
            obj.data[obj.itime, itotal:itotal2, :] = floats[:, 3:]#.copy()
            obj.itotal = itotal2
            obj.ielement = ielement2
            n = ndata
        else:
            fmt = self._endian + (b'2i4s8f' if self.size == 4 else b'2q8s8d')
            s = Struct(fmt)
            nelements = ndata // ntotal  # 11*4
            for unused_i in range(nelements):
                edata = data[n:n+ntotal]
                out = s.unpack(edata)
                (nid_device, eid, fiber, nx, ny, txy, angle, major, minor, tmax, ovm) = out
                nid = nid_device // 10
                fiber = fiber.decode('utf-8').strip()
                assert nid > 0, nid
                self.obj.add_sort1(dt, nid, eid, fiber, nx, ny, txy,
                                   angle, major, minor, tmax, ovm)
                n += ntotal

        assert ndata > 0, ndata
        assert nelements > 0, 'nelements=%r element_type=%s element_name=%r' % (nelements, self.element_type, self.element_name)
        #assert ndata % ntotal == 0, '%s n=%s nwide=%s len=%s ntotal=%s' % (self.element_name, ndata % ntotal, ndata % self.num_wide, ndata, ntotal)
        #assert self.num_wide * 4 * self.factor == ntotal, 'numwide*4=%s ntotal=%s' % (self.num_wide * 4, ntotal)
        assert n > 0, "n = %s result_name=%s" % (n, result_name)
        return n

    def _read_ogs1_table27(self, data, ndata):
        """OGS1 - grid point stresses - volume direct"""
        #is_sort1 = self.is_sort1
        if self.num_wide == 9:  # real/random
            #result_name = 'grid_point_stresses_volume_direct'
            n = self._read_ogs1_table27_numwide9(data, ndata)
        else:
            msg = self.code_information()
            #msg = 'only num_wide=9 is allowed  num_wide=%s' % self.num_wide
            raise RuntimeError(msg)
        return n

    def _read_ogs1_table27_numwide9(self, data, ndata):
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
        result_name = 'grid_point_stresses_volume_direct'
        if self._results.is_not_saved(result_name):
            return ndata

        obj_vector_real = GridPointStressesVolumeDirectArray
        self._results._found_result(result_name)
        slot = getattr(self, result_name)
        n = 0

        #result_name, is_random = self._apply_oes_ato_crm_psd_rms_no(result_name)
        ntotal = 36 * self.factor  # 9 * 4
        nelements = ndata // ntotal
        assert ndata % (nelements * ntotal) == 0, ndata % (nelements * ntotal)
        auto_return, is_vectorized = self._create_oes_object4(
            nelements, result_name, slot, obj_vector_real)
        if auto_return:
            return nelements * ntotal

        obj = self.obj
        dt = self.nonlinear_factor

        if self.use_vector and is_vectorized:
            n = nelements * ntotal
            itotal = obj.ielement
            ielement2 = obj.itotal + nelements
            itotal2 = ielement2

            floats = frombuffer(data, dtype=self.fdtype8).reshape(nelements, 9)#.copy()
            obj._times[obj.itime] = dt
            if obj.itime == 0:
                ints = frombuffer(data, dtype=self.idtype8).reshape(nelements, 9)
                nids = ints[:, 0] // 10
                assert nids.min() > 0, nids.min()
                obj.node[itotal:itotal2] = nids

            #[nid, nx, ny, nz, txy, tyz, txz, pressure, ovm]
            #strings = frombuffer(data, dtype=self._uendian + 'S4').reshape(nelements, 11)[:, 2].copy()
            #obj.location[itotal:itotal2] = strings
            obj.data[obj.itime, itotal:itotal2, :] = floats[:, 1:]#.copy()
            obj.itotal = itotal2
            obj.ielement = ielement2
            n = ndata
        else:
            fmt = mapfmt(self._endian + b'i8f', self.size)
            s = Struct(fmt)
            for unused_i in range(nelements):
                edata = data[n:n+ntotal]
                out = s.unpack(edata)
                (nid_device, nx, ny, nz, txy, tyz, txz, pressure, ovm) = out
                nid = nid_device // 10
                assert nid > 0, nid
                self.obj.add_sort1(dt, nid, nx, ny, nz, txy, tyz, txz, pressure, ovm)
                n += ntotal
        return n


    def _read_ogs1_table35(self, data, ndata):
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
        result_name = 'grid_point_stress_discontinuities'
        if self._results.is_not_saved(result_name):
            return ndata
        self._results._found_result(result_name)
        slot = getattr(self, result_name)
        n = 0

        if self.num_wide == 6:
            obj_vector_real = GridPointStressesSurfaceDiscontinutiesArray

            #result_name, is_random = self._apply_oes_ato_crm_psd_rms_no(result_name)
            ntotal = 6 * 4 * self.factor
            nelements = ndata // ntotal
            assert ndata % (nelements * ntotal) == 0, ndata % (nelements * ntotal)
            auto_return, is_vectorized = self._create_oes_object4(
                nelements, result_name, slot, obj_vector_real)
            if auto_return:
                return nelements * ntotal

            obj = self.obj
            dt = self.nonlinear_factor

            if self.use_vector and is_vectorized:
                n = nelements * ntotal
                itotal = obj.ielement
                ielement2 = obj.itotal + nelements
                itotal2 = ielement2

                floats = frombuffer(data, dtype=self.fdtype).reshape(nelements, 6)#.copy()
                obj._times[obj.itime] = dt
                if obj.itime == 0:
                    ints = frombuffer(data, dtype=self.idtype).reshape(nelements, 6)
                    nids = ints[:, 0] // 10
                    assert nids.min() > 0, nids.min()
                    obj.node[itotal:itotal2] = nids

                #[nid, nx, ny, nz, txy, pressure]
                obj.data[obj.itime, itotal:itotal2, :] = floats[:, 1:]#.copy()
                obj.itotal = itotal2
                obj.ielement = ielement2
                n = ndata
            else:
                s = Struct(mapfmt(self._endian + b'i5f', self.size))
                nelements = ndata // ntotal  # 6*4
                for unused_i in range(nelements):
                    out = s.unpack(data[n:n+ntotal])
                    (nid_device, nx, ny, nz, txy, pressure) = out
                    nid = nid_device // 10
                    assert nid > 0, nid
                    self.obj.add_sort1(dt, nid, nx, ny, nz, txy, pressure)
                    n += ntotal
        else:
            msg = 'only num_wide=11 is allowed  num_wide=%s' % self.num_wide
            raise RuntimeError(msg)
        return n
