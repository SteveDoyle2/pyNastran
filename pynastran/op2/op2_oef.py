#pylint: disable=C0301
"""
Defines the Real/Complex Forces created by:
    FORCE = ALL
"""
from __future__ import print_function
from six import b
from six.moves import range
from struct import Struct
import numpy as np
from numpy import fromstring, vstack, sin, cos, radians, array
from numpy import hstack, zeros

from pyNastran.op2.op2_helper import polar_to_real_imag
from pyNastran.op2.op2_common import OP2Common

from pyNastran.op2.tables.oef_forces.oef_thermal_objects import (
    Real1DHeatFluxArray,
    HeatFlux_2D_3DArray,
    RealChbdyHeatFluxArray,
    RealConvHeatFluxArray,
    RealHeatFluxVUArray,

    RealHeatFluxVUBeamArray,
    RealHeatFluxVU3DArray,

    # TODO: vectorize 2
    #HeatFlux_VUBEAM,
    #HeatFlux_VU_3D,
)
from pyNastran.op2.tables.oef_forces.oef_force_objects import (
    RealRodForceArray, RealViscForceArray,
    RealCBarForceArray, RealCBar100ForceArray,
    RealCBushForceArray,
    RealPlateForceArray,
    RealPlateBilinearForceArray,
    RealSpringForceArray, RealDamperForceArray,
    RealCShearForceArray,
    RealCGapForceArray,
    RealConeAxForceArray,
    RealSolidPressureForceArray,
    RealCBeamForceArray,
    RealBendForceArray,

    # TODO: vectorize 2
    #RealForce_VU_2D, RealForce_VU,
)
from pyNastran.op2.tables.oef_forces.oef_complex_force_objects import (
    ComplexRodForceArray,
    ComplexCBarForceArray,
    ComplexCBeamForceArray,
    ComplexCBushForceArray,
    ComplexCShearForceArray,
    ComplexSpringForceArray,
    ComplexDamperForceArray,
    ComplexViscForceArray,
    ComplexPlateForceArray,
    ComplexPlate2ForceArray,
    ComplexSolidPressureForceArray,
    ComplexCBendForceArray,

    # TODO: vectorize 2
    #ComplexForce_VU_2D, ComplexForce_VU,
)


class OEF(OP2Common):
    """Defines OEFx table reading for element forces/heat flux"""
    def __init__(self):
        OP2Common.__init__(self)

    def _oef_force_code(self):
        """
        Gets the numwide codes for the element to determine if
        the real or complex result should be found.
        The format and sort codes do not always give the right answer...
        """
        real_mapper = {
            1: 3,    # CROD
            2: 1 + (10 - 1) * 11,  # CBEAM
            3: 3,    # CTUBE
            4: 17,   # CSHEAR
            10: 3,    # CONROD
            11: 2,    # CELAS1
            12: 2,    # CELAS2
            13: 2,    # CELAS3
            14: 2,    # CELAS4

            20: 2,    # CDAMP1
            21: 2,    # CDAMP2
            22: 2,    # CDAMP3
            23: 2,    # CDAMP4
            24: 3,    # CVISC
            33: 9,    # CQUAD4
            34: 9,    # CBAR
            35: 7,    # CCONEAX
            38: 9,    # CGAP
            40: 8,    # CBUSH1D ???
            64: 2 + (11 - 2) * 5,  # CQUAD8
            69: 1 + (8 - 1) * 2,  # CBEND
            70: 2 + (11 - 2) * 4,  # CTRIAR
            74: 9,    # CTRIA3
            75: 2 + (11 - 2) * 4,  # CTRIA6


            #76:  16,   # Acoustic Velocity/Pressure CHEXA ???
            76: None,  # dummy so it doesnt go into the real results
            77: 10,   # Acoustic Velocity/Pressure CPENTA
            78: 10,   # Acoustic Velocity/Pressure CTETRA

            82: 2 + (11 - 2) * 5,  # CQUADR
            95: 9,    # composite CQUAD4 ???
            96: 9,    # composite CQUAD8 ???
            97: 9,    # composite CTRIA3 ???
            98: 9,    # composite CTRIA6 ???
            100: 8,    # BARS
            102: 7,    # CBUSH
            144: 2 + (11 - 2) * 5,  # bilinear CQUAD4
            189: 6 + (19 - 6) * 4,  # VUQUAD
            190: 6 + (19 - 6) * 3,  # VUTRIA
            191: 4 + (12 - 4) * 2,  # VUBEAM
            200: 9,    # CWELD
            232: 9,    # composite CQUADR ???
            233: 9,    # composite TRIAR ???
            235: 9,    # punch CQUADR...num_wide in DMAP is wrong...left out first entry...
            236: 8,    # punch CTRIAR
        }
        imag_mapper = {
            1: 5,    # CROD
            2: 1 + (17 - 1) * 11,  # CBEAM
            3: 5,    # CTUBE
            4: 33,   # CSHEAR
            10: 5,    # CONROD

            11: 3,    # CELAS1
            12: 3,    # CELAS2
            13: 3,    # CELAS3
            14: 3,    # CELAS4

            20: 3,    # CDAMP1
            21: 3,    # CDAMP2
            22: 3,    # CDAMP3
            23: 3,    # CDAMP4
            24: 5,    # CVISC
            33: 17,   # CQUAD4
            34: 17,   # CBAR
            35: 7,    # CCONEAX # needed to not crash the code...
            38: 9,    # CGAP
            40: 8,    # CBUSH1D ???
            64: 2 + (19 - 2) * 5,  # CQUAD8
            69: 1 + (14 - 1) * 2,  # CBEND
            70: 2 + (19 - 2) * 4,  # CTRIAR
            74: 17,   # CTRIA3
            75: 2 + (19 - 2) * 4,  # CTRIA6

            76: 16,   # Acoustic Velocity/Pressure CHEXA_PR
            77: 16,   # Acoustic Velocity/Pressure CPENTA_PR
            78: 16,   # Acoustic Velocity/Pressure CTETRA_PR

            82: 2 + (19 - 2) * 5,  # CQUADR
            95: 9,    # composite CQUAD4 ???
            96: 9,    # composite CQUAD8 ???
            97: 9,    # composite CTRIA3 ???
            98: 9,    # composite CTRIA6 ???
            100: 14,   # BARS
            102: 13,   # CBUSH
            144: 2 + (19 - 2) * 5,  # bilinear CQUAD4
            189: 6 + (31 - 6) * 4,  # VUQUAD
            190: 6 + (31 - 6) * 3,  # VUTRIA
            191: 4 + (18 - 4) * 2,  # VUBEAM
            200: 17,   # CWELD
            232: 9,    # composite CQUADR ???
            233: 9,    # composite TRIAR ???
            235: 17,   # punch CQUADR...num_wide in DMAP is wrong...left out first entry...
            236: 16,   # punch CTRIAR
        }
        try:
            real = real_mapper[self.element_type]
        except KeyError:
            real = None

        try:
            imag = imag_mapper[self.element_type]
        except KeyError:
            imag = None
        return (real, imag)

    def _read_oef1_3(self, data, ndata):
        """Table 3 parser for OEF1 table"""
        self._data_factor = 1
        self.words = [
            'aCode', 'tCode', 'element_type', 'isubcase',
            '???', '???', '???', '???',
            'format_code', 'num_wide', 'o_code', '???',
            '???', '???', '???', '???',
            '???', '???', '???', '???',
            '???', '???', '???', '???',
            '???', 'Title', 'subtitle', 'label']

        self.parse_approach_code(data)

        #: element type
        self.element_type = self.add_data_parameter(data, 'element_type', 'i', 3, False)

        # dynamic load set ID/random code
        #self.dLoadID = self.add_data_parameter(data, 'dLoadID', 'i', 8, False)

        #: format code
        self.format_code = self.add_data_parameter(data, 'format_code', 'i', 9, False)

        #: number of words per entry in record
        #: .. note: is this needed for this table ???
        self.num_wide = self.add_data_parameter(data, 'num_wide', 'i', 10, False)

        #: undefined in DMAP...
        self.o_code = self.add_data_parameter(data, 'o_code', 'i', 11, False)

        #: thermal flag; 1 for heat ransfer, 0 otherwise
        self.thermal = self.add_data_parameter(data, 'thermal', 'i', 23, False)

        ## assuming tCode=1
        if self.analysis_code == 1:   # statics
            self.loadID = self.add_data_parameter(data, 'loadID', 'i', 5, False)  # load set ID number
            self.data_names = self.apply_data_code_value('data_names', ['loadID'])
            self.setNullNonlinearFactor()
        elif self.analysis_code == 2:  # normal modes/buckling (real eigenvalues)
            #: mode number
            self.mode = self.add_data_parameter(data, 'mode', 'i', 5)
            #: eigenvalue
            self.eigr = self.add_data_parameter(data, 'eigr', 'f', 6, False)
            self.cycle = 0.
            self.update_mode_cycle('cycle')
            self.data_names = self.apply_data_code_value('data_names', ['mode', 'eigr', 'cycle'])
            # TODO: mode_cycle is not defined?
            #self.data_names = self.apply_data_code_value('data_names', ['mode', 'eigr', 'mode_cycle'])
        elif self.analysis_code == 3:  # differential stiffness 0
            #: load set ID number
            self.loadID = self.add_data_parameter(data, 'loadID', 'i', 5)
            self.data_names = self.apply_data_code_value('data_names', ['loadID'])
        elif self.analysis_code == 4:  # differential stiffness 1
            #: load set ID number
            self.loadID = self.add_data_parameter(data, 'loadID', 'i', 5)
            self.data_names = self.apply_data_code_value('data_names', ['loadID'])
        elif self.analysis_code == 5:   # frequency
            self.freq = self.add_data_parameter(data, 'freq', 'f', 5)  # frequency
            self.data_names = self.apply_data_code_value('data_names', ['freq'])
        elif self.analysis_code == 6:  # transient
            self.time = self.add_data_parameter(data, 'time', 'f', 5)  # time step
            self.data_names = self.apply_data_code_value('data_names', ['time'])
        elif self.analysis_code == 7:  # pre-buckling
            #: load set ID number
            self.loadID = self.add_data_parameter(data, 'loadID', 'i', 5)
            #self.apply_data_code_value('data_names',['lsdvmn'])
            self.data_names = self.apply_data_code_value('data_names', ['loadID'])
        elif self.analysis_code == 8:  # post-buckling
            #: load set ID number
            self.loadID = self.add_data_parameter(data, 'loadID', 'i', 5)
            #: real eigenvalue
            self.eigr = self.add_data_parameter(data, 'eigr', 'f', 6, False)
            self.data_names = self.apply_data_code_value('data_names', ['loadID', 'eigr'])
        elif self.analysis_code == 9:  # complex eigenvalues
            #: mode number
            self.mode = self.add_data_parameter(data, 'mode', 'i', 5)
            #: real eigenvalue
            self.eigr = self.add_data_parameter(data, 'eigr', 'f', 6, False)
            #: imaginary eigenvalue
            self.eigi = self.add_data_parameter(data, 'eigi', 'f', 7, False)
            self.data_names = self.apply_data_code_value('data_names', ['mode', 'eigr', 'eigi'])
        elif self.analysis_code == 10:  # nonlinear statics
            #: load step
            self.load_step = self.add_data_parameter(data, 'load_step', 'f', 5)
            self.data_names = self.apply_data_code_value('data_names', ['load_step'])
        elif self.analysis_code == 11:  # geometric nonlinear statics
            #: load set ID number
            self.loadID = self.add_data_parameter(data, 'loadID', 'i', 5)
            self.data_names = self.apply_data_code_value('data_names', ['loadID'])
        else:
            msg = 'invalid analysis_code...analysis_code=%s' % str(self.analysis_code)
            raise RuntimeError(msg)


        self.fix_format_code()
        self._parse_thermal_code()
        try:
            self.element_name = self.element_mapper[self.element_type]
        except KeyError:
            self.log.error(self.code_information())
            raise
        assert self.element_name != '', self.code_information()

        #self.element_name = self.element_mapper[self.element_type]
        self.data_code['element_name'] = self.element_name

        if self.is_debug_file:
            self.binary_debug.write('  %-14s = %r\n' % ('element_name', self.element_name))
            self.binary_debug.write('  %-14s = %r %s\n' % ('approach_code', self.approach_code,
                                                           self.approach_code_str(self.approach_code)))
            self.binary_debug.write('  %-14s = %r\n' % ('tCode', self.tCode))
            self.binary_debug.write('  %-14s = %r\n' % ('isubcase', self.isubcase))


        self._read_title(data)
        if self.element_type not in self.element_mapper:
            msg = 'element_type = %s' % self.element_type
            return self._not_implemented_or_skip(data, ndata, msg)
        self._write_debug_bits()
        assert self.num_wide != 146, self.code_information()

    def _read_oef2_3(self, data, ndata):
        """Table 3 parser for OEF2 table"""
        pass

    def _read_oef1_4(self, data, ndata):
        """Table 4 parser for OEF1 table"""
        if self.thermal == 0:
            if self.isubcase not in self.case_control_deck.subcases:
                self.subcase = self.case_control_deck.create_new_subcase(self.isubcase)
            self.subcase.add_op2_data(self.data_code, 'FORCE', self.log)
            n = self._read_oef1_loads(data, ndata)
        elif self.thermal == 1:
            n = self._read_oef1_thermal(data, ndata)
        elif self.thermal == 8: # NRL
            n = self._read_oef1_loads(data, ndata)
        else:
            n = self._not_implemented_or_skip(data, ndata, 'thermal=%s' % self.thermal)
        return n

    def _read_oef1_thermal(self, data, ndata):
        """Table 4 parser for OEF1 thermal table"""
        if self._results.is_not_saved('element_forces'):
            return ndata
        n = 0
        #is_magnitude_phase = self.is_magnitude_phase()
        dt = self.nonlinear_factor

        flag = 'element_id'
        if self.element_type in [1, 2, 3, 10, 34, 69]:  # ROD,BEAM,TUBE,CONROD,BAR,BEND
            # 1-CROD
            # 2-CBEAM
            # 3-CTUBE
            # 10-CONROD
            # 34-CBAR
            # 69-CBEND
            obj_vector_real = Real1DHeatFluxArray
            #if self.element_type == 1: # CROD
            if self.element_type == 1:
                result_name = 'crod_thermal_load'
            elif self.element_type == 2:
                result_name = 'cbeam_thermal_load'
            elif self.element_type == 3:
                result_name = 'ctube_thermal_load'
            elif self.element_type == 10:
                result_name = 'conrod_thermal_load'
            elif self.element_type == 34:
                result_name = 'cbar_thermal_load'
            elif self.element_type == 69:
                result_name = 'cbend_thermal_load'
            else:
                raise NotImplementedError('element_type=%s element_name=%s' % (
                    self.element_type, self.element_name))

            if self._results.is_not_saved(result_name):
                return ndata
            self._results._found_result(result_name)
            slot = getattr(self, result_name)

            if self.format_code == 1 and self.num_wide == 9:  # real
                ntotal = 36
                nelements = ndata // ntotal
                auto_return, is_vectorized = self._create_oes_object4(
                    nelements, result_name, slot, obj_vector_real)
                if auto_return:
                    return nelements * self.num_wide * 4
                obj = self.obj
                #if self.is_debug_file:
                    #self.binary_debug.write('  [cap, element1, element2, ..., cap]\n')
                    #self.binary_debug.write('  cap = %i  # assume 1 cap when there could have been multiple\n' % ndata)
                    #self.binary_debug.write('  #elementi = [eid_device, axial, torque]\n')
                    #self.binary_debug.write('  nelements=%i; nnodes=1 # centroid\n' % nelements)

                if self.use_vector and is_vectorized:
                    n = nelements * 4 * self.num_wide
                    itotal = obj.ielement
                    ielement2 = obj.itotal + nelements
                    itotal2 = ielement2

                    floats = fromstring(data, dtype=self.fdtype).reshape(nelements, 9)
                    obj._times[obj.itime] = dt

                    strings = fromstring(data, dtype=self._endian + 'S4').reshape(nelements, 9)
                    s = array([s1+s2 for s1, s2 in zip(strings[:, 1], strings[:, 2])])
                    #print(s)
                    #print('itime = ', obj.itime)
                    #print('---------')
                    if obj.itime == 0:
                        ints = fromstring(data, dtype=self.idtype).reshape(nelements, 9)
                        eids = ints[:, 0] // 10
                        assert eids.min() > 0, eids.min()
                        obj.element[itotal:itotal2] = eids
                        obj.element_data_type[itotal:itotal2] = s
                        #obj.element_type[obj.itime, itotal:itotal2, :] = strings[:, 3:]

                    #[etype, xgrad, ygrad, zgrad, xflux, yflux, zflux]
                    obj.data[obj.itime, itotal:itotal2, :] = floats[:, 3:]
                    obj.itotal = itotal2
                    obj.ielement = ielement2
                else:
                    s = Struct(b(self._endian + 'i8s6f'))
                    for i in range(nelements):
                        edata = data[n:n+ntotal]
                        out = s.unpack(edata)
                        (eid_device, eType, xgrad, ygrad, zgrad, xflux, yflux, zflux) = out
                        eid = eid_device // 10
                        obj.add(dt, eid, eType, xgrad, ygrad, zgrad, xflux, yflux, zflux)
                        n += ntotal
            else:
                msg = self.code_information()
                return self._not_implemented_or_skip(data, ndata, msg)

        elif self.element_type in [33, 53, 64, 74, 75,  # CQUAD4, CTRIAX6, CQUAD8, CTRIA3, CTRIA6
                                   39, 67, 68]:  # TETRA, HEXA, PENTA
            # 33-QUAD4-centroidal
            # 53-TRIAX6
            # 64-QUAD8
            # 74-TRIA3
            # 75-TRIA6

            # 39-TETRA
            # 67-HEXA
            # 68-PENTA
            if self.element_type == 33:
                result_name = 'cquad4_thermal_load'
            elif self.element_type == 53:
                result_name = 'ctriax6_thermal_load'
            elif self.element_type == 64:
                result_name = 'cquad8_thermal_load'
            elif self.element_type == 74:
                result_name = 'ctria3_thermal_load'
            elif self.element_type == 75:
                result_name = 'ctria6_thermal_load'
            elif self.element_type == 39:
                result_name = 'ctetra_thermal_load'
            elif self.element_type == 67:
                result_name = 'cthexa_thermal_load'
            elif self.element_type == 68:
                result_name = 'cpenta_thermal_load'
            else:
                raise NotImplementedError('element_type=%s element_name=%s' % (
                    self.element_type, self.element_name))

            obj_vector_real = HeatFlux_2D_3DArray
            #if self.element_type == 1: # CROD
            #result_name = 'thermalLoad_2D_3D'

            if self._results.is_not_saved(result_name):
                return ndata
            self._results._found_result(result_name)
            slot = getattr(self, result_name)
            if self.format_code == 1 and self.num_wide == 9:  # real - 2D
                # [33, 53, 64, 74, 75]
                ntotal = 36
                nelements = ndata // ntotal
                auto_return, is_vectorized = self._create_oes_object4(
                    nelements, result_name, slot, obj_vector_real)
                if auto_return:
                    return nelements * self.num_wide * 4
                obj = self.obj
                #if self.is_debug_file:
                    #self.binary_debug.write('  [cap, element1, element2, ..., cap]\n')
                    #self.binary_debug.write('  cap = %i  # assume 1 cap when there could have been multiple\n' % ndata)
                    #self.binary_debug.write('  #elementi = [eid_device, axial, torque]\n')
                    #self.binary_debug.write('  nelements=%i; nnodes=1 # centroid\n' % nelements)

                if self.use_vector and is_vectorized:
                    n = nelements * 4 * self.num_wide
                    itotal = obj.ielement
                    ielement2 = obj.itotal + nelements
                    itotal2 = ielement2

                    floats = fromstring(data, dtype=self.fdtype).reshape(nelements, 9)
                    obj._times[obj.itime] = dt
                    #if obj.itime == 0:
                    ints = fromstring(data, dtype=self.idtype).reshape(nelements, 9)
                    eids = ints[:, 0] // 10
                    assert eids.min() > 0, eids.min()
                    obj.element[itotal:itotal2] = eids
                    strings = fromstring(data, dtype=self._endian + 'S4').reshape(nelements, 9)
                    obj.element_data_type[itotal:itotal2] = array([s1+s2 for s1, s2 in zip(strings[:, 1], strings[:, 2])])

                    #[etype, xgrad, ygrad, zgrad, xflux, yflux, zflux]
                    obj.data[obj.itime, itotal:itotal2, :] = floats[:, 3:]
                    obj.itotal = itotal2
                    obj.ielement = ielement2
                else:
                    # no zed on this element for some reason...
                    s = Struct(b(self._endian + 'i8s6f'))
                    for i in range(nelements):
                        edata = data[n:n+ntotal]
                        n += ntotal
                        out = s.unpack(edata)
                        (eid_device, etype, xgrad, ygrad, zgrad, xflux, yflux, zflux) = out
                        eid = eid_device // 10
                        obj.add_sort1(dt, eid, etype, xgrad, ygrad, zgrad, xflux, yflux, zflux)

            elif self.format_code == 1 and self.num_wide == 10:  # real - 3D
                # [39, 67, 68]:  # HEXA,PENTA
                ntotal = 40
                nelements = ndata // ntotal
                auto_return, is_vectorized = self._create_oes_object4(
                    nelements, result_name, slot, obj_vector_real)
                if auto_return:
                    return nelements * self.num_wide * 4
                obj = self.obj
                assert nelements > 0, 'ndata=%s ntotal=%s' % (ndata, ntotal)
                if self.use_vector and is_vectorized:
                    n = nelements * 4 * self.num_wide
                    itotal = obj.ielement
                    ielement2 = obj.itotal + nelements
                    itotal2 = ielement2

                    floats = fromstring(data, dtype=self.fdtype).reshape(nelements, 10)
                    obj._times[obj.itime] = dt
                    if obj.itime == 0:
                        ints = fromstring(data, dtype=self.idtype).reshape(nelements, 10)
                        eids = ints[:, 0] // 10
                        assert eids.min() > 0, eids.min()
                        obj.element[itotal:itotal2] = eids
                        strings = fromstring(data, dtype=self._endian + 'S4').reshape(nelements, 10)
                        obj.element_data_type[itotal:itotal2] = array([s1+s2 for s1, s2 in zip(strings[:, 1], strings[:, 2])])

                    #[etype, xgrad, ygrad, zgrad, xflux, yflux, zflux, zed]
                    obj.data[obj.itime, itotal:itotal2, :] = floats[:, 3:-1]
                    obj.itotal = itotal2
                    obj.ielement = ielement2
                else:
                    s = Struct(b(self._endian + 'i8s6fi'))
                    for i in range(nelements):
                        edata = data[n:n+ntotal]
                        n += ntotal
                        out = s.unpack(edata)
                        (eid_device, etype, xgrad, ygrad, zgrad, xflux, yflux, zflux, zed) = out
                        eid = eid_device // 10
                        obj.add_sort1(dt, eid, etype, xgrad, ygrad, zgrad, xflux, yflux, zflux)
            else:
                raise RuntimeError(self.code_information())
        elif self.element_type in [107, 108, 109]:  # CHBDYE, CHBDYG, CHBDYP
            # 107-CHBDYE
            # 108-CHBDYG
            # 109-CHBDYP
            #if self.table_name in ['OEF1X']:
            if self.element_type == 107:
                result_name = 'chbdye_thermal_load'
            elif self.element_type == 108:
                result_name = 'chbdyg_thermal_load'
            elif self.element_type == 109:
                result_name = 'chbdyp_thermal_load'
            else:
                raise NotImplementedError('element_type=%s element_name=%s' % (
                    self.element_type, self.element_name))
            #elif self.table_name in ['HOEF1']:
                #if self.element_type == 107:
                    #result_name = 'chbdye_thermal_flux'
                #elif self.element_type == 108:
                    #result_name = 'chbdyg_thermal_flux'
                #elif self.element_type == 109:
                    #result_name = 'chbdyp_thermal_flux'
                #else:
                    #raise NotImplementedError('element_type=%s element_name=%s' % (
                        #self.element_type, self.element_name))
            #else:
                #raise NotImplementedError(msg)

            if self.format_code == 1 and self.num_wide == 8:  # real
                #result_name = 'thermalLoad_CHBDY'
                if self._results.is_not_saved(result_name):
                    return ndata
                self._results._found_result(result_name)
                slot = getattr(self, result_name)

                if self.format_code == 1 and self.num_wide == 8:  # real
                    obj_vector_real = RealChbdyHeatFluxArray
                    ntotal = 32
                    nelements = ndata // ntotal
                    auto_return, is_vectorized = self._create_oes_object4(
                        nelements, result_name, slot, obj_vector_real)
                    if auto_return:
                        return nelements * self.num_wide * 4
                    obj = self.obj
                    #if self.is_debug_file:
                        #self.binary_debug.write('  [cap, element1, element2, ..., cap]\n')
                        #self.binary_debug.write('  cap = %i  # assume 1 cap when there could have been multiple\n' % ndata)
                        #self.binary_debug.write('  #elementi = [eid_device, etype, fapplied, free_conv, force_conv, frad, ftotal]\n')
                        #self.binary_debug.write('  nelements=%i; nnodes=1 # centroid\n' % nelements)

                    if self.use_vector and is_vectorized:
                        n = nelements * 4 * self.num_wide
                        itotal = obj.ielement
                        ielement2 = obj.itotal + nelements
                        itotal2 = ielement2

                        floats = fromstring(data, dtype=self.fdtype).reshape(nelements, 8)
                        obj._times[obj.itime] = dt
                        if obj.itime == 0:
                            ints = fromstring(data, dtype=self.idtype).reshape(nelements, 8)
                            eids = ints[:, 0] // 10
                            assert eids.min() > 0, eids.min()
                            obj.element[itotal:itotal2] = eids
                            #obj.element_type[obj.itime, itotal:itotal2, :] = strings[:, 3:]

                        #[fapplied, free_conv, force_conv, frad, ftotal]
                        obj.data[obj.itime, itotal:itotal2, :] = floats[:, 3:]
                        obj.itotal = itotal2
                        obj.ielement = ielement2
                    else:
                        s1 = Struct(b(self._endian + 'i8s5f'))
                        for i in range(nelements):
                            edata = data[n:n+32]
                            n += ntotal
                            out = s1.unpack(edata)
                            (eid_device, etype, fapplied, free_conv, force_conv, frad, ftotal) = out
                            eid = eid_device // 10

                            if self.is_debug_file:
                                self.binary_debug.write('  %s -> [%s, %s, %s, %s, %s, %s, %s]\n'
                                                        % (eid, eid_device, etype, fapplied, free_conv, force_conv, frad, ftotal))
                            obj.add(dt, eid, etype, fapplied, free_conv, force_conv, frad, ftotal)
            else:
                msg = self.code_information()
                return self._not_implemented_or_skip(data, ndata, msg)

        elif self.element_type == 110:
            # 110-CONV
            result_name = 'thermalLoad_CONV'
            if self._results.is_not_saved(result_name):
                return ndata
            self._results._found_result(result_name)
            slot = getattr(self, result_name)

            if self.format_code == 1 and self.num_wide == 4:
                ntotal = 16
                nelements = ndata // ntotal

                if 0:
                    if self.read_mode == 1:
                        return ndata
                    self.create_transient_object(self.thermalLoad_CONV, HeatFlux_CONV)
                else:
                    auto_return, is_vectorized = self._create_oes_object4(
                        nelements, result_name, slot, RealConvHeatFluxArray)
                    if auto_return:
                        return nelements * self.num_wide * 4
                obj = self.obj
                #if self.is_debug_file:
                    #self.binary_debug.write('  [cap, element1, element2, ..., cap]\n')
                    #self.binary_debug.write('  cap = %i  # assume 1 cap when there could have been multiple\n' % ndata)
                    #self.binary_debug.write('  #elementi = [eid_device, etype, fapplied, free_conv, force_conv, frad, ftotal]\n')
                    #self.binary_debug.write('  nelements=%i; nnodes=1 # centroid\n' % nelements)

                if self.use_vector and is_vectorized:
                    n = nelements * 4 * self.num_wide
                    ielement = obj.ielement
                    ielement2 = ielement + nelements

                    floats = fromstring(data, dtype=self.fdtype).reshape(nelements, 4)
                    obj._times[obj.itime] = dt
                    if obj.itime == 0:
                        ints = fromstring(data, dtype=self.idtype).reshape(nelements, 4)
                        eids = ints[:, 0] // 10
                        nids = ints[:, 2]
                        assert eids.min() > 0, eids.min()
                        assert nids.min() >= 0, nids.min()
                        obj.element_node[ielement:ielement2, 0] = eids
                        obj.element_node[ielement:ielement2, 1] = nids

                    #[eid, free_conv, cntl_node, free_conv_k]
                    obj.data[obj.itime, ielement:ielement2, :] = floats[:, [1, 3]]
                    obj.itotal = ielement2
                    obj.ielement = ielement2
                else:
                    s1 = Struct(b(self._endian + 'ifif'))
                    for i in range(nelements):
                        edata = data[n:n+16]
                        n += 16
                        out = s1.unpack(edata)
                        (eid_device, free_conv, cntl_node, free_conv_k) = out
                        eid = eid_device // 10
                        assert cntl_node >= 0, cntl_node
                        obj.add(dt, eid, cntl_node, free_conv, free_conv_k)
            else:
                msg = self.code_information()
                return self._not_implemented_or_skip(data, ndata, msg)

        elif self.element_type in [145, 146, 147]:  # VUHEXA,VUPENTA,VUTETRA
            # TODO: vectorize

            # 145-VUHEXA
            # 146-VUPENTA
            # 147-VUTETRA
            if self.element_type == 147:  # VUTETRA
                nnodes = 4
            elif self.element_type == 146:  # VUPENTA
                nnodes = 6
            elif self.element_type == 145:  # VUHEXA
                nnodes = 8
            else:
                msg = self.code_information()
                return self._not_implemented_or_skip(data, ndata, msg)

            result_name = 'thermalLoad_VU_3D'
            if self._results.is_not_saved(result_name):
                return ndata
            self._results._found_result(result_name)
            slot = getattr(self, result_name)

            numwide_real = 2 + 7 * nnodes
            if self.format_code == 1 and self.num_wide == numwide_real:  # real
                ntotal = 8 + 28 * nnodes
                nelements = ndata // ntotal
                if 1:
                    if self.read_mode == 1:
                        return ndata
                    self.create_transient_object(self.thermalLoad_VU_3D, HeatFlux_VU_3D)
                    is_vectorized = False
                else:
                    nelements = ndata // ntotal
                    auto_return, is_vectorized = self._create_oes_object4(
                        nelements, result_name, slot, RealHeatFluxVU3DArray)
                    if auto_return:
                        self._data_factor = nnodes
                        return nelements * self.num_wide * 4

                obj = self.obj
                if self.use_vector and is_vectorized:
                    n = nelements * 4 * self.num_wide
                    ielement = obj.ielement
                    ielement2 = ielement + nelements
                    itotal = obj.itotal
                    itotal2 = itotal + nelements * nnodes

                    floats = fromstring(data, dtype=self.fdtype).reshape(nelements, numwide_real)
                    floats2 = floats[:, 2:].reshape(nelements * nnodes, 7)
                    obj._times[obj.itime] = dt
                    #if obj.itime == 0:
                    ints = fromstring(data, dtype=self.idtype).reshape(nelements, numwide_real)
                    ints2 = ints[:, 2:].reshape(nelements * nnodes, 7)
                    eids = ints[:, 0] // 10
                    parent = ints[:, 1]
                    assert eids.min() > 0, eids.min()
                    obj.element_parent[ielement:ielement2] = eids
                    obj.element_parent[ielement:ielement2] = parent

                    #[vugrid, xgrad, ygrad, zgrad, xflux, yflux, zflux]
                    obj.vugrid[obj.itime, itotal:itotal2, :] = ints2[:, 0]
                    obj.data[obj.itime, itotal:itotal2, :] = floats2[:, 3:]
                    obj.itotal = itotal2
                    obj.ielement = ielement2
                else:
                    s1 = self.struct_2i
                    s2 = Struct(b(self._endian + 'i6f'))
                    grad_fluxes = []
                    for i in range(nelements):
                        out = s1.unpack(data[n:n+8])
                        n += 8
                        (eid_device, parent) = out
                        eid = eid_device // 10
                        for i in range(nnodes):
                            out = s2.unpack(data[n:n+28])
                            grad_fluxes.append(out)
                            n += 28
                        obj.add_sort1(dt, eid, parent, grad_fluxes)
            else:
                msg = self.code_information()
                return self._not_implemented_or_skip(data, ndata, msg)

        elif self.element_type in [189, 190]:  # VUQUAD,VUTRIA
            # TODO: vectorize
            # 189-VUQUAD
            # 190-VUTRIA
            if self.format_code == 1 and self.num_wide == 27:  # real
                result_name = 'thermalLoad_VU'
                if self._results.is_not_saved(result_name):
                    return ndata
                self._results._found_result(result_name)
                slot = getattr(self, result_name)

                if self.element_type == 189:
                    nnodes = 4
                elif self.element_type == 190:
                    nnodes = 3
                elif self.element_type == 191:
                    nnodes = 2
                else:
                    raise NotImplementedError(self.code_information())

                numwide = 6 + 7 * nnodes
                ntotal = 24 + 28 * nnodes
                assert ntotal == numwide * 4
                nelements = ndata // ntotal
                auto_return, is_vectorized = self._create_oes_object4(
                    nelements, result_name, slot, RealHeatFluxVUArray)
                if auto_return:
                    self._data_factor = nnodes
                    return nelements * self.num_wide * 4

                obj = self.obj
                if self.use_vector and is_vectorized:
                    n = nelements * 4 * self.num_wide
                    itotal = obj.itotal
                    itotal2 = itotal + nelements * nnodes
                    ielement = obj.ielement
                    ielement2 = obj.ielement + nelements

                    ints = fromstring(data, dtype=self.idtype).reshape(nelements, numwide)
                    floats = fromstring(data, dtype=self.fdtype).reshape(nelements, numwide)
                    obj._times[obj.itime] = dt
                    if obj.itime == 0:
                        eids = ints[:, 0] // 10
                        parent = ints[:, 1]
                        coord = ints[:, 2]
                        # icord - 4s
                        theta = ints[:, 4]
                        assert eids.min() > 0, eids.min()
                        obj.element_parent_coord_icord[ielement:ielement2, 0] = eids
                        obj.element_parent_coord_icord[ielement:ielement2, 1] = parent
                        obj.element_parent_coord_icord[ielement:ielement2, 2] = coord
                        obj.element_parent_coord_icord[ielement:ielement2, 3] = theta
                    #obj.data[itotal:itotal2, 0] = floats[:, 4]

                    ints2 = ints[:, 6:].reshape(nelements * nnodes, 7)
                    floats2 = floats[:, 6:].reshape(nelements * nnodes, 7)

                    #[vugrid]
                    #[xgrad, ygrad, zgrad, xflux, yflux, zflux]
                    obj.int_data[obj.itime, itotal:itotal2, 0] = ints2[:, 0]
                    obj.data[obj.itime, itotal:itotal2, :] = floats2[:, 1:]
                    obj.itotal = itotal2
                    obj.ielement = ielement2
                else:
                    s1 = Struct(b(self._endian + '3i4s2i'))
                    s2 = Struct(b(self._endian + 'i6f'))
                    for i in range(nelements):
                        edata = data[n:n+24]  # 6*4
                        n += 24

                        out = s1.unpack(edata)
                        (eid_device, parent, coord, icord, theta, null) = out
                        eid = eid_device // 10
                        data_in = [eid, parent, coord, icord, theta]
                        #self.log.debug('RealHeatFluxVUArray = %s' % data_in)
                        grad_fluxes = []
                        for i in range(nnodes):
                            edata = data[n:n+28]  # 7*4
                            n += 28
                            out = s2.unpack(edata)
                            grad_fluxes.append(out)
                        data_in.append(grad_fluxes)
                        obj.add_sort1(dt, eid, parent, coord, icord, theta, grad_fluxes)
            else:
                msg = self.code_information()
                return self._not_implemented_or_skip(data, ndata, msg)

        elif self.element_type == 191:  # VUBEAM
            # TODO: vectorize
            nnodes = 2
            numwide_real = 4 + 7 * nnodes

            result_name = 'thermalLoad_VUBeam'
            if self._results.is_not_saved(result_name):
                return ndata
            self._results._found_result(result_name)
            slot = getattr(self, result_name)
            if self.format_code == 1 and self.num_wide == numwide_real:  # real
                ntotal = 16 + 28 * nnodes
                nelements = ndata // ntotal

                if 0:
                    if self.read_mode == 1:
                        return ndata
                    #assert self.num_wide==27,self.code_information()
                    self.create_transient_object(self.thermalLoad_VUBeam, HeatFlux_VUBEAM)
                else:
                    auto_return, is_vectorized = self._create_oes_object4(
                        nelements, result_name, slot, RealHeatFluxVUBeamArray)
                    if auto_return:
                        self._data_factor = nnodes
                        return nelements * self.num_wide * 4

                obj = self.obj
                if self.use_vector and is_vectorized:
                    n = nelements * 4 * self.num_wide
                    itotal = obj.itotal
                    itotal2 = itotal + nelements * nnodes
                    ielement = obj.ielement
                    ielement2 = obj.ielement + nelements

                    ints = fromstring(data, dtype=self.idtype).reshape(nelements, numwide_real)
                    floats = fromstring(data, dtype=self.fdtype).reshape(nelements, numwide_real)
                    obj._times[obj.itime] = dt
                    if obj.itime == 0:
                        eids = ints[:, 0] // 10
                        parent = ints[:, 1]
                        coord = ints[:, 2]
                        # icord - 4s
                        assert eids.min() > 0, eids.min()
                        obj.element_parent_coord[ielement:ielement2, 0] = eids
                        obj.element_parent_coord[ielement:ielement2, 1] = parent
                        obj.element_parent_coord[ielement:ielement2, 2] = coord
                    #obj.data[itotal:itotal2, 0] = floats[:, 4]

                    ints2 = ints[:, 4:].reshape(nelements * nnodes, 7)
                    floats2 = floats[:, 4:].reshape(nelements * nnodes, 7)

                    #[vugrid]
                    #[xgrad, ygrad, zgrad, xflux, yflux, zflux]
                    obj.vugrid[obj.itime, itotal:itotal2, 0] = ints2[:, 0]
                    obj.data[obj.itime, itotal:itotal2, :] = floats2[:, 1:]
                    obj.itotal = itotal2
                    obj.ielement = ielement2
                else:
                    s1 = Struct(b(self._endian + 'iii4s'))
                    s2 = Struct(b(self._endian + 'i6f'))
                    for i in range(nelements):
                        edata = data[n:n+16]  # 4*4
                        n += 16

                        out = s1.unpack(edata)
                        (eid_device, parent, coord, icord) = out
                        eid = eid_device // 10
                        data_in = [eid, parent, coord, icord]
                        self.log.debug('VUBeam %s' % data_in)

                        grad_fluxes = []
                        for i in range(nnodes):
                            edata = data[n:n+28]  # 7*4
                            n += 28
                            out = s2.unpack(edata)
                            grad_fluxes.append(out)
                        data_in.append(grad_fluxes)
                        obj.add_sort1(dt, eid, parent, coord, icord, grad_fluxes)

            elif self.element_type == 118:  # CWELDP
                msg = 'OEF sort1 thermal Type=%s num=%s' % (self.element_name, self.element_type)
                return self._not_implemented_or_skip(data, ndata, msg)
            else:
                msg = 'OEF sort1 thermal Type=%s num=%s' % (self.element_name, self.element_type)
                return self._not_implemented_or_skip(data, ndata, msg)
        else:
            msg = 'OEF sort1 thermal Type=%s num=%s' % (self.element_name, self.element_type)
            return self._not_implemented_or_skip(data, ndata, msg)


        assert self.thermal == 1, self.thermal
        assert ndata > 0, ndata
        try:
            assert nelements > 0, 'nelements=%r element_type=%s element_name=%r\n%s' % (nelements, self.element_type, self.element_name, self.code_information())
        except UnboundLocalError:
            raise UnboundLocalError('element_name=%r' % self.element_name)
        #assert ndata % ntotal == 0, '%s n=%s nwide=%s len=%s ntotal=%s' % (self.element_name, ndata % ntotal, ndata % self.num_wide, ndata, ntotal)
        assert self.num_wide * 4 == ntotal, 'numwide*4=%s ntotal=%s' % (self.num_wide*4, ntotal)
        assert n > 0, 'n=%s element_type=%s element_name=%s numwide=%s' % (
            n, self.element_type, self.element_name, self.num_wide)
        return n

    def _print_obj_name_on_crash(func):
        """
        Decorator debugging function to print the object name and an needed parameters
        """
        def new_func(self, data):
            """
            The actual function exec'd by the decorated function.
            """
            try:
                n = func(self, data)
            except:
                raise
                print("----------")
                try:
                    print(self.obj)
                except:
                    print("error printing %r" % self.obj.__class__.__name__)
                print(self.data_code)
                if self.obj is not None:
                    #from pyNastran.utils import object_attributes
                    #print object_attributes(self.obj)
                    print(self.obj.data_code)
                print("----------")
                raise
            return n
        return new_func

    # @_print_obj_name_on_crash
    def _read_oef1_loads(self, data, ndata):
        """
        Reads the OEF1 table; stores the element forces/heat flux.
        """
        if self._results.is_not_saved('element_forces'):
            return ndata
        flag = 'element_id'
        (num_wide_real, num_wide_imag) = self._oef_force_code()
        if self.is_debug_file:
            self.binary_debug.write('  num_wide_real = %r\n' % num_wide_real)
            self.binary_debug.write('  num_wide_imag = %r\n' % num_wide_imag)

        n = 0
        is_magnitude_phase = self.is_magnitude_phase()
        dt = self.nonlinear_factor

        if self.element_type in [1, 3, 10]:  # rods
            #1-CROD
            #3-CTUBE
            #10-CONROD
            obj_real = RealRodForceArray
            obj_complex = ComplexRodForceArray
            if self.element_type == 1: # CROD
                result_name = 'crod_force'
                slot = self.crod_force
            elif self.element_type == 3:  # CTUBE
                result_name = 'ctube_force'
                slot = self.ctube_force
            elif self.element_type == 10:  # CONROD
                result_name = 'conrod_force'
                slot = self.conrod_force
            else:
                msg = 'sort1 Type=%s num=%s' % (self.element_name, self.element_type)
                return self._not_implemented_or_skip(data, ndata, msg)

            if self._results.is_not_saved(result_name):
                return ndata
            self._results._found_result(result_name)

            slot = getattr(self, result_name)
            if self.format_code == 1 and self.num_wide == 3: # real
                ntotal = 3 * 4
                nelements = ndata // ntotal

                auto_return, is_vectorized = self._create_oes_object4(
                    nelements, result_name, slot, obj_real)
                if auto_return:
                    return nelements * self.num_wide * 4
                obj = self.obj
                if self.is_debug_file:
                    self.binary_debug.write('  [cap, element1, element2, ..., cap]\n')
                    self.binary_debug.write('  cap = %i  # assume 1 cap when there could have been multiple\n' % ndata)
                    self.binary_debug.write('  #elementi = [eid_device, axial, torque]\n')
                    self.binary_debug.write('  nelements=%i; nnodes=1 # centroid\n' % nelements)

                if self.use_vector and is_vectorized:
                    n = nelements * 4 * self.num_wide
                    itotal = obj.ielement
                    ielement2 = obj.itotal + nelements
                    itotal2 = ielement2

                    floats = fromstring(data, dtype=self.fdtype).reshape(nelements, 3)
                    obj._times[obj.itime] = dt
                    if obj.itime == 0:
                        ints = fromstring(data, dtype=self.idtype).reshape(nelements, 3)
                        eids = ints[:, 0] // 10
                        assert eids.min() > 0, eids.min()
                        obj.element[itotal:itotal2] = eids

                    #[axial, torsion]
                    obj.data[obj.itime, itotal:itotal2, :] = floats[:, 1:]
                    obj.itotal = itotal2
                    obj.ielement = ielement2
                else:
                    auto_return, is_vectorized = self._create_oes_object4(
                        nelements, result_name, slot, obj_real)
                    if auto_return:
                        return nelements * self.num_wide * 4

                    obj = self.obj
                    s = Struct(b(self._endian + 'iff'))  # 3
                    for i in range(nelements):
                        edata = data[n:n+ntotal]
                        out = s.unpack(edata)
                        (eid_device, axial, torque) = out
                        eid = eid_device // 10
                        if self.is_debug_file:
                            self.binary_debug.write('OEF_Rod - %s\n' % (str(out)))
                        obj.add(dt, eid, axial, torque)
                        n += ntotal

            elif self.format_code in [2, 3] and self.num_wide == 5: # imag
                ntotal = 5 * 4
                nelements = ndata // ntotal

                auto_return, is_vectorized = self._create_oes_object4(
                    nelements, result_name, slot, obj_complex)
                if auto_return:
                    return nelements * self.num_wide * 4

                obj = self.obj
                if self.is_debug_file:
                    self.binary_debug.write('  [cap, element1, element2, ..., cap]\n')
                    self.binary_debug.write('  cap = %i  # assume 1 cap when there could have been multiple\n' % ndata)
                    self.binary_debug.write('  #elementi = [eid_device, axial, torque]\n')
                    self.binary_debug.write('  nelements=%i; nnodes=1 # centroid\n' % nelements)

                if self.use_vector and is_vectorized:
                    n = nelements * 4 * self.num_wide
                    itotal = obj.ielement
                    ielement2 = obj.itotal + nelements
                    itotal2 = ielement2

                    floats = fromstring(data, dtype=self.fdtype).reshape(nelements, 5)
                    obj._times[obj.itime] = dt
                    if obj.itime == 0:
                        ints = fromstring(data, dtype=self.idtype).reshape(nelements, 5)
                        eids = ints[:, 0] // 10
                        assert eids.min() > 0, eids.min()
                        obj.element[itotal:itotal2] = eids

                    #[axial_force, torque]
                    #(eid_device, axial_real, torque_real, axial_imag, torque_imag) = out
                    if is_magnitude_phase:
                        mag = floats[:, [1, 2]]
                        phase = floats[:, [3, 4]]
                        rtheta = radians(phase)
                        real_imag = mag * (cos(rtheta) + 1.j * sin(rtheta))
                    else:
                        real = floats[:, [1, 2]]
                        imag = floats[:, [3, 4]]
                        real_imag = real + 1.j * imag
                    obj.data[obj.itime, itotal:itotal2, :] = real_imag
                    obj.itotal = itotal2
                    obj.ielement = ielement2
                else:
                    s = Struct(b(self._endian + 'i4f'))  # 5
                    for i in range(nelements):
                        edata = data[n:n+20]

                        out = s.unpack(edata)
                        if self.is_debug_file:
                            self.binary_debug.write('OEF_CRod - %s\n' % (str(out)))
                        (eid_device, axial_real, torque_real, axial_imag, torque_imag) = out
                        eid = eid_device // 10
                        if is_magnitude_phase:
                            axial = polar_to_real_imag(axial_real, axial_imag)
                            torque = polar_to_real_imag(torque_real, torque_imag)
                        else:
                            axial = complex(axial_real, axial_imag)
                            torque = complex(torque_real, torque_imag)

                        obj.add_sort1(dt, eid, axial, torque)
                        n += ntotal
            else:
                msg = self.code_information()
                return self._not_implemented_or_skip(data, ndata, msg)

        elif self.element_type == 2:  # cbeam
            #2-CBEAM
            result_name = 'cbeam_force'
            if self._results.is_not_saved(result_name):
                return ndata
            self._results._found_result(result_name)
            slot = getattr(self, result_name)
            if self.format_code == 1 and self.num_wide == 9:  # real centroid ???
                raise RuntimeError('is this used?')
                auto_return, is_vectorized = self._create_oes_object4(
                    nelements, result_name, slot, RealCBeamForceArray)
                if auto_return:
                    return nelements * self.num_wide * 4

                obj = self.obj
                #is_vectorized = False
                if self.use_vector and is_vectorized:
                    n = nelements * 4 * self.num_wide
                    itotal = obj.itotal
                    itotal2 = obj.itotal + nelements
                    ielement2 = obj.ielement + nelements

                    floats = fromstring(data, dtype=self.fdtype).reshape(nelements, 9)[:, 1:]
                    obj._times[obj.itime] = dt
                    if obj.itime == 0:
                        ints = fromstring(data, dtype=self.idtype).reshape(nelements, 9)
                        eids = ints[:, 0] // 10
                        assert eids.min() > 0, eids.min()
                        assert 0 not in eids, eids

                        obj.element[itotal:itotal2] = eids
                        obj.element_node[itotal:itotal2, 0] = eids
                        #obj.element_node[itotal:itotal2, 1] = nids

                    #[sd, bm1, bm2, ts1, ts2, af, ttrq, wtrq]
                    obj.data[obj.itime, itotal:itotal2, :] = floats
                    obj.itotal = itotal2
                    obj.ielement = ielement2
                else:
                    s = Struct(b(self._endian + 'i8f'))  # 36
                    ntotal = 36
                    nelements = ndata // ntotal
                    obj = self.obj
                    for i in range(nelements):
                        edata = data[n:n+36]
                        out = s.unpack(edata)
                        if self.is_debug_file:
                            self.binary_debug.write('OEF_Beam - %s\n' % (str(out)))
                        (eid_device, sd, bm1, bm2, ts1, ts2, af, ttrq, wtrq) = out
                        eid = eid_device // 10
                        n += 36

            elif self.format_code == 1 and self.num_wide == 100:  # real
                ntotal = 400  # 1+(10-1)*11=100 ->100*4 = 400
                nelements = ndata // ntotal
                auto_return, is_vectorized = self._create_oes_object4(
                    nelements, result_name, slot, RealCBeamForceArray)
                if auto_return:
                    self._data_factor = 11
                    return nelements * self.num_wide * 4
                obj = self.obj

                if self.use_vector and is_vectorized:
                    n = nelements * 4 * self.num_wide
                    itotal = obj.itotal
                    itotal2 = obj.itotal + nelements * 11
                    ielement = obj.ielement
                    ielement2 = obj.ielement + nelements

                    floats = fromstring(data, dtype=self.fdtype).reshape(nelements, 100)[:, 1:]
                    obj._times[obj.itime] = dt
                    if obj.itime == 0:
                        ints = fromstring(data, dtype=self.idtype).reshape(nelements, 100)
                        eids = ints[:, 0] // 10
                        assert eids.min() > 0, eids.min()
                        assert 0 not in eids, eids
                        eids2 = np.repeat(eids, 11)

                        ints2 = ints[:, 1:].reshape(nelements * 11, 9)
                        nids = ints2[:, 0]

                        obj.element[itotal:itotal2] = eids2
                        obj.element_node[itotal:itotal2, 0] = eids2
                        obj.element_node[itotal:itotal2, 1] = nids

                    #[nid, sd, bm1, bm2, ts1, ts2, af, ttrq, wtrq]
                    floats2 = floats.reshape(nelements * 11, 9)[:, 1:]
                    #sd = floats2[:, 0]
                    #obj.data[obj.itime, itotal:itotal2, :] = sd
                    obj.data[obj.itime, itotal:itotal2, :] = floats2

                    obj.itotal = itotal2
                    obj.ielement = ielement2
                else:
                    s1 = self.struct_i
                    s2 = Struct(b(self._endian + 'i8f'))  # 36
                    for i in range(nelements):
                        edata = data[n:n+4]
                        eid_device, = s1.unpack(edata)
                        eid = eid_device // 10
                        n += 4

                        for i in range(11):
                            edata = data[n:n+36]
                            out = s2.unpack(edata)
                            if self.is_debug_file:
                                self.binary_debug.write('OEF_Beam - %s\n' % (str(out)))
                            (nid, sd, bm1, bm2, ts1, ts2, af, ttrq, wtrq) = out

                            if i == 0:  # isNewElement
                                obj.add_new_element_sort1(
                                    dt, eid, nid, sd, bm1, bm2, ts1, ts2, af, ttrq, wtrq)
                            elif sd > 0.:
                                obj.add_sort1(
                                    dt, eid, nid, sd, bm1, bm2, ts1, ts2, af, ttrq, wtrq)
                            n += 36
            elif self.format_code in [2, 3] and self.num_wide == 177: # imag
                ntotal = 708  # 3*4
                nelements = ndata // ntotal

                auto_return, is_vectorized = self._create_oes_object4(
                    nelements, result_name, slot, ComplexCBeamForceArray)
                if auto_return:
                    self._data_factor = 11
                    return nelements * self.num_wide * 4

                obj = self.obj
                if self.is_debug_file:
                    self.binary_debug.write('  [cap, element1, element2, ..., cap]\n')
                    self.binary_debug.write('  cap = %i  # assume 1 cap when there could have been multiple\n' % ndata)
                    #self.binary_debug.write('  #elementi = [eid_device, force]\n')
                    #self.binary_debug.write('  nelements=%i; nnodes=1 # centroid\n' % nelements)

                if self.use_vector and is_vectorized:
                    n = nelements * 4 * self.num_wide
                    itotal = obj.itotal
                    itotal2 = obj.itotal + nelements * 11
                    ielement = obj.ielement
                    ielement2 = obj.ielement + nelements

                    floats = fromstring(data, dtype=self.fdtype).reshape(nelements, 177)[:, 1:]
                    obj._times[obj.itime] = dt
                    if obj.itime == 0:
                        ints = fromstring(data, dtype=self.idtype).reshape(nelements, 177)
                        eids = ints[:, 0] // 10
                        assert eids.min() > 0, eids.min()
                        assert 0 not in eids, eids
                        eids2 = np.repeat(eids, 11)

                        ints2 = ints[:, 1:].reshape(nelements * 11, 16)
                        nids = ints2[:, 0]

                        obj.element[itotal:itotal2] = eids2
                        obj.element_node[itotal:itotal2, 0] = eids2
                        obj.element_node[itotal:itotal2, 1] = nids

                    #[nid, sd, bm1r, bm2r, ts1r, ts2r, afr, ttrqr, wtrqr,
                    #          bm1i, bm2i, ts1i, ts2i, afi, ttrqi, wtrqi]
                    floats2 = floats.reshape(nelements * 11, 16)[:, 1:]
                    if is_magnitude_phase:
                        mag = floats2[:, 1:8] # 7
                        phase = floats2[:, 8:]
                        rtheta = radians(phase)
                        real_imag = mag * (cos(rtheta) + 1.j * sin(rtheta))
                    else:
                        real = floats2[:, 1:8]
                        imag = floats2[:, 8:]
                        real_imag = real + 1.j * imag

                    sd = floats2[:, 0]
                    obj.data[obj.itime, itotal:itotal2, 0] = sd
                    obj.data[obj.itime, itotal:itotal2, 1:] = real_imag

                    obj.itotal = itotal2
                    obj.ielement = ielement2
                else:
                    s1 = self.struct_i
                    s2 = Struct(b(self._endian + 'i15f'))
                    ntotal = 708  # (16*11+1)*4 = 177*4
                    nelements = ndata // ntotal
                    for i in range(nelements):
                        edata = data[n:n+4]
                        eid_device, = s1.unpack(edata)
                        eid = eid_device // 10

                        n += 4
                        for i in range(11):
                            edata = data[n:n+64]
                            n += 64
                            out = s2.unpack(edata)
                            if self.is_debug_file:
                                self.binary_debug.write('OEF_Beam - %s\n' % (str(out)))
                            (nid, sd,
                             bm1r, bm2r, ts1r, ts2r, afr, ttrqr, wtrqr,
                             bm1i, bm2i, ts1i, ts2i, afi, ttrqi, wtrqi) = out

                            if is_magnitude_phase:
                                bm1 = polar_to_real_imag(bm1r, bm1i)
                                bm2 = polar_to_real_imag(bm2r, bm2i)
                                ts1 = polar_to_real_imag(ts1r, ts1i)
                                ts2 = polar_to_real_imag(ts2r, ts2i)
                                af = polar_to_real_imag(afr, afi)
                                ttrq = polar_to_real_imag(ttrqr, ttrqi)
                                wtrq = polar_to_real_imag(wtrqr, wtrqi)
                            else:
                                bm1 = complex(bm1r, bm1i)
                                bm2 = complex(bm2r, bm2i)
                                ts1 = complex(ts1r, ts1i)
                                ts2 = complex(ts2r, ts2i)
                                af = complex(afr, afi)
                                ttrq = complex(ttrqr, ttrqi)
                                wtrq = complex(wtrqr, wtrqi)

                            #if i == 0:
                                #obj.add_new_element_sort1(
                                    #dt, eid, nid, sd, bm1, bm2, ts1, ts2,
                                    #af, ttrq, wtrq)
                            #elif sd > 0.:
                            obj.add_sort1(
                                dt, eid, nid, sd, bm1, bm2, ts1, ts2,
                                af, ttrq, wtrq)
                            #else:
                                ## don't add this field
                                #pass
                                #raise RuntimeError('CBEAM error; i=%s sd=%s' % (i, sd))
            else:
                msg = self.code_information()
                return self._not_implemented_or_skip(data, ndata, msg)

        elif self.element_type in [11, 12, 13, 14,   # springs
                                   20, 21, 22, 23]:  # dampers
            # 11-CELAS1
            # 12-CELAS2
            # 13-CELAS3
            # 14-CELAS4

            # 20-CDAMP1
            # 21-CDAMP2
            # 22-CDAMP3
            # 23-CDAMP4
            if self.element_type == 11:
                result_name = 'celas1_force'
                slot = self.celas1_force
                obj_real = RealSpringForceArray
                obj_complex = ComplexSpringForceArray
            elif self.element_type == 12:
                result_name = 'celas2_force'
                slot = self.celas2_force
                obj_real = RealSpringForceArray
                obj_complex = ComplexSpringForceArray
            elif self.element_type == 13:
                result_name = 'celas3_force'
                slot = self.celas3_force
                obj_real = RealSpringForceArray
                obj_complex = ComplexSpringForceArray
            elif self.element_type == 14:
                result_name = 'celas4_force'
                slot = self.celas4_force
                obj_real = RealSpringForceArray
                obj_complex = ComplexSpringForceArray

            elif self.element_type == 20:
                result_name = 'cdamp1_force'
                slot = self.cdamp1_force
                obj_real = RealDamperForceArray
                obj_complex = ComplexDamperForceArray
            elif self.element_type == 21:
                result_name = 'cdamp2_force'
                slot = self.cdamp2_force
                obj_real = RealDamperForceArray
                obj_complex = ComplexDamperForceArray
            elif self.element_type == 22:
                result_name = 'cdamp3_force'
                slot = self.cdamp3_force
                obj_real = RealDamperForceArray
                obj_complex = ComplexDamperForceArray
            elif self.element_type == 23:
                result_name = 'cdamp4_force'
                slot = self.cdamp4_force
                obj_real = RealDamperForceArray
                obj_complex = ComplexDamperForceArray
            else:
                raise NotImplementedError(self.code_information())

            if self._results.is_not_saved(result_name):
                return ndata
            self._results._found_result(result_name)

            if self.format_code == 1 and self.num_wide == 2:  # real
                ntotal = 8 # 2 * 4
                nelements = ndata // ntotal
                auto_return, is_vectorized = self._create_oes_object4(
                    nelements, result_name, slot, obj_real)
                if auto_return:
                    return nelements * self.num_wide * 4

                obj = self.obj
                if self.use_vector and is_vectorized:
                    n = nelements * 4 * self.num_wide
                    itotal = obj.ielement
                    ielement2 = obj.itotal + nelements
                    itotal2 = ielement2

                    floats = fromstring(data, dtype=self.fdtype).reshape(nelements, 2)
                    obj._times[obj.itime] = dt
                    if obj.itime == 0:
                        ints = fromstring(data, dtype=self.idtype).reshape(nelements, 2)
                        eids = ints[:, 0] // 10
                        assert eids.min() > 0, eids.min()
                        obj.element[itotal:itotal2] = eids

                    #(eid_device, force)
                    obj.data[obj.itime, itotal:itotal2, 0] = floats[:, 1]
                    obj.itotal = itotal2
                    obj.ielement = ielement2
                else:
                    s = Struct(b(self._endian + 'if'))  # 2
                    for i in range(nelements):
                        edata = data[n:n + 8]
                        out = s.unpack(edata)
                        if self.is_debug_file:
                            self.binary_debug.write('OEF_SpringDamper - %s\n' % str(out))
                        (eid_device, force) = out
                        eid = eid_device // 10
                        obj.add(dt, eid, force)
                        n += ntotal
            elif self.format_code in [2, 3] and self.num_wide == 3:  # imag
                ntotal = 12  # 3*4
                nelements = ndata // ntotal

                auto_return, is_vectorized = self._create_oes_object4(
                    nelements, result_name, slot, obj_complex)
                if auto_return:
                    return nelements * self.num_wide * 4

                obj = self.obj
                if self.is_debug_file:
                    self.binary_debug.write('  [cap, element1, element2, ..., cap]\n')
                    self.binary_debug.write('  cap = %i  # assume 1 cap when there could have been multiple\n' % ndata)
                    self.binary_debug.write('  #elementi = [eid_device, force]\n')
                    self.binary_debug.write('  nelements=%i; nnodes=1 # centroid\n' % nelements)

                if self.use_vector and is_vectorized:
                    n = nelements * 4 * self.num_wide
                    itotal = obj.ielement
                    ielement2 = obj.itotal + nelements
                    itotal2 = ielement2

                    floats = fromstring(data, dtype=self.fdtype).reshape(nelements, 3)
                    obj._times[obj.itime] = dt
                    if obj.itime == 0:
                        ints = fromstring(data, dtype=self.idtype).reshape(nelements, 3)
                        eids = ints[:, 0] // 10
                        assert eids.min() > 0, eids.min()
                        obj.element[itotal:itotal2] = eids

                    #[spring_force]
                    if is_magnitude_phase:
                        mag = floats[:, 1]
                        phase = floats[:, 2]
                        rtheta = radians(phase)
                        real_imag = mag * (cos(rtheta) + 1.j * sin(rtheta))
                    else:
                        real = floats[:, 1]
                        imag = floats[:, 2]
                        real_imag = real + 1.j * imag
                    obj.data[obj.itime, itotal:itotal2, 0] = real_imag
                    obj.itotal = itotal2
                    obj.ielement = ielement2
                else:
                    s = Struct(b(self._endian + 'i2f'))
                    for i in range(nelements):
                        edata = data[n:n + 12]
                        out = s.unpack(edata)
                        if self.is_debug_file:
                            self.binary_debug.write('OEF_SpringDamper - %s\n' % str(out))
                        (eid_device, force_real, force_imag) = out
                        eid = eid_device // 10
                        if is_magnitude_phase:
                            force = polar_to_real_imag(force_real, force_imag)
                        else:
                            force = complex(force_real, force_imag)
                        obj.add_sort1(dt, eid, force)
                        n += ntotal
            else:
                #msg = 'OEF: element_name=%s element_type=%s' % (self.element_name, self.element_type)
                msg = self.code_information()
                return self._not_implemented_or_skip(data, ndata, msg)

        elif self.element_type == 24:  # CVISC
            result_name = 'cvisc_force'
            if self._results.is_not_saved(result_name):
                return ndata
            self._results._found_result(result_name)

            slot = getattr(self, result_name)
            obj_real = RealViscForceArray

            if self.format_code == 1 and self.num_wide == 3: # real
                ntotal = 12 # 3 * 4
                nelements = ndata // ntotal

                auto_return, is_vectorized = self._create_oes_object4(
                    nelements, result_name, slot, obj_real)
                if auto_return:
                    return nelements * self.num_wide * 4

                obj = self.obj
                if self.use_vector and is_vectorized:
                    n = nelements * 4 * self.num_wide
                    itotal = obj.ielement
                    ielement2 = obj.itotal + nelements
                    itotal2 = ielement2

                    floats = fromstring(data, dtype=self.fdtype).reshape(nelements, 3)
                    obj._times[obj.itime] = dt
                    if obj.itime == 0:
                        ints = fromstring(data, dtype=self.idtype).reshape(nelements, 3)
                        eids = ints[:, 0] // 10
                        assert eids.min() > 0, eids.min()
                        obj.element[itotal:itotal2] = eids

                    #(eid_device, axial, torque)
                    obj.data[obj.itime, itotal:itotal2, :] = floats[:, 1:]
                    obj.itotal = itotal2
                    obj.ielement = ielement2
                else:
                    s = Struct(b(self._endian + 'iff'))
                    for i in range(nelements):
                        edata = data[n:n+12]

                        out = s.unpack(edata)
                        if self.is_debug_file:
                            self.binary_debug.write('OEF_CVisc - %s\n' % (str(out)))
                        (eid_device, axial, torque) = out
                        eid = eid_device // 10
                        obj.add(dt, eid, axial, torque)
                        n += ntotal
            elif self.format_code in [2, 3] and self.num_wide == 5: # complex
                ntotal = 20  # 5*4
                nelements = ndata // ntotal

                auto_return, is_vectorized = self._create_oes_object4(
                    nelements, result_name, slot, ComplexViscForceArray)
                if auto_return:
                    return nelements * self.num_wide * 4

                obj = self.obj
                if self.is_debug_file:
                    self.binary_debug.write('  [cap, element1, element2, ..., cap]\n')
                    self.binary_debug.write('  cap = %i  # assume 1 cap when there could have been multiple\n' % ndata)
                    self.binary_debug.write('  #elementi = [eid_device, axial, torque]\n')
                    self.binary_debug.write('  nelements=%i; nnodes=1 # centroid\n' % nelements)

                if self.use_vector and is_vectorized:
                    n = nelements * 4 * self.num_wide
                    itotal = obj.ielement
                    ielement2 = obj.itotal + nelements
                    itotal2 = ielement2

                    floats = fromstring(data, dtype=self.fdtype).reshape(nelements, 5)
                    obj._times[obj.itime] = dt
                    if obj.itime == 0:
                        ints = fromstring(data, dtype=self.idtype).reshape(nelements, 5)
                        eids = ints[:, 0] // 10
                        assert eids.min() > 0, eids.min()
                        obj.element[itotal:itotal2] = eids

                    #[axial_force, torque]
                    #(eid_device, axial_real, torque_real, axial_imag, torque_imag) = out
                    if is_magnitude_phase:
                        mag = floats[:, [1, 2]]
                        phase = floats[:, [3, 4]]
                        rtheta = radians(phase)
                        real_imag = mag * (cos(rtheta) + 1.j * sin(rtheta))
                    else:
                        real = floats[:, [1, 2]]
                        imag = floats[:, [3, 4]]
                        real_imag = real + 1.j * imag
                    obj.data[obj.itime, itotal:itotal2, :] = real_imag
                    obj.itotal = itotal2
                    obj.ielement = ielement2
                else:
                    s = Struct(b(self._endian + 'i4f'))  # 5
                    for i in range(nelements):
                        edata = data[n:n+20]

                        out = s.unpack(edata)
                        if self.is_debug_file:
                            self.binary_debug.write('OEF_CVisc - %s\n' % (str(out)))
                        (eid_device, axial_real, torque_real, axial_imag, torque_imag) = out
                        eid = eid_device // 10
                        if is_magnitude_phase:
                            axial = polar_to_real_imag(axial_real, axial_imag)
                            torque = polar_to_real_imag(torque_real, torque_imag)
                        else:
                            axial = complex(axial_real, axial_imag)
                            torque = complex(torque_real, torque_imag)

                        obj.add_sort1(dt, eid, axial, torque)
                        n += ntotal
            else:
                msg = self.code_information()
                return self._not_implemented_or_skip(data, ndata, msg)

        elif self.element_type == 34:  # cbar
            # 34-CBAR
            slot = self.cbar_force
            result_name = 'cbar_force'

            obj_real = RealCBarForceArray
            obj_complex = ComplexCBarForceArray

            if self.format_code == 1 and self.num_wide == 9: # real
                ntotal = 36  # 9*4
                nelements = ndata // ntotal
                auto_return, is_vectorized = self._create_oes_object4(
                    nelements, result_name, slot, obj_real)
                if auto_return:
                    return nelements * self.num_wide * 4

                obj = self.obj
                if self.use_vector and is_vectorized:
                    n = nelements * 4 * self.num_wide
                    itotal = obj.ielement
                    ielement2 = obj.itotal + nelements
                    itotal2 = ielement2

                    floats = fromstring(data, dtype=self.fdtype).reshape(nelements, 9)
                    obj._times[obj.itime] = dt
                    if obj.itime == 0:
                        ints = fromstring(data, dtype=self.idtype).reshape(nelements, 9)
                        eids = ints[:, 0] // 10
                        assert eids.min() > 0, eids.min()
                        obj.element[itotal:itotal2] = eids

                    #[bm1a, bm2a, bm1b, bm2b, ts1, ts2, af, trq]
                    obj.data[obj.itime, itotal:itotal2, :] = floats[:, 1:]
                    obj.itotal = itotal2
                    obj.ielement = ielement2
                else:
                    s = Struct(b(self._endian + 'i8f'))  # 9
                    for i in range(nelements):
                        edata = data[n:n + 36]

                        out = s.unpack(edata)
                        if self.is_debug_file:
                            self.binary_debug.write('OEF_CBar - %s\n' % (str(out)))
                        (eid_device, bm1a, bm2a, bm1b, bm2b, ts1, ts2, af, trq) = out
                        eid = eid_device // 10
                        data_in = [eid, bm1a, bm2a, bm1b, bm2b, ts1, ts2, af, trq]
                        obj.add(dt, data_in)
                        n += ntotal
            elif self.format_code in [2, 3] and self.num_wide == 17: # imag
                # TODO: vectorize
                ntotal = 68  # 17*4
                nelements = ndata // ntotal

                auto_return, is_vectorized = self._create_oes_object4(
                    nelements, result_name, slot, obj_complex)
                if auto_return:
                    return nelements * self.num_wide * 4

                obj = self.obj
                s = Struct(b(self._endian + 'i16f'))
                for i in range(nelements):
                    edata = data[n:n + 68]

                    out = s.unpack(edata)
                    (eid_device,
                     bm1ar, bm2ar, bm1br, bm2br, ts1r, ts2r, afr, trqr,
                     bm1ai, bm2ai, bm1bi, bm2bi, ts1i, ts2i, afi, trqi) = out
                    if self.is_debug_file:
                        self.binary_debug.write('OEF_CBar - %s\n' % (str(out)))
                    eid = eid_device // 10
                    if is_magnitude_phase:
                        bm1a = polar_to_real_imag(bm1ar, bm1ai)
                        bm2a = polar_to_real_imag(bm2ar, bm2ai)
                        bm1b = polar_to_real_imag(bm1br, bm1bi)
                        bm2b = polar_to_real_imag(bm2br, bm2bi)
                        ts1 = polar_to_real_imag(ts1r, ts1i)
                        ts2 = polar_to_real_imag(ts2r, ts2i)
                        af = polar_to_real_imag(afr, afi)
                        trq = polar_to_real_imag(trqr, trqi)
                    else:
                        bm1a = complex(bm1ar, bm1ai)
                        bm2a = complex(bm2ar, bm2ai)
                        bm1b = complex(bm1br, bm1bi)
                        bm2b = complex(bm2br, bm2bi)
                        ts1 = complex(ts1r, ts1i)
                        ts2 = complex(ts2r, ts2i)
                        af = complex(afr, afi)
                        trq = complex(trqr, trqi)

                    #data_in = [bm1a, bm2a, bm1b, bm2b, ts1, ts2, af, trq]
                    #print "%s" % (self.get_element_type(self.element_type)), data_in
                    #eid = obj.add_new_eid(out)
                    obj.add_sort1(dt, eid, bm1a, bm2a, bm1b, bm2b, ts1, ts2, af, trq)
                    n += ntotal
            else:
                msg = self.code_information()
                return self._not_implemented_or_skip(data, ndata, msg)
            #print self.barForces

        elif self.element_type == 100:  # cbar
            #100-BARS
            result_name = 'cbar_force_10nodes'
            self._results._found_result(result_name)
            slot = getattr(self, result_name)
            if self.format_code == 1 and self.num_wide == 8:  # real
                ntotal = 32  # 8*4
                nelements = ndata // ntotal
                auto_return, is_vectorized = self._create_oes_object4(
                    nelements, result_name, slot, RealCBar100ForceArray)
                if auto_return:
                    return nelements * self.num_wide * 4

                obj = self.obj
                if self.use_vector and is_vectorized:
                    n = nelements * 4 * self.num_wide
                    itotal = obj.ielement
                    ielement2 = obj.itotal + nelements
                    itotal2 = ielement2

                    floats = fromstring(data, dtype=self.fdtype).reshape(nelements, 8)
                    obj._times[obj.itime] = dt
                    if obj.itime == 0:
                        ints = fromstring(data, dtype=self.idtype).reshape(nelements, 8)
                        eids = ints[:, 0] // 10
                        assert eids.min() > 0, eids.min()
                        obj.element[itotal:itotal2] = eids

                    #[axial, torsion, SMa, SMt]
                    obj.data[obj.itime, itotal:itotal2, :] = floats[:, 1:]
                    obj.itotal = itotal2
                    obj.ielement = ielement2
                else:
                    s = Struct(b(self._endian + 'i7f'))
                    for i in range(nelements):
                        edata = data[n:n+32]

                        out = s.unpack(edata)
                        if self.is_debug_file:
                            self.binary_debug.write('OEF_CBar100 - %s\n' % (str(out)))
                        (eid_device, sd, bm1, bm2, ts1, ts2, af, trq) = out
                        eid = eid_device // 10
                        obj.add_sort1(dt, eid, sd, bm1, bm2, ts1, ts2, af, trq)
                        n += 32
            #elif self.format_code in [2, 3] and self.num_wide == 14:  # imag
            else:
                msg = self.code_information()
                return self._not_implemented_or_skip(data, ndata, msg)
        elif self.element_type in [33, 74]: # centroidal shells
            # 33-CQUAD4
            # 74-CTRIA3
            if self.element_type == 33:
                result_name = 'cquad4_force'
                slot = self.cquad4_force
            elif self.element_type == 74:
                result_name = 'ctria3_force'
                slot = self.ctria3_force
            else:
                msg = 'sort1 Type=%s num=%s' % (self.element_name, self.element_type)
                return self._not_implemented_or_skip(data, ndata, msg)

            assert self._data_factor == 1, self._data_factor
            if self.format_code == 1 and self.num_wide == 9:  # real
                ntotal = 36 # 9*4
                nelements = ndata // ntotal
                auto_return, is_vectorized = self._create_oes_object4(
                    nelements, result_name, slot, RealPlateForceArray)
                if auto_return:
                    return nelements * self.num_wide * 4

                obj = self.obj
                if is_vectorized:
                    n = nelements * 4 * self.num_wide
                    ielement = obj.ielement
                    ielement2 = ielement + nelements

                    floats = fromstring(data, dtype=self.fdtype).reshape(nelements, 9)
                    obj._times[obj.itime] = dt
                    if obj.itime == 0:
                        ints = fromstring(data, dtype=self.idtype).reshape(nelements, 9)
                        eids = ints[:, 0] // 10
                        assert eids.min() > 0, eids.min()
                        obj.element[ielement:ielement2] = eids

                    #[mx, my, mxy, bmx, bmy, bmxy, tx, ty]
                    obj.data[obj.itime, ielement:ielement2, :] = floats[:, 1:]
                    obj.itotal = ielement2
                    obj.ielement = ielement2
                else:
                    s = Struct(b(self._endian + 'i8f'))
                    for i in range(nelements):
                        edata = data[n:n+36]
                        out = s.unpack(edata)
                        if self.is_debug_file:
                            self.binary_debug.write('real_OEF_Plate-%s - %s\n' % (self.element_type, str(out)))
                        (eid_device, mx, my, mxy, bmx, bmy, bmxy, tx, ty) = out
                        eid = eid_device // 10
                        obj.add_sort1(dt, eid, mx, my, mxy, bmx, bmy, bmxy, tx, ty)
                        n += ntotal
            elif self.format_code in [2, 3] and self.num_wide == 17:  # imag
                ntotal = 68
                nelements = ndata // ntotal

                auto_return, is_vectorized = self._create_oes_object4(
                    nelements, result_name, slot, ComplexPlateForceArray)
                if auto_return:
                    return nelements * self.num_wide * 4

                obj = self.obj
                if self.use_vector and is_vectorized:
                    n = nelements * 4 * self.num_wide
                    ielement = obj.ielement
                    ielement2 = ielement + nelements
                    itotal = obj.itotal
                    itotal2 = itotal + nelements

                    floats = fromstring(data, dtype=self.fdtype).reshape(nelements, 17)
                    obj._times[obj.itime] = dt
                    if obj.itime == 0:
                        ints = fromstring(data, dtype=self.idtype).reshape(nelements, 17)
                        eids = ints[:, 0] // 10
                        assert eids.min() > 0, eids.min()
                        obj.element[itotal:itotal2] = eids

                    #[mx, my, mxy, bmx, bmy, bmxy, tx, ty]
                    if is_magnitude_phase:
                        mag = floats[:, [1, 2, 3, 4, 5, 6, 7, 8]]
                        phase = floats[:, [9, 10, 11, 12, 13, 14, 15, 16]]
                        rtheta = radians(phase)
                        real_imag = mag * (cos(rtheta) + 1.j * sin(rtheta))
                    else:
                        real = floats[:, [1, 2, 3, 4, 5, 6, 7, 8]]
                        imag = floats[:, [9, 10, 11, 12, 13, 14, 15, 16]]
                        real_imag = real + 1.j * imag
                    obj.data[obj.itime, itotal:itotal2, :] = real_imag
                    obj.itotal = itotal2
                    obj.ielement = ielement2
                else:
                    s = Struct(b(self._endian + 'i16f'))
                    for i in range(nelements):
                        edata = data[n:n+68]
                        out = s.unpack(edata)
                        (eid_device,
                         mxr, myr, mxyr, bmxr, bmyr, bmxyr, txr, tyr,
                         mxi, myi, mxyi, bmxi, bmyi, bmxyi, txi, tyi) = out
                        eid = eid_device // 10
                        if self.is_debug_file:
                            self.binary_debug.write('complex_OEF_Plate-%s - %s\n' % (self.element_type, str(out)))

                        if is_magnitude_phase:
                            mx = polar_to_real_imag(mxr, mxi)
                            my = polar_to_real_imag(myr, myi)
                            mxy = polar_to_real_imag(mxyr, mxyi)
                            bmx = polar_to_real_imag(bmxr, bmxi)
                            bmy = polar_to_real_imag(bmyr, bmyi)
                            bmxy = polar_to_real_imag(bmxyr, bmxyi)
                            tx = polar_to_real_imag(txr, txi)
                            ty = polar_to_real_imag(tyr, tyi)
                        else:
                            mx = complex(mxr, mxi)
                            my = complex(myr, myi)
                            mxy = complex(mxyr, mxyi)
                            bmx = complex(bmxr, bmxi)
                            bmy = complex(bmyr, bmyi)
                            bmxy = complex(bmxyr, bmxyi)
                            tx = complex(txr, txi)
                            ty = complex(tyr, tyi)
                        obj.add_sort1(dt, eid, mx, my, mxy, bmx, bmy, bmxy, tx, ty)
                        n += ntotal
            else:
                msg = self.code_information()
                return self._not_implemented_or_skip(data, ndata, msg)

        elif self.element_type in [64, 70, 75, 82, 144]: # bilinear shells
            # 64-CQUAD8
            # 70-CTRIAR
            # 75-CTRIA6
            # 82-CQUADR
            # 144-CQUAD4-bilinear
            if self.element_type == 64:
                result_name = 'cquad8_force'
                slot = self.cquad8_force
            elif self.element_type == 70:
                result_name = 'ctriar_force'
                slot = self.ctriar_force
            elif self.element_type == 75:
                result_name = 'ctria6_force'
                slot = self.ctria6_force
            elif self.element_type == 82:
                result_name = 'cquadr_force'
                slot = self.cquadr_force
            elif self.element_type == 144:
                result_name = 'cquad4_force'
                slot = self.cquad4_force
            else:
                msg = self.code_information()
                return self._not_implemented_or_skip(data, ndata, msg)

            if self.element_type in [70, 75]:  # CTRIAR,CTRIA6
                nnodes = 3
            elif self.element_type in [64, 82, 144]:  # CQUAD8,CQUADR,CQUAD4-bilinear
                nnodes = 4
            else:
                msg = 'name=%r type=%r' % (self.element_name, self.element_type)
                return self._not_implemented_or_skip(data, ndata, msg)

            slot = getattr(self, result_name)
            nnodes_all = nnodes + 1
            numwide_real = 2 + nnodes_all * 9 # centroidal node is the + 1
            numwide_imag = 2 + nnodes_all * 17

            if self._results.is_not_saved(result_name):
                return ndata
            self._results._found_result(result_name)

            if self.format_code == 1 and self.num_wide == numwide_real:  # real
                obj_real = RealPlateBilinearForceArray

                ntotal = 8 + nnodes_all * 36 # centroidal node is the + 1
                assert ntotal == self.num_wide * 4, 'ntotal=%s numwide=%s' % (ntotal, self.num_wide * 4)

                nelements = ndata // ntotal
                auto_return, is_vectorized = self._create_oes_object4(
                    nelements, result_name, slot, obj_real)
                if auto_return:
                    self._data_factor = nnodes_all
                    return nelements * self.num_wide * 4

                obj = self.obj
                if self.use_vector and is_vectorized:
                    nlayers = nelements * nnodes_all
                    n = nelements * self.num_wide * 4

                    istart = obj.itotal
                    iend = istart + nlayers
                    obj._times[obj.itime] = dt

                    if obj.itime == 0:
                        ints = fromstring(data, dtype=self.idtype).reshape(nelements, numwide_real)
                        # Nastran makes this a 4 for CQUAD4s instead
                        # of 0 like the bilinear stress element...
                        ints[:, 2] = 0

                        nids = ints[:, 2:].reshape(nlayers, 9)[:, 0]
                        eids = ints[:, 0] // 10
                        eids2 = vstack([eids] * nnodes_all).T.ravel()
                        obj.element_node[istart:iend, 0] = eids2
                        obj.element_node[istart:iend, 1] = nids

                    floats = fromstring(data, dtype=self.fdtype).reshape(nelements, numwide_real)
                    results = floats[:, 2:].reshape(nlayers, 9)[:, 1:]
                    #[mx, my, mxy, bmx, bmy, bmxy, tx, ty]
                    obj.data[obj.itime, istart:iend, :] = results
                else:
                    s1 = Struct(b(self._endian + 'i4si8f'))  # 8+36
                    s2 = Struct(b(self._endian + 'i8f')) # 36

                    for i in range(nelements):
                        edata = data[n:n + 44]

                        out = s1.unpack(edata)
                        if self.is_debug_file:
                            self.binary_debug.write('OEF_Plate2-%s - %s\n' % (self.element_type, str(out)))
                        (eid_device, term, _nid, mx, my, mxy, bmx, bmy, bmxy, tx, ty) = out
                        #term= 'CEN\'
                        nid = 0
                        eid = eid_device // 10
                        obj.add(dt, eid, term, nid, mx, my, mxy, bmx, bmy, bmxy, tx, ty)
                        n += 44
                        for i in range(nnodes):
                            edata = data[n : n + 36]
                            out = s2.unpack(edata)
                            if self.is_debug_file:
                                self.binary_debug.write('    %s\n' % (str(out)))
                            (nid, mx, my, mxy, bmx, bmy, bmxy, tx, ty) = out
                            assert nid > 0, 'nid=%s' % nid
                            obj.add(dt, eid, term, nid, mx, my, mxy, bmx, bmy, bmxy, tx, ty)
                            n += 36
            elif self.format_code in [2, 3] and self.num_wide == numwide_imag: # complex
                ntotal = numwide_imag * 4
                nelements = ndata // ntotal

                auto_return, is_vectorized = self._create_oes_object4(
                    nelements, result_name, slot, ComplexPlate2ForceArray)
                if auto_return:
                    self._data_factor = nnodes_all
                    return nelements * self.num_wide * 4

                obj = self.obj
                if self.use_vector and is_vectorized:
                    n = nelements * 4 * self.num_wide
                    itotal = obj.itotal
                    ielement = obj.ielement
                    ielement2 = obj.ielement + nelements
                    itotal2 = obj.itotal + nelements * nnodes_all

                    floats = fromstring(data, dtype=self.fdtype).reshape(nelements, numwide_imag)
                    obj._times[obj.itime] = dt
                    if obj.itime == 0:
                        ints = fromstring(data, dtype=self.idtype).reshape(nelements, numwide_imag)
                        ints[:, 2] = 0
                        ints2 = ints[:, 2:].reshape(nelements * nnodes_all, 17)

                        eids = ints[:, 0] // 10
                        nids = ints2[:, 0]
                        assert eids.min() > 0, eids.min()
                        eids2 = vstack([eids] * nnodes_all).T.ravel()
                        obj.element[ielement:ielement2] = eids
                        obj.element_node[itotal:itotal2, 0] = eids2
                        obj.element_node[itotal:itotal2, 1] = nids

                    #[mx, my, mxy, bmx, bmy, bmxy, tx, ty]
                    floats2 = floats[:, 2:].reshape(nelements * nnodes_all, 17)
                    if is_magnitude_phase:
                        mag = floats2[:, [1, 2, 3, 4, 5, 6, 7, 8]]
                        phase = floats2[:, [9, 10, 11, 12, 13, 14, 15, 16]]
                        rtheta = radians(phase)
                        real_imag = mag * (cos(rtheta) + 1.j * sin(rtheta))
                    else:
                        real = floats2[:, [1, 2, 3, 4, 5, 6, 7, 8]]
                        imag = floats2[:, [9, 10, 11, 12, 13, 14, 15, 16]]
                        real_imag = real + 1.j * imag
                    obj.data[obj.itime, itotal:itotal2, :] = real_imag
                    obj.itotal = itotal2
                    obj.ielement = ielement2
                else:
                    s1 = Struct(b(self._endian + 'i4s17f'))  # 2+17=19 * 4 = 76
                    s2 = Struct(b(self._endian + 'i16f'))  # 17 * 4 = 68
                    ntotal = 8 + (nnodes + 1) * 68
                    nelements = ndata // ntotal
                    obj = self.obj
                    for i in range(nelements):
                        edata = data[n:n + 76]
                        n += 76

                        out = s1.unpack(edata)
                        if self.is_debug_file:
                            self.binary_debug.write('OEF_Plate2-%s - %s\n' % (self.element_type, str(out)))
                        (eid_device, term, nid,
                         mxr, myr, mxyr, bmxr, bmyr, bmxyr, txr, tyr,
                         mxi, myi, mxyi, bmxi, bmyi, bmxyi, txi, tyi) = out
                        #term = 'CEN\'

                        eid = eid_device // 10
                        if is_magnitude_phase:
                            mx = polar_to_real_imag(mxr, mxi)
                            my = polar_to_real_imag(myr, myi)
                            mxy = polar_to_real_imag(mxyr, mxyi)
                            bmx = polar_to_real_imag(bmxr, bmxi)
                            bmy = polar_to_real_imag(bmyr, bmyi)
                            bmxy = polar_to_real_imag(bmxyr, bmxyi)
                            tx = polar_to_real_imag(txr, txi)
                            ty = polar_to_real_imag(tyr, tyi)
                        else:
                            mx = complex(mxr, mxi)
                            my = complex(myr, myi)
                            mxy = complex(mxyr, mxyi)
                            bmx = complex(bmxr, bmxi)
                            bmy = complex(bmyr, bmyi)
                            bmxy = complex(bmxyr, bmxyi)
                            tx = complex(txr, txi)
                            ty = complex(tyr, tyi)
                        obj.add_new_element_sort1(dt, eid, term, nid, mx, my, mxy, bmx, bmy, bmxy, tx, ty)

                        for i in range(nnodes):  # .. todo:: fix crash...
                            edata = data[n:n+68]
                            n += 68
                            out = s2.unpack(edata)
                            (nid,
                             mxr, myr, mxyr, bmxr, bmyr, bmxyr, txr, tyr,
                             mxi, myi, mxyi, bmxi, bmyi, bmxyi, txi, tyi) = out
                            if is_magnitude_phase:
                                mx = polar_to_real_imag(mxr, mxi)
                                my = polar_to_real_imag(myr, myi)
                                mxy = polar_to_real_imag(mxyr, mxyi)
                                bmx = polar_to_real_imag(bmxr, bmxi)
                                bmy = polar_to_real_imag(bmyr, bmyi)
                                bmxy = polar_to_real_imag(bmxyr, bmxyi)
                                tx = polar_to_real_imag(txr, txi)
                                ty = polar_to_real_imag(tyr, tyi)
                            else:
                                mx = complex(mxr, mxi)
                                my = complex(myr, myi)
                                mxy = complex(mxyr, mxyi)
                                bmx = complex(bmxr, bmxi)
                                bmy = complex(bmyr, bmyi)
                                bmxy = complex(bmxyr, bmxyi)
                                tx = complex(txr, txi)
                                ty = complex(tyr, tyi)
                            if self.is_debug_file:
                                self.binary_debug.write('OEF_Plate2 - eid=%i nid=%s out=%s\n' % (eid, nid, str(out)))
                            obj.add_sort1(dt, eid, nid, mx, my, mxy, bmx, bmy, bmxy, tx, ty)
            else:
                msg = self.code_information()
                return self._not_implemented_or_skip(data, ndata, msg)

        elif self.element_type in [95, 96, 97, 98]: # composites
            # 95 - CQUAD4
            # 96 - CQUAD8
            # 97 - CTRIA3
            # 98 - CTRIA6 (composite)
            if self.format_code == 1 and self.num_wide == 9:  # real
                if self.read_mode == 1:
                    return ndata
                ntotal = 36
                nelements = ndata // ntotal

                s1 = Struct(b('i8sifffif'))
                s2 = Struct(b('i8sifffff'))
                for i in range(nelements):
                    edata = data[n:n+ntotal]  # 4*9
                    out = s1.unpack(edata)

                    # i    8s              i        f
                    (eid, failure_theory, ply_id, failure_index_for_ply,
                     six, seven, flag, nine,
                     #failure_index_for_bonding,
                     #failure_index_for_element,
                     #flag,
                     #direct_stress_or_strain,
                     #interlaminar_stress,
                     #max_of_fb_fp_for_all_plies
                     ) = out
                    if flag != -1:
                        out = s2.unpack(edata)
                    #print(out)
                    n += 36


                ## TODO: add
                #return ndata
                #print self.code_information()
                #self.create_transient_object(self.compositePlateForces, RealCompositePlateForce)  # undefined
                ##return
                #ntotal = 9 * 4
                #nelements = ndata // ntotal
                #if self.is_debug_file:
                    #self.binary_debug.write('  [cap, element1, element2, ..., cap]\n')
                    #self.binary_debug.write('  cap = %i  # assume 1 cap when there could have been multiple\n' % ndata)
                    ##self.binary_debug.write('  #centeri = [eid_device, j, grid, fd1, sx1, sy1, txy1, angle1, major1, minor1, vm1,\n')
                    ##self.binary_debug.write('  #                                fd2, sx2, sy2, txy2, angle2, major2, minor2, vm2,)]\n')
                    ##self.binary_debug.write('  #nodeji = [eid, ilayer, o1, o2, t12, t1z, t2z, angle, major, minor, ovm)]\n')
                    #self.binary_debug.write('  nelements=%i; nnodes=1 # centroid\n' % nelements)

                #eid_old = 0
                #s = Struct(b(self._endian + 'i8si4f4s'))
                #for i in range(nelements):
                    #if i % 10000 == 0:
                        #print 'i = ', i
                    #edata = data[n:n+ntotal]  # 4*9
                    #out = s.unpack(edata)
                    #(eid_device, theory, lamid, failure_index_direct_stress, failure_mode_max_shear,
                             #failure_index_interlaminar_shear, fmax, failure_flag) = out
                    #eid = eid_device // 10
                    #if self.is_debug_file:
                        #if eid > 0:
                            #self.binary_debug.write('  eid=%i; C=[%s]\n' % (', '.join(['%r' % di for di in out]) ))
                        #else:
                            #self.binary_debug.write('      %s  C=[%s]\n' % (' ' * len(str(eid)), ', '.join(['%r' % di for di in out]) ))

                    #if eid > 0:
                        #obj.add_new_eid(eType, dt, eid, o1, o2, t12, t1z, t2z, angle, major, minor, ovm)
                    #else:
                        #obj.add(dt, eid, o1, o2, t12, t1z, t2z, angle, major, minor, ovm)
                    #eid_old = eid
                    #n += ntotal
            else:
                msg = self.code_information()
                return self._not_implemented_or_skip(data, ndata, msg)

        elif self.element_type in [39, 67, 68]: # solids
            # 39-CTETRA
            # 67-CHEXA
            # 68-CPENTA
            if self.read_mode == 1:
                return ndata
            #self._results._found_result('solid_forces')
            raise RuntimeError(self.code_information())
            if self.format_code == 1 and self.num_wide == 0:  # real
                #self.create_transient_object(self.solidForces, RealCSolidForce)
                raise RuntimeError(self.code_information())
            else:
                msg = self.code_information()
                return self._not_implemented_or_skip(data, ndata, msg)

        elif self.element_type == 53:  # ctriax6
            # 53-CTRIAX6
            if self.read_mode == 1:
                return ndata
            self._results._found_result('ctriax_force')
            if self.format_code == 1 and self.num_wide == 0:  # real
                pass
                #self.create_transient_object(self.ctriax_force, RealCTriaxForce)  # undefined
            else:
                msg = self.code_information()
                return self._not_implemented_or_skip(data, ndata, msg)

        elif self.element_type == 4:  # cshear
            # 4-CSHEAR
            result_name = 'cshear_force'
            if self._results.is_not_saved(result_name):
                return ndata
            self._results._found_result(result_name)

            slot = getattr(self, result_name)
            if self.format_code == 1 and self.num_wide == 17:  # real
                ntotal = 68  # 17*4
                nelements = ndata // ntotal

                obj_real = RealCShearForceArray
                auto_return, is_vectorized = self._create_oes_object4(
                    nelements, result_name, slot, obj_real)
                if auto_return:
                    return nelements * self.num_wide * 4

                obj = self.obj
                if self.use_vector and is_vectorized:
                    n = nelements * 4 * self.num_wide
                    itotal = obj.ielement
                    ielement2 = obj.itotal + nelements
                    itotal2 = ielement2

                    floats = fromstring(data, dtype=self.fdtype).reshape(nelements, 17)
                    obj._times[obj.itime] = dt
                    if obj.itime == 0:
                        ints = fromstring(data, dtype=self.idtype).reshape(nelements, 17)
                        eids = ints[:, 0] // 10
                        assert eids.min() > 0, eids.min()
                        obj.element[itotal:itotal2] = eids

                    # [f41, f21, f12, f32, f23, f43, f34, f14, kf1,
                    #  s12, kf2, s23, kf3, s34, kf4, s41]
                    obj.data[obj.itime, itotal:itotal2, :] = floats[:, 1:]
                    obj.itotal = itotal2
                    obj.ielement = ielement2
                else:
                    s = Struct(b(self._endian + 'i16f'))
                    for i in range(nelements):
                        edata = data[n:n+68]

                        out = s.unpack(edata)
                        if self.is_debug_file:
                            self.binary_debug.write('OEF_Shear - %s\n' % (str(out)))
                        (eid_device,
                         f41, f21, f12, f32, f23, f43, f34, f14, kf1,
                         s12, kf2, s23, kf3, s34, kf4, s41) = out
                        eid = eid_device // 10

                        #data_in = [eid,
                                   #f41, f21, f12, f32, f23, f43, f34,
                                   #f14, kf1, s12, kf2, s23, kf3, s34, kf4, s41]
                        #print "%s" % (self.get_element_type(self.element_type)), data_in
                        obj.add(dt, eid,
                                f41, f21, f12, f32, f23, f43, f34,
                                f14, kf1, s12, kf2, s23, kf3, s34, kf4, s41)
                        n += ntotal

            elif self.format_code in [2, 3] and self.num_wide == 33:  # imag
                ntotal = 132  # 33*4
                nelements = ndata // ntotal

                obj_complex = ComplexCShearForceArray
                auto_return, is_vectorized = self._create_oes_object4(
                    nelements, result_name, slot, obj_complex)
                if auto_return:
                    return nelements * self.num_wide * 4

                obj = self.obj
                #if self.is_debug_file:
                    #self.binary_debug.write('  [cap, element1, element2, ..., cap]\n')
                    #self.binary_debug.write('  cap = %i  # assume 1 cap when there could have been multiple\n' % ndata)
                    #self.binary_debug.write('  #elementi = [eid_device, axial, torque]\n')
                    #self.binary_debug.write('  nelements=%i; nnodes=1 # centroid\n' % nelements)

                if self.use_vector and is_vectorized:
                    n = nelements * 4 * self.num_wide
                    itotal = obj.ielement
                    ielement2 = obj.itotal + nelements
                    itotal2 = ielement2

                    floats = fromstring(data, dtype=self.fdtype).reshape(nelements, 33)
                    obj._times[obj.itime] = dt
                    if obj.itime == 0:
                        ints = fromstring(data, dtype=self.idtype).reshape(nelements, 33)
                        eids = ints[:, 0] // 10
                        assert eids.min() > 0, eids.min()
                        obj.element[itotal:itotal2] = eids

                    #[f41r, f21r, f12r, f32r, f23r, f43r, f34r, f14r
                    # kf1r, s12r, kf2r, s23r, kf3r, s34r, kf4r, s41r
                    # f41i, f21i, f12i, f32i, f23i, f43i, f34i, f14i
                    # kf1i, s12i, kf2i, s23i, kf3i, s34i, kf4i, s41i]
                    if is_magnitude_phase:
                        mag = floats[:, [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16]]
                        phase = floats[:, [17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32]]
                        rtheta = radians(phase)
                        real_imag = mag * (cos(rtheta) + 1.j * sin(rtheta))
                    else:
                        real = floats[:, [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16]]
                        imag = floats[:, [17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32]]
                        real_imag = real + 1.j * imag
                    obj.data[obj.itime, itotal:itotal2, :] = real_imag

                    #obj.data[obj.itime, itotal:itotal2, :] = floats[:, 1:]
                    obj.itotal = itotal2
                    obj.ielement = ielement2
                else:
                    #self.create_transient_object(self.cshear_force, ComplexCShearForce)
                    s = Struct(b(self._endian + 'i32f'))

                    for i in range(nelements):
                        edata = data[n:n+132]
                        n += ntotal
                        out = s.unpack(edata)
                        if self.is_debug_file:
                            self.binary_debug.write('OEF_Shear - %s\n' % (str(out)))
                        (eid_device,
                         f41r, f21r, f12r, f32r, f23r, f43r, f34r, f14r, # 8
                         kf1r, s12r, kf2r, s23r, kf3r, s34r, kf4r, s41r, # 16
                         f41i, f21i, f12i, f32i, f23i, f43i, f34i, f14i,
                         kf1i, s12i, kf2i, s23i, kf3i, s34i, kf4i, s41i) = out
                        if is_magnitude_phase:
                            f41r = polar_to_real_imag(f41r, f41i)
                            kf1 = polar_to_real_imag(kf1r, kf1i)
                            f21r = polar_to_real_imag(f21r, f21i)
                            kf2 = polar_to_real_imag(kf2r, kf2i)
                            f12r = polar_to_real_imag(f12r, f12i)
                            kf3 = polar_to_real_imag(kf3r, kf3i)
                            f23r = polar_to_real_imag(f23r, f23i)
                            kf4 = polar_to_real_imag(kf4r, kf4i)
                            f32r = polar_to_real_imag(f32r, f32i)
                            s12 = polar_to_real_imag(s12r, s12i)
                            f43r = polar_to_real_imag(f43r, f43i)
                            s23 = polar_to_real_imag(s23r, s23i)
                            f34r = polar_to_real_imag(f34r, f34i)
                            s34 = polar_to_real_imag(s34r, s34i)
                            f14r = polar_to_real_imag(f14r, f14i)
                            s41 = polar_to_real_imag(s41r, s41i)
                        else:
                            f41 = complex(f41r, f41i)
                            kf1 = complex(kf1r, kf1i)
                            f21 = complex(f21r, f21i)
                            kf2 = complex(kf2r, kf2i)
                            f12 = complex(f12r, f12i)
                            kf3 = complex(kf3r, kf3i)
                            f23 = complex(f23r, f23i)
                            kf4 = complex(kf4r, kf4i)
                            f32 = complex(f32r, f32i)
                            s12 = complex(s12r, s12i)
                            f43 = complex(f43r, f43i)
                            s23 = complex(s23r, s23i)
                            f34 = complex(f34r, f34i)
                            s34 = complex(s34r, s34i)
                            f14 = complex(f14r, f14i)
                            s41 = complex(s41r, s41i)

                        eid = eid_device // 10
                        obj.add_sort1(dt, eid,
                                      f41, f21, f12, f32, f23, f43, f34, f14,
                                      kf1, s12, kf2, s23, kf3, s34, kf4, s41)
            else:
                msg = self.code_information()
                return self._not_implemented_or_skip(data, ndata, msg)

        elif self.element_type == 35:  # coneax
            # 35-CONEAX
            result_name = 'coneax_force'
            if self._results.is_not_saved(result_name):
                return ndata
            self._results._found_result(result_name)

            slot = getattr(self, result_name)
            if self.format_code == 1 and self.num_wide == 7:  # real
                ntotal = 28  # 7*4
                nelements = ndata // ntotal

                auto_return, is_vectorized = self._create_oes_object4(
                    nelements, result_name, slot, RealConeAxForceArray)
                if auto_return:
                    return nelements * self.num_wide * 4

                obj = self.obj
                if self.use_vector and is_vectorized:
                    n = nelements * 4 * self.num_wide
                    itotal = obj.ielement
                    ielement2 = obj.itotal + nelements
                    itotal2 = ielement2

                    floats = fromstring(data, dtype=self.fdtype).reshape(nelements, 7)
                    obj._times[obj.itime] = dt
                    if obj.itime == 0:
                        ints = fromstring(data, dtype=self.idtype).reshape(nelements, 7)
                        eids = ints[:, 0] // 10
                        assert eids.min() > 0, eids.min()
                        obj.element[itotal:itotal2] = eids

                    # [hopa, bmu, bmv, tm, su, sv]
                    obj.data[obj.itime, itotal:itotal2, :] = floats[:, 1:]
                    obj.itotal = itotal2
                    obj.ielement = ielement2
                else:
                    s = Struct(b(self._endian + 'i6f'))
                    for i in range(nelements):
                        edata = data[n:n+ntotal]
                        out = s.unpack(edata)
                        if self.is_debug_file:
                            self.binary_debug.write('OEF_CONEAX-35 - %s\n' % (str(out)))
                        (eid_device, hopa, bmu, bmv, tm, su, sv) = out
                        eid = eid_device // 10
                        obj.add(dt, eid, hopa, bmu, bmv, tm, su, sv)
                        n += ntotal
            else:
                msg = self.code_information()
                return self._not_implemented_or_skip(data, ndata, msg)


        elif self.element_type == 38:  # cgap
            # 38-GAP
            result_name = 'cgap_force'
            self._results._found_result(result_name)
            slot = getattr(self, result_name)
            if self.format_code == 1 and self.num_wide == 9:  # real
                ntotal = 36 # 9*4
                nelements = ndata // ntotal
                obj = self.obj

                auto_return, is_vectorized = self._create_oes_object4(
                    nelements, result_name, slot, RealCGapForceArray)
                if auto_return:
                    return nelements * self.num_wide * 4

                obj = self.obj
                if self.use_vector and is_vectorized:
                    n = nelements * 4 * self.num_wide
                    itotal = obj.ielement
                    ielement2 = obj.itotal + nelements
                    itotal2 = ielement2

                    floats = fromstring(data, dtype=self.fdtype).reshape(nelements, 9)
                    obj._times[obj.itime] = dt
                    if obj.itime == 0:
                        ints = fromstring(data, dtype=self.idtype).reshape(nelements, 9)
                        eids = ints[:, 0] // 10
                        assert eids.min() > 0, eids.min()
                        obj.element[itotal:itotal2] = eids

                    # [fx, sfy, sfz, u, v, w, sv, sw]
                    obj.data[obj.itime, itotal:itotal2, :] = floats[:, 1:]
                    obj.itotal = itotal2
                    obj.ielement = ielement2
                else:
                    s = Struct(b(self._endian + 'i8f'))
                    for i in range(nelements):
                        edata = data[n:n+36]

                        out = s.unpack(edata)
                        if self.is_debug_file:
                            self.binary_debug.write('OEF_CGAP-38 - %s\n' % (str(out)))
                        (eid_device, fx, sfy, sfz, u, v, w, sv, sw) = out
                        eid = eid_device // 10
                        data_in = [eid, fx, sfy, sfz, u, v, w, sv, sw]
                        #print "%s" %(self.get_element_type(self.element_type)), data_in
                        #eid = obj.add_new_eid(out)
                        obj.add(dt, eid, fx, sfy, sfz, u, v, w, sv, sw)
                        n += ntotal
            else:
                msg = self.code_information()
                return self._not_implemented_or_skip(data, ndata, msg)
        elif self.element_type == 69:  # cbend
            # 69-CBEND
            result_name = 'cbend_force'
            self._results._found_result(result_name)
            slot = getattr(self, result_name)
            if self.format_code == 1 and self.num_wide == 15:  # real
                ntotal = 60  # 15*4
                nelements = ndata // ntotal
                auto_return, is_vectorized = self._create_oes_object4(
                    nelements, result_name, slot, RealBendForceArray)
                if auto_return:
                    return nelements * self.num_wide * 4

                obj = self.obj
                if self.use_vector and is_vectorized:
                    n = nelements * 4 * self.num_wide
                    itotal = obj.ielement
                    ielement2 = obj.itotal + nelements
                    itotal2 = ielement2

                    floats = fromstring(data, dtype=self.fdtype).reshape(nelements, 15)
                    obj._times[obj.itime] = dt
                    if obj.itime == 0:
                        ints = fromstring(data, dtype=self.idtype).reshape(nelements, 15)
                        eids = ints[:, 0] // 10
                        nids_a = ints[:, 1]
                        nids_b = ints[:, 8]
                        assert eids.min() > 0, eids.min()
                        obj.element_nodes[itotal:itotal2, 0] = eids
                        obj.element_nodes[itotal:itotal2, 1] = nids_a
                        obj.element_nodes[itotal:itotal2, 2] = nids_b

                    # [nidA, bm1A, bm2A, ts1A, ts2A, afA, trqA,
                    #  nidB, bm1B, bm2B, ts1B, ts2B, afB, trqB]
                    assert floats[:, 2:8].shape[1] == 6, floats[:, 2:8].shape
                    assert floats[:, 9:].shape[1] == 6, floats[:, 9:].shape
                    obj.data[obj.itime, itotal:itotal2, :6] = floats[:, 2:8]
                    obj.data[obj.itime, itotal:itotal2, 6:] = floats[:, 9:]
                    obj.itotal = itotal2
                    obj.ielement = ielement2
                else:
                    s = Struct(b(self._endian + 'i i6fi6f'))
                    for i in range(nelements):
                        edata = data[n:n+ntotal]

                        out = s.unpack(edata)
                        if self.is_debug_file:
                            self.binary_debug.write('OEF_BEND-69 - %s\n' % (str(out)))
                        (eid_device,
                         nidA, bm1A, bm2A, ts1A, ts2A, afA, trqA,
                         nidB, bm1B, bm2B, ts1B, ts2B, afB, trqB) = out
                        eid = eid_device // 10

                        obj.add_sort1(
                            dt, eid,
                            nidA, bm1A, bm2A, ts1A, ts2A, afA, trqA,
                            nidB, bm1B, bm2B, ts1B, ts2B, afB, trqB)
                        n += ntotal
            elif self.format_code in [2, 3] and self.num_wide == 27:  # imag
                # TODO: vectorize
                ntotal = 108  # 27*4
                nelements = ndata // ntotal

                auto_return, is_vectorized = self._create_oes_object4(
                    nelements, result_name, slot, ComplexCBendForceArray)
                if auto_return:
                    return nelements * self.num_wide * 4

                obj = self.obj
                if self.use_vector and is_vectorized:
                    n = nelements * 4 * self.num_wide
                    itotal = obj.ielement
                    ielement2 = obj.itotal + nelements
                    itotal2 = ielement2

                    ireal = [2, 3,  4,  5,  6,  7, 15, 16, 17, 18, 19, 20]
                    iimag = [8, 9, 10, 11, 12, 13, 21, 22, 23, 24, 25, 26]
                    # 0    1
                    # eid, nidA,
                    # 2      3      4       5     6     7
                    # 8      9      10      11    12    13
                    # bm1Ar, bm2Ar, ts1Ar, ts2Ar, afAr, trqAr,
                    # bm1Ai, bm2Ai, ts1Ai, ts2Ai, afAi, trqAi,
                    # 14
                    # nidB
                    # 15     16     17     18     19    20
                    # 21     22     23     24     25    26
                    # bm1Br, bm2Br, ts1Br, ts2Br, afBr, trqBr,
                    # bm1Bi, bm2Bi, ts1Bi, ts2Bi, afBi, trqBi
                    floats = fromstring(data, dtype=self.fdtype).reshape(nelements, 27)
                    obj._times[obj.itime] = dt
                    if obj.itime == 0:
                        ints = fromstring(data, dtype=self.idtype).reshape(nelements, 27)
                        eids = ints[:, 0] // 10
                        nids_a = ints[:, 1]
                        nids_b = ints[:, 14]
                        assert nids_a.min() > 0, nids_a
                        assert nids_b.min() > 0, nids_b
                        assert eids.min() > 0, eids.min()
                        #print(nids_b)
                        obj.element_nodes[itotal:itotal2, 0] = eids
                        obj.element_nodes[itotal:itotal2, 1] = nids_a
                        obj.element_nodes[itotal:itotal2, 2] = nids_b

                    if is_magnitude_phase:
                        mag = floats[:, ireal]
                        phase = floats[:, iimag]
                        rtheta = radians(phase)
                        real_imag = mag * (cos(rtheta) + 1.j * sin(rtheta))
                    else:
                        real = floats[:, ireal]
                        imag = floats[:, iimag]
                        real_imag = real + 1.j * imag
                    obj.data[obj.itime, itotal:itotal2, :] = real_imag
                    #assert floats[:, 1:6].shape[1] == 12, floats[:, 1:6].shape
                    #assert floats[:, 7:].shape[1] == 6, floats[:, 7:].shape
                    #obj.data[obj.itime, itotal:itotal2, :6] = floats[:, 1:6]
                    #obj.data[obj.itime, itotal:itotal2, 6:] = floats[:, 7:]
                    obj.itotal = itotal2
                    obj.ielement = ielement2
                else:
                    s = Struct(b(self._endian + 'i i12f i12f'))
                    for i in range(nelements):
                        edata = data[n:n+108]
                        n += ntotal
                        out = s.unpack(edata)
                        if self.is_debug_file:
                            self.binary_debug.write('OEF_BEND-69 - %s\n' % (str(out)))
                        (eid_device, nidA,
                         bm1Ar, bm2Ar, ts1Ar, ts2Ar, afAr, trqAr,
                         bm1Ai, bm2Ai, ts1Ai, ts2Ai, afAi, trqAi,
                         nidB,
                         bm1Br, bm2Br, ts1Br, ts2Br, afBr, trqBr,
                         bm1Bi, bm2Bi, ts1Bi, ts2Bi, afBi, trqBi) = out
                        eid = eid_device // 10

                        if is_magnitude_phase:
                            bm1A = polar_to_real_imag(bm1Ar, bm1Ai)
                            bm1B = polar_to_real_imag(bm1Br, bm1Bi)
                            bm2A = polar_to_real_imag(bm2Ar, bm2Ai)
                            bm2B = polar_to_real_imag(bm2Br, bm2Bi)
                            ts1A = polar_to_real_imag(ts1Ar, ts1Ai)
                            ts1B = polar_to_real_imag(ts1Br, ts1Bi)
                            ts2A = polar_to_real_imag(ts2Ar, ts2Ai)
                            ts2B = polar_to_real_imag(ts2Br, ts2Bi)
                            afA = polar_to_real_imag(afAr, afAi)
                            afB = polar_to_real_imag(afBr, afBi)
                            trqA = polar_to_real_imag(trqAr, trqAi)
                            trqB = polar_to_real_imag(trqBr, trqBi)
                        else:
                            bm1A = complex(bm1Ar, bm1Ai)
                            bm1B = complex(bm1Br, bm1Bi)
                            bm2A = complex(bm2Ar, bm2Ai)
                            bm2B = complex(bm2Br, bm2Bi)
                            ts1A = complex(ts1Ar, ts1Ai)
                            ts1B = complex(ts1Br, ts1Bi)
                            ts2A = complex(ts2Ar, ts2Ai)
                            ts2B = complex(ts2Br, ts2Bi)
                            afA = complex(afAr, afAi)
                            afB = complex(afBr, afBi)
                            trqA = complex(trqAr, trqAi)
                            trqB = complex(trqBr, trqBi)

                        obj.add_sort1(dt, eid,
                                nidA, bm1A, bm2A, ts1A, ts2A, afA, trqA,
                                nidB, bm1B, bm2B, ts1B, ts2B, afB, trqB)
            else:
                msg = self.code_information()
                return self._not_implemented_or_skip(data, ndata, msg)

        elif self.element_type in [76, 77, 78]:
            # 76-HEXPR
            # 77-PENPR
            # 78-TETPR
            if self.element_type == 76:
                result_name = 'chexa_pressure_force'
                slot = self.chexa_pressure_force
            elif self.element_type == 77:
                result_name = 'cpenta_pressure_force'
                slot = self.cpenta_pressure_force
            elif self.element_type == 77:
                result_name = 'ctetra_pressure_force'
                slot = self.ctetra_pressure_force
            else:
                msg = self.code_information()
                return self._not_implemented_or_skip(data, ndata, msg)

            self._results._found_result(result_name)
            if self.format_code == 1 and self.num_wide == 10:  # real
                ntotal = 40
                nelements = ndata // ntotal
                nelements = ndata // ntotal

                obj_real = RealSolidPressureForceArray
                auto_return, is_vectorized = self._create_oes_object4(
                    nelements, result_name, slot, obj_real)
                if auto_return:
                    return nelements * self.num_wide * 4

                obj = self.obj
                #if self.is_debug_file:
                    #self.binary_debug.write('  [cap, element1, element2, ..., cap]\n')
                    #self.binary_debug.write('  cap = %i  # assume 1 cap when there could have been multiple\n' % ndata)
                    #self.binary_debug.write('  #elementi = [eid_device, axial, torque]\n')
                    #self.binary_debug.write('  nelements=%i; nnodes=1 # centroid\n' % nelements)
                if self.use_vector and is_vectorized:
                    # self.itime = 0
                    # self.ielement = 0
                    # self.itotal = 0
                    #self.ntimes = 0
                    #self.nelements = 0
                    n = nelements * self.num_wide * 4
                    itotal = obj.ielement
                    ielement2 = obj.itotal + nelements
                    itotal2 = ielement2

                    floats = fromstring(data, dtype=self.fdtype).reshape(nelements, 10)
                    obj._times[obj.itime] = dt
                    if obj.itime == 0:
                        ints = fromstring(data, dtype=self.idtype).reshape(nelements, 10)
                        eids = ints[:, 0] // 10
                        assert eids.min() > 0, eids.min()
                        obj.element[itotal:itotal2] = eids

                    #[axial_force, torque]
                    obj.data[obj.itime, itotal:itotal2, :] = floats[:, 3:]
                    obj.itotal = itotal2
                    obj.ielement = ielement2
                else:
                    s = Struct(b(self._endian + 'i8s7f'))
                    for i in range(nelements):
                        edata = data[n : n + 40]
                        n += 40
                        out = s.unpack(edata)
                        if self.is_debug_file:
                            self.binary_debug.write('OEF_PentaPressure-%s %s\n' % (self.element_type, str(out)))
                        (eid_device, ename, ax, ay, az, vx, vy, vz, pressure) = out
                        eid = eid_device // 10
                        obj.add_sort1(dt, eid, ename, ax, ay, az, vx, vy, vz, pressure)

            elif self.format_code in [2, 3] and self.num_wide == 16:  # imag
                ntotal = 64
                nelements = ndata // ntotal
                auto_return, is_vectorized = self._create_oes_object4(
                    nelements, result_name, slot, ComplexSolidPressureForceArray)
                if auto_return:
                    return nelements * self.num_wide * 4

                obj = self.obj
                #if self.is_debug_file:
                    #self.binary_debug.write('  [cap, element1, element2, ..., cap]\n')
                    #self.binary_debug.write('  cap = %i  # assume 1 cap when there could have been multiple\n' % ndata)
                    #self.binary_debug.write('  #elementi = [eid_device, axial, torque]\n')
                    #self.binary_debug.write('  nelements=%i; nnodes=1 # centroid\n' % nelements)

                if self.use_vector and is_vectorized:
                    n = nelements * 4 * self.num_wide
                    itotal = obj.ielement
                    ielement2 = obj.itotal + nelements
                    itotal2 = ielement2

                    floats = fromstring(data, dtype=self.fdtype).reshape(nelements, 16)
                    obj._times[obj.itime] = dt
                    if obj.itime == 0:
                        ints = fromstring(data, dtype=self.idtype).reshape(nelements, 16)
                        eids = ints[:, 0] // 10
                        assert eids.min() > 0, eids.min()
                        obj.element[itotal:itotal2] = eids

                    #[xaccr, yaccr, zaccr, xvelr, yvelr, zvelr, pressure,
                    # xacci, yacci, zacci, xveli, yveli, zveli]
                    if is_magnitude_phase:
                        mag = floats[:, [3, 4, 5, 6, 7, 8, 9]]
                        phase = hstack([
                            floats[:, [10, 11, 12, 13, 14, 15]],
                            zeros((len(floats), 1), dtype='float32')
                        ])
                        rtheta = radians(phase)
                        real_imag = mag * (cos(rtheta) + 1.j * sin(rtheta))
                    else:
                        real = floats[:, [3, 4, 5, 6, 7, 8, 9]]
                        imag = hstack([
                            floats[:, [10, 11, 12, 13, 14, 15]],
                            zeros((len(floats), 1), dtype='float32')
                        ])
                        real_imag = real + 1.j * imag
                    obj.data[obj.itime, itotal:itotal2, :] = real_imag
                    obj.itotal = itotal2
                    obj.ielement = ielement2
                else:
                    s = Struct(b(self._endian + 'i8s13f'))
                    for i in range(nelements):
                        edata = data[n:n+64]
                        n += 64

                        out = s.unpack(edata)
                        if self.is_debug_file:
                            self.binary_debug.write('OEF_PentaPressure-%s %s\n' % (self.element_type, str(out)))
                        (eid_device, ename,
                         axr, ayr, azr, vxr, vyr, vzr, pressure,
                         axi, ayi, azi, vxi, vyi, vzi) = out
                        eid = eid_device // 10
                        ename = ename.decode('utf-8').strip()

                        if is_magnitude_phase:
                            ax = polar_to_real_imag(axr, axi)
                            vx = polar_to_real_imag(vxr, vxi)
                            ay = polar_to_real_imag(ayr, ayi)
                            vy = polar_to_real_imag(vyr, vyi)
                            az = polar_to_real_imag(azr, azi)
                            vz = polar_to_real_imag(vzr, vzi)
                        else:
                            ax = complex(axr, axi)
                            vx = complex(vxr, vxi)
                            ay = complex(ayr, ayi)
                            vy = complex(vyr, vyi)
                            az = complex(azr, azi)
                            vz = complex(vzr, vzi)
                        cpressure = complex(pressure, 0.)
                        obj.add_sort1(dt, eid, ename, ax, ay, az, vx, vy, vz, cpressure)
            else:
                msg = self.code_information()
                return self._not_implemented_or_skip(data, ndata, msg)

        elif self.element_type == 102:  # cbush
            # 102-CBUSH
            self._results._found_result('cbush_force')
            result_name = 'cbush_force'
            slot = getattr(self, result_name)

            numwide_real = 7
            if self.format_code == 1 and self.num_wide == 7:  # real
                ntotal = 28 # 7*4
                nelements = ndata // ntotal

                auto_return, is_vectorized = self._create_oes_object4(
                    nelements, result_name, slot, RealCBushForceArray)
                if auto_return:
                    return nelements * self.num_wide * 4

                obj = self.obj
                if self.use_vector and is_vectorized:
                    # self.itime = 0
                    # self.ielement = 0
                    # self.itotal = 0
                    #self.ntimes = 0
                    #self.nelements = 0
                    n = nelements * self.num_wide * 4

                    istart = obj.itotal
                    iend = istart + nelements
                    obj._times[obj.itime] = dt

                    if obj.itime == 0:
                        ints = fromstring(data, dtype=self.idtype).reshape(nelements, numwide_real)
                        eids = ints[:, 0] // 10
                        obj.element[istart:iend] = eids
                    results = fromstring(data, dtype=self.fdtype).reshape(nelements, numwide_real)

                    #[fx, fy, fz, mx, my, mz]
                    obj.data[obj.itime, istart:iend, :] = results[:, 1:]
                else:
                    s = Struct(b(self._endian + 'i6f'))
                    for i in range(nelements):
                        edata = data[n:n+28]
                        out = s.unpack(edata)
                        if self.is_debug_file:
                            self.binary_debug.write('OEF_CBUSH-102 - %s\n' % (str(out)))
                        (eid_device, fx, fy, fz, mx, my, mz) = out
                        eid = eid_device // 10
                        obj.add(dt, eid, fx, fy, fz, mx, my, mz)
                        n += ntotal
            elif self.format_code in [2, 3] and self.num_wide == 13:  # imag
                # TODO: vectorize
                ntotal = 52  # 13*4
                nelements = ndata // ntotal
                result_name = 'cbush_force'
                auto_return, is_vectorized = self._create_oes_object4(
                    nelements, result_name, slot, ComplexCBushForceArray)
                if auto_return:
                    return nelements * self.num_wide * 4

                s = Struct(b(self._endian + 'i12f'))

                obj = self.obj
                for i in range(nelements):
                    edata = data[n:n + 52]

                    out = s.unpack(edata)
                    if self.is_debug_file:
                        self.binary_debug.write('OEF_CBUSH-102 - %s\n' % (str(out)))
                    (eid_device,
                     fxr, fyr, fzr, mxr, myr, mzr,
                     fxi, fyi, fzi, mxi, myi, mzi) = out
                    eid = eid_device // 10

                    if is_magnitude_phase:
                        fx = polar_to_real_imag(fxr, fxi)
                        mx = polar_to_real_imag(mxr, mxi)
                        fy = polar_to_real_imag(fyr, fyi)
                        my = polar_to_real_imag(myr, myi)
                        fz = polar_to_real_imag(fzr, fzi)
                        mz = polar_to_real_imag(mzr, mzi)
                    else:
                        fx = complex(fxr, fxi)
                        mx = complex(mxr, mxi)
                        fy = complex(fyr, fyi)
                        my = complex(myr, myi)
                        fz = complex(fzr, fzi)
                        mz = complex(mzr, mzi)

                    obj.add_sort1(dt, eid, fx, fy, fz, mx, my, mz)
                    n += ntotal
            else:
                msg = self.code_information()
                return self._not_implemented_or_skip(data, ndata, msg)

        elif self.element_type in [145, 146, 147]:
            # 145-VUHEXA
            # 146-VUPENTA
            # 147-VUTETRA
            if self.read_mode == 1:
                return ndata
            return ndata
        elif self.element_type in [189, 190]:
            # 189-VUQUAD
            # 190-VUTRIA
            if self.read_mode == 1:
                return ndata

            self._results._found_result('force_VU_2D')
            if self.element_type in [189]:  # VUQUAD
                nnodes = 4
                etype = 'VUQUAD4'
            elif self.element_type in [190]:  # VUTRIA
                nnodes = 3
                etype = 'VUTRIA3'
            else:
                raise NotImplementedError(self.code_information())
            numwide_real = 6 + 13 * nnodes
            numwide_imag = 6 + 25 * nnodes

            if self.format_code == 1 and self.num_wide == numwide_real:  # real
                # TODO: vectorize
                self.create_transient_object(self.force_VU_2D, RealForce_VU_2D)

                ntotal = 24 + 52 * nnodes
                nelements = ndata // ntotal
                obj = self.obj

                s1 = Struct(b(self._endian + '3i4s2i'))
                s2 = Struct(b(self._endian + 'i3f3i5fi'))
                for i in range(nelements):
                    edata = data[n:n+24]  # 6*4
                    n += 24

                    out = s1.unpack(edata)
                    if self.is_debug_file:
                        self.binary_debug.write('OEF_Force_%s-%s - %s\n' % (etype, self.element_type, str(out)))
                    (eid_device, parent, coord, icord, theta, _) = out

                    eid = eid_device // 10
                    data_in = [eid, parent, coord, icord, theta]

                    forces = []
                    for i in range(nnodes):
                        edata = data[n:n+52]  # 13*4
                        n += 52
                        out = s2.unpack(edata)
                        if self.is_debug_file:
                            self.binary_debug.write('%s\n' % (str(out)))
                        (vugrid, mfx, mfy, mfxy, ai, bi, ci, bmx, bmy,
                         bmxy, syz, szx, di) = out
                        out2 = (vugrid, mfx, mfy, mfxy, bmx, bmy, bmxy, syz, szx)
                        forces.append(out2)
                    data_in.append(forces)
                    #data_in = [vugrid,mfx,mfy,mfxy,a,b,c,bmx,bmy,bmxy,syz,szx,d]
                    obj.add(nnodes, dt, data_in)

            elif self.format_code in [2, 3] and self.num_wide == numwide_imag:  # imag
                # TODO: vectorize
                self.create_transient_object(self.force_VU_2D, ComplexForce_VU_2D)
                obj = self.obj
                ntotal = 24 + 100 * nnodes
                s1 = Struct(b(self._endian + 'iii4sii'))
                s2 = Struct(b(self._endian + 'i3f3i5fi3f3i5fi'))
                nelements = ndata // ntotal
                for i in range(nelements):
                    edata = data[n:n+24]  # 6*4
                    n += 24

                    out = s1.unpack(edata)
                    if self.is_debug_file:
                        self.binary_debug.write('OEF_Force_%s-%s - %s\n' % (etype, self.element_type, str(out)))
                    (eid_device, parent, coord, icord, theta, _) = out

                    eid = eid_device // 10
                    data_in = [eid, parent, coord, icord, theta]

                    forces = []
                    for i in range(nnodes):
                        edata = data[n:n+100]  # 13*4
                        n += 100
                        out = s2.unpack(edata)
                        if self.is_debug_file:
                            self.binary_debug.write('OEF_Force_%s-%s - %s\n' % (etype, self.element_type, str(out)))
                        [vugrid, mfxr, mfyr, mfxyr, ai, bi, ci, bmxr, bmyr, bmxyr, syzr, szxr, di,
                         mfxi, mfyi, mfxyi, ai, bi, ci, bmxi, bmyi, bmxyi, syzi, szxi, di] = out

                        if is_magnitude_phase:
                            mfx = polar_to_real_imag(mfxr, mfxi)
                            mfy = polar_to_real_imag(mfyr, mfyi)
                            mfxy = polar_to_real_imag(mfxyr, mfxyi)
                            bmx = polar_to_real_imag(bmxr, bmxi)
                            bmy = polar_to_real_imag(bmyr, bmyi)
                            bmxy = polar_to_real_imag(bmxyr, bmxyi)
                            syz = polar_to_real_imag(syzr, syzi)
                            szx = polar_to_real_imag(szxr, szxi)
                        else:
                            mfx = complex(mfxr, mfxi)
                            mfy = complex(mfyr, mfyi)
                            mfxy = complex(mfxyr, mfxyi)
                            bmx = complex(bmxr, bmxi)
                            bmy = complex(bmyr, bmyi)
                            bmxy = complex(bmxyr, bmxyi)
                            syz = complex(syzr, syzi)
                            szx = complex(szxr, szxi)
                        out2 = [vugrid, mfx, mfy, mfxy, bmx, bmy, bmxy, syz, szx]
                        forces.append(out2)

                    data_in.append(forces)
                    #data_in = [vugrid,mfxr,mfyr,mfxyr,bmxr,bmyr,bmxyr,syzr,szxr,
                                     #mfxi,mfyi,mfxyi,bmxi,bmyi,bmxyi,syzi,szxi]
                    obj.add(nnodes, dt, data_in)
            else:
                msg = self.code_information()
                return self._not_implemented_or_skip(data, ndata, msg)


        elif self.element_type == 191:
            # 191-VUBEAM
            if self.read_mode == 1:
                return ndata
            self._results._found_result('force_VU')

            if self.format_code == 1 and self.num_wide == 20:  # real
                # TODO: vectorize
                self.create_transient_object(self.force_VU, RealForce_VU)
                # 20 = 4 + 8 * 2 = 4 = 16
                nnodes = 2
                #ntotal = 16 + 32 * nnodes
                ntotal = self.num_wide * 4
                nelements = ndata // ntotal
                self.create_transient_object(self.force_VU, RealForce_VU)
                obj = self.obj

                s1 = Struct(b(self._endian + 'iii4s'))
                s2 = Struct(b(self._endian + 'i7f'))
                for i in range(nelements):
                    edata = data[n:n+16]  # 8*4
                    n += 16

                    out = s1.unpack(edata)
                    if self.is_debug_file:
                        self.binary_debug.write('OEF_Force_VU-191 - %s\n' % (str(out)))
                    (eid_device, parent, coord, icord) = out

                    eid = eid_device // 10
                    data_in = [eid, parent, coord, icord]

                    forces = []
                    for i in range(nnodes):
                        edata = data[n:n+32]  # 8*4
                        n += 32
                        out = s2.unpack(edata)
                        if self.is_debug_file:
                            self.binary_debug.write('%s\n' % str(out))
                        forces.append(out)
                    data_in.append(forces)

                    #data_in = [vugrid, posit, forceX, shearY, shearZ, torsion, bendY, bendZ]
                    #print "force %s" %(self.get_element_type(self.element_type)), data_in
                    obj.add(nnodes, dt, data_in)
            elif self.format_code == 1 and self.num_wide == 32:  # random
                # TODO: vectorize
                return ndata
            elif self.format_code in [2, 3] and self.num_wide == 32:  # imag
                # TODO: vectorize
                # 32 = 4 + 56/4 * 2 = 4 + 14 * 2 = 4 + 28
                self.create_transient_object(self.force_VU, ComplexForce_VU)
                nnodes = 2
                #ntotal = 16 + 56 * nnodes
                ntotal = self.num_wide * 4
                s1 = Struct(b(self._endian + 'i2i4s'))
                s2 = Struct(b(self._endian + 'i13f'))
                n = 0
                nelements = ndata // ntotal
                obj = self.obj

                for i in range(nelements):
                    edata = data[n:n+16]  # 8*4
                    n += 16

                    out = s1.unpack(edata)
                    if self.is_debug_file:
                        self.binary_debug.write('OEF_Force_191-%s - %s\n' % (self.element_type, str(out)))
                    (eid_device, parent, coord, icord) = out

                    eid = eid_device // 10
                    data_in = [eid, parent, coord, icord]

                    forces = []
                    for i in range(nnodes):
                        edata = data[n:n+56]  # 14*4
                        n += 56
                        out = s2.unpack(edata)
                        if self.is_debug_file:
                            self.binary_debug.write('%s\n' % str(out))
                        [vugrid, posit,
                         force_xr, shear_yr, shear_zr, torsionr, bending_yr, bending_zr,
                         force_xi, shear_yi, shear_zi, torsioni, bending_yi, bending_zi] = out

                        if is_magnitude_phase:
                            force_x = polar_to_real_imag(force_xr, force_xi)
                            shear_y = polar_to_real_imag(shear_yr, shear_yi)
                            shear_z = polar_to_real_imag(shear_zr, shear_zi)
                            torsion = polar_to_real_imag(torsionr, torsioni)
                            bending_y = polar_to_real_imag(bending_yr, bending_yi)
                            bending_z = polar_to_real_imag(bending_zr, bending_zi)
                        else:
                            force_x = complex(force_xr, force_xi)
                            shear_y = complex(shear_yr, shear_yi)
                            shear_z = complex(shear_zr, shear_zi)
                            torsion = complex(torsionr, torsioni)
                            bending_y = complex(bending_yr, bending_yi)
                            bending_z = complex(bending_zr, bending_zi)

                        out2 = [vugrid, posit, force_x, shear_y,
                                shear_z, torsion, bending_z, bending_z]
                        forces.append(out2)
                    data_in.append(forces)

                    #data_in = [vugrid,posit,forceX,shearY,shearZ,torsion,bendY,bendZ]
                    obj.add(nnodes, dt, data_in)
            else:
                msg = self.code_information()
                return self._not_implemented_or_skip(data, ndata, msg)

        elif self.element_type == 228: # CQUADR-NX
            if self.num_wide == 9:
                # real
                return ndata
            elif self.num_wide == 17:
                # complex?
                return ndata
            else:
                raise RuntimeError(self.code_information())

        elif self.element_type in [233, 235]:
            # 233-TRIARLC
            # 235-CQUADR
            if self.read_mode == 1:
                return ndata
            return ndata
        else:
            msg = self.code_information()
            return self._not_implemented_or_skip(data, ndata, msg)

        #assert self.thermal == 0, self.thermal
        assert ndata > 0, ndata
        assert nelements > 0, 'nelements=%r element_type=%s element_name=%r num_wide=%s' % (
            nelements, self.element_type, self.element_name, self.num_wide)
        #assert ndata % ntotal == 0, '%s n=%s nwide=%s len=%s ntotal=%s' % (self.element_name, ndata % ntotal, ndata % self.num_wide, ndata, ntotal)
        assert self.num_wide * 4 == ntotal, 'numwide*4=%s ntotal=%s' % (self.num_wide*4, ntotal)
        assert n > 0, n
        return n

    def _read_oef2_4(self, data):
        """Table 4 parser for OEF2 thermal table"""
        pass
