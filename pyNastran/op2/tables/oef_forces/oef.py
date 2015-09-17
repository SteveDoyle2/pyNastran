#pylint: disable=C0301
from __future__ import print_function
from six.moves import range
from struct import Struct

from pyNastran.op2.op2_helper import polar_to_real_imag
from pyNastran.op2.op2_common import OP2Common

from pyNastran.op2.tables.oef_forces.oef_thermalObjects import (
    # not vectorized
    # TODO: vectorize 7
    HeatFlux_CHBDYx, HeatFlux_2D_3D, HeatFlux_1D, HeatFlux_VU,
    HeatFlux_VUBEAM, HeatFlux_VU_3D, HeatFlux_CONV)
from pyNastran.op2.tables.oef_forces.oef_forceObjects import (
    # vectorized
    RealRodForceArray, RealCRodForce, RealCTubeForce, RealConrodForce,  # rename vectorized classes

    # vectorized (vectorized, unvectorized)
    RealCBarForceArray, RealCBarForce,
    RealPlateForceArray, RealPlateForce,

    # not vectorized
    RealCBar100Force,                                 # TODO: vectorize 1
    RealCBeamForce, RealCShearForce,                  # TODO: vectorize 2
    RealSpringForce, RealDamperForce, RealViscForce,  # TODO: vectorize 3
    RealPlateBilinearForce,                           # TODO: vectorize 1
    RealConeAxForce,                                  # TODO: vectorize 1
    RealCGapForce, RealBendForce,                     # TODO: vectorize 2
    RealPentaPressureForce, RealCBushForce,           # TODO: vectorize 2
    RealForce_VU_2D, RealForce_VU)                    # TODO: vectorize 2
from pyNastran.op2.tables.oef_forces.oef_complexForceObjects import (
    # not vectorized
    # TODO: vectorize 14
    ComplexRodForce, ComplexCBeamForce,
    ComplexCShearForce, ComplexSpringForce,
    ComplexDamperForce, ComplexViscForce,
    ComplexPlateForce, ComplexPlate2Force,
    ComplexBendForce,
    ComplexPentaPressureForce,
    ComplexCBushForce, ComplexForce_VU_2D,
    ComplexCBarForce, ComplexForce_VU)


class OEF(OP2Common):
    """Defines OEFx table reading for element forces/heat flux"""
    def __init__(self):
        OP2Common.__init__(self)

    def OEF_ForceCode(self):
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

    def _read_oef1_3(self, data):
        """Table 3 parser for OEF1 table"""
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
            self.dataNames = self.apply_data_code_value('dataNames', ['loadID'])
            self.setNullNonlinearFactor()
        elif self.analysis_code == 2:  # normal modes/buckling (real eigenvalues)
            #: mode number
            self.mode = self.add_data_parameter(data, 'mode', 'i', 5)
            #: eigenvalue
            self.eigr = self.add_data_parameter(data, 'eigr', 'f', 6, False)
            self.dataNames = self.apply_data_code_value('dataNames', ['mode', 'eigr'])
            # TODO: mode_cycle is not defined?
            #self.dataNames = self.apply_data_code_value('dataNames', ['mode', 'eigr', 'mode_cycle'])
        elif self.analysis_code == 3:  # differential stiffness 0
            #: load set ID number
            self.loadID = self.add_data_parameter(data, 'loadID', 'i', 5)
            self.dataNames = self.apply_data_code_value('dataNames', ['loadID'])
        elif self.analysis_code == 4:  # differential stiffness 1
            #: load set ID number
            self.loadID = self.add_data_parameter(data, 'loadID', 'i', 5)
            self.dataNames = self.apply_data_code_value('dataNames', ['loadID'])
        elif self.analysis_code == 5:   # frequency
            self.freq = self.add_data_parameter(data, 'freq', 'f', 5)  # frequency
            self.dataNames = self.apply_data_code_value('dataNames', ['freq'])
        elif self.analysis_code == 6:  # transient
            self.time = self.add_data_parameter(data, 'time', 'f', 5)  # time step
            self.dataNames = self.apply_data_code_value('dataNames', ['time'])
        elif self.analysis_code == 7:  # pre-buckling
            #: load set ID number
            self.loadID = self.add_data_parameter(data, 'loadID', 'i', 5)
            #self.apply_data_code_value('dataNames',['lsdvmn'])
            self.dataNames = self.apply_data_code_value('dataNames', ['loadID'])
        elif self.analysis_code == 8:  # post-buckling
            #: load set ID number
            self.loadID = self.add_data_parameter(data, 'loadID', 'i', 5)
            #: real eigenvalue
            self.eigr = self.add_data_parameter(data, 'eigr', 'f', 6, False)
            self.dataNames = self.apply_data_code_value('dataNames', ['loadID', 'eigr'])
        elif self.analysis_code == 9:  # complex eigenvalues
            #: mode number
            self.mode = self.add_data_parameter(data, 'mode', 'i', 5)
            #: real eigenvalue
            self.eigr = self.add_data_parameter(data, 'eigr', 'f', 6, False)
            #: imaginary eigenvalue
            self.eigi = self.add_data_parameter(data, 'eigi', 'f', 7, False)
            self.dataNames = self.apply_data_code_value('dataNames', ['mode', 'eigr', 'eigi'])
        elif self.analysis_code == 10:  # nonlinear statics
            #: load step
            self.load_step = self.add_data_parameter(data, 'load_step', 'f', 5)
            self.dataNames = self.apply_data_code_value('dataNames', ['load_step'])
        elif self.analysis_code == 11:  # geometric nonlinear statics
            #: load set ID number
            self.loadID = self.add_data_parameter(data, 'loadID', 'i', 5)
            self.dataNames = self.apply_data_code_value('dataNames', ['loadID'])
        else:
            msg = 'invalid analysis_code...analysis_code=%s' % str(self.analysis_code)
            raise RuntimeError(msg)


        try:
            self.element_name = self.element_mapper[self.element_type]
        except KeyError:
            self.log.error(self.code_information())
            raise

        #self.element_name = self.element_mapper[self.element_type]
        self.data_code['element_name'] = self.element_name

        if self.debug:
            self.binary_debug.write('  element_name  = %r\n' % self.element_name)
            self.binary_debug.write('  approach_code = %r\n' % self.approach_code)
            self.binary_debug.write('  tCode    = %r\n' % self.tCode)
            self.binary_debug.write('  isubcase = %r\n' % self.isubcase)

        self._read_title(data)
        if self.element_type not in self.element_mapper:
            msg = 'element_type = %s' % self.element_type
            return self._not_implemented_or_skip(data, msg)
        self._write_debug_bits()

    def _read_oef2_3(self, data):
        """Table 3 parser for OEF2 table"""
        pass

    def _read_oef1_4(self, data):
        """Table 4 parser for OEF1 table"""
        if self.thermal == 0:
            n = self._read_oef1_loads(data)
        elif self.thermal == 1:
            n = self._read_oef1_thermal(data)
        else:
            n = self._not_implemented_or_skip(data, 'thermal=%s' % self.thermal)
        return n

    def _read_oef1_thermal(self, data):
        """Table 4 parser for OEF1 thermal table"""
        if self.read_mode == 1:
            return len(data)
        if self._results.is_not_saved('element_forces'):
            return len(data)
        n = 0
        #is_magnitude_phase = self.is_magnitude_phase()
        dt = self.nonlinear_factor

        if self.element_type in [1, 2, 3, 10, 34, 69]:  # ROD,BEAM,TUBE,CONROD,BAR,BEND
            # 1-CROD
            # 2-CBEAM
            # 3-CTUBE
            # 10-CONROD
            # 34-CBAR
            # 69-CBEND
            if self.format_code == 1 and self.num_wide == 9:  # real
                self.create_transient_object(self.thermalLoad_1D, HeatFlux_1D)

                ntotal = 36  # 10*4
                s = Struct(self._endian + b'i8s6f')
                nelements = len(data) // ntotal
                for i in range(nelements):
                    edata = data[n:n+ntotal]

                    out = s.unpack(edata)
                    (eid_device, eType, xgrad, ygrad, zgrad, xflux, yflux, zflux) = out
                    eid = (eid_device - self.device_code) // 10

                    data_in = [eid, eType, xgrad, ygrad, zgrad, xflux, yflux, zflux]
                    #print "heatFlux %s" % (self.get_element_type(self.element_type)), data_in
                    #eid = self.obj.add_new_eid(out)
                    self.obj.add(dt, data_in)
                    n += ntotal
            else:
                msg = self.code_information()
                return self._not_implemented_or_skip(data, msg)

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
            if self.format_code == 1 and self.num_wide == 9:  # real - 2D
                # [33, 53, 64, 74, 75]
                self.create_transient_object(self.thermalLoad_2D_3D, HeatFlux_2D_3D)
                # no zed on this element for some reason...
                ntotal = 36
                nelements = len(data) // ntotal

                s = Struct(self._endian + b'i8s6f')
                for i in range(nelements):
                    edata = data[n:n+ntotal]
                    n += ntotal
                    out = s.unpack(edata)
                    (eid_device, etype, xgrad, ygrad, zgrad, xflux, yflux, zflux) = out
                    eid = (eid_device - self.device_code) // 10
                    data_in = [eid, etype, xgrad, ygrad, zgrad, xflux, yflux, zflux]
                    self.obj.add(dt, data_in)

            elif self.format_code == 1 and self.num_wide == 10:  # real - 3D
                # [39, 67, 68]:  # HEXA,PENTA
                self.create_transient_object(self.thermalLoad_2D_3D, HeatFlux_2D_3D)
                ntotal = 40
                nelements = len(data) // ntotal

                s = Struct(self._endian + b'i8s6fi')
                for i in range(nelements):
                    edata = data[n:n+ntotal]
                    n += ntotal
                    out = s.unpack(edata)
                    (eid_device, etype, xgrad, ygrad, zgrad, xflux, yflux, zflux, zed) = out
                    eid = (eid_device - self.device_code) // 10
                    data_in = [eid, etype, xgrad, ygrad, zgrad, xflux, yflux, zflux]
                    self.obj.add(dt, data_in)
            else:
                msg = self.code_information()
                return self._not_implemented_or_skip(data, msg)
        elif self.element_type in [107, 108, 109]:  # CHBDYE, CHBDYG, CHBDYP
            # 107-CHBDYE
            # 108-CHBDYG
            # 109-CHBDYP
            if self.format_code == 1 and self.num_wide == 8:  # real
                self.create_transient_object(self.thermalLoad_CHBDY, HeatFlux_CHBDYx)

                s1 = Struct(self._endian + b'i8s5f')
                ntotal = 32
                nelements = len(data) // ntotal
                for i in range(nelements):
                    edata = data[n:n+32]
                    n += ntotal
                    out = s1.unpack(edata)
                    (eid_device, etype, fapplied, free_conv, force_conv, frad, ftotal) = out
                    eid = (eid_device - self.device_code) // 10

                    data_in = [eid, etype, fapplied, free_conv, force_conv, frad, ftotal]
                    if self.debug4():
                        self.binary_debug.write('  %s -> [%s, %s, %s, %s, %s, %s, %s]'
                                                % (eid, eid_device, etype, fapplied, free_conv, force_conv, frad, ftotal))
                    #print "heatFlux %s" % (self.get_element_type(self.element_type)), data_in
                    self.obj.add(dt, data_in)
            else:
                msg = self.code_information()
                return self._not_implemented_or_skip(data, msg)

        elif self.element_type in [110]:
            # 110-CONV
            if self.format_code == 1 and self.num_wide == 4:
                self.create_transient_object(self.thermalLoad_CONV, HeatFlux_CONV)
                s1 = Struct(self._endian + b'ifif')
                ntotal = 16
                nelements = len(data) // ntotal
                for i in range(nelements):
                    edata = data[n:n+16]
                    n += 16
                    out = s1.unpack(edata)
                    (eid_device, cntl_node, free_conv, free_conv_k) = out
                    eid = (eid_device - self.device_code) // 10
                    data_in = [eid, cntl_node, free_conv, free_conv_k]
                    #print "heatFlux %s" %(self.get_element_type(self.element_type)), data_in
                    self.obj.add(dt, data_in)
            else:
                msg = self.code_information()
                return self._not_implemented_or_skip(data, msg)

        elif self.element_type in [145, 146, 147]:  # VUHEXA,VUPENTA,VUTETRA
            # 145-VUHEXA
            # 146-VUPENTA
            # 147-VUTETRA
            self.create_transient_object(self.thermalLoad_VU_3D, HeatFlux_VU_3D)
            if self.element_type == 147:  # VUTETRA
                nnodes = 4
            elif self.element_type == 146:  # VUPENTA
                nnodes = 6
            elif self.element_type == 145:  # VUHEXA
                nnodes = 8
            else:
                msg = self.code_information()
                return self._not_implemented_or_skip(data, msg)

            numwide_real = 2 + 7 * nnodes
            if self.format_code == 1 and self.num_wide == numwide_real:  # real
                ntotal = 8 + 28 * nnodes
                s1 = Struct(self._endian + b'ii')
                s2 = Struct(self._endian + b'i6f')
                nelements = len(data) // ntotal
                for i in range(nelements):
                    edata = data[n:n+8]  # 2*4
                    n += 8
                    out = s1.unpack(edata)
                    (eid_device, parent) = out
                    eid = (eid_device - self.device_code) // 10
                    data_in = [eid, parent]
                    grad_fluxes = []
                    for i in range(nnodes):
                        edata = data[0:28]
                        n += 28
                        #print "i=%s len(edata)=%s" %(i,len(edata))
                        out = s2.unpack(edata)
                        grad_fluxes.append(out)
                    data_in.append(grad_fluxes)
                    #print "heatFlux %s" %(self.get_element_type(self.element_type)), data_in
                    self.obj.add(nnodes, dt, data_in)
            else:
                msg = self.code_information()
                return self._not_implemented_or_skip(data, msg)

        elif self.element_type in [189, 190]:  # VUQUAD,VUTRIA
            # 189-VUQUAD
            # 190-VUTRIA
            if self.format_code == 1 and self.num_wide == 27:  # real
                self.create_transient_object(self.thermalLoad_VU, HeatFlux_VU)
                if self.element_type == 189:
                    nnodes = 4
                elif self.element_type == 190:
                    nnodes = 3
                elif self.element_type == 191:
                    nnodes = 2
                else:
                    raise NotImplementedError(self.code_information())

                ntotal = 24 + 28 * nnodes
                s1 = Struct(self._endian + b'iii4sii')
                s2 = Struct(self._endian + b'i6f')
                nelements = len(data) // ntotal
                for i in range(nelements):
                    edata = data[n:n+24]  # 6*4
                    n += 24

                    out = s1.unpack(edata)
                    (eid_device, parent, coord, icord, theta, null) = out
                    eid = (eid_device - self.device_code) // 10
                    data_in = [eid, parent, coord, icord, theta]

                    grad_fluxes = []
                    for i in range(nnodes):
                        edata = data[n:n+28]  # 7*4
                        n += 28
                        out = s2.unpack(edata)
                        grad_fluxes.append(out)
                    data_in.append(grad_fluxes)
                    #data_in = [eid, eType, xGrad, yGrad, zGrad, xFlux, yFlux, zFlux]
                    #print "heatFlux %s" %(self.get_element_type(self.element_type)), data_in
                    self.obj.add(nnodes, dt, data_in)
            else:
                msg = self.code_information()
                return self._not_implemented_or_skip(data, msg)

        elif self.element_type in [191]:  # VUBEAM
            #assert self.num_wide==27,self.code_information()
            self.create_transient_object(self.thermalLoad_VUBeam, HeatFlux_VUBEAM)
            nnodes = 2

            numwide_real = 4 + 7 * nnodes
            if self.format_code == 1 and self.num_wide == numwide_real:  # real
                ntotal = 16 + 28 * nnodes
                s1 = Struct(self._endian + b'iii4s')
                s2 = Struct(self._endian + b'i6f')
                nelements = len(data) // ntotal
                for i in range(nelements):
                    edata = data[n:n+16]  # 4*4
                    n += 16

                    out = s1.unpack(edata)
                    (eid_device, parent, coord, icord) = out

                    eid = (eid_device - self.device_code) // 10
                    data_in = [eid, parent, coord, icord]

                    grad_fluxes = []
                    for i in range(nnodes):
                        edata = data[n:n+28]  # 7*4
                        n += 28
                        out = s2.unpack(edata)
                        grad_fluxes.append(out)
                    data_in.append(grad_fluxes)

                    #data_in = [eid, eType, xGrad, yGrad, zGrad, xFlux, yFlux, zFlux]
                    #print "heatFlux %s" %(self.get_element_type(self.element_type)), data_in
                    self.obj.add(nnodes, dt, data_in)
            else:
                msg = 'OEF sort1 thermal Type=%s num=%s' % (self.element_name, self.element_type)
                return self._not_implemented_or_skip(data, msg)
        else:
            msg = 'OEF sort1 thermal Type=%s num=%s' % (self.element_name, self.element_type)
            return self._not_implemented_or_skip(data, msg)


        assert self.thermal == 1, self.thermal
        assert len(data) > 0, len(data)
        assert nelements > 0, 'nelements=%r element_type=%s element_name=%r' % (nelements, self.element_type, self.element_name)
        #assert len(data) % ntotal == 0, '%s n=%s nwide=%s len=%s ntotal=%s' % (self.element_name, len(data) % ntotal, len(data) % self.num_wide, len(data), ntotal)
        assert self.num_wide * 4 == ntotal, 'numwide*4=%s ntotal=%s' % (self.num_wide*4, ntotal)
        assert n > 0, n
        return n

    def print_obj_name_on_crash(func):
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

    @print_obj_name_on_crash
    def _read_oef1_loads(self, data):
        """
        Reads the OEF1 table; stores the element forces/heat flux.
        """
        if self._results.is_not_saved('element_forces'):
            return len(data)
        (num_wide_real, num_wide_imag) = self.OEF_ForceCode()
        if self.debug4():
            self.binary_debug.write('  num_wide_real = %r\n' % num_wide_real)
            self.binary_debug.write('  num_wide_imag = %r\n' % num_wide_imag)

        n = 0
        is_magnitude_phase = self.is_magnitude_phase()
        dt = self.nonlinear_factor

        if self.element_type in [1, 3, 10]:  # rods
            #1-CROD
            #3-CTUBE
            #10-CONROD
            #obj_complex = ComplexRodForce

            obj_vector_real = RealRodForceArray
            obj_vector_complex = None  #ComplexRodForceArray
            if self.element_type == 1: # CROD
                result_vector_name = 'crod_force'
                slot_vector = self.crod_force
                obj_real = RealCRodForce
            elif self.element_type == 3:  # CTUBE
                result_vector_name = 'ctube_force'
                slot_vector = self.ctube_force
                obj_real = RealCTubeForce
            elif self.element_type == 10:  # CONROD
                result_vector_name = 'conrod_force'
                slot_vector = self.conrod_force
                obj_real = RealConrodForce
            else:
                msg = 'sort1 Type=%s num=%s' % (self.element_name, self.element_type)
                return self._not_implemented_or_skip(data, msg)
            result_name = result_vector_name
            slot = slot_vector

            if self.format_code == 1 and self.num_wide == 3: # real
                ntotal = 12 # 3 * 4
                nelements = len(data) // ntotal
                auto_return = self._create_oes_object2(nelements,
                                                       result_name, result_vector_name,
                                                       slot, slot_vector,
                                                       obj_real, obj_vector_real)
                if auto_return:
                    return nelements * self.num_wide * 4

                s = Struct(self._endian + b'iff')  # 3
                for i in range(nelements):
                    edata = data[n:n+ntotal]
                    out = s.unpack(edata)
                    (eid_device, axial, torque) = out
                    eid = (eid_device - self.device_code) // 10
                    if self.debug4():
                        self.binary_debug.write('OEF_Rod - %s\n' % (str(out)))

                    assert eid > 0, eid
                    self.obj.add(dt, eid, axial, torque)
                    n += ntotal

            elif self.format_code in [2, 3] and self.num_wide == 5: # imag
                self.create_transient_object(slot, ComplexRodForce)

                s = Struct(self._endian + b'i4f')
                ntotal = 20 # 5*4
                nelements = len(data) // ntotal
                for i in range(nelements):
                    edata = data[n:n+20]
                    out = s.unpack(edata)
                    if self.debug4():
                        self.binary_debug.write('OEF_Rod - %s\n' % (str(out)))
                    (eid_device, axial_real, torque_real, axial_imag, torque_imag) = out

                    if is_magnitude_phase:
                        axial = polar_to_real_imag(axial_real, axial_imag)
                        torque = polar_to_real_imag(torque_real, torque_imag)
                    else:
                        axial = complex(axial_real, axial_imag)
                        torque = complex(torque_real, torque_imag)
                    eid = (eid_device - self.device_code) // 10
                    assert eid > 0, eid

                    data_in = [eid, axial, torque]
                    #print "%s" % (self.get_element_type(self.element_type)), data_in
                    self.obj.add(dt, data_in)
                    n += ntotal

            else:
                msg = self.code_information()
                return self._not_implemented_or_skip(data, msg)

        elif self.element_type == 2:  # cbeam
            #2-CBEAM
            if self.read_mode == 1:
                return len(data)

            result_name = 'cbeam_force'
            if self._results.is_not_saved(result_name):
                return len(data)
            self._results._found_result(result_name)
            if self.format_code == 1 and self.num_wide == 9:  # real centroid ???
                self.create_transient_object(self.cbeam_force, RealCBeamForce)
                s = Struct(self._endian + b'i8f')  # 36
                ntotal = 36
                nelements = len(data) // ntotal
                for i in range(nelements):
                    edata = data[n:n+36]
                    out = s.unpack(edata)
                    if self.debug4():
                        self.binary_debug.write('OEF_Beam - %s\n' % (str(out)))
                    (eid_device, sd, bm1, bm2, ts1, ts2, af, ttrq, wtrq) = out
                    eid = (eid_device - self.device_code) // 10
                    assert eid > 0, eid
                    n += 36

            elif self.format_code == 1 and self.num_wide == 100:  # real
                self.create_transient_object(self.cbeam_force, RealCBeamForce)
                s1 = Struct(self._endian + b'i')
                s2 = Struct(self._endian + b'i8f')  # 36

                ntotal = 400  # 1+(10-1)*11=100 ->100*4 = 400
                nelements = len(data) // ntotal
                for i in range(nelements):
                    edata = data[n:n+4]
                    eid_device, = s1.unpack(edata)
                    eid = (eid_device - self.device_code) // 10
                    n += 4

                    for i in range(11):
                        edata = data[n:n+36]
                        out = s2.unpack(edata)
                        if self.debug4():
                            self.binary_debug.write('OEF_Beam - %s\n' % (str(out)))
                        (nid, sd, bm1, bm2, ts1, ts2, af, ttrq, wtrq) = out

                        data_in = [eid, nid, sd, bm1, bm2, ts1, ts2, af, ttrq, wtrq]
                        assert eid > 0, eid
                        if i == 0:  # isNewElement
                            self.obj.add_new_element(dt, data_in)
                        elif sd > 0.:
                            self.obj.add(dt, data_in)
                        n += 36
            elif self.format_code in [2, 3] and self.num_wide == 177: # imag
                self.create_transient_object(self.cbeam_force, ComplexCBeamForce)
                s1 = Struct(self._endian + b'i')
                s2 = Struct(self._endian + b'i15f')
                ntotal = 708  # (16*11+1)*4 = 177*4
                nelements = len(data) // ntotal
                for i in range(nelements):
                    edata = data[n:n+4]
                    eid_device, = s1.unpack(edata)
                    eid = (eid_device - self.device_code) // 10
                    assert eid > 0, eid

                    n += 4
                    for i in range(11):
                        edata = data[n:n+64]
                        n += 64

                        out = s2.unpack(edata)
                        if self.debug4():
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
                        #eid = self.obj.add_new_eid(out)
                        if i == 0:  # isNewElement:
                            data_in = [eid, nid, sd, bm1, bm2,
                                       ts1, ts2, af, ttrq, wtrq]
                            #print "%s cNew   " % (self.get_element_type(self.element_type)), data_in
                            self.obj.add_new_element(dt, data_in)
                        elif sd > 0.:
                            data_in = [eid, nid, sd, bm1, bm2,
                                       ts1, ts2, af, ttrq, wtrq]
                            #print "%s cOld   " % (self.get_element_type(self.element_type)), data_in
                            self.obj.add(dt, data_in)
            else:
                msg = self.code_information()
                return self._not_implemented_or_skip(data, msg)
            #print self.beamForces

        elif self.element_type in [11, 12, 13, 14,   # springs
                                   20, 21, 22, 23]:  # dampers
            if self.read_mode == 1:
                return len(data)
            # 11-CELAS1
            # 12-CELAS2
            # 13-CELAS3
            # 14-CELAS4

            # 20-CDAMP1
            # 21-CDAMP2
            # 22-CDAMP3
            # 23-CDAMP4
            if self.element_type == 11:
                slot = self.celas1_force
                real_class_obj = RealSpringForce
                complex_class_obj = ComplexSpringForce
            elif self.element_type == 12:
                slot = self.celas2_force
                real_class_obj = RealSpringForce
                complex_class_obj = ComplexSpringForce
            elif self.element_type == 13:
                slot = self.celas3_force
                real_class_obj = RealSpringForce
                complex_class_obj = ComplexSpringForce
            elif self.element_type == 14:
                slot = self.celas4_force
                real_class_obj = RealSpringForce
                complex_class_obj = ComplexSpringForce

            elif self.element_type == 20:
                slot = self.cdamp1_force
                real_class_obj = RealDamperForce
                complex_class_obj = ComplexDamperForce
            elif self.element_type == 21:
                slot = self.cdamp2_force
                real_class_obj = RealDamperForce
                complex_class_obj = ComplexDamperForce
            elif self.element_type == 22:
                slot = self.cdamp3_force
                real_class_obj = RealDamperForce
                complex_class_obj = ComplexDamperForce
            elif self.element_type == 23:
                slot = self.cdamp4_force
                real_class_obj = RealDamperForce
                complex_class_obj = ComplexDamperForce
            else:
                msg = self.element_type
                return self._not_implemented_or_skip(data, msg)

            if self.format_code == 1 and self.num_wide == 2:  # real
                self.create_transient_object(slot, real_class_obj)

                s = Struct(self._endian + b'if')  # 2
                ntotal = 8  # 2*4
                nelements = len(data) // ntotal
                for i in range(nelements):
                    edata = data[n:n + 8]

                    out = s.unpack(edata)
                    if self.debug4():
                        self.binary_debug.write('OEF_SpringDamper - %s\n' % str(out))
                    (eid_device, force) = out
                    eid = (eid_device - self.device_code) // 10
                    assert eid > 0, eid

                    data_in = [eid, force]
                    #print "%s" % (self.get_element_type(self.element_type)), data_in
                    self.obj.add(dt, data_in)
                    n += ntotal
            elif self.format_code in [2, 3] and self.num_wide == 3:  # imag
                self.create_transient_object(slot, complex_class_obj)

                s = Struct(self._endian + b'i2f')
                ntotal = 12  # 3*4

                nelements = len(data) // ntotal
                for i in range(nelements):
                    edata = data[n:n + 12]
                    out = s.unpack(edata)
                    if self.debug4():
                        self.binary_debug.write('OEF_SpringDamper - %s\n' % str(out))
                    (eid_device, force_real, force_imag) = out
                    eid = (eid_device - self.device_code) // 10
                    assert eid > 0, eid

                    if is_magnitude_phase:
                        force = polar_to_real_imag(force_real, force_imag)
                    else:
                        force = complex(force_real, force_imag)

                    data_in = [eid, force]
                    #print "%s" % (self.get_element_type(self.element_type)), data_in
                    #eid = self.obj.add_new_eid(out)
                    self.obj.add(dt, data_in)
                    n += ntotal
            else:
                #msg = 'OEF: element_name=%s element_type=%s' % (self.element_name, self.element_type)
                msg = self.code_information()
                return self._not_implemented_or_skip(data, msg)

        elif self.element_type == 24:  # CVISC
            if self.read_mode == 1:
                return len(data)
            if self.format_code == 1 and self.num_wide == 3: # real
                self.create_transient_object(self.cvisc_force, RealViscForce)
                s = Struct(self._endian + b'iff')
                ntotal = 12  # 3*4
                nelements = len(data) // 12
                for i in range(nelements):
                    edata = data[n:n+12]

                    out = s.unpack(edata)
                    if self.debug4():
                        self.binary_debug.write('OEF_CVisc - %s\n' % (str(out)))
                    (eid_device, axial, torque) = out
                    eid = (eid_device - self.device_code) // 10
                    assert eid > 0, eid

                    data_in = [eid, axial, torque]
                    #print "%s" % (self.get_element_type(self.element_type)), data_in
                    #eid = self.obj.add_new_eid(out)
                    self.obj.add(dt, data_in)
                    n += ntotal
            elif self.format_code in [2, 3] and self.num_wide == 5: # complex
                self.create_transient_object(self.cvisc_force, ComplexViscForce)
                s = Struct(self._endian + b'i4f')  # 5
                ntotal = 20  # 5*4
                nelements = len(data) // ntotal
                for i in range(nelements):
                    edata = data[n:n+20]

                    out = s.unpack(edata)
                    if self.debug4():
                        self.binary_debug.write('OEF_CVisc - %s\n' % (str(out)))
                    (eid_device, axial_real, torque_real, axial_imag, torque_imag) = out
                    eid = (eid_device - self.device_code) // 10
                    assert eid > 0, eid

                    if is_magnitude_phase:
                        axial = polar_to_real_imag(axial_real, axial_imag)
                        torque = polar_to_real_imag(torque_real, torque_imag)
                    else:
                        axial = complex(axial_real, axial_imag)
                        torque = complex(torque_real, torque_imag)

                    data_in = [eid, axial, torque]
                    #print "%s" % (self.get_element_type(self.element_type)), data_in
                    #eid = self.obj.add_new_eid(out)
                    self.obj.add(dt, data_in)
                    n += ntotal
            else:
                msg = self.code_information()
                return self._not_implemented_or_skip(data, msg)
            #print self.viscForces

        elif self.element_type == 34:  # cbar
            # 34-CBAR
            slot = self.cbar_force
            slot_vector = self.cbar_force
            result_name = 'cbar_force'
            result_vector_name = result_name
            obj_real = RealCBarForce
            obj_vector_real = RealCBarForceArray

            if self.format_code == 1 and self.num_wide == 9: # real
                ntotal = 36  # 9*4
                nelements = len(data) // ntotal
                auto_return = self._create_oes_object2(nelements,
                                                       result_name, result_vector_name,
                                                       slot, slot_vector,
                                                       obj_real, obj_vector_real)
                if auto_return:
                    return nelements * self.num_wide * 4
                s = Struct(self._endian + b'i8f')  # 9
                for i in range(nelements):
                    edata = data[n:n+36]

                    out = s.unpack(edata)
                    if self.debug4():
                        self.binary_debug.write('OEF_CBar - %s\n' % (str(out)))
                    (eid_device, bm1a, bm2a, bm1b, bm2b, ts1, ts2, af, trq) = out
                    eid = (eid_device - self.device_code) // 10
                    assert eid > 0, eid

                    data_in = [eid, bm1a, bm2a, bm1b, bm2b, ts1, ts2, af, trq]
                    #print "%s" % (self.get_element_type(self.element_type)), data_in
                    #eid = self.obj.add_new_eid(out)
                    self.obj.add(dt, data_in)
                    n += ntotal
            elif self.format_code in [2, 3] and self.num_wide == 17: # imag
                if self.read_mode == 1:
                    return len(data)
                self.create_transient_object(self.cbar_force, ComplexCBarForce)
                s = Struct(self._endian + b'i16f')
                ntotal = 68  # 17*4
                nelements = len(data) // ntotal
                for i in range(nelements):
                    edata = data[n:n+68]

                    out = s.unpack(edata)
                    (eid_device,
                     bm1ar, bm2ar, bm1br, bm2br, ts1r, ts2r, afr, trqr,
                     bm1ai, bm2ai, bm1bi, bm2bi, ts1i, ts2i, afi, trqi) = out
                    if self.debug4():
                        self.binary_debug.write('OEF_CBar - %s\n' % (str(out)))
                    eid = (eid_device - self.device_code) // 10
                    assert eid > 0, eid

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

                    data_in = [eid, bm1a, bm2a, bm1b, bm2b, ts1, ts2, af, trq]
                    #print "%s" % (self.get_element_type(self.element_type)), data_in
                    #eid = self.obj.add_new_eid(out)
                    self.obj.add(dt, data_in)
                    n += ntotal
            else:
                msg = self.code_information()
                return self._not_implemented_or_skip(data, msg)
            #print self.barForces

        elif self.element_type == 100:  # cbar
            #100-BARS
            if self.read_mode == 1:
                return len(data)
            if self.format_code == 1 and self.num_wide == 8:  # real
                self.create_transient_object(self.cbar100_force, RealCBar100Force)

                s = Struct(self._endian + b'i7f')
                ntotal = 32  # 8*4
                nelements = len(data) // ntotal
                for i in range(nelements):
                    edata = data[n:n+32]

                    out = s.unpack(edata)
                    if self.debug4():
                        self.binary_debug.write('OEF_CBar100 - %s\n' % (str(out)))
                    (eid_device, sd, bm1, bm2, ts1, ts2, af, trq) = out
                    eid = (eid_device - self.device_code) // 10
                    assert eid > 0, eid

                    data_in = [eid, sd, bm1, bm2, ts1, ts2, af, trq]
                    #print "%s" %(self.get_element_type(self.element_type)), data_in
                    #eid = self.obj.add_new_eid(out)
                    self.obj.add(dt, data_in)
                    n += 32
                #elif self.format_code in [2, 3] and self.num_wide == 14:  # imag
            else:
                msg = self.code_information()
                return self._not_implemented_or_skip(data, msg)

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
                return self._not_implemented_or_skip(data, msg)
            slot_vector = slot
            result_vector_name = result_name

            obj_real = RealPlateForce
            obj_complex = ComplexPlateForce

            ComplexPlateForceArray = None
            obj_vector_real = RealPlateForceArray
            obj_vector_complex = ComplexPlateForceArray

            if self.format_code == 1 and self.num_wide == 9:  # real
                ntotal = 36 # 9*4
                nelements = len(data) // ntotal
                auto_return = self._create_oes_object2(nelements,
                                                       result_name, result_vector_name,
                                                       slot, slot_vector,
                                                       obj_real, obj_vector_real)
                if auto_return:
                    return nelements * self.num_wide * 4

                s = Struct(self._endian + b'i8f')
                for i in range(nelements):
                    edata = data[n:n+36]

                    out = s.unpack(edata)
                    if self.debug4():
                        self.binary_debug.write('real_OEF_Plate-%s - %s\n' % (self.element_type, str(out)))
                    (eid_device, mx, my, mxy, bmx, bmy, bmxy, tx, ty) = out
                    eid = (eid_device - self.device_code) // 10
                    assert eid > 0, 'eid_device=%s eid=%s table_name-%r' % (eid_device, eid, self.table_name)

                    #print "%s" % (self.get_element_type(self.element_type)), data_in
                    #eid = self.obj.add_new_eid(out)
                    self.obj.add(dt, eid, mx, my, mxy, bmx, bmy, bmxy, tx, ty)
                    n += ntotal
            elif self.format_code in [2, 3] and self.num_wide == 17:  # imag
                if self.read_mode == 1:
                    return len(data)
                self.create_transient_object(slot, ComplexPlateForce)
                s = Struct(self._endian + b'i16f')

                ntotal = 68
                nelements = len(data) // ntotal
                for i in range(nelements):
                    edata = data[n:n+68]
                    out = s.unpack(edata)
                    (eid_device,
                     mxr, myr, mxyr, bmxr, bmyr, bmxyr, txr, tyr,
                     mxi, myi, mxyi, bmxi, bmyi, bmxyi, txi, tyi) = out
                    eid = (eid_device - self.device_code) // 10
                    assert eid > 0, eid
                    if self.debug4():
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
                    assert eid > 0, eid
                    #print "%s" % (self.get_element_type(self.element_type)), data_in
                    self.obj.add(dt, eid, mx, my, mxy, bmx, bmy, bmxy, tx, ty)
                    n += ntotal
            else:
                msg = self.code_information()
                return self._not_implemented_or_skip(data, msg)

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
                slot = self.cquad8_force
            elif self.element_type == 144:
                result_name = 'cquad4_force'
                slot = self.cquad4_force
            else:
                msg = self.code_information()
                return self._not_implemented_or_skip(data, msg)
            slot_vector = slot
            result_vector_name = result_name

            if self.element_type in [70, 75]:  # CTRIAR,CTRIA6
                nnodes = 3
            elif self.element_type in [64, 82, 144]:  # CQUAD8,CQUADR,CQUAD4-bilinear
                nnodes = 4
            else:
                msg = 'name=%r type=%r' % (self.element_name, self.element_type)
                return self._not_implemented_or_skip(data, msg)

            numwide_real = 2 + (nnodes + 1) * 9 # centroidal node is the + 1
            numwide_imag = 2 + (nnodes + 1) * 17

            if 0:
                obj_real = RealPlateBilinearForce
                obj_complex = ComplexPlate2Force

                RealPlateForce2Array = None
                ComplexPlateForce2Array = None
                obj_vector_real = RealPlateForce2Array
                obj_vector_complex = ComplexPlateForce2Array

            if self._results.is_not_saved(result_name):
                return len(data)
            self._results._found_result(result_name)
            if self.format_code == 1 and self.num_wide == numwide_real:  # real
                if self.read_mode == 1:
                    return len(data)

                self.create_transient_object(slot, RealPlateBilinearForce)
                ntotal = 8 + (nnodes + 1) * 36 # centroidal node is the + 1
                assert ntotal == self.num_wide * 4, 'ntotal=%s numwide=%s' % (ntotal, self.num_wide * 4)
                nelements = len(data) // ntotal
                #auto_return = self._create_oes_object2(nelements,
                                                       #result_name, result_vector_name,
                                                       #slot, slot_vector,
                                                       #obj_real, obj_vector_real)
                #if auto_return:
                    #return nelements * self.num_wide * 4

                s1 = Struct(self._endian + b'i4si8f')  # 8+36
                s2 = Struct(self._endian + b'i8f') # 36
                nelements = len(data) // ntotal

                for i in range(nelements):
                    edata = data[n:n+44]

                    out = s1.unpack(edata)
                    if self.debug4():
                        self.binary_debug.write('OEF_Plate2-%s - %s\n' % (self.element_type, str(out)))
                    (eid_device, term, nid, mx, my, mxy, bmx, bmy, bmxy, tx, ty) = out
                    #term= 'CEN\'

                    eid = (eid_device - self.device_code) // 10
                    assert eid > 0, eid

                    #print "%s" % (self.get_element_type(self.element_type)), dt, term, nid, mx, my, mxy, bmx, bmy, bmxy, tx, ty
                    self.obj.add_new_element(eid, dt, term, nid, mx, my, mxy, bmx, bmy, bmxy, tx, ty)
                    n += 44
                    for i in range(nnodes):
                        edata = data[n : n + 36]
                        out = s2.unpack(edata)
                        if self.debug4():
                            self.binary_debug.write('    %s\n' % (str(out)))
                        (nid, mx, my, mxy, bmx, bmy, bmxy, tx, ty) = out
                        assert nid > 0, 'nid=%s' % nid
                        data_in = [nid, mx, my, mxy, bmx, bmy, bmxy, tx, ty]
                        #print "***%s    " % (self.get_element_type(self.element_type)), data_in
                        self.obj.add(eid, dt, nid, mx, my, mxy, bmx, bmy, bmxy, tx, ty)
                        n += 36
            elif self.format_code in [2, 3] and self.num_wide == numwide_imag: # complex
                if self.read_mode == 1:
                    return len(data)

                self.create_transient_object(slot, ComplexPlate2Force)
                s1 = Struct(self._endian + b'i4s17f')  # 2+17=19 * 4 = 76
                s2 = Struct(self._endian + b'i16f')  # 17 * 4 = 68
                ntotal = 8 + (nnodes + 1) * 68
                nelements = len(data) // ntotal
                for i in range(nelements):
                    edata = data[n:n + 76]
                    n += 76

                    out = s1.unpack(edata)
                    if self.debug4():
                        self.binary_debug.write('OEF_Plate2-%s - %s\n' % (self.element_type, str(out)))
                    (eid_device, term, nid,
                     mxr, myr, mxyr, bmxr, bmyr, bmxyr, txr, tyr,
                     mxi, myi, mxyi, bmxi, bmyi, bmxyi, txi, tyi) = out
                    #term = 'CEN\'

                    eid = (eid_device - self.device_code) // 10
                    assert eid > 0, eid
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

                    data_in = [term, nid, mx, my, mxy, bmx, bmy, bmxy, tx, ty]
                    #print "%s" % (self.get_element_type(self.element_type)), data_in
                    self.obj.add_new_element(eid, dt, data_in)

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
                        if self.debug4():
                            self.binary_debug.write('OEF_Plate2 - eid=%i nid=%s out=%s\n' % (eid, nid, str(out)))
                        data_in = [nid, mx, my, mxy, bmx, bmy, bmxy, tx, ty]
                        assert eid > 0, eid
                        #print "***%s    " % (self.get_element_type(self.element_type)),data_in
                        self.obj.add(eid, dt, data_in)
            else:
                msg = self.code_information()
                return self._not_implemented_or_skip(data, msg)

        elif self.element_type in [95, 96, 97, 98]: # composites
            # 95 - CQUAD4
            # 96 - CQUAD8
            # 97 - CTRIA3
            # 98 - CTRIA6 (composite)
            if self.read_mode == 1:
                return len(data)
            if self.format_code == 1 and self.num_wide == 9:  # real
                return len(data)
                #print self.code_information()
                #self.create_transient_object(self.compositePlateForces, RealCompositePlateForce)  # undefined
                ##return
                #ntotal = 9 * 4
                #nelements = len(data) // ntotal
                #if self.debug:
                    #self.binary_debug.write('  [cap, element1, element2, ..., cap]\n')
                    #self.binary_debug.write('  cap = %i  # assume 1 cap when there could have been multiple\n' % len(data))
                    ##self.binary_debug.write('  #centeri = [eid_device, j, grid, fd1, sx1, sy1, txy1, angle1, major1, minor1, vm1,\n')
                    ##self.binary_debug.write('  #                                fd2, sx2, sy2, txy2, angle2, major2, minor2, vm2,)]\n')
                    ##self.binary_debug.write('  #nodeji = [eid, iLayer, o1, o2, t12, t1z, t2z, angle, major, minor, ovm)]\n')
                    #self.binary_debug.write('  nelements=%i; nnodes=1 # centroid\n' % nelements)

                #eid_old = 0
                #s = Struct(self._endian + b'i8si4f4s')
                #for i in range(nelements):
                    #if i % 10000 == 0:
                        #print 'i = ', i
                    #edata = data[n:n+ntotal]  # 4*9
                    #out = s.unpack(edata)
                    #(eid_device, theory, lamid, failure_index_direct_stress, failure_mode_max_shear,
                             #failure_index_interlaminar_shear, fmax, failure_flag) = out
                    #eid = (eid_device - self.device_code) // 10
                    #if self.debug4():
                        #if eid > 0:
                            #self.binary_debug.write('  eid=%i; C=[%s]\n' % (', '.join(['%r' % di for di in out]) ))
                        #else:
                            #self.binary_debug.write('      %s  C=[%s]\n' % (' ' * len(str(eid)), ', '.join(['%r' % di for di in out]) ))

                    #if eid > 0:
                        #self.obj.add_new_eid(eType, dt, eid, o1, o2, t12, t1z, t2z, angle, major, minor, ovm)
                    #else:
                        #self.obj.add(dt, eid, o1, o2, t12, t1z, t2z, angle, major, minor, ovm)
                    #eid_old = eid
                    #n += ntotal
            else:
                msg = self.code_information()
                return self._not_implemented_or_skip(data, msg)

        elif self.element_type in [39, 67, 68]: # solids
            # 39-CTETRA
            # 67-CHEXA
            # 68-CPENTA
            if self.read_mode == 1:
                return len(data)
            self._results._found_result('solidForces')
            if self.format_code == 1 and self.num_wide == 0:  # real
                #self.create_transient_object(self.solidForces, RealCSolidForce)
                pass
            else:
                msg = self.code_information()
                return self._not_implemented_or_skip(data, msg)

        elif self.element_type in [53]:  # ctriax6
            # 53-CTRIAX6
            if self.read_mode == 1:
                return len(data)
            self._results._found_result('ctriax_force')
            if self.format_code == 1 and self.num_wide == 0:  # real
                pass
                #self.create_transient_object(self.ctriax_force, RealCTriaxForce)  # undefined
            else:
                msg = self.code_information()
                return self._not_implemented_or_skip(data, msg)

        elif self.element_type == 4:  # cshear
            # 4-CSHEAR
            if self.read_mode == 1:
                return len(data)
            self._results._found_result('cshear_force')
            if self.format_code == 1 and self.num_wide == 17:  # real
                self.create_transient_object(self.cshear_force, RealCShearForce)
                s = Struct(self._endian + b'i16f')
                ntotal = 68  # 17*4
                nelements = len(data) // 68
                for i in range(nelements):
                    edata = data[n:n+68]

                    out = s.unpack(edata)
                    if self.debug4():
                        self.binary_debug.write('OEF_Shear - %s\n' % (str(out)))
                    (eid_device,
                     f41, f21, f12, f32, f23, f43, f34, f14, kf1,
                     s12, kf2, s23, kf3, s34, kf4, s41) = out
                    eid = (eid_device - self.device_code) // 10
                    assert eid > 0, eid

                    data_in = [eid,
                               f41, f21, f12, f32, f23, f43, f34,
                               f14, kf1, s12, kf2, s23, kf3, s34, kf4, s41]
                    #print "%s" % (self.get_element_type(self.element_type)), data_in
                    self.obj.add(dt, data_in)
                    n += ntotal

            elif self.format_code in [2, 3] and self.num_wide == 33:  # imag
                self.create_transient_object(self.cshear_force, ComplexCShearForce)
                s = Struct(self._endian + b'i32f')
                ntotal = 132  # 33*4
                nelements = len(data) // ntotal
                for i in range(nelements):
                    edata = data[n:n+132]
                    n += ntotal

                    out = s.unpack(edata)
                    if self.debug4():
                        self.binary_debug.write('OEF_Shear - %s\n' % (str(out)))
                    (eid_device,
                     f41r, f21r, f12r, f32r, f23r, f43r, f34r, f14r,
                     kf1r, s12r, kf2r, s23r, kf3r, s34r, kf4r, s41r,
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

                    eid = (eid_device - self.device_code) // 10
                    assert eid > 0, eid
                    data_in = [eid,
                               f41, f21, f12, f32, f23, f43, f34, f14,
                               kf1, s12, kf2, s23, kf3, s34, kf4, s41]
                    #print "%s" %(self.get_element_type(self.element_type)), data_in
                    self.obj.add(dt, data_in)
            else:
                msg = self.code_information()
                return self._not_implemented_or_skip(data, msg)

        elif self.element_type in [35]:  # con
            # 35-CON
            if self.read_mode == 1:
                return len(data)

            self._results._found_result('coneax_force')
            if self.format_code == 1 and self.num_wide == 7:  # real
                self.create_transient_object(self.coneax_force, RealConeAxForce)
                ntotal = 28  # 7*4
                s = Struct(self._endian + b'i6f')
                nelements = len(data) // ntotal
                for i in range(nelements):
                    edata = data[n:n+ntotal]
                    out = s.unpack(edata)
                    if self.debug4():
                        self.binary_debug.write('OEF_CONEAX-35 - %s\n' % (str(out)))
                    (eid_device, hopa, bmu, bmv, tm, su, sv) = out
                    eid = (eid_device - self.device_code) // 10
                    assert eid > 0, eid

                    data_in = [eid, hopa, bmu, bmv, tm, su, sv]
                    #print "%s" % (self.get_element_type(self.element_type)), data_in
                    self.obj.add(dt, data_in)
                    n += ntotal
            else:
                msg = self.code_information()
                return self._not_implemented_or_skip(data, msg)


        elif self.element_type == 38:  # cgap
            # 38-GAP
            if self.read_mode == 1:
                return len(data)
            self._results._found_result('cgap_force')
            if self.format_code == 1 and self.num_wide == 9:  # real
                self.create_transient_object(self.cgap_force, RealCGapForce)
                s = Struct(self._endian + b'i8f')
                ntotal = 36 # 9*4
                nelements = len(data) // ntotal
                for i in range(nelements):
                    edata = data[n:n+36]

                    out = s.unpack(edata)
                    if self.debug4():
                        self.binary_debug.write('OEF_CGAP-38 - %s\n' % (str(out)))
                    (eid_device, fx, sfy, sfz, u, v, w, sv, sw) = out
                    eid = (eid_device - self.device_code) // 10
                    assert eid > 0, eid

                    data_in = [eid, fx, sfy, sfz, u, v, w, sv, sw]
                    #print "%s" %(self.get_element_type(self.element_type)),data_in
                    #eid = self.obj.add_new_eid(out)
                    self.obj.add(dt, data_in)
                    n += ntotal
            else:
                msg = self.code_information()
                return self._not_implemented_or_skip(data, msg)
        elif self.element_type == 69:  # cbend
            # 69-CBEND
            if self.read_mode == 1:
                return len(data)
            self._results._found_result('cbend_force')
            if self.format_code == 1 and self.num_wide == 15:  # real
                self.create_transient_object(self.cbend_force, RealBendForce)
                s = Struct(self._endian + b'ii13f')

                ntotal = 60  # 15*4
                nelements = len(data) // ntotal
                for i in range(nelements):
                    edata = data[n:n+ntotal]

                    out = s.unpack(edata)
                    if self.debug4():
                        self.binary_debug.write('OEF_BEND-69 - %s\n' % (str(out)))
                    (eid_device,
                     nidA, bm1A, bm2A, ts1A, ts2A, afA, trqA,
                     nidB, bm1B, bm2B, ts1B, ts2B, afB, trqB) = out
                    eid = (eid_device - self.device_code) // 10
                    assert eid > 0, eid
                    data_in = [eid,
                               nidA, bm1A, bm2A, ts1A, ts2A, afA, trqA,
                               nidB, bm1B, bm2B, ts1B, ts2B, afB, trqB]
                    #print "%s" %(self.get_element_type(self.element_type)), data_in
                    #eid = self.obj.add_new_eid(out)
                    self.obj.add(dt, data_in)
                    n += ntotal
            elif self.format_code in [2, 3] and self.num_wide == 27:  # imag
                self.create_transient_object(self.cbend_force, ComplexBendForce)
                s = Struct(self._endian + b'ii25f')

                ntotal = 108  # 27*4
                nelements = len(data) // ntotal
                for i in range(nelements):
                    edata = data[n:n+108]
                    n += ntotal

                    out = s.unpack(edata)
                    if self.debug4():
                        self.binary_debug.write('OEF_BEND-69 - %s\n' % (str(out)))
                    (eid_device, nidA,
                     bm1Ar, bm2Ar, ts1Ar, ts2Ar, afAr, trqAr,
                     bm1Ai, bm2Ai, ts1Ai, ts2Ai, afAi, trqAi,
                     nidB,
                     bm1Br, bm2Br, ts1Br, ts2Br, afBr, trqBr,
                     bm1Bi, bm2Bi, ts1Bi, ts2Bi, afBi, trqBi) = out
                    eid = (eid_device - self.device_code) // 10
                    assert eid > 0, eid

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

                    data_in = [eid, nidA,
                               bm1A, bm2A, ts1A, ts2A, afA, trqA,
                               nidB,
                               bm1B, bm2B, ts1B, ts2B, afB, trqB]
                    #print "%s" %(self.get_element_type(self.element_type)), data_in
                    self.obj.add(dt, data_in)
            else:
                msg = self.code_information()
                return self._not_implemented_or_skip(data, msg)

        elif self.element_type in [76, 77, 78]:
            # 76-HEXPR
            # 77-PENPR
            # 78-TETPR
            if self.read_mode == 1:
                return len(data)

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
                return self._not_implemented_or_skip(data, msg)

            self._results._found_result(result_name)
            if self.format_code == 1 and self.num_wide == 10:  # real
                ntotal = 40
                nelements = len(data) // ntotal
                self.create_transient_object(slot, RealPentaPressureForce)
                s = Struct(self._endian + b'i8s7f')
                for i in range(nelements):
                    edata = data[n : n + 40]
                    n += 40
                    out = s.unpack(edata)
                    if self.debug4():
                        self.binary_debug.write('OEF_PentaPressure-%s %s\n' % (self.element_type, str(out)))
                    (eid_device, ename, ax, ay, az, vx, vy, vz, pressure) = out
                    eid = (eid_device - self.device_code) // 10
                    assert eid > 0, eid

                    data_in = [eid, ename, ax, ay, az, vx, vy, vz, pressure]
                    #print "%s" % (self.get_element_type(self.element_type)), data_in
                    self.obj.add(dt, data_in)

            elif self.format_code in [2, 3] and self.num_wide == 16:  # imag
                self.create_transient_object(slot, ComplexPentaPressureForce)

                s = Struct(self._endian + b'i8s13f')
                ntotal = 64
                nelements = len(data) // ntotal
                for i in range(nelements):
                    edata = data[n:n+64]
                    n += 64

                    out = s.unpack(edata)
                    if self.debug4():
                        self.binary_debug.write('OEF_PentaPressure-%s %s\n' % (self.element_type, str(out)))
                    (eid_device, eName,
                     axr, ayr, azr, vxr, vyr, vzr, pressure,
                     axi, ayi, azi, vxi, vyi, vzi) = out
                    eid = (eid_device - self.device_code) // 10
                    assert eid > 0, eid
                    ename = eName.decode('utf-8').strip()
                    #print "eType=%s" %(eType)

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

                    data_in = [eid, ename, ax, ay, az, vx, vy, vz, pressure]
                    #print "%s" % (self.get_element_type(self.element_type)) ,data_in
                    self.obj.add(dt, data_in)
            else:
                msg = self.code_information()
                return self._not_implemented_or_skip(data, msg)

        elif self.element_type == 100:  #  bars
            # 100-BARS
            if self.read_mode == 1:
                return len(data)
            self._results._found_result('cbar100_force')
            if self.format_code == 1 and self.num_wide == 8:
                self.create_transient_object(self.cbar100_force, RealCBar100Force)

                s = Struct(self._endian + b'i7f')
                ntotal = 32  # 8*4
                nelements = len(data) // ntotal
                for i in range(nelements):
                    edata = data[n:n+32]
                    n += 32
                    out = s.unpack(edata)
                    if self.debug4():
                        self.binary_debug.write('OEF_CBar100 - %s\n' % (str(out)))
                    (eid_device, sd, bm1, bm2, ts1, ts2, af, trq) = out
                    eid = (eid_device - self.device_code) // 10
                    assert eid > 0, eid

                    data_in = [eid, sd, bm1, bm2, ts1, ts2, af, trq]
                    #print "%s" %(self.get_element_type(self.element_type)), data_in
                    self.obj.add(dt, data_in)
            else:
                msg = self.code_information()
                return self._not_implemented_or_skip(data, msg)
        elif self.element_type == 102:  # cbush
            # 102-CBUSH
            self._results._found_result('cbush_force')
            if self.format_code == 1 and self.num_wide == 7:  # real
                self.create_transient_object(self.cbush_force, RealCBushForce)
                s = Struct(self._endian + b'i6f')
                ntotal = 28 # 7*4
                nelements = len(data) // ntotal
                for i in range(nelements):
                    edata = data[n:n+28]
                    out = s.unpack(edata)
                    if self.debug4():
                        self.binary_debug.write('OEF_CBUSH-102 - %s\n' % (str(out)))
                    (eid_device, fx, fy, fz, mx, my, mz) = out
                    eid = (eid_device - self.device_code) // 10
                    assert eid > 0, eid

                    data_in = [eid, fx, fy, fz, mx, my, mz]
                    #print "%s" % (self.get_element_type(self.element_type)), data_in
                    self.obj.add(dt, data_in)
                    n += ntotal
            elif self.format_code in [2, 3] and self.num_wide == 13:  # imag
                self.create_transient_object(self.cbush_force, ComplexCBushForce)
                s = Struct(self._endian + b'i12f')

                ntotal = 52  # 13*4
                nelements = len(data) // ntotal
                for i in range(nelements):
                    edata = data[n:n+52]

                    out = s.unpack(edata)
                    if self.debug4():
                        self.binary_debug.write('OEF_CBUSH-102 - %s\n' % (str(out)))
                    (eid_device,
                     fxr, fyr, fzr, mxr, myr, mzr,
                     fxi, fyi, fzi, mxi, myi, mzi) = out
                    eid = (eid_device - self.device_code) // 10
                    assert eid > 0, eid
                    #print "eType=%s" % (eType)

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

                    data_in = [eid, fx, fy, fz, mx, my, mz]
                    #print "%s" %(self.get_element_type(self.element_type)), data_in
                    self.obj.add(dt, data_in)
                    n += ntotal
            else:
                msg = self.code_information()
                return self._not_implemented_or_skip(data, msg)

        elif self.element_type in [145, 146, 147]:
            # 145-VUHEXA
            # 146-VUPENTA
            # 147-VUTETRA
            if self.read_mode == 1:
                return len(data)
            return len(data)
        elif self.element_type in [189, 190]:
            # 189-VUQUAD
            # 190-VUTRIA
            if self.read_mode == 1:
                return len(data)

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
                self.create_transient_object(self.force_VU_2D, RealForce_VU_2D)

                ntotal = 24 + 52 * nnodes
                nelements = len(data) // ntotal

                s1 = Struct(self._endian + b'iii4sii')
                s2 = Struct(self._endian + b'i3f3i5fi')
                for i in range(nelements):
                    edata = data[n:n+24]  # 6*4
                    n += 24

                    out = s1.unpack(edata)
                    if self.debug4():
                        self.binary_debug.write('OEF_Force_%s-%s - %s\n' % (etype, self.element_type, str(out)))
                    (eid_device, parent, coord, icord, theta, _) = out

                    eid = (eid_device - self.device_code) // 10
                    data_in = [eid, parent, coord, icord, theta]
                    assert eid > 0, eid

                    forces = []
                    for i in range(nnodes):
                        edata = data[n:n+52]  # 13*4
                        n += 52
                        out = s2.unpack(edata)
                        if self.debug4():
                            self.binary_debug.write('%s\n' % (str(out)))
                        (vugrid, mfx, mfy, mfxy, a, b, c, bmx, bmy,
                         bmxy, syz, szx, d) = out
                        out2 = (vugrid, mfx, mfy, mfxy, bmx, bmy, bmxy, syz, szx)
                        forces.append(out2)
                    data_in.append(forces)
                    #print "etype=%s" % (etype)

                    #data_in = [vugrid,mfx,mfy,mfxy,a,b,c,bmx,bmy,bmxy,syz,szx,d]
                    #print "force %s" %(self.get_element_type(self.element_type)),data_in
                    self.obj.add(nnodes, dt, data_in)

            elif self.format_code in [2, 3] and self.num_wide == numwide_imag:  # imag
                self.create_transient_object(self.force_VU_2D, ComplexForce_VU_2D)
                ntotal = 24 + 100 * nnodes
                s1 = Struct(self._endian + b'iii4sii')
                s2 = Struct(self._endian + b'i3f3i5fi3f3i5fi')
                nelements = len(data) // ntotal
                for i in range(nelements):
                    edata = data[n:n+24]  # 6*4
                    n += 24

                    out = s1.unpack(edata)
                    if self.debug4():
                        self.binary_debug.write('OEF_Force_%s-%s - %s\n' % (etype, self.element_type, str(out)))
                    (eid_device, parent, coord, icord, theta, _) = out

                    eid = (eid_device - self.device_code) // 10
                    assert eid > 0, eid
                    data_in = [eid, parent, coord, icord, theta]

                    forces = []
                    for i in range(nnodes):
                        edata = data[n:n+100]  # 13*4
                        n += 100
                        out = s2.unpack(edata)
                        if self.debug4():
                            self.binary_debug.write('OEF_Force_%s-%s - %s\n' % (etype, self.element_type, str(out)))
                        [vugrid, mfxr, mfyr, mfxyr, a, b, c, bmxr, bmyr, bmxyr, syzr, szxr, d,
                         mfxi, mfyi, mfxyi, a, b, c, bmxi, bmyi, bmxyi, syzi, szxi, d] = out

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
                    #print "eType=%s" %(eType)
                    #data_in = [vugrid,mfxr,mfyr,mfxyr,bmxr,bmyr,bmxyr,syzr,szxr,
                                     #mfxi,mfyi,mfxyi,bmxi,bmyi,bmxyi,syzi,szxi]
                    #print "force %s" %(self.get_element_type(self.element_type)),data_in
                    self.obj.add(nnodes, dt, data_in)
            else:
                msg = self.code_information()
                return self._not_implemented_or_skip(data, msg)


        elif self.element_type in [191]:
            # 191-VUBEAM
            if self.read_mode == 1:
                return len(data)
            self._results._found_result('force_VU')

            if self.format_code == 1 and self.num_wide == 20:  # real
                self.create_transient_object(self.force_VU, RealForce_VU)
                # 20 = 4 + 8 * 2 = 4 = 16
                nnodes = 2
                #ntotal = 16 + 32 * nnodes
                ntotal = self.num_wide * 4
                nelements = len(data) // ntotal
                self.create_transient_object(self.force_VU, RealForce_VU)

                s1 = Struct(self._endian + b'iii4s')
                s2 = Struct(self._endian + b'i7f')
                for i in range(nelements):
                    edata = data[n:n+16]  # 8*4
                    n += 16

                    out = s1.unpack(edata)
                    if self.debug4():
                        self.binary_debug.write('OEF_Force_VU-191 - %s\n' % (str(out)))
                    (eid_device, parent, coord, icord) = out

                    eid = (eid_device - self.device_code) // 10
                    assert eid > 0, eid
                    data_in = [eid, parent, coord, icord]

                    forces = []
                    for i in range(nnodes):
                        edata = data[n:n+32]  # 8*4
                        n += 32
                        #print "i=%s len(data)=%s" % (i, len(edata))
                        out = s2.unpack(edata)
                        if self.debug4():
                            self.binary_debug.write('%s\n' % str(out))
                        forces.append(out)
                    data_in.append(forces)
                    #print "eType=%s" % (eType)

                    #data_in = [vugrid, posit, forceX, shearY, shearZ, torsion, bendY, bendZ]
                    #print "force %s" %(self.get_element_type(self.element_type)), data_in
                    self.obj.add(nnodes, dt, data_in)
            elif self.format_code == 1 and self.num_wide == 32:  # random
                return len(data)
            elif self.format_code in [2, 3] and self.num_wide == 32:  # imag
                # 32 = 4 + 56/4 * 2 = 4 + 14 * 2 = 4 + 28
                self.create_transient_object(self.force_VU, ComplexForce_VU)
                nnodes = 2
                #ntotal = 16 + 56 * nnodes
                ntotal = self.num_wide * 4
                s1 = Struct(self._endian + b'i2i4s')
                s2 = Struct(self._endian + b'i13f')
                n = 0
                nelements = len(data) // ntotal
                for i in range(nelements):
                    edata = data[n:n+16]  # 8*4
                    n += 16

                    out = s1.unpack(edata)
                    if self.debug4():
                        self.binary_debug.write('OEF_Force_%s-%s - %s\n' % (etype, self.element_type, str(out)))
                    (eid_device, parent, coord, icord) = out

                    eid = (eid_device - self.device_code) // 10
                    assert eid > 0, eid
                    data_in = [eid, parent, coord, icord]

                    forces = []
                    for i in range(nnodes):
                        edata = data[n:n+56]  # 14*4
                        n += 56
                        out = s2.unpack(edata)
                        if self.debug4():
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
                    #print "eType=%s" %(eType)

                    #data_in = [vugrid,posit,forceX,shearY,shearZ,torsion,bendY,bendZ]
                    #print "force %s" %(self.get_element_type(self.element_type)),data_in
                    self.obj.add(nnodes, dt, data_in)
            else:
                msg = self.code_information()
                return self._not_implemented_or_skip(data, msg)

        elif self.element_type in [233, 235]:
            # 233-TRIARLC
            # 235-CQUADR
            if self.read_mode == 1:
                return len(data)
            return len(data)
        else:
            msg = self.code_information()
            return self._not_implemented_or_skip(data, msg)

        assert self.thermal == 0, self.thermal
        assert len(data) > 0, len(data)
        assert nelements > 0, 'nelements=%r element_type=%s element_name=%r' % (nelements, self.element_type, self.element_name)
        #assert len(data) % ntotal == 0, '%s n=%s nwide=%s len=%s ntotal=%s' % (self.element_name, len(data) % ntotal, len(data) % self.num_wide, len(data), ntotal)
        assert self.num_wide * 4 == ntotal, 'numwide*4=%s ntotal=%s' % (self.num_wide*4, ntotal)
        assert n > 0, n
        return n

    def _read_oef2_4(self, data):
        """Table 4 parser for OEF2 thermal table"""
        pass
