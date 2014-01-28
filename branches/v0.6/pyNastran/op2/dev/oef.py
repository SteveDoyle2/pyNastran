from struct import Struct, unpack

from pyNastran.op2.op2_helper import polar_to_real_imag
from pyNastran.op2.dev.op2_common import OP2Common

from pyNastran.op2.tables.oef_forces.oef_forceObjects import (
    RealRodForce, RealCBeamForce, RealCShearForce,
    RealSpringForce, RealDamperForce, RealViscForce,
    RealPlateForce, RealConeAxForce, RealPlate2Force,
    RealCBar100Force, RealCGapForce, RealBendForce,
    RealPentaPressureForce, RealCBushForce,
    RealForce_VU_2D, RealCBarForce, RealForce_VU)
from pyNastran.op2.tables.oef_forces.oef_complexForceObjects import (
    ComplexRodForce, ComplexCBeamForce,
    ComplexCShearForce, ComplexSpringForce,
    ComplexDamperForce, ComplexViscForce,
    ComplexPlateForce, ComplexPlate2Force,
    ComplexBendForce,
    ComplexPentaPressureForce,
    ComplexCBushForce, ComplexForce_VU_2D,
    ComplexCBarForce, ComplexForce_VU)
from pyNastran.op2.tables.oef_forces.thermal_elements import ThermalElements

class OEF(OP2Common):
    def __init__(self):
        OP2Common.__init__(self)

    def OEF_ForceCode(self):
        """
        Gets the numwide codes for the element to determine if
        the real or complex result should be found.
        The format and sort codes do not always give the right answer...
        """
        realMapper = {
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
        imagMapper = {
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
            real = realMapper[self.element_type]
        except KeyError:
            real = None

        try:
            imag = imagMapper[self.element_type]
        except KeyError:
            imag = None
        return (real, imag)

    def read_oef1_3(self, data):
        self.words = ['aCode',       'tCode',    'element_type', 'isubcase',
                 '???',         '???',      '???',          '???',
                 'format_code', 'num_wide', 'o_code',       '???',
                 '???',         '???',      '???',          '???',
                 '???',         '???',      '???',          '???',
                 '???',         '???',      '???',          '???',
                 '???', 'Title', 'subtitle', 'label']

        self.parse_approach_code(data)

        #: element type
        self.element_type = self.add_data_parameter( data, 'element_type', 'i', 3, False)

        # dynamic load set ID/random code
        #self.dLoadID = self.add_data_parameter(data, 'dLoadID', 'i', 8, False)

        #: format code
        self.format_code = self.add_data_parameter( data, 'format_code', 'i', 9, False)

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
            self.eign = self.add_data_parameter(data, 'eign', 'f', 6, False)
            self.dataNames = self.apply_data_code_value('dataNames', ['mode', 'eigr', 'mode_cycle'])
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
            self.eigr = self.add_data_parameter( data, 'eigr', 'f', 6, False)
            self.dataNames = self.apply_data_code_value('dataNames', ['lsdvmn', 'eigr'])
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
            raise RuntimeError('invalid analysis_code...analysis_code=%s' % str(self.analysis_code))

        self.element_name = self.element_mapper[self.element_type]
        if self.debug:
            self.binary_debug.write('  element_name = %r\n' % self.element_name)
            self.binary_debug.write('  aCode    = %r\n' % self.aCode)
            self.binary_debug.write('  tCode    = %r\n' % self.tCode)
            self.binary_debug.write('  isubcase = %r\n' % self.isubcase)

        self.read_title(data)
        if self.element_type not in self.element_mapper:
            raise NotImplementedError(self.element_type)
        self.write_debug_bits()

    def read_oef2_3(self, data):
        pass

    def read_oef1_4(self, data):
        if self.thermal == 0:
            self._read_oef1_loads(data)
        elif self.thermal == 1:
            self._read_oef1_thermal(data)
        else:
            self.not_implemented_or_skip('thermal=%s' % self.thermal)

    def _read_oef1_thermal(self, data):
        n = 0
        is_magnitude_phase = self.is_magnitude_phase()
        dt = self.nonlinear_factor

        if self.element_type in [1, 2, 3, 10, 34, 69]:
            # 1-CROD
            # 2-CBEAM
            # 3-CTUBE
            # 10-CONROD
            # 34-CBAR
            # 69-CBEND:
            if self.num_wide == 9:
                ntotal = 36
                format1 = b'i8s6f'  # SORT1

                ntotal = 36  # 10*4
                s = Struct(format1)
                nelements = len(data) // ntotal
                for i in xrange(nelements):
                    edata = data[n:n+ntotal]

                    out = s.unpack(edata)
                    (eid_device, eType, xGrad, yGrad, zGrad, xFlux, yFlux, zFlux) = out
                    eid = (eid_device - self.device_code) // 10

                    data_in = [eid, eType, xGrad, yGrad, zGrad, xFlux, yFlux, zFlux]
                    #print "heatFlux %s" % (self.get_element_type(self.element_type)), data_in
                    eid = self.obj.add_new_eid(out)
                    self.obj.add(dt, data_in)
                    n += ntotal
            else:
                raise NotImplementedError(self.num_wide)

        elif self.element_type in [33, 39, 53, 64, 67, 68, 74, 75]:
            # 33-CQUAD4-centroidal
            # 39-CTETRA
            # 53-CTRIAX6
            # 67-CHEXA
            # 64-QUAD8
            # 74-CTRIA3-centroidal
            # 75-TRIA6
            # 33-CQUAD4-centroidal
            # 68-CPENTA
            return
        elif self.element_type in [107, 108, 109, 110, 145, 146,
                147, 189, 190, 191]:
            # 107-CHBDYE
            # 108-CHBDYG
            # 109-CHBDYP
            # 110-CONV
            # 145-VUHEXA
            # 146-VUPENTA
            # 147-VUTETRA
            # 189-VUQUAD
            # 190-VUTRIA
            # 191-VUBEAM
            return
        else:
            raise NotImplementedError('OEF sort1 thermal Type=%s num=%s' % (self.element_name, self.element_type))

        assert len(data) > 0
        assert nelements > 0, 'nelements=%r element_type=%s element_name=%r' % (nelements, self.element_type, self.element_name)
        assert len(data) % ntotal == 0, '%s n=%s nwide=%s len=%s ntotal=%s' % (self.element_name, len(data) % ntotal, len(data) % self.num_wide, len(data), ntotal)
        assert self.num_wide * 4 == ntotal, 'numwide*4=%s ntotal=%s' % (self.num_wide*4, ntotal)
        assert n > 0, n

    def _read_oef1_loads(self, data):
        (num_wide_real, num_wide_imag) = self.OEF_ForceCode()
        if self.debug4():
            self.binary_debug.write('  num_wide_real = %r\n' % num_wide_real)
            self.binary_debug.write('  num_wide_imag = %r\n' % num_wide_imag)

        n = 0
        is_magnitude_phase = self.is_magnitude_phase()
        dt = self.nonlinear_factor

        if self.element_type in []:
            pass
        elif self.element_type in [1, 3, 10]:
            #1-CROD
            #3-CTUBE
            #10-CONROD
            if self.num_wide == 3: # real
                self.create_transient_object(self.rodForces, RealRodForce)
                format1 = b'iff' # 3
                ntotal = 12 # 3 * 4
                nelements = len(data) // ntotal
                s = Struct(format1)
                for i in xrange(nelements):
                    edata = data[n:n+ntotal]
                    out = s.unpack(edata)
                    (eid_device, axial, torque) = out
                    eid = (eid_device - self.device_code) // 10
                    if self.debug4():
                        self.binary_debug.write('OEF_Rod - %s\n' % (str(out)))

                    data_in = [eid, axial, torque]
                    #print "%s" % (self.get_element_type(self.element_type)), data_in
                    self.obj.add(dt, data_in)
                    n += ntotal

            elif self.num_wide == 5: # imag
                self.create_transient_object(self.rodForces, ComplexRodForce)

                format1 = b'i4f'
                ntotal = 20 # 5*4
                nelements = len(data) // ntotal
                for i in xrange(nelements):
                    edata = data[n:n+20]
                    out = unpack(format1, edata)
                    (eid_device, axial_real, torque_real, axial_imag, torque_imag) = out

                    if is_magnitude_phase:
                        axial = polar_to_real_imag(axial_real, axial_imag)
                        torque = polar_to_real_imag(torque_real, torque_imag)
                    else:
                        axial = complex(axial_real, axial_imag)
                        torque = complex(torque_real, torque_imag)
                    eid = (eid_device - self.device_code) // 10

                    data_in = [eid, axial, torque]
                    #print "%s" % (self.get_element_type(self.element_type)), data_in
                    eid = self.obj.add_new_eid(out)
                    self.obj.add(dt, data_in)
                    n += ntotal
                #print self.rodForces

            else:
                raise NotImplementedError(self.num_wide)
            #print self.rodForces

        elif self.element_type in [2]:
            #2-CBEAM
            if self.num_wide == 9:  # centroid ???
                self.create_transient_object(self.beamForces, RealCBeamForce)
                format1 = b'i8f'  # 36
                ntotal = 36
                nelements = len(data) // ntotal
                for i in xrange(nelements):
                    edata = data[n:n+36]
                    out = unpack(format1, edata)
                    if self.debug4():
                        self.binary_debug.write('OEF_Beam - %s\n' % (str(out)))
                    (eid_device, sd, bm1, bm2, ts1, ts2, af, ttrq, wtrq) = out
                    eid = (eid_device - self.device_code) // 10
                    n += 36

            elif self.num_wide == 100:  # real
                self.create_transient_object(self.beamForces, RealCBeamForce)
                format1 = b'i'
                formatAll = b'i8f'  # 36

                ntotal = 400  # 1+(10-1)*11=100 ->100*4 = 400
                nelements = len(data) // ntotal
                for i in xrange(nelements):
                    edata = data[n:n+4]
                    eid_device, = unpack(format1, edata)
                    eid = (eid_device - self.device_code) // 10
                    n += 4

                    for i in xrange(11):
                        edata = data[n:n+36]
                        out = unpack(formatAll, edata)
                        if self.debug4():
                            self.binary_debug.write('OEF_Beam - %s\n' % (str(out)))
                        (nid, sd, bm1, bm2, ts1, ts2, af, ttrq, wtrq) = out

                        data_in = [eid, nid, sd, bm1, bm2, ts1, ts2, af, ttrq, wtrq]
                        if i == 0:  # isNewElement
                            self.obj.addNewElement(dt, data_in)
                        elif sd > 0.:
                            self.obj.add(dt, data_in)
                        n += 36
            elif self.num_wide == 177: # imag
                self.create_transient_object(self.beamForces, ComplexCBeamForce)
                formatAll = b'i15f'
                format1 = b'i'
                ntotal = 708  # (16*11+1)*4 = 177*4
                nelements = len(data) // ntotal
                for i in xrange(nelements):
                    edata = data[n:n+4]
                    eid_device, = unpack(format1, edata)
                    eid = (eid_device - self.device_code) // 10

                    n += 4
                    for i in xrange(11):
                        edata = data[n:n+64]
                        n += 64

                        out = unpack(formatAll, edata)
                        (nid, sd, bm1r, bm2r, ts1r, ts2r, afr, ttrqr, wtrqr,
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
                        eid = self.obj.add_new_eid(out)
                        if i == 0:  # isNewElement:
                            data_in = [eid, nid, sd, bm1, bm2,
                                       ts1, ts2, af, ttrq, wtrq]
                            #print "%s cNew   " % (self.get_element_type(self.element_type)), data_in
                            self.obj.addNewElement(dt, data_in)
                        elif sd > 0.:
                            data_in = [eid, nid, sd, bm1, bm2,
                                      ts1, ts2, af, ttrq, wtrq]
                            #print "%s cOld   " % (self.get_element_type(self.element_type)), data_in
                            self.obj.add(dt, data_in)
            else:
                raise NotImplementedError(self.num_wide)
            #print self.beamForces

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
            if self.num_wide == 2:
                if self.element_type in [11, 12, 13, 14]:
                    self.create_transient_object(self.springForces, RealSpringForce)
                elif self.element_type in [20, 21, 22, 23]:
                    self.create_transient_object(self.damperForces, RealDamperForce)
                else:
                    raise NotImplementedError(self.element_type)

                format1 = b'if' # 2
                s = Struct(format1)
                ntotal = 8  # 2*4
                nelements = len(data) // ntotal
                for i in xrange(nelements):
                    edata = data[n:n+8]

                    out = s.unpack(edata)
                    if self.debug4():
                        self.binary_debug.write('OEF_Spring - %s\n' % (str(out)))
                    (eid_device, force) = out
                    eid = (eid_device - self.device_code) // 10

                    data_in = [eid, force]
                    #print "%s" % (self.get_element_type(self.element_type)), data_in
                    self.obj.add(dt, data_in)
                    n += ntotal
            elif self.num_wide == 3:
                if self.element_type in [11, 12, 13, 14]:
                    self.create_transient_object(self.springForces, ComplexSpringForce)
                elif self.element_type in [20, 21, 22, 23]:
                    self.create_transient_object(self.damperForces, ComplexDamperForce)
                else:
                    raise NotImplementedError(self.element_type)

                format1 = b'i2f'
                ntotal = 12  # 3*4

                nelements = len(data) // ntotal
                for i in xrange(nelements):
                    edata = data[n:n+12]
                    out = unpack(format1, edata)
                    (eid_device, forceReal, forceImag) = out
                    eid = (eid_device - self.device_code) // 10

                    if is_magnitude_phase:
                        force = polar_to_real_imag(forceReal, forceImag)
                    else:
                        force = complex(forceReal, forceImag)

                    data_in = [eid, force]
                    #print "%s" % (self.get_element_type(self.element_type)), data_in
                    #eid = self.obj.add_new_eid(out)
                    self.obj.add(dt, data_in)
                    n += ntotal
            else:
                raise NotImplementedError()
            #print self.springForces

        elif self.element_type in [24]:  # CVISC
            if self.num_wide == 3: # real
                self.create_transient_object(self.viscForces, RealViscForce)
                format1 = b'iff'
                ntotal = 12  # 3*4
                nelements = len(data) // 12
                for i in xrange(nelements):
                    edata = data[n:n+12]

                    out = unpack(format1, edata)
                    if self.debug4():
                        self.binary_debug.write('OEF_CVisc - %s\n' % (str(out)))
                    (eid_device, axial, torque) = out
                    eid = (eid_device - self.device_code) // 10

                    data_in = [eid, axial, torque]
                    #print "%s" % (self.get_element_type(self.element_type)), data_in
                    #eid = self.obj.add_new_eid(out)
                    self.obj.add(dt, data_in)
                    n += ntotal
            elif self.num_wide == 5: # complex
                self.create_transient_object(self.viscForces, ComplexViscForce)
                format1 = b'i4f'  # 5
                ntotal = 20  # 5*4
                nelements = len(data) // ntotal
                for i in xrange(nelements):
                    edata = data[n:n+20]

                    out = unpack(format1, edata)
                    (eid_device, axial_real, torque_real, axial_imag, torque_imag) = out
                    eid = (eid_device - self.device_code) // 10

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
                raise NotImplementedError(self.num_wide)
            #print self.viscForces

        elif self.element_type in [34]:  # bars
            # 34-CBAR
            if self.num_wide == 9: # real
                self.create_transient_object(self.barForces, RealCBarForce)
                format1 = b'i8f'  # 9
                s = Struct(format1)
                ntotal = 36  # 9*4
                nelements = len(data) // ntotal
                for i in xrange(nelements):
                    edata = data[n:n+36]

                    out = s.unpack(edata)
                    if self.debug4():
                        self.binary_debug.write('OEF_CBar - %s\n' % (str(out)))
                    (eid_device, bm1a, bm2a, bm1b, bm2b, ts1, ts2, af, trq) = out
                    eid = (eid_device - self.device_code) // 10

                    data_in = [eid, bm1a, bm2a, bm1b, bm2b, ts1, ts2, af, trq]
                    #print "%s" % (self.get_element_type(self.element_type)), data_in
                    #eid = self.obj.add_new_eid(out)
                    self.obj.add(dt, data_in)
                    n += ntotal
            elif self.num_wide == 17: # imag
                self.create_transient_object(self.barForces, ComplexCBarForce)
                format1 = b'i16f'
                ntotal = 68  # 17*4
                nelements = len(data) // ntotal
                for i in xrange(nelements):
                    edata = data[n:n+68]

                    out = unpack(format1, edata)
                    (eid_device, bm1ar, bm2ar, bm1br, bm2br, ts1r, ts2r, afr, trqr,
                                 bm1ai, bm2ai, bm1bi, bm2bi, ts1i, ts2i, afi, trqi) = out
                    eid = (eid_device - self.device_code) // 10

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
                raise NotImplementedError(self.num_wide)
            #print self.barForces

        elif self.element_type in [33, 74]: # centroidal shells
            # 33-CQUAD4
            # 74-CTRIA3
            if self.num_wide == 9:
                self.create_transient_object(self.plateForces, RealPlateForce)
                format1 = b'i8f'
                s = Struct(format1)
                ntotal = 36 # 9*4
                nelements = len(data) // ntotal
                for i in xrange(nelements):
                    edata = data[n:n+36]

                    out = s.unpack(edata)
                    if self.debug4():
                        self.binary_debug.write('OEF_Plate-%s - %s\n' % (self.element_type, str(out)))
                    (eid_device, mx, my, mxy, bmx, bmy, bmxy, tx, ty) = out
                    eid = (eid_device - self.device_code) // 10
                    assert eid > 0, 'eid_device=%s eid=%s table_name-%r' % (eid_device, eid, self.table_name)

                    data_in = [eid, mx, my, mxy, bmx, bmy, bmxy, tx, ty]
                    #print "%s" % (self.get_element_type(self.element_type)), data_in
                    #eid = self.obj.add_new_eid(out)
                    self.obj.add(dt, data_in)
                    n += ntotal
            elif self.num_wide == 17:
                self.create_transient_object(self.plateForces, CopmplexPlateForce)  # undefined
                format1 = b'i16f'
                s = Struct(format1)

                ntotal = 68
                nelements = len(data) // ntotal
                for i in xrange(nelements):
                    edata = data[n:n+68]
                    out = s.unpack(edata)
                    (eid_device, mxr, myr, mxyr, bmxr, bmyr, bmxyr, txr, tyr,
                                 mxi, myi, mxyi, bmxi, bmyi, bmxyi, txi, tyi) = out
                    eid = (eid_device - self.device_code) // 10
                    assert eid > 0
                    if self.debug4():
                        self.binary_debug.write('OEF_Plate-%s - %s\n' % (self.element_type, str(out)))

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

                    data_in = [eid, mx, my, mxy, bmx, bmy, bmxy, tx, ty]
                    #print "%s" % (self.get_element_type(self.element_type)), data_in
                    #eid = self.obj.add_new_eid(out)
                    self.obj.add(dt, data_in)
                    n += ntotal
            else:
                raise NotImplementedError(self.num_wide)
            #print self.plateForces

        elif self.element_type in [64, 70, 75, 82, 144]: # bilinear shells
            # 64-CQUAD8
            # 70-CTRIAR
            # 75-CTRIA6,
            # 82-CQUAD8,
            # 144-CQUAD4-bilinear
            if self.element_type in [70, 75]:  # CTRIAR,CTRIA6
                nnodes = 3
            elif self.element_type in [64, 82, 144]:  # CQUAD8,CQUADR,CQUAD4-bilinear
                nnodes = 4
            else:
                raise NotImplementedError('name=%r type=%r' % (self.element_name, self.element_type))

            numwide_real = 2 + (nnodes + 1) * 9 # centroidal node is the + 1
            numwide_imag = 2 + (nnodes + 1) * 17

            if self.num_wide == numwide_real:  # real
                self.create_transient_object(self.plateForces2, RealPlate2Force)
                allFormat = b'i8f' # 36
                format1 = b'i4s'  # 8
                ntotal = 8 + (nnodes+1) * 36 # centroidal node is the + 1
                assert ntotal == self.num_wide * 4, 'ntotal=%s numwide=%s' % (ntotal, self.num_wide * 4)
                nelements = len(data) // ntotal

                for i in xrange(nelements):
                    edata = data[n:n+44]

                    out = unpack(format1 + allFormat, edata)
                    if self.debug4():
                        self.binary_debug.write('OEF_Plate2-%s - %s\n' % (self.element_type, str(out)))
                    (eid_device, term, nid, mx, my, mxy, bmx, bmy, bmxy, tx, ty) = out
                    #term= 'CEN\'

                    eid = (eid_device - self.device_code) // 10
                    assert eid > 0, eid
                    data_in = [term, nid, mx, my, mxy, bmx, bmy, bmxy, tx, ty]
                    #print "%s" % (self.get_element_type(self.element_type)), data_in
                    self.obj.addNewElement(eid, dt, data_in)
                    n += 44
                    for i in xrange(nnodes):
                        edata = data[n : n + 36]
                        out = unpack(allFormat, edata)
                        if self.debug4():
                            self.binary_debug.write('%s\n' % (str(out)))
                        (nid, mx, my, mxy, bmx, bmy, bmxy, tx, ty) = out
                        assert nid > 0, 'nid=%s' % nid
                        #data_in = [nid, mx, my, mxy, bmx, bmy, bmxy, tx, ty]
                        #print "***%s    " % (self.get_element_type(self.element_type)), data_in
                        self.obj.add(eid, dt, out)
                        n += 36
            elif self.num_wide == num_wide_imag: # complex
                self.create_transient_object(self.plateForces2, ComplexPlate2Force)
                format1 = b'i4s'
                allFormat = b'17f'
                ntotal = 8 + (nnodes+1) * 68
                nelements = len(data) // ntotal
                for i in xrange(nelements):
                    edata = data[n:n+76]
                    n += 76

                    out = unpack(format1 + allFormat, edata)
                    (eid_device, term, nid, mxr, myr, mxyr, bmxr, bmyr, bmxyr, txr, tyr,
                                            mxi, myi, mxyi, bmxi, bmyi, bmxyi, txi, tyi) = out
                    #term = 'CEN\'

                    eid = (eid_device - self.device_code) // 10
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
                    self.obj.addNewElement(eid, dt, data_in)

                    for i in xrange(nnodes):  # .. todo:: fix crash...
                        edata = data[n:n+68]
                        n += 68
                        out = unpack(allFormat, edata)

                        (nid, mxr, myr, mxyr, bmxr, bmyr, bmxyr, txr, tyr,
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
                        data_in = [nid, mx, my, mxy, bmx, bmy, bmxy, tx, ty]
                        #print "***%s    " % (self.get_element_type(self.element_type)),data_in
                        self.obj.add(eid, dt, data_in)
            else:
                raise NotImplementedError(self.num_wide)

        elif self.element_type in [95, 96, 97, 98]: # composites
            # 95 - CQUAD4
            # 96 - CQUAD8
            # 97 - CTRIA3
            # 98 - CTRIA6 (composite)
            if self.num_wide == 9:  # real
                print self.code_information()
                self.create_transient_object(self.compositePlateForce, RealCompositePlateForce)  # undefined
                #return
                ntotal = 9 * 4
                nelements = len(data) // ntotal
                if self.debug:
                    self.binary_debug.write('  [cap, element1, element2, ..., cap]\n')
                    self.binary_debug.write('  cap = %i  # assume 1 cap when there could have been multiple\n' % len(data))
                    #self.binary_debug.write('  #centeri = [eid_device, j, grid, fd1, sx1, sy1, txy1, angle1, major1, minor1, vm1,\n')
                    #self.binary_debug.write('  #                                fd2, sx2, sy2, txy2, angle2, major2, minor2, vm2,)]\n')
                    #self.binary_debug.write('  #nodeji = [eid, iLayer, o1, o2, t12, t1z, t2z, angle, major, minor, ovm)]\n')
                    self.binary_debug.write('  nelements=%i; nnodes=1 # centroid\n' % nelements)

                eid_old = 0
                format1 = 'i8si4f4s' # 9
                s = Struct(format1)
                for i in xrange(nelements):
                    if i % 10000 == 0:
                        print 'i = ', i
                    edata = data[n:n+ntotal]  # 4*9
                    out = s.unpack(edata)
                    (eid_device, theory, lamid, failure_index_direct_stress, failure_mode_max_shear,
                             failure_index_interlaminar_shear, fmax, failure_flag) = out
                    eid = (eid_device - self.device_code) // 10
                    if self.debug4():
                        if eid > 0:
                            self.binary_debug.write('  eid=%i; C=[%s]\n' % (', '.join(['%r' % di for di in out]) ))
                        else:
                            self.binary_debug.write('      %s  C=[%s]\n' % (' ' * len(str(eid)), ', '.join(['%r' % di for di in out]) ))

                    if eid > 0:
                        self.obj.add_new_eid(eType, dt, eid, o1, o2, t12, t1z, t2z, angle, major, minor, ovm)
                    else:
                        self.obj.add(dt, eid, o1, o2, t12, t1z, t2z, angle, major, minor, ovm)
                    eid_old = eid
                    n += ntotal
            else:
                raise NotImplementedError(self.num_wide)

        elif self.element_type in [67, 68]: # solids
            if self.num_wide == 0:
                self.create_transient_object(self.shearForces, RealCShearForce)
            else:
                raise NotImplementedError(self.num_wide)

            # 67-CHEXA
            # 68-CPENTA
            return
        elif self.element_type in [53]:
            # 53-CTRIAX6
            if self.num_wide == 0:
                self.create_transient_object(self.ctriaxForce, RealCTriaxForce)  # undefined
            else:
                raise NotImplementedError(self.num_wide)
            return
        elif self.element_type in [4]:
            # 4-CSHEAR
            if self.num_wide == 0:
                self.create_transient_object(self.shearForces, RealCShearForce)
            else:
                raise NotImplementedError(self.num_wide)
            return
        elif self.element_type in [35]:
            # 35-CON
            return
        elif self.element_type in [38]:
            # 38-GAP
            if self.num_wide == 0:
                self.create_transient_object(self.gapForces, RealCGapForce)
            else:
                raise NotImplementedError(self.num_wide)
            return
        elif self.element_type in [69]:
            # 69-CBEND
            if self.num_wide == 0:
                self.create_transient_object(self.bendForces, RealBendForce)
            else:
                raise NotImplementedError(self.num_wide)
            return
        elif self.element_type in [76, 77, 78]:
            # 76-HEXPR
            # 77-PENPR
            # 78-TETPR
            return
        elif self.element_type in [100]:
            # 100-BARS
            if self.num_wide == 0:
                self.create_transient_object(self.bar100Forces, RealCBar100Force)
            else:
                raise NotImplementedError(self.num_wide)
            return
        elif self.element_type in [102]:
            # 102-CBUSH
            if self.num_wide == 0:
                self.create_transient_object(self.bushForces, RealCBushForce)
            else:
                raise NotImplementedError(self.num_wide)
            return
        elif self.element_type in [145, 146, 147]:
            # 145-VUHEXA
            # 146-VUPENTA
            # 147-VUTETRA
            return
        elif self.element_type in [189, 190]:
            # 189-VUQUAD
            # 190-VUTRIA
            return
        elif self.element_type in [191, 233, 235]:
            # 191-VUBEAM
            # 233-TRIARLC
            # 235-CQUADR
            return
        else:
            raise NotImplementedError('OEF sort1 Type=%s num=%s' % (self.element_name, self.element_type))
        assert len(data) > 0
        assert nelements > 0, 'nelements=%r element_type=%s element_name=%r' % (nelements, self.element_type, self.element_name)
        assert len(data) % ntotal == 0, '%s n=%s nwide=%s len=%s ntotal=%s' % (self.element_name, len(data) % ntotal, len(data) % self.num_wide, len(data), ntotal)
        assert self.num_wide * 4 == ntotal, 'numwide*4=%s ntotal=%s' % (self.num_wide*4, ntotal)
        assert n > 0, n

    def read_oef2_4(self, data):
        pass
