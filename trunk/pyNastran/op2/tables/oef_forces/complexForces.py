from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
from struct import Struct

from pyNastran.op2.op2_helper import polar_to_real_imag


class ComplexForces(object):

    def OEF_Rod_alt(self):  # 1-CROD, 3-CTUBE, 10-CONROD
        #device_code = self.device_code
        dt = self.nonlinear_factor
        (format1, extract) = self.getOEF_FormatStart()
        format1 += '4f'
        format1 = bytes(format1)
        is_magnitude_phase = self.is_magnitude_phase()

        n = 0
        s = Struct(format1)
        ntotal = 20 # 5*4
        nelements = len(self.data) // ntotal
        for i in xrange(nelements):
            edata = self.data[n:n+20]
            out = s.unpack(edata)
            (eid_device, axialReal, torqueReal, axialImag, torqueImag) = out

            if is_magnitude_phase:
                axial = polar_to_real_imag(axialReal, axialImag)
                torque = polar_to_real_imag(torqueReal, torqueImag)
            else:
                axial = complex(axialReal, axialImag)
                torque = complex(torqueReal, torqueImag)
            eid = extract(eid_device, dt)

            dataIn = [eid, axial, torque]
            #print "%s" %(self.get_element_type(self.element_type)),dataIn
            self.obj.add(dt, dataIn)
            n += ntotal
        self.data = self.data[n:]
        #print self.rodForces

    def OEF_Beam_alt(self):  # 2-CBEAM
        dt = self.nonlinear_factor
        (format1, extract) = self.getOEF_FormatStart()
        is_magnitude_phase = self.is_magnitude_phase()

        #print self.code_information()
        #nTotal = 16*11+1
        formatAll = 'i15f'
        format1 = bytes(format1)
        formatAll = bytes(formatAll)
        n = 0
        s1 = Struct(format1)
        s2 = Struct(formatAll)
        ntotal = 708  # (16*11+1)*4 = 177*4
        nelements = len(self.data) // ntotal
        for i in xrange(nelements):
            edata = self.data[n:n+4]
            eid_device, = s1.unpack(edata)
            eid = extract(eid_device, dt)

            n += 4
            for i in xrange(11):
                edata = self.data[n:n+64]
                n += 64

                out = s2.unpack(edata)
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
                #print "eidTemp = ",eidTemp
                #print "nid = ",nid
                #print "sd = ",sd

                #eid = self.obj.add_new_eid(out)
                if i == 0:  # isNewElement:
                    dataIn = [eid, nid, sd, bm1, bm2,
                              ts1, ts2, af, ttrq, wtrq]
                    #print "%s cNew   " %(self.get_element_type(self.element_type)),dataIn
                    self.obj.add_new_element(dt, dataIn)
                    #print
                elif sd > 0.:
                    dataIn = [eid, nid, sd, bm1, bm2,
                              ts1, ts2, af, ttrq, wtrq]
                    #print "%s cOld   " %(self.get_element_type(self.element_type)),dataIn
                    self.obj.add(dt, dataIn)
                    #print
                #else: pass
        self.data = self.data[n:]
        #print self.beamForces

    def OEF_Shear_alt(self):  # 4-CSHEAR
        dt = self.nonlinear_factor
        (format1, extract) = self.getOEF_FormatStart()
        format1 += '32f'
        format1 = bytes(format1)
        is_magnitude_phase = self.is_magnitude_phase()

        n = 0
        ntotal = 132  # 33*4
        nelements = len(self.data) // ntotal
        s = Struct(format1)
        for i in xrange(nelements):
            edata = self.data[n:n+132]
            n += ntotal

            out = s.unpack(edata)
            (eid_device,
             f41r, f21r, f12r, f32r, f23r, f43r, f34r, f14r, kf1r, s12r, kf2r, s23r, kf3r, s34r, kf4r, s41r,
             f41i, f21i, f12i, f32i, f23i, f43i, f34i, f14i, kf1i, s12i, kf2i, s23i, kf3i, s34i, kf4i, s41i) = out
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

            eid = extract(eid_device, dt)
            #print "eType=%s" % (eType)

            data_in = [eid, f41, f21, f12, f32, f23, f43, f34, f14,
                            kf1, s12, kf2, s23, kf3, s34, kf4, s41]
            #print "%s" %(self.get_element_type(self.element_type)), data_in
            self.obj.add(dt, data_in)
        #print self.shearForces

    def OEF_Spring_alt(self):  # 11-CELAS1, 12-CELAS2, 13-CELAS3, 14-CELAS4
        dt = self.nonlinear_factor
        (format1, extract) = self.getOEF_FormatStart()
        format1 += 'ff'
        format1 = bytes(format1)
        is_magnitude_phase = self.is_magnitude_phase()

        n = 0
        ntotal = 12
        s = Struct(format1)
        nelements = len(self.data) // 12  # 3*4
        for i in xrange(nelements):
            edata = self.data[n:n+12]

            out = s.unpack(edata)
            (eid_device, forceReal, forceImag) = out
            eid = extract(eid_device, dt)
            #print "eType=%s" %(eType)

            if is_magnitude_phase:
                force = polar_to_real_imag(forceReal, forceImag)
            else:
                force = complex(forceReal, forceImag)

            dataIn = [eid, force]
            #print "%s" %(self.get_element_type(self.element_type)),dataIn
            self.obj.add(dt, dataIn)
            n += ntotal
        self.data = self.data[n:]
        #print self.springForces

    def OEF_CVisc_alt(self):  # 24-CVISC
        dt = self.nonlinear_factor
        (format1, extract) = self.getOEF_FormatStart()
        format1 += '4f'
        format1 = bytes(format1)
        is_magnitude_phase = self.is_magnitude_phase()

        n = 0
        ntotal = 20  # 5*4
        s = Struct(format1)
        nelements = len(self.data) // ntotal
        for i in xrange(nelements):
            edata = self.data[n:n+20]

            out = s.unpack(edata)
            (eid_device, axialReal, torqueReal, axialImag, torqueImag) = out
            eid = extract(eid_device, dt)
            #print "eType=%s" %(eType)

            if is_magnitude_phase:
                axial = polar_to_real_imag(axialReal, axialImag)
                torque = polar_to_real_imag(torqueReal, torqueImag)
            else:
                axial = complex(axialReal, axialImag)
                torque = complex(torqueReal, torqueImag)

            dataIn = [eid, axial, torque]
            #print "%s" %(self.get_element_type(self.element_type)),dataIn
            self.obj.add(dt, dataIn)
            n += ntotal
        self.data = self.data[n:]
        #print self.viscForces

    def OEF_CBar_alt(self):  # 34-CBAR
        dt = self.nonlinear_factor
        (format1, extract) = self.getOEF_FormatStart()
        format1 += '16f'
        format1 = bytes(format1)
        is_magnitude_phase = self.is_magnitude_phase()
        n = 0
        ntotal = 68  # 17*4
        s = Struct(format1)
        nelements = len(self.data) // ntotal
        for i in xrange(nelements):
            edata = self.data[n:n+68]

            out = s.unpack(edata)
            (eid_device, bm1ar, bm2ar, bm1br, bm2br, ts1r, ts2r, afr, trqr,
                         bm1ai, bm2ai, bm1bi, bm2bi, ts1i, ts2i, afi, trqi) = out
            eid = extract(eid_device, dt)
            #print "eType=%s" %(eType)

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

            dataIn = [eid, bm1a, bm2a, bm1b, bm2b, ts1, ts2, af, trq]
            #print "%s" %(self.get_element_type(self.element_type)),dataIn
            self.obj.add(dt, dataIn)
            n += ntotal
        self.data = self.data[n:]
        #print self.barForces

    def OEF_Plate_alt(self):  # 33-CQUAD4,74-CTRIA3
        dt = self.nonlinear_factor
        (format1, extract) = self.getOEF_FormatStart()
        format1 += '16f'
        format1 = bytes(format1)
        s = Struct(format1)
        is_magnitude_phase = self.is_magnitude_phase()

        n = 0
        ntotal = 68
        nelements = len(self.data) // 68
        for i in xrange(nelements):
            edata = self.data[n:n+68]

            out = s.unpack(edata)
            (eid_device, mxr, myr, mxyr, bmxr, bmyr, bmxyr, txr, tyr,
                         mxi, myi, mxyi, bmxi, bmyi, bmxyi, txi, tyi) = out
            eid = extract(eid_device, dt)
            #print "eType=%s" %(eType)

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

            dataIn = [eid, mx, my, mxy, bmx, bmy, bmxy, tx, ty]
            #print "%s" %(self.get_element_type(self.element_type)),dataIn
            self.obj.add(dt, dataIn)
            n += ntotal
        self.data = self.data[n:]
        #print self.plateForces

    def OEF_Plate2_alt(self):  # 64-CQUAD8,70-CTRIAR,75-CTRIA6,82-CQUAD8,144-CQUAD4-bilinear
        dt = self.nonlinear_factor
        (format1, extract) = self.getOEF_FormatStart()
        format1 += '4s'
        is_magnitude_phase = self.is_magnitude_phase()

        if self.element_type in [70, 75]:  # CTRIAR,CTRIA6
            nNodes = 3
        elif self.element_type in [64, 82, 144]:  # CQUAD8,CQUADR,CQUAD4-bilinear
            nNodes = 4
        else:
            raise NotImplementedError(self.code_information())

        allFormat = '17f'
        format1 = bytes(format1)
        allFormat = bytes(allFormat)
        nTotal = 8 + (nNodes+1) * 68
        n = 0
        s1 = Struct(format1 + allFormat)
        s2 = Struct(allFormat)
        nelements = len(self.data) // nTotal
        for i in xrange(nelements):
            edata = self.data[n:n+76]
            n += 76

            out = s1.unpack(edata)
            (eid_device, term, nid, mxr, myr, mxyr, bmxr, bmyr, bmxyr, txr, tyr,
                                    mxi, myi, mxyi, bmxi, bmyi, bmxyi, txi, tyi) = out
            #term = 'CEN\'
            #print "eType=%s" %(eType)

            eid = extract(eid_device, dt)
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

            dataIn = [term, nid, mx, my, mxy, bmx, bmy, bmxy, tx, ty]
            #print "%s" %(self.get_element_type(self.element_type)),dataIn
            self.obj.add_new_element(eid, dt, dataIn)

            for i in xrange(nNodes):  # .. todo:: fix crash...
                edata = self.data[n:n+68]
                n += 68
                out = s2.unpack(edata)

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
                dataIn = [nid, mx, my, mxy, bmx, bmy, bmxy, tx, ty]
                #print "***%s    " %(self.get_element_type(self.element_type)),dataIn
                self.obj.add(eid, dt, dataIn)
        self.data = self.data[n:]
        #print self.plateForces2

    def OEF_Bend_alt(self):  # 69-CBEND
        dt = self.nonlinear_factor
        (format1, extract) = self.getOEF_FormatStart()
        format1 += 'i25f'
        format1 = bytes(format1)
        s = Struct(format1)
        is_magnitude_phase = self.is_magnitude_phase()

        n = 0
        ntotal = 108  # 27*4
        nelements = len(self.data) // ntotal
        for i in xrange(nelements):
            edata = self.data[n:n+108]
            n += ntotal

            out = s.unpack(edata)
            (eid_device,
             nidA, bm1Ar, bm2Ar, ts1Ar, ts2Ar, afAr, trqAr,
             bm1Ai, bm2Ai, ts1Ai, ts2Ai, afAi, trqAi,
             nidB, bm1Br, bm2Br, ts1Br, ts2Br, afBr, trqBr,
             bm1Bi, bm2Bi, ts1Bi, ts2Bi, afBi, trqBi) = out
            eid = extract(eid_device, dt)
            #print "eType=%s" % (eType)

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

            dataIn = [eid,
                      nidA, bm1A, bm2A, ts1A, ts2A, afA, trqA,
                      nidB, bm1B, bm2B, ts1B, ts2B, afB, trqB]
            #print "%s" %(self.get_element_type(self.element_type)),dataIn
            self.obj.add(dt, dataIn)
        #print self.bendForces

    def OEF_PentaPressure_alt(self):  # 76-CHEXA_PR,77-CPENTA_PR,78-CTETRA_PR
        dt = self.nonlinear_factor
        (format1, extract) = self.getOEF_FormatStart()
        format1 += '8s13f'
        format1 = bytes(format1)
        is_magnitude_phase = self.is_magnitude_phase()

        s = Struct(format1)
        n = 0
        ntotal = 64
        nelements = len(self.data) // ntotal
        for i in xrange(nelements):
            eData = self.data[n:n+64]
            n += 64

            out = s.unpack(eData)
            (eid_device, eName, axr, ayr, azr, vxr, vyr, vzr, pressure,
             axi, ayi, azi, vxi, vyi, vzi) = out
            eid = extract(eid_device, dt)
            eName = eName.decode('utf-8').strip()
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

            dataIn = [eid, eName, ax, ay, az, vx, vy, vz, pressure]
            #print "%s" %(self.get_element_type(self.element_type)),dataIn
            self.obj.add(dt, dataIn)
        self.data = self.data[n:]
        #print self.bendForces

    def OEF_CBush_alt(self):  # 102-CBUSH
        dt = self.nonlinear_factor
        (format1, extract) = self.getOEF_FormatStart()
        format1 += '12f'
        format1 = bytes(format1)
        is_magnitude_phase = self.is_magnitude_phase()

        n = 0
        ntotal = 52  # 13*4
        s = Struct(format1)
        nelements = len(self.data) // ntotal
        for i in xrange(nelements):
            edata = self.data[n:n+52]

            out = s.unpack(edata)
            (eid_device, fxr, fyr, fzr, mxr, myr, mzr,
                         fxi, fyi, fzi, mxi, myi, mzi) = out
            eid = extract(eid_device, dt)
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
        self.data = self.data[n:]
        #print self.bushForces

    def OEF_Force_VU_alt(self):  # 191-VUBEAM
        dt = self.nonlinear_factor
        (format1, extract) = self.getOEF_FormatStart()
        format1 += '2i4s'
        format1 = bytes(format1)
        is_magnitude_phase = self.is_magnitude_phase()

        if self.element_type in [191]:
            nNodes = 2
        else:
            raise NotImplementedError(self.code_information())

        formatAll = 'i13f'
        formatAll = bytes(formatAll)
        ntotal = 16 + 56 * nNodes
        s1 = Struct(format1)
        s2 = Struct(formatAll)
        n = 0
        nelements = len(self.data) // ntotal
        for i in xrange(nelements):
            eData = self.data[n:n+16]  # 8*4

            out = s1.unpack(eData)
            (eid_device, parent, coord, icord) = out

            eid = extract(eid_device, dt)
            dataIn = [eid, parent, coord, icord]

            forces = []
            for i in xrange(nNodes):
                eData = self.data[n:n+56]  # 14*4
                out = s2.unpack(eData)
                [vugrid, posit, forceXr, shearYr, shearZr, torsionr, bendingYr, bendingZr,
                 forceXi, shearYi, shearZi, torsioni, bendingYi, bendingZi] = out

                if is_magnitude_phase:
                    forceX = polar_to_real_imag(forceXr, forceXi)
                    shearY = polar_to_real_imag(shearYr, shearYi)
                    shearZ = polar_to_real_imag(shearZr, shearZi)
                    torsion = polar_to_real_imag(torsionr, torsioni)
                    bendingY = polar_to_real_imag(bendingYr, bendingYi)
                    bendingZ = polar_to_real_imag(bendingZr, bendingZi)
                else:
                    forceX = complex(forceXr, forceXi)
                    shearY = complex(shearYr, shearYi)
                    shearZ = complex(shearZr, shearZi)
                    torsion = complex(torsionr, torsioni)
                    bendingY = complex(bendingYr, bendingYi)
                    bendingZ = complex(bendingZr, bendingZi)

                out2 = [vugrid, posit, forceX, shearY,
                        shearZ, torsion, bendingY, bendingZ]
                forces.append(out2)
            dataIn.append(forces)
            #print "eType=%s" %(eType)

            #dataIn = [vugrid,posit,forceX,shearY,shearZ,torsion,bendY,bendZ]
            #print "force %s" %(self.get_element_type(self.element_type)),dataIn
            self.obj.add(nNodes, dt, dataIn)
        self.data = self.data[n:]

    def OEF_Force_VUTRIA_alt(self):  # 189-VUQUAD,190-VUTRIA
        dt = self.nonlinear_factor
        (format1, extract) = self.getOEF_FormatStart()
        format1 += 'ii4sii'
        format1 = bytes(format1)
        is_magnitude_phase = self.is_magnitude_phase()

        if self.element_type in [189]:  # VUQUAD
            nNodes = 4
        elif self.element_type in [190]:  # VUTRIA
            nNodes = 3
        else:
            raise NotImplementedError(self.code_information())

        formatAll = 'i3f3i5fi3f3i5fi'
        format1 = bytes(format1)
        formatAll = bytes(formatAll)
        n = 0
        s1 = Struct(format1)
        s2 = Struct(formatAll)
        ntotal = 24 + 100 * nNodes
        nelements = len(self.data) // ntotal
        while len(self.data) >= ntotal:
            eData = self.data[n:n+24]  # 6*4
            n += 24

            out = s1.unpack(eData)
            (eid_device, parent, coord, icord, theta, _) = out

            eid = extract(eid_device, dt)
            dataIn = [eid, parent, coord, icord, theta]

            forces = []
            for i in xrange(nNodes):
                eData = self.data[n:n+100]  # 13*4
                n += 100
                out = s2.unpack(eData)
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

            dataIn.append(forces)
            #print "eType=%s" %(eType)
            #dataIn = [vugrid,mfxr,mfyr,mfxyr,bmxr,bmyr,bmxyr,syzr,szxr,
                             #mfxi,mfyi,mfxyi,bmxi,bmyi,bmxyi,syzi,szxi]
            #print "force %s" %(self.get_element_type(self.element_type)),dataIn
            self.obj.add(nNodes, dt, dataIn)
        self.data = self.data[n:]
