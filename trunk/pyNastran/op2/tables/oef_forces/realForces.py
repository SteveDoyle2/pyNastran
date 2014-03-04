from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
from struct import Struct


class RealForces(object):

    def OEF_aCode(self):
        if self.analysis_code == 1:   # statics
            format1 = 'i'  # loadID
            #self.add_data_parameter(data,'loadID','i',5,False)   ## load set ID number
        elif self.analysis_code == 2:  # normal modes/buckling (real eigenvalues)
            format1 = 'i'  # mode
            #self.add_data_parameter(data,'mode','i',5)   ## mode number
            #self.add_data_parameter(data,'eign','f',6,False)   ## eigenvalue
        elif self.analysis_code == 3:  # differential stiffness 0
            format1 = 'i'  # loadID
            #self.add_data_parameter(data,'loadID','i',5)   ## load set ID number
        elif self.analysis_code == 4:  # differential stiffness 1
            format1 = 'i'  # loadID
            ## load set ID number
            self.add_data_parameter(data, 'loadID', 'i', 5)
        elif self.analysis_code == 5:   # frequency
            format1 = 'f'  # freq
            #self.add_data_parameter(data,'freq','f',5)   ## frequency

        elif self.analysis_code == 6:  # transient
            format1 = 'f'  # time
            #self.add_data_parameter(data,'time','f',5)   ## time step
            #print "time(5)=%s" %(self.time)
        elif self.analysis_code == 7:  # pre-buckling
            format1 = 'i'  # loadID
            #self.add_data_parameter(data,'loadID','i',5)   ## load set ID number
            #print "loadID(5)=%s" %(self.loadID)
        elif self.analysis_code == 8:  # post-buckling
            format1 = 'i'  # loadID
            #self.add_data_parameter(data,'loadID','i',5)       ## load set ID number
            #self.add_data_parameter(data,'eigr','f',6,False)   ## real eigenvalue
            #print "loadID(5)=%s  eigr(6)=%s" %(self.loadID,self.eigr)
        elif self.analysis_code == 9:  # complex eigenvalues
            format1 = 'i'  # mode
            #self.add_data_parameter(data,'mode','i',5)         ## mode number
            #self.add_data_parameter(data,'eigr','f',6,False)   ## real eigenvalue
            #self.add_data_parameter(data,'eigi','f',7,False)   ## imaginary eigenvalue
            #print "mode(5)=%s  eigr(6)=%s  eigi(7)=%s" %(self.mode,self.eigr,self.eigi)
        elif self.analysis_code == 10:  # nonlinear statics
            format1 = 'f'  # load_step
            ## load step
            self.add_data_parameter(data, 'load_step', 'f', 5)
            #print "load_step(5) = %s" %(self.load_step)
        elif self.analysis_code == 11:  # geometric nonlinear statics
            format1 = 'i'  # loadID
            #self.add_data_parameter(data,'loadID','i',5)   ## load set ID number
            #print "loadID(5)=%s" %(self.loadID)
        else:
            raise RuntimeError('invalid analysis_code...analysis_code=%s' % (str(self.analysis_code) + '\n' + self.code_information()))
        return format1

    def getOEF_FormatStart(self):
        """
        Returns an i or an f depending on if it's SORT2 or not.
        Also returns an extraction function that is called on the first argument
        """
        is_sort1 = self.is_sort1()
        if is_sort1:
            #print "SORT1 - %s" %(self.get_element_type(self.element_type))
            format1 = 'i'  # SORT1
            extract = self.extractSort1
            #dt = self.nonlinear_factor
        else:
            #print "SORT2 - %s" %(self.get_element_type(self.element_type))
            format1 = 'f'  # SORT2
            extract = self.extractSort2
            #eid = self.nonlinear_factor
        return (format1, extract)

    def OEF_Rod(self):  # 1-CROD, 3-CTUBE, 10-CONROD
        dt = self.nonlinear_factor
        (format1, extract) = self.getOEF_FormatStart()
        format1 += 'ff'
        format1 = bytes(format1)

        n = 0
        s = Struct(format1)
        nelements = len(self.data) // 12  # 3*4
        for i in xrange(nelements):
            edata = self.data[n:n+12]
            out = s.unpack(edata)
            if self.make_op2_debug:
                self.op2_debug.write('OEF_Rod - %s\n' % (str(out)))
            (eid_device, axial, torque) = out
            eid = extract(eid_device, dt)
            #print "eType=%s" %(eType)

            #print "%s" %(self.get_element_type(self.element_type)),dataIn
            self.obj.add(dt, eid, axial, torque)
            n += 12
        self.data = self.data[n:]
        #print self.rodForces

    def OEF_CVisc(self):  # 24-CVISC
        dt = self.nonlinear_factor
        (format1, extract) = self.getOEF_FormatStart()
        format1 += 'ff'
        format1 = bytes(format1)
        ntotal = 12  #   # 3*4
        n = 0
        s = Struct(format1)
        nelements = len(self.data) // 12
        for i in xrange(nelements):
            edata = self.data[n:n+12]

            out = s.unpack(edata)
            if self.make_op2_debug:
                self.op2_debug.write('OEF_CVisc - %s\n' % (str(out)))
            (eid_device, axial, torque) = out
            eid = extract(eid_device, dt)
            #print "eType=%s" %(eType)

            dataIn = [eid, axial, torque]
            #print "%s" %(self.get_element_type(self.element_type)),dataIn
            self.obj.add(dt, dataIn)
            n += ntotal
        self.data = self.data[n:]
        #print self.viscForces

    def OEF_Beam(self):  # 2-CBEAM   # TODO is this correct???
        dt = self.nonlinear_factor
        (format1, extract) = self.getOEF_FormatStart()
        #print self.code_information()
        formatAll = 'i8f'
        format1 = bytes(format1)
        formatAll = bytes(formatAll)

        n = 0
        s1 = Struct(format1)
        s2 = Struct(formatAll)
        nelements = len(self.data) // 400  # 1+(10-1)*11=100 ->100*4 = 400
        for i in xrange(nelements):
            #print "eType=%s" %(eType)

            edata = self.data[n:n+4]
            eid_device, = s1.unpack(edata)
            eid = extract(eid_device, dt)
            n += 4

            for i in xrange(11):
                edata = self.data[n:n+36]
                n += 36

                out = s2.unpack(edata)
                if self.make_op2_debug:
                    self.op2_debug.write('OEF_Beam - %s\n' % (str(out)))
                (nid, sd, bm1, bm2, ts1, ts2, af, ttrq, wtrq) = out
                #print "eidTemp = ",eidTemp
                #print "nid = ",nid
                #print "sd = ",sd

                dataIn = [eid, nid, sd, bm1, bm2, ts1, ts2, af, ttrq, wtrq]
                #print "%s        " %(self.get_element_type(self.element_type)),dataIn
                #eid = self.obj.add_new_eid(out)
                if i == 0:  # isNewElement:
                    self.obj.add_new_element(dt, dataIn)
                    #print
                elif sd > 0.:
                    self.obj.add(dt, dataIn)
                #print
            #else: pass
            #print "len(data) = ",len(self.data)
        self.data = self.data[n:]
        #print self.beamForces

    def OEF_Shear(self):  # 4-CSHEAR
        dt = self.nonlinear_factor
        (format1, extract) = self.getOEF_FormatStart()
        format1 += '16f'
        format1 = bytes(format1)
        ntotal = 68  # 17*4
        s = Struct(format1)
        nelements = len(self.data) // 68
        for i in xrange(nelements):
            edata = self.data[n:n+68]

            out = s.unpack(edata)
            if self.make_op2_debug:
                self.op2_debug.write('OEF_Shear - %s\n' % (str(out)))
            (eid_device,
             f41, f21, f12, f32, f23, f43, f34, f14, kf1,
             s12, kf2, s23, kf3, s34, kf4, s41) = out
            eid = extract(eid_device, dt)
            #print "eType=%s" %(eType)

            dataIn = [eid,
                      f41, f21, f12, f32, f23, f43, f34,
                      f14, kf1, s12, kf2, s23, kf3, s34, kf4, s41]
            #print "%s" %(self.get_element_type(self.element_type)),dataIn
            #eid = self.obj.add_new_eid(out)
            self.obj.add(dt, dataIn)
            n += ntotal
        self.data = self.data[n:]
        #print self.shearForces

    def OEF_Spring(self):  # 11-CELAS1, 12-CELAS2, 13-CELAS3, 14-CELAS4
        dt = self.nonlinear_factor
        (format1, extract) = self.getOEF_FormatStart()
        format1 += 'f'
        format1 = bytes(format1)

        n = 0
        s = Struct(format1)
        nelements = len(self.data) // 8 # 2*4
        for i in xrange(nelements):
            edata = self.data[n:n+8]
            n += 8

            out = s.unpack(edata)
            if self.make_op2_debug:
                self.op2_debug.write('OEF_Spring - %s\n' % (str(out)))
            (eid_device, force) = out
            eid = extract(eid_device, dt)
            #print "eType=%s" %(eType)

            dataIn = [eid, force]
            #print "%s" %(self.get_element_type(self.element_type)),dataIn
            self.obj.add(dt, dataIn)
        self.data = self.data[n:]
        #print self.springForces

    def OEF_CBar(self):  # 34-CBAR
        dt = self.nonlinear_factor
        (format1, extract) = self.getOEF_FormatStart()
        format1 += '8f'
        format1 = bytes(format1)

        ntotal = 36  # 9*4
        n = 0
        s = Struct(format1)
        nelements = len(self.data) // ntotal
        for i in xrange(nelements):
            edata = self.data[n:n+36]

            out = s.unpack(edata)
            if self.make_op2_debug:
                self.op2_debug.write('OEF_CBar - %s\n' % (str(out)))
            (eid_device, bm1a, bm2a, bm1b, bm2b, ts1, ts2, af, trq) = out
            eid = extract(eid_device, dt)

            dataIn = [eid, bm1a, bm2a, bm1b, bm2b, ts1, ts2, af, trq]
            #print "%s" %(self.get_element_type(self.element_type)),dataIn
            #eid = self.obj.add_new_eid(out)
            self.obj.add(dt, dataIn)
            n += ntotal
        self.data = self.data[n:]
        #print self.barForces

    def OEF_CBar100(self):  # 100-CBAR
        dt = self.nonlinear_factor
        (format1, extract) = self.getOEF_FormatStart()
        format1 += '7f'
        format1 = bytes(format1)

        n = 0
        s = Struct(format1)
        ntotal = 32  # 8*4
        nelements = len(self.data) // ntotal
        for i in xrange(nelements):
            edata = self.data[n:n+32]

            out = s.unpack(edata)
            if self.make_op2_debug:
                self.op2_debug.write('OEF_CBar100 - %s\n' % (str(out)))
            (eid_device, sd, bm1, bm2, ts1, ts2, af, trq) = out
            eid = extract(eid_device, dt)
            #print "eType=%s" %(eType)

            data_in = [eid, sd, bm1, bm2, ts1, ts2, af, trq]
            #print "%s" %(self.get_element_type(self.element_type)), data_in
            self.obj.add(dt, data_in)
            n += 32
        self.data = self.data[n:]
        #print self.bar100Forces

    def OEF_Plate(self):  # 33-CQUAD4,74-CTRIA3
        dt = self.nonlinear_factor
        (format1, extract) = self.getOEF_FormatStart()
        format1 += '8f'
        format1 = bytes(format1)

        n = 0
        ntotal = 36 # 9*4
        s = Struct(format1)
        nelements = len(self.data) // ntotal
        for i in xrange(nelements):
            edata = self.data[n:n+36]

            out = s.unpack(edata)
            if self.make_op2_debug:
                self.op2_debug.write('OEF_Plate-%s - %s\n' % (self.element_type, str(out)))
            (eid_device, mx, my, mxy, bmx, bmy, bmxy, tx, ty) = out
            eid = extract(eid_device, dt)
            #print "eType=%s" %(eType)

            #print "%s" %(self.get_element_type(self.element_type)),dataIn
            self.obj.add(dt, eid, mx, my, mxy, bmx, bmy, bmxy, tx, ty)
            n += ntotal
        self.data = self.data[n:]
        #print self.plateForces

    def OEF_Plate2(self):  # 64-CQUAD8,70-CTRIAR,75-CTRIA6,82-CQUAD8,144-CQUAD4-bilinear
        dt = self.nonlinear_factor
        (format1, extract) = self.getOEF_FormatStart()
        format1 += '4s'

        if self.element_type in [70, 75]:  # CTRIAR,CTRIA6
            nnodes = 3
        elif self.element_type in [64, 82, 144]:  # CQUAD8,CQUADR,CQUAD4-bilinear
            nnodes = 4
        else:
            raise NotImplementedError(self.code_information())

        allFormat = 'i8f'
        format1 = bytes(format1)
        allFormat = bytes(allFormat)

        n = 0
        ntotal = 44 + nnodes * 36
        s1 = Struct(format1 + allFormat)
        s2 = Struct(allFormat)
        nelements = len(self.data) // ntotal
        for i in xrange(nelements):
            edata = self.data[n:n+44]

            out = s1.unpack(edata)
            if self.make_op2_debug:
                self.op2_debug.write('OEF_Plate2-%s - %s\n' % (self.element_type, str(out)))
            (eid_device, term, nid, mx, my, mxy, bmx, bmy, bmxy, tx, ty) = out
            #term= 'CEN\'
            #print "eType=%s" %(eType)
            eid = extract(eid_device, dt)
            #print "%s" %(self.get_element_type(self.element_type)),dataIn
            self.obj.add_new_element(eid, dt, term, nid, mx, my, mxy, bmx, bmy, bmxy, tx, ty)
            n += 44
            for i in xrange(nnodes):
                edata = self.data[n:n+36]
                out = s2.unpack(edata)
                if self.make_op2_debug:
                    self.op2_debug.write('%s\n' % (str(out)))
                #(nid,mx,my,mxy,bmx,bmy,bmxy,tx,ty) = out
                #dataIn = [nid,mx,my,mxy,bmx,bmy,bmxy,tx,ty]
                #print "***%s    " %(self.get_element_type(self.element_type)),dataIn
                self.obj.add(eid, dt, *out)
                n += 36
        self.data = self.data[n:]
        #print self.plateForces2

    def OEF_CompositePlate_95_96(self):
        pass

    def OEF_ConeAx(self):  # 35-CCONEAX
        dt = self.nonlinear_factor
        (format1, extract) = self.getOEF_FormatStart()
        format1 += '6f'
        format1 = bytes(format1)

        n = 0
        ntotal = 28  # 7*4
        nelements = len(self.data) // ntotal
        s = Struct(format1)
        for i in xrange(nelements):
            edata = self.data[n:n+ntotal]
            out = s.unpack(edata)
            if self.make_op2_debug:
                self.op2_debug.write('OEF_CONEAX-35 - %s\n' % (str(out)))
            (eid_device, hopa, bmu, bmv, tm, su, sv) = out
            eid = extract(eid_device, dt)
            #print "eType=%s" % (eType)

            data_in = [eid, hopa, bmu, bmv, tm, su, sv]
            #print "%s" % (self.get_element_type(self.element_type)), data_in
            #eid = self.obj.add_new_eid(out)
            self.obj.add(dt, data_in)
            n += ntotal
        self.data = self.data[n:]
        #print self.shearForces

    def OEF_CGap(self):  # 38-CGAP
        dt = self.nonlinear_factor
        (format1, extract) = self.getOEF_FormatStart()
        format1 += '8f'
        format1 = bytes(format1)

        n = 0
        ntotal = 36 # 9*4
        nelements = len(self.data) // ntotal
        s = Struct(format1)
        for i in xrange(nelements):
            edata = self.data[n:n+36]

            out = s.unpack(edata)
            if self.make_op2_debug:
                self.op2_debug.write('OEF_CGAP-38 - %s\n' % (str(out)))
            (eid_device, fx, sfy, sfz, u, v, w, sv, sw) = out
            eid = extract(eid_device, dt)
            #print "eType=%s" %(eType)

            data_in = [eid, fx, sfy, sfz, u, v, w, sv, sw]
            #print "%s" %(self.get_element_type(self.element_type)),data_in
            #eid = self.obj.add_new_eid(out)
            self.obj.add(dt, data_in)
            n += ntotal
        self.data = self.data[n:]
        #print self.plateForces

    def OEF_Bend(self):  # 69-CBEND
        dt = self.nonlinear_factor
        (format1, extract) = self.getOEF_FormatStart()
        format1 += 'i13f'
        format1 = bytes(format1)

        n = 0
        ntotal = 60  # 15*4
        s = Struct(format1)
        nelements = len(self.data) // ntotal
        for i in xrange(nelements):
            edata = self.data[n:n+ntotal]

            out = s.unpack(edata)
            if self.make_op2_debug:
                self.op2_debug.write('OEF_BEND-69 - %s\n' % (str(out)))
            (eid_device, nidA, bm1A, bm2A, ts1A, ts2A, afA, trqA,
                         nidB, bm1B, bm2B, ts1B, ts2B, afB, trqB) = out
            eid = extract(eid_device, dt)
            #print "eType=%s" % (eType)

            data_in = [eid,
                       nidA, bm1A, bm2A, ts1A, ts2A, afA, trqA,
                       nidB, bm1B, bm2B, ts1B, ts2B, afB, trqB]
            #print "%s" %(self.get_element_type(self.element_type)), data_in
            self.obj.add(dt, data_in)
            n += ntotal
        self.data = self.data[n:]
        #print self.bendForces

    def OEF_PentaPressure(self):  # 77-CPENTA_PR,78-CTETRA_PR
        dt = self.nonlinear_factor
        (format1, extract) = self.getOEF_FormatStart()
        format1 += '8s7f'
        format1 = bytes(format1)

        n = 0
        s = Struct(format1)
        ntotal = 40
        nelements = len(self.data) // ntotal
        for i in xrange(nelements):
            eData = self.data[n:n+40]
            n += 40
            #print "len(data) = ",len(eData)

            out = s.unpack(eData)
            if self.make_op2_debug:
                self.op2_debug.write('OEF_PentaPressure-%s %s\n' % (self.element_type, str(out)))
            (eid_device, eName, ax, ay, az, vx, vy, vz, pressure) = out
            eid = extract(eid_device, dt)
            #print "eType=%s" %(eType)

            dataIn = [eid, eName, ax, ay, az, vx, vy, vz, pressure]
            #print "%s" %(self.get_element_type(self.element_type)),dataIn
            self.obj.add(dt, dataIn)
        self.data = self.data[n:]
        #print self.pentaPressureForces

    def OEF_CBush(self):  # 102-CBUSH
        dt = self.nonlinear_factor
        (format1, extract) = self.getOEF_FormatStart()
        format1 += '6f'
        format1 = bytes(format1)

        n = 0
        ntotal = 28 # 7*4
        s = Struct(format1)
        nelements = len(self.data) // ntotal
        for i in xrange(nelements):
            edata = self.data[n:n+28]

            out = s.unpack(edata)
            if self.make_op2_debug:
                self.op2_debug.write('OEF_CBUSH-102 - %s\n' % (str(out)))
            (eid_device, fx, fy, fz, mx, my, mz) = out
            eid = extract(eid_device, dt)
            #print "eType=%s" %(eType)

            data_in = [eid, fx, fy, fz, mx, my, mz]
            #print "%s" %(self.get_element_type(self.element_type)), data_in
            #eid = self.obj.add_new_eid(out)
            self.obj.add(dt, data_in)
            n += ntotal
        self.data = self.data[n:]
        #print self.bushForces

    def OEF_Force_VU(self):  # 191-VUBEAM
        dt = self.nonlinear_factor

        (format1, extract) = self.getOEF_FormatStart()
        format1 += 'ii4s'

        if self.element_type in [191]:
            nNodes = 2
        else:
            raise NotImplementedError(self.code_information())

        formatAll = 'i7f'
        format1 = bytes(format1)
        formatAll = bytes(formatAll)

        n = 0
        s1 = Struct(format1)
        s2 = Struct(formatAll)
        ntotal = 16 + 32 * nNodes
        nelements = len(self.data) // ntotal
        for i in xrange(nelements):
            eData = self.data[n:n+16]  # 8*4
            n += 16

            out = s1.unpack(format1, eData)
            if self.make_op2_debug:
                self.op2_debug.write('OEF_Force_VU-191 - %s\n' % (str(out)))
            (eid_device, parent, coord, icord) = out

            eid = extract(eid_device, dt)
            dataIn = [eid, parent, coord, icord]

            forces = []
            for i in xrange(nNodes):
                eData = self.data[n:n+32]  # 8*4
                n += 32
                #print "i=%s len(data)=%s" %(i,len(eData))
                out = s2.unpack(eData)
                if self.make_op2_debug:
                    self.op2_debug.write('%s\n' % str(out))
                forces.append(out)
            dataIn.append(forces)
            #eType = a+b+c+d+e+f+g+h
            #print "eType=%s" %(eType)

            #dataIn = [vugrid,posit,forceX,shearY,shearZ,torsion,bendY,bendZ]
            #print "force %s" %(self.get_element_type(self.element_type)),dataIn
            self.obj.add(nNodes, dt, dataIn)
        self.data = self.data[n:]

    def OEF_Force_VUTRIA(self):  # 189-VUQUAD,190-VUTRIA
        dt = self.nonlinear_factor
        (format1, extract) = self.getOEF_FormatStart()
        format1 += 'ii4sii'

        if self.element_type in [189]:  # VUQUAD
            eType = 'VUQUAD'
            nNodes = 4
        elif self.element_type in [190]:  # VUTRIA
            eType = 'VUTRIA'
            nNodes = 3
        else:
            raise NotImplementedError(self.code_information())

        formatAll = 'i3f3i5fi'
        format1 = bytes(format1)
        formatAll = bytes(formatAll)
        ntotal = 24 + 52 * nNodes
        n = 0
        s1 = Struct(format1)
        s2 = Struct(formatAll)
        nelements = len(self.data) // ntotal
        for i in xrange(nelements):
            eData = self.data[n:n+24]  # 6*4
            n += 24

            out = s1.unpack(eData)
            if self.make_op2_debug:
                self.op2_debug.write('OEF_Force_%s-%s - %s\n' % (eType, self.element_type, str(out)))
            (eid_device, parent, coord, icord, theta, _) = out

            eid = extract(eid_device, dt)
            dataIn = [eid, parent, coord, icord, theta]

            forces = []
            for i in xrange(nNodes):
                eData = self.data[n:n+52]  # 13*4
                n += 52
                #print "i=%s len(data)=%s" %(i,len(eData))
                out = s2.unpack(eData)
                if self.make_op2_debug:
                    self.op2_debug.write('%s\n' % (str(out)))
                (vugrid, mfx, mfy, mfxy, a, b, c, bmx, bmy,
                    bmxy, syz, szx, d) = out
                out2 = (vugrid, mfx, mfy, mfxy, bmx, bmy, bmxy, syz, szx)
                forces.append(out2)
            dataIn.append(forces)
            #print "eType=%s" %(eType)

            #dataIn = [vugrid,mfx,mfy,mfxy,a,b,c,bmx,bmy,bmxy,syz,szx,d]
            #print "force %s" %(self.get_element_type(self.element_type)),dataIn
            self.obj.add(nNodes, dt, dataIn)

        self.data = self.data[n:]
        #print(self.force_VU_2D)
