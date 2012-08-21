from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
import sys
from struct import unpack


class RealForces(object):

    def OEF_aCode(self):
        if self.analysisCode == 1:   # statics
            format1 = 'i'  # loadID
            #self.addDataParameter(data,'loadID','i',5,False)   ## load set ID number
        elif self.analysisCode == 2:  # normal modes/buckling (real eigenvalues)
            format1 = 'i'  # mode
            #self.addDataParameter(data,'mode','i',5)   ## mode number
            #self.addDataParameter(data,'eign','f',6,False)   ## eigenvalue
        elif self.analysisCode == 3:  # differential stiffness 0
            format1 = 'i'  # loadID
            #self.addDataParameter(data,'loadID','i',5)   ## load set ID number
        elif self.analysisCode == 4:  # differential stiffness 1
            format1 = 'i'  # loadID
            ## load set ID number
            self.addDataParameter(data, 'loadID', 'i', 5)
        elif self.analysisCode == 5:   # frequency
            format1 = 'f'  # freq
            #self.addDataParameter(data,'freq','f',5)   ## frequency

        elif self.analysisCode == 6:  # transient
            format1 = 'f'  # time
            #self.addDataParameter(data,'time','f',5)   ## time step
            #print "time(5)=%s" %(self.time)
        elif self.analysisCode == 7:  # pre-buckling
            format1 = 'i'  # loadID
            #self.addDataParameter(data,'loadID','i',5)   ## load set ID number
            #print "loadID(5)=%s" %(self.loadID)
        elif self.analysisCode == 8:  # post-buckling
            format1 = 'i'  # loadID
            #self.addDataParameter(data,'loadID','i',5)       ## load set ID number
            #self.addDataParameter(data,'eigr','f',6,False)   ## real eigenvalue
            #print "loadID(5)=%s  eigr(6)=%s" %(self.loadID,self.eigr)
        elif self.analysisCode == 9:  # complex eigenvalues
            format1 = 'i'  # mode
            #self.addDataParameter(data,'mode','i',5)         ## mode number
            #self.addDataParameter(data,'eigr','f',6,False)   ## real eigenvalue
            #self.addDataParameter(data,'eigi','f',7,False)   ## imaginary eigenvalue
            #print "mode(5)=%s  eigr(6)=%s  eigi(7)=%s" %(self.mode,self.eigr,self.eigi)
        elif self.analysisCode == 10:  # nonlinear statics
            format1 = 'f'  # loadStep
            ## load step
            self.addDataParameter(data, 'loadStep', 'f', 5)
            #print "loadStep(5) = %s" %(self.loadStep)
        elif self.analysisCode == 11:  # geometric nonlinear statics
            format1 = 'i'  # loadID
            #self.addDataParameter(data,'loadID','i',5)   ## load set ID number
            #print "loadID(5)=%s" %(self.loadID)
        else:
            raise RuntimeError('invalid analysisCode...analysisCode=%s' % (str(self.analysisCode) + '\n' + self.codeInformation()))
        return format1

    def getOEF_FormatStart(self):
        """
        Returns an i or an f depending on if it's SORT2 or not.
        Also returns an extraction function that is called on the first argument
        """
        isSort1 = self.isSort1()
        if isSort1:
            #print "SORT1 - %s" %(self.ElementType(self.elementType))
            format1 = 'i'  # SORT1
            extract = self.extractSort1
            #dt = self.nonlinearFactor
        else:
            #print "SORT2 - %s" %(self.ElementType(self.elementType))
            format1 = 'f'  # SORT2
            extract = self.extractSort2
            #eid = self.nonlinearFactor
        return (format1, extract)

    def OEF_Rod(self):  # 1-CROD, 3-CTUBE, 10-CONROD
        dt = self.nonlinearFactor
        (format1, extract) = self.getOEF_FormatStart()
        format1 += 'ff'
        format1 = bytes(format1)

        while len(self.data) >= 12:  # 3*4
            eData = self.data[0:12]
            self.data = self.data[12:]
            #print "len(data) = ",len(eData)

            out = unpack(format1, eData)
            (eid, axial, torque) = out
            eid2 = extract(eid, dt)
            #print "eType=%s" %(eType)

            dataIn = [eid2, axial, torque]
            #print "%s" %(self.ElementType(self.elementType)),dataIn
            #eid = self.obj.addNewEid(out)
            self.obj.add(dt, dataIn)
            #print "len(data) = ",len(self.data)
        #print self.rodForces

    def OEF_CVisc(self):  # 24-CVISC
        dt = self.nonlinearFactor
        (format1, extract) = self.getOEF_FormatStart()
        format1 += 'ff'
        format1 = bytes(format1)

        while len(self.data) >= 12:  # 3*4
            eData = self.data[0:12]
            self.data = self.data[12:]
            #print "len(data) = ",len(eData)

            out = unpack(format1, eData)
            (eid, axial, torque) = out
            eid2 = extract(eid, dt)
            #print "eType=%s" %(eType)

            dataIn = [eid2, axial, torque]
            #print "%s" %(self.ElementType(self.elementType)),dataIn
            #eid = self.obj.addNewEid(out)
            self.obj.add(dt, dataIn)
            #print "len(data) = ",len(self.data)
        #print self.viscForces

    def OEF_Beam(self):  # 2-CBEAM   ## @todo is this correct???
        dt = self.nonlinearFactor
        (format1, extract) = self.getOEF_FormatStart()
        #print self.codeInformation()
        formatAll = 'iffffffff'
        format1 = bytes(format1)
        formatAll = bytes(formatAll)

        while len(self.data) >= 400:  # 1+(10-1)*11=100 ->100*4 = 400
            #print "eType=%s" %(eType)

            eData = self.data[0:4]
            self.data = self.data[4:]
            eid, = unpack(format1, eData)
            eid2 = extract(eid, dt)

            for i in xrange(11):
                eData = self.data[0:36]
                self.data = self.data[36:]
                #print "len(data) = ",len(eData)

                out = unpack(formatAll, eData)
                (nid, sd, bm1, bm2, ts1, ts2, af, ttrq, wtrq) = out
                #print "eidTemp = ",eidTemp
                #print "nid = ",nid
                #print "sd = ",sd

                dataIn = [eid2, nid, sd, bm1, bm2, ts1, ts2, af, ttrq, wtrq]
                #print "%s        " %(self.ElementType(self.elementType)),dataIn
                #eid = self.obj.addNewEid(out)
                if i == 0:  # isNewElement:
                    self.obj.addNewElement(dt, dataIn)
                    #print
                elif sd > 0.:
                    self.obj.add(dt, dataIn)
                #print
            #else: pass
            #print "len(data) = ",len(self.data)
        #print self.beamForces

    def OEF_Shear(self):  # 4-CSHEAR
        dt = self.nonlinearFactor
        (format1, extract) = self.getOEF_FormatStart()
        format1 += 'ffffffffffffffff'
        format1 = bytes(format1)

        while len(self.data) >= 68:  # 17*4
            eData = self.data[0:68]
            self.data = self.data[68:]
            #print "len(data) = ",len(eData)

            out = unpack(format1, eData)
            (eid, f41, f21, f12, f32, f23, f43, f34, f14, kf1,
                s12, kf2, s23, kf3, s34, kf4, s41) = out
            eid2 = extract(eid, dt)
            #print "eType=%s" %(eType)

            dataIn = [eid2, f41, f21, f12, f32, f23, f43, f34,
                      f14, kf1, s12, kf2, s23, kf3, s34, kf4, s41]
            #print "%s" %(self.ElementType(self.elementType)),dataIn
            #eid = self.obj.addNewEid(out)
            self.obj.add(dt, dataIn)
            #print "len(data) = ",len(self.data)
        #print self.shearForces

    def OEF_Spring(self):  # 11-CELAS1, 12-CELAS2, 13-CELAS3, 14-CELAS4
        dt = self.nonlinearFactor
        (format1, extract) = self.getOEF_FormatStart()
        format1 += 'f'
        format1 = bytes(format1)

        while len(self.data) >= 8:  # 2*4
            eData = self.data[0:8]
            self.data = self.data[8:]
            #print "len(data) = ",len(eData)

            out = unpack(format1, eData)
            (eid, force) = out
            eid2 = extract(eid, dt)
            #print "eType=%s" %(eType)

            dataIn = [eid2, force]
            #print "%s" %(self.ElementType(self.elementType)),dataIn
            #eid = self.obj.addNewEid(out)
            self.obj.add(dt, dataIn)
            #print "len(data) = ",len(self.data)
        #print self.springForces

    def OEF_CBar(self):  # 34-CBAR
        dt = self.nonlinearFactor
        (format1, extract) = self.getOEF_FormatStart()
        format1 += 'ffffffff'
        format1 = bytes(format1)

        while len(self.data) >= 36:  # 9*4
            eData = self.data[0:36]
            self.data = self.data[36:]
            #print "len(data) = ",len(eData)

            out = unpack(format1, eData)
            (eid, bm1a, bm2a, bm1b, bm2b, ts1, ts2, af, trq) = out
            eid2 = extract(eid, dt)
            #print "eType=%s" %(eType)

            dataIn = [eid2, bm1a, bm2a, bm1b, bm2b, ts1, ts2, af, trq]
            #print "%s" %(self.ElementType(self.elementType)),dataIn
            #eid = self.obj.addNewEid(out)
            self.obj.add(dt, dataIn)
            #print "len(data) = ",len(self.data)
        #print self.barForces

    def OEF_CBar100(self):  # 100-CBAR
        dt = self.nonlinearFactor
        (format1, extract) = self.getOEF_FormatStart()
        format1 += 'fffffff'
        format1 = bytes(format1)

        while len(self.data) >= 36:  # 9*4
            eData = self.data[0:32]
            self.data = self.data[32:]
            #print "len(data) = ",len(eData)

            out = unpack(format1, eData)
            (eid, sd, bm1, bm2, ts1, ts2, af, trq) = out
            eid2 = extract(eid, dt)
            #print "eType=%s" %(eType)

            dataIn = [eid2, sd, bm1, bm2, ts1, ts2, af, trq]
            #print "%s" %(self.ElementType(self.elementType)),dataIn
            #eid = self.obj.addNewEid(out)
            self.obj.add(dt, dataIn)
            #print "len(data) = ",len(self.data)
        #print self.bar100Forces

    def OEF_Plate(self):  # 33-CQUAD4,74-CTRIA3
        dt = self.nonlinearFactor
        (format1, extract) = self.getOEF_FormatStart()
        format1 += 'ffffffff'
        format1 = bytes(format1)

        while len(self.data) >= 36:  # 9*4
            eData = self.data[0:36]
            self.data = self.data[36:]
            #print "len(data) = ",len(eData)

            out = unpack(format1, eData)
            (eid, mx, my, mxy, bmx, bmy, bmxy, tx, ty) = out
            eid2 = extract(eid, dt)
            #print "eType=%s" %(eType)

            dataIn = [eid2, mx, my, mxy, bmx, bmy, bmxy, tx, ty]
            #print "%s" %(self.ElementType(self.elementType)),dataIn
            #eid = self.obj.addNewEid(out)
            self.obj.add(dt, dataIn)
            #print "len(data) = ",len(self.data)
        #print self.plateForces

    def OEF_Plate2(self):  # 64-CQUAD8,70-CTRIAR,75-CTRIA6,82-CQUAD8,144-CQUAD4-bilinear
        dt = self.nonlinearFactor
        (format1, extract) = self.getOEF_FormatStart()
        format1 += '4s'

        if self.elementType in [70, 75]:  # CTRIAR,CTRIA6
            nNodes = 4
        elif self.elementType in [64, 82, 144]:  # CQUAD8,CQUADR,CQUAD4-bilinear
            nNodes = 5
        else:
            raise NotImplementedError(self.codeInformation())

        allFormat = 'fffffffff'
        format1 = bytes(format1)
        allFormat = bytes(allFormat)

        nTotal = 44 + nNodes * 36
        while len(self.data) >= nTotal:
            eData = self.data[0:44]
            self.data = self.data[44:]
            #print self.printBlock(eData)
            #print "len(data) = ",len(eData)

            out = unpack(format1 + allFormat, eData)
            (eid, term, nid, mx, my, mxy, bmx, bmy, bmxy, tx, ty) = out
            #term= 'CEN\'
            #print "eType=%s" %(eType)

            eid2 = extract(eid, dt)

            dataIn = [term, nid, mx, my, mxy, bmx, bmy, bmxy, tx, ty]
            #print "%s" %(self.ElementType(self.elementType)),dataIn
            self.obj.addNewElement(eid2, dt, dataIn)

            for i in xrange(nNodes - 1):
                eData = self.data[0:36]
                self.data = self.data[36:]
                dataIn = unpack(allFormat, eData)
                #(nid,mx,my,mxy,bmx,bmy,bmxy,tx,ty) = out
                #dataIn = [nid,mx,my,mxy,bmx,bmy,bmxy,tx,ty]
                #print "***%s    " %(self.ElementType(self.elementType)),dataIn

                self.obj.add(eid2, dt, dataIn)
                #print "len(data) = ",len(self.data)

        #print self.plateForces2

    def OEF_ConeAx(self):  # 35-CCONEAX
        dt = self.nonlinearFactor
        (format1, extract) = self.getOEF_FormatStart()
        format1 += '6f'
        format1 = bytes(format1)

        while len(self.data) >= 28:  # 7*4
            eData = self.data[0:28]
            self.data = self.data[28:]
            #print "len(data) = ",len(eData)

            out = unpack(format1, eData)
            (eid, hopa, bmu, bmv, tm, su, sv) = out
            eid2 = extract(eid, dt)
            #print "eType=%s" %(eType)

            dataIn = [eid2, hopa, bmu, bmv, tm, su, sv]
            #print "%s" %(self.ElementType(self.elementType)),dataIn
            #eid = self.obj.addNewEid(out)
            self.obj.add(dt, dataIn)
            #print "len(data) = ",len(self.data)
        #print self.shearForces

    def OEF_CGap(self):  # 38-CGAP
        dt = self.nonlinearFactor
        (format1, extract) = self.getOEF_FormatStart()
        format1 += 'ffffffff'
        format1 = bytes(format1)

        while len(self.data) >= 36:  # 9*4
            eData = self.data[0:36]
            self.data = self.data[36:]
            #print "len(data) = ",len(eData)

            out = unpack(format1, eData)
            (eid, fx, sfy, sfz, u, v, w, sv, sw) = out
            eid2 = extract(eid, dt)
            #print "eType=%s" %(eType)

            dataIn = [eid2, fx, sfy, sfz, u, v, w, sv, sw]
            #print "%s" %(self.ElementType(self.elementType)),dataIn
            #eid = self.obj.addNewEid(out)
            self.obj.add(dt, dataIn)
            #print "len(data) = ",len(self.data)
        #print self.plateForces

    def OEF_Bend(self):  # 69-CBEND
        dt = self.nonlinearFactor
        (format1, extract) = self.getOEF_FormatStart()
        format1 += 'i13f'
        format1 = bytes(format1)

        while len(self.data) >= 60:  # 15*4
            eData = self.data[0:60]
            self.data = self.data[60:]
            #print "len(data) = ",len(eData)

            out = unpack(format1, eData)
            (eid, nidA, bm1A, bm2A, ts1A, ts2A, afA, trqA,
             nidB, bm1B, bm2B, ts1B, ts2B, afB, trqB) = out
            eid2 = extract(eid, dt)
            #print "eType=%s" %(eType)

            dataIn = [eid2, nidA, bm1A, bm2A, ts1A, ts2A, afA, trqA,
                      nidB, bm1B, bm2B, ts1B, ts2B, afB, trqB]
            #print "%s" %(self.ElementType(self.elementType)),dataIn
            #eid = self.obj.addNewEid(out)
            self.obj.add(dt, dataIn)
            #print "len(data) = ",len(self.data)
        #print self.bendForces

    def OEF_PentaPressure(self):  # 77-CPENTA_PR,78-CTETRA_PR
        dt = self.nonlinearFactor
        (format1, extract) = self.getOEF_FormatStart()
        format1 += '8s7f'
        format1 = bytes(format1)

        while len(self.data) >= 40:  # 10*4
            eData = self.data[0:40]
            self.data = self.data[40:]
            #print "len(data) = ",len(eData)

            out = unpack(format1, eData)
            (eid, eName, ax, ay, az, vx, vy, vz, pressure) = out
            eid2 = extract(eid, dt)
            #print "eType=%s" %(eType)

            dataIn = [eid2, eName, ax, ay, az, vx, vy, vz, pressure]
            #print "%s" %(self.ElementType(self.elementType)),dataIn
            #eid = self.obj.addNewEid(out)
            self.obj.add(dt, dataIn)
            #print "len(data) = ",len(self.data)
        #print self.pentaPressureForces

    def OEF_CBush(self):  # 102-CBUSH
        dt = self.nonlinearFactor
        (format1, extract) = self.getOEF_FormatStart()
        format1 += '6f'
        format1 = bytes(format1)

        while len(self.data) >= 28:  # 7*4
            eData = self.data[0:28]
            self.data = self.data[28:]
            #print "len(data) = ",len(eData)

            out = unpack(format1, eData)
            (eid, fx, fy, fz, mx, my, mz) = out
            eid2 = extract(eid, dt)
            #print "eType=%s" %(eType)

            dataIn = [eid2, fx, fy, fz, mx, my, mz]
            #print "%s" %(self.ElementType(self.elementType)),dataIn
            #eid = self.obj.addNewEid(out)
            self.obj.add(dt, dataIn)
            #print "len(data) = ",len(self.data)
        #print self.bushForces

    def OEF_Force_VU(self):  # 191-VUBEAM
        dt = self.nonlinearFactor

        (format1, extract) = self.getOEF_FormatStart()
        format1 += 'ii4s'

        if self.elementType in [191]:
            nNodes = 2
        else:
            raise NotImplementedError(self.codeInformation())

        formatAll = 'i7f'
        format1 = bytes(format1)
        formatAll = bytes(formatAll)

        n = 16 + 32 * nNodes
        while len(self.data) >= n:
            eData = self.data[0:16]  # 8*4
            self.data = self.data[16:]
            #print "len(data) = ",len(eData)

            out = unpack(format1, eData)
            (eid, parent, coord, icord) = out

            eid2 = extract(eid, dt)
            dataIn = [eid2, parent, coord, icord]

            forces = []
            for i in xrange(nNodes):
                eData = self.data[0:32]  # 8*4
                self.data = self.data[32:]
                #print "i=%s len(data)=%s" %(i,len(eData))
                out = unpack(formatAll, eData)
                forces.append(out)
            dataIn.append(forces)
            #eType = a+b+c+d+e+f+g+h
            #print "eType=%s" %(eType)

            #dataIn = [vugrid,posit,forceX,shearY,shearZ,torsion,bendY,bendZ]
            #print "force %s" %(self.ElementType(self.elementType)),dataIn
            #eid = self.obj.addNewEid(out)
            self.obj.add(nNodes, dt, dataIn)
            #print "len(data) = ",len(self.data)
        ###
        if self.makeOp2Debug:
            print("done with OEF_Force_VU")
        #print(self.force_VU)

    def OEF_Force_VUTRIA(self):  # 189-VUQUAD,190-VUTRIA
        dt = self.nonlinearFactor
        (format1, extract) = self.getOEF_FormatStart()
        format1 += 'ii4sii'

        if self.elementType in [189]:  # VUQUAD
            nNodes = 4
        elif self.elementType in [190]:  # VUTRIA
            nNodes = 3
        else:
            raise NotImplementedError(self.codeInformation())

        formatAll = 'ifffiiifffffi'
        format1 = bytes(format1)
        formatAll = bytes(formatAll)
        n = 24 + 52 * nNodes
        while len(self.data) >= n:
            eData = self.data[0:24]  # 6*4
            self.data = self.data[24:]
            #print "len(data) = ",len(eData)

            out = unpack(format1, eData)
            (eid, parent, coord, icord, theta, _) = out

            eid2 = extract(eid, dt)
            dataIn = [eid2, parent, coord, icord, theta]

            forces = []
            for i in xrange(nNodes):
                eData = self.data[0:52]  # 13*4
                self.data = self.data[52:]
                #print "i=%s len(data)=%s" %(i,len(eData))
                out = unpack(formatAll, eData)
                (vugrid, mfx, mfy, mfxy, a, b, c, bmx, bmy,
                    bmxy, syz, szx, d) = out
                out2 = (vugrid, mfx, mfy, mfxy, bmx, bmy, bmxy, syz, szx)
                forces.append(out2)
            dataIn.append(forces)
            #eType = a+b+c+d+e+f+g+h
            #print "eType=%s" %(eType)

            #dataIn = [vugrid,mfx,mfy,mfxy,a,b,c,bmx,bmy,bmxy,syz,szx,d]
            #print "force %s" %(self.ElementType(self.elementType)),dataIn
            #eid = self.obj.addNewEid(out)
            self.obj.add(nNodes, dt, dataIn)
            #print "len(data) = ",len(self.data)
        ###
        #print(self.force_VU_2D)
