from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
import sys
from struct import unpack

from pyNastran.op2.op2_helper import polarToRealImag

#91  -> PENTANL
#2   -> BEAM
#33  -> TUBE
#92  -> CONRODNL


class ComplexElementsStressStrain(object):

    def OES_Rod1_alt(self):
        """
        genericStressReader - works on CROD_1, CELAS2_12
        stress & strain
        formatCode=1 sortCode=1 (eid,axial,axial,torsion,torsion)
        """
        dt = self.nonlinearFactor
        (format1, extract) = self.getOUG_FormatStart()
        nTotal = 12
        format1 += '4f'
        format1 = bytes(format1)
        isMagnitudePhase = self.isMagnitudePhase()

        n = 0
        nEntries = len(self.data) // nTotal
        for i in xrange(nEntries):
            eData = self.data[n:n + nTotal]
            (eid, axialReal, axialImag, torsionReal,
                torsionImag) = unpack(format1, eData)

            if isMagnitudePhase:
                (axial) = polarToRealImag(axialReal, axialImag)
                (torsion) = polarToRealImag(torsionReal, torsionImag)
            else:
                axial = complex(axialReal, axialImag)
                torsion = complex(torsionReal, torsionImag)

            #print "out = ",out
            eid = extract(eid, dt)
            self.obj.addNewEid(dt, eid, axial, torsion)
            n += nTotal
        self.data = self.data[n:]

    def OES_Elas1_alt(self):
        """
        genericStressReader - works on CROD_1, CELAS2_12
        stress & strain
        formatCode=1 sortCode=1 (eid,axial,axial,torsion,torsion)
        """
        dt = self.nonlinearFactor
        (format1, extract) = self.getOUG_FormatStart()
        nTotal = 12
        format1 += '2f'
        format1 = bytes(format1)
        isMagnitudePhase = self.isMagnitudePhase()

        n = 0
        nEntries = len(self.data) // nTotal
        for i in xrange(nEntries):
            eData = self.data[n:n + nTotal]
            (eid, axialReal, axialImag) = unpack(format1, eData)

            if isMagnitudePhase:
                axial = polarToRealImag(axialReal, axialImag)
            else:
                axial = complex(axialReal, axialImag)

            #print "out = ",out
            eid = extract(eid, dt)
            self.obj.addNewEid(dt, eid, axial)
            n += nTotal
        self.data = self.data[n:]

    def OES_CBAR_34_alt(self):
        dt = self.nonlinearFactor
        #print "len(data) = ",len(self.data)
        assert self.numWide == 19, 'invalid numWide...numWide=%s' % (
            self.numWide)

        (format1, extract) = self.getOUG_FormatStart()
        format1 += '18f'
        format1 = bytes(format1)
        isMagnitudePhase = self.isMagnitudePhase()

        while len(self.data) >= 76:
            #self.printBlock(self.data)
            eData = self.data[0:76]
            self.data = self.data[76:]
            #print "len(data) = ",len(eData)

            (eid, s1ar, s2ar, s3ar, s4ar, axialr,
             s1ai, s2ai, s3ai, s4ai, axiali,
             s1br, s2br, s3br, s4br,
             s1bi, s2bi, s3bi, s4bi) = unpack(format1, eData)

            if isMagnitudePhase:
                s1a = polarToRealImag(s1ar, s1ai)
                s1b = polarToRealImag(s1br, s1bi)
                s2a = polarToRealImag(s2ar, s2ai)
                s2b = polarToRealImag(s2br, s2bi)
                s3a = polarToRealImag(s3ar, s3ai)
                s3b = polarToRealImag(s3br, s3bi)
                s4a = polarToRealImag(s4ar, s4ai)
                s4b = polarToRealImag(s4br, s4bi)
                axial = polarToRealImag(axialr, axiali)

            else:
                s1a = complex(s1ar, s1ai)
                s1b = complex(s1br, s1bi)
                s2a = complex(s2ar, s2ai)
                s2b = complex(s2br, s2bi)
                s3a = complex(s3ar, s3ai)
                s3b = complex(s3br, s3bi)
                s4a = complex(s4ar, s4ai)
                s4b = complex(s4br, s4bi)
                axial = complex(axialr, axiali)

            eid2 = extract(eid, dt)
            self.obj.addNewEid('CBAR', dt, eid2, s1a, s2a, s3a, s4a, axial,
                               s1b, s2b, s3b, s4b)

            #print "eid=%i s1=%i s2=%i s3=%i s4=%i axial=%-5i" %(eid,s1a,s2a,s3a,s4a,axial)
            #print "         s1=%i s2=%i s3=%i s4=%i"          %(s1b,s2b,s3b,s4b)
            #print "len(data) = ",len(self.data)

    def OES_CQUAD4_33_alt(self):
        """
        GRID-ID  DISTANCE,NORMAL-X,NORMAL-Y,SHEAR-XY
        """
        dt = self.nonlinearFactor
        (format1, extract) = self.getOUG_FormatStart()
        format1 += '14f'
        format1 = bytes(format1)
        isMagnitudePhase = self.isMagnitudePhase()

        nNodes = 0  # centroid + 4 corner points
        #self.print_section(20)
        #term = data[0:4] CEN/
        #data = data[4:]
        #print "*****"
        #self.printBlock(self.data)
        #print "self.numWide = ",self.numWide
        #print "len(data) = ",len(self.data)

        assert self.numWide == 15, 'invalid numWide...numWide=%s' % (
            self.numWide)
        while len(self.data) >= 60:  # 2+15*5 = 77 -> 77*4 = 308
            #print self.printBlock(self.data[0:100])
            #(eid,) = unpack(b'i',self.data[0:4])
            #print "abcd=",abcd
            #self.data = self.data[8:]  # 2
            eData = self.data[0:60]  # 4*15=60
            self.data = self.data[60:]
            out = unpack(format1, eData)  # 15
            if self.makeOp2Debug:
                self.op2Debug.write('%s\n' % (str(out)))
            (eid, fd1, sx1r, sx1i, sy1r, sy1i, txy1r, txy1i,
             fd2, sx2r, sx2i, sy2r, sy2i, txy2r, txy2i) = out

            if isMagnitudePhase:
                sx1 = polarToRealImag(sx1r, sx1i)
                sx2 = polarToRealImag(sx2r, sx2i)
                sy1 = polarToRealImag(sy1r, sy1i)
                sy2 = polarToRealImag(sy2r, sy2i)
                txy1 = polarToRealImag(txy1r, txy1i)
                txy2 = polarToRealImag(txy2r, txy2i)
            else:
                sx1 = complex(sx1r, sx1i)
                sx2 = complex(sx2r, sx2i)
                sy1 = complex(sy1r, sy1i)
                sy2 = complex(sy2r, sy2i)
                txy1 = complex(txy1r, txy1i)
                txy2 = complex(txy2r, txy2i)

            eid = extract(eid, dt)

            #print "eid=%i grid=%s fd1=%-3.1f sx1=%s sy1=%s txy1=%s" %(eid,'C',fd1,sx1,sy1,txy1)
            #print   "             fd2=%-3.1f sx2=%s sy2=%s txy2=%s\n"       %(fd2,sx2,sy2,txy2)
            #print "nNodes = ",nNodes
            self.obj.addNewEid('CQUAD4', dt, eid, 'C', fd1, sx1, sy1, txy1)
            self.obj.add(dt, eid, 'C', fd2, sx2, sy2, txy2)

            for nodeID in xrange(nNodes):  # nodes pts
                eData = self.data[0:60]  # 4*15=60
                self.data = self.data[60:]
                out = unpack(b'i14f', eData[0:60])
                if self.makeOp2Debug:
                    self.op2Debug.write('%s\n' % (str(out)))
                (grid, fd1, sx1r, sx1i, sy1r, sy1i, txy1r, txy1i,
                 fd2, sx2r, sx2i, sy2r, sy2i, txy2r, txy2i) = out

                if isMagnitudePhase:
                    sx1 = polarToRealImag(sx1r, sx1i)
                    sx2 = polarToRealImag(sx2r, sx2i)
                    sy1 = polarToRealImag(sy1r, sy1i)
                    sy2 = polarToRealImag(sy2r, sy2i)
                    txy1 = polarToRealImag(txy1r, txy1i)
                    txy2 = polarToRealImag(txy2r, txy2i)
                else:
                    sx1 = complex(sx1r, sx1i)
                    sx2 = complex(sx2r, sx2i)
                    sy1 = complex(sy1r, sy1i)
                    sy2 = complex(sy2r, sy2i)
                    txy1 = complex(txy1r, txy1i)
                    txy2 = complex(txy2r, txy2i)

                #print "eid=%i grid=%i fd1=%i sx1=%i sy1=%i txy1=%i\n" %(eid,grid,fd1,sx1,sy1,txy1)
                #print "               fd2=%i sx2=%i sy2=%i txy2=%i\n"          %(fd2,sx2,sy2,txy2)
                #print "len(data) = ",len(self.data)
                #self.printBlock(self.data)
                self.obj.addNewNode(dt, eid, grid, fd1, sx1, sy1, txy1)
                self.obj.add(dt, eid, grid, fd2, sx2, sy2, txy2)

            #print '--------------------'
            #print "len(data) = ",len(self.data)
            #print "tell = ",self.op2.tell()

            #self.print_section(100)
            #self.dn += 348

    def OES_CQUAD4NL_90_alt(self):
        dt = self.nonlinearFactor
        (format1, extract) = self.getOUG_FormatStart()

        assert self.numWide == 25, "numWide=%s not 25" % (self.numWide)
        nTotal = 100  # 4*25
        format1 += '24f'
        format1 = bytes(format1)

        while len(self.data) >= nTotal:
            eData = self.data[0:nTotal]
            self.data = self.data[nTotal:]

            out = unpack(format1, eData)  # numWide=25
            (eid, fd1, sx1, sy1, xxx, txy1, es1, eps1, ecs1, ex1, ey1, xxx, exy1,
                  fd2, sx2, sy2, xxx, txy2, es2, eps2, ecs2, ex2, ey2, xxx, exy2) = out
            eid = extract(eid, dt)

            data = (eid, fd1, sx1, sy1, xxx, txy1, es1, eps1,
                    ecs1, ex1, ey1, xxx, exy1)
            self.obj.addNewEid(self.elementType, dt, data)
            data = (eid, fd2, sx2, sy2, xxx, txy2, es2, eps2,
                    ecs2, ex2, ey2, xxx, exy2)
            self.obj.add(dt, data)

            #print "eid=%s axial=%s equivStress=%s totalStrain=%s effPlasticCreepStrain=%s effCreepStrain=%s linearTorsionalStresss=%s" %(eid,axial,equivStress,totalStrain,effPlasticCreepStrain,effCreepStrain,linearTorsionalStresss)

            if self.makeOp2Debug:
                self.op2Debug.write('%s\n' % (str(out)))

    def OES_CBUSH1D_40_alt(self):
        dt = self.nonlinearFactor
        (format1, extract) = self.getOUG_FormatStart()
        isMagnitudePhase = self.isMagnitudePhase()

        assert self.numWide == 9, "numWide=%s not 9" % (self.numWide)
        nTotal = 36  # 4*9
        format1 += '8f'
        format1 = bytes(format1)

        while len(self.data) >= nTotal:
            eData = self.data[0:nTotal]
            self.data = self.data[nTotal:]

            out = unpack(format1, eData)  # numWide=25
            (eid, fer, uer, aor, aer,
                  fei, uei, aoi, aei) = out
            eid = extract(eid, dt)
            
            if isMagnitudePhase:
                fe = polarToRealImag(fer, fei)
                ue = polarToRealImag(uer, uei)
                ao = polarToRealImag(aor, aoi)
                ae = polarToRealImag(aer, aei)
            else:
                fe = complex(fer, fei)
                ue = complex(uer, uei)
                ao = complex(aor, aoi)
                ae = complex(aer, aei)

            self.obj.addNewEid(self.elementType, dt, eid, fe, ue, ao, ae)

    def OES_CBUSH_102_alt(self):
        dt = self.nonlinearFactor
        (format1, extract) = self.getOUG_FormatStart()
        isMagnitudePhase = self.isMagnitudePhase()

        assert self.numWide == 13, "numWide=%s not 14" % (self.numWide)
        nTotal = 52  # 4*13
        format1 += '12f'
        format1 = bytes(format1)

        while len(self.data) >= nTotal:
            eData = self.data[0:nTotal]
            self.data = self.data[nTotal:]

            out = unpack(format1, eData)  # numWide=25
            (eid, txr,tyr,tzr,rxr,ryr,rzr,
                  txi,tyi,tzi,rxi,ryi,rzi) = out
            eid = extract(eid, dt)
            
            if isMagnitudePhase:
                tx = polarToRealImag(txr, txi)
                ty = polarToRealImag(tyr, tyi)
                tz = polarToRealImag(tzr, tzi)
                rx = polarToRealImag(rxr, rxi)
                ry = polarToRealImag(ryr, ryi)
                rz = polarToRealImag(rzr, rzi)
            else:
                tx = complex(txr, txi)
                ty = complex(tyr, tyi)
                tz = complex(tzr, tzi)
                rx = complex(rxr, rxi)
                ry = complex(ryr, ryi)
                rz = complex(rzr, rzi)

            #data = (eid, tx, ty, tz, rx, ry, rz)
            self.obj.addNewEid(self.elementType, dt, eid, tx, ty, tz, rx, ry, rz)

    def OES_CQUAD4_144_alt(self):
        """
        GRID-ID  DISTANCE,NORMAL-X,NORMAL-Y,SHEAR-XY,ANGLE,MAJOR MINOR,VONMISES
        """
        if self.makeOp2Debug:
            self.op2Debug.write('---CQUAD4_144---\n')

        #self.print_section(20)
        #term = data[0:4] CEN/
        #data = data[4:]
        #print "*****"
        #self.printBlock(self.data)
        #assert self.numWide==87,'invalid numWide...numWide=%s' %(self.numWide)
        #if self.numWide==87: # 2+(17-1)*5 = 87 -> 87*4 = 348

        if self.elementType == 144:  # CQUAD4
            nTotal = 308  # 2+15*5 = 77 -> 87*4 = 308
            nNodes = 4    # centroid + 4 corner points
            eType = 'CQUAD4'
        elif self.elementType == 64:  # CQUAD8 - done
            nTotal = 308  # 2+15*5 = 77 -> 77*4 = 308
            nNodes = 4    # centroid + 4 corner points
            eType = 'CQUAD8'
        elif self.elementType == 82:  # CQUADR
            nTotal = 308  # 2+15*5 = 77 -> 87*4 = 308
            nNodes = 4    # centroid + 4 corner points
            eType = 'CQUAD4'  ## @todo write the word CQUADR

        elif self.elementType == 75:  # CTRIA6
            nTotal = 248  # 2+15*3 = 62 -> 62*4 = 248
            nNodes = 3    # centroid + 3 corner points
            eType = 'CTRIA6'
        elif self.elementType == 70:  # CTRIAR
            nTotal = 248  # 2+15*4 = 62 -> 62*4 = 248
            nNodes = 3    # centroid + 3 corner points
            eType = 'CTRIAR'  ## @todo write the word CTRIAR
        else:
            raise RuntimeError('elementType=%s nTotal not defined...' %
                            (self.elementType))

        assert nTotal == self.numWide * 4, 'eType=%s numWide*4=%s not nTotal=%s' % (self.elementType, self.numWide * 4, nTotal)
        dt = self.nonlinearFactor
        (format1, extract) = self.getOUG_FormatStart()
        format1 += '14f'
        format1 = bytes(format1)
        isMagnitudePhase = self.isMagnitudePhase()

        while len(self.data) >= nTotal:
            (eid, _) = unpack(b'i4s', self.data[0:8])
            self.data = self.data[8:]  # 2
            eid = extract(eid, dt)
            eData = self.data[0:60]  # 4*15
            self.data = self.data[60:]
            out = unpack(format1, eData)  # len=15*4
            if self.makeOp2Debug:
                self.op2Debug.write('%s\n' % (str(out)))
            (grid, fd1, sx1r, sx1i, sy1r, sy1i, txy1r, txy1i,
             fd2, sx2r, sx2i, sy2r, sy2i, txy2r, txy2i) = out
            grid = 'C'

            if isMagnitudePhase:
                sx1 = polarToRealImag(sx1r, sx1i)
                sy1 = polarToRealImag(
                    sy1r, sy1i)
                sx2 = polarToRealImag(sx2r, sx2i)
                sy2 = polarToRealImag(
                    sy2r, sy2i)
                txy1 = polarToRealImag(txy1r, txy1i)
                txy2 = polarToRealImag(
                    txy2r, txy2i)
            else:
                sx1 = complex(sx1r, sx1i)
                sy1 = complex(sy1r, sy1i)
                sx2 = complex(sx2r, sx2i)
                sy2 = complex(sy2r, sy2i)
                txy1 = complex(txy1r, txy1i)
                txy2 = complex(txy2r, txy2i)

            self.obj.addNewEid(eType, dt, eid, grid, fd1, sx1, sy1, txy1)
            self.obj.add(dt, eid, grid, fd2, sx2, sy2, txy2)

            for nodeID in xrange(nNodes):  # nodes pts
                eData = self.data[0:60]  # 4*15=60
                self.data = self.data[60:]
                out = unpack(b'i14f', eData)
                if self.makeOp2Debug:
                    self.op2Debug.write('%s\n' % (str(out)))
                (grid, fd1, sx1r, sx1i, sy1r, sy1i, txy1r, txy1i,
                 fd2, sx2r, sx2i, sy2r, sy2i, txy2r, txy2i) = out

                if isMagnitudePhase:
                    sx1 = polarToRealImag(sx1r, sx1i)
                    sx2 = polarToRealImag(sx2r, sx2i)
                    sy1 = polarToRealImag(sy1r, sy1i)
                    sy2 = polarToRealImag(sy2r, sy2i)
                    txy1 = polarToRealImag(txy1r, txy1i)
                    txy2 = polarToRealImag(txy2r, txy2i)
                else:
                    sx1 = complex(sx1r, sx1i)
                    sx2 = complex(sx2r, sx2i)
                    sy1 = complex(sy1r, sy1i)
                    sy2 = complex(sy2r, sy2i)
                    txy1 = complex(txy1r, txy1i)
                    txy2 = complex(txy2r, txy2i)

                #print "eid=%i grid=%i fd1=%i sx1=%i sy1=%i txy1=%i"   %(eid,grid,fd1,sx1,sy1,txy1)
                #print "               fd2=%i sx2=%i sy2=%i txy2=%i\n"          %(fd2,sx2,sy2,txy2)
                #print "len(data) = ",len(self.data)
                #self.printBlock(self.data)
                self.obj.addNewNode(dt, eid, grid, fd1, sx1, sy1, txy1)
                self.obj.add(dt, eid, grid, fd2, sx2, sy2, txy2)

            #print '--------------------'
            #print "len(data) = ",len(self.data)
            #self.dn += 348

    def OES_CTRIA3_74_alt(self):  # in progress
        """
        DISTANCE,NORMAL-X,NORMAL-Y,SHEAR-XY,ANGLE,MAJOR,MINOR,VONMISES
        stress is extracted at the centroid
        """
        assert self.numWide == 15, 'invalid numWide...numWide=%s' % (
            self.numWide)

        dt = self.nonlinearFactor
        (format1, extract) = self.getOUG_FormatStart()
        format1 += '14f'
        format1 = bytes(format1)
        isMagnitudePhase = self.isMagnitudePhase()

        while len(self.data) >= 60:
            eData = self.data[0:60]  # 4*15=60
            self.data = self.data[60:]
            out = unpack(format1, eData)

            (eid, fd1, sx1r, sx1i, sy1r, sy1i, txy1r, txy1i,
             fd2, sx2r, sx2i, sy2r, sy2i, txy2r, txy2i) = out

            if isMagnitudePhase:
                sx1 = polarToRealImag(sx1r, sx1i)
                sy1 = polarToRealImag(sy1r, sy1i)
                sx2 = polarToRealImag(sx2r, sx2i)
                sy2 = polarToRealImag(sy2r, sy2i)
                txy1 = polarToRealImag(txy1r, txy1i)
                txy2 = polarToRealImag(txy2r, txy2i)
            else:
                sx1 = complex(sx1r, sx1i)
                sy1 = complex(sy1r, sy1i)
                sx2 = complex(sx2r, sx2i)
                sy2 = complex(sy2r, sy2i)
                txy1 = complex(txy1r, txy1i)
                txy2 = complex(txy2r, txy2i)

            eid = extract(eid, dt)
            #print "eid=%i fd1=%i sx1=%i sy1=%i txy1=%i" %(eid,fd1,sx1,sy1,txy1)
            #print  "      fd2=%i sx2=%i sy2=%i txy2=%i\n"   %(fd2,sx2,sy2,txy2)
            self.obj.addNewEid('CTRIA3', dt, eid, 'C', fd1, sx1, sy1, txy1)
            self.obj.add(dt, eid, 'C', fd2, sx2, sy2, txy2)
            if self.makeOp2Debug:
                self.op2Debug.write('%s\n' % str(out))