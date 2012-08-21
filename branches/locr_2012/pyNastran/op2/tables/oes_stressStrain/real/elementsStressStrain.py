#pylint: disable=C0103
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
from struct import unpack

#91  -> PENTANL
#2   -> BEAM
#33  -> TUBE
#92  -> CONRODNL


class RealElementsStressStrain(object):

    def skipOES_Element(self):
        self.log.debug('skipping approach/table/format/sortCode=%s on %s table' % (self.atfsCode, self.tableName))
        #print 'skipping approach/table/format/sortCode=%s on %s table' % (self.atfsCode,self.tableName)
        #print self.codeInformation()
        self.handleResultsBuffer3(self.dummyPass, None, debug=True)

    def dummyPass(self):
        self.data = b''

    def OES_Thermal(self, debug=False):
        #assert self.numWide==5,'invalid numWide...numWide=%s' % (self.numWide)

        dt = self.nonlinearFactor
        (format1, extract) = self.getOUG_FormatStart()
        format1 += '2i5f'  # 1+2+5=8
        format1 = bytes(format1)
        while len(self.data) >= 32:  # 4*8
            #print self.printSection(40)
            eData = self.data[0:32]
            self.data = self.data[32:]
            #print "len(data) = ",len(eData)

            out = unpack(format1, eData)
            (eid, sideID, hbdyID, cnvCoeff, fApplied, fConv,
                fRad, fTotal) = out
            eid = extract(eid, dt)
            #print "eid=%s sideID=%s hbdyID=%s coeff=%s fApplied=%s fConv=%s fRad=%s fTotal=%s" % (eid,sideID,hbdyID,cnvCoeff,fApplied,fConv,fRad,fTotal)
            #self.obj.addNewEid(eid,axial,axialMS,torsion,torsionMS)

            #print "eid=%i axial=%i torsion=%i" % (eid,axial,torsion)
            #print "len(data) = ",len(self.data)
        ###
        #print self.rodStress[self.iSubcase]

    def OES_basicElement(self):
        """
        genericStressReader - works on CROD_1, CELAS2_12
        stress & strain
        formatCode=1 sortCode=0 (eid,axial,axialMS,torsion,torsionMS)
        """
        dt = self.nonlinearFactor
        (format1, extract) = self.getOUG_FormatStart()
        (nTotal, dataFormat) = self.obj.getLength()
        dataFormat = format1 + dataFormat
        #print "nTotal=%s dataFormat=%s len(data)=%s" % (nTotal,dataFormat,len(self.data))
        dataFormat = bytes(dataFormat)
        n = 0
        nEntries = len(self.data) // nTotal
        for i in xrange(nEntries):
            eData = self.data[n:n + nTotal]
            out = unpack(dataFormat, eData)
            #print "out = ",out
            eid = extract(out[0], dt)
            self.obj.addNewEid(dt, eid, out[1:])
            n += nTotal
        self.data = self.data[n:]

    def OES_CBEAM_2(self):
        dt = self.nonlinearFactor
        (formatStart, extract) = self.getOUG_FormatStart()

        nNodes = 10  # 11-1
        nTotal = self.obj.getLengthTotal()
        (n1, format1) = self.obj.getLength1()
        (n2, format2) = self.obj.getLength2()
        format1 = formatStart + format1
        format1 = bytes(format1)
        format2 = bytes(format2)
        while len(self.data) >= nTotal:
            eData = self.data[0:n1]
            self.data = self.data[n1:]
            #print "len(data) = ",len(eData)

            out = unpack(format1, eData)
            #print "outA = ",out
            eid2 = extract(out[0], dt)
            self.obj.addNewEid(dt, eid2, out[1:])

            for iNode in xrange(nNodes):
                eData = self.data[0:n2]
                self.data = self.data[n2:]
                out = unpack(format2, eData)
                #print "outB = ",out
                self.obj.add(dt, eid2, out)

            #print "eid=%i axial=%i torsion=%i" % (eid,axial,torsion)
            #print "len(data) = ",len(self.data)

    def OES_CQUAD4_33(self):
        """
        GRID-ID  DISTANCE,NORMAL-X,NORMAL-Y,SHEAR-XY,ANGLE,MAJOR MINOR,VONMISES
        """
        dt = self.nonlinearFactor
        (format1, extract) = self.getOUG_FormatStart()
        format1 += '16f'
        format1 = bytes(format1)

        nNodes = 0  # centroid + 4 corner points
        #self.printSection(20)
        #term = data[0:4] CEN/
        #data = data[4:]
        #print "*****"
        #self.printBlock(self.data)

        assert self.numWide == 17, ('invalid numWide...numWide=%s'
                                    % (self.numWide))

        while len(self.data) >= 68:  # 2+17*5 = 87 -> 87*4 = 348
            #print self.printBlock(self.data[0:100])
            #(eid,) = unpack(b'i',self.data[0:4])
            #self.data = self.data[8:]  # 2
            eData = self.data[0:68]  # 4*17
            self.data = self.data[68:]
            out = unpack(format1, eData)  # 17
            (eid, fd1, sx1, sy1, txy1, angle1, major1, minor1, maxShear1,
             fd2, sx2, sy2, txy2, angle2, major2, minor2, maxShear2) = out

            eid = extract(eid, dt)

            #print "eid=%i grid=%s fd1=%-3.1f sx1=%i sy1=%i txy1=%i angle1=%i major1=%i minor1=%i vm1=%i" % (eid,'C',fd1,sx1,sy1,txy1,angle1,major1,minor1,maxShear1)
            #print   "             fd2=%-3.1f sx2=%i sy2=%i txy2=%i angle2=%i major2=%i minor2=%i vm2=%i\n"       % (fd2,sx2,sy2,txy2,angle2,major2,minor2,maxShear2)
            #print "nNodes = ",nNodes
            self.obj.addNewEid('CQUAD4', dt, eid, 'C', fd1, sx1, sy1,
                               txy1, angle1, major1, minor1, maxShear1)
            self.obj.add(dt, eid, 'C', fd2, sx2, sy2, txy2,
                         angle2, major2, minor2, maxShear2)

            for nodeID in xrange(nNodes):  # nodes pts
                eData = self.data[0:68]  # 4*17
                self.data = self.data[68:]
                out = unpack(b'i16f', eData[0:68])
                (grid, fd1, sx1, sy1, txy1, angle1, major1, minor1, vm1,
                 fd2, sx2, sy2, txy2, angle2, major2, minor2, vm2,) = out

                #print "eid=%i grid=%i fd1=%i sx1=%i sy1=%i txy1=%i angle1=%i major1=%i minor1=%i vm1=%i" % (eid,grid,fd1,sx1,sy1,txy1,angle1,major1,minor1,vm1)
                #print "               fd2=%i sx2=%i sy2=%i txy2=%i angle2=%i major2=%i minor2=%i vm2=%i\n"          % (fd2,sx2,sy2,txy2,angle2,major2,minor2,vm2)
                #print "len(data) = ",len(self.data)
                #self.printBlock(self.data)
                self.obj.addNewNode(dt, eid, grid, fd1, sx1, sy1,
                                    txy1, angle1, major1, minor1, vm1)
                self.obj.add(dt, eid, grid, fd2, sx2, sy2,
                             txy2, angle2, major2, minor2, vm2)
            #print '--------------------'
            #print "len(data) = ",len(self.data)
            #print "tell = ",self.op2.tell()

            #self.printSection(100)
            #self.dn += 348
        ###

    def OES_CBAR_34(self):
        dt = self.nonlinearFactor
        #print "len(data) = ",len(self.data)
        assert self.numWide == 16, ('invalid numWide...numWide=%s'
                                    % (self.numWide))

        (format1, extract) = self.getOUG_FormatStart()
        format1 += '15f'
        format1 = bytes(format1)

        while len(self.data) >= 64:
            #self.printBlock(self.data)
            eData = self.data[0:64]
            self.data = self.data[64:]
            #print "len(data) = ",len(eData)

            (eid, s1a, s2a, s3a, s4a, axial, smaxa, smina, MSt,
             s1b, s2b, s3b, s4b, smaxb, sminb, MSc) = unpack(format1, eData)
            eid2 = extract(eid, dt)
            self.obj.addNewEid('CBAR', dt, eid2, s1a, s2a, s3a, s4a, axial, smaxa, smina, MSt,
                               s1b, s2b, s3b, s4b, smaxb, sminb, MSc)

            #print "eid=%i s1=%i s2=%i s3=%i s4=%i axial=%-5i smax=%i smin=%i MSt=%i MSc=%i" % (eid,s1a,s2a,s3a,s4a,axial,smaxa,smina,MSt,MSc)
            #print "         s1=%i s2=%i s3=%i s4=%i          smax=%i smin=%i" % (s1b,s2b,s3b,s4b,smaxb,sminb)
            #print "len(data) = ",len(self.data)

    def OES_CSOLID_67(self):
        """
        stress is extracted at the centroid
        CTETRA_39
        CPENTA_67
        CHEXA_68
        """
        dt = self.nonlinearFactor
        (format1, extract) = self.getOUG_FormatStart()
        format1 += "i4si"
        format1 = bytes(format1)

        #nNodes = 5 # 1 centroid + 4 corner points
        #self.printSection(20)
        #term      = self.data[0:4] CEN/
        #self.data = self.data[4:]
        #print "*****"
        ElementType = self.ElementType(self.elementType)
        nNodes = 1  # this is a minimum, it will be reset later
        nNodesExpected = 1
        #assert self.numWide in [109,151,193],'invalid numWide...numWide=%s' % (self.numWide)
        while len(self.data) >= 16 + 84 * nNodesExpected:
            eData = self.data[0:16]
            self.data = self.data[16:]
            #self.printBlock(eData)

            out = unpack(format1, eData)
            (eid, cid, abcd, nNodes) = out
            eid2 = extract(eid, dt)
            #print "abcd = |%s|" % (abcd)
            #print "eid=%s cid=%s nNodes=%s nNodesExpected=%s" % (eid,cid,nNodes,nNodesExpected)

            assert nNodes < 21, self.printBlock(eData)

            if   ElementType == 'TETRA':
                nNodesExpected = 5
            elif ElementType == 'PENTA':
                nNodesExpected = 7
            elif ElementType == 'HEXA':
                nNodesExpected = 9
            else:
                raise NotImplementedError('not supported....EType=%s eType=%s nNodes=%s numWide=%s' % (ElementType, self.elementType, nNodes, self.numWide))

            #print "len(data) = ",len(self.data)
            for nodeID in xrange(nNodesExpected):  # nodes pts, +1 for centroid (???)
                #print "len(data)A = ",len(self.data)
                eData = self.data[0:4 * 21]  # for the stresses
                self.data = self.data[4 * 21:]
                #print "len(data)B = ",len(self.data)
                #self.printBlock(eData)

                #print "self.tableCode = ",self.tableCode
                #print "len(data) = ",len(self.data)

                gridDevice, = unpack(b'i', eData[0:4])
                #print "gridDevice = ",gridDevice
                if gridDevice == 0:
                    grid = 'C'
                else:
                    #grid = (gridDevice - deviceCode) // 10
                    grid = gridDevice

                out = unpack(b'20f', eData[4:4 * 21])
                (sxx, sxy, s1, a1, a2, a3, pressure, svm,
                 syy, syz, s2, b1, b2, b3,
                 szz, sxz, s3, c1, c2, c3) = out

                #print "%s eid=%s cid=%s grid=%s sxx=%-6i txy=%-5i s1=%-6i a1=%i a2=%i a3=%i press=%i vm=%s" % (elementType,eid,cid,grid,sxx,sxy,s1,a1,a2,a3,pressure,svm)
                #print "%s eid=%s cid=%s grid=%s syy=%-6i tyz=%-5i s2=%-6i b1=%i b2=%i b3=%i"                % (elementType,eid,cid,grid,syy,syz,s2,b1,b2,b3)
                #print "%s eid=%s cid=%s grid=%s szz=%-6i txz=%-5i s3=%-6i c1=%i c2=%i c3=%i"                % (elementType,eid,cid,grid,szz,sxz,s3,c1,c2,c3)
                #print ""

                #smax = max(s1,s2,s3)
                #smin = min(s1,s2,s3)

                aCos = []
                bCos = []
                cCos = []
                if nodeID == 0:
                    #print "adding new eid"
                    self.obj.addNewEid(ElementType, cid, dt, eid2, grid, sxx, syy, szz, sxy, syz, sxz, s1, s2, s3, aCos, bCos, cCos, pressure, svm)
                else:
                    self.obj.add(dt, eid2, grid, sxx, syy, szz, sxy, syz, sxz, s1, s2, s3, aCos, bCos, cCos, pressure, svm)
                #self.printBlock(data)
            #print self.solidStress[self.iSubcase]
            #print '--------------------'
            #print "len(data) = ",len(self.data)

            #self.printSection(100)
            #self.printBlock(self.data[0:100])
            #self.printBlock(self.data[1:100])
            #self.printBlock(self.data[2:100])
            #self.printBlock(self.data[3:100])
        ###

    def OES_CTRIA3_74(self):
        """
        DISTANCE,NORMAL-X,NORMAL-Y,SHEAR-XY,ANGLE,MAJOR,MINOR,VONMISES
        stress is extracted at the centroid
        """
        assert self.numWide == 17, ('invalid numWide...numWide=%s'
                                    % (self.numWide))

        dt = self.nonlinearFactor
        (format1, extract) = self.getOUG_FormatStart()
        format1 += '16f'
        format1 = bytes(format1)

        while len(self.data) >= 68:  # 4*17
            eData = self.data[0:68]
            self.data = self.data[68:]
            out = unpack(format1, eData)

            (eid, fd1, sx1, sy1, txy1, angle1, major1, minor1, vm1,
             fd2, sx2, sy2, txy2, angle2, major2, minor2, vm2,) = out
            eid = extract(eid, dt)
            #print "eid=%i fd1=%i sx1=%i sy1=%i txy1=%i angle1=%i major1=%i minor1=%i vm1=%i" % (eid,fd1,sx1,sy1,txy1,angle1,major1,minor1,vm1)
            #print  "      fd2=%i sx2=%i sy2=%i txy2=%i angle2=%i major2=%i minor2=%i vm2=%i\n"   % (fd2,sx2,sy2,txy2,angle2,major2,minor2,vm2)
            self.obj.addNewEid('CTRIA3', dt, eid, 'C', fd1, sx1, sy1,
                               txy1, angle1, major1, minor1, vm1)
            self.obj.add(dt, eid, 'C', fd2, sx2, sy2, txy2,
                         angle2, major2, minor2, vm2)

    def OES_CSOLID_85(self):  # works
        """
        stress is extracted at the centroid
        CTETRA_85
        CPENTA_91 ???
        CHEXA_93  ???
        """
        #print "starting nonlinear solid element..."
        #nNodes = 5 # 1 centroid + 4 corner points
        #self.printSection(20)
        #term      = self.data[0:4] CEN/
        #self.data = self.data[4:]
        #print "*****"
        dt = self.nonlinearFactor
        (format1, extract) = self.getOUG_FormatStart()
        format1 += "4s"
        format1 = bytes(format1)

        nNodes = 4  # this is a minimum, it will be reset later
        nNodesExpected = 1
        assert self.numWide in [82], 'invalid numWide...numWide=%s' % (self.numWide)
        while len(self.data) >= 16 + 84 * nNodesExpected:
            eData = self.data[0:8]
            self.data = self.data[8:]
            #self.printBlock(eData)

            out = unpack(format1, eData)
            (eid, abcd) = out
            #print "abcd = |%s|" % (a+b+c+d)
            #print "eid=%s cid=%s nNodes=%s nNodesExpected=%s" % (eid,cid,nNodes,nNodesExpected)

            assert nNodes < 21, self.printBlock(eData)
            eid = extract(eid, dt)
            if(nNodes in [4, 10]):
                elementType = "CTETRA"
                nNodesExpected = 5
            #elif(nNodes in [6,15]):
            #    elementType = "CPENTA"
            #    nNodesExpected = 7
            #elif(nNodes in [8,20]):
            #    elementType = "CHEXA"
            #    nNodesExpected = 9
            #else:
            #    raise Exception('not supported....nNodes=%s' % (nNodes))

            #print "len(data) = ",len(self.data)
            for nodeID in xrange(nNodesExpected):  # nodes pts, +1 for centroid (???)
                #print "len(data)A = ",len(self.data)
                eData = self.data[0:64]  # 4*16; for the stresses
                self.data = self.data[64:]
                #print "len(data)B = ",len(self.data)
                #self.printBlock(eData)

                #print "self.tableCode = ",self.tableCode
                #print "len(data) = ",len(self.data)
                #gridDevice, = unpack(b'i', eData[0:4])
                #print "gridDevice = ",gridDevice

                out = unpack(b'i15f', eData[:64])  # 18-3 = 15
                (gridGauss, sxx, syy, szz, sxy, syz, sxz, se,
                 eps, ecs, exx, eyy, ezz, exy, eyz, exz) = out

                if gridGauss == 0:
                    gridGauss = 'C'
                else:
                    #grid = (gridDevice - deviceCode) // 10
                    gridGauss = gridGauss

                #print "%s gridGauss=%-5s eid=%s sxx=%g syy=%g szz=%g" % (elementType,gridGauss,eid,sxx,syy,szz)

                #print "%s eid=%s cid=%s grid=%s sxx=%-6i txy=%-5i s1=%-6i a1=%i a2=%i a3=%i press=%i vm=%s" % (elementType,eid,cid,grid,sxx,sxy,s1,a1,a2,a3,pressure,svm)
                #print "%s eid=%s cid=%s grid=%s syy=%-6i tyz=%-5i s2=%-6i b1=%i b2=%i b3=%i"                % (elementType,eid,cid,grid,syy,syz,s2,b1,b2,b3)
                #print "%s eid=%s cid=%s grid=%s szz=%-6i txz=%-5i s3=%-6i c1=%i c2=%i c3=%i"                % (elementType,eid,cid,grid,szz,sxz,s3,c1,c2,c3)
                #print ""

                #smax = max(s1,s2,s3)
                #smin = min(s1,s2,s3)

                #if nodeID==0:
                #    #print "adding new eid"
                #    self.obj.addNewEid(elementType,cid,eid,grid,sxx,syy,szz,sxy,syz,sxz,aCos,bCos,cCos,pressure,svm)
                #else:
                #    self.obj.add(                      eid,grid,sxx,syy,szz,sxy,syz,sxz,aCos,bCos,cCos,pressure,svm)
                #print "eid=%i grid=%i fd1=%i sx1=%i sy1=%i txy1=%i angle1=%i major1=%i minor1=%i vm1=%i" % (eid,grid,fd1,sx1,sy1,txy1,angle1,major1,minor1,vm1)
                #print "               fd2=%i sx2=%i sy2=%i txy2=%i angle2=%i major2=%i minor2=%i vm2=%i\n"          % (fd2,sx2,sy2,txy2,angle2,major2,minor2,vm2)
                #self.printBlock(data)
            #sys.exit('finished a CTETRANL')
            #print self.solidStress[self.iSubcase]
            #print '--------------------'
            #print "len(data) = ",len(self.data)

            #self.printSection(100)
            #self.printBlock(self.data[0:100])
            #self.printBlock(self.data[1:100])
            #self.printBlock(self.data[2:100])
            #self.printBlock(self.data[3:100])
        ###
        #print self.solidStress[self.iSubcase]

    def OES_field1(self):
        if self.isSort1():
            #raise NotImplementedError('SORT1 is not supported')
            return ('i', self.scaleEid)
        elif self.isSort2():
            if self.analysisCode in [1, 2, 3, 4, 7, 8, 9, 11, 12]:  # eid
                return ('i', self.scaleEid)
            elif self.analysisCode == 5:  # freq
                #freq
                return ('f', self.scaleDt)
                #raise NotImplementedError('freq is not supported')
            elif self.analysisCode == 6:  # time
                #time
                return ('f', self.scaleDt)
                #raise NotImplementedError('time is not supported')
            elif self.analysisCode == 10:  # fqts:
                #fqts # freqTime
                return ('f', self.scaleDt)
                #raise NotImplementedError('freqTime is not supported')
            else:
                raise NotImplementedError('invalid SORT2 analysisCode=%s' %
                                          (self.analysisCode))
            ###
        else:
            raise NotImplementedError('invalid SORTx code')
        ###

    def OES_CTRIAX6_53(self):
        #(Format1,scaleValue) = self.OES_field1()
        #Format = Format1+'ifffffff'
        dt = self.nonlinearFactor
        (format1, extract) = self.getOUG_FormatStart()
        format1 += 'i7f'
        format1 = bytes(format1)
        while len(self.data) >= 132:  # (1+8*4)*4 = 33*4 = 132
            eData = self.data[0:36]  # 4*9
            self.data = self.data[36:]
            out = unpack(format1, eData)
            (eid, loc, rs, azs, As, ss, maxp, tmax, octs) = out
            eid = extract(eid, dt)
            #print "eid=%s loc=%s rs=%s azs=%s as=%s ss=%s maxp=%s tmx=%s octs=%s" % (eid,loc,rs,azs,As,ss,maxp,tmax,octs)
            self.obj.addNewEid(dt, eid, loc, rs, azs, As, ss, maxp, tmax, octs)

            for i in xrange(3):
                eData = self.data[0:32]  # 4*8
                self.data = self.data[32:]
                out = unpack(b'i7f', eData)
                (loc, rs, azs, As, ss, maxp, tmax, octs) = out
                #print "eid=%s loc=%s rs=%s azs=%s as=%s ss=%s maxp=%s tmx=%s octs=%s" % (eid,loc,rs,azs,As,ss,maxp,tmax,octs)
                self.obj.add(dt, eid, loc, rs, azs, As, ss, maxp, tmax, octs)

    def OES_RODNL_89_92(self):
        dt = self.nonlinearFactor
        (format1, extract) = self.getOUG_FormatStart()
        format1 += '6f'  # 1+6=7
        format1 = bytes(format1)

        while len(self.data) >= 28:  # len(format1)*4 = 7*4 = 28
            eData = self.data[0:28]
            self.data = self.data[28:]
            out = unpack(format1, eData)

            (eid, axial, equivStress, totalStrain, effPlasticCreepStrain,
                effCreepStrain, linearTorsionalStresss) = out
            eid = extract(eid, dt)
            data = (eid, axial, equivStress, totalStrain, effPlasticCreepStrain, effCreepStrain, linearTorsionalStresss)

            #print "eid=%s axial=%s equivStress=%s totalStrain=%s effPlasticCreepStrain=%s effCreepStrain=%s linearTorsionalStresss=%s" % (eid,axial,equivStress,totalStrain,effPlasticCreepStrain,effCreepStrain,linearTorsionalStresss)
            self.obj.add(self.elementType, dt, data)

    def OES_CQUAD4NL_90(self):
        dt = self.nonlinearFactor
        (format1, extract) = self.getOUG_FormatStart()

        nStep = 52  # 4*13
        format1 += '12f'  # 1+12=13
        format1 = bytes(format1)
        assert 13 == self.numWide, 'numWide=%s not 13' % (self.numWide)

        while len(self.data) >= nStep:
            eData = self.data[0:nStep]
            self.data = self.data[nStep:]

            out = unpack(format1, eData)  # numWide=13
            (eid, fd1, sx1, sy1, sz1, txy1, es1, eps1, ecs1,
                ex1, ey1, ez1, exy1) = out
            eid = extract(eid, dt)
            data = (eid, fd1, sx1, sy1, sz1, txy1, es1, eps1,
                    ecs1, ex1, ey1, ez1, exy1)
            self.obj.addNewEid(self.elementType, dt, data)
            #print "eid=%s axial=%s equivStress=%s totalStrain=%s effPlasticCreepStrain=%s effCreepStrain=%s linearTorsionalStresss=%s" % (eid,axial,equivStress,totalStrain,effPlasticCreepStrain,effCreepStrain,linearTorsionalStresss)

    def OES_CPENTANL_91(self):
        """
        The DMAP manual says fields 3-18 repeat 7 times. but they dont.
        They repeat 6 times.  Other DMAP cards are correct with
        their repeat statements.
        """
        n = 0
        dt = self.nonlinearFactor
        (format1, extract) = self.getOUG_FormatStart()
        format1 += '4s'
        format1 = bytes(format1)

        while len(self.data) >= 456:  # 2+16*7 = 114 -> 114*4 = 456
            eData = self.data[0:8]
            self.data = self.data[8:]
            (eid, a, b, c, d) = unpack(format1, eData)
            eid = extract(eid, dt)
            #out = unpack(b"ii",eData)
            #(eid,cType) = out
            cType = a + b + c + d

            for i in xrange(7):
                #print "len(self.data) = ",len(self.data)
                eData = self.data[0:64]
                self.data = self.data[64:]
                out = unpack(b'i15f', eData)
                assert len(out) == 16
                (grid, sx, sy, sz, sxy, syz, sxz, se,
                 eps, ecs, ex, ey, ez, exy, eyz, exz) = out
                #print "eid=%3s cType=%s sx=%i sy=%i sz=%i sxy=%s syz=%i szx=%i se=%s" % (eid,cType,sx,sy,sz,sxy,syz,sxz,se)
                #print "gid=%3s ecs=%.3g   ex=%.3g ey=%.3g ez=%.3g exy=%.3g eyz=%.3g ezx=%.3g"  % (grid,ecs,ex,ey,ez,exy,eyz,exz)
                #print ""
                assert a == 'G'

            #self.data = self.data[1456:]
        #self.firstPass = True

    def OES_CHEXANL_93(self):
        """
        The DMAP manual says fields 3-18 repeat 9 times. but they dont.
        They repeat 8 times.  Other DMAP cards are correct with
        their repeat statements.
        """
        #print "CHEXANL_93"
        #print "len(self.data) = ",len(self.data)

        n = 0
        dt = self.nonlinearFactor
        (format1, extract) = self.getOUG_FormatStart()
        format1 += '4s'
        format1 = bytes(format1)

        while len(self.data) >= 584:  # 2+16*9 = 146 -> 146*4 = 584
            eData = self.data[0:8]
            self.data = self.data[8:]
            (eid, cType) = unpack(format1, eData)
            eid = extract(eid, dt)
            #out = unpack('2i',eData)
            #(eid,cType) = out

            for i in xrange(9):
                #print "len(self.data) = ",len(self.data)
                eData = self.data[0:64]
                self.data = self.data[64:]
                out = unpack(b'i15f', eData)
                assert len(out) == 16
                (grid, sx, sy, sz, sxy, syz, sxz, se, eps,
                    ecs, ex, ey, ez, exy, eyz, exz) = out
                #print "eid=%3s cType=%s sx=%i sy=%i sz=%i sxy=%s syz=%i szx=%i se=%s" % (eid,cType,sx,sy,sz,sxy,syz,sxz,se)
                #print "gid=%3s ecs=%.3g   ex=%.3g ey=%.3g ez=%.3g exy=%.3g eyz=%.3g ezx=%.3g"  % (grid,ecs,ex,ey,ez,exy,eyz,exz)
                #print ""
                assert a == 'G'

            #self.data = self.data[1456:]
        #self.firstPass = True

    def OES_CBEAM_94(self):
        nNodes = 10  # 11-1

        #nTotal       = self.obj.getLengthTotal()
        #(n1,format1) = self.obj.getLength1()
        #(n2,format2) = self.obj.getLength2()
        nTotal = 2 * 4 + (18 - 3) * 9 * 4
        nTotal = 204

        n1 = 24
        format1 = '4s5f'
        format1 = bytes(format1)
        #print "len(data) = ",len(self.data)
        while len(self.data) >= nTotal:
            eData = self.data[0:8]
            self.data = self.data[8:]
            (eid, gridA) = unpack(b'2i', eData)
            #print "eid=%s gridA=%s" % (eid,gridA)

            for i in xrange(1):
                for j in xrange(4):  # c,d,e,f @ A;    c,d,e,f @ B
                    eData = self.data[0:n1]
                    self.data = self.data[n1:]
                    out = unpack(format1, eData)
                    (loc, nsx, nse, te, epe, ece) = out
                    #print "loc=%s nsx=%s nse=%s te=%s epe=%s ece=%s" % (loc,nsx,nse,te,epe,ece)
                #self.obj.add(eid,out)

    def OES_CQUAD4_95(self):  # works (doesnt handle all stress/strain cases tho)
        """
        GRID-ID  DISTANCE,NORMAL-X,NORMAL-Y,SHEAR-XY,ANGLE,MAJOR MINOR,VONMISES
        composite quad
        """
        eType = self.ElementType(self.elementType)

        #self.printSection(20)
        #term = data[0:4] CEN/
        #data = data[4:]
        assert self.numWide == 11, 'invalid numWide...numWide=%s' % (
            self.numWide)

        dt = self.nonlinearFactor
        (format1, extract) = self.getOUG_FormatStart()
        format1 += 'i9f'
        format1 = bytes(format1)

        while len(self.data) >= 44:  # 2+17*5 = 87 -> 87*4 = 348
            eData = self.data[0:44]  # 4*11
            self.data = self.data[44:]
            out = unpack(format1, eData)
            (eid, iLayer, o1, o2, t12, t1z, t2z, angle, major,
                minor, ovm) = out
            #print "out =",out
            eid = extract(eid, dt)

            if eid != self.eid2:  # originally initialized to None, the buffer doesnt reset it, so it is the old value
                #print "1 - eid=%s iLayer=%i o1=%i o2=%i ovm=%i" % (eid,iLayer,o1,o2,ovm)
                self.obj.addNewEid(eType, dt, eid, o1, o2,
                                   t12, t1z, t2z, angle, major, minor, ovm)
            else:
                #print "2 - eid=%s iLayer=%i o1=%i o2=%i ovm=%i" % (eid,iLayer,o1,o2,ovm)
                self.obj.add(dt, eid, o1, o2, t12, t1z,
                             t2z, angle, major, minor, ovm)
            self.eid2 = eid
            #self.dn += 348
        ###
        #print "3 - eid=%s iLayer=%i o1=%i o2=%i ovm=%i" % (eid,iLayer,o1,o2,ovm)

    def OES_QUAD4FD_139(self):  # hyperelastic
        """
        Hyperelastic Quad
        36+4*7*4 = 148
        """
        #x = 0
        dt = self.nonlinearFactor
        (format1, extract) = self.getOUG_FormatStart()
        format1 += '4si6f'  # 1 + 4+1+6 = 12
        format1 = bytes(format1)
        while len(self.data) >= 148:
            #if x==2:
            #    sys.exit('end of hyperQuad')
            eData = self.data[0:36]  # 4*9
            self.data = self.data[36:]
            out = unpack(format1, eData)

            (eid, Type, ID, sx, sy, sxy, angle, smj, smi) = out
            eid = extract(eid, dt)
            self.obj.addNewEid(dt, [eid, Type, sx, sy, sxy, angle, smj, smi])
            #print "eid=%s Type=%s\n***ID=%s sx=%s sy=%s sxy=%s angle=%s major=%s minor=%s" % (eid,Type,ID,sx,sy,sxy,angle,smj,smi)
            for i in xrange(3):
                eData = self.data[0:28]  # 4*7
                self.data = self.data[28:]
                out = unpack(b'i6f', eData)
                #(ID,sx,sy,sxy,angle,smj,smi) = out
                self.obj.add(dt, eid, out)
                #print "***ID=%s sx=%s sy=%s sxy=%s angle=%s major=%s minor=%s" % (ID,sx,sy,sxy,angle,smj,smi)
            #self.obj.add(data)
            #x+=1
        #sys.exit('end of hyperQuad')

    def OES_CQUADR_82(self):  # not done...
        """
        GRID-ID  DISTANCE,NORMAL-X,NORMAL-Y,SHEAR-XY,ANGLE,MAJOR MINOR,VONMISES
        """
        if self.elementType == 82:  # CQUADR
            nTotal = 348  # 2+17*5 = 87 -> 87*4 = 348
            nNodes = 4    # centroid + 4 corner points
            eType = 'CQUADR'
        else:
            raise RuntimeError('elementType=%s nTotal not defined...'
                               % (self.elementType))

        dt = self.nonlinearFactor
        (format1, extract) = self.getOUG_FormatStart()
        format1 += '16f'  # 1+16 = 17
        format1 = bytes(format1)

        while len(self.data) >= nTotal:
            (eid, _) = unpack(b'i4s', self.data[0:8])
            self.data = self.data[8:]  # 2
            eid = extract(eid, dt)
            eData = self.data[0:68]  # 4*17
            self.data = self.data[68:]
            out = unpack(format1, eData)  # len=17*4
            (grid, fd1, sx1, sy1, txy1, angle1, major1, minor1, vm1,
             fd2, sx2, sy2, txy2, angle2, major2, minor2, vm2,) = out
            grid = 'C'
            self.obj.addNewEid(eType, eid, grid, fd1, sx1, sy1,
                               txy1, angle1, major1, minor1, vm1)
            self.obj.add(eid, grid, fd2, sx2, sy2, txy2,
                         angle2, major2, minor2, vm2)

            for nodeID in xrange(nNodes):  # nodes pts
                eData = self.data[0:68]
                self.data = self.data[68:]
                out = unpack(b'i16f', eData)
                (grid, fd1, sx1, sy1, txy1, angle1, major1, minor1, vm1,
                 fd2, sx2, sy2, txy2, angle2, major2, minor2, vm2,) = out

                #print "eid=%i grid=%i fd1=%i sx1=%i sy1=%i txy1=%i angle1=%i major1=%i minor1=%i vm1=%i" % (eid,grid,fd1,sx1,sy1,txy1,angle1,major1,minor1,vm1)
                #print "               fd2=%i sx2=%i sy2=%i txy2=%i angle2=%i major2=%i minor2=%i vm2=%i\n"          % (fd2,sx2,sy2,txy2,angle2,major2,minor2,vm2)
                #print "len(data) = ",len(self.data)
                #self.printBlock(self.data)
                self.obj.addNewNode(eid, grid, fd1, sx1,
                                    sy1, txy1, angle1, major1, minor1, vm1)
                self.obj.add(eid, grid, fd2, sx2, sy2,
                             txy2, angle2, major2, minor2, vm2)
            ###
        ###

    def OES_CQUAD4_144(self):  # works
        """
        GRID-ID  DISTANCE,NORMAL-X,NORMAL-Y,SHEAR-XY,ANGLE,MAJOR MINOR,VONMISES
        """
        #self.printSection(20)
        #term = data[0:4] CEN/
        #data = data[4:]
        #print "*****"
        #self.printBlock(self.data)
        #assert self.numWide==87,'invalid numWide...numWide=%s' % (self.numWide)
        #if self.numWide==87: # 2+(17-1)*5 = 87 -> 87*4 = 348

        if self.elementType == 144:  # CQUAD4
            nTotal = 348  # 2+17*5 = 87 -> 87*4 = 348
            nNodes = 4    # centroid + 4 corner points
            eType = 'CQUAD4'
        elif self.elementType == 64:  # CQUAD8
            nTotal = 348  # 2+17*5 = 87 -> 87*4 = 348
            nNodes = 4    # centroid + 4 corner points
            eType = 'CQUAD8'
        elif self.elementType == 82:  # CQUADR
            nTotal = 348  # 2+17*5 = 87 -> 87*4 = 348
            nNodes = 4    # centroid + 4 corner points
            eType = 'CQUAD4'  ## @todo write the word CQUADR

        elif self.elementType == 75:  # CTRIA6
            nTotal = 280  # 2+17*3 = 70 -> 70*4 = 280
            nNodes = 3    # centroid + 3 corner points
            eType = 'CTRIA6'
        elif self.elementType == 70:  # CTRIAR
            nTotal = 280  # 2+17*3 = 70 -> 70*4 = 280
            nNodes = 3    # centroid + 3 corner points
            eType = 'CTRIAR'  ## @todo write the word CTRIAR
        else:
            raise NotImplementedError('elementType=%s nTotal not defined...'
                                      % (self.elementType))

        dt = self.nonlinearFactor
        (format1, extract) = self.getOUG_FormatStart()
        format1 += '16f'
        format1 = bytes(format1)

        while len(self.data) >= nTotal:
            (eid, _) = unpack(b'i4s', self.data[0:8])
            self.data = self.data[8:]  # 2
            eid = extract(eid, dt)
            eData = self.data[0:68]
            self.data = self.data[68:]
            out = unpack(format1, eData)  # len=17*4
            (grid, fd1, sx1, sy1, txy1, angle1, major1, minor1, vm1,
             fd2, sx2, sy2, txy2, angle2, major2, minor2, vm2,) = out
            grid = 'C'
            self.obj.addNewEid(eType, dt, eid, grid, fd1, sx1, sy1,
                               txy1, angle1, major1, minor1, vm1)
            self.obj.add(dt, eid, grid, fd2, sx2, sy2, txy2,
                         angle2, major2, minor2, vm2)

            for nodeID in xrange(nNodes):  # nodes pts
                eData = self.data[0:68]
                self.data = self.data[68:]
                out = unpack(b'i16f', eData)
                (grid, fd1, sx1, sy1, txy1, angle1, major1, minor1, vm1,
                 fd2, sx2, sy2, txy2, angle2, major2, minor2, vm2,) = out

                #print "eid=%i grid=%i fd1=%i sx1=%i sy1=%i txy1=%i angle1=%i major1=%i minor1=%i vm1=%i" % (eid,grid,fd1,sx1,sy1,txy1,angle1,major1,minor1,vm1)
                #print "               fd2=%i sx2=%i sy2=%i txy2=%i angle2=%i major2=%i minor2=%i vm2=%i\n"          % (fd2,sx2,sy2,txy2,angle2,major2,minor2,vm2)
                #print "len(data) = ",len(self.data)
                #self.printBlock(self.data)
                self.obj.addNewNode(dt, eid, grid, fd1, sx1, sy1,
                                    txy1, angle1, major1, minor1, vm1)
                self.obj.add(dt, eid, grid, fd2, sx2, sy2,
                             txy2, angle2, major2, minor2, vm2)

            #print '--------------------'
            #print "len(data) = ",len(self.data)
            #print "tell = ",self.op2.tell()

            #self.printSection(100)
            #self.dn += 348
        ###
       # elif self.numWide==77:
       #     while len(self.data) >= 308: # 2+15*5 = 77 -> 77*4 = 308
       #         (eid,_) = unpack(b'i4s',self.data[0:8])
       #         self.data = self.data[8:]  # 2
       #         eData     = self.data[0:62] # 4*15
       #         self.data = self.data[62: ]
       #
       #         out = unpack(b'i14f', eData)  # 15
       #         (grid,fd1,sx1r,sx11,sy1r,sy11,txy1r,txy11,
       #               fd2,sx2r,sx21,sy2r,sy21,txy2r,txy21) = out
       #
       #         for nodeID in xrange(nNodes):   #nodes pts
       #             eData     = self.data[0:4*15]
       #             self.data = self.data[4*15: ]
       #             out = unpack(b'i14f', eData)
       #             (grid,fd1,sx1r,sx11,sy1r,sy11,txy1r,txy11,
       #                   fd2,sx2r,sx21,sy2r,sy21,txy2r,txy21) = out
       # elif self.numWide==47:
       #     while len(self.data) >= 188: # 2+9*5 = 47 -> 47*4 = 188
       #         (eid,_) = unpack(b'i4s',self.data[0:8])
       #         self.data = self.data[8:]  # 2
       #         eData     = self.data[0:36] # 4*9
       #         self.data = self.data[36: ]
       #
       #         out = unpack(b'i8f', eData)  # 9
       #         (grid,fd1,sx1,sy1,txy1,
       #               fd2,sx2,sy2,txy2) = out
       #
       #         for nodeID in xrange(nNodes):   #nodes pts
       #             eData     = self.data[0:36] # 4*9
       #             self.data = self.data[36: ]
       #             out = unpack(b'i8f', eData)
       #             (grid,fd1,sx1,sy1,txy1,
       #                   fd2,sx2,sy2,txy2) = out
       # else:
       #     raise NotImplementedError('invalid numWide')
