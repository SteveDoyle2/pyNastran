#pylint: disable=C0103,C0301,R0914,E1101
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
from struct import unpack
from pyNastran import isRelease

#91  -> PENTANL
#2   -> BEAM
#33  -> TUBE
#92  -> CONRODNL


class RealElementsStressStrain(object):

    def skipOES_Element(self):
        if not isRelease:
           raise NotImplementedError(self.code_information())

        self.log.debug('skipping approach_code=%s, table_code=%s, format_code-%s '
                       'sort_code=%s on %s table' % (self.analysis_code,
                       self.table_code, self.format_code, self.sort_code,
                       self.table_name))
        #print(self.code_information())
        #print("**************skipping**************")
        self.handle_results_buffer(self.dummyPass, None, debug=True)

    def dummyPass(self):
        self.data = b''

    def OES_field1(self):
        if self.is_sort1():
            #raise NotImplementedError('SORT1 is not supported')
            return ('i', self.scaleEid)
        elif self.isSort2():
            if self.analysis_code in [1, 2, 3, 4, 7, 8, 9, 11, 12]:  # eid
                return ('i', self.scaleEid)
            elif self.analysis_code == 5:  # freq
                #freq
                return ('f', self.scaleDt)
                #raise NotImplementedError('freq is not supported')
            elif self.analysis_code == 6:  # time
                #time
                return ('f', self.scaleDt)
                #raise NotImplementedError('time is not supported')
            elif self.analysis_code == 10:  # fqts:
                #fqts # freqTime
                return ('f', self.scaleDt)
                #raise NotImplementedError('freqTime is not supported')
            else:
                raise NotImplementedError('invalid SORT2 analysis_code=%s' %
                                          self.analysis_code)
        else:
            raise NotImplementedError('invalid SORTx code')

    def OES_Thermal(self, debug=False):
        #assert self.num_wide==5,'invalid num_wide...num_wide=%s' % self.num_wide

        dt = self.nonlinear_factor
        (format1, extract) = self.getOUG_FormatStart()
        format1 += '2i5f'  # 1+2+5=8
        format1 = bytes(format1)
        while len(self.data) >= 32:  # 4*8
            #print self.print_section(40)
            eData = self.data[0:32]
            self.data = self.data[32:]
            #print "len(data) = ",len(eData)

            out = unpack(format1, eData)
            if self.make_op2_debug:
                self.op2_debug.write('Thermal - %s\n' % (str(out)))

            (eid, sideID, hbdyID, cnvCoeff, fApplied, fConv, fRad, fTotal) = out
            eid = extract(eid, dt)
            #print('eid=%s sideID=%s hbdyID=%s coeff=%s fApplied=%s fConv=%s '
            #      'fRad=%s fTotal=%s' % (eid,sideID,hbdyID,cnvCoeff,fApplied,
            #                             fConv,fRad,fTotal))
            #self.obj.add_new_eid(eid,axial,axialMS,torsion,torsionMS)

            #print "eid=%i axial=%i torsion=%i" % (eid,axial,torsion)
            #print "len(data) = ",len(self.data)
        #print self.rodStress[self.isubcase]

    def OES_basicElement(self):
        """
        genericStressReader - works on CROD_1, CELAS2_12
        stress & strain
        format_code=1 sort_code=0 (eid,axial,axialMS,torsion,torsionMS)
        """
        dt = self.nonlinear_factor
        (format1, extract) = self.getOUG_FormatStart()
        (nTotal, dataFormat) = self.obj.getLength()
        dataFormat = format1 + dataFormat
        #print "nTotal=%s dataFormat=%s len(data)=%s" % (nTotal,dataFormat,len(self.data))
        dataFormat = bytes(dataFormat)

        nelements = len(self.data) // nTotal
        n = 0
        for i in xrange(nelements):
            eData = self.data[n:n + nTotal]
            out = unpack(dataFormat, eData)
            if self.make_op2_debug:
                self.op2_debug.write('basicElement - %s\n' % (str(out)))
            eid = extract(out[0], dt)
            self.obj.add_new_eid(dt, eid, out[1:])
            n += nTotal
        self.data = self.data[n:]

    def OES_CBEAM_2(self):
        dt = self.nonlinear_factor
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
            if self.make_op2_debug:
                self.op2_debug.write('CBEAM-2 - %s\n' % (str(out)))
            eid2 = extract(out[0], dt)
            self.obj.add_new_eid(dt, eid2, out[1:])

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
        dt = self.nonlinear_factor
        (format1, extract) = self.getOUG_FormatStart()
        format1 += '16f'
        format1 = bytes(format1)

        nNodes = 0  # centroid + 4 corner points
        assert self.num_wide == 17, ('invalid num_wide...num_wide=%s' % self.num_wide)

        while len(self.data) >= 68:  # 2+17*5 = 87 -> 87*4 = 348
            #print self.print_block(self.data[0:100])
            #(eid,) = unpack(b'i',self.data[0:4])
            #self.data = self.data[8:]  # 2
            eData = self.data[0:68]  # 4*17
            self.data = self.data[68:]
            out = unpack(format1, eData)  # 17
            if self.make_op2_debug:
                self.op2_debug.write('CQUAD4-33A - %s\n' % (str(out)))

            (eid, fd1, sx1, sy1, txy1, angle1, major1, minor1, maxShear1,
             fd2, sx2, sy2, txy2, angle2, major2, minor2, maxShear2) = out

            eid = extract(eid, dt)

            #print "eid=%i grid=%s fd1=%-3.1f sx1=%i sy1=%i txy1=%i angle1=%i major1=%i minor1=%i vm1=%i" % (eid,'C',fd1,sx1,sy1,txy1,angle1,major1,minor1,maxShear1)
            #print   "             fd2=%-3.1f sx2=%i sy2=%i txy2=%i angle2=%i major2=%i minor2=%i vm2=%i\n"       % (fd2,sx2,sy2,txy2,angle2,major2,minor2,maxShear2)
            #print "nNodes = ",nNodes
            self.obj.add_new_eid('CQUAD4', dt, eid, 'C', fd1, sx1, sy1,
                               txy1, angle1, major1, minor1, maxShear1)
            self.obj.add(dt, eid, 'C', fd2, sx2, sy2, txy2,
                         angle2, major2, minor2, maxShear2)

            for nodeID in xrange(nNodes):  # nodes pts
                eData = self.data[0:68]  # 4*17
                self.data = self.data[68:]
                out = unpack(b'i16f', eData[0:68])
                if self.make_op2_debug:
                    self.op2_debug.write('CQUAD4-33B - %s\n' % (str(out)))
                (grid, fd1, sx1, sy1, txy1, angle1, major1, minor1, vm1,
                 fd2, sx2, sy2, txy2, angle2, major2, minor2, vm2,) = out

                #print "eid=%i grid=%i fd1=%i sx1=%i sy1=%i txy1=%i angle1=%i major1=%i minor1=%i vm1=%i" % (eid,grid,fd1,sx1,sy1,txy1,angle1,major1,minor1,vm1)
                #print "               fd2=%i sx2=%i sy2=%i txy2=%i angle2=%i major2=%i minor2=%i vm2=%i\n"          % (fd2,sx2,sy2,txy2,angle2,major2,minor2,vm2)
                #print "len(data) = ",len(self.data)
                #self.print_block(self.data)
                self.obj.addNewNode(dt, eid, grid, fd1, sx1, sy1,
                                    txy1, angle1, major1, minor1, vm1)
                self.obj.add(dt, eid, grid, fd2, sx2, sy2,
                             txy2, angle2, major2, minor2, vm2)
            #print '--------------------'
            #print "len(data) = ",len(self.data)
            #print "tell = ",self.op2.tell()

            #self.print_section(100)
            #self.dn += 348

    def OES_CBAR_34(self):
        dt = self.nonlinear_factor
        assert self.num_wide == 16, ('invalid num_wide...num_wide=%s' % self.num_wide)

        (format1, extract) = self.getOUG_FormatStart()
        format1 += '15f'
        format1 = bytes(format1)

        while len(self.data) >= 64:
            eData = self.data[0:64]
            self.data = self.data[64:]

            out = unpack(format1, eData)
            if self.make_op2_debug:
                self.op2_debug.write('CBAR-34 - %s\n' % (str(out)))
            (eid, s1a, s2a, s3a, s4a, axial, smaxa, smina, MSt,
                         s1b, s2b, s3b, s4b, smaxb, sminb, MSc) = out
            eid2 = extract(eid, dt)
            self.obj.add_new_eid('CBAR', dt, eid2, s1a, s2a, s3a, s4a, axial, smaxa, smina, MSt,
                               s1b, s2b, s3b, s4b, smaxb, sminb, MSc)

            #print "eid=%i s1=%i s2=%i s3=%i s4=%i axial=%-5i smax=%i smin=%i MSt=%i MSc=%i" % (eid,s1a,s2a,s3a,s4a,axial,smaxa,smina,MSt,MSc)
            #print "         s1=%i s2=%i s3=%i s4=%i          smax=%i smin=%i" % (s1b,s2b,s3b,s4b,smaxb,sminb)

    def OES_CBUSH1D_40(self):
        dt = self.nonlinear_factor
        (format1, extract) = self.getOUG_FormatStart()

        assert self.num_wide == 8, "num_wide=%s not 8" % self.num_wide
        nTotal = 32  # 4*8
        format1 += '6fi'
        format1 = bytes(format1)

        n = 0
        nEntries = len(self.data) // nTotal
        for i in xrange(nEntries):
            eData = self.data[n:n + nTotal]

            out = unpack(format1, eData)  # num_wide=25
            if self.make_op2_debug:
                self.op2_debug.write('CBUSH1D-40 - %s\n' % (str(out)))
            (eid, fe, ue, ve, ao, ae, ep, fail) = out
            eid = extract(eid, dt)

            # axial_force, axial_displacement, axial_velocity, axial_stress, axial_strain, plastic_strain, is_failed
            self.obj.add_new_eid(self.element_type, dt, eid, fe, ue, ve, ao, ae, ep, fail)
            n += nTotal
        self.data = self.data[n:]

    def OES_CTRIAX6_53(self):
        #(Format1,scaleValue) = self.OES_field1()
        #Format = Format1+'ifffffff'
        dt = self.nonlinear_factor
        (format1, extract) = self.getOUG_FormatStart()
        format1 += 'i7f'
        format1 = bytes(format1)

        nelements = len(self.data) // 132  # (1+8*4)*4 = 33*4 = 132
        ibase = 0
        for i in xrange(nelements):
            out = unpack(format1, self.data[ibase:ibase + 36])
            (eid, loc, rs, azs, As, ss, maxp, tmax, octs) = out
            if self.make_op2_debug:
                self.op2_debugop2_debug.write('CTRIAX6-53A - %s\n' % (str(out)))
            eid = extract(eid, dt)
            #print "eid=%s loc=%s rs=%s azs=%s as=%s ss=%s maxp=%s tmx=%s octs=%s" % (eid,loc,rs,azs,As,ss,maxp,tmax,octs)
            self.obj.add_new_eid(dt, eid, loc, rs, azs, As, ss, maxp, tmax, octs)

            for i in xrange(3):
                out = unpack(b'i7f', self.data[ibase:ibase + 32])
                (loc, rs, azs, As, ss, maxp, tmax, octs) = out
                if self.make_op2_debug:
                    self.op2_debugop2_debug.write('CTRIAX6-53B - %s\n' % (str(out)))
                #print "eid=%s loc=%s rs=%s azs=%s as=%s ss=%s maxp=%s tmx=%s octs=%s" % (eid,loc,rs,azs,As,ss,maxp,tmax,octs)
                self.obj.add(dt, eid, loc, rs, azs, As, ss, maxp, tmax, octs)
                ibase += 32  # 4*8
        n = nelements * 132
        self.data = self.data[n:]

    def OES_CSOLID_39_67_68(self):
        """
        stress is extracted at the centroid
        self.element_type in [39, 67, 68]:   # ctetra/chexa/cpenta (linear)
        CTETRA_39
        CPENTA_67
        CHEXA_68
        """
        dt = self.nonlinear_factor
        (format1, extract) = self.getOUG_FormatStart()
        format1 += "i4si"
        format1 = bytes(format1)

        #nNodes = 5 # 1 centroid + 4 corner points
        element_type = self.get_element_type(self.element_type)
        nNodes = 1  # this is a minimum, it will be reset later
        nnodes_expected = 1

        #assert self.num_wide in [109,151,193],'invalid num_wide...num_wide=%s' % self.num_wide

        # overly complicated way to get the element type
        # so we can figure out how many elements there are
        if self.element_type == 39: # CTETRA
            nnodes_expected = 5
        elif self.element_type == 67: # CPENTA
            nnodes_expected = 7
        elif self.element_type == 68: # CHEXA
            nnodes_expected = 9
        else:
            msg = ('not supported....EType=%s eType=%s nNodes=%s'
                   'num_wide=%s' % (element_type, self.element_type,
                                   nNodes, self.num_wide))
            raise NotImplementedError(msg)

        ntotal = 16 + 84 * nnodes_expected
        nelements = len(self.data) // ntotal

        ibase = 0
        for i in xrange(nelements):
            out = unpack(format1, self.data[ibase:ibase+16])
            if self.make_op2_debug:
                self.op2_debug.write('%s-%sA - %s\n' % (element_type, self.element_type, str(out)))
            (eid, cid, abcd, nNodes) = out
            eid2 = extract(eid, dt)
            #print "abcd = %r" % abcd
            #print "eid=%s cid=%s nNodes=%s nNodesExpected=%s" % (eid,cid,nNodes,nNodesExpected)

            assert nNodes < 21,  self.print_block(self.data[ibase:ibase+16])

            ibase += 16
            for nodeID in xrange(nnodes_expected):  # nodes pts, +1 for centroid (???)
                out = unpack(b'i20f',self.data[ibase:ibase + 84]) # 4*21 = 84
                if self.make_op2_debug:
                    self.op2_debug.write('%s-%sB - %s\n' % (element_type, self.element_type, str(out)))
                (grid_device, sxx, sxy, s1, a1, a2, a3, pressure, svm,
                              syy, syz, s2, b1, b2, b3,
                              szz, sxz, s3, c1, c2, c3) = out

                if grid_device == 0:
                    grid = 'C'
                else:
                    #grid = (grid_device - device_code) // 10
                    grid = grid_device

                #print("%s eid=%s cid=%s grid=%s sxx=%-6i txy=%-5i s1=%-6i a1=%i a2=%i a3=%i press=%i vm=%s" % (element_type,eid,cid,grid,sxx,sxy,s1,a1,a2,a3,pressure,svm))
                #print("%s eid=%s cid=%s grid=%s syy=%-6i tyz=%-5i s2=%-6i b1=%i b2=%i b3=%i"                % (element_type,eid,cid,grid,syy,syz,s2,b1,b2,b3))
                #print("%s eid=%s cid=%s grid=%s szz=%-6i txz=%-5i s3=%-6i c1=%i c2=%i c3=%i"                % (element_type,eid,cid,grid,szz,sxz,s3,c1,c2,c3))
                #print("")

                #smax = max(s1,s2,s3)
                #smin = min(s1,s2,s3)

                aCos = []
                bCos = []
                cCos = []
                if nodeID == 0:
                    self.obj.add_new_eid(element_type, cid, dt, eid2, grid, sxx, syy, szz, sxy, syz, sxz, s1, s2, s3, aCos, bCos, cCos, pressure, svm)
                else:
                    self.obj.add(dt, eid2, grid, sxx, syy, szz, sxy, syz, sxz, s1, s2, s3, aCos, bCos, cCos, pressure, svm)
                ibase += 84
        n = nelements * ntotal
        assert n == ibase, 'n=%s ibase=%s' % (n, ibase)
        self.data = self.data[n:]

    def OES_CTRIA3_74(self):
        """
        DISTANCE,NORMAL-X,NORMAL-Y,SHEAR-XY,ANGLE,MAJOR,MINOR,VONMISES
        stress is extracted at the centroid
        """
        assert self.num_wide == 17, 'invalid num_wide...num_wide=%s' % self.num_wide

        dt = self.nonlinear_factor
        (format1, extract) = self.getOUG_FormatStart()
        format1 += '16f'
        format1 = bytes(format1)

        ntotal = 68  # 4*17
        n = 0
        nelements = len(self.data) // ntotal
        for i in xrange(nelements):
            eData = self.data[n:n + ntotal]
            out = unpack(format1, eData)
            if self.make_op2_debug:
                self.op2_debug.write('CTRIA3-74 - %s\n' % (str(out)))

            (eid, fd1, sx1, sy1, txy1, angle1, major1, minor1, vm1,
             fd2, sx2, sy2, txy2, angle2, major2, minor2, vm2,) = out
            eid = extract(eid, dt)
            #print "eid=%i fd1=%i sx1=%i sy1=%i txy1=%i angle1=%i major1=%i minor1=%i vm1=%i" % (eid,fd1,sx1,sy1,txy1,angle1,major1,minor1,vm1)
            #print  "      fd2=%i sx2=%i sy2=%i txy2=%i angle2=%i major2=%i minor2=%i vm2=%i\n"   % (fd2,sx2,sy2,txy2,angle2,major2,minor2,vm2)
            self.obj.add_new_eid('CTRIA3', dt, eid, 'C', fd1, sx1, sy1,
                               txy1, angle1, major1, minor1, vm1)
            self.obj.add(dt, eid, 'C', fd2, sx2, sy2, txy2,
                         angle2, major2, minor2, vm2)
            n += ntotal
        self.data = self.data[n:]

    def OES_CQUADR_82(self):
        """
        GRID-ID  DISTANCE,NORMAL-X,NORMAL-Y,SHEAR-XY,ANGLE,MAJOR MINOR,VONMISES
        """
        if self.element_type == 82:  # CQUADR
            nTotal = 348  # 2+17*5 = 87 -> 87*4 = 348
            nNodes = 4    # centroid + 4 corner points
            eType = 'CQUADR'
        else:
            raise RuntimeError('element_type=%s nTotal not defined...'
                               % (self.element_type))

        dt = self.nonlinear_factor
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
            if self.make_op2_debug:
                self.op2_debug.write('CQUADR-82A - %s\n' % str(out))
            (grid, fd1, sx1, sy1, txy1, angle1, major1, minor1, vm1,
             fd2, sx2, sy2, txy2, angle2, major2, minor2, vm2,) = out
            grid = 'C'
            self.obj.add_new_eid(eType, eid, grid, fd1, sx1, sy1,
                               txy1, angle1, major1, minor1, vm1)
            self.obj.add(eid, grid, fd2, sx2, sy2, txy2,
                         angle2, major2, minor2, vm2)

            for nodeID in xrange(nNodes):  # nodes pts
                eData = self.data[0:68]
                self.data = self.data[68:]
                out = unpack(b'i16f', eData)
                if self.make_op2_debug:
                    self.op2_debug.write('CQUADR-82B - %s\n' % str(out))
                (grid, fd1, sx1, sy1, txy1, angle1, major1, minor1, vm1,
                 fd2, sx2, sy2, txy2, angle2, major2, minor2, vm2,) = out

                #print "eid=%i grid=%i fd1=%i sx1=%i sy1=%i txy1=%i angle1=%i major1=%i minor1=%i vm1=%i" % (eid,grid,fd1,sx1,sy1,txy1,angle1,major1,minor1,vm1)
                #print "               fd2=%i sx2=%i sy2=%i txy2=%i angle2=%i major2=%i minor2=%i vm2=%i\n"          % (fd2,sx2,sy2,txy2,angle2,major2,minor2,vm2)
                #print "len(data) = ",len(self.data)
                #self.print_block(self.data)
                self.obj.addNewNode(eid, grid, fd1, sx1,
                                    sy1, txy1, angle1, major1, minor1, vm1)
                self.obj.add(eid, grid, fd2, sx2, sy2,
                             txy2, angle2, major2, minor2, vm2)

    def OES_CGAPNL_86(self):
        dt = self.nonlinear_factor
        (format1, extract) = self.getOUG_FormatStart()

        assert self.num_wide == 11, "num_wide=%s not 11" % (self.num_wide)
        ntotal = 44  # 4*11
        format1 += '8f4s4s'
        format1 = bytes(format1)

        n = 0
        nelements = len(self.data) // ntotal
        for i in xrange(nelements):
            eData = self.data[n:n + ntotal]

            out = unpack(format1, eData)  # num_wide=25
            if self.make_op2_debug:
                self.op2_debug.write('CGAPNL-86 - %s\n' % str(out))
            (eid, cpx, shy, shz, au, shv, shw, slv, slp, form1, form2) = out
            eid = extract(eid, dt)

            self.obj.add_new_eid(dt, eid, cpx, shy, shz, au, shv, shw, slv, slp, form1, form2)
            n += ntotal
        n = nelements * ntotal
        self.data = self.data[n:]

    def OES_RODNL_89_92(self):
        dt = self.nonlinear_factor
        (format1, extract) = self.getOUG_FormatStart()
        format1 += '6f'  # 1+6=7
        format1 = bytes(format1)

        #ibase = 0
        nelements = len(self.data) // 28  # len(format1)*4 = 7*4 = 28
        for i in xrange(nelements):
            eData = self.data[0:28]
            self.data = self.data[28:]
            out = unpack(format1, eData)
            if self.make_op2_debug:
                self.op2_debug.write('CRODNL-%s - %s\n' % (self.element_type, str(out)))

            (eid, axial, equivStress, totalStrain, effPlasticCreepStrain,
                effCreepStrain, linearTorsionalStresss) = out
            eid = extract(eid, dt)
            data = (eid, axial, equivStress, totalStrain, effPlasticCreepStrain, effCreepStrain, linearTorsionalStresss)

            #print "eid=%s axial=%s equivStress=%s totalStrain=%s effPlasticCreepStrain=%s effCreepStrain=%s linearTorsionalStresss=%s" % (eid,axial,equivStress,totalStrain,effPlasticCreepStrain,effCreepStrain,linearTorsionalStresss)
            self.obj.add(self.element_type, dt, data)
        #self.data = self.data[ibase:]

    def OES_CQUAD4NL_90(self):
        dt = self.nonlinear_factor
        (format1, extract) = self.getOUG_FormatStart()

        ntotal = 52  # 4*13
        format1 += '12f'  # 1+12=13
        format1 = bytes(format1)
        assert 13 == self.num_wide, 'num_wide=%s not 13' % self.num_wide

        n = 0
        nelements = len(self.data) // ntotal
        for i in xrange(nelements):
            eData = self.data[n:n + ntotal]
            out = unpack(format1, eData)  # num_wide=13
            if self.make_op2_debug:
                self.op2_debug.write('CQUADNL-90 - %s\n' % str(out))

            (eid, fd1, sx1, sy1, sz1, txy1, es1, eps1, ecs1,
                ex1, ey1, ez1, exy1) = out
            eid = extract(eid, dt)
            data = (eid, fd1, sx1, sy1, sz1, txy1, es1, eps1,
                    ecs1, ex1, ey1, ez1, exy1)
            self.obj.add_new_eid(self.element_type, dt, data)
            #print "eid=%s axial=%s equivStress=%s totalStrain=%s effPlasticCreepStrain=%s effCreepStrain=%s linearTorsionalStresss=%s" % (eid,axial,equivStress,totalStrain,effPlasticCreepStrain,effCreepStrain,linearTorsionalStresss)
            n += ntotal
        n = nelements * ntotal
        self.data = self.data[n:]

    def OES_TETRANL_85_PENTANL_91_CHEXANL_93(self):  # TETRANL 85 / PENTANL 91 / HEXANL 93
        """
        The DMAP manual says fields 3-18 repeat 9 times. but they dont.
        They repeat 8 times.  Other DMAP cards are correct with
        their repeat statements.

        stress is extracted at the centroid
        CTETRANL_ - 4 nodes
        CPENTANL_91 - 7 nodes
        CHEXANL_93 - 9 nodes
        """
        dt = self.nonlinear_factor
        (format1, extract) = self.getOUG_FormatStart()
        format1 += '4s'
        format1 = bytes(format1)

        if self.element_type == 85:
            eType = 'CTETRANL'
            nNodes = 5
        elif self.element_type == 91:
            eType = 'CPENTANL'
            nNodes = 7
            #nTotal = 456
        elif self.element_type == 93:
            eType = 'CHEXANL'
            nNodes = 9
            #nTotal = 584
        else:
            raise NotImplementedError(self.element_type)
        ntotal = 8 + 64 * nNodes

        nelements = len(self.data) // ntotal
        for i in xrange(nelements):  # 2+16*9 = 146 -> 146*4 = 584
            eData = self.data[0:8]
            self.data = self.data[8:]
            out = unpack(format1, eData)
            if self.make_op2_debug:
                self.op2_debug.write('%s-%s - %s\n' % (eType, self.element_type, str(out)))
            (eid, cType) = out
            eid = extract(eid, dt)

            for i in xrange(nNodes):
                eData = self.data[0:64]
                self.data = self.data[64:]
                out = unpack(b'i15f', eData)
                if self.make_op2_debug:
                    self.op2_debug.write('%s-%sB - %s\n' % (eType, self.element_type, str(out)))

                assert len(out) == 16
                (grid, sx, sy, sz, sxy, syz, sxz, se, eps,
                  ecs, ex, ey, ez, exy, eyz, exz) = out
                #print "eid=%3s cType=%s sx=%i sy=%i sz=%i sxy=%s syz=%i szx=%i se=%s" % (eid,cType,sx,sy,sz,sxy,syz,sxz,se)
                #print "gid=%3s ecs=%.3g   ex=%.3g ey=%.3g ez=%.3g exy=%.3g eyz=%.3g ezx=%.3g"  % (grid,ecs,ex,ey,ez,exy,eyz,exz)
                #assert cType == 'GRID',cType
                #print("cType =", cType)
        #n = nelements * ntotal
        #self.data = self.data[n:]

    def OES_VUHEXA_145_VUPENTA_146_VUTETRA_147(self):  # VUHEXA 145 /
        """
        VUTETRA - 4 nodes
        VUPENTA 146 - 6 nodes
        VUHEXA 145 - 8 nodes
        """
        dt = self.nonlinear_factor
        (format1, extract) = self.getOUG_FormatStart()
        format1 += 'i'
        format1 = bytes(format1)

        if self.element_type == 147:
            eType = 'VUTETRA'
            nNodes = 4
        elif self.element_type == 146:
            eType = 'VUPENTA'
            nNodes = 6
        elif self.element_type == 145:
            eType = 'VUHEXA'
            nNodes = 8
        else:
            raise NotImplementedError(self.element_type)

        num_wide = 2 + 12 * nNodes
        nTotal = 8 + 48 * nNodes
        assert self.num_wide == num_wide

        while len(self.data) >= nTotal:  # 2+16*9 = 146 -> 146*4 = 584
            eData = self.data[0:8]
            self.data = self.data[8:]
            (eid, parentID) = unpack(format1, eData)
            eid = extract(eid, dt)

            for i in xrange(nNodes):
                eData = self.data[0:48]
                self.data = self.data[48:]
                out = unpack(b'i11f', eData)
                if self.make_op2_debug:
                    self.op2_debug.write('%s-%s - %s\n' % (eType, self.element_type, str(out)))
                assert len(out) == 12
                (grid, xnorm, ynorm, znorm, txy, tyz, txz,
                 prin1, prin2, prin3, mean, vonoRoct) = out
                #print "eid=%3s cType=%s sx=%i sy=%i sz=%i sxy=%s syz=%i szx=%i se=%s" % (eid,cType,sx,sy,sz,sxy,syz,sxz,se)
                #print "gid=%3s ecs=%.3g   ex=%.3g ey=%.3g ez=%.3g exy=%.3g eyz=%.3g ezx=%.3g"  % (grid,ecs,ex,ey,ez,exy,eyz,exz)
                #assert cType == 'GRID',cType
                #print("grid =",grid)

    def OES_CBEAM_94(self):
        #nNodes = 10  # 11-1

        #nTotal       = self.obj.getLengthTotal()
        #(n1,format1) = self.obj.getLength1()
        #(n2,format2) = self.obj.getLength2()
        #ntotal = 2 * 4 + (18 - 3) * 9 * 4
        nTotal = 204

        n1 = 24
        format1 = '4s5f'
        format1 = bytes(format1)

        nelements = len(self.data) // ntotal
        for i in xrange(nelements):
            eData = self.data[0:8]
            self.data = self.data[8:]
            out = unpack(b'2i', eData)
            (eid, gridA) = out
            if self.make_op2_debug:
                self.op2_debug.write('CBEAM-94A - %s\n' % (str(out)))
            #print "eid=%s gridA=%s" % (eid,gridA)

            for i in xrange(1):
                for j in xrange(4):  # c,d,e,f @ A;    c,d,e,f @ B
                    eData = self.data[0:n1]
                    self.data = self.data[n1:]
                    out = unpack(format1, eData)
                    if self.make_op2_debug:
                        self.op2_debug.write('CBEAM-94B - %s\n' % (str(out)))
                    (loc, nsx, nse, te, epe, ece) = out
                    #print "loc=%s nsx=%s nse=%s te=%s epe=%s ece=%s" % (loc,nsx,nse,te,epe,ece)
                #self.obj.add(eid,out)

    def OES_CQUAD4_95(self):  # works (doesnt handle all stress/strain cases tho)
        """
        GRID-ID  DISTANCE,NORMAL-X,NORMAL-Y,SHEAR-XY,ANGLE,MAJOR MINOR,VONMISES
        composite quad

         95 - CQUAD4
         96 - CQUAD8
         97 - CTRIA3
         98 - CTRIA6 (composite)
        """
        eType = self.get_element_type(self.element_type)

        #self.print_section(20)
        #term = data[0:4] CEN/
        #data = data[4:]
        if self.num_wide != 11:
            raise RuntimeError('invalid num_wide; num_wide=%s' % self.num_wide)

        dt = self.nonlinear_factor
        (format1, extract) = self.getOUG_FormatStart()
        format1 += 'i9f'
        format1 = bytes(format1)

        nelements = len(self.data) // 44  # 2+17*5 = 87 -> 87*4 = 348
        ibase = 0
        for i in xrange(nelements):
        #while len(self.data) <= 44:
            eData = self.data[ibase:ibase+44]  # 4*11
            out = unpack(format1, eData)
            if self.make_op2_debug:
                self.op2_debug.write('CQUAD4-95 - %s\n' % (str(out)))
            (eid, iLayer, o1, o2, t12, t1z, t2z, angle, major, minor, ovm) = out
            eid = extract(eid, dt)

            if eid != self.eid2:  # originally initialized to None, the buffer doesnt reset it, so it is the old value
                #print "1 - eid=%s iLayer=%i o1=%i o2=%i ovm=%i" % (eid,iLayer,o1,o2,ovm)
                self.obj.add_new_eid(eType, dt, eid, o1, o2, t12, t1z, t2z, angle, major, minor, ovm)
            else:
                #print "2 - eid=%s iLayer=%i o1=%i o2=%i ovm=%i" % (eid,iLayer,o1,o2,ovm)
                self.obj.add(dt, eid, o1, o2, t12, t1z, t2z, angle, major, minor, ovm)
            self.eid2 = eid
            ibase += 44
            #self.dn += 348
        n = nelements * 44
        self.data = self.data[n:]
        #print "3 - eid=%s iLayer=%i o1=%i o2=%i ovm=%i" % (eid,iLayer,o1,o2,ovm)

    def OES_QUAD4FD_139(self):  # hyperelastic
        """
        Hyperelastic Quad
        36+4*7*4 = 148
        """
        dt = self.nonlinear_factor
        (format1, extract) = self.getOUG_FormatStart()
        format1 += '4si6f'  # 1 + 4+1+6 = 12
        format1 = bytes(format1)

        while len(self.data) >= 148:
            #if x==2:
            #    sys.exit('end of hyperQuad')
            eData = self.data[0:36]  # 4*9
            self.data = self.data[36:]
            out = unpack(format1, eData)
            if self.make_op2_debug:
                self.op2_debug.write('CQUAD4FD-139A - %s\n' % (str(out)))

            (eid, Type, ID, sx, sy, sxy, angle, smj, smi) = out
            eid = extract(eid, dt)
            self.obj.add_new_eid(dt, [eid, Type, sx, sy, sxy, angle, smj, smi])
            #print "eid=%s Type=%s\n***ID=%s sx=%s sy=%s sxy=%s angle=%s major=%s minor=%s" % (eid,Type,ID,sx,sy,sxy,angle,smj,smi)
            for i in xrange(3):
                eData = self.data[0:28]  # 4*7
                self.data = self.data[28:]
                out = unpack(b'i6f', eData)
                if self.make_op2_debug:
                    self.op2_debug.write('CQUAD4FD-139B - %s\n' % (str(out)))
                #(ID,sx,sy,sxy,angle,smj,smi) = out
                self.obj.add(dt, eid, out)
                #print "***ID=%s sx=%s sy=%s sxy=%s angle=%s major=%s minor=%s" % (ID,sx,sy,sxy,angle,smj,smi)
            #self.obj.add(data)
            #x+=1
        #n = nelements * 148
        #self.data = self.data[n:]

    def OES_CBUSH_102(self):
        dt = self.nonlinear_factor
        (format1, extract) = self.getOUG_FormatStart()

        assert self.num_wide == 7, "num_wide=%s not 7" % self.num_wide
        ntotal = 28  # 4*7
        format1 += '6f'
        format1 = bytes(format1)

        n = 0
        nelements = len(self.data) // ntotal
        for i in xrange(nelements):
            eData = self.data[n:n + ntotal]
            out = unpack(format1, eData)  # num_wide=7
            if self.make_op2_debug:
                self.op2_debug.write('CBUSH-102 - %s\n' % str(out))

            (eid, tx, ty, tz, rx, ry, rz) = out
            eid = extract(eid, dt)

            self.obj.add_new_eid(self.element_type, dt, eid, tx, ty, tz, rx, ry, rz)
            n += ntotal
        n = nelements * ntotal
        self.data = self.data[n:]

    def OES_CQUAD4_144(self):
        """
        GRID-ID  DISTANCE,NORMAL-X,NORMAL-Y,SHEAR-XY,ANGLE,MAJOR MINOR,VONMISES
        """
        #term = data[0:4] CEN/
        #data = data[4:]
        #print "*****"
        #self.print_block(self.data)
        #assert self.num_wide==87,'invalid num_wide...num_wide=%s' % self.num_wide
        #if self.num_wide==87: # 2+(17-1)*5 = 87 -> 87*4 = 348

        if self.element_type == 144:  # CQUAD4
            ntotal = 348  # 2+17*5 = 87 -> 87*4 = 348
            nNodes = 4    # centroid + 4 corner points
            eType = 'CQUAD4'
        elif self.element_type == 64:  # CQUAD8
            ntotal = 348  # 2+17*5 = 87 -> 87*4 = 348
            nNodes = 4    # centroid + 4 corner points
            eType = 'CQUAD8'
        elif self.element_type == 82:  # CQUADR
            ntotal = 348  # 2+17*5 = 87 -> 87*4 = 348
            nNodes = 4    # centroid + 4 corner points
            eType = 'CQUAD4'  # TODO write the word CQUADR

        elif self.element_type == 75:  # CTRIA6
            ntotal = 280  # 2+17*3 = 70 -> 70*4 = 280
            nNodes = 3    # centroid + 3 corner points
            eType = 'CTRIA6'
        elif self.element_type == 70:  # CTRIAR
            ntotal = 280  # 2+17*3 = 70 -> 70*4 = 280
            nNodes = 3    # centroid + 3 corner points
            eType = 'CTRIAR'  # TODO write the word CTRIAR
        else:
            raise NotImplementedError('element_type=%s nTotal not defined...'
                                      % (self.element_type))

        dt = self.nonlinear_factor
        (format1, extract) = self.getOUG_FormatStart()
        format1 += '16f'
        format1 = bytes(format1)

        nelements = len(self.data) // ntotal
        nrows = (1 + nNodes) * nelements  # centroid + nnodes
        #print('nrows = %i' % nrows)

        ibase = 0
        gridC = 'C'
        cformat = b'i4s'+format1  # center format
        nformat = b'i16f'         # node format
        for i in xrange(nelements):
            eData = self.data[ibase:ibase + 76]  # 8 + 68

            out = unpack(cformat, eData)  # len=17*4
            if self.make_op2_debug:
                self.op2_debug.write('%s-%sA - %s\n' % (eType, self.element_type, str(out)))
            (eid, _, grid, fd1, sx1, sy1, txy1, angle1, major1, minor1, vm1,
                           fd2, sx2, sy2, txy2, angle2, major2, minor2, vm2,) = out
            #gridC = 'C'
            eid = extract(eid, dt)
            self.obj.add_new_eid(eType, dt, eid, gridC, fd1, sx1, sy1,
                                 txy1, angle1, major1, minor1, vm1)
            self.obj.add(dt, eid, gridC, fd2, sx2, sy2, txy2,
                         angle2, major2, minor2, vm2)

            ibase += 76
            for nodeID in xrange(nNodes):  # nodes pts
                out = unpack(nformat, self.data[ibase:ibase + 68])
                if self.make_op2_debug:
                    self.op2_debug.write('%s-%sB - %s\n' % (eType, self.element_type, str(out)))
                (grid, fd1, sx1, sy1, txy1, angle1, major1, minor1, vm1,
                 fd2, sx2, sy2, txy2, angle2, major2, minor2, vm2,) = out

                #print "eid=%i grid=%i fd1=%i sx1=%i sy1=%i txy1=%i angle1=%i major1=%i minor1=%i vm1=%i" % (eid,grid,fd1,sx1,sy1,txy1,angle1,major1,minor1,vm1)
                #print "               fd2=%i sx2=%i sy2=%i txy2=%i angle2=%i major2=%i minor2=%i vm2=%i\n"          % (fd2,sx2,sy2,txy2,angle2,major2,minor2,vm2)
                self.obj.addNewNode(dt, eid, grid, fd1, sx1, sy1,
                                    txy1, angle1, major1, minor1, vm1)
                self.obj.add(dt, eid, grid, fd2, sx2, sy2,
                             txy2, angle2, major2, minor2, vm2)
                ibase += 68
        n = nelements * ntotal
        self.data = self.data[n:]

    def OES_VUQUAD_189(self):
        if self.element_type == 144:  # CQUAD4
            nTotal = 440  # 6+(33-7)*4 =  -> 110*4 = 440
            nNodes = 4    # 4 corner points
            eType = 'CQUAD4'
        #elif self.element_type == 64:  # CQUAD8
            #nTotal = 348  # 2+17*5 = 87 -> 87*4 = 348
            #nNodes = 4    # centroid + 4 corner points
            #eType = 'CQUAD8'
        #elif self.element_type == 82:  # CQUADR
            #nTotal = 348  # 2+17*5 = 87 -> 87*4 = 348
            #nNodes = 4    # centroid + 4 corner points
            #eType = 'CQUAD4'  # TODO write the word CQUADR
        #elif self.element_type == 75:  # CTRIA6
            #nTotal = 280  # 2+17*3 = 70 -> 70*4 = 280
            #nNodes = 3    # centroid + 3 corner points
            #eType = 'CTRIA6'
        #elif self.element_type == 70:  # CTRIAR
            #nTotal = 280  # 2+17*3 = 70 -> 70*4 = 280
            #nNodes = 3    # centroid + 3 corner points
            #eType = 'CTRIAR'  # TODO write the word CTRIAR
        else:
            raise NotImplementedError('element_type=%s nTotal not defined...'
                                      % (self.element_type))

        dt = self.nonlinear_factor
        (format1, extract) = self.getOUG_FormatStart()
        format1 += '2i4s2i'
        format1 = bytes(format1)

        nelements = len(self.data) // ntotal
        for i in xrange(nelements):
            (eid, parent, coord, icord, theta, itype) = unpack(b'i4s',
                                                               self.data[0:8])
            self.data = self.data[8:]  # 2
            eid = extract(eid, dt)
            eData = self.data[0:68]
            self.data = self.data[68:]

            ### todo...this is never unpacked...
            out = unpack(format1, eData)  # len=17*4
            if self.make_op2_debug:
                self.op2_debug.write('VU%s-%s - %s\n' % (eType, self.element_type, str(out)))

            self.obj.addNewNode(dt, eid, parent, coord, icord, theta, itype)

            self.obj.add_new_eid(eType, dt, eid, parent, coord, icord, theta, itype)
            for nodeID in xrange(nNodes):  # nodes pts
                eData = self.data[0:68]
                self.data = self.data[68:]
                out = unpack(b'i16f', eData)
                if self.make_op2_debug:
                    self.op2_debug.write('%s\n' % (str(out)))

                (vuid, dummy, dummy2, msx, msy, mxy, dummy3, dummy4, dummy5,
                 bcx, bcy, bcxy,tyz,tzx,dummy6,dummy7,dummy8) = out
                self.obj.add(vuid, dummy, dummy2, msx, msy, mxy,
                             dummy3, dummy4, dummy5,
                             bcx, bcy, bcxy,tyz,tzx,
                             dummy6,dummy7,dummy8)
        n = nelements * ntotal
        #self.data = self.data[n:]

    def OES_CELAS_224_225(self):
        dt = self.nonlinear_factor
        element_name = self.data_code['element_name']
        (format1, extract) = self.getOUG_FormatStart()

        assert self.num_wide == 3, "num_wide=%s not 3" % self.num_wide
        ntotal = 12  # 4*3
        format1 += '2f'
        format1 = bytes(format1)

        nelements = len(self.data) // ntotal
        for i in xrange(nelements):
            eData = self.data[0:ntotal]
            self.data = self.data[ntotal:]

            out = unpack(format1, eData)  # num_wide=3
            if self.make_op2_debug:
                self.op2_debug.write('%s-%s - %s\n' % (element_name, self.element_type, str(out)))
            (eid, force, stress) = out
            eid = extract(eid, dt)
            self.obj.add_new_eid(element_name, dt, eid, force, stress)
        #n = nelements * ntotal
        #self.data = self.data[n:]

    def OESRT_CQUAD4_95(self):
        #dt = self.nonlinear_factor
        #element_name = self.data_code['element_name']
        (format1, extract) = self.getOUG_FormatStart()

        assert self.num_wide == 9, "num_wide=%s not 9" % self.num_wide
        nTotal = 36  # 4*9
        format1 += '8si3fi4s'
        format1 = bytes(format1)

        n = 0
        nelements = len(self.data) // ntotal
        for i in xrange(nelements):
            eData = self.data[n:n + nTotal]
            out = unpack(format1, eData)  # num_wide=9
            if self.make_op2_debug:
                self.op2_debug.write('CQUAD4-95 - %s\n' % str(out))
            #print(self.print_block(eData[-4:]))
            #asfd
            #eid,failure, ply, failureIndexPly, failureIndexBonding, failureIndexMax, flag
            # 3,TSAIWU,1,8.5640,0.0,None
            #print('out', out)

            (eid, failure, ply, strengthRatioPly, failureIndexBonding, strengthRatioBonding, flag, flag2) = out
            strengthRatioPly
            #print("eid=%s failure=|%s| ply=%s failureIndexPly=%s  failureIndexBonding=%s strengthRatioBonding=%s flag=%s flag2=%s" % (eid, failure.strip(), ply, failureIndexPly, failureIndexBonding, strengthRatioBonding, flag, flag2))
            print("eid=%s strengthRatioPly=%g failureIndexBonding=%s strengthRatioBonding=%s" % (eid, strengthRatioPly, failureIndexBonding, strengthRatioBonding))

            #self.obj.add_new_eid(element_name, dt, eid, force, stress)
            n += ntotal
        n = nelements * ntotal
        self.data = self.data[n:]
        asd