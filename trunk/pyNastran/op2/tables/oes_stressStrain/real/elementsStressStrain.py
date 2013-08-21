#pylint: disable=C0103,C0301,R0914,E1101
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
from numpy import ones
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
        self.handle_results_buffer(self.dummyPass, None, None, debug=True)

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
                                          (self.analysis_code))
        else:
            raise NotImplementedError('invalid SORTx code')

    def OES_Thermal(self, debug=False):
        #assert self.num_wide==5,'invalid num_wide...num_wide=%s' % (self.num_wide)

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
            (eid, sideID, hbdyID, cnvCoeff, fApplied, fConv,
                fRad, fTotal) = out
            eid = extract(eid, dt)
            #print('eid=%s sideID=%s hbdyID=%s coeff=%s fApplied=%s fConv=%s '
            #      'fRad=%s fTotal=%s' % (eid,sideID,hbdyID,cnvCoeff,fApplied,
            #                             fConv,fRad,fTotal))
            #self.obj.add_new_eid(eid,axial,axialMS,torsion,torsionMS)

            #print "eid=%i axial=%i torsion=%i" % (eid,axial,torsion)
            #print "len(data) = ",len(self.data)
        #print self.rodStress[self.isubcase]

    def OES_CBUSH1D_40(self, name):
        dt = self.nonlinear_factor
        (format1, extract) = self.getOUG_FormatStart()

        assert self.num_wide == 8, "num_wide=%s not 8" % (self.num_wide)
        ntotal = 32  # 4*8
        format1 += '6fi'
        format1 = bytes(format1)

        n = 0
        nelements = len(self.data) // ntotal
        if self.read_mode == 0 or name not in self._selected_names:
            if name not in self._result_names:
                self._result_names.append(name)
            # figure out the shape
            self.obj._increase_size(dt, nelements, nnodes)
            iend = ntotal * nelements
            self.data = self.data[iend:]
            return
        else:
            oes_cbush1d_40__elementStressStrain

        for i in xrange(nelements):
            eData = self.data[n:n + ntotal]

            out = unpack(format1, eData)  # num_wide=25
            (eid, fe, ue, ve, ao, ae, ep, fail) = out
            eid = extract(eid, dt)

            # axial_force, axial_displacement, axial_velocity, axial_stress, axial_strain, plastic_strain, is_failed
            self.obj.add_new_eid(self.element_type, dt, eid, fe, ue, ve, ao, ae, ep, fail)
            n += ntotal
        self.data = self.data[n:]

    def OES_CTRIAX6_53(self, name):
        #(Format1,scaleValue) = self.OES_field1()
        #Format = Format1+'ifffffff'
        dt = self.nonlinear_factor
        (format1, extract) = self.getOUG_FormatStart()
        format1 += 'i7f'
        format1 = bytes(format1)

        nelements = len(self.data) // 132  # (1+8*4)*4 = 33*4 = 132
        ibase = 0
        for i in xrange(nelements):
            (eid, loc, rs, azs, As, ss, maxp, tmax, octs) = unpack(format1, self.data[ibase:ibase + 36])
            eid = extract(eid, dt)
            #print "eid=%s loc=%s rs=%s azs=%s as=%s ss=%s maxp=%s tmx=%s octs=%s" % (eid,loc,rs,azs,As,ss,maxp,tmax,octs)
            self.obj.add_new_eid(dt, eid, loc, rs, azs, As, ss, maxp, tmax, octs)

            ibase += 36  # 4*9
            for i in xrange(3):
                (loc, rs, azs, As, ss, maxp, tmax, octs) = unpack(b'i7f', self.data[ibase:ibase + 32])
                #print "eid=%s loc=%s rs=%s azs=%s as=%s ss=%s maxp=%s tmx=%s octs=%s" % (eid,loc,rs,azs,As,ss,maxp,tmax,octs)
                self.obj.add(dt, eid, loc, rs, azs, As, ss, maxp, tmax, octs)
                ibase += 32  # 4*8
        self.data = self.data[ibase:]

    def OES_CQUADR_82(self, name):
        """
        GRID-ID  DISTANCE,NORMAL-X,NORMAL-Y,SHEAR-XY,ANGLE,MAJOR MINOR,VONMISES
        """
        if self.element_type == 82:  # CQUADR
            ntotal = 348  # 2+17*5 = 87 -> 87*4 = 348
            nNodes = 4    # centroid + 4 corner points
            eType = 'CQUADR'
        else:
            raise RuntimeError('element_type=%s ntotal not defined...'
                               % (self.element_type))

        dt = self.nonlinear_factor
        (format1, extract) = self.getOUG_FormatStart()
        format1 += '16f'  # 1+16 = 17
        format1 = bytes(format1)

        while len(self.data) >= ntotal:
            (eid, _) = unpack(b'i4s', self.data[0:8])
            self.data = self.data[8:]  # 2
            eid = extract(eid, dt)
            eData = self.data[0:68]  # 4*17
            self.data = self.data[68:]
            out = unpack(format1, eData)  # len=17*4
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

    def OES_CGAPNL_86(self, name):
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
            (eid, cpx, shy, shz, au, shv, shw, slv, slp, form1, form2) = out
            eid = extract(eid, dt)

            self.obj.add_new_eid(dt, eid, cpx, shy, shz, au, shv, shw, slv, slp, form1, form2)
            n += ntotal
        self.data = self.data[n:]

    def OES_RODNL_89_92(self, name):
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

            (eid, axial, equivStress, totalStrain, effPlasticCreepStrain,
                effCreepStrain, linearTorsionalStresss) = out
            eid = extract(eid, dt)
            data = (eid, axial, equivStress, totalStrain, effPlasticCreepStrain, effCreepStrain, linearTorsionalStresss)

            #print "eid=%s axial=%s equivStress=%s totalStrain=%s effPlasticCreepStrain=%s effCreepStrain=%s linearTorsionalStresss=%s" % (eid,axial,equivStress,totalStrain,effPlasticCreepStrain,effCreepStrain,linearTorsionalStresss)
            self.obj.add(self.element_type, dt, data)
        #self.data = self.data[ibase:]

    def OES_CQUAD4NL_90(self, name):
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

            (eid, fd1, sx1, sy1, sz1, txy1, es1, eps1, ecs1,
                ex1, ey1, ez1, exy1) = out
            eid = extract(eid, dt)
            data = (eid, fd1, sx1, sy1, sz1, txy1, es1, eps1,
                    ecs1, ex1, ey1, ez1, exy1)
            self.obj.add_new_eid(self.element_type, dt, data)
            #print "eid=%s axial=%s equivStress=%s totalStrain=%s effPlasticCreepStrain=%s effCreepStrain=%s linearTorsionalStresss=%s" % (eid,axial,equivStress,totalStrain,effPlasticCreepStrain,effCreepStrain,linearTorsionalStresss)
            n += ntotal
        self.data = self.data[n:]

    def OES_TETRANL_85_PENTANL_91_CHEXANL_93(self, name):  # TETRANL 85 / PENTANL 91 / HEXANL 93
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
            #ntotal = 456
        elif self.element_type == 93:
            eType = 'CHEXANL'
            nNodes = 9
            #ntotal = 584
        else:
            raise NotImplementedError(self.element_type)
        ntotal = 8 + 64 * nNodes

        nelements = len(self.data) // ntotal
        for i in range(nelements):  # 2+16*9 = 146 -> 146*4 = 584
            eData = self.data[0:8]
            self.data = self.data[8:]
            (eid, cType) = unpack(format1, eData)
            eid = extract(eid, dt)

            for i in xrange(nNodes):
                eData = self.data[0:64]
                self.data = self.data[64:]
                out = unpack(b'i15f', eData)
                assert len(out) == 16
                (grid, sx, sy, sz, sxy, syz, sxz, se, eps,
                  ecs, ex, ey, ez, exy, eyz, exz) = out
                #print "eid=%3s cType=%s sx=%i sy=%i sz=%i sxy=%s syz=%i szx=%i se=%s" % (eid,cType,sx,sy,sz,sxy,syz,sxz,se)
                #print "gid=%3s ecs=%.3g   ex=%.3g ey=%.3g ez=%.3g exy=%.3g eyz=%.3g ezx=%.3g"  % (grid,ecs,ex,ey,ez,exy,eyz,exz)
                #assert cType == 'GRID',cType
                #print("cType =", cType)

    def OES_VUHEXA_145_VUPENTA_146_VUTETRA_147(self, name):  # VUHEXA 145 /
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
        ntotal = 8 + 48 * nNodes
        assert self.num_wide == num_wide

        while len(self.data) >= ntotal:  # 2+16*9 = 146 -> 146*4 = 584
            eData = self.data[0:8]
            self.data = self.data[8:]
            (eid, parentID) = unpack(format1, eData)
            eid = extract(eid, dt)

            for i in xrange(nNodes):
                eData = self.data[0:48]
                self.data = self.data[48:]
                out = unpack(b'i11f', eData)
                assert len(out) == 12
                (grid, xnorm, ynorm, znorm, txy, tyz, txz,
                 prin1, prin2, prin3, mean, vonoRoct) = out
                #print "eid=%3s cType=%s sx=%i sy=%i sz=%i sxy=%s syz=%i szx=%i se=%s" % (eid,cType,sx,sy,sz,sxy,syz,sxz,se)
                #print "gid=%3s ecs=%.3g   ex=%.3g ey=%.3g ez=%.3g exy=%.3g eyz=%.3g ezx=%.3g"  % (grid,ecs,ex,ey,ez,exy,eyz,exz)
                #assert cType == 'GRID',cType
                #print("grid =",grid)

    def OES_CBEAM_94(self, name):
        #nNodes = 10  # 11-1

        #ntotal       = self.obj.getLengthTotal()
        #(n1,format1) = self.obj.getLength1()
        #(n2,format2) = self.obj.getLength2()
        #ntotal = 2 * 4 + (18 - 3) * 9 * 4
        ntotal = 204

        n1 = 24
        format1 = '4s5f'
        format1 = bytes(format1)

        nelements = len(self.data) // ntotal
        for i in xrange(nelements):
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

    def OES_QUAD4FD_139(self, name):  # hyperelastic
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

            (eid, Type, ID, sx, sy, sxy, angle, smj, smi) = out
            eid = extract(eid, dt)
            self.obj.add_new_eid(dt, [eid, Type, sx, sy, sxy, angle, smj, smi])
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

    def OES_VUQUAD_189(self, name):
        if self.element_type == 144:  # CQUAD4
            ntotal = 440  # 6+(33-7)*4 =  -> 110*4 = 440
            nNodes = 4    # 4 corner points
            eType = 'CQUAD4'
        #elif self.element_type == 64:  # CQUAD8
            #ntotal = 348  # 2+17*5 = 87 -> 87*4 = 348
            #nNodes = 4    # centroid + 4 corner points
            #eType = 'CQUAD8'
        #elif self.element_type == 82:  # CQUADR
            #ntotal = 348  # 2+17*5 = 87 -> 87*4 = 348
            #nNodes = 4    # centroid + 4 corner points
            #eType = 'CQUAD4'  # TODO write the word CQUADR
        #elif self.element_type == 75:  # CTRIA6
            #ntotal = 280  # 2+17*3 = 70 -> 70*4 = 280
            #nNodes = 3    # centroid + 3 corner points
            #eType = 'CTRIA6'
        #elif self.element_type == 70:  # CTRIAR
            #ntotal = 280  # 2+17*3 = 70 -> 70*4 = 280
            #nNodes = 3    # centroid + 3 corner points
            #eType = 'CTRIAR'  # TODO write the word CTRIAR
        else:
            raise NotImplementedError('element_type=%s ntotal not defined...'
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
            out = unpack(format1, eData)  # len=17*4
            self.obj.addNewNode(dt, eid, parent, coord, icord, theta, itype)

            self.obj.add_new_eid(eType, dt, eid, parent, coord, icord, theta, itype)
            for nodeID in xrange(nNodes):  # nodes pts
                eData = self.data[0:68]
                self.data = self.data[68:]
                out = unpack(b'i16f', eData)
                (vuid, dummy, dummy2, msx, msy, mxy, dummy3, dummy4, dummy5,
                 bcx, bcy, bcxy,tyz,tzx,dummy6,dummy7,dummy8) = out
                self.obj.add(vuid, dummy, dummy2, msx, msy, mxy,
                             dummy3, dummy4, dummy5,
                             bcx, bcy, bcxy,tyz,tzx,
                             dummy6,dummy7,dummy8)
        #self.data = self.data[ibase:]

    def OES_CELAS_224_225(self, name):
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
            (eid, force, stress) = out
            eid = extract(eid, dt)
            self.obj.add_new_eid(element_name, dt, eid, force, stress)

    def OESRT_CQUAD4_95(self, name):
        #dt = self.nonlinear_factor
        #element_name = self.data_code['element_name']
        (format1, extract) = self.getOUG_FormatStart()

        assert self.num_wide == 9, "num_wide=%s not 9" % self.num_wide
        ntotal = 36  # 4*9
        format1 += '8si3fi4s'
        format1 = bytes(format1)

        n = 0
        nelements = len(self.data) // ntotal
        for i in xrange(nelements):
            eData = self.data[n:n + ntotal]
            out = unpack(format1, eData)  # num_wide=9
            #print(self.print_block(eData[-4:]))
            #asfd
            #eid,failure, ply, failureIndexPly, failureIndexBonding, failureIndexMax, flag
            # 3,TSAIWU,1,8.5640,0.0,None
            #print('out', out)

            (eid, failure, ply, strengthRatioPly, failureIndexBonding, strengthRatioBonding, flag, flag2) = out
            strengthRatioPly
            #print("eid=%s failure=|%s| ply=%s failureIndexPly=%s  failureIndexBonding=%s strengthRatioBonding=%s flag=%s flag2=%s" % (eid, failure.strip(), ply, failureIndexPly, failureIndexBonding, strengthRatioBonding, flag, flag2))
            #print("eid=%s strengthRatioPly=%g failureIndexBonding=%s strengthRatioBonding=%s" % (eid, strengthRatioPly, failureIndexBonding, strengthRatioBonding))

            #self.obj.add_new_eid(element_name, dt, eid, force, stress)
            n += ntotal
        self.data = self.data[n:]
        asd