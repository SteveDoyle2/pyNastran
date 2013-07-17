#pylint: disable=C0103,C0301,R0914,E1101
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
from numpy import ones
from struct import unpack
#from pyNastran import isRelease

#91  -> PENTANL
#2   -> BEAM
#33  -> TUBE
#92  -> CONRODNL

# OES_CBEAM_2
#   CBEAM (2) linear

# OES_CSOLID_39_67_68
#   CTETRA (39) linear
#   CPENTA (67) linear
#   CHEXA  (68) linear

# OES_CBAR_34
#   CBAR (34) linear

# OES_CQUAD4_144
#   CQUAD4 (144) linear 5 nodes
#   CQUAD8 (82)  linear 5 nodes
#   CTRIA6 (75)  linear 4 nodes
#   CTRIAR (70)  linear 4 nodes


class RealElementsStressStrain2(object):

    def OES_CBEAM_2(self):
        dt = self.nonlinear_factor
        (formatStart, extract) = self.getOUG_FormatStart()

        nNodesPerBeam = 10  # 11-1
        ntotal = self.obj.getLengthTotal()
        (n1, format1) = self.obj.getLength1()
        (n2, format2) = self.obj.getLength2()
        format1 = formatStart + format1
        format1 = bytes(format1)
        format2 = bytes(format2)

        nelements = len(self.data) // ntotal
        #nnodes = nelements * (nNodesPerBeam + 1)
        nnodes = nelements * 2
        if self.read_mode == 0:  # figure out the shape
            self.obj._increase_size(dt, nelements, nnodes)
            iend = ntotal * nelements
            self.data = self.data[iend:]
            return
        else:  # read_mode = 1; # we know the shape so we can make a pandas matrix
            print("nelements =",nelements)
        (inode_start, inode_end, ielement_start, ielement_end
            ) = self.obj._preallocate(dt, nnodes, nelements)

        istart = 0
        iend = ntotal

        etypes=['CBEAM'] * nelements
        eids2 = []
        
        eids=[]; nids=[]; xxb=[]
        sxc=[]; sxd=[]; sxe=[]; sxf=[]
        s1=[]; s2=[]; s3=[]; s4=[]; axial=[]; smax=[]; smin=[]
        MSt=[]; MSc=[]

        istart = 0
        iend = n1
        #print('n1=%s n2=%s' % (n1, n2))
        for ielement in xrange(nelements):
            eData = self.data[istart:iend]

            (eType, eid, xxbi, s1i, s2i, s3i, s4i, smaxi, smini, MSti, MSci) = unpack(format1, eData)
            out = unpack(format1, eData)
            #print("outA = ",out)

            #eid2 = extract(eid, dt)

            eids2.append(eid)
            #print('eidA=%s eidB=%s' % (eid, eid2))

            nids.append(0)
            xxb.append(xxbi)
            eids.append(eid)
            s1.append(s1i)
            s2.append(s2i)
            s3.append(s3i)
            s4.append(s4i)
            smax.append(smaxi)
            smin.append(smini)
            MSt.append(MSti)
            MSc.append(MSci)
            
            istart = iend
            iend = istart + n2
            for inode in xrange(nNodesPerBeam-1):
                istart = iend
                iend = istart + n2
            
            #for inode in xrange(nNodesPerBeam):
            if 1:
                #nids.append(inode + 1)
                nids.append(1)
                eData = self.data[istart:iend]
                out = unpack(format2, eData)
                #print("outB = ",out)
                eid2, xxbi, s1i, s2i, s3i, s4i, smaxi, smini, MSti, MSci = out
                #print('eidB=%s' % eid)
                eids.append(eid)
                xxb.append(xxbi)
                s1.append(s1i)
                s2.append(s2i)
                s3.append(s3i)
                s4.append(s4i)
                smax.append(smaxi)
                smin.append(smini)
                MSt.append(MSti)
                MSc.append(MSci)

                istart = iend
                iend = istart + n2
            iend = istart + n1
        self.data = self.data[iend:]

        #print('delta', inode_end - inode_start, len(eids))
        if dt:
            name = self.obj.data_code['name']
            self.obj.data[name][inode_start:inode_end] = ones(inode_end - inode_start) * dt
        self.obj.data['element_id'][inode_start:inode_end] = eids
        self.obj.data['grid'][inode_start:inode_end] = nids
        self.obj.data['xxb'][inode_start:inode_end] = xxb

        #print("inode_start=%r inode_end=%r delta=%r" % (inode_start, inode_end, inode_end-inode_start))
        #print('len(s1) =', len(s1))

        if self.obj.isStress():
            self.obj.data['sxc'][inode_start:inode_end] = s1
            self.obj.data['sxd'][inode_start:inode_end] = s2
            self.obj.data['sxe'][inode_start:inode_end] = s3
            self.obj.data['sxf'][inode_start:inode_end] = s4

            self.obj.data['smax'][inode_start:inode_end] = smax
            self.obj.data['smin'][inode_start:inode_end] = smin
        else:
            print('type', self.obj.__class__.__name__)
            self.obj.data['e1'][inode_start:inode_end] = s1
            self.obj.data['e2'][inode_start:inode_end] = s2
            self.obj.data['e2'][inode_start:inode_end] = s3
            self.obj.data['e3'][inode_start:inode_end] = s4

            self.obj.data['emax'][inode_start:inode_end] = smax
            self.obj.data['emin'][inode_start:inode_end] = smin

        self.obj.data['MS_tension'][inode_start:inode_end] = MSt
        self.obj.data['MS_compression'][inode_start:inode_end] = MSc

        self.obj.element_data['element_id'][ielement_start:ielement_end] = eids2
        self.obj.element_data['element_type'][ielement_start:ielement_end] = etypes

        #print('len(eids) = ', len(self.obj.data['element_id']))
        #print('inode_end // ntotal = ', inode_end)
        if len(self.obj.data['element_id']) == inode_end: # [nodes, elements]
            self.obj._finalize(dt)
        else:
            self.obj._increment(nnodes, nelements)

    def OES_CSOLID_39_67_68(self, name):
        """
        stress is extracted at the centroid
        self.element_type in [39, 67, 68]:   # ctetra/chexa/cpenta (linear)
        CTETRA_39
        CPENTA_67
        CHEXA_68
        """
        #print("***********************************************")
        dt = self.nonlinear_factor
        (format1, extract) = self.getOUG_FormatStart()
        format1 += "i4si"
        format1 = bytes(format1)

        #nnodes = 5 # 1 centroid + 4 corner points
        #self.print_section(20)
        #term      = self.data[0:4] CEN/
        #self.data = self.data[4:]
        ElementType = self.get_element_type(self.element_type)
        #nnodes = 1  # this is a minimum, it will be reset later
        #nnodes_expected = 1
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
                   'num_wide=%s' % (ElementType, self.element_type,
                                   nNodes, self.num_wide))
            raise NotImplementedError(msg)

        ntotal = 16 + 84 * nnodes_expected
        nelements = len(self.data) // ntotal

        nnodes = nelements * nnodes_expected
        if self.read_mode == 0:  # figure out the shape
            self.obj._increase_size(dt, nelements, nnodes)
            iend = ntotal * nelements
            self.data = self.data[iend:]
            return
        else:  # read_mode = 1; # we know the shape so we can make a pandas matrix
            #index_data = (nodeIDs_to_index)
            #index_data = pd.MultiIndex.from_tuples(zip(data))

            #print("nelements =",nelements)
            pass
        (inode_start, inode_end, ielement_start, ielement_end
            ) = self.obj._preallocate(dt, nnodes, nelements)

        ibase = 0

        # element data
        cids=[]; etypes=[]; eids=[]; 

        element_ids=[]; nodes=[]; grids=[];
        oxx=[]; oyy=[]; ozz=[]; txy=[]; txz=[]; tyz=[]
        o1=[]; o2=[]; o3=[]; ovmShear=[]
        for i in xrange(nelements):
            out = unpack(format1, self.data[ibase:ibase+16])
            (eid, cid, abcd, nNodes) = out
            eid2 = extract(eid, dt)
            #print "abcd = |%s|" % (abcd)
            #print "eid=%s cid=%s nNodes=%s nNodesExpected=%s" % (eid,cid,nNodes,nNodesExpected)

            assert nNodes < 21, self.print_block(eData)


            etypes.append(ElementType)
            cids.append(cid)
            eids.append(eid2)

            ibase += 16
            for nodeID in xrange(nnodes_expected):  # nodes pts, +1 for centroid (???)
                out = unpack(b'i20f',self.data[ibase:ibase + 84]) # 4*21 = 84
                (grid_device, sxx, sxy, s1, a1, a2, a3, pressure, svm,
                              syy, syz, s2, b1, b2, b3,
                              szz, sxz, s3, c1, c2, c3) = out

                if grid_device == 0:
                    grid = 'C'
                else:
                    #grid = (grid_device - device_code) // 10
                    grid = grid_device

                element_ids.append(eid2)
                grids.append(grid_device)
                oxx.append(sxx)
                oyy.append(syy)
                ozz.append(szz)
                txy.append(sxy)
                txz.append(sxz)
                tyz.append(syz)
                ovmShear.append(svm)

                #print "%s eid=%s cid=%s grid=%s sxx=%-6i txy=%-5i s1=%-6i a1=%i a2=%i a3=%i press=%i vm=%s" % (element_type,eid,cid,grid,sxx,sxy,s1,a1,a2,a3,pressure,svm)
                #print "%s eid=%s cid=%s grid=%s syy=%-6i tyz=%-5i s2=%-6i b1=%i b2=%i b3=%i"                % (element_type,eid,cid,grid,syy,syz,s2,b1,b2,b3)
                #print "%s eid=%s cid=%s grid=%s szz=%-6i txz=%-5i s3=%-6i c1=%i c2=%i c3=%i"                % (element_type,eid,cid,grid,szz,sxz,s3,c1,c2,c3)
                #print ""
                #smax = max(s1,s2,s3)
                #smin = min(s1,s2,s3)

                aCos = []
                bCos = []
                cCos = []
                #if nodeID == 0:  # center point
                    #self.obj.add_new_eid(ElementType, cid, dt, eid2, grid, sxx, syy, szz, sxy, syz, sxz, s1, s2, s3, aCos, bCos, cCos, pressure, svm)
                #else:
                    #self.obj.add(dt, eid2, grid, sxx, syy, szz, sxy, syz, sxz, s1, s2, s3, aCos, bCos, cCos, pressure, svm)
                ibase += 84
        self.data = self.data[ibase:]

        if dt:
            name = self.obj.data_code['name']
            self.obj.data[name][inode_start:inode_end] = ones(inode_end - inode_start) * dt
        self.obj.data['element_id'][inode_start:inode_end] = element_ids
        self.obj.data['node_id'][inode_start:inode_end] = grids

        #print("inode_start=%r inode_end=%r delta=%r" % (inode_start, inode_end, inode_end-inode_start))
        #print('len(oxx) =', len(oxx))

        if self.obj.isStress():
            self.obj.data['oxx'][inode_start:inode_end] = oxx
            self.obj.data['oyy'][inode_start:inode_end] = oyy
            self.obj.data['ozz'][inode_start:inode_end] = ozz

            self.obj.data['txy'][inode_start:inode_end] = txy
            self.obj.data['txz'][inode_start:inode_end] = txz
            self.obj.data['tyz'][inode_start:inode_end] = tyz

            self.obj.data['o1'][inode_start:inode_end] = o1
            self.obj.data['o2'][inode_start:inode_end] = o2
            self.obj.data['o3'][inode_start:inode_end] = o3
            if self.obj.isVonMises():
                self.obj.data['ovm'][inode_start:inode_end] = ovmShear
            else:
                self.obj.data['max_shear'][inode_start:inode_end] = ovmShear
        else:
            print('type', self.obj.__class__.__name__)
            self.obj.data['exx'][inode_start:inode_end] = oxx
            self.obj.data['eyy'][inode_start:inode_end] = oyy
            self.obj.data['ezz'][inode_start:inode_end] = ozz

            self.obj.data['exy'][inode_start:inode_end] = txy
            self.obj.data['exz'][inode_start:inode_end] = txz
            self.obj.data['eyz'][inode_start:inode_end] = tyz

            self.obj.data['e1'][inode_start:inode_end] = o1
            self.obj.data['e2'][inode_start:inode_end] = o2
            self.obj.data['e3'][inode_start:inode_end] = o3
            if self.obj.isVonMises():
                self.obj.data['evm'][inode_start:inode_end] = ovmShear
            else:
                self.obj.data['max_shear'][inode_start:inode_end] = ovmShear

        self.obj.element_data['element_id'][ielement_start:ielement_end] = eids
        self.obj.element_data['element_type'][ielement_start:ielement_end] = etypes
        self.obj.element_data['cid'][ielement_start:ielement_end] = cids
        # pressure
        # aCos
        # bCos
        # cCos
        #self.obj.data[''][inode_start:inode_end] = translations[:, 5]

        if len(self.obj.data['element_id']) == inode_end // ntotal: # [nodes, elements]
            self.obj._finalize(dt)
        else:
            self.obj._increment(nnodes, nelements)

    def OES_CBAR_34(self):
        dt = self.nonlinear_factor
        assert self.num_wide == 16, ('invalid num_wide...num_wide=%s'
                                    % (self.num_wide))

        (format1, extract) = self.getOUG_FormatStart()
        format1 += '15f'
        format1 = bytes(format1)
        
        ntotal = 64
        nelements = len(self.data) // 64
        nnodes = nelements
        if self.read_mode == 0:  # figure out the shape
            self.obj._increase_size(dt, nelements, nnodes)
            iend = ntotal * nelements
            self.data = self.data[iend:]
            return
        else:  # read_mode = 1; # we know the shape so we can make a pandas matrix
            print("nelements =",nelements)
        (inode_start, inode_end, ielement_start, ielement_end
            ) = self.obj._preallocate(dt, nnodes, nelements)

        istart = 0
        iend = 64
        
        eids=[]; axial=[]; MSt=[]; MSc=[]
        s1a=[]; s2a=[]; s3a=[]; s4a=[]; smaxa=[]; smina=[]
        s1b=[]; s2b=[]; s3b=[]; s4b=[]; smaxb=[]; sminb=[]
        for ielement in xrange(nelements):
            eData = self.data[istart:iend]

            (eid, s1ai, s2ai, s3ai, s4ai, axiali, smaxai, sminai, MSti,
                  s1bi, s2bi, s3bi, s4bi, smaxbi, sminbi, MSci) = unpack(format1, eData)
            eid2 = extract(eid, dt)
            
            eids.append(eid2)
            s1a.append(s1ai)
            s2a.append(s2ai)
            s3a.append(s3ai)
            s4a.append(s4ai)
            smaxa.append(smaxai)
            smina.append(sminai)
            
            s1b.append(s1bi)
            s2b.append(s2bi)
            s3b.append(s3bi)
            s4b.append(s4bi)
            smaxb.append(smaxbi)
            sminb.append(sminbi)
            
            axial.append(axiali)
            MSt.append(MSti)
            MSc.append(MSci)

            #self.obj.add_new_eid('CBAR', dt, eid2, s1a, s2a, s3a, s4a, axial, smaxa, smina, MSt,
            #                   s1b, s2b, s3b, s4b, smaxb, sminb, MSc)

            #print "eid=%i s1=%i s2=%i s3=%i s4=%i axial=%-5i smax=%i smin=%i MSt=%i MSc=%i" % (eid,s1a,s2a,s3a,s4a,axial,smaxa,smina,MSt,MSc)
            #print "         s1=%i s2=%i s3=%i s4=%i          smax=%i smin=%i" % (s1b,s2b,s3b,s4b,smaxb,sminb)
            #print "len(data) = ",len(self.data)
            istart = iend
            iend += 64
        self.data = self.data[iend:]

        #print('delta', inode_end - inode_start, len(eids))
        if dt:
            name = self.obj.data_code['name']
            self.obj.data[name][inode_start:inode_end] = ones(inode_end - inode_start) * dt

        self.obj.data['element_id'][inode_start:inode_end] = eids

        #print("inode_start=%r inode_end=%r delta=%r" % (inode_start, inode_end, inode_end-inode_start))
        #print('len(s1) =', len(s1))

        if self.obj.isStress():
            self.obj.data['s1a'][inode_start:inode_end] = s1a
            self.obj.data['s2a'][inode_start:inode_end] = s2a
            self.obj.data['s3a'][inode_start:inode_end] = s3a
            self.obj.data['s4a'][inode_start:inode_end] = s4a
            self.obj.data['smaxa'][inode_start:inode_end] = smaxa
            self.obj.data['smina'][inode_start:inode_end] = smina

            self.obj.data['s1b'][inode_start:inode_end] = s1b
            self.obj.data['s2b'][inode_start:inode_end] = s2b
            self.obj.data['s3b'][inode_start:inode_end] = s3b
            self.obj.data['s4b'][inode_start:inode_end] = s4b
            self.obj.data['smaxb'][inode_start:inode_end] = smaxb
            self.obj.data['sminb'][inode_start:inode_end] = sminb
        else:
            #print('type', self.obj.__class__.__name__)
            self.obj.data['e1a'][inode_start:inode_end] = s1a
            self.obj.data['e2a'][inode_start:inode_end] = s2a
            self.obj.data['e3a'][inode_start:inode_end] = s3a
            self.obj.data['e4a'][inode_start:inode_end] = s4a
            self.obj.data['emaxa'][inode_start:inode_end] = smaxa
            self.obj.data['emina'][inode_start:inode_end] = smina

            self.obj.data['e1b'][inode_start:inode_end] = s1b
            self.obj.data['e2b'][inode_start:inode_end] = s2b
            self.obj.data['e3b'][inode_start:inode_end] = s3b
            self.obj.data['e4b'][inode_start:inode_end] = s4b
            self.obj.data['emaxb'][inode_start:inode_end] = smaxb
            self.obj.data['eminb'][inode_start:inode_end] = sminb

        self.obj.data['axial'][inode_start:inode_end] = axial
        self.obj.data['MS_tension'][inode_start:inode_end] = MSt
        self.obj.data['MS_compression'][inode_start:inode_end] = MSc

        #self.obj.data['element_id'][ielement_start:ielement_end] = eids
        #self.obj.data['element_type'][ielement_start:ielement_end] = etypes

        print('len(eids) = ', len(self.obj.data['element_id']))
        print('inode_end // ntotal = ', inode_end)
        if len(self.obj.data['element_id']) == inode_end: # [nodes, elements]
            self.obj._finalize(dt)
        else:
            self.obj._increment(nnodes, nelements)

    def OES_CQUAD4_144(self):
        """
        GRID-ID  DISTANCE,NORMAL-X,NORMAL-Y,SHEAR-XY,ANGLE,MAJOR MINOR,VONMISES
        """
        #term = data[0:4] CEN/
        #data = data[4:]
        #print "*****"
        #self.print_block(self.data)
        #assert self.num_wide==87,'invalid num_wide...num_wide=%s' % (self.num_wide)
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
            eType = 'CQUAD4'  #: .. todo:: write the word CQUADR

        elif self.element_type == 75:  # CTRIA6
            ntotal = 280  # 2+17*3 = 70 -> 70*4 = 280
            nNodes = 3    # centroid + 3 corner points
            eType = 'CTRIA6'
        elif self.element_type == 70:  # CTRIAR
            ntotal = 280  # 2+17*3 = 70 -> 70*4 = 280
            nNodes = 3    # centroid + 3 corner points
            eType = 'CTRIAR'  #: .. todo:: write the word CTRIAR
        else:
            raise NotImplementedError('element_type=%s ntotal not defined...'
                                      % (self.element_type))

        dt = self.nonlinear_factor
        (format1, extract) = self.getOUG_FormatStart()
        format1 += '16f'
        format1 = bytes(format1)
        
        nelements = len(self.data) // ntotal
        nrows = (1 + nNodes) * nelements * 2 # 2 layers
        #print('nrows = %i' % nrows)
        
        gridC = 'C'
        cformat = b'i4s'+format1  # center format
        nformat = b'i16f'         # node format

        nnodes = nrows #nelements * nNodes
        if self.read_mode == 0:  # figure out the shape
            self.obj._increase_size(dt, nelements, nnodes)
            iend = ntotal * nelements
            self.data = self.data[iend:]
            return
        else:  # read_mode = 1; # we know the shape so we can make a pandas matrix
            print("nelements =",nelements)
        
        (inode_start, inode_end, ielement_start, ielement_end
            ) = self.obj._preallocate(dt, nnodes, nelements)


        eids=[]; eTypes=[]; eids2=[]; nids=[]; layer=[]
        fd=[]; sx=[]; sy=[]; txy=[]; angle=[]; major=[]; minor=[]; vm=[]
        #['fd1', 'sx1', 'sy1', 'txy1', 'angle1', 'major1', 'minor1', 'vm1',
        # 'fd2', 'sx2', 'sy2', 'txy2', 'angle2', 'major2', 'minor2', 'vm2',]

        ibase = 0
        for i in xrange(nelements):
            eData = self.data[ibase:ibase + 76]  # 8 + 68

            (eid, _, grid, fd1i, sx1i, sy1i, txy1i, angle1i, major1i, minor1i, vm1i,
                           fd2i, sx2i, sy2i, txy2i, angle2i, major2i, minor2i, vm2i,) = unpack(cformat, eData)  # len=17*4
            #gridC = 'C'
            eid = extract(eid, dt)
            eTypes.append(eType)
            eids.append(eid)

            eids2.append(eid)
            eids2.append(eid)

            nids.append(0)
            nids.append(0)

            layer.append(1)
            layer.append(2)

            fd.append(fd1i)
            sx.append(sx1i)
            sy.append(sy1i)
            txy.append(txy1i)
            angle.append(angle1i)
            major.append(major1i)
            minor.append(minor1i)
            vm.append(vm1i)

            fd.append(fd2i)
            sx.append(sx2i)
            sy.append(sy2i)
            txy.append(txy2i)
            angle.append(angle2i)
            major.append(major2i)
            minor.append(minor2i)
            vm.append(vm2i)

            #self.obj.add_new_eid(eType, dt, eid, gridC, fd1i, sx1i, sy1i, txy1i, angle1i, major1i, minor1i, vm1i)
            #self.obj.add(               dt, eid, gridC, fd2i, sx2i, sy2i, txy2i, angle2i, major2i, minor2i, vm2i)

            ibase += 76
            for node_id in xrange(nNodes):  # nodes pts
                out = unpack(nformat, self.data[ibase:ibase + 68])
                (grid, fd1i, sx1i, sy1i, txy1i, angle1i, major1i, minor1i, vm1i,
                       fd2i, sx2i, sy2i, txy2i, angle2i, major2i, minor2i, vm2i,) = out

                eids2.append(eid)
                eids2.append(eid)

                nids.append(grid)
                nids.append(grid)

                layer.append(1)
                layer.append(2)

                fd.append(fd1i)
                sx.append(sx1i)
                sy.append(sy1i)
                txy.append(txy1i)
                angle.append(angle1i)
                major.append(major1i)
                minor.append(minor1i)
                vm.append(vm1i)
                
                fd.append(fd2i)
                sx.append(sx2i)
                sy.append(sy2i)
                txy.append(txy2i)
                angle.append(angle2i)
                major.append(major2i)
                minor.append(minor2i)
                vm.append(vm2i)
                #print "eid=%i grid=%i fd1=%i sx1=%i sy1=%i txy1=%i angle1=%i major1=%i minor1=%i vm1=%i" % (eid,grid,fd1,sx1,sy1,txy1,angle1,major1,minor1,vm1)
                #print "               fd2=%i sx2=%i sy2=%i txy2=%i angle2=%i major2=%i minor2=%i vm2=%i\n"          % (fd2,sx2,sy2,txy2,angle2,major2,minor2,vm2)
                #self.obj.addNewNode(dt, eidi, gridi, fd1i, sx1i, sy1i, txy1i, angle1i, major1i, minor1i, vm1i)
                #self.obj.add(       dt, eidi, gridi, fd2i, sx2i, sy2i, txy2i, angle2i, major2i, minor2i, vm2i)
                ibase += 68
        self.data = self.data[ibase:]

        print('delta', inode_end - inode_start, len(eids))
        if dt:
            name = self.obj.data_code['name']
            self.obj.data[name][inode_start:inode_end] = ones(inode_end - inode_start) * dt

        self.obj.data['element_id'][inode_start:inode_end] = eids

        #print("inode_start=%r inode_end=%r delta=%r" % (inode_start, inode_end, inode_end-inode_start))
        #print('len(s1) =', len(s1))


        #self.obj.element_data['element_id'][ielement_start:ielement_end] = eids
        #self.obj.element_data['element_type'][ielement_start:ielement_end] = eTypes

        #print('nids',nids)
        self.obj.data['element_id'][inode_start:inode_end] = eids2
        self.obj.data['node_id'][inode_start:inode_end] = nids
        self.obj.data['layer'][inode_start:inode_end] = layer

        if self.obj.isStress():
            self.obj.data['oxx'][inode_start:inode_end] = sx
            self.obj.data['oyy'][inode_start:inode_end] = sy
            self.obj.data['txy'][inode_start:inode_end] = txy
            self.obj.data['omax'][inode_start:inode_end] = major
            self.obj.data['omin'][inode_start:inode_end] = minor
            if self.obj.isVonMises():
                key = 'ovm'
            else:
                key = 'max_shear'
        else:
            print('type', self.obj.__class__.__name__)
            self.obj.data['exx'][inode_start:inode_end] = sx
            self.obj.data['eyy'][inode_start:inode_end] = sy
            self.obj.data['exy'][inode_start:inode_end] = txy
            self.obj.data['emax'][inode_start:inode_end] = major
            self.obj.data['emin'][inode_start:inode_end] = minor
            if self.obj.isVonMises():
                key = 'evm'
            else:
                key = 'max_shear'

        self.obj.data['angle'][inode_start:inode_end] = angle
        self.obj.data[key][inode_start:inode_end] = vm

        #self.obj.data['element_id'][ielement_start:ielement_end] = eids
        #self.obj.data['element_type'][ielement_start:ielement_end] = etypes

        #print('len(eids) = ', len(self.obj.data['element_id']))
        #print('inode_end // ntotal = ', inode_end)
        if len(self.obj.data['element_id']) == inode_end: # [nodes, elements]
            self.obj._finalize(dt)
        else:
            self.obj._increment(nnodes, nelements)
