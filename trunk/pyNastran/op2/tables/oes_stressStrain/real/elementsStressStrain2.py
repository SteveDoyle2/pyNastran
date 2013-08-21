"""
 OES_basicElement - TODO
   CROD   (1)   linear (centroid)
   CTUBE  (3)   linear (centroid)
   CONROD (10)  linear (centroid)
   CELAS1 (11)  linear (centroid)
   CELAS2 (12)  linear (centroid)
   CELAS3 (13)  linear (centroid)
   CELAS4 (???) linear (centroid)

   CDAMP1 (???) linear (centroid)
   CDAMP2 (???) linear (centroid)
   CDAMP3 (???) linear (centroid)
   CDAMP4 (???) linear (centroid)

 OES_CBEAM_2 - beamStress / beamStrain
   CBEAM (2) linear

 OES_CSHEAR_4 - shearStress / shearStrain - TODO
   CSHEAR (4) linear

 OES_CSOLID_39_67_68 - solidStress / solidStrain
   CTETRA (39) linear (centroid)
   CPENTA (67) linear (centroid)  ## HEXA???
   CHEXA  (68) linear (centroid)  ## PENTA??

 OES_CBAR_34 - barStress / barStrain
   CBAR (34) linear

 OES_CTRIA3_74 - plateStress / plateStrain
   CTRIA3 (74)  linear 1 node (centroid)

 OES_CQUAD4_33 - plateStress / plateStrain
   CQUAD4 (33)  linear 1 node (centroid)

 OES_CQUAD4_144 - plateStress / plateStrain
   CQUAD4 (144) linear 5 nodes (centroid + 4 corners)
   CQUAD8 (82)  linear 5 nodes (centroid + 4 corners)
   CTRIA6 (75)  linear 4 nodes (centroid + 3 corners)
   CTRIAR (70)  linear 4 nodes (centroid + 3 corners)

 OES_CQUAD4_95 - compositePlateStress/compositePlateStrain - TODO
   CQUAD4 (95)  composite ? nodes (centroid + 4 corners)
   CQUAD8 (96)  composite ? nodes (centroid + 4 corners)
   CTRIA6 (97)  composite ? nodes (centroid + 3 corners)
   CTRIAR (98)  composite ? nodes (centroid + 3 corners)

OES_CBUSH_102 - bushStress / bushStrain
   CBUSH (102) linear 1 node (centroid)

 OES_CELAS_224_225 - nonlinearSpringStress - TODO
   CELAS1 (224???)  nonlinear (centroid)
   CELAS2 (???)     nonlinear (centroid)
   CELAS3 (225???)  nonlinear (centroid)

"""
#pylint: disable=C0103,C0301,R0914,E1101
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
from numpy import ones
from struct import unpack
#from pyNastran import isRelease

class RealElementsStressStrain2(object):

    def OES_basicElement(self, name):
        """
        genericStressReader - works on CROD_1, CELAS2_12
        stress & strain
        format_code=1 sort_code=0 (eid,axial,axialMS,torsion,torsionMS)
        """
        dt = self.nonlinear_factor
        (format1, extract) = self.getOUG_FormatStart()
        (ntotal, dataFormat) = self.obj.getLength()
        dataFormat = format1 + dataFormat
        #print("ntotal=%s dataFormat=%s len(data)=%s" % (ntotal,dataFormat,len(self.data)))
        dataFormat = bytes(dataFormat)

        nelements = len(self.data) // ntotal
        if self.read_mode == 0 or name not in self._selected_names:
            if name not in self._result_names:
                self._result_names.append(name)
            # figure out the shape
            self.obj._increase_size(dt, nelements)
            iend = ntotal * nelements
            self.data = self.data[iend:]
            return

        n = 0
        for i in xrange(nelements):
            eData = self.data[n:n + ntotal]
            out = unpack(dataFormat, eData)
            #print "out = ",out
            eid = extract(out[0], dt)
            self.obj.add_new_eid(dt, eid, out[1:])
            n += ntotal
        self.data = self.data[n:]

    #-------------------------------------------------------------------------
    # beamStress / beamStrain
    #-------------------------------------------------------------------------
    def OES_CBEAM_2(self, name):
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
        if self.read_mode == 0 or name not in self._selected_names:
            if name not in self._result_names:
                self._result_names.append(name)
            # figure out the shape
            self.obj._increase_size(dt, nelements, nnodes)
            iend = ntotal * nelements
            self.data = self.data[iend:]
            return
        else:  # read_mode = 1; # we know the shape so we can make a pandas matrix
            #print("nelements =",nelements)
            pass
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
            #print('type', self.obj.__class__.__name__)
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

    #-------------------------------------------------------------------------
    # solidStress / solidStrain
    #-------------------------------------------------------------------------
    def OES_CSOLID_39_67_68(self, name):
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

        #nnodes = 5 # 1 centroid + 4 corner points
        #self.print_section(20)
        #term      = self.data[0:4] CEN/
        #self.data = self.data[4:]
        element_type = self.get_element_type(self.element_type)
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
            msg = ('not supported....EType=%s eType=%s '
                   'num_wide=%s' % (ElementType, self.element_type, self.num_wide))
            raise NotImplementedError(msg)

        ntotal = 16 + 84 * nnodes_expected
        nelements = len(self.data) // ntotal
        overflow = len(self.data) % ntotal

        nnodes = nelements * nnodes_expected
        if self.read_mode == 0 or name not in self._selected_names:
            if name not in self._result_names:
                self._result_names.append(name)

            # figure out the shape
            self.obj._increase_size(dt, nelements, nnodes)
            iend = ntotal * nelements
            self.data = self.data[iend:]
            return
        else:
            # read_mode = 1; # we know the shape so we can make a pandas matrix
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

            assert nNodes < 21, 'nNodes=%s' % nNodes #, self.print_block(eData)

            etypes.append(element_type)
            cids.append(cid)
            eids.append(eid2)

            ibase += 16
            for nodeID in xrange(nnodes_expected):  # nodes pts, +1 for centroid (???)
                out = unpack(b'i20f',self.data[ibase:ibase + 84]) # 4*21 = 84
                (grid_device, sxx, sxy, s1, a1, a2, a3, pressure, svm,
                              syy, syz, s2, b1, b2, b3,
                              szz, sxz, s3, c1, c2, c3) = out
                #print("grid_device=%s" % grid_device)
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
                o1.append(s1)
                o2.append(s3)
                o3.append(s3)
                ovmShear.append(svm)

                #print("%s eid=%s cid=%s grid=%3s sxx=%-6i txy=%-5i s1=%-6i a1=%i a2=%i a3=%i press=%i vm=%s" % (element_type,eid2,cid,grid,sxx,sxy,s1,a1,a2,a3,pressure,svm))
                #print("%s eid=%s cid=%s grid=%3s syy=%-6i tyz=%-5i s2=%-6i b1=%i b2=%i b3=%i"                % (element_type,eid2,cid,grid,syy,syz,s2,b1,b2,b3))
                #print("%s eid=%s cid=%s grid=%3s szz=%-6i txz=%-5i s3=%-6i c1=%i c2=%i c3=%i"                % (element_type,eid2,cid,grid,szz,sxz,s3,c1,c2,c3))
                #print("")
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

        [koxx, koyy, kozz, ktxy, ktxz, ktyz, ko1, ko2, ko3,  kovm] = self.obj._get_headers()
        assert len(oxx) == inode_end - inode_start, 'inode_start=%s inode_end=%s len(oxx)=%s' % (inode_start, inode_end, len(oxx))
        self.obj.data[koxx][inode_start:inode_end] = oxx
        self.obj.data[koyy][inode_start:inode_end] = oyy
        self.obj.data[kozz][inode_start:inode_end] = ozz

        self.obj.data[ktxy][inode_start:inode_end] = txy
        self.obj.data[ktxz][inode_start:inode_end] = txz
        self.obj.data[ktyz][inode_start:inode_end] = tyz

        self.obj.data[ko1][inode_start:inode_end] = o1
        self.obj.data[ko2][inode_start:inode_end] = o2
        self.obj.data[ko3][inode_start:inode_end] = o3
        self.obj.data[kovm][inode_start:inode_end] = ovmShear

        assert len(oxx) == inode_end - inode_start, 'ielement_start=%s ielement_end=%s len(eids)=%s' % (ielement_start, ielement_end, len(eids))
        self.obj.element_data['element_id'][ielement_start:ielement_end] = eids
        self.obj.element_data['element_type'][ielement_start:ielement_end] = etypes
        self.obj.element_data['cid'][ielement_start:ielement_end] = cids
        self.obj.element_data['nnodes'][ielement_start:ielement_end] = ones(nelements) * nnodes_expected

        # pressure
        # aCos
        # bCos
        # cCos
        #self.obj.data[''][inode_start:inode_end] = translations[:, 5]

        if len(self.obj.data['element_id']) == inode_end: # [nodes, elements]
            #print("self.element_data =\n", self.obj.element_data)
            self.obj._finalize(dt)
            #print("self.element_data =\n", self.obj.element_data)
        else:
            #print('increment...', overflow, ntotal)
            self.obj._increment(nnodes, nelements)

    #-------------------------------------------------------------------------
    # barStress / barStrain
    #-------------------------------------------------------------------------
    def OES_CBAR_34(self, name):
        dt = self.nonlinear_factor
        assert self.num_wide == 16, ('invalid num_wide...num_wide=%s'
                                    % (self.num_wide))

        (format1, extract) = self.getOUG_FormatStart()
        format1 += '15f'
        format1 = bytes(format1)

        ntotal = 64
        nelements = len(self.data) // 64
        nnodes = nelements
        if self.read_mode == 0 or name not in self._selected_names:
            if name not in self._result_names:
                self._result_names.append(name)
            # figure out the shape
            self.obj._increase_size(dt, nelements, nnodes)
            iend = ntotal * nelements
            self.data = self.data[iend:]
            return
        else:  # read_mode = 1; # we know the shape so we can make a pandas matrix
            #print("nelements =",nelements)
            pass

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

        #print('len(eids) = ', len(self.obj.data['element_id']))
        #print('inode_end // ntotal = ', inode_end)
        if len(self.obj.data['element_id']) == inode_end: # [nodes, elements]
            self.obj._finalize(dt)
        else:
            self.obj._increment(nnodes, nelements)

    #-------------------------------------------------------------------------
    # bushStress / bushStrain
    #-------------------------------------------------------------------------
    def OES_CBUSH_102(self, name):
        dt = self.nonlinear_factor
        (format1, extract) = self.getOUG_FormatStart()

        assert self.num_wide == 7, "num_wide=%s not 7" % self.num_wide
        ntotal = 28  # 4*7
        format1 += '6f'
        format1 = bytes(format1)

        n = 0
        nelements = len(self.data) // ntotal

        if self.read_mode == 0 or name not in self._selected_names:
            if name not in self._result_names:
                self._result_names.append(name)
            # figure out the shape
            self.obj._increase_size(dt, nelements)
            iend = ntotal * nelements
            self.data = self.data[iend:]
            return
        else:  # read_mode = 1; # we know the shape so we can make a pandas matrix
            #print("nelements =",nelements)
            pass

        (ielement_start, ielement_end) = self.obj._preallocate(dt, nelements)

        istart = 0
        iend = ntotal
        etypes=['CBUSH'] * nelements

        eids = []
        T1 = []; T2=[]; T3=[]; R1=[]; R2=[]; R3=[]
        for i in xrange(nelements):
            eData = self.data[n:n + ntotal]
            out = unpack(format1, eData)  # num_wide=7

            (eid, tx, ty, tz, rx, ry, rz) = out
            eid = extract(eid, dt)

            eids.append(eid)
            T1.append(tx)
            T2.append(ty)
            T3.append(tz)
            R1.append(rx)
            R2.append(ry)
            R3.append(rz)

            #self.obj.add_new_eid(self.element_type, dt, eid, tx, ty, tz, rx, ry, rz)
            n += ntotal
        self.data = self.data[n:]

        #print('delta', inode_end - inode_start, len(eids))
        if dt:
            name = self.obj.data_code['name']
            self.obj.data[name][ielement_start:ielement_end] = ones(ielement_end - ielement_start) * dt

        #print("inode_start=%r inode_end=%r delta=%r" % (inode_start, inode_end, inode_end-inode_start))
        #print('len(s1) =', len(s1))

        self.obj.data['element_id'][ielement_start:ielement_end] = eids
        self.obj.data['element_type'][ielement_start:ielement_end] = etypes
        #if self.obj.isStress():
        self.obj.data['T1'][ielement_start:ielement_end] = T1
        self.obj.data['T2'][ielement_start:ielement_end] = T2
        self.obj.data['T3'][ielement_start:ielement_end] = T3
        self.obj.data['R1'][ielement_start:ielement_end] = R1
        self.obj.data['R2'][ielement_start:ielement_end] = R2
        self.obj.data['R3'][ielement_start:ielement_end] = R3

        #print('len(eids) = ', len(self.obj.data['element_id']))
        if len(self.obj.data['element_id']) == ielement_end: # [nodes, elements]
            self.obj._finalize(dt)
        else:
            self.obj._increment(nelements)

    #-------------------------------------------------------------------------
    # compositePlateStress / compositePlateStrain
    #-------------------------------------------------------------------------
    def OES_CQUAD4_95(self, name):  # works (doesnt handle all stress/strain cases tho)
        """
        GRID-ID  DISTANCE,NORMAL-X,NORMAL-Y,SHEAR-XY,ANGLE,MAJOR MINOR,VONMISES
        composite quad

         95 - CQUAD4
         96 - CQUAD8
         97 - CTRIA3
         98 - CTRIA6 (composite)
        """
        dt = self.nonlinear_factor
        nloops = len(self.data) // 44  # 2+17*5 = 87 -> 87*4 = 348
        #print("nloops = ", nloops)
        ntotal = 44
        if self.read_mode == 0 or name not in self._selected_names:

            if name not in self._result_names:
                self._result_names.append(name)
                
                ibase = 0
                #nlayers = 0
                eid_nlayer = {}
                (format1, extract) = self.getOUG_FormatStart()
                format1 += 'i'
                format1 = bytes(format1)
                #print('format1 = %r' % format1)
                for i in xrange(nloops):
                    eData = self.data[ibase:ibase+8]  # 4*11
                    (eid, iLayer) = unpack(format1, eData)
                    eid = extract(eid, dt)
                    if eid != self.eid2:  # originally initialized to None, the buffer doesnt reset it, so it is the old value
                        eid_nlayer[eid] = 1
                    else:
                        eid_nlayer[eid] += 1
                    #nlayers += 1
                    #self.eid2 = eid
                    #ibase += 44
                #print('eid_layer =', eid_layer)
                #print('nlayers =', nlayers)
                #self.data = self.data[ibase:]

                # figure out the shape
                #print("nloops=%s" % nloops)
                self.obj._increase_size(dt, nloops, eid_nlayer)
                iend = ntotal * nloops
                self.data = self.data[iend:]
                return
            else:  # read_mode = 1; # we know the shape so we can make a pandas matrix
                #print("nelements =", nelements)
                pass

        (inode_start, inode_end, ielement_start, ielement_end
            ) = self.obj._preallocate(dt, nloops)

        eType = self.get_element_type(self.element_type)

        #self.print_section(20)
        #term = data[0:4] CEN/
        #data = data[4:]
        if self.num_wide != 11:
            raise RuntimeError('invalid num_wide; num_wide=%s' % self.num_wide)

        (format1, extract) = self.getOUG_FormatStart()
        format1 += 'i9f'
        format1 = bytes(format1)

        ibase = 0
        eids = []
        eids_elements = []
        ilayers = []
        etypes = [eType] * nloops
        o1s = []
        o2s = []
        
        for i in xrange(nloops):
        #while len(self.data) <= 44:
            eData = self.data[ibase:ibase+44]  # 4*11
            (eid, ilayer, o1, o2, t12, t1z, t2z, angle, major, minor, ovm) = unpack(format1, eData)
            eid = extract(eid, dt)
            
            #eids_elements.append(eid)

            eids.append(eid)
            ilayers.append(ilayer)
            o1s.append(o1)
            o2s.append(o2)

            if eid != self.eid2:  # originally initialized to None, the buffer doesnt reset it, so it is the old value
                #print "1 - eid=%s iLayer=%i o1=%i o2=%i ovm=%i" % (eid,iLayer,o1,o2,ovm)
                self.obj.add_new_eid(eType, dt, eid, o1, o2, t12, t1z, t2z, angle, major, minor, ovm)
            else:
                #print "2 - eid=%s iLayer=%i o1=%i o2=%i ovm=%i" % (eid,iLayer,o1,o2,ovm)
                self.obj.add(dt, eid, o1, o2, t12, t1z, t2z, angle, major, minor, ovm)
            self.eid2 = eid
            ibase += 44
            #self.dn += 348
        self.data = self.data[ibase:]
        #print "3 - eid=%s iLayer=%i o1=%i o2=%i ovm=%i" % (eid,iLayer,o1,o2,ovm)
        #----------------------------------------------------------------------
        self.obj.element_data['element_id'][ielement_start:ielement_end] = eids_elements
        #print('self.obj.element_data\n', self.obj.element_data)
        self.obj.element_data['element_type'][ielement_start:ielement_end] = etypes
        #self.obj.element_data['nlayer'][ielement_start:ielement_end] = nnodes_list
        #print('self.obj.element_data\n', self.obj.element_data)

        self.obj.data['element_id'][inode_start:inode_end] = eids
        self.obj.data['layer'][inode_start:inode_end] = layers
        self.obj.data['o1'][inode_start:inode_end] = o1s
        self.obj.data['o2'][inode_start:inode_end] = o2s

        #headers = self.obj._get_headers()
        #(kfd, koxx, koyy, ktxy, komax, komin, kovm) = headers

        assert  inode_end - inode_start == len(fd)
        if dt:
            name = self.obj.data_code['name']
            self.obj.data[name][inode_start:inode_end] = ones(inode_end - inode_start) * dt
        self.obj.data['element_id'][inode_start:inode_end] = eids

        self.obj.data[kfd][inode_start:inode_end] = fd
        #print(self.obj.data.keys())
        #print('self.obj.data\n', self.obj.data)
        self.obj.data[koxx][inode_start:inode_end] = sx
        self.obj.data[koyy][inode_start:inode_end] = sy
        self.obj.data[ktxy][inode_start:inode_end] = txy
        self.obj.data[komax][inode_start:inode_end] = major
        self.obj.data[komin][inode_start:inode_end] = minor
        self.obj.data['angle'][inode_start:inode_end] = angle
        self.obj.data[kovm][inode_start:inode_end] = vm

        #self.obj.data['element_id'][ielement_start:ielement_end] = eids
        #self.obj.data['element_type'][ielement_start:ielement_end] = etypes

        #print('len(eids) = ', len(self.obj.data['element_id']))
        #print('inode_end // ntotal = ', inode_end)
        if self.obj._is_full(nnodes, nelements):
        #if len(self.obj.data['element_id']) == inode_end: # [nodes, elements]
            self.obj._finalize(dt)

    #-------------------------------------------------------------------------
    # plateStress / plateStrain
    #-------------------------------------------------------------------------
    def OES_CTRIA3_74(self, name):
        """
        DISTANCE,NORMAL-X,NORMAL-Y,SHEAR-XY,ANGLE,MAJOR,MINOR,VONMISES
        stress is extracted at the centroid
        """
        assert self.num_wide == 17, ('invalid num_wide...num_wide=%s'
                                    % self.num_wide)

        dt = self.nonlinear_factor
        (format1, extract) = self.getOUG_FormatStart()
        format1 += '16f'
        format1 = bytes(format1)

        ntotal = 68  # 4*17
        n = 0
        nelements = len(self.data) // ntotal
        nnodes = nelements * 2 # 2 layers at the centroid
        #nnodes = 1
        #nrows = nelements * 2 # 2 layers
        
        if self.read_mode == 0 or name not in self._selected_names:
            if name not in self._result_names:
                self._result_names.append(name)

            # figure out the shape
            #print("nnodes=%s nelements=%s" % (nnodes, nelements))
            self.obj._increase_size(dt, nnodes, nelements)
            iend = ntotal * nelements
            self.data = self.data[iend:]
            return
        else:  # read_mode = 1; # we know the shape so we can make a pandas matrix
            #print("nelements =", nelements)
            pass

        (inode_start, inode_end, ielement_start, ielement_end
            ) = self.obj._preallocate(dt, nnodes, nelements)
        #print("inode_start=%s inode_end=%s ielement_start=%s ielement_end=%s" % (inode_start, inode_end, ielement_start, ielement_end))

        eids=[]; eTypes=[]; eids2=[]; nids=[]; layer=[]
        fd=[]; sx=[]; sy=[]; txy=[]; angle=[]; major=[]; minor=[]; vm=[]
        eTypes = ['CTRIA3'] * nelements
        nnodes_list = ones(nelements)
        for i in xrange(nelements):
            eData = self.data[n:n + ntotal]
            out = unpack(format1, eData)

            (eid, fd1i, sx1i, sy1i, txy1i, angle1i, major1i, minor1i, vm1i,
                  fd2i, sx2i, sy2i, txy2i, angle2i, major2i, minor2i, vm2i) = out
            eid = extract(eid, dt)

            #eTypes.append(eType)
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

            #print("eid=%2i fd1=%i sx1=%i sy1=%i txy1=%i angle1=%i major1=%i minor1=%i vm1=%i" % (eid,fd1i,sx1i,sy1i,txy1i,angle1i,major1i,minor1i,vm1i))
            #print("        fd2=%i sx2=%i sy2=%i txy2=%i angle2=%i major2=%i minor2=%i vm2=%i\n"   % (fd2i,sx2i,sy2i,txy2i,angle2i,major2i,minor2i,vm2i))
            n += ntotal
        self.data = self.data[n:]
        #-----------------------------------------------------------------------
        #print('delta', inode_end - inode_start, len(eids))
        #istart = self.obj._size_start
        #iend = istart + nnodes

        #print('---------------------------------------------------------')
        #print("inode_start=%r inode_end=%r delta=%r" % (inode_start, inode_end, inode_end-inode_start))
        #print('len(s1) =', len(s1))

        self.obj.element_data['element_id'][ielement_start:ielement_end] = eids
        #print('self.obj.element_data\n', self.obj.element_data)
        self.obj.element_data['element_type'][ielement_start:ielement_end] = eTypes
        self.obj.element_data['nnodes'][ielement_start:ielement_end] = nnodes_list
        #print('self.obj.element_data\n', self.obj.element_data)
        #print('nids',nids)
        assert  inode_end - inode_start == len(eids2), 'inode_start=%s inode_end=%s len(eids2)=%s' % (inode_start, inode_end, len(eids2))
        self.obj.data['element_id'][inode_start:inode_end] = eids2
        self.obj.data['node_id'][inode_start:inode_end] = nids
        self.obj.data['layer'][inode_start:inode_end] = layer

        headers = self.obj._get_headers()
        (kfd, koxx, koyy, ktxy, komax, komin, kovm) = headers

        assert  inode_end - inode_start == len(fd)
        if dt:
            name = self.obj.data_code['name']
            self.obj.data[name][inode_start:inode_end] = ones(inode_end - inode_start) * dt
        self.obj.data['element_id'][inode_start:inode_end] = eids

        self.obj.data[kfd][inode_start:inode_end] = fd
        #print(self.obj.data.keys())
        #print('self.obj.data\n', self.obj.data)
        self.obj.data[koxx][inode_start:inode_end] = sx
        self.obj.data[koyy][inode_start:inode_end] = sy
        self.obj.data[ktxy][inode_start:inode_end] = txy
        self.obj.data[komax][inode_start:inode_end] = major
        self.obj.data[komin][inode_start:inode_end] = minor
        self.obj.data['angle'][inode_start:inode_end] = angle
        self.obj.data[kovm][inode_start:inode_end] = vm

        #self.obj.data['element_id'][ielement_start:ielement_end] = eids
        #self.obj.data['element_type'][ielement_start:ielement_end] = etypes

        #print('len(eids) = ', len(self.obj.data['element_id']))
        #print('inode_end // ntotal = ', inode_end)
        if self.obj._is_full(nnodes, nelements):
        #if len(self.obj.data['element_id']) == inode_end: # [nodes, elements]
            self.obj._finalize(dt)
        #else:
            #self.obj._increment(nnodes, nelements)

    def OES_CQUAD4_33(self, name):
        """
        GRID-ID  DISTANCE,NORMAL-X,NORMAL-Y,SHEAR-XY,ANGLE,MAJOR MINOR,VONMISES
        """
        dt = self.nonlinear_factor
        (format1, extract) = self.getOUG_FormatStart()
        format1 += '16f'
        format1 = bytes(format1)

        nNodes = 0  # centroid + 4 corner points
        #self.print_section(20)
        #term = data[0:4] CEN/
        #data = data[4:]
        #print "*****"
        #self.print_block(self.data)

        assert self.num_wide == 17, ('invalid num_wide...num_wide=%s'
                                    % self.num_wide)

        istart = 0
        iend = 68
        ntotal = 68 # 4*17
        nelements = len(self.data) // 68

        nnodes = nelements * 2 # 2 layers
        #nnodes = nrows #nelements * nNodes
        if self.read_mode == 0 or name not in self._selected_names:
            if name not in self._result_names:
                self._result_names.append(name)

            # figure out the shape
            self.obj._increase_size(dt, nnodes, nelements)
            iend = ntotal * nelements
            self.data = self.data[iend:]
            return
        else:  # read_mode = 1; # we know the shape so we can make a pandas matrix
            #print("nelements =", nelements)
            pass
        
        (inode_start, inode_end, ielement_start, ielement_end
            ) = self.obj._preallocate(dt, nnodes, nelements)
        #print('***dt=%s nnodes=%s nelements=%s' % (dt, nnodes, nelements))

        eids=[]; eTypes=[]; eids2=[]; nids=[]; layer=[]; nnodes_list=[]
        fd=[]; sx=[]; sy=[]; sxy=[]; angle=[]; major=[]; minor=[]; svm=[]

        layer += [1, 2] * nelements
        #eTypes = ['CQUAD4'] * nelements
        eTypes = []

        for i in xrange(nelements):
            #print self.print_block(self.data[0:100])
            #(eid,) = unpack(b'i',self.data[0:4])
            #self.data = self.data[8:]  # 2
            eData = self.data[istart:iend]
            #print("len(eData) = ", len(eData))
            out = unpack(format1, eData)  # 17
            (eid, fd1, sx1, sy1, txy1, angle1, major1, minor1, ovm1,
                  fd2, sx2, sy2, txy2, angle2, major2, minor2, ovm2) = out

            eid = extract(eid, dt)

            #print("eid=%i grid=%s fd1=%-3.1f sx1=%i sy1=%i txy1=%i angle1=%i major1=%i minor1=%i vm1=%i" % (eid,'C',fd1,sx1,sy1,txy1,angle1,major1,minor1,maxShear1))
            #print(  "             fd2=%-3.1f sx2=%i sy2=%i txy2=%i angle2=%i major2=%i minor2=%i vm2=%i\n"       % (fd2,sx2,sy2,txy2,angle2,major2,minor2,maxShear2))
            #print "nNodes = ",nNodes
            #self.obj.add_new_eid('CQUAD4', dt, eid, 'C', fd1, sx1, sy1,
            #                   txy1, angle1, major1, minor1, maxShear1)
            #self.obj.add(dt, eid, 'C', fd2, sx2, sy2, txy2,
            #             angle2, major2, minor2, maxShear2)

            istart = iend
            iend += ntotal
            eids2.append(eid)

            eids.append(eid)
            eids.append(eid)
            fd.append(fd1)
            fd.append(fd2)
            sx.append(sx1)
            sx.append(sx2)
            sy.append(sy1)
            sy.append(sy2)
            angle.append(angle1)
            angle.append(angle2)
            minor.append(minor1)
            minor.append(minor2)
            major.append(major1)
            major.append(major2)
            svm.append(ovm1)
            svm.append(ovm2)
            eTypes.append('CQUAD4_33')
        
        from numpy import array
        eTypes = array(eTypes)
        self.data = self.data[istart:]
        #-----------------------------------------------------------------------
        #print('---------------------------------------------------------')
        inode_start = self.obj._inode_start
        inode_end = inode_start + nnodes

        ielement_start = self.obj._ielement_start
        ielement_end = ielement_start + nelements

        #print('delta', inode_end - inode_start, len(eids))
        
        #nnodes = inode_end - inode_start
        #nnodes = inode_end - inode_start

        ndt, nelements_size, nnodes_size, dts = self.obj._get_shape()
        #print("ndt=%s nelements_size=%s nnodes_size=%s dts=%s" % (ndt, nelements_size, nnodes_size, str(dts)))

        #print("inode_start=%r inode_end=%r delta=%r" % (inode_start, inode_end, inode_end-inode_start))
        #print("ielement_start=%r ielement_end=%r delta=%r" % (ielement_start, ielement_end, ielement_end-ielement_start))
        #print('len(svm) =', len(svm))
        assert nelements == len(eids2), 'nelements=%s len(eids)=%s' % (nelements, len(eids2))
        assert len(eids2) == ielement_end - ielement_start, 'ielement_start=%s ielement_end=%s len(eids)=%s' % (ielement_start, ielement_end, len(eids2))
        #print('len(element_data)=', len(self.obj.element_data))
        self.obj.element_data['element_id'][ielement_start:ielement_end] = eids2
        #print(self.obj.element_data)
        self.obj.element_data['element_type'][ielement_start:ielement_end] = eTypes
        self.obj.element_data['nnodes'][ielement_start:ielement_end] = ones(nelements)

        #print('nids',nids)
        assert len(eids) == 2*nelements, '2*nelements=%s len(eids)=%s' % (2*nelements, len(eids))
        assert len(eids) == inode_end - inode_start, 'inode_start=%s inode_end=%s delta=%s len(eids)=%s' % (inode_start, inode_end,
            inode_end - inode_start, len(eids))

        if dt:
            name = self.obj.data_code['name']
            self.obj.data[name][inode_start:inode_end] = ones(nnodes) * dt
        self.obj.data['element_id'][inode_start:inode_end] = eids
        #self.obj.data['element_id'][inode_start:inode_end] = eids
        self.obj.data['node_id'][inode_start:inode_end] = nids
        self.obj.data['layer'][inode_start:inode_end] = layer
        #print(self.obj.element_data)

        headers = self.obj._get_headers()
        (kfd, koxx, koyy, ktxy, komax, komin, kovm) = headers

        self.obj.data[kfd][inode_start:inode_end] = fd
        self.obj.data[koxx][inode_start:inode_end] = sx
        self.obj.data[koyy][inode_start:inode_end] = sy
        self.obj.data[ktxy][inode_start:inode_end] = sxy
        self.obj.data[komax][inode_start:inode_end] = major
        self.obj.data[komin][inode_start:inode_end] = minor
        self.obj.data['angle'][inode_start:inode_end] = angle
        self.obj.data[kovm][inode_start:inode_end] = svm

        #self.obj.data['element_id'][ielement_start:ielement_end] = eids
        #self.obj.data['element_type'][ielement_start:ielement_end] = etypes

        #print('len(eids) = ', len(self.obj.data['element_id']))
        #print('inode_end // ntotal = ', inode_end)
        #if len(self.obj.data['element_id']) == inode_end: # [nodes, elements]
        if self.obj._is_full(nnodes, nelements):
            #print('finalize')
            #sys.exit('finalize')
            self.obj._finalize(dt)
        #else:
            #print('increment')
            #self.obj._increment(nnodes, nelements)

    def OES_CQUAD4_144(self, name):
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
                                      % self.element_type)

        dt = self.nonlinear_factor
        (format1, extract) = self.getOUG_FormatStart()
        format1 += '16f'
        format1 = bytes(format1)

        nelements = len(self.data) // ntotal
        nnodes = (1 + nNodes) * nelements * 2 # 2 layers
        #print('nrows = %i' % nrows)

        gridC = 'C'
        cformat = b'i4s'+format1  # center format
        nformat = b'i16f'         # node format

        #nnodes = nrows #nelements * nNodes
        if self.read_mode == 0 or name not in self._selected_names:
            if name not in self._result_names:
                self._result_names.append(name)

            # figure out the shape
            self.obj._increase_size(dt, nnodes, nelements)
            iend = ntotal * nelements
            self.data = self.data[iend:]
            return
        else:  # read_mode = 1; # we know the shape so we can make a pandas matrix
            #print("nelements =", nelements)
            pass

        (inode_start, inode_end, ielement_start, ielement_end
            ) = self.obj._preallocate(dt, nnodes, nelements)


        eids=[]; eTypes=[]; eids2=[]; nids=[]; layer=[] #; nnodes_list=[]
        fd=[]; sx=[]; sy=[]; sxy=[]; angle=[]; major=[]; minor=[]; vm=[]
        #['fd1', 'sx1', 'sy1', 'txy1', 'angle1', 'major1', 'minor1', 'vm1',
        # 'fd2', 'sx2', 'sy2', 'txy2', 'angle2', 'major2', 'minor2', 'vm2',]

        ibase = 0
        nnodes_temp = 1 + nNodes
        
        nnodes_list = ones(nelements) * nnodes_temp
        layer = [1, 2] * (nelements * nnodes_temp)

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

            fd.append(fd1i)
            sx.append(sx1i)
            sy.append(sy1i)
            sxy.append(txy1i)
            angle.append(angle1i)
            major.append(major1i)
            minor.append(minor1i)
            vm.append(vm1i)

            fd.append(fd2i)
            sx.append(sx2i)
            sy.append(sy2i)
            sxy.append(txy2i)
            angle.append(angle2i)
            major.append(major2i)
            minor.append(minor2i)
            vm.append(vm2i)

            #print("eid=%i grid=%i fd1=%i sx1=%i sy1=%i txy1=%i angle1=%i major1=%i minor1=%i vm1=%i" % (eid,grid,fd1,sx1,sy1,txy1,angle1,major1,minor1,vm1))
            #print("               fd2=%i sx2=%i sy2=%i txy2=%i angle2=%i major2=%i minor2=%i vm2=%i\n"        % (fd2,sx2,sy2,txy2,angle2,major2,minor2,vm2))
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

                #layer.append(1)
                #layer.append(2)

                fd.append(fd1i)
                sx.append(sx1i)
                sy.append(sy1i)
                sxy.append(txy1i)
                angle.append(angle1i)
                major.append(major1i)
                minor.append(minor1i)
                vm.append(vm1i)

                fd.append(fd2i)
                sx.append(sx2i)
                sy.append(sy2i)
                sxy.append(txy2i)
                angle.append(angle2i)
                major.append(major2i)
                minor.append(minor2i)
                vm.append(vm2i)
                #print("eid=%i grid=%i fd1=%i sx1=%i sy1=%i txy1=%i angle1=%i major1=%i minor1=%i vm1=%i" % (eid,grid,fd1,sx1,sy1,txy1,angle1,major1,minor1,vm1))
                #print("               fd2=%i sx2=%i sy2=%i txy2=%i angle2=%i major2=%i minor2=%i vm2=%i\n"          % (fd2,sx2,sy2,txy2,angle2,major2,minor2,vm2))
                #self.obj.addNewNode(dt, eidi, gridi, fd1i, sx1i, sy1i, txy1i, angle1i, major1i, minor1i, vm1i)
                #self.obj.add(       dt, eidi, gridi, fd2i, sx2i, sy2i, txy2i, angle2i, major2i, minor2i, vm2i)
                ibase += 68
        self.data = self.data[ibase:]

        #-----------------------------------------------------------------------
        #print('delta', inode_end - inode_start, len(eids))
        if dt:
            name = self.obj.data_code['name']
            self.obj.data[name][inode_start:inode_end] = ones(inode_end - inode_start) * dt

        self.obj.data['element_id'][inode_start:inode_end] = eids

        #print("inode_start=%r inode_end=%r delta=%r" % (inode_start, inode_end, inode_end-inode_start))
        #print('len(s1) =', len(s1))

        self.obj.element_data['element_id'][ielement_start:ielement_end] = eids
        self.obj.element_data['element_type'][ielement_start:ielement_end] = eTypes
        self.obj.element_data['nnodes'][ielement_start:ielement_end] = nnodes_list

        #print('nids',nids)
        self.obj.data['element_id'][inode_start:inode_end] = eids2
        self.obj.data['node_id'][inode_start:inode_end] = nids
        self.obj.data['layer'][inode_start:inode_end] = layer

        headers = self.obj._get_headers()
        (kfd, oxx, oyy, txy, omax, omin, ovm) = headers

        self.obj.data[kfd][inode_start:inode_end] = fd
        self.obj.data[oxx][inode_start:inode_end] = sx
        self.obj.data[oyy][inode_start:inode_end] = sy
        self.obj.data[txy][inode_start:inode_end] = sxy
        self.obj.data[omax][inode_start:inode_end] = major
        self.obj.data[omin][inode_start:inode_end] = minor
        self.obj.data['angle'][inode_start:inode_end] = angle
        self.obj.data[ovm][inode_start:inode_end] = vm

        #self.obj.data['element_id'][ielement_start:ielement_end] = eids
        #self.obj.data['element_type'][ielement_start:ielement_end] = etypes

        #print('len(eids) = ', len(self.obj.data['element_id']))
        #print('inode_end // ntotal = ', inode_end)
        #if len(self.obj.data['element_id']) == inode_end: # [nodes, elements]
        if self.obj._is_full(nnodes, nelements):
            self.obj._finalize(dt)
        #else:
            #self.obj._increment(nnodes, nelements)
