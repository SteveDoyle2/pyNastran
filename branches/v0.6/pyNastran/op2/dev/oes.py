from struct import Struct, unpack

from pyNastran.op2.dev.op2_common import OP2Common
from pyNastran.op2.op2_helper import polar_to_real_imag


class OES(OP2Common):
    def __init__(self):
        OP2Common.__init__(self)

    def read_oes1_3(self, data):
        self.words = ['aCode',       'tCode',    'element_type', 'isubcase',
                 '???',         '???',      '???',          'load_set'
                 'format_code', 'num_wide', 's_code',       '???',
                 '???',         '???',      '???',          '???',
                 '???',         '???',      '???',          '???',
                 '???',         '???',      '???',          '???',
                 '???', 'Title', 'subtitle', 'label']

        self.parse_approach_code(data)  # 3

        ## element type
        self.element_type = self.add_data_parameter(data, 'element_type', 'i', 3, False)

        ## load set ID
        self.load_set = self.add_data_parameter(data, 'load_set', 'i', 8, False)

        ## format code
        self.format_code = self.add_data_parameter(data, 'format_code', 'i', 9, False)

        ## number of words per entry in record
        ## .. note:: is this needed for this table ???
        self.num_wide = self.add_data_parameter(data, 'num_wide', 'i', 10, False)

        ## stress/strain codes
        self.s_code = self.add_data_parameter(data, 's_code', 'i', 11, False)

        ## thermal flag; 1 for heat ransfer, 0 otherwise
        self.thermal = self.add_data_parameter(data, 'thermal', 'i', 23, False)

        ## assuming tCode=1
        if self.analysis_code == 1:   # statics / displacement / heat flux
            ## load set number
            self.lsdvmn = self.add_data_parameter(data, 'lsdvmn', 'i', 5, False)
            self.dataNames = self.apply_data_code_value('dataNames', ['lsdvmn'])
            self.setNullNonlinearFactor()
        elif self.analysis_code == 2:  # real eigenvalues
            ## mode number
            self.mode = self.add_data_parameter(data, 'mode', 'i', 5)
            ## real eigenvalue
            #self.eign = self.add_data_parameter(data, 'eign', 'f', 6, False)
            ## mode or cycle TODO confused on the type - F1???
            self.mode_cycle = self.add_data_parameter(data, 'mode_cycle', 'f', 7, False)
            self.dataNames = self.apply_data_code_value('dataNames', ['mode', 'eigr', 'mode_cycle'])
        #elif self.analysis_code==3: # differential stiffness
            #self.lsdvmn = self.get_values(data,'i',5) ## load set number
            #self.data_code['lsdvmn'] = self.lsdvmn
        #elif self.analysis_code==4: # differential stiffness
        #    self.lsdvmn = self.get_values(data,'i',5) ## load set number
        elif self.analysis_code == 5:   # frequency
            ## frequency
            self.freq = self.add_data_parameter(data, 'freq', 'f', 5)
            self.dataNames = self.apply_data_code_value('dataNames', ['freq'])
        elif self.analysis_code == 6:  # transient
            ## time step
            self.dt = self.add_data_parameter(data, 'dt', 'f', 5)
            self.dataNames = self.apply_data_code_value('dataNames', ['dt'])
        elif self.analysis_code == 7:  # pre-buckling
            ## load set
            self.lsdvmn = self.add_data_parameter(data, 'lsdvmn', 'i', 5)
            self.dataNames = self.apply_data_code_value('dataNames', ['lsdvmn'])
        elif self.analysis_code == 8:  # post-buckling
            ## mode number
            self.lsdvmn = self.add_data_parameter(data, 'lsdvmn', 'i', 5)
            self.eigr = self.add_data_parameter(data, 'eigr', 'f', 6, False)  # real eigenvalue
            self.dataNames = self.apply_data_code_value('dataNames', ['lsdvmn', 'eigr'])
        elif self.analysis_code == 9:  # complex eigenvalues
            ## mode number
            self.mode = self.add_data_parameter(data, 'mode', 'i', 5)
            ## real eigenvalue
            self.eigr = self.add_data_parameter(data, 'eigr', 'f', 6, False)
            ## imaginary eigenvalue
            self.eigi = self.add_data_parameter(data, 'eigi', 'f', 7, False)
            self.dataNames = self.apply_data_code_value('dataNames', ['mode', 'eigr', 'eigi'])
        elif self.analysis_code == 10:  # nonlinear statics
            ## load step
            self.lftsfq = self.add_data_parameter(data, 'lftsfq', 'f', 5)
            self.dataNames = self.apply_data_code_value('dataNames', ['lftsfq'])
        elif self.analysis_code == 11:  # old geometric nonlinear statics
            ## load set number
            self.lsdvmn = self.add_data_parameter(data, 'lsdvmn', 'i', 5)
            self.dataNames = self.apply_data_code_value('dataNames', ['lsdvmn'])
        elif self.analysis_code == 12:  # contran ? (may appear as aCode=6)  --> straight from DMAP...grrr...
            ## Time step ??? --> straight from DMAP
            self.dt = self.add_data_parameter(data, 'dt', 'f', 5)
            self.dataNames = self.apply_data_code_value('dataNames', ['dt'])
        else:
            raise RuntimeError('invalid analysis_code...analysis_code=%s' %
                               self.analysis_code)
        # tCode=2
        #if self.analysis_code==2: # sort2
        #    self.lsdvmn = self.get_values(data,'i',5)

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

    def read_oes1_4(self, data):
        #assert self.isubtable == -4, self.isubtable
        #if self.debug:
            #self.binary_debug.write('  element_name = %r\n' % self.element_name)
        #print "element_name =", self.element_name

        if self.is_sort1():
            self.read_oes1_4_sort1(data)
        else:
            raise NotImplementedError('sort2 Type=%s num=%s' % (self.element_name, self.element_type))
        if self.debug:
            self.binary_debug.write('*'* 20 + '\n\n')

    def read_oes1_4_sort1(self, data):
        if self.thermal == 0:
            self.read_oes_loads(data)
        elif self.thermal == 1:
            self.read_oes_thermal(data)
        else:
            raise NotImplementedError(self.thermal)

    def read_oes_thermal(self, data):
        n = 0

    def read_oes_loads(self, data):
        n = 0
        is_magnitude_phase = self.is_magnitude_phase()
        if self.element_type in [1, 3, 10]: # CROD, CTUBE, CONROD
            if self.num_wide == 5:  # real
                ntotal = 5 * 4
                nelements = len(data) // ntotal
                if self.debug:
                    self.binary_debug.write('  [cap, element1, element2, ..., cap]\n')
                    self.binary_debug.write('  cap = %i  # assume 1 cap when there could have been multiple\n' % len(data))
                    self.binary_debug.write('  #elementi = [eid_device, axial, axial_margin, torsion, torsion_margin]\n')
                    self.binary_debug.write('  nelements=%i; nnodes=1 # centroid\n' % nelements)
                format1 = b'i4f'
                s = Struct(format1)
                for i in xrange(nelements):
                    edata = data[n:n+ntotal]
                    out = s.unpack(edata)
                    (eid_device, axial, axial_margin, torsion, torsion_margin) = out
                    eid = (eid_device - self.device_code) // 10
                    if self.debug4():
                        self.binary_debug.write('  ----------\n')
                        self.binary_debug.write('  eid = %i\n' % eid)
                        self.binary_debug.write('  C = [%s]\n' % ', '.join(['%r' % di for di in out]) )

                    #print "eid =", eid
                    n += ntotal
            elif self.num_wide == 5: # real/imaginary
                ntotal = 20
                format1 = b'i4f'
                nelements = len(data) // ntotal
                for i in xrange(nelements):
                    edata = data[n:n + ntotal]
                    (eid_device, axialReal, axialImag, torsionReal,
                        torsionImag) = unpack(format1, edata)
                    eid = (eid_device - self.device_code) // 10

                    if is_magnitude_phase:
                        axial = polar_to_real_imag(axialReal, axialImag)
                        torsion = polar_to_real_imag(torsionReal, torsionImag)
                    else:
                        axial = complex(axialReal, axialImag)
                        torsion = complex(torsionReal, torsionImag)
                    assert eid > 0, eid

                    #print "out = ",out
                    #self.obj.add_new_eid(dt, eid, axial, torsion)
                    n += ntotal
            elif self.num_wide == 8:  # ???
                ntotal = 32
                format1 = b'i'
                nelements = len(data) // ntotal
                for i in xrange(nelements):
                    edata = data[n:n + 4]
                    eid_device, = unpack(format1, edata)
                    eid = (eid_device - self.device_code) // 10
                    assert eid > 0, eid
                    n += ntotal
            else:
                raise NotImplementedError(self.num_wide)

        elif self.element_type == 2: # CBEAM
            ntotal = 444 # 44 + 10*40  (11 nodes)
            #def getLengthTotal(self):
                #return 444  # 44+10*40   (11 nodes)

            #def getLength1(self):
                #return (44, 'ifffffffff')

            #def getLength2(self):
                #return (40, 'ifffffffff')

            nelements = len(data) // ntotal
            s = Struct(b'i')
            for i in xrange(nelements):
                edata = data[n:n+4]
                eid, = s.unpack(edata)
                n += ntotal

        elif self.element_type == 34: # CBAR
            if self.num_wide == 16:
                ntotal = 16 * 4
                nelements = len(data) // ntotal
                if self.debug:
                    self.binary_debug.write('  [cap, element1, element2, ..., cap]\n')
                    #self.binary_debug.write('  cap = %i  # assume 1 cap when there could have been multiple\n' % len(data))
                    self.binary_debug.write('  #elementi = [eid_device, s1a, s2a, s3a, s4a, axial, smaxa, smina, MSt,\n')
                    self.binary_debug.write('                           s1b, s2b, s3b, s4b, smaxb, sminb,        MSc]\n')
                    self.binary_debug.write('  nelements=%i; nnodes=1 # centroid\n' % nelements)

                s = Struct(b'i15f')
                for i in xrange(nelements):
                    edata = data[n:n+ntotal]
                    out = s.unpack(edata)
                    (eid_device, s1a, s2a, s3a, s4a, axial, smaxa, smina, MSt,
                                 s1b, s2b, s3b, s4b, smaxb, sminb, MSc) = out
                    eid = (eid_device - self.device_code) // 10
                    if self.debug4():
                        self.binary_debug.write('  ----------\n')
                        self.binary_debug.write('  eid = %i\n' % eid)
                        self.binary_debug.write('  C%i = [%s]\n' % (i, ', '.join(['%r' % di for di in out]) ))
                    n += ntotal
                    #print "eid =", eid
            elif self.num_wide == 19:
                format1 = b'i18f'

                ntotal = 76
                nelements = len(data) // ntotal
                for i in xrange(nelements):
                    edata = data[n:n+ntotal]
                    n += ntotal

                    (eid_device, s1ar, s2ar, s3ar, s4ar, axialr,
                     s1ai, s2ai, s3ai, s4ai, axiali,
                     s1br, s2br, s3br, s4br,
                     s1bi, s2bi, s3bi, s4bi) = unpack(format1, edata)

                    if is_magnitude_phase:
                        s1a = polar_to_real_imag(s1ar, s1ai)
                        s1b = polar_to_real_imag(s1br, s1bi)
                        s2a = polar_to_real_imag(s2ar, s2ai)
                        s2b = polar_to_real_imag(s2br, s2bi)
                        s3a = polar_to_real_imag(s3ar, s3ai)
                        s3b = polar_to_real_imag(s3br, s3bi)
                        s4a = polar_to_real_imag(s4ar, s4ai)
                        s4b = polar_to_real_imag(s4br, s4bi)
                        axial = polar_to_real_imag(axialr, axiali)

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

                    eid = (eid_device - self.device_code) // 10
                    #self.obj.add_new_eid('CBAR', dt, eid, s1a, s2a, s3a, s4a, axial,
                    #                                      s1b, s2b, s3b, s4b)
                    #print "eid=%i s1=%i s2=%i s3=%i s4=%i axial=%-5i" %(eid,s1a,s2a,s3a,s4a,axial)
                    #print "         s1=%i s2=%i s3=%i s4=%i"          %(s1b,s2b,s3b,s4b)
            else:
                raise NotImplementedError(self.num_wide)

        elif self.element_type in [11, 12, 13]: # CELAS1, CELAS2, CELAS3
            if self.num_wide == 2:  # real
                ntotal = 2 * 4
                nelements = len(data) // ntotal
                s = Struct(b'if')
                for i in xrange(nelements):
                    edata = data[n:n+ntotal]
                    out = s.unpack(edata)
                    (eid_device, ox) = out
                    eid = (eid_device - self.device_code) // 10
                    if self.debug4():
                        self.binary_debug.write('  ----------\n')
                        self.binary_debug.write('  eid = %i\n' % eid)
                        self.binary_debug.write('  result%i = [%i, %f]\n' % (i, eid_device, ox) )
                    n += ntotal
                    #print "eid =", eid
            elif self.num_wide == 3:
                ntotal = 12
                format1 = b'i2f'

                nelements = len(data) // ntotal
                for i in xrange(nelements):
                    edata = data[n:n + ntotal]
                    (eid_device, axialReal, axialImag) = unpack(format1, edata)
                    if is_magnitude_phase:
                        axial = polar_to_real_imag(axialReal, axialImag)
                    else:
                        axial = complex(axialReal, axialImag)

                    #print "out = ",out
                    eid = (eid_device - self.device_code) // 10
                    assert eid > 0
                    #self.obj.add_new_eid(dt, eid, axial)
                    n += ntotal
            else:
                raise NotImplementedError(self.num_wide)

        elif self.element_type in [39, 67, 68]: # TETRA, HEXA, PENTA
            # 67-CHEXA
            return

            if self.element_type == 39: # CTETRA
                nnodes_expected = 5  # 1 centroid + 4 corner points
            elif self.element_type == 67:
                nnodes_expected = 9
            elif self.element_type == 68:
                nnodes_expected = 7
            else:
                raise NotImplementedError('sort1 Type=%s num=%s' % (self.element_name, self.element_type))

            numwide_real = 4 + 21 * nnodes_expected

            if self.num_wide == numwide_real:
                ntotal = 16 + 84 * nnodes_expected
                #self.show_data(data[:ntotal])
                nelements = len(data) // ntotal
                #print "nelements =", nelements
                s = Struct(b'i')
                for i in xrange(nelements):
                    edata = data[n:n+ntotal]
                    eid_device, = s.unpack(data[n:n+4])
                    eid = (eid_device - self.device_code) // 10

                    if self.debug4():
                        self.binary_debug.write('  ----------\n')
                        self.binary_debug.write('  eid = %i\n' % eid)
                        #self.binary_debug.write('  C = [%s]\n' % ', '.join(['%r' % di for di in out]) )
                    #out = unpack(b'ii4si', edata)
                    n += ntotal
                    #print "eid =", eid
            else:
                raise NotImplementedError(self.num_wide)

        #=========================
        # plates
        elif self.element_type in [33]: # QUAD4-centroidal
            numwide_real = 17
            numwide_imag = 15
            if self.num_wide == numwide_real:
                return
                ntotal = 68  # 4*17
                format1 = b'i16f'
                s = Struct(format1)
                nelements = len(data) // ntotal
                for i in xrange(nelements):
                    edata = data[n:n+ntotal]
                    out = s.unpack(edata)
                    #if self.debug4():
                        #self.binary_debug.write('CQUAD4-33A - %s\n' % (str(out)))

                    (eid_device, fd1, sx1, sy1, txy1, angle1, major1, minor1, maxShear1,
                                 fd2, sx2, sy2, txy2, angle2, major2, minor2, maxShear2) = out

                    eid = (eid_device - self.device_code) // 10
                    if self.debug4():
                        self.binary_debug.write('  ----------\n')
                        self.binary_debug.write('  eid = %i\n' % eid)
                        self.binary_debug.write('  C = [%s]\n' % ', '.join(['%r' % di for di in out]) )

                    #print "eid=%i grid=%s fd1=%-3.1f sx1=%i sy1=%i txy1=%i angle1=%i major1=%i minor1=%i vm1=%i" % (eid,'C',fd1,sx1,sy1,txy1,angle1,major1,minor1,maxShear1)
                    #print   "             fd2=%-3.1f sx2=%i sy2=%i txy2=%i angle2=%i major2=%i minor2=%i vm2=%i\n"       % (fd2,sx2,sy2,txy2,angle2,major2,minor2,maxShear2)
                    #print "nNodes = ",nNodes
                    #self.obj.add_new_eid('CQUAD4', dt, eid, 'C', fd1, sx1, sy1,
                    #                   txy1, angle1, major1, minor1, maxShear1)
                    #self.obj.add(dt, eid, 'C', fd2, sx2, sy2, txy2,
                    #             angle2, major2, minor2, maxShear2)
                    #print "eid =", eid
                    n += ntotal
            elif self.num_wide == numwide_imag:
                format1 = b'i14f'
                nnodes = 0  # centroid + 4 corner points

                ntotal = 4 * (15 * (nnodes + 1))
                nelements = len(data) // ntotal
                for i in xrange(nelements):
                    edata = data[n:n+60]  # 4*15=60
                    n += 60
                    out = unpack(format1, edata)  # 15
                    if self.debug4():
                        self.binary_debug.write('  %s\n' % str(out))
                    (eid_device, fd1, sx1r, sx1i, sy1r, sy1i, txy1r, txy1i,
                                 fd2, sx2r, sx2i, sy2r, sy2i, txy2r, txy2i) = out

                    if is_magnitude_phase:
                        sx1 = polar_to_real_imag(sx1r, sx1i)
                        sx2 = polar_to_real_imag(sx2r, sx2i)
                        sy1 = polar_to_real_imag(sy1r, sy1i)
                        sy2 = polar_to_real_imag(sy2r, sy2i)
                        txy1 = polar_to_real_imag(txy1r, txy1i)
                        txy2 = polar_to_real_imag(txy2r, txy2i)
                    else:
                        sx1 = complex(sx1r, sx1i)
                        sx2 = complex(sx2r, sx2i)
                        sy1 = complex(sy1r, sy1i)
                        sy2 = complex(sy2r, sy2i)
                        txy1 = complex(txy1r, txy1i)
                        txy2 = complex(txy2r, txy2i)

                    eid = (eid_device - self.device_code) // 10

                    #print "eid=%i grid=%s fd1=%-3.1f sx1=%s sy1=%s txy1=%s" %(eid,'C',fd1,sx1,sy1,txy1)
                    #print   "             fd2=%-3.1f sx2=%s sy2=%s txy2=%s\n"       %(fd2,sx2,sy2,txy2)
                    #self.obj.add_new_eid('CQUAD4', dt, eid, 'C', fd1, sx1, sy1, txy1)
                    #self.obj.add(dt, eid, 'C', fd2, sx2, sy2, txy2)

                    for nodeID in xrange(nnodes):  # nodes pts
                        edata = self.data[n:n+60]  # 4*15=60
                        n += 60
                        out = unpack(b'i14f', edata)
                        if self.debug4():
                            self.binary_debug.write('  %s\n' % str(out))
                        (grid, fd1, sx1r, sx1i, sy1r, sy1i, txy1r, txy1i,
                               fd2, sx2r, sx2i, sy2r, sy2i, txy2r, txy2i) = out

                        if is_magnitude_phase:
                            sx1 = polar_to_real_imag(sx1r, sx1i)
                            sx2 = polar_to_real_imag(sx2r, sx2i)
                            sy1 = polar_to_real_imag(sy1r, sy1i)
                            sy2 = polar_to_real_imag(sy2r, sy2i)
                            txy1 = polar_to_real_imag(txy1r, txy1i)
                            txy2 = polar_to_real_imag(txy2r, txy2i)
                        else:
                            sx1 = complex(sx1r, sx1i)
                            sx2 = complex(sx2r, sx2i)
                            sy1 = complex(sy1r, sy1i)
                            sy2 = complex(sy2r, sy2i)
                            txy1 = complex(txy1r, txy1i)
                            txy2 = complex(txy2r, txy2i)

                        #print "eid=%i grid=%i fd1=%i sx1=%i sy1=%i txy1=%i\n" %(eid,grid,fd1,sx1,sy1,txy1)
                        #print "               fd2=%i sx2=%i sy2=%i txy2=%i\n"          %(fd2,sx2,sy2,txy2)
                        #self.obj.addNewNode(dt, eid, grid, fd1, sx1, sy1, txy1)
                        #self.obj.add(dt, eid, grid, fd2, sx2, sy2, txy2)

            else:
                raise NotImplementedError(self.num_wide)

        elif self.element_type in [74]: # TRIA3
            numwide_real = 17
            numwide_imag = 15
            if self.num_wide == numwide_real:
                return
                ntotal = 68  # 4*17
                #format1 = 'i16f'
                nelements = len(data) // ntotal
                s = Struct(b'i16f')
                if self.debug:
                    self.binary_debug.write('  [cap, element1, element2, ..., cap]\n')
                    self.binary_debug.write('  cap = %i  # assume 1 cap when there could have been multiple\n' % len(data))
                    self.binary_debug.write('  #elementi = [eid_device, fd1, sx1, sy1, txy1, angle1, major1, minor1, vm1,\n')
                    self.binary_debug.write('  #                        fd2, sx2, sy2, txy2, angle2, major2, minor2, vm2,]\n')
                    self.binary_debug.write('  nelements=%i; nnodes=1 # centroid\n' % nelements)

                for i in xrange(nelements):
                    edata = data[n:n + ntotal]
                    out = s.unpack(edata)
                    #if self.debug4():
                        #self.binary_debug.write('CTRIA3-74 - %s\n' % (str(out)))

                    (eid_device, fd1, sx1, sy1, txy1, angle1, major1, minor1, vm1,
                                 fd2, sx2, sy2, txy2, angle2, major2, minor2, vm2,) = out
                    eid = (eid_device - self.device_code) // 10

                    if self.debug4():
                        self.binary_debug.write('  ----------\n')
                        self.binary_debug.write('  eid = %i\n' % eid)
                        self.binary_debug.write('  centroid = [%s]\n' % ', '.join(['%r' % di for di in out]) )

                    #print "eid =", eid
                    n += ntotal
            if self.num_wide == numwide_imag:
                raise NotImplementedError(self.num_wide)
            else:
                raise NotImplementedError(self.num_wide)

        elif self.element_type in [64, 144, 70, 75, 82]:  # 64-cquad8/cquad4/70-ctriar/ctria6/cquadr
            # 64-CQUAD8
            # 70-CTRIAR
            # 75-CTRIA6
            # 82-CQUADR
            # 144-CQUAD4-bilinear
            if self.element_type in [64, 82, 144]:
                nnodes = 4 # + 1 centroid
            elif self.element_type in [70, 75]:
                nnodes = 3 # + 1 centroid
            else:
                raise NotImplementedError('sort1 Type=%s num=%s' % (self.element_name, self.element_type))

            numwide_real = 2 + 17 * (nnodes + 1)
            numwide_imag = 2 + 15 * (nnodes + 1)

            if self.num_wide == numwide_real:
                ntotal = 4 * (2 + 17 * (nnodes + 1))
                #assert ntotal == 348, ntotal
                center_format = b'i4si16f'
                node_format = b'i16f'
                cs = Struct(center_format)
                ns = Struct(node_format)

                nelements = len(data) // ntotal
                if self.debug:
                    self.binary_debug.write('  [cap, element1, element2, ..., cap]\n')
                    self.binary_debug.write('  cap = %i  # assume 1 cap when there could have been multiple\n' % len(data))
                    self.binary_debug.write('  #elementi = [centeri, node1i, node2i, node3i, node4i]\n')
                    self.binary_debug.write('  #centeri = [eid_device, j, grid, fd1, sx1, sy1, txy1, angle1, major1, minor1, vm1,\n')
                    self.binary_debug.write('  #                                fd2, sx2, sy2, txy2, angle2, major2, minor2, vm2,)]\n')
                    self.binary_debug.write('  #nodeji = [grid, fd1, sx1, sy1, txy1, angle1, major1, minor1, vm1,\n')
                    self.binary_debug.write('  #                fd2, sx2, sy2, txy2, angle2, major2, minor2, vm2,)]\n')
                    self.binary_debug.write('  nelements=%i; nnodes=%i # +1 centroid\n' % (nelements, nnodes))

                for i in xrange(nelements):
                    edata = data[n:n+76]

                    out = cs.unpack(edata)  # len=17*4
                    ## j is the number of nodes, so CQUAD4 -> 4, but we don't need to save it...
                    (eid_device, j, grid, fd1, sx1, sy1, txy1, angle1, major1, minor1, vm1,
                                          fd2, sx2, sy2, txy2, angle2, major2, minor2, vm2,) = out
                    eid = (eid_device - self.device_code) // 10

                    if self.debug4():
                        d = [eid_device, j, grid,
                             fd1, sx1, sy1, txy1, angle1, major1, minor1, vm1,
                             fd2, sx2, sy2, txy2, angle2, major2, minor2, vm2]
                        self.binary_debug.write('  ----------\n')
                        self.binary_debug.write('  eid = %i\n' % eid)
                        self.binary_debug.write('  C = [%s]\n' % ', '.join(['%r' % di for di in d]) )

                    #print "eid=%i grid=%s fd1=%i sx1=%i sy1=%i txy1=%i angle1=%i major1=%i minor1=%i vm1=%i" % (eid,grid,fd1,sx1,sy1,txy1,angle1,major1,minor1,vm1)
                    #print "               fd2=%i sx2=%i sy2=%i txy2=%i angle2=%i major2=%i minor2=%i vm2=%i\n"        % (fd2,sx2,sy2,txy2,angle2,major2,minor2,vm2)
                    n += 76
                    for inode in range(nnodes):
                        out = ns.unpack(data[n:n + 68])
                        if self.debug4():
                            d = tuple([grid,
                                      fd1, sx1, sy1, txy1, angle1, major1, minor1, vm1,
                                      fd2, sx2, sy2, txy2, angle2, major2, minor2, vm2])
                            self.binary_debug.write('  node%i = [%s]\n' % (inode+1, ', '.join(['%r' % di for di in d])))
                        (grid, fd1, sx1, sy1, txy1, angle1, major1, minor1, vm1,
                               fd2, sx2, sy2, txy2, angle2, major2, minor2, vm2,) = out
                        assert isinstance(grid, int), grid
                        assert grid > 0, grid

                        #print "eid=%i grid=%i fd1=%i sx1=%i sy1=%i txy1=%i angle1=%i major1=%i minor1=%i vm1=%i" % (eid,grid,fd1,sx1,sy1,txy1,angle1,major1,minor1,vm1)
                        #print "               fd2=%i sx2=%i sy2=%i txy2=%i angle2=%i major2=%i minor2=%i vm2=%i\n"        % (fd2,sx2,sy2,txy2,angle2,major2,minor2,vm2)
                        n += 68
            elif self.num_wide == numwide_imag:
                format1 = b'i14f' # 15
                ntotal = numwide_imag * 4
                assert self.num_wide * 4 == ntotal, 'numwide*4=%s ntotal=%s' % (self.num_wide*4, ntotal)
                nelements = len(data) // ntotal
                for i in xrange(nelements):
                    (eid_device, _) = unpack(b'i4s', data[n:n+8])
                    n += 8

                    eid = (eid_device - self.device_code) // 10
                    assert eid > 0
                    edata = data[n:n+60]  # 4*15
                    n += 60
                    out = unpack(format1, edata)  # len=15*4
                    if self.debug4():
                        self.binary_debug.write('%s\n' % (str(out)))
                    (grid, fd1, sx1r, sx1i, sy1r, sy1i, txy1r, txy1i,
                           fd2, sx2r, sx2i, sy2r, sy2i, txy2r, txy2i) = out
                    grid = 'C'

                    if is_magnitude_phase:
                        sx1 = polar_to_real_imag(sx1r, sx1i)
                        sy1 = polar_to_real_imag(sy1r, sy1i)
                        sx2 = polar_to_real_imag(sx2r, sx2i)
                        sy2 = polar_to_real_imag(sy2r, sy2i)
                        txy1 = polar_to_real_imag(txy1r, txy1i)
                        txy2 = polar_to_real_imag(txy2r, txy2i)
                    else:
                        sx1 = complex(sx1r, sx1i)
                        sy1 = complex(sy1r, sy1i)
                        sx2 = complex(sx2r, sx2i)
                        sy2 = complex(sy2r, sy2i)
                        txy1 = complex(txy1r, txy1i)
                        txy2 = complex(txy2r, txy2i)

                    #self.obj.add_new_eid(eType, dt, eid, grid, fd1, sx1, sy1, txy1)
                    #self.obj.add(dt, eid, grid, fd2, sx2, sy2, txy2)

                    for nodeID in xrange(nnodes):  # nodes pts
                        edata = data[n:n+60]  # 4*15=60
                        n += 60
                        out = unpack(b'i14f', edata)
                        if self.debug4():
                            self.binary_debug.write('%s\n' % (str(out)))
                        (grid, fd1, sx1r, sx1i, sy1r, sy1i, txy1r, txy1i,
                               fd2, sx2r, sx2i, sy2r, sy2i, txy2r, txy2i) = out

                        if is_magnitude_phase:
                            sx1 = polar_to_real_imag(sx1r, sx1i)
                            sx2 = polar_to_real_imag(sx2r, sx2i)
                            sy1 = polar_to_real_imag(sy1r, sy1i)
                            sy2 = polar_to_real_imag(sy2r, sy2i)
                            txy1 = polar_to_real_imag(txy1r, txy1i)
                            txy2 = polar_to_real_imag(txy2r, txy2i)
                        else:
                            sx1 = complex(sx1r, sx1i)
                            sx2 = complex(sx2r, sx2i)
                            sy1 = complex(sy1r, sy1i)
                            sy2 = complex(sy2r, sy2i)
                            txy1 = complex(txy1r, txy1i)
                            txy2 = complex(txy2r, txy2i)

                        #print "eid=%i grid=%i fd1=%i sx1=%i sy1=%i txy1=%i"   %(eid,grid,fd1,sx1,sy1,txy1)
                        #print "               fd2=%i sx2=%i sy2=%i txy2=%i\n"          %(fd2,sx2,sy2,txy2)
                        #self.obj.addNewNode(dt, eid, grid, fd1, sx1, sy1, txy1)
                        #self.obj.add(dt, eid, grid, fd2, sx2, sy2, txy2)
            else:
                raise NotImplementedError(self.num_wide)

        elif self.element_type in [88, 90]: # nonlinear shells
            # 88-CTRIA3NL
            # 90-CQUAD4NL

            if self.num_wide == 13:  # real
                ntotal = 52  # 4*13
                format1 = b'i12f'  # 1+12=13
                s = Struct(format1)
                nelements = len(data) // ntotal
                for i in xrange(nelements):
                    edata = data[n:n + ntotal]
                    out = s.unpack(edata)
                    if self.debug4():
                        self.binary_debug.write('CQUADNL-90 - %s\n' % str(out))

                    (eid_device, fd1, sx1, sy1, sz1, txy1, es1, eps1, ecs1,
                                      ex1, ey1, ez1, exy1) = out
                    eid = (eid_device - self.device_code) // 10
                    indata = (eid, fd1, sx1, sy1, sz1, txy1, es1, eps1, ecs1,
                                        ex1, ey1, ez1, exy1)
                    #self.obj.add_new_eid(self.element_type, dt, indata)
                    #print "eid=%s axial=%s equivStress=%s totalStrain=%s effPlasticCreepStrain=%s effCreepStrain=%s linearTorsionalStresss=%s" % (eid,axial,equivStress,totalStrain,effPlasticCreepStrain,effCreepStrain,linearTorsionalStresss)
                    n += ntotal
            elif self.num_wide == 25:  # complex
                ntotal = 100  # 4*25
                format1 = b'i24f'  # 1+24=25
                s = Struct(format1)
                nelements = len(data) // ntotal
                for i in xrange(nelements):
                    edata = data[n:n + ntotal]
                    out = s.unpack(edata)
                    if self.debug4():
                        self.binary_debug.write('CQUADNL-90 - %s\n' % str(out))
                    (eid_device, fd1, sx1, sy1, undef1, txy1, es1, eps1, ecs1, ex1, ey1, undef2, etxy1,
                                 fd2, sx2, sy2, undef3, txy2, es2, eps2, ecs2, ex2, ey2, undef4, etxy2) = out
                    eid = (eid_device - self.device_code) // 10
                    indata = (eid, fd1, sx1, sy1, undef1, txy1, es1, eps1, ecs1, ex1, ey1, undef2, etxy1,
                                   fd2, sx2, sy2, undef3, txy2, es2, eps2, ecs2, ex2, ey2, undef4, etxy2)
                    #self.obj.add_new_eid(self.element_type, dt, indata)
                    n += ntotal
            else:
                raise NotImplementedError(self.num_wide)

        elif self.element_type in [95, 96, 97, 98]: # composite shell
             # 95 - CQUAD4
             # 96 - CQUAD8
             # 97 - CTRIA3
             # 98 - CTRIA6 (composite)
            numwide_real = 11
            if self.num_wide == numwide_real:
                ntotal = 44
                nelements = len(data) // ntotal
                format1 = 'ii9f'
                s = Struct(format1)
                if self.debug:
                    self.binary_debug.write('  [cap, element1, element2, ..., cap]\n')
                    self.binary_debug.write('  cap = %i  # assume 1 cap when there could have been multiple\n' % len(data))
                    self.binary_debug.write('  #centeri = [eid_device, j, grid, fd1, sx1, sy1, txy1, angle1, major1, minor1, vm1,\n')
                    self.binary_debug.write('  #                                fd2, sx2, sy2, txy2, angle2, major2, minor2, vm2,)]\n')
                    self.binary_debug.write('  #nodeji = [eid, iLayer, o1, o2, t12, t1z, t2z, angle, major, minor, ovm)]\n')
                    self.binary_debug.write('  nelements=%i; nnodes=1 # centroid\n' % nelements)
                return

                eid_old = 0
                for i in xrange(nelements):
                    if i % 10000 == 0:
                        print 'i = ', i
                    edata = data[n:n+44]  # 4*11
                    out = s.unpack(edata)
                    (eid_device, iLayer, o1, o2, t12, t1z, t2z, angle, major, minor, ovm) = out
                    eid = (eid_device - self.device_code) // 10

                    if self.debug4():
                        self.binary_debug.write('  ----------\n')
                        self.binary_debug.write('  eid = %i\n' % eid)
                        self.binary_debug.write('  C = [%s]\n' % ', '.join(['%r' % di for di in out]) )

                    if eid != eid_old:  # originally initialized to None, the buffer doesnt reset it, so it is the old value
                        #print "1 - eid=%s iLayer=%i o1=%i o2=%i ovm=%i" % (eid,iLayer,o1,o2,ovm)
                        #self.obj.add_new_eid(eType, dt, eid, o1, o2, t12, t1z, t2z, angle, major, minor, ovm)
                        pass
                    else:
                        pass
                        #print "2 - eid=%s iLayer=%i o1=%i o2=%i ovm=%i" % (eid,iLayer,o1,o2,ovm)
                        #self.obj.add(dt, eid, o1, o2, t12, t1z, t2z, angle, major, minor, ovm)
                    eid_old = eid
                    n += 44
            else:
                raise NotImplementedError(self.num_wide)

        #=========================
        elif self.element_type in [53]: # axial plates
            # 53 - CTRIAX6
            if self.num_wide == 33: # real
                format1 = b'2i7f'  # 36

                ntotal = 132  # (1+8*4)*4 = 33*4 = 132
                nelements = len(data) // ntotal
                n = 0
                s = Struct(format1)
                ns = Struct(b'i7f')
                for i in xrange(nelements):
                    out = s.unpack(data[n:n + 36])
                    (eid_device, loc, rs, azs, As, ss, maxp, tmax, octs) = out
                    if self.debug4():
                        self.binary_debug.write('CTRIAX6-53A - %s\n' % (str(out)))
                    eid = (eid_device - self.device_code) // 10
                    #print "eid=%s loc=%s rs=%s azs=%s as=%s ss=%s maxp=%s tmx=%s octs=%s" % (eid,loc,rs,azs,As,ss,maxp,tmax,octs)
                    #self.obj.add_new_eid(dt, eid, loc, rs, azs, As, ss, maxp, tmax, octs)
                    n += 36
                    for i in xrange(3):
                        out = ns.unpack(data[n:n + 32])
                        (loc, rs, azs, As, ss, maxp, tmax, octs) = out
                        if self.debug4():
                            self.binary_debug.write('CTRIAX6-53B - %s\n' % (str(out)))
                        #print "eid=%s loc=%s rs=%s azs=%s as=%s ss=%s maxp=%s tmx=%s octs=%s" % (eid,loc,rs,azs,As,ss,maxp,tmax,octs)
                        #self.obj.add(dt, eid, loc, rs, azs, As, ss, maxp, tmax, octs)
                        n += 32  # 4*8
            elif self.num_wide == 37: # imag
                raise NotImplementedError(self.num_wide)
            else:
                raise NotImplementedError(self.num_wide)
        elif self.element_type in [102]: # bush
            # 102-CBUSH
            numwide_real = 7
            if self.num_wide == numwide_real:
                return
                assert self.num_wide == 7, "num_wide=%s not 7" % self.num_wide
                ntotal = 28  # 4*7
                format1 = b'i6f'

                n = 0
                nelements = len(data) // ntotal
                s = Struct(format1)
                for i in xrange(nelements):
                    edata = data[n:n + ntotal]
                    out = s.unpack(edata)  # num_wide=7
                    if self.debug4():
                        self.binary_debug.write('CBUSH-102 - %s\n' % str(out))

                    (eid_device, tx, ty, tz, rx, ry, rz) = out
                    eid = (eid_device - self.device_code) // 10

                    #self.obj.add_new_eid(self.element_type, dt, eid, tx, ty, tz, rx, ry, rz)
                    n += ntotal
            elif self.num_wide == 13:  # imag
                ntotal = 52  # 4*13
                format1 = b'i12f'

                n = 0
                nelements = len(data) // ntotal
                s = Struct(format1)
                for i in xrange(nelements):
                    edata = data[n:n + ntotal]
                    out = s.unpack(edata)  # num_wide=7
                    if self.debug4():
                        self.binary_debug.write('CBUSH-102 - %s\n' % str(out))

                    (eid_device, txr, tyr, tzr, rxr, ryr, rzr,
                                 txi, tyi, tzi, rxi ,ryi, rzi) = out
                    eid = (eid_device - self.device_code) // 10

                    if is_magnitude_phase:
                        tx = polar_to_real_imag(txr, txi)
                        ty = polar_to_real_imag(tyr, tyi)
                        tz = polar_to_real_imag(tzr, tzi)
                        rx = polar_to_real_imag(rxr, rxi)
                        ry = polar_to_real_imag(ryr, ryi)
                        rz = polar_to_real_imag(rzr, rzi)
                    else:
                        tx = complex(txr, txi)
                        ty = complex(tyr, tyi)
                        tz = complex(tzr, tzi)
                        rx = complex(rxr, rxi)
                        ry = complex(ryr, ryi)
                        rz = complex(rzr, rzi)
                    #self.obj.add_new_eid(self.element_type, dt, eid, tx, ty, tz, rx, ry, rz)
                    n += ntotal
            else:
                raise NotImplementedError(self.num_wide)

        elif self.element_type in [40]:  # bush
            # 40-CBUSH1D
            numwide_real = 8
            if self.num_wide == numwide_real:
                assert self.num_wide == 8, "num_wide=%s not 8" % self.num_wide
                ntotal = 32  # 4*8
                format1 = b'i6fi'

                n = 0
                s = Struct(format1)
                nelements = len(data) // ntotal
                for i in xrange(nelements):
                    edata = data[n:n + ntotal]

                    out = s.unpack(edata)  # num_wide=25
                    if self.debug4():
                        self.binary_debug.write('CBUSH1D-40 - %s\n' % (str(out)))
                    (eid_device, fe, ue, ve, ao, ae, ep, fail) = out
                    eid = (eid_device - self.device_code) // 10

                    # axial_force, axial_displacement, axial_velocity, axial_stress, axial_strain, plastic_strain, is_failed
                    #self.obj.add_new_eid(self.element_type, dt, eid, fe, ue, ve, ao, ae, ep, fail)
                    n += ntotal
        elif self.element_type in [82]:  # cquadr
            # 82-CQUADR
            #GRID-ID  DISTANCE,NORMAL-X,NORMAL-Y,SHEAR-XY,ANGLE,MAJOR MINOR,VONMISES

            if self.element_type == 82:  # CQUADR
                ntotal = 348  # 2+17*5 = 87 -> 87*4 = 348
                nNodes = 4    # centroid + 4 corner points
                eType = 'CQUADR'
            else:
                raise RuntimeError('element_type=%s ntotal not defined...' % self.element_type)

            numwide_real = 2 + 17*(nnodes + 1)
            if self.num_wide == numwide_real:
                ntotal = numwide_real * 4
                format1 = b'i16f'  # 1+16 = 17 * 4 = 68
                s = Struct(format1)
                n = 0
                nelements = len(data) // ntotal
                for i in xrange(nelements):
                    hdata = data[n:n+8]
                    n += 8
                    (eid_device, _) = unpack(b'i4s', hdata)
                    eid = (eid_device - self.device_code) // 10

                    edata = data[n:n+68]  # 4*17
                    n += 68
                    out = s.unpack(edata)  # len=17*4
                    if self.debug4():
                        self.binary_debug.write('CQUADR-82A - %s\n' % str(out))
                    (grid, fd1, sx1, sy1, txy1, angle1, major1, minor1, vm1,
                           fd2, sx2, sy2, txy2, angle2, major2, minor2, vm2,) = out
                    grid = 'C'
                    #self.obj.add_new_eid(eType, eid, grid, fd1, sx1, sy1,
                    #                   txy1, angle1, major1, minor1, vm1)
                    #self.obj.add(eid, grid, fd2, sx2, sy2, txy2,
                    #             angle2, major2, minor2, vm2)

                    for nodeID in xrange(nNodes):  # nodes pts
                        edata = data[n:n+68]
                        out = s.unpack(edata)
                        if self.debug4():
                            self.binary_debug.write('CQUADR-82B - %s\n' % str(out))
                        (grid, fd1, sx1, sy1, txy1, angle1, major1, minor1, vm1,
                               fd2, sx2, sy2, txy2, angle2, major2, minor2, vm2,) = out

                        #print "eid=%i grid=%i fd1=%i sx1=%i sy1=%i txy1=%i angle1=%i major1=%i minor1=%i vm1=%i" % (eid,grid,fd1,sx1,sy1,txy1,angle1,major1,minor1,vm1)
                        #print "               fd2=%i sx2=%i sy2=%i txy2=%i angle2=%i major2=%i minor2=%i vm2=%i\n"          % (fd2,sx2,sy2,txy2,angle2,major2,minor2,vm2)
                        #self.obj.addNewNode(eid, grid, fd1, sx1,
                        #                    sy1, txy1, angle1, major1, minor1, vm1)
                        #self.obj.add(eid, grid, fd2, sx2, sy2,
                        #             txy2, angle2, major2, minor2, vm2)
                        n+= 68
            else:
                raise NotImplementedError(self.num_wide)

        elif self.element_type in [87, 89, 92]:  # nonlinear rods
            # 87-CTUBENL
            # 89-RODNL
            # 92-CONRODNL
            #resultName = self.makeOES_Object(self.nonlinearRodStress, NonlinearRodObject, 'nonlinearRodStress',
            #                                 self.nonlinearRodStrain, NonlinearRodObject, 'nonlinearRodStrain')
            numwide_real = 7
            if self.num_wide == numwide_real:
                format1 = b'i6f'  # 1+6=7
                s = Struct(format1)
                ntotal = 28  # len(format1)*4 = 7*4 = 28
                nelements = len(data) // ntotal
                for i in xrange(nelements):
                    edata = data[n:n+ntotal]
                    out = s.unpack(edata)

                    (eid_device, axial, equivStress, totalStrain, effPlasticCreepStrain,
                        effCreepStrain, linearTorsionalStresss) = out
                    eid = (eid_device - self.device_code) // 10
                    if self.debug4():
                        self.binary_debug.write('CRODNL-%s - %s\n' % (self.element_type, str(out)))
                    indata = (eid, axial, equivStress, totalStrain, effPlasticCreepStrain, effCreepStrain, linearTorsionalStresss)
                    #print "eid=%s axial=%s equivStress=%s totalStrain=%s effPlasticCreepStrain=%s effCreepStrain=%s linearTorsionalStresss=%s" % (eid,axial,equivStress,totalStrain,effPlasticCreepStrain,effCreepStrain,linearTorsionalStresss)
                    #self.obj.add(self.element_type, dt, indata)
                    n += ntotal
            else:
                raise NotImplementedError(self.num_wide)

        elif self.element_type in [224, 225]: # nonlinear spring
            # 224-CELAS1
            # 225-CELAS3
            # nonlinearSpringStress
            numwide_real = 3
            if self.num_wide == numwide_real:
                assert self.num_wide == 3, "num_wide=%s not 3" % self.num_wide
                format1 = b'i2f'
                ntotal = 12  # 4*3
                nelements = len(data) // ntotal
                s = Struct(format1)
                for i in xrange(nelements):
                    edata = data[0:ntotal]
                    out = s.unpack(edata)  # num_wide=3
                    (eid_device, force, stress) = out
                    eid = (eid_device - self.device_code) // 10
                    if self.debug4():
                        self.binary_debug.write('%s-%s - %s\n' % (self.element_name, self.element_type, str(out)))
                    #self.obj.add_new_eid(element_name, dt, eid, force, stress)
            else:
                raise NotImplementedError(self.num_wide)

        elif self.element_type in [4]: # shear
            # 4-CSHEAR
            if self.num_wide == 1:
                asdf
            elif self.num_wide == 5:
                ntotal = 20  # 4*5
                format1 = b'i4f'

                n = 0
                nelements = len(data) // ntotal
                s = Struct(format1)
                for i in xrange(nelements):
                    edata = data[n:n + ntotal]
                    out = s.unpack(edata)  # num_wide=5
                    if self.debug4():
                        self.binary_debug.write('CSHEAR-4 - %s\n' % str(out))

                    (eid_device, etmaxr, etmaxi, etavgr, etavgi) = out
                    eid = (eid_device - self.device_code) // 10

                    if is_magnitude_phase:
                        etmax = polar_to_real_imag(etmaxr, etmaxi)
                        etavg = polar_to_real_imag(etavgr, etavgi)
                    else:
                        etmax = complex(etmaxr, etmaxi)
                        etavg = complex(etavgr, etavgi)
                    #self.obj.add_new_eid(self.element_type, dt, eid, etmax, etavg)
                    n += ntotal
            else:
                raise NotImplementedError(self.num_wide)
        elif self.element_type in [35]:
            # 35-CON
            return
        elif self.element_type in [60, 61]:
            # 60-DUM8
            # 61-DUM9
            return
        elif self.element_type in [69]:
            # 69-CBEND
            return
        elif self.element_type in [86]:
            # 86-GAPNL
            return
        elif self.element_type in [94]:
            # 94-BEAMNL
            return
        elif self.element_type in [85, 91, 93]:
            # 85-TETRANL
            # 91-PENTANL
            # 93-HEXANL
            return
        elif self.element_type in [100]:  # bars
            # 100-BARS
            return
        elif self.element_type in [101]:
            # 101-AABSF
            return
        elif self.element_type in [140, 145, 146, 147, 201]:
            # 140-HEXA8FD, 201-QUAD4FD
            # 145-VUHEXA
            # 146-VUPENTA
            # 147-VUTETRA
            return
        elif self.element_type in [47, 48, 139, 189, 190]:
            # 47-AXIF2
            # 48-AXIF3
            # 139-QUAD4FD
            # 189-VUQUAD
            # 190-VUTRIA
            return
        elif self.element_type in [191]:
            # 191-VUBEAM
            return
        elif self.element_type in [50, 51, 203]:
            # 203-SLIF1D?
            # 50-SLOT3
            # 51-SLOT4
            return
        elif self.element_type in [160, 161, 162, 163, 164, 165, 166, 167, 168,
                                   169, 170, 171, 172, 202,
                                   204, 218, 211, 213, 214, 
                                   216, 217, 219, 220, 221, 222, 223,
                                   226, 232, 233, 235]:
            # 160-PENTA6FD
            # 161-TETRA4FD
            # 162-TRIA3FD
            # 163-HEXAFD
            # 164-QUADFD
            # 165-PENTAFD
            # 166-TETRAFD
            # 167-TRIAFD
            # 168-TRIAX3FD
            # 169-TRIAXFD
            # 170-QUADX4FD
            # 171-QUADXFD
            # 172-QUADRNL
            # 202-HEXA8FD
            # 204-PENTA6FD
            # 211-TRIAFD
            # 213-TRIAXFD
            # 214-QUADX4FD
            # 216-TETRA4FD
            # 217-TRIA3FD
            # 218-HEXAFD
            # 219-QUADFD
            # 220-PENTAFD
            # 221-TETRAFD
            # 223-QUADXFD
            # 222-TRIAX3FD
            # 226-BUSH
            # 232-QUADRLC
            # 233-TRIARLC
            # 235-CQUADR
            return
        else:
            raise NotImplementedError('sort1 Type=%s num=%s' % (self.element_name, self.element_type))

        assert len(data) > 0
        assert nelements > 0, 'nelements=%r element_type=%s element_name=%r' % (nelements, self.element_type, self.element_name)
        assert len(data) % ntotal == 0, '%s n=%s len=%s ntotal=%s' % (self.element_name, len(data) % ntotal, len(data), ntotal)
        assert self.num_wide * 4 == ntotal, 'numwide*4=%s ntotal=%s' % (self.num_wide*4, ntotal)
        assert self.thermal == 0, self.thermal