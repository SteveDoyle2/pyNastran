from struct import unpack

class OES(object):
    def __init__(self):
        self.oes_mapper = {
            1: 'ROD',
            2 : 'CBEAM',
            10: 'CONROD',
            11: 'CELAS1',
            12: 'CELAS2',
            13: 'CELAS3',
            33: 'QUAD4', # 1 node
            34 : 'CBAR',
            39: 'TETRA',
            67: 'HEXA',
            68: 'PENTA',
            74: 'TRIA3', # 1 node
            144: 'QUAD144', # 5 nodes
        }
        pass

    def read_oes1x1_4(self, data):
        #assert self.isubtable == -4, self.isubtable
        self.element_name = self.oes_mapper[self.element_type]

        if self.debug:
            self.binary_debug.write('  Table3\n')
            self.binary_debug.write('  element_name = %r\n' % self.element_name)

        if self.is_sort1():
            self.read_oes1x1_4_sort1(data)
        else:
            raise NotImplementedError('sort2 Type=%s num=%s' % (self.element_name, self.element_type))
        print "element_name =", self.element_name

        if self.debug:
            self.binary_debug.write('\n')

        print "----------------"
        #sys.exit('stopping...')

    def read_oes1x1_4_sort1(self, data):
        print "element_name =", self.element_name
        n = 0
        if self.element_type in [1, 10]: # CROD
            ntotal = 5 * 4
            nelements = len(data) // ntotal
            if self.debug:
                self.binary_debug.write('[cap, element1, element2, ..., cap]\n')
                self.binary_debug.write('  cap = %i  # assume 1 cap when there could have been multiple\n' % len(data))
                self.binary_debug.write('  #elementi = [eid_device, axial, axial_margin, torsion, torsion_margin]\n')
                self.binary_debug.write('  nelements=%i; nnodes=1 # centroid\n' % nelements)
            for i in xrange(nelements):
                edata = data[n:n+ntotal]
                out = unpack('i4f', edata)
                (eid_device, axial, axial_margin, torsion, torsion_margin) = out
                eid = (eid_device - self.device_code) // 10
                if self.debug:
                    self.binary_debug.write('  ----------\n')
                    self.binary_debug.write('  eid = %i\n' % eid)
                    self.binary_debug.write('  C = [%s]\n' % ', '.join(['%r' % di for di in out]) )

                print "eid =", eid
                n += ntotal

        elif self.element_type == 2: # CBEAM
            ntotal = 444 # 44 + 10*40  (11 nodes)
            #def getLengthTotal(self):
                #return 444  # 44+10*40   (11 nodes)

            #def getLength1(self):
                #return (44, 'ifffffffff')

            #def getLength2(self):
                #return (40, 'ifffffffff')

            nelements = len(data) // ntotal
            for i in xrange(nelements):
                edata = data[n:n+4]
                eid, = unpack('i', edata)
                n += ntotal

        elif self.element_type == 34: # CBAR
            ntotal = 16 * 4
            nelements = len(data) // ntotal
            if self.debug:
                self.binary_debug.write('[cap, element1, element2, ..., cap]\n')
                #self.binary_debug.write('  cap = %i  # assume 1 cap when there could have been multiple\n' % len(data))
                self.binary_debug.write('  #elementi = [eid_device, s1a, s2a, s3a, s4a, axial, smaxa, smina, MSt,\n')
                self.binary_debug.write('                           s1b, s2b, s3b, s4b, smaxb, sminb,        MSc]\n')
                self.binary_debug.write('  nelements=%i; nnodes=1 # centroid\n' % nelements)
            for i in xrange(nelements):
                edata = data[n:n+ntotal]
                out = unpack('i15f', edata)
                (eid_device, s1a, s2a, s3a, s4a, axial, smaxa, smina, MSt,
                             s1b, s2b, s3b, s4b, smaxb, sminb, MSc) = out
                eid = (eid_device - self.device_code) // 10
                if self.debug:
                    self.binary_debug.write('  ----------\n')
                    self.binary_debug.write('  eid = %i\n' % eid)
                    self.binary_debug.write('  C%i = [%s]\n' % (i, ', '.join(['%r' % di for di in out]) ))
                n += ntotal
                print "eid =", eid
        elif self.element_type in [39, 67, 68]: # TETRA, HEXA, PENTA
            if self.element_type == 39: # CTETRA
                nnodes_expected = 5  # 1 centroid + 4 corner points
            elif self.element_type == 67:
                nnodes_expected = 9
            elif self.element_type == 68:
                nnodes_expected = 7
            else:
                raise NotImplementedError('sort1 Type=%s num=%s' % (self.element_name, self.element_type))

            ntotal = 16 + 84 * nnodes_expected
            #self.show_data(data[:ntotal])
            nelements = len(data) // ntotal
            print "nelements =", nelements
            for i in xrange(nelements):
                edata = data[n:n+ntotal]
                eid_device, = unpack('i', data[n:n+4])
                eid = (eid_device - self.device_code) // 10

                if self.debug:
                    self.binary_debug.write('  ----------\n')
                    self.binary_debug.write('  eid = %i\n' % eid)
                    #self.binary_debug.write('  C = [%s]\n' % ', '.join(['%r' % di for di in out]) )
                #out = unpack('ii4si', edata)
                n += ntotal
                print "eid =", eid

        #=========================
        elif self.element_type in [33]: # QUAD4-centroidal
            ntotal = 68  # 4*17
            format1 = 'i16f'
            nelements = len(data) // ntotal
            for i in xrange(nelements):
                edata = data[n:n+ntotal]
                out = unpack(format1, edata)
                #if self.make_op2_debug:
                    #self.op2_debug.write('CQUAD4-33A - %s\n' % (str(out)))

                (eid_device, fd1, sx1, sy1, txy1, angle1, major1, minor1, maxShear1,
                             fd2, sx2, sy2, txy2, angle2, major2, minor2, maxShear2) = out

                eid = (eid_device - self.device_code) // 10
                if self.debug:
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

        elif self.element_type in [74]: # TRIA3
            ntotal = 68  # 4*17
            format1 = 'i16f'
            nelements = len(data) // ntotal
            if self.debug:
                self.binary_debug.write('[cap, element1, element2, ..., cap]\n')
                self.binary_debug.write('  cap = %i  # assume 1 cap when there could have been multiple\n' % len(data))
                self.binary_debug.write('  #elementi = [eid_device, fd1, sx1, sy1, txy1, angle1, major1, minor1, vm1,\n')
                self.binary_debug.write('                           fd2, sx2, sy2, txy2, angle2, major2, minor2, vm2,]\n')
                self.binary_debug.write('  nelements=%i; nnodes=1 # centroid\n' % (nelements))

            for i in xrange(nelements):
                edata = data[n:n + ntotal]
                out = unpack(format1, edata)
                #if self.make_op2_debug:
                    #self.op2_debug.write('CTRIA3-74 - %s\n' % (str(out)))

                (eid_device, fd1, sx1, sy1, txy1, angle1, major1, minor1, vm1,
                             fd2, sx2, sy2, txy2, angle2, major2, minor2, vm2,) = out
                eid = (eid_device - self.device_code) // 10

                if self.debug:
                    #d = [eid_device, j, grid,
                         #fd1, sx1, sy1, txy1, angle1, major1, minor1, vm1,
                         #fd2, sx2, sy2, txy2, angle2, major2, minor2, vm2]
                    #print "d =", d
                    self.binary_debug.write('  ----------\n')
                    self.binary_debug.write('  eid = %i\n' % eid)
                    self.binary_debug.write('  centroid = [%s]\n' % ', '.join(['%r' % di for di in out]) )

                print "eid =", eid
                n += ntotal
        elif self.element_type in [144]: # CQUAD4-bilinear
            if self.element_type == 144:
                nnodes = 4 # + 1 centroid
            else:
                raise NotImplementedError('sort1 Type=%s num=%s' % (self.element_name, self.element_type))
            #ntotal = 16 + 84 * nnodes_expected
            ntotal = 4 * (2 + 17 * (nnodes + 1))
            assert ntotal == 348, ntotal
            nelements = len(data) // ntotal
            center_format = 'i4si16f'
            node_format = 'i16f'
            if self.debug:
                self.binary_debug.write('[cap, element1, element2, ..., cap]\n')
                self.binary_debug.write('  cap = %i  # assume 1 cap when there could have been multiple\n' % len(data))
                self.binary_debug.write('  #elementi = [centeri, node1i, node2i, node3i, node4i]\n')
                self.binary_debug.write('  #centeri = [eid_device, j, grid, fd1, sx1, sy1, txy1, angle1, major1, minor1, vm1,\n')
                self.binary_debug.write('  #                                fd2, sx2, sy2, txy2, angle2, major2, minor2, vm2,)]\n')
                self.binary_debug.write('  #nodeji = [grid, fd1, sx1, sy1, txy1, angle1, major1, minor1, vm1,\n')
                self.binary_debug.write('  #                fd2, sx2, sy2, txy2, angle2, major2, minor2, vm2,)]\n')
                self.binary_debug.write('  nelements=%i; nnodes=%i # +1 centroid\n' % (nelements, nnodes))

            for i in xrange(nelements):
                edata = data[n:n+76]
                #eid_device, = unpack('i', data[n:n+4])

                out = unpack(center_format, edata)  # len=17*4
                ## j is the number of nodes, so CQUAD4 -> 4, but we don't need to save it...
                (eid_device, j, grid, fd1, sx1, sy1, txy1, angle1, major1, minor1, vm1,
                                      fd2, sx2, sy2, txy2, angle2, major2, minor2, vm2,) = out
                eid = (eid_device - self.device_code) // 10

                if self.debug:
                    d = [eid_device, j, grid,
                         fd1, sx1, sy1, txy1, angle1, major1, minor1, vm1,
                         fd2, sx2, sy2, txy2, angle2, major2, minor2, vm2]
                    print "d =", d
                    self.binary_debug.write('  ----------\n')
                    self.binary_debug.write('  eid = %i\n' % eid)
                    self.binary_debug.write('  C = [%s]\n' % ', '.join(['%r' % di for di in d]) )

                #print "eid=%i grid=%s fd1=%i sx1=%i sy1=%i txy1=%i angle1=%i major1=%i minor1=%i vm1=%i" % (eid,grid,fd1,sx1,sy1,txy1,angle1,major1,minor1,vm1)
                #print "               fd2=%i sx2=%i sy2=%i txy2=%i angle2=%i major2=%i minor2=%i vm2=%i\n"        % (fd2,sx2,sy2,txy2,angle2,major2,minor2,vm2)
                n += 76
                for inode in range(nnodes):
                    out = unpack(node_format, data[n:n + 68])
                    if self.debug:
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
        else:
            raise NotImplementedError('sort1 Type=%s num=%s' % (self.element_name, self.element_type))
        assert nelements > 0, nelements
        assert len(data) % ntotal == 0, '%s n=%s len=%s ntotal=%s' % (self.element_name, len(data) % ntotal, len(data), ntotal)


    def apply_data_code_value(self, name, value):
        pass
    def setNullNonlinearFactor(self):
        pass

    def is_sort1(self):
        if self.sort_bits[0] == 0:
            return True
        return False

    def is_sort2(self):
        return not(self.is_sort1())

    def read_oes1x1_3(self, data):
        self.words = ['aCode',       'tCode',    'element_type', 'isubcase',
                 '???',         '???',      '???',          'load_set'
                 'format_code', 'num_wide', 's_code',       '???',
                 '???',         '???',      '???',          '???',
                 '???',         '???',      '???',          '???',
                 '???',         '???',      '???',          '???',
                 '???', 'Title', 'subtitle', 'label']

        #self.show_data(data)
        self.parse_approach_code(data)  # 3
        ## element type
        self.add_data_parameter(data, 'element_type', 'i', 3, False)
        ## load set ID
        self.add_data_parameter(data, 'load_set', 'i', 8, False)
        ## format code
        self.add_data_parameter(data, 'format_code', 'i', 9, False)
        ## number of words per entry in record
        ## .. note:: is this needed for this table ???
        self.add_data_parameter(data, 'num_wide', 'i', 10, False)
        ## stress/strain codes
        self.add_data_parameter(data, 's_code', 'i', 11, False)
        ## thermal flag; 1 for heat ransfer, 0 otherwise
        self.add_data_parameter(data, 'thermal', 'i', 23, False)

        ## assuming tCode=1
        if self.analysis_code == 1:   # statics / displacement / heat flux
            ## load set number
            self.add_data_parameter(data, 'lsdvmn', 'i', 5, False)
            self.apply_data_code_value('dataNames', ['lsdvmn'])
            self.setNullNonlinearFactor()
        elif self.analysis_code == 2:  # real eigenvalues
            ## mode number
            self.add_data_parameter(data, 'mode', 'i', 5)
            ## real eigenvalue            self.add_data_parameter(data, 'eign', 'f', 6, False)
            ## mode or cycle TODO confused on the type - F1???
            self.add_data_parameter(data, 'mode_cycle', 'f', 7, False)
            self.apply_data_code_value('dataNames', ['mode', 'eigr', 'mode_cycle'])
        #elif self.analysis_code==3: # differential stiffness
            #self.lsdvmn = self.get_values(data,'i',5) ## load set number
            #self.data_code['lsdvmn'] = self.lsdvmn
        #elif self.analysis_code==4: # differential stiffness
        #    self.lsdvmn = self.get_values(data,'i',5) ## load set number
        elif self.analysis_code == 5:   # frequency
            ## frequency
            self.add_data_parameter(data, 'freq', 'f', 5)
            self.apply_data_code_value('dataNames', ['freq'])
        elif self.analysis_code == 6:  # transient
            ## time step
            self.add_data_parameter(data, 'dt', 'f', 5)
            self.apply_data_code_value('dataNames', ['dt'])
        elif self.analysis_code == 7:  # pre-buckling
            ## load set
            self.add_data_parameter(data, 'lsdvmn', 'i', 5)
            self.apply_data_code_value('dataNames', ['lsdvmn'])
        elif self.analysis_code == 8:  # post-buckling
            ## mode number
            self.add_data_parameter(data, 'lsdvmn', 'i', 5)
            self.add_data_parameter(data, 'eigr', 'f', 6, False)  # real eigenvalue
            self.apply_data_code_value('dataNames', ['lsdvmn', 'eigr'])
        elif self.analysis_code == 9:  # complex eigenvalues
            ## mode number
            self.add_data_parameter(data, 'mode', 'i', 5)
            ## real eigenvalue
            self.add_data_parameter(data, 'eigr', 'f', 6, False)
            ## imaginary eigenvalue            self.add_data_parameter(data, 'eigi', 'f', 7, False)
            self.apply_data_code_value('dataNames', ['mode', 'eigr', 'eigi'])
        elif self.analysis_code == 10:  # nonlinear statics
            ## load step
            self.add_data_parameter(data, 'lftsfq', 'f', 5)
            self.apply_data_code_value('dataNames', ['lftsfq'])
        elif self.analysis_code == 11:  # old geometric nonlinear statics
            ## load set number
            self.add_data_parameter(data, 'lsdvmn', 'i', 5)
            self.apply_data_code_value('dataNames', ['lsdvmn'])
        elif self.analysis_code == 12:  # contran ? (may appear as aCode=6)  --> straight from DMAP...grrr...
            ## Time step ??? --> straight from DMAP
            self.add_data_parameter(data, 'dt', 'f', 5)
            self.apply_data_code_value('dataNames', ['dt'])
        else:
            raise RuntimeError('invalid analysis_code...analysis_code=%s' %
                               self.analysis_code)
        # tCode=2
        #if self.analysis_code==2: # sort2
        #    self.lsdvmn = self.get_values(data,'i',5)

        if not self.is_sort1() and not isRelease:
            raise NotImplementedError('sort2...')

        if self.debug:
            self.binary_debug.write('  element_name = %r\n' % self.oes_mapper[self.element_type])
            self.binary_debug.write('  aCode    = %r\n' % self.aCode)
            self.binary_debug.write('  tCode    = %r\n' % self.tCode)
            self.binary_debug.write('  isubcase = %r\n' % self.isubcase)

        #self.show_data(data[200:])
        self.read_title(data)
        if self.element_type not in self.oes_mapper:
            raise NotImplementedError(self.element_type)

        if self.debug:
            msg = ''
            for i, param in enumerate(self.words):
                if param == '???':
                    param = 0
                msg += '%s, ' % param
                if i % 5 == 4:
                    msg += '\n             '

            self.binary_debug.write('  recordi = [%s]\n\n' % msg)
        #sys.exit()

    def read_title(self, data):
        assert len(data) == 584, len(data)
        Title, subtitle, label = unpack('128s128s128s', data[200:])  # titleSubtitleLabel

        self.Title = Title.strip()
        #: the subtitle of the subcase
        self.subtitle = subtitle.strip()
        #: the label of the subcase
        self.label = label.strip()
        if self.debug:
            self.binary_debug.write('  Title    = %r\n' % self.Title)
            self.binary_debug.write('  subtitle = %r\n' % self.subtitle)
            self.binary_debug.write('  label    = %r\n' % self.label)

    def finish_oes(self):
        del self.element_type
        del self.load_set
        del self.format_code
        del self.num_wide
        del self.s_code
        del self.thermal
