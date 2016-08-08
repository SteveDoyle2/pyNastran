from six.moves import zip, range
from numpy import array

from pyNastran.op2.result_objects.table_object import RealTableArray, ComplexTableArray
from pyNastran.f06.f06_formatting import write_floats_13e


class ComplexEigenvectorArray(ComplexTableArray):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        ComplexTableArray.__init__(self, data_code, is_sort1, isubcase, dt)

    def write_f06(self, f, header=None, page_stamp='PAGE %s', page_num=1, is_mag_phase=False, is_sort1=True):
        if header is None:
            header = []

        words = []
        #freq = self.eigrs[i]
        #freq = 0.0
        #words.append('%16s = %12E\n' % ('EIGENVALUE', freq))
        #words.append('%%16s = %%12E\n')
        has_cycle = False
        if hasattr(self, 'cycles'):
            has_cycle = True
        #if has_cycle:
        #words.append('%s = %s          C O M P L E X   E I G E N V E C T O R   N O. %s\n' % ('%16s', '%12E', '%10i'))
            #msg.append('%16s = %12E          C O M P L E X   E I G E N V E C T O R   N O . %10i\n \n' % ('CYCLES', self.mode_cycle, imode))
        #else:
        words.append('      COMPLEX EIGENVALUE = %s, %s\n' % ('%12E', '%12E'))
        words.append('                                       C O M P L E X   E I G E N V E C T O R   NO. %s\n' % '%10i')
        #msg.append('                                         C O M P L E X   E I G E N V E C T O R   N O . %10i\n \n' % (imode))

        #msg.append('      POINT ID.   TYPE          T1             T2             T3             R1             R2             R3\n')
        #words += self.get_table_marker()
        return self._write_f06_transient_block(words, header, page_stamp, page_num, f, is_mag_phase, is_sort1)


class RealEigenvectorArray(RealTableArray):
    def __init__(self, data_code, is_sort1, isubcase, dt, f06_flag=False):
        RealTableArray.__init__(self, data_code, is_sort1, isubcase, dt)

        if f06_flag:
            self.itime = -1
            self._times = []
            self.data = []
            self.node_gridtype = []

            self.data_code = data_code
            self.apply_data_code()
            self.set_data_members()

    def build_f06_vectorization(self):
        self.data = array(self.data, dtype='float32')
        self.node_gridtype = array(self.node_gridtype)
        self._times = array(self._times)

    def add_new_f06_mode(self, imode):
        self._times.append(imode)
        self.itime += 1
        self.ntimes += 1
        self.data.append([])
        self.node_gridtype = []

    def update_f06_mode(self, data_code, imode):
        """
        this method is called if the object
        already exits and a new time step is found
        """
        #assert mode >= 0.
        self.data_code = data_code
        self.apply_data_code()
        #self.caseVal = imode
        self.add_new_f06_mode(imode)
        self.set_data_members()

    def read_f06_data(self, data_code, lines):
        """
        it is assumed that all data coming in is correct, so...

           [node_id, grid_type, t1, t2, t3, r1, r2, r3]
           [100, 'G', 1.0, 2.0, 3.0, 4.0, 5.0, 6.0] #  valid

           [101, 'S', 1.0, 2.0] #  invalid
           [101, 'S', 1.0, 0.0, 0.0, 0.0, 0.0, 0.0] #  valid
           [102, 'S', 2.0, 0.0, 0.0, 0.0, 0.0, 0.0] #  valid
        """
        imode = data_code['mode']
        if imode not in self._times:
            self.update_f06_mode(data_code, imode)

        itime = self.itime
        datai = self.data[itime]

        for line in lines:
            if len(line) != 8:
                msg = 'invalid length; even spoints must be in \n'
                msg += '[nid, type, t1, t2, t3, r1, r2, t3] format\nline=%s\n' % line
                msg += 'expected length=8; length=%s' % len(line)
                raise RuntimeError(msg)
            (node_id, grid_type, t1, t2, t3, r1, r2, r3) = line
            self.node_gridtype.append([node_id, self.cast_grid_type(grid_type)])
            datai.append([t1, t2, t3, r1, r2, r3])
        assert self.eigrs[-1] == data_code['eigr'], 'eigrs=%s\ndata_code[eigrs]=%s' %(self.eigrs, data_code['eigr'])

    def write_f06(self, f, header=None, page_stamp='PAGE %s', page_num=1, is_mag_phase=False, is_sort1=True):
        if header is None:
            header = []
        #if self.nonlinear_factor is not None:
            #return self._write_f06_transient(header, page_stamp, page_num, f, is_mag_phase=is_mag_phase, is_sort1=is_sort1)
        # modes get added
        words = '                                         R E A L   E I G E N V E C T O R   N O . %10i\n \n' \
                '      POINT ID.   TYPE          T1             T2             T3             R1             R2             R3\n'

        msg = []
        #if not len(header) >= 3:
            #header.append('')
        for itime in range(self.ntimes):
            node = self.node_gridtype[:, 0]
            gridtype = self.node_gridtype[:, 1]
            t1 = self.data[itime, :, 0]
            t2 = self.data[itime, :, 1]
            t3 = self.data[itime, :, 2]
            r1 = self.data[itime, :, 3]
            r2 = self.data[itime, :, 4]
            r3 = self.data[itime, :, 5]

            dt = self._times[itime]
            #if isinstance(dt, float):
                #header[1] = ' %s = %10.4E\n' % (self.data_code['name'], dt)
            #else:
                #header[1] = ' %s = %10i\n' % (self.data_code['name'], dt)
            msg += header + [words % dt]
            for node_id, gridtypei, t1i, t2i, t3i, r1i, r2i, r3i in zip(node, gridtype, t1, t2, t3, r1, r2, r3):
                sgridtype = self.recast_gridtype_as_string(gridtypei)
                vals = [t1i, t2i, t3i, r1i, r2i, r3i]
                vals2 = write_floats_13e(vals)
                (dx, dy, dz, rx, ry, rz) = vals2
                msg.append('%14i %6s     %-13s  %-13s  %-13s  %-13s  %-13s  %s\n' % (node_id, sgridtype, dx, dy, dz, rx, ry, rz))

            msg.append(page_stamp % page_num)
            f.write(''.join(msg))
            msg = ['']
            page_num += 1
        return page_num - 1


#class Eigenvector(RealTableObject):  # approach_code=2, sort_code=0, thermal=0
    #"""
    #::

      #EIGENVALUE =  6.158494E+07
          #CYCLES =  1.248985E+03         R E A L   E I G E N V E C T O R   N O .          1

      #POINT ID.   TYPE          T1             T2             T3             R1             R2             R3
             #1      G      2.547245E-17  -6.388945E-16   2.292728E+00  -1.076928E-15   2.579163E-17   0.0
          #2002      G     -6.382321E-17  -1.556607E-15   3.242408E+00  -6.530917E-16   1.747180E-17   0.0
          #2003      G      0.0            0.0            0.0            0.0            0.0            0.0
    #"""
    #def __init__(self, data_code, is_sort1, isubcase, imode):
        #RealTableObject.__init__(self, data_code, is_sort1, isubcase, imode)
        ##self.caseVal = mode
        #self.update_dt = self.update_mode
        ##print "mode = %s" %(mode)
        ##print "data_code = ",self.data_code

        ##assert mode>=0.
        #self.gridTypes = {}
        #self.translations = {imode: {}}
        #self.rotations = {imode: {}}

    #def read_f06_data(self, data_code, data):
        #"""
        #it is now assumed that all data coming in is correct, so...


        #so...
           #[node_id, grid_type, t1, t2, t3, r1, r2, r3]
           #[100, 'G', 1.0, 2.0, 3.0, 4.0, 5.0, 6.0] #  valid

           #[101, 'S', 1.0, 2.0] #  invalid
           #[101, 'S', 1.0, 0.0, 0.0, 0.0, 0.0, 0.0] #  valid
           #[102, 'S', 2.0, 0.0, 0.0, 0.0, 0.0, 0.0] #  valid
        #is no longer valid
        #"""
        #imode = data_code['mode']
        #if imode not in self.translations:
            #self.update_mode(data_code, imode)

        #for line in data:
            #if len(line) != 8:
                #msg = 'invalid length; even spoints must be in \n'
                #msg += '[nid, type, t1, t2, t3, r1, r2, t3] format\nline=%s\n' % line
                #msg += 'expected length=8; length=%s' % len(line)
                #raise RuntimeError(msg)
            #(node_id, grid_type, t1, t2, t3, r1, r2, r3) = line
            #self.gridTypes[node_id] = grid_type
            #self.translations[imode][node_id] = array([t1, t2, t3], dtype='float32')
            #self.rotations[imode][node_id] = array([r1, r2, r3], dtype='float32')
        #assert self.eigrs[-1] == data_code['eigr'], 'eigrs=%s\ndata_code[eigrs]=%s' %(self.eigrs, data_code['eigr'])

    #def update_mode(self, data_code, imode):
        #"""
        #this method is called if the object
        #already exits and a new time step is found
        #"""
        ##assert mode>=0.
        #self.data_code = data_code
        #self.apply_data_code()
        ##self.caseVal = imode
        ##print "mode = %s" %(str(mode))
        #self.add_new_mode(imode)
        #self.set_data_members()

    #def add_new_mode(self, imode):
        #self.translations[imode] = {}
        #self.rotations[imode] = {}

    #def eigenvalues(self):
        #return self.eigrs

    #def write_f06(self, f, header=None, page_stamp='PAGE %s', page_num=1, is_mag_phase=False, is_sort1=True):
        #if header is None:
            #header = []
        #"""
        #::

          #EIGENVALUE =  6.158494E+07
              #CYCLES =  1.248985E+03         R E A L   E I G E N V E C T O R   N O .          1

          #POINT ID.   TYPE          T1             T2             T3             R1             R2             R3
                 #1      G      2.547245E-17  -6.388945E-16   2.292728E+00  -1.076928E-15   2.579163E-17   0.0
              #2002      G     -6.382321E-17  -1.556607E-15   3.242408E+00  -6.530917E-16   1.747180E-17   0.0
              #2003      G      0.0            0.0            0.0            0.0            0.0            0.0
        #"""
        #msg = []
        #hasCycle = hasattr(self, 'mode_cycle')

        ##print "self.eigrs =", self.eigrs
        ##print "dir",dir(self)
        #for i, (imode, eig_vals) in enumerate(sorted(iteritems(self.translations))):
            #msg += header
            #freq = self.eigrs[i]
            #msg.append('%16s = %13E\n' % ('EIGENVALUE', freq))

            #if hasCycle:
                #cycle = sqrt(abs(freq))/(2. * pi)
                ##if isinstance(self.mode_cycle, float):
                    ##msg.append('%16s = %13E         R E A L   E I G E N V E C T O R   N O . %10i\n \n' % ('CYCLES', self.mode_cycle, imode))
                ##else:
                    ##msg.append('%16s = %13i         R E A L   E I G E N V E C T O R   N O . %10i\n \n' % ('CYCLES', self.mode_cycle, imode))
                #msg.append('%16s = %13E         R E A L   E I G E N V E C T O R   N O . %10i\n \n' % ('CYCLES', cycle, imode))
                ##msg.append('%16s = %13E         R E A L   E I G E N V E C T O R   N O . %10i\n \n' % ('CYCLES', self.mode_cycle, imode))
            #else:
                #msg.append('                                         R E A L   E I G E N V E C T O R   N O . %10i\n \n' % imode)

            #msg.append('      POINT ID.   TYPE          T1             T2             T3             R1             R2             R3\n')
            #for node_id, displacement in sorted(iteritems(eig_vals)):
                #rotation = self.rotations[imode][node_id]
                #grid_type = self.gridTypes[node_id]
                #(dx, dy, dz) = displacement
                #(rx, ry, rz) = rotation

                #vals = [dx, dy, dz, rx, ry, rz]
                #vals2 = write_floats_13e(vals)
                #[dx, dy, dz, rx, ry, rz] = vals2
                #msg.append('%14i %6s     %-13s  %-13s  %-13s  %-13s  %-13s  %s\n' % (node_id, grid_type, dx, dy, dz, rx, ry, rz))
            #msg.append(page_stamp % page_num)
            #f.write(''.join(msg))
            #msg = ['']
            #page_num += 1
        #return page_num - 1


#class RealEigenvector(ScalarObject):  # approach_code=2, sort_code=0, thermal=0
    #"""
    #::

                                           #R E A L   E I G E N V E C T O R   N O .          1
        #POINT ID.   TYPE          T1             T2             T3             R1             R2             R3
               #1      G      0.0            0.0            0.0            0.0            1.260264E-01   0.0
               #7      G      9.999849E-01   0.0            6.728968E-03   0.0            8.021386E-03   0.0
    #"""
    #def __init__(self, data_code, isubcase, imode):
        #ScalarObject.__init__(self, data_code, isubcase)
        ##self.caseVal = mode
        ##print "mode = %s" %(iMode)
        #self.caseVal = self.getUnsteadyValue()

        ##assert mode>=0.
        #self.gridTypes = {}
        #self.translations = {imode: {}}
        #self.rotations = {imode: {}}

    #def add_new_mode(self, imode):
        #self.translations[imode] = {}
        #self.rotations[imode] = {}

    #def update_dt(self, data_code, dt):
        ##print " self.data_code = ",self.data_code
        #self.data_code = data_code
        #self.apply_data_code()
        #self.set_data_members()
        #self.caseVal = dt

        ##print "*self.data_code = ",self.data_code
        #self.translations[self.caseVal] = {}
        #self.rotations[self.caseVal] = {}

    #def delete_transient(self, dt):
        #del self.translations[dt]
        #del self.rotations[dt]

    #def get_transients(self):
        #k = self.translations.keys()
        #k.sort()
        #return k

    #def add(self, node_id, grid_type, v1, v2, v3, v4, v5, v6):
        #msg = "node_id=%s v1=%s v2=%s v3=%s" % (node_id, v1, v2, v3)
        #msg += "           v4=%s v5=%s v6=%s" % (v4, v5, v6)
        ##print(msg)
        #assert 0 < node_id < 1000000000, msg
        ##assert nodeID not in self.translations

        #if grid_type == 1:
            #Type = 'G'
        #elif grid_type == 2:
            #Type = 'S'
        #elif grid_type == 7:
            #Type = 'L'
        #else:
            #raise ValueError('invalid grid type...grid_type=%s' % grid_type)

        #self.gridTypes[node_id] = Type
        ##print 'self.caseVal = %s' %(self.caseVal),type(self.caseVal)
        ##print "d = ",self.translations
        #self.translations[self.caseVal][node_id] = [v1, v2, v3]
        #self.rotations[self.caseVal][node_id] = [v4, v5, v6]

    #def modes(self):
        #return sorted(self.translations.keys())

    #def eigenvalues(self):
        #return self.eigrs

    #def write_f06(self, f, header=None, page_stamp='PAGE %s', page_num=1, is_mag_phase=False, is_sort1=True):
        #if header is None:
            #header = []
        #"""
        #::

          #EIGENVALUE =  6.158494E+07
                                             #R E A L   E I G E N V E C T O R   N O .          1

          #POINT ID.   TYPE          T1             T2             T3             R1             R2             R3
                 #1      G      2.547245E-17  -6.388945E-16   2.292728E+00  -1.076928E-15   2.579163E-17   0.0
              #2002      G     -6.382321E-17  -1.556607E-15   3.242408E+00  -6.530917E-16   1.747180E-17   0.0
              #2003      G      0.0            0.0            0.0            0.0            0.0            0.0
        #"""
        #msg = []
        ##print self.data_code
        #for i, (imode, eig_vals) in enumerate(sorted(iteritems(self.translations))):
            #msg += header
            #freq = self.eigrs[i]
            #msg.append('%16s = %12E\n' % ('EIGENVALUE', freq))
            #msg.append('                                         R E A L   E I G E N V E C T O R   N O . %10i\n \n' % (imode))
            #msg.append('      POINT ID.   TYPE          T1             T2             T3             R1             R2             R3\n')
            #for node_id, translation in sorted(iteritems(eig_vals)):
                #rotation = self.rotations[imode][node_id]
                #grid_type = self.gridTypes[node_id]
                #(dx, dy, dz) = translation
                #(rx, ry, rz) = rotation

                #vals = [dx, dy, dz, rx, ry, rz]
                #vals2 = write_floats_13e(vals)
                #[dx, dy, dz, rx, ry, rz] = vals2
                #if grid_type == 'S':
                    #msg.append('%14i %6s     %s\n' % (node_id, grid_type, dx))
                #else:
                    #msg.append('%14i %6s     %-13s  %-13s  %-13s  %-13s  %-13s  %s\n' % (node_id, grid_type, dx, dy, dz, rx, ry, rz))
            #msg.append(page_stamp % page_num)
            #f.write(''.join(msg))
            #msg = ['']
            #page_num += 1
        #return page_num - 1


#class ComplexEigenvector(ComplexTableObject):  # approach_code=2, sort_code=0, thermal=0
    #def __init__(self, data_code, is_sort1, isubcase, iMode, apply_data_code=True):
        #ComplexTableObject.__init__(self, data_code, is_sort1, isubcase, iMode, apply_data_code=apply_data_code)
        #self.caseVal = iMode
        #self.update_dt = self.update_mode

        ##print "mode = %s" %(mode)

        ##assert mode>=0.
        ##self.gridTypes = {}
        ##self.translations = {iMode: {}}
        ##self.rotations    = {iMode: {}}

    #def update_mode(self, data_code, imode):
        #"""
        #this method is called if the object
        #already exits and a new time step is found
        #"""
        ##assert mode>=0.
        #self.caseVal = imode
        #self.data_code = data_code
        #self.apply_data_code()
        #self.set_data_members()
        #self.add_new_mode(imode)

    #def add_new_mode(self, imode):
        #self.translations[imode] = {}
        #self.rotations[imode] = {}

    #def eigenvalues(self):
        #return sorted(self.translations.keys())

    #def read_f06_data(self, data_code, data):
        #"""
        #it is now assumed that all data coming in is correct, so...


        #so...
           #[node_id, grid_type, t1, t2, t3, r1, r2, r3]
           #[100, 'G', 1.0, 2.0, 3.0, 4.0, 5.0, 6.0] #  valid

           #[101, 'S', 1.0, 2.0] #  invalid
           #[101, 'S', 1.0, 0.0, 0.0, 0.0, 0.0, 0.0] #  valid
           #[102, 'S', 2.0, 0.0, 0.0, 0.0, 0.0, 0.0] #  valid
        #is no longer valid
        #"""
        #imode = data_code['mode']
        #if imode not in self.translations:
            #self.update_mode(data_code, imode)

        #for line in data:
            #if len(line) != 8:
                #msg = 'invalid length; even spoints must be in \n'
                #msg += '[nid, type, t1, t2, t3, r1, r2, t3] format\nline=%s\n' % line
                #msg += 'expected length=8; length=%s' % len(line)
                #raise RuntimeError(msg)
            #(node_id, grid_type, t1, t2, t3, r1, r2, r3) = line
            #self.gridTypes[node_id] = grid_type
            #self.translations[imode][node_id] = array([t1, t2, t3], dtype='complex64')
            #self.rotations[imode][node_id] = array([r1, r2, r3], dtype='complex64')
        #assert self.eigrs[-1] == data_code['eigr'], 'eigrs=%s\ndata_code[eigrs]=%s' %(self.eigrs, data_code['eigr'])
        #assert self.eigis[-1] == data_code['eigi'], 'eigis=%s\ndata_code[eigis]=%s' %(self.eigis, data_code['eigi'])

    #def write_f06(self, f, header=None, page_stamp='PAGE %s', page_num=1, is_mag_phase=False, is_sort1=True):
        #if header is None:
            #header = []
        #msg = []
        ##print self.data_code
        #hasCycle = hasattr(self, 'mode_cycle')
        #assert len(self.eigis) == len(self.translations)
        #for i, (imode, eig_vals) in enumerate(sorted(iteritems(self.translations))):
            #msg += header
            #eigr = self.eigrs[i]
            #eigi = self.eigis[i]
            #eigr = 0. if eigr == -0. else eigr
            #eigi = 0. if eigi == -0. else eigi

            #f.write(''.join(header))
            #f.write('      %16s = %12E, %12E\n' % ('COMPLEX EIGENVALUE', eigr, eigi))
            #if hasCycle:
                #f.write('%16s = %12E          C O M P L E X   E I G E N V E C T O R   N O . %10i\n \n' % ('CYCLES', self.mode_cycle, imode))
            #else:
                ## TODO: should there be an extra ' \n'???
                #f.write('                                       C O M P L E X   E I G E N V E C T O R   NO. %10i\n' % (imode))

            #if is_mag_phase:
                #f.write('                                                         (MAGNITUDE/PHASE)\n \n')
            #else:
                #f.write('                                                          (REAL/IMAGINARY)\n \n')
            #f.write('      POINT ID.   TYPE          T1             T2             T3             R1             R2             R3\n')
            #for node_id, displacement in sorted(iteritems(eig_vals)):
                #rotation = self.rotations[imode][node_id]
                #grid_type = self.gridTypes[node_id]
                #(dx, dy, dz) = displacement
                #(rx, ry, rz) = rotation

                #vals = [dx, dy, dz, rx, ry, rz]
                #vals2 = write_imag_floats_13e(vals, is_mag_phase)
                #[dxr, dyr, dzr, rxr, ryr, rzr,
                 #dxi, dyi, dzi, rxi, ryi, rzi] = vals2
                #if grid_type == 'S':
                    #f.write('0%13i %6s     %-13s\n' % (node_id, grid_type, dxr))
                    #f.write(' %13s %6s     %-13s\n' % ('', '', dxi))
                #else:
                    #f.write('0%13i %6s     %-13s  %-13s  %-13s  %-13s  %-13s  %s\n' % (node_id, grid_type, dxr, dyr, dzr, rxr, ryr, rzr))
                    #f.write(' %13s %6s     %-13s  %-13s  %-13s  %-13s  %-13s  %s\n' % ('', '', dxi, dyi, dzi, rxi, ryi, rzi))

            #f.write(page_stamp % page_num)
            #msg = ['']
            #page_num += 1
        #return page_num - 1
