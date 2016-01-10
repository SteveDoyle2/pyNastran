from pyNastran.op2.resultObjects.tableObject import RealTableArray, ComplexTableArray


class RealLoadVectorArray(RealTableArray):  # table_code=2, sort_code=0, thermal=0

    def __init__(self, data_code, is_sort1, isubcase, dt):
        RealTableArray.__init__(self, data_code, is_sort1, isubcase, dt)

    def write_f06(self, header, page_stamp, page_num=1, f=None, is_mag_phase=False, is_sort1=True):
        words = ['                                                     L O A D   V E C T O R\n', ]
        #words += self.get_table_marker()
        write_words = True
        if self.nonlinear_factor is not None:
            return self._write_f06_transient_block(
                words, header, page_stamp, page_num, f, write_words,
                is_mag_phase=is_mag_phase, is_sort1=is_sort1)
        return self._write_f06_block(
            words, header, page_stamp, page_num, f, write_words,
            is_mag_phase=False, is_sort1=True
        )


class ComplexLoadVectorArray(ComplexTableArray):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        ComplexTableArray.__init__(self, data_code, is_sort1, isubcase, dt)

    def write_f06(self, header, page_stamp, page_num=1, f=None, is_mag_phase=False, is_sort1=True):
        words = ['                                               C O M P L E X   L O A D   V E C T O R\n', ]
        return self._write_f06_transient_block(
            words, header, page_stamp, page_num, f, is_mag_phase, is_sort1)


class RealTemperatureVectorArray(RealTableArray):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        RealTableArray.__init__(self, data_code, is_sort1, isubcase, dt)

    def write_f06(self, header, page_stamp, page_num=1, f=None, is_mag_phase=False, is_sort1=True):
        words = [
            '                                              T E M P E R A T U R E   V E C T O R\n',
            ' \n',
            '      POINT ID.   TYPE      ID   VALUE     ID+1 VALUE     ID+2 VALUE     ID+3 VALUE     ID+4 VALUE     ID+5 VALUE\n'
        ]
        #words += self.get_table_marker()
        write_words = False
        if self.nonlinear_factor is not None:
            return self._write_f06_transient_block(
                words, header, page_stamp, page_num, f, write_words,
                is_mag_phase=is_mag_phase, is_sort1=is_sort1)
        return self._write_f06_block(words, header, page_stamp, page_num, f, write_words,
                                     is_mag_phase=is_mag_phase, is_sort1=is_sort1)


#class RealThermalVector(RealTableObject):
    #def __init__(self, data_code, is_sort1, isubcase, dt):
        #RealTableObject.__init__(self, data_code, is_sort1, isubcase, dt)

    #def write_f06(self, header, page_stamp, page_num=1, f=None,
                  #is_mag_phase=False, is_sort1=True):
        #if self.nonlinear_factor is not None:
            #return self._write_f06_transient(header, page_stamp, page_num, f,
                                             #is_mag_phase=is_mag_phase, is_sort1=is_sort1)
        #msg = header + ['                                              T E M P E R A T U R E   V E C T O R\n',
                        #' \n',
                        #'      POINT ID.   TYPE      ID   VALUE     ID+1 VALUE     ID+2 VALUE     ID+3 VALUE     ID+4 VALUE     ID+5 VALUE\n']
        #f.write(''.join(msg))
        #for nodeID, translation in sorted(iteritems(self.translations)):
            #rotation = self.rotations[nodeID]
            #grid_type = self.gridTypes[nodeID]

            #(dx, dy, dz) = translation
            #(rx, ry, rz) = rotation
            #vals = [dx, dy, dz, rx, ry, rz]
            #vals2 = write_floats_13e(vals)
            ##if not is_all_zeros:
            #[dx, dy, dz, rx, ry, rz] = vals2
            #f.write('%14i %6s     %-13s  %-13s  %-13s  %-13s  %-13s  %s\n' % (nodeID, grid_type, dx, dy, dz, rx, ry, rz))

        #f.write(page_stamp % page_num)
        #return page_num

    #def _write_f06_transient(self, header, page_stamp, page_num=1, f=None,
                             #is_mag_phase=False, is_sort1=True):
        #words = ['                                              T E M P E R A T U R E   V E C T O R\n',
                 #' \n',
                 #'      POINT ID.   TYPE      ID   VALUE     ID+1 VALUE     ID+2 VALUE     ID+3 VALUE     ID+4 VALUE     ID+5 VALUE\n']

        #for dt, translations in sorted(iteritems(self.translations)):
            #header[1] = ' %s = %10.4E\n' % (self.data_code['name'], dt)
            #f.write(''.join(header + words))
            #for nodeID, translation in sorted(iteritems(translations)):
                #rotation = self.rotations[dt][nodeID]
                #grid_type = self.gridTypes[nodeID]

                #(dx, dy, dz) = translation
                #(rx, ry, rz) = rotation

                #vals = [dx, dy, dz, rx, ry, rz]
                #vals2 = write_floats_13e(vals)
                #if not is_all_zeros:
                    #[dx, dy, dz, rx, ry, rz] = vals2
                    #f.write('%14i %6s     %-13s  %-13s  %-13s  %-13s  %-13s  %s\n' % (nodeID, grid_type, dx, dy, dz, rx, ry, rz))

            #f.write(page_stamp % page_num)
            #page_num += 1
        #return page_num - 1


#class RealThermalLoadVector(RealThermalVector):     # table_code=2, thermal=1
    #def __init__(self, data_code, is_sort1, isubcase, dt):
        #RealThermalVector.__init__(self, data_code, is_sort1, isubcase, dt)


#class RealThermalVelocityVector(RealThermalVector):  # table_code=10, thermal=1
    #def __init__(self, data_code, is_sort1, isubcase, dt):
        #RealThermalVector.__init__(self, data_code, is_sort1, isubcase, dt)

class RealThermalVelocityVectorArray(RealTableArray):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        RealTableArray.__init__(self, data_code, is_sort1, isubcase, dt)

    def write_f06(self, header, page_stamp, page_num=1, f=None, is_mag_phase=False, is_sort1=True):
        words = [
            '                                              THERMAL VELOCITY   V E C T O R\n',
            ' \n',
            '      POINT ID.   TYPE      ID   VALUE     ID+1 VALUE     ID+2 VALUE     ID+3 VALUE     ID+4 VALUE     ID+5 VALUE\n'
        ]
        #words += self.get_table_marker()
        write_words = False
        if self.nonlinear_factor is not None:
            return self._write_f06_transient_block(
                words, header, page_stamp, page_num, f, write_words,
                is_mag_phase=is_mag_phase, is_sort1=is_sort1)
        return self._write_f06_block(
            words, header, page_stamp, page_num, f, write_words,
            is_mag_phase=is_mag_phase, is_sort1=is_sort1)

