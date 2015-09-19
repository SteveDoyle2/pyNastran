from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
from six import iteritems
from pyNastran.op2.tables.oes_stressStrain.real.oes_objects import StressObject, StrainObject
from pyNastran.f06.f06_formatting import writeFloats13E


class RealBushStress(StressObject):

    def __init__(self, data_code, is_sort1, isubcase, dt):
        StressObject.__init__(self, data_code, isubcase)
        self.eType = {}

        self.code = [self.format_code, self.sort_code, self.s_code]
        self.translations = {}
        self.rotations = {}

        self.dt = dt
        if is_sort1:
            if dt is not None:
                self.add_new_eid = self.add_new_eid_sort1
        else:
            assert dt is not None
            self.add_new_eid = self.add_new_eid_sort2

    def get_stats(self):
        nelements = len(self.eType)

        msg = self.get_data_code()
        if self.dt is not None:  # transient
            ntimes = len(self.translations)
            msg.append('  type=%s ntimes=%s nelements=%s\n'
                       % (self.__class__.__name__, ntimes, nelements))
        else:
            msg.append('  imaginary type=%s nelements=%s\n' % (self.__class__.__name__,
                                                               nelements))
        msg.append('  eType, translations, rotations\n')
        return msg

    def add_f06_data(self, data, transient):
        if transient is None:
            for line in data:
                (eType, eid, tx, ty, tz, rx, ry, rz) = line
                self.eType[eid] = 'CBUSH'
                self.translations[eid] = [tx, ty, tz]
                self.rotations[eid] = [rx, ry, rz]
            return

        (dtName, dt) = transient
        self.data_code['name'] = dtName

        if dt not in self.translations:
            self.update_dt(self.data_code, dt)

        for line in data:
            (eType, eid, tx, ty, tz, rx, ry, rz) = line
            self.eType[eid] = 'CBUSH'
            self.translations[dt][eid] = [tx, ty, tz]
            self.rotations[dt][eid] = [rx, ry, rz]

    def delete_transient(self, dt):
        del self.translations[dt]
        del self.rotations[dt]

    def get_transients(self):
        k = self.translations.keys()
        k.sort()
        return k

    def add_new_transient(self, dt):
        """
        initializes the transient variables
        """
        self.dt = dt
        self.translations[dt] = {}
        self.rotations[dt] = {}

    def add_new_eid(self, eType, dt, eid, tx, ty, tz, rx, ry, rz):
        self.eType[eid] = eType
        self.translations[eid] = [tx, ty, tz]
        self.rotations[eid] = [rx, ry, rz]

    def add_new_eid_sort1(self, eType, dt, eid, tx, ty, tz, rx, ry, rz):
        if dt not in self.translations:
            self.add_new_transient(dt)
        self.eType[eid] = eType
        self.translations[dt][eid] = [tx, ty, tz]
        self.rotations[dt][eid] = [rx, ry, rz]

    def write_f06(self, header, page_stamp, page_num=1, f=None, is_mag_phase=False, is_sort1=True):
        if self.nonlinear_factor is not None:
            return self._write_f06_transient(header, page_stamp, page_num, f, is_mag_phase=is_mag_phase, is_sort1=is_sort1)

        msg = header + [
            '                                  S T R E S S E S   I N   B U S H   E L E M E N T S        ( C B U S H )\n \n',
            '                  ELEMENT-ID        STRESS-TX     STRESS-TY     STRESS-TZ    STRESS-RX     STRESS-RY     STRESS-RZ \n',
        ]

        for eid, (tx, ty, tz) in sorted(iteritems(self.translations)):
            eType = self.eType[eid]
            (rx, ry, rz) = self.rotations[eid]

            vals = [tx, ty, tz, rx, ry, rz]
            (vals2, is_all_zeros) = writeFloats13E(vals)
            [tx, ty, tz, rx, ry, rz] = vals2
            msg.append('0                   %8i     %-13s %-13s %-13s %-13s %-13s %s\n' % (eid, tx, ty, tz, rx, ry, rz))

        msg.append(page_stamp % page_num)
        f.write(''.join(msg))
        return page_num

    def _write_f06_transient(self, header, page_stamp, page_num=1, f=None, is_mag_phase=False, is_sort1=True):
        words = [
            '                                  S T R E S S E S   I N   B U S H   E L E M E N T S        ( C B U S H )\n \n',
            '                  ELEMENT-ID        STRESS-TX     STRESS-TY     STRESS-TZ    STRESS-RX     STRESS-RY     STRESS-RZ \n',
        ]
        msg = []
        for dt, translations in sorted(iteritems(self.translations)):
            header[1] = ' %s = %10.4E\n' % (self.data_code['name'], dt)
            msg += header + words
            for eid, (tx, ty, tz) in sorted(iteritems(translations)):
                eType = self.eType[eid]
                (rx, ry, rz) = self.rotations[dt][eid]

                vals = [tx, ty, tz, rx, ry, rz]
                (vals2, is_all_zeros) = writeFloats13E(vals)
                [tx, ty, tz, rx, ry, rz] = vals2
                msg.append('0%8i   %-13s  %-13s  %-13s  %-13s  %-13s  %s\n' % (eid, tx, ty, tz, rx, ry, rz))

            msg.append(page_stamp % page_num)
            f.write(''.join(msg))
            page_num += 1
        return page_num - 1


class RealBushStrain(StrainObject):
    """
    """
    def __init__(self, data_code, is_sort1, isubcase, dt):
        StrainObject.__init__(self, data_code, isubcase)
        self.eType = {}

        self.code = [self.format_code, self.sort_code, self.s_code]
        self.translations = {}
        self.rotations = {}

        if is_sort1:
            if dt is not None:
                self.add_new_eid = self.add_new_eid_sort1
        else:
            assert dt is not None
            #self.add_new_eid = self.add_new_eid_sort2

    def get_stats(self):
        nelements = len(self.eType)

        msg = self.get_data_code()
        if self.dt is not None:  # transient
            ntimes = len(self.translations)
            msg.append('  type=%s ntimes=%s nelements=%s\n' % (self.__class__.__name__, ntimes, nelements))
        else:
            msg.append('  imaginary type=%s nelements=%s\n' % (self.__class__.__name__, nelements))
        msg.append('  eType, translations, rotations\n')
        return msg

    def add_f06_data(self, data, transient):
        if transient is None:
            for line in data:
                (eType, eid, tx, ty, tz, rx, ry, rz) = line
                self.eType[eid] = 'CBUSH'
                self.translations[eid] = [tx, ty, tz]
                self.rotations[eid] = [rx, ry, rz]
            return

        (dtName, dt) = transient
        self.data_code['name'] = dtName

        if dt not in self.translations:
            self.update_dt(self.data_code, dt)

        for line in data:
            (eType, eid, tx, ty, tz, rx, ry, rz) = line
            self.eType[eid] = 'CBUSH'
            self.translations[dt][eid] = [tx, ty, tz]
            self.rotations[dt][eid] = [rx, ry, rz]

    def delete_transient(self, dt):
        del self.translations[dt]
        del self.rotations[dt]

    def get_transients(self):
        k = self.translations.keys()
        k.sort()
        return k

    def add_new_transient(self, dt):
        """
        initializes the transient variables
        """
        self.dt = dt
        self.translations[dt] = {}
        self.rotations[dt] = {}

    def add_new_eid(self, eType, dt, eid, tx, ty, tz, rx, ry, rz):
        self.eType[eid] = eType
        self.translations[eid] = [tx, ty, tz]
        self.rotations[eid] = [rx, ry, rz]

    def add_new_eid_sort1(self, eType, dt, eid, tx, ty, tz, rx, ry, rz):
        if dt not in self.translations:
            self.add_new_transient(dt)
        self.eType[eid] = eType
        self.translations[dt][eid] = [tx, ty, tz]
        self.rotations[dt][eid] = [rx, ry, rz]

    def write_f06(self, header, page_stamp, page_num=1, f=None, is_mag_phase=False, is_sort1=True):
        if self.nonlinear_factor is not None:
            return self._write_f06_transient(header, page_stamp, page_num, f, is_mag_phase=is_mag_phase, is_sort1=is_sort1)

        msg = header + [
            '                                    S T R A I N S   I N   B U S H   E L E M E N T S        ( C B U S H )\n'
            ' \n'
            '                  ELEMENT-ID        STRAIN-TX     STRAIN-TY     STRAIN-TZ    STRAIN-RX     STRAIN-RY     STRAIN-RZ \n'
        ]
        for eid, (tx, ty, tz) in sorted(iteritems(self.translations)):
            eType = self.eType[eid]
            (rx, ry, rz) = self.rotations[eid]
            vals = [tx, ty, tz, rx, ry, rz]
            (vals2, is_all_zeros) = writeFloats13E(vals)
            [tx, ty, tz, rx, ry, rz] = vals2
            msg.append('0                   %8i     %-13s %-13s %-13s %-13s %-13s %s\n' % (eid, tx, ty, tz, rx, ry, rz))

        msg.append(page_stamp % page_num)
        f.write(''.join(msg))
        return page_num

    def __write_f06_transient(self, header, page_stamp, page_num=1, f=None, is_mag_phase=False, is_sort1=True):
        f.write('%s _write_f06_transient not implemented...\n' % self.__class__.__name__)
        raise NotImplementedError('CBUSH')
        return page_num - 1
