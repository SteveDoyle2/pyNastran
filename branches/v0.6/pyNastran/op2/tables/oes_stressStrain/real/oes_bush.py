from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)

from .oes_objects import StressObject, StrainObject
from pyNastran.f06.f06_formatting import writeFloats13E


class BushStressObject(StressObject):

    def __init__(self, data_code, is_sort1, isubcase, dt=None):
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

    def write_f06(self, header, pageStamp, pageNum=1, f=None, is_mag_phase=False):
        if self.nonlinear_factor is not None:
            return self._write_f06_transient(header, pageStamp, pageNum, f, is_mag_phase)

        msg = header + [
            '                                  S T R E S S E S   I N   B U S H   E L E M E N T S        ( C B U S H )\n\n',
            '                  ELEMENT-ID        STRESS-TX     STRESS-TY     STRESS-TZ    STRESS-RX     STRESS-RY     STRESS-RZ \n',
        ]

        for eid, (tx, ty, tz) in sorted(self.translations.iteritems()):
            eType = self.eType[eid]
            (rx, ry, rz) = self.rotations[eid]

            vals = [tx, ty, tz, rx, ry, rz]
            (vals2, isAllZeros) = writeFloats13E(vals)
            [tx, ty, tz, rx, ry, rz] = vals2
            msg.append('0%8i   %13s  %13s  %13s  %13s  %13s  %-s\n' % (eid, tx, ty, tz, rx, ry, rz.rstrip()))

        msg.append(pageStamp % pageNum)
        f.write(''.join(msg))
        return pageNum

    def _write_f06_transient(self, header, pageStamp, pageNum=1, f=None, is_mag_phase=False):
        words = [
            '                                  S T R E S S E S   I N   B U S H   E L E M E N T S        ( C B U S H )\n\n',
            '                  ELEMENT-ID        STRESS-TX     STRESS-TY     STRESS-TZ    STRESS-RX     STRESS-RY     STRESS-RZ \n',
        ]
        msg = []
        for dt, translations in sorted(self.translations.iteritems()):
            header[1] = ' %s = %10.4E\n' % (self.data_code['name'], dt)
            msg += header + words
            for eid, (tx, ty, tz) in sorted(translations.iteritems()):
                eType = self.eType[eid]
                (rx, ry, rz) = self.rotations[dt][eid]

                vals = [tx, ty, tz, rx, ry, rz]
                (vals2, isAllZeros) = writeFloats13E(vals)
                [tx, ty, tz, rx, ry, rz] = vals2
                msg.append('0%8i   %13s  %13s  %13s  %13s  %13s  %-s\n' % (eid, tx, ty, tz, rx, ry, rz.rstrip()))

            msg.append(pageStamp % pageNum)
            f.write(''.join(msg))
            pageNum += 1
        return pageNum - 1


class BushStrainObject(StrainObject):
    """
    """
    def __init__(self, data_code, is_sort1, isubcase, dt=None):
        raise NotImplementedError('is this used?')
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
            self.add_new_eid = self.add_new_eid_sort2

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


    def write_f06(self, header, pageStamp, pageNum=1, f=None, is_mag_phase=False):
        raise NotImplementedError('CBUSH')
        return 'BushStress write_f06 not implemented...', pageNum
        if self.nonlinear_factor is not None:
            return self._write_f06_transient(header, pageStamp, pageNum, f, is_mag_phase)

        msg = header + [
            '                                  S T R A I N S    I N   B A R   E L E M E N T S          ( C B A R )\n',
            '  ELEMENT        SA1            SA2            SA3            SA4           AXIAL          SA-MAX         SA-MIN     M.S.-T\n',
            '    ID.          SB1            SB2            SB3            SB4           STRAIN         SB-MAX         SB-MIN     M.S.-C\n',
        ]
        for eid, E1s in sorted(self.e1.iteritems()):
            eType = self.eType[eid]
            axial = self.axial[eid]

            e1 = self.e1[eid]
            e2 = self.e2[eid]
            e3 = self.e3[eid]
            e4 = self.e4[eid]
            vals = [e1[0], e2[0], e3[0], e4[0], axial,
                    e1[1], e2[1], e3[1], e4[1]]
            (vals2, isAllZeros) = writeFloats13E(vals)
            [e10, e20, e30, e40, axial,
             e11, e21, e31, e41] = vals2

            msg.append('0%8i   %13s  %13s  %13s  %13s  %13s  %13s  %13s %-s\n' % (eid, e10, e20, e30, e40, axial.rstrip()))
            msg.append(' %8s   %13s  %13s  %13s  %13s  %13s  %13s  %13s %-s\n' % ('', e11, e21, e31, e41.rstrip()))
        msg.append(pageStamp % pageNum)
        f.write(''.join(msg))
        return pageNum

    def _write_f06_transient(self, header, pageStamp, pageNum=1, f=None, is_mag_phase=False):
        raise NotImplementedError('CBUSH')
        words = [
            '                                  S T R A I N S    I N   B A R   E L E M E N T S           ( C B A R )\n',
            '  ELEMENT        SA1            SA2            SA3            SA4           AXIAL          SA-MAX         SA-MIN     M.S.-T\n',
            '    ID.          SB1            SB2            SB3            SB4           STRAIN         SB-MAX         SB-MIN     M.S.-C\n',
        ]
        msg = []
        for dt, E1s in sorted(self.e1.iteritems()):
            header[1] = ' %s = %10.4E\n' % (self.data_code['name'], dt)
            msg += header + words
            for eid, e1s in sorted(E1s.iteritems()):
                eType = self.eType[eid]
                axial = self.axial[eid]

                e1 = self.e1[eid]
                e2 = self.e2[eid]
                e3 = self.e3[eid]
                e4 = self.e4[eid]
                vals = [e1[0], e2[0], e3[0], e4[0], axial,
                        e1[1], e2[1], e3[1], e4[1]]
                (vals2, isAllZeros) = writeFloats13E(vals)
                [e10, e20, e30, e40,
                 e11, e21, e31, e41] = vals2

                msg.append('0%8i   %13s  %13s  %13s  %13s  %13s  %13s  %13s %-s\n' % (eid, e10, e20, e30, e40, axial.rstrip()))
                msg.append(' %8s   %13s  %13s  %13s  %13s  %13s  %13s  %13s %-s\n' % ('', e11, e21, e31, e41.rstrip()))
            msg.append(pageStamp % pageNum)
            f.write(''.join(msg))
            msg = ['']
            pageNum += 1
        return pageNum - 1