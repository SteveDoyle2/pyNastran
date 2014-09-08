from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
#import sys

from ..real.oes_objects import StressObject, StrainObject


class ComplexRodDamper(StressObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        StressObject.__init__(self, data_code, isubcase)
        self.eType = 'CBUSH'

        self.code = [self.format_code, self.sort_code, self.s_code]
        self.axial = {}
        self.torsion = {}


class ComplexRodStress(StressObject):
    """
    """
    def __init__(self, data_code, is_sort1, isubcase, dt):
        StressObject.__init__(self, data_code, isubcase)
        self.eType = 'CROD'

        self.code = [self.format_code, self.sort_code, self.s_code]
        self.axial = {}
        self.torsion = {}

        self.dt = dt
        if is_sort1:
            if dt is not None:
                #self.add = self.add_sort1
                self.add_new_eid = self.add_new_eid_sort1
        else:
            assert dt is not None
            #self.add = self.addSort2
            self.add_new_eid = self.add_new_eid_sort2

    def getLength(self):
        return (20, '4f')

    def add_f06_data(self, data, transient):
        if transient is None:
            for line in data:
                (eid, axial, torsion) = line
                self.axial[eid] = axial
                self.torsion[eid] = torsion
            return

        (dtName, dt) = transient
        self.dt = dt
        self.data_code['name'] = dtName
        if dt not in self.s1:
            self.update_dt(self.data_code, dt)

        for line in data:
            (eid, axial, torsion) = line
            self.axial[dt][eid] = axial
            self.torsion[dt][eid] = torsion

    def delete_transient(self, dt):
        del self.axial[dt]
        del self.torsion[dt]

    def get_transients(self):
        k = self.axial.keys()
        k.sort()
        return k

    def add_new_transient(self, dt):
        """
        initializes the transient variables
        """
        self.dt = dt
        self.axial[dt] = {}
        self.torsion[dt] = {}

    def add_new_eid(self, dt, eid, axial, torsion):
        #print "Rod Stress add..."
        #(axial, torsion) = out
        assert isinstance(eid, int)
        self.axial[eid] = axial
        self.torsion[eid] = torsion

    def add_new_eid_sort1(self, dt, eid, axial, torsion):
        #(axial, torsion) = out

        if dt not in self.axial:
            self.add_new_transient(dt)
        self.axial[dt][eid] = axial
        self.torsion[dt][eid] = torsion

    def add_new_eid_sort2(self, eid, dt, axial, torsion):
        #(axial, torsion) = out

        if dt not in self.axial:
            self.add_new_transient(dt)
        self.axial[dt][eid] = axial
        self.torsion[dt][eid] = torsion

    def _write_f06(self, header, pageStamp, page_num=1, f=None, is_mag_phase=False):
        f.write('ComplexRodStress write_f06 not implemented...\n')
        return page_num
        raise NotImplementedError()
        if self.nonlinear_factor is not None:
            return self._write_f06_transient(header, pageStamp, page_num, f, is_mag_phase)

        msg = header + ['                                     S T R E S S E S   I N   R O D   E L E M E N T S      ( C R O D )\n',
                        '       ELEMENT       AXIAL       SAFETY      TORSIONAL     SAFETY       ELEMENT       AXIAL       SAFETY      TORSIONAL     SAFETY\n',
                        '         ID.        STRESS       MARGIN        STRESS      MARGIN         ID.        STRESS       MARGIN        STRESS      MARGIN\n']
        out = []
        for eid in sorted(self.axial):
            axial = self.axial[eid]
            torsion = self.torsion[eid]
            (vals2, is_all_zeros) = writeFloatsImag13E([axial, torsion])
            (axial, torsion) = vals2
            out.append([eid, axial, torsion])

        nOut = len(out)
        nWrite = nOut
        if nOut % 2 == 1:
            nWrite = nOut - 1
        for i in xrange(0, nWrite, 2):
            #print i,out[i:]
            outLine = '      %8i   %13s  %10.4E %13s  %10.4E   %8i   %13s  %10.4E %13s  %10.4E\n' % (tuple(out[i] + out[i + 1]))
            msg.append(outLine)

        if nOut % 2 == 1:
            outLine = '      %8i   %13s  %10.4E %13s  %10.4E\n' % (
                tuple(out[-1]))
            msg.append(outLine)
        msg.append(pageStamp % page_num)
        f.write(''.join(msg))
        return page_num

    def _write_f06_transient(self, header, pageStamp, page_num=1, f=None, is_mag_phase=False):
        return 'ComplexRodStress _write_f06_transient not implemented...', page_num
        raise NotImplementedError()
        words = ['                                     S T R E S S E S   I N   R O D   E L E M E N T S      ( C R O D )\n',
                 '       ELEMENT       AXIAL       SAFETY      TORSIONAL     SAFETY       ELEMENT       AXIAL       SAFETY      TORSIONAL     SAFETY\n',
                 '         ID.        STRESS       MARGIN        STRESS      MARGIN         ID.        STRESS       MARGIN        STRESS      MARGIN\n']
        msg = []
        for dt, axials in sorted(self.axial.iteritems()):
            dtLine = '%14s = %12.5E\n' % (self.data_code['name'], dt)
            header[2] = dtLine
            msg += header + words
            out = []
            for eid in sorted(axials):
                axial = self.axial[dt][eid]
                torsion = self.torsion[dt][eid]

                (vals2, is_all_zeros) = writeFloatsImag13E([axial, torsion])
                (axial, torsion) = vals2
                out.append([eid, axial, MSa, torsion, MSt])

            nOut = len(out)
            nWrite = nOut
            if nOut % 2 == 1:
                nWrite = nOut - 1
            for i in xrange(0, nWrite, 2):
                outLine = '      %8i   %13s  %10.4E %13s  %10.4E   %8i   %13s  %10.4E %13s  %10.4E\n' % (tuple(out[i] + out[i + 1]))
                msg.append(outLine)

            if nOut % 2 == 1:
                outLine = '      %8i   %13s  %10.4E %13s  %10.4E\n' % (tuple(out[-1]))
                msg.append(outLine)
            msg.append(pageStamp % page_num)
            f.write(''.join(msg))
            page_num += 1
        return page_num - 1


class ComplexRodStrain(StrainObject):
    """
    """
    def __init__(self, data_code, is_sort1, isubcase, dt):
        StrainObject.__init__(self, data_code, isubcase)
        self.eType = 'CROD'  # {} # 'CROD/CONROD/CTUBE'

        self.code = [self.format_code, self.sort_code, self.s_code]

        self.axial = {}
        self.torsion = {}

        self.dt = dt
        if is_sort1:
            if dt is not None:
                #self.add = self.add_sort1
                self.add_new_eid = self.add_new_eid_sort1
        else:
            assert dt is not None
            #self.add = self.addSort2
            self.add_new_eid = self.add_new_eid_sort2

    def add_f06_data(self, data, transient):
        if transient is None:
            for line in data:
                (eid, axial, torsion) = line
                self.axial[eid] = axial
                self.torsion[eid] = torsion
            return

        (dtName, dt) = transient
        self.dt = dt
        self.data_code['name'] = dtName
        if dt not in self.s1:
            self.update_dt(self.data_code, dt)

        for line in data:
            (eid, axial, torsion) = line
            self.axial[dt][eid] = axial
            self.torsion[dt][eid] = torsion

    def delete_transient(self, dt):
        del self.axial[dt]
        del self.torsion[dt]

    def get_transients(self):
        k = self.axial.keys()
        k.sort()
        return k

    def add_new_transient(self, dt):
        """
        initializes the transient variables
        """
        self.dt = dt
        self.axial[self.dt] = {}
        self.torsion[self.dt] = {}

    def add_new_eid(self, dt, eid, axial, torsion):
        #(axial, torsion) = out
        assert eid >= 0
        #self.eType = self.eType
        self.axial[eid] = axial
        self.torsion[eid] = torsion

    def add_new_eid_sort1(self, dt, eid, axial, torsion):
        #(axial, torsion) = out
        assert eid >= 0
        #self.eType[eid] = self.element_type
        if dt not in self.axial:
            self.add_new_transient(dt)
        self.axial[dt][eid] = axial
        self.torsion[dt][eid] = torsion

    def add_new_eid_sort2(self, eid, dt, out):
        (axial, torsion) = out
        assert eid >= 0
        #self.eType[eid] = self.element_type
        if dt not in self.axial:
            self.add_new_transient(dt)
        self.axial[dt][eid] = axial
        self.torsion[dt][eid] = torsion

    def _write_f06(self, header, pageStamp, page_num=1, f=None, is_mag_phase=False):
        return 'ComplexRodStrain write_f06 not implemented...', page_num
        raise NotImplementedError()
        if self.dt is not None:
            return self._write_f06_transient(header, pageStamp, page_num, f, is_mag_phase)

        msg = header + ['                                       S T R A I N S   I N   R O D   E L E M E N T S      ( C R O D )\n',
                        '       ELEMENT       AXIAL       SAFETY      TORSIONAL     SAFETY       ELEMENT       AXIAL       SAFETY      TORSIONAL     SAFETY\n',
                        '         ID.        STRAIN       MARGIN        STRAIN      MARGIN         ID.        STRAIN       MARGIN        STRAIN      MARGIN\n']
        out = []
        for eid in sorted(self.axial):
            axial = self.axial[eid]
            torsion = self.torsion[eid]
            (vals2, is_all_zeros) = writeFloatsImag13E([axial, torsion])
            (axial, torsion) = vals2
            out.append([eid, axial, torsion])

        nOut = len(out)
        nWrite = nOut
        if nOut % 2 == 1:
            nWrite = nOut - 1
        for i in xrange(0, nWrite, 2):
            outLine = '      %8i   %13s  %10.4E %13s  %10.4E   %8i   %13s  %10.4E %13s  %10.4E\n' % (tuple(out[i] + out[i + 1]))
            msg.append(outLine)

        if nOut % 2 == 1:
            outLine = '      %8i   %13s  %10.4E %13s  %10.4E\n' % (
                tuple(out[-1]))
            msg.append(outLine)
        msg.append(pageStamp % page_num)
        f.write(''.join(msg))
        return page_num

    def _write_f06_transient(self, header, pageStamp, page_num=1, f=None, is_mag_phase=False):
        return 'ComplexRodStress _write_f06_transient not implemented...', page_num
        raise NotImplementedError()
        words = ['                                       S T R A I N S   I N   R O D   E L E M E N T S      ( C R O D )\n',
                 '       ELEMENT       AXIAL       SAFETY      TORSIONAL     SAFETY       ELEMENT       AXIAL       SAFETY      TORSIONAL     SAFETY\n',
                 '         ID.        STRAIN       MARGIN        STRAIN      MARGIN         ID.        STRAIN       MARGIN        STRAIN      MARGIN\n']
        msg = []
        for dt, axials in sorted(self.axial.iteritems()):
            dtLine = '%14s = %12.5E\n' % (self.data_code['name'], dt)
            header[2] = dtLine
            msg += header + words
            out = []
            for eid in sorted(axials):
                axial = self.axial[dt][eid]
                torsion = self.torsion[dt][eid]

                out.append([eid, axial, torsion])

            nOut = len(out)
            nWrite = nOut
            if nOut % 2 == 1:
                nWrite = nOut - 1
            for i in xrange(0, nWrite, 2):
                outLine = '      %8i   %13.6E  %10.4E %13.6E  %10.4E   %8i   %13.6E  %10.4E %13.6E  %10.4E\n' % (tuple(out[i] + out[i + 1]))
                msg.append(outLine)

            if nOut % 2 == 1:
                outLine = '      %8i   %13.6E  %10.4E %13.6E  %10.4E\n' % (
                    tuple(out[-1]))
                msg.append(outLine)
            msg.append(pageStamp % page_num)
            page_num += 1
            f.write(''.join(msg))
        return page_num - 1