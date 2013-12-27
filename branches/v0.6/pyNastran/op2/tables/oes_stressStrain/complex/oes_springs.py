from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)

from ..real.oes_objects import StressObject, StrainObject  # ,array
#from oes_complexObjects import complexStressObject,complexStrainObject

from pyNastran.f06.f06_formatting import writeImagFloats13E # writeFloats13E

complexStressObject = StressObject
complexStrainObject = StrainObject


class ComplexCelasStressObject(complexStressObject):
    """
    ::

                                S T R E S S E S   I N   S C A L A R   S P R I N G S        ( C E L A S 2 )
          TIME         STRESS              TIME         STRESS              TIME         STRESS              TIME         STRESS
      0.0            0.0               5.000000E-02   0.0               1.000000E-01   0.0               1.500000E-01   0.0
      2.000000E-01   0.0               2.500000E-01   0.0               3.000000E-01   0.0               3.500000E-01   0.0
    """
    def __init__(self, data_code, is_sort1, isubcase, dt=None):
        StressObject.__init__(self, data_code, isubcase)
        self.eType = {}
        self.element_name = self.data_code['element_name']

        self.code = [self.format_code, self.sort_code, self.s_code]
        self.stress = {}

        self.dt = dt
        if is_sort1:
            if dt is not None:
                #self.add = self.add_sort1
                self.add_new_eid = self.add_new_eid_sort1
        else:
            assert dt is not None
            #self.add = self.addSort2
            self.add_new_eid = self.add_new_eid_sort2

    def get_stats(self):
        nelements = len(self.eType)

        msg = self.get_data_code()
        if self.dt is not None:  # transient
            ntimes = len(self.stress)
            msg.append('  type=%s ntimes=%s nelements=%s\n'
                       % (self.__class__.__name__, ntimes, nelements))
        else:
            msg.append('  type=%s nelements=%s\n' % (self.__class__.__name__,
                                                     nelements))
        msg.append('  eType, stress\n')
        return msg

    def delete_transient(self, dt):
        del self.stress[dt]

    def get_transients(self):
        k = self.stress.keys()
        k.sort()
        return k

    def add_new_transient(self, dt):
        """initializes the transient variables"""
        self.element_name = self.data_code['element_name']
        self.dt = dt
        self.stress[dt] = {}

    def add_new_eid(self, dt, eid, stress):
        self.eType[eid] = self.element_name
        self.stress[eid] = stress

    def add_new_eid_sort1(self, dt, eid, stress):
        if dt not in self.stress:
            self.add_new_transient(dt)
        self.eType[eid] = self.element_name
        self.stress[dt][eid] = stress

    def add_new_eid_sort2(self, eid, dt, stress):
        if dt not in self.stress:
            self.add_new_transient(dt)
        self.eType[eid] = self.element_name
        self.stress[dt][eid] = stress

    def write_f06(self, header, pageStamp, pageNum=1, f=None, is_mag_phase=False):
        """
        .. todo:: doesnt write...
        """
        if self.nonlinear_factor is not None:
            return self._write_f06_transient(header, pageStamp, pageNum, f, is_mag_phase)
        return 'ComplexCelasStressObject write_f06 not implemented...\n'
        #raise NotImplementedError()

    def _write_f06_transient(self, header, pageStamp, pageNum=1, f=None, is_mag_phase=False):
        """
        .. todo:: improve formatting
        """
        words = ['                         C O M P L E X   F O R C E S   I N   S C A L A R   S P R I N G S   ( C E L A S 1 )\n',
                 '                                                          (REAL/IMAGINARY)\n',
                 ' \n',
                 '                ELEMENT                                                   ELEMENT\n',
                 '                  ID.                    FORCE                              ID.                    FORCE\n']
#                   1001       1.537879E+01 /  0.0                            1002       1.537879E+01 /  0.0
#                   1003       1.537879E+01 /  0.0                            1004       1.537879E+01 /  0.0
#                   1005       1.537879E+01 /  0.0                            1006       1.537879E+01 /  0.0
#                   1007       7.689395E+00 /  0.0                            1008       7.689395E+00 /  0.0
#                   1009       7.689395E+00 /  0.0                            1010       7.689395E+00 /  0.0
        msg = []
        is_mag_phase = False
        for dt, Stress in sorted(self.stress.iteritems()):
            if isinstance(dt, float):  # fix
                header[1] = ' %s = %10.4E float %s\n' % (self.data_code[
                    'name'], dt, self.analysis_code)
            else:
                header[1] = ' %s = %10i integer %s\n' % (self.data_code[
                    'name'], dt, self.analysis_code)
            msg += header + words

            i = 0
            for elementID, stress in sorted(Stress.iteritems()):

                if is_mag_phase:
                    stressr = abs(stressr)
                    stressi = angle(stress, deg=True)
                else:
                    stressr = stress.real
                    stressi = stress.imag

                (vals2, isAllZeros) = writeImagFloats13E([stress], is_mag_phase)
                if i == 0:
                    elementID1 = elementID
                    [stress1Real, stress1Imag] = vals2
                if i == 1:
                    elementID2 = elementID
                    [stress2Real, stress2Imag] = vals2
                    msg.append('%14i %13s / %13s  %14i %13s / %-s\n' % (elementID1, stress1Real, stress1Imag, elementID2, stress2Real, stress2Imag.rstrip()))
                    i = -1
                i += 1
            msg.append(pageStamp + str(pageNum) + '\n')
            f.write(''.join(msg))
            msg = ['']
            pageNum += 1
        return pageNum - 1

    def __reprTransient__(self):
        msg = '---CELASx STRESSES---\n'
        msg += '%-6s %6s ' % ('EID', 'eType')
        headers = ['stress']
        for header in headers:
            msg += '%10s ' % header
        msg += '\n'

        for dt, stress in sorted(self.stress.iteritems()):
            msg += '%s = %g\n' % (self.data_code['name'], dt)
            for eid, istress in sorted(stress.iteritems()):
                msg += '%-6g %6s ' % (eid, self.eType[eid])
                if abs(istress) < 1e-6:
                    msg += '%10s ' % '0'
                else:
                    msg += '%10g ' % istress
                msg += '\n'
        return msg


class ComplexCelasStrainObject(complexStrainObject):
    def __init__(self, data_code, is_sort1, isubcase, dt=None):
        StrainObject.__init__(self, data_code, isubcase)
        self.eType = {}
        self.element_name = self.data_code['element_name']

        self.code = [self.format_code, self.sort_code, self.s_code]

        self.strain = {}

        self.dt = dt
        if is_sort1:
            if dt is not None:
                #self.add = self.add_sort1
                self.add_new_eid = self.add_new_eid_sort1
        else:
            assert dt is not None
            #self.add = self.addSort2
            self.add_new_eid = self.add_new_eid_sort2

    def delete_transient(self, dt):
        del self.strain[dt]

    def get_transients(self):
        k = self.strain.keys()
        k.sort()
        return k

    def add_new_transient(self, dt):
        """
        initializes the transient variables
        """
        self.strain[dt] = {}

    def add_new_eid(self, dt, eid, strain):
        assert eid >= 0
        #self.eType = self.eType
        self.eType[eid] = self.element_name
        self.strain[eid] = strain

    def add_new_eid_sort1(self, dt, eid, strain):
        assert eid >= 0

        self.eType[eid] = self.element_type
        self.strain[dt][eid] = strain