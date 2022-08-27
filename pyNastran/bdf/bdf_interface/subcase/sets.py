from .utils import write_set
from .subcase_base import CaseControlCard


class SET(CaseControlCard):
    type = 'SET'
    def __init__(self, set_id, values):
        super(SET, self).__init__()
        self.set_id = int(set_id)

        #values2 = expand_thru_case_control(values)
        self.value = values

    @property
    def key(self):
        """temporary method to emulate the old key attribute"""
        return '%s %s' % (self.type, self.set_id)

    def __iter__(self):
        """temporary method to emulate the old list access style"""
        value = self
        options = None
        param_type = 'OBJ-type'
        return iter([value, options, param_type])

    @classmethod
    def add_from_case_control(cls, line_upper, lines, i):
        """add method used by the CaseControl class"""
        line = lines[i]
        sline = line_upper.split('=')
        assert len(sline) == 2, sline

        key, value = sline
        try:
            (key, set_id) = key.split()
        except Exception:
            raise RuntimeError(key)

        assert key.upper() == key, key
        unused_options = int(set_id)

        #if self.debug:
            #self.log.debug('SET-type key=%r set_id=%r' % (key, set_id))
        fivalues = value.rstrip(' ,').split(',')  # float/int values

        #: .. todo:: should be more efficient multiline reader...
        # read more lines....
        if line[-1].strip() == ',':
            i += 1
            #print("rawSETLine = %r" % (lines[i]))
            while 1:
                if lines[i].strip()[-1] == ',':
                    fivalues += lines[i][:-1].split(',')
                else:  # last case
                    fivalues += lines[i].split(',')
                    #print("fivalues last = i=%s %r" % (i, lines[i]))
                    i += 1
                    break
                i += 1
        #print("len(fivalues) = %s" % len(fivalues))
        return cls(set_id, fivalues)

    def write(self, spaces):
        """
        writes
        SET 80 = 3926, 3927, 3928, 4141, 4142, 4143, 4356, 4357, 4358, 4571,
             4572, 4573, 3323 THRU 3462, 3464 THRU 3603, 3605 THRU 3683,
             3910 THRU 3921, 4125 THRU 4136, 4340 THRU 4351

        """
        return write_set(self.set_id, self.value, spaces=spaces)

    def __repr__(self):
        """see `write`"""
        return write_set(self.set_id, self.value)

class SETMC(SET):
    """
    SETMC 121 = ACCE/99(T3),1200(T1),1399(R2)
    SETMC 222 = STRESS/134(22)
    SETMC 343 = ACCE/99(T3),1200(T1),1399(R2),STRESS/134(22)
    SETMC 122 = DISP/45(T1) 45(T2) 45(T3),
                 38(T1) 38(T2) 38(T3),
            VELO/45(T1) 45(T2) 45(T3),
                 38(T1) 38(T2) 38(T3),
            ACCE/45(T1) 45(T2) 45(T3),
                 38(T1) 38(T2) 38(T3)

    """
    type = 'SETMC'
    def __init__(self, set_id, values):
        super(SETMC, self).__init__(set_id, values)
