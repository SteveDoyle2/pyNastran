from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
from pyNastran.bdf.cards.baseCard import BaseCard



def ssq(*listA):
    """
    sum of squares
    .. note:: used for DEQATN
    """
    out = 0.
    for x in listA:
        out += x * x
    return out


def sum2(*listA):
    """
    sum of listA
    .. note:: used for DEQATN
    """
    return sum(listA)


def mod(x, y):
    """
    x%y
    .. note:: used for DEQATN
    """
    return x % y


def logx(x, y):
    """
    log base x of y
    .. note:: used for DEQATN
    """
    log(y, x)


def dim(x, y):
    """
    .. note:: used for DEQATN
    """
    return x - min(x, y)


def db(p, pref):
    """
    sound pressure in decibels
    would capitalize it, but you wouldnt be able to call the function...
    """
    return 20. * log(p / pref)


class DEQATN(BaseCard):  # needs work...
    type = 'DEQATN'

    def __init__(self, card=None, data=None, comment=''):
        if comment:
            self._comment = comment

        self.dtable = None
        new_card = ''
        found_none = False
        #print(card)
        line0 = card[0]
        name_eqid = line0[:16]
        #print('name_eqid = %r' % name_eqid)
        assert ',' not in name_eqid, name_eqid

        try:
            name, eq_id = name_eqid.split()
        except ValueError:
            msg = 'cannot split %r\n' % name_eqid
            msg += "Expected data of the form 'DEQATN  100'\n"
            msg += 'card=%s' % card
            raise ValueError(msg)

        self.equation_id = int(eq_id)

        line0_eq = line0[16:]
        eqs = [line0_eq] + card[1:]
        self.eqs = []
        neqs = len(eqs)
        is_join = False
        for i, eq in enumerate(eqs):
            if is_join:
                eq = eqi.rstrip() + eq.lstrip()
            eqi = eq.strip()
            if i == 0 and eqi == '':
                #self.eqs.append(eqi)
                continue

            if i == 0:
                # first line
                if eqi.endswith(';'):
                    eqi = eqi[:-1]
                    assert not eqi.endswith(';'), eq
                else:
                    is_join = True
                assert len(eqi) <= 56, eqi
            elif i != neqs-1:
                # mid line
                assert len(eqi) <= 64, eqi
                if eqi.endswith(';'):
                    eqi = eqi[:-1]
                    is_join = False
                    assert not eqi.endswith(';'), eq
                else:
                    is_join = True


            else:
                # last line
                pass
                is_join = False
            if not is_join:
                if '=' not in eqi:
                    raise SyntaxError('line=%r expected an equal sign' % eqi)
                self.eqs.append(eqi)
            #print(i, eqi)
        assert not is_join
        #assert len(eqs) <= 8, 'len(eqID)==%s' % (len(self.eqID))
        #self._setup_equation()

    def _setup_equation(self):
        default_values = {}
        if self.dtable is not None:
            self.dtable.default_values = {}
        func_name, nargs, func_str = fortran_to_python(self.eqs, default_values)
        self.func_name = func_name
        #print('**************', func_str)
        exec func_str
        print(locals().keys())
        func = locals()[func_name]
        setattr(self, func_name, func)
        #print(func)
        self.func = func
        self.nargs = nargs

    def cross_reference(self, model):
        # TODO: get deafults from DTABLE
        # TODO: get limits from DCONSTR
        self._setup_equation()
        self.dtable = model.dtable

    def uncross_reference(self):
        del self.func
        # del getattr(self, self.func_name)
        del self.func_name
        del self.nargs
        del self.dtable

    def evaluate(self, *args):
        #args2 = args[:self.nargs]
        #print('args =', args2)
        assert len(args) <= self.nargs, 'nargs=%s len(args)=%s; func_name=%s' % (self.nargs, len(args), self.func_name)
        return self.func(*args)
        #self.func(*args)

    def raw_fields(self):
        return [self.write_card()]
    def repr_fields(self):
        return self.raw_fields()

    def write_card(self, size=8, is_double=False):
        #self.evaluate(1, 2)
        equation_line0 = self.eqs[0]
        assert len(equation_line0) < 56, equation_line0
        msg = 'DEQATN  %8i%56s' % (self.equation_id, equation_line0)
        assert len(equation_line0) < 56, equation_line0
        for eq in self.eqs[1:]:
            msg += '        %64s\n' % eq
            assert len(eq) < 64, eq
        return msg

def fortran_to_python(lines, default_values):
    #print(lines)
    msg = ''
    #line0 = lines[0].lower()
    #print('line0=%r' % line0)
    #nlines = len(lines)
    for i, line in enumerate(lines):
        print('line=%r' % line)
        # line = line.upper()
        try:
            f, eq = line.split('=')
        except:
            raise SyntaxError('= not found in %r' % (line))
        f = f.strip()
        eq = eq.strip()

        if i == 0:
            #print('eq = %r' % eq)
            try:
                float(eq)
                is_float = True
            except ValueError:
                is_float = False

            if is_float:
                func_name, arguments = f.strip('(,)').split('(')
                func_name = func_name.strip(' ')
                variables = arguments.split(',')
                #print('func_name=%r' % func_name)
                val = float(eq)
                vals = []
                for var in variables:
                    if var in default_values:
                        vals.append('%s=%s' % (var, default_values[var]))
                    else:
                        vals.append('%s=%s' % (var, val))
                vals2 = ', '.join(vals)
                msg += 'def %s(%s):\n' % (func_name, vals2)
                msg += '    try:\n'
                for var in variables:
                    #msg += "    assert isinstance(%s, float), '%s is not a float; type(%s)=%s' % (%s)")
                    #msg += '        %s = float(%s)\n' % (var, var)
                    msg += '        %s = float(%s)\n' % (var, var)
                msg += '    except:\n'
                msg += '        print(locals())\n'
                msg += '        raise\n'
            else:
                func_name, arguments = f.strip('(,)').split('(')
                func_name = func_name.strip(' ')
                variables = arguments.split(',')
                msg += 'def %s:\n' % f
                for var in variables:
                    msg += '    %s = float(%s)\n' % (var, var)

                out = eq
                # if nlines > 1:
                    # msg += '    try:\n'
                    # msg += '        pass\n'
        else:
            out = f
            msg += '    %s = %s\n' % (out, eq)
        #print('  i=%s f=%r eq=%r' % (i, f, eq))
    #if nlines > 1:
        # msg += '    except:\n'
        # msg += '        print(locals())\n'
        # msg += '        raise\n'
    msg += '    return %s\n' % out

    print(msg)
    nargs = len(variables)
    return func_name, nargs, msg
