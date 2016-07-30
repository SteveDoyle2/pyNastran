"""
Defines the DEQATN class and sub-functions.

The capitalization of the sub-functions is important.
"""
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
from six import exec_
from copy import deepcopy

from pyNastran.bdf.cards.base_card import BaseCard

from numpy import cos, sin, tan
from numpy import arcsin as asin, arccos as acos, arctan as atan, arctan2 as atan2
from numpy import arcsinh as asinh, arccosh as acosh, arctanh as atanh
# atan2h
from numpy import mean, exp, sqrt, square, sum
from numpy import log, log10, mod, abs
from numpy.linalg import norm

def rss(*args):  # good
    """2-norm; generalized magnitude of vector for N components"""
    return norm(args)

def avg(*args):
    """average"""
    return mean(args)

def ssq(*args):
    """sum of squares"""
    return square(args).sum()

def logx(x, y):
    """log base_x(y)"""
    return log(y**x) / log(x)

def dim(x, y):
    """positive difference"""
    return x - min(x, y)

def db(p, pref):
    """sound pressure in decibels"""
    return 20. * log(p / pref)

def atan2h(db, pref):
    raise NotImplementedError()

def invdb(db, pref):
    raise NotImplementedError()

def dba(p, pref, f):
    raise NotImplementedError()

def invdba(dba, pref, f):
    raise NotImplementedError()

class DEQATN(BaseCard):  # needs work...
    type = 'DEQATN'

    def __init__(self, name, equation_id, eqs, comment=''):
        if comment:
            self._comment = comment
        self.dtable = None
        self.func = None
        self.name = name
        self.equation_id = equation_id
        self.eqs = eqs

    @classmethod
    def add_card(cls, card, comment=''):
        new_card = ''
        found_none = False
        #print(card)
        line0 = card[0]
        if '\t' in line0:
            line0 = line0.expandtabs()

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

        equation_id = int(eq_id)

        line0_eq = line0[16:]
        eqs_temp = [line0_eq] + card[1:]
        eqs = []
        neqs = len(eqs_temp)
        is_join = False
        for i, eq in enumerate(eqs_temp):
            #print('i=%s join=%s eq=%r' % (i, is_join, eq))
            if is_join:
                eq = eqi.rstrip() + eq.lstrip()
            eqi = eq.strip().replace(' ', '')
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
                assert len(eqi) <= 64, 'len(eqi)=%s eq=%r' % (len(eqi), eqi)
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
                eqs.append(eqi)
            #print(i, eqi)
        #assert not is_join
        if is_join:
            eqs.append(eqi)
        _raw_eqs = deepcopy(eqs)  # TODO: temporary
        assert len(eqs) > 0, eqs
        #assert len(eqs) <= 8, 'len(eqID)==%s' % (len(eqID))
        #_setup_equation()
        return DEQATN(name, equation_id, eqs, comment=comment)

    def _setup_equation(self):
        default_values = {}
        if self.dtable is not None:
            self.dtable_ref.default_values = {}
        func_name, nargs, func_str = fortran_to_python(self.eqs, default_values)
        self.func_name = func_name
        #print('**************', func_str)
        exec_(func_str)
        #print(locals().keys())
        func = locals()[func_name]
        setattr(self, func_name, func)
        #print(func)
        self.func = func
        self.nargs = nargs

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        # TODO: get deafults from DTABLE
        # TODO: get limits from DCONSTR
        self.dtable = model.dtable
        self.dtable_ref = self.dtable
        self._setup_equation()

    def uncross_reference(self):
        del self.func
        del self.f
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
        eqs = split_equations(self.eqs)
        equation_line0 = eqs[0]
        #assert len(equation_line0) <= 56, equation_line0
        msg = 'DEQATN  %-8i%-56s\n' % (self.equation_id, equation_line0)
        assert len(equation_line0) <= 56, equation_line0
        for eq in eqs[1:]:
            msg += '        %-64s\n' % eq
            assert len(eq) <= 64, eq
        #print(msg)
        return msg

def split_equations(lines):
    # first line must be < 56
    # second line may be < 64
    lines2 = []
    for i, line in enumerate(lines):
        #print('-------------------------')
        # we'll add ; to the end of each line
        if i == 0:
            lines2 += split_equation([], line.strip() + ';', 56)
        else:
            lines2 += split_equation([], line.strip() + ';', 64)

    # remove the trailing semicolon
    lines2[-1] = lines2[-1][:-1]
    return lines2

def split_equation(lines_out, line, n, isplit=0):
    #print('n=%s -> line=%r len=%s' % (n, line, len(line)))
    if len(line) <= n:
        lines_out.append(line.strip())
        return lines_out
    # equation must be split
    line0 = line[:n][::-1].replace('**', '^')
    # fore, aft = line0.split('+-()*', 1)
    #print('line0 = %r; len=%s' % (str(line0[::-1]), len(line0)))
    out = {}
    for operator in ('+', '*', '^', '-', ')', ',', '='):
        if operator in line0:
            i = line0.index(operator)
            out[i] = operator

    try:
        imin = min(out)
    except ValueError:
        msg = "Couldn't find an operator ()+-/*= in %r\n" % line[n:]
        msg += 'line = %r' % line
        raise ValueError(msg)

    operator = out[imin]
    #print('operator = %r' % operator)
    fore, aft = line0.split(operator, 1)
    i = len(aft) + 1

    line_out = line[:i]
    #print('appending %r; len=%s' % (line_out, len(line_out)))
    #print('fore = %r' % fore[::-1])
    #print('aft  = %r' % aft[::-1])
    lines_out.append(line_out.replace('^', '**').strip())
    isplit += 1
    if isplit > 10:
        raise RuntimeError()
    lines_out = split_equation(lines_out, line[i:], n, isplit+1)
    return lines_out

def fortran_to_python_short(line, default_values):
    func_str = 'def func(args):\n'
    func_str += '    return %s(args)\n' % line.strip()
    d = {}
    print(func_str)
    exec_(func_str, globals(), d)
    return d['func']

def fortran_to_python(lines, default_values):
    #print(lines)
    msg = ''
    #line0 = lines[0].lower()
    #print('line0=%r' % line0)
    #nlines = len(lines)
    assert len(lines) > 0, lines
    is_eq_defined = False
    for i, line in enumerate(lines):
        #print('line=%r' % line)
        # line = line.upper()
        line = line.lower()
        try:
            f, eq = line.split('=')
        except:
            raise SyntaxError('= not found in %r' % (line))
        f = f.strip()
        eq = eq.strip()
        #print('f=%r eq=%r' % (f, eq))

        if i == 0:
            #print('eq = %r' % eq)
            try:
                float(eq)
                is_float = True
            except ValueError:
                is_float = False

            #print('is_float =', is_float)
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
                #print(msg)
                is_eq_defined = True
                out = eq
                # if nlines > 1:
                    # msg += '    try:\n'
                    # msg += '        pass\n'
        else:
            out = f
            msg += '    %s = %s\n' % (out, eq)
        #print('  i=%s f=%r eq=%r' % (i, f, eq))
    #if is_
    #if nlines > 1:
        # msg += '    except:\n'
        # msg += '        print(locals())\n'
        # msg += '        raise\n'
    if is_eq_defined:
        msg += '    return %s\n' % out
    else:
        msg = 'def %s(%s):\n' % (func_name, vals2)
        try:
            msg += '    return %s\n' % float(eq)
        except ValueError:
            msg += '    return %s\n' % str(eq)
            raise NotImplementedError(msg)

    #print(msg)
    nargs = len(variables)
    return func_name, nargs, msg
