# coding: utf-8
"""
Defines the DEQATN class and sub-functions.

The capitalization of the sub-functions is important.

"""
from __future__ import annotations
from typing import TYPE_CHECKING
import numpy as np
from numpy import (
    cos, sin, tan, log, log10, mean, exp, sqrt, square, mod, abs, sum,
    arcsin as asin, arccos as acos, arctan as atan, arctan2 as atan2,
    arcsinh as asinh, arccosh as acosh, arctanh as atanh)
# atan2h
from numpy.linalg import norm  # type: ignore

from pyNastran.bdf.cards.base_card import BaseCard
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf import BDF

def pi(num):
    """weird way to multiply π by a number"""
    return np.pi * num

def rss(*args):  # good
    """2-norm; generalized magnitude of vector for N components"""
    return norm(args)

def avg(*args):
    """average"""
    return np.mean(args)

def ssq(*args):
    """sum of squares"""
    return np.square(args).sum()

def logx(x, y):
    """log base_x(y)"""
    return log(y**x) / log(x)

def dim(x, y):
    """positive difference"""
    return x - min(x, y)

def db(p, pref):
    """sound pressure in decibels"""
    return 20. * log(p / pref)

#def _Log(z):

def atan2h(x, y):
    """
    Hyperbolic arctangent
    >>> arctanh(z) = 1/2 Log((1+z)/(1-z))

    real:
        the following must be true: |x1| > |x2| and x2 ≠ 0.

    complex:
        x1 = a + bi
        x2 = b + di
        a = b = 0 and (sign of c) = (sign of d):
             the result is 0.
        a = b = 0 and (sign of c) ≠ (sign of d):
             the result is π.
        c = d = 0 and (sign of a) = (sign of b):
             the result is π/2.
        c = d = 0 and (sign of a) ≠ (sign of b):
             the result is −π/2
    """
    #integer_float_types = (int, np.int32, float, np.float32)
    #if isinstance(x, integer_float_types):
        #assert x >= 0, 'x=%s y=%s' % (x, y)
        #assert y > 0, 'x=%s y=%s' % (x, y)
    #return np.arctanh(x, y)
    raise NotImplementedError()

def invdb(dbi: float, pref: float) -> float:
    """inverse Db"""
    return 10. ** (dbi / 20. + log(pref))

def dba(p: float, pref: float, f: float) -> float:
    """
    sound pressure in decibels (perceived)

    Parameters
    ----------
    p : float
        structural responses or acoustic pressure
    f : float
        forcing frequency
    pref : float
        reference pressure

    Returns
    -------
    dbi : float
        acoustic pressure in Decibels

    """
    ta1, ta2 = _get_ta(f)
    return 20. * log(p / pref) + 10 * log(ta1) + 10. * log(ta2)

def invdba(dbai: float, pref: float, f: float) -> float:
    """
    Inverse Dba

    Parameters
    ----------
    dbai : float
        acoustic pressure in Decibels (perceived)
    f : float
        forcing frequency
    pref : float
        reference pressure

    Returns
    -------
    p : float
        structural responses or acoustic pressure

    """
    ta1, ta2 = _get_ta(f)
    #dbai = dba(p, pref, f)
    return 10. ** ((dbai - 10. * log(ta1) - 10. * log(ta2))/20)

def _get_ta(f: float) -> float:
    """gets the factors for dba, invdba"""
    k1 = 2.242882e16
    k3 = 1.562339
    p1 = 20.598997
    p2 = 107.65265
    p3 = 737.86223
    p4 = 12194.22
    ta1 = k3 * f**4 / ((f**2 + p2**2) * (f**2 + p3**2))
    ta2 = k1 * f**4 / ((f**2 + p1**2)**2 * (f**2 + p4**2)**2)
    return ta1, ta2

# we'll add _ to the beginning of these variables
BUILTINS = ['del', 'eval', 'yield', 'async', 'await', 'property',
            'slice', 'filter', 'map']


class DEQATN(BaseCard):  # needs work...
    """
    Design Equation Defintion
    Defines one or more equations for use in design sensitivity analysis.

    +--------+------+-----+-----+-----+-----+-------+-----+
    |    1   |   2  |  3  |  4  |  5  |  6  |    7  |  8  |
    +========+======+=====+=====+=====+=====+=======+=====+
    | DEQATN | EQID |            EQUATION                 |
    +--------+------+-------------------------------------+
    |        |                EQUATION (cont.)            |
    +--------+--------------------------------------------+
    """
    type = 'DEQATN'
    _properties = ['dtable']

    def __init__(self, equation_id, eqs, comment=''):
        """
        Creates a DEQATN card

        Parameters
        ----------
        equation_id : int
            the id of the equation
        eqs : List[str]
            the equations, which may overbound the field
            split them by a semicolon (;)
        comment : str; default=''
            a comment for the card

        DEQATN  41      F1(A,B,C,D,R) = A+B *C–(D**3 + 10.0) + sin(PI(1) * R)
                        + A**2 / (B - C); F = A + B - F1 * D

        def F1(A, B, C, D, R):
            F1 = A+B *C-(D**3 + 10.0) + sin(PI(1) * R) + A**2 / (B – C)
            F = A + B - F1 * D
            return F

        eqs = [
            'F1(A,B,C,D,R) = A+B *C–(D**3 + 10.0) + sin(PI(1) * R) + A**2 / (B – C)',
            'F = A + B – F1 * D',
        ]
        >>> deqatn = DEQATN(41, eq, comment='')

        """
        if comment:
            self.comment = comment
        self.dtable = None
        self.func = None
        self.equation_id = equation_id
        self.eqs = eqs
        self.func_str = ''

    @classmethod
    def _init_from_empty(cls):
        equation_id = 1
        eqs = []
        return DEQATN(equation_id, eqs, comment='')

    @classmethod
    def add_card(cls, card: List[str], comment: str=''):
        """
        Adds a DEQATN card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : List[str]
            this card is special and is not a ``BDFCard`` like other cards
        comment : str; default=''
            a comment for the card

        """
        #print(card)
        line0 = card[0]
        if '\t' in line0:
            line0 = line0.expandtabs()

        name_eqid = line0[:16]
        #print('name_eqid = %r' % name_eqid)
        assert ',' not in name_eqid, name_eqid

        try:
            name, eq_id = name_eqid.split()
            assert name.strip().upper() == 'DEQATN', card
        except ValueError:
            msg = 'cannot split %r\n' % name_eqid
            msg += "Expected data of the form 'DEQATN  100'\n"
            msg += 'card=%s' % card
            raise ValueError(msg)

        equation_id = int(eq_id)

        # combine the equations into a single organized block
        line0_eq = line0[16:]
        eqs_temp = [line0_eq] + card[1:]
        #eqs_temp2 = [line.replace(';;', ';') for line in eqs_temp]
        #for line in eqs_temp2:
            #print(line)
        eqs = lines_to_eqs(eqs_temp)
        return DEQATN(equation_id, eqs, comment=comment)

    def _setup_equation(self) -> None:
        """
        creates an executable equation object from self.eqs

        x = 10.
        >>> deqatn.func(x)
        42.0

        >>> deqatn.func_str
        def stress(x):
            x = float(x)
            return x + 32.

        """
        default_values = {}
        if self.dtable is not None:
            default_values = self.dtable_ref.default_values
        func_name, nargs, func_str = fortran_to_python(
            self.equation_id, self.eqs, default_values, str(self))
        self.func_str = func_str
        self.func_name = func_name
        try:
            exec(func_str)
        except SyntaxError:
            print(func_str)
            raise
        #print(locals().keys())
        func = locals()[func_name]
        setattr(self, func_name, func)
        #print(func)
        self.func = func
        self.nargs = nargs

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        # TODO: get defaults from DTABLE
        # TODO: get limits from DCONSTR
        self.dtable = model.dtable
        self.dtable_ref = self.dtable
        self._setup_equation()

    def safe_cross_reference(self, model: BDF) -> None:
        self.cross_reference(model)

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        del self.func
        #del self.f
        #del getattr(self, self.func_name)
        setattr(self, self.func_name, None)
        del self.func_name
        del self.nargs
        del self.dtable, self.dtable_ref

    def _verify(self, xref: bool) -> None:
        pass

    def evaluate(self, *args) -> float:
        """Makes a call to self.func"""
        #args2 = args[:self.nargs]
        #print('args =', args2)
        if len(args) > self.nargs:
            msg = 'len(args) > nargs\n'
            msg += 'nargs=%s len(args)=%s; func_name=%s' % (
                self.nargs, len(args), self.func_name)
            raise RuntimeError(msg)
        return self.func(*args)
        #self.func(*args)

    def raw_fields(self) -> List[str]:
        return [self.write_card()]

    def repr_fields(self) -> List[str]:
        return self.raw_fields()

    def write_card(self, size: int=8, is_double: bool=False) -> str:
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

def lines_to_eqs(eqs_in: List[str]) -> List[str]:
    """splits the equations"""
    eqs_wrapped = _split_equations_by_semicolon(eqs_in)
    eqs = _join_wrapped_equation_lines(eqs_in, eqs_wrapped)
    assert len(eqs) > 0, eqs
    return eqs

def _split_equations_by_semicolon(eqs_in: List[str]) -> List[str]:
    """helper for ``lines_to_eqs``"""
    eqs_temp_out = []
    nchars = 72 - 16
    for iline, eq in enumerate(eqs_in):
        if iline == 0:
            eq2 = eq[:nchars].strip(' \t\n')
        else:
            eq2 = eq.expandtabs()[8:nchars].strip(' \t\n')

        semicolon = ';' if eq2.rstrip().endswith(';') else ''
        eq2 = eq2.rstrip(' \t;')
        #nline = len(eq.rstrip('; \n')) + 16

        #print('eq2=%r' % eq2)
        if ';' in eq2:
            eq2s = eq2.split(';')
            eq_tempi = [eqi.strip() + ';' for eqi in eq2s if eqi.strip()]
            #for check_line in eq2s:
                #print(check_line)
                #_check_for_valid_line(check_line, eq)

            #print('eq_tempi = %r' % eq_tempi)
            #eq_tempi[-1] += semicolon
            eqs_temp_out += eq_tempi
        else:
            check_line = eq2 + semicolon
            #print('check_line = %r' % (check_line))
            #_check_for_valid_line(check_line, eq)
            eqs_temp_out.append(check_line)
        nchars = 72
    return eqs_temp_out

#def _check_for_valid_line(check_line, full_line):
    #if '=' not in check_line:
        #msg = (
            #'expected an equal sign (the first 8 characters are removed)\n'
            #'line     =%r\n'
            #'full_line=%r' % (check_line, full_line))
        #raise SyntaxError(msg)

def _join_wrapped_equation_lines(unused_eqs_temp_in, eqs_temp: List[str]) -> List[str]:
    """helper for ``lines_to_eqs``"""
    eqs = []
    neqs = len(eqs_temp)
    is_join = False
    eqi = ''
    for i, eq in enumerate(eqs_temp):
        #print(f'i={i} join={is_join} eq={eq!r}')
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
            #assert len(eqi) <= 64, 'len(eqi)=%s eq=%r' % (len(eqi), eqi)
            if eqi.endswith(';'):
                eqi = eqi[:-1]
                is_join = False
                assert not eqi.endswith(';'), eq
            else:
                is_join = True
        else:
            # last line
            is_join = False
        if not is_join:
            if '=' not in eqi:
                raise SyntaxError('line=%r expected an equal sign' % eqi)
            eqs.append(eqi)
        #print(i, eqi)
    #assert not is_join
    if is_join:
        eqs.append(eqi)
    return eqs

def split_equations(lines: List[str]) -> List[str]:
    """takes an overbounded DEQATN card and shortens it"""
    # first line must be < 56
    # second line may be < 64
    lines2 = []
    for i, line in enumerate(lines):
        #print('-------------------------')
        # we'll add ; to the end of each line
        if i == 0:
            lines2 += _split_equation([], line.strip() + ';', 56)
        else:
            lines2 += _split_equation([], line.strip() + ';', 64)

    # remove the trailing semicolon
    lines2[-1] = lines2[-1][:-1]
    return lines2

def _split_equation(lines_out: List[str], line: str, n: int,
                    isplit: int=0) -> List[str]:
    """
    Takes an overbounded DEQATN line and shortens it using recursion

    Parameters
    ----------
    lines_out : List[str]
        len(lines) = 0 : first iteration
        len(lines) = 1 : second iteration
    line : str
        the line to split
    n : int
        the maximum number of characters allowed
        the first line of the DEQATN has a different number of fields
        allowed vs. subsequent lines
    isplit : int; default=0
        the number of levels deep in the recursive function we are

    Returns
    -------
    lines_out : List[str]
        the long line broken into shorter lines

    """
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
    unused_fore, aft = line0.split(operator, 1)
    i = len(aft) + 1

    line_out = line[:i]
    #print('appending %r; len=%s' % (line_out, len(line_out)))
    #print('fore = %r' % fore[::-1])
    #print('aft  = %r' % aft[::-1])
    lines_out.append(line_out.replace('^', '**').strip())
    isplit += 1
    if isplit > 360:
        raise RuntimeError('Recursion issue splitting line; isplit=%i' %  isplit)
    lines_out = _split_equation(lines_out, line[i:], n, isplit+1)
    return lines_out

def fortran_to_python_short(line: str, unused_default_values: Any) -> Any:
    """the function used by the DRESP2"""
    func_str = 'def func(args):\n'
    func_str += '    return %s(args)\n' % line.strip()
    local_dict = {}
    exec(func_str, globals(), local_dict)
    return local_dict['func']

def split_to_equations(lines: List[str]) -> List[str]:
    """
    Splits a line like::

        b = a + z; c = 42

    into::

        b = a + z
        c = 42
    """
    equation_lines = []
    for line in lines:
        line = line.rstrip(' ;')
        if ';' in line:
            lines2 = line.split(';')
            equation_lines.extend(lines2)
        else:
            equation_lines.append(line)
    return equation_lines

def fortran_to_python(deqatn_id: int,
                      lines: List[str],
                      default_values: Dict[str, Union[float, np.ndarray]],
                      comment: str='') -> Tuple[str, int, str]:
    """
    Creates the python function

    Parameters
    ----------
    lines : List[str]
        the equations to write broken up by statement
    default_values : dict[name] = value
        the default values from the DTABLE card

    Returns
    -------
    func_name : str
        the name of the function
    nargs : int
        the number of variables to the function
    func_msg : str
        the python function

    def f(x, y=10.):
        '''
        $ deqatn
        DEQATN  1000    f(x,y) = x+y
        '''
        try:
            if isinstance(x, (int, float, str)):
                x = float(x)
            if isinstance(y, (int, float, str)):
                y = float(y)
        except:
            print(locals())
            raise
        f = x + y
        return f

    """
    func_msg = ''
    variables = []
    assert len(lines) > 0, lines

    equation_lines = split_to_equations(lines)
    for i, line in enumerate(equation_lines):
        #print('--------------------')
        line = line.lower()
        #func_msg += '#i=%s\n' % i
        assert ';' not in line, line
        try:
            # f(x, y) = 10.
            # f(x, y) = abs(x) + y
            # f = 42.
            f, eq = line.split('=')
        except ValueError:
            if '=' not in line:
                raise SyntaxError('= not found in %r' % (line))
            msg = 'only 1 = sign may be found a line\n'
            msg += 'line = %r\n' % line
            if len(lines) > 1:
                msg += 'lines:\n%s' % '\n'.join(lines)
            raise SyntaxError(msg)
        f = f.strip()
        eq = eq.strip().rstrip(';')

        #print('f=%r eq=%r' % (f, eq))
        for builtin in BUILTINS:
            if builtin == f:
                f = f.replace(builtin, builtin + '_')
            if builtin == eq:
                eq = eq.replace(builtin, builtin + '_')

        if i == 0:
            func_name, func_msg, variables = write_function_header(
                f, eq, default_values, comment)
            f = func_name  # return the value...
            func_msg += '    # i=0 write_function_header\n'
            #print(func_msg)
        else:
            out = f
            func_msg += '    %s = %s\n' % (out, eq)
            #print('out = %r' % out)
            #print('eq = %r' % eq)

    func_msg += '    return %s' % f
    #print(func_msg)
    if func_name in variables:
        raise RuntimeError(f'The function name {func_name!r} for DEQATN,{deqatn_id:d} '
                           f'must not also be an argument; arguments={variables}')
    nargs = len(variables)
    return func_name, nargs, func_msg


def write_function_header(func_header: str, eq: str,
                          default_values: Dict[str, float],
                          comment: str='') -> Tuple[str, str, List[str]]:
    """
    initializes the python function

    def f(x, y=10.):
        '''
        $ deqatn
        DEQATN  1000    f(x,y) = x+y
        '''
        try:
            if isinstance(x, (int, float, str)):
                x = float(x)
            if isinstance(y, (int, float, str)):
                y = float(y)
        except:
            print(locals())
            raise

    Parameters
    ----------
    f : str
        the function header
        f(a, b, c)
    eq : str
        the value on the other side of the equals sign (f=eq)
        1.
        max(a, b, c)
    default_values : dict[name] = value
        the default values from the DTABLE card

    Returns
    -------
    func_name : str
        the name of the function ``f``
    msg : str
        see above
    variables : List[str]
        the variables used by the equation header
        a, b, c

    """
    msg = ''

    try:
        float(eq)
        is_float = True
    except ValueError:
        is_float = False

    func_name, arguments = func_header.strip('(,)').split('(')
    func_name = func_name.strip(' ')
    variables = arguments.split(',')
    variables = ['_' + var if var in BUILTINS else var
                 for var in variables]

    if func_name in BUILTINS:
        func_name = '_' + func_name

    if is_float:
        # f(a,b,c) = 1.
        #
        # means
        #
        # def f(a,b,c):
        #     f = 1.
        #
        func_line = _write_function_line(func_name, variables, default_values)
    else:
        # f(a,b,c) = min(a,b,c)
        #
        # means
        #
        # def f(a,b,c):
        #     f = min(a,b,c)
        #
        func_line = _write_function_line(func_name, variables, default_values)
    msg += func_line
    msg += _write_comment(comment)
    msg += _write_variables(variables)
    for builtin in BUILTINS:
        ubuiltin = '_' + builtin
        if builtin in func_line and ubuiltin not in func_line:
            raise RuntimeError('cannot have an equation with %r\n%s' % (builtin, func_line))
        if builtin in variables and ubuiltin not in variables:
            raise RuntimeError('cannot have an equation with %r\n%s' % (builtin, variables))
            #import re
            #eq = 'YIELD_A_YIELD'
            #eq = '/YIELD'
            #p = re.compile(r"\byield\b", flags=re.IGNORECASE)
            #p2 = p.sub(eq,'_yield')
            #print('P2 = %r' % p2)
            #y = re.search(r"\byield\b", eq, flags=re.IGNORECASE)
            #if y is not None:
                #print('groups= ', y.groups())
                #for group in y.groups():
                    #print('group = %r' % group)
                #print(y.group(0))
                #print('***eq = %r' % eq)
    for builtin in BUILTINS:
        if builtin in eq and '_' + builtin not in eq:
            eq = eq.replace(builtin, '_'+builtin)

    msg += '    %s = %s\n' % (func_name, eq)
    return func_name, msg, variables

def _write_function_line(func_name: str, variables: List[str],
                         default_values: Dict[str, float]) -> str:
    """writes the ``def f(x, y, z=1.):`` part of the function"""
    vals = []
    is_default = False
    for var in variables:
        if var in BUILTINS:
            var += '_'

        if var in default_values:
            vals.append('%s=%s' % (var, default_values[var]))
            is_default = True
        else:
            vals.append('%s' % (var))
            if is_default:
                msg = 'default variables must be set at the end of the function\n'
                msg += 'variables = %s\n' % variables
                msg += 'default_values = %s' % default_values
                raise RuntimeError(msg)
    vals2 = ', '.join(vals)
    msg = 'def %s(%s):\n' % (func_name, vals2)
    return msg

def _write_comment(comment: str) -> str:
    """writes the deqatn to the comment block"""
    lines = comment.split('\n')
    msgi = '\n    '.join(lines)
    msg = '    """\n    %s"""\n' % msgi
    return msg

def _write_variables(variables: List[str]) -> str:
    """type checks the inputs"""
    msg = '    try:\n'
    for var in variables:
        if var in BUILTINS:
            var += '_'
        #msg += "    assert isinstance(%s, float), '%s is not a float; type(%s)=%s' % (%s)")
        #msg += '        %s = float(%s)\n' % (var, var)
        msg += '        if isinstance(%s, (int, str)):\n' % var
        msg += '            %s = float(%s)\n' % (var, var)
    msg += '    except:\n'
    msg += '        print(locals())\n'
    msg += '        raise\n'
    return msg
