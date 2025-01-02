from __future__ import annotations
from typing import TYPE_CHECKING
from copy import deepcopy

from pyNastran.bdf.cards.deqatn import lines_to_eqs, write_deqatn, _setup_deqatn
from pyNastran.bdf.cards.base_card import BaseCard
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.dev.bdf_vectorized3.bdf import BDF


class DEQATN(BaseCard):  # needs work...
    """
    Design Equation Definition
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

    def __init__(self, equation_id: int, eqs: list[str],
                 ifile: int, comment: str=''):
        """
        Creates a DEQATN card

        Parameters
        ----------
        equation_id : int
            the id of the equation
        eqs : list[str]
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
        self.dtable_ref = None
        self.equation_id = equation_id
        self.eqs = eqs

        self.func = None
        self.func_str = ''
        self.func_name = ''
        self.nargs = 0

    #@classmethod
    #def _init_from_empty(cls):
        #equation_id = 1
        #eqs = []
        #return DEQATN(equation_id, eqs, comment='')

    def __deepcopy__(self, memo: dict[str, Any]):
        copy = type(self)(1, [], 0)
        memo[id(self)] = copy
        if self.comment:
            copy.comment = self.comment
        copy.equation_id = deepcopy(self.equation_id, memo)
        copy.eqs = deepcopy(self.eqs, memo)
        copy.func = deepcopy(self.func, memo)
        copy.func_str = deepcopy(self.func_str, memo)
        copy.func_name = deepcopy(self.func_name, memo)
        copy.nargs = deepcopy(self.nargs, memo)
        return copy

    @classmethod
    def add_card(cls, card: list[str], ifile: int, comment: str=''):
        """
        Adds a DEQATN card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : list[str]
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
        return DEQATN(equation_id, eqs, ifile, comment=comment)

    def _setup_equation(self, model: BDF) -> None:
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
        func_name, nargs, func_str = _setup_deqatn(
            self.equation_id, self.eqs,
            self.dtable_ref,
            str(self),
        )
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
        #self.dtable = model.dtable
        #self.dtable_ref = self.dtable
        self._setup_equation(model)

    #def safe_cross_reference(self, model: BDF) -> None:
        #self.cross_reference(model)

    #def uncross_reference(self) -> None:
        #"""Removes cross-reference links"""
        #del self.func
        ##del self.f
        ##del getattr(self, self.func_name)
        #setattr(self, self.func_name, None)
        #del self.func_name
        #del self.nargs
        #del self.dtable, self.dtable_ref

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

    def raw_fields(self) -> list[str]:
        return [self.write_card()]

    def repr_fields(self) -> list[str]:
        return self.raw_fields()

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        msg = write_deqatn(self.equation_id, self.eqs,
                           size=size, is_double=is_double)
        return msg
