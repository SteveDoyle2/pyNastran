from __future__ import annotations
from typing import TYPE_CHECKING
import numpy as np

from pyNastran.bdf.cards.bdf_tables import TABLED1, TABDMP1, read_table_lax, read_table
from pyNastran.bdf.field_writer_8 import set_blank_if_default, print_card_8
from pyNastran.bdf.bdf_interface.assign_type import (
    integer, integer_or_blank,
    string_or_blank,
)
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf import BDF
    from pyNastran.bdf.bdf_interface.bdf_card import BDFCard


class TABLED1_ZONA(TABLED1):
    """
    Dynamic Load Tabular Function, Form 1
    Defines a tabular function for use in generating frequency-dependent and
    time-dependent dynamic loads.

    +---------+------+-------+-------+--------+-----+-----+------+------+
    |    1    |   2  |   3   |   4   |    5   |  6  |  7  |  8   |   9  |
    +=========+======+=======+=======+========+=====+=====+======+======+
    | TABLED1 |  TID | XAXIS | YAXIS | EXTRAP |     |     |      |      |
    +---------+------+-------+-------+--------+-----+-----+------+------+
    |         |  x1  |  y1   |   x2  |   y2   | x3  | y3  | etc. | ENDT |
    +---------+------+-------+-------+--------+-----+-----+------+------+
    | TABLED1 |  32  |       |       |        |     |     |      |      |
    +---------+------+-------+-------+--------+-----+-----+------+------+
    |         | -3.0 |  6.9  |  2.0  |   5.6  | 3.0 | 5.6 | ENDT |      |
    +---------+------+-------+-------+--------+-----+-----+------+------+

    ..note:: EXTRAP is NX specific

    """
    type = 'TABLED1'

    @classmethod
    def _init_from_empty(cls):
        tid = 1
        x = [0., 1.]
        y = [0., 1.]
        return TABLED1(tid, x, y, xaxis='LINEAR', yaxis='LINEAR', extrap=0, comment='')

    def __init__(self, tid: int, x: np.ndarray, y: np.ndarray,
                 xaxis: str='LINEAR', yaxis: str='LINEAR',
                 extrap: int=0, comment: str=''):
        """
        Creates a TABLED1, which is a dynamic load card that is applied
        by the DAREA card

        Parameters
        ----------
        tid : int
            table id
        x : list[float]
            nvalues
        y : list[float]
            nvalues
        xaxis : str
            LINEAR, LOG
        yaxis : str
            LINEAR, LOG
        extrap : int; default=0
            Extrapolation method:
                0 : linear
                1 : constant
            .. note:: this is NX specific
        comment : str; default=''
            a comment for the card

        """
        Table.__init__(self)
        if comment:
            self.comment = comment
        self.tid = tid
        self.extrap = extrap
        self.x = np.asarray(x, dtype='float64')
        self.y = np.asarray(y, dtype='float64')
        self.xaxis = xaxis
        self.yaxis = yaxis
        assert self.xaxis in ['LINEAR', 'LOG'], f'xaxis={self.xaxis!r}'
        assert self.yaxis in ['LINEAR', 'LOG'], f'yaxis={self.yaxis!r}'

    @classmethod
    def add_card(cls, card: BDFCard, comment: str=''):
        """
        Adds a TABLED1 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        table_id = integer(card, 1, 'tid')
        xaxis = string_or_blank(card, 2, 'xaxis', default='LINEAR')
        yaxis = string_or_blank(card, 3, 'yaxis', default='LINEAR')
        extrap = integer_or_blank(card, 4, 'extrap', default=0)

        x, y = read_table(card, table_id, 'TABLED1')
        return TABLED1(table_id, x, y, xaxis=xaxis, yaxis=yaxis, extrap=extrap, comment=comment)

    @classmethod
    def add_card_lax(cls, card: BDFCard, comment: str=''):
        """
        Adds a TABLED1 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        table_id = integer(card, 1, 'tid')
        xaxis = string_or_blank(card, 2, 'xaxis', default='LINEAR')
        yaxis = string_or_blank(card, 3, 'yaxis', default='LINEAR')
        extrap = integer_or_blank(card, 4, 'yaxis', default=0)

        x, y = read_table_lax(card, table_id, 'TABLED1')
        return TABLED1(table_id, x, y, xaxis=xaxis, yaxis=yaxis, extrap=extrap, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment: str=''):
        table_id, extrap, xcode, ycode, xs, ys = data
        xaxis = _map_axis(xcode)
        yaxis = _map_axis(ycode)
        x = np.array(xs, dtype='float64')
        y = np.array(ys, dtype='float64')
        return TABLED1(table_id, x, y, xaxis=xaxis, yaxis=yaxis, extrap=extrap, comment=comment)

    def raw_fields(self):
        xy = []
        for xi, yi in zip(self.x, self.y):
            xy.extend([xi, yi])
        list_fields = ['TABLED1', self.tid, self.xaxis, self.yaxis, self.extrap,
                       None, None, None, None] + xy + ['ENDT']
        return list_fields

    def repr_fields(self):
        extrap = set_blank_if_default(self.extrap, 0)
        xy = []
        for xi, yi in zip(self.x, self.y):
            xy.extend([xi, yi])
        list_fields = ['TABLED1', self.tid, self.xaxis, self.yaxis, extrap,
                       None, None, None, None] + xy + ['ENDT']
        return list_fields

    def interpolate(self, x):
        if isinstance(x, float):
            x = [x]
        x = np.asarray(x)
        #nx = x.size
        #ny = self.y.size

        # xj follow xi
        i = np.searchsorted(self.x, x, side='left') - 1
        j = i + 1
        #k = np.where(j == ny)[0]

        # TODO: handle out of range errors
        xi = self.x[i]
        yi = self.y[i]
        try:
            xj = self.x[j]
            yj = self.y[j]
        except IndexError:
            #print('table.x = %s' % self.x)
            #print('table.y = %s' % self.y)
            #print('x = %s' % x)
            #print('yi = %s' % yi)
            return yi

        # TODO: could probably speed this up with log rules
        if self.xaxis == 'LINEAR' and self.yaxis == 'LINEAR':
            dx = xj - xi
            y = (xj - x) / dx * yi + (x - xi) / dx * yj
        elif self.xaxis == 'LOG' and self.yaxis == 'LINEAR':
            dx = np.log(xj / xi)
            y = np.log(xj / x) / dx * yi + np.log(x / xi) / dx * yj
        elif self.xaxis == 'LINEAR' and self.yaxis == 'LOG':
            dx = xj - xi
            lny = (xj - x) / dx * np.log(yi) + (x - xi) / dx * np.log(yj)
            y = np.exp(lny)
        elif self.xaxis == 'LOG' and self.yaxis == 'LOG':
            dx = np.log(xj / xi)
            lny = (xj - x) / dx * np.log(yi) + (x - xi) / dx * np.log(yj)
            y = np.exp(lny)
        else:
            raise NotImplementedError('xaxis=%r yaxis=%r' % (self.xaxis, self.yaxis))
        return y


class TABDMP1_ZONA(TABDMP1):
    @classmethod
    def add_card(cls, card: BDFCard, comment: str=''):
        """
        Adds a TABDMP1 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        table_id = integer(card, 1, 'tid')
        Type = string_or_blank(card, 2, 'Type', default='G')
        x, y = read_table_lax(card, table_id, 'TABDMP1')
        return TABDMP1(table_id, x, y, Type=Type, comment=comment)

