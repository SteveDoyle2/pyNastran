"""
defines:
 - CBAR
 - CBARAO
 - BAROR
 - CBEAM3
 - CBEND

"""
# pylint: disable=R0904,R0902,E1101,E1103,C0111,C0302,C0103,W0101
from __future__ import annotations
from typing import Tuple, TYPE_CHECKING

import numpy as np
from numpy.linalg import norm

from pyNastran.utils.numpy_utils import integer_types
from pyNastran.bdf.field_writer_8 import set_blank_if_default
from pyNastran.bdf.cards.base_card import BaseCard, Element
from pyNastran.bdf.bdf_interface.assign_type import (
    integer, integer_or_blank, integer_double_or_blank, double_or_blank,
    integer_string_or_blank, string, integer_or_double,
    double)
from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.field_writer_16 import print_card_16
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.nptyping import NDArray3float, NDArray33float
    from cpylog import SimpleLogger
    from pyNastran.bdf.bdf import BDF


class LineElement(Element):  # CBAR, CBEAM, CBEAM3, CBEND
    def __init__(self):
        Element.__init__(self)
        self.pid_ref = None  # type: Optional[Any]
        #self.nodes_ref = None

    def C(self):
        """torsional constant"""
        if self.pid_ref is None:
            raise RuntimeError('Element eid=%i has not been '
                               'cross referenced.\n%s' % (self.eid, str(self)))
        return self.pid_ref.C()

    def Area(self):
        """returns the area of the element face"""
        raise NotImplementedError('implement self.Area() for %s' % self.type)

    def E(self):
        """returns the Young's Modulus, :math:`E`"""
        if self.pid_ref is None:
            raise RuntimeError('Element eid=%i has not been '
                               'cross referenced.\n%s' % (self.eid, str(self)))
        return self.pid_ref.mid_ref.E()

    def G(self):
        """returns the Shear Modulus, :math:`G`"""
        if self.pid_ref is None:
            raise RuntimeError('Element eid=%i has not been '
                               'cross referenced.\n%s' % (self.eid, str(self)))
        return self.pid_ref.mid_ref.G()

    def J(self):
        """returns the Polar Moment of Inertia, :math:`J`"""
        if self.pid_ref is None:
            raise RuntimeError('Element eid=%i has not been '
                               'cross referenced.\n%s' % (self.eid, str(self)))
        return self.pid_ref.J()

    def I11(self):
        """returns the Moment of Inertia, :math:`I_{11}`"""
        if self.pid_ref is None:
            raise RuntimeError('Element eid=%i has not been '
                               'cross referenced.\n%s' % (self.eid, str(self)))
        return self.pid_ref.I11()

    def I22(self):
        """returns the Moment of Inertia, :math:`I_{22}`"""
        if self.pid_ref is None:
            raise RuntimeError('Element eid=%i has not been '
                               'cross referenced.\n%s' % (self.eid, str(self)))
        return self.pid_ref.I22()

    def I12(self):
        """returns the Moment of Inertia, :math:`I_{12}`"""
        if self.pid_ref is None:
            raise RuntimeError('Element eid=%i has not been '
                               'cross referenced.\n%s' % (self.eid, str(self)))
        return self.pid_ref.I12()

    def Nu(self):
        """Get Poisson's Ratio, :math:`\nu`"""
        if self.pid_ref is None:
            raise RuntimeError('Element eid=%i has not been '
                               'cross referenced.\n%s' % (self.eid, str(self)))
        return self.pid_ref.mid_ref.nu

    def Rho(self):
        """Get the material density, :math:`\rho`"""
        #print(str(self.pid), type(self.pid))
        #raise NotImplementedError('implement self.Rho() for %s' % self.type)
        if self.pid_ref is None:
            raise RuntimeError('Element eid=%i has not been '
                               'cross referenced.\n%s' % (self.eid, str(self)))
        return self.pid_ref.mid_ref.rho

    def Nsm(self):
        """Placeholder method for the non-structural mass, :math:`nsm`"""
        raise NotImplementedError('implement self.Area() for %s' % self.type)

    def MassPerLength(self):
        """Get the mass per unit length, :math:`\frac{m}{L}`"""
        if self.pid_ref is None:
            raise RuntimeError('Element eid=%i has not been '
                               'cross referenced.\n%s' % (self.eid, str(self)))
        return self.pid_ref.MassPerLength()

    def Mass(self):
        r"""
        Get the mass of the element.

        .. math:: m = \left( \rho A + nsm \right) L

        """
        L = self.Length()
        mass = L * self.MassPerLength()
        #try:
            #mass = (self.Rho() * self.Area() + self.Nsm()) * L
        #except TypeError:
            #msg = 'TypeError on eid=%s pid=%s:\n' % (self.eid, self.Pid())
            #msg += 'rho = %s\narea = %s\nnsm = %s\nL = %s' % (self.Rho(),
            #                                                  self.Area(),
            #                                                  self.Nsm(), L)
            #raise TypeError(msg)

        return mass

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.nodes = self.node_ids
        self.pid = self.Pid()
        self.nodes_ref = None
        self.pid_ref = None

    def Length(self):
        r"""
        Gets the length, :math:`L`, of the element.

        .. math:: L = \sqrt{  (n_{x2}-n_{x1})^2+(n_{y2}-n_{y1})^2+(n_{z2}-n_{z1})^2  }

        """
        L = norm(self.nodes_ref[1].get_position() - self.nodes_ref[0].get_position())
        return L

    def get_edge_ids(self):
        """
        Return the edge IDs
        """
        node_ids = self.node_ids
        return [(node_ids[0], node_ids[1])]


class BAROR(BaseCard):
    """
    +-------+---+-----+---+---+-------+-----+-------+------+
    |   1   | 2 |  3  | 4 | 5 |   6   |  7  |   8   |  9   |
    +=======+===+=====+===+===+=======+=====+=======+======+
    | BAROR |   | PID |   |   | G0/X1 |  X2 |  X3   | OFFT |
    +-------+---+-----+---+---+-------+-----+-------+------+
    | BAROR |   | 39  |   |   |  0.6  | 2.9 | -5.87 | GOG  |
    +-------+---+-----+---+---+-------+-----+-------+------+

    """
    type = 'BAROR'

    @classmethod
    def _init_from_empty(cls):
        pid = 1
        is_g0 = True
        g0 = 1
        x = None
        return BAROR(pid, is_g0, g0, x, offt='GGG', comment='')

    def __init__(self, pid, is_g0, g0, x, offt='GGG', comment=''):
        BaseCard.__init__(self)
        if comment:
            self.comment = comment
        if x is None:
            x = np.array([0., 0., 0.])
        self.n = 0
        self.pid = pid
        self.is_g0 = is_g0
        self.g0 = g0
        self.x = x
        self.offt = offt
        #if isinstance(offt, integer_types):
            #raise NotImplementedError('the integer form of offt is not supported; offt=%s' % offt)

    @classmethod
    def add_card(cls, card, comment=''):
        pid = integer_or_blank(card, 2, 'pid')

        # x / g0
        field5 = integer_double_or_blank(card, 5, 'g0_x1', 0.)
        if isinstance(field5, integer_types):
            is_g0 = True
            g0 = field5
            x = [0., 0., 0.]
        elif isinstance(field5, float):
            is_g0 = False
            g0 = None
            x = np.array([field5,
                          double_or_blank(card, 6, 'x2', 0.),
                          double_or_blank(card, 7, 'x3', 0.)],
                         dtype='float64')
        else:
            raise NotImplementedError('BAROR field5 = %r' % field5)
        offt = integer_string_or_blank(card, 8, 'offt', 'GGG')
        assert len(card) <= 9, 'len(BAROR card) = %i\ncard=%s' % (len(card), card)
        return BAROR(pid, is_g0, g0, x, offt=offt, comment=comment)

    def raw_fields(self):
        """
        Gets the fields of the card in their full form
        """
        list_fields = ['BAROR', None, None] + self.x.tolist() + [self.offt]
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)

class CBARAO(BaseCard):
    """
    Per MSC 2016.1
    +--------+------+-------+------+-----+--------+-----+----+----+
    |   1    |  2   |   3   |  4   |  5  |    6   |  7  | 8  |  9 |
    +========+======+=======+======+=====+========+=====+====+====+
    | CBARAO | EID  | SCALE |  X1  | X2  |  X3    | X4  | X5 | X6 |
    +--------+------+-------+------+-----+--------+-----+----+----+
    | CBARAO | 1065 |  FR   | 0.2  | 0.4 |  0.6   | 0.8 |    |    |
    +--------+------+-------+------+-----+--------+-----+----+----+

    Alternate form (not supported):
    +--------+------+-------+------+-----+--------+-----+----+----+
    |   1    |  2   |   3   |  4   |  5  |    6   |  7  | 8  |  9 |
    +========+======+=======+======+=====+========+=====+====+====+
    | CBARAO | EID  | SCALE | NPTS | X1  | DELTAX |     |    |    |
    +--------+------+-------+------+-----+--------+-----+----+----+
    | CBARAO | 1065 |  FR   |  4   | 0.2 |  0.2   |     |    |    |
    +--------+------+-------+------+-----+--------+-----+----+----+

    """
    type = 'CBARAO'

    @classmethod
    def _init_from_empty(cls):
        eid = 1
        scale = 'FR'
        x = [0.5]
        return CBARAO(eid, scale, x, comment='')

    def __init__(self, eid, scale, x, comment=''):
        """
        Creates a CBARAO card, which defines additional output locations
        for the CBAR card.

        It also changes the OP2 element type from a CBAR-34 to a CBAR-100.
        However, it is ignored if there are no PLOAD1s in the model.
        Furthermore, the type is changed for the whole deck, regardless of
        whether there are PLOAD1s in the other load cases.

        Parameters
        ----------
        eid : int
            element id
        scale : str
            defines what x means
            LE : x is in absolute coordinates along the bar
            FR : x is in fractional
        x : List[float]
            the additional output locations (doesn't include the end points)
            len(x) <= 6
        comment : str; default=''
            a comment for the card

        MSC only

        """
        if comment:
            self.comment = comment
        self.eid = eid
        self.scale = scale
        self.x = np.unique(x).tolist()

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a CBARAO card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        eid = integer(card, 1, 'eid')
        scale = string(card, 2, 'scale')
        x1_npoints = integer_or_double(card, 3, 'x1/npoints')
        if isinstance(x1_npoints, integer_types):
            npoints = x1_npoints
            assert 0 < npoints < 7, 'CBARAO npoints=%r must be 1-6' % npoints
            x1 = double(card, 4, 'x1')
            delta_x = double(card, 5, 'delta_x')
            x = np.linspace(x1, x1 + delta_x * (npoints-1), num=npoints)
            assert len(x) == npoints, x
        else:
            x = [
                x1_npoints,
                double_or_blank(card, 4, 'x2'),
                double_or_blank(card, 5, 'x3'),
                double_or_blank(card, 6, 'x4'),
                double_or_blank(card, 7, 'x5'),
                double_or_blank(card, 8, 'x6'),
            ]
            x = [xi for xi in x if xi is not None]
        assert len(card) <= 9, 'len(CBARAO card) = %i\ncard=%s' % (len(card), card)
        return CBARAO(eid, scale, x, comment=comment)

    def _verify(self, xref):
        pass

    def raw_fields(self):
        list_fields = ['CBARAO', self.eid, self.scale] + self.x
        return list_fields

    def repr_fields(self):
        return self.raw_fields()

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)

class CBAR(LineElement):
    """
    +-------+-----+-----+-----+-----+-----+-----+-----+------+
    |   1   |  2  |  3  |  4  |  5  |  6  |  7  |  8  |   9  |
    +=======+=====+=====+=====+=====+=====+=====+=====+======+
    | CBAR  | EID | PID | GA  | GB  | X1  | X2  | X3  | OFFT |
    +-------+-----+-----+-----+-----+-----+-----+-----+------+
    |       | PA  | PB  | W1A | W2A | W3A | W1B | W2B | W3B  |
    +-------+-----+-----+-----+-----+-----+-----+-----+------+

    or

    +-------+-----+-----+-----+-----+-----+-----+-----+------+
    |   1   |  2  |  3  |  4  |  5  |  6  |  7  |  8  |   9  |
    +=======+=====+=====+=====+=====+=====+=====+=====+======+
    | CBAR  | EID | PID | GA  | GB  | G0  |     |     | OFFT |
    +-------+-----+-----+-----+-----+-----+-----+-----+------+
    |       | PA  | PB  | W1A | W2A | W3A | W1B | W2B | W3B  |
    +-------+-----+-----+-----+-----+-----+-----+-----+------+

    +-------+-------+-----+-------+-------+--------+-------+-------+-------+
    |   1   |   2   |  3  |   4   |   5   |    6   |   7   |   8   |   9   |
    +=======+=======+=====+=======+=======+========+=======+=======+=======+
    |  CBAR |   2   |  39 |   7   |   6   |  105   |       |       |  GGG  |
    +-------+-------+-----+-------+-------+--------+-------+-------+-------+
    |       |       | 513 |  0.0  |  0.0  |    -9. |  0.0  |  0.0  |   -9. |
    +-------+-------+-----+-------+-------+--------+-------+-------+-------+

    """
    type = 'CBAR'
    _field_map = {
        1: 'eid', 2:'pid', 3:'ga', 4:'gb',
        8:'offt', 9:'pa', 10:'pb',
        #, 'W1A', 'W2A', 'W3A', 'W1B', 'W2B', 'W3B'
    }
    def update_by_cp_name(self, cp_name, value):
        if cp_name == 'W1A':
            self.wa[0] = value
        elif cp_name == 'W2A':
            self.wa[1] = value
        elif cp_name == 'W3A':
            self.wa[2] = value
        elif cp_name == 'W1B':
            self.wb[0] = value
        elif cp_name == 'W2B':
            self.wb[1] = value
        elif cp_name == 'W3B':
            self.wb[2] = value
        elif cp_name == 'X1':
            self.x[0] = value
        elif cp_name == 'X2':
            self.x[1] = value
        elif cp_name == 'X3':
            self.x[2] = value
        else:  # pragma: no cover
            msg = 'CBAR: cp_name=%r must be added to update_by_cp_name' % cp_name
            raise NotImplementedError(msg)

    def _update_field_helper(self, n, value):
        if n == 11:
            self.wa[0] = value
        elif n == 12:
            self.wa[1] = value
        elif n == 13:
            self.wa[2] = value
        elif n == 14:
            self.wb[0] = value
        elif n == 15:
            self.wb[1] = value
        elif n == 16:
            self.wb[2] = value
        else:
            if self.g0 is not None:
                if n == 5:
                    self.g0 = value
                else:
                    raise KeyError('Field %r=%r is an invalid %s entry.' % (n, value, self.type))
            else:
                if n == 5:
                    self.x[0] = value
                elif n == 6:
                    self.x[1] = value
                elif n == 7:
                    self.x[2] = value
                else:
                    raise KeyError('Field %r=%r is an invalid %s entry.' % (n, value, self.type))

    @classmethod
    def export_to_hdf5(cls, h5_file, model, eids):
        """exports the elements in a vectorized way"""
        #comments = []
        pids = []
        nodes = []
        x = []
        g0 = []
        offt = []
        unused_bit = []
        pa = []
        pb = []
        wa = []
        wb = []
        nan = np.full(3, np.nan)
        encoding = model._encoding
        for eid in eids:
            element = model.elements[eid]
            #comments.append(element.comment)
            pids.append(element.pid)
            nodes.append(element.nodes)
            if element.g0 is None:
                x.append(element.x)
                g0.append(-1)
            else:
                x.append(nan)
                g0.append(element.g0)

            offti = element.offt
            if isinstance(offti, integer_types):
                offti = str(offti)
            offt.append(offti.encode(encoding))
            pa.append(element.pa)
            pb.append(element.pb)
            wa.append(element.wa)
            wb.append(element.wb)
        #h5_file.create_dataset('_comment', data=comments)
        h5_file.create_dataset('eid', data=eids)
        h5_file.create_dataset('nodes', data=nodes)
        h5_file.create_dataset('pid', data=pids)
        #print('x =', x)
        #print('g0 =', g0)
        h5_file.create_dataset('x', data=x)
        h5_file.create_dataset('g0', data=g0)
        h5_file.create_dataset('offt', data=offt)

        h5_file.create_dataset('pa', data=pa)
        h5_file.create_dataset('pb', data=pb)

        h5_file.create_dataset('wa', data=wa)
        h5_file.create_dataset('wb', data=wb)

    def __init__(self, eid, pid, nids,
                 x, g0, offt='GGG',
                 pa=0, pb=0, wa=None, wb=None, comment=''):
        """
        Adds a CBAR card

        Parameters
        ----------
        pid : int
            property id
        mid : int
            material id
        nids : List[int, int]
            node ids; connected grid points at ends A and B
        x : List[float, float, float]
            Components of orientation vector, from GA, in the displacement
            coordinate system at GA (default), or in the basic coordinate system
        g0 : int
            Alternate method to supply the orientation vector using grid
            point G0. Direction of is from GA to G0. is then transferred
            to End A
        offt : str; default='GGG'
            Offset vector interpretation flag
        pa / pb : int; default=0
            Pin Flag at End A/B.  Releases the specified DOFs
        wa / wb : List[float, float, float]
            Components of offset vectors from the grid points to the end
            points of the axis of the shear center
        comment : str; default=''
            a comment for the card

        """
        LineElement.__init__(self)
        if comment:
            self.comment = comment
        if wa is None:
            wa = np.zeros(3, dtype='float64')
        else:
            wa = np.asarray(wa)
        if wb is None:
            wb = np.zeros(3, dtype='float64')
        else:
            wb = np.asarray(wb)

        if x is not None:
            x = np.asarray(x)
        if isinstance(offt, str):
            offt = offt.replace('E', 'O')
            offt = int(offt) if offt.isdigit() else offt

        self.eid = eid
        self.pid = pid
        self.x = x
        self.g0 = g0
        self.ga = nids[0]
        self.gb = nids[1]
        self.offt = offt
        self.pa = pa
        self.pb = pb
        self.wa = wa
        self.wb = wb
        self.pid_ref = None
        self.ga_ref = None
        self.gb_ref = None
        self.g0_ref = None
        self.g0_vector = None

    def validate(self):
        msg = ''
        if self.x is None:
            if not isinstance(self.g0, integer_types):
                msg += 'CBAR eid=%s: x is None, so g0=%s must be an integer' % (self.eid, self.g0)
        else:
            if not isinstance(self.x, (list, np.ndarray)):
                msg += 'CBAR eid=%s: x=%s and g0=%s, so x must be a list; type(x)=%s' % (
                    self.eid, self.x, self.g0, type(self.x))
        if msg:
            raise ValueError(msg)

        if self.g0 in [self.ga, self.gb]:
            msg = 'G0=%s cannot be GA=%s or GB=%s' % (self.g0, self.ga, self.gb)
            raise RuntimeError(msg)

        if isinstance(self.offt, integer_types):
            assert self.offt in [1, 2, 21, 22, 41, 42], 'invalid offt; offt=%i' % self.offt
            #raise NotImplementedError('invalid offt; offt=%i' % self.offt)
        elif not isinstance(self.offt, str):
            raise SyntaxError('invalid offt expected a string of length 3 '
                              'offt=%r; Type=%s' % (self.offt, type(self.offt)))

        if not isinstance(self.offt, integer_types):
            check_offt(self)

    @classmethod
    def add_card(cls, card, baror=None, comment=''):
        """
        Adds a CBAR card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        beamor : BAROR() or None
            defines the defaults
        comment : str; default=''
            a comment for the card

        """
        eid = integer(card, 1, 'eid')
        pid_default = eid
        x1_default, x2_default, x3_default = 0., 0., 0.
        offt_default = 'GGG'
        if baror is not None:
            if baror.pid is not None:
                pid_default = baror.pid
            if baror.x is None:
                x1_default = baror.g0
                x2_default = None
                x3_default = None
            else:
                x1_default, x2_default, x3_default = baror.x
            offt_default = baror.offt

        pid = integer_or_blank(card, 2, 'pid', pid_default)
        ga = integer(card, 3, 'ga')
        gb = integer(card, 4, 'gb')
        x, g0 = init_x_g0(card, eid, x1_default, x2_default, x3_default)

        # doesn't exist in NX nastran
        offt = integer_string_or_blank(card, 8, 'offt', offt_default)
        #print('cls.offt = %r' % (cls.offt))

        pa = integer_or_blank(card, 9, 'pa', 0)
        pb = integer_or_blank(card, 10, 'pb', 0)

        wa = np.array([double_or_blank(card, 11, 'w1a', 0.0),
                       double_or_blank(card, 12, 'w2a', 0.0),
                       double_or_blank(card, 13, 'w3a', 0.0)], dtype='float64')

        wb = np.array([double_or_blank(card, 14, 'w1b', 0.0),
                       double_or_blank(card, 15, 'w2b', 0.0),
                       double_or_blank(card, 16, 'w3b', 0.0)], dtype='float64')
        assert len(card) <= 17, 'len(CBAR card) = %i\ncard=%s' % (len(card), card)
        return CBAR(eid, pid, [ga, gb], x, g0,
                    offt, pa, pb, wa, wb, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        #: .. todo:: verify
        #data = [[eid,pid,ga,gb,pa,pb,w1a,w2a,w3a,w1b,w2b,w3b],[f,g0]]
        #data = [[eid,pid,ga,gb,pa,pb,w1a,w2a,w3a,w1b,w2b,w3b],[f,x1,x2,x3]]

        main = data[0]
        flag = data[1][0]
        if flag in [0, 1]:
            g0 = None
            x = np.array([data[1][1],
                          data[1][2],
                          data[1][3]], dtype='float64')
        else:
            g0 = data[1][1]
            x = None

        eid = main[0]
        pid = main[1]
        ga = main[2]
        gb = main[3]
        #self.offt = str(data[4]) # GGG
        offt = 'GGG'  #: .. todo:: offt can be an integer; translate to char
        pa = main[4]
        pb = main[5]

        wa = np.array([main[6], main[7], main[8]], dtype='float64')
        wb = np.array([main[9], main[10], main[11]], dtype='float64')
        return CBAR(eid, pid, [ga, gb], x, g0,
                    offt, pa, pb, wa, wb, comment=comment)

    def _verify(self, xref):
        eid = self.eid
        unused_pid = self.Pid()
        unused_edges = self.get_edge_ids()
        if xref:  # True
            assert self.pid_ref.type in ['PBAR', 'PBARL', 'PBRSECT'], '%s%s' % (self, self.pid_ref)
            mid = self.Mid()
            A = self.Area()
            nsm = self.Nsm()
            mpl = self.MassPerLength()
            L = self.Length()
            mass = self.Mass()
            assert isinstance(mid, int), 'mid=%r' % mid
            assert isinstance(nsm, float), 'nsm=%r' % nsm
            assert isinstance(A, float), 'eid=%s A=%r' % (eid, A)
            assert isinstance(L, float), 'eid=%s L=%r' % (eid, L)
            assert isinstance(mpl, float), 'eid=%s mass_per_length=%r' % (eid, mpl)
            assert isinstance(mass, float), 'eid=%s mass=%r' % (eid, mass)
            assert L > 0.0, 'eid=%s L=%s' % (eid, L)

    def Mid(self):
        if self.pid_ref is None:
            msg = 'Element eid=%i has not been cross referenced.\n%s' % (self.eid, str(self))
            raise RuntimeError(msg)
        return self.pid_ref.Mid()

    def Area(self):
        if self.pid_ref is None:
            msg = 'Element eid=%i has not been cross referenced.\n%s' % (self.eid, str(self))
            raise RuntimeError(msg)
        A = self.pid_ref.Area()
        assert isinstance(A, float)
        return A

    def J(self):
        if self.pid_ref is None:
            msg = 'Element eid=%i has not been cross referenced.\n%s' % (self.eid, str(self))
            raise RuntimeError(msg)
        j = self.pid_ref.J()
        if not isinstance(j, float):
            msg = 'J=%r must be a float; CBAR eid=%s pid=%s pidType=%s' % (
                j, self.eid, self.pid_ref.pid, self.pid_ref.type)
            raise TypeError(msg)
        return j

    def Length(self):
        # TODO: consider w1a and w1b in the length formulation
        L = norm(self.gb_ref.get_position() - self.ga_ref.get_position())
        assert isinstance(L, float)
        return L

    def Nsm(self):
        if self.pid_ref is None:
            msg = 'Element eid=%i has not been cross referenced.\n%s' % (self.eid, str(self))
            raise RuntimeError(msg)
        nsm = self.pid_ref.Nsm()
        assert isinstance(nsm, float)
        return nsm

    def I1(self):
        if self.pid_ref is None:
            msg = 'Element eid=%i has not been cross referenced.\n%s' % (self.eid, str(self))
            raise RuntimeError(msg)
        return self.pid_ref.I1()

    def I2(self):
        if self.pid_ref is None:
            msg = 'Element eid=%i has not been cross referenced.\n%s' % (self.eid, str(self))
            raise RuntimeError(msg)
        return self.pid_ref.I2()

    def Centroid(self):
        if self.pid_ref is None:
            msg = 'Element eid=%i has not been cross referenced.\n%s' % (self.eid, str(self))
            raise RuntimeError(msg)
        return (self.ga_ref.get_position() + self.gb_ref.get_position()) / 2.

    def center_of_mass(self):
        return self.Centroid()

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        #if self.g0:
        #    self.x = nodes[self.g0].get_position() - nodes[self.ga].get_position()
        msg = ', which is required by CBAR eid=%s' % (self.eid)
        self.ga_ref = model.Node(self.ga, msg=msg)
        self.gb_ref = model.Node(self.gb, msg=msg)
        self.pid_ref = model.Property(self.pid, msg=msg)
        if model.is_nx:
            assert self.offt == 'GGG', 'NX only support offt=GGG; offt=%r' % self.offt

        if self.g0:
            self.g0_ref = model.nodes[self.g0]
            self.g0_vector = self.g0_ref.get_position() - self.ga_ref.get_position()
        else:
            self.g0_vector = self.x

    def safe_cross_reference(self, model, xref_errors):
        msg = ', which is required by CBAR eid=%s' % (self.eid)
        self.ga_ref = model.Node(self.ga, msg=msg)
        self.gb_ref = model.Node(self.gb, msg=msg)
        self.nodes_ref = [self.ga_ref, self.gb_ref]
        self.pid_ref = model.safe_property(self.pid, self.eid, xref_errors, msg=msg)

        if self.g0:
            try:
                self.g0_ref = model.nodes[self.g0]
                self.g0_vector = self.g0_ref.get_position() - self.ga_ref.get_position()
            except KeyError:
                model.log.warning('Node=%s%s' % (self.g0, msg))
        else:
            self.g0_vector = self.x

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.pid = self.Pid()
        self.ga = self.Ga()
        self.gb = self.Gb()
        self.g0 = self.G0()
        self.ga_ref = None
        self.gb_ref = None
        self.g0_ref = None
        self.pid_ref = None

    def Ga(self):
        """gets Ga/G1"""
        if self.ga_ref is None:
            return self.ga
        return self.ga_ref.nid

    def Gb(self):
        """gets Gb/G2"""
        if self.gb_ref is None:
            return self.gb
        return self.gb_ref.nid

    def G0(self):
        """gets G0"""
        if self.g0_ref is None:
            return self.g0
        return self.g0_ref.nid

    def get_x_g0_defaults(self):
        """
        X and G0 compete for the same fields, so the method exists to
        make it easier to write the card

        Returns
        -------
        x_g0 : varies
            g0 : List[int, None, None]
            x : List[float, float, float]

        Notes
        -----
        Used by CBAR and CBEAM

        """
        if self.g0 is not None:
            return self.G0(), None, None
        #print('x =', self.x)
        #print('g0 =', self.g0)
        #x1 = set_blank_if_default(self.x[0], 0.0)
        #x2 = set_blank_if_default(self.x[1], 0.0)
        #x3 = set_blank_if_default(self.x[2], 0.0)
        return list(self.x)

    def get_axes(self, model):
        """
        Gets the axes of a CBAR/CBEAM, while respecting the OFFT flag.

        Notes
        -----
        :func:`pyNastran.bdf.cards.elements.bars.rotate_v_wa_wb` for a
        description of the OFFT flag.

        is_passed: bool
        out: (wa, wb, ihat, jhat, khat)
        """
        is_failed = True

        check_offt(self)
        is_failed = True
        ihat = None
        yhat = None
        zhat = None

        eid = self.eid
        (nid1, nid2) = self.node_ids
        node1 = model.nodes[nid1]
        node2 = model.nodes[nid2]
        xyz1 = node1.get_position()
        xyz2 = node2.get_position()

        elem = model.elements[eid]
        pid_ref = elem.pid_ref
        if pid_ref is None:
            pid_ref = model.Property(elem.pid)
        assert not isinstance(pid_ref, integer_types), elem

        is_failed, (wa, wb, ihat, yhat, zhat) = self.get_axes_by_nodes(
            model, pid_ref, node1, node2, xyz1, xyz2, model.log)
        return is_failed, (wa, wb, ihat, yhat, zhat)

    def get_axes_by_nodes(self, model, pid_ref, node1, node2, xyz1, xyz2, log):
        """
        Gets the axes of a CBAR/CBEAM, while respecting the OFFT flag.

        Notes
        -----
        :func:`pyNastran.bdf.cards.elements.bars.rotate_v_wa_wb` for a
        description of the OFFT flag.

        """
        #TODO: not integrated with CBAR yet...

        is_failed = True
        eid = self.eid
        #centroid = (n1 + n2) / 2.
        #i = n2 - n1
        #Li = norm(i)
        #ihat = i / Li

        elem = self
        #(nid1, nid2) = elem.node_ids
        #node1 = model.nodes[nid1]
        #node2 = model.nodes[nid2]
        #xyz1 = node1.get_position()
        #xyz2 = node2.get_position()

        # wa/wb are not considered in i_offset
        # they are considered in ihat
        i = xyz2 - xyz1
        Li = norm(i)
        if Li == 0.:
            msg = 'xyz1=%s xyz2=%s\n%s' % (xyz1, xyz2, self)
            raise ValueError(msg)
        i_offset = i / Li

        unused_v, wa, wb, xform = rotate_v_wa_wb(
            model, elem,
            xyz1, xyz2, node1, node2,
            i_offset, i, eid, Li, log)
        if wb is None:
            # one or more of v, wa, wb are bad

            # xform is xform_offset...assuming None
            ihat = None
            yhat = None
            zhat = None
            return is_failed, (wa, wb, ihat, yhat, zhat)

        ihat = xform[0, :]
        yhat = xform[1, :]
        zhat = xform[2, :]

        is_failed = False
        return is_failed, (wa, wb, ihat, yhat, zhat)

    def get_orientation_vector(self, xyz):
        """
        Element offsets are defined in a Cartesian system located at the
        connecting grid point. The components of the offsets are always
        defined in units of translation, even if the displacement
        coordinate system is cylindrical or spherical.

        For example, in Figure 11-11, the grid point displacement
        coordinate system is cylindrical, and the offset vector is
        defined using Cartesian coordinates u1, u2, and u3 in units of
        translation.
        """
        if self.g0:
            v = xyz[self.g0] - xyz[self.Ga()]
        else:
            v = self.x
        assert self.offt == 'GGG', self.offt
        return v

    @property
    def node_ids(self):
        return [self.Ga(), self.Gb()]

    def get_edge_ids(self):
        return [tuple(sorted(self.node_ids))]

    @property
    def nodes(self):
        return [self.ga, self.gb]

    @nodes.setter
    def nodes(self, values):
        self.ga = values[0]
        self.gb = values[1]

    @property
    def nodes_ref(self):
        return [self.ga_ref, self.gb_ref]

    @nodes_ref.setter
    def nodes_ref(self, values):
        assert values is not None, values
        self.ga_ref = values[0]
        self.gb_ref = values[1]

    def raw_fields(self):
        """Gets the fields of the card in their full form"""
        (x1, x2, x3) = self.get_x_g0_defaults()

        # offt doesn't exist in NX nastran
        offt = set_blank_if_default(self.offt, 'GGG')

        list_fields = ['CBAR', self.eid, self.Pid(), self.Ga(), self.Gb(), x1, x2,
                       x3, offt, self.pa, self.pb] + list(self.wa) + list(self.wb)
        return list_fields

    def repr_fields(self):
        """Gets the fields of the card in their reduced form"""
        pa = set_blank_if_default(self.pa, 0)
        pb = set_blank_if_default(self.pb, 0)

        w1a = set_blank_if_default(self.wa[0], 0.0)
        w2a = set_blank_if_default(self.wa[1], 0.0)
        w3a = set_blank_if_default(self.wa[2], 0.0)

        w1b = set_blank_if_default(self.wb[0], 0.0)
        w2b = set_blank_if_default(self.wb[1], 0.0)
        w3b = set_blank_if_default(self.wb[2], 0.0)
        x1, x2, x3 = self.get_x_g0_defaults()

        # offt doesn't exist in NX nastran
        offt = set_blank_if_default(self.offt, 'GGG')

        list_fields = ['CBAR', self.eid, self.Pid(), self.Ga(), self.Gb(), x1, x2,
                       x3, offt, pa, pb, w1a, w2a, w3a, w1b, w2b, w3b]
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)

    def write_card_16(self, is_double=False):
        card = self.repr_fields()
        return self.comment + print_card_16(card)


class CBEAM3(LineElement):  # was CBAR
    _properties = ['node_ids', 'nodes']
    """
    Defines a three-node beam element

    """
    @classmethod
    def _init_from_empty(cls):
        eid = 1
        pid = 1
        x = None
        nids = [1, 2, 3]
        g0 = 4
        return CBEAM3(eid, pid, nids, x, g0,
                      wa=None, wb=None, wc=None,
                      tw=None, s=None, comment='')

    type = 'CBEAM3'
    def __init__(self, eid, pid, nids, x, g0,
                 wa=None, wb=None, wc=None, tw=None, s=None, comment=''):
        LineElement.__init__(self)

        if wa is None:
            wa = np.zeros(3, dtype='float64')
        if wb is None:
            wb = np.zeros(3, dtype='float64')
        if wc is None:
            wc = np.zeros(3, dtype='float64')
        if tw is None:
            tw = np.zeros(3, dtype='float64')
        if s is None:
            s = np.zeros(3, dtype='int32')

        if comment:
            self.comment = comment
        self.eid = eid
        self.pid = pid
        self.ga = nids[0]
        self.gb = nids[1]
        self.gc = nids[2]
        self.x = x
        self.g0 = g0
        self.wa = wa
        self.wb = wb
        self.wc = wc
        self.tw = tw
        self.s = s
        self.ga_ref = None
        self.gb_ref = None
        self.gc_ref = None
        self.g0_ref = None
        self.pid_ref = None

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a CBEAM3 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        eid = integer(card, 1, 'eid')
        pid = integer_or_blank(card, 2, 'pid', eid)
        ga = integer(card, 3, 'ga')
        gb = integer(card, 4, 'gb')
        gc = integer_or_blank(card, 5, 'gc')

        # card, eid, x1_default, x2_default, x3_default
        x, g0 = init_x_g0_cbeam3(card, eid, 0., 0., 0.)
        wa = np.array([double_or_blank(card, 9, 'w1a', 0.0),
                       double_or_blank(card, 10, 'w2a', 0.0),
                       double_or_blank(card, 11, 'w3a', 0.0)], dtype='float64')

        wb = np.array([double_or_blank(card, 12, 'w1b', 0.0),
                       double_or_blank(card, 13, 'w2b', 0.0),
                       double_or_blank(card, 14, 'w3b', 0.0)], dtype='float64')

        wc = np.array([double_or_blank(card, 15, 'w1c', 0.0),
                       double_or_blank(card, 16, 'w2c', 0.0),
                       double_or_blank(card, 17, 'w3c', 0.0)], dtype='float64')

        tw = np.array([double_or_blank(card, 18, 'twa', 0.),
                       double_or_blank(card, 19, 'twb', 0.),
                       double_or_blank(card, 20, 'twc', 0.)], dtype='float64')

        # TODO: what are the defaults?
        s = np.array([integer_or_blank(card, 21, 'sa', -1),
                      integer_or_blank(card, 22, 'sb', -1),
                      integer_or_blank(card, 23, 'sc', -1)], dtype='int32')
        assert len(card) <= 24, 'len(CBEAM3 card) = %i\ncard=%s' % (len(card), card)
        return CBEAM3(eid, pid, [ga, gb, gc], x, g0,
                      wa=wa, wb=wb, wc=wc, tw=tw, s=s, comment=comment)

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        msg = ', which is required by CBEAM3 eid=%s' % (self.eid)
        self.ga_ref = model.Node(self.ga, msg=msg)
        self.gb_ref = model.Node(self.gb, msg=msg)
        if self.gc is not None:
            self.gc_ref = model.Node(self.gc, msg=msg)
        self.pid_ref = model.Property(self.pid, msg=msg)
        if self.g0:
            self.g0_ref = model.Node(self.g0, msg=msg)

    def safe_cross_reference(self, model, xref_errors):
        msg = ', which is required by CBEAM3 eid=%s' % (self.eid)
        self.ga_ref = model.Node(self.ga, msg=msg)
        self.gb_ref = model.Node(self.gb, msg=msg)
        if self.gc is not None:
            self.gc_ref = model.Node(self.gc, msg=msg)
        self.pid_ref = model.safe_property(self.pid, self.eid, xref_errors, msg=msg)
        if self.g0:
            self.g0_ref = model.Node(self.g0, msg=msg)

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.ga = self.Ga()
        self.gb = self.Gb()
        self.gc = self.Gc()
        if self.g0_ref is not None:
            self.g0 = self.G0()
            self.g0_ref = None
        self.pid = self.Pid()
        self.ga_ref = None
        self.gb_ref = None
        self.gc_ref = None
        self.pid_ref = None

    def center_of_mass(self):
        return np.zeros(3)
    def Centroid(self):
        return np.zeros(3)
    def MassPerLength(self):
        return 0.

    def Length(self):
        r"""
        We'll fit a 2nd order polynomial to the x, y, and z coefficients.
        We know GA (t=0), GB(t=1), and GC (t=0.5 assumed).
        This gives us A, for:
         - y1 = a1*t^2 = b1*t + c1

        where:
         - yi is for the x, y, and z terms

        [xa, xb, xc] = [A][a1, b1, c1].T
        [ya, yb, yc] = [A][a2, b2, c2].T
        [za, zb, zc] = [A][a3, b3, c3].T

        [a1, b1, c1] = [A^-1][[xa, xb, xc].T
        A = [0.  , 0.  , 1.  ]
            [1.  , 1.  , 1.  ]
            [0.25, 0.5 , 1.  ]
        Ainv = [ 2.,  2., -4.]
               [-3., -1.,  4.]
               [ 1.,  0.,  0.]

        """
        xyza = self.ga_ref.get_position() + self.wa
        xyzb = self.gb_ref.get_position() + self.wb
        if self.gc is not None:
            xyzc = self.gc_ref.get_position() + self.wc
            xa, ya, za = xyza
            length = self._integrate(
                xyza - xa,
                xyzb - ya,
                xyzc - za,
            )
        else:
            length = np.linalg.norm(xyzb - xyza)
        return length

    def _integrate(self, xabc, yabc, zabc):
        """
        We integrate:
         y = sqrt(x'(t)^2 + y'(t)^2 + z'(t)^2)*dt from 0 to 1
         y = sqrt(r)

        /where:
         - x(t) = a*t^2 + b*t + c
         - x'(t) = 2*a*t + b
         - x'(t)^2 = 4*(a*t)^2 + 2*a*b*t + b^2

        expanding terms:
         - x'(t)^2 = 4*(a1*t)^2 + 2*a1*b1*t + b1^2
         - y'(t)^2 = 4*(a2*t)^2 + 2*a2*b2*t + b2^2
         - z'(t)^2 = 4*(a3*t)^2 + 2*a3*b3*t + b3^2

        grouping terms:
         - a = 4 * (a1 ** 2 + a2 ** 2 + a3 ** 2) * t^2
         - b = 2 * (a1 * b1 + a2 * b2 + a3 * b3) * t
         - c = b1 ** 2 + b2 ** 2 + b3 ** 2
         - y = integrate(sqrt(a*t^2 + b*t + t), t, 0., 1.)
        Looking up integral formulas, we get a really complicated integral.

        for a = 0 (and rewriting):
         - y = integrate(sqrt(a*t + b), t, 0., 1.)
         - y = (2*b/3a + 2*t/3) * sqrt(at + b)
        or:
         - y = (2*c/3b + 2*t/3) * sqrt(b*t + c)
        """
        #print(xabc, yabc, zabc)
        Ainv = np.array([
            [2., 2., -4.],
            [-3., -1., 4.],
            [1., 0., 0.],
        ])
        a1, b1, unused_c1 = Ainv @ xabc.reshape(3, 1)
        a2, b2, unused_c2 = Ainv @ yabc.reshape(3, 1)
        a3, b3, unused_c3 = Ainv @ zabc.reshape(3, 1)
        #print(x, y, z)
        #print('---------------------')
        #print(a1, a2, a3)
        a = 4 * (a1 ** 2 + a2 ** 2 + a3 ** 2)[0]
        b = 2 * (a1 * b1 + a2 * b2 + a3 * b3)[0]
        c = (b1 ** 2 + b2 ** 2 + b3 ** 2)[0]
        #print(a, b, c)
        if np.allclose(a, 0) and np.allclose(b, 0.):
            length = c ** 0.5
        elif np.allclose(a, 0):
            #print(self.ga_ref)
            #print(self.gb_ref)
            #print(self.gc_ref)
            # dx = self.gb_ref.get_position() - self.ga_ref.get_position()
            #coeff = ((2 * b) / (3 * a) + (2 * t) / 3)
            #length = coeff * np.sqrt(a*t + b)
            dx = [xabc[1], yabc[1], zabc[1]]
            length = np.linalg.norm(dx)
            #print('dx=', dx, length)
        else:
            t = 1.
            at2btc = a * t **2 + b * t + c
            length = (
                (b + 2 * a * t) / (4 * a) * np.sqrt(at2btc) +
                (4 * a * c - b ** 2) / (8 * a ** 1.5) * np.log(
                    2*a*t + b + 2*np.sqrt(a*(at2btc)))
            )
        return length

    def Area(self):
        if isinstance(self.pid_ref, integer_types):
            msg = 'Element eid=%i has not been cross referenced.\n%s' % (self.eid, str(self))
            raise RuntimeError(msg)
        A = self.pid_ref.Area()

        xyza = self.ga_ref.get_position() + self.wa
        xyzb = self.gb_ref.get_position() + self.wb
        Aa, Ab, Ac = A # self.pid_ref.area
        L = self.Length()
        if self.gc is not None:
            xa, ya, za = xyza
            xb, yb, zb = xyzb
            xc, yc, zc = self.gc_ref.get_position() + self.wc
            area = self._integrate(
                np.array([Aa * xa, Ab * xb, Ac * xc]) / L,
                np.array([Aa * ya, Ab * yb, Ac * yc]) / L,
                np.array([Aa * za, Ab * zb, Ac * zc]) / L,
            )
        else:
            area = np.linalg.norm(xyzb*Ab - xyza*Aa) / L
        assert isinstance(area, float), area
        return area

    def Volume(self):
        if isinstance(self.pid_ref, integer_types):
            msg = 'Element eid=%i has not been cross referenced.\n%s' % (self.eid, str(self))
            raise RuntimeError(msg)
        A = self.pid_ref.Area()

        xyza = self.ga_ref.get_position() + self.wa
        xyzb = self.gb_ref.get_position() + self.wb
        Aa, Ab, Ac = A # self.pid_ref.area
        if self.gc is not None:
            xa, ya, za = xyza
            xb, yb, zb = xyzb
            xc, yc, zc = self.gc_ref.get_position() + self.wc
            volume = self._integrate(
                np.array([Aa * xa, Ab * xb, Ac * xc]),
                np.array([Aa * ya, Ab * yb, Ac * yc]),
                np.array([Aa * za, Ab * zb, Ac * zc]),
            )
        else:
            volume = np.linalg.norm(xyzb*Ab - xyza*Aa)
        assert isinstance(volume, float), volume
        return volume

    def Nsm(self):
        if isinstance(self.pid_ref, integer_types):
            msg = 'Element eid=%i has not been cross referenced.\n%s' % (self.eid, str(self))
            raise RuntimeError(msg)
        nsm = self.pid_ref.Nsm()
        assert isinstance(nsm, float), nsm

        xyza = self.ga_ref.get_position() + self.wa
        xyzb = self.gb_ref.get_position() + self.wb
        if self.gc is not None:
            xa, ya, za = xyza
            xb, yb, zb = xyzb
            xc, yc, zc = self.gc_ref.get_position() + self.wc
            nsma, nsmb, nsmc = self.pid_ref.nsm
            L = self.Length()
            nsm = self._integrate(
                np.array([nsma * xa, nsmb * xb, nsmc * xc]) / L,
                np.array([nsma * ya, nsmb * yb, nsmc * yc]) / L,
                np.array([nsma * za, nsmb * zb, nsmc * zc]) / L,
            )
        else:
            nsm = np.linalg.norm(xyzb*nsmb - xyza*nsmb) / L
        assert isinstance(nsm, float), nsm
        return nsm

    def Ga(self):
        """gets node 1"""
        if self.ga_ref is None:
            return self.ga
        return self.ga_ref.nid

    def Gb(self):
        """gets node 2"""
        if self.gb_ref is None:
            return self.gb
        return self.gb_ref.nid

    def Gc(self):
        """gets the node between node 1 and 2"""
        if self.gc_ref is None:
            return self.gc
        return self.gc_ref.nid

    def G0(self):
        """gets the orientation vector node"""
        if self.g0_ref is None:
            return self.g0
        return self.g0_ref.nid

    @property
    def nodes(self):
        return [self.ga, self.gb, self.gc]

    @property
    def node_ids(self):
        return [self.Ga(), self.Gb(), self.Gc()]


    def get_edge_ids(self):
        nids = [self.Ga(), self.Gb()]
        return [tuple(sorted(nids))]

    def get_x_g0_defaults(self):
        """
        X and G0 compete for the same fields, so the method exists to
        make it easier to write the card

        Returns
        -------
        x_g0 : varies
            g0 : List[int, None, None]
            x : List[float, float, float]

        Notes
        -----
        Used by CBAR, CBEAM, and CBEAM3

        """
        if self.g0 is not None:
            return self.G0(), None, None
        return list(self.x)

    def raw_fields(self):
        x1, x2, x3 = self.get_x_g0_defaults()
        ga, gb, gc = self.node_ids
        list_fields = ['CBEAM3', self.eid, self.Pid(), ga, gb, gc, x1, x2, x3] + \
                  list(self.wa) + list(self.wb) + list(self.wc) + list(self.tw) + list(self.s)
        return list_fields

    def repr_fields(self):
        w1a = set_blank_if_default(self.wa[0], 0.0)
        w2a = set_blank_if_default(self.wa[1], 0.0)
        w3a = set_blank_if_default(self.wa[2], 0.0)
        w1b = set_blank_if_default(self.wb[0], 0.0)
        w2b = set_blank_if_default(self.wb[1], 0.0)
        w3b = set_blank_if_default(self.wb[2], 0.0)
        w1c = set_blank_if_default(self.wc[0], 0.0)
        w2c = set_blank_if_default(self.wc[1], 0.0)
        w3c = set_blank_if_default(self.wc[2], 0.0)

        twa = set_blank_if_default(self.tw[0], 0.0)
        twb = set_blank_if_default(self.tw[1], 0.0)
        twc = set_blank_if_default(self.tw[2], 0.0)

        x1, x2, x3 = self.get_x_g0_defaults()
        ga, gb, gc = self.node_ids
        list_fields = ['CBEAM3', self.eid, self.Pid(), ga, gb, gc, x1, x2, x3,
                       w1a, w2a, w3a, w1b, w2b, w3b, w1c, w2c, w3c,
                       twa, twb, twc, self.s[0], self.s[1], self.s[2]]
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)

    def _verify(self, xref):
        unused_edges = self.get_edge_ids()

def init_x_g0_cbeam3(card, eid, x1_default, x2_default, x3_default):
    """reads the x/g0 field for the CBEAM3"""
    field6 = integer_double_or_blank(card, 6, 'g0_x1', x1_default)

    if isinstance(field6, integer_types):
        g0 = field6
        x = None
    elif isinstance(field6, float):
        g0 = None
        x = np.array([field6,
                      double_or_blank(card, 7, 'x2', x2_default),
                      double_or_blank(card, 8, 'x3', x3_default)], dtype='float64')
        if norm(x) == 0.0:
            msg = 'G0 vector defining plane 1 is not defined.\n'
            msg += 'G0 = %s\n' % g0
            msg += 'X  = %s\n' % x
            raise RuntimeError(msg)
    else:
        msg = ('field5 on %s (G0/X1) is the wrong type...id=%s field5=%s '
               'type=%s' % (card.field(0), eid, field6, type(field6)))
        raise RuntimeError(msg)
    return x, g0


class CBEND(LineElement):
    type = 'CBEND'
    _field_map = {
        1: 'eid', 2:'pid', 3:'ga', 4:'gb', 8:'geom',
    }
    _properties = ['node_ids']

    def _update_field_helper(self, n, value):
        if self.g0 is not None:
            if n == 5:
                self.g0 = value
            else:
                raise KeyError('Field %r=%r is an invalid %s entry.' % (n, value, self.type))
        else:
            if n == 5:
                self.x[0] = value
            elif n == 6:
                self.x[1] = value
            elif n == 7:
                self.x[2] = value
            else:
                raise KeyError('Field %r=%r is an invalid %s entry.' % (n, value, self.type))

    @classmethod
    def _init_from_empty(cls):
        eid = 1
        pid = 1
        nids = [1, 2]
        g0 = 4
        x = None
        geom = 1
        return CBEND(eid, pid, nids, g0, x, geom, comment='')

    def __init__(self, eid, pid, nids, g0, x, geom, comment=''):
        """
        Creates a CBEND card

        Parameters
        ----------
        eid : int
            element id
        pid : int
            property id (PBEND)
        nids : List[int, int]
            node ids; connected grid points at ends A and B
        g0 : int
            ???
        x : List[float, float, float]
            ???
        geom : int
            1 : The center of curvature lies on the line AO (or its extension) or vector v.
            2 : The tangent of centroid arc at end A is parallel to line AO or vector v.
                Point O (or vector v) and the arc must be on the same side of the chord AB.
            3 : The bend radius (RB) is specified on the PBEND entry:
                Points A, B, and O (or vector v) define a plane parallel or coincident
                with the plane of the element arc. Point O (or vector v) lies on the
                opposite side of line AB from the center of the curvature.
            4 : THETAB is specified on the PBEND entry. Points A, B, and O (or vector v)
                define a plane parallel or coincident with the plane of the element arc.
                Point O (or vector v) lies on the opposite side of line AB from the center
                of curvature.
        comment : str; default=''
            a comment for the card

        """
        LineElement.__init__(self)
        if comment:
            self.comment = comment
        self.eid = eid
        self.pid = pid
        self.ga = nids[0]
        self.gb = nids[1]

        if g0 is None:
            assert x is not None, 'g0=%s x=%s; one must not be None' % (g0, x)
        self.g0 = g0
        self.x = x
        self.geom = geom
        assert self.geom in [1, 2, 3, 4], 'geom is invalid geom=%r' % self.geom
        self.ga_ref = None
        self.gb_ref = None
        self.pid_ref = None
        if self.g0 is not None:
            assert isinstance(self.g0, integer_types), self.get_stats()

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a CBEND card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        eid = integer(card, 1, 'eid')
        pid = integer_or_blank(card, 2, 'pid', eid)
        ga = integer(card, 3, 'ga')
        gb = integer(card, 4, 'gb')
        x1_g0 = integer_double_or_blank(card, 5, 'x1_g0', 0.0)
        if isinstance(x1_g0, integer_types):
            g0 = x1_g0
            x = None
        elif isinstance(x1_g0, float):
            g0 = None
            x = np.array([double_or_blank(card, 5, 'x1', 0.0),
                          double_or_blank(card, 6, 'x2', 0.0),
                          double_or_blank(card, 7, 'x3', 0.0)], dtype='float64')
            if norm(x) == 0.0:
                msg = 'G0 vector defining plane 1 is not defined.\n'
                msg += 'G0 = %s\n' % g0
                msg += 'X  = %s\n' % x
                raise RuntimeError(msg)
        else:
            raise ValueError('invalid x1Go=%r on CBEND' % x1_g0)
        geom = integer(card, 8, 'geom')

        assert len(card) == 9, 'len(CBEND card) = %i\ncard=%s' % (len(card), card)
        return CBEND(eid, pid, [ga, gb], g0, x, geom, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        #data = [[eid, pid, ga, gb, geom], [f, x1, x2, x3]]
        #data = [[eid, pid, ga, gb, geom], [f, g0]]

        main = data[0]
        flag = data[1][0]
        if flag in [0, 1]:
            g0 = None
            x = np.array([data[1][1],
                          data[1][2],
                          data[1][3]], dtype='float64')
        else:
            g0 = data[1][1]
            x = None

        eid = main[0]
        pid = main[1]
        ga = main[2]
        gb = main[3]
        geom = main[4]
        return CBEND(eid, pid, [ga, gb], g0, x, geom, comment=comment)

    def get_x_g0_defaults(self):
        if self.g0 is not None:
            return self.g0, None, None
        #print('x =', self.x)
        #print('g0 =', self.g0)
        #x1 = set_blank_if_default(self.x[0], 0.0)
        #x2 = set_blank_if_default(self.x[1], 0.0)
        #x3 = set_blank_if_default(self.x[2], 0.0)
        return list(self.x)

    def Length(self):
        # TODO: consider w1a and w1b in the length formulation
        L = norm(self.gb_ref.get_position() - self.ga_ref.get_position())
        assert isinstance(L, float)
        return L

        #prop = self.pid_ref
        #bend_radius = prop.rb
        #theta_bend = prop.thetab
        #length_oa = None
        #if self.geom == 1:
            #The center of curvature lies on the line AO
            #(or its extension) or vector .
            #pass
        #elif self.geom == 2:
            # The tangent of centroid arc at end A is
            # parallel to line AO or vector . Point O (or
            # vector) and the arc must be on the
            # same side of the chord .
            #pass
        #elif self.geom == 3:
            # The bend radius (RB) is specified on the
            # PBEND entry: Points A, B, and O (or
            # vector ) define a plane parallel or
            # coincident with the plane of the element
            # arc. Point O (or vector ) lies on the
            # opposite side of line AB from the center of
            # the curvature.
            #pass
        #elif self.geom == 4:
            # THETAB is specified on the PBEND entry.
            # Points A, B, and O (or vector ) define a
            # plane parallel or coincident with the plane
            # of the element arc. Point O (or vector )
            # lies on the opposite side of line AB from the
            # center of curvature.
            #pass
        #else:
            #raise RuntimeError('geom=%r is not supported on the CBEND' % self.geom)
        #return L

    def validate(self):
        if self.g0 is not None:
            assert isinstance(self.g0, integer_types), 'g0=%s must be an integer' % self.g0
        if self.g0 in [self.ga, self.gb]:
            msg = 'G0=%s cannot be GA=%s or GB=%s' % (self.g0, self.ga, self.gb)
            raise RuntimeError(msg)
        #BEND ELEMENT %1 BEND RADIUS OR ARC ANGLE INCONSISTENT
        #WITH GEOM OPTION
        #RB is nonzero on PBEND entry when GEOM option on CBEND entry is 1,
        #2, or 4 or RB is zero when GEOM option is 3 or AB is nonzero when
        #when GEOM option is 1, 2, or 3 or B is <= 0. or > 180, when
        #GEOM option is 4.

    @property
    def node_ids(self):
        return [self.Ga(), self.Gb()]

    @property
    def nodes(self):
        return [self.ga, self.gb]

    @nodes.setter
    def nodes(self, values):
        self.ga = values[0]
        self.gb = values[1]

    def Ga(self):
        if self.ga_ref is None:
            return self.ga
        return self.ga_ref.nid

    def Gb(self):
        if self.gb_ref is None:
            return self.gb
        return self.gb_ref.nid

    #def get_edge_ids(self):
        #return [tuple(sorted(self.node_ids))]

    @property
    def nodes_ref(self):
        return [self.ga_ref, self.gb_ref]

    @nodes_ref.setter
    def nodes_ref(self, values):
        assert values is not None, values
        self.ga_ref = values[0]
        self.gb_ref = values[1]

    def Area(self):
        if isinstance(self.pid, integer_types):
            msg = 'Element eid=%i has not been cross referenced.\n%s' % (self.eid, str(self))
            raise RuntimeError(msg)
        return self.pid_ref.Area()

    def _verify(self, xref):
        unused_edges = self.get_edge_ids()

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        msg = ', which is required by CBEND eid=%s' % (self.eid)
        #self.g0 = model.nodes[self.g0]
        self.ga_ref = model.Node(self.ga, msg=msg)
        self.gb_ref = model.Node(self.gb, msg=msg)
        self.pid_ref = model.Property(self.pid, msg=msg)

    def safe_cross_reference(self, model, xref_errors):
        msg = ', which is required by CBEND eid=%s' % (self.eid)
        self.ga_ref = model.Node(self.ga, msg=msg)
        self.gb_ref = model.Node(self.gb, msg=msg)
        self.pid_ref = model.safe_property(self.pid, self.eid, xref_errors, msg=msg)

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        node_ids = self.node_ids
        self.ga = node_ids[0]
        self.gb = node_ids[1]
        self.pid = self.Pid()
        self.ga_ref = None
        self.gb_ref = None
        self.pid_ref = None

    def Centroid(self):
        if self.pid_ref is None:
            msg = 'Element eid=%i has not been cross referenced.\n%s' % (self.eid, str(self))
            raise RuntimeError(msg)
        return (self.ga_ref.get_position() + self.gb_ref.get_position()) / 2.

    def center_of_mass(self):
        return self.Centroid()

    def raw_fields(self):
        (x1, x2, x3) = self.get_x_g0_defaults()
        list_fields = ['CBEND', self.eid, self.Pid(), self.Ga(), self.Gb(),
                       x1, x2, x3, self.geom]
        return list_fields

    def repr_fields(self):
        return self.raw_fields()

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)

def init_x_g0(card, eid, x1_default, x2_default, x3_default):
    """common method to read the x/g0 field for the CBAR, CBEAM, CBEAM3"""
    field5 = integer_double_or_blank(card, 5, 'g0_x1', x1_default)

    if isinstance(field5, integer_types):
        g0 = field5
        x = None
    elif isinstance(field5, float):
        g0 = None
        x = np.array([field5,
                      double_or_blank(card, 6, 'x2', x2_default),
                      double_or_blank(card, 7, 'x3', x3_default)], dtype='float64')
        if norm(x) == 0.0:
            msg = 'G0 vector defining plane 1 is not defined.\n'
            msg += 'G0 = %s\n' % g0
            msg += 'X  = %s\n' % x
            raise RuntimeError(msg)
    else:
        msg = ('field5 on %s (G0/X1) is the wrong type...id=%s field5=%s '
               'type=%s' % (card.field(0), eid, field5, type(field5)))
        raise RuntimeError(msg)
    return x, g0

def get_bar_vector(model, elem, node1, node2, xyz1):
    """helper method for ``rotate_v_wa_wb``"""
    cd1 = node1.Cd()
    cd2 = node2.Cd()
    if model is None:
        cd1_ref = node1.cd_ref
        cd2_ref = node2.cd_ref

        # get the vector v, which defines the projection on to the elemental
        # coordinate frame
        if elem.g0:
            #msg = 'which is required by %s eid=%s\n%s' % (elem.type, elem.g0, str(elem))
            g0_ref = elem.g0_ref
            n0 = g0_ref.get_position()
            v = n0 - xyz1
        else:
            v = cd1_ref.transform_node_to_global(elem.x)

    else:
        msg = ', which is required by %s=%s' % (elem.type, elem.eid)
        cd1_ref = model.Coord(cd1)
        cd2_ref = model.Coord(cd2)

        # get the vector v, which defines the projection on to the elemental
        # coordinate frame
        if elem.g0:
            #msg = 'which is required by %s eid=%s\n%s' % (elem.type, elem.g0, str(elem))
            g0_ref = model.Node(elem.g0, msg=msg)
            n0 = g0_ref.get_position()
            v = n0 - xyz1
        else:
            v = cd1_ref.transform_node_to_global(elem.x)
            cd1_ref = model.Coord(cd1)
            cd2_ref = model.Coord(cd2)

    return v, cd1, cd1_ref, cd2, cd2_ref

def rotate_v_wa_wb(model: BDF, elem, xyz1, xyz2, node1, node2, ihat_offset, i_offset, eid,
                   Li_offset,
                   log: SimpleLogger) -> Tuple[NDArray3float, NDArray3float, NDArray3float, NDArray33float]:
    """
    Rotates v, wa, wb

    Parameters
    ----------
    model : BDF()
        BDF : assume the model isn't xref'd
        None : use the xref'd values
    elem : CBAR() / CBEAM()
       the CBAR/CBEAM
    xyz1 / xyz2 : (3, ) float ndarray
        the xyz locations for node 1 / 2
    node1 / node2 : GRID()
        the xyz object for node 1 / 2
    ihat_offset : (3, ) float ndarray
        the normalized x-axis (not including the CBEAM offset)
    i_offset : (3, ) float ndarray
        the unnormalized x-axis (not including the CBEAM offset)
    eid : int
        the element id
    Li_offset : float
        the length of the CBAR/CBEAM (not including the CBEAM offset)
    log : Log()
        a logging object or None

    Returns
    -------
    v : List[float, float, float]
        the projection vector that defines the y-axis (jhat)
    wa : List[float, float, float]
       the offset vector at A
    wb : List[float, float, float]
       the offset vector at B
    xform : (3, 3) float ndarray
        a vstack of the [ihat, jhat, khat] axes

    Notes
    -----
    This section details the OFFT flag.

    ABC or A-B-C (an example is G-G-G or B-G-G)
    while the slots are:
     - A -> orientation; values=[G, B]
     - B -> End A; values=[G, O]
     - C -> End B; values=[G, O]

    and the values for A,B,C mean:
     - B -> basic
     - G -> global
     - O -> orientation

    so for example G-G-G, that's global for all terms.
    BOG means basic orientation, orientation end A, global end B

    so now we're left with what does basic/global/orientation mean?
    - basic -> the global coordinate system defined by cid=0
    - global -> the local coordinate system defined by the
                CD field on the GRID card, but referenced by
                the CBAR/CBEAM
    - orientation -> wa/wb are defined in the xform_offset (yz) frame;
                     this is likely the easiest frame for a user

    """
    check_offt(elem)
    v, cd1, cd1_ref, cd2, cd2_ref = get_bar_vector(model, elem, node1, node2, xyz1)

    #--------------------------------------------------------------------------
    offt_vector, offt_end_a, offt_end_b = elem.offt

    # rotate v
    if offt_vector == 'G':
        # end A
        # global - cid != 0
        if cd1 != 0:
            v = cd1_ref.transform_node_to_global_assuming_rectangular(v)
    elif offt_vector == 'B':
        # basic - cid = 0
        pass
    else:
        msg = 'offt_vector=%r is not supported; offt=%s' % (offt_vector, elem.offt)
        return None, None, None, None

    yhat_offset, zhat_offset = get_bar_yz_transform(
        v, ihat_offset, eid, xyz1, xyz2, node1.nid, node2.nid,
        i_offset, Li_offset)
    xform_offset = np.vstack([ihat_offset, yhat_offset, zhat_offset]) # 3x3 unit matrix

    #--------------------------------------------------------------------------
    # rotate wa
    # wa defines the offset at end A
    wa = elem.wa
    #ia = n1
    if offt_end_a == 'G':
        if cd1 != 0:
            wa = cd1_ref.transform_node_to_global_assuming_rectangular(wa)
    elif offt_end_a == 'B':
        pass
    elif offt_end_a == 'O':
        # rotate point wa from the local frame to the global frame
        wa = wa @ xform_offset
        #ia = n1 + wa
    else:
        msg = 'offt_end_a=%r is not supported; offt=%s' % (offt_end_a, elem.offt)
        log.error(msg)
        return v, None, None, xform_offset

    #--------------------------------------------------------------------------
    # rotate wb
    # wb defines the offset at end B
    wb = elem.wb
    #ib = n2
    if offt_end_b == 'G':
        if cd2 != 0:
            # MasterModelTaxi
            wb = cd2_ref.transform_node_to_global_assuming_rectangular(wb)
    elif offt_end_b == 'B':
        pass
    elif offt_end_b == 'O':
        # rotate point wb from the local frame to the global frame
        wb = wb @ xform_offset
        #ib = n2 + wb
    else:
        msg = 'offt_end_b=%r is not supported; offt=%s' % (offt_end_b, elem.offt)
        model.log.error(msg)
        return v, wa, None, xform_offset

    #--------------------------------------------------------------------------
    #i = ib - ia # (xyz2 + wb) - (xyz1 + wa)
    #i = (xyz2 + wb) - (xyz1 + wa)
    i = i_offset
    Li = norm(i)
    ihat = i / Li
    yhat, zhat = get_bar_yz_transform(v, ihat, eid, xyz1, xyz2, node1.nid, node2.nid, i, Li)

    #print('  n1=%s n2=%s' % (n1, n2))
    #print('  ib=%s ia=%s' % (ib, ia))
    #print('  wa=%s wb=%s' % (wa, wb))
    #print('  ioffset=%s i=%s' % (i_offset, i))
    #print('  ihat=%s' % (ihat))
    #print('  yhat=%s' % (yhat))
    #print('  zhat=%s' % (zhat))
    #print("")

    xform = np.vstack([ihat, yhat, zhat]) # 3x3 unit matrix

    return v, wa, wb, xform

def get_bar_yz_transform(v, ihat, eid, xyz1, xyz2, nid1, nid2, i, Li):
    """
    helper method for ``_get_bar_yz_arrays``

    Parameters
    ----------
    v : List[float, float, float]
        the projection vector that defines the y-axis (jhat)
    ihat : (3, ) float ndarray
        the normalized x-axis (not including the CBEAM offset)
    eid : int
        the element id
    xyz1 / xyz2 : (3, ) float ndarray
        the xyz locations for node 1 / 2
    nid1 / nid2  : int
        the node ids for xyz1 / xyz2
    i : (3, ) float ndarray
        the unnormalized x-axis (not including the CBEAM offset)
    Li : float
        the length of the CBAR/CBEAM (not including the CBEAM offset)

    Returns
    -------
    yhat (3, ) float ndarray
       the CBAR/CBEAM's y-axis
    zhat (3, ) float ndarray
       the CBAR/CBEAM's z-axis

    """
    vhat = v / norm(v) # j
    try:
        z = np.cross(ihat, vhat) # k
    except ValueError:
        msg = 'Invalid vector length\n'
        msg += 'xyz1=%s\n' % str(xyz1)
        msg += 'xyz2=%s\n' % str(xyz2)
        msg += 'nid1=%s\n' % str(nid1)
        msg += 'nid2=%s\n' % str(nid2)
        msg += 'i   =%s\n' % str(i)
        msg += 'Li  =%s\n' % str(Li)
        msg += 'ihat=%s\n' % str(ihat)
        msg += 'v   =%s\n' % str(v)
        msg += 'vhat=%s\n' % str(vhat)
        msg += 'z=cross(ihat, vhat)'
        print(msg)
        raise ValueError(msg)

    zhat = z / norm(z)
    yhat = np.cross(zhat, ihat) # j

    if norm(ihat) == 0.0 or norm(yhat) == 0.0 or norm(z) == 0.0:
        print('  invalid_orientation - eid=%s yhat=%s zhat=%s v=%s i=%s n%s=%s n%s=%s' % (
            eid, yhat, zhat, v, i, nid1, xyz1, nid2, xyz2))
    elif not np.allclose(norm(yhat), 1.0) or not np.allclose(norm(zhat), 1.0) or Li == 0.0:
        print('  length_error        - eid=%s Li=%s Lyhat=%s Lzhat=%s'
              ' v=%s i=%s n%s=%s n%s=%s' % (
                  eid, Li, norm(yhat), norm(zhat), v, i, nid1, xyz1, nid2, xyz2))
    return yhat, zhat

def check_offt(element):
    """
    B,G,O
    Note: The character 'O' in the table replaces the obsolete character 'E'

    """
    if isinstance(element.offt, integer_types):
        raise SyntaxError('invalid offt expected a string of length 3; '
                          'offt=%r; Type=%s\n%s' % (element.offt, type(element.offt), str(element)))
    msg = 'invalid offt parameter of %s...offt=%s' % (element.type, element.offt)
    assert element.offt[0] in ['G', 'B'], msg
    assert element.offt[1] in ['G', 'O', 'E'], msg
    assert element.offt[2] in ['G', 'O', 'E'], msg
