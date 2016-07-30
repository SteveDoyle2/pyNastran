# pylint: disable=R0904,R0902,E1101,E1103,C0111,C0302,C0103,W0101
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
from six.moves import range

from numpy.linalg import norm

from pyNastran.utils import integer_types
from pyNastran.bdf.field_writer_8 import set_blank_if_default
from pyNastran.bdf.cards.base_card import Element #, Mid
from pyNastran.bdf.bdf_interface.assign_type import (
    integer, integer_or_blank, double, double_or_blank)
from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.field_writer_16 import print_card_16


class RodElement(Element):  # CROD, CONROD, CTUBE

    def __init__(self):
        Element.__init__(self)

    def cross_reference(self, model):
        msg = ' which is required by %s eid=%s' % (self.type, self.eid)
        self.nodes = model.Nodes(self.nodes, msg=msg)
        self.pid = model.Property(self.pid, msg=msg)
        self.nodes_ref = self.nodes
        self.pid_ref = self.pid

    def uncross_reference(self):
        self.nodes = self.node_ids
        self.pid = self.Pid()
        del self.nodes_ref, self.pid_ref

    def nodeIDs(self):
        self.deprecated('self.nodeIDs()', 'self.node_ids', '0.8')
        return self.node_ids

    @property
    def node_ids(self):
        return self._nodeIDs(allow_empty_nodes=False)

    def get_edge_ids(self):
        return [tuple(sorted(self.node_ids))]

    def Rho(self):
        r"""returns the material density  \f$ \rho \f$"""
        return self.pid_ref.mid_ref.rho

    def Length(self):
        r"""
        Gets the length of the element.

        .. math:: L = \sqrt{  (n_{x2}-n_{x1})^2+(n_{y2}-n_{y1})^2+(n_{z2}-n_{z1})^2  }
        """
        L = norm(self.nodes_ref[1].get_position() - self.nodes_ref[0].get_position())
        return L

    def Mass(self):
        r"""
        get the mass of the element.

        .. math:: m = \left( \rho A + nsm \right) L
        """
        L = self.Length()
        mass = (self.Rho() * self.Area() + self.Nsm()) * L
        return mass


class CROD(RodElement):
    type = 'CROD'
    _field_map = {
        1: 'eid', 2:'pid',
    }

    def _update_field_helper(self, n, value):
        if n == 3:
            self.nodes[0] = value
        elif n == 4:
            self.nodes[1] = value
        else:
            raise KeyError('Field %r=%r is an invalid %s entry.' % (n, value, self.type))

    def __init__(self, eid, pid, nids, comment=''):
        """
        +------+-----+-----+----+----+
        | CROD | EID | PID | N1 | N2 |
        +------+-----+-----+----+----+
        """
        RodElement.__init__(self)
        if comment:
            self._comment = comment
        self.eid = eid
        self.pid = pid
        self.prepare_node_ids(nids)
        assert len(self.nodes) == 2

    @classmethod
    def add_card(cls, card, comment=''):
        eid = integer(card, 1, 'eid')
        pid = integer_or_blank(card, 2, 'pid', eid)
        nids = [integer(card, 3, 'n1'),
                integer(card, 4, 'n2')]
        assert len(card) == 5, 'len(CROD card) = %i\ncard=%s' % (len(card), str(card))
        return CROD(eid, pid, nids, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        eid = data[0]
        pid = data[1]
        nids = data[2:4]
        return CROD(eid, pid, nids, comment=comment)

    def _verify(self, xref=False):
        eid = self.Eid()
        pid = self.Pid()
        edges = self.get_edge_ids()
        assert isinstance(pid, int), 'pid=%r' % pid
        if xref:  # True
            mid = self.Mid()
            L = self.Length()
            A = self.Area()
            nsm = self.Nsm()
            mpa = self.MassPerLength()
            mass = self.Mass()
            assert isinstance(mid, integer_types), 'mid=%r' % mid
            assert isinstance(L, float), 'L=%r' % L
            assert isinstance(A, float), 'A=%r' % A
            assert isinstance(nsm, float), 'nsm=%r' % nsm
            assert isinstance(mpa, float), 'mass_per_length=%r' % mpa
            assert isinstance(mass, float), 'mass=%r' % mass

        c = self.Centroid()
        for i in range(3):
            assert isinstance(c[i], float), 'centroid[%i]=%r' % (i, c[i])

    def Centroid(self):
        return (self.nodes_ref[0].get_position() + self.nodes_ref[1].get_position()) / 2.

    def Eid(self):
        return self.eid

    def Mid(self):
        if isinstance(self.pid, integer_types):
            msg = 'Element eid=%i has not been cross referenced.\n%s' % (self.eid, str(self))
            raise RuntimeError(msg)
        return self.pid_ref.Mid()

    def Area(self):
        if isinstance(self.pid, integer_types):
            msg = 'Element eid=%i has not been cross referenced.\n%s' % (self.eid, str(self))
            raise RuntimeError(msg)
        return self.pid_ref.A

    def Nsm(self):
        if isinstance(self.pid, integer_types):
            msg = 'Element eid=%i has not been cross referenced.\n%s' % (self.eid, str(self))
            raise RuntimeError(msg)
        return self.pid_ref.nsm

    def E(self):
        if isinstance(self.pid, integer_types):
            msg = 'Element eid=%i has not been cross referenced.\n%s' % (self.eid, str(self))
            raise RuntimeError(msg)
        return self.pid_ref.mid_ref.E()

    def G(self):
        if isinstance(self.pid, integer_types):
            msg = 'Element eid=%i has not been cross referenced.\n%s' % (self.eid, str(self))
            raise RuntimeError(msg)
        return self.pid_ref.mid_ref.G()

    def J(self):
        if isinstance(self.pid, integer_types):
            msg = 'Element eid=%i has not been cross referenced.\n%s' % (self.eid, str(self))
            raise RuntimeError(msg)
        return self.pid_ref.J()

    def C(self):
        if isinstance(self.pid, integer_types):
            msg = 'Element eid=%i has not been cross referenced.\n%s' % (self.eid, str(self))
            raise RuntimeError(msg)
        return self.pid_ref.c

    def MassPerLength(self):
        if isinstance(self.pid, integer_types):
            msg = 'Element eid=%i has not been cross referenced.\n%s' % (self.eid, str(self))
            raise RuntimeError(msg)
        massPerLength = self.pid_ref.mid_ref.rho * self.pid_ref.A + self.pid_ref.nsm
        return massPerLength

    def raw_fields(self):
        list_fields = ['CROD', self.eid, self.Pid()] + self.node_ids
        return list_fields

    def repr_fields(self):
        return self.raw_fields()

    def write_card(self, size=8, is_double=False):
        card = self.raw_fields()
        return self.comment + print_card_8(card)


class CTUBE(RodElement):
    type = 'CTUBE'
    _field_map = {
        1: 'eid', 2:'pid',
    }

    def _update_field_helper(self, n, value):
        if n == 3:
            self.nodes[0] = value
        elif n == 4:
            self.nodes[1] = value
        else:
            raise KeyError('Field %r=%r is an invalid %s entry.' % (n, value, self.type))

    def __init__(self, eid, pid, nids, comment=''):
        RodElement.__init__(self)
        if comment:
            self._comment = comment
        self.eid = eid
        self.pid = pid
        self.prepare_node_ids(nids)
        assert len(self.nodes) == 2

    @classmethod
    def add_card(cls, card, comment=''):
        eid = integer(card, 1, 'eid')
        pid = integer_or_blank(card, 2, 'pid', eid)
        nids = [integer(card, 3, 'n1'),
                integer(card, 4, 'n2')]
        assert len(card) == 5, 'len(CTUBE card) = %i\ncard=%s' % (len(card), card)
        return CTUBE(eid, pid, nids, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        eid = data[0]
        pid = data[1]
        nids = data[2:4]
        return CTUBE(eid, pid, nids, comment=comment)

    def _verify(self, xref=False):
        pid = self.Pid()
        A = self.Area()
        edges = self.get_edge_ids()
        assert isinstance(pid, int), 'pid=%r' % pid
        assert isinstance(A, float), 'A=%r' % A
        if xref:
            L = self.Length()
            nsm = self.Nsm()
            assert isinstance(L, float), 'L=%r' % L
            assert isinstance(nsm, float), 'nsm=%r' % nsm
            if self.pid_ref.mid_ref.type == 'MAT1':
                mpa = self.pid_ref.MassPerLength()
                mass = self.Mass()
                assert isinstance(mpa, float), 'mass_per_length=%r' % mpa
                assert isinstance(mass, float), 'mass=%r' % mass
            elif self.pid_ref.mid_ref.type == 'MAT4':
                pass
            else:
                raise NotImplementedError('_verify does not support self.pid_ref.mid_ref.type=%s' % self.pid_ref.mid_ref.type)

        c = self.Centroid()
        for i in range(3):
            assert isinstance(c[i], float), 'centroid[%i]=%r' % (i, c[i])

    def Eid(self):
        return self.eid

    def Mid(self):
        return self.pid_ref.Mid()

    def Mass(self):
        return self.pid_ref.MassPerLength() * self.Length()

    def Nsm(self):
        if isinstance(self.pid, integer_types):
            raise RuntimeError('Element eid=%i has not been cross referenced.\n%s' % (self.eid, str(self)))
        return self.pid_ref.Nsm()

    def Area(self):
        if isinstance(self.pid, integer_types):
            raise RuntimeError('Element eid=%i has not been cross referenced.\n%s' % (self.eid, str(self)))
        return self.pid_ref.Area()

    def E(self):
        if isinstance(self.pid, integer_types):
            raise RuntimeError('Element eid=%i has not been cross referenced.\n%s' % (self.eid, str(self)))
        return self.pid_ref.mid_ref.E()

    def G(self):
        if isinstance(self.pid, integer_types):
            raise RuntimeError('Element eid=%i has not been cross referenced.\n%s' % (self.eid, str(self)))
        return self.pid_ref.mid_ref.G()

    def J(self):
        if isinstance(self.pid, integer_types):
            raise RuntimeError('Element eid=%i has not been cross referenced.\n%s' % (self.eid, str(self)))
        return self.pid_ref.J()

    def Centroid(self):
        return (self.nodes_ref[0].get_position() + self.nodes_ref[1].get_position()) / 2.

    def raw_fields(self):
        list_fields = ['CTUBE', self.eid, self.Pid()] + self.node_ids
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class CONROD(RodElement):
    type = 'CONROD'
    _field_map = {
        1: 'eid', 4:'mid', 5:'A', 6:'j', 7:'c', 8:'nsm',
    }

    def _update_field_helper(self, n, value):
        if n == 2:
            self.nodes[0] = value
        elif n == 3:
            self.nodes[1] = value
        else:
            raise KeyError('Field %r=%r is an invalid %s entry.' % (n, value, self.type))

    def __init__(self, eid, mid, nids, A, j=0.0, c=0.0, nsm=0.0, comment=''):
        RodElement.__init__(self)
        if comment:
            self._comment = comment
        self.eid = eid
        self.mid = mid
        self.A = A
        self.j = j
        self.c = c
        self.nsm = nsm
        self.prepare_node_ids(nids)
        assert len(self.nodes) == 2

    @classmethod
    def add_card(cls, card, comment=''):
        eid = integer(card, 1, 'eid')
        nids = [integer(card, 2, 'n1'),
                integer(card, 3, 'n2')]
        mid = integer(card, 4, 'mid')
        A = double(card, 5, 'A')
        j = double_or_blank(card, 6, 'j', 0.0)
        c = double_or_blank(card, 7, 'c', 0.0)
        nsm = double_or_blank(card, 8, 'nsm', 0.0)
        return CONROD(eid, mid, nids, A, j, c, nsm, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        eid = data[0]
        nids = data[1:3]
        mid = data[3]
        A = data[4]
        j = data[5]
        c = data[6]
        nsm = data[7]
        return CONROD(eid, mid, nids, A, j, c, nsm, comment=comment)

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ' which is required by %s eid=%s' % (self.type, self.eid)
        self.nodes = model.Nodes(self.nodes, msg=msg)
        self.nodes_ref = self.nodes
        self.mid = model.Material(self.mid, msg=msg)
        self.mid_ref = self.mid

    def uncross_reference(self):
        self.nodes = self.node_ids
        self.mid = self.Mid()
        del self.nodes_ref, self.mid_ref

    def _verify(self, xref=False):
        pid = self.Pid()
        assert pid is None, 'pid=%r' % pid
        edges = self.get_edge_ids()
        if xref:  # True
            mid = self.Mid()
            L = self.Length()
            A = self.Area()
            nsm = self.Nsm()
            mpa = self.MassPerLength()
            mass = self.Mass()
            assert isinstance(mid, integer_types), 'mid=%r' % mid
            assert isinstance(L, float), 'L=%r' % L
            assert isinstance(A, float), 'A=%r' % A
            assert isinstance(nsm, float), 'nsm=%r' % nsm
            assert isinstance(mpa, float), 'mass_per_length=%r' % mpa
            assert isinstance(mass, float), 'mass=%r' % mass

            c = self.Centroid()
            for i in range(3):
                assert isinstance(c[i], float), 'centroid[%i]=%r' % (i, c[i])

    def Centroid(self):
        return (self.nodes_ref[0].get_position() + self.nodes_ref[1].get_position()) / 2.

    def Mid(self):
        if isinstance(self.mid, integer_types):
            return self.mid
        #elif self.mid is None:
            #print ("No material defined for element ", self.eid)
            #return None
        else:
            return self.mid_ref.mid

    def Eid(self):
        return self.eid

    def Pid(self):
        return None

    def MassPerLength(self):
        if isinstance(self.mid, integer_types):
            raise RuntimeError('Element eid=%i has not been cross referenced.\n%s' % (self.eid, str(self)))
        massPerLength = self.mid_ref.rho * self.A + self.nsm
        return massPerLength

    def C(self):
        """torsional constant"""
        return self.c

    def Area(self):
        return self.A

    def J(self):
        r"""returns the Polar Moment of Inertia, :math:`J`"""
        return self.j

    def Nsm(self):
        """Placeholder method for the non-structural mass"""
        return self.nsm

    def E(self):
        r"""returns the Young's Modulus, :math:`E`$"""
        if isinstance(self.mid, integer_types):
            raise RuntimeError('Element eid=%i has not been cross referenced.\n%s' % (self.eid, str(self)))
        return self.mid_ref.E()

    def G(self):
        r"""returns the Shear Modulus, :math:`G`"""
        if isinstance(self.mid, integer_types):
            raise RuntimeError('Element eid=%i has not been cross referenced.\n%s' % (self.eid, str(self)))
        return self.mid_ref.G()

    def Rho(self):
        r"""returns the material density, :math:`\rho`"""
        if isinstance(self.mid, integer_types):
            raise RuntimeError('Element eid=%i has not been cross referenced.\n%s' % (self.eid, str(self)))
        return self.mid_ref.rho

    #def writeCodeAster(self):
        #msg = ''
        #msg += "    POUTRE=_F(GROUP_MA='CONROD_%s',\n" % self.eid
        #msg += "              SECTION='CERCLE',  # circular section\n"
        #if self.Thickness():
            #msg += "              CARA=('R','EP'),   # radius, thickness\n"
            #msg += "              VALE=(%g,%g),\n" % (
                #self.Radius(), self.Thickness())
        #else:
            #msg += "              CARA=('R')   # radius\n"
            #msg += "              VALE=(%g),\n" % self.Radius()
        #return msg

    def raw_fields(self):
        list_fields = [
            'CONROD', self.eid] + self.node_ids + [
                self.Mid(), self.A, self.j, self.c, self.nsm]
        return list_fields

    def repr_fields(self):
        j = set_blank_if_default(self.j, 0.0)
        c = set_blank_if_default(self.c, 0.0)
        nsm = set_blank_if_default(self.nsm, 0.0)
        list_fields = [
            'CONROD', self.eid] + self.node_ids + [self.Mid(), self.A, j, c, nsm]
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)
