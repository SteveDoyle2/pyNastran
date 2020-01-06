# pylint: disable=C0103
from __future__ import annotations
from typing import TYPE_CHECKING
from numpy.linalg import norm  # type: ignore

from pyNastran.utils.numpy_utils import integer_types
from pyNastran.bdf.field_writer_8 import set_blank_if_default
from pyNastran.bdf.cards.base_card import Element #, Mid
from pyNastran.bdf.bdf_interface.assign_type import (
    integer, integer_or_blank, double_or_blank)
from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.field_writer_16 import print_card_16
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf import BDF


class RodElement(Element):  # CROD, CONROD, CTUBE

    def __init__(self):
        Element.__init__(self)

    @property
    def node_ids(self):
        return self._node_ids(nodes=self.nodes_ref, allow_empty_nodes=False)

    def get_edge_ids(self):
        return [tuple(sorted(self.node_ids))]

    def Mass(self):
        r"""
        get the mass of the element.

        .. math:: m = \left( \rho A + nsm \right) L
        """
        L = self.Length()
        mass = (self.Rho() * self.Area() + self.Nsm()) * L
        return mass


class CROD(RodElement):
    """
    +------+-----+-----+----+----+
    |   1  |  2  |  3  |  4 |  5 |
    +======+=====+=====+====+====+
    | CROD | EID | PID | N1 | N2 |
    +------+-----+-----+----+----+
    """
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

    @classmethod
    def export_to_hdf5(cls, h5_file, model, eids):
        """exports the elements in a vectorized way"""
        #comments = []
        pids = []
        nodes = []
        for eid in eids:
            element = model.elements[eid]
            #comments.append(element.comment)
            pids.append(element.pid)
            nodes.append(element.nodes)
        #h5_file.create_dataset('_comment', data=comments)
        h5_file.create_dataset('eid', data=eids)
        h5_file.create_dataset('pid', data=pids)
        h5_file.create_dataset('nodes', data=nodes)

    def __init__(self, eid, pid, nids, comment=''):
        """
        Creates a CROD card

        Parameters
        ----------
        eid : int
            element id
        pid : int
            property id (PROD)
        nids : List[int, int]
            node ids
        comment : str; default=''
            a comment for the card
        """
        RodElement.__init__(self)
        if comment:
            self.comment = comment
        self.eid = eid
        self.pid = pid
        self.nodes = self.prepare_node_ids(nids)
        assert len(self.nodes) == 2
        self.nodes_ref = None
        self.pid_ref = None

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a CROD card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        eid = integer(card, 1, 'eid')
        pid = integer_or_blank(card, 2, 'pid', eid)
        nids = [integer(card, 3, 'n1'),
                integer(card, 4, 'n2')]
        assert len(card) == 5, 'len(CROD card) = %i\ncard=%s' % (len(card), str(card))
        return CROD(eid, pid, nids, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        """
        Adds a CROD card from the OP2

        Parameters
        ----------
        data : List[varies]
            a list of fields defined in OP2 format
        comment : str; default=''
            a comment for the card
        """
        eid = data[0]
        pid = data[1]
        nids = data[2:4]
        return CROD(eid, pid, nids, comment=comment)


    def cross_reference(self, model: BDF) -> None:
        msg = ', which is required by CROD eid=%s' % (self.eid)
        self.nodes_ref = model.Nodes(self.nodes, msg=msg)
        self.pid_ref = model.Property(self.pid, msg=msg)

    def safe_cross_reference(self, model, xref_errors):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ', which is required by CROD eid=%s' % self.eid
        self.nodes_ref = model.Nodes(self.node_ids, msg=msg)
        self.pid_ref = model.safe_property(self.pid, self.eid, xref_errors, msg=msg)

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.nodes = self.node_ids
        self.pid = self.Pid()
        self.nodes_ref = None
        self.pid_ref = None

    def _verify(self, xref):
        eid = self.eid
        pid = self.Pid()
        unused_edges = self.get_edge_ids()
        assert isinstance(eid, int), 'eid=%r' % eid
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

    def Length(self):
        r"""
        Gets the length of the element.

        .. math:: L = \sqrt{  (n_{x2}-n_{x1})^2+(n_{y2}-n_{y1})^2+(n_{z2}-n_{z1})^2  }
        """
        if self.pid_ref is None:
            msg = 'Element eid=%i has not been cross referenced.\n%s' % (self.eid, str(self))
            raise RuntimeError(msg)
        L = norm(self.nodes_ref[1].get_position() - self.nodes_ref[0].get_position())
        return L

    def Rho(self):
        r"""returns the material density  \f$ \rho \f$"""
        if self.pid_ref is None:
            msg = 'Element eid=%i has not been cross referenced.\n%s' % (self.eid, str(self))
            raise RuntimeError(msg)
        return self.pid_ref.mid_ref.rho

    def Centroid(self):
        return (self.nodes_ref[0].get_position() + self.nodes_ref[1].get_position()) / 2.

    def center_of_mass(self):
        return self.Centroid()

    def Mid(self):
        if self.pid_ref is None:
            msg = 'Element eid=%i has not been cross referenced.\n%s' % (self.eid, str(self))
            raise RuntimeError(msg)
        return self.pid_ref.Mid()

    def Area(self):
        if self.pid_ref is None:
            msg = 'Element eid=%i has not been cross referenced.\n%s' % (self.eid, str(self))
            raise RuntimeError(msg)
        return self.pid_ref.A

    def Nsm(self):
        if self.pid_ref is None:
            msg = 'Element eid=%i has not been cross referenced.\n%s' % (self.eid, str(self))
            raise RuntimeError(msg)
        return self.pid_ref.nsm

    def E(self):
        if self.pid_ref is None:
            msg = 'Element eid=%i has not been cross referenced.\n%s' % (self.eid, str(self))
            raise RuntimeError(msg)
        return self.pid_ref.mid_ref.E()

    def G(self):
        if self.pid_ref is None:
            msg = 'Element eid=%i has not been cross referenced.\n%s' % (self.eid, str(self))
            raise RuntimeError(msg)
        return self.pid_ref.mid_ref.G()

    def J(self):
        if self.pid_ref is None:
            msg = 'Element eid=%i has not been cross referenced.\n%s' % (self.eid, str(self))
            raise RuntimeError(msg)
        return self.pid_ref.J()

    def C(self):
        if self.pid_ref is None:
            msg = 'Element eid=%i has not been cross referenced.\n%s' % (self.eid, str(self))
            raise RuntimeError(msg)
        return self.pid_ref.c

    def MassPerLength(self):
        if self.pid_ref is None:
            msg = 'Element eid=%i has not been cross referenced.\n%s' % (self.eid, str(self))
            raise RuntimeError(msg)
        mass_per_length = self.pid_ref.mid_ref.rho * self.pid_ref.A + self.pid_ref.nsm
        return mass_per_length

    def raw_fields(self):
        list_fields = ['CROD', self.eid, self.Pid()] + self.node_ids
        return list_fields

    def repr_fields(self):
        return self.raw_fields()

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.raw_fields()
        return self.comment + print_card_8(card)

    def write_card_16(self, is_double=False):
        card = self.raw_fields()
        return self.comment + print_card_16(card)

class CTUBE(RodElement):
    """
    +-------+-----+-----+----+----+
    |   1   |  2  |  3  |  4 |  5 |
    +=======+=====+=====+====+====+
    | CTUBE | EID | PID | N1 | N2 |
    +-------+-----+-----+----+----+
    """
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

    @classmethod
    def export_to_hdf5(cls, h5_file, model, eids):
        """exports the elements in a vectorized way"""
        #comments = []
        pids = []
        nodes = []
        for eid in eids:
            element = model.elements[eid]
            #comments.append(element.comment)
            pids.append(element.pid)
            nodes.append(element.nodes)
        #h5_file.create_dataset('_comment', data=comments)
        h5_file.create_dataset('eid', data=eids)
        h5_file.create_dataset('pid', data=pids)
        h5_file.create_dataset('nodes', data=nodes)

    def __init__(self, eid, pid, nids, comment=''):
        """
        Creates a CTUBE card

        Parameters
        ----------
        eid : int
            element id
        pid : int
            property id
        nids : List[int, int]
            node ids
        comment : str; default=''
            a comment for the card
        """
        RodElement.__init__(self)
        if comment:
            self.comment = comment
        self.eid = eid
        self.pid = pid
        self.nodes = self.prepare_node_ids(nids)
        assert len(self.nodes) == 2
        self.nodes_ref = None
        self.pid_ref = None

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a CTUBE card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        eid = integer(card, 1, 'eid')
        pid = integer_or_blank(card, 2, 'pid', eid)
        nids = [integer(card, 3, 'n1'),
                integer(card, 4, 'n2')]
        assert len(card) == 5, 'len(CTUBE card) = %i\ncard=%s' % (len(card), card)
        return CTUBE(eid, pid, nids, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        """
        Adds a CTUBE card from the OP2

        Parameters
        ----------
        data : List[varies]
            a list of fields defined in OP2 format
        comment : str; default=''
            a comment for the card
        """
        eid = data[0]
        pid = data[1]
        nids = data[2:4]
        return CTUBE(eid, pid, nids, comment=comment)

    def cross_reference(self, model: BDF) -> None:
        msg = ', which is required by CTUBE eid=%s' % (self.eid)
        self.nodes_ref = model.Nodes(self.nodes, msg=msg)
        self.pid_ref = model.Property(self.pid, msg=msg)

    def safe_cross_reference(self, model, xref_errors):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ', which is required by CTUBE eid=%s' % self.eid
        self.nodes_ref = model.Nodes(self.node_ids, msg=msg)
        self.pid_ref = model.safe_property(self.pid, self.eid, xref_errors, msg=msg)
        ## TODO: xref coord

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.nodes = self.node_ids
        self.pid = self.Pid()
        self.nodes_ref = None
        self.pid_ref = None

    def _verify(self, xref):
        pid = self.Pid()
        unused_edges = self.get_edge_ids()
        assert isinstance(pid, int), 'pid=%r' % pid
        if xref:
            A = self.Area()
            L = self.Length()
            nsm = self.Nsm()
            assert isinstance(A, float), 'A=%r' % A
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
                msg = '_verify does not support self.pid_ref.mid_ref.type=%s' % (
                    self.pid_ref.mid_ref.type)
                raise NotImplementedError(msg)

            c = self.Centroid()
            for i in range(3):
                assert isinstance(c[i], float), 'centroid[%i]=%r' % (i, c[i])

    def Length(self):
        r"""
        Gets the length of the element.

        .. math:: L = \sqrt{  (n_{x2}-n_{x1})^2+(n_{y2}-n_{y1})^2+(n_{z2}-n_{z1})^2  }
        """
        if self.pid_ref is None:
            msg = 'Element eid=%i has not been cross referenced.\n%s' % (self.eid, str(self))
            raise RuntimeError(msg)
        L = norm(self.nodes_ref[1].get_position() - self.nodes_ref[0].get_position())
        return L

    def Rho(self):
        r"""returns the material density  \f$ \rho \f$"""
        if self.pid_ref is None:
            msg = 'Element eid=%i has not been cross referenced.\n%s' % (self.eid, str(self))
            raise RuntimeError(msg)
        return self.pid_ref.mid_ref.rho

    def Mid(self):
        if self.pid_ref is None:
            msg = 'Element eid=%i has not been cross referenced.\n%s' % (self.eid, str(self))
            raise RuntimeError(msg)
        return self.pid_ref.Mid()

    def Mass(self):
        if self.pid_ref is None:
            msg = 'Element eid=%i has not been cross referenced.\n%s' % (self.eid, str(self))
            raise RuntimeError(msg)
        return self.pid_ref.MassPerLength() * self.Length()

    def Nsm(self):
        if self.pid_ref is None:
            msg = 'Element eid=%i has not been cross referenced.\n%s' % (self.eid, str(self))
            raise RuntimeError(msg)
        return self.pid_ref.Nsm()

    def Area(self):
        if self.pid_ref is None:
            msg = 'Element eid=%i has not been cross referenced.\n%s' % (self.eid, str(self))
            raise RuntimeError(msg)
        return self.pid_ref.Area()

    def E(self):
        if self.pid_ref is None:
            msg = 'Element eid=%i has not been cross referenced.\n%s' % (self.eid, str(self))
            raise RuntimeError(msg)
        return self.pid_ref.mid_ref.E()

    def G(self):
        if self.pid_ref is None:
            msg = 'Element eid=%i has not been cross referenced.\n%s' % (self.eid, str(self))
            raise RuntimeError(msg)
        return self.pid_ref.mid_ref.G()

    def J(self):
        if self.pid_ref is None:
            msg = 'Element eid=%i has not been cross referenced.\n%s' % (self.eid, str(self))
            raise RuntimeError(msg)
        return self.pid_ref.J()

    def Centroid(self):
        return (self.nodes_ref[0].get_position() + self.nodes_ref[1].get_position()) / 2.

    def center_of_mass(self):
        return self.Centroid()

    def raw_fields(self):
        list_fields = ['CTUBE', self.eid, self.Pid()] + self.node_ids
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class CONROD(RodElement):
    """
    +--------+-----+-----+----+-----+---+---+---+-----+
    |   1    |  2  |  3  |  4 |  5  | 6 | 7 | 8 |  9  |
    +========+=====+=====+====+=====+===+===+===+=====+
    | CONROD | EID | N1  | N2 | MID | A | J | C | NSM |
    +--------+-----+-----+----+-----+---+---+---+-----+
    """
    type = 'CONROD'
    pid = -10 # 10 is the element type per DMAP
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

    @classmethod
    def export_to_hdf5(cls, h5_file, model, eids):
        """exports the elements in a vectorized way"""
        #comments = []
        nodes = []
        mids = []
        A = []
        J = []
        c = []
        nsm = []
        for eid in eids:
            element = model.elements[eid]
            #comments.append(element.comment)
            mids.append(element.mid)
            nodes.append(element.nodes)
            A.append(element.A)
            J.append(element.j)
            c.append(element.c)
            nsm.append(element.nsm)
        #h5_file.create_dataset('_comment', data=comments)
        h5_file.create_dataset('eid', data=eids)
        h5_file.create_dataset('nodes', data=nodes)
        h5_file.create_dataset('mid', data=mids)
        h5_file.create_dataset('A', data=A)
        h5_file.create_dataset('J', data=J)
        h5_file.create_dataset('c', data=c)
        h5_file.create_dataset('nsm', data=nsm)

    def __init__(self, eid, mid, nids, A=0.0, j=0.0, c=0.0, nsm=0.0, comment=''):
        """
        Creates a CONROD card

        Parameters
        ----------
        eid : int
            element id
        mid : int
            material id
        nids : List[int, int]
            node ids
        A : float
            area
        j : float; default=0.
            polar moment of inertia
        c : float; default=0.
            stress factor
        nsm : float; default=0.
            non-structural mass per unit length
        comment : str; default=''
            a comment for the card
        """
        RodElement.__init__(self)
        if comment:
            self.comment = comment
        self.eid = eid
        self.mid = mid
        self.A = A
        self.j = j
        self.c = c
        self.nsm = nsm
        self.nodes = self.prepare_node_ids(nids)
        assert len(self.nodes) == 2
        self.nodes_ref = None
        self.mid_ref = None

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a CONROD card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        eid = integer(card, 1, 'eid')
        nids = [integer(card, 2, 'n1'),
                integer(card, 3, 'n2')]
        mid = integer(card, 4, 'mid')
        A = double_or_blank(card, 5, 'A', 0.0)
        j = double_or_blank(card, 6, 'j', 0.0)
        c = double_or_blank(card, 7, 'c', 0.0)
        nsm = double_or_blank(card, 8, 'nsm', 0.0)
        assert len(card) <= 9, 'len(CONROD card) = %i\ncard=%s' % (len(card), str(card))
        return CONROD(eid, mid, nids, A, j, c, nsm, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        """
        Adds a CONROD card from the OP2

        Parameters
        ----------
        data : List[varies]
            a list of fields defined in OP2 format
        comment : str; default=''
            a comment for the card
        """
        eid = data[0]
        nids = data[1:3]
        mid = data[3]
        A = data[4]
        j = data[5]
        c = data[6]
        nsm = data[7]
        return CONROD(eid, mid, nids, A, j, c, nsm, comment=comment)

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ', which is required by CONROD eid=%s' % (self.eid)
        self.nodes_ref = model.Nodes(self.nodes, msg=msg)
        self.mid_ref = model.Material(self.mid, msg=msg)

    def safe_cross_reference(self, model, xref_errors):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ', which is required by CONROD eid=%s' % self.eid
        self.nodes_ref = model.Nodes(self.node_ids, msg=msg)
        self.mid_ref = model.safe_material(self.mid, self.eid, xref_errors, msg=msg)

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.nodes = self.node_ids
        self.mid = self.Mid()
        self.nodes_ref = None
        self.mid_ref = None

    def _verify(self, xref):
        pid = self.Pid()
        assert pid == -10, 'pid=%r' % pid
        unused_edges = self.get_edge_ids()
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

    def Length(self):
        r"""
        Gets the length of the element.

        .. math:: L = \sqrt{  (n_{x2}-n_{x1})^2+(n_{y2}-n_{y1})^2+(n_{z2}-n_{z1})^2  }
        """
        if self.mid_ref is None:
            msg = 'Element eid=%i has not been cross referenced.\n%s' % (self.eid, str(self))
            raise RuntimeError(msg)
        L = norm(self.nodes_ref[1].get_position() - self.nodes_ref[0].get_position())
        return L

    def Rho(self):
        r"""returns the material density  \f$ \rho \f$"""
        if self.mid_ref is None:
            msg = 'Element eid=%i has not been cross referenced.\n%s' % (self.eid, str(self))
            raise RuntimeError(msg)
        return self.mid_ref.rho

    def Centroid(self):
        """Get the centroid of the element (save as the center of mass for the CONROD)"""
        if self.mid_ref is None:
            msg = 'Element eid=%i has not been cross referenced.\n%s' % (self.eid, str(self))
            raise RuntimeError(msg)
        return (self.nodes_ref[0].get_position() + self.nodes_ref[1].get_position()) / 2.

    def center_of_mass(self):
        """Get the center of mass of the element (save as the centroid for the CONROD)"""
        return self.Centroid()

    def Mid(self):
        if self.mid_ref is None:
            return self.mid
        #elif self.mid is None:
            #print ("No material defined for element ", self.eid)
            #return None
        return self.mid_ref.mid

    def Pid(self):
        """Spoofs the property id for the CONROD"""
        return self.pid

    def MassPerLength(self):
        """Gets the mass per length of the CONROD"""
        if self.mid_ref is None:
            msg = 'Element eid=%i has not been cross referenced.\n%s' % (self.eid, str(self))
            raise RuntimeError(msg)
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
        if self.mid_ref is None:
            msg = 'Element eid=%i has not been cross referenced.\n%s' % (self.eid, str(self))
            raise RuntimeError(msg)
        return self.mid_ref.E()

    def G(self):
        r"""returns the Shear Modulus, :math:`G`"""
        if self.mid_ref is None:
            msg = 'Element eid=%i has not been cross referenced.\n%s' % (self.eid, str(self))
            raise RuntimeError(msg)
        return self.mid_ref.G()

    #def write_code_aster(self):
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

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)
