"""
All nodes are defined in this file.  This includes:

 * Node
   * XPoint
     * EPOINT
     * SPOINT
   * XPoints
     * EPOINTs
     * SPOINTs
   * GRID
   * GRDSET
   * GRIDB
 * POINT
 * Ring
   * RINGAX

All ungrouped elements are Node objects.

The EPOINT/SPOINT classes refer to a single EPOINT/SPOINT.  The
EPOINTs/SPOINTs classes are for multiple degrees of freedom
(e.g. an SPOINT card).
"""
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
import numpy as np
from six import string_types, PY2
if PY2:
    u = unicode
else:
    u = str

from pyNastran.utils import integer_types
from pyNastran.bdf.field_writer_8 import set_string8_blank_if_default
from pyNastran.bdf.field_writer_16 import set_string16_blank_if_default

from pyNastran.bdf.field_writer_8 import set_blank_if_default
from pyNastran.bdf.cards.base_card import BaseCard, expand_thru
from pyNastran.bdf.cards.collpase_card import collapse_thru_packs
from pyNastran.bdf.bdf_interface.assign_type import (
    integer, integer_or_blank, double, double_or_blank, blank, integer_or_string)
from pyNastran.bdf.field_writer_8 import print_card_8, print_float_8, print_int_card
from pyNastran.bdf.field_writer_16 import print_float_16, print_card_16
from pyNastran.bdf.field_writer_double import print_scientific_double, print_card_double

class Ring(BaseCard):
    """Generic Ring base class"""
    def __init__(self):
        pass

class Node(BaseCard):
    """Generic Node base class"""
    def __init__(self):
        pass

class RINGAX(Ring):
    """
    Defines a ring for conical shell problems.

    +-------+-----+-----+-----+----+-----+-----+------+-----+
    |   1   |  2  |  3  |  4  |  5 |  6  |  7  |  8   |  9  |
    +=======+=====+=====+=====+====+=====+=====+======+=====+
    |RINGAX | MID |     |  R  |  Z |     |     | PS   |     |
    +-------+-----+-----+-----+----+-----+-----+------+-----+
    """
    type = 'RINGAX'

    #: allows the get_field method and update_field methods to be used
    _field_map = {1: 'mid', 3:'R', 4:'z', 7:'ps'}

    def __init__(self, card=None, data=None, comment=''):  # this card has missing fields
        """
        Creates the RINGAX card

        :param card:
          a BDFCard object
        :type card:
          BDFCard
        :param data:
          a list with the RINGAX fields defined in OP2 format
        :type data:
          LIST
        :param comment:
          a comment for the card
        :type comment:
          string
        """
        Ring.__init__(self)
        if comment:
            self._comment = comment
        if card:
            #: Node ID
            self.nid = integer(card, 1, 'nid')
            blank(card, 2, 'blank')

            #: Radius
            self.R = double(card, 3, 'R')
            self.z = double(card, 4, 'z')
            blank(card, 5, 'blank')
            blank(card, 6, 'blank')

            #: local SPC constraint
            self.ps = integer_or_blank(card, 7, 'ps')
            assert len(card) <= 8, 'len(RINGAX card) = %i\ncard=%s' % (len(card), card)
        else:  # hasn't been validated
            self.nid = data[0]
            self.R = data[1]
            self.z = data[2]
            self.ps = data[3]
            assert len(data) == 4, data

    def raw_fields(self):
        """
        Gets the fields in their unmodified form

        :returns fields:
          the fields that define the card
        :type fields:
          LIST
        """
        list_fields = ['RINGAX', self.nid, None, self.R, self.z, None,
                       None, self.ps]
        return list_fields

    def write_card(self, size=8, is_double=False):
        """
        The writer method used by BDF.write_card

        :param size:
          the size of the card (8/16)
        :type size:
          int
        """
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)

class XPoint(Node):
    """common class for EPOINT/SPOINT"""
    def __init__(self, nid, comment):
        Node.__init__(self)
        if comment:
            self._comment = comment
        self.nid = nid
        assert isinstance(nid, integer_types), nid

    @property
    def type(self):
        """dummy method for EPOINT/SPOINT classes"""
        raise NotImplementedError('This method should be overwritten by the parent class')

    def raw_fields(self):
        """
        Gets the fields in their unmodified form

        :returns fields:
          the fields that define the card
        :type fields:
          LIST
        """
        lists_fields = []
        if isinstance(self.nid, integer_types):
            list_fields = [self.type, self.nid]
            lists_fields.append(list_fields)
        else:
            singles, doubles = collapse_thru_packs(self.nid)
            if singles:
                list_fields = [self.type] + singles
            if doubles:
                for spoint_double in doubles:
                    list_fields = [self.type] + spoint_double
                    lists_fields.append(list_fields)
        return lists_fields

    def write_card(self, size=8, is_double=False):
        """
        The writer method used by BDF.write_card

        :param size:   unused
        :param is_double: unused
        """
        msg = self.comment
        lists_fields = self.repr_fields()
        for list_fields in lists_fields:
            if 'THRU' not in list_fields:
                msg += print_int_card(list_fields)
            else:
                msg += print_card_8(list_fields)
        return msg

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        pass

class SPOINT(XPoint):
    """defines the SPOINT class"""
    type = 'SPOINT'

    def __init__(self, nid, comment=''):
        """
        Creates the SPOINT card

        :param card:
          a BDFCard object
        :type card:
          BDFCard
        :param data:
          a list with the SPOINT fields defined in OP2 format
        :type data:
          LIST
        :param comment:
          a comment for the card
        :type comment:
          string
        """
        XPoint.__init__(self, nid, comment)


class EPOINT(XPoint):
    """defines the EPOINT class"""
    type = 'EPOINT'

    def __init__(self, nid, comment=''):
        """
        Creates the EPOINT card

        :param card:
          a BDFCard object
        :type card:
          BDFCard
        :param data:
          a list with the SPOINT fields defined in OP2 format
        :type data:
          LIST
        :param comment:
          a comment for the card
        :type comment:
          string
        """
        XPoint.__init__(self, nid, comment)


def _get_compressed_xpoints(Type, xpoints):
    """
    Gets the spoints in sorted, short form.

      uncompresed:  SPOINT,1,3,5
      compressed:   SPOINT,1,3,5

      uncompresed:  SPOINT,1,2,3,4,5
      compressed:   SPOINT,1,THRU,5

      uncompresed:  SPOINT,1,2,3,4,5,7
      compressed:   SPOINT,7
                    SPOINT,1,THRU,5

    Type = 'SPOINT'
    spoints = [1, 2, 3, 4, 5]
    fields = _get_compressed_xpoints(Type, spoints)
    >>> fields
    ['SPOINT', 1, 'THRU', 5]
    """
    spoints = list(xpoints)
    spoints.sort()

    singles, doubles = collapse_thru_packs(spoints)

    lists_fields = []
    if singles:
        list_fields = [Type] + singles
        lists_fields.append(list_fields)
    if doubles:
        for spoint_double in doubles:
            list_fields = [Type] + spoint_double
            lists_fields.append(list_fields)
    return lists_fields

class XPoints(Node):
    """common class for EPOINTs and SPOINTs"""

    @property
    def type(self):
        """dummy method for EPOINTs/SPOINTs classes"""
        raise NotImplementedError('This method should be overwritten by the parent class')

    def __init__(self, ids, comment=''):
        Node.__init__(self)
        if comment:
            self._comment = comment
        self.points = ids

    @classmethod
    def add_card(cls, card, comment=''):
        fields = []
        for i in range(1, len(card)):
            field = integer_or_string(card, i, 'ID%i' % i)
            fields.append(field)
        points = set(expand_thru(fields))
        return cls(points, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        fields = data
        assert isinstance(data, list), data
        assert isinstance(data[0], int), data
        points = set(expand_thru(fields))
        return cls(points, comment=comment)

    def __len__(self):
        """
        Returns the number of degrees of freedom for the EPOINTs/SPOINTs class

        Returns
        -------
        ndofs : int
            the number of degrees of freedom
        """
        return len(self.points)

    def get_ndof(self):  # TODO: deprecate?
        return self.__len__()

    def add_points(self, sList):
        """
        Adds more EPOINTs/SPOINTs to this object
        """
        self.points = self.points.union(set(sList))

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        pass

    def raw_fields(self):
        """
        Gets the fields in their unmodified form

        :returns fields:
          the fields that define the card
        :type fields:
          LIST
        """
        points = list(self.points)
        points.sort()
        return [self.type] + points

    def write_card(self, size=8, is_double=False):
        """
        The writer method used by BDF.write_card

        :param size:   unused
        :param is_double: unused
        """
        lists_fields = _get_compressed_xpoints(self.type, self.points)
        msg = self.comment
        for list_fields in lists_fields:
            if 'THRU' not in list_fields:
                msg += print_int_card(list_fields)
            else:
                msg += print_card_8(list_fields)
        return msg

class SPOINTs(XPoints):
    """
    +--------+-----+------+-----+-----+-----+-----+-----+-----+
    |   1    |  2  |  3   |  4  |  5  |  6  |  7  |  8  |  9  |
    +========+=====+======+=====+=====+=====+=====+=====+=====+
    | SPOINT | ID1 | THRU | ID2 |     |     |     |     |     |
    +--------+-----+------+-----+-----+-----+-----+-----+-----+
    | SPOINT | ID1 | ID1  | ID3 | ID4 | ID5 | ID6 | ID7 | ID8 |
    +--------+-----+------+-----+-----+-----+-----+-----+-----+
    |        | ID8 | etc. |     |     |     |     |     |     |
    +--------+-----+------+-----+-----+-----+-----+-----+-----+
    """
    type = 'SPOINT'

    def __init__(self, ids, comment=''):
        """
        Creates the SPOINTs card that contains many SPOINTs

        Parameters
        ----------
        ids : List[int]
            SPOINT ids
        comment : str
            a comment for the card
        """
        XPoints.__init__(self, ids, comment=comment)

    def create_spointi(self):
        """
        Creates individal SPOINT objects
        """
        spoints = []
        for nid in self.points:
            spoints.append(SPOINT(nid))
        return spoints


class EPOINTs(XPoints):
    """
    +--------+-----+------+-----+-----+-----+-----+-----+-----+
    |   1    |  2  |  3   |  4  |  5  |  6  |  7  |  8  |  9  |
    +========+=====+======+=====+=====+=====+=====+=====+=====+
    | EPOINT | ID1 | THRU | ID2 |     |     |     |     |     |
    +--------+-----+------+-----+-----+-----+-----+-----+-----+
    | EPOINT | ID1 | ID1  | ID3 | ID4 | ID5 | ID6 | ID7 | ID8 |
    +--------+-----+------+-----+-----+-----+-----+-----+-----+
    |        | ID8 | etc. |     |     |     |     |     |     |
    +--------+-----+------+-----+-----+-----+-----+-----+-----+
    """
    type = 'EPOINT'

    def __init__(self, ids, comment=''):
        """
        Creates the EPOINTs card that contains many EPOINTs

        Parameters
        ----------
        ids : List[int]
            EPOINT ids
        comment : str
            a comment for the card
        """
        XPoints.__init__(self, ids, comment=comment)

    def create_epointi(self):
        """
        Creates individal EPOINT objects
        """
        points = []
        for nid in self.points:
            points.append(EPOINT(nid))
        return points


class GRDSET(Node):
    """
    Defines default options for fields 3, 7, 8, and 9 of all GRID entries.

    +--------+-----+----+----+----+----+----+----+------+
    |    1   |  2  | 3  | 4  | 5  | 6  |  7 | 8  |  9   |
    +========+=====+====+====+====+====+====+====+======+
    | GRDSET |     | CP |    |    |    | CD | PS | SEID |
    +--------+-----+----+----+----+----+----+----+------+
    """
    type = 'GRDSET'

    #: allows the get_field method and update_field methods to be used
    _field_map = {1: 'nid', 2:'cp', 6:'cd', 7:'ps', 8:'seid'}

    def __init__(self, cp, cd, ps, seid, comment=''):
        """
        Creates the GRDSET card
        """
        if comment:
            self._comment = comment

        #: Output Coordinate System
        self.cp = cp

        #: Analysis coordinate system
        self.cd = cd

        #: Default SPC constraint on undefined nodes
        self.ps = ps

        #: Superelement ID
        self.seid = seid

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Parameters
        ----------
        card : BDFCard()
           a BDFCard object
        comment : str
           a comment for the card
        """
        #: Grid point coordinate system
        blank(card, 1, 'blank')

        cp = integer_or_blank(card, 2, 'cp', 0)
        blank(card, 3, 'blank')
        blank(card, 4, 'blank')
        blank(card, 5, 'blank')
        cd = integer_or_blank(card, 6, 'cd', 0)

        ps = str(integer_or_blank(card, 7, 'ps', ''))
        seid = integer_or_blank(card, 8, 'seid', 0)
        assert len(card) <= 9, 'len(GRDSET card) = %i\ncard=%s' % (len(card), card)
        return GRDSET(cp, cd, ps, seid, comment=comment)

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ' which is required by the GRDSET'
        self.cp = model.Coord(self.cp, msg=msg)
        self.cp_ref = self.cp
        self.cd = model.Coord(self.cd, msg=msg)
        self.cd_ref = self.cd
        #self.seid = model.SuperElement(self.seid, msg)
        #self.seid_ref = self.seid

    def uncross_reference(self):
        self.cp = self.Cp()
        self.cd = self.Cd()
        del self.cp_ref, self.cd_ref

    def Cd(self):
        """
        Gets the output coordinate system

        Returns
        -------
        cd : int
            the output coordinate system
        """
        if isinstance(self.cd, integer_types):
            return self.cd
        else:
            return self.cd.cid

    def Cp(self):
        """
        Gets the analysis coordinate system

        Returns
        -------
        cp : int
            the analysis coordinate system
        """
        if isinstance(self.cp, integer_types):
            return self.cp
        else:
            return self.cp.cid

    def Ps(self):
        """
        Gets the GRID-based SPC

        Returns
        -------
        ps : str
            the GRID-based SPC
        """
        return self.ps

    def SEid(self):
        """
        Gets the Superelement ID

        Returns
        -------
        seid : int
            the Superelement ID
        """
        if isinstance(self.seid, integer_types):
            return self.seid
        else:
            return self.seid.seid

    def _verify(self, xref):
        """
        Verifies all methods for this object work

        :param xref: has this model been cross referenced
        :type xref:  bool
        """
        cp = self.Cp()
        seid = self.SEid()
        cd = self.Cd()
        ps = self.Ps()
        assert isinstance(cp, int), 'cp=%r' % cp
        assert isinstance(cd, int), 'cd=%r' % cd
        assert isinstance(ps, string_types), 'ps=%r' % ps
        assert isinstance(seid, int), 'seid=%r' % seid

    def raw_fields(self):
        """
        Gets the fields in their unmodified form

        :returns fields:
          the fields that define the card
        :type fields:
          LIST
        """
        list_fields = ['GRDSET', None, self.Cp(), None, None, None,
                       self.Cd(), self.ps, self.SEid()]
        return list_fields

    def repr_fields(self):
        """
        Gets the fields in their simplified form

        :returns fields:
          the fields that define the card
        :type fields:
          LIST
        """
        cp = set_blank_if_default(self.Cp(), 0)
        cd = set_blank_if_default(self.Cd(), 0)
        ps = set_blank_if_default(self.ps, 0)
        seid = set_blank_if_default(self.SEid(), 0)
        list_fields = ['GRDSET', None, cp, None, None, None, cd, ps, seid]
        return list_fields

    def write_card(self, f, size, is_double):
        """
        The writer method used by BDF.write_card

        :param size:
          the size of the card (8/16)
        :type size:
          int
        """
        card = self.repr_fields()
        f.write(print_card_8(card))


class GRIDB(Node):
    """defines the GRIDB class"""
    type = 'GRIDB'

    #: allows the get_field method and update_field methods to be used
    _field_map = {1: 'nid', 4:'phi', 6:'cd', 7:'ps', 8:'idf'}

    def __init__(self, nid, phi, cd, ps, idf, comment=''):
        """
        Creates the GRIDB card
        """
        if comment:
            self._comment = comment
        Node.__init__(self)

        #: node ID
        self.nid = nid
        self.phi = phi
        # analysis coordinate system
        self.cd = cd
        #: local SPC constraint
        self.ps = ps
        self.idf = idf

        assert self.nid > 0, 'nid=%s' % self.nid
        assert self.phi >= 0, 'phi=%s' % self.phi
        assert self.cd >= 0, 'cd=%s' % self.cd
        assert self.ps >= 0, 'ps=%s' % self.ps
        assert self.idf >= 0, 'idf=%s' % self.idf

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Parameters
        ----------
        card : BDFCard()
           the BDFCard object
        comment : str
          a comment for the card
        """
        nid = integer(card, 1, 'nid')
        phi = double(card, 4, 'phi')
        cd = integer(card, 6, 'cd')
        ps = integer(card, 7, 'ps')
        idf = integer(card, 8, 'idf')
        return GRIDB(nid, phi, cd, ps, idf, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        """
        Parameters
        ----------
        data : List[varies]
            a list with the GRIDB fields defined in OP2 format
        """
        nid = data[0]
        phi = data[1]
        cd = data[2]
        ps = data[3]
        idf = data[4]
        return GRIDB(nid, phi, cd, ps, idf, comment=comment)

    def _verify(self, xref):
        """
        Verifies all methods for this object work

        Returns
        -------
        xref : bool
            has this model been cross referenced
        """
        pass

    def Cd(self):
        """
        Gets the output coordinate system

        Returns
        -------
        cd : int
            the output coordinate system
        """
        if isinstance(self.cd, integer_types):
            return self.cd
        else:
            return self.cd.cid

    def raw_fields(self):
        """
        Gets the fields in their unmodified form

        :returns fields:
          the fields that define the card
        :type fields:
          LIST
        """
        list_fields = ['GRIDB', self.nid, None, None, self.phi, None,
                       self.Cd(), self.ps, self.idf]
        return list_fields

    def repr_fields(self):
        """
        Gets the fields in their simplified form

        :returns fields:
          the fields that define the card
        :type fields:
          LIST
        """
        #phi = set_blank_if_default(self.phi, 0.0)
        cd = set_blank_if_default(self.Cd(), 0)
        ps = set_blank_if_default(self.ps, 0)
        idf = set_blank_if_default(self.idf, 0)
        list_fields = ['GRIDB', self.nid, None, None, self.phi, None, cd, ps,
                       idf]
        return list_fields

    def write_card(self, size=8, is_double=False):
        """
        The writer method used by BDF.write_card

        Parameters
        ----------
        size : int; default=8
            the size of the card (8/16)
        is_double : bool; default=False
            should this card be written with double precision

        Returns
        -------
        msg : str
            the card as a string
        """
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)


class GRID(Node):
    """
    +------+-----+----+----+----+----+----+----+------+
    |   1  |  2  | 3  | 4  | 5  | 6  |  7 | 8  |  9   |
    +======+=====+====+====+====+====+====+====+======+
    | GRID | NID | CP | X1 | X2 | X3 | CD | PS | SEID |
    +------+-----+----+----+----+----+----+----+------+
    """
    type = 'GRID'

    #: allows the get_field method and update_field methods to be used
    _field_map = {1: 'nid', 2:'cp', 6:'cd', 7:'ps', 8:'seid'}

    def _get_field_helper(self, n):
        """
        Gets complicated parameters on the GRID card

        Parameters
        ----------
        n : int
            the field number to update

        Returns
        -------
        value : float
            the value for the appropriate field
        """
        if n == 3:
            value = self.xyz[0]
        elif n == 4:
            value = self.xyz[1]
        elif n == 5:
            value = self.xyz[2]
        else:
            raise KeyError('Field %r is an invalid %s entry.' % (n, self.type))
        return value

    def _update_field_helper(self, n, value):
        """
        Updates complicated parameters on the GRID card

        Parameters
        ----------
        n : int
            the field number to update
        value : float
            the value for the appropriate field
        """
        if n == 3:
            self.xyz[0] = value
        elif n == 4:
            self.xyz[1] = value
        elif n == 5:
            self.xyz[2] = value
        else:
            raise KeyError('Field %r=%r is an invalid %s entry.' % (n, value, self.type))

    def __init__(self, nid, cp=0, xyz=None, cd=0, ps='', seid=0, comment=''):
        """
        Creates the GRID card in a functional way

        Parameters
        ----------
        nid : int
            node id
        cp : int; default=0
            the xyz coordinate frame
        xyz : (3, ) float ndarray; default=None -> [0., 0., 0.]
            the xyz/r-theta-z/rho-theta-phi values
        cd : int; default=0
            the analysis coordinate frame
        ps : str; default=''
            Additional SPCs in the analysis coordinate frame (e.g. '123').
            This corresponds to DOF set ``SG``.
        seid : int; default=0
            ???
        """
        Node.__init__(self)
        if comment:
            self._comment = comment
        self.nid = nid
        self.cp = cp
        if xyz is None:
            xyz = [0., 0., 0.]
        self.xyz = np.asarray(xyz, dtype='float64')
        assert self.xyz.size == 3, self.xyz.shape
        self.cd = cd
        self.ps = ps
        self.seid = seid
        self._validate_input()

    @classmethod
    def add_op2_data(cls, data, comment=''):
        nid = data[0]
        cp = data[1]
        xyz = np.array(data[2:5])
        cd = data[5]
        ps = data[6]
        seid = data[7]
        if ps == 0:
            ps = ''
        return GRID(nid, cp, xyz, cd, ps, seid, comment=comment)

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Creates the GRID card

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        data : List[int/float]; default=None
            a list with the GRID fields defined in OP2 format
        comment : str; default=''
          a comment for the card
        """
        nfields = len(card)
        #: Node ID
        nid = integer(card, 1, 'nid')

        #: Grid point coordinate system
        cp = integer_or_blank(card, 2, 'cp', 0)

        #: node location in local frame
        xyz = np.array([
            double_or_blank(card, 3, 'x1', 0.),
            double_or_blank(card, 4, 'x2', 0.),
            double_or_blank(card, 5, 'x3', 0.)], dtype='float64')

        if nfields > 6:
            #: Analysis coordinate system
            cd = integer_or_blank(card, 6, 'cd', 0)

            #: SPC constraint
            ps = u(integer_or_blank(card, 7, 'ps', ''))

            #: Superelement ID
            seid = integer_or_blank(card, 8, 'seid', 0)
            assert len(card) <= 9, 'len(GRID card) = %i\ncard=%s' % (len(card), card)
        else:
            cd = 0
            ps = ''
            seid = 0
        return GRID(nid, cp, xyz, cd, ps, seid, comment=comment)

    def _validate_input(self):
        assert self.nid > 0, 'nid=%s' % (self.nid)
        assert self.cp >= 0, 'cp=%s' % (self.cp)
        assert self.cd >= -1, 'cd=%s' % (self.cd)
        assert self.seid >= 0, 'seid=%s' % (self.seid)
        assert len(self.xyz) == 3

    def Position(self, debug=False):
        self.deprecated('Position()', 'get_position()', '0.8')
        return self.get_position()

    def PositionWRT(self, model, cid):
        self.deprecated('PositionWRT(self, model, cid)',
                        'get_position_wrt(model, cid)', '0.8')
        return self.get_position_wrt(model, cid)

    def UpdatePosition(self, model, xyz, cid=0):
        self.deprecated('UpdatePosition(self, model, xyz, cid',
                        'set_position(self, model, xyz, cid=cid)', '0.8')
        return self.set_position(model, xyz, cid=cid)

    def Nid(self):
        """
        Gets the GRID ID

        Returns
        -------
        nid : int
            node ID
        """
        return self.nid

    def Ps(self):
        """
        Gets the GRID-based SPC

        Returns
        -------
        ps : int
            the GRID-based SPC
        """
        return self.ps

    def Cd(self):
        """
        Gets the output coordinate system

        Returns
        -------
        cd : int
            the output coordinate system
        """
        if isinstance(self.cd, integer_types):
            return self.cd
        else:
            return self.cd.cid

    def Cp(self):
        """
        Gets the analysis coordinate system

        Returns
        -------
        cp : int
            the analysis coordinate system
        """
        if isinstance(self.cp, integer_types):
            return self.cp
        else:
            return self.cp.cid

    def SEid(self):
        """
        Gets the Superelement ID

        Returns
        -------
        seid : int
            the Superelement ID
        """
        if isinstance(self.seid, integer_types):
            return self.seid
        else:
            return self.seid.seid

    def _verify(self, xref):
        """
        Verifies all methods for this object work

        Parameters
        ----------
        xref : bool
            has this model been cross referenced
        """
        nid = self.Nid()
        cp = self.Cp()
        cd = self.Cd()
        xyz = self.xyz
        ps = self.Ps()
        seid = self.SEid()
        assert isinstance(nid, int), 'nid=%r' % nid
        assert isinstance(cp, int), 'cp=%r' % cp
        assert isinstance(cd, int), 'cd=%r' % cd
        assert isinstance(ps, string_types), 'ps=%r' % ps
        assert isinstance(seid, int), 'seid=%r' % seid
        if xref:
            pos_xyz = self.get_position()

    def get_ndof(self):
        """
        Gets the number of degrees of freedom for the GRID

        Returns
        -------
        six : int
            the value 6
        """
        return 6

    def set_position(self, model, xyz, cid=0):
        """
        Updates the GRID location

        Parameters
        ----------
        xyz : (3, ) float ndarray
            the location of the node.
        cp : int; default=0 (global)
            the analysis coordinate system
        """
        self.xyz = xyz
        msg = ' which is required by GRID nid=%s' % self.nid
        self.cp = model.Coord(cid, msg=msg)

    def get_position(self):
        """
        Gets the point in the global XYZ coordinate system.

        Returns
        -------
        xyz : (3, ) float ndarray
            the position of the GRID in the global coordinate system
        """
        try:
            xyz = self.cp_ref.transform_node_to_global(self.xyz)
        except AttributeError:
            if self.cp == 0:
                return self.xyz
            raise
        return xyz

    def get_position_wrt(self, model, cid):
        """
        Gets the location of the GRID which started in some arbitrary
        system and returns it in the desired coordinate system

        Parameters
        ----------
        model : BDF()
            the BDF object
        cid : int
            the desired coordinate ID

        Returns
        -------
        xyz : (3, ) float ndarray
            the position of the GRID in an arbitrary coordinate system
        """
        if cid == self.Cp(): # same coordinate system
            return self.xyz

        # converting the xyz point arbitrary->global
        p = self.cp_ref.transform_node_to_global(self.xyz)

        # a matrix global->local matrix is found
        msg = ' which is required by %s nid=%s' % (self.type, self.nid)
        coord_b = model.Coord(cid, msg=msg)
        xyz = coord_b.transform_node_to_local(p)
        return xyz

    def cross_reference(self, model, grdset=None):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        grdset : GRDSET / None; default=None
            a GRDSET if available (default=None)

        .. note::  The gridset object will only update the fields that
                   have not been set
        """
        if grdset:
            # update using a gridset object
            if not self.cp:
                self.cp = grdset.cp
                self.cp_ref = self.cp
            if not self.cd:
                self.cd = grdset.cd
                self.cd_ref = self.cd
            if not self.ps:
                self.ps = grdset.ps
                self.ps_ref = self.ps
            if not self.seid:
                self.seid = grdset.seid
                self.seid_ref = self.seid
        msg = ' which is required by %s nid=%s' % (self.type, self.nid)
        self.cp = model.Coord(self.cp, msg=msg)
        self.cp_ref = self.cp
        if self.cd != -1:
            self.cd = model.Coord(self.cd, msg=msg)
            self.cd_ref = self.cd

    def uncross_reference(self):
        self.cp = self.Cp()
        self.cd = self.Cd()
        del self.cp_ref, self.cd_ref
        if hasattr(self, 'elements'):
            del self.elements

    def raw_fields(self):
        """
        Gets the fields in their unmodified form

        Returns
        -------
        fields : List[int/float/str]
            the fields that define the card
        """
        list_fields = ['GRID', self.nid, self.Cp()] + list(self.xyz) + \
                      [self.Cd(), self.ps, self.SEid()]
        return list_fields

    def repr_fields(self):
        """
        Gets the fields in their simplified form

        Returns
        -------
        fields : List[int/float/str]
            the fields that define the card
        """
        cp = set_blank_if_default(self.Cp(), 0)
        cd = set_blank_if_default(self.Cd(), 0)
        seid = set_blank_if_default(self.SEid(), 0)
        list_fields = ['GRID', self.nid, cp] + list(self.xyz) + [cd, self.ps,
                                                                 seid]
        return list_fields

    def write_card(self, size=8, is_double=False):
        """
        The writer method used by BDF.write_card

        Parameters
        ----------
        size : int; default=8
            the size of the card (8/16)
        is_double : bool; default=False
            should this card be written with double precision

        Returns
        -------
        msg : str
            the card as a string
        """
        if size == 8:
            return self.write_card_8()
        else:
            return self.write_card_16(is_double)

    def write_card_8(self):
        """
        Writes a GRID card in 8-field format
        """
        xyz = self.xyz
        cp = self.Cp()
        cd = self.Cd()

        cps = set_string8_blank_if_default(cp, 0)
        if [cd, self.ps, self.seid] == [0, '', 0]:
            # default
            msg = 'GRID    %8i%8s%s%s%s\n' % (
                self.nid, cps,
                print_float_8(xyz[0]),
                print_float_8(xyz[1]),
                print_float_8(xyz[2]))
            return self.comment + msg
        else:
            cds = set_string8_blank_if_default(cd, 0)
            seid = set_string8_blank_if_default(self.SEid(), 0)
            msg = 'GRID    %8i%8s%s%s%s%s%8s%s\n' % (
                self.nid, cps,
                print_float_8(xyz[0]),
                print_float_8(xyz[1]),
                print_float_8(xyz[2]),
                cds, self.ps, seid)
            return self.comment + msg

    def write_card_16(self, is_double=False):
        """
        Writes a GRID card in 16-field format
        """
        xyz = self.xyz
        cp = set_string16_blank_if_default(self.Cp(), 0)
        cd = set_string16_blank_if_default(self.Cd(), 0)
        seid = set_string16_blank_if_default(self.SEid(), 0)

        if is_double:
            if [cd, self.ps, self.seid] == [0, '', 0]:
                msg = ('GRID*   %16i%16s%16s%16s\n'
                       '*       %16s\n' % (
                           self.nid,
                           cp,
                           print_scientific_double(xyz[0]),
                           print_scientific_double(xyz[1]),
                           print_scientific_double(xyz[2])))
            else:
                msg = ('GRID*   %16i%16s%16s%16s\n'
                       '*       %16s%16s%16s%16s\n' % (
                           self.nid,
                           cp,
                           print_scientific_double(xyz[0]),
                           print_scientific_double(xyz[1]),
                           print_scientific_double(xyz[2]),
                           cd, self.ps, seid))
        else:
            if [cd, self.ps, self.seid] == [0, '', 0]:
                msg = ('GRID*   %16i%16s%16s%16s\n'
                       '*       %16s\n' % (
                           self.nid,
                           cp,
                           print_float_16(xyz[0]),
                           print_float_16(xyz[1]),
                           print_float_16(xyz[2])))
            else:
                msg = ('GRID*   %16i%16s%16s%16s\n'
                       '*       %16s%16s%16s%16s\n' % (
                           self.nid,
                           cp,
                           print_float_16(xyz[0]),
                           print_float_16(xyz[1]),
                           print_float_16(xyz[2]),
                           cd, self.ps, seid))
        return self.comment + msg


class POINT(Node):
    """
    +-------+-----+----+----+----+----+----+----+-----+
    |   1   |  2  | 3  | 4  | 5  | 6  |  7 | 8  |  9  |
    +=======+=====+====+====+====+====+====+====+=====+
    | POINT | NID | CP | X1 | X2 | X3 |    |    |     |
    +-------+-----+----+----+----+----+----+----+-----+
    """
    type = 'POINT'
    _field_map = {1: 'nid', 2:'cp'}

    def _get_field_helper(self, n):
        """
        Gets complicated parameters on the POINT card

        Parameters
        ----------
        n : int
            the field number to update

        Returns
        -------
        value : varies
            the value for the appropriate field
        """
        if n == 3:
            value = self.xyz[0]
        elif n == 4:
            value = self.xyz[1]
        elif n == 5:
            value = self.xyz[2]
        else:
            raise KeyError('Field %r is an invalid %s entry.' % (n, self.type))
        return value

    def _update_field_helper(self, n, value):
        """
        Updates complicated parameters on the POINT card

        Parameters
        ----------
        n : int
            the field number to update
        value : varies
            the value for the appropriate field
        """
        if n == 3:
            self.xyz[0] = value
        elif n == 4:
            self.xyz[1] = value
        elif n == 5:
            self.xyz[2] = value
        else:
            raise KeyError('Field %r=%r is an invalid %s entry.' % (n, value, self.type))

    def __init__(self, nid, cp, xyz, cd, ps, seid, comment=''):
        """
        Creates the POINT card
        """
        if comment:
            self._comment = comment
        Node.__init__(self)

        #: Node ID
        self.nid = nid

        #: Grid point coordinate system
        self.cp = cp

        #: node location in local frame
        self.xyz = xyz

        #: Analysis coordinate system
        self.cd = cd

        #: SPC constraint
        self.ps = ps

        #: Superelement ID
        self.seid = seid
        assert self.nid > 0, 'nid=%s' % (self.nid)
        assert self.cp >= 0, 'cp=%s' % (self.cp)
        assert len(xyz) == 3

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Parameters
        ----------
        card : BDFCard()
           the BDFCard object
        comment : str
          a comment for the card
        """
        nid = integer(card, 1, 'nid')
        cp = integer_or_blank(card, 2, 'cp', 0)

        xyz = np.array([
            double_or_blank(card, 3, 'x1', 0.),
            double_or_blank(card, 4, 'x2', 0.),
            double_or_blank(card, 5, 'x3', 0.)], dtype='float64')

        cd = blank(card, 6, 'cd', 0)
        ps = blank(card, 7, 'ps', '')

        seid = blank(card, 8, 'seid', 0)
        assert len(card) <= 9, 'len(POINT card) = %i\ncard=%s' % (len(card), card)
        return POINT(nid, cp, xyz, cd, ps, seid, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        nid = data[0]
        cp = data[1]
        xyz = np.array(data[2:5])
        ps = ''
        seid = 0
        cd = 0
        return POINT(nid, cp, xyz, cd, ps, seid, comment=comment)

    def Position(self):
        self.deprecated('Position()', 'get_position()', '0.8')
        return self.get_position()

    def PositionWRT(self, model, cid):
        self.deprecated('Position()', 'get_position_wrt()', '0.8')
        return self.get_position_wrt(model, cid)

    def UpdatePosition(self, model, xyz, cid=0):
        self.deprecated('UpdatePosition(self, model, xyz, cid)',
                        'set_position(model, xyz, cid)', '0.8')
        return self.set_position(model, xyz, cid=cid)

    def set_position(self, model, xyz, cid=0):
        """
        Updates the POINT location

        Parameters
        ----------
        xyz : (3,) float ndarray
            the location of the node
        cp : int; default=0 (global)
            the analysis coordinate system
        """
        self.xyz = xyz
        msg = ' which is required by POINT nid=%s' % self.nid
        self.cp = model.Coord(cid, msg=msg)

    def get_position(self):
        """
        Gets the point in the global XYZ coordinate system.

        Returns
        -------
        position : (3,) float ndarray
            the position of the POINT in the globaly coordinate system
        """
        p = self.cp.transform_node_to_global(self.xyz)
        return p

    def get_position_wrt(self, model, cid):
        """
        Gets the location of the POINT which started in some arbitrary
        system and returns it in the desired coordinate system

        Parameters
        ----------
        model : BDF()
            the BDF model object
        cid : int
            the desired coordinate ID

        Returns
        -------
        xyz : (3,) ndarray
            the position of the POINT in an arbitrary coordinate system
        """
        if cid == self.Cp(): # same coordinate system
            return self.xyz

        # converting the xyz point arbitrary->global
        p = self.cp.transform_node_to_global(self.xyz)

        # a matrix global->local matrix is found
        msg = ' which is required by %s nid=%s' % (self.type, self.nid)
        coord_b = model.Coord(cid, msg=msg)
        xyz = coord_b.transform_node_to_local(p)
        return xyz

    def Cp(self):
        """
        Gets the analysis coordinate system

        Returns
        -------
        cp : int
            the analysis coordinate system
        """
        if isinstance(self.cp, integer_types):
            return self.cp
        else:
            return self.cp_ref.cid

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        self.cp = model.Coord(self.cp)
        self.cp_ref = self.cp

    def uncross_reference(self):
        self.cp = self.Cp()
        del self.cp_ref

    def raw_fields(self):
        """
        Gets the fields in their unmodified form

        Returns
        -------
        fields : list[varies]
            the fields that define the card
        """
        list_fields = ['POINT', self.nid, self.Cp()] + list(self.xyz)
        return list_fields

    def repr_fields(self):
        """
        Gets the fields in their simplified form

        :returns fields:
          the fields that define the card
        :type fields:
          LIST
        """
        cp = set_blank_if_default(self.Cp(), 0)
        list_fields = ['POINT', self.nid, cp] + list(self.xyz)
        return list_fields

    def write_card(self, size=8, is_double=False):
        """
        The writer method used by BDF.write_card

        Parameters
        ----------
        size : int
            the size of the card (8/16)
        """
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        if is_double:
            return self.comment + print_card_double(card)
        return self.comment + print_card_16(card)
