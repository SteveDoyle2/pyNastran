# pylint: disable=C0103,R0902,R0904,R0914
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
from numpy import array

from pyNastran.bdf.deprecated import GridDeprecated, PointDeprecated, SPOINTsDeprecated

from pyNastran.bdf.field_writer_8 import set_string8_blank_if_default
from pyNastran.bdf.field_writer_16 import set_string16_blank_if_default

from pyNastran.bdf.field_writer_8 import set_blank_if_default
from pyNastran.bdf.cards.baseCard import BaseCard, expand_thru, collapse_thru_packs
from pyNastran.bdf.bdfInterface.assign_type import (integer, integer_or_blank,
    double, double_or_blank, blank, integer_or_string)
from pyNastran.bdf.field_writer_8 import print_card_8, print_float_8, print_int_card
from pyNastran.bdf.field_writer_16 import print_float_16, print_card_16
from pyNastran.bdf.field_writer_double import print_scientific_double, print_card_double

class Ring(BaseCard):
    """Generic Ring base class"""
    def __init__(self, card, data):
        assert card is None or data is None


class Node(BaseCard):
    """Generic Node base class"""
    def __init__(self, card, data):
        assert card is None or data is None


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
        :param self:
          the RINGAX object pointer
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
        Ring.__init__(self, card, data)
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
            assert len(card) <= 8, 'len(RINGAX card) = %i' % len(card)
        else:  # hasn't been validated
            self.nid = data[0]
            self.R = data[1]
            self.z = data[2]
            self.ps = data[3]
            assert len(data) == 4, data

    def raw_fields(self):
        """
        Gets the fields in their unmodified form

        :param self:
          the RINGAX object pointer
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

        :param self:
          the RINGAX object pointer
        :param size:
          the size of the card (8/16)
        :type size:
          int
        """
        card = self.repr_fields()
        if size == 8:
            return self.comment() + print_card_8(card)
        return self.comment() + print_card_16(card)


class SPOINT(Node):
    type = 'SPOINT'

    def __init__(self, nid, comment=''):
        """
        Creates the SPOINT card

        :param self:
          the SPOINT object pointer
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
        Node.__init__(self, card=None, data=None)
        if comment:
            self._comment = comment
        self.nid = nid
        assert isinstance(nid, int), nid

    def cross_reference(self, model):
        """
        Cross links the card

        :param self:
          the SPOINT object pointer
        :param model:
          the BDF object
        :type model:
          BDF
        """
        pass

    def raw_fields(self):
        """
        Gets the fields in their unmodified form

        :param self:
          the SPOINT object pointer
        :returns fields:
          the fields that define the card
        :type fields:
          LIST
        """
        lists_fields = []
        if isinstance(self.nid, int):
            list_fields = ['SPOINT', self.nid]
            lists_fields.append(list_fields)
        else:
            singles, doubles = collapse_thru_packs(self.nid)
            if singles:
                list_fields = ['SPOINT'] + singles
            if doubles:
                for spoint_double in doubles:
                    list_fields = ['SPOINT'] + spoint_double
                    lists_fields.append(list_fields)
        return lists_fields

    def write_card(self, size=8, is_double=False):
        """
        The writer method used by BDF.write_card

        :param self:   the SPOINT object pointer
        :param size:   unused
        :param is_double: unused
        """
        lists_fields = self.repr_fields()
        msg = self.comment()
        for list_fields in lists_fields:
            if 'THRU' not in list_fields:
                msg += print_int_card(list_fields)
            else:
                msg += print_card_8(list_fields)
        return msg


class SPOINTs(Node):
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

    def __init__(self, card=None, data=None, comment=''):
        """
        Creates the SPOINTs card that contains many SPOINTs
        :param self:
          the SPOINTs object pointer
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
        if comment:
            self._comment = comment
        Node.__init__(self, card, data)

        if card:
            fields = []
            for i in range(1, len(card)):
                field = integer_or_string(card, i, 'ID%i' % i)
                fields.append(field)
        else:
            fields = data
            assert isinstance(data, list), data
            assert isinstance(data[0], int), data
        self.spoints = set(expand_thru(fields))

    def get_ndof(self):
        """
        Returns the number of degrees of freedom for the SPOINTs class

        :param self:
          the SPOINT object pointer
        :returns ndofs:
          the number of degrees of freedom
        :type ndofs:
          int
        """
        return len(self.spoints)

    def addSPoints(self, sList):
        """
        Adds more SPOINTs to this object

        :param self:
          the SPOINT object pointer
        """
        self.spoints = self.spoints.union(set(sList))

    def cross_reference(self, model):
        """
        Cross links the card

        :param self:
          the SPOINT object pointer
        :param model:
          the BDF object
        :type model:
          BDF
        """
        pass

    def createSPOINTi(self):
        """
        Creates individal SPOINT objects

        :param self:
          the SPOINT object pointer
        """
        spoints = []
        for nid in self.spoints:
            spoints.append(SPOINT(nid))
        return spoints

    def raw_fields(self):
        """
        Gets the fields in their unmodified form

        :param self:
          the SPOINT object pointer
        :returns fields:
          the fields that define the card
        :type fields:
          LIST
        """
        spoints = list(self.spoints)
        spoints.sort()
        return ['SPOINT'] + spoints

    def _get_compressed_spoints(self):
        """
        Gets the spoints in sorted, short form.

          uncompresed:  SPOINT,1,3,5
          compressed:   SPOINT,1,3,5

          uncompresed:  SPOINT,1,2,3,4,5
          compressed:   SPOINT,1,THRU,5

          uncompresed:  SPOINT,1,2,3,4,5,7
          compressed:   SPOINT,7
                        SPOINT,1,THRU,5
        """
        spoints = list(self.spoints)
        spoints.sort()

        singles, doubles = collapse_thru_packs(spoints)

        lists_fields = []
        if singles:
            list_fields = ['SPOINT'] + singles
            lists_fields.append(list_fields)
        if doubles:
            for spoint_double in doubles:
                list_fields = ['SPOINT'] + spoint_double
                lists_fields.append(list_fields)
        return lists_fields

    def write_card(self, size=8, is_double=False):
        """
        The writer method used by BDF.write_card

        :param self:   the SPOINT object pointer
        :param size:   unused
        :param is_double: unused
        """
        lists_fields = self._get_compressed_spoints()
        msg = self.comment()
        for list_fields in lists_fields:
            if 'THRU' not in list_fields:
                msg += print_int_card(list_fields)
            else:
                msg += print_card_8(list_fields)
        return msg


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

    def __init__(self, card=None, data=None, comment=''):
        """
        Creates the GRDSET card

        :param self:
          the GRDSET object pointer
        :param card:
          a BDFCard object
        :type card:
          BDFCard
        :param data:
          a list with the GRDSET fields defined in OP2 format
        :type data:
          LIST
        :param comment:
          a comment for the card
        :type comment:
          string
        """
        if comment:
            self._comment = comment

        #: Grid point coordinate system
        blank(card, 1, 'blank')

        #: Output Coordinate System
        self.cp = integer_or_blank(card, 2, 'cp', 0)
        blank(card, 3, 'blank')
        blank(card, 4, 'blank')
        blank(card, 5, 'blank')

        #: Analysis coordinate system
        self.cd = integer_or_blank(card, 6, 'cd', 0)

        #: Default SPC constraint on undefined nodes
        self.ps = str(integer_or_blank(card, 7, 'ps', ''))

        #: Superelement ID
        self.seid = integer_or_blank(card, 8, 'seid', 0)
        assert len(card) <= 9, 'len(GRDSET card) = %i' % len(card)

    def cross_reference(self, model):
        """
        Cross links the card

        :param self:
          the SPOINT object pointer
        :param model:
          the BDF object
        :type model:
          BDF
        """
        msg = ' which is required by the GRDSET'
        self.cp = model.Coord(self.cp, msg=msg)
        self.cd = model.Coord(self.cd, msg=msg)
        #self.seid = model.SuperElement(self.seid, msg)

    def Cd(self):
        """
        Gets the output coordinate system

        :param self: the GRDSET object pointer
        :returns cd: the output coordinate system
        :type cd:    int
        """
        if isinstance(self.cd, int):
            return self.cd
        else:
            return self.cd.cid

    def Cp(self):
        """
        Gets the analysis coordinate system

        :param self: the GRDSET object pointer
        :returns cp: the analysis coordinate system
        :type cp:    int
        """
        if isinstance(self.cp, int):
            return self.cp
        else:
            return self.cp.cid

    def Ps(self):
        """
        Gets the GRID-based SPC

        :param self: the GRID object pointer
        :returns ps: the GRID-based SPC
        """
        return self.ps

    def SEid(self):
        """
        Gets the Superelement ID

        :param self:   the GRDSET object pointer
        :returns seid: the Superelement ID
        :type seid:    int
        """
        if isinstance(self.seid, int):
            return self.seid
        else:
            return self.seid.seid

    def _verify(self, xref):
        """
        Verifies all methods for this object work

        :param self: the GRDSET object pointer
        :param xref: has this model been cross referenced
        :type xref:  bool
        """
        cp = self.Cp()
        seid = self.SEid()
        cd = self.Cd()
        ps = self.Ps()
        assert isinstance(cp, int), 'cp=%r' % cp
        assert isinstance(cd, int), 'cd=%r' % cd
        assert isinstance(ps, str), 'ps=%r' % ps
        assert isinstance(seid, int), 'seid=%r' % seid

    def raw_fields(self):
        """
        Gets the fields in their unmodified form

        :param self:
          the GRDSET object pointer
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

        :param self:
          the GRDSET object pointer
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

        :param self:
          the SPOINT object pointer
        :param size:
          the size of the card (8/16)
        :type size:
          int
        """
        card = self.repr_fields()
        f.write(print_card_8(card))


class GRIDB(Node):
    type = 'GRIDB'

    #: allows the get_field method and update_field methods to be used
    _field_map = {1: 'nid', 4:'phi', 6:'cd', 7:'ps', 8:'idf'}

    def __init__(self, card=None, data=None, comment=''):
        """
        Creates the GRIDB card

        :param self:
          the GRIDB object pointer
        :param card:
          a BDFCard object
        :type card:
          BDFCard
        :param data:
          a list with the GRIDB fields defined in OP2 format
        :type data:
          LIST
        :param comment:
          a comment for the card
        :type comment:
          string
        """
        if comment:
            self._comment = comment
        Node.__init__(self, card, data)

        if card:
            #: node ID
            self.nid = integer(card, 1, 'nid')

            self.phi = double(card, 4, 'phi')

            # analysis coordinate system
            self.cd = integer(card, 6, 'cd')

            #: local SPC constraint
            self.ps = integer(card, 7, 'ps')

            self.idf = integer(card, 8, 'idf')
        else:
            self.nid = data[0]
            self.phi = data[1]
            self.cd = data[2]
            self.ps = data[3]
            self.idf = data[4]

        assert self.nid > 0, 'nid=%s' % self.nid
        assert self.phi >= 0, 'phi=%s' % self.phi
        assert self.cd >= 0, 'cd=%s' % self.cd
        assert self.ps >= 0, 'ps=%s' % self.ps
        assert self.idf >= 0, 'idf=%s' % self.idf

    def _verify(self, xref):
        """
        Verifies all methods for this object work

        :param self: the GRIDB object pointer
        :param xref: has this model been cross referenced
        :type xref:  bool
        """
        pass

    def Cd(self):
        """
        Gets the output coordinate system

        :param self: the GRIDB object pointer
        :returns cd: the output coordinate system
        :type cd:    int
        """
        if isinstance(self.cd, int):
            return self.cd
        else:
            return self.cd.cid

    def raw_fields(self):
        """
        Gets the fields in their unmodified form

        :param self:
          the GRIDB object pointer
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

        :param self:
          the GRIDB object pointer
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

        :param self:
          the GRIDB object pointer
        :param size:
          the size of the card (8/16)
        :type size:
          int
        """
        card = self.repr_fields()
        if size == 8:
            return self.comment() + print_card_8(card)
        return self.comment() + print_card_16(card)


class GRID(Node, GridDeprecated):
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

        :param self:  the GRID object pointer
        :param n:     the field number to update
        :type n:      int
        :param value: the value for the appropriate field
        :type field:  varies
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

        :param self:  the GRID object pointer
        :param n:     the field number to update
        :type n:      int
        :param value: the value for the appropriate field
        :type field:  varies
        """
        if n == 3:
            self.xyz[0] = value
        elif n == 4:
            self.xyz[1] = value
        elif n == 5:
            self.xyz[2] = value
        else:
            raise KeyError('Field %r=%r is an invalid %s entry.' % (n, value, self.type))

    def __init__(self, card=None, data=None, comment=''):
        """
        Creates the GRID card

        :param self:
          the GRID object pointer
        :param card:
          a BDFCard object
        :type card:
          BDFCard
        :param data:
          a list with the GRID fields defined in OP2 format
        :type data:
          LIST
        :param comment:
          a comment for the card
        :type comment:
          string
        """
        if comment:
            self._comment = comment
        Node.__init__(self, card, data)

        if card:
            #: Node ID
            self.nid = integer(card, 1, 'nid')

            #: Grid point coordinate system
            self.cp = integer_or_blank(card, 2, 'cp', 0)

            #: node location in local frame
            self.xyz = array([
                double_or_blank(card, 3, 'x1', 0.),
                double_or_blank(card, 4, 'x2', 0.),
                double_or_blank(card, 5, 'x3', 0.)], dtype='float64')

            #: Analysis coordinate system
            self.cd = integer_or_blank(card, 6, 'cd', 0)

            #: SPC constraint
            self.ps = str(integer_or_blank(card, 7, 'ps', ''))

            #: Superelement ID
            self.seid = integer_or_blank(card, 8, 'seid', 0)
            assert len(card) <= 9, 'len(GRID card) = %i' % len(card)
        else:
            self.nid = data[0]
            self.cp = data[1]
            self.xyz = array(data[2:5])
            self.cd = data[5]
            self.ps = data[6]
            self.seid = data[7]
            if self.ps == 0:
                self.ps = ''
            assert len(self.xyz) == 3

        assert self.nid > 0, 'nid=%s' % (self.nid)
        assert self.cp >= 0, 'cp=%s' % (self.cp)
        assert self.cd >= -1, 'cd=%s' % (self.cd)
        assert self.seid >= 0, 'seid=%s' % (self.seid)

    def Nid(self):
        """
        Gets the GRID ID

        :param self: the GRID object pointer
        :returns nid: node ID
        :type nid:    int
        """
        return self.nid

    def Ps(self):
        """
        Gets the GRID-based SPC

        :param self: the GRID object pointer
        :returns ps: the GRID-based SPC
        """
        return self.ps

    def Cd(self):
        """
        Gets the output coordinate system

        :param self: the GRID object pointer

        :returns cd: the output coordinate system
        :type cd:    int
        """
        if isinstance(self.cd, int):
            return self.cd
        else:
            return self.cd.cid

    def Cp(self):
        """
        Gets the analysis coordinate system

        :param self: the GRID object pointer
        :returns cp: the analysis coordinate system
        :type cp:    int
        """
        if isinstance(self.cp, int):
            return self.cp
        else:
            return self.cp.cid

    def SEid(self):
        """
        Gets the Superelement ID

        :param self:   the GRID object pointer
        :returns seid: the Superelement ID
        :type seid:    int
        """
        if isinstance(self.seid, int):
            return self.seid
        else:
            return self.seid.seid

    def _verify(self, xref):
        """
        Verifies all methods for this object work

        :param self: the GRID object pointer
        :param xref: has this model been cross referenced
        :type xref:  bool
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
        assert isinstance(ps, str), 'ps=%r' % ps
        assert isinstance(seid, int), 'seid=%r' % seid
        if xref:
            pos_xyz = self.Position()

    def get_ndof(self):
        """
        Gets the number of degrees of freedom for the GRID

        :param self:  the GRID object pointer
        :returns six: the value 6
        :type six:    int
        """
        return 6

    def set_position(self, model, xyz, cid=0):
        """
        Updates the GRID location

        :param self: the GRID object pointer
        :param xyz:  the location of the node.
        :type xyz:   TYPE = NDARRAY.  SIZE=(3,)
        :param cp:   the analysis coordinate system.  (default=0; global)
        :type cp:    int
        """
        self.xyz = xyz
        msg = ' which is required by GRID nid=%s' % self.nid
        self.cp = model.Coord(cid, msg=msg)

    def get_position(self, debug=False):
        """
        Gets the point in the global XYZ coordinate system.

        :param self:   the GRID object pointer
        :param debug:  developer debug (default=False)
        :type debug:   bool
        :returns xyz:  the position of the GRID in the globaly
                       coordinate system
        :type xyz:     TYPE = NDARRAY.  SIZE=(3,)
        """
        xyz = self.cp.transform_node_to_global(self.xyz, debug=debug)
        return xyz

    def get_position_wrt(self, model, cid, debug=False):
        """
        Gets the location of the GRID which started in some arbitrary
        system and returns it in the desired coordinate system

        :param self:  the object pointer
        :param model: the BDF model object
        :type model:  BDF()
        :param cid:   the desired coordinate ID
        :type cid:    int
        :param debug: developer debug (default=False)
        :type debug:  bool
        :returns xyz: the position of the GRID in an arbitrary
                      coordinate system
        :type xyz:    TYPE = NDARRAY.  SIZE=(3,)
        """
        if cid == self.Cp(): # same coordinate system
            return self.xyz

        # converting the xyz point arbitrary->global
        p = self.cp.transform_node_to_global(self.xyz, debug=debug)

        # a matrix global->local matrix is found
        msg = ' which is required by %s nid=%s' % (self.type, self.nid)
        coordB = model.Coord(cid, msg=msg)
        xyz = coordB.transform_node_to_local(p, debug=debug)
        return xyz

    def cross_reference(self, model, grdset=None):
        """
        Cross links the card

        :param self:   the GRID object pointer
        :param model:  the BDF object
        :type model:   BDF()
        :param grdset: a GRDSET if available (default=None)
        :type grdset:  GRDSET() or None

        .. note::  The gridset object will only update the fields that
                   have not been set
        """
        if grdset:  # update using a gridset object
            if not self.cp:
                self.cp = grdset.cp
            if not self.cd:
                self.cd = grdset.cd
            if not self.ps:
                self.ps = grdset.ps
            if not self.seid:
                self.seid = grdset.seid
        msg = ' which is required by %s nid=%s' % (self.type, self.nid)
        self.cp = model.Coord(self.cp, msg=msg)
        if self.cd != -1:
            self.cd = model.Coord(self.cd, msg=msg)

    def raw_fields(self):
        """
        Gets the fields in their unmodified form

        :param self:
          the GRID object pointer
        :returns fields:
          the fields that define the card
        :type fields:
          LIST
        """
        list_fields = ['GRID', self.nid, self.Cp()] + list(self.xyz) + \
                      [self.Cd(), self.ps, self.SEid()]
        return list_fields

    def repr_fields(self):
        """
        Gets the fields in their simplified form

        :param self:
          the GRID object pointer
        :returns fields:
          the fields that define the card
        :type fields:
          LIST
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

        :param self:
          the GRID object pointer
        :param size:
          the size of the card (8/16)
        :type size:
          int
        :param isdouble:
          should this card be written with double precision (default=False)
        :type is_double:
          bool
        """
        xyz = self.xyz
        if size == 8:
            cp = set_string8_blank_if_default(self.Cp(), 0)
            cd = set_string8_blank_if_default(self.Cd(), 0)
            seid = set_string8_blank_if_default(self.SEid(), 0)
            msg = '%-8s%8i%8s%s%s%s%s%8s%s\n' % (
                'GRID', self.nid, cp,
                print_float_8(xyz[0]),
                print_float_8(xyz[1]),
                print_float_8(xyz[2]),
                cd, self.ps, seid)
        else:
            cp = set_string16_blank_if_default(self.Cp(), 0)
            cd = set_string16_blank_if_default(self.Cd(), 0)
            seid = set_string16_blank_if_default(self.SEid(), 0)
            if is_double:
                msg = ('%-8s%16i%16s%16s%16s\n'
                       '%-8s%16s%16s%16s%16s\n' % (
                           'GRID*', self.nid,
                           cp,
                           print_scientific_double(xyz[0]),
                           print_scientific_double(xyz[1]),
                           '*',
                           print_scientific_double(xyz[2]),
                           cd, self.ps, seid))
            else:
                msg = ('%-8s%16i%16s%16s%16s\n'
                       '%-8s%16s%16s%16s%16s\n' % (
                           'GRID*', self.nid,
                           cp,
                           print_float_16(xyz[0]),
                           print_float_16(xyz[1]),
                           '*',
                           print_float_16(xyz[2]),
                           cd, self.ps, seid))
        return self.comment() + msg.rstrip() + '\n'


class POINT(Node, PointDeprecated):
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

        :param self:  the POINT object pointer
        :param n:     the field number to update
        :type n:      int
        :param value: the value for the appropriate field
        :type field:  varies
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

        :param self:  the POINT object pointer
        :param n:     the field number to update
        :type n:      int
        :param value: the value for the appropriate field
        :type field:  varies
        """
        if n == 3:
            self.xyz[0] = value
        elif n == 4:
            self.xyz[1] = value
        elif n == 5:
            self.xyz[2] = value
        else:
            raise KeyError('Field %r=%r is an invalid %s entry.' % (n, value, self.type))

    def __init__(self, card=None, data=None, comment=''):
        """
        Creates the POINT card

        :param self:
          the POINT object pointer
        :param card:
          a BDFCard object
        :type card:
          BDFCard
        :param data:
          a list with the POINT fields defined in OP2 format
        :type data:
          LIST
        :param comment:
          a comment for the card
        :type comment:
          string
        """
        if comment:
            self._comment = comment
        Node.__init__(self, card, data)

        if card:
            #: Node ID
            self.nid = integer(card, 1, 'nid')

            #: Grid point coordinate system
            self.cp = integer_or_blank(card, 2, 'cp', 0)

            #: node location in local frame
            self.xyz = array([
                double_or_blank(card, 3, 'x1', 0.),
                double_or_blank(card, 4, 'x2', 0.),
                double_or_blank(card, 5, 'x3', 0.)], dtype='float64')

            #: Analysis coordinate system
            self.cd = blank(card, 6, 'cd', 0)

            #: SPC constraint
            self.ps = blank(card, 7, 'ps', '')

            #: Superelement ID
            self.seid = blank(card, 8, 'seid', 0)
            assert len(card) <= 9, 'len(POINT card) = %i' % len(card)
        else:
            self.nid = data[0]
            self.cp = data[1]
            self.xyz = array(data[2:5])
            assert len(self.xyz) == 3
            self.ps = ''
            self.seid = 0
            self.cd = 0

        assert self.nid > 0, 'nid=%s' % (self.nid)
        assert self.cp >= 0, 'cp=%s' % (self.cp)

    def set_position(self, model, xyz, cid=0):
        """
        Updates the POINT location

        :param self: the POINT object pointer
        :param xyz:  the location of the node.
        :type xyz:   TYPE = NDARRAY.  SIZE=(3,)
        :param cp:   the analysis coordinate system.  (default=0; global)
        :type cp:    int
        """
        self.xyz = xyz
        msg = ' which is required by POINT nid=%s' % self.nid
        self.cp = model.Coord(cid, msg=msg)

    def get_position(self, debug=False):
        """
        Gets the point in the global XYZ coordinate system.

        :param self:  the POINT object pointer
        :param debug: developer debug (default=False)
        :returns position: the position of the POINT in the globaly
                           coordinate system
        """
        p = self.cp.transform_node_to_global(self.xyz, debug=debug)
        return p

    def get_position_wrt(self, model, cid, debug=False):
        """
        Gets the location of the POINT which started in some arbitrary
        system and returns it in the desired coordinate system

        :param self:  the POINT object pointer
        :param model: the BDF model object
        :type model:  BDF()
        :param cid:   the desired coordinate ID
        :type cid:    int
        :param debug: debug (default=False)
        :type debug:  bool
        :returns xyz: the position of the POINT in an arbitrary
                      coordinate system
        :type xyz:    TYPE = NDARRAY.  SIZE=(3,)
        """
        if cid == self.Cp(): # same coordinate system
            return self.xyz

        # converting the xyz point arbitrary->global
        p = self.cp.transform_node_to_global(self.xyz, debug=debug)

        # a matrix global->local matrix is found
        msg = ' which is required by %s nid=%s' % (self.type, self.nid)
        coordB = model.Coord(cid, msg=msg)
        xyz = coordB.transform_node_to_local(p, debug=debug)
        return xyz

    def Cp(self):
        """
        Gets the analysis coordinate system

        :param self: the POINT object pointer
        :returns cp: the analysis coordinate system
        :type cp:    int
        """
        if isinstance(self.cp, int):
            return self.cp
        else:
            return self.cp.cid

    def cross_reference(self, model):
        """
        Cross links the card

        :param self:   the GRID object pointer
        :param model:  the BDF object
        :type model:   BDF()
        """
        self.cp = model.Coord(self.cp)

    def raw_fields(self):
        """
        Gets the fields in their unmodified form

        :param self:
          the GRID object pointer
        :returns fields:
          the fields that define the card
        :type fields:
          LIST
        """
        list_fields = ['POINT', self.nid, self.Cp()] + list(self.xyz)
        return list_fields

    def repr_fields(self):
        """
        Gets the fields in their simplified form

        :param self:
          the GRID object pointer
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

        :param self:
          the GRID object pointer
        :param size:
          the size of the card (8/16)
        :type size:
          int
        """
        card = self.repr_fields()
        if size == 8:
            return self.comment() + print_card_8(card)
        if is_double:
            return self.comment() + print_card_double(card)
        return self.comment() + print_card_16(card)
