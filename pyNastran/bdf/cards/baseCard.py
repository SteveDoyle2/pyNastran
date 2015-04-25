# pylint: disable=R0904,R0902,C0111,C0103
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
from six import string_types, integer_types
from six.moves import zip, range

from pyNastran.bdf.fieldWriter import print_card
from pyNastran.bdf.field_writer_8 import is_same
from pyNastran.bdf.bdfInterface.assign_type import interpret_value
from pyNastran.bdf.deprecated import BaseCardDeprecated, ElementDeprecated


class BaseCard(BaseCardDeprecated):
    def __init__(self):
        pass

    def comment(self):
        if hasattr(self, '_comment'):
            return '%s' % self._comment
        return ''

    def _test_update_fields(self):
        n = 1
        while 1:
            try:
                self.update_field(n, 1.0)  # dummy updating the field
            except IndexError:
                return
            except KeyError:
                return

    def update_field(self, n, value):
        """
        Updates a field based on it's field number.

        :param self: the BaseCard object
        :param n: the integer field number
        :param value: the value to update the field to

        .. note::
            This is dynamic if the card length changes.

        .. code-block:: python

          nid = 1
          node = model.nodes[nid]
          # ['GRID', nid, cp, x, y, z]
          node.update_field(3, 0.1) # change the z coordinate
        """
        try:
            key_name = self._field_map[n]
            setattr(self, key_name, value)
        except KeyError:
            self._update_field_helper(n, value)

    def _update_field_helper(self, n, value):
        raise NotImplementedError()

    def _get_field_helper(self, n):
        raise NotImplementedError()

    def get_field(self, n):
        """
        Gets a field based on it's field number

        :param self: the BaseCard object
        :param n: the integer field number
        :returns value: the value of the field

        .. code-block:: python

          nid = 1
          node = model.nodes[nid]
          # ['GRID', nid, cp, x, y, z]
          z = node.get_field(5)

        """
        try:
            key_name = self._field_map[n]
            value = getattr(self, key_name)
        except KeyError:
            value = self._get_field_helper(n)
        return value

    def write_code_aster(self):
        """default method for Code Aster"""
        return ('# skipping %s because write_code_aster is not implemented\n'
                % self.type)

    def _verify(self, xref):
        """
        Verifies all methods for this object work

        :param self: the object pointer
        :param xref: has this model been cross referenced
        :type xref:  bool
        """
        print('# skipping _verify (type=%s) because _verify is '
              'not implemented\n' % self.type)

    def _is_same_fields(self, fields1, fields2):
        for (field1, field2) in zip(fields1, fields2):
            if not is_same(field1, field2):
                return False
        return True

    def _is_same_card(self, card):
        if self.type != card.type:
            return False
        fields1 = self.raw_fields()
        fields2 = card.raw_fields()
        return self._is_same_fields(fields1, fields2)

    def print_raw_card(self, size=8, is_double=False):  ## TODO: needs to be fixed
        """A card's raw fields include all defaults for all fields"""
        list_fields = self.raw_fields()
        return self.comment() + print_card(list_fields, size=size, is_double=is_double)

    def repr_fields(self):
        """gets the card's fields"""
        return self.raw_fields()

    def print_card(self, size=8, is_double=False):
        list_fields = self.repr_fields()
        return self.comment() + print_card(list_fields, size=size, is_double=is_double)

    def print_repr_card(self, size=8, is_double=False):
        list_fields = self.repr_fields()
        return self.comment() + print_card(list_fields, size=size, is_double=is_double)

    def __repr__(self):
        """
        Prints a card in the simplest way possible
        (default values are left blank).
        """
        try:
            return self.print_card(size=8)
        except:
            try:
                return self.print_card(size=16)
            except:
                print('problem printing %s card' % self.type)
                fields = self.repr_fields()
                print("fields = ", fields)
                raise


class Property(BaseCard):
    def __init__(self, card, data):
        assert card is None or data is None

    def Pid(self):
        """
        returns the property ID of an property

        :param self:  the Property pointer
        :returns pid: the Property ID
        :type pid:    int
        """
        return self.pid

    def Mid(self):
        """
        returns the material ID of an element

        :param self:  the Property pointer
        :returns mid: the Material ID
        :type mid:    int
        """
        if isinstance(self.mid, int):
            return self.mid
        else:
            return self.mid.mid

    def cross_reference(self, model):
        """dummy cross reference method for a Property"""
        msg = ' which is required by %s pid=%s' % (self.type, self.pid)
        self.mid = model.Material(self.mid, msg)


class Material(BaseCard):
    """Base Material Class"""
    def __init__(self, card, data):
        BaseCard.__init__(self)

    def cross_reference(self, model):
        """dummy cross reference method for a Material"""
        pass

    def Mid(self):
        """
        returns the material ID of an element

        :param self:  the Property pointer
        :returns mid: the Material ID
        :type mid:    int
        """
        return self.mid


class Element(BaseCard, ElementDeprecated):
    """defines the Element class"""
    pid = 0  # CONM2, rigid

    def __init__(self, card, data):
        BaseCard.__init__(self)
        assert card is None or data is None
        #: the list of node IDs for an element (default=None)
        self.nodes = None

    def Pid(self):
        """
        returns the property ID of an element
        :param self:  the Element pointer
        :returns pid: the Property ID
        :type pid:    int
        """
        if isinstance(self.pid, int):
            return self.pid
        else:
            return self.pid.pid

    def get_node_positions(self, nodes=None):
        """returns the positions of multiple node objects"""
        if not nodes:
            nodes = self.nodes

        positions = []
        for node in nodes:
            if node is not None:
                positions.append(node.Position())
            else:
                positions.append(None)
        return positions

    def _nodeIDs(self, nodes=None, allowEmptyNodes=False, msg=''):
        """returns nodeIDs for repr functions"""
        try:
            if not nodes:
                nodes = self.nodes

            if allowEmptyNodes:
                nodes2 = []
                for i, node in enumerate(nodes):
                    #print("node=%r type=%r" % (node, type(node)))
                    if node == 0 or node is None:
                        nodes2.append(None)
                    elif isinstance(node, int):
                        nodes2.append(node)
                    else:
                        nodes2.append(node.nid)
                return nodes2
            else:
                try:
                    nodeIDs = []
                    for i, node in enumerate(nodes):
                        #print("node=%r type=%r" % (node, type(node)))
                        if isinstance(node, int):
                            nodeIDs.append(node)
                        else:
                            nodeIDs.append(node.nid)

                    #if isinstance(nodes[0], int):
                        #nodeIDs = [node for node in nodes]
                    #else:
                        #nodeIDs = [node.nid for node in nodes]
                except:
                    print('type=%s nodes=%s allowEmptyNodes=%s\nmsg=%s' % (
                        self.type, nodes, allowEmptyNodes, msg))
                    raise
                assert 0 not in nodeIDs, 'nodeIDs = %s' % nodeIDs
                return nodeIDs
        except:
            print('type=%s nodes=%s allowEmptyNodes=%s\nmsg=%s' % (
                self.type, nodes, allowEmptyNodes, msg))
            raise

    def prepare_node_ids(self, nids, allow_empty_nodes=False):
        """Verifies all node IDs exist and that they're integers"""
        self.nodes = []
        for nid in nids:
            if isinstance(nid, integer_types):
                self.nodes.append(nid)
            elif nid is None and allow_empty_nodes:
                self.nodes.append(None)
            else:  # string???
                #self.nodes.append(int(nid))
                msg = 'this element may have missing nodes...\n'
                msg += 'nids=%s allow_empty_nodes=False;\ntype(nid)=%s' % (nids, type(nid))
                raise RuntimeError(msg)

    @property  # I think this means you can just call it as an attribute...
    def faces(self):
        """
        Gets the faces of the element

        :returns:  dictionary with face number as the keys and
                   a list of nodes (integer pointers) as the values.
        .. note::  The order of the nodes are consistent with ANSYS numbering.

        Example
        =======
        >>> print element.faces
        """
        faces = {}
        nodes = self.node_ids
        if self.type.startswith('CHEX'):  # CHEXA
            faces[1] = [nodes[0], nodes[1], nodes[2], nodes[3]]
            faces[2] = [nodes[0], nodes[1], nodes[5], nodes[4]]
            faces[3] = [nodes[1], nodes[2], nodes[6], nodes[5]]
            faces[4] = [nodes[2], nodes[3], nodes[7], nodes[6]]
            faces[5] = [nodes[3], nodes[0], nodes[4], nodes[7]]
            faces[6] = [nodes[4], nodes[5], nodes[6], nodes[7]]
        elif self.type.startswith('CTET'):  # CTETRA
            faces[1] = [nodes[0], nodes[1], nodes[2]]
            faces[2] = [nodes[0], nodes[1], nodes[3]]
            faces[3] = [nodes[1], nodes[2], nodes[3]]
            faces[4] = [nodes[2], nodes[0], nodes[3]]
        elif self.type.startswith('CPYR'):  # PYRAMID
            faces[1] = [nodes[0], nodes[1], nodes[2], nodes[3]]
            faces[2] = [nodes[0], nodes[1], nodes[4]]
            faces[3] = [nodes[1], nodes[2], nodes[4]]
            faces[4] = [nodes[2], nodes[3], nodes[4]]
            faces[5] = [nodes[3], nodes[0], nodes[4]]
        elif self.type.startswith('CPEN'):  # CPENTA
            faces[1] = [nodes[0], nodes[1], nodes[2]]
            faces[2] = [nodes[3], nodes[4], nodes[5]]
            faces[3] = [nodes[0], nodes[1], nodes[4], nodes[3]]
            faces[4] = [nodes[1], nodes[2], nodes[5], nodes[4]]
            faces[5] = [nodes[2], nodes[0], nodes[3], nodes[5]]
        elif self.type.startswith('CQUAD'):  # CQUADx
            # both sides
            faces[1] = [nodes[0], nodes[1], nodes[2], nodes[3]]
            faces[2] = [nodes[1], nodes[0], nodes[3], nodes[2]]
        elif self.type.startswith('CTRI'):  # CTRIAx
            # both sides
            faces[1] = [nodes[0], nodes[1], nodes[2]]
            faces[2] = [nodes[1], nodes[0], nodes[2]]
        else:
            faces = None
        return faces

    def nodes2face(self, nodes):
        """
        returns the face number that matches the list of nodes input

        :param nodes:   list of nodes
        :returns faces: the face number as an integer

        .. warning:: It's assumed you have the nodes in the proper order.
        """
        assert isinstance(nodes, list), 'nodes=%s' % str(nodes)
        face = None
        for i in self.faces.keys():
            chck = self.faces[i]
            if nodes == chck:
                face = i
        return face


def expand_thru(fields, set_fields=True, sort_fields=False):
    """
    Expands a list of values of the form [1,5,THRU,9,13]
    to be [1,5,6,7,8,9,13]

    :param fields:      the fields to expand
    :param set_fields:  should the fields be converted to a set and then back
                        to a list? This is useful for [2, 'THRU' 5, 1] (default=True)
    :param sort_fields: should the fields be sorted at the end? (default=False)
    """
    # ..todo:  should this be removed...is the field capitalized when read in?
    fields = [field.upper() if isinstance(field, string_types) else field for field in fields]

    if isinstance(fields, int):
        return [fields]
    if len(fields) == 1:
        return [int(fields[0])]
    out = []
    nfields = len(fields)
    i = 0
    while i < nfields:
        if fields[i] == 'THRU':
            istart = int(fields[i - 1])
            iend = int(fields[i + 1])
            for j in range(istart, iend + 1): # adding 1 to iend for the range offset
                out.append(j)
            i += 2
        else:
            out.append(int(fields[i]))
            i += 1

    if set_fields:
        out = list(set(out))
    if sort_fields:
        out.sort()
    return out


def expand_thru_by(fields, set_fields=True, sort_fields=False):
    """
    Expands a list of values of the form [1,5,THRU,9,BY,2,13]
    to be [1,5,7,9,13]

    :param fields:      the fields to expand
    :param set_fields:  should the fields be converted to a set and then back
                        to a list? This is useful for [2, 'THRU' 5, 1] (default=True)
    :param sort_fields: should the fields be sorted at the end? (default=False)

    .. todo:: not tested
    .. note:: used for QBDY3 and what else ???
    """
    # ..todo:  should this be removed...is the field capitalized when read in?
    fields = [field.upper()
              if isinstance(field, string_types) else field for field in fields]

    if len(fields) == 1:
        return [interpret_value(fields[0])]
    out = []
    nfields = len(fields)
    i = 0
    by = 1
    while i < nfields:
        if fields[i] == 'THRU':
            by = 1
            byCase = False
            if i + 2 < nfields and fields[i + 2] == 'BY':
                by = interpret_value(fields[i + 3])
            else:
                by = 1
                byCase = True
            min_value = interpret_value(fields[i - 1])
            max_value = interpret_value(fields[i + 1])
            maxR = int((max_value - min_value) // by + 1)  # max range value

            for j in range(0, maxR):  # +1 is to include final point
                value = min_value + by * j
                out.append(value)

            if byCase:  # null/standard case
                i += 2
            else:     # BY case
                i += 3
        else:
            out.append(interpret_value(fields[i]))
            i += 1

    if set_fields:
        out = list(set(out))
    if sort_fields:
        out.sort()
    return out


def expand_thru_exclude(self, fields):
    """
    Expands a list of values of the form [1,5,THRU,11,EXCEPT,7,8,13]
    to be [1,5,6,9,10,11,13]

    .. warning:: hasnt been tested
    """
    # ..todo:  should this be removed...is the field capitalized when read in?
    fields = [interpret_value(field.upper())
              if isinstance(field, string_types) else field for field in fields]

    fields_out = []
    nfields = len(fields)
    for i in range(nfields):
        if fields[i] == 'THRU':
            sorted_list = []
            for j in range(fields[i - 1], fields[i + 1]):
                sorted_list.append(fields[j])

        elif fields[i] == 'EXCLUDE':
            stored_set = set(sorted_list)
            while fields[i] < max(sorted_list):
                stored_set.remove(fields[i])
            sorted_list = list(stored_set)
        else:
            if sorted_list:
                fields_out += sorted_list
            fields_out.append(fields[i])
    return fields_out


def collapse_thru_by(fields, get_packs=False):
    """
    :param fields:    the list of fields to collapse
    :param get_packs: get the list of packs so "special" formatting can be done

    fields              packs
    [1, 2, 3...150]  -> [1, 150, 1]
    [1, 3, 5...150]  -> [1, 150, 2]
    """
    assert 'THRU' not in fields, fields
    fields.sort()
    packs = condense(fields)
    if get_packs:
        return packs
    fields2 = build_thru(packs)
    #assert fields == expand_thru_by(fields2)  # why doesn't this work?
    return fields2


def collapse_thru_by_float(fields):
    assert 'THRU' not in fields, fields
    fields.sort()
    packs = condense(fields)
    fields2 = build_thru_float(packs)
    #assert fields == expand_thru_by(fields2)  # why doesn't this work?
    return fields2


def collapse_thru(fields):
    assert 'THRU' not in fields, fields
    fields.sort()
    packs = condense(fields)
    fields2 = build_thru(packs, max_dv=1)
    #assert fields == expand_thru_by(fields2), fields2  # why doesn't this work?
    return fields2


def collapse_thru_packs(fields):
    assert 'THRU' not in fields, fields
    fields.sort()
    packs = condense(fields)
    singles, doubles = build_thru_packs(packs, max_dv=1)

    #assert fields == expand_thru_by(fields2), fields2  # why doesn't this work?
    return singles, doubles


def condense(value_list):
    """
    Builds a list of packs (list of 3 values representing the first, last,
    and delta values for condensing a SET card.

    .. seealso:: build_thru
    """
    if len(value_list) == 0:
        return []
    if len(value_list) == 1:
        return [[value_list[0], value_list[0], 1]]
    value_list.sort()
    packs = []

    dv_old = None
    firstVal = value_list[0]
    lastVal = firstVal

    for val in value_list[1:]:
        try:
            dv = val - lastVal
        except TypeError:
            print("lastVal=%r val=%r" % (lastVal, val))
            print("valueList=%r" % value_list)
            raise

        # sets up the first item of the pack
        if dv_old is None:
            dv_old = dv

        # fill up the pack
        if dv_old == dv:
            lastVal = val
        else:
            packs.append([firstVal, lastVal, dv_old])
            lastVal = val
            dv_old = None
            firstVal = val

    # fills the last pack
    if dv_old == dv:
        packs.append([firstVal, val, dv])
    else:
        packs.append([firstVal, val, dv_old])
    return packs


def build_thru_packs(packs, max_dv=1):
    """
    # invalid
    SET1,4000, 1, 3, THRU, 10, 20, THRU, 30

    # valid
    SET1,4000, 1
    SET1,4000, 3,  THRU, 10
    SET1,4000, 20, THRU, 30

    returns
      singles = [1]
      doubles = [[3, 'THRU', 10], [20, 'THRU', 30]]
    """
    singles = []
    doubles = []
    for (first_val, last_val, by) in packs:
        if first_val == last_val:
            singles.append(first_val)
        else:
            if by == 1:
                if last_val - first_val < 3: # dont make extra THRU cards
                    for val in range(first_val, last_val + 1):
                        singles.append(val)
                else:
                    double = [first_val, 'THRU', last_val]
                    doubles.append(double)
            else:
                for val in range(first_val, last_val + 1, by):
                    singles.append(val)
    return singles, doubles


def build_thru(packs, max_dv=None):
    """
    Takes a pack [1,7,2] and converts it into fields used by a SET card.
    The values correspond to the first value, last value, and delta in the
    list.  This means that [1,1001,2] represents 500 values.
    [1,1001,1] represents 1001 values and will be written as [1,THRU,1001]..

    :param packs: list of packs (list of 3 values: [first, last, delta] )
    :param maxDV: integer defining the max allowable delta between two values
            (default=None; no limit)
    """
    fields = []
    for (first_val, last_val, dv) in packs:
        if first_val == last_val:
            fields.append(first_val)
        elif dv == 1:
            if last_val - first_val > 2:
                fields.append(first_val)
                fields.append('THRU')
                fields.append(last_val)
            elif last_val - first_val == 2:
                fields.append(first_val)
                fields.append(first_val + 1)
                fields.append(last_val)
            else:
                fields.append(first_val)
                fields.append(last_val)
        else:
            if max_dv is None:
                if last_val - first_val > 4 * dv:
                    fields.append(first_val)
                    fields.append('THRU')
                    fields.append(last_val)
                    fields.append('BY')
                    fields.append(dv)
                else:
                    for v in range(first_val, last_val + dv, dv):
                        fields.append(v)
            else:
                for v in range(first_val, last_val + dv, dv):
                    fields.append(v)
    return fields


def build_thru_float(packs, max_dv=None):
    """
    Takes a pack [1,7,2] and converts it into fields used by a SET card.
    The values correspond to the first value, last value, and delta in the
    list.  This means that [1,1001,2] represents 500 values.
    [1,1001,1] represents 1001 values and will be written as [1,THRU,1001]..

    :param packs: list of packs (list of 3 values: [first, last, delta] )
    :param max_dv: integer defining the max allowable delta between two values
                   (default=None; no limit)
    """
    fields = []
    for (first_val, last_val, dv) in packs:
        if last_val - first_val > 4 * dv:
            fields.append(first_val)
            fields.append('THRU')
            fields.append(last_val)
            fields.append('BY')
            fields.append(dv)
        else:
            nv = int(round((last_val - first_val) / dv)) + 1
            for i in range(nv):
                v = first_val + i * dv
                fields.append(v)
    return fields
