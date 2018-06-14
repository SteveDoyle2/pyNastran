from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
from typing import List, Dict, Union, Optional, Any
from six import string_types, PY2
from six.moves import zip, range

import numpy as np
#from numpy import nan, empty, unique

from pyNastran.bdf.bdf_interface.bdf_card import BDFCard
from pyNastran.utils import object_attributes, object_methods, integer_types
from pyNastran.bdf.field_writer import print_card
from pyNastran.bdf.field_writer_8 import is_same
from pyNastran.bdf.bdf_interface.assign_type import interpret_value
from pyNastran.bdf.utils import deprecated


if PY2:
    def long_range(start, stop=None, step=1):
        if isinstance(stop, long):
            out = []
            i = start
            while i < stop:
                out.append(i)
                i += step
            return out
        elif stop is None:
            i = 0
            out = []
            while i < start:
                out.append(i)
                i += step
            return out
        else:
            return xrange(start, stop, step)
    range = long_range
else:
    pass


class BaseCard(object):
    """
    Defines a series of base methods for every card class
    (e.g., GRID, CTRIA3) including:

     - deepcopy()
     - get_stats()
     - validate()
     - object_attributes(mode='public', keys_to_skip=None)
     - object_methods(self, mode='public', keys_to_skip=None)
     - comment
     - update_field(self, n, value)

    """
    def __init__(self):
        pass

    def __deepcopy__(self, memo_dict):
        #raw_fields = self.repr_fields()
        raw_fields = self.raw_fields()
        card = BDFCard(raw_fields)
        return self.add_card(card)

    def get_stats(self):
        # type: () -> str
        """Prints out an easy to read summary of the card"""
        msg = '---%s---\n' % self.type
        for name in sorted(self.object_attributes()):
            value = getattr(self, name)
            msg += '  %-6s : %r\n' % (name, value)
        return msg

    def deprecated(self, old_name, new_name, deprecated_version):
        # type: (str, str, str) -> None
        """deprecates methods"""
        deprecated(old_name, new_name, deprecated_version, levels=[0, 1, 2])

    def validate(self):
        """card checking method that should be overwritten"""
        pass

    def object_attributes(self, mode='public', keys_to_skip=None):
        """.. seealso:: `pyNastran.utils.object_attributes(...)`"""
        if keys_to_skip is None:
            keys_to_skip = []

        my_keys_to_skip = [
        ]
        return object_attributes(self, mode=mode, keys_to_skip=keys_to_skip+my_keys_to_skip)

    def object_methods(self, mode='public', keys_to_skip=None):
        """.. seealso:: `pyNastran.utils.object_methods(...)`"""
        if keys_to_skip is None:
            keys_to_skip = []
        my_keys_to_skip = []

        my_keys_to_skip = [
        ]
        return object_methods(self, mode=mode, keys_to_skip=keys_to_skip+my_keys_to_skip)

    @property
    def comment(self):
        # type: () -> str
        """accesses the comment"""
        # just for testing
        #self.deprecated('comment()', 'comment2()', '0.7')
        if hasattr(self, '_comment'):
            return '%s' % self._comment
        return ''

    @comment.setter
    def comment(self, new_comment):
    # type: (str) -> None
        """sets a comment"""
        #comment = new_comment.rstrip()
        #self._comment = comment + '\n' if comment else ''
        self._comment = _format_comment(new_comment)

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
        # type: (int, Optional[Union[int, float, str]]) -> None
        """
        Updates a field based on it's field number.

        Parameters
        ----------
        n : int
            the field number
        value : int/float/str/None
           the value to update the field to

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
        """
        dynamic method for non-standard attributes
        (e.g., node.update_field(3, 0.1) to update z)

        """
        msg = '%s has not overwritten _update_field_helper; out of range' % self.__class__.__name__
        raise IndexError(msg)

    def _get_field_helper(self, n):
        """dynamic method for non-standard attributes (e.g., node.get_field(3, 0.1) to get z)"""
        msg = '%s has not overwritten _get_field_helper; out of range' % self.__class__.__name__
        raise IndexError(msg)

    def get_field(self, n):
        # type: (int) -> Optional[Union[int, float, str]]
        """
        Gets a field based on it's field number

        Parameters
        ----------
        n : int
           the field number

        Returns
        -------
        value : int/float/str/None
            the value of the field

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

    def _verify(self, xref):
        # type: (bool) -> None
        """
        Verifies all methods for this object work

        Parameters
        ----------
        xref : bool
            has this model been cross referenced

        """
        print('# skipping _verify (type=%s) because _verify is '
              'not implemented' % self.type)

    def __eq__(self, card):
        # type: (Any) -> bool
        """
        Enables functions like:

        .. code-block:: python

          >>> GRID(nid=1, ...) === GRID(nid=1, ...)
          True
          >>> GRID(nid=1, ...) === GRID(nid=2, ...)
          False
          >>> GRID(nid=1, ...) === CQUAD4(eid=1, ...)
          False

        """
        if not isinstance(card, self.__class__):
            return False
        if self.type != card.type:
            return False
        fields1 = self.raw_fields()
        fields2 = card.raw_fields()
        return self._is_same_fields(fields1, fields2)

    def _is_same_fields(self, fields1, fields2):
        # type: (List[Optional[Union[int, float, str]]], List[Optional[Union[int, float, str]]]) -> bool
        for (field1, field2) in zip(fields1, fields2):
            if not is_same(field1, field2):
                return False
        return True

    def print_raw_card(self, size=8, is_double=False):
        # type: (int, bool) -> str
        """A card's raw fields include all defaults for all fields"""
        list_fields = self.raw_fields()
        return self.comment + print_card(list_fields, size=size, is_double=is_double)

    def repr_fields(self):
    # type: () -> List[Optional[Union[int, float, str]]]
        """
        Gets the fields in their simplified form

        Returns
        -------
        fields : List[varies]
            the fields that define the card

        """
        return self.raw_fields()

    def print_card(self, size=8, is_double=False):
        # type: (int, bool) -> str
        """prints the card in 8/16/16-double format"""
        list_fields = self.repr_fields()
        return self.comment + print_card(list_fields, size=size, is_double=is_double)

    def print_repr_card(self, size=8, is_double=False):
        # type: (int, bool) -> str
        """prints the card in 8/16/16-double format"""
        list_fields = self.repr_fields()
        return self.comment + print_card(list_fields, size=size, is_double=is_double)

    def __repr__(self):
        # type: () -> str
        """
        Prints a card in the simplest way possible
        (default values are left blank).

        """
        comment = self.comment
        list_fields = self.repr_fields()
        try:
            return comment + print_card(list_fields, size=8)
        except:
            try:
                return comment + print_card(list_fields, size=16)
            except:
                print('problem printing %s card' % self.type)
                print("list_fields = ", list_fields)
                raise

    def rstrip(self):
        # type: () -> str
        try:
            msg = '%s' % str(self)
        except UnicodeEncodeError:
            comment = self.comment
            self.comment = ''
            msg = '$ dropped comment due to unicode error\n%s' % str(self)
            self.comment = comment
        return msg.rstrip()

    def write_card(self, size=8, is_double=False):
        # type: (int, bool) -> str
        """
        Writes the card with the specified width and precision

        Parameters
        ----------
        size : int (default=8)
            size of the field; {8, 16}
        is_double : bool (default=False)
            is this card double precision

        Returns
        -------
        msg : str
            the string representation of the card

        """
        raise NotImplementedError('%s has not overwritten write_card' % self.__class__.__name__)

    def write_card_16(self, is_double=False):
        fields = self.repr_fields()
        return print_card(fields, size=16, is_double=False)

class Property(BaseCard):
    """Base Property Class"""
    def __init__(self):
        # type: () -> None
        """dummy init"""
        pass

    def Pid(self):
        # type: () -> int
        """
        returns the property ID of an property

        Returns
        -------
        pid : int
            the Property ID

        """
        return self.pid

    def Mid(self):
        # type: () -> int
        """
        returns the material ID of an element

        Returns
        -------
        mid : int
            the Material ID

        """
        if self.mid_ref is None:
            return self.mid
        return self.mid_ref.mid

    def write_card_8(self):
        # type: () -> str
        return self.write_card()

    def write_card_16(self, is_double=False):
        # type: (bool) -> str
        return self.write_card()


class Material(BaseCard):
    """Base Material Class"""
    def __init__(self):
        # type: () -> None
        """dummy init"""
        BaseCard.__init__(self)

    @property
    def TRef(self):
        # type: () -> float
        if not hasattr(self, 'tref'):
            raise AttributeError('%r object has no attribute tref' % self.type)
        return self.tref

    @TRef.setter
    def TRef(self, tref):
        # type: (float) -> None
        """sets the self.Tref attributes"""
        if not hasattr(self, 'tref'):
            raise AttributeError('%r object has no attribute tref' % self.type)
        self.tref = tref

    def cross_reference(self, model):
        # type: (Any) -> None
        """dummy cross reference method for a Material"""
        pass

    def Mid(self):
        # type: () -> Any
        """
        returns the material ID of an element

        Returns
        -------
        mid : int
            the Material ID
        """
        return self.mid


class Element(BaseCard):
    """defines the Element class"""
    pid = 0  # CONM2, rigid

    def __init__(self):
        """dummy init"""
        BaseCard.__init__(self)
        #: the list of node IDs for an element (default=None)
        #self.nodes = None

    def verify_unique_node_ids(self):
        node_ids = self.node_ids
        self._verify_unique_node_ids(node_ids)

    def _verify_unique_node_ids(self, required_node_ids, non_required_node_ids=None):
        if required_node_ids:
            if non_required_node_ids:
                raise NotImplementedError('only required nodes implemented')
            else:
                urnids = np.unique(required_node_ids)
                n_unique_node_ids = len(urnids)
                n_node_ids = len(required_node_ids)
                if n_unique_node_ids != n_node_ids:
                    msg = 'nunique_node_ids=%s nnode_ids=%s' % (n_unique_node_ids, n_node_ids)
                    raise RuntimeError(msg)
        else:
            raise NotImplementedError('only required nodes implemented')

    def Pid(self):
        # type: () -> int
        """
        Gets the Property ID of an element

        Returns
        -------
        pid : int
            the Property ID

        """
        if self.pid_ref is None:
            return self.pid
        return self.pid_ref.pid

    def get_node_positions(self, nodes=None):
        # type: (Any) -> np.ndarray
        """returns the positions of multiple node objects"""
        if nodes is None:
            nodes = self.nodes_ref

        nnodes = len(nodes)
        positions = np.empty((nnodes, 3), dtype='float64')
        positions.fill(np.nan)
        for i, node in enumerate(nodes):
            assert not isinstance(node, int), self.type
            if node is not None:
                positions[i, :] = node.get_position()
        return positions

    def get_node_positions_no_xref(self, model, nodes=None):
        # type: (Any, Any) -> np.ndarray
        """returns the positions of multiple node objects"""
        if not nodes:
            nodes = self.nodes

        nnodes = len(nodes)
        positions = np.empty((nnodes, 3), dtype='float64')
        positions.fill(np.nan)
        for i, nid in enumerate(nodes):
            if nid is not None:
                node = model.Node(nid)
                positions[i, :] = node.get_position_no_xref(model)
        return positions

    def _node_ids(self, nodes=None, allow_empty_nodes=False, msg=''):
        # type: (Optional[List[Any]], bool, str) -> List[int]
        """returns nodeIDs for repr functions"""
        return _node_ids(self, nodes, allow_empty_nodes, msg)

    def prepare_node_ids(self, nids, allow_empty_nodes=False):
        # type: (List[int], bool) -> None
        """Verifies all node IDs exist and that they're integers"""
        self.nodes = nids
        self.validate_node_ids(allow_empty_nodes)

    def validate_node_ids(self, allow_empty_nodes=False):
        # type: (bool) -> None
        if allow_empty_nodes:
            # only put valid nodes in here
            nids2 = [nid for nid in self.nodes ]
            if len(nids2) == 0:
                msg = '%s requires at least one node id be specified; node_ids=%s' % (
                    self.type, nids2)
                raise ValueError(msg)

            #unique_nodes = unique(nids2)
            #if len(nids2) != len(unique_nodes):
                #msg = '%s requires that all node ids be unique; node_ids=%s' % (self.type, nids2)
                #raise IndexError(msg)

            # remove 0 nodes
            nodes = [nid if nid is not 0 else None
                     for nid in self.nodes]
        else:
            nodes = self.nodes
            #unique_nodes = unique(self.nodes)
            #if len(self.nodes) != len(unique_nodes):
                #msg = '%s requires that all node ids be unique; node_ids=%s' % (
                    #self.type, self.nodes)
                #raise IndexError(msg)


        nodes2 = []
        for nid in nodes:
            if isinstance(nid, integer_types):
                nodes2.append(nid)
            elif nid is None and allow_empty_nodes:
                nodes2.append(None)
            else:  # string???
                #nodes.append(int(nid))
                msg = 'this element may have missing nodes...\n'
                msg += 'nids=%s allow_empty_nodes=False;\ntype(nid)=%s' % (self.nodes, type(nid))
                raise RuntimeError(msg)
        self.nodes = nodes2

    # this doesn't belong in Element as CBARs don't have faces
    #@property
    #def faces(self):
        ## () -> Dict[int, List[int]]
        #"""
        #Gets the faces of the element

        #Returns
        #-------
        #faces : Dict[int] = [face1, face2, ...]
            #key = face number
            #value = a list of nodes (integer pointers) as the values.

        #Notes
        #-----
        #The order of the nodes are consistent with normals that point outwards
        #The face numbering is meaningless

        #.. note::  The order of the nodes are consistent with ANSYS numbering; is this current?
        #.. warning:: higher order element ids not verified with ANSYS; is this current?

        #Examples
        #--------
        #>>> print(element.faces)

        #"""
        #faces = {}
        #try:
            #nodes = self.node_ids
        #except AttributeError:
            #return None
        #if self.type.startswith('CQUAD'): # CQUADx
            ## both sides
            #faces[1] = [nodes[0], nodes[1], nodes[2], nodes[3]]  # CQUAD8/9?
            #faces[2] = [nodes[1], nodes[0], nodes[3], nodes[2]]
        #elif self.type.startswith('CTRI'):  # CTRIAx
            ## both sides
            #faces[1] = [nodes[0], nodes[1], nodes[2]]  # CTRIA6?
            #faces[2] = [nodes[1], nodes[0], nodes[2]]
        #else:
            #faces = None
        #return faces

    #def nodes2face(self, nodes):
        #"""
        #returns the face number that matches the list of nodes input

        #Parameters
        #----------
        #nodes : List[node]
            #a series of nodes

        #Returns
        #-------
        #face_id : int
            #the face number as an integer

        #.. warning:: It's assumed you have the nodes in the proper order.

        #"""
        #assert isinstance(nodes, list), 'nodes=%s' % str(nodes)
        #face = None
        #for i in self.faces.keys():
            #chck = self.faces[i]
            #if nodes == chck:
                #face = i
        #return face

def _format_comment(comment):
    # type: (str) -> str
    r"""Format a card comment to precede the card using
    nastran-compatible comment character $. The comment
    string can have multiple lines specified as linebreaks.

    Empty comments or just spaces are returned as an empty string.

    Examples
    --------
    >>> _format_comment('a comment\ntaking two lines')
    $a comment
    $taking two lines

    >>> _format_comment('')
    <empty string>

    >>> _format_comment('       ')
    <empty string>

    >>> _format_comment('$ a comment within a comment looks weird')
    $$ a comment within a comment looks weird

    """
    if comment.strip() == '':  # deals with a bunch of spaces
        return ''
    else:
        return ''.join([u'${}\n'.format(comment_line)
                        for comment_line in comment.rstrip().split('\n')])

def _node_ids(card, nodes=None, allow_empty_nodes=False, msg=''):
    try:
        if not nodes:
            nodes = card.nodes
            assert nodes is not None, card.__dict__

        if allow_empty_nodes:
            nodes2 = []
            for node in nodes:
                if node == 0 or node is None:
                    nodes2.append(None)
                elif isinstance(node, integer_types):
                    nodes2.append(node)
                else:
                    nodes2.append(node.nid)
            assert nodes2 is not None, str(card)
            return nodes2
        else:
            try:
                node_ids = []
                for node in nodes:
                    if isinstance(node, integer_types):
                        node_ids.append(node)
                    else:
                        node_ids.append(node.nid)

                #if isinstance(nodes[0], integer_types):
                    #node_ids = [node for node in nodes]
                #else:
                    #node_ids = [node.nid for node in nodes]
            except:
                print('type=%s nodes=%s allow_empty_nodes=%s\nmsg=%s' % (
                    card.type, nodes, allow_empty_nodes, msg))
                raise
            assert 0 not in node_ids, 'node_ids = %s' % node_ids
            assert node_ids is not None, str(card)
            return node_ids
    except:
        print('type=%s nodes=%s allow_empty_nodes=%s\nmsg=%s' % (
            card.type, nodes, allow_empty_nodes, msg))
        raise
    raise RuntimeError('huh...')

def expand_thru(fields, set_fields=True, sort_fields=False):
    # type: (List[str], bool, bool) -> List[int]
    """
    Expands a list of values of the form [1,5,THRU,9,13]
    to be [1,5,6,7,8,9,13]

    Parameters
    ----------
    fields : List[int/str]
        the fields to expand
    set_fields : bool; default=True
        Should the fields be converted to a set and then back to a list?
        This is useful for [2, 'THRU' 5, 1]
    sort_fields : bool; default=False
        Should the fields be sorted at the end?

    """
    # ..todo:  should this be removed...is the field capitalized when read in?
    if isinstance(fields, integer_types):
        return [fields]
    #elif isinstance(fields[0], integer_types):  # don't use this [1, 'THRU', 10]
        #return fields
    elif len(fields) == 1:
        return [int(fields[0])]

    fields = [field.upper()
              if isinstance(field, string_types) else field for field in fields]

    out = []
    nfields = len(fields)
    i = 0
    while i < nfields:
        if isinstance(fields[i], string_types) and fields[i] == 'THRU':
            istart = int(fields[i - 1])
            iend = int(fields[i + 1])

            # adding 1 to iend for the range offset
            for j in range(istart, iend + 1):
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


def expand_thru_by(fields, set_fields=True, sort_fields=True,
                   require_int=True, allow_blanks=False):
    # type: (List[str], bool, bool, bool, bool) -> List[int]
    """
    Expands a list of values of the form [1,5,THRU,9,BY,2,13]
    to be [1,5,7,9,13]

    Parameters
    ----------
    fields : List[int/str]
        the fields to expand
    set_fields : bool; default=True
        Should the fields be converted to a set and then back to a list
        to remove duplicates?
        This is useful for [2, 'THRU' 5, 1]
    sort_fields : bool; default=False
        Should the fields be sorted at the end?
    require_int : bool; default=True
        True : all data must be integers
        False : floats are allowed (e.g., DDVAL)
    allow_blanks : bool; default=Fals
        True : blank/Nones are ignored (e.g., NSM1/NSML1)
        False : crash

    .. todo:: not tested

    Notes
    -----
    used for QBDY3 and what else ???

    """
    if require_int:
        func = int
    else:
        func = interpret_value

    # ..todo:  should this be removed...is the field capitalized when read in?
    fields = [field.upper()
              if isinstance(field, string_types) else field for field in fields]

    if len(fields) == 1:
        return [func(fields[0])]
    out = []
    nfields = len(fields)
    i = 0
    by = 1
    while i < nfields:
        #print('fields[i]=%r' % fields[i])
        is_blank = (
            allow_blanks and (
                (isinstance(fields[i], string_types) and fields[i].strip() == '') or
                fields[i] is None)
        )
        if is_blank:
            #print('blank=%s' % fields[i])
            i += 1
            continue
        if fields[i] == 'THRU':
            by = 1
            by_case = False
            if i + 2 < nfields and fields[i + 2] == 'BY':
                by = func(fields[i + 3])
            else:
                by = 1
                by_case = True
            min_value = func(fields[i - 1])
            max_value = func(fields[i + 1])
            max_range = int((max_value - min_value) // by + 1)  # max range value

            for j in range(0, max_range):  # +1 is to include final point
                value = min_value + by * j
                out.append(value)
            out.append(max_value)

            if by_case:  # null/standard case
                # A thru B
                i += 2
            else:     # BY case
                # A thru B by C
                i += 4
        else:
            out.append(func(fields[i]))
            i += 1

    if set_fields:
        out = list(set(out))
    if sort_fields:
        out.sort()
    return out


def expand_thru_exclude(fields):
    # type: (List[str]) -> List[int]
    """
    Expands a list of values of the form [1,5,THRU,11,EXCEPT,7,8,13]
    to be [1,5,6,9,10,11,13]

    .. warning:: hasn't been tested

    """
    # ..todo:  should this be removed...is the field capitalized when read in?
    fields = [interpret_value(field.upper())
              if isinstance(field, string_types) else field for field in fields]

    fields_out = []  # type: List[int]
    nfields = len(fields)
    for i in range(nfields):
        #print('fields[%i] = %r' % (i, fields[i]))
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

def break_word_by_trailing_integer(pname_fid):
    """
    Splits a word

    Examples
    --------
    >>> break_word_by_trailing_integer('T11')
    ('T', '11')
    >>> break_word_by_trailing_integer('THETA11')
    ('THETA', '11')

    """
    nums = []
    i = 0
    for i, letter in enumerate(reversed(pname_fid)):
        if letter.isdigit():
            nums.append(letter)
        else:
            break

    num = ''.join(nums[::-1])
    if not num:
        msg = ("pname_fid=%r does not follow the form 'T1', 'T11', 'THETA42' "
               "(letters and a number)" % pname_fid)
        raise SyntaxError(msg)
    word = pname_fid[:-i]
    assert len(word)+len(num) == len(pname_fid), 'word=%r num=%r pname_fid=%r' % (word, num, pname_fid)
    return word, num

def break_word_by_trailing_parentheses_integer_ab(pname_fid):
    """
    Splits a word

    Examples
    --------
    >>> break_word_by_trailing_parentheses_integer('A(11)')
    ('A', '11')
    >>> break_word_by_trailing_parentheses_integer('NSM(11)')
    ('NSM', '11')
    >>> break_word_by_trailing_parentheses_integer('NSM(B)')
    ('NSM', 'B')

    """
    assert pname_fid.endswith(')'), pname_fid
    word, num = pname_fid[:-1].split('(')
    if num not in ['A', 'B']:
        num = int(num)
    return word, num
