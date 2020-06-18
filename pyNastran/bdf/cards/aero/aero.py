# coding: utf-8
# pylint: disable=R0902,R0904,R0914,C0302,C0111,C0103,R0913
"""
All aero cards are defined in this file.  This includes:

 * AECOMP
 * AEFACT
 * AELINK
 * AELIST
 * AEPARM
 * AESURF / AESURFS
 * CAERO1 / CAERO2 / CAERO3 / CAERO4 / CAERO5
 * PAERO1 / PAERO2 / PAERO3 / PAERO4 / PAERO5
 * SPLINE1 / SPLINE2 / SPLINE3 / SPLINE4 / SPLINE5
 * MONPNT1 / MONPNT2 / MONPNT3

All cards are BaseCard objects.

"""
from __future__ import annotations
from itertools import count
import math
from typing import List, Union, Any, TYPE_CHECKING

import numpy as np

from pyNastran.utils.numpy_utils import integer_types
from pyNastran.utils import object_attributes
from pyNastran.bdf.field_writer_8 import set_blank_if_default, print_card_8, print_float_8
from pyNastran.bdf.cards.base_card import BaseCard, expand_thru
from pyNastran.bdf.bdf_interface.assign_type import (
    fields, integer, integer_or_blank, double, double_or_blank, string,
    string_or_blank, integer_or_string,
    interpret_value, parse_components, components_or_blank)
from pyNastran.bdf.cards.utils import wipe_empty_fields
from pyNastran.bdf.cards.aero.utils import (
    points_elements_from_quad_points, create_axisymmetric_body)
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf import BDF, BDFCard


class AECOMP(BaseCard):
    """
    Defines a component for use in monitor point definition or external splines.

    +--------+-------+----------+-------+-------+-------+-------+-------+-------+
    |   1    |   2   |    3     |   4   |   5   |   6   |   7   |   8   |   9   |
    +========+=======+==========+=======+=======+=======+=======+=======+=======+
    | AECOMP | NAME  | LISTTYPE | LIST1 | LIST2 | LIST3 | LIST4 | LIST5 | LIST6 |
    +--------+-------+----------+-------+-------+-------+-------+-------+-------+
    |        | LIST7 |   etc.   |       |       |       |       |       |       |
    +--------+-------+----------+-------+-------+-------+-------+-------+-------+
    | AECOMP | WING  |  AELIST  | 1001  | 1002  |       |       |       |       |
    +--------+-------+----------+-------+-------+-------+-------+-------+-------+

    Attributes
    ----------
    name : str
        The name.
    list_type : str
        {'SET1', 'AELIST', 'CAEROx'}
    lists : list[int]
        list of values of AECOMP lists

    """
    type = 'AECOMP'
    allowed_list_types = ['SET1', 'AELIST', 'CAERO']

    @classmethod
    def _init_from_empty(cls):
        name = 'name'
        list_type = 'CAERO'
        lists = [1]
        return AECOMP(name, list_type, lists, comment='')

    def __init__(self, name: str, list_type: List[str],
                 lists: Union[int, List[int]],
                 comment: str='') -> None:
        """
        Creates an AECOMP card

        Parameters
        ----------
        name : str
            the name of the component
        list_type : str
            One of CAERO, AELIST or CMPID for aerodynamic components and
            SET1 for structural components. Aerodynamic components are
            defined on the aerodynamic ks-set mesh while the structural
            components are defined on the g-set mesh.
        lists : List[int, int, ...]; int
            The identification number of either SET1, AELIST or CAEROi
            entries that define the set of grid points that comprise
            the component
        comment : str; default=''
            a comment for the card

        """
        BaseCard.__init__(self)
        if comment:
            self.comment = comment
        if isinstance(lists, integer_types):
            lists = [lists]
        elif not isinstance(lists, (list, tuple)):
            raise TypeError('AECOMP; type(lists)=%s and must be a list/tuple' % type(lists))

        self.name = name
        self.list_type = list_type
        self.lists = lists
        self.lists_ref = None

    def validate(self) -> None:
        if not self.list_type in ['SET1', 'AELIST', 'CAERO', 'CMPID']:
            msg = 'list_type=%r not in [SET1, AELIST, CAERO, CMPID]' % self.list_type
            raise RuntimeError(msg)

    @classmethod
    def add_card(cls, card: Any, comment: str='') -> AECOMP:
        """
        Adds an AECOMP card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        name = string(card, 1, 'name')
        list_type = string(card, 2, 'list_type')
        j = 1
        lists = []
        for i in range(3, len(card)):
            list_i = integer(card, i, '%s_%i' % (list_type, j))
            lists.append(list_i)
            j += 1
        return AECOMP(name, list_type, lists, comment=comment)

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        msg = ', which is required by AECOMP name=%r' % self.name
        if self.list_type == 'SET1':
            self.lists_ref = [model.SET1(key, msg) for key in self.lists]
        elif self.list_type == 'AELIST':
            self.lists_ref = [model.AELIST(key, msg) for key in self.lists]
        elif self.list_type == 'CAERO':
            self.lists_ref = [model.CAero(key, msg) for key in self.lists]
        #elif self.list_type == 'CMPID':
            # AEQUAD4,/AETRIA3
        else:
            raise NotImplementedError(self.list_type)

    def safe_cross_reference(self, model: BDF):
        msg = ', which is required by AECOMP name=%r' % self.name
        #return
        lists_ref = []
        if self.list_type == 'SET1':
            for key in self.lists:
                try:
                    ref = model.SET1(key, msg)
                except KeyError:
                    ref = None
                lists_ref.append(ref)
        elif self.list_type == 'AELIST':
            for key in self.lists:
                try:
                    ref = model.AELIST(key, msg)
                except KeyError:
                    ref = None
                lists_ref.append(ref)
        elif self.list_type == 'CAERO':
            for key in self.lists:
                try:
                    ref = model.CAero(key, msg)
                except KeyError:
                    ref = None
                lists_ref.append(ref)
        #elif self.list_type == 'CMPID':
            # AEQUAD4,/AETRIA3
        else:
            raise NotImplementedError(self.list_type)
        self.lists_ref = lists_ref

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.lists = self.get_lists()
        self.lists_ref = None

    def get_lists(self):
        if self.lists_ref is None:
            return self.lists
        if self.list_type == 'SET1':
            lists = [set1 if isinstance(set1, integer_types)
                     else set1.sid for set1 in self.lists_ref]
        elif self.list_type == 'AELIST':
            lists = [aelist if isinstance(aelist, integer_types)
                     else aelist.sid for aelist in self.lists_ref]
        elif self.list_type == 'CAERO':
            lists = [caero if isinstance(caero, integer_types)
                     else caero.eid for caero in self.lists_ref]
        #elif self.list_type == 'CMPID':
            # AEQUAD4,/AETRIA3
        else:
            raise NotImplementedError(self.list_type)
        return lists

    def raw_fields(self):
        list_fields = ['AECOMP', self.name, self.list_type] + self.get_lists()
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        """
        The writer method used by BDF.write_card()

        Parameters
        -----------
        size : int; default=8
            the size of the card (8/16)

        """
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class AECOMPL(BaseCard):
    """
    +---------+--------+--------+--------+---------+--------+--------+--------+--------+
    |    1    |    2   |    3   |    4   |    5    |    6   |    7   |    8   |    9   |
    +=========+========+========+========+=========+========+========+========+========+
    | AECOMPL |  NAME  | LABEL1 | LABEL2 | LABEL3  | LABEL4 | LABEL5 | LABEL6 | LABEL7 |
    +---------+--------+--------+--------+---------+--------+--------+--------+--------+
    |         | LABEL8 |  etc.  |        |         |        |        |        |        |
    +---------+--------+--------+--------+---------+--------+--------+--------+--------+
    | AECOMPL |  HORIZ |  STAB  |  ELEV  | BALANCE |        |        |        |        |
    +---------+--------+--------+--------+---------+--------+--------+--------+--------+
    """
    type = 'AECOMPL'

    @classmethod
    def _init_from_empty(cls):
        name = 'HORIZ'
        labels = 'ELEV'
        return AECOMPL(name, labels, comment='')

    def __init__(self, name: str,
                 labels: List[str],
                 comment: str='') -> None:
        """
        Creates an AECOMPL card

        Parameters
        ----------
        name : str
            the name of the component
        labels : List[str, str, ...]; str
            A string of 8 characters referring to the names of other components
            defined by either AECOMP or other AECOMPL entries.
        comment : str; default=''
            a comment for the card

        """
        BaseCard.__init__(self)
        if comment:
            self.comment = comment
        if isinstance(labels, str):
            labels = [labels]
        elif not isinstance(labels, (list, tuple)):
            raise TypeError('AECOMPL; type(labels)=%s and must be a list/tuple' % type(labels))

        self.name = name
        self.labels = labels
        #self.labels_ref = None

    @classmethod
    def add_card(cls, card: BDFCard, comment: str='') -> AECOMP:
        """
        Adds an AECOMPL card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        name = string(card, 1, 'name')
        labels = []
        j = 1
        for i in range(2, len(card)):
            label = string(card, i, 'label_%i' % j)
            labels.append(label)
            j += 1
        return AECOMPL(name, labels, comment=comment)

    def cross_reference(self, model: BDF) -> None:
        pass

    def safe_cross_reference(self, model):
        pass

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        pass

    def raw_fields(self):
        list_fields = ['AECOMPL', self.name] + self.labels
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        """
        The writer method used by BDF.write_card()

        Parameters
        -----------
        size : int; default=8
            the size of the card (8/16)

        """
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class AEFACT(BaseCard):
    """
    Defines real numbers for aeroelastic analysis.

    +--------+-----+----+--------+-----+----+----+----+----+
    |   1    |   2 |  3 |    4   |  5  |  6 |  7 | 8  |  9 |
    +========+=====+====+========+=====+====+====+====+====+
    | AEFACT | SID | D1 |   D2   | D3  | D4 | D5 | D6 | D7 |
    +--------+-----+----+--------+-----+----+----+----+----+
    |        | D8  | D9 |  etc.  |     |    |    |    |    |
    +--------+-----+----+--------+-----+----+----+----+----+
    | AEFACT | 97  |.3  |  0.7   | 1.0 |    |    |    |    |
    +--------+-----+----+--------+-----+----+----+----+----+

    TODO: Are these defined in percentages and thus,
          should they be normalized if they are not?

    """
    type = 'AEFACT'

    @classmethod
    def _init_from_empty(cls):
        sid = 1
        fractions = [0., 1.,]
        return AEFACT(sid, fractions, comment='')

    def __init__(self, sid, fractions, comment=''):
        """
        Creates an AEFACT card, which is used by the CAEROx / PAEROx card
        to adjust the spacing of the sub-paneleing (and grid point
        paneling in the case of the CAERO3).

        Parameters
        ----------
        sid : int
            unique id
        fractions : List[float, ..., float]
            list of percentages
        comment : str; default=''
            a comment for the card

        """
        super(AEFACT, self).__init__()
        if comment:
            self.comment = comment
        #: Set identification number. (Unique Integer > 0)
        self.sid = sid
        #: Number (float)
        self.fractions = np.asarray(fractions, dtype='float64')

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds an AEFACT card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        sid = integer(card, 1, 'sid')

        fractions = []
        for i in range(2, len(card)):
            fraction = double(card, i, 'factor_%i' % (i - 1))
            fractions.append(fraction)
        assert len(card) > 2, 'len(AEFACT card) = %i\n%s' % (len(card), card)
        return AEFACT(sid, fractions, comment=comment)

    #def cross_reference(self, model: BDF) -> None:
        #pass

    #def uncross_reference(self):
        #"""Removes cross-reference links"""
        #pass

    def raw_fields(self):
        """
        Gets the fields in their unmodified form

        Returns
        -------
        fields : List[int/float/str]
            the fields that define the card

        """
        list_fields = ['AEFACT', self.sid] + list(self.fractions)
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        """
        The writer method used by BDF.write_card()

        Parameters
        -----------
        size : int; default=8
            the size of the card (8/16)

        """
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class AELINK(BaseCard):
    r"""
    Defines relationships between or among AESTAT and AESURF entries, such
    that:

    .. math:: u^D + \Sigma_{i=1}^n C_i u_i^I = 0.0

    +--------+-------+-------+--------+------+-------+----+-------+----+
    |   1    |   2   |   3   |   4    |   5  |   6   |  7 |   8   |  9 |
    +========+=======+=======+========+======+=======+====+=======+====+
    | AELINK |  ID   | LABLD | LABL1  |  C1  | LABL2 | C2 | LABL3 | C3 |
    +--------+-------+-------+--------+------+-------+----+-------+----+
    |        | LABL4 |  C4   |  etc.  |      |       |    |       |    |
    +--------+-------+-------+--------+------+-------+----+-------+----+
    | AELINK |  10   | INBDA |  OTBDA | -2.0 |       |    |       |    |
    +--------+-------+-------+--------+------+-------+----+-------+----+
    """
    type = 'AELINK'

    @classmethod
    def _init_from_empty(cls):
        aelink_id = 1
        label = 'ELEV'
        independent_labels = ['ELEV1', 'ELEV2']
        linking_coefficents = [1., 2.]
        return AELINK(aelink_id, label, independent_labels, linking_coefficents, comment='')

    def __init__(self, aelink_id: Union[int, str],
                 label: str, independent_labels: List[str],
                 linking_coefficents: List[float],
                 comment: str='') -> None:
        """
        Creates an AELINK card, which defines an equation linking
        AESTAT and AESURF cards

        Parameters
        ----------
        aelink_id : int/str
            unique id
        label : str
            name of the dependent AESURF card
        independent_labels : List[str, ..., str]
            name for the independent variables (AESTATs)
        linking_coefficents : List[float]
            linking coefficients
        comment : str; default=''
            a comment for the card

        """
        BaseCard.__init__(self)
        if comment:
            self.comment = comment
        #: defines the dependent variable name (string)
        self.label = label

        #: defines the independent variable name (string)
        self.independent_labels = independent_labels

        #: linking coefficients (real)
        self.linking_coefficents = linking_coefficents

        if isinstance(aelink_id, str):
            if aelink_id != 'ALWAYS':
                raise RuntimeError("The only valid ID that is a string is 'ALWAYS'")
            aelink_id = 0
        #: an ID=0 is applicable to the global subcase, ID=1 only subcase 1
        self.aelink_id = aelink_id

    @property
    def Cis(self):
        self.deprecated('Cis', 'linking_coefficents', '1.2')
        return self.linking_coefficents
    @Cis.setter
    def Cis(self, linking_coefficents):
        self.deprecated('Cis', 'linking_coefficents', '1.2')
        self.linking_coefficents = linking_coefficents

    def object_attributes(self, mode='public', keys_to_skip=None,
                          filter_properties=False):
        """.. seealso:: `pyNastran.utils.object_methods(...)`"""
        if keys_to_skip is None:
            keys_to_skip = []
        my_keys_to_skip = ['Cis']
        return BaseCard.object_attributes(self, mode=mode,
                                          keys_to_skip=keys_to_skip+my_keys_to_skip,
                                          filter_properties=filter_properties)

    def object_methods(self, mode='public', keys_to_skip=None):
        """.. seealso:: `pyNastran.utils.object_methods(...)`"""
        if keys_to_skip is None:
            keys_to_skip = []
        my_keys_to_skip = ['Cis']
        return BaseCard.object_methods(self, mode=mode, keys_to_skip=keys_to_skip+my_keys_to_skip)

    def validate(self):
        if len(self.independent_labels) != len(self.linking_coefficents):
            msg = 'nlabels=%s nci=%s\nindependent_labels=%s linking_coefficents=%s\n%s' % (
                len(self.independent_labels), len(self.linking_coefficents),
                self.independent_labels, self.linking_coefficents, str(self))
            raise RuntimeError(msg)
        if len(self.independent_labels) == 0:
            msg = 'nlabels=%s nci=%s\nindependent_labels=%s linking_coefficents=%s\n%s' % (
                len(self.independent_labels), len(self.linking_coefficents),
                self.independent_labels, self.linking_coefficents, str(self))
            raise RuntimeError(msg)

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds an AELINK card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        aelink_id = integer_or_string(card, 1, 'ID')
        label = string(card, 2, 'label')
        independent_labels = []
        linking_coefficents = []

        list_fields = [interpret_value(field, card) for field in card[3:]]
        assert len(list_fields) % 2 == 0, 'list_fields=%s' % list_fields
        for i in range(0, len(list_fields), 2):
            independent_label = list_fields[i]
            linking_coefficent = list_fields[i + 1]
            independent_labels.append(independent_label)
            linking_coefficents.append(linking_coefficent)
        return AELINK(aelink_id, label, independent_labels, linking_coefficents,
                      comment=comment)

    #def uncross_reference(self) -> None:
        #"""Removes cross-reference links"""
        #pass

    def raw_fields(self):
        """
        Gets the fields in their unmodified form

        Returns
        -------
        list_fields : List[int/float/str]
            the fields that define the card

        """
        list_fields = ['AELINK', self.aelink_id, self.label]
        for (ivar, ival) in zip(self.independent_labels, self.linking_coefficents):
            list_fields += [ivar, ival]
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        """
        The writer method used by BDF.write_card()

        Parameters
        -----------
        size : int; default=8
            the size of the card (8/16)

        """
        card = self.raw_fields()
        return self.comment + print_card_8(card)


class AELIST(BaseCard):
    """
    Defines a list of aerodynamic elements to undergo the motion prescribed
    with the AESURF Bulk Data entry for static aeroelasticity.

    +---------+------+------+------+------+------+------+------+------+
    |    1    |   2  |   3  |  4   |   5  |   6  |  7   |  8   |  9   |
    +=========+======+======+======+======+======+======+======+======+
    |  AELIST |  SID |  E1  |  E2  |  E3  |  E4  |  E5  |  E6  |  E7  |
    +---------+------+------+------+------+------+------+------+------+
    |         |  E8  | etc. |      |      |      |      |      |      |
    +---------+------+------+------+------+------+------+------+------+
    |  AELIST |  75  | 1001 | THRU | 1075 | 1101 | THRU | 1109 | 1201 |
    +---------+------+------+------+------+------+------+------+------+
    |         | 1202 |      |      |      |      |      |      |      |
    +---------+------+------+------+------+------+------+------+------+

    Notes
    -----
    1. These entries are referenced by the AESURF entry.
    2. When the THRU option is used, all intermediate grid points must exist.
       The word THRU may not appear in field 3 or 9 (2 or 9 for continuations).
    3. Intervening blank fields are not allowed.

    """
    type = 'AELIST'

    @classmethod
    def _init_from_empty(cls):
        return AELIST(1, [1], comment='')

    def __init__(self, sid: int, elements: List[int], comment: str='') -> None:
        """
        Creates an AELIST card, which defines the aero boxes for
        an AESURF/SPLINEx.

        Parameters
        ----------
        sid : int
            unique id
        elements : List[int, ..., int]; int
            list of box ids
        comment : str; default=''
            a comment for the card

        """
        BaseCard.__init__(self)
        if comment:
            self.comment = comment

        if isinstance(elements, integer_types):
            elements = [elements]

        if isinstance(elements, np.ndarray):
            assert len(elements.shape) == 1, elements.shape
            elements = elements.tolist()
        if not isinstance(elements, (list, tuple)):
            raise TypeError('AELIST; type(elements)=%s and must be a list/tuple' % type(elements))

        #: Set identification number. (Integer > 0)
        self.sid = sid
        #: List of aerodynamic boxes generated by CAERO1 entries to define a
        #: surface. (Integer > 0 or 'THRU')
        self.elements = expand_thru(elements)
        self.elements.sort()

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds an AELIST card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        sid = integer(card, 1, 'sid')
        elements = fields(integer_or_string, card, 'eid', i=2, j=len(card))
        return AELIST(sid, elements, comment=comment)

    def cross_reference(self, model: BDF) -> None:
        pass

    def safe_cross_reference(self, model):
        pass

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        pass

    def clean_ids(self):
        self.elements = list(set(self.elements))
        self.elements.sort()

    def raw_fields(self):
        """
        Gets the fields in their unmodified form

        Returns
        -------
        fields : List[int/float/str]
            the fields that define the card

        """
        list_fields = ['AELIST', self.sid] + self.elements
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class AEPARM(BaseCard):
    """
    Defines a general aerodynamic trim variable degree-of-freedom (aerodynamic
    extra point). The forces associated with this controller will be derived
    from AEDW, AEFORCE and AEPRESS input data.

    +--------+----+--------+-------+
    |    1   | 2  |    3   |   4   |
    +========+====+========+=======+
    | AEPARM | ID | LABEL  | UNITS |
    +--------+----+--------+-------+
    | AEPARM | 5  | THRUST | LBS   |
    +--------+----+--------+-------+
    """
    type = 'AEPARM'
    _field_map = {
        1: 'id', 2:'label', 3:'units'
    }

    @classmethod
    def _init_from_empty(cls):
        aeparm_id = 1
        label = 'name'
        units = ''
        return AEPARM(aeparm_id, label, units, comment='')

    def __init__(self, aeparm_id: int, label: str, units: str, comment: str='') -> None:
        """
        Creates an AEPARM card, which defines a new trim variable.

        Parameters
        ----------
        id : int
            the unique id
        label : str
            the variable name
        units : str
            unused by Nastran
        comment : str; default=''
            a comment for the card

        """
        BaseCard.__init__(self)
        if comment:
            self.comment = comment
        self.aeparm_id = aeparm_id
        self.label = label
        self.units = units

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds an AEPARM card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        aeparm_id = integer(card, 1, 'aeparm_id')
        label = string(card, 2, 'label')
        units = card.field(3)
        units = '' if units is None else units

        assert len(card) <= 4, 'len(AEPARM card) = %i\ncard=%s' % (len(card), card)
        return AEPARM(aeparm_id, label, units, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        """
        Adds an AEPARM card from the OP2

        Parameters
        ----------
        data : List[varies]
            a list of fields defined in OP2 format
        comment : str; default=''
            a comment for the card

        """
        aeparm_id = data[0]
        label = data[1]
        units = data[2]
        assert len(data) == 3, 'data = %s' % data
        return AEPARM(aeparm_id, label, units, comment=comment)

    def cross_reference(self, model: BDF) -> None:
        pass

    def safe_cross_reference(self, model):
        pass

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        pass

    def raw_fields(self):
        """
        Gets the fields in their unmodified form

        Returns
        -------
        fields : List[int/float/str]
            the fields that define the card

        """
        list_fields = ['AEPARM', self.aeparm_id, self.label, self.units]
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.raw_fields()
        return self.comment + print_card_8(card)


class AESURF(BaseCard):
    """
    Specifies an aerodynamic control surface as a member of the set of
    aerodynamic extra points. The forces associated with this controller will
    be derived from rigid rotation of the aerodynamic model about the hinge
    line(s) and from AEDW, AEFORCE and AEPRESS input data. The mass properties
    of the control surface can be specified using an AESURFS entry.

    +--------+--------+-------+-------+-------+--------+--------+--------+--------+
    |    1   |   2    |   3   |   4   |   5   |   6    |    7   |   8    |   9    |
    +========+========+=======+=======+=======+========+========+========+========+
    | AESURF |   ID   | LABEL | CID1  | ALID1 |  CID2  | ALID2  |  EFF   |  LDW   |
    +--------+--------+-------+-------+-------+--------+--------+--------+--------+
    |        |  CREFC | CREFS | PLLIM | PULIM | HMLLIM | HMULIM | TQLLIM | TQULIM |
    +--------+--------+-------+-------+-------+--------+--------+--------+--------+
    """
    type = 'AESURF'
    _field_map = {
        1: 'aesid', 2:'label', 3:'cid1', 4:'alid1', 5:'cid2', 6:'alid2',
        7:'eff', 8:'ldw', 9:'crefc', 10:'crefs', 11:'pllim', 12:'pulim',
        13:'hmllim', 14:'hmulim', 15:'tqllim', '16':'tqulim',
    }

    @classmethod
    def _init_from_empty(cls):
        aesid = 1
        label = 'name'
        cid1 = 1
        alid1 = 1
        return AESURF(aesid, label, cid1, alid1,
                      cid2=None, alid2=None, eff=1.0, ldw='LDW',
                      crefc=1.0, crefs=1.0, pllim=-np.pi/2., pulim=np.pi/2.,
                      hmllim=None, hmulim=None, tqllim=None, tqulim=None, comment='')

    def __init__(self, aesid: int, label: str, cid1: int, alid1: int,
                 cid2: Optional[int]=None, alid2: Optional[int]=None,
                 eff: float=1.0, ldw: str='LDW',
                 crefc: float=1.0, crefs: float=1.0,
                 pllim: float=-np.pi/2., pulim: float=np.pi/2.,
                  # hinge moment lower/upper limits
                 hmllim: Optional[int]=None, hmulim: Optional[int]=None,
                  # TABLEDi deflection limits vs. dynamic pressure
                 tqllim: Optional[int]=None, tqulim: Optional[int]=None,
                 comment='') -> None:
        """
        Creates an AESURF card, which defines a control surface

        Parameters
        ----------
        aesid : int
            controller number
        label : str
            controller name
        cid1 / cid2 : int / None
            coordinate system id for primary/secondary control surface
        alid1 / alid2 : int / None
            AELIST id for primary/secondary control surface
        eff : float; default=1.0
            Control surface effectiveness
        ldw : str; default='LDW'
            Linear downwash flag;  ['LDW', 'NODLW']
        crefc : float; default=1.0
            reference chord for the control surface
        crefs : float; default=1.0
            reference area for the control surface
        pllim / pulim : float; default=-pi/2 / pi/2
            Lower/Upper deflection limits for the control surface in radians
        hmllim / hmulim : float; default=None
            Lower/Upper hinge moment limits for the control surface in
            force-length units
        tqllim / tqulim : int; default=None
            Set identification numbers of TABLEDi entries that provide the
            lower/upper deflection limits for the control surface as a
            function of the dynamic pressure
        comment : str; default=''
            a comment for the card

        """
        BaseCard.__init__(self)
        if comment:
            self.comment = comment
        #: Controller identification number
        self.aesid = aesid

        #: Controller name.
        self.label = label

        #: Identification number of a rectangular coordinate system with a
        #: y-axis that defines the hinge line of the control surface
        #: component.
        self.cid1 = cid1
        #: Identification of an AELIST Bulk Data entry that identifies all
        #: aerodynamic elements that make up the control surface
        #: component. (Integer > 0)
        self.alid1 = alid1

        self.cid2 = cid2
        self.alid2 = alid2

        #: Control surface effectiveness. See Remark 4. (Real != 0.0;
        #: Default=1.0)
        self.eff = eff
        #: Linear downwash flag. See Remark 2.
        #: (Character, one of LDW or NOLDW; Default=LDW).
        self.ldw = ldw
        #: Reference chord length for the control surface. (Real>0.0;
        #: Default=1.0)
        self.crefc = crefc
        #: Reference surface area for the control surface. (Real>0.0;
        #: Default=1.0)
        self.crefs = crefs

        #: Lower and upper deflection limits for the control surface in
        #: radians. (Real, Default = +/- pi/2)
        self.pllim = pllim
        self.pulim = pulim

        #: Lower and upper hinge moment limits for the control surface in
        #: force-length units. (Real, Default = no limit) -> 1e8
        self.hmllim = hmllim
        self.hmulim = hmulim

        #: Set identification numbers of TABLEDi entries that provide the
        #: lower and upper deflection limits for the control surface as a
        #: function of the dynamic pressure. (Integer>0, Default = no limit)
        self.tqllim = tqllim
        self.tqulim = tqulim
        self.cid1_ref = None
        self.cid2_ref = None
        self.alid1_ref = None
        self.alid2_ref = None
        self.tqllim_ref = None
        self.tqulim_ref = None
        assert self.ldw in {'LDW', 'NOLDW'}, self.ldw

    #def validate(self):
        #assert self.ldw in {'LDW', 'NOLDW'}, self.ldw

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds an AESURF card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        aesid = integer(card, 1, 'aesid')
        label = string(card, 2, 'label')

        cid1 = integer(card, 3, 'cid1')
        alid1 = integer(card, 4, 'alid1')

        cid2 = integer_or_blank(card, 5, 'cid2')
        alid2 = integer_or_blank(card, 6, 'alid2')

        eff = double_or_blank(card, 7, 'eff', 1.0)
        ldw = string_or_blank(card, 8, 'ldw', 'LDW')
        crefc = double_or_blank(card, 9, 'crefc', 1.0)
        crefs = double_or_blank(card, 10, 'crefs', 1.0)

        pllim = double_or_blank(card, 11, 'pllim', -np.pi / 2.)
        pulim = double_or_blank(card, 12, 'pulim', np.pi / 2.)

        hmllim = double_or_blank(card, 13, 'hmllim')
        hmulim = double_or_blank(card, 14, 'hmulim')
        tqllim = integer_or_blank(card, 15, 'tqllim')
        tqulim = integer_or_blank(card, 16, 'tqulim')
        assert len(card) <= 17, 'len(AESURF card) = %i\ncard=%s' % (len(card), card)
        return AESURF(aesid, label, cid1, alid1, cid2, alid2, eff, ldw,
                      crefc, crefs, pllim, pulim, hmllim, hmulim,
                      tqllim, tqulim, comment=comment)

    def Cid1(self):
        if self.cid1_ref is not None:
            return self.cid1_ref.cid
        return self.cid1

    def Cid2(self):
        if self.cid2_ref is not None:
            return self.cid2_ref.cid
        return self.cid2

    def aelist_id1(self):
        if self.alid1_ref is not None:
            return self.alid1_ref.sid
        return self.alid1

    def aelist_id2(self):
        if self.alid2_ref is not None:
            return self.alid2_ref.sid
        return self.alid2

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        msg = ', which is required by AESURF eid=%s' % self.label
        self.cid1_ref = model.Coord(self.cid1, msg=msg)
        if self.cid2 is not None:
            self.cid2_ref = model.Coord(self.cid2)
        self.alid1_ref = model.AELIST(self.alid1)
        if self.alid2:
            self.alid2_ref = model.AELIST(self.alid2)
        if self.tqllim is not None:
            self.tqllim_ref = model.TableD(self.tqllim)
        if self.tqulim is not None:
            self.tqulim_ref = model.TableD(self.tqulim)

    def safe_cross_reference(self, model, xref_errors):
        msg = ', which is required by AESURF aesid=%s' % self.aesid
        self.cid1_ref = model.safe_coord(self.cid1, self.aesid, xref_errors, msg=msg)
        if self.cid2 is not None:
            self.cid2_ref = model.safe_coord(self.cid2, self.aesid, xref_errors, msg=msg)

        self.alid1_ref = model.safe_aelist(self.alid1, self.aesid, xref_errors, msg=msg)
        if self.alid2:
            self.alid2_ref = model.safe_aelist(self.alid2, self.aesid, xref_errors, msg=msg)

        if self.tqllim is not None:
            self.tqllim_ref = model.safe_tabled(self.tqllim, self.aesid, xref_errors, msg=msg)
        if self.tqulim is not None:
            self.tqulim_ref = model.safe_tabled(self.tqulim, self.aesid, xref_errors, msg=msg)

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.cid1 = self.Cid1()
        self.cid2 = self.Cid2()
        self.cid1_ref = None
        self.cid2_ref = None

        self.alid1 = self.aelist_id1()
        self.alid2 = self.aelist_id2()
        self.alid1_ref = None
        self.alid2_ref = None
        #self.tqulim
        #self.tqllim
        self.tqllim_ref = None
        self.tqulim_ref = None

    def update(self, unused_model, maps):
        coord_map = maps['coord']
        aelist_map = maps['aelist']
        self.cid1 = coord_map[self.cid1]
        if self.cid2:
            self.cid2 = coord_map[self.cid2]

        self.alid1 = aelist_map[self.alid1]
        if self.alid2:
            self.alid2 = aelist_map[self.alid2]

    def raw_fields(self):
        """
        Gets the fields in their unmodified form

        Returns
        -------
        fieldsreset_camera[int/float/str]
            the fields that define the card

        """
        list_fields = ['AESURF', self.aesid, self.label, self.Cid1(), self.aelist_id1(),
                       self.Cid2(), self.aelist_id2(), self.eff, self.ldw,
                       self.crefc, self.crefs, self.pllim, self.pulim, self.hmllim,
                       self.hmulim, self.tqllim, self.tqulim]
        return list_fields

    def repr_fields(self):
        """
        Gets the fields in their simplified form

        Returns
        -------
        fields : List[int/float/str]
            the fields that define the card

        """
        eff = set_blank_if_default(self.eff, 1.0)
        ldw = set_blank_if_default(self.ldw, 'LDW')
        crefc = set_blank_if_default(self.crefc, 1.0)
        crefs = set_blank_if_default(self.crefs, 1.0)

        pllim = set_blank_if_default(self.pllim, -np.pi / 2.)
        pulim = set_blank_if_default(self.pulim, np.pi / 2.)

        list_fields = ['AESURF', self.aesid, self.label, self.Cid1(), self.aelist_id1(),
                       self.Cid2(), self.aelist_id2(), eff, ldw, crefc, crefs,
                       pllim, pulim, self.hmllim, self.hmulim, self.tqllim,
                       self.tqulim]
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
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
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class AESURFS(BaseCard):
    """
    Optional specification of the structural nodes associated with an
    aerodynamic control surface that has been defined on an AESURF entry. The
    mass associated with these structural nodes define the control surface
    moment(s) of inertia about the hinge line(s).
    Specifies rigid body motions to be used as trim variables in static
    aeroelasticity.

    +---------+------+-------+---+-------+---+-------+
    |    1    |  2   |   3   | 4 |   5   | 6 |   7   |
    +=========+======+=======+===+=======+===+=======+
    | AESURFS |  ID  | LABEL |   | LIST1 |   | LIST2 |
    +---------+------+-------+---+-------+---+-------+
    | AESURFS | 6001 | ELEV  |   | 6002  |   | 6003  |
    +---------+------+-------+---+-------+---+-------+
    """
    type = 'AESURFS'

    @classmethod
    def _init_from_empty(cls):
        aesid = 1
        label = 'name'
        list1 = 1
        list2 = 2
        return AESURFS(aesid, label, list1, list2, comment='')

    def __init__(self, aesid: int, label: str,
                 list1: List[int], list2: List[int],
                 comment: str='') -> None:
        """
        Creates an AESURFS card

        Parameters
        ----------
        aesid : int
            the unique id
        label : str
            the AESURF name
        list1 / list2 : int / None
            the list (SET1) of node ids for the primary/secondary
            control surface(s) on the AESURF card
        comment : str; default=''
            a comment for the card

        """
        BaseCard.__init__(self)
        if comment:
            self.comment = comment
        self.aesid = aesid
        self.label = label
        self.list1 = list1
        self.list2 = list2
        self.list1_ref = None
        self.list2_ref = None

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds an AESURFS card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        aesid = integer(card, 1, 'ID')
        label = string(card, 2, 'label')
        list1 = integer(card, 4, 'list1')
        list2 = integer(card, 6, 'list2')
        assert len(card) <= 7, 'len(AESURFS card) = %i\ncard=%s' % (len(card), card)
        return AESURFS(aesid, label, list1, list2, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        aesid = data[0]
        label = data[1]
        list1 = data[2]
        list2 = data[3]
        assert len(data) == 4, 'data = %s' % data
        return AESURFS(aesid, label, list1, list2, comment=comment)

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        msg = ', which is required by AESURFS aesid=%s' % self.aesid
        self.list1_ref = model.Set(self.list1, msg)
        self.list1_ref.cross_reference_set(model, 'Node', msg)

        self.list2_ref = model.Set(self.list1, msg=msg)
        self.list2_ref.cross_reference_set(model, 'Node', msg)

    def safe_cross_reference(self, model):
        msg = ', which is required by AESURFS aesid=%s' % self.aesid
        try:
            self.list1_ref = model.Set(self.list1, msg=msg)
            self.list1_ref.cross_reference_set(model, 'Node', msg)
        except KeyError:
            pass

        try:
            self.list2_ref = model.Set(self.list1, msg=msg)
            self.list2_ref.cross_reference_set(model, 'Node', msg)
        except KeyError:
            pass

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.list1 = self.List1()
        self.list2 = self.List2()
        self.list1_ref = None
        self.list2_ref = None

    def List1(self):
        if self.list1_ref is not None:
            return self.list1_ref.sid
        return self.list1

    def List2(self):
        if self.list2_ref is not None:
            return self.list2_ref.sid
        return self.list2

    def raw_fields(self):
        """
        Gets the fields in their unmodified form

        Returns
        -------
        fields : List[int/float/str]
            the fields that define the card

        """
        list_fields = ['AESURFS', self.aesid, self.label, None, self.List1(), None,
                       self.List2()]
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.raw_fields()
        return self.comment + print_card_8(card)


class CAERO1(BaseCard):
    """
    Defines an aerodynamic macro element (panel) in terms of two leading edge
    locations and side chords. This is used for Doublet-Lattice theory for
    subsonic aerodynamics and the ZONA51 theory for supersonic aerodynamics.

    +--------+-----+-----+----+-------+--------+--------+--------+------+
    |   1    |  2  |  3  | 4  |   5   |   6    |    7   |   8    |   9  |
    +========+=====+=====+====+=======+========+========+========+======+
    | CAERO1 | EID | PID | CP | NSPAN | NCHORD |  LSPAN | LCHORD | IGID |
    +--------+-----+-----+----+-------+--------+--------+--------+------+
    |        |  X1 | Y1  | Z1 |  X12  |   X4   |   Y4   |   Z4   | X43  |
    +--------+-----+-----+----+-------+--------+--------+--------+------+

    ::

      1
      | \
      |   \
      |     \
      |      4
      |      |
      |      |
      2------3

    Attributes
    ----------
    eid : int
        element id
    pid : int, PAERO1
        int : PAERO1 ID
        PAERO1 : PAERO1 object (xref)
    igroup : int
        Group number
    p1 : (1, 3) ndarray float
        xyz location of point 1 (leading edge; inboard)
    p4 : (1, 3) ndarray float
        xyz location of point 4 (leading edge; outboard)
    x12 : float
        distance along the flow direction from node 1 to node 2; (typically x, root chord)
    x43 : float
        distance along the flow direction from node 4 to node 3; (typically x, tip chord)
    cp : int, CORDx
        int : coordinate system
        CORDx : Coordinate object (xref)
    nspan : int
        int > 0 : N spanwise boxes distributed evenly
        int = 0 : use lchord
    nchord : int
        int > 0 : N chordwise boxes distributed evenly
        int = 0 : use lchord
    lspan : int, AEFACT
        int > 0 : AEFACT reference for non-uniform nspan
        int = 0 : use nspan
        AEFACT : AEFACT object  (xref)
    lchord : int, AEFACT
        int > 0 : AEFACT reference for non-uniform nchord
        int = 0 : use nchord
        AEFACT : AEFACT object  (xref)
    comment : str; default=''
         a comment for the card

    """
    type = 'CAERO1'
    _field_map = {
        1: 'sid', 2:'pid', 3:'cp', 4:'nspan', 5:'nchord',
        6:'lspan', 7:'lchord', 8:'igroup', 12:'x12', 16:'x43',
    }
    _properties = ['_field_map', 'shape', 'xy', 'min_max_eid', 'npanels']
    def _get_field_helper(self, n):
        """
        Gets complicated parameters on the CAERO1 card

        Parameters
        ----------
        n : int
            the field number to update
        value : int/float
            the value for the appropriate field

        """
        if n == 9:
            out = self.p1[0]
        elif n == 10:
            out = self.p1[1]
        elif n == 11:
            out = self.p1[2]

        elif n == 13:
            out = self.p4[0]
        elif n == 14:
            out = self.p4[1]
        elif n == 15:
            out = self.p4[2]
        else:
            raise KeyError('Field %r is an invalid CAERO1 entry.' % n)
        return out

    def _update_field_helper(self, n, value):
        """
        Updates complicated parameters on the CAERO1 card

        Parameters
        ----------
        n : int
            the field number to update
        value : int/float
            the value for the appropriate field

        """
        if n == 9:
            self.p1[0] = value
        elif n == 10:
            self.p1[1] = value
        elif n == 11:
            self.p1[2] = value

        elif n == 13:
            self.p4[0] = value
        elif n == 14:
            self.p4[1] = value
        elif n == 15:
            self.p4[2] = value
        else:
            raise KeyError('Field %r=%r is an invalid CAERO1 entry.' % (n, value))

    @classmethod
    def _init_from_empty(cls):
        eid = 1
        pid = 1
        igroup = 1
        p1 = [0., 0., 0.]
        x12 = 1.
        p4 = [0., 10., 0.]
        x43 = 0.5
        return CAERO1(eid, pid, igroup, p1, x12, p4, x43,
                      cp=0, nspan=0, lspan=0, nchord=0, lchord=0, comment='')

    def _finalize_hdf5(self, encoding):
        """hdf5 helper function"""
        self.p1 = np.asarray(self.p1)
        self.p4 = np.asarray(self.p4)

    def __init__(self, eid, pid, igroup, p1, x12, p4, x43,
                 cp=0, nspan=0, lspan=0, nchord=0, lchord=0, comment=''):
        """
        Defines a CAERO1 card, which defines a simplified lifting surface
        (e.g., wing/tail).

        Parameters
        ----------
        eid : int
            element id
        pid : int, PAERO1
            int : PAERO1 ID
            PAERO1 : PAERO1 object (xref)
        igroup : int
            Group number
        p1 : (1, 3) ndarray float
            xyz location of point 1 (leading edge; inboard)
        p4 : (1, 3) ndarray float
            xyz location of point 4 (leading edge; outboard)
        x12 : float
            distance along the flow direction from node 1 to node 2; (typically x, root chord)
        x43 : float
            distance along the flow direction from node 4 to node 3; (typically x, tip chord)
        cp : int, CORDx; default=0
            int : coordinate system
            CORDx : Coordinate object (xref)
        nspan : int; default=0
            int > 0 : N spanwise boxes distributed evenly
            int = 0 : use lchord
        nchord : int; default=0
            int > 0 : N chordwise boxes distributed evenly
            int = 0 : use lchord
        lspan : int, AEFACT; default=0
            int > 0 : AEFACT reference for non-uniform nspan
            int = 0 : use nspan
            AEFACT : AEFACT object  (xref)
        lchord : int, AEFACT; default=0
            int > 0 : AEFACT reference for non-uniform nchord
            int = 0 : use nchord
            AEFACT : AEFACT object  (xref)
        comment : str; default=''
             a comment for the card

        """
        BaseCard.__init__(self)
        if cp is None:
            cp = 0
        if lspan is None:
            lspan = 0
        if nspan is None:
            nspan = 0
        if nchord is None:
            nchord = 0
        if lchord is None:
            lchord = 0
        p1 = np.asarray(p1)
        p4 = np.asarray(p4)

        if comment:
            self.comment = comment
        #: Element identification number
        self.eid = eid

        #: Property identification number of a PAERO2 entry.
        self.pid = pid

        #: Coordinate system for locating point 1.
        self.cp = cp
        self.nspan = nspan
        self.lspan = lspan
        self.nchord = nchord
        self.lchord = lchord
        self.igroup = igroup
        self.p1 = p1
        self.x12 = x12
        self.p4 = p4
        self.x43 = x43
        self.pid_ref = None
        self.cp_ref = None
        self.lchord_ref = None
        self.lspan_ref = None
        self.ascid_ref = None
        self.box_ids = None
        #self._init_ids() #TODO: make this work here?

    def validate(self):
        msg = ''
        is_failed = False
        if not isinstance(self.p1, np.ndarray):
            msg += 'p1=%s and must be a numpy array\n' % (self.p1)
            is_failed = True
        if not isinstance(self.p4, np.ndarray):
            msg += 'p1=%s and must be a numpy array\n' % (self.p1)
            is_failed = True

        if self.x12 <= 0.:
            msg += 'X12=%s and must be greater than or equal to 0\n' % (self.x12)
            is_failed = True
        #if self.x43 <= 0.:
            #msg += 'X43=%s and must be greater than or equal to 0\n' % (self.x43)
            #is_failed = True

        if self.nspan == 0 and self.lspan == 0:
            msg += 'NSPAN or LSPAN must be greater than 0; nspan=%r nlspan=%s\n' % (
                self.nspan, self.lspan)
            is_failed = True
        if self.nspan != 0 and self.lspan != 0:
            msg += 'Either NSPAN or LSPAN must 0; nspan=%r nlspan=%s\n' % (
                self.nspan, self.lspan)
            is_failed = True

        if self.nchord == 0 and self.lchord == 0:
            msg += 'NCHORD or LCHORD must be greater than 0; nchord=%r lchord=%s\n' % (
                self.nchord, self.lchord)
            is_failed = True
        if self.nchord != 0 and self.lchord != 0:
            msg += 'Either NCHORD or LCHORD must 0; nchord=%r lchord=%s\n' % (
                self.nchord, self.lchord)
            is_failed = True
        if is_failed:
            msg += str(self)
            raise ValueError(msg)
        assert len(self.p1) == 3, 'p1=%s' % self.p1
        assert len(self.p4) == 3, 'p4=%s' % self.p4

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a CAERO1 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        eid = integer(card, 1, 'eid')
        pid = integer(card, 2, 'pid')
        cp = integer_or_blank(card, 3, 'cp', 0)
        nspan = integer_or_blank(card, 4, 'nspan', 0)
        nchord = integer_or_blank(card, 5, 'nchord', 0)
        lspan = integer_or_blank(card, 6, 'lspan', 0)
        lchord = integer_or_blank(card, 7, 'lchord', 0)
        igroup = integer(card, 8, 'igid')

        p1 = np.array([
            double_or_blank(card, 9, 'x1', 0.0),
            double_or_blank(card, 10, 'y1', 0.0),
            double_or_blank(card, 11, 'z1', 0.0)])
        x12 = double_or_blank(card, 12, 'x12', 0.)

        p4 = np.array([
            double_or_blank(card, 13, 'x4', 0.0),
            double_or_blank(card, 14, 'y4', 0.0),
            double_or_blank(card, 15, 'z4', 0.0)])
        x43 = double_or_blank(card, 16, 'x43', 0.)

        assert len(card) <= 17, 'len(CAERO1 card) = %i\ncard=%s' % (len(card), card)
        return CAERO1(eid, pid, igroup, p1, x12, p4, x43,
                      cp=cp, nspan=nspan, lspan=lspan, nchord=nchord, lchord=lchord,
                      comment=comment)

    @classmethod
    def add_quad(cls, eid, pid, span, chord, igroup,
                 p1, p2, p3, p4, cp=0, spanwise='y', comment=''):
        r"""
        ::

          1
          | \
          |   \
          |     \
          |      4
          |      |
          |      |
          2------3

        TODO: CP not handled correctly

        """
        x12 = p2[0] - p1[0]
        x43 = p3[0] - p4[0]
        nspan = 0
        lspan = 0
        nchord = 0
        lchord = 0

        if spanwise.lower() == 'y':
            y41 = p4[1] - p1[1]
            y32 = p3[1] - p2[1]
            dspan = max(y41, y32)
        elif spanwise.lower() == 'z':
            y41 = p4[2] - p1[2]
            y32 = p3[2] - p2[2]
            dspan = max(y41, y32)
        else:
            raise NotImplementedError('spanwise=%r; expected=[y, z]' % spanwise.lower())

        dx = max(x12, x43)
        if isinstance(span, integer_types):
            nspan = span
        elif isinstance(span, AEFACT):
            lspan = span.sid
        elif isinstance(span, float):
            nspan = int(math.ceil(dspan / span))
            if nspan <= 0:
                msg = 'y41=%s y32=%s; dspan=%s span=%s nspan=%s; nspan must be greater than 0' % (
                    y41, y32, dspan, span, nspan)
                raise ValueError(msg)
        else:
            raise TypeError(span)

        if isinstance(chord, integer_types):
            nchord = chord
        elif isinstance(chord, AEFACT):
            lchord = chord.sid
        elif isinstance(chord, float):
            nchord = int(math.ceil(dx / chord))
            if nchord <= 0:
                msg = 'x12=%s x43=%s; dx=%s chord=%s nchord=%s; nchord must be greater than 0' % (
                    x12, x43, dx, chord, nchord)
                raise ValueError(msg)
        else:
            raise TypeError(chord)

        return CAERO1(eid, pid, igroup, p1, x12, p4, x43,
                      cp=cp, nspan=nspan, lspan=lspan, nchord=nchord, lchord=lchord,
                      comment=comment)

    def flip_normal(self):
        """flips the CAERO1 normal vector"""
        self.p1, self.p4 = self.p4, self.p1
        self.x12, self.x43 = self.x43, self.x12

    def _init_ids(self, dtype='int32'):
        """
        Fill `self.box_ids` with the sub-box ids. Shape is (nchord, nspan)

        """
        nchord, nspan = self.shape
        assert nchord >= 1, 'nchord=%s' % nchord
        assert nspan >= 1, 'nspan=%s' % nspan
        self.box_ids = np.zeros((nchord, nspan), dtype=dtype)

        npanels = nchord * nspan
        try:
            self.box_ids = np.arange(self.eid, self.eid + npanels,
                                     dtype=dtype).reshape(nspan, nchord).T
        except OverflowError:
            if dtype == 'int64':
                # we already tried int64
                msg = 'eid=%s lchord=%s lspan=%s nchord=%s' % (
                    self.eid, self.lchord, self.lspan, nchord)
                raise OverflowError(msg)
            self._init_ids(dtype='int64')

    def Cp(self):
        if self.cp_ref is not None:
            return self.cp_ref.cid
        return self.cp

    def Pid(self):
        if self.pid_ref is not None:
            return self.pid_ref.pid
        return self.pid

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        msg = ', which is required by CAERO1 eid=%s' % self.eid
        self.pid_ref = model.PAero(self.pid, msg=msg)
        self.cp_ref = model.Coord(self.cp, msg=msg)
        if model.sol in [144, 145, 146, 200]:
            self.ascid_ref = model.Acsid(msg=msg)
        else:
            self.ascid_ref = model.safe_acsid(msg=msg)

        if self.nchord == 0:
            assert isinstance(self.lchord, integer_types), self.lchord
            self.lchord_ref = model.AEFact(self.lchord, msg)
        if self.nspan == 0:
            assert isinstance(self.lspan, integer_types), self.lspan
            self.lspan_ref = model.AEFact(self.lspan, msg)
        self._init_ids()

    def safe_cross_reference(self, model, xref_errors):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        msg = ', which is required by CAERO1 eid=%s' % self.eid
        try:
            self.pid_ref = model.PAero(self.pid, msg=msg)
        except KeyError:
            pass

        self.cp_ref = model.safe_coord(self.cp, self.eid, xref_errors, msg=msg)
        self.ascid_ref = model.safe_acsid(msg=msg)

        if self.nchord == 0:
            assert isinstance(self.lchord, integer_types), self.lchord
            self.lchord_ref = model.safe_aefact(self.lchord, self.eid, xref_errors, msg)

        if self.nspan == 0:
            assert isinstance(self.lspan, integer_types), self.lspan
            self.lspan_ref = model.safe_aefact(self.lspan, self.eid, xref_errors, msg)

        self._init_ids()

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.pid = self.Pid()
        self.cp = self.Cp()
        self.lchord = self.get_LChord()
        self.lspan = self.get_LSpan()
        self.pid_ref = None
        self.cp_ref = None
        self.lchord_ref = None
        self.lspan_ref = None
        self.ascid_ref = None

    def update(self, maps):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        #msg = ', which is required by CAERO1 eid=%s' % self.eid
        paero_map = maps['paero']
        coord_map = maps['coord']
        aefact_map = maps['aefact']
        self.pid = paero_map[self.pid]
        self.cp = coord_map[self.cp]
        #self.acsid = coord_map[self.acsid]  # AERO/AEROS card

        if self.nchord == 0:
            self.lchord = aefact_map[self.lchord]
        if self.nspan == 0:
            self.lspan = aefact_map[self.lspan]
        #self._init_ids(model)

    @property
    def min_max_eid(self):
        """
        Gets the min and max element ids of the CAERO card

        Returns
        -------
        min_max_eid : (2, ) list
            The [min_eid, max_eid]

        """
        nchord, nspan = self.shape
        return [self.eid, self.eid + nchord * nspan]

    def get_leading_edge_points(self):
        """gets the leading edge points"""
        if self.cp_ref is None and self.cp == 0:
            p1 = self.p1
            p4 = self.p4
        else:
            p1 = self.cp_ref.transform_node_to_global(self.p1)
            p4 = self.cp_ref.transform_node_to_global(self.p4)
        return p1, p4

    def get_points(self):
        """
        Get the 4 corner points for the CAERO card

        Returns
        -------
        p1234 : (4, 3) list
             List of 4 corner points in the global frame

        """
        p1, p4 = self.get_leading_edge_points()
        if self.ascid_ref is None:
            # yes, this really does list + array addition
            p2 = p1 + np.array([self.x12, 0., 0.])
            p3 = p4 + np.array([self.x43, 0., 0.])
        else:
            p2 = p1 + self.ascid_ref.transform_vector_to_global(np.array([self.x12, 0., 0.]))
            p3 = p4 + self.ascid_ref.transform_vector_to_global(np.array([self.x43, 0., 0.]))
        return [p1, p2, p3, p4]

    def get_box_index(self, box_id):
        """
        Get the index of ``self.box_ids`` that coresponds to the given box id.

        Parameters
        -----------
        box_id : int
            Box id to ge tthe index of.

        Returns
        --------
        index : tuple
            Index of ``self.box_ids`` that coresponds to the given box id.

        """
        if box_id not in self.box_ids:
            self._box_id_error(box_id)
        index = np.where(self.box_ids == box_id)
        index = (index[0][0], index[1][0])
        return index

    def get_box_quarter_chord_center(self, box_id):
        """
        The the location of the quarter chord of the box along the centerline.

        Parameters
        -----------
        box_id : int
            Box id.

        Returns
        --------
        xyz_quarter_chord : ndarray
            Location of box quater chord in global.

        """
        return self._get_box_x_chord_center(box_id, 0.25)

    def get_box_mid_chord_center(self, box_id):
        """
        The the location of the mid chord of the box along the centerline.

        Parameters
        -----------
        box_id : int
            Box id.

        Returns
        --------
        xyz_mid_chord : ndarray
            Location of box mid chord in global.

        """
        return self._get_box_x_chord_center(box_id, 0.5)

    def _get_box_x_chord_center(self, box_id, x_chord):
        """
        The the location of the x_chord of the box along the centerline.
        """
        if self.lchord != 0 or self.lspan != 0:
            raise NotImplementedError()
        ichord, ispan = self.get_box_index(box_id)

        le_vector = self.p4 - self.p1
        delta_xyz = le_vector * ((ispan + 0.5)/self.nspan)
        yz = delta_xyz[1:3] + self.p1[1:3]
        chord = ((ispan + 0.5)/self.nspan) * (self.x43 - self.x12) + self.x12
        x = (ichord + x_chord)/self.nchord * chord + self.p1[0] + delta_xyz[0]
        return np.array([x, yz[0], yz[1]])

    def _box_id_error(self, box_id):
        """
        Raise box_id IndexError.
        """
        msg = '%i not in range of aero box ids\nRange: %i to %i' % (box_id, self.box_ids[0, 0],
                                                                    self.box_ids[-1, -1])
        raise IndexError(msg)

    @property
    def npanels(self):
        nchord, nspan = self.shape
        return nchord * nspan

    @property
    def shape(self):
        """returns (nelements_nchord, nelements_span)"""
        if self.nchord == 0:
            x = self.lchord_ref.fractions
            nchord = len(x) - 1
        else:
            nchord = self.nchord

        if self.nspan == 0:
            y = self.lspan_ref.fractions
            nspan = len(y) - 1
        else:
            nspan = self.nspan
        if nchord < 1 or nspan < 1:
            msg = 'CAERO1 eid=%s nchord=%s nspan=%s lchord=%s lspan=%s' % (
                self.eid, self.nchord, self.nspan, self.lchord, self.lspan)
            raise RuntimeError(msg)
        return nchord, nspan

    def get_npanel_points_elements(self):
        """
        Gets the number of sub-points and sub-elements for the CAERO card

        Returns
        -------
        npoints : int
            The number of nodes for the CAERO
        nelmements : int
            The number of elements for the CAERO

        """
        nchord, nspan = self.shape
        nelements = nchord * nspan
        npoints = (nchord + 1) * (nspan + 1)
        return npoints, nelements

    @property
    def xy(self):
        """
        Returns
        -------
        x : (nchord,) ndarray
            The percentage x location in the chord-wise direction of each panel
        y : (nspan,) ndarray
            The percentage y location in the span-wise direction of each panel

        """
        if self.nchord == 0:
            x = self.lchord_ref.fractions
            nchord = len(x) - 1
        else:
            nchord = self.nchord
            x = np.linspace(0., 1., nchord + 1)

        if self.nspan == 0:
            y = self.lspan_ref.fractions
            nspan = len(y) - 1
        else:
            nspan = self.nspan
            y = np.linspace(0., 1., nspan + 1)

        if nchord < 1 or nspan < 1:
            msg = 'CAERO1 eid=%s nchord=%s nspan=%s lchord=%s lspan=%s' % (
                self.eid, self.nchord, self.nspan, self.lchord, self.lspan)
            raise RuntimeError(msg)
        return x, y

    def panel_points_elements(self):
        """
        Gets the sub-points and sub-elements for the CAERO1 card

        Returns
        -------
        points : (nnodes,3) ndarray of floats
            the array of points
        elements : (nelements,4) ndarray of integers
            the array of point ids

        """
        p1, p2, p3, p4 = self.get_points()
        x, y = self.xy
        # We're reordering the points so we get the node ids and element ids
        # to be consistent with Nastran.  This is only useful if you're plotting
        # aero panel forces
        #
        # this gives us chordwise panels and chordwise nodes
        return points_elements_from_quad_points(p1, p4, p3, p2, y, x, dtype='int32')

    def set_points(self, points):
        self.p1 = points[0]
        p2 = points[1]
        p3 = points[2]
        self.p4 = points[3]
        self.x12 = p2[0] - self.p1[0]
        self.x43 = p3[0] - self.p4[0]
        assert self.x12 >= 0., 'p1=%s p2=%s' % (self.p1, p2)
        assert self.x43 >= 0., 'p4=%s p3=%s' % (self.p4, p3)
        assert self.x12 > 0. or self.x43 > 0., 'points=%s' % (points)
        self.p1 = np.asarray(self.p1)
        self.p4 = np.asarray(self.p4)

    def shift(self, dxyz):
        """shifts the aero panel"""
        self.p1 += dxyz
        self.p4 += dxyz

    def raw_fields(self):
        """
        Gets the fields in their unmodified form

        Returns
        -------
        fields : list
          the fields that define the card

        """
        lchord = self.get_LChord()
        lspan = self.get_LSpan()
        list_fields = (['CAERO1', self.eid, self.Pid(), self.Cp(), self.nspan,
                        self.nchord, lspan, lchord, self.igroup, ] +
                       list(self.p1) + [self.x12] + list(self.p4) + [self.x43])
        return list_fields

    def get_LChord(self):
        if self.lchord_ref is not None:
            return self.lchord_ref.sid
        return self.lchord

    def get_LSpan(self):
        if self.lspan_ref is not None:
            return self.lspan_ref.sid
        return self.lspan

    def repr_fields(self):
        """
        Gets the fields in their simplified form

        Returns
        -------
        fields : LIST
            The fields that define the card

        """
        cp = set_blank_if_default(self.Cp(), 0)
        nspan = set_blank_if_default(self.nspan, 0)
        nchord = set_blank_if_default(self.nchord, 0)
        lchord = set_blank_if_default(self.get_LChord(), 0)
        lspan = set_blank_if_default(self.get_LSpan(), 0)
        list_fields = (['CAERO1', self.eid, self.Pid(), cp, nspan, nchord,
                        lspan, lchord, self.igroup] + list(self.p1) +
                       [self.x12] + list(self.p4) + [self.x43])
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class CAERO2(BaseCard):
    """
    Aerodynamic Body Connection
    Defines aerodynamic slender body and interference elements for
    Doublet-Lattice aerodynamics.

    +--------+-----+-----+----+-----+------+-----+------+------+
    |    1   |  2  |  3  |  4 |  5  |   6  |   7 |   8  |  9   |
    +========+=====+=====+====+=====+======+=====+======+======+
    | CAERO2 | EID | PID | CP | NSB | NINT | LSB | LINT | IGID |
    +--------+-----+-----+----+-----+------+-----+------+------+
    |        | X1  |  Y1 | Z1 | X12 |      |     |      |      |
    +--------+-----+-----+----+-----+------+-----+------+------+
    """
    type = 'CAERO2'
    _field_map = {
        1: 'sid', 2:'pid', 3:'cp', 4:'nsb', 5:'lsb',
        6:'nint', 7:'lint', 8:'igroup', 12:'x12',
    }
    _properties = ['nboxes']

    def _get_field_helper(self, n):
        """
        Gets complicated parameters on the CAERO2 card

        Parameters
        ----------
        n : int
            The field number to update

        Returns
        -------
        value : int, float, None
            The value for the appropriate field

        """
        if n == 9:
            out = self.p1[0]
        elif n == 10:
            out = self.p1[1]
        elif n == 11:
            out = self.p1[2]
        else:
            raise KeyError('Field %r is an invalid CAERO2 entry.' % n)
        return out

    def _update_field_helper(self, n, value):
        """
        Updates complicated parameters on the CAERO2 card

        Parameters
        ----------
        n : int
            The field number to update
        value : int, float, None
            The value for the appropriate field

        """
        if n == 9:
            self.p1[0] = value
        elif n == 10:
            self.p1[1] = value
        elif n == 11:
            self.p1[2] = value
        else:
            raise KeyError('Field %r=%r is an invalid CAERO2 entry.' % (n, value))

    @classmethod
    def _init_from_empty(cls):
        eid = 1
        pid = 1
        igroup = 1
        p1 = [0., 0., 0.]
        x12 = 10.
        return CAERO2(eid, pid, igroup, p1, x12, cp=0, nsb=0, nint=0, lsb=0, lint=0, comment='')

    def __init__(self, eid, pid, igroup, p1, x12,
                 cp=0, nsb=0, nint=0, lsb=0, lint=0, comment=''):
        """
        Defines a CAERO2 card, which defines a slender body
        (e.g., fuselage/wingtip tank).

        Parameters
        ----------
        eid : int
            element id
        pid : int, PAERO2
            int : PAERO2 ID
            PAERO2 : PAERO2 object (xref)
        igroup : int
            Group number
        p1 : (1, 3) ndarray float
            xyz location of point 1 (forward position)
        x12 : float
            length of the CAERO2
        cp : int, CORDx; default=0
            int : coordinate system
            CORDx : Coordinate object (xref)
        nsb : int; default=0
            Number of slender body elements
        lsb : int; default=0
            AEFACT id for defining the location of the slender body elements
        nint : int; default=0
            Number of interference elements
        lint : int; default=0
            AEFACT id for defining the location of interference elements
        comment : str; default=''
            a comment for the card

        """
        BaseCard.__init__(self)
        if lsb is None:
            lsb = 0
        if lint is None:
            lint = 0
        if nint is None:
            nint = 0
        if nsb is None:
            nsb = 0

        if comment:
            self.comment = comment
        p1 = np.asarray(p1)

        #: Element identification number
        self.eid = eid

        #: Property identification number of a PAERO2 entry.
        self.pid = pid

        #: Coordinate system for locating point 1.
        self.cp = cp

        #: Number of slender body elements. If NSB > 0, then NSB equal
        #: divisions are assumed; if zero or blank, specify a list of
        #: divisions in LSB. (Integer >= 0)
        self.nsb = nsb

        #: Number of interference elements. If NINT > 0, then NINT equal
        #: divisions are assumed; if zero or blank, specify a list of
        #: divisions in LINT. (Integer >= 0)
        self.nint = nint

        #: ID of an AEFACT Bulk Data entry for slender body division
        #: points; used only if NSB is zero or blank. (Integer >= 0)
        self.lsb = lsb

        #: ID of an AEFACT data entry containing a list of division
        #: points for interference elements; used only if NINT is zero
        #: or blank. (Integer > 0)
        self.lint = lint

        #: Interference group identification. Aerodynamic elements with
        #: different IGIDs are uncoupled. (Integer >= 0)
        self.igroup = igroup

        #: Location of point 1 in coordinate system CP
        self.p1 = p1

        #: Length of body in the x-direction of the aerodynamic coordinate
        #: system.  (Real > 0)
        self.x12 = x12

        self.pid_ref = None
        self.cp_ref = None
        self.lint_ref = None
        self.lsb_ref = None
        self.ascid_ref = None

    def validate(self):
        #print('nsb=%s lsb=%s' % (self.nsb, self.lsb))
        #print('nint=%s lint=%s' % (self.nint, self.lint))
        assert isinstance(self.lsb, integer_types), self.lsb
        assert isinstance(self.lint, integer_types), self.lint
        assert len(self.p1) == 3, 'CAERO2: p1=%s' % self.p1
        if self.nsb == 0 and self.lsb == 0:
            msg = 'CAERO2: nsb=%s lsb=%s; nsb or lsb must be > 0' % (self.nsb, self.lsb)
            raise ValueError(msg)
        if self.nint == 0 and self.lint == 0:
            msg = 'CAERO2: nint=%s lint=%s; nint or lint must be > 0' % (self.nint, self.lint)
            raise ValueError(msg)
        assert len(self.p1) == 3, 'CAERO2: p1=%s' % self.p1

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a CAERO2 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        eid = integer(card, 1, 'eid')
        pid = integer(card, 2, 'pid')
        cp = integer_or_blank(card, 3, 'cp', 0)
        nsb = integer_or_blank(card, 4, 'nsb', 0)
        nint = integer_or_blank(card, 5, 'nint', 0)

        lsb = integer_or_blank(card, 6, 'nsb=%s lsb' % nsb, 0)
        lint = integer_or_blank(card, 7, 'nint=%s lint' % nint, 0)
        igid = integer(card, 8, 'igid')

        p1 = np.array([
            double_or_blank(card, 9, 'x1', 0.0),
            double_or_blank(card, 10, 'y1', 0.0),
            double_or_blank(card, 11, 'z1', 0.0)])
        x12 = double_or_blank(card, 12, 'x12', 0.)
        assert len(card) <= 13, 'len(CAERO2 card) = %i\ncard=%s' % (len(card), card)
        return CAERO2(eid, pid, igid, p1, x12,
                      cp=cp, nsb=nsb, nint=nint, lsb=lsb, lint=lint,
                      comment=comment)

    def Cp(self):
        if self.cp_ref is not None:
            return self.cp_ref.cid
        return self.cp

    def Pid(self):
        if self.pid_ref is not None:
            return self.pid_ref.pid
        return self.pid

    def Lsb(self):  # AEFACT
        if self.lsb_ref is not None:
            return self.lsb_ref.sid
        return self.lsb

    def Lint(self):  # AEFACT
        if self.lint_ref is not None:
            return self.lint_ref.sid
        return self.lint

    @property
    def nboxes(self):
        if self.nsb > 0:
            return self.nsb
        return len(self.lsb_ref.fractions) # AEFACT

    def _init_ids(self, dtype='int32'):
        self.box_ids = np.arange(0, self.nboxes, dtype=dtype)

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        msg = ', which is required by CAERO2 eid=%s' % self.eid
        self.pid_ref = model.PAero(self.pid, msg=msg)  # links to PAERO2
        self.cp_ref = model.Coord(self.cp, msg=msg)
        if self.nsb == 0:
            self.lsb_ref = model.AEFact(self.lsb, msg=msg)
        if self.nint == 0:
            self.lint_ref = model.AEFact(self.lint, msg=msg)
        self.ascid_ref = model.Acsid(msg=msg)
        self._init_ids()

    def safe_cross_reference(self, model, xref_errors):
        msg = ', which is required by CAERO2 eid=%s' % self.eid
        self.pid_ref = model.safe_paero(self.pid, self.eid, xref_errors, msg=msg)  # links to PAERO2

        self.cp_ref = model.safe_coord(self.cp, self.eid, xref_errors, msg=msg)

        if self.nsb == 0:
            self.lsb_ref = model.safe_aefact(self.lsb, self.eid, xref_errors, msg=msg)
        if self.nint == 0:
            self.lint_ref = model.safe_aefact(self.lint, self.eid, xref_errors, msg=msg)
        self.ascid_ref = model.safe_acsid(msg=msg)
        self._init_ids()

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.pid = self.Pid()
        self.cp = self.Cp()
        if self.nsb == 0:
            self.lsb = self.Lsb()
        if self.nint == 0:
            self.lint = self.Lint()
        self.pid_ref = None
        self.cp_ref = None
        self.lint_ref = None
        self.lsb_ref = None
        self.ascid_ref = None

    def get_points(self):
        """creates a 1D representation of the CAERO2"""
        p1 = self.cp_ref.transform_node_to_global(self.p1)
        p2 = p1 + self.ascid_ref.transform_vector_to_global(np.array([self.x12, 0., 0.]))

        #print("x12 = %s" % self.x12)
        #print("pcaero[%s] = %s" % (self.eid, [p1,p2]))
        return [p1, p2]

    def get_points_elements_3d(self):
        """
        Gets the points/elements in 3d space as CQUAD4s
        The idea is that this is used by the GUI to display CAERO panels.

        TODO: doesn't support the aero coordinate system

        """
        paero2 = self.pid_ref

        if self.nsb == 0:
            xstation = self.lsb_ref.fractions
            nx = len(xstation) - 1
            #print('xstation = ', xstation)
        else:
            nx = self.nsb
            station = np.linspace(0., nx, num=nx+1) # *dx?
        assert nx > 0, 'nx=%s' % nx


        #print('paero2 - pid=%s lrsb=%s lrib=%s' % (paero2.pid, paero2.lrsb, paero2.lrib))
        if paero2.lrsb in [0, None]:
            radii_slender = np.ones(nx + 1) * paero2.width
        else:
            radii_slender = paero2.lrsb_ref.fractions

        # TODO: not suppported
        if paero2.lrib in [0, None]:
            unused_radii_interference = np.ones(nx + 1) * paero2.width
        else:
            #print('lrib = ', paero2.lrib)
            unused_radii_interference = paero2.lrib_ref.fractions
        radii = radii_slender

        # TODO: not suppported
        #theta_interference1 = paero2.theta1
        #theta_interference2 = paero2.theta2

        if self.nsb != 0:
            p1, p2 = self.get_points()
            L = p2 - p1
            #print('L=%s nx=%s' % (L, nx))
            dxyz = L / nx
            #print('dxyz\n%s' % (dxyz))
            dx, dy, dz = dxyz
            xstation = station * dx
            ystation = station * dy
            zstation = station * dz
        else:
            p1, p2 = self.get_points()
            L = p2 - p1
            dxi = xstation.max() - xstation.min()

            #print('L=%s nx=%s dxi=%s' % (L, nx, dxi))
            xratio = xstation / dxi
            #print('xstation/dxi=%s' % xratio)
            dxyz = np.zeros((nx+1, 3))
            for i, xr in enumerate(xratio):
                dxyz[i, :] = xr * L
            ystation = dxyz[:, 1]
            zstation = dxyz[:, 2]

        # I think this just lets you know what directions it can pivot in
        # and therefore doesn't affect visualization
        #assert paero2.orient == 'ZY', paero2.orient
        aspect_ratio = paero2.AR

        assert len(radii) == (nx + 1), 'len(radii)=%s nx=%s' % (len(radii), nx)
        if len(xstation) != (nx + 1):
            msg = 'len(xstation)=%s nx=%s\nxstation=%s\n%s' % (
                len(xstation), nx, xstation, str(self))
            raise RuntimeError(msg)

        xyz, elems = create_axisymmetric_body(
            xstation, ystation, zstation, radii, aspect_ratio,
            p1)
        assert xyz is not None, str(self)
        return xyz, elems

    def set_points(self, points):
        self.p1 = np.asarray(points[0])
        p2 = np.asarray(points[1])
        x12 = p2 - self.p1
        self.x12 = x12[0]

    def shift(self, dxyz):
        """shifts the aero panel"""
        self.p1 += dxyz

    def raw_fields(self):
        """
        Gets the fields in their unmodified form

        Returns
        -------
        fields : list
            The fields that define the card

        """
        list_fields = (['CAERO2', self.eid, self.Pid(), self.Cp(), self.nsb,
                        self.nint, self.Lsb(), self.Lint(), self.igroup, ] + list(self.p1)
                       + [self.x12])
        return list_fields

    def repr_fields(self):
        """
        Gets the fields in their simplified form

        Returns
        -------
        fields : list
            The fields that define the card

        """
        cp = set_blank_if_default(self.Cp(), 0)
        nint = set_blank_if_default(self.nint, 0)
        lsb = set_blank_if_default(self.Lsb(), 0)
        lint = set_blank_if_default(self.Lint(), 0)
        list_fields = (['CAERO2', self.eid, self.Pid(), cp, self.nsb, nint,
                        lsb, lint, self.igroup, ] + list(self.p1) +
                       [self.x12])
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class CAERO3(BaseCard):
    type = 'CAERO3'
    _properties = ['shape', 'xy']

    @classmethod
    def _init_from_empty(cls):
        eid = 1
        pid = 1
        list_w = 1
        p1 = [0., 0., 0.]
        p4 = [0., 10., 0.]
        x12 = 10.
        x43 = 10.
        return CAERO3(eid, pid, list_w, p1, x12, p4, x43,
                      cp=0, list_c1=None, list_c2=None, comment='')

    def __init__(self, eid, pid, list_w,
                 p1, x12, p4, x43,
                 cp=0, list_c1=None, list_c2=None,
                 comment=''):
        """
        Creates a CAERO2 card, which defines a wing with a wing break/cant.

        Parameters
        ----------
        eid : int
            element id
        pid : int
            PAERO3 property id
        p1 : (3,) float ndarray
            ???
        x12 : float
            ???
        p4 : (3,) float ndarray
            ???
        x43 : float
            ???
        cp : int; default=0
            coordinate system for locating point 1
        list_w : int
            ???
        list_c1 : int; default=None
            defines an AEFACT for ???
        list_c2 : int; default=None
            defines an AEFACT for ???
        comment : str; default=''
            a comment for the card

        """
        assert cp != 100
        BaseCard.__init__(self)
        if comment:
            self.comment = comment
        p1 = np.asarray(p1)
        p4 = np.asarray(p4)

        #: Element identification number
        self.eid = eid

        #: Property identification number of a PAERO3 entry.
        self.pid = pid

        #: Coordinate system for locating point 1.
        self.cp = cp
        self.list_w = list_w
        self.list_c1 = list_c1
        self.list_c2 = list_c2

        self.p1 = p1
        self.x12 = x12
        self.p4 = p4
        self.x43 = x43
        self.pid_ref = None
        self.cp_ref = None
        self.ascid_ref = None
        self.list_w_ref = None
        self.list_c1_ref = None
        self.list_c2_ref = None

    def validate(self):
        assert len(self.p1) == 3, 'p1=%s' % self.p1
        assert len(self.p4) == 3, 'p4=%s' % self.p4
        assert self.x12 > 0., 'x12=%s' % self.x12
        assert self.x43 >= 0., 'x43=%s' % self.x43
        assert isinstance(self.cp, int), 'cp=%r' % self.cp

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a CAERO3 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        eid = integer(card, 1, 'eid')
        pid = integer(card, 2, 'pid')
        cp = integer_or_blank(card, 3, 'cp', 0)
        list_w = integer(card, 4, 'list_w')
        list_c1 = integer_or_blank(card, 5, 'list_c1')
        list_c2 = integer_or_blank(card, 6, 'list_c2')
        p1 = np.array([
            double_or_blank(card, 9, 'x1', 0.0),
            double_or_blank(card, 10, 'y1', 0.0),
            double_or_blank(card, 11, 'z1', 0.0)])
        x12 = double(card, 12, 'x12')
        p4 = np.array([
            double_or_blank(card, 13, 'x4', 0.0),
            double_or_blank(card, 14, 'y4', 0.0),
            double_or_blank(card, 15, 'z4', 0.0)])
        x43 = double_or_blank(card, 16, 'x43', 0.0)
        assert len(card) <= 17, 'len(CAERO3 card) = %i\ncard=%s' % (len(card), card)
        return CAERO3(eid, pid, list_w, p1, x12, p4, x43,
                      cp=cp, list_c1=list_c1, list_c2=list_c2, comment=comment)

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        msg = ', which is required by CAERO3 eid=%s' % self.eid
        self.pid_ref = model.PAero(self.pid, msg=msg)  # links to PAERO3
        self.cp_ref = model.Coord(self.cp, msg=msg)
        if self.list_w is not None:
            self.list_w_ref = model.AEFact(self.list_w, msg=msg)
        if self.list_c1_ref is not None:
            self.list_c1_ref = model.AEFact(self.list_c1, msg=msg)
        if self.list_c2 is not None:
            self.list_c2_ref = model.AEFact(self.list_c2, msg=msg)
        self.ascid_ref = model.Acsid(msg=msg)

    def safe_cross_reference(self, model, xref_errors):
        msg = ', which is required by CAERO3 eid=%s' % self.eid
        self.pid_ref = model.safe_paero(self.pid, self.eid, xref_errors, msg=msg)  # links to PAERO3
        self.cp_ref = model.safe_coord(self.cp, self.eid, xref_errors, msg=msg)

        if self.list_w is not None:
            self.list_w_ref = model.safe_aefact(self.list_w, self.eid, xref_errors, msg=msg)

        if self.list_c1 is not None:
            self.list_c1_ref = model.safe_aefact(self.list_c1, self.eid, xref_errors, msg=msg)

        if self.list_c2 is not None:
            self.list_c2_ref = model.safe_aefact(self.list_c2, self.eid, xref_errors, msg=msg)
        try:
            self.ascid_ref = model.Acsid(msg=msg)
        except KeyError:
            model.log.warning('cannot find an aero coordinate system for %s' % msg)

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.pid = self.Pid()
        self.cp = self.Cp()
        if self.list_w != self.List_w():
            self.list_w = self.List_w()
        if self.list_c1 != self.List_c1():
            self.list_c1 = self.List_c1()
        if self.list_c2 != self.List_c2():
            self.list_c2 = self.List_c2()

        self.pid_ref = None
        self.cp_ref = None
        self.ascid_ref = None
        self.list_w_ref = None
        self.list_c1_ref = None
        self.list_c2_ref = None

    def get_points(self):
        """
        Get the 4 corner points for the CAERO card

        Returns
        -------
        p1234 : (4, 3) list
             List of 4 corner points in the global frame

        """
        p1 = self.cp_ref.transform_node_to_global(self.p1)
        p4 = self.cp_ref.transform_node_to_global(self.p4)
        p2 = p1 + self.ascid_ref.transform_vector_to_global(np.array([self.x12, 0., 0.]))
        p3 = p4 + self.ascid_ref.transform_vector_to_global(np.array([self.x43, 0., 0.]))
        return [p1, p2, p3, p4]

    def panel_points_elements(self):
        """
        Gets the sub-points and sub-elements for the CAERO card

        Returns
        -------
        points : (nnodes,3) ndarray of floats
            the array of points
        elements : (nelements,4) ndarray of integers
            the array of point ids

        """
        p1, p2, p3, p4 = self.get_points()
        x, y = self.xy
        return points_elements_from_quad_points(p1, p2, p3, p4, x, y, dtype='int32')

    def get_npanel_points_elements(self):
        """
        Gets the number of sub-points and sub-elements for the CAERO card

        Returns
        -------
        npoints : int
            The number of nodes for the CAERO
        nelmements : int
            The number of elements for the CAERO

        """
        nchord, nspan = self.shape
        nelements = nchord * nspan
        npoints = (nchord + 1) * (nspan + 1)
        return npoints, nelements

    @property
    def shape(self):
        """returns (nelements_nchord, nelements_span)"""
        nchord = 2
        nspan = self.pid_ref.nbox
        return nchord, nspan

    @property
    def xy(self):
        """
        Returns
        -------
        x : (nchord,) ndarray
            The percentage x location in the chord-wise direction of each panel
        y : (nspan,) ndarray
            The percentage y location in the span-wise direction of each panel

        """
        nchord, nspan = self.shape
        x = np.linspace(0., 1., nchord + 1)
        y = np.linspace(0., 1., nspan + 1)

        if nchord < 1 or nspan < 1:
            msg = 'CAERO3 eid=%s nchord=%s nspan=%s' % (
                self.eid, nchord, nspan)
            raise RuntimeError(msg)
        return x, y

    #def get_points_elements_3d(self):
        #"""
        #Gets the points/elements in 3d space as CQUAD4s
        #The idea is that this is used by the GUI to display CAERO panels.

        #TODO: doesn't support the aero coordinate system
        #"""
        #paero2 = self.pid_ref

    def Cp(self):
        if self.cp_ref is not None:
            return self.cp_ref.cid
        return self.cp

    def Pid(self):
        if self.pid_ref is not None:
            return self.pid_ref.pid
        return self.pid

    def List_w(self):
        if self.list_w_ref is not None:
            return self.list_w_ref.sid
        return self.list_w

    def List_c1(self):
        if self.list_c1_ref is not None:
            return self.list_c1_ref.sid
        return self.list_c1

    def List_c2(self):
        if self.list_c2_ref is not None:
            return self.list_c2_ref.sid
        return self.list_c2

    def raw_fields(self):
        """
        Gets the fields in their unmodified form

        Returns
        -------
        fields : list
            The fields that define the card

        """
        list_fields = (['CAERO3', self.eid, self.Pid(), self.Cp(), self.List_w(),
                        self.List_c1(), self.List_c2(), None, None] + list(self.p1) + [self.x12] +
                       list(self.p4) + [self.x43])
        return list_fields

    def repr_fields(self):
        """
        Gets the fields in their simplified form

        Returns
        -------
        fields : list
            The fields that define the card

        """
        cp = set_blank_if_default(self.Cp(), 0)
        list_fields = (['CAERO3', self.eid, self.Pid(), cp, self.List_w(),
                        self.List_c1(), self.List_c2(), None, None] + list(self.p1) + [self.x12] +
                       list(self.p4) + [self.x43])
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class CAERO4(BaseCard):
    """
    Aerodynamic Macro-Strip Element Connection
    Defines an aerodynamic macro element for Strip theory.

    +--------+-----+-----+----+-------+--------+--------+--------+------+
    |   1    |  2  |  3  | 4  |   5   |   6    |    7   |   8    |   9  |
    +========+=====+=====+====+=======+========+========+========+======+
    | CAERO4 | EID | PID | CP | NSPAN | NCHORD |        |        |      |
    +--------+-----+-----+----+-------+--------+--------+--------+------+
    |        |  X1 | Y1  | Z1 |  X12  |   X4   |   Y4   |   Z4   | X43  |
    +--------+-----+-----+----+-------+--------+--------+--------+------+
    """
    type = 'CAERO4'
    _properties = ['shape', 'xy']

    @classmethod
    def _init_from_empty(cls):
        eid = 1
        pid = 1
        p1 = [0., 0., 0.]
        p4 = [0., 10., 0.]
        x12 = 10.
        x43 = 10.
        return CAERO4(eid, pid, p1, x12, p4, x43, cp=0, nspan=0, lspan=0, comment='')

    def __init__(self, eid, pid, p1, x12, p4, x43,
                 cp=0, nspan=0, lspan=0, comment=''):
        """
        Defines a CAERO4 card, which defines a strip theory surface.

        Parameters
        ----------
        eid : int
            element id
        pid : int, PAERO4
            int : PAERO4 ID
            PAERO4 : PAERO4 object (xref)
        p1 : (1, 3) ndarray float
            xyz location of point 1 (leading edge; inboard)
        p4 : (1, 3) ndarray float
            xyz location of point 4 (leading edge; outboard)
        x12 : float
            distance along the flow direction from node 1 to node 2
            (typically x, root chord)
        x43 : float
            distance along the flow direction from node 4 to node 3
            (typically x, tip chord)
        cp : int, CORDx; default=0
            int : coordinate system
            CORDx : Coordinate object (xref)
        nspan : int; default=0
            int > 0 : N spanwise boxes distributed evenly
            int = 0 : use lchord
        lspan : int, AEFACT; default=0
            int > 0 : AEFACT reference for non-uniform nspan
            int = 0 : use nspan
        comment : str; default=''
             a comment for the card

        """
        BaseCard.__init__(self)
        if comment:
            self.comment = comment

        #: Element identification number
        self.eid = eid

        #: Property identification number of a PAERO4 entry.
        self.pid = pid

        #: Coordinate system for locating point 1.
        self.cp = cp
        self.nspan = nspan
        self.lspan = lspan
        self.p1 = np.asarray(p1)
        self.x12 = x12
        self.p4 = np.asarray(p4)
        self.x43 = x43
        self.pid_ref = None
        self.cp_ref = None
        self.lspan_ref = None
        self.box_ids = None

    def validate(self):
        if self.nspan == 0 and self.lspan == 0:
            msg = 'NSPAN or LSPAN must be greater than 0; nspan=%r nlspan=%s\n' % (
                self.nspan, self.lspan)
            raise RuntimeError(msg)
        assert len(self.p1) == 3, 'p1=%s' % self.p1
        assert len(self.p4) == 3, 'p4=%s' % self.p4
        assert self.x12 > 0., 'x12=%s' % self.x12
        assert self.x43 >= 0., 'x43=%s' % self.x43

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a CAERO4 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        eid = integer(card, 1, 'eid')
        pid = integer(card, 2, 'pid')
        cp = integer_or_blank(card, 3, 'cp', 0)
        nspan = integer_or_blank(card, 4, 'nspan', 0)
        lspan = integer_or_blank(card, 5, 'lspan', 0)

        p1 = np.array([
            double_or_blank(card, 9, 'x1', 0.0),
            double_or_blank(card, 10, 'y1', 0.0),
            double_or_blank(card, 11, 'z1', 0.0)])
        x12 = double_or_blank(card, 12, 'x12', 0.)

        p4 = np.array([
            double_or_blank(card, 13, 'x4', 0.0),
            double_or_blank(card, 14, 'y4', 0.0),
            double_or_blank(card, 15, 'z4', 0.0)])
        x43 = double_or_blank(card, 16, 'x43', 0.)
        assert len(card) <= 17, 'len(CAERO4 card) = %i\ncard=%s' % (len(card), card)
        return CAERO4(eid, pid, p1, x12, p4, x43,
                      cp=cp, nspan=nspan, lspan=lspan, comment=comment)

    def get_points(self):
        p1 = self.cp_ref.transform_node_to_global(self.p1)
        p4 = self.cp_ref.transform_node_to_global(self.p4)
        p2 = p1 + np.array([self.x12, 0., 0.])
        p3 = p4 + np.array([self.x43, 0., 0.])
        return [p1, p2, p3, p4]

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        msg = ', which is required by CAERO4 eid=%s' % self.eid
        self.cp_ref = model.Coord(self.cp, msg=msg)
        self.pid_ref = model.PAero(self.pid, msg=msg)  # links to PAERO4 (not added)

        if self.nspan == 0:
            assert isinstance(self.lspan, integer_types), self.lspan
            self.lspan_ref = model.AEFact(self.lspan, msg)
        self._init_ids()

    def safe_cross_reference(self, model, xref_errors):
        msg = ', which is required by CAERO4 eid=%s' % self.eid
        self.pid_ref = model.safe_paero(self.pid, self.eid, xref_errors, msg=msg)  # links to PAERO4 (not added)
        self.cp_ref = model.safe_coord(self.cp, self.eid, xref_errors, msg=msg)

        if self.nspan == 0:
            assert isinstance(self.lspan, integer_types), self.lspan
            self.lspan_ref = model.safe_aefact(self.lspan, self.eid, xref_errors, msg)
        self._init_ids()

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.pid = self.Pid()
        self.cp = self.Cp()
        if self.nspan == 0:
            self.lspan = self.get_LSpan()
        self.pid_ref = None
        self.cp_ref = None
        self.lspan_ref = None

    def Cp(self):
        if self.cp_ref is not None:
            return self.cp_ref.cid
        return self.cp

    def Pid(self):
        if self.pid_ref is not None:
            return self.pid_ref.pid
        return self.pid

    def _init_ids(self, dtype='int32'):
        """
        Fill `self.box_ids` with the sub-box ids. Shape is (nchord, nspan)
        """
        nchord, nspan = self.shape
        assert nchord >= 1, 'nchord=%s' % nchord
        assert nspan >= 1, 'nspan=%s' % nspan
        self.box_ids = np.zeros((nchord, nspan), dtype=dtype)

        try:
            for ichord in range(nchord):
                for ispan in range(nspan):
                    self.box_ids[ichord, ispan] = self.eid + ichord + ispan * nchord
        except OverflowError:
            if dtype == 'int64':
                msg = 'eid=%s ichord=%s ispan=%s nchord=%s' % (
                    self.eid, ichord, ispan, nchord)
                raise OverflowError(msg)
            self._init_ids(dtype='int64')

    @property
    def shape(self):
        """returns (nelements_nchord, nelements_span)"""
        nchord = 1

        if self.nspan == 0:
            y = self.lspan_ref.fractions
            nspan = len(y) - 1
        else:
            nspan = self.nspan
        if nspan < 1:
            msg = 'CAERO4 eid=%s nspan=%s lspan=%s' % (
                self.eid, self.nspan, self.lspan)
            raise RuntimeError(msg)
        return nchord, nspan

    def get_npanel_points_elements(self):
        """
        Gets the number of sub-points and sub-elements for the CAERO card

        Returns
        -------
        npoints : int
            The number of nodes for the CAERO
        nelmements : int
            The number of elements for the CAERO

        """
        nchord, nspan = self.shape
        nelements = nchord * nspan
        npoints = (nchord + 1) * (nspan + 1)
        return npoints, nelements

    @property
    def xy(self):
        """
        Returns
        -------
        x : (nchord,) ndarray
            The percentage x location in the chord-wise direction of each panel
        y : (nspan,) ndarray
            The percentage y location in the span-wise direction of each panel

        """
        x = np.linspace(0., 1., num=2)  # nchord=1

        if self.nspan == 0:
            y = self.lspan_ref.fractions
            nspan = len(y) - 1
        else:
            nspan = self.nspan
            y = np.linspace(0., 1., nspan + 1)

        if nspan < 1:
            msg = 'CAERO4 eid=%s nspan=%s lspan=%s' % (
                self.eid, self.nspan, self.lspan)
            raise RuntimeError(msg)
        return x, y

    def panel_points_elements(self):
        """
        Gets the sub-points and sub-elements for the CAERO card

        Returns
        -------
        points : (nnodes,3) ndarray of floats
            the array of points
        elements : (nelements,4) ndarray of integers
            the array of point ids

        """
        p1, p2, p3, p4 = self.get_points()
        x, y = self.xy
        return points_elements_from_quad_points(p1, p2, p3, p4, x, y, dtype='int32')

    def get_LSpan(self):
        if isinstance(self.lspan, integer_types):
            return self.lspan
        return self.lspan_ref.sid

    def raw_fields(self):
        """
        Gets the fields in their unmodified form

        Returns
        -------
        fields : list
            The fields that define the card

        """
        list_fields = (['CAERO4', self.eid, self.Pid(), self.Cp(), self.nspan,
                        self.get_LSpan(), None, None, None,] + list(self.p1) + [self.x12] +
                       list(self.p4) + [self.x43])
        return list_fields

    def repr_fields(self):
        """
        Gets the fields in their simplified form

        Returns
        -------
        fields : list
            The fields that define the card

        """
        cp = set_blank_if_default(self.Cp(), 0)

        nspan = set_blank_if_default(self.nspan, 0)
        lspan = set_blank_if_default(self.get_LSpan(), 0)

        list_fields = (['CAERO4', self.eid, self.Pid(), cp, nspan,
                        lspan, None, None, None,] + list(self.p1) + [self.x12] +
                       list(self.p4) + [self.x43])
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class CAERO5(BaseCard):
    """
    Defines an aerodynamic macro element for Piston theory.

    +--------+------+------+-----+-------+-------+-------+--------+-------+
    |   1    |  2   |   3  |  4  |   5   |   6   |   7   |   8    |   9   |
    +========+======+======+=====+=======+=======+=======+========+=======+
    | CAERO5 | EID  | PID  | CP  | NSPAN | LSPAN | NTHRY | NTHICK |       |
    +--------+------+------+-----+-------+-------+-------+--------+-------+
    |        | X1   |  Y1  | Z1  |  X12  |  X4   |  Y4   |   Z4   |  X43  |
    +--------+------+------+-----+-------+-------+-------+--------+-------+
    | CAERO5 | 6000 | 6001 | 100 |       |  315  |   0   |   0    |       |
    +--------+------+------+-----+-------+-------+-------+--------+-------+
    |        | 0.0  |  0.0 | 0.0 |  1.0  |  0.2  |  1.0  |   0.   |  0.8  |
    +--------+------+------+-----+-------+-------+-------+--------+-------+
    """
    type = 'CAERO5'

    @classmethod
    def _init_from_empty(cls):
        eid = 1
        pid = 1
        p1 = [0., 0., 0.]
        p4 = [0., 10., 0.]
        x12 = 1.
        x43 = 0.5
        nspan = 5
        return CAERO5(eid, pid, p1, x12, p4, x43,
                      cp=0, nspan=nspan, lspan=0, ntheory=0, nthick=0, comment='')

    def __init__(self, eid, pid, p1, x12, p4, x43,
                 cp=0, nspan=0, lspan=0, ntheory=0, nthick=0,
                 comment=''):
        """
        Defines a CAERO5 card, which defines elements for Piston theory
        (high supersonic flow where the normal Mach is less than 1).

        Parameters
        ----------
        eid : int
            element id
        pid : int
            PAERO5 ID
        p1 : (1, 3) ndarray float
            xyz location of point 1 (leading edge; inboard)
        p4 : (1, 3) ndarray float
            xyz location of point 4 (leading edge; outboard)
        x12 : float
            distance along the flow direction from node 1 to node 2; (typically x, root chord)
        x43 : float
            distance along the flow direction from node 4 to node 3; (typically x, tip chord)
        cp : int, CORDx; default=0
            int : coordinate system
        nspan : int; default=0
            int > 0 : N spanwise boxes distributed evenly
            int = 0 : use lchord
        lspan : int, AEFACT; default=0
            int > 0 : AEFACT reference for non-uniform nspan
            int = 0 : use nspan
        ntheory : int; default=0
            ???
            valid_theory = {0, 1, 2}
        nthick : int; default=0
            ???
        comment : str; default=''
             a comment for the card

        """
        BaseCard.__init__(self)
        if comment:
            self.comment = comment

        #: Element identification number
        self.eid = eid

        #: Property identification number of a PAERO5 entry.
        self.pid = pid

        #: Coordinate system for locating point 1.
        self.cp = cp
        self.nspan = nspan
        self.lspan = lspan
        self.ntheory = ntheory
        self.nthick = nthick
        self.p1 = np.asarray(p1, dtype='float64')
        self.x12 = x12
        self.p4 = np.asarray(p4, dtype='float64')
        self.x43 = x43
        self.pid_ref = None
        self.cp_ref = None
        self.lspan_ref = None

        assert self.x12 > 0., 'x12=%s' % self.x12
        if not (self.nspan > 0 or self.lspan > 0):
            msg = 'nspan=%r or lspan=%r must be > 0' % (self.nspan, self.lspan)
            raise ValueError(msg)
        if not (self.x12 > 0.0 or self.x43 > 0.0):
            msg = 'x12=%r or x43=%r must be > 0.0' % (self.x12, self.x43)
            raise ValueError(msg)

    def validate(self):
        assert self.ntheory in [0, 1, 2], 'ntheory=%r' % self.ntheory

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a CAERO5 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        eid = integer(card, 1, 'eid')
        pid = integer(card, 2, 'pid')
        cp = integer_or_blank(card, 3, 'cp', 0)
        nspan = integer_or_blank(card, 4, 'nspan', 0)
        lspan = integer_or_blank(card, 5, 'lspan', 0)
        ntheory = integer_or_blank(card, 6, 'ntheory', 0)
        nthick = integer_or_blank(card, 7, 'nthick')
        # 8 - blank
        p1 = np.array([
            double_or_blank(card, 9, 'x1', 0.0),
            double_or_blank(card, 10, 'y1', 0.0),
            double_or_blank(card, 11, 'z1', 0.0)])
        x12 = double(card, 12, 'x12')
        p4 = np.array([
            double_or_blank(card, 13, 'x4', 0.0),
            double_or_blank(card, 14, 'y4', 0.0),
            double_or_blank(card, 15, 'z4', 0.0)])
        x43 = double_or_blank(card, 16, 'x43', 0.0)
        assert len(card) <= 17, 'len(CAERO3 card) = %i\ncard=%s' % (len(card), card)
        return CAERO5(eid, pid, p1, x12, p4, x43,
                      cp=cp, nspan=nspan, lspan=lspan, ntheory=ntheory, nthick=nthick,
                      comment=comment)

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        msg = ', which is required by CAERO5 eid=%s' % self.eid
        self.pid_ref = model.PAero(self.pid, msg=msg)
        self.cp_ref = model.Coord(self.cp, msg=msg)
        if self.nspan == 0:
            self.lspan_ref = model.AEFact(self.lspan, msg=msg)

    def safe_cross_reference(self, model, xref_errors):
        xref_errors = {}
        msg = ', which is required by CAERO5 eid=%s' % self.eid
        self.pid_ref = model.safe_paero(self.pid, self.eid, xref_errors, msg=msg)
        self.cp_ref = model.safe_coord(self.cp, self.eid, xref_errors, msg=msg)
        if self.nspan == 0:
            self.lspan_ref = model.safe_aefact(self.lspan, self.eid, xref_errors, msg=msg)

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.pid = self.Pid()
        self.cp = self.Cp()
        if self.nspan == 0:
            self.lspan = self.LSpan()
        self.pid_ref = None
        self.cp_ref = None
        self.lspan_ref = None

    def get_points(self):
        p1 = self.cp_ref.transform_node_to_global(self.p1)
        p4 = self.cp_ref.transform_node_to_global(self.p4)
        p2 = p1 + np.array([self.x12, 0., 0.])
        p3 = p4 + np.array([self.x43, 0., 0.])
        return [p1, p2, p3, p4]

    def get_npanel_points_elements(self):
        msg = 'CAERO5 eid=%s nspan=%s lspan=%s' % (
            self.eid, self.nspan, self.lspan)
        if self.nspan == 0:
            y = self.lspan_ref.fractions
            nspan = len(y) - 1
        else:
            nspan = self.nspan
        assert nspan >= 1, msg

        nchord = 1
        nelements = nchord * nspan
        npoints = (nchord + 1) * (nspan + 1)
        return npoints, nelements

    def panel_points_elements(self):
        p1, p2, p3, p4 = self.get_points()

        msg = 'CAERO5 eid=%s nspan=%s lspan=%s' % (
            self.eid, self.nspan, self.lspan)
        if self.nspan == 0:
            y = self.lspan_ref.fractions
            nspan = len(y) - 1
        else:
            nspan = self.nspan
            y = np.linspace(0., 1., nspan + 1)
        assert nspan >= 1, msg

        x = np.array([0., 1.], dtype='float64')
        assert nspan >= 1, msg

        return points_elements_from_quad_points(p1, p2, p3, p4, x, y, dtype='int32')

    def c1_c2(self, mach):
        p1, unused_p2, unused_p3, p4 = self.get_points()
        #i = p2 - p1
        #ihat = i / norm(i)
        #k = cross(ihat, p4-p1)
        #khat = k / norm(k)
        #jhat = cross(khat, ihat)
        #b = self.p4 - self.p1
        L = np.linalg.norm(p4 - p1)

        # b = L * cos(Lambda)
        # ci = L * sin(Lambda)

        if self.ntheory == 0:
            # piston theory
            pass
        elif self.ntheory == 1:
            raise NotImplementedError('ntheory=%s' % self.ntheory)
            #gamma = 1.4
            #lambda_sweep = 0.
            #c1 = 1.
            #secL = 1 / np.cos(lambda_sweep)
            #secL2 = secL ** 2
            #ma2_secL2 = mach ** 2 - secL2
            #c1 = mach / ma2_secL2 ** 0.5
            #c2 = (mach ** 4 * (gamma + 1) - 4 * secL2 * ma2_secL2) / (4 * ma2_secL2 ** 2)
        else:
            gamma = 1.4

            # the advance in x
            ci = p4[0] - p1[0]

            # sweep angle
            lambda_sweep = np.arcsin(ci / L)

            sec_lambda = 1 / np.cos(lambda_sweep)
            sec_lambda2 = sec_lambda ** 2
            ma2_sec_lambda2 = mach ** 2 - sec_lambda2
            c1 = mach / ma2_sec_lambda2 ** 0.5
            c2 = (
                (mach ** 4 * (gamma + 1) - 4 * sec_lambda2 * ma2_sec_lambda2) /
                (4 * ma2_sec_lambda2 ** 2)
            )
        return c1, c2

    def Cp(self):
        if self.cp_ref is not None:
            return self.cp_ref.cid
        return self.cp

    def Pid(self):
        if self.pid_ref is not None:
            return self.pid_ref.pid
        return self.pid

    def LSpan(self):
        if self.lspan_ref is not None:
            return self.lspan_ref.sid
        return self.lspan

    def repr_fields(self):
        """
        Gets the fields in their simplified form

        Returns
        -------
        fields : list
            The fields that define the card
        """
        nspan = self.nspan
        lspan = self.LSpan()
        ntheory = self.ntheory
        cp = set_blank_if_default(self.Cp(), 0)
        list_fields = (['CAERO5', self.eid, self.Pid(), cp, nspan, lspan,
                        ntheory, self.nthick, None,] + list(self.p1) + [self.x12] +
                       list(self.p4) + [self.x43])
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class PAERO5(BaseCard):
    type = 'PAERO5'
    _properties = ['ltaus_id', 'lxis_id']

    @classmethod
    def _init_from_empty(cls):
        pid = 1
        caoci = [0., 0., 0.]
        return PAERO5(pid, caoci, nalpha=0, lalpha=0, nxis=0, lxis=0, ntaus=0, ltaus=0, comment='')

    def __init__(self, pid, caoci,
                 nalpha=0, lalpha=0,
                 nxis=0, lxis=0,
                 ntaus=0, ltaus=0, comment=''):
        """
        +--------+-------+--------+--------+---------+-------+-------+-------+
        |   1    |   2   |    3   |   4    |    5    |   6   |   7   |   8   |
        +========+=======+========+========+=========+=======+=======+=======+
        | PAERO5 |  PID  | NALPHA | LALPHA |  NXIS   | LXIS  | NTAUS | LTAUS |
        +--------+-------+--------+--------+---------+-------+-------+-------+
        |        | CAOC1 | CAOC2  | CAOC3  |  CAOC4  | CAOC5 |       |       |
        +--------+-------+--------+--------+---------+-------+-------+-------+
        | PAERO5 | 7001  |   1    |  702   |    1    | 701   |   1   |  700  |
        +--------+-------+--------+--------+---------+-------+-------+-------+
        |        |  0.0  |  0.0   |  5.25  | 3.99375 |  0.0  |       |       |
        +--------+-------+--------+--------+---------+-------+-------+-------+
        """
        BaseCard.__init__(self)
        if comment:
            self.comment = comment
        self.pid = pid
        self.nalpha = nalpha
        self.lalpha = lalpha

        # number of dimensionless chord coordinates in zeta ()
        self.nxis = nxis
        # ID of AEFACT that lists zeta
        self.lxis = lxis

        # number of dimensionless thickess coordinates in tau
        self.ntaus = ntaus
        # ID of AEFACT that lists thickness ratios (t/c)
        self.ltaus = ltaus

        # ca/c - control surface chord / strip chord
        self.caoci = np.array(caoci, dtype='float64')
        self.lxis_ref = None
        self.ltaus_ref = None

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a PAERO5 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        pid = integer(card, 1, 'property_id')
        nalpha = integer_or_blank(card, 2, 'nalpha', default=0)
        lalpha = integer_or_blank(card, 3, 'lalpha', default=0)

        nxis = integer_or_blank(card, 4, 'nxis', default=0)
        lxis = integer_or_blank(card, 5, 'lxis', default=0)

        ntaus = integer_or_blank(card, 6, 'ntaus', default=0)
        ltaus = integer_or_blank(card, 7, 'ltaus', default=0)

        caoci = []
        for n, i in enumerate(range(9, len(card))):
            ca = double(card, i, 'ca/ci_%i' % (n+1))
            caoci.append(ca)
        return PAERO5(pid, caoci,
                      nalpha=nalpha, lalpha=lalpha, nxis=nxis, lxis=lxis,
                      ntaus=ntaus, ltaus=ltaus,
                      comment=comment)
    @property
    def lxis_id(self):
        if self.lxis_ref is not None:
            return  self.lxis_ref.sid
        return self.lxis

    @property
    def ltaus_id(self):
        if self.ltaus_ref is not None:
            return self.ltaus_ref.sid
        return self.ltaus

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        msg = ', which is required by PAERO5 eid=%s' % self.pid
        if self.lxis != 0:
            self.lxis_ref = model.AEFact(self.lxis_id, msg=msg)
        if self.ltaus != 0:
            self.ltaus_ref = model.AEFact(self.ltaus_id, msg=msg)

    def safe_cross_reference(self, model, xref_errors):
        msg = ', which is required by PAERO5 eid=%s' % self.pid
        if self.lxis != 0:
            self.lxis_ref = model.safe_aefact(self.lxis_id, self.pid, xref_errors, msg=msg)
        if self.ltaus != 0:
            self.ltaus_ref = model.safe_aefact(self.ltaus_id, self.pid, xref_errors, msg=msg)

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.lxis = self.lxis_id
        self.ltaus = self.ltaus_id
        self.lxis_ref = None
        self.ltaus_ref = None

    def raw_fields(self):
        list_fields = ['PAERO5', self.pid, self.nalpha, self.lalpha, self.nxis,
                       self.lxis_id, self.ntaus, self.ltaus_id] + list(self.caoci)
        return list_fields

    def repr_fields(self):
        list_fields = self.raw_fields()
        list_fields.insert(8, None)
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        return self.comment + print_card_8(card)

    #def integrals(self):
        ## chord location
        #x = self.lxis.fractions

        ## thickness
        #y = self.ltaus.fractions

        ## slope of airfoil semi-thickness
        #yp = derivative1(y/2, x)

        ## x hinge
        #for xh in self.caoci:
            #I1 = integrate(yp, x, 0., 1.)
            #I2 = integrate(x * yp, x, 0., 1.)
            #I3 = integrate(x**2*yp, x, 0., 1.)
            #I4 = integrate(yp**2, x, 0., 1.)
            #I5 = integrate(x**2 * yp**2, x, 0., 1.)

            #J1 = integrate(yp, x, xh, 1.)
            #J2 = integrate(x * yp, x, xh, 1.)
            #J3 = integrate(x**2*yp, x, xh, 1.)
            #J4 = integrate(yp**2, x, xh, 1.)
            #J5 = integrate(x**2 * yp**2, x, xh, 1.)

        #return(I1, I2, I3, I4, I5,
               #J1, J2, J3, J4, J5)


class MONPNT1(BaseCard):
    """
    +---------+---------+------+-----+-----+-------+------+----+----+
    |    1    |    2    |  3   |  4  |  5  |   6   |   7  | 8  | 9  |
    +=========+=========+======+=====+=====+=======+======+====+====+
    | MONPNT1 |  NAME   |                   LABEL                   |
    +---------+---------+------+-----+-----+-------+------+----+----+
    |         |  AXES   | COMP | CP  |  X  |   Y   |   Z  | CD |    |
    +---------+---------+------+-----+-----+-------+------+----+----+
    | MONPNT1 | WING155 |    Wing Integrated Load to Butline 155    |
    +---------+---------+------+-----+-----+-------+------+----+----+
    |         |    34   | WING |     | 0.0 | 155.0 | 15.0 |    |    |
    +---------+---------+------+-----+-----+-------+------+----+----+
    """
    type = 'MONPNT1'

    @classmethod
    def _init_from_empty(cls):
        name = 'WING'
        label = 'Wing Integrated Load to Butline'
        axes = '6'
        aecomp_name = 'FLAP'
        xyz = [0., 1., 2.]
        return MONPNT1(name, label, axes, aecomp_name, xyz, cp=0, cd=None, comment='')

    def __init__(self, name, label, axes, aecomp_name, xyz, cp=0, cd=None, comment=''):
        """
        Creates a MONPNT1 card

        Parameters
        ----------
        name : str
            Character string of up to 8 characters identifying the
            monitor point
        label : str
            A string comprising no more than 56 characters
            that identifies and labels the monitor point.
        axes : str
            components {1,2,3,4,5,6}
        aecomp_name : str
            name of the AECOMP/AECOMPL entry
        xyz : List[float, float, float]; default=None
            The coordinates in the CP coordinate system about which the
            loads are to be monitored.
            None : [0., 0., 0.]
        cp : int, CORDx; default=0
           coordinate system of XYZ
        cd : int; default=None -> cp
            the coordinate system for load outputs
        comment : str; default=''
            a comment for the card

        Notes
        -----
        CD - MSC specific field

        """
        BaseCard.__init__(self)
        if comment:
            self.comment = comment
        if cd is None:
            cd = cp
        xyz = np.asarray(xyz)

        self.name = name
        self.label = label
        self.axes = axes
        self.comp = aecomp_name
        self.cp = cp
        self.xyz = xyz
        self.cd = cd
        assert len(xyz) == 3, xyz
        self.cp_ref = None
        self.cd_ref = None

    @classmethod
    def add_card(cls, card, comment=''):
        name = string(card, 1, 'name')

        label_fields = [labeli for labeli in card[2:8] if labeli is not None]
        label = ''.join(label_fields).strip()
        assert len(label) <= 56, label

        axes = parse_components(card, 9, 'axes')
        comp = string(card, 10, 'comp')
        cp = integer_or_blank(card, 11, 'cp', 0)
        xyz = [
            double_or_blank(card, 12, 'x', default=0.0),
            double_or_blank(card, 13, 'y', default=0.0),
            double_or_blank(card, 14, 'z', default=0.0),
        ]
        cd = integer_or_blank(card, 15, 'cd', cp)
        return MONPNT1(name, label, axes, comp, xyz, cp=cp, cd=cd, comment=comment)

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        msg = ', which is required by MONPNT1 name=%s' % self.name
        self.cp_ref = model.Coord(self.cp, msg=msg)
        self.cd_ref = model.Coord(self.cd, msg=msg)

    def safe_cross_reference(self, model, xref_errors):
        msg = ', which is required by MONPNT1 name=%s' % self.name
        self.cp_ref = model.safe_coord(self.cp, self.name, xref_errors, msg=msg)
        self.cd_ref = model.safe_coord(self.cd, self.name, xref_errors, msg=msg)

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.cp = self.Cp()
        self.cd = self.Cd()
        self.cp_ref = None
        self.cd_ref = None

    def Cp(self):
        if self.cp_ref is not None:
            return self.cp_ref.cid
        return self.cp

    def Cd(self):
        if self.cd_ref is not None:
            return self.cd_ref.cid
        return self.cd

    def raw_fields(self):
        list_fields = [
            'MONPNT1', self.name, self.label.strip(), self.axes, self.comp,
            self.Cp(),] + list(self.xyz) + [self.Cd()]
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        cp = self.Cp()
        x, y, z = self.xyz
        cd = self.Cd()

        # Default = the coordinate system specified by the CP field
        if cd == cp:
            cd = ''

        msg = 'MONPNT1 %-8s%s\n' % (self.name, self.label)
        msg += '        %-8s%-8s%-8s%-8s%-8s%-8s%-8s\n' % (
            self.axes, self.comp, cp,
            print_float_8(x), print_float_8(y), print_float_8(z),
            cd)
        #card = self.repr_fields()
        return self.comment + msg

    def __repr__(self):
        return self.write_card()


class MONPNT2(BaseCard):
    """MSC Nastran specific card"""
    type = 'MONPNT2'

    @classmethod
    def _init_from_empty(cls):
        name = 'WING'
        label = 'Wing Integrated Load to Butline'
        table = 'MYTABLE'
        Type = 'CAT'
        nddl_item = 42
        eid = 2
        return MONPNT2(name, label, table, Type, nddl_item, eid, comment='')

    def __init__(self, name, label, table, Type, nddl_item, eid, comment=''):
        BaseCard.__init__(self)
        if comment:
            self.comment = comment
        self.name = name
        self.label = label
        self.table = table
        self.Type = Type
        self.nddl_item = nddl_item
        self.eid = eid

    def validate(self):
        assert self.table in ['STRESS', 'FORCE', 'STRAIN'], self.table

    @classmethod
    def add_card(cls, card, comment=''):
        name = string(card, 1, 'name')

        label_fields = [labeli for labeli in card[2:8] if labeli is not None]
        label = ''.join(label_fields).strip()
        assert len(label) <= 56, label

        table = string(card, 9, 'table')
        Type = string(card, 10, 'type')
        nddl_item = integer_or_string(card, 11, 'nddl_item')
        #nddl_item = integer_or_blank(card, 11, 'nddl_item')
        eid = integer_or_blank(card, 12, 'eid')
        return MONPNT2(name, label, table, Type, nddl_item, eid, comment=comment)

    def cross_reference(self, model: BDF) -> None:
        pass

    def safe_cross_reference(self, model, unused_xref_errors):
        self.cross_reference(model)

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        pass

    def raw_fields(self):
        list_fields = [
            'MONPNT2', self.name, self.label.strip(),
            self.table, self.Type, self.nddl_item, self.eid]
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        msg = 'MONPNT2 %-8s%s\n' % (self.name, self.label)
        msg += ('        %-8s%-8s%-8s%-8s\n' % (
            self.table, self.Type, self.nddl_item, self.eid
        ))
        #card = self.repr_fields()
        return self.comment + msg.rstrip() + '\n'

    def __repr__(self):
        return self.write_card()


class MONPNT3(BaseCard):
    """MSC Nastran specific card"""
    type = 'MONPNT3'

    @classmethod
    def _init_from_empty(cls):
        name = 'WING'
        label = 'Wing Integrated Load to Butline'
        axes = '6'
        grid_set = 10
        elem_set = 11
        xyz = [0., 1., 2.]
        return MONPNT3(name, label, axes, grid_set, elem_set, xyz,
                       cp=0, cd=None, xflag=None, comment='')

    def __init__(self, name, label, axes, grid_set, elem_set, xyz,
                 cp=0, cd=None, xflag=None, comment=''):
        BaseCard.__init__(self)
        if comment:
            self.comment = comment
        if cd is None:
            cd = cp
        xyz = np.asarray(xyz)

        self.name = name
        self.label = label
        self.axes = axes
        #self.comp = comp
        self.grid_set = grid_set
        self.elem_set = elem_set
        self.xyz = xyz
        self.xflag = xflag
        self.cp = cp
        self.cd = cd
        self.cp_ref = None
        self.cd_ref = None

    @classmethod
    def add_card(cls, card, comment=''):
        name = string(card, 1, 'name')

        label_fields = [labeli for labeli in card[2:8] if labeli is not None]
        label = ''.join(label_fields).strip()
        assert len(label) <= 56, label

        axes = parse_components(card, 9, 'axes')
        grid_set = integer(card, 10, 'grid_set')
        elem_set = integer_or_blank(card, 11, 'elem_set')
        cp = integer_or_blank(card, 12, 'cp', 0)
        xyz = [
            double_or_blank(card, 13, 'x', default=0.0),
            double_or_blank(card, 14, 'y', default=0.0),
            double_or_blank(card, 15, 'z', default=0.0),
        ]
        xflag = string_or_blank(card, 16, 'xflag')
        cd = integer_or_blank(card, 17, 'cd', cp)
        return MONPNT3(name, label, axes, grid_set, elem_set, xyz,
                       cp=cp, cd=cd, xflag=xflag, comment=comment)

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        msg = ', which is required by MONPNT3 name=%s' % self.name
        self.cp_ref = model.Coord(self.cp, msg=msg)
        self.cd_ref = model.Coord(self.cd, msg=msg)

    def safe_cross_reference(self, model, xref_errors):
        msg = ', which is required by MONPNT3 name=%s' % self.name
        self.cp_ref = model.safe_coord(self.cp, self.name, xref_errors, msg=msg)
        self.cd_ref = model.safe_coord(self.cd, self.name, xref_errors, msg=msg)

    def uncross_reference(self) -> None:
        self.cp = self.Cp()
        self.cd = self.Cd()
        self.cp_ref = None
        self.cd_ref = None

    def Cp(self):
        if self.cp_ref is not None:
            return self.cp_ref.cid
        return self.cp

    def Cd(self):
        if self.cd_ref is not None:
            return self.cd_ref.cid
        return self.cd

    def raw_fields(self):
        list_fields = [
            'MONPNT3', self.name, self.label.strip(),
            self.axes, self.grid_set, self.elem_set, self.Cp()
            ] + list(self.xyz) + [self.xflag, self.Cd()]
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        cp = self.Cp()
        cd = self.Cd()
        if cp == cd:
            cd = ''
        xflag = self.xflag
        if xflag is None:
            xflag = ''

        x, y, z = self.xyz
        msg = 'MONPNT3 %-8s%s\n' % (self.name, self.label)
        msg += ('        %-8s%-8s%-8s%-8s%-8s%-8s%-8s%-8s\n'
                '         %-8s' % (
                    self.axes, self.grid_set, self.elem_set, cp,
                    print_float_8(x), print_float_8(y), print_float_8(z),
                    xflag, cd
                ))
        #card = self.repr_fields()
        return self.comment + msg.rstrip() + '\n'

    def __repr__(self):
        return self.write_card()

class MONDSP1(BaseCard):
    """
    +---------+---------+------+-----+-----+-------+------+----+----+
    |    1    |    2    |  3   |  4  |  5  |   6   |   7  | 8  | 9  |
    +=========+=========+======+=====+=====+=======+======+====+====+
    | MONPNT1 |  NAME   |                   LABEL                   |
    +---------+---------+------+-----+-----+-------+------+----+----+
    |         |  AXES   | COMP | CP  |  X  |   Y   |   Z  | CD |    |
    +---------+---------+------+-----+-----+-------+------+----+----+
    | MONPNT1 | WING155 |    Wing Integrated Load to Butline 155    |
    +---------+---------+------+-----+-----+-------+------+----+----+
    |         |    34   | WING |     | 0.0 | 155.0 | 15.0 |    |    |
    +---------+---------+------+-----+-----+-------+------+----+----+
    """
    type = 'MONDSP1'

    @classmethod
    def _init_from_empty(cls):
        name = 'WING'
        label = 'Wing Integrated Load to Butline'
        axes = '6'
        aecomp_name = 'FLAP'
        xyz = [0., 1., 2.]
        return MONDSP1(name, label, axes, aecomp_name, xyz, cp=0, cd=None, ind_dof='123', comment='')

    def __init__(self, name, label, axes, aecomp_name, xyz, cp=0, cd=None, ind_dof='123', comment=''):
        """
        Creates a MONDSP1 card

        Parameters
        ----------
        name : str
            Character string of up to 8 characters identifying the
            monitor point
        label : str
            A string comprising no more than 56 characters
            that identifies and labels the monitor point.
        axes : str
            components {1,2,3,4,5,6}
        aecomp_name : str
            name of the AECOMP/AECOMPL entry
        xyz : List[float, float, float]; default=None
            The coordinates in the CP coordinate system about which the
            loads are to be monitored.
            None : [0., 0., 0.]
        cp : int, CORDx; default=0
           coordinate system of XYZ
        cd : int; default=None -> cp
            the coordinate system for load outputs
        ind_dof : str; default='123'
            the dofs to map
        comment : str; default=''
            a comment for the card

        Notes
        -----
        MSC specific card

        """
        BaseCard.__init__(self)
        if comment:
            self.comment = comment
        if cd is None:
            cd = cp
        xyz = np.asarray(xyz)

        self.name = name
        self.label = label
        self.axes = axes
        self.comp = aecomp_name
        self.cp = cp
        self.xyz = xyz
        self.cd = cd
        self.ind_dof = ind_dof
        assert len(xyz) == 3, xyz
        self.cp_ref = None
        self.cd_ref = None

    @classmethod
    def add_card(cls, card, comment=''):
        name = string(card, 1, 'name')

        label_fields = [labeli for labeli in card[2:8] if labeli is not None]
        label = ''.join(label_fields).strip()
        assert len(label) <= 56, label

        axes = parse_components(card, 9, 'axes')
        comp = string(card, 10, 'comp')
        cp = integer_or_blank(card, 11, 'cp', 0)
        xyz = [
            double_or_blank(card, 12, 'x', default=0.0),
            double_or_blank(card, 13, 'y', default=0.0),
            double_or_blank(card, 14, 'z', default=0.0),
        ]
        cd = integer_or_blank(card, 15, 'cd', cp)
        ind_dof = components_or_blank(card, 16, 'ind_dof', '123')
        return MONDSP1(name, label, axes, comp, xyz, cp=cp, cd=cd, ind_dof=ind_dof, comment=comment)

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        msg = ', which is required by MONDSP1 name=%s' % self.name
        self.cp_ref = model.Coord(self.cp, msg=msg)
        self.cd_ref = model.Coord(self.cd, msg=msg)

    def safe_cross_reference(self, model, xref_errors):
        msg = ', which is required by MONDSP1 name=%s' % self.name
        self.cp_ref = model.safe_coord(self.cp, self.name, xref_errors, msg=msg)
        self.cd_ref = model.safe_coord(self.cd, self.name, xref_errors, msg=msg)

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.cp = self.Cp()
        self.cd = self.Cd()
        self.cp_ref = None
        self.cd_ref = None

    def Cp(self):
        if self.cp_ref is not None:
            return self.cp_ref.cid
        return self.cp

    def Cd(self):
        if self.cd_ref is not None:
            return self.cd_ref.cid
        return self.cd

    def raw_fields(self):
        list_fields = [
            'MONDSP1', self.name, self.label.strip(), self.axes, self.comp,
            self.Cp(),] + list(self.xyz) + [self.Cd(), self.ind_dof]
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        cp = self.Cp()
        x, y, z = self.xyz
        cd = self.Cd()

        # Default = the coordinate system specified by the CP field
        if cd == cp:
            cd = ''
        msg = 'MONDSP1 %-8s%s\n' % (self.name, self.label)
        msg += '        %-8s%-8s%-8s%-8s%-8s%-8s%-8s%-8s\n' % (
            self.axes, self.comp, cp,
            print_float_8(x), print_float_8(y), print_float_8(z),
            cd, self.ind_dof)
        #card = self.repr_fields()
        return self.comment + msg

    def __repr__(self):
        return self.write_card()


class PAERO1(BaseCard):
    """
    Defines associated bodies for the panels in the Doublet-Lattice method.

    +--------+-----+----+----+----+----+----+----+
    |    1   |  2  |  3 |  4 |  5 |  6 |  7 |  8 |
    +========+=====+====+====+====+====+====+====+
    | PAERO1 | PID | B1 | B2 | B3 | B4 | B5 | B6 |
    +--------+-----+----+----+----+----+----+----+

    """
    type = 'PAERO1'
    _field_map = {1: 'pid'}
    _properties = ['_field_map']

    def _get_field_helper(self, n):
        """
        Gets complicated parameters on the PAERO1 card

        Parameters
        ----------
        n : int
            the field number to update

        Returns
        -------
        value : int
            the value for the appropriate field

        """
        return self.Bi[n - 1]

    def _update_field_helper(self, n, value):
        """
        Updates complicated parameters on the PAERO1 card

        Parameters
        ----------
        n : int
            the field number to update
        value : varies
            the value for the appropriate field

        """
        self.caero_body_ids[n - 1] = value

    def object_attributes(self, mode='public', keys_to_skip=None, filter_properties=False):
        """
        List the names of attributes of a class as strings. Returns public
        attributes as default.

        Parameters
        ----------
        mode : str
            defines what kind of attributes will be listed
            * 'public' - names that do not begin with underscore
            * 'private' - names that begin with single underscore
            * 'both' - private and public
            * 'all' - all attributes that are defined for the object
        keys_to_skip : List[str]; default=None -> []
            names to not consider to avoid deprecation warnings

        Returns
        -------
        attribute_names : List[str]
            sorted list of the names of attributes of a given type or None
            if the mode is wrong

        """
        if keys_to_skip is None:
            keys_to_skip = []
        my_keys_to_skip = ['Bi']
        return object_attributes(self, mode=mode, keys_to_skip=keys_to_skip+my_keys_to_skip,
                                 filter_properties=filter_properties)

    @property
    def Bi(self):
        self.deprecated('Bi', 'caero_body_ids', '1.2')
        return self.caero_body_ids

    @Bi.setter
    def Bi(self, Bi):
        self.deprecated('Bi', 'caero_body_ids', '1.2')
        self.caero_body_ids = Bi

    @classmethod
    def _init_from_empty(cls):
        pid = 1
        return PAERO1(pid, caero_body_ids=None, comment='')

    def __init__(self, pid, caero_body_ids=None, comment=''):
        """
        Creates a PAERO1 card, which defines associated bodies for the
        panels in the Doublet-Lattice method.

        Parameters
        ----------
        pid : int
            PAERO1 id
        caero_body_ids : List[int]; default=None
            CAERO2 ids that are within the same IGROUP group
        comment : str; default=''
            a comment for the card

        """
        BaseCard.__init__(self)
        if comment:
            self.comment = comment
        self.pid = pid
        if caero_body_ids is None:
            caero_body_ids = []
        self.caero_body_ids = caero_body_ids

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a PAERO1 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        pid = integer(card, 1, 'pid')
        caero_body_ids = [interpret_value(field, card) for field in card[2:]]
        caero_body_ids2 = []

        for caero_body_id in caero_body_ids:
            if isinstance(caero_body_id, integer_types) and caero_body_id >= 0:
                caero_body_ids2.append(caero_body_id)
            elif caero_body_id is not None:
                msg = 'invalid caero_body_id value on PAERO1; caero_body_id=%r' % (caero_body_id)
                raise RuntimeError(msg)
            #else:
                #pass
        return PAERO1(pid, caero_body_ids, comment=comment)

    def cross_reference(self, model: BDF) -> None:
        pass

    def safe_cross_reference(self, model, xref_errors):
        pass

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        pass

    def Bodies(self):
        self.deprecated('Bodies', 'N/A', '1.2')
        return self.Bi

    def raw_fields(self):
        """
        Gets the fields in their unmodified form

        Returns
        -------
        fields : list[varies]
            the fields that define the card

        """
        list_fields = ['PAERO1', self.pid] + self.caero_body_ids
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.raw_fields()
        return self.comment + print_card_8(card)


class PAERO2(BaseCard):
    """
    Defines the cross-sectional properties of aerodynamic bodies.

    +--------+------+--------+-------+------+------+------+------+------+
    |    1   |  2   |   3    |    4  |   5  |   6  |   7  |   8  |   9  |
    +========+======+========+=======+======+======+======+======+======+
    | PAERO2 | PID  | ORIENT | WIDTH |  AR  | LRSB | LRIB | LTH1 | LTH2 |
    +--------+------+--------+-------+------+------+------+------+------+
    | THI1   | THN1 |  THI2  |  THN2 | THI3 | THN3 |      |      |      |
    +--------+------+--------+-------+------+------+------+------+------+

    """
    type = 'PAERO2'
    _field_map = {
        1: 'pid', 2:'orient', 3:'width', 4:'AR', 5:'lrsb', 6:'lrib',
        #7: 'lth1', 8:'lth2',
    }
    _properties = ['_field_map', 'lth1', 'lth2', ]

    def _get_field_helper(self, n):
        """
        Gets complicated parameters on the PAERO2 card

        Parameters
        ----------
        n : int
            the field number to update

        Returns
        -------
        value : varies
            the value for the appropriate field

        """
        nnew = n - 8
        spot = nnew // 2
        i = nnew % 2
        if i == 0:
            value = self.thi[spot]
        else:
            value = self.thn[spot]
        return value

    def _update_field_helper(self, n, value):
        """
        Updates complicated parameters on the PAERO2 card

        Parameters
        ----------
        n : int
            the field number to update
        value : varies
            the value for the appropriate field

        """
        nnew = n - 8
        spot = nnew // 2
        i = nnew % 2
        if i == 0:
            self.thi[spot] = value
        else:
            self.thn[spot] = value

    @classmethod
    def _init_from_empty(cls):
        pid = 1
        width = 10.
        AR = 1.
        thi = [None]
        thn = [None]
        orient = 'ZY'
        return PAERO2(pid, orient, width, AR, thi, thn, lrsb=None, lrib=None, lth=None, comment='')

    #def _finalize_hdf5(self, encoding):
        #"""hdf5 helper function"""
        #pass
        #print(self.get_stats())

    def __init__(self, pid, orient, width, AR,
                 thi, thn, lrsb=None, lrib=None, lth=None,
                 comment=''):
        """
        Creates a PAERO2 card, which defines additional cross-sectional
        properties for the CAERO2 geometry.

        Parameters
        ----------
        pid : int
            PAERO1 id
        orient : str
            Orientation flag. Type of motion allowed for bodies. Refers
            to the aerodynamic coordinate system of ACSID. See AERO entry.
            valid_orientations = {Z, Y, ZY}
        width : float
            Reference half-width of body and the width of the constant
            width interference tube
        AR : float
            Aspect ratio of the interference tube (height/width)
        thi / thn : List[int]
            The first (thi) and last (thn) interference element of a body
            to use the theta1/theta2 array
        lrsb : int; default=None
            int : AEFACT id containing a list of slender body half-widths
                  at the end points of the slender body elements
            None : use width
        lrib : int; default=None
            int : AEFACT id containing a list of interference body
                  half-widths at the end points of the interference elements
            None : use width
        lth : List[int, int]; default=None
            AEFACT id for defining theta arrays for interference calculations
            for theta1/theta2
        comment : str; default=''
            a comment for the card

        """
        BaseCard.__init__(self)
        if comment:
            self.comment = comment

        #: Property identification number. (Integer > 0)
        self.pid = pid

        #: Orientation flag. Type of motion allowed for bodies. Refers to
        #: the aerodynamic coordinate system of ACSID. See AERO entry.
        #: (Character = 'Z', 'Y', or 'ZY')
        self.orient = orient

        #: Reference half-width of body and the width of the constant width
        #: interference tube. (Real > 0.0)
        self.width = width

        #: Aspect ratio of the interference tube (height/width). float>0.
        self.AR = AR

        #: Identification number of an AEFACT entry containing a list of
        #: slender body half-widths at the end points of the
        #: slender body elements. If blank, the value of WIDTH will be used.
        #: (Integer > 0 or blank)
        self.lrsb = lrsb

        #: Identification number of an AEFACT entry containing a list of
        #: slender body half-widths at the end points of the
        #: interference elements. If blank, the value of WIDTH will be used.
        #: (Integer > 0 or blank)
        self.lrib = lrib

        #: Identification number of AEFACT entries for defining theta arrays for
        #: interference calculations. (Integer >= 0)
        if lth is None:
            lth = [None, None]
        self.lth = lth
        self.thi = thi
        self.thn = thn
        if self.lrsb == 0:
            self.lrsb = None
        if self.lrib == 0:
            self.lrib = None
        self.lrsb_ref = None
        self.lrib_ref = None

    @property
    def lth1(self):
        return self.lth[0]
    @property
    def lth2(self):
        return self.lth[1]

    @lth1.setter
    def lth1(self, lth1):
        self.lth[0] = lth1

    @lth2.setter
    def lth2(self, lth2):
        self.lth[1] = lth2

    def validate(self):
        assert self.orient in ['Z', 'Y', 'ZY'], 'PAERO2: orient=%r' % self.orient
        assert isinstance(self.AR, float), 'PAERO2: AR=%r type=%s' % (self.AR, type(self.AR))
        assert isinstance(self.width, float), 'PAERO2: width=%r type=%s' % (self.width, type(self.width))
        assert isinstance(self.thi, list), 'PAERO2: thi=%s type=%s' % (self.thi, type(self.thi))
        assert isinstance(self.thn, list), 'PAERO2: thn=%s type=%s' % (self.thn, type(self.thn))

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a PAERO2 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        pid = integer(card, 1, 'pid')
        orient = string(card, 2, 'orient')
        width = double(card, 3, 'width')
        AR = double(card, 4, 'AR')
        lrsb = integer_or_blank(card, 5, 'lrsb')
        lrib = integer_or_blank(card, 6, 'lrib')
        lth1 = integer_or_blank(card, 7, 'lth1')
        lth2 = integer_or_blank(card, 8, 'lth2')
        thi = []
        thn = []
        list_fields = [interpret_value(field, card) for field in card[9:]]
        nfields = len(list_fields)
        lth = [lth1, lth2]
        for i in range(9, 9 + nfields, 2):
            thi.append(integer(card, i, 'lth'))
            thn.append(integer(card, i + 1, 'thn'))
        return PAERO2(pid, orient, width, AR, thi, thn,
                      lrsb=lrsb, lrib=lrib, lth=lth,
                      comment=comment)

    def cross_reference(self, model: BDF) -> None:
        msg = ', which is required by PAERO2 eid=%s' % self.pid
        if self.lrsb is not None and isinstance(self.lrsb, integer_types):
            self.lrsb_ref = model.AEFact(self.lrsb, msg=msg)
        if self.lrib is not None and isinstance(self.lrib, integer_types):
            self.lrib_ref = model.AEFact(self.lrib, msg=msg)

    def safe_cross_reference(self, model, xref_errors):
        msg = ', which is required by PAERO2 eid=%s' % self.pid
        if self.lrsb is not None and isinstance(self.lrsb, integer_types):
            self.lrsb_ref = model.safe_aefact(self.lrsb, self.pid, xref_errors, msg=msg)
        if self.lrib is not None and isinstance(self.lrib, integer_types):
            self.lrib_ref = model.safe_aefact(self.lrib, self.pid, xref_errors, msg=msg)

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.lrsb = self.Lrsb()
        self.lrib = self.Lrib()
        self.lrsb_ref = None
        self.lrib_ref = None

    def Lrsb(self):
        """AEFACT id"""
        if self.lrsb_ref is not None:
            return self.lrsb_ref.sid
        return self.lrsb

    def Lrib(self):
        """AEFACT id"""
        if self.lrib_ref is not None:
            return self.lrib_ref.sid
        return self.lrib

    def raw_fields(self):
        """
        Gets the fields in their unmodified form

        Returns
        -------
        fields : list[varies]
            the fields that define the card

        """
        list_fields = ['PAERO2', self.pid, self.orient, self.width,
                       self.AR, self.Lrsb(), self.Lrib()] + self.lth
        for (thi, thn) in zip(self.thi, self.thn):
            list_fields += [thi, thn]
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class PAERO3(BaseCard):
    """
    Defines the number of Mach boxes in the flow direction and the
    location of cranks and control surfaces of a Mach box lifting
    surface.

    +--------+------+------+-------+------+-----+------+------+------+
    |    1   |   2  |   3  |   4   |   5  |  6  |   7  |   8  |  9   |
    +========+======+======+=======+======+=====+======+======+======+
    | PAERO3 |  PID | NBOX | NCTRL |      |  X5 |  Y5  |  X6  |  Y6  |
    +--------+------+------+-------+------+-----+------+------+------+
    |        |  X7  |  Y7  |   X8  |  Y8  |  X9 |  Y9  |  X10 |  Y10 |
    +--------+------+------+-------+------+-----+------+------+------+
    |        |  X11 |  Y11 |  X12  |  Y12 |     |      |      |      |
    +--------+------+------+-------+------+-----+------+------+------+
    | PAERO3 | 2001 |  15  |   1   |      | 0.  |  65. |      |      |
    +--------+------+------+-------+------+-----+------+------+------+
    |        |  78. |  65. |  108. |  65. | 82. | 97.5 | 112. | 97.5 |
    +--------+------+------+-------+------+-----+------+------+------+
    |        |  86. | 130. |  116. | 130. |     |      |      |      |
    +--------+------+------+-------+------+-----+------+------+------+

    """
    type = 'PAERO3'
    _field_map = {
        1: 'pid', 2:'orient', 3:'width', 4:'AR',
    }
    _properties = ['npoints']

    def _get_field_helper(self, n):
        """
        Gets complicated parameters on the PAERO3 card

        Parameters
        ----------
        n : int
            the field number to update

        Returns
        -------
        value : varies
            the value for the appropriate field

        """
        nnew = n - 6
        if nnew < 0:
            raise RuntimeError('field n=%i on PAERO3 is invalid' % n)
        spot = nnew // 2
        i = nnew % 2
        if i == 0:
            value = self.x[spot]
        else:
            value = self.y[spot]
        return value

    def _update_field_helper(self, n, value):
        """
        Updates complicated parameters on the PAERO3 card

        Parameters
        ----------
        n : int
            the field number to update
        value :varies
            the value for the appropriate field

        """
        nnew = n - 6
        if nnew < 0:
            raise RuntimeError('field n=%i on PAERO3 is invalid' % n)
        spot = nnew // 2
        i = nnew % 2
        if i == 0:
            self.x[spot] = value
        else:
            self.y[spot] = value

    @classmethod
    def _init_from_empty(cls):
        pid = 1
        nbox = 1
        ncontrol_surfaces = 1
        x = [0., 0., 0.]
        y = [0., 10., 0.]
        return PAERO3(pid, nbox, ncontrol_surfaces, x, y, comment='')

    def __init__(self, pid, nbox, ncontrol_surfaces, x, y, comment=''):
        """
        Creates a PAERO3 card, which defines the number of Mach boxes
        in the flow direction and the location of cranks and control
        surfaces of a Mach box lifting surface.

        Parameters
        ----------
        pid : int
            PAERO1 id
        nbox : int
            Number of Mach boxes in the flow direction; 0 < nbox < 50
        ncontrol_surfaces : int
            Number of control surfaces. (0, 1, or 2)
        x / y : List[float, None]
            float : locations of points 5 through 12, which are in the
            aerodynamic coordinate system, to define the cranks and
            control surface geometry.
        comment : str; default=''
            a comment for the card

        """
        BaseCard.__init__(self)
        if comment:
            self.comment = comment
        #: Property identification number. (Integer > 0)
        self.pid = pid
        self.nbox = nbox
        self.ncontrol_surfaces = ncontrol_surfaces
        self.x = x
        self.y = y

    def validate(self):
        assert len(self.x) == len(self.y), 'nx=%s ny=%s' % (len(self.x), len(self.y))
        assert len(self.x) <= 8, 'nx=%s'  % len(self.x)
        for i, xi, yi in zip(count(), self.x, self.y):
            if xi is None or yi is None:
                assert xi == yi, 'x%i=%s y%i=%s must be None or floats' % (i, xi, i, yi)

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a PAERO3 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        pid = integer(card, 1, 'pid')
        nbox = integer(card, 2, 'nbox')
        ncontrol_surfaces = integer(card, 3, 'ncontrol_surfaces')
        x = []
        y = []
        nfields = card.nfields

        j = 5
        for i in range(5, nfields, 2):
            xi = double_or_blank(card, i, 'x%i' % j)
            yi = double_or_blank(card, i + 1, 'y%i' % j)
            x.append(xi)
            y.append(yi)
            j += 1
        return PAERO3(pid, nbox, ncontrol_surfaces, x, y, comment=comment)

    def cross_reference(self, model: BDF) -> None:
        pass

    def safe_cross_reference(self, model, unused_xref_errors):
        return self.cross_reference(model)

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        pass

    @property
    def npoints(self):
        return len(self.x)

    def raw_fields(self):
        """
        Gets the fields in their unmodified form

        Returns
        -------
        fields : list[varies]
            the fields that define the card

        """
        list_fields = ['PAERO3', self.pid, self.nbox, self.ncontrol_surfaces, None]
        for (x, y) in zip(self.x, self.y):
            list_fields += [x, y]
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        return self.comment + print_card_8(card)

class PAERO4(BaseCard):
    """
    Defines properties of each strip element for Strip theory.

    +--------+------+-------+--------+-------+-------+--------+--------+--------+
    |    1   |   2  |   3   |   4    |   5   |   6   |    7   |   8    |    9   |
    +========+======+=======+========+=======+=======+========+========+========+
    | PAERO4 | PID  | CLA   |  LCLA  |  CIRC | LCIRC |  DOC1  |  CAOC1 | GAPOC1 |
    +--------+------+-------+--------+-------+-------+--------+--------+--------+
    |        | DOC2 | CAOC2 | GAPOC2 |  DOC3 | CAOC3 | GAPOC3 |  etc.  |        |
    +--------+------+-------+--------+-------+-------+--------+--------+--------+
    | PAERO4 | 6001 |   1   |   501  |   0   |   0   |   0.0  |   0.0  |   0.0  |
    +--------+------+-------+--------+-------+-------+--------+--------+--------+
    |        | 0.50 |  0.25 |  0.02  |  0.53 |  0.24 |   0.0  |        |        |
    +--------+------+-------+--------+-------+-------+--------+--------+--------+
    ## TODO: what happens for DOC4?

    """
    type = 'PAERO4'
    _field_map = {
        1: 'pid', #2:'orient', 3:'width', 4:'AR',
    }

    #def _get_field_helper(self, n):
        #"""
        #Gets complicated parameters on the PAERO3 card

        #Parameters
        #----------
        #n : int
            #the field number to update

        #Returns
        #-------
        #value : varies
            #the value for the appropriate field
        #"""
        #nnew = n - 6
        #if nnew < 0:
            #raise RuntimeError('field n=%i on PAERO3 is invalid' % n)
        #spot = nnew // 2
        #i = nnew % 2
        #if i == 0:
            #value = self.x[spot]
        #else:
            #value = self.y[spot]
        #return value

    #def _update_field_helper(self, n, value):
        #"""
        #Updates complicated parameters on the PAERO3 card

        #Parameters
        #----------
        #n : int
            #the field number to update
        #value :varies
            #the value for the appropriate field
        #"""
        #nnew = n - 6
        #if nnew < 0:
            #raise RuntimeError('field n=%i on PAERO3 is invalid' % n)
        #spot = nnew // 2
        #i = nnew % 2
        #if i == 0:
            #self.x[spot] = value
        #else:
            #self.y[spot] = value

    @classmethod
    def _init_from_empty(cls):
        pid = 1
        docs = [1, 2]
        caocs = [1, 2]
        gapocs = [1, 2]
        return PAERO4(pid, docs, caocs, gapocs,
                      cla=0, lcla=0, circ=0, lcirc=0, comment='')

    def __init__(self, pid, docs, caocs, gapocs,
                 cla=0, lcla=0, circ=0, lcirc=0, comment=''):
        BaseCard.__init__(self)
        if comment:
            self.comment = comment
        #: Property identification number. (Integer > 0)
        self.pid = pid
        self.cla = cla
        self.lcla = lcla
        self.circ = circ
        self.lcirc = lcirc
        self.docs = docs
        self.caocs = caocs
        self.gapocs = gapocs
        assert isinstance(docs, list), docs
        assert isinstance(caocs, list), caocs
        assert isinstance(gapocs, list), gapocs

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a PAERO4 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        pid = integer(card, 1, 'pid')
        cla = integer_or_blank(card, 2, 'cla', 0)
        lcla = integer_or_blank(card, 3, 'lcla', 0) # ???

        circ = integer_or_blank(card, 4, 'circ', 0)
        lcirc = integer_or_blank(card, 5, 'lcirc', 0) # ???
        nfields = card.nfields

        j = 0
        docs = []
        caocs = []
        gapocs = []
        for i in range(6, nfields, 3):
            doc = double(card, i, 'doc_%i' % j)
            caoc = double(card, i + 1, 'caoc_%i' % j)
            gapoc = double(card, i + 2, 'gapoc_%i' % j)
            docs.append(doc)
            caocs.append(caoc)
            gapocs.append(gapoc)
            j += 1
        return PAERO4(pid, docs, caocs, gapocs,
                      cla=cla, lcla=lcla, circ=circ, lcirc=lcirc, comment=comment)

    def cross_reference(self, model: BDF) -> None:
        pass

    def safe_cross_reference(self, model: BDF, xref_errors):
        pass

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        pass

    def raw_fields(self):
        """
        Gets the fields in their unmodified form

        Returns
        -------
        fields : list[varies]
            the fields that define the card

        """
        list_fields = ['PAERO4', self.pid, self.cla, self.lcla, self.circ, self.lcirc]
        for doc, caoc, gapoc in zip(self.docs, self.caocs, self.gapocs):
            list_fields += [doc, caoc, gapoc]
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class Spline(BaseCard):
    def __init__(self):
        BaseCard.__init__(self)


class SPLINE1(Spline):
    """
    Surface Spline Methods
    Defines a surface spline for interpolating motion and/or forces for
    aeroelastic problems on aerodynamic geometries defined by regular
    arrays of aerodynamic points.

    +---------+-------+-------+------+------+------+----+------+-------+
    |    1    |   2   |    3  |   4  |   5  |   6  |  7 |   8  |   9   |
    +=========+=======+=======+======+======+======+====+======+=======+
    | SPLINE1 | EID   | CAERO | BOX1 | BOX2 | SETG | DZ | METH | USAGE |
    +---------+-------+-------+------+------+------+----+------+-------+
    |         | NELEM | MELEM |      |      |      |    |      |       |
    +---------+-------+-------+------+------+------+----+------+-------+
    | SPLINE1 |   3   |  111  | 115  | 122  |  14  | 0. |      |       |
    +---------+-------+-------+------+------+------+----+------+-------+

    """
    type = 'SPLINE1'
    _field_map = {
        1: 'eid', 2:'caero', 3:'box1', 4:'box2', 5:'setg', 6:'dz',
        7: 'method', 8:'usage', 9:'nelements', 10:'melements',
    }
    _properties = ['aero_element_ids']

    @classmethod
    def _init_from_empty(cls):
        eid = 1
        caero = 1
        box1 = 1
        box2 = 2
        setg = 3
        return SPLINE1(eid, caero, box1, box2, setg,
                       dz=0., method='IPS', usage='BOTH',
                       nelements=10, melements=10, comment='')

    def __init__(self, eid: int, caero: int, box1: int, box2: int, setg: int,
                 dz: float=0., method: str='IPS',
                 usage: str='BOTH', nelements: int=10,
                 melements: int=10, comment: str=''):
        """
        Creates a SPLINE1, which defines a surface spline.

        Parameters
        ----------
        eid : int
            spline id
        caero : int
            CAEROx id that defines the plane of the spline
        box1 / box2 : int
            First/last box id that is used by the spline
        setg : int
            SETx id that defines the list of GRID points that are used
            by the surface spline
        dz : float; default=0.0
            linear attachment flexibility
            dz = 0.; spline passes through all grid points
        method : str; default=IPS
            method for spline fit
            valid_methods = {IPS, TPS, FPS}
            IPS : Harder-Desmarais Infinite Plate Spline
            TPS : Thin Plate Spline
            FPS : Finite Plate Spline
        usage : str; default=BOTH
            Spline usage flag to determine whether this spline applies
            to the force transformation, displacement transformation, or
            both
            valid_usage = {FORCE, DISP, BOTH}
        nelements : int; default=10
            The number of FE elements along the local spline x-axis if
            using the FPS option
        melements : int; default=10
            The number of FE elements along the local spline y-axis if
            using the FPS option
        comment : str; default=''
            a comment for the card

        """
        Spline.__init__(self)
        if comment:
            self.comment = comment
        self.eid = eid
        self.caero = caero
        self.box1 = box1
        self.box2 = box2
        self.setg = setg
        self.dz = dz
        self.method = method
        self.usage = usage
        self.nelements = nelements
        self.melements = melements
        self.caero_ref = None
        self.setg_ref = None

    def validate(self):
        assert self.method in ['IPS', 'TPS', 'FPS'], 'method = %s' % self.method
        if self.method == 'FPS':
            assert self.nelements > 0, f'nelements={self.nelements} method={self.method}'
            assert self.melements > 0, f'melements={self.melements} method={self.method}'
        assert self.box2 >= self.box1, 'box1=%s box2=%s' % (self.box1, self.box2)
        assert self.usage in ['FORCE', 'DISP', 'BOTH'], 'usage = %s' % self.usage

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a SPLINE1 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        eid = integer(card, 1, 'eid')
        caero = integer(card, 2, 'caero')
        box1 = integer(card, 3, 'box1')
        box2 = integer(card, 4, 'box2')
        setg = integer(card, 5, 'setg')
        dz = double_or_blank(card, 6, 'dz', 0.0)
        method = string_or_blank(card, 7, 'method', 'IPS')
        usage = string_or_blank(card, 8, 'usage', 'BOTH')
        nelements = integer_or_blank(card, 9, 'nelements', 10)
        melements = integer_or_blank(card, 10, 'melements', 10)
        assert len(card) <= 11, 'len(SPLINE1 card) = %i\ncard=%s' % (len(card), card)
        return SPLINE1(eid, caero, box1, box2, setg, dz, method, usage,
                       nelements, melements, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        eid = data[0]
        caero = data[1]
        box1 = data[2]
        box2 = data[3]
        setg = data[4]
        dz = data[5]
        method = data[6]
        usage = data[7]
        nelements = data[8]
        melements = data[9]
        assert len(data) == 10, 'data = %s' % data
        return SPLINE1(eid, caero, box1, box2, setg, dz, method, usage,
                       nelements, melements, comment=comment)

    @property
    def aero_element_ids(self):
        return np.arange(self.box1, self.box2 + 1)

    def CAero(self):
        if self.caero_ref is not None:
            return self.caero_ref.eid
        return self.caero

    def Set(self):
        if self.setg_ref is not None:
            return self.setg_ref.sid
        return self.setg

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        msg = ', which is required by SPLINE1 eid=%s' % self.eid
        self.caero_ref = model.CAero(self.caero, msg=msg)
        self.setg_ref = model.Set(self.setg, msg=msg)

        if self.setg_ref.type == 'SET2':
            self.setg_ref.cross_reference_set(model, 'MACRO', msg=msg)
        else:
            self.setg_ref.cross_reference_set(model, 'Node', msg=msg)
            nnodes = len(self.setg_ref.ids)
            if nnodes < 3:
                msg = 'SPLINE1 requires at least 3 nodes; nnodes=%s\n' % (nnodes)
                msg += str(self)
                msg += str(self.setg_ref)
                raise RuntimeError(msg)

    def safe_cross_reference(self, model, xref_errors):
        msg = ', which is required by SPLINE1 eid=%s' % self.eid
        self.caero_ref = model.safe_caero(self.caero, self.eid, xref_errors, msg=msg)
        #self.setg_ref = model.safe_set(self, self.setg, self.eid, xref_errors, msg=msg)
    #def safe_set(self, setg, set_type, self.eid, xref_errors, msg=''):
        try:
            self.setg_ref = model.Set(self.setg, msg=msg)
            try:
                self.setg_ref.safe_cross_reference(model, 'Node', msg=msg)
            except:
                print(self.setg_ref)
                raise

            nnodes = len(self.setg_ref.ids)
            if nnodes < 3:
                msg = 'SPLINE1 requires at least 3 nodes; nnodes=%s\n' % (nnodes)
                msg += str(self)
                msg += str(self.setg_ref)
                model.log.warning(msg)
                msg = ''
        except KeyError:
            model.log.warning('failed to find SETx set_id=%s%s; allowed_sets=%s' % (
                self.setg, msg, np.unique(list(model.sets.keys()))))

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.caero = self.CAero()
        self.setg = self.Set()
        self.caero_ref = None
        self.setg_ref = None

    def raw_fields(self):
        """
        Gets the fields in their unmodified form

        Returns
        -------
        fields : list[varies]
            the fields that define the card

        """
        list_fields = ['SPLINE1', self.eid, self.CAero(), self.box1, self.box2,
                       self.Set(), self.dz, self.method, self.usage, self.nelements,
                       self.melements]
        return list_fields

    def repr_fields(self):
        dz = set_blank_if_default(self.dz, 0.)
        method = set_blank_if_default(self.method, 'IPS')
        usage = set_blank_if_default(self.usage, 'BOTH')
        nelements = set_blank_if_default(self.nelements, 10)
        melements = set_blank_if_default(self.melements, 10)

        list_fields = ['SPLINE1', self.eid, self.CAero(), self.box1, self.box2,
                       self.Set(), dz, method, usage, nelements, melements]
        list_fields = wipe_empty_fields(list_fields)
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class SPLINE2(Spline):
    """
    Linear Spline
    Defines a beam spline for interpolating motion and/or forces for
    aeroelastic problems on aerodynamic geometries defined by regular
    arrays of aerodynamic points.

    +---------+------+-------+-------+-------+------+----+------+-----+
    |    1    |   2  |   3   |   4   |   5   |  6   |  7 |   8  |  9  |
    +=========+======+=======+=======+=======+======+====+======+=====+
    | SPLINE2 | EID  | CAERO |  ID1  |  ID2  | SETG | DZ | DTOR | CID |
    +---------+------+-------+-------+-------+------+----+------+-----+
    |         | DTHX | DTHY  | None  | USAGE |      |    |      |     |
    +---------+------+-------+-------+-------+------+----+------+-----+
    | SPLINE2 |   5  |   8   |  12   | 24    | 60   | 0. | 1.0  |  3  |
    +---------+------+-------+-------+-------+------+----+------+-----+
    |         |  1.  |       |       |       |      |    |      |     |
    +---------+------+-------+-------+-------+------+----+------+-----+

    """
    type = 'SPLINE2'
    _field_map = {
        1: 'eid', 2:'caero', 3:'id1', 4:'id2', 5:'setg', 6:'dz',
        7: 'dtor', 8:'cid', 9:'dthx', 10:'dthy',
    }
    _properties = ['aero_element_ids']

    @classmethod
    def _init_from_empty(cls):
        eid = 1
        caero = 1
        box1 = 1
        box2 = 2
        setg = 1
        return SPLINE2(eid, caero, box1, box2, setg,
                       dz=0.0, dtor=1.0, cid=0, dthx=0., dthy=0., usage='BOTH', comment='')

    def __init__(self, eid: int, caero: int,
                 box1: int, box2: int, setg: int,
                 dz: float=0.0, dtor: float=1.0,
                 cid: int=0,
                 dthx: float=0.0, dthy: float=0.0,
                 usage: str='BOTH', comment: str=''):
        """
        Creates a SPLINE2 card, which defines a beam spline.

        Parameters
        ----------
        eid : int
            spline id
        caero : int
            CAEROx id that defines the plane of the spline
        box1 / box2 : int
            First/last box/body id that is used by the spline
        setg : int
            SETx id that defines the list of GRID points that are used
            by the beam spline
        dz : float; default=0.0
            linear attachment flexibility
            dz = 0.; spline passes through all grid points
        dtor : float; default=1.0
            Torsional flexibility ratio (EI/GJ).
            Use 1.0 for bodies (CAERO2).
        cid : int; default=0
            Rectangular coordinate system for which the y-axis defines the
            axis of the spline. Not used for bodies, CAERO2
        dthx : float; default=0.
            Rotational attachment flexibility.
            DTHX : Used for rotation about the spline's x-axis (in-plane
                   bending rotations).  It is not used for bodies (CAERO2).
            DTHY : Used for rotation about the spline's y-axis (torsion).
                   It is used for slope of bodies.
        usage : str; default=BOTH
            Spline usage flag to determine whether this spline applies
            to the force transformation, displacement transformation, or
            both
            valid_usage = {FORCE, DISP, BOTH}
        comment : str; default=''
            a comment for the card

        """
        Spline.__init__(self)
        if comment:
            self.comment = comment

        self.eid = eid
        self.caero = caero
        self.box1 = box1
        self.box2 = box2
        self.setg = setg
        self.dz = dz
        self.dtor = dtor
        self.cid = cid
        self.dthx = dthx
        self.dthy = dthy
        self.usage = usage
        self.cid_ref = None
        self.caero_ref = None
        self.setg_ref = None

    def validate(self):
        assert self.box2 >= self.box1, 'box2=%s box1=%s' % (self.box2, self.box1)

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a SPLINE2 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        eid = integer(card, 1, 'eid')
        caero = integer(card, 2, 'caero')
        id1 = integer(card, 3, 'id1')
        id2 = integer(card, 4, 'id2')
        setg = integer(card, 5, 'setg')
        dz = double_or_blank(card, 6, 'dz', 0.0)
        dtor = double_or_blank(card, 7, 'dtor', 1.0)
        cid = integer_or_blank(card, 8, 'cid', 0)
        dthx = double_or_blank(card, 9, 'dthx', 0.)
        dthy = double_or_blank(card, 10, 'dthy', 0.)

        usage = string_or_blank(card, 12, 'usage', 'BOTH')
        assert len(card) <= 13, 'len(SPLINE2 card) = %i\ncard=%s' % (len(card), card)
        return SPLINE2(eid, caero, id1, id2, setg, dz, dtor, cid,
                       dthx, dthy, usage, comment=comment)

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        msg = ', which is required by SPLINE2 eid=%s' % self.eid
        self.cid_ref = model.Coord(self.Cid(), msg=msg)
        self.caero_ref = model.CAero(self.CAero(), msg=msg)
        self.setg_ref = model.Set(self.Set(), msg=msg)

        if self.setg_ref.type == 'SET2':
            self.setg_ref.cross_reference_set(model, 'MACRO', msg=msg)
        else:
            self.setg_ref.cross_reference_set(model, 'Node', msg=msg)
            nnodes = len(self.setg_ref.ids)
            if nnodes < 2:
                msg = 'SPLINE2 requires at least 2 nodes; nnodes=%s\n' % (nnodes)
                msg += str(self)
                msg += str(self.setg_ref)
                raise RuntimeError(msg)

    def safe_cross_reference(self, model, xref_errors):
        msg = ', which is required by SPLINE2 eid=%s' % self.eid
        self.cid_ref = model.safe_coord(self.Cid(), self.eid, xref_errors, msg=msg)

        try:
            self.caero_ref = model.CAero(self.CAero(), msg=msg)
        except KeyError:
            pass

        try:
            self.setg_ref = model.Set(self.Set(), msg=msg)
            if self.setg_ref.type == 'SET2':
                self.setg_ref.cross_reference_set(model, 'MACRO', msg=msg)
            else:
                self.setg_ref.cross_reference_set(model, 'Node', msg=msg)

                nnodes = len(self.setg_ref.ids)
                if nnodes < 2:
                    msg = 'SPLINE2 requires at least 2 nodes; nnodes=%s\n' % (nnodes)
                    msg += str(self)
                    msg += str(self.setg_ref)
                    raise RuntimeError(msg)
        except KeyError:
            pass

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.cid = self.Cid()
        self.caero = self.CAero()
        self.setg = self.Set()
        self.cid_ref = None
        self.caero_ref = None
        self.setg_ref = None

    @property
    def aero_element_ids(self):
        return np.arange(self.box1, self.box2 + 1)

    def Cid(self):
        if self.setg_ref is not None:
            return self.cid_ref.cid
        return self.cid

    def CAero(self):
        if self.setg_ref is not None:
            return self.caero_ref.eid
        return self.caero

    def Set(self):
        if self.setg_ref is not None:
            return self.setg_ref.sid
        return self.setg

    def raw_fields(self):
        """
        Gets the fields in their unmodified form

        Returns
        -------
        fields : list[varies]
            the fields that define the card

        """
        list_fields = ['SPLINE2', self.eid, self.CAero(), self.box1, self.box2,
                       self.Set(), self.dz, self.dtor, self.Cid(), self.dthx,
                       self.dthy, None, self.usage]
        return list_fields

    def repr_fields(self):
        dz = set_blank_if_default(self.dz, 0.)
        usage = set_blank_if_default(self.usage, 'BOTH')
        list_fields = ['SPLINE2', self.eid, self.CAero(), self.box1, self.box2,
                       self.Set(), dz, self.dtor, self.Cid(), self.dthx, self.dthy,
                       None, usage]
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class SPLINE3(Spline):
    """
    Defines a constraint equation for aeroelastic problems.
    Useful for control surface constraints.

    +---------+------+-------+-------+------+----+----+-----+-------+
    |    1    |  2   |   3   |   4   |  5   |  6 |  7 |  8  |   9   |
    +=========+======+=======+=======+======+====+====+=====+=======+
    | SPLINE3 | EID  | CAERO | BOXID | COMP | G1 | C1 | A1  | USAGE |
    +---------+------+-------+-------+------+----+----+-----+-------+
    |         |  G2  |  C2   |  A2   |      | G3 | C3 | A2  |       |
    +---------+------+-------+-------+------+----+----+-----+-------+
    |         |  G4  |  C4   |  A4   | etc. |    |    |     |       |
    +---------+------+-------+-------+------+----+----+-----+-------+
    | SPLINE3 | 7000 |  107  |  109  |  6   | 5  | 3  | 1.0 | BOTH  |
    +---------+------+-------+-------+------+----+----+-----+-------+
    |         |  43  |   5   | -1.0  |      |    |    |     |       |
    +---------+------+-------+-------+------+----+----+-----+-------+

    """
    type = 'SPLINE3'
    _properties = ['node_ids']
    _field_map = {
        1: 'eid', 2:'caero', 3:'box_id',
        7: 'a1', 8:'usage',
    }
        #5:'g1', 6:'c1',
        #9:G2,C2,A2...

    @classmethod
    def _init_from_empty(cls):
        eid = 1
        caero = 1
        box_id = 1
        components = 1
        nodes = [2, 3]
        displacement_components = [2, 3]
        coeffs = [2., 3.]
        return SPLINE3(eid, caero, box_id, components, nodes, displacement_components, coeffs,
                       usage='BOTH')

    def __init__(self, eid: int, caero: int, box_id: int,
                 components: int,
                 nodes: List[int],
                 displacement_components: List[int],
                 coeffs: List[float],
                 usage: str='BOTH', comment: str=''):
        """
        Creates a SPLINE3 card, which is useful for control surface
        constraints.

        Parameters
        ----------
        eid : int
            spline id
        caero : int
            CAEROx id that defines the plane of the spline
        box_id : int
           Identification number of the aerodynamic box number.
        components : int
           The component of motion to be interpolated.
           3, 5          (CAERO1)
           2, 3, 5, 6    (CAERO2)
           3             (CAERO3)
           3, 5, 6       (CAERO4)
           3, 5, 6       (CAERO5)
           1, 2, 3, 5, 6 (3D Geometry)
           2-lateral displacement
           3-transverse displacement
           5-pitch angle
           6-relative control angle for CAERO4/5; yaw angle for CAERO2
        nodes : List[int]
           Grid point identification number of the independent grid point.
        displacement_components :  : List[int]
           Component numbers in the displacement coordinate system.
           1-6 (GRIDs)
           0 (SPOINTs)
        coeffs : List[float]
           Coefficient of the constraint relationship.
        usage : str; default=BOTH
            Spline usage flag to determine whether this spline applies
            to the force transformation, displacement transformation, or
            both
            valid_usage = {FORCE, DISP, BOTH}
        comment : str; default=''
            a comment for the card

        """
        Spline.__init__(self)
        if comment:
            self.comment = comment

        if isinstance(nodes, integer_types):
            nodes = [nodes]
        if isinstance(displacement_components, integer_types):
            displacement_components = [displacement_components]
        if isinstance(coeffs, float):
            coeffs = [coeffs]

        self.eid = eid
        self.caero = caero
        self.box_id = box_id
        #if isinstance(components, integer_types):
            #components = str(components)
        self.components = components
        self.usage = usage

        self.nodes = nodes
        self.displacement_components = displacement_components
        self.coeffs = coeffs
        self.nodes_ref = None
        self.caero_ref = None

    def validate(self):
        msg = ''
        if self.components not in [0, 1, 2, 3, 4, 5, 6]:
            msg += 'components=%r must be [0, 1, 2, 3, 4, 5, 6]\n' % (
                self.components)

        if not len(self.nodes) == len(self.displacement_components):
            msg += 'nnodes=%s ndisplacement_components=%s must be equal\n' % (
                len(self.nodes), len(self.displacement_components))
        if not len(self.nodes) == len(self.coeffs):
            msg += 'nnodes=%s ncoeffs=%s must be equal\n' % (
                len(self.nodes), len(self.coeffs))

        for i, disp_component  in enumerate(self.displacement_components):
            if disp_component not in [0, 1, 2, 3, 4, 5, 6]:
                msg += 'i=%s displacement_component=%s must be [0, 1, 2, 3, 4, 5, 6]\n' % (
                    i, disp_component)

        if self.usage not in ['FORCE', 'DISP', 'BOTH']:
            msg += 'usage=%r must be in [FORCE, DISP, BOTH]\n' % self.usage

        if msg:
            msg += str(self)
            raise RuntimeError(msg)

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a SPLINE3 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        eid = integer(card, 1, 'eid')
        caero = integer(card, 2, 'caero')
        box_id = integer(card, 3, 'box_id')
        components = integer(card, 4, 'comp')
        node = integer(card, 5, 'G1')
        coeff = integer(card, 6, 'C1')
        displacement_component = double(card, 7, 'A1')
        usage = string_or_blank(card, 8, 'usage', 'BOTH')

        nfields = len(card) - 1
        nrows = nfields // 8
        if nfields % 8:
            nrows += 1

        nodes = [node]
        coeffs = [coeff]
        displacement_components = [displacement_component]
        i = 2
        for irow in range(1, nrows):
            #print('G%i' % i)
            j = 1 + irow * 8
            node = integer(card, j, 'G%i' % i)
            coeff = integer(card, j + 1, 'C%i' % i)
            displacement_component = double(card, j + 2, 'A%i' % i)
            nodes.append(node)
            coeffs.append(coeff)
            displacement_components.append(displacement_component)
            i += 1
            if card.field(j + 4) or card.field(j + 5) or card.field(j + 6):
                node = integer(card, j + 4, 'G%i' % i)
                coeff = parse_components(card, j + 5, 'C%i' % i)
                displacement_component = double(card, j + 6, 'A%i' % i)
                nodes.append(node)
                coeffs.append(coeff)
                displacement_components.append(displacement_component)
                i += 1
        return SPLINE3(eid, caero, box_id, components,
                       nodes, coeffs, displacement_components, usage,
                       comment=comment)

    def cross_reference(self, model: BDF) -> None:
        msg = ', which is required by SPLINE3 eid=%s' % self.eid
        self.nodes_ref = model.Nodes(self.nodes, msg=msg)
        self.caero_ref = model.CAero(self.caero, msg=msg)

    def safe_cross_reference(self, model, xref_errors):
        msg = ', which is required by SPLINE3 eid=%s' % self.eid
        self.nodes_ref = model.Nodes(self.nodes, msg=msg)
        self.caero_ref = model.safe_caero(self.caero, self.eid, xref_errors, msg=msg)

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.caero = self.CAero()
        self.nodes = self.node_ids
        self.nodes_ref = None
        self.caero_ref = None

    def CAero(self):
        if self.caero_ref is not None:
            return self.caero_ref.eid
        return self.caero

    @property
    def node_ids(self):
        if self.nodes_ref is None:
            return self.nodes
        return [node.nid for node in self.nodes_ref]

    def raw_fields(self):
        """
        +---------+------+-------+-------+------+----+----+-----+-------+
        |    1    |  2   |   3   |   4   |  5   |  6 |  7 |  8  |   9   |
        +=========+======+=======+=======+======+====+====+=====+=======+
        | SPLINE3 | EID  | CAERO | BOXID | COMP | G1 | C1 | A1  | USAGE |
        +---------+------+-------+-------+------+----+----+-----+-------+
        |         |  G2  |  C2   |  A2   | ---- | G3 | C3 | A2  |  ---  |
        +---------+------+-------+-------+------+----+----+-----+-------+
        |         |  G4  |  C4   |  A4   | etc. |    |    |     |       |
        +---------+------+-------+-------+------+----+----+-----+-------+

        """
        list_fields = [
            'SPLINE3', self.eid, self.CAero(), self.box_id, self.components,
            self.nodes[0], self.displacement_components[0], self.coeffs[0], self.usage]
        for nid, disp_c, coeff in zip(self.nodes[1:], self.displacement_components[1:],
                                      self.coeffs[1:]):
            list_fields += [nid, disp_c, coeff, None]
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class SPLINE4(Spline):
    """
    Surface Spline Methods
    Defines a curved surface spline for interpolating motion and/or forces for
    aeroelastic problems on general aerodynamic geometries using either the
    Infinite Plate, Thin Plate or Finite Plate splining method.

    +---------+-------+-------+--------+-----+------+----+------+-------+
    |    1    |   2   |   3   |    4   |  5  |   6  |  7 |   8  |   9   |
    +=========+=======+=======+========+=====+======+====+======+=======+
    | SPLINE4 |  EID  | CAERO | AELIST |     | SETG | DZ | METH | USAGE |
    +---------+-------+-------+--------+-----+------+----+------+-------+
    |         | NELEM | MELEM |        |     |      |    |      |       |
    +---------+-------+-------+--------+-----+------+----+------+-------+
    | SPLINE4 |   3   | 111   |   115  |     |  14  | 0. | IPS  |       |
    +---------+-------+-------+--------+-----+------+----+------+-------+

    """
    type = 'SPLINE4'
    _properties = ['aero_element_ids']
    _field_map = {
        1: 'eid', 2:'caero', 3:'aelist', 5:'setg', 6:'dz',
        7: 'method', 8:'usage', 9:'nelements', 10:'melements',
    }

    @classmethod
    def _init_from_empty(cls):
        eid = 1
        caero = 1
        aelist = 1
        setg = 1
        dz = 1.
        method = 'IPS'
        usage = 'BOTH'
        nelements = 1
        melements = 1
        return SPLINE4(eid, caero, aelist, setg, dz, method, usage, nelements, melements,
                       comment='')

    def __init__(self, eid: int, caero: int, aelist: int, setg: int,
                 dz: float, method: str, usage: str,
                 nelements: int, melements: int, comment: str=''):
        """
        Creates a SPLINE4 card, which defines a curved Infinite Plate,
        Thin Plate, or Finite Plate Spline.

        Parameters
        ----------
        eid : int
            spline id
        caero : int
            CAEROx id that defines the plane of the spline
        box1 / box2 : int
            First/last box id that is used by the spline
        setg : int
            SETx id that defines the list of GRID points that are used
            by the surface spline
        dz : float; default=0.0
            linear attachment flexibility
            dz = 0.; spline passes through all grid points
        method : str; default=IPS
            method for spline fit
            valid_methods = {IPS, TPS, FPS}
            IPS : Harder-Desmarais Infinite Plate Spline
            TPS : Thin Plate Spline
            FPS : Finite Plate Spline
        usage : str; default=BOTH
            Spline usage flag to determine whether this spline applies
            to the force transformation, displacement transformation, or
            both
            valid_usage = {FORCE, DISP, BOTH}
        nelements / melements : int; default=10
            The number of FE elements along the local spline x/y-axis if
            using the FPS option
        comment : str; default=''
            a comment for the card

        """
        Spline.__init__(self)
        if comment:
            self.comment = comment
        self.eid = eid
        self.caero = caero
        self.aelist = aelist
        self.setg = setg
        self.dz = dz
        self.method = method
        self.usage = usage
        self.nelements = nelements
        self.melements = melements
        self.caero_ref = None
        self.setg_ref = None
        self.aelist_ref = None

    def validate(self):
        assert self.method in ['IPS', 'TPS', 'FPS'], 'method = %s' % self.method
        assert self.usage in ['FORCE', 'DISP', 'BOTH'], 'uasge = %s' % self.usage

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a SPLINE4 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        eid = integer(card, 1, 'eid')
        caero = integer(card, 2, 'caero')
        aelist = integer(card, 3, 'aelist')
        # None
        setg = integer(card, 5, 'setg')
        dz = double_or_blank(card, 6, 'dz', 0.0)
        method = string_or_blank(card, 7, 'method', 'IPS')
        usage = string_or_blank(card, 8, 'usage', 'BOTH')
        nelements = integer_or_blank(card, 9, 'nelements', 10)
        melements = integer_or_blank(card, 10, 'melements', 10)
        assert len(card) <= 11, 'len(SPLINE4 card) = %i\ncard=%s' % (len(card), card)
        return SPLINE4(eid, caero, aelist, setg, dz, method, usage,
                       nelements, melements, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        eid = data[0]
        caero = data[1]
        aelist = data[2]
        setg = data[3]
        dz = data[4]
        method = data[5]
        usage = data[6]
        nelements = data[7]
        melements = data[8]
        assert len(data) == 9, 'data = %s' % (data)
        return SPLINE4(eid, caero, aelist, setg, dz, method, usage,
                       nelements, melements, comment=comment)

    def CAero(self) -> int:
        if self.caero_ref is not None:
            return self.caero_ref.eid
        return self.caero

    def AEList(self) -> int:
        if self.aelist_ref is not None:
            return self.aelist_ref.sid
        return self.aelist

    def Set(self) -> int:
        if self.setg_ref is not None:
            return self.setg_ref.sid
        return self.setg

    @property
    def aero_element_ids(self) -> List[int]:
        return self.aelist_ref.elements

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        msg = ', which is required by SPLINE4 eid=%s' % self.eid
        self.caero_ref = model.CAero(self.CAero(), msg=msg)
        self.setg_ref = model.Set(self.Set(), msg=msg)
        self.aelist_ref = model.AEList(self.aelist, msg=msg)

        if self.setg_ref.type == 'SET2':
            self.setg_ref.cross_reference_set(model, 'MACRO', msg=msg)
        else:
            self.setg_ref.cross_reference_set(model, 'Node', msg=msg)
            nnodes = len(self.setg_ref.ids)
            if nnodes < 3:
                msg = 'SPLINE4 requires at least 3 nodes; nnodes=%s\n' % (nnodes)
                msg += str(self)
                msg += str(self.setg_ref)
                raise RuntimeError(msg)


    def safe_cross_reference(self, model: BDF, xref_errors):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        msg = ', which is required by SPLINE4 eid=%s' % self.eid
        self.caero_ref = model.safe_caero(self.caero, self.eid, xref_errors, msg=msg)
        self.aelist_ref = model.safe_aelist(self.aelist, self.eid, xref_errors, msg=msg)
        self.setg_ref = model.Set(self.Set(), msg=msg)

        if self.setg_ref.type == 'SET2':
            self.setg_ref.cross_reference_set(model, 'MACRO', msg=msg)
        else:
            self.setg_ref.cross_reference_set(model, 'Node', msg)

            nnodes = len(self.setg_ref.ids)
            if nnodes < 3:
                msg = 'SPLINE4 requires at least 3 nodes; nnodes=%s\n' % (nnodes)
                msg += str(self)
                msg += str(self.setg_ref)
                raise ValueError(msg)

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.caero = self.CAero()
        self.setg = self.Set()
        self.aelist = self.AEList()
        self.caero_ref = None
        self.setg_ref = None
        self.aelist_ref = None

    def raw_fields(self):
        """
        Gets the fields in their unmodified form

        Returns
        -------
        fields : list[varies]
            the fields that define the card

        """
        list_fields = ['SPLINE4', self.eid, self.CAero(), self.AEList(), None,
                       self.Set(), self.dz, self.method, self.usage, self.nelements,
                       self.melements]
        return list_fields

    def repr_fields(self):
        dz = set_blank_if_default(self.dz, 0.)
        method = set_blank_if_default(self.method, 'IPS')
        usage = set_blank_if_default(self.usage, 'BOTH')
        nelements = set_blank_if_default(self.nelements, 10)
        melements = set_blank_if_default(self.melements, 10)

        list_fields = ['SPLINE4', self.eid, self.CAero(), self.AEList(), None,
                       self.Set(), dz, method, usage, nelements, melements]
        list_fields = wipe_empty_fields(list_fields)
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class SPLINE5(Spline):
    """
    Linear Spline
    Defines a 1D beam spline for interpolating motion and/or forces for
    aeroelastic problems on aerodynamic geometries defined by irregular arrays
    of aerodynamic points. The interpolating beam supports axial rotation and
    bending in the yz-plane.

    +=========+======+=======+========+=======+======+====+=======+=======+
    |    1    |  2   |    3  |    4   |   5   |   6  |  7 |   8   |   9   |
    +=========+======+=======+========+=======+======+====+=======+=======+
    | SPLINE5 | EID  | CAERO | AELIST |       | SETG | DZ | DTOR  |  CID  |
    +---------+------+-------+--------+-------+------+----+-------+-------+
    |         | DTHX | DTHY  |        | USAGE | METH |    | FTYPE | RCORE |
    +---------+------+-------+--------+-------+------+----+-------+-------+

    METH, FTYPE, RCORE are in 2012+ (not MSC.2005r2 or NX.10)

    """
    type = 'SPLINE5'
    _properties = ['aero_element_ids']
    _field_map = {
        1: 'eid', 2:'caero', 3:'aelist', 5:'setg', 6:'dz',
        7: 'dtor', 8:'cid', 9:'thx', 10:'thy', 12:'usage',
        13 : 'meth', 15 : 'ftype', 16 : 'rcore',
    }

    @classmethod
    def _init_from_empty(cls):
        eid = 1
        caero = 1
        aelist = 1
        setg = 1
        thx = 1.
        thy = 1.
        return SPLINE5(eid, caero, aelist, setg, thx, thy,
                       dz=0., dtor=1.0, cid=0, usage='BOTH', method='BEAM',
                       ftype='WF2', rcore=None, comment='')

    def __init__(self, eid, caero, aelist, setg, thx, thy, dz=0., dtor=1.0, cid=0,
                 usage='BOTH', method='BEAM', ftype='WF2', rcore=None, comment=''):
        Spline.__init__(self)
        if comment:
            self.comment = comment

        self.eid = eid
        self.caero = caero
        self.aelist = aelist
        self.setg = setg
        self.dz = dz
        self.dtor = dtor
        self.cid = cid
        self.thx = thx
        self.thy = thy
        self.usage = usage
        self.method = method
        self.ftype = ftype
        self.rcore = rcore
        self.cid_ref = None
        self.caero_ref = None
        self.setg_ref = None
        self.aelist_ref = None

    def validate(self):
        assert isinstance(self.cid, integer_types), self.cid
        assert self.method in ['BEAM', 'RIS'], self.method
        assert self.ftype in ['WF0', 'WF2'], self.ftype

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a SPLINE5 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        eid = integer(card, 1, 'eid')
        caero = integer(card, 2, 'caero')
        aelist = integer(card, 3, 'aelist')
        # None
        setg = integer(card, 5, 'setq')
        dz = double_or_blank(card, 6, 'dz', 0.0)
        dtor = double_or_blank(card, 7, 'dtor', 1.0)
        cid = integer_or_blank(card, 8, 'cid', 0)
        thx = double(card, 9, 'thx')
        thy = double(card, 10, 'thy')
        usage = string_or_blank(card, 12, 'usage', 'BOTH')
        # per nast/tpl/fmondsp.dat, METH can be a double(0.0) ???
        method = string_or_blank(card, 13, 'meth', 'BEAM')
        ftype = string_or_blank(card, 15, 'ftype', 'WF2')
        rcore = double_or_blank(card, 16, 'rcore')

        usage = string_or_blank(card, 12, 'usage', 'BOTH')
        assert len(card) <= 16, 'len(SPLINE5 card) = %i\n%s' % (len(card), card)
        return SPLINE5(eid, caero, aelist, setg, thx, thy, dz=dz, dtor=dtor, cid=cid,
                       usage=usage, method=method, ftype=ftype, rcore=rcore, comment=comment)
    @property
    def aero_element_ids(self):
        return self.aelist_ref.elements


    def Cid(self):
        if self.cid_ref is not None:
            return self.cid_ref.cid
        return self.cid

    def CAero(self):
        if self.caero_ref is not None:
            return self.caero_ref.eid
        return self.caero

    def AEList(self):
        if self.aelist_ref is not None:
            return self.aelist_ref.sid
        return self.aelist

    def Set(self):
        if self.setg_ref is not None:
            return self.setg_ref.sid
        return self.setg

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        msg = ', which is required by SPLINE5 eid=%s' % self.eid
        self.cid_ref = model.Coord(self.cid, msg=msg)
        self.caero_ref = model.CAero(self.caero, msg=msg)
        self.setg_ref = model.Set(self.setg, msg=msg)
        self.aelist_ref = model.AEList(self.aelist, msg=msg)

        if self.setg_ref.type == 'SET2':
            self.setg_ref.cross_reference_set(model, 'MACRO', msg=msg)
        else:
            self.setg_ref.cross_reference_set(model, 'Node', msg=msg)
            nnodes = len(self.setg_ref.ids)
            if nnodes < 3:
                msg = 'SPLINE1 requires at least 3 nodes; nnodes=%s\n' % (nnodes)
                msg += str(self)
                msg += str(self.setg_ref)
                raise RuntimeError(msg)

    def safe_cross_reference(self, model, xref_errors):
        msg = ', which is required by SPLINE5 eid=%s' % self.eid
        self.cid_ref = model.safe_coord(self.cid, self.eid, xref_errors, msg=msg)
        self.caero_ref = model.safe_caero(self.caero, self.eid, xref_errors, msg=msg)

        try:
            self.setg_ref = model.Set(self.setg, msg=msg)
            nnodes = len(self.setg_ref.ids)
            if nnodes < 3:
                msg = 'SPLINE5 requires at least 3 nodes; nnodes=%s\n' % (nnodes)
                msg += str(self)
                msg += str(self.setg_ref)
                raise RuntimeError(msg)
        except KeyError:
            pass

        try:
            self.setg_ref.cross_reference_set(model, 'Node', msg)
            self.aelist_ref = model.AEList(self.aelist, msg=msg)
        except KeyError:
            pass

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.cid = self.Cid()
        self.caero = self.CAero()
        self.setg = self.Set()
        self.aelist = self.AEList()
        self.caero_ref = None
        self.setg_ref = None
        self.aelist_ref = None

    def raw_fields(self):
        """
        Gets the fields in their unmodified form

        Returns
        -------
        fields : list[varies]
            the fields that define the card

        """
        list_fields = ['SPLINE5', self.eid, self.CAero(), self.AEList(), None,
                       self.Set(), self.dz, self.dtor, self.Cid(), self.thx, self.thy,
                       None, self.usage, self.method, None, self.ftype, self.rcore]
        return list_fields

    def repr_fields(self):
        dz = set_blank_if_default(self.dz, 0.)
        usage = set_blank_if_default(self.usage, 'BOTH')

        list_fields = ['SPLINE5', self.eid, self.CAero(), self.AEList(), None,
                       self.Set(), dz, self.dtor, self.Cid(), self.thx, self.thy,
                       None, usage, self.method, None, self.ftype, self.rcore]
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        return self.comment + print_card_8(card)
