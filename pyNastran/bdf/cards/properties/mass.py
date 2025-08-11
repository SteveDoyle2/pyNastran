# pylint: disable=C0103
"""
All mass properties are defined in this file.  This includes:

 * NSM
 * PMASS

All mass properties are PointProperty and Property objects.

"""
from __future__ import annotations
import warnings
from typing import TYPE_CHECKING

import numpy as np

from pyNastran.bdf.cards.base_card import expand_thru_by, expand_thru, BaseCard, Property
from pyNastran.bdf.bdf_interface.assign_type import (
    integer, integer_or_string, double, double_or_blank, string)
from pyNastran.bdf.bdf_interface.assign_type_force import (
    force_double_or_blank)
from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.field_writer_16 import print_card_16
from pyNastran.utils.numpy_utils import integer_types, float_types
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf import BDF
    from pyNastran.bdf.bdf_interface.bdf_card import BDFCard


class NSMx(Property):
    """Common class for NSM and NSML"""
    _field_map = {
        1: 'sid', 2:'Type', 3:'id', 4:'value'
    }

    #: Set points to either Property entries or Element entries.
    #: Properties are:
    valid_properties = [
        'PSHELL', 'PCOMP', 'PCOMPG', 'PBAR', 'PBARL', 'PBEAM', 'PBEAML', 'PBCOMP',
        'PROD', 'CONROD', 'PBEND', 'PSHEAR', 'PTUBE', 'PRAC2D', # 'PCONEAX',
        'ELEMENT', 'PDUM8',
    ]

    def __init__(self, sid: int, nsm_type: str, pid_eid: int, value: float, comment: str=''):
        """
        Creates an NSM/NSM1 card

        Parameters
        ----------
        sid : int
            Case control NSM id
        nsm_type : str
            Type of card the NSM is applied to
            valid_properties = {
                PSHELL, PCOMP, PBAR, PBARL, PBEAM, PBEAML, PBCOMP,
                PROD, CONROD, PBEND, PSHEAR, PTUBE, PCONEAX, PRAC2D,
                ELEMENT
            }
        pid_eid : int
            property id or element id depending on nsm_type
        value : float
            the non-structural pass per unit length/area
        comment : str; default=''
            a comment for the card
        """
        Property.__init__(self)
        if comment:
            self.comment = comment
        self.sid = sid
        self.nsm_type = nsm_type
        self.id = pid_eid
        self.value = value
        assert isinstance(pid_eid, integer_types), 'pid_eid=%s type=%s' % (pid_eid, type(pid_eid))
        assert isinstance(value, float_types), 'value=%s type=%s' % (value, type(value))

        if self.nsm_type not in self.valid_properties:
            raise TypeError('nsm_type=%r must be in [%s]' % (
                self.nsm_type, ', '.join(self.valid_properties)))

    @property
    def ids(self):
        return [self.id]

    @classmethod
    def add_card(cls, card, icard=0, comment=''):
        """
        Adds an NSM card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        icard : int; default=0
            the index of the card that's being parsed
        comment : str; default=''
            a comment for the card
        """
        noffset = 2 * icard
        sid = integer(card, 1, 'sid')
        nsm_type = string(card, 2, 'Type')
        pid_eid = integer(card, 3 + noffset, 'pid/eid')
        value = double(card, 4 + noffset, 'value')
        return cls(sid, nsm_type, pid_eid, value, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        """
        Adds an NSM/NSM1 card from the OP2

        Parameters
        ----------
        data : list[varies]
            a list of fields defined in OP2 format
        comment : str; default=''
            a comment for the card
        """
        sid = data[0]
        #sid=9  prop_set=PBEA ID=538976333 value=0.0
        #sid=10 prop_set=PDUM ID=538976312 value=2.80259692865e-45
        #sid=10 prop_set=ELEM ID=542395973 value=0.0
        nsm_type = data[1]
        pid_eid = data[2]
        value = data[3]
        return cls(sid, nsm_type, pid_eid, value, comment=comment)

    def cross_reference(self, model: BDF) -> None:
        pass



    def raw_fields(self):
        #nodes = self.node_ids
        list_fields = [self.type, self.sid, self.nsm_type, self.id, self.value]
        return list_fields

    def repr_fields(self):
        return self.raw_fields()

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)


class NSM1x(Property):
    """Common class for NSM1 and NSML1"""
    valid_properties = [
        'PSHELL', 'PCOMP', 'PCOMPG', 'PBAR', 'PBARL', 'PBEAM', 'PBEAML', 'PBCOMP',
        'PROD', 'CONROD', 'PBEND', 'PSHEAR', 'PTUBE', 'PRAC2D', # 'PCONEAX',
        'PBUSH',
        'ELEMENT',
    ]

    def __init__(self, sid, nsm_type, value, ids, comment=''):
        """
        Creates an NSM1/NSML1 card

        Parameters
        ----------
        sid : int
            Case control NSM id
        nsm_type : str
            Type of card the NSM is applied to
            valid_properties = {
                PSHELL, PCOMP, PBAR, PBARL, PBEAM, PBEAML, PBCOMP,
                PROD, CONROD, PBEND, PSHEAR, PTUBE, PCONEAX, PRAC2D,
                ELEMENT
            }
        value : float
            NSM1:  the non-structural pass per unit length/area
            NSML1: the total non-structural pass per unit length/area;
                   the nsm will be broken down based on a weighted area/length
        ids : list[int]
            property ids or element ids depending on nsm_type
        comment : str; default=''
            a comment for the card
        """
        Property.__init__(self)
        if comment:
            self.comment = comment

        if isinstance(ids, integer_types):
            ids = [ids]
        elif isinstance(ids, str):
            assert ids == 'ALL', f'ids={ids!r} is not ALL'
            ids = [ids]
        elif isinstance(ids, list):
            pass
        elif isinstance(ids, np.ndarray):
            ids = ids.tolist()
        else:  # pragma: no cover
            raise TypeError((ids, type(ids)))

        if ids == ['ALL']:
            pass
        else:
            # With the 'THRU' and 'THRU', 'BY' forms, blanks fields are
            # allowed for readability. Any combination of a list of IDs
            # and 'THRU' and 'THRU', 'BY' is allowed. The "THRU" and "BY"
            # lists may have missing IDs. That is the list of IDs in a
            # THRU range need not be continuous.
            ids = expand_thru_by(ids, allow_blanks=True)
        assert len(ids) > 0, ids
        self.sid = sid
        self.nsm_type = nsm_type

        # TODO: expand the node ids
        self.ids = ids
        self.value = value
        #print(str(self))
        if self.nsm_type not in self.valid_properties:
            raise TypeError('nsm_type=%r must be in [%s]' % (
                self.nsm_type, ', '.join(self.valid_properties)))
        assert isinstance(ids, list), 'ids=%r is not a list' % (ids)

    @property
    def Type(self):
        """gets the nsm_type"""
        return self.nsm_type

    @Type.setter
    def Type(self, nsm_type):
        """sets the nsm_type"""
        self.nsm_type = nsm_type

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds an NSM1/NSML1 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        sid = integer(card, 1, 'sid')
        nsm_type = string(card, 2, 'Type')
        value = double(card, 3, 'value')

        # TODO: doesn't support 1 THRU 11 BY 2
        ids = []
        _id = 1
        nfields = len(card)
        if nfields == 5:
            id1 = integer_or_string(card, 4, 'ID_1')
            if id1 != 'ALL' and not isinstance(id1, int):
                msg = ('*ID_1 = %r (field #4) on card must be an integer or ALL.\n'
                       'card=%s' % (id1, card))
                raise SyntaxError(msg)
            ids = id1
        else:
            # we'll handle expansion in the init
            ids = card[4:]
        return cls(sid, nsm_type, value, ids, comment=comment)

    def cross_reference(self, model: BDF) -> None:
        pass

    def raw_fields(self):
        list_fields = [self.type, self.sid, self.nsm_type, self.value] + self.ids
        return list_fields

    def repr_fields(self):
        return self.raw_fields()

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)


class NSM1(NSM1x):
    """
    Defines a set of non structural mass.

    +------+-----+------+-------+-----+----+----+----+----+
    |  1   |  2  |  3   |   4   |  5  | 6  | 7  | 8  | 9  |
    +======+=====+======+=======+=====+====+====+====+====+
    | NSM1 | SID | TYPE | VALUE | ID  | ID | ID | ID | ID |
    +------+-----+------+-------+-----+----+----+----+----+
    |      |  ID |  ID  |  ID   | etc |    |    |    |    |
    +------+-----+------+-------+-----+----+----+----+----+
    """
    type = 'NSM1'

    @classmethod
    def _init_from_empty(cls):
        sid = 1
        nsm_type = 'PSHELL'
        pid_eid = 42
        value = 1.
        return NSM1(sid, nsm_type, value, pid_eid, comment='')

    def __init__(self, sid, nsm_type, value, pid_eid, comment=''):
        """See ``NSM1x``"""
        assert isinstance(value, float), f'NSM1; value={value!r} and must be a float'
        NSM1x.__init__(self, sid, nsm_type, value, pid_eid, comment='')


class NSM(NSMx):
    """
    Defines a set of non structural mass.

    +-----+-----+------+----+-------+----+-------+----+-------+
    |  1  |  2  |  3   |  4 |   5   | 6  |   7   | 8  |   9   |
    +=====+=====+======+====+=======+====+=======+====+=======+
    | NSM | SID | TYPE | ID | VALUE | ID | VALUE | ID | VALUE |
    +-----+-----+------+----+-------+----+-------+----+-------+
    """
    type = 'NSM'
    _properties = ['ids']

    @classmethod
    def _init_from_empty(cls):
        sid = 1
        nsm_type = 'PSHELL'
        pid_eid = 42
        value = 1.
        return NSM(sid, nsm_type, pid_eid, value, comment='')

    def __init__(self, sid, nsm_type, pid_eid, value, comment=''):
        """See ``NSMx``"""
        NSMx.__init__(self, sid, nsm_type, pid_eid, value, comment='')


class NSML(NSMx):
    """
    Defines a set of lumped non structural mass.

    +------+-----+------+----+-------+----+-------+----+-------+
    |  1   |  2  |  3   |  4 |   5   | 6  |   7   | 8  |   9   |
    +======+=====+======+====+=======+====+=======+====+=======+
    | NSML | SID | TYPE | ID | VALUE | ID | VALUE | ID | VALUE |
    +------+-----+------+----+-------+----+-------+----+-------+
    """
    type = 'NSML'
    _properties = ['ids']

    @classmethod
    def _init_from_empty(cls):
        sid = 1
        nsm_type = 'PSHELL'
        pid_eid = 42
        value = 1.
        return NSML(sid, nsm_type, pid_eid, value, comment='')

    def __init__(self, sid, nsm_type, pid_eid, value, comment=''):
        """
        Creates an NSML card, which defines lumped non-structural mass

        Parameters
        ----------
        sid : int
            Case control NSM id
        nsm_type : str
            Type of card the NSM is applied to
            valid_properties = {
                PSHELL, PCOMP, PBAR, PBARL, PBEAM, PBEAML, PBCOMP,
                PROD, CONROD, PBEND, PSHEAR, PTUBE, PCONEAX, PRAC2D,
                ELEMENT
            }
        pid_eid : int
            property id or element id depending on nsm_type
        value : float
            the non-structural pass per unit length/area
        comment : str; default=''
            a comment for the card
        """
        NSMx.__init__(self, sid, nsm_type, pid_eid, value, comment=comment)


def get_id_id0(ids: list[int] | list[str],
               mydict: dict[int, BaseCard]) -> tuple[list[int], int]:
    if ids == ['ALL']:
        all_ids = list(mydict)
        id0 = all_ids[0]
    else:
        id0 = ids[0]
        all_ids = ids
        assert isinstance(all_ids, list), all_ids
        assert isinstance(all_ids[0], integer_types), all_ids
    assert isinstance(id0, integer_types), (self.ids, all_ids)
    return all_ids, id0


class NSML1(NSM1x):
    """
    Defines lumped non structural mass entries by VALUE,ID list.

    +-------+-------+---------+-------+-------+------+-------+------+------+
    |   1   |   2   |    3    |   4   |   5   |  6   |   7   |  8   |  9   |
    +=======+=======+=========+=======+=======+======+=======+======+======+
    | NSML1 |  SID  |   TYPE  | VALUE |  ID   |  ID  |   ID  |  ID  |  ID  |
    +-------+-------+---------+-------+-------+------+-------+------+------+
    |       |   ID  |    ID   |  ID   |  etc  |      |       |      |      |
    +-------+-------+---------+-------+-------+------+-------+------+------+
    | NSML1 |   3   | ELEMENT |  .044 |  1240 | 1500 |       |      |      |
    +-------+-------+---------+-------+-------+------+-------+------+------+
    | NSML1 |  SID  |  TYPE   | VALUE |   ID  | THRU |   ID  |  ID  | THRU |
    +-------+-------+---------+-------+-------+------+-------+------+------+
    |       |   ID  |    ID   | THRU  |   ID  |  ID  |  THRU |  ID  |  ID  |
    +-------+-------+---------+-------+-------+------+-------+------+------+
    |       |  THRU |    ID   |  ...  |       |      |       |      |      |
    +-------+-------+---------+-------+-------+------+-------+------+------+
    | NSML1 |   15  |  PSHELL |  .067 |  1240 | THRU |  1760 |      |      |
    +-------+-------+---------+-------+-------+------+-------+------+------+
    |       |  2567 |   THRU  |  2568 | 35689 | THRU | 40998 |      |      |
    +-------+-------+---------+-------+-------+------+-------+------+------+
    |       |   76  |   THRU  |  300  |       |      |       |      |      |
    +-------+-------+---------+-------+-------+------+-------+------+------+
    | NSML1 |  SID  |   TYPE  | VALUE |   ID  | THRU |   ID  |  BY  |   N  |
    +-------+-------+---------+-------+-------+------+-------+------+------+
    |       |   ID  |   THRU  |   ID  |   BY  |   N  |       |      |      |
    +-------+-------+---------+-------+-------+------+-------+------+------+
    | NSML1 |   3   |  PSHELL |  .067 |  1240 | THRU |  1760 | 1763 | 1764 |
    +-------+-------+---------+-------+-------+------+-------+------+------+
    |       |  2567 |   THRU  |  2568 | 35689 |  TO  | 40999 |  BY  |   2  |
    +-------+-------+---------+-------+-------+------+-------+------+------+
    |       | 76666 |  76668  | 79834 |       |      |       |      |      |
    +-------+-------+---------+-------+-------+------+-------+------+------+
    """
    type = 'NSML1'

    @classmethod
    def _init_from_empty(cls):
        sid = 1
        nsm_type = 'PSHELL'
        ids = [1, 2]
        value = 1.
        return NSML1(sid, nsm_type, value, ids, comment='')

    def __init__(self, sid, nsm_type, value, ids, comment=''):
        """
        Creates an NSML card, which defines lumped non-structural mass

        Parameters
        ----------
        sid : int
            Case control NSM id
        nsm_type : str
            Type of card the NSM is applied to
            valid_properties = {
                PSHELL, PCOMP, PBAR, PBARL, PBEAM, PBEAML, PBCOMP,
                PROD, CONROD, PBEND, PSHEAR, PTUBE, PCONEAX, PRAC2D,
                ELEMENT
            }
        value : float
            the non-structural pass per unit length/area
        ids : list[int]
            property ids or element ids depending on nsm_type
        comment : str; default=''
            a comment for the card
        """
        assert isinstance(value, float), f'NSML1; value={value!r} and must be a float'
        NSM1x.__init__(self, sid, nsm_type, value, ids, comment=comment)

    def get_eid_mass_cg_by_element(self, model: BDF) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """lumped mass to be distributed

        - Area element calculation (for example, CQUAD4):
          Mass for a particular element = (Element area) x (VALUE / ΣElement area for all IDs)
        - Line element calculation (for example, CBEAM):
          Mass for a particular element = (Element length) x (VALUE / ΣElement length for all IDs)

        """
        value = self.value
        divide_by_sum = True

        use_length = False
        use_area = False

        line_elements = {
            'CBAR', 'CBEAM', 'CROD', 'CTUBE', 'CBUSH', 'CONROD'}
        shell_elements = {
            'CTRIA3', 'CTRIA6', 'CTRIAR',
            'CQUAD4', 'CQUADR', 'CQUAD8', 'CQUAD'}

        line_properties = {
            'PROD', 'PTUBE',
            'PBAR', 'PBARL', 'PBEAM', 'PBEAML',
            'PBUSH',
        }

        if self.nsm_type in {'ELEMENT', 'CONROD'}:
            all_ids, id0 = get_id_id0(self.ids, model.elements)

            if self.nsm_type == 'CONROD':
                all_ids = [eid for eid in all_ids
                           if model.elements[eid].type == 'CONROD']
                id0 = all_ids[0]

            ncards = len(all_ids)
            eids = np.zeros(ncards, dtype='int32')
            cg = np.zeros((ncards, 3), dtype='float64')

            elem0 = model.elements[id0]
            etype = elem0.type
            if etype in line_elements:
                use_length = True
                length = np.zeros(ncards, dtype='float64')
                for i, eid in enumerate(all_ids):
                    eids[i] = eid
                    elem = model.elements[eid]
                    assert elem.type in line_elements, elem.type
                    length[i] = elem.Length()
                    cg[i, :] = elem.Centroid()
            elif etype in shell_elements:
                use_area = True
                area = np.zeros(ncards, dtype='float64')
                for i, eid in enumerate(all_ids):
                    eids[i] = eid
                    elem = model.elements[eid]
                    assert elem.type in shell_elements, elem.type
                    area[i] = elem.Area()
                    cg[i, :] = elem.Centroid()
                #raise NotImplementedError((self.nsm_type, 'shell'))
            else:
                raise NotImplementedError((self.nsm_type, 'shell', etype))

        elif self.nsm_type in line_properties:
            pline_properties = {self.nsm_type}
            if self.nsm_type in {'PBAR', 'PBARL'}:
                pline_properties = {'PBAR', 'PBARL'}
                allowed_etype = 'CBAR'
            elif self.nsm_type in {'PBEAM', 'PBEAML'}:
                pline_properties = {'PBEAM', 'PBEAML'}
                allowed_etype = 'CBEAM'
            elif self.nsm_type == 'PROD':
                allowed_etype = 'CROD'
            elif self.nsm_type == 'PTUBE':
                allowed_etype = 'CTUBE'
            elif self.nsm_type == 'PBUSH':
                allowed_etype = 'CBUSH'
            else:
                raise NotImplementedError((self.nsm_type, 'line'))

            if self.ids == ['ALL']:
                all_ids = [pid for pid, prop in model.properties.items()
                           if prop.type in pline_properties]
                id0 = all_ids[0]
            else:
                all_ids, id0 = get_id_id0(self.ids, model.properties)

            use_length = True
            prop = model.properties[id0]
            ptype = prop.type
            assert ptype in pline_properties, f'allowed={pline_properties}\n{prop}'
            elems = [(eid, elem) for eid, elem in sorted(model.elements.items())
                     if elem.type == allowed_etype and elem.pid in self.ids]

            ncards = len(elems)
            eids = np.zeros(ncards, dtype='int32')
            length = np.zeros(ncards, dtype='float64')
            cg = np.zeros((ncards, 3), dtype='float64')
            for i, (eid, elem) in enumerate(elems):
                eids[i] = eid
                elem = model.elements[eid]
                assert prop.type in pline_properties, (elem, prop)
                length[i] = elem.Length()
                cg[i, :] = elem.Centroid()

        elif self.nsm_type in {'PSHELL', 'PCOMP', 'PCOMPG'}:
            use_area = True

            if self.nsm_type == 'PSHELL':
                shell_properties = {'PSHELL'}
            elif self.nsm_type in {'PCOMP', 'PCOMPG'}:
                shell_properties = {'PCOMP', 'PCOMPG'}
            else:
                raise NotImplementedError((self.nsm_type, 'shell'))

            if self.ids == ['ALL']:
                all_ids = [pid for pid, prop in model.properties.items()
                           if prop.type in shell_properties]
            else:
                all_ids = self.ids
            id0 = all_ids[0]

            prop = model.properties[id0]
            ptype = prop.type
            assert ptype in shell_properties, prop
            elems = [(eid, elem) for eid, elem in sorted(model.elements.items())
                     if elem.type in shell_elements and elem.pid in self.ids]

            ncards = len(elems)
            eids = np.zeros(ncards, dtype='int32')
            area = np.zeros(ncards, dtype='float64')
            cg = np.zeros((ncards, 3), dtype='float64')
            for i, (eid, elem) in enumerate(elems):
                eids[i] = eid
                elem = model.elements[eid]
                assert prop.type in shell_properties, (elem, prop)
                area[i] = elem.Area()
                cg[i, :] = elem.Centroid()
        else:  # pragma: no cover
            raise NotImplementedError(self.nsm_type)

        # here's where we do the distribution
        if use_area:
            mass = value * area / area.sum()
        elif use_length:
            mass = value * length / length.sum()
        else:
            raise RuntimeError((use_area, use_length))
        return eids, mass, cg


class NSMADD(BaseCard):
    """
    Defines non structural mass as the sum of the sets listed.

    +--------+----+----+-----+
    |    1   | 2  |  3 |  4  |
    +========+====+====+=====+
    | NSMADD | 2  |  1 |  3  |
    +--------+----+----+-----+
    """
    type = 'NSMADD'
    _properties = ['ids', 'nsm_ids', 'conid']

    @classmethod
    def _init_from_empty(cls):
        sid = 1
        sets = [1, 2]
        return NSMADD(sid, sets, comment='')

    def __init__(self, sid, sets, comment=''):
        """
        Creates an NSMADD card, which sum NSM sets

        Parameters
        ----------
        sid : int
            the NSM Case Control value
        sets : list[int]
            the NSM, NSM1, NSML, NSML1 values
        comment : str; default=''
            a comment for the card
        """
        BaseCard.__init__(self)
        if comment:
            self.comment = comment
        self.sid = sid
        self.sets = expand_thru(sets)
        self.sets.sort()
        self.sets_ref = None

    def validate(self):
        usets = np.unique(self.sets)
        if not len(self.sets) == len(usets):
            warnings.warn(f'NSMADD sid={self.sid} has duplicate set ids\n{str(self.sets)}')
        self.sets_ref = None

    @classmethod
    def add_card(cls, card: BDFCard, comment: str=''):
        """
        Adds an NSMADD card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        sid = integer(card, 1, 'sid')
        sets = card.fields(2)
        return NSMADD(sid, sets, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        """
        Adds an NSMADD card from the OP2

        Parameters
        ----------
        data : list[varies]
            a list of fields defined in OP2 format
        comment : str; default=''
            a comment for the card
        """
        sid = data[0]
        sets = list(data[1:-1])
        return NSMADD(sid, sets, comment=comment)

    @property
    def nsm_ids(self):
        """gets the nonstructural-mass ids"""
        if self.sets_ref is None:
            return self.sets

        nsm_ids = []
        for nsm in self.sets_ref:
            if isinstance(nsm, integer_types):
                nsm_ids.append(nsm)
            elif isinstance(nsm, list):
                nsm_ids.append(nsm[0].sid)
            else:
                raise TypeError('type=%s; nsm=\n%s' % (type(nsm), nsm))
        return nsm_ids

    @property
    def ids(self):
        return self.nsm_ids

    @property
    def conid(self) -> int:
        return self.sid

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ', which is required by NSMADD=%s' % self.sid
        nsms, failed_nsms = self._nsm_failed(model)
        if failed_nsms:
            nsm_list = list(model.nsms)
            raise KeyError(f'failed_nsms={failed_nsms}; nsms={nsm_list}{msg}')
        self.sets_ref = nsms

    def _nsm_failed(self, model: BDF) -> tuple[list[NSM], list[int]]:
        msg = ', which is required by NSMADD=%s' % self.sid
        nsms = []
        failed_nsms = []
        for nsm_id in self.sets:
            try:
                nsm = model.NSM(nsm_id, msg=msg)
                nsms.append(nsm)
            except KeyError:
                failed_nsms.append(nsm_id)
        return nsms, failed_nsms

    def safe_cross_reference(self, model: BDF, debug: bool=True):
        msg = ', which is required by NSMADD=%s' % self.sid
        nsms, failed_nsms = self._nsm_failed(model)
        if failed_nsms:
            nsm_list = list(model.nsms)
            model.log.warning(f'failed_nsms={failed_nsms}; nsms={nsm_list}{msg}')
        self.sets_ref = nsms

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.sets = self.nsm_ids
        self.sets_ref = None

    def raw_fields(self):
        fields = ['NSMADD', self.sid] + self.nsm_ids
        return fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.raw_fields()
        return self.comment + print_card_8(card)

    def write_card_16(self, is_double=False):
        card = self.raw_fields()
        return self.comment + print_card_16(card)


class PMASS(Property):
    """
    Scalar Mass Property
    Specifies the mass value of a scalar mass element (CMASS1 or CMASS3 entries).

    +-------+------+------+------+------+------+----+------+----+
    |   1   |   2  |   3  |   4  |   5  |   6  |  7 |   8  |  9 |
    +=======+======+======+======+======+======+====+======+====+
    | PMASS | PID1 |  M1  | PID2 |  M2  | PID3 | M3 | PID4 | M4 |
    +-------+------+------+------+------+------+----+------+----+
    | PMASS |   7  | 4.29 |   6  | 13.2 |      |    |      |    |
    +-------+------+------+------+------+------+----+------+----+
    """
    type = 'PMASS'
    _field_map = {
        1: 'pid', 2:'mass',
    }
    pname_fid_map = {
        # 1-based
        3 : 'mass',
    }
    _properties = ['node_ids']

    @classmethod
    def _init_from_empty(cls):
        pid = 1
        mass = 1.
        return PMASS(pid, mass, comment='')

    def __init__(self, pid: int, mass: float, comment: str=''):
        """
        Creates an PMASS card, which defines a mass applied to a single DOF

        Parameters
        ----------
        pid : int
            Property id used by a CMASS1/CMASS3 card
        mass : float
            the mass to apply
        comment : str; default=''
            a comment for the card
        """
        Property.__init__(self)
        if comment:
            self.comment = comment
        self.pid = pid
        self.mass = mass

    @classmethod
    def add_card(cls, card: BDFCard, icard: int=0, comment: str=''):
        """
        Adds a PMASS card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        icard : int; default=0
            the index of the card that's being parsed
        comment : str; default=''
            a comment for the card
        """
        icard *= 2
        #: Property ID
        pid = integer(card, 1 + icard, 'pid')
        mass = double_or_blank(card, 2 + icard, 'mass', default=0.)
        return PMASS(pid, mass, comment=comment)

    @classmethod
    def add_card_lax(cls, card: BDFCard, icard: int=0, comment: str=''):
        """
        Adds a PMASS card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        icard : int; default=0
            the index of the card that's being parsed
        comment : str; default=''
            a comment for the card
        """
        icard *= 2
        #: Property ID
        pid = integer(card, 1 + icard, 'pid')
        mass = force_double_or_blank(card, 2 + icard, 'mass', default=0.)
        return PMASS(pid, mass, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        """
        Adds a PMASS card from the OP2

        Parameters
        ----------
        data : list[varies]
            a list of fields defined in OP2 format
        comment : str; default=''
            a comment for the card
        """
        pid = data[0]
        mass = data[1]
        return PMASS(pid, mass, comment=comment)

    def cross_reference(self, model: BDF) -> None:
        pass

    def safe_cross_reference(self, model: BDF, xref_errors) -> None:
        pass

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        pass

    def _verify(self, xref):
        pid = self.Pid()
        mass = self.Mass()
        assert isinstance(pid, integer_types), 'pid=%r' % pid
        assert isinstance(mass, float), 'mass=%r' % mass

    def Mass(self) -> float:
        return self.mass

    def raw_fields(self):
        list_fields = ['PMASS', self.pid, self.mass]
        return list_fields

    def repr_fields(self):
        return self.raw_fields()

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)

NSMs = NSM | NSM1 | NSML | NSML1
