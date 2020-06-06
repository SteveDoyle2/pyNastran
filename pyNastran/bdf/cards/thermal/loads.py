# pylint: disable=R0902,R0904,R0914,C0111
from __future__ import annotations
from typing import TYPE_CHECKING
from pyNastran.utils.numpy_utils import integer_types
from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.field_writer_16 import print_card_16
from pyNastran.bdf.field_writer_double import print_card_double
from pyNastran.bdf.cards.utils import wipe_empty_fields
from pyNastran.bdf.cards.thermal.thermal import ThermalCard
from pyNastran.bdf.field_writer_8 import set_blank_if_default
from pyNastran.bdf.cards.base_card import expand_thru, expand_thru_by, BaseCard
from pyNastran.bdf.cards.collpase_card import collapse_thru_by
from pyNastran.bdf.bdf_interface.assign_type import (
    integer, integer_or_blank, double, double_or_blank, integer_or_string,
    integer_double_or_blank, string, fields)
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf import BDF


#class ThermalLoadDefault(ThermalCard):
    #def __init__(self, card, data):
        #ThermalCard.__init__(self)


class ThermalLoad(ThermalCard):
    def __init__(self):
        ThermalCard.__init__(self)

class QVOL(ThermalLoad):
    """
    Defines a rate of volumetric heat addition in a conduction element.

    +------+------+------+---------+------+------+------+------+------+
    |  1   |   2  |   3  |    4    |   5  |   6  |   7  |   8  |   9  |
    +======+======+======+=========+======+======+======+======+======+
    | QVOL | SID  | QVOL | CNTRLND | EID1 | EID2 | EID3 | EID4 | EID5 |
    +------+------+------+---------+------+------+------+------+------+
    |      | EID6 | etc. |         |      |      |      |      |      |
    +------+------+------+---------+------+------+------+------+------+

    """
    type = 'QVOL'
    _properties = ['element_ids']

    @classmethod
    def _init_from_empty(cls):
        sid = 1
        #flag = 'LINE'
        qvol = 1.0
        control_point = 10
        elements = [1, 2]
        return QVOL(sid, qvol, control_point, elements, comment='')

    def __init__(self, sid, qvol, control_point, elements, comment=''):
        ThermalLoad.__init__(self)
        if comment:
            self.comment = comment
        #: Load set identification number. (Integer > 0)
        self.sid = sid
        self.qvol = qvol
        self.control_point = control_point
        if isinstance(elements, integer_types):
            elements = [elements]
        self.elements = elements
        self.elements_ref = None

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a QVOL card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        sid = integer(card, 1, 'sid')
        qvol = double(card, 2, 'qvol')
        control_point = integer_or_blank(card, 3, 'control_id', 0)

        i = 1
        eids = []
        for ifield in range(4, len(card)):
            eid = integer_or_string(card, ifield, 'eid_%i' % i)
            eids.append(eid)
            i += 1
        elements = expand_thru_by(eids)
        return QVOL(sid, qvol, control_point, elements, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        """
        Adds a QVOL card from the OP2

        Parameters
        ----------
        data : List[varies]
            a list of fields defined in OP2 format
        comment : str; default=''
            a comment for the card

        """
        sid, qvol, control_point, eid = data
        return QVOL(sid, qvol, control_point, eid, comment=comment)

    def get_loads(self):
        return [self]

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        msg = ', which is required by QVOL sid=%s' % self.sid
        self.elements_ref = model.Elements(self.elements, msg=msg)

    def safe_cross_reference(self, model, xref_errors):
        try:
            return self.cross_reference(model)
        except KeyError:
            model.log.warning('failed cross-referencing\n%s' % str(self))
            raise

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.elements = self.element_ids
        self.elements_ref = None

    @property
    def element_ids(self):
        if self.elements_ref is None:
            return self.elements
        eids = []
        for eid_ref in self.elements_ref:
            eids.append(eid_ref.eid)
        return eids

    def Eids(self):
        return self.element_ids

    def raw_fields(self):
        list_fields = ['QVOL', self.sid, self.qvol, self.control_point] + self.element_ids
        return list_fields

    def repr_fields(self):
        eids = collapse_thru_by(self.element_ids)
        list_fields = ['QVOL', self.sid, self.qvol, self.control_point] + eids
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        if is_double:
            return self.comment + print_card_double(card)
        return self.comment + print_card_16(card)


class QVECT(ThermalLoad):
    """
    Thermal Vector Flux Load

    Defines thermal vector flux from a distant source into a face of one
    or more CHBDYi boundary condition surface elements.

    +-------+------+------+-------+-----+---------+---------+---------+---------+
    |   1   |   2  |   3  |   4   |  5  |     6   |    7    |    8    |    9    |
    +=======+======+======+=======+=====+=========+=========+=========+=========+
    | QVECT | SID  |  Q0  | TSOUR | CE  | E1/TID1 | E2/TID2 | E3/TID3 | CNTRLND |
    +-------+------+------+-------+-----+---------+---------+---------+---------+
    |       | EID1 | EID2 |  etc. |     |         |         |         |         |
    +-------+------+------+-------+-----+---------+---------+---------+---------+

    """
    type = 'QVECT'
    _properties = ['element_ids']

    @classmethod
    def _init_from_empty(cls):
        sid = 1
        q0 = 1.0
        eids = [1, 2]
        return QVECT(sid, q0, eids, t_source=None, ce=0,
                     vector_tableds=None, control_id=0, comment='')

    def __init__(self, sid, q0, eids, t_source=None,
                 ce=0, vector_tableds=None, control_id=0, comment=''):
        """
        Creates a QVECT card

        Parameters
        ----------
        sid : int
            Load set identification number. (Integer > 0)
        q0 : float; default=None
            Magnitude of thermal flux vector into face
        t_source : float; default=None
            Temperature of the radiant source
        ce : int; default=0
            Coordinate system identification number for thermal vector flux
        vector_tableds : List[int/float, int/float, int/float]
            vector : float; default=None
                directional cosines in coordinate system CE) of
                the thermal vector flux
                None : [0.0, 0.0, 0.0]
            tabled : int
                TABLEDi entry identification numbers defining the
                components as a function of time
        control_id : int; default=0
            Control point
        eids : List[int] or THRU
            Element identification number of a CHBDYE, CHBDYG, or
            CHBDYP entry
        comment : str; default=''
            a comment for the card

        """
        ThermalLoad.__init__(self)
        if comment:
            self.comment = comment
        #: Load set identification number. (Integer > 0)
        self.sid = sid
        self.q0 = q0
        self.t_source = t_source
        self.ce = ce
        self.control_id = control_id

        if vector_tableds is None:
            self.vector_tableds = [0., 0., 0.]
        else:
            self.vector_tableds = vector_tableds
        self.eids = eids
        self.eids_ref = None

    def validate(self):
        assert isinstance(self.eids, list), 'type(eids)=%s' % type(self.eids)

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a QVECT card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        sid = integer(card, 1, 'sid')
        q0 = double(card, 2, 'q0')
        t_source = double_or_blank(card, 3, 't_source')
        ce = integer_or_blank(card, 4, 'ce', 0)
        vector_tableds = [
            integer_double_or_blank(card, 5, 'e1_tabled1', 0.0),
            integer_double_or_blank(card, 6, 'e2_tabled2', 0.0),
            integer_double_or_blank(card, 7, 'e3_tabled3', 0.0),
        ]
        control_id = integer_or_blank(card, 8, 'control_id', 0)

        i = 1
        eids = []
        for ifield in range(9, len(card)):
            eid = integer_or_string(card, ifield, 'eid_%i' % i)
            eids.append(eid)
            assert eid != 0, card
            i += 1
        elements = expand_thru_by(eids)
        return QVECT(sid, q0, elements, t_source=t_source,
                     ce=ce, vector_tableds=vector_tableds, control_id=control_id,
                     comment=comment)

    def get_loads(self):
        return [self]

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        msg = ', which is required by QVECT sid=%s' % self.sid
        self.eids_ref = model.Elements(self.eids, msg=msg)

    def safe_cross_reference(self, model, xref_errors):
        try:
            return self.cross_reference(model)
        except KeyError:
            model.log.warning('failed cross-referencing\n%s' % str(self))

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.eids = self.element_ids
        self.eids_ref = None

    @property
    def element_ids(self):
        if self.eids_ref is None:
            return self.eids
        eids = []
        for eid_ref in self.eids_ref:
            eids.append(eid_ref.eid)
        return eids

    def Eids(self):
        return self.element_ids

    def raw_fields(self):
        list_fields = [
            'QVECT', self.sid, self.q0, self.t_source, self.ce
            ] + self.vector_tableds + [self.control_id] + self.element_ids
        return list_fields

    def repr_fields(self):
        eids = collapse_thru_by(self.element_ids)
        list_fields = [
            'QVECT', self.sid, self.q0, self.t_source, self.ce
            ] + self.vector_tableds + [self.control_id] + eids
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        if is_double:
            return self.comment + print_card_double(card)
        return self.comment + print_card_16(card)


class QBDY1(ThermalLoad):
    """
    Defines a uniform heat flux into CHBDYj elements.

    """
    type = 'QBDY1'
    _properties = ['element_ids']

    @classmethod
    def _init_from_empty(cls):
        sid = 1
        qflux = 1.0
        eids = [1, 2]
        return QBDY1(sid, qflux, eids, comment='')

    def __init__(self, sid, qflux, eids, comment=''):
        ThermalLoad.__init__(self)
        if comment:
            self.comment = comment
        #: Load set identification number. (Integer > 0)
        self.sid = sid
        #: Heat flux into element (FLOAT)
        self.qflux = qflux
        #: CHBDYj element identification numbers (Integer)
        #: .. todo:: use expand_thru_by ???
        assert len(eids) > 0
        self.eids = expand_thru(eids)
        self.eids_ref = None

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a QBDY1 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        sid = integer(card, 1, 'sid')
        qflux = double(card, 2, 'qflux')
        eids = []
        j = 1
        for i in range(3, len(card)):
            eid = integer_or_string(card, i, 'eid%i' % j)
            eids.append(eid)
            j += 1
        return QBDY1(sid, qflux, eids, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        """
        Adds a QBDY1 card from the OP2

        Parameters
        ----------
        data : List[varies]
            a list of fields defined in OP2 format
        comment : str; default=''
            a comment for the card

        """
        sid = data[0]
        qflux = data[1]
        eids = data[2:]
        return QBDY1(sid, qflux, eids, comment=comment)

    def get_loads(self):
        return [self]

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        msg = ', which is required by QBDY1 sid=%s' % self.sid
        self.eids_ref = model.Elements(self.eids, msg=msg)

    def safe_cross_reference(self, model, xref_errors):
        try:
            return self.cross_reference(model)
        except KeyError:
            model.log.warning('failed cross-referencing\n%s' % str(self))

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.eids = self.element_ids
        self.eids_ref = None

    def nQFluxTerms(self):
        return len(self.qflux)

    @property
    def element_ids(self):
        if self.eids_ref is None:
            return self.eids
        eids = []
        for eid_ref in self.eids_ref:
            eids.append(eid_ref.eid)
        return eids

    def Eids(self):
        return self.element_ids

    def raw_fields(self):
        list_fields = ['QBDY1', self.sid, self.qflux] + self.element_ids
        return list_fields

    def repr_fields(self):
        eids = collapse_thru_by(self.element_ids)
        list_fields = ['QBDY1', self.sid, self.qflux] + eids
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        if is_double:
            return self.comment + print_card_double(card)
        return self.comment + print_card_16(card)


class QBDY2(ThermalLoad):  # not tested
    """
    Defines a uniform heat flux load for a boundary surface.

    """
    type = 'QBDY2'

    @classmethod
    def _init_from_empty(cls):
        sid = 1
        eids = [1]
        qflux = [1.0]
        return QBDY2(sid, qflux, eids, comment='')

    def __init__(self, sid, eid, qfluxs, comment=''):
        ThermalLoad.__init__(self)
        if comment:
            self.comment = comment

        #: Load set identification number. (Integer > 0)
        self.sid = sid
        #: Identification number of an CHBDYj element. (Integer > 0)
        self.eid = eid
        #: Heat flux at the i-th grid point on the referenced CHBDYj
        #: element. (Real or blank)
        if isinstance(qfluxs, float):
            qfluxs = [qfluxs]
        self.qfluxs = qfluxs
        self.eid_ref = None

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a QBDY2 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        sid = integer(card, 1, 'sid')
        eid = integer(card, 2, 'eid')

        qfluxs = []
        j = 1
        for i in range(3, len(card)):
            q = double_or_blank(card, i, 'qFlux%i' % j)
            qfluxs.append(q)
            j += 1

        assert len(qfluxs) > 0, qfluxs
        qfluxs = wipe_empty_fields(qfluxs)
        return QBDY2(sid, eid, qfluxs, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        """
        Adds a QBDY2 card from the OP2

        Parameters
        ----------
        data : List[varies]
            a list of fields defined in OP2 format
        comment : str; default=''
            a comment for the card

        """
        sid = data[0]
        eid = data[1]
        qfluxs = [data[2]]
        return QBDY2(sid, eid, qfluxs, comment=comment)

    def get_loads(self):
        return [self]

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        msg = ', which is required by QBDY2 sid=%s' % self.sid
        self.eid_ref = model.Element(self.eid, msg=msg)

    def safe_cross_reference(self, model, xref_errors):
        try:
            return self.cross_reference(model)
        except KeyError:
            model.log.warning('failed cross-referencing\n%s' % str(self))

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.eid = self.Eid()
        self.eid_ref = None

    def Eid(self):
        if self.eid_ref is not None:
            return self.eid_ref.eid
        return self.eid

    def nQFluxTerms(self):
        return len(self.qfluxs)

    def raw_fields(self):
        list_fields = ['QBDY2', self.sid, self.Eid()] + self.qfluxs
        return list_fields

    def repr_fields(self):
        return self.raw_fields()

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        if is_double:
            return self.comment + print_card_double(card)
        return self.comment + print_card_16(card)


class QBDY3(ThermalLoad):
    """
    Defines a uniform heat flux load for a boundary surface.

    """
    type = 'QBDY3'
    _properties = ['element_ids']

    @classmethod
    def _init_from_empty(cls):
        sid = 1
        q0 = 1.0
        cntrlnd = 10
        eids = [1, 2]
        return QBDY3(sid, q0, cntrlnd, eids, comment='')

    def __init__(self, sid, q0, cntrlnd, eids, comment=''):
        """
        Creates a QBDY3 card

        Parameters
        ----------
        sid : int
            Load set identification number. (Integer > 0)
        q0 : float; default=None
            Magnitude of thermal flux vector into face
        control_id : int; default=0
            Control point
        eids : List[int] or THRU
            Element identification number of a CHBDYE, CHBDYG, or
            CHBDYP entry
        comment : str; default=''
            a comment for the card

        """
        ThermalLoad.__init__(self)
        if comment:
            self.comment = comment

        #: Load set identification number. (Integer > 0)
        self.sid = sid

        #: Heat flux into element
        self.q0 = q0

        #: Control point for thermal flux load. (Integer > 0; Default = 0)
        self.cntrlnd = cntrlnd

        #: CHBDYj element identification numbers
        self.eids = expand_thru_by(eids)
        self.eids_ref = None

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a QBDY3 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        sid = integer(card, 1, 'sid')
        q0 = double(card, 2, 'q0')
        cntrlnd = integer_or_blank(card, 3, 'cntrlnd', 0)

        nfields = card.nfields
        eids = fields(integer_or_string, card, 'eid', i=4, j=nfields)
        return QBDY3(sid, q0, cntrlnd, eids, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        """
        Adds a QBDY3 card from the OP2

        Parameters
        ----------
        data : List[varies]
            a list of fields defined in OP2 format
        comment : str; default=''
            a comment for the card

        """
        sid = data[0]
        q0 = data[1]
        cntrlnd = data[2]
        eids = list(data[3:])
        return QBDY3(sid, q0, cntrlnd, eids, comment=comment)

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        msg = ', which is required by QBDY3 sid=%s' % self.sid
        eids = []
        for eid in self.eids:
            eids.append(model.Element(eid, msg=msg))
        self.eids_ref = eids

    def safe_cross_reference(self, model, xref_errors):
        try:
            return self.cross_reference(model)
        except KeyError:
            model.log.warning('failed cross-referencing\n%s' % str(self))

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.eids = self.element_ids
        self.eids_ref = None

    @property
    def element_ids(self):
        if self.eids_ref is None:
            return self.eids
        eids = []
        for eid_ref in self.eids_ref:
            eids.append(eid_ref.eid)
        return eids

    def Eids(self):
        return self.element_ids

    def raw_fields(self):
        eids = self.element_ids
        eids.sort()
        list_fields = (['QBDY3', self.sid, self.q0, self.cntrlnd] +
                       collapse_thru_by(eids))
        return list_fields

    def repr_fields(self):
        cntrlnd = set_blank_if_default(self.cntrlnd, 0)
        eids = self.element_ids
        eids.sort()
        list_fields = ['QBDY3', self.sid, self.q0, cntrlnd] + collapse_thru_by(eids)
        return list_fields

    def get_loads(self):
        return [self]

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        if is_double:
            return self.comment + print_card_double(card)
        return self.comment + print_card_16(card)


class QHBDY(ThermalLoad):
    """
    Defines a uniform heat flux into a set of grid points.

    """
    type = 'QHBDY'
    _properties = ['flag_to_nnodes']
    flag_to_nnodes = {
        'POINT' : (1, 1),
        'LINE' : (2, 2),
        'REV' : (2, 2),
        'AREA3' : (3, 3),
        'AREA4' : (4, 4),
        'AREA6' : (4, 6), # 4-6
        'AREA8' : (5, 8), # 5-8
    }

    @classmethod
    def _init_from_empty(cls):
        sid = 1
        flag = 'LINE'
        q0 = 1.0
        grids = [1, 2]
        return QHBDY(sid, flag, q0, grids, af=None, comment='')

    def __init__(self, sid, flag, q0, grids, af=None, comment=''):
        """
        Creates a QHBDY card

        Parameters
        ----------
        sid : int
            load id
        flag : str
            valid_flags = {POINT, LINE, REV, AREA3, AREA4, AREA6, AREA8}
        q0 : float
            Magnitude of thermal flux into face. Q0 is positive for heat
            into the surface
        af : float; default=None
            Area factor depends on type
        grids : List[int]
            Grid point identification of connected grid points
        comment : str; default=''
            a comment for the card

        """
        ThermalLoad.__init__(self)
        if comment:
            self.comment = comment

        #: Load set identification number. (Integer > 0)
        self.sid = sid
        self.flag = flag

        #: Magnitude of thermal flux into face. Q0 is positive for heat
        #: into the surface. (Real)
        self.q0 = q0

        #: Area factor depends on type. (Real > 0.0 or blank)
        self.af = af

        #: Grid point identification of connected grid points.
        #: (Integer > 0 or blank)
        self.grids = grids

        assert flag in ['POINT', 'LINE', 'REV', 'AREA3', 'AREA4', 'AREA6', 'AREA8'], str(self)

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a QHBDY card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        sid = integer(card, 1, 'eid')
        flag = string(card, 2, 'flag')

        q0 = double(card, 3, 'q0')
        af = double_or_blank(card, 4, 'af')
        nnodes_required, nnodes_max = cls.flag_to_nnodes[flag]

        grids = []
        if nnodes_required == nnodes_max:
            for i in range(nnodes_required):
                grid = integer(card, 5 + i, 'grid%i' % (i + 1))
                grids.append(grid)
        else:
            int_node_count = 0
            for i in range(nnodes_max):
                grid = integer_or_blank(card, 5 + i, 'grid%i' % (i + 1))
                if grid is not None:
                    int_node_count += 1
                grids.append(grid)
            if int_node_count < nnodes_required:
                msg = 'int_node_count=%s nnodes_required=%s' % (int_node_count, nnodes_required)
                raise RuntimeError(msg)
        return QHBDY(sid, flag, q0, grids, af=af, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        """
        Adds a QHBDY card from the OP2

        Parameters
        ----------
        data : List[varies]
            a list of fields defined in OP2 format
        comment : str; default=''
            a comment for the card

        """
        sid = data[0]
        flag_int = data[1]
        q0 = data[2]
        af = data[3]
        grids = list(data[4:])

        if flag_int == 1:
            flag = 'POINT'
            nnodes = 1
        elif flag_int == 2:
            flag = 'LINE'
            nnodes = 2
        elif flag_int == 3:
            flag = 'REV'
            nnodes = 2
        elif flag_int == 5:
            flag = 'AREA4'
            nnodes = 4
        elif flag_int == 9:
            flag = 'AREA8'
            nnodes = 8
        else:  # pragma: no cover
            raise NotImplementedError(f'QHBDY sid={sid} flag_int={flag_int} data={data[2:]}')
        grids2 = grids[:nnodes]
        #print(sid, flag_int, flag, q0, af, grids, grids2)
        return QHBDY(sid, flag, q0, grids2, af=af, comment=comment)

    def get_loads(self):
        return [self]

    def cross_reference(self, model: BDF) -> None:
        pass

    def safe_cross_reference(self, model, xref_errors):
        try:
            return self.cross_reference(model)
        except KeyError:
            model.log.warning('failed cross-referencing\n%s' % str(self))

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        pass

    def raw_fields(self):
        list_fields = ['QHBDY', self.sid, self.flag, self.q0, self.af] + self.grids
        return list_fields

    def repr_fields(self):
        return self.raw_fields()

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        if is_double:
            return self.comment + print_card_double(card)
        return self.comment + print_card_16(card)


class TEMPRB(ThermalLoad):
    """

    +--------+-----+----+-------+----+-------+----+----+
    |    1   |  2  |  3 |   4   |  5 |   6   |  7 |  8 |
    +========+=====+====+=======+====+=======+====+====+
    | TEMPRB | SID | G1 |  T1   | G2 |  T2   | G3 | T3 |
    +--------+-----+----+-------+----+-------+----+----+

    +--------+------+------+------+------+------+------+------+------+
    |    1   |   2  |   3  |   4  |   5  |   6  |   7  |   8  |  9   |
    +========+======+======+======+======+======+======+======+======+
    | TEMPRB | SID  | EID1 |   TA |  TB  | TP1A | TP1B | TP2A | TP2B |
    +--------+------+------+------+------+------+------+------+------+
    |        | TCA  |  TDA |  TEA |  TFA |  TCB |  TDB |  TEB | TFB  |
    +--------+------+------+------+------+------+------+------+------+
    |        | EID2 | EID3 | EID4 | EID5 | EID6 | EID7 | etc. |      |
    +--------+------+------+------+------+------+------+------+------+
    | TEMPRB | 200  |   1  | 68.0 | 23.0 | 0.0  | 28.0 | 2.5  |      |
    +--------+------+------+------+------+------+------+------+------+
    |        | 68.0 | 91.0 | 45.0 | 48.0 | 80.0 | 20.0 |      |      |
    +--------+------+------+------+------+------+------+------+------+
    |        |   9  |  10  |      |      |      |      |      |      |
    +--------+------+------+------+------+------+------+------+------+

    """

    type = 'TEMPRB'

    @classmethod
    def _init_from_empty(cls):
        sid = 1
        eids = [1]
        ta = 1.0
        tb = 1.0
        tp1 = [1.0, 1.0]
        tp2 = [1.0, 1.0]
        tai = [1.0, 1.0, 1.0, 1.0]
        tbi = [1.0, 1.0, 1.0, 1.0]
        return TEMPRB(sid, ta, tb, tp1, tp2, tai, tbi, eids, comment='')

    #@classmethod
    #def export_to_hdf5(cls, h5_file, model, loads):
        #"""exports the loads in a vectorized way"""
        ##encoding = model._encoding
        ##comments = []
        #unused_sid = loads[0].sid
        ##h5_file.create_dataset('sid', data=sid)
        #for i, load in enumerate(loads):
            #sub_group = h5_file.create_group(str(i))
            ##comments.append(loads.comment)

            #node = []
            #temp = []
            #for nid, tempi in load.temperatures.items():
                #node.append(nid)
                #temp.append(tempi)
            #sub_group.create_dataset('node', data=node)
            #sub_group.create_dataset('temperature', data=temp)
        #h5_file.create_dataset('_comment', data=comments)

    def __init__(self, sid, ta, tb, tp1, tp2,
                 tai, tbi, eids, comment=''):
        """
        Creates a TEMPRB card

        Parameters
        ----------
        sid : int
            Load set identification number
        comment : str; default=''
            a comment for the card

        """
        ThermalLoad.__init__(self)
        if comment:
            self.comment = comment
        #: Load set identification number. (Integer > 0)
        self.sid = sid

        # temperatures at end A/B
        self.ta = ta
        self.tb = tb

        # gradients at A/B
        self.tp1 = tp1
        self.tp2 = tp2

        # temperatures at points C/D/E/F at A/B
        self.tai = tai
        self.tbi = tbi
        self.eids = eids

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a TEMPRB card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        sid = integer(card, 1, 'sid')

        #TEMPRB SID  EID1   TA   TB TP1A TP1B TP2A TP2B
               #TCA   TDA  TEA  TFA  TCB  TDB  TEB TFB
               #EID2 EID3 EID4 EID5 EID6 EID7 -etc.-
        #TEMPRB 200 1 68.0 23.0 0.0 28.0 2.5
               #68.0 91.0 45.0 48.0 80.0 20.0
               #9 10

        ##+========+======+======+======+======+======+======+======+======+
        ##| TEMPRB | SID  | EID1 |   TA |  TB  | TP1A | TP1B | TP2A | TP2B |
        ##|        | TCA  |  TDA |  TEA |  TFA |  TCB |  TDB |  TEB | TFB  |
        ##|        | EID2 | EID3 | EID4 | EID5 | EID6 | EID7 | etc. |      |
        ##+--------+------+------+------+------+------+------+------+------+

        eid1 = integer(card, 2, 'EID1')
        # TODO: the default for MSC is 0.0
        ta = double(card, 3, 'T(A)')
        tb = double(card, 4, 'T(B)')
        tp1a = double_or_blank(card, 5, 'TP1(A)')
        tp1b = double_or_blank(card, 6, 'TP1(B)')
        tp2a = double_or_blank(card, 7, 'TP2(A)')
        tp2b = double_or_blank(card, 8, 'TP2(B)')
        tp1 = [tp1a, tp1b]
        tp2 = [tp2a, tp2b]

        # TODO: the default for MSC is TA
        tca = double_or_blank(card, 9, 'TC(A)')
        tda = double_or_blank(card, 10, 'TD(A)')
        tea = double_or_blank(card, 11, 'TE(A)')
        tfa = double_or_blank(card, 12, 'TF(A)')

        # TODO: the default for MSC is TB
        tcb = double_or_blank(card, 13, 'TC(B)')
        tdb = double_or_blank(card, 14, 'TD(B)')
        teb = double_or_blank(card, 15, 'TE(B)')
        tfb = double_or_blank(card, 16, 'TF(B)')
        eids = [eid1]

        nfields = len(card)
        #assert nfields <= 16, 'len(card)=%i card=%s' % (len(card), card)
        for ifield in range(17, nfields):
            eid = integer(card, ifield, f'eid{ifield-15}')
            eids.append(eid)

        tai = [tca, tda, tea, tfa]
        tbi = [tcb, tdb, teb, tfb]
        return TEMPRB(sid, ta, tb, tp1, tp2,
                      tai, tbi, eids, comment=comment)

    def cross_reference(self, model: BDF) -> None:
        pass

    def safe_cross_reference(self, model, debug=True):
        pass

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        pass

    def raw_fields(self):
        """Writes the TEMPRB card"""
        eid = self.eids[0]
        list_fields = [
            'TEMPRB', self.sid, eid, self.ta, self.tb,
            self.tp1[0], self.tp1[1],
            self.tp2[0], self.tp2[1],
        ] + self.tai + self.tbi + self.eids[1:]
        return list_fields

    def repr_fields(self):
        """Writes the TEMPRB card"""
        return self.raw_fields()

    def get_loads(self):
        return [self]

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        if is_double:
            return self.comment + print_card_double(card)
        return self.comment + print_card_16(card)

class TEMP(ThermalLoad):
    """
    Defines temperature at grid points for determination of thermal loading,
    temperature-dependent material properties, or stress recovery.

    +------+-----+----+-------+----+-------+----+----+
    |   1  |  2  |  3 |   4   |  5 |   6   |  7 |  8 |
    +======+=====+====+=======+====+=======+====+====+
    | TEMP | SID | G1 |  T1   | G2 |  T2   | G3 | T3 |
    +------+-----+----+-------+----+-------+----+----+
    | TEMP |  3  | 94 | 316.2 | 49 | 219.8 |    |    |
    +------+-----+----+-------+----+-------+----+----+

    """
    type = 'TEMP'

    @classmethod
    def _init_from_empty(cls):
        sid = 1
        temperatures = {1 : 1.0}
        return TEMP(sid, temperatures, comment='')

    @classmethod
    def export_to_hdf5(cls, h5_file, model, loads):
        """exports the loads in a vectorized way"""
        #encoding = model._encoding
        #comments = []
        unused_sid = loads[0].sid
        #h5_file.create_dataset('sid', data=sid)
        for i, load in enumerate(loads):
            sub_group = h5_file.create_group(str(i))
            #comments.append(loads.comment)

            node = []
            temp = []
            for nid, tempi in load.temperatures.items():
                node.append(nid)
                temp.append(tempi)
            sub_group.create_dataset('node', data=node)
            sub_group.create_dataset('temperature', data=temp)
        #h5_file.create_dataset('_comment', data=comments)

    def __init__(self, sid, temperatures, comment=''):
        """
        Creates a TEMP card

        Parameters
        ----------
        sid : int
            Load set identification number
        temperatures : dict[nid] : temperature
            nid : int
                node id
            temperature : float
                the nodal temperature
        comment : str; default=''
            a comment for the card

        """
        ThermalLoad.__init__(self)
        if comment:
            self.comment = comment
        #: Load set identification number. (Integer > 0)
        self.sid = sid

        #: dictionary of temperatures where the key is the grid ID (Gi)
        #: and the value is the temperature (Ti)
        self.temperatures = temperatures

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a TEMP card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        sid = integer(card, 1, 'sid')

        nfields = len(card)
        assert nfields <= 8, 'len(card)=%i card=%s' % (len(card), card)

        ntemps = (nfields -2) // 2
        assert nfields % 2 == 0, card
        assert nfields // 2 > 1, card

        temperatures = {}
        for i in range(ntemps):
            n = i * 2 + 2
            gi = integer(card, n, 'g' + str(i))
            temperaturei = double(card, n + 1, 'T' + str(i))
            temperatures[gi] = temperaturei
        return TEMP(sid, temperatures, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        """
        Adds a TEMP card from the OP2

        Parameters
        ----------
        data : List[varies]
            a list of fields defined in OP2 format
        comment : str; default=''
            a comment for the card

        """
        sid = data[0]
        temperatures = {data[1]: data[2]}
        return TEMP(sid, temperatures, comment=comment)

    def add(self, temp_obj):
        assert self.sid == temp_obj.sid
        for (gid, temp) in temp_obj.temperatures.items():
            self.temperatures[gid] = temp

    def cross_reference(self, model: BDF) -> None:
        pass

    def safe_cross_reference(self, model, debug=True):
        pass

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        pass

    def raw_fields(self):
        """Writes the TEMP card"""
        list_fields = ['TEMP', self.sid]
        ntemps = len(self.temperatures) - 1
        for i, (gid, temp) in enumerate(sorted(self.temperatures.items())):
            list_fields += [gid, temp]
            if i % 3 == 2 and ntemps > i:  # start a new TEMP card
                list_fields += [None, 'TEMP', self.sid]
        return list_fields

    def repr_fields(self):
        """Writes the TEMP card"""
        return self.raw_fields()

    def get_loads(self):
        return [self]

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        if is_double:
            return self.comment + print_card_double(card)
        return self.comment + print_card_16(card)


class TEMPP1(BaseCard):
    """
    +--------+------+------+------+--------+------+------+------+
    |   1    |   2  |   3  |   4  |    5   |   6  |   7  |   8  |
    +========+======+======+======+========+======+======+======+
    | TEMPP1 | SID  | EID1 | TBAR | TPRIME |  T1  |  T2  |      |
    +--------+------+------+------+--------+------+------+------+
    |        | EID2 | EID3 | EID4 | EID5   | EID6 | EID7 | etc. |
    +--------+------+------+------+--------+------+------+------+

    +--------+------+------+------+--------+------+------+------+
    | TEMPP1 |  2   |  24  | 62.0 |  10.0  | 57.0 | 67.0 |      |
    +--------+------+------+------+--------+------+------+------+
    |        |  26  |  21  |  19  |   30   |      |      |      |
    +--------+------+------+------+--------+------+------+------+

    Alternate Form
    +--------+------+------+------+--------+------+------+------+
    |        | EID2 | THRU | EIDi |  EIDj  | THRU | EIDk |      |
    +--------+------+------+------+--------+------+------+------+
    |        |  1   | THRU |  10  |   30   | THRU |  61  |      |
    +--------+------+------+------+--------+------+------+------+

    """
    type = 'TEMPP1'

    @classmethod
    def _init_from_empty(cls):
        sid = 1
        eid = 1
        tbar = 2.0
        tprime = 1.0
        t_stress = 10.
        return TEMPP1(sid, eid, tbar, tprime, t_stress, comment='')

    def __init__(self, sid, eid, tbar, tprime, t_stress, comment=''):
        BaseCard.__init__(self)
        self.comment = comment
        self.sid = sid
        self.eid = eid
        self.tbar = tbar
        self.tprime = tprime
        self.t_stress = t_stress

    @classmethod
    def add_op2_data(cls, data, comment=''):
        """
        Adds a TEMPP1 card from the OP2

        Parameters
        ----------
        data : List[varies]
            a list of fields defined in OP2 format
        comment : str; default=''
            a comment for the card
        """
        sid, eid, t, tprime, ts1, ts2 = data
        return TEMPP1(sid, eid, t, tprime, [ts1, ts2], comment=comment)

    def raw_fields(self):
        """Writes the TEMP card"""
        list_fields = ['TEMPP1', self.sid, self.eid, self.tbar] + self.t_stress
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        list_fields = ['TEMPP1', self.sid, self.eid, self.tbar] + self.t_stress
        return print_card_8(list_fields)


class TEMPB3(BaseCard):
    """
    +--------+--------+--------+--------+-------+-------+--------+--------+--------+
    |   1    |    2   |    3   |    4   |   5   |   6   |    7   |    8   |   9    |
    +========+========+========+========+=======+=======+========+========+========+
    | TEMPB3 |   SID  |   EID  |  T(A)  |  T(B) |  T(C) | TPY(A) | TPZ(A) | TPY(B) |
    +--------+--------+--------+--------+-------+-------+--------+--------+--------+
    |        | TPZ(B) | TPY(C) | TPZ(C) | TC(A) | TD(A) |  TE(A) |  TF(A) |  TC(B) |
    +--------+--------+--------+--------+-------+-------+--------+--------+--------+
    |        |  TD(B) |  TE(B) |  TF(B) | TC(C) | TD(C) |  TE(C) |  TF(C) |        |
    +--------+--------+--------+--------+-------+-------+--------+--------+--------+
    |        | Element ID List
    +--------+--------+--------+--------+-------+-------+--------+--------+--------+
    """
    type = 'TEMPB3'

    @classmethod
    def _init_from_empty(cls):
        sid = 1
        eid = 1
        t = [0., 0., 0.]
        tpy = [0., 0., 0.]
        tpz = [0., 0., 0.]
        tc = [0., 0., 0.]
        td = [0., 0., 0.]
        te = [0., 0., 0.]
        tf = [0., 0., 0.]
        eids = [1, 2]
        return TEMPB3(sid, eid, t, tpy, tpz, tc, td, te, tf, eids, comment='')

    def __init__(self, sid, eid, t, tpy, tpz, tc, td, te, tf, eids, comment=''):
        BaseCard.__init__(self)
        self.comment = comment
        self.sid = sid
        self.eid = eid
        self.t = t
        self.tpy = tpy
        self.tpz = tpz
        self.tc = tc
        self.td = td
        self.te = te
        self.tf = tf
        self.eids = eids
        assert len(eids) > 0, str(self)
        self.eid_ref = None
        self.eids_ref = None

    @classmethod
    def add_card(cls, card, icard=0, comment=''):
        """
        Adds a TEMPD card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        icard : int; default=0
            sid to be parsed
        comment : str; default=''
            a comment for the card

        """
        sid = integer(card, 1, 'sid')
        eid = integer(card, 2, 'eid')

        t = []
        ifield = 3
        for var in ['A', 'B', 'C']:
            ti = double_or_blank(card, ifield, 'T(%s)' % var, default=0.0)
            t.append(ti)
            ifield += 1

        tpy = []
        tpz = []
        for var in ['A', 'B', 'C']:
            tpyi = double_or_blank(card, ifield, 'TPY(%s)' % var, default=0.0)
            tpzi = double_or_blank(card, ifield + 1, 'TPZ(%s)' % var, default=0.0)
            tpy.append(tpyi)
            tpz.append(tpzi)
            ifield += 2
        tc = []
        td = []
        te = []
        tf = []
        for var in ['A', 'B', 'C']:
            tci = double_or_blank(card, ifield, 'TC(%s)' % var, default=0.0)
            tdi = double_or_blank(card, ifield + 1, 'TD(%s)' % var, default=0.0)
            tei = double_or_blank(card, ifield + 2, 'TE(%s)' % var, default=0.0)
            tfi = double_or_blank(card, ifield + 3, 'TF(%s)' % var, default=0.0)
            tc.append(tci)
            td.append(tdi)
            te.append(tei)
            tf.append(tfi)
            ifield += 4
        ifield += 1 # skip the None

        list_fields = card[ifield:]
        eids = expand_thru_by(list_fields, set_fields=True, sort_fields=True,
                              require_int=True, allow_blanks=False)
        return TEMPB3(sid, eid, t, tpy, tpz, tc, td, te, tf, eids, comment=comment)

    def cross_reference(self, model: BDF) -> None:
        self.eid_ref = model.Element(self.eid)
        self.eids_ref = model.Elements(self.eids)

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.eid_ref = None
        self.eids_ref = None

    def get_loads(self):
        return [self]

    def raw_fields(self):
        """Writes the TEMP card"""
        list_fields = ['TEMPB3', self.sid, self.eid, ] + self.t
        for tpyi, tpzi in zip(self.tpy, self.tpz):
            list_fields.append(tpyi)
            list_fields.append(tpzi)
        for tci, tdi, tei, tfi in zip(self.tc, self.td, self.tf, self.tf):
            list_fields.append(tci)
            list_fields.append(tdi)
            list_fields.append(tei)
            list_fields.append(tfi)
        list_fields.append(None)
        list_fields += self.eids
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        list_fields = self.raw_fields()
        return print_card_8(list_fields)

# Loads
#-------------------------------------------------------
# Default Loads


class TEMPD(BaseCard):
    """
    Defines a temperature value for all grid points of the structural model
    that have not been given a temperature on a TEMP entry

    +-------+------+----+------+----+------+----+------+----+
    |   1   |  2   | 3  |  4   |  5 |  6   | 7  |  8   | 9  |
    +=======+======+====+======+====+======+====+======+====+
    | TEMPD | SID1 | T1 | SID2 | T2 | SID3 | T3 | SID4 | T4 |
    +-------+------+----+------+----+------+----+------+----+

    """
    type = 'TEMPD'

    @classmethod
    def _init_from_empty(cls):
        sid = 1
        temperatures = {1 : 1.0}
        return TEMPD(sid, temperatures, comment='')

    def __init__(self, sid, temperature, comment=''):
        """
        Creates a TEMPD card

        Parameters
        ----------
        sid : int
            Load set identification number. (Integer > 0)
        temperature : float
            default temperature
        comment : str; default=''
            a comment for the card

        """
        BaseCard.__init__(self)
        if comment:
            self.comment = comment
        self.sid = sid
        self.temperature = temperature

    @classmethod
    def add_card(cls, card, icard=0, comment=''):
        """
        Adds a TEMPD card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        icard : int; default=0
            sid to be parsed
        comment : str; default=''
            a comment for the card

        """
        nfields = len(card) - 1
        assert nfields % 2 == 0, 'card=%s' % card
        i = 2 * icard
        sid = integer(card, i + 1, 'sid')
        temperature = double(card, i + 2, 'temp')
        return TEMPD(sid, temperature, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        """
        Adds a TEMPD card from the OP2

        Parameters
        ----------
        data : List[varies]
            a list of fields defined in OP2 format
        comment : str; default=''
            a comment for the card

        """
        sid = data[0]
        temperature = data[1]
        return TEMPD(sid, temperature, comment=comment)

    def add(self, tempd_obj):
        for (lid, tempd) in tempd_obj.temperature.items():
            self.temperature[lid] = tempd

    def cross_reference(self, model: BDF) -> None:
        pass

    def safe_cross_reference(self, model, xref_errors):
        pass

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        pass

    def raw_fields(self):
        """Writes the TEMPD card"""
        list_fields = ['TEMPD', self.sid, self.temperature]
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.raw_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        if is_double:
            return self.comment + print_card_double(card)
        return self.comment + print_card_16(card)
