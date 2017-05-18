# pylint: disable=C0103,R0902,R0904,R0914,C0111
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
from six import  iteritems
from six.moves import range

from pyNastran.utils import integer_types
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


class ThermalLoadDefault(ThermalCard):
    def __init__(self, card, data):
        pass


class ThermalLoad(ThermalCard):
    def __init__(self):
        pass

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

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ' which is required by QVOL sid=%s' % self.sid
        self.elements = model.Elements(self.elements, msg=msg)
        self.elements_ref = self.elements

    def safe_cross_reference(self, model):
        try:
            return self.cross_reference(model)
        except KeyError:
            model.log.warning('failed cross-referencing\n%s' % str(self))
            raise

    def uncross_reference(self):
        self.elements = self.element_ids
        del self.elements_ref

    def _eid(self, eid):
        if isinstance(eid, integer_types):
            return eid
        return eid.eid

    @property
    def element_ids(self):
        return self.Eids()

    def Eids(self):
        eids = []
        for eid in self.elements:
            eids.append(self._eid(eid))
        return eids

    def raw_fields(self):
        list_fields = ['QVOL', self.sid, self.qvol, self.control_point] + self.element_ids
        return list_fields

    def repr_fields(self):
        eids = collapse_thru_by(self.element_ids)
        list_fields = ['QVOL', self.sid, self.qvol, self.control_point] + eids
        return list_fields

    def write_card(self, size=8, is_double=False):
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

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ' which is required by QVECT sid=%s' % self.sid
        self.eids = model.Elements(self.eids, msg=msg)
        self.eidss_ref = self.eids

    def safe_cross_reference(self, model):
        try:
            return self.cross_reference(model)
        except KeyError:
            model.log.warning('failed cross-referencing\n%s' % str(self))

    def uncross_reference(self):
        self.eids = self.element_ids
        del self.elements_ref

    def _eid(self, eid):
        if isinstance(eid, integer_types):
            return eid
        return eid.eid

    @property
    def element_ids(self):
        return self.Eids()

    def Eids(self):
        eids = []
        for eid in self.eids:
            eids.append(self._eid(eid))
        return eids

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

    def write_card(self, size=8, is_double=False):
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

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ' which is required by QBDY1 sid=%s' % self.sid
        self.eids = model.Elements(self.eids, msg=msg)
        self.eids_ref = self.eids

    def safe_cross_reference(self, model):
        try:
            return self.cross_reference(model)
        except KeyError:
            model.log.warning('failed cross-referencing\n%s' % str(self))

    def uncross_reference(self):
        self.eids = self.element_ids
        del self.eids_ref

    def _eid(self, eid):
        if isinstance(eid, integer_types):
            return eid
        return eid.eid

    def nQFluxTerms(self):
        return len(self.qflux)

    @property
    def element_ids(self):
        return self.Eids()

    def Eids(self):
        eids = []
        for eid in self.eids:
            eids.append(self._eid(eid))
        return eids

    def raw_fields(self):
        list_fields = ['QBDY1', self.sid, self.qflux] + self.element_ids
        return list_fields

    def repr_fields(self):
        eids = collapse_thru_by(self.element_ids)
        list_fields = ['QBDY1', self.sid, self.qflux] + eids
        return list_fields

    def write_card(self, size=8, is_double=False):
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

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ' which is required by QBDY2 sid=%s' % self.sid
        self.eid = model.Element(self.eid, msg=msg)
        self.eid_ref = self.eid

    def safe_cross_reference(self, model):
        try:
            return self.cross_reference(model)
        except KeyError:
            model.log.warning('failed cross-referencing\n%s' % str(self))

    def uncross_reference(self):
        self.eid = self.Eid()
        del self.eid_ref

    def Eid(self):
        if isinstance(self.eid, integer_types):
            return self.eid
        return self.eid.eid

    def nQFluxTerms(self):
        return len(self.qfluxs)

    def raw_fields(self):
        list_fields = ['QBDY2', self.sid, self.Eid()] + self.qfluxs
        return list_fields

    def repr_fields(self):
        return self.raw_fields()

    def write_card(self, size=8, is_double=False):
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

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ' which is required by QBDY3 sid=%s' % self.sid
        for i, eid in enumerate(self.eids):
            self.eids[i] = model.Element(eid, msg=msg)
        self.eids_ref = self.eids

    def safe_cross_reference(self, model):
        try:
            return self.cross_reference(model)
        except KeyError:
            model.log.warning('failed cross-referencing\n%s' % str(self))

    def uncross_reference(self):
        self.eids = self.element_ids
        del self.eids_ref

    def Eids(self):
        eids = []
        for eid in self.eids:
            eids.append(self.Eid(eid))
        return eids

    def Eid(self, eid):
        if isinstance(eid, integer_types):
            return eid
        return eid.eid

    def raw_fields(self):
        eids = self.Eids()
        eids.sort()
        list_fields = (['QBDY3', self.sid, self.q0, self.cntrlnd] +
                       collapse_thru_by(eids))
        return list_fields

    def repr_fields(self):
        cntrlnd = set_blank_if_default(self.cntrlnd, 0)
        eids = self.Eids()
        eids.sort()
        list_fields = ['QBDY3', self.sid, self.q0, cntrlnd] + collapse_thru_by(eids)
        return list_fields

    def get_loads(self):
        """
        .. todo:: return loads
        """
        return []

    def write_card(self, size=8, is_double=False):
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
    flag_to_nnodes = {
        'POINT' : (1, 1),
        'LINE' : (2, 2),
        'REV' : (2, 2),
        'AREA3' : (3, 3),
        'AREA4' : (4, 4),
        'AREA6' : (4, 6), # 4-6
        'AREA8' : (5, 8), # 5-8
    }

    def __init__(self, sid, flag, q0, grids, af=None, comment=''):
        ThermalLoad.__init__(self)
        if comment:
            self.comment = comment

        #: Load set identification number. (Integer > 0)
        self.sid = sid
        self.flag = flag
        assert flag in ['POINT', 'LINE', 'REV', 'AREA3', 'AREA4', 'AREA6', 'AREA8']

        #: Magnitude of thermal flux into face. Q0 is positive for heat
        #: into the surface. (Real)
        self.q0 = q0

        #: Area factor depends on type. (Real > 0.0 or blank)
        self.af = af

        #: Grid point identification of connected grid points.
        #: (Integer > 0 or blank)
        self.grids = grids

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
        assert flag in ['POINT', 'LINE', 'REV', 'AREA3', 'AREA4',
                        'AREA6', 'AREA8']

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
        flag = data[1]
        q0 = data[2]
        af = data[3]
        grids = data[4:]
        return QHBDY(sid, flag, q0, grids, af=af, comment=comment)

    def get_loads(self):
        return [self]

    def cross_reference(self, model):
        pass

    def safe_cross_reference(self, model):
        try:
            return self.cross_reference(model)
        except KeyError:
            model.log.warning('failed cross-referencing\n%s' % str(self))

    def uncross_reference(self):
        pass

    def raw_fields(self):
        list_fields = ['QHBDY', self.sid, self.flag, self.q0, self.af] + self.grids
        return list_fields

    def repr_fields(self):
        return self.raw_fields()

    def write_card(self, size=8, is_double=False):
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

        temperatures = {}
        for i in range(ntemps):
            n = i * 2 + 2
            gi = integer(card, n, 'g' + str(i))
            Ti = double(card, n + 1, 'T' + str(i))
            temperatures[gi] = Ti
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
        for (gid, temp) in iteritems(temp_obj.temperatures):
            self.temperatures[gid] = temp

    def cross_reference(self, model):
        pass

    def safe_cross_reference(self, model, debug=True):
        pass

    def uncross_reference(self):
        pass

    def raw_fields(self):
        """Writes the TEMP card"""
        list_fields = ['TEMP', self.sid]
        ntemps = len(self.temperatures) - 1
        for i, (gid, temp) in enumerate(sorted(iteritems(self.temperatures))):
            list_fields += [gid, temp]
            if i % 3 == 2 and ntemps > i:  # start a new TEMP card
                list_fields += [None, 'TEMP', self.sid]
        return list_fields

    def repr_fields(self):
        """Writes the TEMP card"""
        return self.raw_fields()

    def get_loads(self):
        """
        .. todo:: return loads
        """
        return []

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        if is_double:
            return self.comment + print_card_double(card)
        return self.comment + print_card_16(card)


class TEMPP1(BaseCard):
    type = 'TEMPP1'
    def __init__(self, sid, eid, tbar, tprime, t_stress, comment=''):
        self.comment = comment
        self.sid = sid
        self.eid = eid
        self.tbar = tbar
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
        return TEMPP1(sid, eid, t, tprime, [ts1, ts2])

    def write_card(self, size=8, is_double=False):
        list_fields = ['TEMPP1', self.sid, self.eid, self.tbar] + self.t_stress
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
        for (lid, tempd) in iteritems(tempd_obj.temperatures):
            self.temperatures[lid] = tempd

    def cross_reference(self, model):
        pass

    def uncross_reference(self):
        pass

    def raw_fields(self):
        """Writes the TEMPD card"""
        list_fields = ['TEMPD', self.sid, self.temperature]
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.raw_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        if is_double:
            return self.comment + print_card_double(card)
        return self.comment + print_card_16(card)
