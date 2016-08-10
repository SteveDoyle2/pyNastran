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
    string, fields)


class ThermalLoadDefault(ThermalCard):
    def __init__(self, card, data):
        pass


class ThermalLoad(ThermalCard):
    def __init__(self):
        pass

class QVOL(ThermalLoad):
    """
    Defines a rate of volumetric heat addition in a conduction element.
    """
    type = 'QVOL'

    def __init__(self, sid, qvol, control_point, elements, comment=''):
        ThermalLoad.__init__(self)
        if comment:
            self._comment = comment
        #: Load set identification number. (Integer > 0)
        self.sid = sid
        self.qvol = qvol
        self.control_point = control_point
        self.elements = elements

    @classmethod
    def add_card(cls, card, comment=''):
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

    def getLoads(self):
        self.deprecated('getLoads()', 'get_loads()', '0.8')
        return self.get_loads()

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


class QBDY1(ThermalLoad):
    """
    Defines a uniform heat flux into CHBDYj elements.
    """
    type = 'QBDY1'

    def __init__(self, sid, qFlux, eids, comment=''):
        ThermalLoad.__init__(self)
        if comment:
            self._comment = comment
        #: Load set identification number. (Integer > 0)
        self.sid = sid
        #: Heat flux into element (FLOAT)
        self.qFlux = qFlux
        #: CHBDYj element identification numbers (Integer)
        #: .. todo:: use expand_thru_by ???
        assert len(eids) > 0
        self.eids = expand_thru(eids)

    @classmethod
    def add_card(cls, card, comment=''):
        sid = integer(card, 1, 'sid')
        qFlux = double(card, 2, 'qFlux')
        eids = []
        j = 1
        for i in range(3, len(card)):
            eid = integer_or_string(card, i, 'eid%i' % j)
            eids.append(eid)
            j += 1
        return QBDY1(sid, qFlux, eids, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        sid = data[0]
        qFlux = data[1]
        eids = data[2:]
        return QBDY1(sid, qFlux, eids, comment=comment)

    def getLoads(self):
        self.deprecated('getLoads()', 'get_loads()', '0.8')
        return self.get_loads()

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
        return len(self.qFlux)

    @property
    def element_ids(self):
        return self.Eids()

    def Eids(self):
        eids = []
        for eid in self.eids:
            eids.append(self._eid(eid))
        return eids

    def raw_fields(self):
        list_fields = ['QBDY1', self.sid, self.qFlux] + self.element_ids
        return list_fields

    def repr_fields(self):
        eids = collapse_thru_by(self.element_ids)
        list_fields = ['QBDY1', self.sid, self.qFlux] + eids
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
            self._comment = comment

        #: Load set identification number. (Integer > 0)
        self.sid = sid
        #: Identification number of an CHBDYj element. (Integer > 0)
        self.eid = eid
        #: Heat flux at the i-th grid point on the referenced CHBDYj
        #: element. (Real or blank)
        self.qFlux = qfluxs

    @classmethod
    def add_card(cls, card, comment=''):
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
        sid = data[0]
        eid = data[1]
        qfluxs = [data[2]]
        return QBDY2(sid, eid, qfluxs, comment=comment)

    def getLoads(self):
        self.deprecated('getLoads()', 'get_loads()', '0.8')
        return self.get_loads()

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
        return len(self.qFlux)

    def raw_fields(self):
        list_fields = ['QBDY2', self.sid, self.Eid()] + self.qFlux
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

    def __init__(self, sid, Q0, cntrlnd, eids, comment=''):
        ThermalLoad.__init__(self)
        if comment:
            self._comment = comment

        #: Load set identification number. (Integer > 0)
        self.sid = sid

        #: Heat flux into element
        self.Q0 = Q0

        #: Control point for thermal flux load. (Integer > 0; Default = 0)
        self.cntrlnd = cntrlnd

        #: CHBDYj element identification numbers
        self.eids = expand_thru_by(eids)

    @classmethod
    def add_card(cls, card, comment=''):
        sid = integer(card, 1, 'sid')
        Q0 = double(card, 2, 'Q0')
        cntrlnd = integer_or_blank(card, 3, 'cntrlnd', 0)

        nfields = card.nfields
        eids = fields(integer_or_string, card, 'eid', i=4, j=nfields)
        return QBDY3(sid, Q0, cntrlnd, eids, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        sid = data[0]
        Q0 = data[1]
        cntrlnd = data[2]
        eids = list(data[3:])
        return QBDY3(sid, Q0, cntrlnd, eids, comment=comment)

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
        list_fields = (['QBDY3', self.sid, self.Q0, self.cntrlnd] +
                       collapse_thru_by(eids))
        return list_fields

    def repr_fields(self):
        cntrlnd = set_blank_if_default(self.cntrlnd, 0)
        eids = self.Eids()
        eids.sort()
        list_fields = ['QBDY3', self.sid, self.Q0, cntrlnd] + collapse_thru_by(eids)
        return list_fields

    def getLoads(self):
        self.deprecated('getLoads()', 'get_loads()', '0.8')
        return self.get_loads()

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
        'POINT' : 1,
        'LINE' : 2,
        'REV' : 2,
        'AREA3' : 3,
        'AREA4' : 4,
        'AREA6' : 6,
        'AREA8' : 8,
    }

    def __init__(self, sid, flag, Q0, af, grids, comment=''):
        ThermalLoad.__init__(self)
        if comment:
            self._comment = comment

        #: Load set identification number. (Integer > 0)
        self.sid = sid
        self.flag = flag
        assert flag in ['POINT', 'LINE', 'REV', 'AREA3', 'AREA4',
                        'AREA6', 'AREA8']

        #: Magnitude of thermal flux into face. Q0 is positive for heat
        #: into the surface. (Real)
        self.Q0 = Q0

        #: Area factor depends on type. (Real > 0.0 or blank)
        self.af = af

        #: Grid point identification of connected grid points.
        #: (Integer > 0 or blank)
        self.grids = grids

    @classmethod
    def add_card(cls, card, comment=''):
        sid = integer(card, 1, 'eid')

        flag = string(card, 2, 'flag')
        assert flag in ['POINT', 'LINE', 'REV', 'AREA3', 'AREA4',
                        'AREA6', 'AREA8']

        Q0 = double(card, 3, 'Q0')
        af = double_or_blank(card, 4, 'af')
        nnodes = cls.flag_to_nnodes[flag]

        grids = []
        for i in range(nnodes):
            grid = integer(card, 5 + i, 'grid%i' % (i + 1))
            grids.append(grid)
        return QHBDY(sid, flag, Q0, af, grids, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        sid = data[0]
        flag = data[1]
        Q0 = data[2]
        af = data[3]
        grids = data[4:]
        return QHBDY(sid, flag, Q0, af, grids, comment=comment)

    def getLoads(self):
        self.deprecated('getLoads()', 'get_loads()', '0.8')
        return self.get_loads()

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
        list_fields = ['QHBDY', self.sid, self.flag, self.Q0, self.af] + self.grids
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
    | TEMP | SID | G1 |  T1   | G2 |  T2   | G3 | T3 |
    +------+-----+----+-------+----+-------+----+----+
    | TEMP |  3  | 94 | 316.2 | 49 | 219.8 |    |    |
    +------+-----+----+-------+----+-------+----+----+
    """
    type = 'TEMP'

    def __init__(self, sid, temperatures, comment=''):
        ThermalLoad.__init__(self)
        if comment:
            self._comment = comment
        #: Load set identification number. (Integer > 0)
        self.sid = sid

        #: dictionary of temperatures where the key is the grid ID (Gi)
        #: and the value is the temperature (Ti)
        self.temperatures = temperatures

    @classmethod
    def add_card(cls, card, comment=''):
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

    def getLoads(self):
        self.deprecated('getLoads()', 'get_loads()', '0.8')
        return self.get_loads()

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
        BaseCard.__init__(self)
        if comment:
            self._comment = comment
        self.sid = sid
        self.temperature = temperature

    @classmethod
    def add_card(cls, card, icard=0, comment=''):
        nfields = len(card) - 1
        assert nfields % 2 == 0, 'card=%s' % card
        i = 2 * icard
        sid = integer(card, i + 1, 'sid')
        temperature = double(card, i + 2, 'temp')
        return TEMPD(sid, temperature, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
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
