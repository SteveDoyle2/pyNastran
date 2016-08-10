# coding: utf-8
"""
All dynamic loads are defined in this file.  This includes:

 * ACSRCE
 * DLOAD
 * TLOAD1
 * TLOAD2
 * RLOAD1
 * RLOAD2
"""
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
from six.moves import zip, range
import numpy as np

from pyNastran.utils import integer_types, integer_float_types
from pyNastran.bdf.field_writer_8 import set_blank_if_default
from pyNastran.bdf.bdf_interface.assign_type import (
    integer, integer_or_blank, double_or_blank, integer_string_or_blank,
    integer_double_or_blank, double)
from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.field_writer_16 import print_card_16
from pyNastran.bdf.field_writer_double import print_card_double
from pyNastran.bdf.cards.loads.loads import TabularLoad, LoadCombination, BaseCard


class ACSRCE(BaseCard):
    r"""
    Defines acoustic source as a function of power vs. frequency.

    +-----+----------+---------------+-----------------+-------+-----+---+
    |  1  |    2     |       3       |        4        |   5   |  6  | 7 |
    +=====+==========+===============+=================+=======+=====+===+
    | SID | EXCITEID | DELAYI/DELAYR | DPHASEI/DPHASER | TP/RP | RHO | B |
    +-----+----------+---------------+-----------------+-------+-----+---+

    ..math ::
      C = \sqrt(B ⁄ ρ)
      Source Strength = {A} * 1/(2πf)  * \sqrt( 8πC P(f) / ρ) ^ (ei(θ + 2πfτ))

    """
    type = 'ACSRCE'

    def __init__(self, sid, excite_id, delay, dphase, power, rho, b, comment=''):
        if comment:
            self._comment = comment
        self.sid = sid
        self.excite_id = excite_id
        self.delay = delay
        self.dphase = dphase
        self.power = power
        self.rho = rho
        self.b = b

    @classmethod
    def add_card(cls, card, comment=''):
        sid = integer(card, 1, 'sid')
        excite_id = integer(card, 2, 'excite_id') # DAREA, FBALOAD, SLOAD
        delay = integer_double_or_blank(card, 3, 'delay', 0) # DELAY, FBADLAY
        dphase = integer_double_or_blank(card, 4, 'dphase', 0) # DPHASE, FBAPHAS
        power = integer_double_or_blank(card, 5, 'power/tp/rp', 0) # TABLEDi/power
        rho = double(card, 6, 'rho')
        b = double(card, 7, 'bulk modulus')

        assert len(card) <= 8, 'len(ACSRCE card) = %i\n%s' % (len(card), card)
        return ACSRCE(sid, excite_id, delay, dphase, power, rho, b, comment=comment)

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ' which is required by %s=%s' % (self.type, self.sid)
        if isinstance(self.delay, integer_types) and self.delay > 0:
            self.delay = model.DELAY(self.delay, msg=msg)
            self.delay_ref = self.delay
        if isinstance(self.dphase, integer_types) and self.dphase > 0:
            self.dphase = model.Table(self.dphase, msg=msg)
            self.dphase_ref = self.dphase
        if isinstance(self.power, integer_types):
            print('power')
            self.power = model.Table(self.power, msg=msg)
            self.power_ref = self.power

        #load_ids2 = []
        #for load_id in self.load_ids:
            #load_id2 = model.get_dload_entries(load_id, msg=msg)
            #load_ids2.append(load_id2)
        #self.load_ids = load_ids2
        #self.load_ids_ref = self.load_ids

    def safe_cross_reference(self, model):
        return self.cross_reference(model)

    #def uncross_reference(self):
        #self.load_ids = [self.LoadID(load) for load in self.load_ids]
        #del self.load_ids_ref

    def Delay(self):
        #if self.delay in [0, 0.0]:
            #return 0
        if isinstance(self.delay, integer_float_types):
            return self.delay
        #print('self.delay =', self.delay)
        return self.delay_ref.tid

    def DPhase(self):
        #if self.dphase in [0, 0.0]:
            #return 0
        if isinstance(self.dphase, integer_float_types):
            return self.dphase
        #print('self.dphase =', self.dphase)
        return self.dphase_ref.tid

    def Power(self):
        #if self.power in [0, 0.0]:
            #return 0
        if isinstance(self.power, integer_float_types):
            return self.power
        #print('self.power =', self.power)
        return self.power_ref.tid

    def get_load_at_freq(self, freq):
        r"""
        ..math ::
          C = \sqrt(B ⁄ ρ)
          Source_strength = {A} * 1/(2πf)  * \sqrt( 8πC P(f) / ρ) ^ (ei(θ + 2πfτ))
        """
        C = np.sqrt(self.b / self.rho)
        ei = np.exp(1) * 1.j
        A = 0.0
        pi = np.pi
        if self.delay in [0, 0.]:
            tau = 0.
        else:
            print('delay\n', self.delay)
            tau = self.delay_ref.value
        Pf = self.power_ref.interpolate(freq)
        if self.dphase in [0, 0.]:
            theta = 0.
        else:
            print('dphase\n', self.dphase)
            theta = self.dphase_ref.interpolate(freq)
        strength = A / (2.* pi * freq) * np.sqrt(8*pi*C*Pf / self.rho) ** (ei*(theta + 2*pi*freq*tau))

        return 0.0

    def raw_fields(self):
        list_fields = ['ACSRCE', self.sid, self.excite_id, self.Delay(), self.DPhase(),
                       self.Power(), self.rho, self.b]
        return list_fields

    def repr_fields(self):
        return self.raw_fields()

    def write_card(self, size=8, is_double=False):
        card = self.raw_fields()
        if size == 16:
            return self.comment + print_card_16(card)
        return self.comment + print_card_8(card)


class DLOAD(LoadCombination):
    type = 'DLOAD'

    def __init__(self, sid, scale, scale_factors, load_ids, comment=''):
        LoadCombination.__init__(self, sid, scale, scale_factors, load_ids,
                                 comment=comment)

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        load_ids2 = []
        msg = ' which is required by %s=%s' % (self.type, self.sid)
        for load_id in self.load_ids:
            load_id2 = model.get_dload_entries(load_id, msg=msg)
            load_ids2.append(load_id2)
        self.load_ids = load_ids2
        self.load_ids_ref = self.load_ids

    def uncross_reference(self):
        self.load_ids = [self.LoadID(load) for load in self.load_ids]
        del self.load_ids_ref

    def raw_fields(self):
        list_fields = ['DLOAD', self.sid, self.scale]
        for (scale_factor, load_id) in zip(self.scale_factors, self.load_ids):
            list_fields += [scale_factor, self.LoadID(load_id)]
        return list_fields

    def repr_fields(self):
        return self.raw_fields()

    def write_card(self, size=8, is_double=False):
        card = self.raw_fields()
        if size == 16:
            return self.comment + print_card_16(card)
        return self.comment + print_card_8(card)


class RLOAD1(TabularLoad):
    r"""
    Defines a frequency-dependent dynamic load of the form
    for use in frequency response problems.

    .. math::
      \left\{ P(f)  \right\}  = \left\{A\right\} [ C(f)+iD(f)]
         e^{  i \left\{\theta - 2 \pi f \tau \right\} }

    +--------+-----+----------+-------+--------+----+----+------+
    | RLOAD1 | SID | EXCITEID | DELAY | DPHASE | TC | TD | TYPE |
    +--------+-----+----------+-------+--------+----+----+------+
    | RLOAD1 | 5   | 3        |       |        | 1  |    |      |
    +--------+-----+----------+-------+--------+----+----+------+
    """
    type = 'RLOAD1'

    def __init__(self, sid, excite_id, delay, dphase, tc, td, Type, comment=''):
        TabularLoad.__init__(self)
        if comment:
            self._comment = comment
        self.sid = sid
        self.excite_id = excite_id
        self.delay = delay
        self.dphase = dphase
        self.tc = tc
        self.td = td
        self.Type = Type

    def validate(self):
        msg = ''
        is_failed = False
        if self.tc > 0 or self.td > 0:
            msg += 'either RLOAD TC or TD > 0; tc=%s td=%s\n' % (self.tc, self.td)

        if self.Type in [0, 'L', 'LO', 'LOA', 'LOAD']:
            self.Type = 'LOAD'
        elif self.Type in [1, 'D', 'DI', 'DIS', 'DISP']:
            self.Type = 'DISP'
        elif self.Type in [2, 'V', 'VE', 'VEL', 'VELO']:
            self.Type = 'VELO'
        elif self.Type in [3, 'A', 'AC', 'ACC', 'ACCE']:
            self.Type = 'ACCE'
        else:
            msg += 'invalid RLOAD1 type  Type=%r\n' % self.Type
            is_failed = True

        if is_failed:
            msg += str(self)
            raise RuntimeError(msg)

    @classmethod
    def add_card(cls, card, comment=''):
        sid = integer(card, 1, 'sid')
        excite_id = integer(card, 2, 'excite_id')
        delay = integer_double_or_blank(card, 3, 'delay', 0)
        dphase = integer_double_or_blank(card, 4, 'dphase', 0)
        tc = integer_double_or_blank(card, 5, 'tc', 0)
        td = integer_double_or_blank(card, 6, 'td', 0)
        Type = integer_string_or_blank(card, 7, 'Type', 'LOAD')

        assert len(card) <= 8, 'len(RLOAD1 card) = %i\ncard=%s' % (len(card), card)
        return RLOAD1(sid, excite_id, delay, dphase, tc, td, Type, comment=comment)

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ' which is required by RLOAD1 sid=%s' % (self.sid)
        if self.tc > 0:
            self.tc = model.Table(self.tc, msg=msg)
            self.tc_ref = self.tc
        if self.td > 0:
            self.td = model.Table(self.td, msg=msg)
            self.td_ref = self.td
        if isinstance(self.delay, integer_types) and self.delay > 0:
            self.delay = model.DELAY(self.delay_id, msg=msg)
            self.delay_ref = self.delay
        if isinstance(self.dphase, integer_types) and self.dphase > 0:
            self.dphase = model.DPHASE(self.dphase, msg=msg)
            self.dphase_ref = self.dphase

    def safe_cross_reference(self, model):
        msg = ' which is required by RLOAD1 sid=%s' % (self.sid)
        if self.tc > 0:
            self.tc = model.Table(self.tc, msg=msg)
            self.tc_ref = self.tc
        if self.td > 0:
            self.td = model.Table(self.td, msg=msg)
            self.td_ref = self.td
        if isinstance(self.delay, integer_types) and self.delay > 0:
            self.delay = model.DELAY(self.delay_id, msg=msg)
            self.delay_ref = self.delay
        if isinstance(self.dphase, integer_types) and self.dphase > 0:
            self.dphase = model.DPHASE(self.dphase, msg=msg)
            self.dphase_ref = self.dphase

    def uncross_reference(self):
        self.tc = self.Tc()
        self.td = self.Td()
        self.delay = self.delay_id
        if self.tc > 0:
            del self.tc_ref
        if self.td > 0:
            del self.td_ref
        if isinstance(self.delay, integer_types) and self.delay > 0:
            del self.delay_ref

    def getLoads(self):
        self.deprecated('getLoads()', 'get_loads()', '0.8')
        return self.get_loads()

    def get_loads(self):
        return [self]

    def Tc(self):
        if self.tc in [0, 0.0]:
            return 0
        elif isinstance(self.tc, integer_types):
            return self.tc
        return self.tc.tid

    def Td(self):
        if self.td in [0, 0.0]:
            return 0
        elif isinstance(self.td, integer_types):
            return self.td
        return self.td.tid

    @property
    def delay_id(self):
        if self.delay in [0, 0.0]:
            return 0
        elif isinstance(self.delay, integer_float_types):
            return self.delay
        return self.delay_ref.sid

    @property
    def dphase_id(self):
        #if self.dphase in [0, 0.0]:
            #return 0
        if isinstance(self.dphase, integer_types):
            return self.dphase
        return self.dphase_ref.sid

    def get_load_at_freq(self, freq, scale=1.):
        # A = 1. # points to DAREA or SPCD
        if isinstance(self.tc, float):
            c = float(self.tc)
        elif self.tc == 0:
            c = 0.
        else:
            c = self.tc.interpolate(freq)

        if isinstance(self.td, float):
            d = float(self.td)
        elif self.td == 0:
            d = 0.
        else:
            d = self.td.interpolate(freq)

        if isinstance(self.dphase, float):
            dphase = self.dphase
        elif self.dphase == 0 or self.dphase is None:
            dphase = 0.0
        else:
            #print('DPHASE is not supported; type=%s' % type(self.dphase))
            #print('dphase = ', self.dphase)
            nids, comps, dphases = self.dphase_ref.get_dphase_at_freq(freq)
            assert len(dphases) == 1, dphases
            dphase = dphases[0]

        if isinstance(self.delay, float):
            tau = self.delay
        elif self.delay == 0 or self.dphase is None:
            tau = 0.0
        else:
            tau = self.delay.get_delay_at_freq(freq)

        out = (c + 1.j * d) * np.exp(dphase - 2 * np.pi * freq * tau)
        return out

    def raw_fields(self):
        list_fields = ['RLOAD1', self.sid, self.excite_id, self.delay_id, self.dphase_id,
                       self.Tc(), self.Td(), self.Type]
        return list_fields

    def repr_fields(self):
        Type = set_blank_if_default(self.Type, 'LOAD')
        list_fields = ['RLOAD1', self.sid, self.excite_id, self.delay_id, self.dphase_id,
                       self.Tc(), self.Td(), Type]
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        if is_double:
            return self.comment + print_card_double(card)
        return self.comment + print_card_16(card)


class RLOAD2(TabularLoad):
    r"""
    Defines a frequency-dependent dynamic load of the form
    for use in frequency response problems.

    .. math:: \left\{ P(f)  \right\}  = \left\{A\right\} * B(f)
        e^{  i \left\{ \phi(f) + \theta - 2 \pi f \tau \right\} }

    +--------+-----+----------+-------+--------+----+----+------+
    | RLOAD2 | SID | EXCITEID | DELAY | DPHASE | TB | TP | TYPE |
    +--------+-----+----------+-------+--------+----+----+------+
    | RLOAD2 | 5   | 3        |       |        | 1  |    |      |
    +--------+-----+----------+-------+--------+----+----+------+
    """
    type = 'RLOAD2'

    def __init__(self, sid, excite_id, delay, dphase, tb, tp, Type, comment=''):
        TabularLoad.__init__(self)
        if comment:
            self._comment = comment
        self.sid = sid
        self.excite_id = excite_id
        self.delay = delay
        self.dphase = dphase
        self.tb = tb
        self.tp = tp
        self.Type = Type

    @classmethod
    def add_card(cls, card, comment=''):
        sid = integer(card, 1, 'sid')
        excite_id = integer(card, 2, 'excite_id')
        delay = integer_double_or_blank(card, 3, 'delay', 0)
        dphase = integer_double_or_blank(card, 4, 'dphase', 0)
        tb = integer_or_blank(card, 5, 'tb', 0)
        tp = integer_or_blank(card, 6, 'tp', 0)
        Type = integer_string_or_blank(card, 7, 'Type', 'LOAD')

        if Type in [0, 'L', 'LO', 'LOA', 'LOAD']:
            Type = 'LOAD'
        elif Type in [1, 'D', 'DI', 'DIS', 'DISP']:
            Type = 'DISP'
        elif Type in [2, 'V', 'VE', 'VEL', 'VELO']:
            Type = 'VELO'
        elif Type in [3, 'A', 'AC', 'ACC', 'ACCE']:
            Type = 'ACCE'
        else:
            msg = 'invalid RLOAD2 type  Type=%r' % Type
            raise RuntimeError(msg)
        assert len(card) <= 8, 'len(RLOAD2 card) = %i\ncard=%s' % (len(card), card)
        return RLOAD2(sid, excite_id, delay, dphase, tb, tp, Type, comment=comment)

    def get_load_at_freq(self, freq, scale=1.):
        # A = 1. # points to DAREA or SPCD
        if isinstance(self.tb, float):
            b = self.tb
        elif self.tb == 0:
            b = 0.0
        else:
            b = self.tb_ref.interpolate(freq)

        if isinstance(self.tp, float):
            p = self.tp
        elif self.tp == 0:
            p = 0.0
        else:
            p = self.tp_ref.interpolate(freq)

        if isinstance(self.dphase, float):
            dphase = self.dphase
        elif self.dphase == 0 or self.dphase is None:
            dphase = 0.0
        else:
            print('self.dphase = ', self.dphase)
            raise NotImplementedError('DPHASE is not supported')
            dphase = self.dphase_ref.get_dphase_at_freq(freq)

        if isinstance(self.delay, float):
            tau = self.delay
        elif self.delay == 0:
            tau = 0.0
        else:
            #raise NotImplementedError('DELAY is not supported')
            nids, comps, taus = self.delay_ref.get_delay_at_freq(freq)
            assert len(taus) == 1, taus
            tau = taus[0]

        try:
            out = b * np.exp(1.j * p + dphase - 2 * np.pi * freq * tau)
        except TypeError:
            print('b =', b)
            print('p =', p)
            print('dphase =', dphase)
            print('freq =', freq)
            print('tau =', tau)
            raise
        return out

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ' which is required by RLOAD2=%s' % (self.sid)
        if self.tb:
            self.tb = model.Table(self.tb, msg=msg)
            self.tb_ref = self.tb
        if self.tp:
            self.tp = model.Table(self.tp, msg=msg)
            self.tp_ref = self.tp
        if isinstance(self.delay, integer_types) and self.delay > 0:
            self.delay = model.DELAY(self.delay, msg=msg)
            self.delay_ref = self.delay
        if isinstance(self.dphase, integer_types) and self.dphase > 0:
            self.dphase = model.DPHASE(self.dphase, msg=msg)
            self.dphase_ref = self.dphase

    def safe_cross_reference(self, model):
        msg = ' which is required by RLOAD2=%s' % (self.sid)
        if self.tb:
            self.tb = model.Table(self.tb, msg=msg)
            self.tb_ref = self.tb
        if self.tp:
            self.tp = model.Table(self.tp, msg=msg)
            self.tp_ref = self.tp
        if self.delay > 0:
            self.delay = model.DELAY(self.delay, msg=msg)
            self.delay_ref = self.delay
        if self.dphase:
            self.dphase = model.DPHASE(self.dphase, msg=msg)
            self.dphase_ref = self.dphase

    def uncross_reference(self):
        self.tb = self.Tb()
        self.tp = self.Tp()
        self.delay = self.delay_id
        if self.tb > 0:
            del self.tb_ref
        if self.tp > 0:
            del self.tp_ref
        if self.delay > 0:
            del self.delay_ref

    def getLoads(self):
        self.deprecated('getLoads()', 'get_loads()', '0.8')
        return self.get_loads()

    def get_loads(self):
        return [self]

    def LoadID(self):
        return self.sid

    def Tb(self):
        if self.tb == 0:
            return 0
        elif isinstance(self.tb, integer_types):
            return self.tb
        return self.tb_ref.tid

    def Tp(self):
        if self.tp == 0:
            return 0
        elif isinstance(self.tp, integer_types):
            return self.tp
        return self.tp_ref.tid

    @property
    def delay_id(self):
        if self.delay in [0, 0.]:
            return 0
        elif isinstance(self.delay, integer_types):
            return self.delay
        return self.delay.sid

    @property
    def dphase_id(self):
        #if self.dphase in [0, 0.0]:
            #return 0
        print(self.dphase)
        if isinstance(self.dphase, integer_types):
            return self.dphase
        return self.dphase_ref.sid

    def raw_fields(self):
        list_fields = ['RLOAD2', self.sid, self.excite_id, self.delay_id, self.dphase_id,
                       self.Tb(), self.Tp(), self.Type]
        return list_fields

    def repr_fields(self):
        Type = set_blank_if_default(self.Type, 0.0)
        list_fields = ['RLOAD2', self.sid, self.excite_id, self.delay_id, self.dphase_id,
                       self.Tb(), self.Tp(), Type]
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        if is_double:
            return self.comment + print_card_double(card)
        return self.comment + print_card_16(card)


class TLOAD1(TabularLoad):
    r"""
    Transient Response Dynamic Excitation, Form 1

    Defines a time-dependent dynamic load or enforced motion of the form:

    .. math::
      \left\{ P(t) \right\} = \left\{ A \right\} \cdot F(t-\tau)

    for use in transient response analysis.
    """
    type = 'TLOAD1'

    def __init__(self, sid, excite_id, delay, Type, tid, us0, vs0, comment=''):
        TabularLoad.__init__(self)
        if comment:
            self._comment = comment
        #: load ID
        self.sid = sid

        #: Identification number of DAREA or SPCD entry set or a thermal load
        #: set (in heat transfer analysis) that defines {A}. (Integer > 0)
        self.excite_id = excite_id

        #: If it is a non-zero integer, it represents the
        #: identification number of DELAY Bulk Data entry that defines .
        #: If it is real, then it directly defines the value of that will
        #: be used for all degrees-of-freedom that are excited by this
        #: dynamic load entry.  See also Remark 9. (Integer >= 0,
        #: real or blank)
        self.delay = delay

        #: Defines the type of the dynamic excitation. (LOAD,DISP, VELO, ACCE)
        self.Type = Type

        #: Identification number of TABLEDi entry that gives F(t). (Integer > 0)
        self.tid = tid

        #: Factor for initial displacements of the enforced degrees-of-freedom.
        #: (Real; Default = 0.0)
        self.us0 = us0
        #: Factor for initial velocities of the enforced degrees-of-freedom.
        #: (Real; Default = 0.0)
        self.vs0 = vs0

    @classmethod
    def add_card(cls, card, comment=''):
        sid = integer(card, 1, 'sid')
        excite_id = integer(card, 2, 'excite_id')
        delay = integer_double_or_blank(card, 3, 'delay', 0)
        Type = integer_string_or_blank(card, 4, 'Type', 'LOAD')
        tid = integer(card, 5, 'tid')
        us0 = double_or_blank(card, 6, 'us0', 0.0)
        vs0 = double_or_blank(card, 7, 'vs0', 0.0)

        if Type in [0, 'L', 'LO', 'LOA', 'LOAD']:
            Type = 'LOAD'
        elif Type in [1, 'D', 'DI', 'DIS', 'DISP']:
            Type = 'DISP'
        elif Type in [2, 'V', 'VE', 'VEL', 'VELO']:
            Type = 'VELO'
        elif Type in [3, 'A', 'AC', 'ACC', 'ACCE']:
            Type = 'ACCE'
        else:
            msg = 'invalid TLOAD1 type  Type=%r' % Type
            raise RuntimeError(msg)
        assert len(card) <= 8, 'len(TLOAD1 card) = %i\ncard=%s' % (len(card), card)
        return TLOAD1(sid, excite_id, delay, Type, tid, us0, vs0, comment=comment)

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
        msg = ' which is required by %s=%s' % (self.type, self.sid)
        if self.tid:
            self.tid = model.Table(self.tid, msg=msg)
            self.tid_ref = self.tid
        if isinstance(self.delay, integer_types) and self.delay > 0:
            self.delay = model.DELAY(self.delay, msg=msg)
            self.delay_ref = self.delay

    def safe_cross_reference(self, model, debug=True):
        msg = ' which is required by %s=%s' % (self.type, self.sid)
        if self.tid:
            #try:
            self.tid = model.Table(self.tid, msg=msg)
            self.tid_ref = self.tid
            #except
        if isinstance(self.delay, integer_types) and self.delay > 0:
            self.delay = model.DELAY(self.delay_id, msg=msg)
            self.delay_ref = self.delay

    def uncross_reference(self):
        self.tid = self.Tid()
        self.delay = self.delay_id
        if self.tid > 0:
            del self.tid_ref
        if self.delay > 0:
            del self.delay_ref

    def Tid(self):
        if self.tid == 0:
            return 0
        elif isinstance(self.tid, integer_types):
            return self.tid
        return self.tid_ref.tid

    @property
    def delay_id(self):
        if self.delay in [0, 0.]:
            return 0
        elif isinstance(self.delay, integer_float_types):
            return self.delay
        return self.delay_ref.sid

    def get_load_at_time(self, time, scale=1.):
        # A = 1. # points to DAREA or SPCD
        if isinstance(time, float):
            time = np.array([time])
        else:
            time = np.asarray(time)

        if isinstance(self.delay, float):
            tau = self.delay
        elif self.delay == 0 or self.delay is None:
            tau = 0.0
        else:
            #raise NotImplementedError('DELAY is not supported')
            tau = self.delay.get_delay_at_time(time)

        i = np.where(time - tau > 0)
        time2 = time[i]
        resp = self.tid_ref.interpolate(time2)
        is_spcd = False
        if self.Type == 'VELO' and is_spcd:
            resp[0] = self.us0
        if self.Type == 'ACCE' and is_spcd:
            resp[0] = self.vs0
        return resp

    def raw_fields(self):
        list_fields = ['TLOAD1', self.sid, self.excite_id, self.delay_id, self.Type,
                       self.Tid(), self.us0, self.vs0]
        return list_fields

    def repr_fields(self):
        us0 = set_blank_if_default(self.us0, 0.0)
        vs0 = set_blank_if_default(self.vs0, 0.0)
        list_fields = ['TLOAD1', self.sid, self.excite_id, self.delay_id, self.Type,
                       self.Tid(), us0, vs0]
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        if is_double:
            return self.comment + print_card_double(card)
        return self.comment + print_card_16(card)


class TLOAD2(TabularLoad):
    r"""
    Transient Response Dynamic Excitation, Form 1

    Defines a time-dependent dynamic load or enforced motion of the form:

    .. math::
      \left\{ P(t) \right\} = \left\{ A \right\} \cdot F(t-\tau)

    for use in transient response analysis.
    """
    type = 'TLOAD2'

    def __init__(self, sid, excite_id, delay, Type, T1, T2,
                 frequency, phase, c, b, us0, vs0, comment=''):
        TabularLoad.__init__(self)
        if comment:
            self._comment = comment
        #: load ID
        #: SID must be unique for all TLOAD1, TLOAD2, RLOAD1, RLOAD2, and ACSRCE entries.
        self.sid = sid
        self.excite_id = excite_id
        self.delay = delay

        #: Defines the type of the dynamic excitation. (Integer; character
        #: or blank; Default = 0)
        self.Type = Type

        #: Time constant. (Real >= 0.0)
        self.T1 = T1
        #: Time constant. (Real; T2 > T1)
        self.T2 = T2
        #: Frequency in cycles per unit time. (Real >= 0.0; Default = 0.0)
        self.frequency = frequency
        #: Phase angle in degrees. (Real; Default = 0.0)
        self.phase = phase
        #: Exponential coefficient. (Real; Default = 0.0)
        self.c = c
        #: Growth coefficient. (Real; Default = 0.0)
        self.b = b
        #: Factor for initial displacements of the enforced degrees-of-freedom.
        #: (Real; Default = 0.0)
        self.us0 = us0
        #: Factor for initial velocities of the enforced degrees-of-freedom
        #: (Real; Default = 0.0)
        self.vs0 = vs0

    @classmethod
    def add_card(cls, card, comment=''):
        sid = integer(card, 1, 'sid')
        excite_id = integer(card, 2, 'excite_id')
        delay = integer_or_blank(card, 3, 'delay', 0)
        Type = integer_string_or_blank(card, 4, 'Type', 'LOAD')

        T1 = double_or_blank(card, 5, 'T1', 0.0)
        T2 = double_or_blank(card, 6, 'T2', T1)
        frequency = double_or_blank(card, 7, 'frequency', 0.)
        phase = double_or_blank(card, 8, 'phase', 0.)
        c = double_or_blank(card, 9, 'c', 0.)
        b = double_or_blank(card, 10, 'b', 0.)
        us0 = double_or_blank(card, 11, 'us0', 0.)
        vs0 = double_or_blank(card, 12, 'vs0', 0.)

        if Type in [0, 'L', 'LO', 'LOA', 'LOAD']:
            Type = 'LOAD'
        elif Type in [1, 'D', 'DI', 'DIS', 'DISP']:
            Type = 'DISP'
        elif Type in [2, 'V', 'VE', 'VEL', 'VELO']:
            Type = 'VELO'
        elif Type in [3, 'A', 'AC', 'ACC', 'ACCE']:
            Type = 'ACCE'
        elif Type in [5, 6, 7, 12, 13]:
            pass
        else:
            msg = 'invalid TLOAD2 type  Type=%r' % Type
            raise RuntimeError(msg)
        assert len(card) <= 13, 'len(TLOAD2 card) = %i\ncard=%s' % (len(card), card)
        return TLOAD2(sid, excite_id, delay, Type, T1, T2, frequency, phase,
                      c, b, us0, vs0, comment=comment)

    def get_load_at_time(self, time, scale=1.):
        if isinstance(time, float):
            time = np.array([time])
        else:
            time = np.asarray(time)

        # A = 1. # points to DAREA or SPCD
        #xy = array(self.tid.table.table)
        #x = xy[:, 0]
        #y = xy[:, 1]
        #assert x.shape == y.shape, 'x.shape=%s y.shape=%s' % (str(x.shape), str(y.shape))
        #f = interp1d(x, y)

        if isinstance(self.delay, float):
            tau = self.delay
        elif self.delay == 0 or self.delay is None:
            tau = 0.0
        else:
            tau = self.delay_ref.get_delay_at_time(time)

        t1 = self.T1 + tau
        t2 = self.T2 + tau
        f = self.frequency
        p = self.phase
        f = np.zeros(time.shape, dtype=time.dtype)

        i = np.where(t1 <= time)[0]
        j = np.where(time[i] <= t2)[0]
        i = i[j]
        f[i] = scale * time[i] ** self.b * np.exp(self.c * time[i]) * np.cos(2 * np.pi * f * time[i] + p)

        is_spcd = False
        resp = f
        if self.Type == 'VELO' and is_spcd:
            f[0] = self.us0
        if self.Type == 'ACCE' and is_spcd:
            f[0] = self.vs0
        return f

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
        msg = ' which is required by TLOAD2 sid=%s' % (self.sid)
        if isinstance(self.delay, integer_types) and self.delay > 0:
            self.delay = model.DELAY(self.delay_id, msg=msg)
            self.delay_ref = self.delay
        # TODO: excite_id

    def safe_cross_reference(self, model, debug=True):
        msg = ' which is required by TLOAD2 sid=%s' % (self.sid)
        if isinstance(self.delay, integer_types) and self.delay > 0:
            self.delay = model.DELAY(self.delay_id, msg=msg)
            self.delay_ref = self.delay
        # TODO: excite_id

    def uncross_reference(self):
        self.delay = self.delay_id
        if isinstance(self.delay, integer_types) and self.delay > 0:
            del self.delay_ref

    @property
    def delay_id(self):
        if self.delay == 0:
            return 0
        elif isinstance(self.delay, integer_types):
            return self.delay
        return self.delay_ref.sid

    def raw_fields(self):
        list_fields = ['TLOAD2', self.sid, self.excite_id, self.delay_id, self.Type,
                       self.T1, self.T2, self.frequency, self.phase, self.c, self.b,
                       self.us0, self.vs0]
        return list_fields

    def repr_fields(self):
        frequency = set_blank_if_default(self.frequency, 0.0)
        phase = set_blank_if_default(self.phase, 0.0)
        c = set_blank_if_default(self.c, 0.0)
        b = set_blank_if_default(self.b, 0.0)

        us0 = set_blank_if_default(self.us0, 0.0)
        vs0 = set_blank_if_default(self.vs0, 0.0)
        list_fields = ['TLOAD2', self.sid, self.excite_id, self.delay_id, self.Type,
                       self.T1, self.T2, frequency, phase, c, b, us0, vs0]
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        if is_double:
            return self.comment + print_card_double(card)
        return self.comment + print_card_16(card)
