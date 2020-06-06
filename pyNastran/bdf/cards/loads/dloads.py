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
from __future__ import annotations
from typing import TYPE_CHECKING
import numpy as np

from pyNastran.utils.numpy_utils import integer_types
from pyNastran.bdf.field_writer_8 import set_blank_if_default
from pyNastran.bdf.bdf_interface.assign_type import (
    integer, double_or_blank, integer_string_or_blank,
    integer_double_or_blank, double)
from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.field_writer_16 import print_card_16
from pyNastran.bdf.field_writer_double import print_card_double
from pyNastran.bdf.cards.loads.loads import DynamicLoad, LoadCombination, BaseCard
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf import BDF


class ACSRCE(BaseCard):
    r"""
    Defines acoustic source as a function of power vs. frequency.

    +--------+-----+----------+---------------+-----------------+-------+-----+---+
    |   1    |  2  |    3     |       4       |        5        |   6   |  7  | 8 |
    +========+=====+==========+===============+=================+=======+=====+===+
    | ACSRCE | SID | EXCITEID | DELAYI/DELAYR | DPHASEI/DPHASER | TP/RP | RHO | B |
    +--------+-----+----------+---------------+-----------------+-------+-----+---+

    ..math ::
      C = \sqrt(B ⁄ ρ)
      Source Strength = {A} * 1/(2πf)  * \sqrt( 8πC P(f) / ρ) ^ (ei(θ + 2πfτ))

    """
    type = 'ACSRCE'

    @classmethod
    def _init_from_empty(cls):
        sid = 1
        excite_id = 2
        rho = 3.
        b = 5.
        return ACSRCE(sid, excite_id, rho, b,
                      delay=0, dphase=0, power=0, comment='')

    def __init__(self, sid, excite_id, rho, b,
                 delay=0, dphase=0, power=0, comment=''):
        """
        Creates an ACSRCE card

        Parameters
        ----------
        sid : int
            load set id number (referenced by DLOAD)
        excite_id : int
            Identification number of a DAREA or SLOAD entry that lists
            each degree of freedom to apply the excitation and the
            corresponding scale factor, A, for the excitation
        rho : float
            Density of the fluid
        b : float
            Bulk modulus of the fluid
        delay : int; default=0
            Time delay, τ.
        dphase : int / float; default=0
            the dphase; if it's 0/blank there is no phase lag
            float : delay in units of time
            int : delay id
        power : int; default=0
            Power as a function of frequency, P(f).
            float : value of P(f) used over all frequencies for all
                    degrees of freedom in EXCITEID entry.
            int : TABLEDi entry that defines P(f) for all degrees of
                  freedom in EXCITEID entry.
        comment : str; default=''
            a comment for the card
        """
        if comment:
            self.comment = comment
        self.sid = sid
        self.excite_id = excite_id
        self.delay = delay
        self.dphase = dphase
        self.power = power
        self.rho = rho
        self.b = b
        self.power_ref = None
        self.sloads_ref = None
        self.delay_ref = None
        self.dphase_ref = None
        #self.dphases_ref = None
        #self.delays_ref = None

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a ACSRCE card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        sid = integer(card, 1, 'sid')
        excite_id = integer(card, 2, 'excite_id') # DAREA, FBALOAD, SLOAD
        delay = integer_double_or_blank(card, 3, 'delay', 0) # DELAY, FBADLAY
        dphase = integer_double_or_blank(card, 4, 'dphase', 0) # DPHASE, FBAPHAS
        power = integer_double_or_blank(card, 5, 'power/tp/rp', 0) # TABLEDi/power
        rho = double(card, 6, 'rho')
        b = double(card, 7, 'bulk modulus')

        assert len(card) <= 8, 'len(ACSRCE card) = %i\n%s' % (len(card), card)
        return ACSRCE(sid, excite_id, rho, b,
                      delay=delay, dphase=dphase, power=power, comment=comment)

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        cmsg = ', which is required by ACSRCE=%s' % (self.sid)

        # TODO: excite_id = DAREA, FBALOAD, SLOAD
        sloads_ref = {}
        lseqs_ref = {}
        for load_id, loads in model.loads.items():
            for load in loads:
                if load.type == 'SLOAD':
                    #if load_id not in sloads_ref:
                        #sloads_ref[load_id] = []
                    for nid in load.node_ids:
                        sloads_ref[(load_id, nid, 0)] = load
                elif load.type == 'LSEQ':
                    load_idi = load.lid_ref[0].sid
                    #print(load)
                    #print(load.lid)
                    excite_idi = load.excite_id
                    #print('load_idi  = %s' % load_idi)
                    #print('excite_id = %s' % excite_idi)
                    assert load_idi not in lseqs_ref
                    lseqs_ref[load_idi] = load
        if sloads_ref:
            self.sloads_ref = sloads_ref
            sload_keys = list(sloads_ref.keys())
            #print('sload_keys =', sload_keys)
        else:
            sload_keys = []

        if self.excite_id not in model.dareas and self.excite_id not in lseqs_ref:
            darea_keys = list(model.dareas.keys())
            dphase_keys = list(model.dphases.keys())
            delay_keys = list(model.delays.keys())
            msg = 'excite_id=%s delay=%s dphase=%s\n' % (
                self.excite_id, self.delay, self.dphase)
            msg += '  darea_keys=%s\n' % darea_keys
            msg += '  sloads(load_id, nid, comp)=%s\n' % sload_keys
            msg += '  dphases(sid)=%s\n' % dphase_keys
            msg += '  delays(delay_id)=%s\n' % delay_keys
            #raise RuntimeError(msg)
            #print(msg)

        if isinstance(self.delay, integer_types) and self.delay > 0:
            delays_ref = {}
            for sload_key in sload_keys:
                nid = sload_key[1]
                delay_key = (self.delay, nid, 0)
                delays_ref[sload_key] = model.DELAY(self.delay, msg=cmsg)
            if delays_ref:
                self.delay_ref = delays_ref

        if isinstance(self.dphase, integer_types) and self.dphase > 0:
            dphases_ref = {}
            for sload_key in sload_keys:
                nid = sload_key[1]
                dphase_key = (self.dphase, nid, 0)
                dphases_ref[sload_key] = model.DPHASE(self.dphase, msg=cmsg)
            if dphases_ref:
                self.dphase_ref = dphases_ref

        if isinstance(self.power, integer_types) and self.power > 0:
            self.power_ref = model.TableD(self.power, msg=cmsg)

        #load_ids2 = []
        #for load_id in self.load_ids:
            #load_id2 = model.DLoad(load_id, consider_dload_combinations=False, msg=msg)
            #load_ids2.append(load_id2)
        #self.load_ids = load_ids2
        #self.load_ids_ref = self.load_ids

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.power = self.Power()
        self.dphase = self.DPhase()
        self.delay = self.Delay()
        #self.sloads = self.

        #self.tb = self.Tb()
        #self.tp = self.Tp()
        #self.delay = self.delay_id
        #if self.tb > 0:
            #del self.tb_ref
        #if self.tp > 0:
            #del self.tp_ref
        self.power_ref = None
        self.sloads_ref = None
        self.delay_ref = None
        self.dphase_ref = None
        #self.dphases_ref = None
        #self.delays_ref = None

    def safe_cross_reference(self, model, xref_errors):
        return self.cross_reference(model)

    #def uncross_reference(self) -> None:
        #self.load_ids = [self.LoadID(load) for load in self.load_ids]
        #del self.load_ids_ref

    def Delay(self):
        if self.delay_ref is not None:
            return next(self.delay_ref.values()).sid
        elif self.delay in [0, 0.0]:
            return 0
        else:
            return self.delay

    def DPhase(self):
        if self.dphase_ref is not None:
            return next(self.delay_ref.values()).tid
        elif self.dphase in [0, 0.0]:
            return 0
        else:
            return self.dphase

    def Power(self):
        if self.power_ref is not None:
            return self.power_ref.tid
        return self.power

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
            #print('delay\n', self.delay_ref)
            tau = self.delay_ref.value
        Pf = self.power_ref.interpolate(freq)
        if self.dphase in [0, 0.]:
            theta = 0.
        else:
            #print('dphase\n', self.dphase_ref)
            theta = self.dphase_ref.interpolate(freq)
        strength = A / (2.* pi * freq) * np.sqrt(8*pi*C*Pf / self.rho) ** (ei*(theta + 2*pi*freq*tau))

        return 0.0

    def raw_fields(self):
        list_fields = ['ACSRCE', self.sid, self.excite_id, self.Delay(), self.DPhase(),
                       self.Power(), self.rho, self.b]
        return list_fields

    def repr_fields(self):
        return self.raw_fields()

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.raw_fields()
        if size == 16:
            return self.comment + print_card_16(card)
        return self.comment + print_card_8(card)


class DLOAD(LoadCombination):
    """
    +-------+-----+----+------+----+----+----+----+----+
    |   1   |  2  |  3 |   4  |  5 |  6 |  7 |  8 |  9 |
    +=======+=====+====+======+====+====+====+====+====+
    | DLOAD | SID |  S |  S1  | L1 | S2 | L2 | S3 | L3 |
    +-------+-----+----+------+----+----+----+----+----+
    |       | S4  | L4 | etc. |    |    |    |    |    |
    +-------+-----+----+------+----+----+----+----+----+
    """
    type = 'DLOAD'

    @classmethod
    def _init_from_empty(cls):
        sid = 1
        scale = 1.
        scale_factors = [1., 2.]
        load_ids = [1, 2]
        return DLOAD(sid, scale, scale_factors, load_ids, comment='')

    def __init__(self, sid, scale, scale_factors, load_ids, comment=''):
        """
        Creates a DLOAD card

        Parameters
        ----------
        sid : int
            Load set identification number. See Remarks 1. and 4. (Integer > 0)
        scale : float
            Scale factor. See Remarks 2. and 8. (Real)
        Si : List[float]
            Scale factors. See Remarks 2., 7. and 8. (Real)
        load_ids : List[int]
            Load set identification numbers of RLOAD1, RLOAD2, TLOAD1,
            TLOAD2, and ACSRCE entries. See Remarks 3 and 7. (Integer > 0)
        comment : str; default=''
            a comment for the card

        """
        LoadCombination.__init__(self, sid, scale, scale_factors, load_ids,
                                 comment=comment)

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        dload_ids2 = []
        msg = ', which is required by DLOAD=%s' % (self.sid)
        for dload_id in self.load_ids:
            dload_id2 = model.DLoad(dload_id, consider_dload_combinations=False, msg=msg)
            dload_ids2.append(dload_id2)
        self.load_ids_ref = dload_ids2

    def safe_cross_reference(self, model, xref_errors, debug=True):
        dload_ids2 = []
        msg = ', which is required by DLOAD=%s' % (self.sid)
        for dload_id in self.load_ids:
            try:
                dload_id2 = model.DLoad(dload_id, consider_dload_combinations=False, msg=msg)
            except KeyError:
                if debug:
                    msg = 'Couldnt find dload_id=%i, which is required by %s=%s' % (
                        dload_id, self.type, self.sid)
                    model.log.warning(msg)
                continue
            dload_ids2.append(dload_id2)
        self.load_ids_ref = dload_ids2

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.load_ids = [self.LoadID(dload) for dload in self.get_load_ids()]
        self.load_ids_ref = None

    def raw_fields(self):
        list_fields = ['DLOAD', self.sid, self.scale]
        for (scale_factor, load_id) in zip(self.scale_factors, self.get_load_ids()):
            list_fields += [scale_factor, self.LoadID(load_id)]
        return list_fields

    def repr_fields(self):
        return self.raw_fields()

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.raw_fields()
        if size == 16:
            return self.comment + print_card_16(card)
        return self.comment + print_card_8(card)


class RLOAD1(DynamicLoad):
    r"""
    Defines a frequency-dependent dynamic load of the form
    for use in frequency response problems.

    .. math::
      \left\{ P(f)  \right\}  = \left\{A\right\} [ C(f)+iD(f)]
         e^{  i \left\{\theta - 2 \pi f \tau \right\} }

    +--------+-----+----------+-------+--------+----+----+------+
    |   1    |  2  |     3    |   4   |   5    |  6 |  7 |   8  |
    +========+=====+==========+=======+========+====+====+======+
    | RLOAD1 | SID | EXCITEID | DELAY | DPHASE | TC | TD | TYPE |
    +--------+-----+----------+-------+--------+----+----+------+
    | RLOAD1 |  5  |    3     |       |        | 1  |    |      |
    +--------+-----+----------+-------+--------+----+----+------+

    NX allows DELAY and DPHASE to be floats
    """
    type = 'RLOAD1'
    _properties = ['delay_id', 'dphase_id']

    @classmethod
    def _init_from_empty(cls):
        sid = 1
        excite_id = 1
        return RLOAD1(sid, excite_id, delay=0, dphase=0, tc=0, td=0, Type='LOAD', comment='')

    def __init__(self, sid, excite_id, delay=0, dphase=0, tc=0, td=0, Type='LOAD', comment=''):
        """
        Creates an RLOAD1 card, which defienes a frequency-dependent load
        based on TABLEDs.

        Parameters
        ----------
        sid : int
            load id
        excite_id : int
            node id where the load is applied
        delay : int/float; default=None
            the delay; if it's 0/blank there is no delay
            float : delay in units of time
            int : delay id
        dphase : int/float; default=None
            the dphase; if it's 0/blank there is no phase lag
            float : delay in units of time
            int : delay id
        tc : int/float; default=0
            TABLEDi id that defines C(f) for all degrees of freedom in
            EXCITEID entry
        td : int/float; default=0
            TABLEDi id that defines D(f) for all degrees of freedom in
            EXCITEID entry
        Type : int/str; default='LOAD'
            the type of load
            0/LOAD
            1/DISP
            2/VELO
            3/ACCE
            4, 5, 6, 7, 12, 13 - MSC only
        comment : str; default=''
            a comment for the card

        """
        DynamicLoad.__init__(self)
        if comment:
            self.comment = comment
        Type = update_loadtype(Type)

        self.sid = sid
        self.excite_id = excite_id
        self.delay = delay
        self.dphase = dphase
        self.tc = tc
        self.td = td
        self.Type = Type
        assert sid > 0, self
        self.tc_ref = None
        self.td_ref = None
        self.delay_ref = None
        self.dphase_ref = None

    def validate(self):
        msg = ''
        is_failed = False
        if self.tc > 0 or self.td > 0:
            msg += 'either RLOAD1 TC or TD > 0; tc=%s td=%s\n' % (self.tc, self.td)

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
        """
        Adds a RLOAD1 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        sid = integer(card, 1, 'sid')
        excite_id = integer(card, 2, 'excite_id')
        delay = integer_double_or_blank(card, 3, 'delay', 0)
        dphase = integer_double_or_blank(card, 4, 'dphase', 0)
        tc = integer_double_or_blank(card, 5, 'tc', 0)
        td = integer_double_or_blank(card, 6, 'td', 0)
        Type = integer_string_or_blank(card, 7, 'Type', 'LOAD')

        assert len(card) <= 8, 'len(RLOAD1 card) = %i\ncard=%s' % (len(card), card)
        return RLOAD1(sid, excite_id, delay, dphase, tc, td, Type, comment=comment)

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ', which is required by RLOAD1 sid=%s' % (self.sid)
        _cross_reference_excite_id(self, model, msg)
        if isinstance(self.tc, integer_types) and self.tc:
            self.tc_ref = model.TableD(self.tc, msg=msg)
        if isinstance(self.td, integer_types) and self.td:
            self.td_ref = model.TableD(self.td, msg=msg)
        if isinstance(self.delay, integer_types) and self.delay > 0:
            self.delay_ref = model.DELAY(self.delay_id, msg=msg)
        if isinstance(self.dphase, integer_types) and self.dphase > 0:
            self.dphase_ref = model.DPHASE(self.dphase, msg=msg)

    def safe_cross_reference(self, model, xref_errors, ):
        msg = ', which is required by RLOAD1 sid=%s' % (self.sid)
        _cross_reference_excite_id(self, model, msg)
        if isinstance(self.tc, integer_types) and self.tc:
            self.tc_ref = model.TableD(self.tc, msg=msg)
        if isinstance(self.td, integer_types) and self.td:
            self.td_ref = model.TableD(self.td, msg=msg)
        if isinstance(self.delay, integer_types) and self.delay > 0:
            self.delay_ref = model.DELAY(self.delay_id, msg=msg)
        if isinstance(self.dphase, integer_types) and self.dphase > 0:
            self.dphase_ref = model.DPHASE(self.dphase, msg=msg)

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.tc = self.Tc()
        self.td = self.Td()
        self.delay = self.delay_id
        self.dphase = self.dphase_id
        self.tc_ref = None
        self.td_ref = None
        self.delay_ref = None
        self.dphase_ref = None

    def get_loads(self):
        return [self]

    def Tc(self):
        if self.tc_ref is not None:
            return self.tc_ref.tid
        elif self.tc in [0, 0.0]:
            return 0
        return self.tc

    def Td(self):
        if self.td_ref is not None:
            return self.td_ref.tid
        elif self.td in [0, 0.0]:
            return 0
        return self.td

    @property
    def delay_id(self):
        if self.delay_ref is not None:
            return self.delay_ref.sid
        elif self.delay in [0, 0.]:
            return 0
        return self.delay

    @property
    def dphase_id(self):
        if self.dphase_ref is not None:
            return self.dphase_ref.sid
        elif self.dphase in [0, 0.0]:
            return 0
        return self.dphase

    def get_load_at_freq(self, freq, scale=1.):
        # A = 1. # points to DAREA or SPCD
        if isinstance(freq, float):
            freq = np.array([freq])
        else:
            freq = np.asarray(freq)

        if isinstance(self.tc, float):
            c = float(self.tc)
        elif self.tc == 0:
            c = 0.
        else:
            c = self.tc_ref.interpolate(freq)

        if isinstance(self.td, float):
            d = float(self.td)
        elif self.td == 0:
            d = 0.
        else:
            d = self.td_ref.interpolate(freq)

        if isinstance(self.dphase, float):
            dphase = self.dphase
        elif self.dphase == 0:
            dphase = 0.0
        else:
            nids, comps, dphases = self.dphase_ref.get_dphase_at_freq(freq)
            assert len(dphases) == 1, 'dphases=%s\n%s' % (dphases, self.dphase_ref)
            dphase = dphases[0]

        if isinstance(self.delay, float):
            tau = self.delay
        elif self.delay == 0:
            tau = 0.0
        else:
            nids, comps, taus = self.delay_ref.get_delay_at_freq(freq)
            assert len(taus) == 1, 'taus=%s\n%s' % (taus, self.delay_ref)
            tau = taus[0]

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

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        if is_double:
            return self.comment + print_card_double(card)
        return self.comment + print_card_16(card)


def _cross_reference_excite_id_backup(self, model, msg):  # pragma: no cover
    """not quite done...not sure how to handle the very odd xref

    EXCITEID may refer to one or more static load entries (FORCE, PLOADi, GRAV, etc.).
    """
    excite_id_ref = []
    case_control = model.case_control_deck
    if case_control is not None:
        #print('cc = %r' % case_control)
        for key, subcase in sorted(model.case_control_deck.subcases.items()):
            #print(subcase, type(subcase))
            #if 'LOADSET' in subcase:
                #lseq_id = subcase['LOADSET'][0]
                #lseq = model.Load(lseq_id, consider_load_combinations=False, msg=msg)[0]
                #self.excite_id_ref = lseq
                ##self.dload_id = lseq.
            #if 'DLOAD' in subcase:
                if self.excite_id in model.loads:
                    # FORCE, FORCE1, FORCE2, PLOAD4, GRAV
                    # changes the magnitudes of the load, not the direction
                    model.log.debug('excite_id load = %s' % self.excite_id)
                    #print('  dloads =', list(model.dloads.keys()))
                    #print('  dareas =', list(model.dareas.keys()))
                    excite_id_ref += model.loads[self.excite_id]
                if self.excite_id in model.dareas:
                    model.log.debug('excite_id darea = %s' % self.excite_id)
                    darea_ref = model.DAREA(self.excite_id, msg=msg)
                    excite_id_ref.append(darea_ref)
                if self.excite_id in model.dload_entries:
                    # this is probably wrong...
                    # it was added to pass TestLoads.test_loads_nonlinear_thermal1, but
                    # I think QVECT should be in self.loads, not self.dload_entries...
                    model.log.debug('excite_id dload_entries = %s' % self.excite_id)
                    excite_id_ref += model.dload_entries
                #  what about TEMPBC?
            #else:
                #msg = ('LOADSET and DLOAD are not found in the case control deck\n%s' %
                       #str(model.case_control_deck))
                #raise RuntimeError(msg)
    #else:
        #model.log.warning('could not find excite_id=%i for\n%s' % (self.excite_id, str(self)))
        #self.excite_id_ref = model.DAREA(self.excite_id, msg=msg)
    if len(excite_id_ref) == 0:
        print('excite_id = %s' % self.excite_id)
        print('  loads  =', list(model.loads.keys()))
        print('  dareas =', list(model.dareas.keys()))
        print('  dloads =', list(model.dloads.keys()))
        print('  dload_entries =', list(model.dload_entries.keys()))
        model.log.warning('could not find excite_id=%i for\n%s' % (self.excite_id, str(self)))
        raise RuntimeError('could not find excite_id=%i for\n%s' % (self.excite_id, str(self)))

def get_lseqs_by_excite_id(model, excite_id):
    from collections import defaultdict

    # get the lseqs that correspond to the correct EXCITE_ID id
    lseq_sids = defaultdict(list)
    for sid, loads in model.load_combinations.items():
        for load in loads:
            if load.type == 'LSEQ':
                if excite_id == load.excite_id:
                    #print(load)
                    lseq_sids[sid].append(load)
    #for sid, loads in lseqs.items():
        #print(sid, loads)
    return lseq_sids

def _cross_reference_excite_id(self, model, msg):
    """not quite done...not sure how to handle the very odd xref

    EXCITEID may refer to one or more static load entries (FORCE, PLOADi, GRAV, etc.).
    """
    #print('*' * 80)
    lseq_sids = get_lseqs_by_excite_id(model, self.excite_id)

    # find all the LOADSETs in the model
    # LOADSETs reference LSEQs by sid
    valid_lseqs = []
    if lseq_sids:
        # get the sid for the LSEQ
        case_control = model.case_control_deck
        if case_control is not None:
            #print('cc = %r' % case_control)
            for key, subcase in sorted(model.case_control_deck.subcases.items()):
                if 'LOADSET' in subcase:
                    lseq_sid = subcase['LOADSET'][0]
                    if lseq_sid in lseq_sids:
                        model.log.debug('adding LOADSET = %i' % lseq_sid)
                        valid_lseqs.append(lseq_sid)
        if valid_lseqs:
            valid_lseqs = list(set(valid_lseqs))
            valid_lseqs.sort()
            #assert len(valid_lseqs) == 1, 'valid_lseqs=%s' % valid_lseqs
    #print('valid_lseqs =', valid_lseqs)
    #  can Case Control LOADSET be substituded for Case Control DLOAD id?

    excite_id_ref = []
    if self.excite_id in model.loads:
        # FORCE, FORCE1, FORCE2, PLOAD4, GRAV
        # changes the magnitudes of the load, not the direction
        model.log.debug('excite_id load = %s' % self.excite_id)
        #print('  dloads =', list(model.dloads.keys()))
        #print('  dareas =', list(model.dareas.keys()))
        excite_id_ref += model.loads[self.excite_id]

    if self.excite_id in model.dareas:
        model.log.debug('excite_id darea = %s' % self.excite_id)
        darea_ref = model.DAREA(self.excite_id, msg=msg)
        excite_id_ref.append(darea_ref)

    if self.excite_id in model.bcs:
        # CONV, TEMPBC
        model.log.debug('excite_id bcs = %s' % self.excite_id)
        excite_id_ref = model.bcs[self.excite_id]

    if self.excite_id in model.dload_entries:  #  this is probably wrong...
        # this is probably wrong...
        # it was added to pass TestLoads.test_loads_nonlinear_thermal1, but
        # I think QVECT should be in self.loads, not self.dload_entries...
        model.log.debug('excite_id dload_entries = %s' % self.excite_id)
        excite_id_ref += model.dload_entries

    if self.excite_id in model.load_combinations:  #  this should be right...
        # C:\NASA\m4\formats\git\examples\move_tpl\nlstrs2.op2
        model.log.debug('excite_id load_combinations = %s' % self.excite_id)
        excite_id_ref = model.load_combinations[self.excite_id]

    #  handles LSEQ
    if valid_lseqs:
        for lseq_sid in valid_lseqs:
            excite_id_ref += lseq_sids[lseq_sid]

    #  what about SPCD?

    if len(excite_id_ref) == 0:
        print(model.get_bdf_stats())
        print('excite_id = %s' % self.excite_id)
        print('  loads  =', list(model.loads.keys()))
        print('  dareas =', list(model.dareas.keys()))
        print('  bcs =', list(model.bcs.keys()))
        print('  dloads =', list(model.dloads.keys()))
        print('  dload_entries =', list(model.dload_entries.keys()))
        print('  load_combinations =', list(model.load_combinations.keys())) #  what about LSEQ
        if lseq_sids:
            sids = list(lseq_sids.keys())
            print('  lseq_excite_ids=%s; lseq_sids=%s; valid_lseqs=%s' % (
                self.excite_id, sids, valid_lseqs))
        else:
            print('  lseq_sids = []')
        model.log.warning('could not find excite_id=%i for\n%s' % (self.excite_id, str(self)))
        raise RuntimeError('could not find excite_id=%i for\n%s' % (self.excite_id, str(self)))


class RLOAD2(DynamicLoad):
    r"""
    Defines a frequency-dependent dynamic load of the form
    for use in frequency response problems.

    .. math:: \left\{ P(f)  \right\}  = \left\{A\right\} * B(f)
        e^{  i \left\{ \phi(f) + \theta - 2 \pi f \tau \right\} }

    +--------+-----+----------+-------+--------+----+----+------+
    |   1    |  2  |     3    |   4   |    5   |  6 |  7 |  8   |
    +========+=====+==========+=======+========+====+====+======+
    | RLOAD2 | SID | EXCITEID | DELAY | DPHASE | TB | TP | TYPE |
    +--------+-----+----------+-------+--------+----+----+------+
    | RLOAD2 |  5  |    3     |       |        | 1  |    |      |
    +--------+-----+----------+-------+--------+----+----+------+

    NX allows DELAY and DPHASE to be floats
    """
    type = 'RLOAD2'
    _properties = ['delay_id', 'dphase_id']

    @classmethod
    def _init_from_empty(cls):
        sid = 1
        excite_id = 1
        return RLOAD2(sid, excite_id, delay=0, dphase=0, tb=0, tp=0, Type='LOAD', comment='')

    # P(f) = {A} * B(f) * e^(i*phi(f), + theta - 2*pi*f*tau)
    def __init__(self, sid, excite_id, delay=0, dphase=0, tb=0, tp=0, Type='LOAD', comment=''):
        """
        Creates a nRLOAD2 card, which defienes a frequency-dependent load
        based on TABLEDs.

        Parameters
        ----------
        sid : int
            load id
        excite_id : int
            node id where the load is applied
        delay : int/float; default=None
            the delay; if it's 0/blank there is no delay
            float : delay in units of time
            int : delay id
        dphase : int/float; default=None
            the dphase; if it's 0/blank there is no phase lag
            float : delay in units of time
            int : delay id
        tb : int/float; default=0
            TABLEDi id that defines B(f) for all degrees of freedom in
            EXCITEID entry
        tc : int/float; default=0
            TABLEDi id that defines C(f) for all degrees of freedom in
            EXCITEID entry
        td : int/float; default=0
            TABLEDi id that defines D(f) for all degrees of freedom in
            EXCITEID entry
        tp : int/float; default=0
            TABLEDi id that defines phi(f) for all degrees of freedom in
            EXCITEID entry
        Type : int/str; default='LOAD'
            the type of load
            0/LOAD
            1/DISP
            2/VELO
            3/ACCE
            4, 5, 6, 7, 12, 13 - MSC only
        comment : str; default=''
            a comment for the card

        """
        DynamicLoad.__init__(self)
        if comment:
            self.comment = comment
        Type = update_loadtype(Type)

        self.sid = sid
        self.excite_id = excite_id
        self.delay = delay
        self.dphase = dphase
        self.tb = tb
        self.tp = tp
        self.Type = Type
        self.tb_ref = None
        self.tp_ref = None
        self.delay_ref = None
        self.dphase_ref = None

    #@property
    #def Type(self):
        #"""gets the load_type"""
        #return self.load_type
    #@Type.setter
    #def Type(self, load_type):
        #"""sets the load_type"""
        #self.load_type = load_type

    def validate(self):
        msg = ''
        is_failed = False
        if self.tb > 0 or self.tp > 0:
            msg += 'either RLOAD2 TB or TP > 0; tb=%s tp=%s\n' % (self.tb, self.tp)

        if self.Type in [0, 'L', 'LO', 'LOA', 'LOAD']:
            self.Type = 'LOAD'
        elif self.Type in [1, 'D', 'DI', 'DIS', 'DISP']:
            self.Type = 'DISP'
        elif self.Type in [2, 'V', 'VE', 'VEL', 'VELO']:
            self.Type = 'VELO'
        elif self.Type in [3, 'A', 'AC', 'ACC', 'ACCE']:
            self.Type = 'ACCE'
        else:
            msg += 'invalid RLOAD2 type  Type=%r\n' % self.Type
            is_failed = True

        if is_failed:
            msg += str(self)
            raise RuntimeError(msg)

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a RLOAD2 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        sid = integer(card, 1, 'sid')
        excite_id = integer(card, 2, 'excite_id')
        delay = integer_double_or_blank(card, 3, 'delay', 0)
        dphase = integer_double_or_blank(card, 4, 'dphase', 0)
        tb = integer_double_or_blank(card, 5, 'tb', 0)
        tp = integer_double_or_blank(card, 6, 'tp', 0)
        Type = integer_string_or_blank(card, 7, 'Type', 'LOAD')

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
            nids, comps, dphases = self.dphase_ref.get_dphase_at_freq(freq)
            assert len(dphases) == 1, dphases
            dphase = dphases[0]

        if isinstance(self.delay, float):
            tau = self.delay
        elif self.delay == 0:
            tau = 0.0
        else:
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

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ', which is required by RLOAD2=%s' % (self.sid)
        _cross_reference_excite_id(self, model, msg)
        if isinstance(self.tb, integer_types) and self.tb:
            self.tb_ref = model.TableD(self.tb, msg=msg)
        if isinstance(self.tp, integer_types) and self.tp:
            self.tp_ref = model.TableD(self.tp, msg=msg)
        if isinstance(self.delay, integer_types) and self.delay > 0:
            self.delay_ref = model.DELAY(self.delay, msg=msg)
        if isinstance(self.dphase, integer_types) and self.dphase > 0:
            self.dphase_ref = model.DPHASE(self.dphase, msg=msg)

    def safe_cross_reference(self, model, xref_errors, ):
        msg = ', which is required by RLOAD2=%s' % (self.sid)
        _cross_reference_excite_id(self, model, msg)
        if isinstance(self.tb, integer_types) and self.tb:
            self.tb_ref = model.TableD(self.tb, msg=msg)
        if isinstance(self.tp, integer_types) and self.tp:
            self.tp_ref = model.TableD(self.tp, msg=msg)
        if isinstance(self.delay, integer_types) and self.delay > 0:
            self.delay_ref = model.DELAY(self.delay, msg=msg)
        if isinstance(self.dphase, integer_types) and self.dphase > 0:
            self.dphase_ref = model.DPHASE(self.dphase, msg=msg)

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.tb = self.Tb()
        self.tp = self.Tp()
        self.delay = self.delay_id
        self.dphase = self.dphase_id
        self.tb_ref = None
        self.tp_ref = None
        self.delay_ref = None
        self.dphase_ref = None

    def get_loads(self):
        return [self]

    def LoadID(self):
        return self.sid

    def Tb(self):
        if self.tb_ref is not None:
            return self.tb_ref.tid
        elif self.tb == 0:
            return 0
        return self.tb

    def Tp(self):
        if self.tp_ref is not None:
            return self.tp_ref.tid
        elif self.tp == 0:
            return 0
        return self.tp

    @property
    def delay_id(self):
        if self.delay_ref is not None:
            return self.delay_ref.sid
        elif self.delay == 0:
            return 0
        return self.delay

    @property
    def dphase_id(self):
        if self.dphase_ref is not None:
            return self.dphase_ref.sid
        elif self.dphase == 0:
            return 0
        return self.dphase

    def raw_fields(self):
        list_fields = ['RLOAD2', self.sid, self.excite_id, self.delay_id, self.dphase_id,
                       self.Tb(), self.Tp(), self.Type]
        return list_fields

    def repr_fields(self):
        Type = set_blank_if_default(self.Type, 0.0)
        list_fields = ['RLOAD2', self.sid, self.excite_id, self.delay_id, self.dphase_id,
                       self.Tb(), self.Tp(), Type]
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        if is_double:
            return self.comment + print_card_double(card)
        return self.comment + print_card_16(card)


class TLOAD1(DynamicLoad):
    r"""
    Transient Response Dynamic Excitation, Form 1

    Defines a time-dependent dynamic load or enforced motion of the form:

    .. math::
      \left\{ P(t) \right\} = \left\{ A \right\} \cdot F(t-\tau)

    for use in transient response analysis.

    MSC 20005.2
    +--------+-----+----------+-------+------+-----+-----+-----+
    |    1   |  2  |     3    |   4   |   5  |  6  |  7  |  8  |
    +========+=====+==========+=======+======+=====+=====+=====+
    | TLOAD1 | SID | EXCITEID | DELAY | TYPE | TID | US0 | VS0 |
    +--------+-----+----------+-------+------+-----+-----+-----+

    NX 11
    +--------+-----+----------+-------+------+-----+
    |    1   |  2  |     3    |   4   |   5  |  6  |
    +========+=====+==========+=======+======+=====+
    | TLOAD1 | SID | EXCITEID | DELAY | TYPE | TID |
    +--------+-----+----------+-------+------+-----+
    """
    type = 'TLOAD1'
    _properties = ['delay_id']

    @classmethod
    def _init_from_empty(cls):
        sid = 1
        excite_id = 1
        tid = 1
        return TLOAD1(sid, excite_id, tid, delay=0, Type='LOAD', us0=0.0, vs0=0.0, comment='')

    def __init__(self, sid, excite_id, tid, delay=0, Type='LOAD',
                 us0=0.0, vs0=0.0, comment=''):
        """
        Creates a TLOAD1 card, which defienes a time-dependent load
        based on a DTABLE.

        Parameters
        ----------
        sid : int
            load id
        excite_id : int
            node id where the load is applied
        tid : int
            TABLEDi id that defines F(t) for all degrees of freedom in
            EXCITEID entry
            float : MSC not supported
        delay : int/float; default=None
            the delay; if it's 0/blank there is no delay
            float : delay in units of time
            int : delay id
        Type : int/str; default='LOAD'
            the type of load
            0/LOAD
            1/DISP
            2/VELO
            3/ACCE
            4, 5, 6, 7, 12, 13 - MSC only
        us0 : float; default=0.
            Factor for initial displacements of the enforced degrees-of-freedom
            MSC only
        vs0 : float; default=0.
            Factor for initial velocities of the enforced degrees-of-freedom
            MSC only
        comment : str; default=''
            a comment for the card
        """
        DynamicLoad.__init__(self)
        if delay is None:
            delay = 0
        Type = update_loadtype(Type)

        if comment:
            self.comment = comment

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

        self.tid_ref = None
        self.delay_ref = None

    def validate(self):
        if self.Type in [0, 'L', 'LO', 'LOA', 'LOAD']:
            self.Type = 'LOAD'
        elif self.Type in [1, 'D', 'DI', 'DIS', 'DISP']:
            self.Type = 'DISP'
        elif self.Type in [2, 'V', 'VE', 'VEL', 'VELO']:
            self.Type = 'VELO'
        elif self.Type in [3, 'A', 'AC', 'ACC', 'ACCE']:
            self.Type = 'ACCE'
        elif self.Type in [4, 5, 6, 7, 12, 13]:  # MSC-only
            pass
        else:
            msg = 'invalid TLOAD1 type  Type=%r' % self.Type
            raise RuntimeError(msg)

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a TLOAD1 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        sid = integer(card, 1, 'sid')
        excite_id = integer(card, 2, 'excite_id')
        delay = integer_double_or_blank(card, 3, 'delay', 0)
        Type = integer_string_or_blank(card, 4, 'Type', 'LOAD')
        tid = integer(card, 5, 'tid')
        us0 = double_or_blank(card, 6, 'us0', 0.0)
        vs0 = double_or_blank(card, 7, 'vs0', 0.0)

        assert len(card) <= 8, 'len(TLOAD1 card) = %i\ncard=%s' % (len(card), card)
        return TLOAD1(sid, excite_id, tid, delay=delay, Type=Type, us0=us0, vs0=vs0, comment=comment)

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
        msg = ', which is required by TLOAD1=%s' % (self.sid)
        _cross_reference_excite_id(self, model, msg)
        if self.tid:
            self.tid_ref = model.TableD(self.tid, msg=msg)
        if isinstance(self.delay, integer_types) and self.delay > 0:
            self.delay_ref = model.DELAY(self.delay, msg=msg)

    def safe_cross_reference(self, model, debug=True):
        msg = ', which is required by TLOAD1=%s' % (self.sid)
        _cross_reference_excite_id(self, model, msg)
        if self.tid:
            #try:
            self.tid_ref = model.TableD(self.tid, msg=msg)
            #except
        if isinstance(self.delay, integer_types) and self.delay > 0:
            self.delay_ref = model.DELAY(self.delay_id, msg=msg)

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.tid = self.Tid()
        self.delay = self.delay_id
        self.tid_ref = None
        self.delay_ref = None

    def Tid(self):
        if self.tid_ref is not None:
            return self.tid_ref.tid
        elif self.tid == 0:
            return 0
        else:
            return self.tid

    @property
    def delay_id(self):
        if self.delay_ref is not None:
            return self.delay_ref.sid
        elif self.delay == 0:
            return 0
        return self.delay

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
            tau = self.delay_ref.get_delay_at_time(time)

        i = np.where(time - tau > 0)
        time2 = time[i]
        resp = self.tid_ref.interpolate(time2)
        is_spcd = False
        if self.Type == 'VELO' and is_spcd:
            resp[0] = self.us0
        if self.Type == 'ACCE' and is_spcd:
            resp[0] = self.vs0
        return resp * scale

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

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        if is_double:
            return self.comment + print_card_double(card)
        return self.comment + print_card_16(card)


class TLOAD2(DynamicLoad):
    r"""
    Transient Response Dynamic Excitation, Form 1

    Defines a time-dependent dynamic load or enforced motion of the form:

    .. math::
    \left\{ P(t) \right\} = \left\{ A \right\} e^(C*t) cos(2 \pi f t + \phi)

      P(t) = 0                                            (t<T1+tau or t >  T2+tau)
      P(t) = {A} * t^b * e^(C*t) * cos(2*pi*f*t + phase)  (T1+tau <=   t <= T2+tau)

    for use in transient response analysis.

    MSC 2016.1
    +--------+-----+----------+-------+------+-----+-----+--------+---------+
    |    1   |  2  |     3    |   4   |   5  |  6  |  7  |    8   |    9    |
    +========+=====+==========+=======+======+=====+=====+========+=========+
    | TLOAD2 | SID | EXCITEID | DELAY | TYPE | T1  | T2  |  FREQ  |  PHASE  |
    +--------+-----+----------+-------+------+-----+-----+--------+---------+
    |        |  C  |     B    |  US0  |  VS0 |     |     |        |         |
    +--------+-----+----------+-------+------+-----+-----+--------+---------+

    NX 11
    +--------+-----+----------+-------+------+-----+-----+--------+---------+
    |    1   |  2  |     3    |   4   |   5  |  6  |  7  |    8   |    9    |
    +========+=====+==========+=======+======+=====+=====+========+=========+
    | TLOAD2 | SID | EXCITEID | DELAY | TYPE | T1  | T2  |  FREQ  |  PHASE  |
    +--------+-----+----------+-------+------+-----+-----+--------+---------+
    |        |  C  |     B    |       |      |     |     |        |         |
    +--------+-----+----------+-------+------+-----+-----+--------+---------+

    """
    type = 'TLOAD2'
    _properties = ['delay_id']

    @classmethod
    def _init_from_empty(cls):
        sid = 1
        excite_id = 1
        return TLOAD2(sid, excite_id, delay=0, Type='LOAD', T1=0., T2=None,
                      frequency=0., phase=0., c=0., b=0., us0=0., vs0=0., comment='')

    def __init__(self, sid, excite_id, delay=0, Type='LOAD', T1=0., T2=None,
                 frequency=0., phase=0., c=0., b=0., us0=0., vs0=0., comment=''):
        """
        Creates a TLOAD2 card, which defines a exponential time dependent
        load based on constants.

        Parameters
        ----------
        sid : int
            load id
        excite_id : int
            node id where the load is applied
        delay : int/float; default=None
            the delay; if it's 0/blank there is no delay
            float : delay in units of time
            int : delay id
        Type : int/str; default='LOAD'
            the type of load
            0/LOAD
            1/DISP
            2/VELO
            3/ACCE
            4, 5, 6, 7, 12, 13 - MSC only
        T1 : float; default=0.
            time constant (t1 > 0.0)
            times below this are ignored
        T2 : float; default=None
            time constant (t2 > t1)
            times above this are ignored
        frequency : float; default=0.
            Frequency in cycles per unit time.
        phase : float; default=0.
            Phase angle in degrees.
        c : float; default=0.
            Exponential coefficient.
        b : float; default=0.
            Growth coefficient.
        us0 : float; default=0.
            Factor for initial displacements of the enforced degrees-of-freedom
            MSC only
        vs0 : float; default=0.
            Factor for initial velocities of the enforced degrees-of-freedom
            MSC only
        comment : str; default=''
            a comment for the card
        """
        DynamicLoad.__init__(self)
        if comment:
            self.comment = comment
        if T2 is None:
            T2 = T1

        Type = update_loadtype(Type)

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

        self.delay_ref = None

    def validate(self):
        if self.Type in [0, 'L', 'LO', 'LOA', 'LOAD']:
            self.Type = 'LOAD'
        elif self.Type in [1, 'D', 'DI', 'DIS', 'DISP']:
            self.Type = 'DISP'
        elif self.Type in [2, 'V', 'VE', 'VEL', 'VELO']:
            self.Type = 'VELO'
        elif self.Type in [3, 'A', 'AC', 'ACC', 'ACCE']:
            self.Type = 'ACCE'
        elif self.Type in [5, 6, 7, 12, 13]: # MSC only
            pass
        else:
            msg = 'invalid TLOAD2 type  Type=%r' % self.Type
            raise RuntimeError(msg)

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a TLOAD2 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        sid = integer(card, 1, 'sid')
        excite_id = integer(card, 2, 'excite_id')
        delay = integer_double_or_blank(card, 3, 'delay', 0)
        Type = integer_string_or_blank(card, 4, 'Type', 'LOAD')

        T1 = double_or_blank(card, 5, 'T1', 0.0)
        T2 = double_or_blank(card, 6, 'T2', T1)
        frequency = double_or_blank(card, 7, 'frequency', 0.)
        phase = double_or_blank(card, 8, 'phase', 0.)
        c = double_or_blank(card, 9, 'c', 0.)
        b = double_or_blank(card, 10, 'b', 0.)
        us0 = double_or_blank(card, 11, 'us0', 0.)
        vs0 = double_or_blank(card, 12, 'vs0', 0.)

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
        #resp = f
        if self.Type == 'VELO' and is_spcd:
            f[0] = self.us0
        if self.Type == 'ACCE' and is_spcd:
            f[0] = self.vs0
        return f

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
        msg = ', which is required by TLOAD2 sid=%s' % (self.sid)
        _cross_reference_excite_id(self, model, msg)
        if isinstance(self.delay, integer_types) and self.delay > 0:
            self.delay_ref = model.DELAY(self.delay_id, msg=msg)
        # TODO: excite_id

    def safe_cross_reference(self, model, xref_errors, debug=True):
        msg = ', which is required by TLOAD2 sid=%s' % (self.sid)
        _cross_reference_excite_id(self, model, msg)
        if isinstance(self.delay, integer_types) and self.delay > 0:
            self.delay_ref = model.DELAY(self.delay_id, msg=msg)
        # TODO: excite_id

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.delay = self.delay_id
        self.delay_ref = None

    @property
    def delay_id(self):
        if self.delay_ref is not None:
            return self.delay_ref.sid
        elif self.delay == 0:
            return 0
        return self.delay

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

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        if is_double:
            return self.comment + print_card_double(card)
        return self.comment + print_card_16(card)

def update_loadtype(load_type):
    if load_type in [0, 'L', 'LO', 'LOA', 'LOAD']:
        load_type = 'LOAD'
    elif load_type in [1, 'D', 'DI', 'DIS', 'DISP']:
        load_type = 'DISP'
    elif load_type in [2, 'V', 'VE', 'VEL', 'VELO']:
        load_type = 'VELO'
    elif load_type in [3, 'A', 'AC', 'ACC', 'ACCE']:
        load_type = 'ACCE'
    return load_type
