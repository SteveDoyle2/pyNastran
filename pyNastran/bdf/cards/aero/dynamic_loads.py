# coding: utf-8
"""
All aero cards are defined in this file.  This includes:

 * AERO
 * FLFACT
 * FLUTTER
 * GUST
 * MKAERO1 / MKAERO2

All cards are BaseCard objects.

"""
from __future__ import annotations
from itertools import count
from typing import TYPE_CHECKING
import numpy as np

from pyNastran.utils.numpy_utils import integer_types
from pyNastran.bdf.field_writer_8 import set_blank_if_default, print_card_8
from pyNastran.bdf.field_writer_16 import print_card_16
from pyNastran.bdf.cards.base_card import BaseCard
from pyNastran.utils.atmosphere import (
    make_flfacts_eas_sweep, make_flfacts_alt_sweep, make_flfacts_mach_sweep,
    atm_density, _velocity_factor)
from pyNastran.bdf.bdf_interface.assign_type import (
    integer, integer_or_blank, double, double_or_blank, string,
    fields, string_or_blank, double_string_or_blank, interpret_value)
from pyNastran.bdf.cards.utils import wipe_empty_fields
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf import BDF


class Aero(BaseCard):
    """Base class for AERO and AEROS cards."""
    def __init__(self):
        """
        Common class for AERO, AEROS

        Attributes
        ----------
        acsid : int; default=0
            aerodyanmic coordinate system
            defines the direction of the wind
        sym_xz : int; default=0
            xz symmetry flag (+1=symmetry; -1=antisymmetric)
        sym_xy : int; default=0
            xy symmetry flag (+1=symmetry; -1=antisymmetric)

        """
        BaseCard.__init__(self)
        self.sym_xy = None
        self.sym_xz = None
        self.acsid = None
        self.acsid_ref = None

    def Acsid(self):
        try:
            return self.acsid_ref.cid
        except AttributeError:
            return self.acsid

    @property
    def is_symmetric_xy(self):
        if self.sym_xy == 1:
            return True
        return False

    @property
    def is_symmetric_xz(self):
        if self.sym_xz == 1:
            return True
        return False

    @property
    def is_anti_symmetric_xy(self):
        if self.sym_xy == -1:
            return True
        return False

    @property
    def is_anti_symmetric_xz(self):
        if self.sym_xz == -1:
            return True
        return False

    def set_ground_effect(self, enable):  # TODO: verify
        if enable:
            self.sym_xy = -1
        else:
            self.sym_xy = 1


class AERO(Aero):
    """
    Gives basic aerodynamic parameters for unsteady aerodynamics.

    +------+-------+----------+------+--------+-------+-------+
    |   1  |   2   |    3     |   4  |   5    |   6   |   7   |
    +======+=======+==========+======+========+=======+=======+
    | AERO | ACSID | VELOCITY | REFC | RHOREF | SYMXZ | SYMXY |
    +------+-------+----------+------+--------+-------+-------+
    | AERO |   3   |   1.3+   | 100. |  1.-5  |   1   |  -1   |
    +------+-------+----------+------+--------+-------+-------+
    """
    type = 'AERO'
    _properties = ['is_anti_symmetric_xy', 'is_anti_symmetric_xz',
                   'is_symmetric_xy', 'is_symmetric_xz']
    _field_map = {
        1: 'acsid', 2:'velocity', 3:'cRef', 4:'rhoRef', 5:'symXZ',
        6:'symXY',
    }

    @classmethod
    def _init_from_empty(cls):
        velocity = 1.
        cref = 1.
        rho_ref = 1.
        return AERO(velocity, cref, rho_ref, acsid=0, sym_xz=0, sym_xy=0, comment='')

    def __init__(self, velocity, cref, rho_ref, acsid=0, sym_xz=0, sym_xy=0, comment=''):
        """
        Creates an AERO card

        Parameters
        ----------
        velocity : float
            the airspeed
        cref : float
            the aerodynamic chord
        rho_ref : float
            FLFACT density scaling factor
        acsid : int; default=0
            aerodyanmic coordinate system
            defines the direction of the wind
        sym_xz : int; default=0
            xz symmetry flag (+1=symmetry; -1=antisymmetric)
        sym_xy : int; default=0
            xy symmetry flag (+1=symmetry; -1=antisymmetric)
        comment : str; default=''
            a comment for the card

        """
        Aero.__init__(self)
        if comment:
            self.comment = comment

        #: Aerodynamic coordinate system identification
        if acsid is None:
            acsid = 0
        self.acsid = acsid

        #: Velocity for aerodynamic force data recovery and to calculate the BOV
        #: parameter
        self.velocity = velocity

        #: Reference length for reduced frequency
        self.cref = cref

        #: Reference density
        self.rho_ref = rho_ref

        #: Symmetry key for the aero coordinate x-z plane. See Remark 6.
        #: (Integer = +1 for symmetry, 0 for no symmetry, and -1 for antisymmetry;
        #: Default = 0)
        self.sym_xz = sym_xz

        #: The symmetry key for the aero coordinate x-y plane can be used to
        #: simulate ground effect. (Integer = -1 for symmetry, 0 for no symmetry,
        #: and +1 for antisymmetry; Default = 0)
        self.sym_xy = sym_xy

    def validate(self):
        msg = ''
        if not isinstance(self.acsid, integer_types):
            msg += 'acsid=%r must be an integer; type=%s' % (
                self.acsid, type(self.acsid))
        if not isinstance(self.sym_xz, integer_types):
            msg = 'sym_xz=%r must be an integer; type=%s' % (
                self.sym_xz, type(self.sym_xz))
        if not isinstance(self.sym_xy, integer_types):
            msg = 'sym_xy=%r must be an integer; type=%s' % (
                self.sym_xy, type(self.sym_xy))
        if msg:
            raise TypeError(msg + str(self))

    def cross_reference(self, model: BDF) -> None:
        """
        Cross refernece aerodynamic coordinate system.

        Parameters
        ----------
        model : BDF
            The BDF object.

        """
        msg = ', which is required by AERO'
        self.acsid_ref = model.Coord(self.acsid, msg=msg)

    def safe_cross_reference(self, model, xref_errors):
        """
        Safe cross refernece aerodynamic coordinate system.

        Parameters
        ----------
        model : BDF
            The BDF object.

        """
        msg = ', which is required by AERO'
        self.acsid_ref = model.safe_coord(self.acsid, None, xref_errors, msg=msg)

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds an AERO card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        acsid = integer_or_blank(card, 1, 'acsid', 0)
        velocity = double_or_blank(card, 2, 'velocity')
        cref = double(card, 3, 'cRef')
        rho_ref = double(card, 4, 'rho_ref')
        sym_xz = integer_or_blank(card, 5, 'symXZ', 0)
        sym_xy = integer_or_blank(card, 6, 'symXY', 0)
        assert len(card) <= 7, 'len(AERO card) = %i\ncard=%s' % (len(card), card)
        return AERO(velocity, cref, rho_ref, acsid=acsid, sym_xz=sym_xz, sym_xy=sym_xy,
                    comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        acsid = data[0]
        velocity = data[1]
        cref = data[2]
        rho_ref = data[3]
        sym_xz = data[4]
        sym_xy = data[5]
        assert len(data) == 6, 'data = %s' % data
        return AERO(acsid, velocity, cref, rho_ref, sym_xz, sym_xy,
                    comment=comment)

        # T is the tabular function
        #angle = self.wg*self.t*(t-(x-self.x0)/self.V)

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.acsid_ref = None

    def update(self, maps):
        """
        maps = {
            'coord' : cid_map,
        }

        """
        cid_map = maps['coord']
        self.acsid = cid_map[self.acsid]

    def raw_fields(self):
        """
        Gets the fields in their unmodified form

        Returns
        -------
        fields : List[int/float/str]
           the fields that define the card

        """
        list_fields = ['AERO', self.Acsid(), self.velocity, self.cref,
                       self.rho_ref, self.sym_xz, self.sym_xy]
        return list_fields

    def repr_fields(self):
        """
        Gets the fields in their simplified form

        Returns
        -------
        fields : List[varies]
          the fields that define the card

        """
        sym_xz = set_blank_if_default(self.sym_xz, 0)
        sym_xy = set_blank_if_default(self.sym_xy, 0)
        list_fields = ['AERO', self.Acsid(), self.velocity, self.cref,
                       self.rho_ref, sym_xz, sym_xy]
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


class FLFACT(BaseCard):
    """
    +--------+-----+----+------+-----+----+----+----+----+
    |   1    |  2  |  3 |   4  |  5  | 6  | 7  | 8  | 9  |
    +========+=====+====+======+=====+====+====+====+====+
    | FLFACT | SID | F1 | F2   | F3  | F4 | F5 | F6 | F7 |
    +--------+-----+----+------+-----+----+----+----+----+
    |        | F8  | F9 | etc. |     |    |    |    |    |
    +--------+-----+----+------+-----+----+----+----+----+
    | FLFACT | 97  | .3 |  .7  | 3.5 |    |    |    |    |
    +--------+-----+----+------+-----+----+----+----+----+

    # delta quantity approach

    +--------+-----+-------+------+-------+----+--------+
    |   1    |  2  |  3    |   4  |   5   | 6  |   7    |
    +========+=====+=======+======+=======+====+========+
    | FLFACT | SID | F1    | THRU | FNF   | NF |  FMID  |
    +--------+-----+-------+------+-------+----+--------+
    | FLFACT | 201 | 0.200 | THRU | 0.100 | 11 | 0.1333 |
    +--------+-----+-------+------+-------+----+--------+
    """
    type = 'FLFACT'

    @classmethod
    def _init_from_empty(cls):
        sid = 1
        factors = [1.]
        return FLFACT(sid, factors, comment='')

    def __init__(self, sid, factors, comment=''):
        """
        Creates an FLFACT card, which defines factors used for flutter
        analysis.  These factors define either:
         - density
         - mach
         - velocity
         - reduced frequency
        depending on the FLUTTER method chosen (e.g., PK, PKNL, PKNLS)

        Parameters
        ----------
        sid : int
            the id of a density, reduced_frequency, mach, or velocity table
            the FLUTTER card defines the meaning
        factors : varies
            values : List[float, ..., float]
                list of factors
            List[f1, THRU, fnf, nf, fmid]
                f1 : float
                    first value
                THRU : str
                    the word THRU
                fnf : float
                    second value
                nf : int
                    number of values
                fmid : float; default=(f1 + fnf) / 2.
                    the mid point to bias the array
                TODO: does f1 need be be greater than f2/fnf???
        comment : str; default=''
            a comment for the card

        """
        BaseCard.__init__(self)
        if comment:
            self.comment = comment
        self.sid = sid
        #self.f1 = f1
        #self.fnf = fnf
        #self.nf = nf
        #self.fmid = fmid

        # the dumb string_types thing is because we also get floats
        if len(factors) > 1 and isinstance(factors[1], str) and factors[1] == 'THRU':
            #msg = 'embedded THRUs not supported yet on FLFACT card\n'
            nfactors = len(factors)
            if nfactors == 4:
                (f1, _thru, fnf, nf) = factors
                fmid = (f1 + fnf) / 2.
            elif nfactors == 5:
                (f1, _thru, fnf, nf, fmid) = factors
                #assert _thru.upper() == 'THRU', 'factors=%s' % str(factors)
            else:
                raise RuntimeError('factors must be length 4/5; factors=%s' % factors)
            i = np.linspace(0, nf, nf, endpoint=False) + 1
            factors = (
                (f1*(fnf - fmid) * (nf-i) + fnf * (fmid - f1) * (i-1)) /
                (   (fnf - fmid) * (nf-i) +       (fmid - f1) * (i-1))
            )
        self.factors = np.asarray(factors)

    def validate(self):
        if len(self.factors) == 0:
            raise ValueError('FLFACT sid=%s is empty; factors=%s' % (self.sid, str(self.factors)))

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds an FLFACT card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        sid = integer(card, 1, 'sid')
        assert len(card) > 2, 'len(FLFACT card)=%s; card=%s' % (len(card), card)
        field3 = double_string_or_blank(card, 3, 'THRU')
        if field3 is None:
            f1 = double(card, 2, 'f1')
            factors = [f1]
            assert len(card) == 3, 'len(FLFACT card)=%s; card=%s' % (len(card), card)
        elif isinstance(field3, float):
            factors = fields(double, card, 'factors', i=2, j=len(card))
        elif isinstance(field3, str) and field3 == 'THRU':
            f1 = double(card, 2, 'f1')
            fnf = double(card, 4, 'fnf')
            nf = integer(card, 5, 'nf')
            fmid_default = (f1 + fnf) / 2.
            fmid = double_or_blank(card, 6, 'fmid', fmid_default)
            assert len(card) <= 7, 'len(FLFACT card)=%s; card=%s' % (len(card), card)
            factors = [f1, 'THRU', fnf, nf, fmid]
        else:
            raise SyntaxError('expected a float or string for FLFACT field 3; value=%r' % field3)
        return FLFACT(sid, factors, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        sid = data[0]
        factors = data[1:]
        return FLFACT(sid, factors, comment=comment)

    def max(self):
        return self.factors.max()

    def min(self):
        return self.factors.min()

    #def uncross_reference(self) -> None:
        #pass

    def raw_fields(self):
        """
        Gets the fields in their unmodified form

        Returns
        -------
        fields : list[varies]
            the fields that define the card

        """
        list_fields = ['FLFACT', self.sid] + list(self.factors)
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)


class FLUTTER(BaseCard):
    """
    Defines data needed to perform flutter analysis.

    +---------+-----+--------+------+------+-------+-------+-------------+------+
    |    1    |  2  |   3    |  4   |  5   |   6   |   7   |      8      |  9   |
    +=========+=====+========+======+======+=======+=======+=============+======+
    | FLUTTER | SID | METHOD | DENS | MACH | RFREQ | IMETH | NVALUE/OMAX | EPS  |
    +---------+-----+--------+------+------+-------+-------+-------------+------+
    | FLUTTER | 19  |   K    | 119  | 219  | 319   |   S   |      5      | 1.-4 |
    +---------+-----+--------+------+------+-------+-------+-------------+------+
    """
    type = 'FLUTTER'
    _field_map = {
        1: 'sid', 2:'method', 3:'density', 4:'mach', 5:'reduced_freq_velocity', 6:'imethod',
        8:'epsilon',
    }
    _properties = ['_field_map', 'headers', ]

    @classmethod
    def _init_from_empty(cls):
        sid = 1
        method = 'PKNL'
        density = 1
        mach = 1
        reduced_freq_velocity = 1
        return FLUTTER(sid, method, density, mach, reduced_freq_velocity,
                       imethod='L', nvalue=None, omax=None, epsilon=1.0e-3, comment='')

    def _get_field_helper(self, n):
        """
        Gets complicated parameters on the FLUTTER card

        Parameters
        ----------
        n : int
            the field number to update

        Returns
        -------
        value : int/float/str
            the value for the appropriate field

        """
        if n == 7:
            if self.method in ['K', 'KE']:
                value = self.nvalue
            elif self.method in ['PKS', 'PKNLS']:
                value = self.omax
            else:
                value = self.nvalue
            return value
        else:
            raise KeyError('Field %r is an invalid FLUTTER entry.' % (n))

    def _update_field_helper(self, n, value):
        """
        Updates complicated parameters on the FLUTTER card

        Parameters
        ----------
        n : int
            the field number to update
        value : int/float/str
            the value for the appropriate field

        """
        if n == 7:
            if self.method in ['K', 'KE']:
                self.nvalue = value
            elif self.method in ['PKS', 'PKNLS']:
                self.omax = value
            else:
                self.nvalue = value
        else:
            raise KeyError('Field %r=%r is an invalid FLUTTER entry.' % (n, value))

    def __init__(self, sid, method, density, mach, reduced_freq_velocity,
                 imethod='L', nvalue=None, omax=None, epsilon=1.0e-3, comment=''):
        """
        Creates a FLUTTER card, which is required for a flutter (SOL 145)
        analysis.

        Parameters
        ----------
        sid : int
            flutter id
        method : str
            valid methods = [K, KE,
                             PKS, PKNLS, PKNL, PKE]
        density : int
            defines a series of air densities in units of mass/volume
            PARAM,WTMASS does not affect this
            AERO affects this
            references an FLFACT id
        mach : int
            defines a series of the mach numbers
            references an FLFACT id
        reduced_freq_velocity : int
            Defines a series of either:
               1) reduced frequencies - K, KE
               2) velocities - PK, PKNL, PKS, PKNLS
            depending on the method chosen.
            references an FLFACT id
        imethod : str; default='L'
            Choice of interpolation method for aerodynamic matrix interpolation.
            imethods :
               1) L - linear
               2) S - surface
               3) TCUB - termwise cubic
        nvalue : int
            Number of eigenvalues beginning with the first eigenvalue for
            output and plots
        omax : float
            For the PKS and PKNLS methods, OMAX specifies the maximum frequency, in
            Hz., to be used in he flutter sweep.
            MSC only.
        epsilon : float; default=1.0e-3
            Convergence parameter for k. Used in the PK and PKNL methods only
        comment : str; default=''
            a comment for the card

        """
        BaseCard.__init__(self)
        if comment:
            self.comment = comment
        self.sid = sid
        if method in ['PK', 'PKNL', 'PKNLS']:
            imethod = 'L'
        #else:
            #assert imethod in ['S', 'L', None], imethod
        self.method = method
        self.density = density
        self.mach = mach

        # KFREQ - K, KE
        # VEL - PK, PKNL, PKS, PKNLS
        self.reduced_freq_velocity = reduced_freq_velocity

        #
        self.imethod = imethod
        self.nvalue = nvalue
        self.omax = omax
        self.epsilon = epsilon
        self.density_ref = None
        self.mach_ref = None
        self.reduced_freq_velocity_ref = None

    def validate(self):
        msg = ''
        if self.method not in ['K', 'KE', 'PK', 'PKNL', 'PKS', 'PKNLS']:
            msg += 'method = %r; allowed=[K, KE, PKS, PKNLS, PKNL, PK]\n' % self.method
        if self.imethod not in ['L', 'S', 'TCUB']:
            msg += 'imethod = %r; allowed=[L, S, TCUB]\n' % self.imethod
        if msg:
            raise ValueError(msg + str(self))

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a FLUTTER card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        sid = integer(card, 1, 'sid')
        method = string(card, 2, 'method (K, KE, PKS, PKNLS, PKNL, PK)')
        density_id = integer(card, 3, 'density')
        mach_id = integer(card, 4, 'mach')
        reduced_freq_velocity_id = integer(card, 5, 'reduced_freq_velocity')

        if method in ['K', 'KE']:
            imethod = string_or_blank(card, 6, 'imethod', 'L')
            nvalue = integer_or_blank(card, 7, 'nvalue')
            omax = None
            assert imethod in ['L', 'S', 'TCUB'], 'imethod = %s' % imethod  # linear-surface
        elif method in ['PKS', 'PKNLS']:
            imethod = None
            nvalue = None
            omax = double_or_blank(card, 7, 'omax')
        elif method == 'PKNL':
            nvalue = integer_or_blank(card, 7, 'nvalue')
            omax = None
            imethod = None
        elif method == 'PK':
            nvalue = integer_or_blank(card, 7, 'nvalue')
            omax = None
            imethod = None
        else:
            raise NotImplementedError('FLUTTER method=%r' % method)

        assert method in ['K', 'KE', 'PK', 'PKS', 'PKNL', 'PKNLS', None], method
        epsilon = double_or_blank(card, 8, 'epsilon', 1e-3)  # not defined in QRG
        assert len(card) <= 9, 'len(FLUTTER card) = %i\ncard=%s' % (len(card), card)
        return FLUTTER(sid, method, density_id, mach_id, reduced_freq_velocity_id,
                       imethod=imethod, nvalue=nvalue, omax=omax,
                       epsilon=epsilon, comment=comment)

    def make_flfacts_eas_sweep(self, model: BDF,
                               alt: float, eass: List[float],
                               alt_units: str='m',
                               velocity_units: str='m/s',
                               density_units: str='kg/m^3',
                               eas_units: str='m/s') -> Tuple[Any, Any, Any]:
        """
        Makes a sweep across equivalent airspeed for a constant altitude.

        Parameters
        ----------
        model : BDF
            the BDF model object
        alt : float
            Altitude in alt_units
        eass : List[float]
            Equivalent airspeed in eas_units
        alt_units : str; default='m'
            the altitude units; ft, kft, m
        velocity_units : str; default='m/s'
            the velocity units; ft/s, m/s, in/s, knots
        density_units : str; default='kg/m^3'
            the density units; slug/ft^3, slinch/in^3, kg/m^3
        eas_units : str; default='m/s'
            the equivalent airspeed units; ft/s, m/s, in/s, knots

        """
        eass.sort()
        rho, mach, velocity = make_flfacts_eas_sweep(
            alt, eass,
            alt_units=alt_units, velocity_units=velocity_units,
            density_units=density_units, eas_units=eas_units)
        flfact_rho = self.sid + 1
        flfact_mach = self.sid + 2
        flfact_velocity = self.sid + 3
        flfact_eas = self.sid + 4

        comment = ' density: min=%.3e max=%.3e %s' % (
            rho.min(), rho.max(), density_units,
        )
        model.add_flfact(flfact_rho, rho, comment=comment)
        model.add_flfact(flfact_mach, mach, comment=' Mach: %s' % mach.min())
        comment = ' velocity: min=%.3f max=%.3f %s' % (
            velocity.min(), velocity.max(), velocity_units)
        model.add_flfact(flfact_velocity, velocity, comment=comment)

        # eas in velocity units
        comment = ' EAS: min=%.3f max=%.3f %s' % (
            eass.min(), eass.max(), eas_units)
        model.add_flfact(flfact_eas, eass, comment=comment)

    def make_flfacts_alt_sweep(self, model: BDF, mach, alts, eas_limit: float=1000.,
                               alt_units: str='m',
                               velocity_units: str='m/s',
                               density_units: str='kg/m^3',
                               eas_units: str='m/s') -> Tuple[Any, Any, Any]:
        """makes an altitude sweep"""
        alts.sort()
        alts = alts[::-1]
        rho, mach, velocity = make_flfacts_alt_sweep(
            mach, alts, eas_limit=eas_limit,
            alt_units=alt_units,
            velocity_units=velocity_units,
            density_units=density_units,
            eas_units=eas_units)
        flfact_rho = self.sid + 1
        flfact_mach = self.sid + 2
        flfact_velocity = self.sid + 3
        flfact_eas = self.sid + 4
        flfact_alt = self.sid + 5

        alts2 = alts[:len(rho)]
        assert len(rho) == len(alts2)
        comment = ' density: min=%.3e max=%.3e %s; alt min=%.0f max=%.0f %s' % (
            rho.min(), rho.max(), density_units,
            alts2.min(), alts2.max(), alt_units,
        )
        model.add_flfact(flfact_rho, rho, comment=comment)
        model.add_flfact(flfact_mach, mach, comment=' Mach: %s' % mach.min())
        comment = ' velocity: min=%.3f max=%.3f %s' % (
            velocity.min(), velocity.max(), velocity_units)
        model.add_flfact(flfact_velocity, velocity, comment=comment)

        # eas in velocity units
        rho0 = atm_density(0., alt_units=alt_units, density_units=density_units)
        eas = velocity * np.sqrt(rho / rho0)
        kvel = _velocity_factor(velocity_units, eas_units)

        eas_in_eas_units = eas * kvel
        comment = ' EAS: min=%.3f max=%.3f %s' % (
            eas_in_eas_units.min(), eas_in_eas_units.max(), eas_units)
        model.add_flfact(flfact_eas, eas_in_eas_units, comment=comment)

        comment = ' Alt: min=%.3f max=%.3f %s' % (alts2.min(), alts2.max(), alt_units)
        model.add_flfact(flfact_alt, alts2, comment=comment)

    def make_flfacts_mach_sweep(self, model, alt, machs, eas_limit=1000., alt_units='m',
                                velocity_units='m/s',
                                density_units='kg/m^3',
                                eas_units='m/s'):
        """makes a mach sweep"""
        machs.sort()
        machs = machs[::-1]
        rho, mach, velocity = make_flfacts_mach_sweep(
            alt, machs, eas_limit=eas_limit,
            alt_units=alt_units,
            velocity_units=velocity_units,
            density_units=density_units,
            eas_units=eas_units)

        machs2 = machs[:len(rho)]
        assert len(rho) == len(machs2)

        flfact_rho = self.sid + 1
        flfact_mach = self.sid + 2
        flfact_velocity = self.sid + 3
        flfact_eas = self.sid + 4

        comment = ' density: min=%.3e max=%.3e %s; alt %.0f %s' % (
            rho.min(), rho.max(), density_units,
            alt, alt_units,
        )
        model.add_flfact(flfact_rho, rho, comment=comment)
        comment = ' Mach: min=%s max=%s' % (mach.min(), mach.max())
        model.add_flfact(flfact_mach, mach, comment=comment)
        comment = ' velocity: min=%.3f max=%.3f %s' % (
            velocity.min(), velocity.max(), velocity_units)
        model.add_flfact(flfact_velocity, velocity, comment=comment)

        # eas in velocity units
        rho0 = atm_density(0., alt_units=alt_units, density_units=density_units)
        eas = velocity * np.sqrt(rho / rho0)
        kvel = _velocity_factor(velocity_units, eas_units)

        eas_in_eas_units = eas * kvel
        comment = ' EAS: min=%.3f max=%.3f %s' % (
            eas_in_eas_units.min(), eas_in_eas_units.max(), eas_units)
        model.add_flfact(flfact_eas, eas_in_eas_units, comment=comment)

    @property
    def headers(self):
        headers = ['density', 'mach']
        if self.method in ['PK', 'PKS', 'PKNL', 'PKNLS']:
            headers.append('velocity')
        elif self.method in ['K', 'KE']:
            headers.append('reduced_frequency')
        else:
            raise NotImplementedError('FLUTTER method=%r' % self.method)
        return headers

    @classmethod
    def add_op2_data(cls, data, comment=''):
        assert len(data) == 8, 'FLUTTER = %s' % data
        sid = data[0]
        method = data[1]
        density = data[2]
        mach = data[3]
        reduced_freq_velocity = data[4]
        method = data[5]
        imethod = data[6]
        nvalue = data[7]
        omax = data[8]
        epsilon = None
        return FLUTTER(sid, method, density, mach, reduced_freq_velocity,
                       imethod, nvalue, omax,
                       epsilon, comment=comment)

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        msg = ', which is required by FLUTTER sid=%s' % self.sid
        self.density_ref = model.FLFACT(self.density, msg=msg)
        self.mach_ref = model.FLFACT(self.mach, msg=msg)
        self.reduced_freq_velocity_ref = model.FLFACT(self.reduced_freq_velocity, msg=msg)

    def safe_cross_reference(self, model):
        msg = ', which is required by FLUTTER sid=%s' % self.sid
        try:
            self.density_ref = model.FLFACT(self.density, msg=msg)
        except KeyError:
            pass
        try:
            self.mach_ref = model.FLFACT(self.mach, msg=msg)
        except KeyError:
            pass
        try:
            self.reduced_freq_velocity_ref = model.FLFACT(self.reduced_freq_velocity, msg=msg)
        except KeyError:
            pass

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.density = self.get_density()
        self.mach = self.get_mach()
        self.reduced_freq_velocity = self.get_rfreq_vel()
        self.density_ref = None
        self.mach_ref = None
        self.reduced_freq_velocity_ref = None

    def get_density(self):
        if self.density_ref is not None:
            return self.density_ref.sid
        return self.density

    def get_mach(self):
        if self.mach_ref is not None:
            return self.mach_ref.sid
        return self.mach

    def get_rfreq_vel(self):
        if self.reduced_freq_velocity_ref is not None:
            return self.reduced_freq_velocity_ref.sid
        return self.reduced_freq_velocity

    def _get_raw_nvalue_omax(self):
        if self.method in ['K', 'KE']:
            #assert self.imethod in ['L', 'S'], 'imethod = %s' % self.imethod
            return self.imethod, self.nvalue
        elif self.method in ['PKS', 'PKNLS']:
            return self.imethod, self.omax
        # PK, PKNL
        return self.imethod, self.nvalue

    def _get_repr_nvalue_omax(self):
        if self.method in ['K', 'KE']:
            imethod = set_blank_if_default(self.imethod, 'L')
            #assert self.imethod in ['L', 'S'], 'imethod = %s' % self.imethods
            return imethod, self.nvalue
        elif self.method in ['PKS', 'PKNLS']:
            return self.imethod, self.omax
        # PK, PKNL
        return self.imethod, self.nvalue

    def raw_fields(self):
        """
        Gets the fields in their unmodified form

        Returns
        -------
        fields : list[varies]
            the fields that define the card

        """
        (imethod, nvalue) = self._get_raw_nvalue_omax()
        list_fields = ['FLUTTER', self.sid, self.method, self.get_density(),
                       self.get_mach(), self.get_rfreq_vel(), imethod, nvalue, self.epsilon]
        return list_fields

    def repr_fields(self):
        (imethod, nvalue) = self._get_repr_nvalue_omax()
        epsilon = set_blank_if_default(self.epsilon, 0.001)
        list_fields = ['FLUTTER', self.sid, self.method, self.get_density(), self.get_mach(),
                       self.get_rfreq_vel(), imethod, nvalue, epsilon]
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class GUST(BaseCard):
    """
    Defines a stationary vertical gust for use in aeroelastic response
    analysis.

    +------+-----+-------+-----+-----+------+
    |   1  |  2  |   3   |  4  |  5  |  6   |
    +======+=====+=======+=====+=====+======+
    | GUST | SID | DLOAD | WG  | X0  |  V   |
    +------+-----+-------+-----+-----+------+
    | GUST | 133 |   61  | 1.0 | 0.  | 1.+4 |
    +------+-----+-------+-----+-----+------+
    """
    type = 'GUST'
    _field_map = {
        1: 'sid', 2:'dload', 3:'wg', 4:'x0', 5:'V',
    }

    @classmethod
    def _init_from_empty(cls):
        sid = 1
        dload = 1
        wg = 1.
        x0 = 0.
        return GUST(sid, dload, wg, x0, V=None, comment='')

    def __init__(self, sid, dload, wg, x0, V=None, comment=''):
        """
        Creates a GUST card, which defines a stationary vertical gust
        for use in aeroelastic response analysis.

        Parameters
        ----------
        sid : int
            gust load id
        dload : int
            TLOADx or RLOADx entry that defines the time/frequency
            dependence
        wg : float
            Scale factor (gust velocity/forward velocity) for gust
            velocity
        x0 : float
            Streamwise location in the aerodynamic coordinate system of
            the gust reference point.
        V : float; default=None
            float : velocity of the vehicle (must be the same as the
                    velocity on the AERO card)
            None : ???
        comment : str; default=''
            a comment for the card

        """
        BaseCard.__init__(self)
        if comment:
            self.comment = comment
        self.sid = sid
        self.dload = dload
        self.wg = wg
        self.x0 = x0
        self.V = V

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a GUST card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        sid = integer(card, 1, 'sid')
        dload = integer(card, 2, 'dload')
        wg = double(card, 3, 'wg')
        x0 = double(card, 4, 'x0')
        V = double_or_blank(card, 5, 'V')
        assert len(card) <= 6, 'len(GUST card) = %i\ncard=%s' % (len(card), card)
        return GUST(sid, dload, wg, x0, V=V, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        sid = data[0]
        dload = data[1]
        wg = data[2]
        x0 = data[3]
        V = data[4]
        assert len(data) == 5, 'data = %s' % data
        return GUST(sid, dload, wg, x0, V, comment=comment)

    #def Angle(self):
        #angle = self.wg * self.t * (t-(x-self.x0) / self.V) # T is the tabular
        #return angle

    #def uncross_reference(self) -> None:
        #pass

    def _verify(self, model, xref):
        if model.aero:
            pass
            #assert model.aero.V == self.V

    def raw_fields(self):
        """
        Gets the fields in their unmodified form

        Returns
        -------
        fields : list[varies]
            the fields that define the card

        """
        list_fields = ['GUST', self.sid, self.dload, self.wg, self.x0, self.V]
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class MKAERO1(BaseCard):
    """
    Provides a table of Mach numbers (m) and reduced frequencies (k) for
    aerodynamic matrix calculation.

    +---------+----+----+----+----+----+----+----+----+
    |    1    |  2 | 3  |  4 | 5  | 6  | 7  | 8  | 9  |
    +=========+====+====+====+====+====+====+====+====+
    | MKAERO1 | m1 | m2 | m3 | m4 | m5 | m6 | m7 | m8 |
    +---------+----+----+----+----+----+----+----+----+
    |         | k1 | k2 | k3 | k4 | k5 | k6 | k7 | k8 |
    +---------+----+----+----+----+----+----+----+----+
    """
    type = 'MKAERO1'

    @classmethod
    def _init_from_empty(cls):
        machs = [1.]
        reduced_freqs = [1.]
        return MKAERO1(machs, reduced_freqs, comment='')

    def __init__(self, machs, reduced_freqs, comment=''):
        """
        Creates an MKAERO1 card, which defines a set of mach and
        reduced frequencies.

        Parameters
        ----------
        machs : List[float]
            series of Mach numbers
        reduced_freqs : List[float]
            series of reduced frequencies
        comment : str; default=''
            a comment for the card

        """
        BaseCard.__init__(self)
        if comment:
            self.comment = comment
        self.machs = np.unique(machs)
        self.reduced_freqs = np.unique(reduced_freqs)

    def validate(self):
        msg = ''
        if None in self.machs:
            msg += 'MKAERO1; None in machs=%s\n' % (self.machs)
        if None in self.reduced_freqs:
            msg += 'MKAERO1; None in rfreqs=%s\n' % (self.reduced_freqs)
        if len(self.machs) == 0:
            msg += 'MKAERO1; nmachs=%s machs=%s\n' % (len(self.machs), self.machs)
        if len(self.reduced_freqs) == 0:
            msg += 'MKAERO1; nrfreqs=%s rfreqs=%s' % (len(self.reduced_freqs), self.reduced_freqs)
        if msg:
            raise ValueError(msg.rstrip())

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds an MKAERO1 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        list_fields = [interpret_value(field, card) for field in card[1:]]
        nfields = len(list_fields) - 8
        machs = []
        reduced_freqs = []
        for i in range(1, 1 + nfields):
            machs.append(double_or_blank(card, i, 'mach'))
            reduced_freqs.append(double_or_blank(card, i + 8, 'rFreq'))
        machs = wipe_empty_fields(machs)
        reduced_freqs = wipe_empty_fields(reduced_freqs)
        return MKAERO1(machs, reduced_freqs, comment=comment)

    def mklist(self):
        mklist = []
        for mach in self.machs:
            for kfreq in self.reduced_freqs:
                mklist.append([mach, kfreq])
        return mklist

    def raw_fields(self):
        """
        Gets the fields in their unmodified form

        Returns
        -------
        fields : list[varies]
            the fields that define the card

        """
        #list_fields = ['MKAERO1']
        #for (i, mach, rfreq) in zip(count(), self.machs, self.reduced_freqs):
        #    list_fields += [mach, rfreq]

        # kind of a hack because there isn't a good way to do this for
        # duplicately-defined MKAERO1s
        machs = [None] * max(8, len(self.machs))
        freqs = [None] * max(8, len(self.reduced_freqs))
        for i, mach in enumerate(self.machs):
            machs[i] = mach
        for i, freq in enumerate(self.reduced_freqs):
            freqs[i] = freq
        list_fields = ['MKAERO1'] + machs + freqs
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        nmachs = len(self.machs)
        nreduced_freqs = len(self.reduced_freqs)
        if nmachs > 8 or nreduced_freqs > 8:
            mach_sets = []
            rfreq_sets = []
            imach = 0
            ifreq = 0
            while imach < nmachs:
                mach_sets.append(self.machs[imach:imach+8])
                imach += 8
            while ifreq < nreduced_freqs:
                rfreq_sets.append(self.reduced_freqs[ifreq:ifreq+8])
                ifreq += 8
            msg = self.comment

            #print('mach_sets = %s' % mach_sets)
            #print('rfreq_sets = %s' % rfreq_sets)
            for mach_set in mach_sets:
                for rfreq_set in rfreq_sets:
                    msg += MKAERO1(mach_set, rfreq_set).write_card(
                        size=size, is_double=is_double)
            return msg

        machs = [None] * 8
        reduced_freqs = [None] * 8
        if not 0 < len(self.machs) <= 8:
            msg = 'MKAERO1; nmachs=%s machs=%s' % (len(self.machs), self.machs)
            raise ValueError(msg)
        if not 0 < len(self.reduced_freqs) <= 8:
            msg = 'MKAERO1; nrfreqs=%s rfreqs=%s' % (len(self.reduced_freqs), self.reduced_freqs)
            raise ValueError(msg)

        for i, mach in zip(count(), self.machs):
            machs[i] = mach
        for i, rfreq in zip(count(), self.reduced_freqs):
            reduced_freqs[i] = rfreq
        return self.comment + print_card_8(['MKAERO1'] + machs + reduced_freqs)

    def __repr__(self):
        return self.write_card()


class MKAERO2(BaseCard):
    """
    Provides a table of Mach numbers (m) and reduced frequencies (k) for
    aerodynamic matrix calculation.

    +---------+----+----+----+----+----+----+----+----+
    |    1    | 2  | 3  | 4  | 5  | 6  | 7  | 8  | 9  |
    +=========+====+====+====+====+====+====+====+====+
    | MKAERO2 | m1 | k1 | m2 | k2 | m3 | k3 | m4 | k4 |
    +---------+----+----+----+----+----+----+----+----+
    """
    type = 'MKAERO2'

    @classmethod
    def _init_from_empty(cls):
        machs = [1.]
        reduced_freqs = [1.]
        return MKAERO2(machs, reduced_freqs, comment='')

    def __init__(self, machs, reduced_freqs, comment=''):
        """
        Creates an MKAERO2 card, which defines a set of mach and
        reduced frequency pairs.

        Parameters
        ----------
        machs : List[float]
            series of Mach numbers
        reduced_freqs : List[float]
            series of reduced frequencies
        comment : str; default=''
            a comment for the card

        """
        BaseCard.__init__(self)
        if comment:
            self.comment = comment
        self.machs = machs
        self.reduced_freqs = reduced_freqs

    def validate(self):
        if len(self.machs) == 0:
            msg = 'MKAERO2; nmachs=%s machs=%s' % (len(self.machs), self.machs)
            raise ValueError(msg)
        if len(self.reduced_freqs) == 0:
            msg = 'MKAERO2; nrfreqs=%s rfreqs=%s' % (len(self.reduced_freqs), self.reduced_freqs)
            raise ValueError(msg)
        if len(self.machs) != len(self.reduced_freqs):
            msg = 'MKAERO2; len(machs)=%s len(rfreqs)=%s; should be the same' % (
                len(self.machs), len(self.reduced_freqs))
            raise ValueError(msg)

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds an MKAERO2 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        list_fields = card.fields(1)
        nfields = len(list_fields)
        machs = []
        reduced_freqs = []
        for i in range(1, 1 + nfields, 2):
            machs.append(double(card, i, 'mach'))
            reduced_freqs.append(double(card, i + 1, 'rFreq'))
        return MKAERO2(machs, reduced_freqs, comment=comment)


    def mklist(self):
        mklist = []
        for mach, kfreq in zip(self.machs, self.reduced_freqs):
            mklist.append([mach, kfreq])
        return mklist

    def raw_fields(self):
        """
        Gets the fields in their unmodified form

        Returns
        -------
        fields : list[varies]
            the fields that define the card

        """
        list_fields = ['MKAERO2']
        for (mach, rfreq) in zip(self.machs, self.reduced_freqs):
            list_fields += [mach, rfreq]
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        cards = []
        list_fields = ['MKAERO2']
        nvalues = 0
        for mach, rfreq in zip(self.machs, self.reduced_freqs):
            list_fields += [mach, rfreq]
            nvalues += 1
            if nvalues == 4:
                cards.append(print_card_8(list_fields))
                list_fields = ['MKAERO2']
                nvalues = 0
        if nvalues:
            cards.append(print_card_8(list_fields))
        else:
            if len(self.machs) != len(self.reduced_freqs) or len(self.machs) == 0:
                msg = 'MKAERO2: len(machs)=%s len(reduced_freqs)=%s' % (
                    len(self.machs), len(self.reduced_freqs))
                raise ValueError(msg)

        return self.comment + ''.join(cards)

    def __repr__(self):
        return self.write_card()
