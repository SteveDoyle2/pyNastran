from __future__ import annotations
from typing import Optional, TYPE_CHECKING
import numpy as np

from pyNastran.bdf.field_writer_8 import set_blank_if_default
from pyNastran.bdf.bdf_interface.assign_type import (
    integer, integer_or_blank,
    double, double_or_blank,
    double_string_or_blank, string_or_blank, fields,)

from pyNastran.bdf.cards.base_card import expand_thru

from pyNastran.dev.bdf_vectorized3.cards.base_card import (
    VectorizedBaseCard, make_idim, hslice_by_idim,
    parse_check, save_ifile_comment)
from pyNastran.dev.bdf_vectorized3.cards.write_utils import (
    array_str,
    array_default_int, array_default_float, array_default_str,
    array_float_nan, get_print_card_size)
from pyNastran.bdf.cards.aero.dynamic_loads import von_karman_psd, dryden_psd

if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf_interface.bdf_card import BDFCard
    #from pyNastran.dev.bdf_vectorized3.bdf import BDF
    from pyNastran.dev.bdf_vectorized3.types import TextIOLike


class FLFACT(VectorizedBaseCard):
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
    _id_name = 'flfact_id'
    _skip_equality_check = True  # assume unequal
    @VectorizedBaseCard.clear_check
    def clear(self) -> None:
        self.flfact_id = np.array([], dtype='int32')
        self.nfactors = np.array([], dtype='int32')
        self.factors = np.array([], dtype='float64')

    def slice_card_by_flfact_id(self, flfact_id: np.ndarray) -> FLFACT:
        assert self.n > 0, self.n
        assert len(self.flfact_id) > 0, self.flfact_id
        i = self.index(flfact_id)
        cls_obj = self.slice_card_by_index(i)
        assert cls_obj.n > 0, cls_obj
        return cls_obj

    def geom_check(self, missing: dict[str, np.ndarray]):
        pass

    def add(self, sid: int, factors: list[float],
            ifile: int=0, comment: str='') -> int:
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
            values : list[float, ..., float]
                list of factors
            list[f1, THRU, fnf, nf, fmid]
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
        self.cards.append((sid, factors, ifile, comment))
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, ifile: int, comment: str='') -> int:
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
        #return FLFACT(sid, factors, comment=comment)
        self.cards.append((sid, factors, ifile, comment))
        self.n += 1
        return self.n - 1

    @VectorizedBaseCard.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        ifile = np.zeros(ncards, dtype='int32')
        flfact_id = np.zeros(ncards, dtype='int32')
        nfactors = np.zeros(ncards, dtype='int32')
        #self.elements = np.array([], dtype='int32')
        comment = {}
        all_factors = []
        for icard, card in enumerate(self.cards):
            sid, factors, ifilei, commenti = card
            factors = expand_thru(factors, set_fields=False, sort_fields=False)
            ifile[icard] = ifilei
            if commenti:
                comment[sid] = commenti
            flfact_id[icard] = sid
            nfactors[icard] = len(factors)
            all_factors.extend(factors)
        factors = np.array(all_factors, dtype='float64')
        self._save(flfact_id, nfactors, factors,
                   ifile=ifile, comment=comment)
        self.cards = []

    def _save(self, flfact_id, nfactors, factors,
              ifile=None, comment=None):
        ncards = len(flfact_id)
        if ifile is None:
            ifile = np.zeros(ncards, dtype='int32')
        if len(self.flfact_id):
            ifile = np.hstack([self.ifile, ifile])
            flfact_id = np.hstack([self.flfact_id, flfact_id])
            factors = np.hstack([self.factors, factors])
            nfactors = np.hstack([self.nfactors, nfactors])
        save_ifile_comment(self, ifile, comment)
        self.flfact_id = flfact_id
        self.nfactors = nfactors
        self.factors = factors
        self.n = len(flfact_id)

    def sort(self) -> None:
        uflfact = np.unique(self.flfact_id)
        if np.array_equal(uflfact, self.flfact_id):
            return
        i = np.argsort(self.flfact_id)
        self.__apply_slice__(self, i)

    def __apply_slice__(self, flfact: FLFACT, i: np.ndarray) -> None:
        flfact.flfact_id = self.flfact_id[i]
        ifactor = self.ifactor
        flfact.factors = hslice_by_idim(i, ifactor, self.factors)
        flfact.nfactors = self.nfactors[i]
        flfact.n = len(i)

    @property
    def ifactor(self) -> np.ndarray:
        return make_idim(self.n, self.nfactors)

    @property
    def max_id(self) -> int:
        return self.flfact_id.max()

    @parse_check
    def write_file(self, bdf_file: TextIOLike, size: int=8,
                   is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card, size = get_print_card_size(size, self.max_id)

        for sid, (ifactor0, ifactor1) in zip(self.flfact_id, self.ifactor):
            factors = self.factors[ifactor0:ifactor1].tolist()
            list_fields = ['FLFACT', sid] + factors
            bdf_file.write(print_card(list_fields))
        return


class GUST(VectorizedBaseCard):
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
    _id_name = 'gust_id'
    @VectorizedBaseCard.clear_check
    def clear(self) -> None:
        self.gust_id = np.array([], dtype='int32')

    #def __len__(self) -> int:
        #return len(self.name)

    def add(self, sid: int, dload: int, wg: float, x0: float,
            V: Optional[float]=None,
            ifile: int=0, comment: str='') -> int:
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
        self.cards.append((sid, dload, wg, x0, V, ifile, comment))
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, ifile: int, comment: str='') -> int:
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
        assert len(card) <= 6, f'len(GUST card) = {len(card):d}\ncard={card}'
        #return GUST(sid, dload, wg, x0, V=V, comment=comment)
        self.cards.append((sid, dload, wg, x0, V, ifile, comment))
        self.n += 1
        return self.n - 1

    @VectorizedBaseCard.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        ifile = np.zeros(ncards, dtype='int32')
        gust_id = np.zeros(ncards, dtype='int32')
        dload_id = np.zeros(ncards, dtype='int32')

        wg = np.zeros(ncards, dtype='float64')
        x0 = np.zeros(ncards, dtype='float64')
        V = np.zeros(ncards, dtype='float64')
        comment = {}
        for icard, card in enumerate(self.cards):
            (gust_idi, dloadi, wgi, x0i, Vi, ifilei, commenti) = card
            ifile[icard] = ifilei
            if commenti:
                comment[gust_idi] = commenti
            gust_id[icard] = gust_idi
            dload_id[icard] = dloadi
            wg[icard] = wgi
            x0[icard] = x0i
            V[icard] = Vi
        self._save(gust_id, dload_id, wg, x0, V,
                   ifile=ifile, comment=comment)
        #self.sort()
        self.cards = []

    def _save(self, gust_id, dload_id, wg, x0, V,
              ifile=None, comment=None):
        ncards = len(gust_id)
        if ifile is None:
            ifile = np.zeros(ncards, dtype='int32')
        if len(self.gust_id):
            ifile = np.stack([self.ifile, ifile])
            gust_id = np.hstack([self.gust_id, gust_id])
            dload_id = np.hstack([self.dload_id, dload_id])
            wg = np.hstack([self.wg, wg])
            x0 = np.hstack([self.x0, x0])
            V = np.hstack([self.V, V])
        save_ifile_comment(self, ifile, comment)
        self.gust_id = gust_id
        self.dload_id = dload_id
        self.wg = wg
        self.x0 = x0
        self.V = V

    def __apply_slice__(self, gust: GUST, i: np.ndarray) -> None:
        gust.gust_id = self.gust_id[i]
        gust.dload_id = self.dload_id[i]
        gust.wg = self.wg[i]
        gust.x0 = self.x0[i]
        gust.V = self.V[i]
        gust.n = len(i)

    def convert(self, xyz_scale: float=1.0,
                velocity_scale: float=1.0, **kwargs) -> None:
        self.x0 *= xyz_scale
        self.V *= velocity_scale

    def geom_check(self, missing: dict[str, np.ndarray]):
        #mids = hstack_msg([prop.material_id for prop in self.allowed_materials],
                          #msg=f'no materials for {self.type}')
        #mids.sort()
        coords = self.model.coord.coord_id
        #all_aecomp_names = self.model.aecomp.name
        #aecomp_names = np.unique(self.comp)
        #ucoords = np.unique(np.hstack([self.cp, self.cd]))
        #geom_check(self,
                   #missing,
                   #coord=(coords, ucoords),
                   #aecomp=(all_aecomp_names, aecomp_names))

    @property
    def max_id(self) -> int:
        return max(self.gust_id.max(), self.dload_id.max())

    @parse_check
    def write_file(self, bdf_file: TextIOLike, size: int=8,
                   is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card, size = get_print_card_size(size, self.max_id)

        gust_ids = array_str(self.gust_id, size=size)
        dload_ids = array_str(self.dload_id, size=size)
        for gust_id, dload, wg, x0, V in zip(gust_ids, dload_ids, self.wg, self.x0, self.V):
            list_fields = ['GUST', gust_id, dload, wg, x0, V]
            bdf_file.write(print_card(list_fields))
        return

    @staticmethod
    def von_karman_psd(omega: np.ndarray, sigma: float, L: float,
                       V: float) -> np.ndarray:
        """Von Karman vertical turbulence PSD.

        Parameters
        ----------
        omega : (N,) ndarray
            Circular frequency [rad/s].
        sigma : float
            RMS gust velocity.
        L : float
            Turbulence scale length.
        V : float
            True airspeed.

        Returns
        -------
        (N,) ndarray
            PSD [velocity^2 / (rad/s)].
        """
        return von_karman_psd(omega, sigma, L, V)

    @staticmethod
    def dryden_psd(omega: np.ndarray, sigma: float, L: float,
                   V: float) -> np.ndarray:
        """Dryden vertical turbulence PSD.

        Parameters
        ----------
        omega : (N,) ndarray
            Circular frequency [rad/s].
        sigma : float
            RMS gust velocity.
        L : float
            Turbulence scale length.
        V : float
            True airspeed.

        Returns
        -------
        (N,) ndarray
            PSD [velocity^2 / (rad/s)].
        """
        return dryden_psd(omega, sigma, L, V)

    @staticmethod
    def gust_psd(omega: np.ndarray, sigma: float, L: float, V: float,
                 gust_model: str = 'von_karman') -> np.ndarray:
        """Compute gust PSD using the specified turbulence model.

        Parameters
        ----------
        omega : (N,) ndarray
            Circular frequency [rad/s].
        sigma : float
            RMS gust velocity.
        L : float
            Turbulence scale length.
        V : float
            True airspeed.
        gust_model : str
            'von_karman' or 'dryden'.

        Returns
        -------
        (N,) ndarray
            PSD [velocity^2 / (rad/s)].
        """
        if gust_model == 'von_karman':
            return von_karman_psd(omega, sigma, L, V)
        elif gust_model == 'dryden':
            return dryden_psd(omega, sigma, L, V)
        raise ValueError(f"Unknown gust model: '{gust_model}'. "
                         f"Use 'von_karman' or 'dryden'.")


FLUTTER_METHODS = {'K', 'KE', 'PK', 'PKS', 'PKNL', 'PKNLS'}
class FLUTTER(VectorizedBaseCard):
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
    _id_name = 'flutter_id'
    @VectorizedBaseCard.clear_check
    def clear(self) -> None:
        self.flutter_id = np.array([], dtype='int32')
        self.method = np.array([], dtype='|U8')
        self.density_flfact_id = np.array([], dtype='int32')
        self.mach_flfact_id = np.array([], dtype='int32')
        self.rfreq_velocity_flfact_id = np.array([], dtype='int32')
        self.imethod = np.array([], dtype='|U8')
        self.nvalue = np.array([], dtype='int32')
        self.omax = np.array([], dtype='float64')
        self.eps = np.array([], dtype='float64')

    #def __len__(self) -> int:
        #return len(self.name)

    def add(self, flutter_id: int, method: str,
            density: int, mach: int, reduced_freq_velocity: int,
            imethod: str='L',
            nvalue=None, omax=None,
            epsilon: float=1.0e-3,
            ifile: int=0, comment: str='', validate: bool=False) -> int:
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
        #(flutter_idi, methodi,
         #density_flfact_idi, mach_flfact_idi, rfreq_flfact_idi,
         #nvaluei, omaxi, epsiloni, comment) = card
        assert method in FLUTTER_METHODS, f'method={method} allowed={FLUTTER_METHODS}'
        self.cards.append((flutter_id, method, imethod,
                           density, mach, reduced_freq_velocity,
                           nvalue, omax, epsilon, ifile, comment))
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, ifile: int, comment: str='') -> int:
        """
        Adds a GUST card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        flutter_id = integer(card, 1, 'sid')
        method = string_or_blank(card, 2, 'method (K, KE, PKS, PKNLS, PKNL, PK)', default='L')
        density_id = integer(card, 3, 'density')
        mach_id = integer(card, 4, 'mach')
        reduced_freq_velocity_id = integer(card, 5, 'reduced_freq_velocity')

        omax = None
        imethod = string_or_blank(card, 6, 'imethod', default='L')
        if method in ['K', 'KE']:
            nvalue = integer_or_blank(card, 7, 'nvalue')
            assert imethod in ['L', 'S', 'TCUB'], 'imethod = %s' % imethod  # linear-surface
        elif method in ['PKS', 'PKNLS']:
            nvalue = None
            omax = double_or_blank(card, 7, 'omax')
        elif method == 'PKNL':
            nvalue = integer_or_blank(card, 7, 'nvalue')
        elif method == 'PK':
            nvalue = integer_or_blank(card, 7, 'nvalue')
        else:  # pragma: no cover
            raise NotImplementedError('FLUTTER method=%r' % method)

        assert method in FLUTTER_METHODS, f'method={method} allowed={FLUTTER_METHODS}'
        epsilon = double_or_blank(card, 8, 'epsilon', default=1e-3)  # not defined in QRG
        assert len(card) <= 9, f'len(FLUTTER card) = {len(card):d}\ncard={card}'
        #return FLUTTER(sid, method, density_id, mach_id, reduced_freq_velocity_id,
                       #imethod=imethod, nvalue=nvalue, omax=omax,
                       #epsilon=epsilon, comment=comment)
        self.cards.append((flutter_id, method, imethod,
                           density_id, mach_id, reduced_freq_velocity_id,
                           nvalue, omax, epsilon, ifile, comment))
        self.n += 1
        return self.n - 1

    @VectorizedBaseCard.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        ifile = np.zeros(ncards, dtype='int32')
        flutter_id = np.zeros(ncards, dtype='int32')
        method = np.zeros(ncards, dtype='|U8')
        density_flfact_id = np.zeros(ncards, dtype='int32')
        mach_flfact_id = np.zeros(ncards, dtype='int32')
        rfreq_velocity_flfact_id = np.zeros(ncards, dtype='int32')
        imethod = np.zeros(ncards, dtype='|U8')
        nvalue = np.zeros(ncards, dtype='int32')
        omax = np.zeros(ncards, dtype='float64')
        epsilon = np.zeros(ncards, dtype='float64')
        comment = {}
        for icard, card in enumerate(self.cards):
            (flutter_idi, methodi, imethodi,
             density_flfact_idi, mach_flfact_idi, rfreq_velocity_flfact_idi,
             nvaluei, omaxi, epsiloni, ifilei, commenti) = card
            if nvaluei is None:
                nvaluei = 0
            if omaxi is None:
                omaxi = 0
            ifile[icard] = ifilei
            if commenti:
                comment[flutter_idi] = commenti
            flutter_id[icard] = flutter_idi
            method[icard] = methodi
            imethod[icard] = imethodi
            density_flfact_id[icard] = density_flfact_idi
            mach_flfact_id[icard] = mach_flfact_idi
            rfreq_velocity_flfact_id[icard] = rfreq_velocity_flfact_idi
            nvalue[icard] = nvaluei
            omax[icard] = omaxi
            epsilon[icard] = epsiloni
        self._save(flutter_id, method, imethod,
                   density_flfact_id, mach_flfact_id, rfreq_velocity_flfact_id,
                   nvalue, omax, epsilon,
                   ifile=ifile, comment=comment)
        self.sort()
        self.cards = []

    def _save(self, flutter_id, method, imethod,
              density_flfact_id, mach_flfact_id, rfreq_velocity_flfact_id,
              nvalue, omax, epsilon,
              ifile=None, comment=None):
        ncards = len(flutter_id)
        if ifile is None:
            ifile = np.zeros(ncards, dtype='int32')
        if len(self.flutter_id):
            ifile = np.stack([self.ifile, ifile])
            flutter_id = np.hstack([self.flutter_id, flutter_id])
            method = np.hstack([self.method, method])
            imethod = np.hstack([self.imethod, imethod])
            density_flfact_id = np.hstack([self.density_flfact_id, density_flfact_id])
            mach_flfact_id = np.hstack([self.mach_flfact_id, mach_flfact_id])
            rfreq_velocity_flfact_id = np.hstack([
                self.rfreq_velocity_flfact_id,
                rfreq_velocity_flfact_id])
            nvalue = np.hstack([self.nvalue, nvalue])
            omax = np.hstack([self.omax, omax])
            epsilon = np.hstack([self.epsilon, epsilon])
        save_ifile_comment(self, ifile, comment)
        self.flutter_id = flutter_id
        self.method = method
        self.imethod = imethod
        self.density_flfact_id = density_flfact_id
        self.mach_flfact_id = mach_flfact_id
        self.rfreq_velocity_flfact_id = rfreq_velocity_flfact_id
        self.nvalue = nvalue
        self.omax = omax
        self.epsilon = epsilon

    def __apply_slice__(self, flutter: FLUTTER, i: np.ndarray) -> None:
        flutter.flutter_id = self.flutter_id[i]
        flutter.method = self.method[i]
        flutter.density_flfact_id = self.density_flfact_id[i]
        flutter.mach_flfact_id = self.mach_flfact_id[i]
        flutter.rfreq_velocity_flfact_id = self.rfreq_velocity_flfact_id[i]
        flutter.imethod = self.imethod[i]
        flutter.nvalue = self.nvalue[i]
        flutter.omax = self.omax[i]
        flutter.epsilon = self.epsilon[i]
        flutter.n = len(i)

    def geom_check(self, missing: dict[str, np.ndarray]):
        #mids = hstack_msg([prop.material_id for prop in self.allowed_materials],
                          #msg=f'no materials for {self.type}')
        #mids.sort()
        coords = self.model.coord.coord_id
        #all_aecomp_names = self.model.aecomp.name
        #aecomp_names = np.unique(self.comp)
        #ucoords = np.unique(np.hstack([self.cp, self.cd]))
        #geom_check(self,
                   #missing,
                   #coord=(coords, ucoords),
                   #aecomp=(all_aecomp_names, aecomp_names))

    def _get_repr_nvalue_omax(self):
        if self.method in ['K', 'KE']:
            imethod = set_blank_if_default(self.imethod, 'L')
            #assert self.imethod in ['L', 'S'], 'imethod = %s' % self.imethods
            return imethod, self.nvalue
        elif self.method in ['PKS', 'PKNLS']:
            return self.imethod, self.omax
        # PK, PKNL
        return self.imethod, self.nvalue

    @property
    def max_id(self) -> int:
        return max(self.flutter_id.max(), self.density_flfact_id.max(),
                   self.mach_flfact_id.max(), self.rfreq_velocity_flfact_id.max())

    @parse_check
    def write_file(self, bdf_file: TextIOLike, size: int=8,
                   is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card, size = get_print_card_size(size, self.max_id)

        flutter_ids = array_str(self.flutter_id, size=size)
        densities = array_str(self.density_flfact_id, size=size)
        machs = array_str(self.mach_flfact_id, size=size)
        rfreqs = array_str(self.rfreq_velocity_flfact_id, size=size)
        epsilons = array_default_float(self.epsilon, default=1e-3, size=size, is_double=False)
        imethods = array_default_str(self.imethod, default='L', size=size)

        # mixed field based on method flag
        iomax = (self.method == 'PKS') & (self.method == 'PKNLS')
        invalue = ~iomax
        omaxs = array_float_nan(self.omax, size=size, is_double=False)
        nvalues = array_default_int(self.nvalue, default=0, size=size)
        omaxs[invalue] = nvalues[invalue]

        for flutter_id, method, density, mach, rfreq, \
            imethod, nvalue_omax, epsilon in zip(
                flutter_ids, self.method, densities, machs, rfreqs,
                imethods, omaxs, epsilons):
            #(imethod, nvalue) = self._get_repr_nvalue_omax()
            #epsilon = set_blank_if_default(self.epsilon, 0.001)
            list_fields = ['FLUTTER', flutter_id, method, density, mach,
                           rfreq, imethod, nvalue_omax, epsilon]
            bdf_file.write(print_card(list_fields))
        return



def get_mklist(mkaeros: list) -> np.ndarray:
    """Gets the (nmk, 2) array of [mach, reduced_freq] pairs from MKAERO1/MKAERO2 cards.

    Parameters
    ----------
    mkaeros : list[MKAERO1 | MKAERO2]
        list of MKAERO1/MKAERO2 card objects

    Returns
    -------
    mk_array : (nmk, 2) float ndarray
        sorted unique [mach, reduced_freq] pairs;
        empty (0,) array if no MKAERO cards

    """
    mklist = []
    for mkaero in mkaeros:
        mklist += mkaero.mklist()
    if not mklist:
        return np.array([])
    mk_array = np.array(mklist, dtype='float64')
    mk_array = np.unique(mk_array, axis=0)
    return mk_array
