from __future__ import annotations
from itertools import count
from typing import Optional, TYPE_CHECKING
import numpy as np

from pyNastran.utils.numpy_utils import integer_types

from .utils import split_filename_dollar
from pyNastran.bdf.cards.aero.dynamic_loads import Aero
from pyNastran.bdf.field_writer_8 import (
    set_blank_if_default, print_card_8, print_float_8)
from pyNastran.bdf.cards.base_card import BaseCard
from pyNastran.bdf.bdf_interface.assign_type import (
    integer, integer_or_blank, double, double_or_blank, string,
    string_or_blank, integer_or_string, integer_or_double, blank,
    integer_string_or_blank,
    # string_multifield_or_blank,
    string_multifield_dollar_int_or_blank,
)
from pyNastran.bdf.cards.aero.aero import AELINK
from pyNastran.bdf.cards.aero.static_loads import AEROS
from pyNastran.bdf.cards.aero.dynamic_loads import AERO  # MKAERO1,
from pyNastran.bdf.cards.coordinate_systems import CoordBase
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf import BDF
    from pyNastran.bdf.bdf_interface.bdf_card import BDFCard


class ACOORD(CoordBase):  # not done
    """
    Defines a general coordinate system using three rotational angles as
    functions of coordinate values in the reference coordinate system.
    The CORD3G entry is used with the MAT9 entry to orient material principal
    axes for 3-D composite analysis.

    +--------+-----+--------+--------+--------+-------+-------+--------+
    |    1   |  2  |    3   |    4   |    5   |   6   |   7   |    8   |
    +========+=====+========+========+========+=======+=======+========+
    | ACOORD |  ID | XORIGN | YORIGN | ZORIGN | DELTA | THETA |        |
    +--------+-----+--------+--------+--------+-------+-------+--------+
    | ACOORD |  10 |  250.0 |  52.5  |  15.0  |  0.0  |  0.0  |        |
    +--------+-----+--------+--------+--------+-------+-------+--------+
    """
    type = 'ACOORD'
    Type = 'R'

    @property
    def rid(self):
        return None

    def __init__(self, cid: int,
                 origin: np.ndarray,
                 delta: float,
                 theta: float,
                 comment: str=''):
        """
        Defines the ACOORD card

        Parameters
        ----------
        cid : int
            coordinate system id
        origin : list[float]
            the xyz origin
        delta : float
            pitch angle; enforced using a boundary condition
        theta : float
            roll angle; physically rolls the model
        comment : str; default=''
            a comment for the card

        """
        CoordBase.__init__(self)
        if comment:
            self.comment = comment
        self.cid = cid
        self.origin = origin
        self.delta = delta
        self.theta = theta

    @classmethod
    def add_card(cls, card: BDFCard, comment: str=''):
        """
        Adds a ACOORD card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        cid = integer(card, 1, 'cid')
        origin_x = double(card, 2, 'origin_x')
        origin_y = double(card, 3, 'origin_y')
        origin_z = double(card, 4, 'origin_z')
        origin = np.array([
            origin_x, origin_y, origin_z])
        delta = double(card, 5, 'delta (pitch)')
        theta = double(card, 6, 'theta (roll)')
        dunno = double_or_blank(card, 7, 'zero', default=0.0)
        assert dunno == 0.0, blank(card, 7, 'unknown')
        assert len(card) <= 8, f'len(ACOORD card) = {len(card):d}\ncard={card}'
        return ACOORD(cid, origin, delta, theta, comment=comment)

    def setup(self):
        self.i = np.array([1., 0., 0.])
        self.j = np.array([0., 1., 0.])
        self.k = np.array([0., 0., 1.])

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        pass

    def safe_cross_reference(self, model: BDF, xref_errors) -> None:
        self.cross_reference(model)

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        pass

    def coord_to_xyz(self, p):
        return p
        #return self.acoord_transform_to_global(p)

    def acoord_transform_to_global(self, p):
        """
        Parameters
        ----------
        p : (3,) float ndarray
            the point to transform

        .. warning:: not done, just setting up how you'd do this
        .. note::    per http://en.wikipedia.org/wiki/Euler_angles

         "This means for example that a convention named (YXZ) is the result
         of performing first an intrinsic Z rotation, followed by X and
         Y rotations, in the moving axes (Note: the order of multiplication
         of matrices is the opposite of the order in which they're
         applied to a vector)."

        """
        ct = np.cos(np.radians(self.theta))
        st = np.sin(np.radians(self.theta))
        #if rotation == 1:
            #p = self.rotation_x(ct, st) @ p
        #elif rotation == 2:
        p = self.rotation_y(ct, st) @ p
        #elif rotation == 3:
            #p = self.rotation_z(ct, st) @ p
        #else:
            #raise RuntimeError('rotation=%s rotations=%s' % (rotation, rotations))
        return p

    def rotation_x(self, ct, st):
        matrix = np.array([[1., 0., 0.],
                           [ct, 0., -st],
                           [-st, 0., ct]])
        return matrix

    def rotation_y(self, ct, st):
        matrix = np.array([[ct, 0., st],
                           [0., 1., 0.],
                           [-st, 0., ct]])
        return matrix

    def rotation_z(self, ct, st):
        matrix = np.array([[ct, st, 0.],
                           [-st, ct, 0.],
                           [0., 0., 1.]])
        return matrix

    def raw_fields(self):
        list_fields = ['ACOORD', self.cid] + list(self.origin) + [self.delta, self.theta]
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class AEROZ(Aero):
    """
    Gives basic aerodynamic parameters for unsteady aerodynamics.

    +-------+-------+-------+------+------+-------+-------+-------+
    |   1   |   2   |   3   |  4   |  5   |   6   |   7   |   8   |
    +=======+=======+=======+======+======+=======+=======+=======+
    | AEROS | ACSID | RCSID | REFC | REFB | REFS  | SYMXZ | SYMXY |
    +-------+-------+-------+------+------+-------+-------+-------+
    | AEROS |   10  |   20  | 10.  | 100. | 1000. |   1   |       |
    +-------+-------+-------+------+------+-------+-------+-------+
    """
    type = 'AEROZ'
    #_field_map = {
        #1: 'acsid', 2:'rcsid', 3:'cRef', 4:'bRef', 5:'Sref',
        #6:'symXZ', 7:'symXY',
    #}

    def __init__(self, fm_mass_unit: str, fm_length_unit: str,
                 cref: float, bref: float, sref: float,
                 flip: str='NO', acsid: int=0, rcsid: int=0,
                 sym_xz: str='NO',
                 xyz_ref: Optional[list[float]]=None,
                 comment: str=''):
        """
        Creates an AEROZ card

        Parameters
        ----------
        cref : float
            the aerodynamic chord
        bref : float
            the wing span
            for a half model, this should be the full span
            for a full model, this should be the full span
        sref : float
            the wing area
            for a half model, this should be the half area
            for a full model, this should be the full area
        acsid : int; default=0
            aerodyanmic coordinate system
            defines the direction of the wind
        rcsid : int; default=0
            coordinate system for rigid body motions
        sym_xz : str; default='NO'
            xz symmetry flag (+1=symmetry; -1=antisymmetric)
        comment : str; default=''
            a comment for the card

        """
        Aero.__init__(self)
        if comment:
            self.comment = comment

        if xyz_ref is None:
            xyz_ref = [0., 0., 0.]
        self.fm_mass_unit = fm_mass_unit
        self.fm_length_unit = fm_length_unit
        self.flip = flip

        #: Aerodynamic coordinate system identification.
        self.acsid = acsid

        #: Reference coordinate system identification for rigid body motions.
        self.rcsid = rcsid

        #: Reference chord length
        self.cref = cref

        #: Reference span
        self.bref = bref

        #: Reference wing area
        self.sref = sref

        #: Symmetry key for the aero coordinate x-z plane. See Remark 6.
        #: (Integer = +1 for symmetry, 0 for no symmetry, and -1 for antisymmetry;
        #: Default = 0)
        self.sym_xz = sym_xz

        self.xyz_ref = np.asarray(xyz_ref, dtype='float64')

        if self.acsid is None:
            self.acsid = 0
        if self.rcsid is None:
            self.rcsid = 0
        # if self.sym_xz is None:
        #     self.sym_xz = 0
        # if self.sym_xy is None:
        #     self.sym_xy = 0
        self.acsid_ref = None
        self.rcsid_ref = None

        # YES-aero=half,structure=half
        # NO-aero=full; structure=full
        # H2F-aero=full; structure=half
        assert sym_xz in ['YES', 'NO', 'H2F'], 'sym_xz=%r' % flip

        # YES-structure=left,aero=right
        assert flip in ['YES', 'NO'], f'flip={flip!r}'
        assert fm_mass_unit in ['SLIN', 'LBM', 'SLUG', 'NONE'], f'fm_mass_unit={fm_mass_unit!r}'
        assert fm_length_unit in ['IN', 'FT', 'M', 'CM', 'MM', 'NONE'], f'fm_length_unit={fm_length_unit!r}'

    @classmethod
    def add_card(cls, card: BDFCard, comment: str=''):
        """
        Adds an AEROZ card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card


        $       ACSID XZSYM FLIP FMMUNIT FMLUNIT REFC   REFB   REFS
        $+ABC   REFX  REFY  REFZ
        AEROZ   0     YES   NO   SLIN    IN       22.73 59.394 1175.8
                59.53 0.0   0.0

        """
        acsid = integer_or_blank(card, 1, 'acsid', 0)
        sym_xz = string(card, 2, 'sym_xz')
        flip = string(card, 3, 'flip')
        fm_mass_unit = string_or_blank(card, 4, 'fm_mass_unit', default='NONE')
        fm_length_unit = string(card, 5, 'fm_length_unit')

        #rcsid = integer_or_blank(card, 2, 'rcsid', 0)

        cref = double_or_blank(card, 6, 'cRef', default=1.)
        bref = double_or_blank(card, 7, 'bRef', default=1.)
        sref = double_or_blank(card, 8, 'Sref', default=1.)

        xref = double_or_blank(card, 9, 'xRef', default=0.)
        yref = double_or_blank(card, 10, 'yRef', default=0.)
        zref = double_or_blank(card, 11, 'zref', default=0.)
        xyz_ref = [xref, yref, zref]

        assert len(card) <= 12, f'len(AEROZ card) = {len(card):d}\ncard={card}'

        # faking data to not change gui
        rcsid = 0
        #sym_xy = 0
        return AEROZ(fm_mass_unit, fm_length_unit,
                     cref, bref, sref, acsid=acsid, rcsid=rcsid,
                     sym_xz=sym_xz, flip=flip, xyz_ref=xyz_ref,
                     comment=comment)

    @property
    def weight_unit(self) -> str:
        if self.fm_mass_unit == 'NONE':
            weight_unit = 'NONE'
        elif self.fm_length_unit in ['IN', 'FT']:
            weight_unit = 'LBF'
        elif self.fm_length_unit in ['MM', 'CM', 'M', 'KM']:
            weight_unit ='N'
        else:
            raise RuntimeError(f'weight_unit={self.fm_length_unit!r}')
        # UNIT_MAP = {
        #     ('SLIN', 'IN'): 'LBF',
        #     ('SLUG', 'FT'): 'LBF',
        #     ('LBF', 'IN'): 'LBF',
        # }
        # key = (self.fm_mass_unit, self.fm_length_unit)
        # return UNIT_MAP.get(key, '???')
        return weight_unit

    def Acsid(self) -> int:
        try:
            return self.acsid_ref.cid
        except AttributeError:
            return self.acsid

    def Rcsid(self) -> int:
        try:
            return self.rcsid_ref.cid
        except AttributeError:
            return self.rcsid

    # def validate(self):
    #     msg = ''
    #     if not isinstance(self.acsid, integer_types):
    #         msg += 'acsid=%s must be an integer; type=%s\n' % (self.acsid, type(self.acsid))
    #     if not isinstance(self.rcsid, integer_types):
    #         msg += 'rcsid=%s must be an integer; type=%s\n' % (self.rcsid, type(self.rcsid))
    #     if not isinstance(self.cref, float):
    #         msg += 'cref=%s must be an float; type=%s\n' % (self.cref, type(self.cref))
    #     if not isinstance(self.bref, float):
    #         msg += 'bref=%s must be an float; type=%s\n' % (self.bref, type(self.bref))
    #     if not isinstance(self.sref, float):
    #         msg += 'sref=%s must be an float; type=%s\n' % (self.sref, type(self.sref))
    #     if not isinstance(self.sym_xz, integer_types):
    #         msg += 'sym_xz=%s must be an integer; type=%s\n' % (self.sym_xz, type(self.sym_xz))
    #     if not isinstance(self.sym_xy, integer_types):
    #         msg += 'sym_xy=%s must be an integer; type=%s\n' % (self.sym_xy, type(self.sym_xy))
    #     if msg:
    #         raise TypeError('There are errors on the AEROS card:\n%s%s' % (msg, self))

    def cross_reference(self, model: BDF) -> None:
        """
        Cross reference aerodynamic coordinate system.

        Parameters
        ----------
        model : BDF
            The BDF object.

        """
        msg = ', which is required by AEROZ'
        self.acsid_ref = model.Coord(self.acsid, msg=msg)
        self.rcsid_ref = model.Coord(self.rcsid, msg=msg)

    def safe_cross_reference(self, model: BDF, xref_errors):
        """
        Safe cross reference aerodynamic coordinate system.

        Parameters
        ----------
        model : BDF
            The BDF object.

        """
        msg = ', which is required by AEROZ'
        self.acsid_ref = model.safe_coord(self.acsid, None, xref_errors, msg=msg)
        self.rcsid_ref = model.safe_coord(self.rcsid, None, xref_errors, msg=msg)

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.acsid_ref = None
        self.rcsid_ref = None

    def convert_to_zona(self, unused_model):
        #$       ACSID XZSYM FLIP FMMUNIT FMLUNIT REFC   REFB   REFS
        #$+ABC   REFX  REFY  REFZ
        #AEROZ   0     YES   NO   SLIN    IN       22.73 59.394 1175.8
                #59.53 0.0   0.0
        cref = self.cref
        bref = self.bref
        sref = self.sref
        acsid = self.acsid
        rho_ref = 1.0
        if self.sym_xz == 'NO':
            sym_xz = 0
        elif self.sym_xz == 'YES':
            sym_xz = 1
        else:
            raise NotImplementedError(self.sym_xz)
        assert sym_xz in [0, 1], sym_xz
        aeros = AEROS(cref, bref, sref, acsid=acsid, rcsid=0, sym_xz=sym_xz, sym_xy=0,
                      comment=str(self))

        velocity = 1.
        aero = AERO(velocity, cref, rho_ref, acsid=acsid, sym_xz=sym_xz, sym_xy=0,
                    comment='')
        return aeros, aero

    def raw_fields(self):
        """
        Gets the fields in their unmodified form

        Returns
        -------
        fields : list[varies]
            the fields that define the card

        """
        raise NotImplementedError()
        #list_fields = ['AEROS', self.Acsid(), self.Rcsid(), self.cref,
                       #self.bref, self.sref, self.sym_xz, self.sym_xy]
        #return list_fields

    def repr_fields(self):
        """
        Gets the fields in their simplified form

        Returns
        -------
        fields : list[varies]
          the fields that define the card

        """
        unused_sym_xz = set_blank_if_default(self.sym_xz, 0)
        unused_sym_xy = set_blank_if_default(self.sym_xy, 0)
        #$       ACSID XZSYM FLIP FMMUNIT FMLUNIT REFC   REFB   REFS
        #$+ABC   REFX  REFY  REFZ
        #AEROZ   0     YES   NO   SLIN    IN       22.73 59.394 1175.8
                #59.53 0.0   0.0

        list_fields = ['AEROZ', self.Acsid(), self.sym_xz, self.flip,
                       self.fm_mass_unit, self.fm_length_unit,
                       self.cref, self.bref, self.sref] + list(self.xyz_ref)
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class ATTACH(BaseCard):
    """
    Defines a rigid body connection between aerodynamic boxes and
    structural finite element grid points.
    """
    type = 'ATTACH'

    def __init__(self, attach_id: int, model: str,
                 setk: int, refgrid: int, comment: str=''):
        BaseCard.__init__(self)

        if comment:
            self.comment = comment
        self.attach_id = attach_id
        self.model = model
        self.setk = setk
        self.refgrid = refgrid

    @classmethod
    def add_card(cls, card: BDFCard, comment: str=''):
        # ['ATTACH', '80', '80', '3023']
        # ['ATTACH', '80', 'FFIN#1', '80', '3023']
        attach_id = integer(card, 1, 'attach_id')
        model = integer_string_or_blank(card, 2, 'model', default='NONE') # not used
        if isinstance(model, integer_types):
            setk = integer(card, 3, 'setk')
            # model = 'NONE'
            refgrid = integer(card, 3, 'refgrid')
        else:
            setk = integer(card, 3, 'setk')  # PANLST1/2/3
            refgrid = integer(card, 4, 'refgrid')
        assert len(card) <= 5, f'len(ATTACH card) = {len(card):d}\ncard={card}'
        return ATTACH(attach_id, model, setk, refgrid, comment=comment)

    def cross_reference(self, model: BDF) -> None:
        msg = f', which is required by ATTACH={self.attach_id}'
        # self.setk_ref = model.Set(self.setk, msg)

    def safe_cross_reference(self, model: BDF, xref_errors) -> None:
        self.cross_reference(model)

    def uncross_reference(self) -> None:
        return

    def repr_fields(self):
        """
        Gets the fields in their simplified form

        Returns
        -------
        fields : list[varies]
          the fields that define the card

        """
        list_fields = ['ATTACH', self.attach_id, self.model, self.setk, self.refgrid]
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class MLDPRNT(BaseCard):
    # _field_map = {
    #     1: 'sid', 2: 'mach', 3: 'q', 8: 'aeqr',
    # }

    def __init__(self, mldprnt_id: int, filename: str,
                 form: str,
                 tspnt: float, tepnt: float,
                 labels: list[str], ikeys: list[int],
                 psd_time: str='TIME', sof='NO', comment: str=''):
        BaseCard.__init__(self)
        if comment:
            self.comment = comment

        self.mldprnt_id = mldprnt_id
        self.filename = filename
        self.form = form
        self.psd_time = psd_time
        self.tspnt = tspnt
        self.tepnt = tepnt
        self.sof = sof
        self.labels = labels
        self.ikeys = ikeys
        assert form in {'TABLE', 'IDEAS', 'FEMAP', 'ESA'}, f'form={form!r}'
        assert psd_time in {'PSD', 'TIME'}, f'psd_time={psd_time!r}'
        assert sof in {'YES', 'NO', 'RFA'}, f'sof={sof!r}'

    @classmethod
    def add_card(cls, card: BDFCard, comment: str=''):
        """
        Adds a TRIM card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        # MLDPRNT IDPRNT  ---FILENM---- FORM    PSD    TSPNT TEPNT SOF
        #         LABEL1  IKEY1  LABEL2 IKEY2
        # MLDPRNT 10      HA144MLD.PLT  TABLE                      YES
        #         STATE   X      MODALX 1       EXTOUT 100
        mldprnt_id = integer(card, 1, 'mldprnt_id')
        filename = string_multifield_dollar_int_or_blank(
            card, (2, 3), 'filename', default='')
        form = string_or_blank(card, 4, 'form', default='TABLE')
        psd_time = string_or_blank(card, 5, 'psd/time', default='TIME')
        tspnt = double_or_blank(card, 6, 'TSPNT')
        tepnt = double_or_blank(card, 7, 'TEPNT')
        sof = string_or_blank(card, 8, 'SOF', default='NO')

        nfields_left = len(card) - 9
        assert nfields_left % 2 == 0, nfields_left
        assert nfields_left // 2 > 0, nfields_left

        labels = []
        ikeys = []
        j = 1
        for ifield in range(9, len(card), 2):
            label = string(card, ifield, f'label{j}')
            ikey = integer_or_string(card, ifield+1, f'ikey{j}')
            labels.append(label)
            ikeys.append(ikey)
            j += 1
        assert len(card) > 9, f'len(MLDPRNT card) = {len(card):d}\ncard={card}'
        return MLDPRNT(mldprnt_id, filename, form,
                       tspnt, tepnt, labels, ikeys,
                       psd_time=psd_time, sof=sof, comment=comment)

    # def validate(self):
    #     assert self.true_g in ['TRUE', 'G'], 'true_g=%r' % self.true_g

    def cross_reference(self, model: BDF) -> None:
        pass

    def safe_cross_reference(self, model: BDF, xref_errors):
        self.cross_reference(model)

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
        filenamea, filenameb = split_filename_dollar(self.filename)
        list_fields = [
            'MLDPRNT', self.mldprnt_id, filenamea, filenameb,
            self.form, self.psd_time, self.tspnt, self.tepnt, self.sof]
        for label, ikey in zip(self.labels, self.ikeys):
            list_fields.append(label)
            list_fields.append(ikey)
        return list_fields

    def repr_fields(self):
        list_fields = self.raw_fields()
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class MLDSTAT(BaseCard):
    type = 'MLDSTAT'
    # _field_map = {
    #     1: 'sid', 2: 'mach', 3: 'q', 8: 'aeqr',
    # }
    def __init__(self, mldstat_id: int, mldtrim_id: int,
                 transform,
                 filename: str,
                 states: list[str], values: list[float],
                 dx_tox: str='YES',
                 state_space_arr: str='',
                 state_space_brr: str='',
                 comment: str=''):
        BaseCard.__init__(self)
        if comment:
            self.comment = comment

        self.mldstat_id = mldstat_id
        self.mldtrim_id = mldtrim_id
        self.transform = transform
        self.filename = filename
        self.states = states
        self.values = values
        self.dx_tox = dx_tox
        self.state_space_arr = state_space_arr
        self.state_space_brr = state_space_brr
        self.mldtrim_ref = None
        self.ssa_ref = None
        self.ssb_ref = None

    @classmethod
    def add_card(cls, card: BDFCard, comment: str=''):
        """
        Adds a MLDSTAT card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        # MLDSTAT IDSTAT IDTRIM TRNSFM ARR BRR DXTOX FILENM CONT
        #         STATE1 INITIAL
        mldstat_id = integer(card, 1, 'mldstat_id')
        mldtrim_id = integer(card, 2, 'mldtrim_id')
        field3 = card.field(3)
        transform = str(field3) if field3 is not None else ''
        field4 = card.field(4)
        state_space_arr = str(field4) if field4 is not None else ''
        state_space_brr = ''
        field5 = card.field(5)
        if field5 is None or field5 == '':
            dx_tox = 'YES'
        elif isinstance(field5, int):
            dx_tox = 'NO' if field5 == 0 else 'YES'
        else:
            dx_tox = str(field5).strip()
        if dx_tox not in {'YES', 'NO'}:
            dx_tox = 'NO' if dx_tox == '0' else 'YES'

        # filename = ''
        filename = string_multifield_dollar_int_or_blank(
            card, (6, 7), 'filename', default='')
        i = 9
        j = 1
        states = []
        values = []
        while i < len(card):
            field_state = card.field(i)
            if field_state is None:
                i += 2
                j += 1
                continue
            state = str(field_state)
            field_val = card.field(i + 1)
            if field_val is None:
                value = 0.0
            elif isinstance(field_val, (int, float)):
                value = float(field_val)
            else:
                try:
                    value = float(field_val)
                except (ValueError, TypeError):
                    value = 0.0
            states.append(state)
            values.append(value)
            i += 2
            j += 1
        return MLDSTAT(mldstat_id, mldtrim_id, transform, filename,
                       states, values, dx_tox=dx_tox,
                       state_space_arr=state_space_arr,
                       state_space_brr=state_space_brr,
                       comment=comment)

    # def validate(self):
    #     assert self.true_g in ['TRUE', 'G'], 'true_g=%r' % self.true_g

    def cross_reference(self, model: BDF) -> None:
        # msg = f', which is required by MLDSTAT={self.mldstat_id}\n{str(self)}'
        if self.mldtrim_id:
            self.mldtrim_ref = model.zaero.mldtrim[self.mldtrim_id]
        if self.state_space_arr:
            self.ssa_ref = model.dmi[self.state_space_arr]
        if self.state_space_arr:
            self.ssb_ref = model.dmi[self.state_space_brr]

    def safe_cross_reference(self, model: BDF, xref_errors):
        self.cross_reference(model)

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.mldtrim_ref = None
        self.ssa_ref = None
        self.ssb_ref = None

    def raw_fields(self):
        """
        Gets the fields in their unmodified form

        Returns
        -------
        fields : list[varies]
            the fields that define the card

        """
        filename_a, filename_b = split_filename_dollar(self.filename)
        list_fields = [
            'MLDSTAT', self.mldstat_id, self.mldtrim_id, self.transform,
            self.state_space_arr, self.state_space_brr, self.dx_tox, filename_a, filename_b,]
        for state, value in zip(self.states, self.values):
            list_fields.append(state)
            list_fields.append(value)
        return list_fields

    def repr_fields(self):
        list_fields = self.raw_fields()
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        msg = (
            f'MLDSTAT {self.mldstat_id:<8d}{self.mldtrim_id:<8d}{self.transform:<8s}'
            f'{self.state_space_arr:<8s}{self.state_space_brr:<8s}'
            f'{self.dx_tox:<8s}{self.filename:>16s}\n        '
        )
        for istate, state, value in zip(count(), self.states, self.values):
            msg += f'{state:>8s}{print_float_8(value)}'
            if istate > 0 and istate % 4 == 0:
                msg += '\n        '
        return self.comment + msg.rstrip() + '\n'


class MINSTAT(BaseCard):
    type = 'MINSTAT'
    # _field_map = {
    #     1: 'sid', 2: 'mach', 3: 'q', 8: 'aeqr',
    # }

    def __init__(self, minstat_id, aerolag_id, itmax, apcid,
                 pweight_id, dinit_id, klist_id, min_inp,
                 print_flag, save_flag, filename,
                 msmod, aerolag_gust_id=0, dinit_gust_id=0,
                 comment: str=''):
        BaseCard.__init__(self)
        if comment:
            self.comment = comment

        self.minstat_id = minstat_id
        self.aerolag_id = aerolag_id
        self.itmax = itmax
        self.apcid = apcid
        self.pweight_id = pweight_id
        self.dinit_id = dinit_id
        self.klist_id = klist_id
        self.min_inp = min_inp
        self.print_flag = print_flag
        self.save_flag = save_flag
        self.filename = filename
        self.msmod = msmod
        self.aerolag_gust_id = aerolag_gust_id
        self.dinit_gust_id = dinit_gust_id

    @classmethod
    def add_card(cls, card: BDFCard, comment: str=''):
        """
        Adds a MINSTAT card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        # MINSTAT SID   LAGID ITMAX  APCID APWID DMATID KLIST   MININP
        #         PRINT SAVE  FILENM       MSMOD LAGIDG DMATIDG
        # MINSTAT 20    10    100    21    22    23     0       3
        #         -2    SAVE  APPROX.DAT
        minstat_id = integer(card, 1, 'mldstat_id')
        aerolag_id = integer(card, 2, 'LAGID, aerolag_id')
        itmax = integer_or_blank(card, 3, 'itmax', default=None)
        apcid = integer_or_blank(card, 4, 'apcid', default=0)
        pweight_id = integer_or_blank(card, 5, 'apwid, pweight_id', default=0)
        dinit_id = integer_or_blank(card, 6, 'DINIT, dinit_id', default=0)
        klist_id = integer_or_blank(card, 7, 'klist_id', default=0)
        min_inp = integer_or_blank(card, 8, 'min_inp', default=0)
        print_flag = integer_or_blank(card, 9, 'print_flag', default=0)
        save_flag = string_or_blank(card, 10, 'save_flag', default='')
        filename = string_multifield_dollar_int_or_blank(
            card, (11, 12), 'filename', default='')
        msmod = integer_or_blank(card, 13, 'msmod', default=0)
        aerolag_gust_id = integer_or_blank(card, 14, 'LAGIDG, aerolag2_id', default=0)
        dinit_gust_id = integer_or_blank(card, 15, 'DMATIDG, dinit_id', default=0)

        return MINSTAT(minstat_id, aerolag_id, itmax, apcid,
                       pweight_id, dinit_id, klist_id, min_inp,
                       print_flag, save_flag, filename,
                       msmod, aerolag_gust_id, dinit_gust_id,
                       comment=comment)

    # def validate(self):
    #     assert self.true_g in ['TRUE', 'G'], 'true_g=%r' % self.true_g

    def cross_reference(self, model: BDF) -> None:
        pass

    def safe_cross_reference(self, model: BDF, xref_errors):
        self.cross_reference(model)

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
        list_fields = [
            'MINSTAT',
            self.minstat_id, self.aerolag_id, self.itmax, self.apcid,
            self.pweight_id, self.dinit_id, self.klist_id, self.min_inp,
            self.print_flag, self.save_flag, self.filename,
            self.msmod, self.aerolag_gust_id, self.dinit_gust_id,]
        return list_fields

    def repr_fields(self):
        list_fields = self.raw_fields()
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class MLDCOMD(BaseCard):
    type = 'MLDCOMD'

    # _field_map = {
    #     1: 'sid', 2: 'mach', 3: 'q', 8: 'aeqr',
    # }
    def __init__(self, mldcomd_id: int,
                 extinp_ids: list[int], table_ids: list[int],
                 comment: str=''):
        BaseCard.__init__(self)
        if comment:
            self.comment = comment

        self.mldcomd_id = mldcomd_id
        self.extinp_ids = extinp_ids
        self.table_ids = table_ids
        self.extinps_ref = None
        self.tables_ref = None

    @classmethod
    def add_card(cls, card: BDFCard, comment: str=''):
        """
        Adds a TRIM card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        # MLDCOMD SETID CONT
        #         EXTINP1 IDTAB1 EXTINP2 IDTAB2
        # MLDTRIM GRAVITY NZ THKCAM MODAL
        mldcomd_id = integer(card, 1, 'mldcomd_id')

        i = 9
        j = 1
        extinp_ids = []
        table_ids = []
        while i < len(card):
            extinp_id = integer_or_blank(card, i, f'extinp_id{j}')
            if extinp_id is None:
                i += 2
                j += 1
                continue
            table_id = integer_or_blank(card, i+1, f'table_id{j}', default=0)
            extinp_ids.append(extinp_id)
            table_ids.append(table_id)
            i += 2
            j += 1
        return MLDCOMD(mldcomd_id, extinp_ids, table_ids, comment=comment)

    # def validate(self):
    #     assert self.true_g in ['TRUE', 'G'], 'true_g=%r' % self.true_g

    def cross_reference(self, model: BDF) -> None:
        extinps_ref = []
        tables_ref = []
        zaero = model.zaero
        msg = f', which is required by MLDCOMD={self.mldcomd_id}\n{str(self)}'
        for extinp_id, tabled_id in zip(self.extinp_ids, self.table_ids):
            extinp = zaero.extinp[extinp_id]
            extinp.cross_reference(model, self)
            tabled = model.TableD(tabled_id, msg)
            extinps_ref.append(extinp)
            tables_ref.append(tabled)
        self.extinps_ref = extinps_ref
        self.tables_ref = tables_ref

    def safe_cross_reference(self, model: BDF, xref_errors):
        self.cross_reference(model)

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.extinps_ref = []
        self.tables_ref = []

    def raw_fields(self):
        """
        Gets the fields in their unmodified form

        Returns
        -------
        fields : list[varies]
            the fields that define the card

        """
        list_fields = ['MLDCOMD', self.mldcomd_id, None, None, None, None, None, None, None]
        for key, value in zip(self.extinp_ids, self.table_ids):
            list_fields.append(key)
            list_fields.append(value)
        return list_fields

    def repr_fields(self):
        list_fields = self.raw_fields()
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class MLDTIME(BaseCard):
    type = 'MLDTIME'

    # _field_map = {
    #     1: 'sid', 2: 'mach', 3: 'q', 8: 'aeqr',
    # }
    def __init__(self, mldtime_id: int, tstart: float, tend: float, dt: float,
                 out_dt: int, print_flag: int, method: str,
                 comment: str=''):
        BaseCard.__init__(self)
        if comment:
            self.comment = comment

        self.mldtime_id = mldtime_id
        self.tstart = tstart
        self.tend = tend
        self.dt = dt
        self.out_dt = out_dt
        self.print_flag = print_flag
        self.method = method

    @classmethod
    def add_card(cls, card: BDFCard, comment: str=''):
        """
        Adds a TRIM card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        # MLDTIME IDTIME TSTART TEND DT   OUTDT PRINT METHOD
        # MLDTIME 1      0.0    2.0  0.01 5     -1
        mldtime_id = integer(card, 1, 'mldcomd_id')
        tstart = double(card, 2, 'tstart')
        tend = double(card, 3, 'tend')
        dt = double(card, 4, 'dt')
        out_dt = integer_or_blank(card, 5, 'out_dt', default=1)
        print_flag = integer_or_blank(card, 6, 'print_flag', default=-1)
        method = integer_or_blank(card, 7, 'method')

        return MLDTIME(mldtime_id, tstart, tend, dt, out_dt,
                       print_flag, method, comment=comment)

    # def validate(self):
    #     assert self.true_g in ['TRUE', 'G'], 'true_g=%r' % self.true_g

    def cross_reference(self, model: BDF) -> None:
        pass

    def safe_cross_reference(self, model: BDF, xref_errors):
        self.cross_reference(model)

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
        list_fields = [
            'MLDTIME', self.mldtime_id, self.tstart, self.tend,
            self.dt, self.out_dt,
            self.print_flag, self.method]
        return list_fields

    def repr_fields(self):
        list_fields = self.raw_fields()
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class MLDTRIM(BaseCard):
    type = 'MLDTRIM'

    # _field_map = {
    #     1: 'sid', 2: 'mach', 3: 'q', 8: 'aeqr',
    # }
    def __init__(self, mldtrim_id, gravity: float, nz: float,
                 thkcam: str,
                 trim_vars: list[str], values: list[float],
                 modal_dmi: str='', comment: str=''):
        BaseCard.__init__(self)
        if comment:
            self.comment = comment

        self.mldtrim_id = mldtrim_id
        self.gravity = gravity
        self.nz = nz
        self.thkcam = thkcam
        self.modal_dmi = modal_dmi
        self.trim_vars = trim_vars
        self.values = values

    @classmethod
    def add_card(cls, card: BDFCard, comment: str=''):
        """
        Adds a TRIM card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        # MLDTRIM GRAVITY NZ THKCAM MODAL
        mldtrim_id = integer(card, 1, 'mldtrim_id')
        gravity = double(card, 2, 'gravity')
        nz = double(card, 3, 'nz')
        thkcam = string(card, 4, 'thkcam')
        modal_dmi = string_or_blank(card, 5, 'modal_dmi', default='')

        i = 9
        j = 1
        trim_vars = []
        values = []
        while i < len(card):
            trim_var = string(card, i, f'trim_var{j}')
            value = double(card, i+1, f'value{j}')
            trim_vars.append(trim_var)
            values.append(value)
            i += 2
            j += 1

        return MLDTRIM(mldtrim_id, gravity, nz, thkcam,
                       trim_vars, values, modal_dmi=modal_dmi, comment=comment)

    # def validate(self):
    #     assert self.true_g in ['TRUE', 'G'], 'true_g=%r' % self.true_g

    def cross_reference(self, model: BDF) -> None:
        pass

    def safe_cross_reference(self, model: BDF, xref_errors):
        self.cross_reference(model)

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
        list_fields = [
            'MLDTRIM', self.mldtrim_id,  self.gravity,
            self.nz, self.thkcam, self.modal_dmi, None, None, None]
        for trim_var, value in zip(self.trim_vars, self.values):
            list_fields.append(trim_var)
            list_fields.append(value)
        return list_fields

    def repr_fields(self):
        list_fields = self.raw_fields()
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class EXTFILE(BaseCard):
    type = 'EXTFILE'

    def __init__(self, extfile_id: int, filename: str, comment: str=''):
        BaseCard.__init__(self)

        if comment:
            self.comment = comment
        self.extfile_id = extfile_id
        self.filename = filename

    @classmethod
    def add_card(cls, card: BDFCard, comment: str=''):
        extfile_id = integer(card, 1, 'extfile_id')
        assert len(card) == 3, f'len(EXTFILE card) = {len(card):d}\ncard={card}'
        filename = card.field(2)
        assert len(filename) > 0, card
        return EXTFILE(extfile_id, filename, comment=comment)

    def cross_reference(self, model: BDF) -> None:
        pass

    def safe_cross_reference(self, model: BDF, xref_errors) -> None:
        self.cross_reference(model)

    def uncross_reference(self) -> None:
        return

    def repr_fields(self):
        """
        Gets the fields in their simplified form

        Returns
        -------
        fields : list[varies]
          the fields that define the card

        """
        list_fields = ['EXTFILE', self.extfile_id, self.filename]
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        assert len(str(self.extfile_id)) < 7, self.extfile_id
        assert len(self.filename) < 55, self.filename
        msg =  f'EXTFILE {self.extfile_id:<8d}{self.filename}\n'
        card = self.repr_fields()
        return self.comment + msg

    def __repr__(self):
        return self.write_card(size=8)


class SPLINE0(BaseCard):
    """Defines a rigid body spline (infinite plate spline with no elasticity).

    SPLINE0 EID     MODEL   CP      SETK
    """
    type = 'SPLINE0'

    def __init__(self, eid: int, model_name: str, cp: int,
                 setk: int, comment: str = ''):
        BaseCard.__init__(self)
        if comment:
            self.comment = comment
        self.eid = eid
        self.model_name = model_name
        self.cp = cp
        self.setk = setk

    @classmethod
    def add_card(cls, card, comment: str = ''):
        eid = integer(card, 1, 'eid')
        model_name = string_or_blank(card, 2, 'model', default='')
        cp = integer_or_blank(card, 3, 'cp', default=0)
        setk = integer_or_blank(card, 4, 'setk', default=0)
        return SPLINE0(eid, model_name, cp, setk, comment=comment)

    def raw_fields(self):
        return ['SPLINE0', self.eid, self.model_name, self.cp, self.setk]

    def write_card(self, size: int = 8, is_double: bool = False) -> str:
        card = self.raw_fields()
        return self.comment + print_card_8(card)

    def __repr__(self):
        return self.write_card(size=8)


class PBODY7(BaseCard):
    """Defines properties for a BODY7 aerodynamic body element.

    PBODY7 PID     IPBODY  ...
    """
    type = 'PBODY7'

    def __init__(self, pid: int, ipbody: int, fields: list, comment: str = ''):
        BaseCard.__init__(self)
        if comment:
            self.comment = comment
        self.pid = pid
        self.ipbody = ipbody
        self.fields = fields

    @classmethod
    def add_card(cls, card, comment: str = ''):
        pid = integer(card, 1, 'pid')
        ipbody = integer_or_blank(card, 2, 'ipbody', default=0)
        fields = []
        for i in range(3, len(card)):
            val = card.field(i)
            fields.append(val)
        return PBODY7(pid, ipbody, fields, comment=comment)

    def raw_fields(self):
        return ['PBODY7', self.pid, self.ipbody] + self.fields

    def write_card(self, size: int = 8, is_double: bool = False) -> str:
        card = self.raw_fields()
        return self.comment + print_card_8(card)

    def __repr__(self):
        return self.write_card(size=8)


class TRIMFLT(BaseCard):
    """Defines trim conditions for flutter analysis on a deformed shape.

    TRIMFLT IDFLT   TITLA   ALPHA   ...
    """
    type = 'TRIMFLT'

    def __init__(self, trimflt_id: int, title: str, alpha: float,
                 fields: list, comment: str = ''):
        BaseCard.__init__(self)
        if comment:
            self.comment = comment
        self.trimflt_id = trimflt_id
        self.title = title
        self.alpha = alpha
        self.fields = fields

    @classmethod
    def add_card(cls, card, comment: str = ''):
        trimflt_id = integer(card, 1, 'trimflt_id')
        # Field 2 can be a string title or an integer reference ID
        field2 = card.field(2)
        if field2 is None or field2 == '':
            title = ''
        elif isinstance(field2, int):
            title = str(field2)
        else:
            title = str(field2)
        alpha = double_or_blank(card, 3, 'alpha', default=0.0)
        fields = []
        for i in range(4, len(card)):
            fields.append(card.field(i))
        return TRIMFLT(trimflt_id, title, alpha, fields, comment=comment)

    def raw_fields(self):
        return ['TRIMFLT', self.trimflt_id, self.title, self.alpha] + self.fields

    def write_card(self, size: int = 8, is_double: bool = False) -> str:
        card = self.raw_fields()
        return self.comment + print_card_8(card)

    def __repr__(self):
        return self.write_card(size=8)


class CONMLST(BaseCard):
    """Defines concentrated mass additions for a flutter/gust subcase.

    CONMLST SID                                                         +C1
    +C1     FACTOR1 CID1    FACTOR2 CID2    ...

    Each (FACTOR, CID) pair scales the referenced CONM1/CONM2 card.
    """
    type = 'CONMLST'

    def __init__(self, conmlst_id: int, factors: list[float],
                 conm_ids: list[int], comment: str = ''):
        BaseCard.__init__(self)
        if comment:
            self.comment = comment
        self.conmlst_id = conmlst_id
        self.factors = factors
        self.conm_ids = conm_ids

    @classmethod
    def add_card(cls, card, comment: str = ''):
        conmlst_id = integer(card, 1, 'conmlst_id')
        # Continuation fields: pairs of (factor, conm_id)
        factors = []
        conm_ids = []
        i = 9  # first continuation field
        j = 1
        while i < len(card):
            factor = double_or_blank(card, i, f'factor{j}')
            if factor is None:
                break
            conm_id = integer_or_blank(card, i + 1, f'conm_id{j}', default=0)
            factors.append(factor)
            conm_ids.append(conm_id)
            i += 2
            j += 1
        return CONMLST(conmlst_id, factors, conm_ids, comment=comment)

    def raw_fields(self):
        fields = ['CONMLST', self.conmlst_id, None, None, None, None, None, None, None]
        for factor, cid in zip(self.factors, self.conm_ids):
            fields.extend([factor, cid])
        return fields

    def write_card(self, size: int = 8, is_double: bool = False) -> str:
        card = self.raw_fields()
        return self.comment + print_card_8(card)

    def __repr__(self):
        return self.write_card(size=8)


class CPFACT(BaseCard):
    """Defines a pressure scaling factor for aerodynamic panels.

    CPFACT  SID     IDMK    SYM     COMP    TYPE    PANLST  FACTOR1 FACTOR2 +CP
    +CP     STRIP1  ...
    """
    type = 'CPFACT'

    def __init__(self, cpfact_id: int, idmk: int, sym: str,
                 comp: str, cptype: str, panlst: str,
                 factor1: float, factor2: float,
                 strips: list[int], comment: str = ''):
        BaseCard.__init__(self)
        if comment:
            self.comment = comment
        self.cpfact_id = cpfact_id
        self.idmk = idmk
        self.sym = sym
        self.comp = comp
        self.cptype = cptype
        self.panlst = panlst
        self.factor1 = factor1
        self.factor2 = factor2
        self.strips = strips

    @classmethod
    def add_card(cls, card, comment: str = ''):
        cpfact_id = integer(card, 1, 'cpfact_id')
        idmk = integer(card, 2, 'idmk')
        sym = string_or_blank(card, 3, 'sym', default='BOTH')
        comp = string_or_blank(card, 4, 'comp', default='ALL')
        cptype = string_or_blank(card, 5, 'cptype', default='FEM')
        panlst = string_or_blank(card, 6, 'panlst', default='ALL')
        factor1 = double_or_blank(card, 7, 'factor1', default=1.0)
        factor2 = double_or_blank(card, 8, 'factor2', default=0.0)
        strips = []
        for i in range(9, len(card)):
            val = integer_or_blank(card, i, f'strip{i - 8}')
            if val is not None:
                strips.append(val)
        return CPFACT(cpfact_id, idmk, sym, comp, cptype, panlst,
                      factor1, factor2, strips, comment=comment)

    def raw_fields(self):
        fields = ['CPFACT', self.cpfact_id, self.idmk, self.sym,
                  self.comp, self.cptype, self.panlst,
                  self.factor1, self.factor2] + self.strips
        return fields

    def write_card(self, size: int = 8, is_double: bool = False) -> str:
        card = self.raw_fields()
        return self.comment + print_card_8(card)

    def __repr__(self):
        return self.write_card(size=8)


def _make_generic_zaero_card_class(card_type_name: str):
    """Factory that creates a card class for a specific ZAERO card type."""

    class _Card(BaseCard):
        type = card_type_name

        def __init__(self, card_id, fields, comment=''):
            BaseCard.__init__(self)
            if comment:
                self.comment = comment
            self.card_id = card_id
            self._fields = fields

        @classmethod
        def add_card(cls, card, comment: str = ''):
            card_id = 0
            raw = card.field(1)
            if isinstance(raw, int):
                card_id = raw
            elif isinstance(raw, str):
                try:
                    card_id = int(raw)
                except (ValueError, TypeError):
                    pass
            fields = [card.field(i) for i in range(1, len(card))]
            return cls(card_id, fields, comment=comment)

        def raw_fields(self):
            return [self.type] + self._fields

        def write_card(self, size: int = 8, is_double: bool = False) -> str:
            card = self.raw_fields()
            return self.comment + print_card_8(card)

        def __repr__(self):
            return self.write_card(size=8)

    _Card.__name__ = card_type_name
    _Card.__qualname__ = card_type_name
    return _Card


# Generate card classes for all remaining ZAERO cards
ASEOUT = _make_generic_zaero_card_class('ASEOUT')
class CMARGIN(BaseCard):
    """Defines gain/phase margin analysis parameters for ASE.

    CMARGIN SID     GMHI    GMLO    PMHI    PMLO                        +CM
    +CM     DF      PRINT   NROOT
    """

    type = 'CMARGIN'

    def __init__(self, cmargin_id: int, gm_high: float, gm_low: float,
                 pm_high: float, pm_low: float, df: float = 1e-4,
                 print_flag: int = 0, nroot: int = 0, comment: str = ''):
        BaseCard.__init__(self)
        if comment:
            self.comment = comment
        self.cmargin_id = cmargin_id
        self.gm_high = gm_high
        self.gm_low = gm_low
        self.pm_high = pm_high
        self.pm_low = pm_low
        self.df = df
        self.print_flag = print_flag
        self.nroot = nroot

    @classmethod
    def add_card(cls, card, comment: str = ''):
        cmargin_id = integer(card, 1, 'cmargin_id')
        gm_high = double_or_blank(card, 2, 'gm_high', default=50.0)
        gm_low = double_or_blank(card, 3, 'gm_low', default=-50.0)
        pm_high = double_or_blank(card, 4, 'pm_high', default=60.0)
        pm_low = double_or_blank(card, 5, 'pm_low', default=-60.0)
        df = double_or_blank(card, 9, 'df', default=1e-4)
        print_flag = integer_or_blank(card, 10, 'print_flag', default=0)
        nroot = integer_or_blank(card, 11, 'nroot', default=0)
        return CMARGIN(cmargin_id, gm_high, gm_low, pm_high, pm_low,
                       df, print_flag, nroot, comment=comment)

    def raw_fields(self):
        return ['CMARGIN', self.cmargin_id, self.gm_high, self.gm_low,
                self.pm_high, self.pm_low, None, None, None,
                self.df, self.print_flag, self.nroot]

    def write_card(self, size: int = 8, is_double: bool = False) -> str:
        card = self.raw_fields()
        return self.comment + print_card_8(card)

    def __repr__(self):
        return self.write_card(size=8)
OUTPUT4 = _make_generic_zaero_card_class('OUTPUT4')
CROSPSD = _make_generic_zaero_card_class('CROSPSD')
FOILSEC = _make_generic_zaero_card_class('FOILSEC')
GENGUST = _make_generic_zaero_card_class('GENGUST')
MFTGUST = _make_generic_zaero_card_class('MFTGUST')
DMIS = _make_generic_zaero_card_class('DMIS')
CELLWNG = _make_generic_zaero_card_class('CELLWNG')
CELLBOX = _make_generic_zaero_card_class('CELLBOX')
INPCFD = _make_generic_zaero_card_class('INPCFD')
OMITCFD = _make_generic_zaero_card_class('OMITCFD')
WT1AJJ = _make_generic_zaero_card_class('WT1AJJ')
WT1FRC = _make_generic_zaero_card_class('WT1FRC')
WT2AJJ = _make_generic_zaero_card_class('WT2AJJ')
WTUCP = _make_generic_zaero_card_class('WTUCP')
TRIMOBJ = _make_generic_zaero_card_class('TRIMOBJ')
TRIMCON = _make_generic_zaero_card_class('TRIMCON')
APCNSND = _make_generic_zaero_card_class('APCNSND')
APCNSCP = _make_generic_zaero_card_class('APCNSCP')


class APCONST(BaseCard):
    """Defines constraints for rational aerodynamic approximation (RFA).

    APCONST SID     DA0     DA1     DA2     NRP     NCP     FR(1)   FC(1)
    """
    type = 'APCONST'

    def __init__(self, sid: int, da0: int, da1: int, da2: int,
                 nrp: int, ncp: int,
                 fr_values: list[float], fc_values: list[float],
                 comment: str = ''):
        BaseCard.__init__(self)
        if comment:
            self.comment = comment
        self.sid = sid
        self.da0 = da0
        self.da1 = da1
        self.da2 = da2
        self.nrp = nrp
        self.ncp = ncp
        self.fr_values = fr_values
        self.fc_values = fc_values

    @classmethod
    def add_card(cls, card, comment: str = ''):
        sid = integer(card, 1, 'sid')
        da0 = integer_or_blank(card, 2, 'da0', default=0)
        da1 = integer_or_blank(card, 3, 'da1', default=0)
        da2 = integer_or_blank(card, 4, 'da2', default=0)
        nrp = integer_or_blank(card, 5, 'nrp', default=0)
        ncp = integer_or_blank(card, 6, 'ncp', default=0)
        fr_values = []
        fc_values = []
        i = 7
        for j in range(max(nrp, 1)):
            if i >= len(card):
                break
            fr = integer_or_double(card, i, f'fr{j + 1}')
            fr_values.append(fr)
            i += 1
        for j in range(max(ncp, 1)):
            if i >= len(card):
                break
            fc = integer_or_double(card, i, f'fc{j + 1}')
            fc_values.append(fc)
            i += 1
        return APCONST(sid, da0, da1, da2, nrp, ncp,
                       fr_values, fc_values, comment=comment)

    def raw_fields(self):
        fields = ['APCONST', self.sid, self.da0, self.da1, self.da2,
                  self.nrp, self.ncp] + self.fr_values + self.fc_values
        return fields

    def write_card(self, size: int = 8, is_double: bool = False) -> str:
        card = self.raw_fields()
        return self.comment + print_card_8(card)

    def __repr__(self):
        return self.write_card(size=8)
