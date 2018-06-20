# coding: utf-8
"""
All ZONA aero cards are defined in this file.  This includes:
 * TRIM

All cards are BaseCard objects.
"""
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
from itertools import count
from six import string_types
import numpy as np

from pyNastran.bdf.cards.aero.dynamic_loads import Aero
from pyNastran.utils import integer_types
from pyNastran.bdf.field_writer_8 import set_blank_if_default, print_card_8
from pyNastran.bdf.cards.base_card import BaseCard
from pyNastran.bdf.bdf_interface.assign_type import (
    integer, integer_or_blank, double, double_or_blank, string,
    string_or_blank, integer_or_string, double_or_string, blank,
)
from pyNastran.bdf.cards.aero.aero import Spline
from pyNastran.bdf.cards.aero.utils import elements_from_quad, points_elements_from_quad_points
from pyNastran.bdf.cards.coordinate_systems import Coord

class ACOORD(Coord):  # not done
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

    def __init__(self, cid, origin, delta, theta, comment=''):
        """
        Defines the CORD3G card

        Parameters
        ----------
        cid : int
            coordinate system id
        origin : List[float]
            the xyz origin
        delta : float
            pitch angle
        theta : float
            roll angle
        comment : str; default=''
            a comment for the card

        """
        Coord.__init__(self)
        if comment:
            self.comment = comment
        self.cid = cid
        self.origin = origin
        self.delta = delta
        self.theta = theta

    @classmethod
    def add_card(cls, card, comment=''):
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
        origin = [origin_x, origin_y, origin_z]
        delta = double(card, 5, 'delta')
        theta = double(card, 6, 'theta')
        assert len(card) <= 7, 'len(ACOORD card) = %i\ncard=%s' % (len(card), card)
        return ACOORD(cid, origin, delta, theta, comment=comment)

    def setup(self):
        self.i = np.array([1., 0., 0.])
        self.j = np.array([0., 1., 0.])
        self.k = np.array([0., 0., 1.])

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        pass

    def uncross_reference(self):
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
            #p = np.dot(self.rotation_x(ct, st), p)
        #elif rotation == 2:
        p = np.dot(self.rotation_y(ct, st), p)
        #elif rotation == 3:
            #p = np.dot(self.rotation_z(ct, st), p)
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
        list_fields = ['ACOORD', self.cid] + self.origin + [self.delta, self.theta]
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class AESURFZ(BaseCard):
    """
    Specifies an aerodynamic control surface for aeroservoelastic, static
    aeroelastic/trim analysis, or the transient response analysis.

    +---------+--------+-------+-------+-------+--------+--------+--------+--------+
    |    1    |   2    |   3   |   4   |   5   |   6    |    7   |   8    |   9    |
    +=========+========+=======+=======+=======+========+========+========+========+
    | AESURFZ | LABEL  |  TYPE |  CID  |  SETK |  SETG  |  ACTID |        |        |
    +---------+--------+-------+-------+-------+--------+--------+--------+--------+
    | AESURFZ | RUDDER |  ASYM |   1   |   10  |   20   |    0   |        |        |
    +---------+--------+-------+-------+-------+--------+--------+--------+--------+
    """
    type = 'AESURFZ'

    @property
    def aesid(self):
        return self.label
    @property
    def alid1_ref(self):
        return None

    def __init__(self, label, surface_type, cid, setk, setg, actuator_tf,
                 comment=''):
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

        #: Controller name.
        self.label = label
        self.surface_type = surface_type
        self.setk = setk
        self.setg = setg
        self.actuator_tf = actuator_tf
        setk, setg, actuator_tf

        #: Identification number of a rectangular coordinate system with a
        #: y-axis that defines the hinge line of the control surface
        #: component.
        self.cid = cid

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
        label = string(card, 1, 'label')
        surface_type = string(card, 2, 'TYPE')
        cid = integer(card, 3, 'CID')
        setk = integer(card, 4, 'SETK') # PANLST1, PANLST2, PANLST3
        setg = integer(card, 5, 'SETG') # SET1, SETADD
        actuator_tf = integer_or_blank(card, 6, 'ACTID') # ACTU card
        assert len(card) <= 6, 'len(AESURFZ card) = %i\ncard=%s' % (len(card), card)
        assert surface_type in ['SYM', 'ANTISYM', 'ASYM']
        return AESURFZ(label, surface_type, cid, setk, setg, actuator_tf, comment=comment)

    def Cid(self):
        if self.cid_ref is not None:
            return self.cid_ref.cid
        return self.cid

    def aelist_id1(self):
        if self.alid1_ref is not None:
            return self.alid1_ref.sid
        return self.alid1

    def aelist_id2(self):
        if self.alid2_ref is not None:
            return self.alid2_ref.sid
        return self.alid2

    def AELIST_id1(self):
        self.deprecated('AESURF.AELIST_id1()', 'AESURF.aelist_id1()', '1.1')
        return self.aelist_id1()

    def AELIST_id2(self):
        self.deprecated('AESURF.AELIST_id2()', 'AESURF.aelist_id2()', '1.1')
        return self.aelist_id2()

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        self.cid_ref = model.Coord(self.cid)

        #self.alid1_ref = model.AELIST(self.alid1)
        #if self.alid2:
            #self.alid2_ref = model.AELIST(self.alid2)
        #if self.tqllim is not None:
            #self.tqllim_ref = model.TableD(self.tqllim)
        #if self.tqulim is not None:
            #self.tqulim_ref = model.TableD(self.tqulim)
        self.setk_ref = model.zona.panlsts[self.setk]
        self.setk_ref.cross_reference(model)
        self.aero_element_ids = self.setk_ref.aero_element_ids

    def safe_cross_reference(self, model, xref_errors):
        msg = ', which is required by AESURF aesid=%s' % self.aesid
        self.cid_ref = model.safe_coord(self.cid, self.aesid, xref_errors, msg=msg)
        #if self.cid2 is not None:
            #self.cid2_ref = model.safe_coord(self.cid2, self.aesid, xref_errors, msg=msg)
        #try:
            #self.alid1_ref = model.AELIST(self.alid1)
        #except KeyError:
            #pass
        #if self.alid2:
            #try:
                #self.alid2_ref = model.AELIST(self.alid2)
            #except KeyError:
                #pass
        #if self.tqllim is not None:
            #try:
                #self.tqllim_ref = model.TableD(self.tqllim)
            #except KeyError:
                #pass
        #if self.tqulim is not None:
            #try:
                #self.tqulim_ref = model.TableD(self.tqulim)
            #except KeyError:
                #pass
        self.setk_ref = model.zona.panlsts[self.setk]
        self.setk_ref.cross_reference(model)
        self.aero_element_ids = self.setk_ref.aero_element_ids

    def uncross_reference(self):
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
        list_fields = ['AESURFZ', self.label, self.surface_type, self.cid,
                       self.setk, self.setg, self.actuator_tf]
        return list_fields

    def repr_fields(self):
        """
        Gets the fields in their simplified form

        Returns
        -------
        fields : List[int/float/str]
            the fields that define the card

        """
        return self.raw_fields()

    def write_card(self, size=8, is_double=False):
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
    type = 'AEROS'
    _field_map = {
        1: 'acsid', 2:'rcsid', 3:'cRef', 4:'bRef', 5:'Sref',
        6:'symXZ', 7:'symXY',
    }

    def __init__(self, fm_mass_unit, fm_length_unit,
                 cref, bref, sref,
                 flip='NO', acsid=0, rcsid=0, sym_xz=0, xyz_ref=None, comment=''):
        """
        Creates an AEROS card

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
        sym_xz : int; default=0
            xz symmetry flag (+1=symmetry; -1=antisymmetric)
        comment : str; default=''
            a comment for the card

        """
        Aero.__init__(self)
        if comment:
            self.comment = comment

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

        self.xyz_ref = xyz_ref

        if self.acsid is None:
            self.acsid = 0
        if self.rcsid is None:
            self.rcsid = 0
        if self.sym_xz is None:
            self.sym_xz = 0
        if self.sym_xy is None:
            self.sym_xy = 0
        self.acsid_ref = None
        self.rcsid_ref = None

    def Acsid(self):
        try:
            return self.acsid_ref.cid
        except AttributeError:
            return self.acsid

    def Rcsid(self):
        try:
            return self.rcsid_ref.cid
        except AttributeError:
            return self.rcsid

    #def validate(self):
        #msg = ''
        #if not isinstance(self.acsid, integer_types):
            #msg += 'acsid=%s must be an integer; type=%s\n' % (self.acsid, type(self.acsid))
        #if not isinstance(self.rcsid, integer_types):
            #msg += 'rcsid=%s must be an integer; type=%s\n' % (self.rcsid, type(self.rcsid))
        #if not isinstance(self.cref, float):
            #msg += 'cref=%s must be an float; type=%s\n' % (self.cref, type(self.cref))
        #if not isinstance(self.bref, float):
            #msg += 'bref=%s must be an float; type=%s\n' % (self.bref, type(self.bref))
        #if not isinstance(self.sref, float):
            #msg += 'sref=%s must be an float; type=%s\n' % (self.sref, type(self.sref))
        #if not isinstance(self.sym_xz, integer_types):
            #msg += 'sym_xz=%s must be an integer; type=%s\n' % (self.sym_xz, type(self.sym_xz))
        #if not isinstance(self.sym_xy, integer_types):
            #msg += 'sym_xy=%s must be an integer; type=%s\n' % (self.sym_xy, type(self.sym_xy))
        #if msg:
            #raise TypeError('There are errors on the AEROS card:\n%s%s' % (msg, self))

    def cross_reference(self, model):
        """
        Cross refernece aerodynamic coordinate system.

        Parameters
        ----------
        model : BDF
            The BDF object.

        """
        msg = ', which is required by AEROS'
        self.acsid_ref = model.Coord(self.acsid, msg=msg)
        self.rcsid_ref = model.Coord(self.rcsid, msg=msg)

    def safe_cross_reference(self, model, xref_errors):
        """
        Safe cross refernece aerodynamic coordinate system.

        Parameters
        ----------
        model : BDF
            The BDF object.

        """
        msg = ', which is required by AEROS'
        self.acsid_ref = model.safe_coord(self.acsid, None, xref_errors, msg=msg)
        self.rcsid_ref = model.safe_coord(self.rcsid, None, xref_errors, msg=msg)

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds an AEROS card from ``BDF.add_card(...)``

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
        fm_mass_unit = string(card, 4, 'fm_mass_unit')
        fm_length_unit = string(card, 5, 'fm_length_unit')

        # YES-aero=half,structure=half
        # NO-aero=full; structure=full
        # H2F-aero=full; structure=half
        assert sym_xz in ['YES', 'NO', 'H2F'], 'sym_xz=%r' % flip


        # YES-structure=left,aero=right
        assert flip in ['YES', 'NO'], 'flip=%r' % flip
        assert fm_mass_unit in ['SLIN', 'LBM'], 'fm_mass_unit=%r' % fm_mass_unit
        assert fm_length_unit in ['IN'], 'fm_length_unit=%r' % fm_length_unit

        #rcsid = integer_or_blank(card, 2, 'rcsid', 0)

        cref = double(card, 6, 'cRef')
        bref = double(card, 7, 'bRef')
        sref = double(card, 8, 'Sref')

        xref = double(card, 9, 'xRef')
        yref = double(card, 10, 'yRef')
        zref = double(card, 11, 'zref')
        xyz_ref = [xref, yref, zref]

        assert len(card) <= 12, 'len(AEROZ card) = %i\ncard=%s' % (len(card), card)
        rcsid = 0
        sym_xy = 0
        return AEROZ(fm_mass_unit, fm_length_unit,
                     cref, bref, sref, acsid=acsid, rcsid=rcsid,
                     sym_xz=sym_xz, flip=flip, xyz_ref=xyz_ref,
                     comment=comment)

    def uncross_reference(self):
        self.acsid_ref = None
        self.rcsid_ref = None

    def update(self, maps):
        """
        maps = {
            'coord' : cid_map,
        }

        """
        cid_map = maps['coord']
        self.acsid = cid_map[self.acsid]
        self.rcsid = cid_map[self.rcsid]

    def raw_fields(self):
        """
        Gets the fields in their unmodified form

        Returns
        -------
        fields : list[varies]
            the fields that define the card

        """
        asdf
        #list_fields = ['AEROS', self.Acsid(), self.Rcsid(), self.cref,
                       #self.bref, self.sref, self.sym_xz, self.sym_xy]
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
        #$       ACSID XZSYM FLIP FMMUNIT FMLUNIT REFC   REFB   REFS
        #$+ABC   REFX  REFY  REFZ
        #AEROZ   0     YES   NO   SLIN    IN       22.73 59.394 1175.8
                #59.53 0.0   0.0

        list_fields = ['AEROZ', self.Acsid(), self.sym_xz, self.flip,
                       self.fm_mass_unit, self.fm_length_unit,
                       self.cref, self.bref, self.sref] + list(self.xyz_ref)
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        return self.comment + print_card_8(card)

class PANLST3(Spline):
    """
    Defines a set of aerodynamic boxes by the LABEL entry in CAERO7 or BODY7
    bulk data cards.

    +---------+------+-------+-------+------+------+----+-----+-------+
    |    1    |  2   |   3   |   4   |  5   |   6  |  7 |  8  |   9   |
    +=========+======+=======+=======+======+======+====+=====+=======+
    | SPLINE1 | EID  | MODEL |  CP   | SETK | SETG | DZ | EPS |       |
    +---------+------+-------+-------+------+------+----+-----+-------+
    | SPLINE1 | 100  |       |       |  1   |  10  | 0. |     |       |
    +---------+------+-------+-------+------+------+----+-----+-------+
    PANLST3 SETID LABEL1 LABEL2 LABEL3 … -etc- …
    PANLST3 100 WING HTAIL

    PANLST3 is referred to by SPLINEi, ATTACH, LOADMOD, CPFACT, JETFRC, and/or
    AESURFZ bulk data card.
    """
    type = 'PANLST3'

    def __init__(self, eid, panel_groups, comment=''):
        """
        Creates a PANLST3 card

        Parameters
        ----------
        eid : int
            spline id
        comment : str; default=''
            a comment for the card

        """
        # https://www.zonatech.com/Documentation/ZAERO_9.2_Users_3rd_Ed.pdf
        Spline.__init__(self)
        if comment:
            self.comment = comment

        self.eid = eid
        self.panel_groups = panel_groups # points to CAERO7 / BODY7
        self.aero_element_ids = []

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
        group_id = 1
        panel_groups = []
        for ifield in range(2, len(card)):
            name = string(card, ifield, 'group_%i'%  (group_id))
            panel_groups.append(name)
        return PANLST3(eid, panel_groups, comment=comment)

    def cross_reference(self, model):
        msg = ', which is required by PANLST3 eid=%s' % self.eid
        #self.nodes_ref = model.Nodes(self.nodes, msg=msg)
        caero_refs = []
        aero_element_ids = []
        for caero_label in self.panel_groups:
            caero_eid = model.zona.caero_to_name_map[caero_label]
            caero_ref = model.CAero(caero_eid, msg=msg)
            caero_refs.append(caero_ref)
            eid = caero_ref.eid
            npanels = caero_ref.npanels
            if npanels == 0:
                model.log.warning('skipping PANLST3 because there are 0 panels in:\n%r' % caero_ref)
                continue
            aero_element_ids2 = range(eid, eid + npanels)
            assert len(aero_element_ids2) == npanels, npanels
            aero_element_ids += aero_element_ids2
        self.caero_refs = caero_refs
        self.aero_element_ids = aero_element_ids

    def safe_cross_reference(self, model, xref_errors):
        self.cross_reference(model)

    def raw_fields(self):
        list_fields = ['PANLST3', self.eid] + self.panel_groups
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        return self.comment + print_card_8(card)

class BODY7(BaseCard):
    """
    Defines an aerodynamic body macroelement of a body-like component.

    +--------+-----+-----+----+-----+------+-----+------+------+
    |    1   |  2  |  3  |  4 |  5  |   6  |   7 |   8  |  9   |
    +========+=====+=====+====+=====+======+=====+======+======+
    | CAERO2 | EID | PID | CP | NSB | NINT | LSB | LINT | IGID |
    +--------+-----+-----+----+-----+------+-----+------+------+
    |        | X1  |  Y1 | Z1 | X12 |      |     |      |      |
    +--------+-----+-----+----+-----+------+-----+------+------+
    BODY7 BID LABEL IPBODY7 ACOORD NSEG IDMESH1 IDMESH2 IDMESH3 CONT
    CONT IDMESH4 -etc
    BODY7 4 BODY 2 8 4 20 21 22 +BC
    +BC 23
    """
    type = 'BODY7'

    def __init__(self, eid, label, pid, nseg, idmeshes, acoord=0, comment=''):
        """
        Defines a BODY7 card, which defines a slender body
        (e.g., fuselage/wingtip tank).

        Parameters
        ----------
        eid : int
            body id
        label : str
            An arbitrary character string used to define the body.
        pid  : int; default=0
            Identification number of PBODY7 bulk data card
            (specifying body wake and/or inlet aerodynamic boxes)
        acoord : int; default=0
            Identification number of ACOORD bulk data card
            (specifying body center line location and orientation)
        nseg : int
            Number of body segments
        idmeshes : List[int]
            Identification number of SEGMESH bulk data card (specifying body segments).
        comment : str; default=''
            a comment for the card

        """
        BaseCard.__init__(self)

        if comment:
            self.comment = comment

        #: Element identification number
        self.eid = eid
        self.label = label

        #: Property identification number of a PAERO7 entry.
        self.pid = pid
        self.nseg = nseg
        self.idmeshes = idmeshes
        self.acoord = acoord
        self.pid_ref = None
        self.acoord_ref = None

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a BODY7 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        eid = integer(card, 1, 'eid')
        label = string(card, 2, 'label')
        pid = integer_or_blank(card, 3, 'pid')
        acoord = integer_or_blank(card, 4, 'acoord', 0)
        nseg = integer_or_blank(card, 5, 'nseg')

        idmeshes = []
        for i, ifield in enumerate(range(6, len(card))):
            segmesh = integer(card, ifield, 'idmesh_%i' % (i+1))
            idmeshes.append(segmesh)
        assert len(card) <= 13, 'len(BODY7 card) = %i\ncard=%s' % (len(card), card)
        return BODY7(eid, label, pid, nseg, idmeshes, acoord=acoord, comment=comment)

    def ACoord(self):
        if self.acoord_ref is not None:
            return self.acoord_ref.cid
        return self.acoord

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

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        msg = ', which is required by BODY7 eid=%s' % self.eid
        self.segmesh_refs = []
        for segmesh_id in self.idmeshes:
            segmesh_ref = model.PAero(segmesh_id, msg=msg)  # links to SEGMESH/PAERO7
            self.segmesh_refs.append(segmesh_ref)

        #if self.pid is not None:
            #self.pid_ref = model.PAero(self.pid, msg=msg)  # links to PAERO7
        #self.acoord_ref = model.Coord(self.acoord, msg=msg)
        #if self.nsb == 0:
            #self.lsb_ref = model.AEFact(self.lsb, msg=msg)
        #if self.nint == 0:
            #self.lint_ref = model.AEFact(self.lint, msg=msg)
        if self.acoord:
            self.acoord_ref = model.Coord(self.acoord, msg=msg)
        #self.ascid_ref = model.Acsid(msg=msg)
        self.ascid_ref = model.Coord(0, msg=msg)

    def safe_cross_reference(self, model, xref_errors, debug=False):
        self.cross_reference(model)

    def uncross_reference(self):
        self.pid = self.Pid()
        self.acoord = self.ACoord()
        self.pid_ref = None
        self.acoord_ref = None

    def get_points(self):
        """creates a 1D representation of the CAERO2"""
        p1 = self.cp_ref.transform_node_to_global(self.p1)
        p2 = p1 + self.ascid_ref.transform_vector_to_global(np.array([self.x12, 0., 0.]))

        #print("x12 = %s" % self.x12)
        #print("pcaero[%s] = %s" % (self.eid, [p1,p2]))
        return [p1, p2]

    @property
    def npanels(self):
        npanels = 0
        for segmesh in self.segmesh_refs:
            for itype in segmesh.itypes:
                if itype not in [3]:
                    return 0
            nx = segmesh.naxial
            ny = segmesh.nradial
            npanels += nx * ny
        return npanels

    def get_points_elements_3d(self):
        """
        Gets the points/elements in 3d space as CQUAD4s
        The idea is that this is used by the GUI to display CAERO panels.

        TODO: doesn't support the aero coordinate system

        """
        paero2 = self.pid_ref
        xyz = []
        element = []
        nelements = 0
        for segmesh in self.segmesh_refs:
            #print(segmesh)
            try:
                xyzi, elementi = self._get_points_elements_3di(segmesh)
            except NotImplementedError:
                return None, None
            xyz.append(xyzi)
            element.append(elementi + nelements)
            nelements += elementi.shape[0]

        xyzs = np.vstack(xyz)
        elements = np.vstack(element)
        return xyzs, elements

    def _get_points_elements_3di(self, segmesh):
        """
        points (nchord, nspan) float ndarray; might be backwards???
            the points
        elements (nquads, 4) int ndarray
            series of quad elements
            nquads = (nchord-1) * (nspan-1)
        """
        lengths_y = []
        lengths_z = []
        #print(segmesh.get_stats())
        nx = segmesh.naxial
        ny = segmesh.nradial

        xs = []
        ys = []
        zs = []
        for itype, x, y, z, camber, idy_ref, idz_ref in zip(segmesh.itypes,
                                                            segmesh.xs, segmesh.ys, segmesh.zs,
                                                            segmesh.cambers,
                                                            segmesh.idys_ref, segmesh.idzs_ref):
            if itype not in [3]:
                msg = ('only itype=3 is supported; not %r (1=body of revolution, '
                       '2=elliptical body, 3=arbitrary body)' % itype)
                raise NotImplementedError(msg)
            assert camber == 0.0, 'camber is not supported'
            #print(x, y, z, idy_ref, idz_ref)
            #print(idy_ref.get_stats())
            ynodes = idy_ref.fractions
            znodes = idz_ref.fractions
            assert len(ynodes) == len(znodes), 'len(ynodes)=%s len(znodes)=%s' % (len(ynodes), len(znodes))
            nnodes = len(ynodes)

            origin_x, origin_y, origin_z = self.acoord_ref.origin
            x_offset = origin_x + x
            y_offset = origin_y + y
            z_offset = origin_z + z
            xs.append([x_offset] * nnodes)
            ys.append(y_offset + ynodes)
            zs.append(z_offset + znodes)

        xyz = np.vstack([
            np.hstack(xs),
            np.hstack(ys),
            np.hstack(zs),
        ]).T
        #print(xyz)
        #print('ynodes = %s' % ynodes)
        #print('znodes = %s' % znodes)
        elements = elements_from_quad(nx, ny, dtype='int32')
        return xyz, elements

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
        list_fields = ['BODY7', self.eid, self.label, self.Pid(), self.ACoord(),
                        self.nseg] + self.idmeshes
        return list_fields

    def repr_fields(self):
        """
        Gets the fields in their simplified form

        Returns
        -------
        fields : list
            The fields that define the card

        """
        return self.raw_fields()

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        return self.comment + print_card_8(card)

class SEGMESH(BaseCard):
    """
    Defines a grid system for a body segment; referenced by the BODY7 bulk data card.

    +---------+--------+-------+------+---------+-------+------+------+-----+
    |    1    |    2   |   3   |   4  |    5    |   6   |   7  |   8  |  9  |
    +=========+========+=======+======+=========+=======+======+======+=====+
    | SEGMESH | IDMESH | NAXIS | NRAD | NOSERAD | IAXIS |      |      |     |
    +---------+--------+-------+------+---------+-------+------+------+-----+
    |         | ITYPE1 |   X1  | CAM1 |   YR1   |  ZR1  | IDY1 | IDZ1 |     |
    +---------+--------+-------+------+---------+-------+------+------+-----+
    |         | ITYPE2 |   X2  | CAM2 |   YR2   |  ZR2  | IDY2 | IDZ2 |     |
    +---------+--------+-------+------+---------+-------+------+------+-----+
    |         | ITYPE3 |   X3  | CAM3 |   YR3   |  ZR3  | IDY3 | IDZ3 |     |
    +---------+--------+-------+------+---------+-------+------+------+-----+
    | SEGMESH |    2   |  3    |  6   |         |       |      |      |     |
    +---------+--------+-------+------+---------+-------+------+------+-----+
    |         |    1   |  0.0  |  0.0 | 0.0     |       |      |      |     |
    +---------+--------+-------+------+---------+-------+------+------+-----+
    |         |    1   |  1.0  |  0.0 | 0.5     |       |      |      |     |
    +---------+--------+-------+------+---------+-------+------+------+-----+
    |         |    3   |  2.0  |      |         |       |  103 | 104  |     |
    +---------+--------+-------+------+---------+-------+------+------+-----+

    """
    type = 'SEGMESH'

    @property
    def pid(self):
        return self.segmesh_id

    def __init__(self, segmesh_id, naxial, nradial, nose_radius, iaxis,
                 itypes, xs, cambers, ys, zs, idys, idzs, comment=''):
        """
        Defines a SEGMESH card, which is similar to a PARO2 and defines additional
        slender body parameters.

        Parameters
        ----------
        DMESH : int
            Body segment mesh identification number.
        NAXIS : int
            Number of axial stations (i.e., divisions) of the segment. (min=2).
        NRAD : int
            Number of circumferential points of the segment (min=3).
        NOSERAD : float
            Nose radius of blunt body.
            NOSERAD is active only if ZONA7U (Hypersonic Aerodynamic Method)
            is used (the METHOD entry of the MKAEROZ Bulk Data equals 2 or –2).
            Furthermore, NOSERAD is used only if the SEGMESH bulk data card is
            the first segment defined in the BODY7 bulk data card.
        IAXIS : int
            The index of the axial station where the blunt nose ends.
            IAXIS is active only if ZONA7U (Hypersonic Aerodynamic
            Method) is used.
        ITYPEi : int
            Type of input used to define the circumferential box cuts
            - 1 body of revolution
            - 2 elliptical body
            - 3 arbitrary body
        Xi : List[float]
            X-location of the axial station; Xi must be in ascending
            order. (i.e., Xi+1 > Xi)
        cambers : List[float]
            Body camber at the Xi axial station. (Real)
        YRi : List[float]
            Body cross-sectional radius if ITYPEi = 1 or the semi-axis length
            of the elliptical body parallel to the Y-axis if ITYPEi=2.
        ZRi : List[float]
            The semi-axis length of the elliptical body parallel to the Z-axis.
            Used only if ITYPEi=2. (Real)
        IDYi : int
            Identification number of AEFACT bulk data card that specifies
            NRAD number of the Y-coordinate locations of the circumferential
            points at the Xi axial station. Use only if ITYPEi=3.
        IDZi : int
            Identification number of AEFACT bulk data card that specifies
            NRAD number of the Z-coordinate locations of the circumferential
            points at the Xi axial station. Use only if ITYPEi=3.
        eid : int
            body id
        label : str
            An arbitrary character string used to define the body.
        pid  : int; default=0
            Identification number of PBODY7 bulk data card
            (specifying body wake and/or inlet aerodynamic boxes)
        acoord : int; default=0
            Identification number of ACOORD bulk data card
            (specifying body center line location and orientation)
        nseg : int
            Number of body segments
        IDMESHi : List[int]
            Identification number of SEGMESH bulk data card (specifying body segments).
        comment : str; default=''
            a comment for the card

        """
        BaseCard.__init__(self)
        self.segmesh_id = segmesh_id
        self.naxial = naxial
        self.nradial = nradial
        self.nose_radius = nose_radius
        self.iaxis = iaxis
        self.itypes = itypes
        self.cambers = cambers
        self.xs = xs
        self.ys = ys
        self.zs = zs
        self.idys = idys
        self.idzs = idzs

        self.idys_ref = None
        self.idzs_ref = None

        #self.cp_ref = None
        #self.lint_ref = None
        #self.lsb_ref = None
        #self.ascid_ref = None

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a SEGMESH card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        segmesh_id = integer(card, 1, 'segmesh_id')
        naxial = integer(card, 2, 'naxial')
        nradial = integer(card, 3, 'nradial')
        nose_radius = double_or_blank(card, 4, 'nose_radius')
        iaxis = integer_or_blank(card, 5, 'iaxis')

        itypes = []
        xs = []
        ys = []
        zs = []
        cambers = []
        idys = []
        idzs = []
        for counter, ifield in enumerate(range(9, len(card), 8)):
            itype = integer(card, ifield, 'itype%i' % (counter+1))
            x = double_or_blank(card, ifield+1, 'itype%i' % (counter+1), 0.)
            camber = double_or_blank(card, ifield+2, 'camber%i' % (counter+1), 0.)
            y = double_or_blank(card, ifield+3, 'y%i' % (counter+1), 0.)
            z = double_or_blank(card, ifield+4, 'z%i' % (counter+1), 0.)
            idy = integer_or_blank(card, ifield+5, 'idy%i' % (counter+1))
            idz = integer_or_blank(card, ifield+6, 'idz%i' % (counter+1))
            itypes.append(itype)
            xs.append(x)
            ys.append(y)
            zs.append(z)
            cambers.append(camber)
            idys.append(idy)
            idzs.append(idz)
        assert len(itypes) == naxial, 'naxial=%s nradial=%s len(itypes)=%s' % (naxial, nradial, len(itypes))
        return SEGMESH(segmesh_id, naxial, nradial, nose_radius, iaxis,
                       itypes, xs, cambers, ys, zs, idys, idzs)

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

    def cross_reference(self, model):
        msg = ', which is required by SEGMESH eid=%s' % self.pid
        idys_ref = []
        idzs_ref = []
        for idy in self.idys:
            idy_ref = None
            if idy is not None and isinstance(idy, integer_types):
                idy_ref = model.AEFact(idy, msg=msg)
            idys_ref.append(idy_ref)

        for idz in self.idzs:
            idz_ref = None
            if idz is not None and isinstance(idz, integer_types):
                idz_ref = model.AEFact(idz, msg=msg)
            idzs_ref.append(idz_ref)
        self.idys_ref = idys_ref
        self.idzs_ref = idzs_ref
        #print(self.idys_ref)

    def safe_cross_reference(self, model, xref_errors, debug=False):
        return self.cross_reference(model)

    def uncross_reference(self):
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

        #Rs = []
        assert len(radii) == (nx + 1), 'len(radii)=%s nx=%s' % (len(radii), nx)
        if len(xstation) != (nx + 1):
            msg = 'len(xstation)=%s nx=%s\nxstation=%s\n%s' % (
                len(xstation), nx, xstation, str(self))
            raise RuntimeError(msg)

        xs = []
        ys = []
        zs = []
        yzs = []
        for i, xi, yi, zi, radius in zip(count(), xstation, ystation, zstation, radii):
            #print('  station=%s xi=%.4f radius=%s' % (i, xi, radius))
            yz = self.create_ellipse(aspect_ratio, radius)
            yzs.append(yz)
            try:
                y = yz[:, 0] + yi
                z = yz[:, 1] + zi
            except ValueError:
                print('yz = %s' % yz)
                print('yz.shape = %s' % str(yz.shape))
                print('dy = %s' % dy)
                print('dz = %s' % dz)
                raise
            ntheta = yz.shape[0]
            x = np.ones(ntheta) * xi
            xs.append(x)
            ys.append(y)
            zs.append(z)
            #Rs.append(np.sqrt(y**2 + z**2))
        #print('yz.shape=%s xs.shape=%s' % (str(np.array(yzs).shape), str(np.array(xs).shape)))
        #xyz = np.hstack([yzs, xs])
        xs = np.array(xs)
        ys = np.array(ys)
        zs = np.array(zs)
        try:
            xyz = np.vstack([
                np.hstack(xs),
                np.hstack(ys),
                np.hstack(zs),
            ]).T + p1
        except:
            print('xs =', xs.shape)
            print('ys =', ys.shape)
            print('zs =', zs.shape)
            raise

        #R = np.hstack(Rs)
        #print('xyz.shape =', xyz.shape)
        #print('xyz =', xyz)
        #print('R =', R)

        ny = ntheta
        elems = elements_from_quad(nx+1, ny)
        #print('elems =\n', elems)
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
        list_fields = [
            'SEGMESH', self.segmesh_id, self.naxial, self.nradial, self.nose_radius,
            self.iaxis, None, None, None]
        for itype, x, camber, y, z, idy, idz in zip(self.itypes, self.xs, self.cambers,
                                                    self.ys, self.zs, self.idys, self.idzs):
            list_fields += [itype, x, camber, y, z, idy, idz, None]
        return list_fields

    def repr_fields(self):
        """
        Gets the fields in their simplified form

        Returns
        -------
        fields : list
            The fields that define the card

        """
        return self.raw_fields()

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        return self.comment + print_card_8(card)

class CAERO7(BaseCard):
    """
    Totally wrong...

    Defines an aerodynamic macro element (panel) in terms of two leading edge
    locations and side chords. This is used for Doublet-Lattice theory for
    subsonic aerodynamics and the ZONA51 theory for supersonic aerodynamics.

    +--------+-----+-------+--------+-------+--------+--------+--------+---------+
    |   1    |  2  |   3   |   4    |   5   |   6    |    7   |   8    |    9    |
    +========+=====+=======+========+=======+========+========+========+=========+
    | CAERO7 | WID | LABEL | ACOORD | NSPAN | NCHORD |  LSPAN |  ZTAIC | PAFOIL7 |
    +--------+-----+-------+--------+-------+--------+--------+--------+---------+
    |        | XRL |  YRL  |   ZRL  |  RCH  |  LRCHD | ATTCHR | ACORDR |         |
    +--------+-----+-------+--------+-------+--------+--------+--------+---------+
    |        | XTL |  YTL  |   ZTL  |  TCH  |  LTCHD | ATTCHT | ACORDT |         |
    +--------+-----+-------+--------+-------+--------+--------+--------+---------+

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
    igid : int
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
    type = 'CAERO7'
    _field_map = {
        1: 'sid', 2:'pid', 3:'cp', 4:'nspan', 5:'nchord',
        6:'lspan', 7:'lchord', 8:'igid', 12:'x12', 16:'x43',
    }

    def __init__(self, eid, label, p1, x12, p4, x43,
                 cp=0, nspan=0, nchord=0, lspan=0, p_airfoil=None, ztaic=None, comment=''):
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
        if nspan is None:
            nspan = 0
        if nchord is None:
            nchord = 0
        p1 = np.asarray(p1)
        p4 = np.asarray(p4)

        if comment:
            self.comment = comment
        #: Element identification number
        self.eid = eid

        #: Property identification number of a PAERO2 entry.
        self.label = label

        #: Coordinate system for locating point 1.
        self.cp = cp
        self.nspan = nspan
        self.nchord = nchord
        self.lspan = lspan
        self.p_airfoil = p_airfoil
        self.p1 = p1
        self.x12 = x12
        self.p4 = p4
        self.x43 = x43
        self.ztaic = None

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
        if self.nspan == 0 and self.lspan == 0:
            msg += 'NSPAN or LSPAN must be greater than 0; nspan=%r nlspan=%s\n' % (
                self.nspan, self.lspan)
            is_failed = True
        if self.nspan <= 0:
            msg += 'NSPAN must be greater than 0; nspan=%r\n' % (
                self.nspan)
            is_failed = True

        if self.nchord <= 0:
            msg += 'NCHORD must be greater than 0; nchord=%r\n' % (
                self.nchord)
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

        CAERO7  100101  RIBW1           2       25
                998.904  39.821 230.687 1298.159        310001
                1121.821 61.134 236.560 1175.243
        CAERO7  100201  RIBW2           2       25
                1121.821 61.134 236.560 1175.243
                1244.258 84.704 243.625 1052.805

        """
        eid = integer(card, 1, 'eid')
        name = string(card, 2, 'name')
        cp = integer_or_blank(card, 3, 'cp', 0)
        nspan = integer_or_blank(card, 4, 'nspan', 0)
        nchord = integer_or_blank(card, 5, 'nchord', 0)
        lspan = integer_or_blank(card, 6, 'aefact_lchord', 0)
        if lspan:
            lspan = 0
        ztaic = integer_or_blank(card, 7, 'ztaic')
        p_airfoil = integer_or_blank(card, 8, 'aefact')
        #assert cp == 0
        #igroup = integer(card, 8, 'igid')

        x1 = double_or_blank(card, 9, 'x1', 0.0)
        y1 = double_or_blank(card, 10, 'y1', 0.0)
        z1 = double_or_blank(card, 11, 'z1', 0.0)
        p1 = np.array([x1, y1, z1])
        x12 = double_or_blank(card, 12, 'x12', 0.)
        lchord_root = integer_or_blank(card, 13, 'lchord_root')
        attach_root = integer_or_blank(card, 14, 'attach_root')
        achord_root = integer_or_blank(card, 15, 'achord_root')

        x4 = double_or_blank(card, 17, 'x4', 0.0)
        y4 = double_or_blank(card, 18, 'y4', 0.0)
        z4 = double_or_blank(card, 19, 'z4', 0.0)
        p4 = np.array([x4, y4, z4])
        x43 = double_or_blank(card, 20, 'x43', 0.)
        lchord_tip = integer_or_blank(card, 21, 'lchord_tip')
        attach_tip = integer_or_blank(card, 22, 'attach_tip')
        achord_tip = integer_or_blank(card, 23, 'achord_tip')

        assert len(card) <= 23, 'len(CAERO7 card) = %i\ncard=%s' % (len(card), card)
        return CAERO7(eid, name, p1, x12, p4, x43,
                      cp=cp, nspan=nspan, nchord=nchord, lspan=lspan, p_airfoil=p_airfoil,
                      comment=comment)

    def flip_normal(self):
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
        i = 0
        try:
            self.box_ids = np.arange(self.eid, self.eid + npanels, dtype=dtype).reshape(nspan, nchord).T
        except OverflowError:
            if dtype == 'int64':
                # we already tried int64
                msg = 'eid=%s ichord=%s ispan=%s nchord=%s' % (
                    self.eid, ichord, ispan, nchord)
                raise OverflowError(msg)
            self._init_ids(dtype='int64')

    def Cp(self):
        if self.cp_ref is not None:
            return self.cp_ref.cid
        return self.cp

    #def Pid(self):
        #if self.pid_ref is not None:
            #return self.pid_ref.pid
        #return self.pid

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        msg = ', which is required by CAERO1 eid=%s' % self.eid
        #self.pid_ref = model.PAero(self.pid, msg=msg)
        self.cp_ref = model.Coord(self.cp, msg=msg)
        self.ascid_ref = model.Acsid(msg=msg)

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
        #try:
            #self.pid_ref = model.PAero(self.pid, msg=msg)
        #except KeyError:
            #pass

        self.cp_ref = model.safe_coord(self.cp, self.eid, xref_errors, msg=msg)
        self.ascid_ref = model.safe_acsid(msg=msg)

        if self.nchord == 0:
            assert isinstance(self.lchord, integer_types), self.lchord
            self.lchord_ref = model.safe_aefact(self.lchord, self.eid, xref_errors, msg)

        if self.nspan == 0:
            assert isinstance(self.lspan, integer_types), self.lspan
            self.lspan_ref = model.safe_aefact(self.lspan, self.eid, xref_errors, msg)

        self._init_ids()

    def uncross_reference(self):
        #self.pid = self.Pid()
        self.cp = self.Cp()
        self.lchord = self.get_LChord()
        self.lspan = self.get_LSpan()
        #self.pid_ref = None
        self.cp_ref = None
        self.lchord_ref = None
        self.lspan_ref = None
        self.ascid_ref = None

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

    def get_points(self):
        """
        Get the 4 corner points for the CAERO card

        Returns
        -------
        p1234 : (4, 3) list
             List of 4 corner points in the global frame

        """
        if self.cp_ref is None and self.cp == 0:
            p1 = self.p1
            p4 = self.p4
        else:
            p1 = self.cp_ref.transform_node_to_global(self.p1)
            p4 = self.cp_ref.transform_node_to_global(self.p4)

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
        list_fields = (['CAERO1', self.eid, self.name, self.Cp(), self.nspan,
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
        #lchord = set_blank_if_default(self.get_LChord(), 0)
        lspan = set_blank_if_default(self.get_LSpan(), 0)
        lchord = 0
        lspan = 0
        list_fields = (
            ['CAERO7', self.eid, self.label, cp, nspan, nchord, lspan, self.ztaic, self.p_airfoil,] +
            list(self.p1) + [self.x12, None, None, None, None] +
            list(self.p4) + [self.x43, None, None, None, None])
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class TRIM(BaseCard):
    """
    Specifies constraints for aeroelastic trim variables.

    +------+--------+------+--------+--------+-----+--------+-----+----------+
    |   1  |   2    |   3  |    4   |    5   |  6  |    7   |  8  |     9    |
    +======+========+======+========+========+=====+========+=====+==========+
    | TRIM |   ID   | MACH |    Q   | LABEL1 | UX1 | LABEL2 | UX2 | IS_RIGID |
    +------+--------+------+--------+--------+-----+--------+-----+----------+
    |      | LABEL3 |  UX3 | LABEL4 |   UX4  | ... |        |     |          |
    +------+--------+------+--------+--------+-----+--------+-----+----------+
    """
    type = 'TRIM_ZONA'
    _field_map = {
        1: 'sid', 2:'mach', 3:'q', 8:'aeqr',
    }

    def __init__(self, sid, mkaeroz, q, cg, true_g, nxyz, pqr, loadset,
                 labels, uxs, comment=''):
        """
        Creates a TRIM card for a static aero (144) analysis.

        Parameters
        ----------
        sid : int
            the trim id; referenced by the Case Control TRIM field
        q : float
            dynamic pressure
        labels : List[str]
            names of the fixed variables
        uxs : List[float]
            values corresponding to labels
        comment : str; default=''
            a comment for the card

        """
        BaseCard.__init__(self)
        if comment:
            self.comment = comment
        #: Trim set identification number. (Integer > 0)
        self.sid = sid
        self.mkaeroz = mkaeroz
        #: Dynamic pressure. (Real > 0.0)
        self.q = q

        self.cg = cg
        self.nxyz = nxyz
        self.true_g = true_g
        self.pqr = pqr
        self.loadset = loadset

        #: The label identifying aerodynamic trim variables defined on an
        #: AESTAT or AESURF entry.
        self.labels = labels

        #: The magnitude of the aerodynamic extra point degree-of-freedom.
        #: (Real)
        self.uxs = uxs

    def validate(self):
        assert self.q > 0.0, 'q=%s' % self.q
        if len(set(self.labels)) != len(self.labels):
            msg = 'not all labels are unique; labels=%s' % str(self.labels)
            raise RuntimeError(msg)
        if len(self.labels) != len(self.uxs):
            msg = 'nlabels=%s != nux=%s; labels=%s uxs=%s' % (
                len(self.labels), len(self.uxs), str(self.labels), str(self.uxs))
            raise RuntimeError(msg)

    def cross_reference(self, model):
        pass
        #self.suport = model.suport
        #self.suport1 = model.suport1
        #self.aestats = model.aestats
        #self.aelinks = model.aelinks
        #self.aesurf = model.aesurf

    def safe_cross_reference(self, model):
        pass

    def uncross_reference(self):
        pass

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a TRIM card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        sid = integer(card, 1, 'sid')
        mkaeroz = integer(card, 2, 'mkaeroz')
        q = double(card, 3, 'q')
        # 5
        # 6
        cg = [
            double(card, 7, 'cg-x'),
            double(card, 8, 'cg-y'),
            double(card, 9, 'cg-z'),
        ]

        wtmass = double(card, 9, 'wtmass')
        weight = double(card, 10, 'weight')
        inertia = [
            double(card, 11, 'Ixx'),
            double(card, 12, 'Ixy'),
            double(card, 13, 'Iyy'),
            double(card, 14, 'Ixz'),
            double(card, 15, 'Iyz'),
            double(card, 16, 'Izz'),
        ]
        #  TRUE/G  NX      NY      NZ      PDOT    QDOT    RDOT    LOADSET
        true_g = string(card, 17, 'TRUE/G')
        nx = double_or_string(card, 18, 'NX')
        ny = double_or_string(card, 19, 'NY')
        nz = double_or_string(card, 20, 'NZ')
        nxyz = [nx, ny, nz]

        p = double_or_string(card, 21, 'P')
        q = double_or_string(card, 22, 'Q')
        r = double_or_string(card, 23, 'R')
        pqr = [p, q, r]
        loadset = blank(card, 24, 'loadset')

        labels = []
        uxs = []

        i = 25
        n = 1
        while i < len(card):
            label = integer(card, i, 'label%i' % n)
            ux = double_or_string(card, i + 1, 'ux%i' % n)
            if isinstance(ux, string_types):
                assert ux == 'FREE', 'ux=%r' % ux
            print('  label=%s ux=%s' % (label, ux))
            labels.append(label)
            uxs.append(ux)
            i += 2
            n += 1
        return TRIM(sid, mkaeroz, q, cg, true_g, nxyz, pqr, loadset,
                    labels, uxs, comment=comment)

    def raw_fields(self):
        """
        Gets the fields in their unmodified form

        Returns
        -------
        fields : list[varies]
            the fields that define the card

        """
        mach = 1.0
        aeqr = 0.0
        list_fields = ['TRIM', self.sid, mach, self.q]
        nlabels = len(self.labels)
        assert nlabels > 0, self.labels
        for (i, label, ux) in zip(count(), self.labels, self.uxs):
            list_fields += [label, ux]
            if i == 1:
                list_fields += [aeqr]
        if nlabels == 1:
            list_fields += [None, None, aeqr]
        return list_fields

    def repr_fields(self):
        # fixes a Nastran bug
        #aeqr = set_blank_if_default(self.aeqr, 1.0)
        aeqr = 0.

        list_fields = self.raw_fields()
        list_fields[8] = aeqr
        return list_fields

    def write_card(self, size=8, is_double=False):
        return ''
        #card = self.repr_fields()
        #return self.comment + print_card_8(card)

class SPLINE1(Spline):
    """
    Defines an infinite plate spline method for displacements and loads
    transferal between CAERO7 macroelement and structural grid points.

    +---------+------+-------+-------+------+------+----+-----+-------+
    |    1    |  2   |   3   |   4   |  5   |   6  |  7 |  8  |   9   |
    +=========+======+=======+=======+======+======+====+=====+=======+
    | SPLINE1 | EID  | MODEL |  CP   | SETK | SETG | DZ | EPS |       |
    +---------+------+-------+-------+------+------+----+-----+-------+
    | SPLINE1 | 100  |       |       |  1   |  10  | 0. |     |       |
    +---------+------+-------+-------+------+------+----+-----+-------+
    """
    type = 'SPLINE1_ZONA'

    def __init__(self, eid, setk, setg, model=None, cp=None, dz=None, eps=0.01, comment=''):
        """
        Creates a SPLINE1 card, which is useful for control surface
        constraints.

        Parameters
        ----------
        eid : int
            spline id
        comment : str; default=''
            a comment for the card

        """
        # https://www.zonatech.com/Documentation/ZAERO_9.2_Users_3rd_Ed.pdf
        Spline.__init__(self)
        if comment:
            self.comment = comment

        self.eid = eid
        self.model = model
        self.cp = cp
        self.setk = setk
        self.setg = setg
        self.dz = dz
        self.eps = eps
        self.setk_ref = None
        self.setg_ref = None
        self.aero_element_ids = []

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
        model = blank(card, 2, 'model')
        cp = blank(card, 3, 'cp')

        setk = integer(card, 4, 'setk')
        setg = integer(card, 5, 'setg')
        dz = blank(card, 6, 'dz')
        eps = double_or_blank(card, 6, 'eps', 0.01)
        return SPLINE1(eid, setk, setg, model=model, cp=cp, dz=dz, eps=eps,
                       comment=comment)

    def cross_reference(self, model):
        msg = ', which is required by SPLINE1 eid=%s' % self.eid
        self.setg_ref = model.Set(self.setg, msg=msg)
        self.setg_ref.cross_reference_set(model, 'Node', msg=msg)

        self.setk_ref = model.zona.panlsts[self.setk]
        self.setk_ref.cross_reference(model)
        self.aero_element_ids = self.setk_ref.aero_element_ids

    def safe_cross_reference(self, model, xref_errors):
        msg = ', which is required by SPLINE1 eid=%s' % self.eid
        try:
            self.setg_ref = model.Set(self.setg, msg=msg)
            self.setg_ref.safe_cross_reference(model, 'Node', msg=msg)
        except KeyError:
            model.log.warning('failed to find SETx set_id=%s%s; allowed_sets=%s' % (
                self.setg, msg, np.unique(list(model.sets.keys()))))

        try:
            self.setk_ref = model.zona.panlsts[self.setk]
            self.setk_ref.safe_cross_reference(model, xref_errors)
            self.aero_element_ids = self.setk_ref.aero_element_ids
        except KeyError:
            pass

    def uncross_reference(self):
        self.setk_ref = None
        self.setg_ref = None

    def raw_fields(self):
        list_fields = ['SPLINE1', self.eid, self.model, self.cp, self.setk, self.setg,
                       self.dz, self.eps]
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        return self.comment + print_card_8(card)

class SPLINE2(Spline):
    """
    Defines an infinite plate spline method for displacements and loads
    transferal between CAERO7 macroelement and structural grid points.

    +---------+------+-------+------+------+----+-----+-------+-------+
    |    1    |  2   |   3   |  5   |   6  |  6 |  7  |   8   |   9   |
    +=========+======+=======+======+======+====+=====+=======+=======+
    | SPLINE2 | EID  | MODEL | SETK | SETG | DZ | EPS |  CP   | CURV  |
    +---------+------+-------+------+------+----+-----+-------+-------+
    | SPLINE2 | 100  |       |  1   |  10  | 0. |     |       |       |
    +---------+------+-------+------+------+----+-----+-------+-------+
    """
    type = 'SPLINE2_ZONA'

    def __init__(self, eid, setk, setg, model=None, dz=None, eps=0.01, cp=None, curvature=None, comment=''):
        """
        Creates a SPLINE1 card, which is useful for control surface
        constraints.

        Parameters
        ----------
        eid : int
            spline id
        comment : str; default=''
            a comment for the card

        """
        # https://www.zonatech.com/Documentation/ZAERO_9.2_Users_3rd_Ed.pdf
        Spline.__init__(self)
        if comment:
            self.comment = comment

        self.eid = eid
        self.model = model
        self.cp = cp
        self.setk = setk
        self.setg = setg
        self.dz = dz
        self.eps = eps
        self.curvature = curvature
        self.setk_ref = None
        self.setg_ref = None

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
        model = string_or_blank(card, 2, 'model')
        setk = integer(card, 3, 'setk')
        setg = integer(card, 4, 'setg')
        dz = blank(card, 5, 'dz')
        eps = double_or_blank(card, 6, 'eps', 0.01)
        cp = integer_or_blank(card, 7, 'cp', 0)
        curvature = double_or_blank(card, 8, 'curvature', 1.0)
        return SPLINE2(eid, setk, setg, model=model, cp=cp, dz=dz, eps=eps,
                       curvature=curvature, comment=comment)

    def cross_reference(self, model):
        msg = ', which is required by SPLINE1 eid=%s' % self.eid
        self.setg_ref = model.Set(self.setg, msg=msg)
        self.setg_ref.cross_reference_set(model, 'Node', msg=msg)
        #self.nodes_ref = model.Nodes(self.nodes, msg=msg)
        #self.caero_ref = model.CAero(self.caero, msg=msg)
        self.setk_ref = model.zona.panlsts[self.setk]
        self.setk_ref.cross_reference(model)

    def safe_cross_reference(self, model, xref_errors):
        try:
            msg = ', which is required by SPLINE1 eid=%s' % self.eid
            self.setg_ref = model.Set(self.setg, msg=msg)
            self.setg_ref.cross_reference_set(model, 'Node', msg=msg)
        except:
            pass
        #self.nodes_ref = model.Nodes(self.nodes, msg=msg)
        #self.caero_ref = model.CAero(self.caero, msg=msg)
        self.setk_ref = model.zona.panlsts[self.setk]
        self.setk_ref.safe_cross_reference(model, xref_errors)

    def uncross_reference(self):
        self.setk_ref = None
        self.setg_ref = None

    def raw_fields(self):
        list_fields = ['SPLINE2', self.eid, self.model, self.setk, self.setg,
                       self.dz, self.eps, self.cp, self.curvature]
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        return self.comment + print_card_8(card)

class SPLINE3(Spline):
    """
    Defines a 3-D spline for the BODY7 and CAERO7 macroelement.

    +---------+------+-------+-------+------+------+----+-----+-------+
    |    1    |  2   |   3   |   4   |  5   |   6  |  7 |  8  |   9   |
    +=========+======+=======+=======+======+======+====+=====+=======+
    | SPLINE3 | EID  | MODEL |  CP   | SETK | SETG | DZ | EPS |       |
    +---------+------+-------+-------+------+------+----+-----+-------+
    | SPLINE3 | 100  |       |  N/A  |  1   |  10  | 0. |     |       |
    +---------+------+-------+-------+------+------+----+-----+-------+
    """
    type = 'SPLINE3_ZONA'

    def __init__(self, eid, setk, setg, model=None, cp=None, dz=None, eps=0.01, comment=''):
        """
        Creates a SPLINE3 card, which is useful for control surface
        constraints.

        Parameters
        ----------
        eid : int
            spline id
        comment : str; default=''
            a comment for the card

        """
        # https://www.zonatech.com/Documentation/ZAERO_9.2_Users_3rd_Ed.pdf
        Spline.__init__(self)
        if comment:
            self.comment = comment

        self.eid = eid
        self.model = model
        self.cp = cp
        self.setk = setk
        self.setg = setg
        self.dz = dz
        self.eps = eps
        self.setk_ref = None
        self.setg_ref = None
        self.aero_element_ids = []

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
        model = blank(card, 2, 'model')
        cp = blank(card, 3, 'cp')

        setk = integer(card, 4, 'setk')
        setg = integer(card, 5, 'setg')
        dz = blank(card, 6, 'dz')
        eps = double_or_blank(card, 6, 'eps', 0.01)
        return SPLINE3(eid, setk, setg, model=model, cp=cp, dz=dz, eps=eps,
                       comment=comment)

    def cross_reference(self, model):
        msg = ', which is required by SPLINE3 eid=%s' % self.eid
        self.setg_ref = model.Set(self.setg, msg=msg)
        self.setg_ref.cross_reference_set(model, 'Node', msg=msg)
        #self.nodes_ref = model.Nodes(self.nodes, msg=msg)
        #self.caero_ref = model.CAero(self.caero, msg=msg)
        self.setk_ref = model.zona.panlsts[self.setk]
        self.setk_ref.cross_reference(model)
        self.aero_element_ids = self.setk_ref.aero_element_ids

    def safe_cross_reference(self, model, xref_errors):
        msg = ', which is required by SPLINE3 eid=%s' % self.eid
        try:
            self.setg_ref = model.Set(self.setg, msg=msg)
            self.setg_ref.cross_reference_set(model, 'Node', msg=msg)
        except:
            pass
        self.setk_ref = model.zona.panlsts[self.setk]
        self.setk_ref.safe_cross_reference(model, xref_errors)
        self.aero_element_ids = self.setk_ref.aero_element_ids

    def uncross_reference(self):
        self.setk_ref = None
        self.setg_ref = None

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
        list_fields = ['SPLINE3', self.eid, self.model, self.cp, self.setk, self.setg,
                       self.dz, self.eps]
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        return self.comment + print_card_8(card)
