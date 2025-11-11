"""
Defines the following contact cards:

nx contact:
 - BCONP
 - BLSEG
 - BCPARA
 - BCRPARA
 - BCTPARA
 - BCTADD
 - BCTSET
 - BSURF
 - BSURFS
 - BFRIC

msc contact:
 - BCAUTOP (todo)
 - BCBDPRP (todo)
 - BCBMRAD (todo)
 - BCBODY
 - BCBODY1 (todo)
 - BCBZIER (todo)
 - BCGRID  (todo)
 - BCHANGE (todo)
 - BCMATL  (todo)
 - BCMOVE  (todo)
 - BCNURB2 (todo)
 - BCONECT (todo)
 - BCONP
 - BCONPRG (todo)
 - BCONPRP (todo)

 - BCONUDS (todo)
 - BCPARA
 - BCPROP (todo)
 - BCRIGID (todo)
 - BCRGSRF (todo)
 - BCSCAP (todo)
 - BCSEG (todo)
 - BCTABLE (todo)
 - BCTABL1 (todo)
 - BCTRIM (todo)
 - BFRlC (todo)
 - BOUTPUT (todo)
 - BSURF
 - BWIDTH (todo)
 - DYPARAM,CONTACT (todo)

glue:
 - BGADD
 - BGSET

"""
from __future__ import annotations
import warnings
from typing import Optional, Any, TYPE_CHECKING

from pyNastran.utils.numpy_utils import integer_types
from pyNastran.bdf.cards.base_card import BaseCard, expand_thru_by, _node_ids
from pyNastran.bdf.bdf_interface.internal_get import coord_id
from pyNastran.bdf.bdf_interface.assign_type import (
    integer, integer_or_blank, integer_string_or_blank, double_or_blank,
    integer_double_or_blank, string, string_or_blank, string_choice_or_blank,
    double, blank)
from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.field_writer_16 import print_card_16
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf import BDF
    from pyNastran.bdf.bdf_interface.bdf_card import BDFCard


def get_table3d(model: BDF, fric: int) -> TABLE3D:
    return model.table3d[fric]


class BFRIC(BaseCard):
    """
    Slideline Contact Friction
    Defines frictional properties between two bodies in contact.

    +-------+------+-------+-------+-------+
    |   1   |   2  |   3   |   4   |   5   |
    +=======+======+=======+=======+=======+
    | BFRIC | FID  | FSTIF |       |  MU1  |
    +-------+------+-------+-------+-------+

    """
    type = 'BFRIC'

    @classmethod
    def _init_from_empty(cls):
        friction_id = 1
        mu1 = 0.2
        return BFRIC(friction_id, mu1)

    def __init__(self, friction_id: int, mu1: float, fstiff=None, comment=''):
        """
        Creates a BFRIC card, which defines a frictional contact.

        Parameters
        ----------
        friction_id : int
            Friction identification number.
        mu1 : float
            Coefficient of static friction.
        fstiff : float; default=None
            Frictional stiffness in stick. See Remarks 2 and 3
            Default=automatically selected by the program.

        """
        BaseCard.__init__(self)
        if comment:
            self.comment = comment
        self.friction_id = friction_id
        self.fstiff = fstiff
        self.mu1 = mu1

    @classmethod
    def add_card(cls, card: BDFCard, comment: str=''):
        friction_id = integer(card, 1, 'friction_id')
        fstiff = double_or_blank(card, 2, 'fstiff')
        #
        mu1 = double(card, 4, 'mu1')
        return BFRIC(friction_id, mu1, fstiff=fstiff, comment='')

    def raw_fields(self):
        list_fields = [
            'BFRIC', self.friction_id, self.fstiff, None, self.mu1]
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        return self.comment + print_card_8(card)

class BLSEG(BaseCard):
    """
    3D Contact Region Definition by Shell Elements (SOLs 101, 601 and 701)

    Defines a 3D contact region by shell element IDs.

    +=======+====+====+======+====+====+=====+====+====+
    |   1   |  2 |  3 |   4  |  5 |  6 |  7  |  8 |  9 |
    +-------+----+----+------+----+----+-----+----+----+
    | BLSEG | ID | G1 |  G2  | G3 | G4 | G5  | G6 | G7 |
    +-------+----+----+------+----+----+-----+----+----+
    | BLSEG | ID | G1 | THRU | G2 | BY | INC |    |    |
    +-------+----+----+------+----+----+-----+----+----+

    """
    type = 'BLSEG'

    @classmethod
    def _init_from_empty(cls):
        line_id = 1
        nodes = [1]
        return BLSEG(line_id, nodes, comment='')

    def __init__(self, line_id, nodes, comment=''):
        if comment:
            self.comment = comment

        self.line_id = line_id
        self.nodes = expand_thru_by(nodes)
        self.nodes_ref = None

    @classmethod
    def add_card(cls, card: BDFCard, comment: str=''):
        """
        Adds a BLSEG card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        line_id = integer(card, 1, 'line_id')
        #: Number (float)
        nfields = card.nfields
        i = 2
        nodes = []
        while i < nfields:
            d = integer_string_or_blank(card, i, 'field_%s' % i)
            if d is not None:
                nodes.append(d)
            i += 1
        return BLSEG(line_id, nodes, comment=comment)

    def cross_reference(self, model: BDF) -> None:
        msg = f', which is required by BLSEG line_id={self.line_id}'
        self.nodes_ref = model.Nodes(self.nodes, msg=msg)

    def uncross_reference(self) -> None:
        #msg = f', which is required by BLSEG line_id={self.line_id}'
        self.nodes = self.node_ids
        self.nodes_ref = None

    @property
    def node_ids(self):
        """returns nodeIDs for repr functions"""
        return _node_ids(self, nodes=self.nodes_ref, allow_empty_nodes=False, msg='')

    def raw_fields(self):
        list_fields = ['BLSEG', self.line_id] + self.node_ids
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class BCBODY(BaseCard):
    """
    +--------+---------+--------+--------+----------+---------+---------+---------+---------+
    |    1   |    2    |    3   |   4    |    5     |     6   |    7    |    8    |    9    |
    +========+=========+========+========+==========+=========+=========+=========+=========+
    | BCBODY | BID     | DIM    | BEHAV  |   BSID   |   ISTYP | FRIC    |  IDSPL  | CONTROL |
    +--------+---------+--------+--------+----------+---------+---------+---------+---------+
    |        | NLOAD   | ANGVEL | DCOS1  |   DCOS2  |   DCOS3 | VELRB1  |  VELRB2 | VELRB3  |
    +--------+---------+--------+--------+----------+---------+---------+---------+---------+
    |        | ADVANCE | SANGLE | COPTB  |   USER   |         |         |         |         |
    +--------+---------+--------+--------+----------+---------+---------+---------+---------+
    |        | CTYPE   | ISMALL | ITYPE  |   IAUG   |  PENALT | AUGDIST |         |         |
    +--------+---------+--------+--------+----------+---------+---------+---------+---------+
    |        | RIGID   | CGID   | NENT   | --- Rigid Body Name ---      |         |         |
    +--------+---------+--------+--------+----------+---------+---------+---------+---------+
    |        | APPROV  | A      |  N1    | N2       |    N3   |    V1   |    V2   |    V3   |
    +--------+---------+--------+--------+----------+---------+---------+---------+---------+
    |        | RTEMP   | G(temp)|  Tempr | T(Tempr) |         |         |         |         |
    +--------+---------+--------+--------+----------+---------+---------+---------+---------+
    |        | SINK    | G(sink)|  Tsink | T(Tsink) |         |         |         |         |
    +--------+---------+--------+--------+----------+---------+---------+---------+---------+
    |        | GROW    | GF1    |  GF2   |   GF3    | TAB-GF1 | TAB-GF2 | TAB-GF3 |         |
    +--------+---------+--------+--------+----------+---------+---------+---------+---------+
    |        | HEAT    | CFILM  |  TSINK |   CHEAT  | TBODY   | HCV     | HNC     |  ITYPE  |
    +--------+---------+--------+--------+----------+---------+---------+---------+---------+
    |        | BNC     | EMISS  |  HBL   |          |         |         |         |         |
    +--------+---------+--------+--------+----------+---------+---------+---------+---------+
    """
    type = 'BCBODY'
    @classmethod
    def _init_from_empty(cls):
        contact_id = 1
        bsid = 1
        dim = '3D'
        behav = 'DEFORM'
        istype = 4
        fric = 5
        idispl = 0
        word_dict = {}
        return BCBODY(contact_id, bsid, word_dict,
                      dim=dim, behav=behav, istype=istype, fric=fric, idispl=idispl, comment='')

    def __init__(self, contact_id: int, bsid: int,
                 word_dict: dict[str, Any],
                 dim: str='3D', behav: str='DEFORM',
                 istype: int=0, fric: int | float=0, idispl: int=0,
                 heat=None, grow=None,
                 comment: str=''):
        if comment:
            self.comment = comment

        self.contact_id = contact_id
        self.dim = dim
        self.behav = behav
        self.bsid = bsid
        self.istype = istype
        self.fric = fric
        self.idispl = idispl
        self.word_dict = word_dict
        assert behav in ['DEFORM', 'RIGID', 'SYMM', 'ACOUS', 'WORK', 'HEAT'], behav
        for key, values in word_dict.items():
            if key == 'GROW':
                assert len(values) == 6, f'ngrow={len(values)}; grow={values}'
            elif key == 'HEAT':
                # heat = [cfilm, tsink, cheat, tbody, hcv, hnc, itype,
                #         bnc, emiss, hbl, hnl, bnl, hnle,
                #         hnce, bnce, cmb, cms]
                assert len(values) == 17, f'nheat={len(values)}; heat={values}'
            else:
                raise NotImplementedError((key, values))

    @classmethod
    def add_card(cls, card: BDFCard, comment: str=''):
        """
        Adds a BCBODY card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        BID : int (4,1)
           Contact body identification number referenced by
           BCTABLE, BCHANGE, or BCMOVE. (Integer > 0; Required)
        dim : str; default='3D'
           Dimension of body.
           DIM=2D planar body in x-y plane of the basic coordinate system,
                  composed of 2D elements or curves.
           DIM=3D any 3D body composed of rigid surfaces, shell elements
                  or solid elements.
        behav (4,8); default=DEFORM
           Behavior of curve or surface (Character)
            - DEFORM body is deformable
            - RIGID body is rigid
            - SYMM body is a symmetry body
            - ACOUS indicates an acoustic body
            - WORK indicates body is a workpiece
            - HEAT indicates body is a heat-rigid body.
           See Remark 3. for Rigid Bodies..
        BSID : int
            Identification number of a BSURF, BCBOX, BCPROP
            or BCMATL entry if BEHAV=DEFORM. (Integer > 0)
        ISTYP : int (4,3); default=0
           Check of contact conditions. (Integer > 0)
        ISTYP : int
           is not supported in segment-to-segment contact.
           For a deformable body:
           =0 symmetric penetration, double sided contact.
           =1 unsymmetric penetration, single sided contact. (Integer > 0)
           =2 double-sided contact with automatic optimization of contact constraint
              equations (this option is known as “optimized contact”).
              Notes: single-sided contact (ISTYP=1) with the contact bodies arranged properly
              using the contact table frequently runs much faster than ISTYP=2.
              For a rigid body:
           =0 no symmetry condition on rigid body.
           =1 rigid body is a symmetry plane.
        FRIC : int/float (6,7); default=0
            Friction coefficient. (Real > 0 or integer)
            If the value is an integer it represents the ID of a TABL3Di.
        IDSPL : int (4,5); default=0
            Set IDSPL=1 to activate the SPLINE (analytical contact)
            option for a deformable body and for a rigid contact surface.
            Set it to zero or leave blank to not have analytical contact. (Integer)
        NLOAD : int or None
            Enter a positive number if "load controlled" and rotations are allowed (Integer). The
            positive number is the grid number where the moments or rotations are applied. The
            rotations are specified using SPCD at grid ID NLOAD and can be specified using dof's
            1-3 (for rotation about x, y, z respectively), or by dof's 4-6 (for rotation about x, y, z
            respectively).
            Note: This rotation takes the position of the grid point defined in CGID field as the
            center of rotation.
        ANGVEL : int/float; default=0.0
            Angular velocity or angular position about local axis through center of rotation. If the
            value is an integer it represents the ID of a TABLED1, TABLED2 or TABL3D, i.e., a
            time-dependent or multi-dimensional table; however, no log scales, only linear scales.
            (Real or Integer; Default = 0.0)
        DCOSi : int/float; default=0.0
            Components of direction cosine of local axis if ANGVEL is nonzero. If the value is an
            integer, it represents the ID of a TABLED1, TABLED2 or TABL3D, i.e., a time-dependent
            or multi-dimensional table; however, no log scales, only linear scales. (Real
            or Integer; Default=0.0) In 2D contact only DCOS3 is used and the Default is 1.0.
        VELRBi : int/float; default=0.0
            Translation velocity or final position (depending on the value of CONTROL) of rigid
            body at the grid point defined in CGID filed. For velocity control only, if the value is
            an integer, it represents the ID of TABLED1, TABLED2 or TABL3D, i.e., a time-dependent
            or multi-dimensional table; however, no log scales, only linear scales. Only
            VELRB1 and VELRB2 are used in 2D contact. (Real or Integer; Default = 0.0)

        """
        contact_id = integer(card, 1, 'contact_id')
        dim = string_choice_or_blank(card, 2, 'dim',
                                     ('2D', '3D'),
                                     default='3D')

        behav = string_choice_or_blank(card, 3, 'behav',
                                       ('RIGID', 'DEFORM', 'SYMM', 'ACOUS', 'WORK', 'HEAT'),
                                       default='DEFORM')
        if behav == 'DEFORM':
            bsid = integer(card, 4, 'bsid')
        else:
            bsid = integer_double_or_blank(card, 4, 'bsid')

        istype = integer_or_blank(card, 5, 'istype', default=0)
        fric = integer_double_or_blank(card, 6, 'fric', default=0)
        idispl = integer_or_blank(card, 7, 'idispl', default=0)
        control = integer_or_blank(card, 8, 'control', default=0)

        # NLOAD   | ANGVEL | DCOS1  | DCOS2|  DCOS3 | VELRB1  | VELRB2 | VELRB3
        word_nload = integer_string_or_blank(card, 9, 'nload (int) / word (str)', default=None)
        i = 9
        if word_nload is None or isinstance(word_nload, int):
            nload = integer_or_blank(card, 9, 'nload')
            ang_vel = double_or_blank(card, 10, 'ang_vel', default=0.0)
            dcos = [
                integer_double_or_blank(card, 11, 'dcos1', default=0.0),
                integer_double_or_blank(card, 12, 'dcos2', default=0.0),
                integer_double_or_blank(card, 13, 'dcos3', default=0.0),
            ]
            vel_rb = [
                integer_double_or_blank(card, 14, 'vel_rb1', default=0.0),
                integer_double_or_blank(card, 15, 'vel_rb2', default=0.0),
                integer_double_or_blank(card, 16, 'vel_rb3', default=0.0),
            ]
            i += 8

        heat = []
        grow = []
        word_dict = {}
        old_word = None
        while i < len(card):
            word = string_or_blank(card, i, 'word (str)', default=None)
            if word is None:
                raise RuntimeError(f'should be broken by {old_word}')
            assert word in ['ADVANCE', 'HEAT', 'GROW', 'NURBS', 'PATCH3D',
                            'RIGID', 'BEZIER'], word
            #print('*', word)
            if word == 'ADVANCE':
                # | ADVANCE | SANGLE | COPTB | | MIDNO |
                sangle = double_or_blank(card, i+1, 'sangle', default=60.)
                coptb = integer(card, i+2, 'coptb')
                blank(card, i+3, 'user')
                blank(card, i+4, 'min_nod')
                # “ADVANCE”
                #     The entries for this continuation line are for advanced options starting with
                #     MD Nastran R2.
                # SANGLE
                #     Threshold for automatic discontinuity detection in degrees. (Real; Default = 60.0)
                #     Used for SPLINE option in SOL 400 only. SANGLE is not used when IDSPL ≥ 0.
                # COPTB
                #     Flag to indicate how body surfaces may contact. See Remark 9. on the BCTABLE entry.
                #     (Integer; Default = 0)
                # MIDNOD
                #     Mid-side node projection flag. (Integer > 0; Default = 0)
                #     When MIDNOD > 0 and IDSPL 0, the mid-side grid of quadratic elements are
                #     projected onto the selected spline surfaces. This operation is performed before the
                #     contact process starts and it may change the location of grids in contact bodies. It may
                #     operate in combination with the initial stress-free contact.
                i += 8
            elif word == 'HEAT':
                # “HEAT”
                #    The entries of this continuation line(s) are for contact in heat transfer in a pure thermal
                #    analysis or in a coupled thermal/structural analysis. In a pure structural analysis they are
                #    ignored.
                # CFILM (9,1)/(10,1)
                #     Heat transfer coefficient (film) to environment. (Real or Integer, Default = 0.0) If Real,
                #     the value entered is the film coefficient. If Integer, the value entered is the ID of a
                #     TABLEM1 or TABLEM2 entry specifying the heat transfer coefficient vs temperature
                #     or a TABL3D entry specifying the film coefficient vs temperature and possibly other
                #     variables.
                # TSINK (9,2)/(10,2)
                #     Environment sink temperature. (Real or Integer, Default = 0.0). If Real, the value
                #     entered is the sink temperature. If Integer, the value entered is the ID of a TABLED1
                #     or TABLED2 entry specifying temperature vs time or a TABL3D entry specifying the
                #     sink temperature vs time and possibly other variables. When entered as a negative
                #     integer its absolute value is a scalar point identification number. If a scalar point is
                #     specified on this entry it need not be defined on an SPOINT entry.
                # CHEAT (9,3)/(10,3)
                #     Contact heat transfer coefficient. (Real or Integer; Default = 0.0). If Real, the value
                #     entered is the contact heat transfer coefficient. If Integer, the value entered is the ID of
                #     a TABLEM1 or TABLEM2 entry specifying the contact heat transfer coefficient vs
                #     temperature or a TABL3D entry specifying the contact heat transfer coefficient vs
                #     temperature and possibly other variables.
                # TBODY (9,4)/(10,4)
                #     Body temperature. (Real or Integer; Default = 0.0). If Real, the value entered is the body
                #     temperature. If Integer, the value entered is the ID of a TABLED1 or TABLED2 entry
                #     specifying the body temperature vs time or a TABL3D entry specifying the body
                #     temperature vs time and possibly other variables. When entered as a negative integer its
                #     absolute value is a scalar point identification number. If a scalar point is specified on
                #     this entry it need not be defined on an SPOINT entry.
                # HCV (9,5)/(10,5)
                #     Convection coefficient for near field behavior (Real or Integer; Default = 0.0). If Real
                #     the value entered is the near field convection coefficient. If Integer, the value entered is
                #     the ID of a TABLEM1 or TABLEM2 entry specifying the near field convection
                #     coefficient vs temperature or a TABL3D entry specifying the
                # HEAT CFILM TSINK CHEAT TBODY HCV HNC  ITYPE
                #      BNC   EMISS HBL   HNL   BNL HNLE BNLE
                #      HNCE  BNCE  CMB   CMS
                assert len(heat) == 0, heat
                cfilm = double(card, i+1, 'cfilm')
                tsink = double(card, i+2, 'tsink')
                cheat = double(card, i+3, 'cheat')
                tbody = double(card, i+4, 'tbody')
                hcv = double(card, i+5, 'hcv')
                hnc = double(card, i+6, 'hnc')
                itype = integer_or_blank(card, i+7, 'itype')
                i += 8

                bnc = double_or_blank(card, i+1, 'bnc', default=1.)
                emiss = double_or_blank(card, i+2, 'emiss', default=0.)
                hbl = double_or_blank(card, i+3, 'hbl', default=0.)
                hnl = integer_double_or_blank(card, i+4, 'hnl', default=0.)
                bnl = integer_double_or_blank(card, i+5, 'bnl', default=1.)
                hnle = integer_double_or_blank(card, i+6, 'hnle', default=0.)
                bnle = integer_double_or_blank(card, i+7, 'bnle', default=1.)
                i += 8

                hnce = integer_double_or_blank(card, i+1, 'hnce', default=0.)
                bnce = integer_double_or_blank(card, i+2, 'bnce', default=1.)
                cmb = double_or_blank(card, i+3, 'cmb', default=0.)
                cms = double_or_blank(card, i+4, 'cms', default=0.)
                heat = [cfilm, tsink, cheat, tbody, hcv, hnc, itype,
                        bnc, emiss, hbl, hnl, bnl, hnle,
                        hnce, bnce, cmb, cms]
                word_dict[word] = heat
                i += 8
            elif word == 'GROW':
                assert grow is None, grow
                #'GROW' GF1 GF2 GF3 TAB-GF1 TAB-GF2 TAB-GF3
                gf1 = double_or_blank(card, i+1, 'GF1', default=1.0)
                gf2 = double_or_blank(card, i+2, 'GF2', default=1.0)
                gf3 = double_or_blank(card, i+3, 'GF3', default=1.0)
                tab_gf1 = integer_or_blank(card, i+4, 'tab_GF1')
                tab_gf2 = integer_or_blank(card, i+5, 'tab_GF2')
                tab_gf3 = integer_or_blank(card, i+6, 'tab_GF3')
                #blank = blank(card, i+7, 'GROW blank')
                #print('grow values =', [gf1, gf2, gf3, tab_gf1, tab_gf2, tab_gf3])
                grow = [gf1, gf2, gf3, tab_gf1, tab_gf2, tab_gf3]
                word_dict[word] = grow
                i += 8
                #GF1 GF2 GF3 TAB-GF1 TAB-GF2 TAB-GF3
            elif word == 'NURBS':
                i, values = _get_bcbody_section_values(card, i, word)
                #print('end of NURBS -> ', valuei)
            elif word == 'PATCH3D':
                #'PATCH3D' NPATCH
                #          IDP G1 G2 G3 G4
                #          IDP G1 G2 G3 G4
                i, values = _get_bcbody_section_values(card, i, word)
            elif word == 'RIGID':
                i += 8
            elif word == 'BEZIER':
                i, values = _get_bcbody_section_values(card, i, word)
                #print(word, values)
            else:
                raise NotImplementedError(word)
            old_word = word

        return BCBODY(contact_id, bsid, word_dict,
                      dim=dim, behav=behav, istype=istype,
                      fric=fric, idispl=idispl,
                      comment=comment)

    def cross_reference(self, model: BDF) -> None:
        if isinstance(self.fric, integer_types) and self.fric > 0:
            fric_ref = get_table3d(model, fric)

    def raw_fields(self):
        #     | BCBODY | BID     | DIM    | BEHAV  |   BSID   |   ISTYP | FRIC    |  IDSPL  | CONTROL |
        #     |   ***  | NLOAD   | ANGVEL | DCOS1  |   DCOS2  |   DCOS3 | VELRB1  |  VELRB2 | VELRB3  |
        #     |        | ADVANCE | SANGLE | COPTB  |   USER   |         |         |         |         |
        #     |        | CTYPE   | ISMALL | ITYPE  |   IAUG   |  PENALT | AUGDIST |         |         |
        #     |   ***  | RIGID   | CGID   | NENT   | --- Rigid Body Name ---      |         |         |
        #     |        | APPROV  | A      |  N1    | N2       |    N3   |    V1   |    V2   |    V3   |
        #     |        | RTEMP   | G(temp)|  Tempr | T(Tempr) |         |         |         |         |
        #     |        | SINK    | G(sink)|  Tsink | T(Tsink) |         |         |         |         |
        #     |   ***  | GROW    | GF1    |  GF2   |   GF3    | TAB-GF1 | TAB-GF2 | TAB-GF3 |         |
        #     |   ***  | HEAT    | CFILM  |  TSINK |   CHEAT  | TBODY   | HCV     | HNC     |  ITYPE  |
        #     |        | BNC     | EMISS  |  HBL   |          |         |         |         |         |
        list_fields = [
            # line 1
            'BCBODY', self.contact_id, self.dim, self.behav, self.bsid, self.istype, self.fric, self.idispl, None,
        ]
        if self.behav == 'HEAT':
            # 17+1 = 18; 24-18 = 6
            list_fields += ['HEAT'] + self.heat + [None] * 6
        elif self.behav == 'GROW':
            # 6+1 = 7; 8-7=1
            list_fields += ['GROW'] + self.grow + [None]
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        return self.comment + print_card_8(card)

def _get_bcbody_section_values(card, i: int, word: str) -> tuple[int, list[Any]]:
    """gets all the values of a multi-line section"""
    i0 = i
    values = []
    i += 8
    valuei = isinstance(string_or_blank(card, i, f'{word}_word'), str)
    #print('i', i, 1)
    while not isinstance(valuei, str) and i < len(card):
        i += 8
        valuei = string_or_blank(card, i, f'{word}_word')
        #print('i', i, valuei)
    values = card[i0:i]
    return i, values


class BCONP(BaseCard):
    """
    3D Contact Region Definition by Shell Elements (SOLs 101, 601 and 701)

    Defines a 3D contact region by shell element IDs.

    +-------+----+-------+--------+-----+------+--------+-------+-----+
    |   1   |  2 |   3   |   4    |  5  |   6  |   7    |   8   |  9  |
    +=======+====+=======+========+=====+======+========+=======+=====+
    | BCONP | ID | SLAVE | MASTER |     | SFAC | FRICID | PTYPE | CID |
    +-------+----+-------+--------+-----+------+--------+-------+-----+
    | BCONP | 95 |   10  |   15   |     |  1.0 |   33   |   1   |     |
    +-------+----+-------+--------+-----+------+--------+-------+-----+

    """
    type = 'BCONP'
    @classmethod
    def _init_from_empty(cls):
        contact_id = 1
        slave = 2
        master = 3
        sfac = 4
        friction_id = 5
        ptype = 'cat'
        cid = 0
        return BCONP(contact_id, slave, master, sfac, friction_id, ptype, cid, comment='')

    def __init__(self, contact_id, slave, master, sfac, friction_id, ptype, cid, comment=''):
        if comment:
            self.comment = comment

        self.contact_id = contact_id
        self.slave = slave
        self.master = master
        self.sfac = sfac
        self.friction_id = friction_id
        self.ptype = ptype
        self.cid = cid
        self.cid_ref = None
        self.friction_id_ref = None
        self.slave_ref = None
        self.master_ref = None

    @classmethod
    def add_card(cls, card: BDFCard, comment: str=''):
        """
        Adds a BCONP card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        contact_id = integer(card, 1, 'contact_id')
        slave = integer(card, 2, 'slave')
        master = integer(card, 3, 'master')
        sfac = double_or_blank(card, 5, 'sfac', default=1.0)
        friction_id = integer_or_blank(card, 6, 'fric_id')
        ptype = integer_or_blank(card, 7, 'ptype', default=1)
        cid = integer_or_blank(card, 8, 'cid', default=0)
        return BCONP(contact_id, slave, master, sfac, friction_id, ptype, cid,
                     comment=comment)

    def cross_reference(self, model: BDF) -> None:
        msg = f', which is required by BCONP line_id={self.contact_id}'
        #self.nodes_ref = model.Nodes(self.nodes, msg=msg)
        self.cid_ref = model.Coord(self.cid, msg=msg)
        if self.friction_id is not None:
            self.friction_id_ref = model.bfric[self.friction_id]
        self.slave_ref = model.blseg[self.slave]
        self.master_ref = model.blseg[self.master]

    def uncross_reference(self) -> None:
        self.cid = self.Cid()
        self.friction_id = self.FrictionId()
        self.cid_ref = None
        self.friction_id_ref = None
        self.slave_ref = None
        self.master_ref = None

    def Cid(self) -> int:
        return coord_id(self.cid_ref, self.cid)

    def FrictionId(self) -> int:
        if self.friction_id_ref is not None:
            return self.friction_id_ref.friction_id
        return self.friction_id

    def Slave(self) -> int:
        if self.slave_ref is not None:
            return self.slave_ref.line_id
        return self.slave

    def Master(self) -> int:
        if self.master_ref is not None:
            return self.master_ref.line_id
        return self.master

    def raw_fields(self):
        list_fields = [
            'BCONP', self.contact_id, self.Slave(), self.Master(), None, self.sfac,
            self.FrictionId(), self.ptype, self.Cid()]
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class BSURF(BaseCard):
    """
    3D Contact Region Definition by Shell Elements (SOLs 101, 601 and 701)

    Defines a 3D contact region by shell element IDs.

    +-------+------+------+-------+-------+--------+------+------+------+
    |   1   |   2  |   3  |   4   |   5   |    6   |  7   |   8  |   9  |
    +=======+======+======+=======+=======+========+======+======+======+
    | BSURF |  ID  | EID1 | EID2  | EID3  |  EID4  | EID5 | EID6 | EID7 |
    +-------+------+------+-------+-------+--------+------+------+------+
    |       | EID8 | EID9 | EID10 |  etc. |        |      |      |      |
    +-------+------+------+-------+-------+--------+------+------+------+
    | BSURF |  ID  | EID1 |  THRU | EID2  |   BY   | INC  |      |      |
    +-------+------+------+-------+-------+--------+------+------+------+
    |       | EID8 | EID9 | EID10 | EID11 |  etc.  |      |      |      |
    +-------+------+------+-------+-------+--------+------+------+------+
    |       | EID8 | THRU | EID9  |  BY   |  INC   |      |      |      |
    +-------+------+------+-------+-------+--------+------+------+------+
    | BSURF |  15  |  5   | THRU  |  21   |   BY   |  4   |      |      |
    +-------+------+------+-------+-------+--------+------+------+------+
    |       |  27  |  30  |  32   |  33   |        |      |      |      |
    +-------+------+------+-------+-------+--------+------+------+------+
    |       |  35  | THRU |  44   |       |        |      |      |      |
    +-------+------+------+-------+-------+--------+------+------+------+
    |       |  67  |  68  |  70   |  85   |   92   |      |      |      |
    +-------+------+------+-------+-------+--------+------+------+------+

    """
    type = 'BSURF'

    @classmethod
    def _init_from_empty(cls):
        sid = 1
        eids = [1]
        return BSURF(sid, eids, comment='')

    def __init__(self, sid, eids, comment=''):
        if comment:
            self.comment = comment
        #: Set identification number. (Unique Integer > 0)
        self.sid = sid
        #: Element identification numbers of shell elements. (Integer > 0)
        self.eids = eids

    @classmethod
    def add_card(cls, card: BDFCard, comment: str=''):
        """
        Adds a BSURF card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        sid = integer(card, 1, 'sid')
        #: Number (float)
        nfields = card.nfields
        i = 2
        eid_data = []
        while i < nfields:
            d = integer_string_or_blank(card, i, 'field_%s' % i)
            if d is not None:
                eid_data.append(d)
            i += 1
        eids = expand_thru_by(eid_data)
        return BSURF(sid, eids, comment=comment)

    def raw_fields(self):
        fields = ['BSURF', self.sid]
        return fields + list(self.eids)

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class BSURFS(BaseCard):
    """
    Defines a 3D contact region by the faces of the CHEXA, CPENTA or CTETRA
    elements.

    Notes
    -----
    1. The continuation field is optional.
    2. BSURFS is a collection of one or more element faces on solid elements.
       BSURFS defines a contact region which may act as a contact source
       (contactor) or target.
    3. The ID must be unique with respect to all other BSURFS, BSURF, and
       BCPROP entries.

    """
    type = 'BSURFS'

    @classmethod
    def _init_from_empty(cls):
        bsurfs_id = 1
        eids = [1]
        g1s = [1]
        g2s = [1]
        g3s = [1]
        return BSURFS(bsurfs_id, eids, g1s, g2s, g3s, comment='')

    def __init__(self, bsurfs_id: int,
                 eids: list[int],
                 g1s: list[int],
                 g2s: list[int],
                 g3s: list[int],
                 comment: str=''):
        if comment:
            self.comment = comment
        #: Identification number of a contact region. See Remarks 2 and 3.
        #: (Integer > 0)
        self.id = bsurfs_id

        #: Element identification numbers of solid elements. (Integer > 0)
        self.eids = eids

        #: Identification numbers of 3 corner grid points on the face (triangular
        #: or quadrilateral) of the solid element. (Integer > 0)
        self.g1s = g1s
        self.g2s = g2s
        self.g3s = g3s

    @classmethod
    def add_card(cls, card: BDFCard, comment: str=''):
        """
        Adds a BSURFS card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        bsurfs_id = integer(card, 1, 'id')
        eids = []
        g1s = []
        g2s = []
        g3s = []

        n = card.nfields - 5
        i = 0
        j = 1
        while i < n:
            eid = integer(card, 5 + i, 'eid%s' % j)
            g1 = integer(card, 5 + i + 1, 'g3_%s' % j)
            g2 = integer(card, 5 + i + 2, 'g2_%s' % j)
            g3 = integer(card, 5 + i + 3, 'g1_%s' % j)
            j += 1
            i += 4
            eids.append(eid)
            g1s.append(g1)
            g2s.append(g2)
            g3s.append(g3)
        return BSURFS(bsurfs_id, eids, g1s, g2s, g3s, comment=comment)

    def raw_fields(self):
        fields = ['BSURFS', self.id, None, None, None]
        for eid, g1, g2, g3 in zip(self.eids, self.g1s, self.g2s, self.g3s):
            fields += [eid, g1, g2, g3]
        return fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class BCTSET(BaseCard):
    """
    3D Contact Set Definition (SOLs 101, 601 and 701 only)
    Defines contact pairs of a 3D contact set.

    +--------+-------+------+-------+-------+-------+-------+
    |   1    |   2   | 3    |  4    |   5   |   6   |   7   |
    +========+=======+======+=======+=======+=======+=======+
    | BCTSET | CSID  | SID1 | TID1  | FRIC1 | MIND1 | MAXD1 |
    +--------+-------+------+-------+-------+-------+-------+
    |        |       | SID2 | TID2  | FRIC2 | MIND2 | MAXD2 |
    +--------+-------+------+-------+-------+-------+-------+
    |        |  etc. |      |       |       |       |       |
    +--------+-------+------+-------+-------+-------+-------+

    """
    type = 'BCTSET'

    @classmethod
    def _init_from_empty(cls):
        csid = 1
        sids = [1]
        tids = [1]
        frictions = [0.01]
        min_distances = [0.1]
        max_distances = [1.]
        return BCTSET(csid, sids, tids, frictions, min_distances, max_distances, comment='', sol=101)

    def __init__(self, csid, sids, tids, frictions, min_distances, max_distances,
                 comment='', sol=101):
        if comment:
            self.comment = comment
        #: CSID Contact set identification number. (Integer > 0)
        self.csid = csid
        #: SIDi Source region (contactor) identification number for contact pair i.
        #: (Integer > 0)
        self.sids = sids

        #: TIDi Target region identification number for contact pair i. (Integer > 0)
        self.tids = tids

        #: FRICi Static coefficient of friction for contact pair i. (Real; Default=0.0)
        self.frictions = frictions

        #: MINDi Minimum search distance for contact. (Real) (Sol 101 only)
        self.min_distances = min_distances

        #: MAXDi Maximum search distance for contact. (Real) (Sol 101 only)
        self.max_distances = max_distances

    @classmethod
    def add_card(cls, card, comment='', sol=101):
        csid = integer(card, 1, 'csid')
        sids = []
        tids = []
        frictions = []
        min_distances = []
        max_distances = []

        nfields = card.nfields
        i = 2
        j = 1
        while i < nfields:
            sids.append(integer(card, i, 'sid%s' % j))
            tids.append(integer(card, i + 1, 'tid%s' % j))
            frictions.append(double_or_blank(card, i + 2, 'fric%s' % j, default=0.0))

            # we ignore what the NX QRG says about the min/max distance for SOL 401
            # if you don't, Nastran will crash
            #if sol == 101:
            min_distances.append(double_or_blank(card, i + 3, 'mind%s' % j, default=0.0))
            max_distances.append(double_or_blank(card, i + 4, 'maxd%s' % j, default=0.0))
            #else:
                #min_distances.append(None)
                #max_distances.append(None)

            i += 8
            j += 1
        return BCTSET(csid, sids, tids, frictions, min_distances,
                      max_distances, comment=comment,
                      sol=sol)

    def raw_fields(self):
        fields = ['BCTSET', self.csid]
        for sid, tid, fric, mind, maxd in zip(self.sids, self.tids, self.frictions,
                                              self.min_distances, self.max_distances):
            fields += [sid, tid, fric, mind, maxd, None, None, None]
        return fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)


class BCRPARA(BaseCard):
    """
    +---------+------+------+--------+------+-----+---+---+---+----+
    |    1    |   2  |   3  |   4    |   5  |  6  | 7 | 8 | 9 | 10 |
    +=========+======+======+========+======+=====+===+===+===+====+
    | BCRPARA | CRID | SURF | OFFSET | TYPE | GP  |   |   |   |    |
    +---------+------+------+--------+------+-----+---+---+---+----+
    """
    type = 'BCRPARA'

    @classmethod
    def _init_from_empty(cls):
        crid = 1
        return BCRPARA(crid, offset=None, surf='TOP', Type='FLEX', grid_point=0, comment='')

    def __init__(self, crid: int, offset: Optional[float]=None,
                 surf: str='TOP', Type: str='FLEX', grid_point: int=0,
                 comment: str=''):
        """
        Creates a BCRPARA card

        Parameters
        ----------
        crid : int
            CRID Contact region ID.
        offset : float; default=None
            Offset distance for the contact region (Real > 0.0).
            None : OFFSET value in BCTPARA entry
        surf : str; default='TOP'
            SURF Indicates the contact side. See Remark 1.  {'TOP', 'BOT'; )
        Type : str; default='FLEX'
            Indicates whether a contact region is a rigid surface if it
            is used as a target region. {'RIGID', 'FLEX'}.
            This is not supported for SOL 101.
        grid_point : int; default=0
            Control grid point for a target contact region with TYPE=RIGID
            or when the rigid-target algorithm is used.  The grid point
            may be used to control the motion of a rigid surface.
            (Integer > 0).  This is not supported for SOL 101.
        comment : str; default=''
            a comment for the card

        """
        if comment:
            self.comment = comment

        #: CRID Contact region ID. (Integer > 0)
        self.crid = crid

        #: SURF Indicates the contact side. See Remark 1. (Character = "TOP" or
        #: "BOT"; Default = "TOP")
        self.surf = surf

        #: Offset distance for the contact region. See Remark 2. (Real > 0.0,
        #: Default =OFFSET value in BCTPARA entry)
        self.offset = offset

        #: Indicates whether a contact region is a rigid surface if it is used as a
        #: target region. See Remarks 3 and 4. (Character = "RIGID" or "FLEX",
        #: Default = "FLEX"). This is not supported for SOL 101.
        self.Type = Type

        #: Control grid point for a target contact region with TYPE=RIGID or
        #: when the rigid-target algorithm is used. The grid point may be
        #: used to control the motion of a rigid surface. (Integer > 0)
        #: This is not supported for SOL 101.
        self.grid_point = grid_point

    @classmethod
    def add_card(cls, card: BDFCard, comment: str=''):
        """
        Adds a BCRPARA card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        crid = integer(card, 1, 'crid')
        surf = string_or_blank(card, 2, 'surf', default='TOP')
        offset = double_or_blank(card, 3, 'offset', default=None)
        Type = string_or_blank(card, 4, 'type', default='FLEX')
        grid_point = integer_or_blank(card, 5, 'grid_point', default=0)
        return BCRPARA(crid, surf=surf, offset=offset, Type=Type,
                       grid_point=grid_point, comment=comment)

    def raw_fields(self):
        fields = ['BCRPARA', self.crid, self.surf, self.offset, self.Type, self.grid_point]
        return fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)


class BCPARA(BaseCard):
    """
    Defines contact parameters used in SOL 600.

    +--------+---------+--------+--------+--------+--------+---------+--------+
    |   1    |    2    |    3   |   4    |   5    |   6    |    7    |    8   |
    +========+=========+========+========+========+========+=========+========+
    | BCPARA |  CSID   | Param1 | Value1 | Param2 | Value2 | Param3  | Value3 |
    +--------+---------+--------+--------+--------+--------+---------+--------+
    |        | Param4  | Value4 | Param5 | Value5 |  etc.  |         |        |
    +--------+---------+--------+--------+--------+--------+---------+--------+
    | BCPARA | NBODIES |   4    |  BIAS  |   0.5  |        |         |        |
    +--------+---------+--------+--------+--------+--------+---------+--------+

    """
    type = 'BCPARA'

    @classmethod
    def _init_from_empty(cls):
        csid = 1
        params = {'NBODIES' : 4}
        return BCTPARM(csid, params, comment='')

    def _finalize_hdf5(self, encoding):
        """hdf5 helper function"""
        keys, values = self.params
        self.params = {key : value for key, value in zip(keys, values)}

    def __init__(self, csid: int, params: dict[str, Optional[int | float]],
                 comment: str=''):
        """
        Creates a BCPARA card

        Parameters
        ----------
        csid : int
            ID is not used and should be set to zero. Only one BCPARA should be
            entered and it applies to all subcases.
        csid : int
            Contact set ID. Parameters defined in this command apply to
            contact set CSID defined by a BCTSET entry. (Integer > 0)
        params : dict[key] : int/float
            the optional parameters
        comment : str; default=''
            a comment for the card

        """
        if comment:
            self.comment = comment

        #: Contact set ID. Parameters defined in this command apply to
        #: contact set CSID defined by a BCTSET entry. (Integer > 0)
        self.csid = csid
        self.params = params

    @classmethod
    def add_card(cls, card: BDFCard, comment: str=''):
        """
        Adds a BCPARA card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        csid = integer_or_blank(card, 1, 'csid', default=0)
        i = 2
        j = 1
        params = {}
        ivalue_line = 0
        allowed_params = [
            'BIAS', 'DDULMT', 'ERRBAS', 'ERROR', 'FKIND',
            'FNTOL', 'FTYPE', 'IBSEP', 'ISPLIT', 'ICSEP',
            'LINQUAD', 'LINCNT', 'MAXNOD', 'METHOD', 'MAXENT',
            'NBODIES', 'NLGLUE', 'SLDLMT', 'SEGSYM', 'THKOFF',
        ]
        while i < card.nfields:
            param = string_or_blank(card, i, f'param{j}')
            i += 1

            if param is None:
                blank(card, i, f'blank_value{j}')
                i += 1
                j += 1
                continue
            #print('param', param)
            assert param in allowed_params, f'param={param!r} is not supported; allowed={allowed_params}'
            if param == 'BIAS':
                # Contact tolerance bias factor. (0.0 ≤ Real ≤ 1.0.;
                # Default = 0.9 for IGLUE=0, if field left blank or 0.0 (to obtain a near zero value,
                # enter 1.0E-16).
                # Default = 0.0 for IGLUE <>0. Note that when IGLUE<>0, BIAS can only be given
                # by the BCTABLE or BCONPRG.
                # Default = 0.0 for BEHAVE=SYMM on BCBODY, if field left blank or 0.0.
                value = double_or_blank(card, i, f'{param}_value{j}', default=0.0)
            elif param == 'DDULMT':
                # Maximum value of DDU in a Newton-Raphson cycle. (Real ≥ 0.0;
                # Default = 0.0, no limitation)
                value = double_or_blank(card, i, f'{param}_value{j}', default=0.0)
            elif param == 'ERRBAS':
                # Error computation option.
                # Integer
                # 0 = compute error globally or
                # 1 = calculate error based on each pair of slave-master;
                # Default = 0)
                value = integer_or_blank(card, i, f'{param}_value{j}', default=0)
                assert value in [0, 1], f'ERRBAS must be [0, 1]; ERRBAS={value}'
            elif param == 'ERROR':
                # Distance below which a node is considered touching a body. (Real; Default = blank).
                # Automatically calculated if left blank. If left blank, the code calculates ERROR as the
                # smallest value resulting from: Either dividing the smallest nonzero element
                # dimension (plates or solids) in the contact body by 20. Or dividing the thinnest shell
                # thickness in the contact body by 4. This value is then used for all contact pairs.
                value = double_or_blank(card, i, f'{param}_value{j}')
            elif param == 'FKIND':
                # FKIND
                # (2,5)
                # Friction kind. (Integer 0 or 1)
                # 0 Friction based on nodal stress.
                # 1 Default if friction is present and FKIND is not entered.
                # Friction based on nodal force.
                value = integer_or_blank(card, i, f'value{j}', 1)
                assert value in [0, 1], f'FKIND must be [6]; FKIND={value}'
            elif param == 'FNTOL':
                # FNTOL Separation force (or stress if separation is controlled by stress as determined by
                # IBSEP) above which a node separates from a body. Automatically calculated if left
                # blank. (Real; Default = blank)
                value = double_or_blank(card, i, f'{param}_value{j}')

            elif param == 'FTYPE':
                # FTYPE Friction type. See Remark 5. (Integer)
                # 0 No friction. (Default)
                # 6 Bilinear Coulomb friction.
                # 7 Bilinear Shear friction.
                #
                # FTYPE
                # Friction type. (Integer)
                # 0 No friction. (Default)
                # 1 Shear friction.
                # 2 Coulomb Friction.
                # 3 Shear friction for rolling.
                # 4 Coulomb friction for rolling.
                # 5 Stick-slip Coulomb friction.
                # 6 Bilinear Coulomb friction. (Default if friction is present and FTYPE is not entered.)
                # 7 Bilinear Shear friction.
                value = integer_or_blank(card, i, f'{param}_value{j}', default=0)
                assert value in [0, 1, 2, 3, 4, 5, 6], f'FTYPE must be [0, 1, 2, 3, 4, 5, 6]; FTYPE={value}'
            elif param == 'IBSEP':
                # Flag for separation based on stresses or forces. (Integer > 0; Default = 0)
                # 0 Separation based on forces.
                # 1 Separation based on absolute stresses (force/area)
                # 2 Separation based on absolute stress (extrapolating integration point stresses)
                # 3 Relative nodal stress (force/area)
                # 4 Separation based on relative stress (extrapolating integration point stresses)
                # Only option 2 and 4 can be used with mid-side node elements where the mid-side
                # nodes contact (LINQUAD=-1). For segment to segment contact, the program will
                # set IBSEP to 2 internally. See Remarks 6. and 10.
                value = integer_or_blank(card, i, f'{param}_value{j}', default=0)
                assert value in [0, 1, 2, 3, 4], f'IBSEP must be [0, 1, 2, 3, 4]; IBSEP={value}'
            elif param == 'ISPLIT':
                # ISPLIT (2,7) Flag for increment splitting procedure. (Integer > 0; Default = 3 for
                # statics and 0 for dynamics)
                # 0 Uses increment splitting procedures for the fixed time step
                # procedures.
                # BCPARA 1115
                # Contact Parameters in SOL 600
                # Main Index
                # 1 Suppresses splitting for the fixed time step procedures. Any
                # penetration that occurred in the previous increment is
                # adjusted for equilibrium at the start of the next increment.
                # This method may require smaller time steps then the other
                # methods
                # 2 Suppresses splitting for adaptive time step procedures. Any
                # penetration that occurred in the previous increment is
                # adjusted for equilibrium at the start of the next increment.
                # This method may require smaller time steps then the other
                # methods.
                # 3 Uses contact procedure which does not require increment
                # splitting (3 is not available for dynamics). If a run does not
                # converge due to an excessive number of “iterative
                # penetration checking” messages, ISPLIT=2 may help,
                # however the time steps may need to be set to smaller values.
                value = integer_or_blank(card, i, f'{param}_value{j}', default=3)
                assert value in [0, 1, 2, 3], f'ISPLIT must be [0, 1, 2, 3]; ISPLIT={value}'
            elif param == 'ICSEP':
                # ICSEP Flag to control separation. Not used for segment-to-segment contact. (Integer > 0;
                # Default = 0)
                # 0 The node separates and an iteration occurs if the force on the node is
                # greater than the separation force.
                # 1 If a node which was in contact at the end of the previous increment
                # has a force greater than the separation force, the node does NOT
                # separate in this increment, but separates at the beginning of the next
                # increment.
                # 2 If a new node comes into contact during this increment, it is not
                # allowed to separate during this increment, prevents chattering.
                # 3 Both 1 and 2 are in effect.
                value = integer_or_blank(card, i, f'{param}_value{j}', default=0)
                assert value in [0, 1, 2, 3], f'ICSEP must be [0, 1, 2, 3]; ICSEP={value}'
            elif param == 'LINQUAD':
                # LINQUAD (2,14)
                # Higher order element contact flag (Integer; Default = 1).
                # 1 The outer boundary of a contact body is described by the
                # corner nodes only and mid-side nodes can’t come into
                # contact.
                # -1 The outer boundary is described by a quadratic field and
                # both corner and mid-side nodes are considered in contact.
                # If this flag is set to -1 and IBSEP is blank, IBSEP will be reset
                # to 2. This option is only available with Marc 2003 and
                # subsequent releases.
                value = integer_or_blank(card, i, f'{param}_value{j}', default=0)
                assert value in [-1, 1], f'LINQUAD must be [-1, 1]; LINQUAD={value}'
            elif param == 'LINCNT':
                value = integer_or_blank(card, i, f'{param}_value{j}', default=0)
                assert value in [-1, 0, 1], f'LINCNT must be [-1, 0, 1]; LINQUAD={value}'
            elif param == 'METHOD':
                #METHOD Flag to select Contact methods. (Character)
                #NODESURF Node to segment contact. (Default)
                #SEGTOSEG Segment to segment contact.
                value = string_choice_or_blank(card, i, f'{param}_value{j}',
                                       ('NODESURF', 'SEGTOSEG', 'SEGSMALL', 'SEGLARGE'),
                                        default='NODESURF')
            elif param == 'MAXENT':
                # MAXENT (2,2) Maximum number of entities created for any contact body. (Integer >
                # 0 or blank; default is max element number or 1.5 times the number of
                # nodes whichever is smaller)
                value = integer_or_blank(card, i, f'{param}_value{j}')
            elif param == 'MAXNOD':
                # MAXNOD (2,3) Maximum number of nodes that lie on the periphery of any
                # deformable contact body. (Integer > 0 or blank; default is the number
                # of nodes)
                value = integer_or_blank(card, i, f'{param}_value{j}')
            elif param == 'NBODIES':
                # NBODIES (2,1) Number of contact bodies defined in the analysis. (Integer > 0 or blank)
                value = integer(card, i, f'{param}_value{j}')
            elif param == 'NLGLUE':
                #(SOLs 101 and 400 only)
                #If all slave's for the BCTABLE or BCONPRG corresponding to the first loadcase
                #(first subcase and first step) contain IGLUE >0, permanent glued contact with small
                #rotation condition will be used for all SLAVE entries in all subcases and all steps
                #unless BCPARA,0,,1 is specified. If IGLUE < 0 exists, permanent glued
                value = integer_or_blank(card, i, f'{param}_value{j}', default=1)
                assert value in [0, 1], f'NLGLUE must be [0, 1]; NLGLUE={value}'
            elif param == 'SLDLMT':
                #SLDLMT Maximum allowed sliding distance, beyond it the contact segments are
                # redefined, for segment to segment contact analysis with large deformation.
                # (Real ≥ 0.0; Default = 0.0) See remark 9.
                value = double_or_blank(card, i, f'{param}_value{j}', default=0.0)
            elif param == 'SEGSYM':
                #Specify symmetric or non-symmetric friction matrix in segment to segment contact
                #analysis. (Integer 0 = symmetric matrix or 1 = non-symmetric matrix; Default = 0)
                value = integer_or_blank(card, i, f'{param}_value{j}', default=1)
                assert value in [0, 1], f'SEGSYM must be [0, 1]; SEGSYM={value}'
            elif param == 'THKOFF':
                #Ignore thickness from the tolerance used by ISEARCH=2 in node-to-surface contact
                #or from the characteristic length (for PENALT and AUGDIST) in segment-tosegment
                #contact. (Integer 0 = do not ignore thickness or 1 = remove thickness;
                #Default = 0)
                value = integer_or_blank(card, i, f'{param}_value{j}', default=1)
                assert value in [0, 1], f'THKOFF must be [0, 1]; THKOFF={value}'
            else:
                raise NotImplementedError(f'param={param} card={card}')

            params[param] = value
            i += 1
            j += 1

            ivalue_line += 1
            if ivalue_line == 3:
                if i == 8:
                    ivalue_line = 0
                    #print('*1', i, j, ivalue_line)
                    i += 1
                else:
                    #print('*2', i, j, ivalue_line)
                    i += 2
        return BCPARA(csid, params, comment=comment)

    def raw_fields(self):
        fields = ['BCPARA', self.csid]
        ivalue = 1
        ivalue_line = 1
        for key, value in sorted(self.params.items()):
            #print(ivalue, key, value)
            fields.append(key)
            fields.append(value)
            if ivalue_line == 3:
                if ivalue == 3:
                    #print('*')
                    fields.append(None)
                else:
                    #print('**')
                    fields.append(None)
                    fields.append(None)
                ivalue_line = 0
            ivalue += 1
            ivalue_line += 1
        return fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)


class BCTPARM(BaseCard):
    """
    Contact Parameters (SOLs 101, 103, 111, 112, and 401).
    Control parameters for the contact algorithm.

    +---------+--------+--------+--------+--------+--------+---------+--------+
    |    1    |   2    |    3   |   4    |   5    |   6    |    7    |    8   |
    +=========+========+========+========+========+========+=========+========+
    | BCTPARM | CSID   | Param1 | Value1 | Param2 | Value2 | Param3  | Value3 |
    +---------+--------+--------+--------+--------+--------+---------+--------+
    |         | Param4 | Value4 | Param5 | Value5 |  etc.  |         |        |
    +---------+--------+--------+--------+--------+--------+---------+--------+
    | BCTPARM |   1    | PENN   |  10.0  |  PENT  |  0.5   |   CTOL  | 0.001  |
    +---------+--------+--------+--------+--------+--------+---------+--------+
    |         | SHLTHK |   1    |        |        |        |         |        |
    +---------+--------+--------+--------+--------+--------+---------+--------+

    """
    type = 'BCTPARM'

    @classmethod
    def _init_from_empty(cls):
        csid = 1
        params = {'CSTIFF' : 1}
        return BCTPARM(csid, params, comment='')

    def _finalize_hdf5(self, encoding: str):
        """hdf5 helper function"""
        keys, values = self.params
        self.params = {key : value for key, value in zip(keys, values)}

    def __init__(self, csid: int,
                 params: dict[str, int | float],
    comment: str=''):
        """
        Creates a BCTPARM card

        Parameters
        ----------
        csid : int
            Contact set ID. Parameters defined in this command apply to
            contact set CSID defined by a BCTSET entry. (Integer > 0)
        params : dict[key] : value
            the optional parameters
        comment : str; default=''
            a comment for the card

        """
        if comment:
            self.comment = comment

        #: Contact set ID. Parameters defined in this command apply to
        #: contact set CSID defined by a BCTSET entry. (Integer > 0)
        self.csid = csid
        self.params = params

    def _finalize_hdf5(self, encoding: str) -> None:
        keys = self.params[0]
        values = self.params[1]
        self.params = {}
        for key, value in zip(keys, values):
            if isinstance(value, bytes):
                value = value.decode(encoding)
            self.params[key] = value

    @classmethod
    def add_card(cls, card: BDFCard, comment: str=''):
        """
        Adds a BCTPARM card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        csid = integer(card, 1, 'csid')
        params = read_param_dict(card, cls._get_param_value_from_card)
        return BCTPARM(csid, params, comment=comment)

    def _get_param_value_from_card(
            card: BDFCard,
            i: int, j: int) -> tuple[str, int | float]:
        param = string(card, i, 'param%s' % j)
        if param == 'TYPE' and 0:
            value = integer_or_blank(card, i+1, f'value{j:d}', default=0)
            assert value in [0, 1, 2], 'TYPE must be [0, 1, 2]; TYPE=%r' % value
        elif param == 'PENN':
            # PENN 10.0
            value = double(card, i+1, 'value%s' % j)
        elif param == 'PENT':
            # PENT 0.5
            value = double(card, i+1, 'value%s' % j)
        elif param == 'CTOL':
            # CTOL 10.0
            value = double(card, i+1, 'value%s' % j)
        elif param == 'SHLTHK':
            # SHLTHK 1
            value = integer(card, i+1, 'value%s' % j)
        # elif param == 'TYPE': # NX
        # value = string_or_blank(card, i, 'value%s' % j, default='FLEX').upper()
        # assert value in ['FLEX', 'RIGID', 'COATING'], 'TYPE must be [FLEX, RIGID, COATING.]; CSTIFF=%r' % value

        # elif param == 'NSIDE':
        # value = integer_or_blank(card, i+1, 'value%s' % j, default=1)
        # assert value in [1, 2], 'NSIDE must be [1, 2]; NSIDE=%r' % value
        # elif param == 'TBIRTH':
        # value = double_or_blank(card, i+1, 'value%s' % j, default=0.0)
        # elif param == 'TDEATH':
        # value = double_or_blank(card, i+1, 'value%s' % j, default=0.0)
        # elif param == 'INIPENE':
        # value = integer_or_blank(card, i+1, 'value%s' % j, default=0)
        # assert value in [0, 1, 2, 3], 'INIPENE must be [0, 1, 2]; INIPENE=%r' % value
        # elif param == 'PDEPTH':
        # value = double_or_blank(card, i+1, 'value%s' % j, default=0.0)
        # elif param == 'SEGNORM':
        # value = integer_or_blank(card, i+1, 'value%s' % j, default=0)
        # assert value in [-1, 0, 1], 'SEGNORM must be [-1, 0, 1]; SEGNORM=%r' % value
        # elif param == 'OFFTYPE':
        # value = integer_or_blank(card, i+1, 'value%s' % j, default=0)
        # assert value in [0, 1, 2], 'OFFTYPE must be [0, 1, 2]; OFFTYPE=%r' % value
        # elif param == 'OFFSET':
        # value = double_or_blank(card, i+1, 'value%s' % j, default=0.0)
        # elif param == 'TZPENE':
        # value = double_or_blank(card, i+1, 'value%s' % j, default=0.0)

        # elif param == 'CSTIFF':
        # value = integer_or_blank(card, i+1, 'value%s' % j, default=0)
        # assert value in [0, 1], 'CSTIFF must be [0, 1]; CSTIFF=%r' % value
        # elif param == 'TIED':
        # value = integer_or_blank(card, i+1, 'value%s' % j, default=0)
        # assert value in [0, 1], 'TIED must be [0, 1]; TIED=%r' % value
        # elif param == 'TIEDTOL':
        # value = double_or_blank(card, i+1, 'value%s' % j, default=0.0)
        # elif param == 'EXTFAC':
        # value = double_or_blank(card, i+1, 'value%s' % j, default=0.001)
        # assert 1.0E-6 <= value <= 0.1, 'EXTFAC must be 1.0E-6 < EXTFAC < 0.1; EXTFAC=%r' % value
        else:
            # FRICMOD, FPARA1/2/3/4/5, EPSN, EPST, CFACTOR1, PENETOL
            # NCMOD, TCMOD, RFORCE, LFORCE, RTPCHECK, RTPMAX, XTYPE
            # ...
            value = integer_double_or_blank(card, i+1, 'value%s' % j)
            assert value is not None, '%s%i must not be None' % (param, j)
        return param, value

    def raw_fields(self) -> list[str | int | float]:
        fields = ['BCTPARM', self.csid]
        fields.extend(params_listify(self.params))
        return fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)

class BCTPARA(BaseCard):
    """
    Defines parameters for a surface-to-surface contact region.

    +---------+--------+--------+--------+--------+--------+---------+--------+
    |    1    |   2    |    3   |   4    |   5    |   6    |    7    |    8   |
    +=========+========+========+========+========+========+=========+========+
    | BCTPARA | CSID   | Param1 | Value1 | Param2 | Value2 | Param3  | Value3 |
    +---------+--------+--------+--------+--------+--------+---------+--------+
    |         | Param4 | Value4 | Param5 | Value5 |  etc.  |         |        |
    +---------+--------+--------+--------+--------+--------+---------+--------+
    | BCTPARA |   1    | TYPE   |   0    | NSIDE  |   2    | SEGNORM |  -1    |
    +---------+--------+--------+--------+--------+--------+---------+--------+
    |         | CSTIFF |   1    | OFFSET | 0.015  |        |         |        |
    +---------+--------+--------+--------+--------+--------+---------+--------+

    """
    type = 'BCTPARA'

    @classmethod
    def _init_from_empty(cls):
        csid = 1
        params = {'CSTIFF' : 1}
        return BCTPARA(csid, params, comment='')

    def _finalize_hdf5(self, encoding):
        """hdf5 helper function"""
        keys, values = self.params
        self.params = {key : value for key, value in zip(keys, values)}

    def __init__(self, csid: int,
                 params: dict[str, Any], comment: str=''):
        """
        Creates a BCTPARA card

        Parameters
        ----------
        csid : int
            Contact set ID. Parameters defined in this command apply to
            contact set CSID defined by a BCTSET entry. (Integer > 0)
        params : dict[key] : value
            the optional parameters
        comment : str; default=''
            a comment for the card

        """
        if comment:
            self.comment = comment

        #: Contact set ID. Parameters defined in this command apply to
        #: contact set CSID defined by a BCTSET entry. (Integer > 0)
        self.csid = csid
        self.params = params

    @classmethod
    def add_card(cls, card: BDFCard, comment: str=''):
        """
        Adds a BCTPARA card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        csid = integer(card, 1, 'csid')

        all_params = [
            'TYPE', 'NSIDE', 'TBIRTH', 'TDEATH', 'INIPENE',
            'PDEPTH', 'SEGNORM', 'OFFTYPE', 'OFFSET',
            'TZPENE', 'CSTIFF', 'TIED', 'TIEDTOL', 'EXTFAC']
        all_params.sort()
        params = read_param_dict(card, cls._get_param_value_from_card)
        for param in params:
            if param not in all_params:
                warnings.warn(f'param={param!r} is not allowed; allowed={all_params}\n{str(card)}')

        return BCTPARA(csid, params, comment=comment)

    @staticmethod
    def _get_param_value_from_card(card: BDFCard,
                                   i: int, j: int) -> tuple[str, int | float]:
        param = string(card, i, f'param{j}')
        # assert param in all_params, f'param={param!r} is not allowed; allowed={all_params}'
        i += 1
        if param == 'TYPE':
            value = integer_or_blank(card, i, 'value%s' % j, default=0)
            assert value in [0, 1, 2], 'TYPE must be [0, 1, 2]; TYPE=%r' % value
        #elif param == 'TYPE': # NX
            #value = string_or_blank(card, i, 'value%s' % j, 'FLEX').upper()
            #assert value in ['FLEX', 'RIGID', 'COATING'], 'TYPE must be [FLEX, RIGID, COATING.]; CSTIFF=%r' % value

        elif param == 'NSIDE':
            value = integer_or_blank(card, i, 'value%s' % j, default=1)
            assert value in [1, 2], 'NSIDE must be [1, 2]; NSIDE=%r' % value
        elif param == 'TBIRTH':
            value = double_or_blank(card, i, 'value%s' % j, default=0.0)
        elif param == 'TDEATH':
            value = double_or_blank(card, i, 'value%s' % j, default=0.0)
        elif param == 'INIPENE':
            value = integer_or_blank(card, i, 'value%s' % j, default=0)
            assert value in [0, 1, 2, 3], 'INIPENE must be [0, 1, 2]; INIPENE=%r' % value
        elif param == 'PDEPTH':
            value = double_or_blank(card, i, 'value%s' % j, default=0.0)
        elif param == 'SEGNORM':
            value = integer_or_blank(card, i, 'value%s' % j, default=0)
            assert value in [-1, 0, 1], 'SEGNORM must be [-1, 0, 1]; SEGNORM=%r' % value
        elif param == 'OFFTYPE':
            value = integer_or_blank(card, i, 'value%s' % j, default=0)
            assert value in [0, 1, 2], 'OFFTYPE must be [0, 1, 2]; OFFTYPE=%r' % value
        elif param == 'OFFSET':
            value = double_or_blank(card, i, 'value%s' % j, default=0.0)
        elif param == 'TZPENE':
            value = double_or_blank(card, i, 'value%s' % j, default=0.0)

        elif param == 'CSTIFF':
            value = integer_or_blank(card, i, 'value%s' % j, default=0)
            assert value in [0, 1], 'CSTIFF must be [0, 1]; CSTIFF=%r' % value
        elif param == 'TIED':
            value = integer_or_blank(card, i, 'value%s' % j, default=0)
            assert value in [0, 1], 'TIED must be [0, 1]; TIED=%r' % value
        elif param == 'TIEDTOL':
            value = double_or_blank(card, i, 'value%s' % j, default=0.0)
        elif param == 'EXTFAC':
            value = double_or_blank(card, i, 'value%s' % j, default=0.001)
            assert 1.0E-6 <= value <= 0.1, 'EXTFAC must be 1.0E-6 < EXTFAC < 0.1; EXTFAC=%r' % value
        else:
            # FRICMOD, FPARA1/2/3/4/5, EPSN, EPST, CFACTOR1, PENETOL
            # NCMOD, TCMOD, RFORCE, LFORCE, RTPCHECK, RTPMAX, XTYPE
            # ...
            value = integer_double_or_blank(card, i, 'value%s' % j)
            assert value is not None, '%s%i must not be None' % (param, j)
        return param, value

    def raw_fields(self) -> list:
        fields = ['BCTPARA', self.csid]
        fields.extend(params_listify(self.params))
        return fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)


class BCTADD(BaseCard):
    """
    +--------+------+----+-------+----+----+----+----+----+
    |   1    |  2   | 3  |   4   |  5 | 6  |  7 | 8  |  9 |
    +========+======+====+=======+====+====+====+====+====+
    | BCTADD | CSID | SI |  S2   | S3 | S4 | S5 | S6 | S7 |
    +--------+------+----+-------+----+----+----+----+----+
    |        |  S8  | S9 |  etc. |    |    |    |    |    |
    +--------+------+----+-------+----+----+----+----+----+

    Remarks:
    1. To include several contact sets defined via BCTSET entries in a model,
       BCTADD must be used to combine the contact sets. CSID in BCTADD is
       then selected with the Case Control command BCSET.
    2. Si must be unique and may not be the identification of this or any other
       BCTADD entry.

    """
    type = 'BCTADD'

    @classmethod
    def _init_from_empty(cls):
        csid = 1
        contact_sets = [1, 2]
        return BCTADD(csid, contact_sets, comment='')

    def __init__(self, csid, contact_sets, comment=''):
        if comment:
            self.comment = comment
        #: Contact set identification number. (Integer > 0)
        self.csid = csid

        #: Identification numbers of contact sets defined via BCTSET entries.
        #: (Integer > 0)
        self.contact_sets = contact_sets

    @classmethod
    def add_card(cls, card: BDFCard, comment: str=''):
        """
        Adds a BCTADD card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        csid = integer(card, 1, 'csid')
        contact_sets = []

        i = 1
        j = 1
        while i < card.nfields:
            contact_set = integer(card, i, 'S%i' % j)
            contact_sets.append(contact_set)
            i += 1
            j += 1
        return BCTADD(csid, contact_sets, comment=comment)

    def raw_fields(self):
        fields = ['BCTADD'] + self.contact_sets
        return fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)

class BGADD(BaseCard):
    """
    +-------+------+----+-------+----+----+----+----+----+
    |   1   |  2   | 3  |   4   |  5 | 6  |  7 | 8  |  9 |
    +=======+======+====+=======+====+====+====+====+====+
    | BGADD | GSID | SI |  S2   | S3 | S4 | S5 | S6 | S7 |
    +-------+------+----+-------+----+----+----+----+----+
    |       |  S8  | S9 |  etc. |    |    |    |    |    |
    +-------+------+----+-------+----+----+----+----+----+

    """
    type = 'BGADD'

    @classmethod
    def _init_from_empty(cls):
        glue_id = 1
        contact_sets = [1, 2]
        return BGADD(glue_id, contact_sets, comment='')

    def __init__(self, glue_id, contact_sets, comment=''):
        if comment:
            self.comment = comment
        #: Glue identification number. (Integer > 0)
        self.glue_id = glue_id

        #: Identification numbers of contact sets defined via BCTSET entries.
        #: (Integer > 0)
        self.contact_sets = contact_sets

    @classmethod
    def add_card(cls, card: BDFCard, comment: str=''):
        """
        Adds a BGADD card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        glue_id = integer(card, 1, 'glue_id')
        contact_sets = []

        i = 1
        j = 1
        while i < card.nfields:
            contact_set = integer(card, i, 'S%i' % j)
            contact_sets.append(contact_set)
            i += 1
            j += 1
        return BGADD(glue_id, contact_sets, comment=comment)

    def raw_fields(self):
        fields = ['BGADD'] + self.contact_sets
        return fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)

class BGSET(BaseCard):
    """
    +-------+------+------+------+---------+----+------+------+----+
    |   1   |  2   |  3   |   4  |    5    | 6  |  7   |   8  |  9 |
    +=======+======+======+======+=========+====+======+======+====+
    | BGSET | GSID | SID1 | TID1 | SDIST1  |    | EXT1 |      |    |
    +-------+------+------+------+---------+----+------+------+----+
    |       |      | SID2 | TID2 | SDIST2  |    | EXT2 |      |    |
    +-------+------+------+------+---------+----+------+------+----+
    """
    type = 'BGSET'

    @classmethod
    def _init_from_empty(cls):
        glue_id = 1
        sids = [1]
        tids = [1]
        sdists = [0.01]
        exts = [1.]
        return BGSET(glue_id, sids, tids, sdists, exts, comment='', sol=101)

    def __init__(self, glue_id: int,
                 sids: list[int],
                 tids: list[int],
                 sdists: list[float],
                 exts: list[float],
                 comment: str='', sol: int=101):
        if comment:
            self.comment = comment
        #: GSID Glue set identification number. (Integer > 0)
        self.glue_id = glue_id
        #: SIDi Source region (contactor) identification number for contact pair i.
        #: (Integer > 0)
        self.sids = sids

        #: TIDi Target region identification number for contact pair i. (Integer > 0)
        self.tids = tids

        #: SDISTi Search distance for glue regions (Real); (Default=10.0)
        self.sdists = sdists

        #: EXTi Extension factor for target region (SOLs 402 and 601 only).
        self.exts = exts

    @classmethod
    def add_card(cls, card: BDFCard, comment: str='', sol: int=101):
        glue_id = integer(card, 1, 'glue_id')
        sids = []
        tids = []
        sdists = []
        exts = []

        nfields = card.nfields
        i = 2
        j = 1
        while i < nfields:
            #SIDi Source region identification number for glue pair i. (Integer > 0)
            #TIDi Target region identification number for glue pair i. (Integer > 0)
            #SDISTi Search distance for glue regions (Real); (Default=10.0)
            #EXTi Extension factor for target region (SOLs 402 and 601 only).

            sids.append(integer(card, i, 'sid%s' % j))
            tids.append(integer(card, i + 1, 'tid%s' % j))
            sdists.append(double_or_blank(card, i + 2, 'fric%s' % j, default=0.0))
            #if sol == 101:
            exts.append(double_or_blank(card, i + 4, 'mind%s' % j, default=0.0))
            #else:
                #exts.append(None)
            i += 8
            j += 1
        return BGSET(glue_id, sids, tids, sdists, exts,
                     comment=comment, sol=sol)

    def raw_fields(self):
        fields = ['BGSET', self.glue_id]
        for sid, tid, sdist, ext in zip(self.sids, self.tids, self.sdists,
                                        self.exts):
            fields += [sid, tid, sdist, None, ext, None, None, None]
        return fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)


def params_listify(params: dict[str, int | float],
                   ) -> list[str | int | float]:
    """
    Annoying function for contact params.

    3 sets of params are allowed per line, but:
     - 1 blank on the first line
     - 2 blanon the 2nd/3rd/etc. line

    """
    i = 0
    fields = []
    for j, (key, value) in enumerate(sorted(params.items())):
        if j == 3:
            fields.append(None)
            i = 0
        elif j > 0 and i == 3:
            fields.append(None)
            fields.append(None)
            i = 0
        fields.append(key)
        fields.append(value)
        i += 1
    return fields


def read_param_dict(card: BDFCard,
                    read_func) -> dict[str, int | float]:
    """reading version of params_listify"""
    params = {}
    iline = 0
    ifield = 2
    iparam = 0
    while ifield < card.nfields:
        param, value = read_func(card, ifield, iparam+1)
        params[param] = value
        iparam += 1
        if iparam % 3 == 0:
            # every 3 params, reset the line
            ifield = 9 + 8 * iline
            iline += 1
        else:
            ifield += 2
    assert len(params) > 0, card
    return params
