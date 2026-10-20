"""
Defines the following contact cards:

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
 - BCONP (see contact.py)
 - BCONPRG (todo)
 - BCONPRP (todo)

 - BCONUDS (todo)
 - BCPARA (see contact.py)
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
 - BSURF (see contact.py)
 - BWIDTH (todo)
 - DYPARAM,CONTACT (todo)

"""
from __future__ import annotations
from typing import Any, TYPE_CHECKING

from pyNastran.utils.numpy_utils import integer_types
from pyNastran.bdf.cards.base_card import BaseCard
from pyNastran.bdf.bdf_interface.assign_type import (
    integer, integer_or_blank, integer_string_or_blank, double_or_blank,
    integer_double_or_blank, string_or_blank, string_choice_or_blank,
    double, blank)
from pyNastran.bdf.field_writer_8 import print_card_8
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf import BDF
    from pyNastran.bdf.bdf_interface.bdf_card import BDFCard


def get_table3d(model: BDF, fric: int) -> TABLE3D:
    return model.table3d[fric]


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
                      dim=dim, behav=behav, istype=istype, fric=fric, idispl=idispl)

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
            fric_ref = get_table3d(model, self.fric)

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
    val = string_or_blank(card, i, f'{word}_word')
    valuei = isinstance(val, str)
    #print('i', i, 1)
    while not isinstance(valuei, str) and i < len(card):
        i += 8
        valuei = string_or_blank(card, i, f'{word}_word')
        #print('i', i, valuei)
    values = card[i0:i]
    return i, values
