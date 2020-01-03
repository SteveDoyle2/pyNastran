"""
All set cards are defined in this file.  This includes:

* bcs
  * RADM, RADBC
* views
  * VIEW, VIEW3D

"""
from __future__ import annotations
from typing import TYPE_CHECKING
import warnings

from pyNastran.utils.numpy_utils import integer_types, float_types
from pyNastran.bdf.field_writer_8 import set_blank_if_default, print_card_8
from pyNastran.bdf.field_writer_16 import print_card_16
from pyNastran.bdf.cards.base_card import BaseCard, expand_thru_by
from pyNastran.bdf.cards.collpase_card import collapse_thru_by
from pyNastran.bdf.bdf_interface.assign_type import (
    fields, integer, double, integer_or_blank, double_or_blank,
    integer_or_string, string_or_blank)

from pyNastran.bdf.cards.thermal.thermal import ThermalBC
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf import BDF


class RADM(ThermalBC):
    """
    Defines the radiation properties of a boundary element for heat transfer
    analysis

    """
    type = 'RADM'

    @classmethod
    def _init_from_empty(cls):
        radmid = 2
        absorb = 2.0
        emissivity = 1.0
        return RADM(radmid, absorb, emissivity, comment='')

    def __init__(self, radmid, absorb, emissivity, comment=''):
        ThermalBC.__init__(self)
        if comment:
            self.comment = comment

        #: Material identification number
        self.radmid = radmid

        self.absorb = absorb
        if isinstance(emissivity, float_types):
            self.emissivity = [emissivity]
        else:
            self.emissivity = emissivity

    def validate(self):
        assert self.radmid > 0, str(self)
        msg = ''
        if self.absorb is not None:
            if not 0. <= self.absorb <= 1.0:
                msg += 'absorb=%s not in range 0.0 <= absorb <= 1.0\n' % (self.absorb)
        for i, emissivityi in enumerate(self.emissivity):
            if not 0. <= emissivityi <= 1.0:
                msg += 'emissivity[%i]=%s\n' % (i, emissivityi)
        if msg:
            warnings.warn(msg + str(self))
            #raise RuntimeError(msg + str(self))

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a RADM card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        nfields = card.nfields
        radmid = integer(card, 1, 'radmid')
        absorb = double_or_blank(card, 2, 'absorb')
        emissivity = fields(double, card, 'emissivity', i=3, j=nfields)
        return RADM(radmid, absorb, emissivity, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        """
        Adds a RADM card from the OP2

        Parameters
        ----------
        data : List[varies]
            a list of fields defined in OP2 format
        comment : str; default=''
            a comment for the card

        """
        radmid, absorb = data[:2]
        emissivity = data[2:]
        return RADM(radmid, absorb, emissivity, comment=comment)

    #def cross_reference(self, model: BDF) -> None:
        #pass

    def raw_fields(self):
        list_fields = ['RADM', self.radmid, self.absorb] + self.emissivity
        return list_fields

    def repr_fields(self):
        list_fields = ['RADM', self.radmid, self.absorb] + self.emissivity
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        """
        The writer method used by BDF.write_card()

        Parameters
        -----------
        size : int; default=8
            the size of the card (8/16)

        """
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)


class RADBC(ThermalBC):
    """
    Specifies an CHBDYi element face for application of radiation boundary
    conditions

    """
    type = 'RADBC'

    @classmethod
    def _init_from_empty(cls):
        nodamb = 1
        famb = 1.0
        cntrlnd = 10
        eids = [1, 2]
        return RADBC(nodamb, famb, cntrlnd, eids, comment='')

    def __init__(self, nodamb, famb, cntrlnd, eids, comment=''):
        ThermalBC.__init__(self)
        if comment:
            self.comment = comment

        #: NODAMB Ambient point for radiation exchange. (Integer > 0)
        self.nodamb = nodamb

        #: Radiation view factor between the face and the ambient point.
        #: (Real > 0.0)
        self.famb = famb

        #: Control point for thermal flux load. (Integer > 0; Default = 0)
        self.cntrlnd = cntrlnd

        #: CHBDYi element identification number
        if isinstance(eids, int):
            eids = [eids]
        self.eids = expand_thru_by(eids)

        assert self.nodamb > 0
        assert self.famb > 0.0
        assert self.cntrlnd >= 0
        self.eids_ref = None

    def validate(self):
        min_eid = min(self.eids)
        if min_eid < 1:
            msg = 'min(eids)=%i' % min_eid
            warnings.warn(msg)

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a RADBC card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        nodamb = integer(card, 1, 'nodamb')
        famb = double(card, 2, 'famb')
        cntrlnd = integer_or_blank(card, 3, 'cntrlnd', 0)

        nfields = card.nfields
        eids = fields(integer_or_string, card, 'eid', i=4, j=nfields)
        return RADBC(nodamb, famb, cntrlnd, eids, comment=comment)

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        msg = ', which is required by RADBC pid=%s' % self.nodamb
        elems = []
        for eid in self.eids:
            elem = model.Element(eid, msg=msg)
            elems.append(elem)
        self.eids_ref = elems

    def Eids(self):
        if self.eids_ref is None:
            return self.eids
        eids = []
        for eid_ref in self.eids_ref:
            eids.append(eid_ref.eid)
        return eids

    def raw_fields(self):
        list_fields = (['RADBC', self.nodamb, self.famb, self.cntrlnd] +
                       self.Eids())
        return list_fields

    def repr_fields(self):
        cntrlnd = set_blank_if_default(self.cntrlnd, 0)
        eids = collapse_thru_by(self.Eids())
        list_fields = ['RADBC', self.nodamb, self.famb, cntrlnd] + eids
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        """
        The writer method used by BDF.write_card()

        Parameters
        -----------
        size : int; default=8
            the size of the card (8/16)

        """
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)

# Boundary Conditions
#-------------------------------------------------------
# View Factors
class VIEW(BaseCard):
    """
    Defines radiation cavity and shadowing for radiation
    view factor calculations.

    +------+-------+---------+-------+----+----+--------+
    |   1  |   2   |    3    |   4   | 5  |  6 |    7   |
    +======+=======+=========+=======+====+====+========+
    | VIEW | IVIEW | ICAVITY | SHADE | NB | NG | DISLIN |
    +------+-------+---------+-------+----+----+--------+
    | VIEW |   1   |    1    | BOTH  | 2  | 3  |  0.25  |
    +------+-------+---------+-------+----+----+--------+

    """
    type = 'VIEW'

    @classmethod
    def _init_from_empty(cls):
        iview = 2
        icavity = 2
        return VIEW(iview, icavity, shade='BOTH', nbeta=1, ngamma=1, dislin=0.0, comment='')

    def __init__(self, iview, icavity, shade='BOTH', nbeta=1, ngamma=1,
                 dislin=0.0, comment=''):
        """
        Creates a VIEW, which defines a 2D view factor

        Parameters
        ----------
        iview : int
            Identification number
        icavity : int
            Cavity identification number for grouping the radiant exchange faces of
            CHBDYi elements
        shade : str; default='BOTH'
            Shadowing flag for the face of CHBDYi element
            - NONE means the face can neither shade nor be shaded by other faces
            - KSHD means the face can shade other faces
            - KBSHD means the face can be shaded by other faces
            - BOTH means the face can both shade and be shaded by other faces
        nbeta / ngamma : int; default=1 / 1
            Subelement mesh size in the beta/gamma direction. (Integer > 0)
        dislin : float; default=0.0
            The displacement of a surface perpendicular to the surface
        comment : str; default=''
            a comment for the card

        """
        BaseCard.__init__(self)
        if comment:
            self.comment = comment

        #: Material identification number
        self.iview = iview
        self.icavity = icavity
        self.shade = shade
        self.nbeta = nbeta
        self.ngamma = ngamma
        self.dislin = dislin

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a VIEW card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        iview = integer(card, 1, 'iview')
        icavity = integer(card, 2, 'icavity')
        shade = string_or_blank(card, 3, 'shade', default='BOTH')
        nbeta = integer_or_blank(card, 4, 'nbeta', default=1)
        ngamma = integer_or_blank(card, 5, 'ngamma', default=1)
        dislin = double_or_blank(card, 6, 'dislin', default=0.0)
        return VIEW(iview, icavity, shade=shade, nbeta=nbeta, ngamma=ngamma,
                    dislin=dislin, comment=comment)

    #def cross_reference(self, model: BDF) -> None:
        #pass

    def raw_fields(self):
        list_fields = ['VIEW', self.iview, self.icavity, self.shade,
                       self.nbeta, self.ngamma, self.dislin]
        return list_fields

    def repr_fields(self):
        list_fields = ['VIEW', self.iview, self.icavity, self.shade,
                       self.nbeta, self.ngamma, self.dislin]
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        """
        The writer method used by BDF.write_card()

        Parameters
        -----------
        size : int; default=8
            the size of the card (8/16)

        """
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)


class VIEW3D(BaseCard):
    """
    View Factor Definition - Gaussian Integration Method

    Defines parameters to control and/or request the Gaussian Integration
    method of view factor calculation for a specified cavity.

    +--------+---------+------+------+------+------+--------+------+--------+
    |    1   |    2    |   3  |  4   |   5  |   6  |    7   |   8  |    9   |
    +========+=========+======+======+======+======+========+======+========+
    | VIEW3D | ICAVITY | GITB | GIPS | CIER | ETOL |  ZTOL  | WTOL | RADCHK |
    +--------+---------+------+------+------+------+--------+------+--------+
    | VIEW3D |    1    |   2  |   2  |   4  |      | 1.0E-6 |      |        |
    +--------+---------+------+------+------+------+--------+------+--------+

    """
    type = 'VIEW3D'

    @classmethod
    def _init_from_empty(cls):
        icavity = 2
        return VIEW3D(icavity, gitb=4, gips=4, cier=4, error_tol=0.1,
                      zero_tol=1e-10, warp_tol=0.01, rad_check=3, comment='')

    def __init__(self, icavity, gitb=4, gips=4, cier=4,
                 error_tol=0.1, zero_tol=1e-10, warp_tol=0.01,
                 rad_check=3, comment=''):
        """
        Creates a VIEW3D, which defines a 3D view factor

        Parameters
        ----------
        icavity : int
            Radiant cavity identification number on RADCAV entry. (Integer > 0)
        gitb : int; default=4
            Gaussian integration order to be implemented in calculating net
            effective view factors in the presence of third-body shadowing.
            (Integer 2, 3, 4, 5, 6 or 10)
        gips : int; default=4
            Gaussian integration order to be implemented in calculating net
            effective view factors in the presence of self-shadowing.
            (Integer 2, 3, 4, 5, 6 or 10)
        cier : int; default=4
            Discretization level used in the semi-analytic contour integration
            method. (1 < Integer < 20)
        error_tol : float; default=0.1
            Error estimate above which a corrected view factor is calculated
            using the semi-analytic contour integration method. (Real > 0.0)
        zero_tol : float; default=1e-10
            Assumed level of calculation below which the numbers are considered
            to be zero. (Real > 0.0)
        warp_tol : float; default=0.01
            Assumed degree of warpage above which the actual value of will be
            calculated. (0.0 < Real < 1.0)
        rad_check : int; default=3
             Type of diagnostic output desired for the radiation exchange surfaces.
        comment : str; default=''
            a comment for the card

        """
        BaseCard.__init__(self)
        if comment:
            self.comment = comment

        #: Material identification number
        self.icavity = icavity
        self.gitb = gitb
        self.gips = gips
        self.cier = cier
        self.error_tol = error_tol
        self.zero_tol = zero_tol
        self.warp_tol = warp_tol
        self.rad_check = rad_check

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a VIEW3D card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        icavity = integer(card, 1, 'icavity')
        gitb = integer_or_blank(card, 2, 'gitb', 4)
        gips = integer_or_blank(card, 3, 'gips', 4)
        cier = integer_or_blank(card, 4, 'cier', 4)
        error_tol = double_or_blank(card, 5, 'error_tol', 0.1)
        zero_tol = double_or_blank(card, 6, 'zero_tol', 1e-10)
        warp_tol = double_or_blank(card, 7, 'warp_tol', 0.01)
        rad_check = integer_or_blank(card, 8, 'rad_check', 3)
        return VIEW3D(icavity, gitb=gitb, gips=gips, cier=cier,
                      error_tol=error_tol, zero_tol=zero_tol, warp_tol=warp_tol,
                      rad_check=rad_check, comment=comment)

    #def cross_reference(self, model: BDF) -> None:
        #pass

    def raw_fields(self):
        list_fields = ['VIEW3D', self.icavity, self.gitb, self.gips, self.cier,
                       self.error_tol, self.zero_tol, self.warp_tol, self.rad_check]
        return list_fields

    def repr_fields(self):
        list_fields = ['VIEW3D', self.icavity, self.gitb, self.gips, self.cier,
                       self.error_tol, self.zero_tol, self.warp_tol, self.rad_check]
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        """
        The writer method used by BDF.write_card()

        Parameters
        -----------
        size : int; default=8
            the size of the card (8/16)

        """
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)

#--------------------------------------------------------
class RADCAV(ThermalBC):
    """
    Identifies the characteristics of each radiant enclosure.

    +--------+---------+--------+------+---------+--------+-------+------+--------+
    |    1   |    2    |    3   |  4   |     5   |    6   |   7   |   8  |    9   |
    +========+=========+========+========+=======+========+=======+======+========+
    | RADCAV | ICAVITY | ELEAMB | SHADOW | SCALE | PRTPCH | NFECI | RMAX |        |
    +--------+---------+--------+------+---------+--------+-------+------+--------+
    |        |  SET11  |  SET12 |  SET21 | SET22 |  SET31 | SET32 | etc. |        |
    +--------+---------+--------+------+---------+--------+-------+------+--------+
    | RADCAV |    1    |    1   |        |       |        |       | .99  |        |
    +--------+---------+--------+------+---------+--------+-------+------+--------+
    |        |    3    |    5   |    4   |   5   |    7   |   5   |      |        |
    +--------+---------+--------+------+---------+--------+-------+------+--------+

    """
    type = 'RADCAV'

    @classmethod
    def _init_from_empty(cls):
        icavity = 2
        sets = [1, 2, 3]
        return RADCAV(icavity, sets, ele_amb=None, shadow='YES', scale=0.0,
                      prtpch=None, nefci=None, rmax=0.1, ncomp=32, comment='')

    def __init__(self, icavity, sets, ele_amb=None,
                 shadow='YES', scale=0.0, prtpch=None,
                 nefci=None, rmax=0.1, ncomp=32, comment=''):
        ThermalBC.__init__(self)
        if comment:
            self.comment = comment

        self.icavity = icavity
        self.ele_amb = ele_amb
        self.shadow = shadow
        self.scale = scale
        self.prtpch = prtpch
        self.nefci = nefci
        self.rmax = rmax
        self.ncomp = ncomp

        if isinstance(sets, integer_types):
            sets = [sets]
        self.sets = sets

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a RADCAV card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        icavity = integer(card, 1, 'icavity')
        ele_amb = integer_or_blank(card, 2, 'ele_amb')
        shadow = string_or_blank(card, 3, 'shadow', 'YES')
        scale = double_or_blank(card, 4, 'scale', 0.0)
        prtpch = integer_or_blank(card, 5, 'prtpch')
        nefci = string_or_blank(card, 6, 'nefci')
        rmax = double_or_blank(card, 7, 'rmax', 1.0)
        ncomp = integer_or_blank(card, 8, 'ncomp', 32)

        sets = fields(integer, card, 'set', i=9, j=card.nfields)
        return RADCAV(icavity, sets, ele_amb=ele_amb,
                      shadow=shadow, scale=scale, prtpch=prtpch,
                      nefci=nefci, rmax=rmax, ncomp=ncomp, comment=comment)

    #@classmethod
    #def add_op2_data(cls, data, comment=''):
        #"""
        #Adds a RADM card from the OP2

        #Parameters
        #----------
        #data : List[varies]
            #a list of fields defined in OP2 format
        #comment : str; default=''
            #a comment for the card
        #"""
        #radmid, absorb = data[:2]
        #emissivity = data[2:]

    #def cross_reference(self, model: BDF) -> None:
        #pass

    def raw_fields(self):
        list_fields = ['RADCAV', self.icavity, self.ele_amb, self.shadow, self.scale,
                       self.prtpch, self.nefci, self.rmax, self.ncomp] + self.sets
        return list_fields

    def repr_fields(self):
        list_fields = ['RADCAV', self.icavity, self.ele_amb, self.shadow, self.scale,
                       self.prtpch, self.nefci, self.rmax, self.ncomp] + self.sets
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        """
        The writer method used by BDF.write_card()

        Parameters
        -----------
        size : int; default=8
            the size of the card (8/16)

        """
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)

class RADLST(ThermalBC):
    """
    Identifies the characteristics of each radiant enclosure.

    +--------+---------+--------+------+-------+--------+-------+------+--------+
    |    1   |    2    |    3   |  4   |   5   |    6   |   7   |   8  |    9   |
    +========+=========+========+======+=======+========+=======+======+========+
    | RADLST | ICAVITY | MTXTYP | EID1 |  EID2 |  EID3  |  EID4 | EID5 |  EID6  |
    +--------+---------+--------+------+-------+--------+-------+------+--------+
    |        |   EID7  |  etc.  |      |       |        |       |      |        |
    +--------+---------+--------+------+-------+--------+-------+------+--------+
    | RADLST |    3    |    5   |  4   |   5   |    7   |   5   |      |        |
    +--------+---------+--------+------+-------+--------+-------+------+--------+

    """
    type = 'RADCAV'

    @classmethod
    def _init_from_empty(cls):
        icavity = 2
        eids = [1, 2, 3]
        return RADLST(icavity, eids, matrix_type=1, comment='')

    def __init__(self, icavity, eids, matrix_type=1, comment=''):
        ThermalBC.__init__(self)
        if comment:
            self.comment = comment

        self.icavity = icavity
        self.matrix_type = matrix_type
        assert isinstance(matrix_type, integer_types), matrix_type
        if isinstance(eids, integer_types):
            eids = [eids]
        self.eids = eids

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a RADLST card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        icavity = integer(card, 1, 'icavity')
        matrix_type = integer_or_blank(card, 2, 'matrix_type', 1)

        eids = fields(integer, card, 'eid', i=3, j=card.nfields)
        return RADLST(icavity, eids, matrix_type=matrix_type, comment=comment)

    def raw_fields(self):
        list_fields = ['RADLST', self.icavity, self.matrix_type] + self.eids
        return list_fields

    #def repr_fields(self):
        #list_fields = ['RADCAV', self.icavity, self.ele_amb, self.shadow, self.scale,
                       #self.prtpch, self.nefci, self.rmax, self.ncomp] + self.sets
        #return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        """
        The writer method used by BDF.write_card()

        Parameters
        -----------
        size : int; default=8
            the size of the card (8/16)

        """
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)

class RADMTX(ThermalBC):
    """
    Provides the Fji=Aj*fji exchange factors for all the faces of a
    radiation enclosure specified in the corresponding RADLST entry.

    +--------+---------+--------+--------+--------+--------+--------+--------+--------+
    |    1   |    2    |    3   |   4    |    5   |    6   |    7   |    8   |    9   |
    +========+=========+========+========+========+========+========+========+========+
    | RADMTX | ICAVITY | INDEX  |  Fi,j  | Fi+1,j | Fi+2,j | Fi+3,j | Fi+4,j | Fi+5,j |
    +--------+---------+--------+--------+--------+--------+--------+--------+--------+
    |        | Fi+6,j  |  etc.  |        |        |        |        |        |        |
    +--------+---------+--------+--------+--------+--------+--------+--------+--------+
    | RADMTX |    2    |    1   |  0.0   |  0.1   |   0.2  |   0.2  |   0.3  |   0.2  |
    +--------+---------+--------+--------+--------+--------+--------+--------+--------+

    """
    type = 'RADMTX'

    @classmethod
    def _init_from_empty(cls):
        icavity = 1
        index = 2
        exchange_factors = [1., 2.]
        return RADMTX(icavity, index, exchange_factors, comment='')

    def __init__(self, icavity, index, exchange_factors, comment=''):
        ThermalBC.__init__(self)
        if comment:
            self.comment = comment

        self.icavity = icavity
        assert isinstance(index, integer_types), type(index)
        self.index = index

        if isinstance(exchange_factors, float_types):
            exchange_factors = [exchange_factors]
        self.exchange_factors = exchange_factors

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a RADMTX card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        icavity = integer(card, 1, 'icavity')
        index = integer(card, 2, 'index')
        exchange_factors = fields(double, card, 'eid', i=3, j=card.nfields)
        return RADMTX(icavity, index, exchange_factors, comment=comment)

    def raw_fields(self):
        assert isinstance(self.index, integer_types), type(self.index)
        list_fields = ['RADMTX', self.icavity, self.index] + self.exchange_factors
        return list_fields

    #def repr_fields(self):
        #list_fields = ['RADCAV', self.icavity, self.ele_amb, self.shadow, self.scale,
                       #self.prtpch, self.nefci, self.rmax, self.ncomp] + self.sets
        #return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        """
        The writer method used by BDF.write_card()

        Parameters
        -----------
        size : int; default=8
            the size of the card (8/16)

        """
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)
