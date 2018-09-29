#pylint: disable=C0103,C0111,R0902
"""
All ungrouped properties are defined in this file.  This includes:
 * PCOMPS (SolidProperty)
 * PLSOLID (SolidProperty)
 * PSOLID (SolidProperty)
"""
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
import numpy as np

from pyNastran.utils.numpy_utils import integer_types
from pyNastran.bdf.field_writer_8 import set_blank_if_default
from pyNastran.bdf.cards.base_card import Property
from pyNastran.bdf.bdf_interface.assign_type import (
    integer, integer_or_blank, double, double_or_blank, string_or_blank,
    integer_string_or_blank)
from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.field_writer_16 import print_card_16

class SolidProperty(Property):
    def __init__(self):
        Property.__init__(self)
        self.mid_ref = None

    def Rho(self):
        return self.mid_ref.rho


class PLSOLID(SolidProperty):
    """
    Defines a fully nonlinear (i.e., large strain and large rotation)
    hyperelastic solid element.

    +---------+-----+-----+-----+
    |    1    |  2  |  3  |  4  |
    +=========+=====+=====+=====+
    | PLSOLID | PID | MID | STR |
    +---------+-----+-----+-----+
    | PLSOLID |  20 |  21 |     |
    +---------+-----+-----+-----+
    """
    type = 'PLSOLID'
    _field_map = {
        1: 'pid', 2:'mid', 3:'str',
    }

    def __init__(self, pid, mid, stress_strain='GRID', ge=0., comment=''):
        """
        Creates a PLSOLID card

        Parameters
        ----------
        pid : int
            property id
        mid : int
            material id
        stress_strain : str
            Location of stress and strain output
            valid types = {GRID, GAUSS}
        ge : float; default=0.
            damping coefficient
        comment : str; default=''
            a comment for the card
        """
        SolidProperty.__init__(self)
        if stress_strain == 'GAUS':
            stress_strain = 'GAUSS'

        if comment:
            self.comment = comment
        #: Property ID
        self.pid = pid
        #: Material ID
        self.mid = mid
        #: Location of stress and strain output
        self.stress_strain = stress_strain

        self.ge = ge
        assert isinstance(pid, integer_types), type(pid)
        assert isinstance(mid, integer_types), type(mid)
        self.mid_ref = None

    def validate(self):
        if self.stress_strain not in ['GRID', 'GAUSS']:
            raise RuntimeError('STR="%s" doesnt have a valid stress/strain '
                               'output value set; valid=["GRID", "GAUSS"]\n'
                               % self.stress_strain)

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a PLSOLID card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        pid = integer(card, 1, 'pid')
        mid = integer(card, 2, 'mid')
        stress_strain = string_or_blank(card, 3, 'stress_strain', 'GRID')
        assert len(card) <= 4, 'len(PLSOLID card) = %i\ncard=%s' % (len(card), card)
        return PLSOLID(pid, mid, stress_strain, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        """
        Adds a PLSOLID card from the OP2

        Parameters
        ----------
        data : List[varies]
            a list of fields defined in OP2 format
        comment : str; default=''
            a comment for the card
        """
        pid = data[0]
        mid = data[1]
        ge = data[2]
        stress_strain = data[3]
        return PLSOLID(pid, mid, stress_strain, ge, comment=comment)

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ', which is required by PLSOLID pid=%s' % self.pid
        self.mid_ref = model.HyperelasticMaterial(self.mid, msg) # MATHP, MATHE

    def uncross_reference(self):
        self.mid = self.Mid()
        self.mid_ref = None

    def Rho(self):
        """
        Returns the density
        """
        return self.mid_ref.rho

    def _verify(self, xref):
        pass

    def raw_fields(self):
        stress_strain = set_blank_if_default(self.stress_strain, 'GRID')
        fields = ['PLSOLID', self.pid, self.Mid(), stress_strain]
        return fields

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)


class PCOMPS(SolidProperty):
    type = 'PCOMPS'
    _field_map = {
        #1: 'pid', 2:'mid', 3:'cordm', 4:'integ', 5:'stress',
        #6:'isop', 7:'fctn',
    }

    def __init__(self, pid, global_ply_ids, mids, thicknesses, thetas,
                 cordm=0, psdir=13, sb=None, nb=None, tref=0.0, ge=0.0,
                 failure_theories=None, interlaminar_failure_theories=None,
                 souts=None, comment=''):
        if comment:
            self.comment = comment
        nplies = len(mids)
        if failure_theories is None:
            failure_theories = [None] * nplies
        if interlaminar_failure_theories is None:
            interlaminar_failure_theories = [None] * nplies
        if souts is None:
            souts = ['NO'] * nplies
        SolidProperty.__init__(self)
        self.pid = pid
        self.cordm = cordm
        self.psdir = psdir
        self.sb = sb
        self.nb = nb
        self.tref = tref
        self.ge = ge
        assert psdir in [12, 13, 21, 23, 31, 32], psdir
        self.global_ply_ids = global_ply_ids
        self.mids = mids
        self.thicknesses = np.array(thicknesses, dtype='float64')
        self.thetas = np.array(thetas, dtype='float64')
        self.failure_theories = failure_theories
        self.interlaminar_failure_theories = interlaminar_failure_theories
        self.souts = souts
        self.mids_ref = None

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a PCOMPS card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        pid = integer(card, 1, 'pid')
        cordm = integer_or_blank(card, 2, 'cordm', 0)
        psdir = integer_or_blank(card, 3, 'psdir', 13)
        sb = double_or_blank(card, 4, 'sb')
        nb = double_or_blank(card, 5, 'nb')
        tref = double_or_blank(card, 6, 'tref', 0.0)
        ge = double_or_blank(card, 7, 'ge', 0.0)
        nfields = len(card) - 1
        #nrows =
        ifield = 9
        global_ply_ids = []
        mids = []
        thicknesses = []
        thetas = []
        failure_theories = []
        interlaminar_failure_theories = []
        souts = []
        iply = 1
        while ifield < nfields:
            global_ply_id = integer(card, ifield, 'global_ply_id_%i' % iply)
            mid = integer(card, ifield + 1, 'mid_%i' % iply)
            t = double(card, ifield + 2, 'thickness_%i' % iply)
            theta = double(card, ifield + 3, 'theta_%i' % iply)
            ft = string_or_blank(card, ifield + 4, 'failure_theory_%i' % iply)
            ift = string_or_blank(card, ifield + 5, 'interlaminar_failure_theory_%i' % iply)
            sout = string_or_blank(card, ifield + 6, 'sout_%i' % iply, 'NO')
            global_ply_ids.append(global_ply_id)
            mids.append(mid)
            thicknesses.append(t)
            thetas.append(theta)
            failure_theories.append(ft)
            interlaminar_failure_theories.append(ift)
            souts.append(sout)
            iply += 1
            ifield += 8
        assert len(card) <= ifield, 'len(PCOMPS card) = %i\ncard=%s' % (len(card), card)
        return PCOMPS(pid, global_ply_ids, mids, thicknesses, thetas,
                      cordm, psdir, sb, nb, tref, ge,
                      failure_theories, interlaminar_failure_theories, souts,
                      comment=comment)

    def _verify(self, xref):
        """
        Verifies all methods for this object work

        Parameters
        ----------
        xref : bool
            has this model been cross referenced
        """
        pass

    def Mid(self):
        return self.mids[0]

    def Rho(self):
        """
        Returns the density
        """
        rho = np.zeros(self.thicknesses.size, dtype='float64')
        for i, mid in enumerate(self.mids_ref):
            rho[i] = mid.rho
        rhot = rho * self.thicknesses / self.thicknesses.sum()
        #print('rhot =', rhot)
        return rhot.mean()

    def Mids(self):
        return self.mids

    def cross_reference(self, model):
        msg = ', which is required by PSOLID pid=%s' % self.pid
        self.mids_ref = []
        for mid in self.mids:
            mid_ref = model.Material(mid, msg=msg)
            self.mids_ref.append(mid_ref)

    def uncross_reference(self):
        self.mids_ref = None

    def raw_fields(self):
        fields = ['PCOMPS', self.pid, self.cordm, self.psdir, self.sb,
                  self.nb, self.tref, self.ge, None]
        for glply, mid, t, theta, ft, ift, sout in zip(self.global_ply_ids,
                                                       self.mids, self.thicknesses, self.thetas,
                                                       self.failure_theories,
                                                       self.interlaminar_failure_theories,
                                                       self.souts):
            fields += [glply, mid, t, theta, ft, ift, sout, None]
        return fields

    def repr_fields(self):
        #cordm = set_blank_if_default(self.cordm, 0)
        #fctn = set_blank_if_default(self.fctn, 'SMECH')
        fields = ['PCOMPS', self.pid, self.cordm, self.psdir, self.sb,
                  self.nb, self.tref, self.ge, None]
        for glply, mid, t, theta, ft, ift, sout in zip(self.global_ply_ids,
                                                       self.mids, self.thicknesses, self.thetas,
                                                       self.failure_theories,
                                                       self.interlaminar_failure_theories,
                                                       self.souts):
            fields += [glply, mid, t, theta, ft, ift, sout, None]
        return fields

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        # this card has integers & strings, so it uses...
        return self.comment + print_card_8(card)

class PSOLID(SolidProperty):
    """
    +--------+-----+-----+-------+-----+--------+---------+------+
    |    1   |  2  |  3  |   4   |  5  |    6   |    7    |   8  |
    +========+=====+=====+=======+=====+========+=========+======+
    | PSOLID | PID | MID | CORDM | IN  | STRESS |   ISOP  | FCTN |
    +--------+-----+-----+-------+-----+--------+---------+------+
    | PSOLID |  1  |     |   1   | 0   |        |         |      |
    +--------+-----+-----+-------+-----+--------+---------+------+
    | PSOLID |  2  | 100 |   6   | TWO |  GRID  | REDUCED |      |
    +--------+-----+-----+-------+-----+--------+---------+------+
    """
    type = 'PSOLID'
    _field_map = {
        1: 'pid', 2:'mid', 3:'cordm', 4:'integ', 5:'stress',
        6:'isop', 7:'fctn',
    }

    def __init__(self, pid, mid, cordm=0, integ=None, stress=None, isop=None,
                 fctn='SMECH', comment=''):
        """
        Creates a PSOLID card

        Parameters
        ----------
        pid : int
            property id
        mid : int
            material id
        cordm : int; default=0
            material coordinate system
        integ : int; default=None
            None-varies depending on element type
            0, 'BUBBLE'
            1, 'GAUSS'
            2, 'TWO'
            3, 'THREE'
            REDUCED
            FULL
        stress : int/str; default=None
            None/GRID, 1-GAUSS
        isop : int/str; default=None
            0-REDUCED
            1-FULL
        fctn : str; default='SMECH'
            PFLUID/SMECH
        comment : str; default=''
            a comment for the card
        """
        SolidProperty.__init__(self)
        if comment:
            self.comment = comment
        #: Property ID
        self.pid = pid
        #: Material ID
        self.mid = mid
        self.cordm = cordm
        #valid_integration = ['THREE', 'TWO', 'FULL', 'BUBBLE',
        #                     2, 3, None, 'REDUCED']

        # note that None is supposed to vary depending on element type
        # 0-BUBBLE (not for midside nodes)
        # 1-GAUSS
        # 2-TWO
        # 3-THREE
        self.integ = integ

        # blank/GRID
        # 1-GAUSS
        self.stress = stress

        # note that None is supposed to vary depending on element type
        # 0-REDUCED
        # 1-FULL
        self.isop = isop

        # PFLUID
        # SMECH
        self.fctn = fctn
        self.mid_ref = None

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a PSOLID card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        pid = integer(card, 1, 'pid')
        mid = integer(card, 2, 'mid')
        cordm = integer_or_blank(card, 3, 'cordm', 0)
        integ = integer_string_or_blank(card, 4, 'integ')
        stress = integer_string_or_blank(card, 5, 'stress')
        isop = integer_string_or_blank(card, 6, 'isop')
        fctn = string_or_blank(card, 7, 'fctn', 'SMECH')
        assert len(card) <= 8, 'len(PSOLID card) = %i\ncard=%s' % (len(card), card)
        return cls(pid, mid, cordm, integ, stress, isop,
                   fctn, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        """
        Adds a PSOLID card from the OP2

        Parameters
        ----------
        data : List[varies]
            a list of fields defined in OP2 format
        comment : str; default=''
            a comment for the card
        """
        pid = data[0]
        mid = data[1]
        cordm = data[2]
        integ = data[3]
        stress = data[4]
        isop = data[5]
        fctn = data[6].decode('latin1')

        if fctn == 'SMEC':
            fctn = 'SMECH'
        elif fctn == 'PFLU':
            fctn = 'PFLUID'
        else:
            raise NotImplementedError('PSOLID; fctn=%r' % fctn)
        return PSOLID(pid, mid, cordm, integ, stress, isop,
                      fctn, comment=comment)

    def cross_reference(self, model):
        # type: (Any) -> None
        """cross reference method for a PSOLID"""
        msg = ', which is required by PSOLID pid=%s' % (self.pid)
        self.mid_ref = model.Material(self.mid, msg)

    def uncross_reference(self):
        # type: () -> None
        self.mid = self.Mid()
        self.mid_ref = None

    def E(self):
        return self.mid_ref.E()

    def G(self):
        return self.mid_ref.G()

    def Nu(self):
        return self.mid_ref.Nu()

    def materials(self):
        return [self.mid_ref]

    def Rho(self):
        """
        Returns the density
        """
        return self.mid_ref.rho

    def _verify(self, xref):
        pid = self.Pid()
        mid = self.Mid()
        assert isinstance(pid, integer_types), 'pid=%r' % pid
        assert isinstance(mid, integer_types), 'mid=%r' % mid

        if xref:
            if self.mid_ref.type not in ['MAT1', 'MAT4', 'MAT5', 'MAT9', 'MAT10', 'MAT11']:
                msg = 'mid=%i self.mid_ref.type=%s' % (mid, self.mid_ref.type)
                raise TypeError(msg)

    def raw_fields(self):
        fields = ['PSOLID', self.pid, self.Mid(), self.cordm, self.integ,
                  self.stress, self.isop, self.fctn]
        return fields

    def repr_fields(self):
        cordm = set_blank_if_default(self.cordm, 0)
        fctn = set_blank_if_default(self.fctn, 'SMECH')
        fields = ['PSOLID', self.pid, self.Mid(), cordm, self.integ,
                  self.stress, self.isop, fctn]
        return fields

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        # this card has integers & strings, so it uses...
        return self.comment + print_card_8(card)


class PIHEX(PSOLID):
    type = 'PIHEX'
    def __init__(self, pid, mid, cordm=0, integ=None, stress=None, isop=None,
                 fctn='SMECH', comment=''):
        PSOLID.__init__(self, pid, mid, cordm, integ, stress, isop,
                        fctn, comment=comment)

    def raw_fields(self):
        fields = ['PIHEX', self.pid, self.Mid(), self.cordm, self.integ,
                  self.stress, self.isop, self.fctn]
        return fields

    def repr_fields(self):
        cordm = set_blank_if_default(self.cordm, 0)
        fctn = set_blank_if_default(self.fctn, 'SMECH')
        fields = ['PIHEX', self.pid, self.Mid(), cordm, self.integ,
                  self.stress, self.isop, fctn]
        return fields
