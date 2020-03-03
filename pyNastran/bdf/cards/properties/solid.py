#pylint: disable=C0103,C0111,R0902
"""
All solid properties are defined in this file.  This includes:
 * PCOMPS
 * PLSOLID
 * PSOLID

 All solid properties are Property objects.

"""
from __future__ import annotations
from typing import Dict, Any, TYPE_CHECKING
import numpy as np

from pyNastran.utils.numpy_utils import integer_types
from pyNastran.bdf.field_writer_8 import set_blank_if_default
from pyNastran.bdf.cards.base_card import Property
from pyNastran.bdf.bdf_interface.assign_type import (
    integer, integer_or_blank, double, double_or_blank, string_or_blank,
    integer_string_or_blank)
from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.field_writer_16 import print_card_16
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf import BDF


class PLSOLID(Property):
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

    @classmethod
    def _init_from_empty(cls):
        pid = 1
        mid = 1
        return PLSOLID(pid, mid, stress_strain='GRID', ge=0., comment='')

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
        Property.__init__(self)
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

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ', which is required by PLSOLID pid=%s' % self.pid
        self.mid_ref = model.HyperelasticMaterial(self.mid, msg) # MATHP, MATHE

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.mid = self.Mid()
        self.mid_ref = None

    def Rho(self):
        """Returns the density"""
        return self.mid_ref.rho

    def _verify(self, xref):
        pass

    def raw_fields(self):
        stress_strain = set_blank_if_default(self.stress_strain, 'GRID')
        fields = ['PLSOLID', self.pid, self.Mid(), stress_strain]
        return fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)


class PCOMPS(Property):
    type = 'PCOMPS'
    _field_map = {
        #1: 'pid', 2:'mid', 3:'cordm', 4:'integ', 5:'stress',
        #6:'isop', 7:'fctn',
    }  # type: Dict[int,str]

    @classmethod
    def _init_from_empty(cls):
        pid = 1
        global_ply_ids = [1, 2]
        mids = [10, 20]
        thicknesses = [0.1, 0.2]
        thetas = [30., 60.]
        return PCOMPS(pid, global_ply_ids, mids, thicknesses, thetas,
                      cordm=0, psdir=13, sb=None, nb=None, tref=0.0,
                      ge=0.0, failure_theories=None,
                      interlaminar_failure_theories=None, souts=None, comment='')

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
        Property.__init__(self)
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

    def cross_reference(self, model: BDF) -> None:
        msg = ', which is required by PSOLID pid=%s' % self.pid
        self.mids_ref = []
        for mid in self.mids:
            mid_ref = model.Material(mid, msg=msg)
            self.mids_ref.append(mid_ref)

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.mids = self.material_ids
        self.mids_ref = None

    @property
    def material_ids(self):
        if self.mids_ref is None:
            return self.mids
        mids = []
        for mid_ref in self.mids_ref:
            mids.append(mid_ref.mid)
        return mids

    def raw_fields(self):
        fields = ['PCOMPS', self.pid, self.cordm, self.psdir, self.sb,
                  self.nb, self.tref, self.ge, None]
        mids = self.material_ids
        for glply, mid, t, theta, ft, ift, sout in zip(self.global_ply_ids,
                                                       mids, self.thicknesses, self.thetas,
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
        mids = self.material_ids
        for glply, mid, t, theta, ft, ift, sout in zip(self.global_ply_ids,
                                                       mids, self.thicknesses, self.thetas,
                                                       self.failure_theories,
                                                       self.interlaminar_failure_theories,
                                                       self.souts):
            fields += [glply, mid, t, theta, ft, ift, sout, None]
        return fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        # this card has integers & strings, so it uses...
        return self.comment + print_card_8(card)


class PSOLID(Property):
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

    @classmethod
    def _init_from_empty(cls):
        pid = 1
        mid = 1
        return PSOLID(pid, mid, cordm=0, integ=None, stress=None,
                      isop=None, fctn='SMECH', comment='')

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
        Property.__init__(self)
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
        if integ == 0:
            integ = 'BUBBLE'
        elif integ == 1:
            integ = 'GAUSS'
        elif integ == 2:
            integ = 'TWO'
        elif integ == 3:
            integ = 'THREE'

        self.integ = integ

        # blank/GRID
        # 1-GAUSS
        if stress == 0:
            stress = 'GRID'
        elif stress == 1:
            stress = 'GAUSS'
        self.stress = stress

        # note that None is supposed to vary depending on element type
        # 0-REDUCED
        # 1-FULL
        if isop == 0:
            isop = 'REDUCED'
        elif isop == 1:
            isop = 'FULL'
        elif isop == 2:
            pass
        self.isop = isop

        # PFLUID
        # SMECH
        if fctn == 'SMEC':
            fctn = 'SMECH'
        elif fctn == 'PFLU':
            fctn = 'PFLUID'
        self.fctn = fctn
        self.mid_ref = None

    @classmethod
    def export_to_hdf5(cls, h5_file, model, pids):
        """exports the properties in a vectorized way"""
        encoding = model._encoding
        #comments = []
        mid = []
        cordm = []
        integ = []
        stress = []
        isop = []
        fctn = []
        assert len(pids) > 0, pids
        for pid in pids:
            prop = model.properties[pid]
            #comments.append(prop.comment)
            mid.append(prop.mid)
            cordm.append(prop.cordm)
            if prop.integ is None:
                integ.append(b'')
            elif prop.integ == 0:
                integ.append(b'BUBBLE')
            elif prop.integ == 1:
                integ.append(b'GAUSS')
            elif prop.integ == 2:
                integ.append(b'TWO')
            elif prop.integ == 3:
                integ.append(b'THREE')
            else:
                #0, 'BUBBLE'
                #1, 'GAUSS'
                #2, 'TWO'
                #3, 'THREE'
                #REDUCED
                #FULL
                integ.append(prop.integ.encode(encoding))

            if prop.stress is None:
                stress.append(b'')
            elif prop.stress == 1:
                stress.append(b'GAUSS')
            else:
                # None/GRID, 1-GAUSS
                stress.append(prop.stress.encode(encoding))

            if prop.isop is None:
                isop.append(b'')
            elif prop.isop == 0:
                isop.append(b'REDUCED')
            elif prop.isop == 1:
                isop.append(b'FULL')
            elif prop.isop == 2:
                isop.append('2')
            else:
                isop.append(prop.isop.encode(encoding))

            # 'SMECH'
            fctn.append(prop.fctn.encode(encoding))
        #fctn = np.array(fctn, dtype='<8U')
        #h5_file.create_dataset('_comment', data=comments)
        h5_file.create_dataset('pid', data=pids)
        h5_file.create_dataset('mid', data=mid)
        h5_file.create_dataset('cordm', data=cordm)
        h5_file.create_dataset('integ', data=integ)
        h5_file.create_dataset('stress', data=stress)
        h5_file.create_dataset('isop', data=isop)
        h5_file.create_dataset('fctn', data=fctn)

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

        #integ : int; default=None
        #None-varies depending on element type
        #0, 'BUBBLE'
        #1, 'GAUSS'
        #2, 'TWO'
        #3, 'THREE'
        #REDUCED
        #FULL

        # stress : int, string, or blank
        #    blank/GRID
        #    1-GAUSS
        if stress == 0:
            stress = 'GRID'
        elif stress == 1:
            stress = 'GAUSS'
        else:  # pragma: no cover
            raise NotImplementedError('stress=%s and must be [0, 1]' % stress)

        if isop == 0:
            isop = 'REDUCED'
        elif isop == 1:
            isop = 'FULL'
        elif isop == 2:
            pass
        else:  # pragma: no cover
            raise NotImplementedError('isop=%s and must be [0, 1, 2]' % isop)

        fctn = fctn.rstrip()
        if fctn == 'SMEC':
            fctn = 'SMECH'
        elif fctn == 'PFLU':
            fctn = 'PFLUID'
        else:  # pragma: no cover
            raise NotImplementedError('PSOLID; fctn=%r' % fctn)
        return PSOLID(pid, mid, cordm, integ, stress, isop,
                      fctn, comment=comment)

    def cross_reference(self, model: BDF) -> None:
        """cross reference method for a PSOLID"""
        msg = ', which is required by PSOLID pid=%s' % (self.pid)
        self.mid_ref = model.Material(self.mid, msg)

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
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

    def write_card(self, size: int=8, is_double: bool=False) -> str:
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
