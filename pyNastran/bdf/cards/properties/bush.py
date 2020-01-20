"""
All bush properties are defined in this file.  This includes:

 *   PBUSH
 *   PBUSH1D
 *   PBUSH2D (not implemented)
 *   PBUSHT

All bush properties are BushingProperty and Property objects.

"""
from __future__ import annotations
from typing import TYPE_CHECKING
from pyNastran.utils.numpy_utils import integer_types
from pyNastran.bdf.cards.base_card import Property
from pyNastran.bdf.bdf_interface.assign_type import (
    integer, integer_or_blank, double, double_or_blank, string,
    string_or_blank, blank, fields)
from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.field_writer_16 import print_card_16
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf import BDF


class BushingProperty(Property):
    type = 'BushingProperty'

    def __init__(self):
        Property.__init__(self)

    def cross_reference(self, model: BDF) -> None:
        pass

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        pass


class PBUSH(BushingProperty):
    """
    Generalized Spring-and-Damper Property
    Defines the nominal property values for a generalized spring-and-damper
    structural element.

    +-------+-----+------+-----+-----+-----+-----+-----+----+
    |   1   |  2  |   3  |  4  |  5  |  6  |  7  |  8  | 9  |
    +=======+=====+======+=====+=====+=====+=====+=====+====+
    | PBUSH | PID |   K  |  K1 |  K2 |  K3 |  K4 |  K5 | K6 |
    +-------+-----+------+-----+-----+-----+-----+-----+----+
    |       |  B  |  B1  |  B2 |  B3 |  B4 |  B5 |  B6 |    |
    +-------+-----+------+-----+-----+-----+-----+-----+----+
    |       | GE  |  GE1 | GE2 | GE3 | GE4 | GE5 | GE6 |    |
    +-------+-----+------+-----+-----+-----+-----+-----+----+
    |       | RCV |  SA  |  ST |  EA |  ET |     |     |    |
    +-------+-----+------+-----+-----+-----+-----+-----+----+
    |       |  M  | MASS |     |     |     |     |     |    |
    +-------+-----+------+-----+-----+-----+-----+-----+----+
    """
    type = 'PBUSH'
    _field_map = {
        1: 'pid',
    }
    pname_map = {
        -2 : 'K1', -3 : 'K2', -4 : 'K3', -5 : 'K4', -6 : 'K5', -7 : 'K6',
        -8 : 'B1', -9 : 'B2', -10 : 'B3', -11 : 'B4', -12 : 'B5', -13 : 'B6',
        -14 : 'GE1', -15 : 'GE2', -16 : 'GE3', -17 : 'GE4', -18 : 'GE5', -19 : 'GE6',
        -20 : 'SA', -21 : 'ST', -22 : 'EA', -23 : 'ET',
    }
    def update_by_pname_fid(self, name, value):
        if name == 'B1':
            self.Bi[0] = value
        elif name == 'B2':
            self.Bi[1] = value
        elif name == 'B3':
            self.Bi[2] = value
        elif name == 'B4':
            self.Bi[3] = value
        elif name == 'B5':
            self.Bi[4] = value
        elif name == 'B6':
            self.Bi[5] = value

        elif name == 'K1':
            self.Ki[0] = value
        elif name == 'K2':
            self.Ki[1] = value
        elif name == 'K3':
            self.Ki[2] = value
        elif name == 'K4':
            self.Ki[3] = value
        elif name == 'K5':
            self.Ki[4] = value
        elif name == 'K6':
            self.Ki[5] = value

        elif name == 'GE1':
            self.GEi[0] = value
        elif name == 'GE2':
            self.GEi[1] = value
        elif name == 'GE3':
            self.GEi[2] = value
        elif name == 'GE4':
            self.GEi[3] = value
        elif name == 'GE5':
            self.GEi[4] = value
        elif name == 'GE6':
            self.GEi[5] = value
        #elif name == 'M':
            #self.mass
        elif isinstance(name, int) and name in self.pname_map:
            name2 = self.pname_map[name]
            self.update_by_pname_fid(name2, value)
        else:
            raise NotImplementedError('property_type=%r has not implemented %r in pname_map' % (
                self.type, name))
    #pname_fid_map = {
        #4 : 'A', 'A' : 'A',
        #5 : 'i1', 'I1' : 'i1',
        #6 : 'i2', 'I2' : 'i2',
        #7 : 'i12', 'I12' : 'i12',
        #5 : 'j', 'J' : 'j',
    #}

    def __init__(self, pid, k, b, ge, rcv=None, mass=None, comment=''):
        """
        Creates a PBUSH card, which defines a property for a CBUSH

        Parameters
        ----------
        pid : int
            property id
        k : List[float]
            Nominal stiffness values in directions 1 through 6.
            len(k) = 6
        b : List[float]
            Nominal damping coefficients in direction 1 through 6 in units of
            force per unit velocity
            len(b) = 6
        ge : List[float]
            Nominal structural damping constant in directions 1 through 6.
            len(ge) = 6
        rcv : List[float]; default=None -> (None, None, None, None)
            [sa, st, ea, et] = rcv
            length(mass_fields) = 4
        mass : float; default=None
            lumped mass of the CBUSH
            This is an MSC only parameter.
        comment : str; default=''
            a comment for the card

        """
        BushingProperty.__init__(self)
        if comment:
            self.comment = comment

        #: Property ID
        self.pid = pid
        self.vars = []

        # K parameter
        self.Ki = k
        if k:
            nk = len(k)
            if nk < 6:
                k.extend([0.] * (6 - nk))
            self.vars.append('K')
        # B parameter
        self.Bi = b
        if b:
            nb = len(b)
            if nb < 6:
                b.extend([0.] * (6 - nb))
            self.vars.append('B')

        # GE parameter
        self.GEi = ge
        if ge:
            nge = len(ge)
            if nge < 6:
                ge.extend([0.] * (6 - nge))
            self.vars.append('GE')

        # RCV parameters
        if rcv is None:
            sa, st, ea, et = None, None, None, None
        else:
            sa, st, ea, et = rcv
        if sa is not None or st is not None or ea is not None or et is not None:
            self.vars.append('RCV')
        self.sa = sa
        self.st = st
        self.ea = ea
        self.et = et

        # M parameter (MSC only; in 2016, not in 2005)
        self.mass = mass
        if mass:
            self.vars.append('M')

    @classmethod
    def _init_from_empty(cls):
        pid = 1
        k = [1.]
        b = [1.]
        ge = [1.]
        return PBUSH(pid, k, b, ge, rcv=None, mass=None, comment='')

    def validate(self):
        assert isinstance(self.Ki, list), 'PBUSH: pid=%i type(Ki)=%s Ki=%s' % (self.pid, type(self.Ki), self.Ki)
        assert isinstance(self.Bi, list), 'PBUSH: pid=%i type(Bi)=%s Bi=%s' % (self.pid, type(self.Bi), self.Bi)
        assert isinstance(self.GEi, list), 'PBUSH: pid=%i type(GEi)=%s GEi=%s' % (self.pid, type(self.GEi), self.GEi)

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a PBUSH card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        k_fields = []
        b_fields = []
        ge_fields = []
        rcv_fields = [None, None, None, None]
        mass = None

        pid = integer(card, 1, 'pid')

        nfields = card.nfields
        istart = 2
        while istart < nfields:
            pname = string(card, istart, 'pname')

            if pname == 'K':
                # Flag indicating that the next 1 to 6 fields are stiffness values in
                # the element coordinate system.
                #self.k = string(card, istart, 'k')

                #: Nominal stiffness values in directions 1 through 6.
                #: See Remarks 2 and 3. (Real; Default = 0.0)
                k_fields = cls._read_var(card, 'Ki', istart + 1, istart + 7)

            elif pname == 'B':
                # Flag indicating that the next 1 to 6 fields are force-per-velocity
                # damping.
                #self.b = string(card, istart, 'b')

                #: Force per unit velocity (Real)
                #: Nominal damping coefficients in direction 1 through 6 in units of
                #: force per unit velocity. See Remarks 2, 3, and 9. (Real; Default=0.)
                b_fields = cls._read_var(card, 'Bi', istart + 1, istart + 7)

            elif pname == 'GE':
                # Flag indicating that the next fields, 1 through 6 are structural
                # damping constants. See Remark 7. (Character)
                #self.ge = string(card, istart, 'ge')

                #: Nominal structural damping constant in directions 1 through 6. See
                #: Remarks 2. and 3. (Real; Default = 0.0)
                ge_fields = cls._read_var(card, 'GEi', istart + 1, istart + 7)
            elif pname == 'RCV':
                rcv_fields = cls._read_rcv(card, istart)
            elif pname == 'M':
                # Lumped mass of the cbush; default=0.0
                mass = double_or_blank(card, istart + 1, 'mass', 0.)
            else:
                raise RuntimeError('unsupported PBUSH type; pname=%r' % pname)
                #break #  old version...
            istart += 8
        return PBUSH(pid, k_fields, b_fields, ge_fields, rcv_fields, mass,
                     comment=comment)

    @classmethod
    def _read_var(cls, card, var_prefix, istart, iend):
        Ki = fields(double_or_blank, card, var_prefix, istart, iend)
        return Ki

    @classmethod
    def add_op2_data(cls, data, comment=''):
        """
        Adds a PBUSH card from the OP2

        Parameters
        ----------
        data : List[varies]
            a list of fields defined in OP2 format
        comment : str; default=''
            a comment for the card

        """
        (pid, k1, k2, k3, k4, k5, k6, b1, b2, b3, b4, b5, b6,
         g1, g2, g3, g4, g5, g6, sa, st, ea, et) = data
        k_fields = [k1, k2, k3, k4, k5, k6]
        b_fields = [b1, b2, b3, b4, b5, b6]
        ge_fields = [g1, g2, g3, g4, g5, g6]
        rcv_fields = [sa, st, ea, et]
        mass = 0.
        return PBUSH(pid, k_fields, b_fields, ge_fields, rcv_fields, mass,
                     comment=comment)

    def _verify(self, xref):
        pid = self.Pid()
        assert isinstance(pid, integer_types), 'pid=%r' % pid

    @classmethod
    def _read_rcv(cls, card, istart):
        # Flag indicating that the next 1 to 4 fields are stress or strain
        # coefficients. (Character)
        #self.rcv = string(card, istart, 'rcv')
        sa = double_or_blank(card, istart + 1, 'sa', 1.)
        st = double_or_blank(card, istart + 2, 'st', 1.)
        ea = double_or_blank(card, istart + 3, 'ea', 1.)
        et = double_or_blank(card, istart + 4, 'et', 1.)
        return [sa, st, ea, et]

    def raw_fields(self):
        list_fields = ['PBUSH', self.pid]
        for var in self.vars:
            if var == 'K':
                list_fields += ['K'] + self.Ki
            elif var == 'B':
                list_fields += ['B'] + self.Bi
            elif var == 'GE':
                list_fields += ['GE'] + self.GEi
            elif var == 'RCV':
                list_fields += ['RCV', self.sa, self.st, self.ea, self.et]
            elif var == 'M':
                list_fields += ['M', self.mass]
            else:
                raise RuntimeError('not supported PBUSH field...')
            nspaces = 8 - (len(list_fields) - 1) % 8

            if nspaces == 8:
                list_fields += [None]
            elif nspaces < 8:
                list_fields += [None] * (nspaces + 1)
        return list_fields

    def repr_fields(self):
        return self.raw_fields()

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)


class PBUSH1D(BushingProperty):
    """
    +---------+--------+-------+--------+--------+-------+-------+-------+
    |   1     |    2   |   3   |    4   |    5   |   6   |   7   |   8   |
    +=========+========+=======+========+========+=======+=======+=======+
    | PBUSH1D |   PID  |   K   |    C   |    M   |       |   SA  |   SE  |
    +---------+--------+-------+--------+--------+-------+-------+-------+
    |         | SHOCKA | TYPE  |   CVT  |   CVC  | EXPVT | EXPVC |  IDTS |
    +---------+--------+-------+--------+--------+-------+-------+-------+
    |         | IDETS  | IDECS | IDETSD | IDECSD |       |       |       |
    +---------+--------+-------+--------+--------+-------+-------+-------+
    |         | SPRING |  TYPE |   IDT  |   IDC  | IDTDU | IDCDU |       |
    +---------+--------+-------+--------+--------+-------+-------+-------+
    |         | DAMPER |  TYPE |   IDT  |   IDC  | IDTDV | IDCDV |       |
    +---------+--------+-------+--------+--------+-------+-------+-------+
    |         | GENER  |  IDT  |   IDC  |  IDTDU | IDCDU | IDTDV | IDCDV |
    +---------+--------+-------+--------+--------+-------+-------+-------+
    """
    type = 'PBUSH1D'
    _properties = ['pname_fid_map']
    pname_fid_map = {
        'K' : 'k',
        'C' : 'c',
        'M' : 'm',
    }
    @classmethod
    def _init_from_empty(cls):
        pid = 1
        return PBUSH1D(pid, k=0., c=0., m=0., sa=0., se=0., optional_vars=None, comment='')

    def __init__(self, pid: int,
                 k: float=0., c: float=0., m: float=0.,
                 sa: float=0., se: float=0., optional_vars=None,
                 comment: str=''):
        """
        Creates a PBUSH1D card

        Parameters
        ----------
        pid : int
            property id
        k : float
           stiffness
        c : float
            Viscous damping
        m : float
            mass
        sa : float
            Stress recovery coefficient [1/area].
        se : float
            Strain recovery coefficient [1/length].
        optional_vars : dict[name] : value; default=None
            name : str
                SHOCKA, SPRING, DAMPER, GENER
            values : List[varies]
                the values
            SHOCKA:
                Coefficients of the following force versus
                velocity/displacement relationship
                F(u, v) = Cv * S(u) * sign(v) * |v|^EXPV
                TYPE CVT CVC EXPVT EXPVC
                IDTS IDETS IDECS IDETSD IDECSD
                CVT/CVC : int
                    Viscous damping coefficient CV for tension v > 0, force
                    per unit velocity.
                EXPVT/EXPVC : int
                    Exponent of velocity EXPV for tension v > 0 (or compression v < 0).
                IDTS : int
                    Identification number of a TABLEDi entry for tension and
                    compression if TYPE=TABLE. The TABLEDi entry
                    defines the scale factor S, versus displacement u.
                IDETS/IDECS : int
                    Identification number of a DEQATN entry for tension if
                    TYPE=EQUAT. The DEQATN entry defines the scale
                    factor S, versus displacement u, for tension u > 0
                    (or compression v < 0).
                IDETSD/IDECSD : int
                    Identification number of a DEQATN entry for tension if
                    TYPE=EQUAT. The DEQATN entry defines the defines the scale
                    factor S, versus displacement u, for tension u > 0
                    (or compression v < 0).
            SPRING:
                Nonlinear elastic spring element in terms of a force versus
                displacement relationship
                TYPE IDT IDC IDTDU IDCDU
            DAMPER:
                Nonlinear viscous element in terms of a force versus
                velocity relationship.
                TYPE IDT IDC IDTDV IDCDV
            GENER:
                General nonlinear elastic spring and viscous damper
                element in terms of a force versus displacement and
                velocity relationship.  For this element, the relationship
                can only be defined with TYPE=EQUAT (and it's implicit).
                IDT IDC IDTDU IDCDU IDTDV IDCDV

            TYPE : int
               the type of the result; {TABLE, EQUAT}
            IDT/IDC : int
                tension/compression table/equation
            IDTDU/IDCDU : int
                du/dt tension/compression table/eq
            IDTDV/IDCDV : int
                dv/dt tension/compression table/eq

        """
        BushingProperty.__init__(self)
        if comment:
            self.comment = comment

        #: Property ID
        self.pid = pid
        self.k = k
        self.c = c
        self.m = m

        self.sa = sa
        self.se = se

        # SPRING parameters
        self.spring_type = None
        self.spring_idt = None
        self.spring_idc = None
        self.spring_idtdu = None
        self.spring_idcdu = None

        # DAMPER parameters
        self.damper_type = None
        self.damper_idt = None
        self.damper_idc = None
        self.damper_idtdv = None
        self.damper_idcdv = None

        # GENER parameters
        #self.gener_idt = None
        #self.gener_idc = None
        #self.gener_idtdu = None
        #self.gener_idcdu = None
        #self.gener_idtdv = None
        #self.gener_idcdv = None

        # SHOCK parameters
        self.shock_type = None
        self.shock_cvt = None
        self.shock_cvc = None
        self.shock_exp_vt = None
        self.shock_exp_vc = None
        self.shock_idts = None

        self.shock_idets = None
        self.shock_idecs = None
        self.shock_idetsd = None
        self.shock_idecsd = None
        if optional_vars:
            for key, values in optional_vars.items():
                if key == 'SHOCKA':
                    (shock_type, shock_cvt, shock_cvc, shock_exp_vt, shock_exp_vc,
                     shock_idts, shock_idets, shock_idecs, shock_idetsd, shock_idecsd
                    ) = values
                    self.shock_type = shock_type
                    self.shock_cvt = shock_cvt
                    self.shock_cvc = shock_cvc
                    self.shock_exp_vt = shock_exp_vt
                    self.shock_exp_vc = shock_exp_vc
                    self.shock_idts = shock_idts
                    self.shock_idets = shock_idets
                    self.shock_idecs = shock_idecs
                    self.shock_idetsd = shock_idetsd
                    self.shock_idecsd = shock_idecsd
                    assert isinstance(self.shock_type, str), f'shock_type={shock_type}'

                elif key == 'SPRING':
                    (spring_type, spring_idt, spring_idc, spring_idtdu,
                     spring_idcdu) = values
                    self.spring_type = spring_type
                    self.spring_idt = spring_idt
                    self.spring_idtc = spring_idc
                    self.spring_idc = spring_idc
                    self.spring_idtdu = spring_idtdu
                    self.spring_idcdu = spring_idcdu
                    assert isinstance(self.spring_type, str), f'spring_type={spring_type}'

                elif key == 'DAMPER':
                    (damper_type, damper_idt, damper_idc, damper_idtdv,
                     damper_idcdv) = values
                    self.damper_type = damper_type
                    self.damper_idt = damper_idt
                    self.damper_idc = damper_idc
                    self.damper_idtdv = damper_idtdv
                    self.damper_idcdv = damper_idcdv
                    assert isinstance(self.damper_type, str), f'damper_type={damper_type}'

                elif key == 'GENER':
                    (
                        gener_idt, gener_idc,
                        gener_idtdu, gener_idcdu,
                        gener_idtdv, gener_idcdv
                    ) = values

                    self.gener_idt = gener_idt
                    self.gener_idc = gener_idc
                    self.gener_idtdu = gener_idtdu
                    self.gener_idcdu = gener_idcdu
                    self.gener_idtdv = gener_idtdv
                    self.gener_idcdv = gener_idcdv
                else:
                    msg = ('PBUSH1D: pid=%s key=%r and must be '
                           '{SHOCKA, SPRING, DAPMER, GENER}' % self.pid)
                    raise RuntimeError(msg)
        if optional_vars is None:
            self.vars = []
        else:
            self.vars = list(optional_vars.keys())
            self.vars.sort()

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a PBUSH1D card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        pid = integer(card, 1, 'pid')
        k = double_or_blank(card, 2, 'k', 0.0)
        c = double_or_blank(card, 3, 'c', 0.0)
        m = double_or_blank(card, 4, 'm', 0.0)

        sa = double_or_blank(card, 6, 'sa', 0.0)
        se = double_or_blank(card, 7, 'se', 0.0)

        nfields = card.nfields
        optional_vars = {}
        istart = 9
        while istart < nfields:
            pname = string(card, istart, 'pname')
            if pname == 'SHOCKA':
                istart, out = cls._read_shock(card, istart)
                optional_vars['SHOCKA'] = out

            elif pname == 'SPRING':
                out = cls._read_spring(card, istart)
                optional_vars['SPRING'] = out

            elif pname == 'DAMPER':
                out = cls._read_damper(card, istart)
                optional_vars['DAMPER'] = out

            elif pname == 'GENER':
                out = cls._read_gener(card, istart)
                optional_vars['GENER'] = out
            else:
                break
            istart += 8

        #self.pid = pid
        #self.k = k
        #self.c = c
        #self.m = m

        #self.sa = sa
        #self.se = se
        return PBUSH1D(pid, k=k, c=c, m=m, sa=sa, se=se,
                       optional_vars=optional_vars, comment=comment)

    #@classmethod
    #def add_op2_data(cls, data, comment=''):
        #pid = data[0]
        #b = data[1]
        #raise NotImplementedError('PBUSH1D data...')

    def _verify(self, xref):
        pid = self.Pid()
        assert isinstance(pid, int), 'pid=%r' % pid

    @staticmethod
    def _read_shock(card, istart):
        """
        F(u, v) = Cv * S(u) * sign(v) * |v|^ev
        """
        shock_type = string_or_blank(card, istart + 1, 'shockType')
        shock_cvt = double(card, istart + 2, 'shockCVT')
        shock_cvc = double_or_blank(card, istart + 3, 'shockCVC')
        shock_exp_vt = double_or_blank(card, istart + 4, 'shockExpVT', 1.0)
        shock_exp_vc = double_or_blank(card, istart + 5,
                                       'shockExpVC', shock_exp_vt)

        if shock_type == 'TABLE':
            shock_idts = None
            shock_idets = None
            shock_idecs = None
            shock_idetsd = None
            shock_idecsd = None
            # shock_idts = integer(card, istart + 6, 'shockIDTS')
            # shock_idets = blank(card, istart + 9, 'shockIDETS')
            # shock_idecs = blank(card, istart + 10, 'shockIDECS')
            # shock_idetsd = blank(card, istart + 11, 'shockIDETSD')
            # shock_idecsd = blank(card, istart + 12, 'shockIDECSD')
        elif shock_type == 'EQUAT':
            shock_idts = blank(card, istart + 6, 'shockIDTS')
            shock_idets = integer(card, istart + 9, 'shockIDETS')
            shock_idecs = integer_or_blank(card, istart + 10,
                                           'shockIDECS', shock_idets)
            shock_idetsd = integer(card, istart + 11, 'shockIDETSD')
            shock_idecsd = integer_or_blank(card, istart + 11,
                                            'shockIDECSD', shock_idetsd)
            #def DEquation(self):
                #if isinstance(self.dequation, int):
                    #return self.dequation
                #return self.dequation.equation_id
        else:
            msg = 'Invalid shockType=%r on card\n%s' % (shock_type, card)
            raise RuntimeError(msg)

        out = (
            shock_type, shock_cvt, shock_cvc, shock_exp_vt, shock_exp_vc,
            shock_idts, shock_idets, shock_idecs, shock_idetsd, shock_idecsd
        )
        istart += 8
        return istart, out

    @staticmethod
    def _read_spring(card, istart):
        """
        F(u) = Ft(u)
        """
        spring_type = string_or_blank(card, istart + 1, 'springType')
        spring_idt = integer(card, istart + 2, 'springIDT')

        if spring_type == 'TABLE':
            spring_idc = blank(card, istart + 3, 'springIDC')
            spring_idtdu = blank(card, istart + 4, 'springIDTDU')
            spring_idcdu = blank(card, istart + 5, 'springIDCDU')
        elif spring_type == 'EQUAT':
            spring_idc = integer_or_blank(card, istart + 3,
                                          'springIDC', spring_idt)
            spring_idtdu = integer(card, istart + 4, 'springIDTDU')
            spring_idcdu = integer_or_blank(card, istart + 5,
                                            'springIDCDU', spring_idtdu)
        else:
            msg = 'Invalid springType=%r on card\n%s' % (spring_type, card)
            raise RuntimeError(msg)

        return spring_type, spring_idt, spring_idc, spring_idtdu, spring_idcdu

    @staticmethod
    def _read_damper(card, istart):
        """
        F(v) = Ft(u)
        """
        damper_type = string_or_blank(card, istart + 1, 'damperType')
        damper_idt = integer(card, istart + 2, 'damperIDT')
        if damper_type == 'TABLE':
            damper_idc = blank(card, istart + 3, 'damperIDC')
            damper_idtdv = blank(card, istart + 4, 'damperIDTDV')
            damper_idcdv = blank(card, istart + 5, 'damperIDCDV')
        elif damper_type == 'EQUAT':
            damper_idc = integer_or_blank(card, istart + 3, 'damperIDC')
            damper_idtdv = integer(card, istart + 4, 'damperIDTDV')
            damper_idcdv = integer_or_blank(card, istart + 5, 'damperIDCDV', damper_idtdv)
        else:
            msg = 'Invalid damper_type=%r on card\n%s' % (damper_type, card)
            raise RuntimeError(msg)

        return damper_type, damper_idt, damper_idc, damper_idtdv, damper_idcdv

    @staticmethod
    def _read_gener(card, istart):
        """
        F(u, v) = Ft(u, v)
        """
        gener_idt = integer(card, istart + 2, 'generIDT')
        gener_idc = integer_or_blank(card, istart + 3,
                                     'generIDC', gener_idt)
        gener_idtdu = integer(card, istart + 4, 'generIDTDU')
        gener_idcdu = integer_or_blank(card, istart + 5,
                                       'generIDCDU', gener_idtdu)
        gener_idtdv = integer(card, istart + 6, 'generIDTDV')
        gener_idcdv = integer_or_blank(card, istart + 7,
                                       'generIDCDV', gener_idtdv)
        out = (
            gener_idt, gener_idc,
            gener_idtdu, gener_idcdu,
            gener_idtdv, gener_idcdv
        )
        return out

    def _shock_fields(self):
        list_fields = ['SHOCKA', self.shock_type, self.shock_cvt, self.shock_cvc,
                       self.shock_exp_vt, self.shock_exp_vc, self.shock_idts, None, None,
                       self.shock_idets, self.shock_idecs, self.shock_idetsd,
                       self.shock_idecsd]
        return list_fields

    def _spring_fields(self):
        list_fields = ['SPRING', self.spring_type, self.spring_idt,
                       self.spring_idc, self.spring_idtdu, self.spring_idcdu]
        return list_fields

    def _damper_fields(self):
        list_fields = ['DAMPER', self.damper_type, self.damper_idt,
                       self.damper_idc, self.damper_idtdv, self.damper_idcdv]
        return list_fields

    def _gener_fields(self):
        list_fields = ['GENER', None, self.gener_idt, self.gener_idc,
                       self.gener_idtdu, self.gener_idcdu, self.gener_idtdv,
                       self.gener_idcdv]
        return list_fields

    def raw_fields(self):
        list_fields = ['PBUSH1D', self.pid, self.k, self.c, self.m, None,
                       self.sa, self.se, None]
        for var in self.vars:
            if var == 'SHOCKA':
                list_fields += self._shock_fields()
            elif var == 'SPRING':
                list_fields += self._spring_fields()
            elif var == 'DAMPER':
                list_fields += self._damper_fields()
            elif var == 'GENER':
                list_fields += self._gener_fields()
            else:
                msg = 'var=%s not supported PBUSH1D field...' % var
                raise RuntimeError(msg)
            nspaces = 8 - (len(list_fields) - 1) % 8

            if nspaces < 8:
                list_fields += [None] * (nspaces)
        return list_fields

    def repr_fields(self):
        return self.raw_fields()

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)


class PBUSH2D(BushingProperty):
    type = 'PBUSH2D'

    def __init__(self, card=None, comment=''):
        BushingProperty.__init__(self, card)
        if comment:
            self.comment = comment
        raise NotImplementedError()

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
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)

def _append_nones(list_obj, nrequired):
    """this function has side effects"""
    n_none = nrequired - len(list_obj)
    list_obj.extend([None] * n_none)

class PBUSHT(BushingProperty):
    type = 'PBUSHT'

    @classmethod
    def _init_from_empty(cls):
        pid = 1
        k_tables = []
        b_tables = []
        ge_tables = []
        kn_tables = []
        return PBUSHT(pid, k_tables, b_tables, ge_tables, kn_tables, comment='')

    def update_by_pname_fid(self, name, value):
        if name == 'TGEID1':
            self.ge_tables[0] = value
        elif name == 'TGEID2':
            self.ge_tables[1] = value
        else:
            raise NotImplementedError('%r has not implemented update_by_pname_fid for %r' % (self.type, name))

    def __init__(self, pid, k_tables, b_tables,
                 ge_tables, kn_tables, comment=''):
        BushingProperty.__init__(self)
        if comment:
            self.comment = comment
        self.pid = pid

        _append_nones(k_tables, 6)
        _append_nones(b_tables, 6)
        _append_nones(ge_tables, 6)
        _append_nones(kn_tables, 6)

        self.k_tables = k_tables
        self.b_tables = b_tables
        self.ge_tables = ge_tables
        self.kn_tables = kn_tables

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a PBUSHT card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        k_tables = []
        b_tables = []
        ge_tables = []
        kn_tables = []

        pid = integer(card, 1, 'pid')
        nfields = len(card) - 1
        nrows = nfields // 8
        if nfields % 8 != 0:
            nrows += 1

        for irow in range(nrows):
            ifield = 1 + irow * 8
            param = string(card, ifield + 1, 'param_type')
            table = []
            for j in range(6):
                table_value = integer_or_blank(card, ifield + j + 2, param + '%i' % (j+1))
                table.append(table_value)
            if param == 'K':
                k_tables = table
            elif param == 'B':
                b_tables = table
            elif param == 'GE':
                ge_tables = table
            elif param == 'KN':
                kn_tables = table
            else:
                raise ValueError(param)
        return PBUSHT(pid, k_tables, b_tables, ge_tables, kn_tables,
                      comment=comment)

    def raw_fields(self):
        list_fields = ['PBUSHT', self.pid]
        if self.k_tables:
            list_fields += ['K'] + self.k_tables + [None]
        if self.b_tables:
            list_fields += ['B'] + self.b_tables + [None]
        if self.ge_tables:
            list_fields += ['GE'] + self.ge_tables + [None]
        if self.kn_tables:
            list_fields += ['KN'] + self.kn_tables
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)
