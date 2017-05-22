# pylint: disable=C0103,R0902,R0904,R0914
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
from six import iteritems, string_types
from six.moves import zip, range
import numpy as np

from pyNastran.utils import integer_types
from pyNastran.bdf.field_writer_8 import set_blank_if_default
from pyNastran.bdf.cards.base_card import (BaseCard, expand_thru_by)
from pyNastran.bdf.cards.deqatn import fortran_to_python_short
    #collapse_thru_by_float, condense, build_thru_float)
from pyNastran.bdf.bdf_interface.assign_type import (
    integer, integer_or_blank, integer_or_string, integer_string_or_blank,
    double, double_or_blank, string, string_or_blank,
    integer_double_or_blank, integer_double_string_or_blank,
    double_string_or_blank, interpret_value)
from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.field_writer_16 import print_card_16
from pyNastran.bdf.field_writer_double import print_card_double
from pyNastran.bdf.cards.utils import build_table_lines


def validate_dvcrel(validate, Type, cp_name):
    """
    Valdiates the DVCREL1/2

    Notes
    -----
    1.  words that start with integers (e.g., 12I/T**3) doesn't support
        strings
    """
    if validate:
        msg = 'DVCRELx: Type=%r cp_name=%r is invalid' % (Type, cp_name)
        if Type in ['CQUAD4']:
            assert cp_name in ['T1', 'T2', 'T3', 'T4'], msg # 'ZOFFS',
        elif Type in ['CTRIA3']:
            assert cp_name in [], msg
        elif Type in ['CONM2']:
            assert cp_name in ['M', 'X1', 'X2', 'X3'], msg
        elif Type in ['CBAR']:
            assert cp_name in ['X1', 'X2', 'X3'], msg
        elif Type in ['CBEAM']:
            assert cp_name in ['X1', 'X2', 'X3'], msg
        elif Type in ['CELAS1']:
            assert cp_name in [], msg
        elif Type in ['CBUSH']:
            assert cp_name in ['X1', 'X2', 'X3', 'S', 'S1', 'S2', 'S3'], msg
        else:
            raise NotImplementedError(msg)


def validate_dvmrel(validate, mat_type, mp_name):
    """
    Valdiates the DVMREL1/2

    Notes
    -----
    1.  words that start with integers (e.g., 12I/T**3) doesn't support
        strings
    """
    if validate:
        msg = 'DVMRELx: mat_type=%r mp_name=%r is invalid' % (mat_type, mp_name)
        if mat_type in ['MAT1']:
            assert mp_name in ['E', 'G', 'NU', 'GE', 'RHO', 'A', 'TREF'], msg
        elif mat_type == 'MAT2':
            assert mp_name in ['G11', 'G12', 'G22', 'G33', 'GE', 'RHO',
                               'A1', 'A2', 'A3', 'TREF'], msg
        elif mat_type == 'MAT3':
            assert mp_name in ['EX', 'ETH', 'EZ', 'NUTH', 'NUXTH', 'NUTHZ', 'NUZX', 'RHO'], msg
        elif mat_type == 'MAT8':
            assert mp_name in ['E1', 'G1Z', 'NU12'], msg
        elif mat_type == 'MAT9':
            assert mp_name in ['G11', 'G22', 'G33', 'G44', 'G55', 'G66', 'RHO'], msg
        elif mat_type == 'MAT10':
            assert mp_name in ['BULK', 'RHO'], msg
        elif mat_type == 'MAT11':
            assert mp_name in ['E1', 'E2', 'E3', 'G12', 'G13', 'G23', 'RHO'], msg
        else:
            raise NotImplementedError(msg)


def validate_dvprel(prop_type, pname_fid, validate):
    """
    Valdiates the DVPREL1/2

    Notes
    -----
    1.  words that start with integers (e.g., 12I/T**3) doesn't support
        strings
    2.  FID > 0 --> references the Property Card
    3.  FID < 0 --> references the EPT card
    """
    if validate:
        msg = 'DVPREL1: prop_type=%r pname_fid=%r is invalid' % (prop_type, pname_fid)
        #if prop_type == 'CELAS2':
            #assert pname_fid in ['K', 'GE', 'S'], msg
        #elif prop_type == 'CELAS4':
            #assert pname_fid in ['K'], msg
        if prop_type == 'PELAS':
            if pname_fid in ['K1', 3]:
                pname_fid = 'K1'
            elif pname_fid in ['GE1', 4]:
                pname_fid = 'GE1'
            else:
                raise NotImplementedError('PELAST pname_fid=%r is invalid' % pname_fid)
            #assert pname_fid in [3, 4, 'K1', 'GE1'], msg
        elif prop_type == 'PELAST':
            if pname_fid in ['TKID', 3]:
                pname_fid = 'TKID'
            else:
                raise NotImplementedError('PELAST pname_fid=%r is invalid' % pname_fid)
            assert pname_fid in [3, 4, 'TKID'], msg

        elif prop_type == 'PROD':
            if pname_fid in ['A', 4]:
                pname_fid = 'A'
            elif pname_fid in ['J', 5]:
                pname_fid = 'J'
            #elif pname_fid in ['C', 6]:
                #pname_fid = 'C'
            else:
                raise NotImplementedError('PROD pname_fid=%r is invalid' % pname_fid)
            assert pname_fid in [4, 'A', 5, 'J'], msg
        elif prop_type == 'PTUBE':
            assert pname_fid in [4, 5], msg

        #elif prop_type == 'CBAR':
            #assert pname_fid in ['X1', 'X2'], msg
        elif prop_type == 'PBAR':
            assert pname_fid in [4, 5, 6, 7, 12, 13, 14, 15, 16, 17, 18, 19, 'A', 'I1', 'J'], msg
        elif prop_type == 'PBARL':
            assert pname_fid in [12, 13, 14, 15, 16, 17, 'DIM1', 'DIM2', 'DIM3', 'DIM4'], msg

        #elif prop_type == 'CBEAM':
            #assert pname_fid in ['X1', 'X2', 'X3', 'W1A', 'W2A', 'W3A', 'W1B', 'W2B', 'W3B'], msg
        elif prop_type == 'PBEAM':
            assert pname_fid in ['I1', 'I2', 'A', 'J',
                                 'I1(B)', 'I2(B)',
                                 '-8', '-9', '-10', '-14', '-15', '-16', '-17', '-18', '-19', '-20', '-21'], msg
        elif prop_type == 'PBEAML':
            assert pname_fid in ['DIM1', 'DIM2', 'DIM3', 'DIM4', 'DIM5', 'DIM6',
                                 'DIM1(A)',
                                 'DIM1(B)', 'DIM2(B)', 'I1(B)', 'I2(B)',
                                 'NSM'], msg # 'DIM(B)'

        #elif prop_type == 'CQUAD4':
            #assert pname_fid in ['T1', 'T2', 'T3', 'T4'], msg
        elif prop_type == 'PSHELL':
            if pname_fid in ['T', 4]:
                pname_fid = 'T'
            elif pname_fid in [6]:  # 12I/T**3 doesn't support strings
                pass
            else:
                raise NotImplementedError('PSHELL pname_fid=%r is invalid' % pname_fid)
            #if cp_name in '12I/T**3':
                #cp_name =
            #assert pname_fid in ['T', 4, 6], msg
        elif prop_type == 'PCOMP':
            if isinstance(pname_fid, str):
                if pname_fid in ['Z0', 'SB']:
                    pass
                else:
                    word, num = break_word_by_trailing_integer(pname_fid)
                    if word not in ['T', 'THETA']:
                        raise RuntimeError('word=%r\n%s' % (word, msg))
            else:
                assert pname_fid in [#'Z0',
                                     3, #3-z0
                                     # 13-t1, 14-theta1, 17-t2, 18-theta2
                                     13, 14, 17, 18,
                                     23, 24, 27, 28,
                                     33, 34, 37, 38,
                                     43, 44, 47, 48], msg
        elif prop_type == 'PCOMPG':
            #if pname_fid in ['T', 4]:
                #pname_fid = 'T'
            #elif pname_fid in [6]:  # 12I/T**3 doesn't support strings
                #pass
            #else:
                #raise NotImplementedError('PSHELL pname_fid=%r is invalid' % pname_fid)
            #if cp_name in '12I/T**3':
            assert pname_fid in ['Z0', 'SB',
                                 14, 24, 34, 44,
                                 15, 25, 75, 85], msg

        #elif prop_type == 'CBUSH':
            #assert pname_fid in ['X1', 'X2', 'X3', 'S', 'S1'], msg
        elif prop_type == 'PBUSH':
            assert pname_fid in [18,
                                 'K1', 'K2', 'K3', 'K4', 'K5', 'K6',
                                 'B2',
                                 'GE1', 'GE3', 'GE4', 'GE5', 'GE6',
                                 '-13'], msg
        elif prop_type == 'PBUSH1D':
            assert pname_fid in ['K', 'C'], msg
        elif prop_type == 'PBUSHT':
            assert pname_fid in ['TBID1', 'TGEID1', 'TGEID2'], msg

        # CGAP
        elif prop_type == 'PGAP':
            assert pname_fid in [5], msg
        elif prop_type == 'PVISC':
            assert pname_fid in ['CE1'], msg

        #elif prop_type == 'CDAMP2':
            #assert pname_fid in ['B'], msg
        elif prop_type == 'PDAMP':
            assert pname_fid in [3, 'B1'], msg

        #elif prop_type == 'CMASS2':
            #assert pname_fid in ['M'], msg
        #elif prop_type == 'CMASS4':
            #assert pname_fid in ['M'], msg
        elif prop_type == 'PMASS':
            assert pname_fid in [3], msg

        #elif prop_type == 'CONM2':
            #assert pname_fid in ['M', 'X1', 'X2', 'I11', 'I22'], msg

        elif prop_type == 'PSHEAR':
            if pname_fid in ['T', 4]:
                pname_fid = 'T'
            else:
                raise NotImplementedError('PSHEAR pname_fid=%r is invalid' % pname_fid)

        elif prop_type == 'PWELD':
            assert pname_fid in ['D'], msg

        elif prop_type == 'PBEND':
            raise RuntimeError('Nastran does not support the PBEND')
        else:
            raise NotImplementedError(msg)
    return pname_fid


class OptConstraint(BaseCard):
    def __init__(self):
        pass


class DCONSTR(OptConstraint):
    """
    +---------+------+-----+------------+------------+-------+--------+
    |   1     |   2  |  3  |     4      |      5     |   6   |   7    |
    +=========+======+=====+============+============+=======+========+
    | DCONSTR | DCID | RID | LALLOW/LID | UALLOW/UID | LOWFQ | HIGHFQ |
    +---------+------+-----+------------+------------+-------+--------+
    """
    type = 'DCONSTR'
    def __init__(self, oid, dresp_id, lid=-1.e20, uid=1.e20,
                 lowfq=0., highfq=1.e20, comment=''):
        """

        """
        if comment:
            self.comment = comment
        self.oid = oid
        # DRESP entry
        self.dresp_id = dresp_id
        # lower bound
        self.lid = lid
        # upper bound
        self.uid = uid
        # low end of frequency range (Hz)
        self.lowfq = lowfq
        # high end of frequency range (Hz)
        self.highfq = highfq

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a DCONSTR card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        oid = integer(card, 1, 'oid')
        dresp_id = integer(card, 2, 'dresp_id')
        lid = integer_double_or_blank(card, 3, 'lid', -1e20)
        uid = integer_double_or_blank(card, 4, 'uid', 1e20)
        lowfq = double_or_blank(card, 5, 'lowfq', 0.0)
        highfq = double_or_blank(card, 6, 'highfq', 1e20)
        assert len(card) <= 7, 'len(DCONSTR card) = %i\ncard=%s' % (len(card), card)
        return DCONSTR(oid, dresp_id, lid, uid, lowfq, highfq, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        """
        Adds a DCONSTR card from the OP2

        Parameters
        ----------
        data : List[varies]
            a list of fields defined in OP2 format
        comment : str; default=''
            a comment for the card
        """
        oid = data[0]
        rid = data[1]
        lid = data[2]
        uid = data[3]
        lowfq = data[4]
        highfq = data[5]
        return DCONSTR(oid, rid, lid, uid, lowfq, highfq, comment=comment)

    def DRespID(self):
        if isinstance(self.dresp_id, integer_types):
            return self.dresp_id
        else:
            return self.dresp_id_ref.dresp_id

    def Rid(self):
        return self.DRespID()

    def Lid(self):
        if isinstance(self.lid, (integer_types, float)):
            return self.lid
        else:
            return self.lid_ref.tid

    def Uid(self):
        if isinstance(self.uid, (integer_types, float)):
            return self.uid
        else:
            return self.uid_ref.tid

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ' which is required by DCONSTR oid=%s' % (self.oid)
        self.dresp_id = model.DResp(self.DRespID(), msg)
        self.dresp_id_ref = self.dresp_id
        if isinstance(self.lid, integer_types):
            self.lid = model.TableD(self.lid, msg)
            self.lid_ref = self.lid
        if isinstance(self.uid, integer_types):
            self.uid = model.TableD(self.uid, msg)
            self.uid_ref = self.uid

    def uncross_reference(self):
        self.dresp_id = self.DRespID()
        self.lid = self.Lid()
        self.uid = self.Uid()

        if isinstance(self.lid, integer_types):
            del self.lid
        if isinstance(self.uid, integer_types):
            del self.uid_ref
        del self.dresp_id_ref

    def raw_fields(self):
        list_fields = ['DCONSTR', self.oid, self.DRespID(), self.Lid(),
                       self.Uid(), self.lowfq, self.highfq]
        return list_fields

    def repr_fields(self):
        lid = set_blank_if_default(self.Lid(), -1e20)
        uid = set_blank_if_default(self.Uid(), 1e20)
        lowfq = set_blank_if_default(self.lowfq, 0.0)
        highfq = set_blank_if_default(self.highfq, 1e20)
        list_fields = ['DCONSTR', self.oid, self.DRespID(), lid, uid, lowfq, highfq]
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        if is_double:
            return self.comment + print_card_double(card)
        return self.comment + print_card_16(card)


class DESVAR(OptConstraint):
    type = 'DESVAR'
    """
    +--------+-----+-------+-------+-----+-----+-------+-------+
    |    1   |  2  |    3  |   4   |  5  |  6  |    7  |   8   |
    +========+=====+=======+=======+=====+=====+=======+=======+
    | DESVAR | OID | LABEL | XINIT | XLB | XUB | DELXV | DDVAL |
    +--------+-----+-------+-------+-----+-----+-------+-------+
    """
    def __init__(self, desvar_id, label, xinit, xlb=-1e20, xub=1e20,
                 delx=None, ddval=None, comment=''):
        """
        Creates a DESVAR card

        Parameters
        ----------
        desvar_id : int
            design variable id
        label : str
            name of the design variable
        xinit : float
            the starting point value for the variable
        xlb : float; default=-1.e20
            the lower bound
        xub : float; default=1.e20
            the lower bound
        delx : float; default=1.e20
            fractional change allowed for design variables during
            approximate optimization
            NX  if blank : take from DOPTPRM; otherwise 1.0
            MSC if blank : take from DOPTPRM; otherwise 0.5
        ddval : int; default=None
            int : DDVAL id
                  allows you to set discrete values
            None : continuous
        comment : str; default=''
            a comment for the card
        """
        if comment:
            self.comment = comment
        self.desvar_id = desvar_id
        #: user-defined name for printing purposes
        self.label = label
        self.xinit = xinit
        self.xlb = xlb
        self.xub = xub
        assert xlb <= xub, 'desvar_id=%s xlb=%s xub=%s' % (desvar_id, xlb, xub)
        assert xinit >= xlb, 'desvar_id=%s xlb=%s xub=%s' % (desvar_id, xlb, xub)
        assert xinit <= xub, 'desvar_id=%s xlb=%s xub=%s' % (desvar_id, xlb, xub)

        # controls change for a single optimization cycle
        # taken from DOPTPRM if None; else default=1.
        self.delx = delx
        # DDVAL id if you want discrete values
        self.ddval = ddval

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a DESVAR card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        desvar_id = integer(card, 1, 'desvar_id')
        label = string(card, 2, 'label')
        xinit = double(card, 3, 'xinit')
        xlb = double_or_blank(card, 4, 'xlb', -1e20)
        xub = double_or_blank(card, 5, 'xub', 1e20)
        delx = double_or_blank(card, 6, 'delx')
        ddval = integer_or_blank(card, 7, 'ddval')
        assert len(card) <= 8, 'len(DESVAR card) = %i\ncard=%s' % (len(card), card)
        return DESVAR(desvar_id, label, xinit, xlb=xlb, xub=xub,
                      delx=delx, ddval=ddval, comment=comment)

    def raw_fields(self):
        list_fields = ['DESVAR', self.desvar_id, self.label, self.xinit, self.xlb,
                       self.xub, self.delx, self.ddval]
        return list_fields

    def repr_fields(self):
        """
        Gets the fields in their simplified form

        Returns
        -------
        fields : List[varies]
          the fields that define the card
        """
        xlb = set_blank_if_default(self.xlb, -1e20)
        xub = set_blank_if_default(self.xub, 1e20)

        label = self.label.strip()
        if len(label) <= 6:
            label = ' %6s ' % label
        list_fields = ['DESVAR', self.desvar_id, label, self.xinit, xlb,
                       xub, self.delx, self.ddval]
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)


class DDVAL(OptConstraint):
    """
    +-------+-----+-------+-------+-------+-------+-------+-------+-------+
    |   1   |  2  |   3   |   4   |   5   |   6   |   7   |   8   |   9   |
    +=======+=====+=======+=======+=======+=======+=======+=======+=======+
    | DDVAL | ID  | DVAL1 | DVAL2 | DVAL3 | DVAL4 | DVAL5 | DVAL6 | DVAL7 |
    +-------+-----+-------+-------+-------+-------+-------+-------+-------+
    | DDVAL | ID  | DVAL1 | THRU  | DVAL2 | BY    |  INC  |       |       |
    +-------+-----+-------+-------+-------+-------+-------+-------+-------+

    +-------+-----+-------+-------+-------+-------+-------+-------+-------+
    | DDVAL | 110 |  0.1  |  0.2  |  0.3  |  0.5  |  0.6  |  0.4  |       |
    +-------+-----+-------+-------+-------+-------+-------+-------+-------+
    |       | .7  | THRU  |  1.0  |  BY   | 0.05  |       |       |       |
    +-------+-----+-------+-------+-------+-------+-------+-------+-------+
    |       | 1.5 |  2.0  |       |       |       |       |       |       |
    +-------+-----+-------+-------+-------+-------+-------+-------+-------+
    """
    type = 'DDVAL'

    def __init__(self, oid, ddvals, comment=''):
        if comment:
            self.comment = comment
        if isinstance(ddvals, float):
            ddvals = [ddvals]
        else:
            ddvals = expand_thru_by(ddvals)
            ddvals.sort()
        self.oid = oid
        self.ddvals = ddvals

    def validate(self):
        if not isinstance(self.ddvals, list):
            raise TypeError('ddvals=%s must be a list' % self.ddvals)

        msg = ''
        for iddval, ddval in enumerate(self.ddvals):
            if not isinstance(ddval, float):
                msg += 'ddval[%i]=%s is not a float; type=%s\n' % (iddval, ddval, type(ddval))
        if msg:
            raise TypeError(msg)

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a DDVAL card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        oid = integer(card, 1, 'oid')
        n = 1
        ddvals = []
        for i in range(2, len(card)):
            ddval = double_string_or_blank(card, i, 'DDVAL%s' % n)
            if ddval is not None:
                ddvals.append(ddval)
        return DDVAL(oid, ddvals, comment=comment)

    def raw_fields(self):
        self.ddvals.sort()
        list_fields = ['DDVAL', self.oid] + self.ddvals
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        if is_double:
            return self.comment + print_card_double(card)
        return self.comment + print_card_16(card)


class DOPTPRM(OptConstraint):
    """causes a Nastran core dump if FSDMAX is nonzero and there is no stress case"""
    type = 'DOPTPRM'
    defaults = {
        'APRCOD' : 2,
        'AUTOSE' : 0,
        'CONV1' : 0.001,
        'CONV2' : 1e-20,
        'CONVDV' : 0.001,  # 0.0001 for topology optimization
        'CONVPR' : 0.001,
        'CT' : -0.03,
        'CTMIN' : 0.003,
        'DELB' : 0.0001,
        'DELP' : 0.2,
        'DELX' : 0.5,  # 0.2 for topology optimization
        'DLXESL' : 0.5,
        'DESMAX' : 5,  # 30 for topology optimization
        'DISCOD' : 1,
        'DISBEG' : 0,
        'DPMAX' : 0.5,
        'DPMIN' : 0.01,
        'DRATIO' : 0.1,
        'DSMXESL' : 20,
        'DXMAX' : 1.0,
        'DXMIN' : 0.05,  # 1e-5 for topology optimization
        'ETA1' : 0.01,
        'ETA2' : 0.25,
        'ETA3' : 0.7,
        'FSDALP' : 0.9,
        'FSDMAX' : 0,
        'GMAX' : 0.005,
        'GSCAL' : 0.001,
        'IGMAX' : 0,
        'IPRINT' : 0,
        'ISCAL' : 0,
        'METHOD' : 0,
        'NASPRO' : 0,
        'OBJMOD' : 0,
        'OPTCOD' : 0,
        'P1' : 0,
        'P2' : 1,

        'P2CALL' : None,  #: .. todo:: DEFAULT PRINTS ALL???
        'P2CBL' : None,   #: .. todo:: DEFAULT PRINTS ALL???
        'P2CC' : None,    #: .. todo:: DEFAULT PRINTS ALL???
        'P2CDDV' : None,  #: .. todo:: DEFAULT PRINTS ALL???
        'P2CM' : None,    #: .. todo:: DEFAULT PRINTS ALL???
        'P2CP' : None,    #: .. todo:: DEFAULT PRINTS ALL???
        'P2CR' : None,    #: .. todo:: DEFAULT PRINTS ALL???
        'P2RSET' : None,  #: .. todo:: DEFAULT PRINTS ALL???
        'PENAL' : 0.0,
        'PLVIOL' : 0,
        'PTOL' : 1e35,
        'STPSCL' : 1.0,
        'TCHECK' : -1,
        'TDMIN' : None,  #: .. todo:: ???
        'TREGION' : 0,
        'UPDFAC1' : 2.0,
        'UPDFAC2' : 0.5,
        }
    def __init__(self, params, comment=''):
        """
        Design Optimization Parameters
        Overrides default values of parameters used in design optimization

        +---------+--------+------+--------+------+--------+------+--------+------+
        |    1    |    2   |   3  |   4    |  5   |    6   |   7  |    8   |   9  |
        +=========+========+======+========+======+========+======+========+======+
        | DOPTPRM | PARAM1 | VAL1 | PARAM2 | VAL2 | PARAM3 | VAL3 | PARAM4 | VAL4 |
        +---------+--------+------+--------+------+--------+------+--------+------+
        |         | PARAM5 | VAL5 | -etc.- |      |        |      |        |      |
        +---------+--------+------+--------+------+--------+------+--------+------+
        """
        if comment:
            self.comment = comment
        self.params = params

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a DOPTPRM card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        nfields = len(card) - 1
        params = {}
        for i in range(0, nfields, 2):
            param = string_or_blank(card, i + 1, 'param')
            default_value = None
            if param is None:
                continue
            if param in cls.defaults:
                default_value = cls.defaults[param]
            val = integer_double_or_blank(card, i + 2, '%s_value' % param, default_value)
            params[param] = val
        return DOPTPRM(params, comment=comment)

    def raw_fields(self):
        list_fields = ['DOPTPRM']
        for param, val in sorted(iteritems(self.params)):
            list_fields += [param, val]
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        if is_double:
            return self.comment + print_card_double(card)
        return self.comment + print_card_16(card)


class DLINK(OptConstraint):
    type = 'DLINK'

    def __init__(self, oid, ddvid, IDv, Ci, c0=0., cmult=1., comment=''):
        """
        Multiple Design Variable Linking
        Relates one design variable to one or more other design variables.

        +-------+------+-------+--------+-------+------+----+------+----+
        |   1   |   2  |   3   |   4    |   5   |   6  |  7 |   8  | 9  |
        +=======+======+=======+========+=======+======+====+======+====+
        | DLINK |  ID  | DDVID |   C0   | CMULT | IDV1 | C1 | IDV2 | C2 |
        +-------+------+-------+--------+-------+------+----+------+----+
        |       | IDV3 |   C3  | -etc.- |       |      |    |      |    |
        +-------+------+-------+--------+-------+------+----+------+----+
        """
        if comment:
            self.comment = comment
        self.oid = oid
        self.ddvid = ddvid
        self.c0 = c0
        self.cmult = cmult
        self.IDv = IDv
        self.Ci = Ci

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a DLINK card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        oid = integer(card, 1, 'oid')
        ddvid = integer(card, 2, 'ddvid')
        c0 = double_or_blank(card, 3, 'c0', 0.)
        cmult = double_or_blank(card, 4, 'cmult', 1.)

        nfields = len(card) - 4
        n = nfields // 2
        idvs = []
        Ci = []

        for i in range(n):
            j = 2 * i + 5
            idv = integer(card, j, 'IDv' + str(i))
            ci = double(card, j + 1, 'Ci' + str(i))
            idvs.append(idv)
            Ci.append(ci)
        return DLINK(oid, ddvid, idvs, Ci, c0=c0, cmult=cmult, comment=comment)

    def raw_fields(self):
        list_fields = ['DLINK', self.oid, self.ddvid, self.c0, self.cmult]
        for (idv, ci) in zip(self.IDv, self.Ci):
            list_fields += [idv, ci]
        return list_fields

    def repr_fields(self):
        c0 = set_blank_if_default(self.c0, 0.)
        cmult = set_blank_if_default(self.cmult, 1.)
        list_fields = ['DLINK', self.oid, self.ddvid, c0, cmult]
        for (idv, ci) in zip(self.IDv, self.Ci):
            list_fields += [idv, ci]
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        if is_double:
            return self.comment + print_card_double(card)
        return self.comment + print_card_16(card)


def validate_dresp1(property_type, response_type, atta, attb, atti):
    property_types = [
        'ELEM',
        #'PELAS',
        'PROD',
        # 'PTUBE',
        'PBAR', 'PBARL',
        'PBEAM', #'PBEAML',
        #'PSHEAR',
        'PSHELL', 'PCOMP',
        #'PSOLID',
        #'PKNL',
        None,
    ]
    stress_types = ['PBARL']

    response_types = [
        #'WEIGHT', 'VOLUME',

        # static
        #'DISP', 'STRESS', 'COMP', 'SPCFORCE', 'ESE', 'GPFORCP',

        # aero
        #'FLUTTER', 'STABDER', 'TRIM', 'DIVERG',

        # ???
        #'CEIG', 'CSTRAT',

        # modal
        #'EIGN', 'LAMA',

        # time
        #'TDISP', 'TVELO', 'TACCL', 'TFORC', 'TSPCF',

        # frequency
        #'FREQ',
        #'FRDISP', 'FRVELO', 'FRACCL',
        #'FRSPCF', 'FRMASS', 'FRFORC', 'FRSTRE',

        # random
        #'PSDDISP', 'PSDACCL',
        #'RMSDISP', 'RMSVELO', 'RMSACCL',

        #'CFAILURE', 'TOTSE',
    ]
    msg = 'DRESP1 ptype=%s rtype=%s atta=%s attb=%s atti=%s' % (
        property_type, response_type, atta, attb, atti)
    #print(msg)
    if property_type not in property_types:
        if property_type == 'ELEM':
            if response_type == 'CFAILURE':
                pass
            elif response_type == 'STRESS':
                pass
            else:
                msg = 'DRESP1 ptype=%s rtype=%s atta=%s attb=%s atti=%s' % (
                    property_type, response_type, atta, attb, atti)
                raise RuntimeError(msg)
        #elif property_type is stress_properties:
            #msg = 'DRESP1 ptype=%s rtype=%s atta=%s attb=%s atti=%s' % (
                #property_type, response_type, atta, attb, atti)
            #raise RuntimeError(msg)
        #else:
            #msg = 'DRESP1 ptype=%s rtype=%s atta=%s attb=%s atti=%s' % (
                #property_type, response_type, atta, attb, atti)
            raise RuntimeError(msg)

    if response_type == 'FLUTTER':
        #print(msg)
        assert property_type in [None, 'PKNL'], 'DRESP1 ptype=%r rtype=%r atta=%r attb=%r atti=%r' % (
            property_type, response_type, atta, attb, atti)
        assert atta is None, 'DRESP1 ptype=%s rtype=%s atta=%s attb=%s atti=%s' % (
            property_type, response_type, atta, attb, atti)
        assert attb is None, 'DRESP1 ptype=%s rtype=%s atta=%s attb=%s atti=%s' % (
            property_type, response_type, atta, attb, atti)
        assert len(atti) == 4, 'DRESP1 ptype=%s rtype=%s atta=%s attb=%s atti=%s' % (
            property_type, response_type, atta, attb, atti)
    elif property_type is None:
        if response_type == 'WEIGHT':
            assert atta in [1, 2, 3, 4, 5, 6, None], 'DRESP1 ptype=%s rtype=%s atta=%s attb=%s atti=%s' % (
                property_type, response_type, atta, attb, atti)
            assert attb in [1, 2, 3, 4, 5, 6, None], 'DRESP1 ptype=%s rtype=%s atta=%s attb=%s atti=%s' % (
                property_type, response_type, atta, attb, atti)
            if len(atti) == 0:
                atti = ['ALL']
            for attii in atti:
                if attii != 'ALL':
                    assert isinstance(attii, integer_types), 'DRESP1 ptype=%s rtype=%s atta=%s attb=%s atti=%s' % (
                        property_type, response_type, atta, attb, atti)
                    assert attii >= 0, 'DRESP1 ptype=%s rtype=%s atta=%s attb=%s atti=%s' % (
                        property_type, response_type, atta, attb, atti)

        elif response_type == 'VOLUME':
            assert atta is None, 'DRESP1 ptype=%s rtype=%s atta=%s attb=%s atti=%s' % (
                property_type, response_type, atta, attb, atti)
            assert attb is None, 'DRESP1 ptype=%s rtype=%s atta=%s attb=%s atti=%s' % (
                property_type, response_type, atta, attb, atti)
            if len(atti) == 0:
                atti = ['ALL']
            for attii in atti:
                assert attii in [0, 'ALL'], 'DRESP1 ptype=%s rtype=%s atta=%s attb=%s atti=%s' % (
                    property_type, response_type, atta, attb, atti)

        elif response_type == 'DISP':
            atta = str(atta)
            for attai in atta:
                assert atta in '123456', 'DRESP1 ptype=%s rtype=%s atta=%s attb=%s atti=%s' % (
                    property_type, response_type, atta, attb, atti)  # 8???
            assert attb is None, 'DRESP1 ptype=%s rtype=%s atta=%s attb=%s atti=%s' % (
                property_type, response_type, atta, attb, atti)
            assert len(atti) == 1, 'DRESP1 ptype=%s rtype=%s atta=%s attb=%s atti=%s' % (
                property_type, response_type, atta, attb, atti)

        elif response_type in ['FRDISP', 'FRVELO', 'FRACCL', 'FRSPCF']:  # frequency displacement
            assert atta in [1, 2, 3, 4, 5, 7, 8, 9], 'DRESP1 ptype=%s rtype=%s atta=%s attb=%s atti=%s' % (
                property_type, response_type, atta, attb, atti)

            if attb is None or isinstance(attb, float):
                # blank is all frequencies
                # float is a specific frequency
                pass
            else:
                assert attb in ['SUM', 'AVG', 'SSQ', 'RSS', 'MAX', 'MIN', 'AVE'], 'DRESP1 ptype=%s rtype=%s atta=%s attb=%s atti=%s' % (
                    property_type, response_type, atta, attb, atti)  # remove AVE?

            assert len(atti) >= 1, 'DRESP1 ptype=%s rtype=%s atta=%s attb=%s atti=%s' % (
                property_type, response_type, atta, attb, atti)
            for attii in atti:
                assert isinstance(attii, integer_types), 'DRESP1 ptype=%s rtype=%s atta=%s attb=%s atti=%s' % (
                    property_type, response_type, atta, attb, atti)
                assert attii > 0, 'DRESP1 ptype=%s rtype=%s atta=%s attb=%s atti=%s' % (
                    property_type, response_type, atta, attb, atti)

        elif response_type in ['TDISP', 'TVELO', 'TACCL', 'TSPCF']:  # time displacement
            atta = str(atta)
            for attai in atta:
                assert atta in '123456', 'DRESP1 ptype=%s rtype=%s atta=%s attb=%s atti=%s' % (
                    property_type, response_type, atta, attb, atti)  # 8???

            if attb is None or isinstance(attb, float):
                # blank is all times
                # float is a specific times
                pass
            else:
                assert attb in ['SUM', 'AVG', 'SSQ', 'RSS', 'MAX', 'MIN'], 'DRESP1 ptype=%s rtype=%s atta=%s attb=%s atti=%s' % (
                    property_type, response_type, atta, attb, atti)

            assert len(atti) >= 1, 'DRESP1 ptype=%s rtype=%s atta=%s attb=%s atti=%s' % (
                property_type, response_type, atta, attb, atti)
            for attii in atti:
                assert isinstance(attii, integer_types), 'DRESP1 ptype=%s rtype=%s atta=%s attb=%s atti=%s' % (
                    property_type, response_type, atta, attb, atti)
                assert attii > 0, 'DRESP1 ptype=%s rtype=%s atta=%s attb=%s atti=%s' % (
                    property_type, response_type, atta, attb, atti)

        elif response_type in ['RMSDISP', 'RMSACCL']:
            assert atta in [1, 2, 3], 'DRESP1 ptype=%s rtype=%s atta=%s attb=%s atti=%s' % (
                property_type, response_type, atta, attb, atti)
            assert attb in [140, 1000], 'DRESP1 ptype=%s rtype=%s atta=%s attb=%s atti=%s' % (
                property_type, response_type, atta, attb, atti)
            assert len(atti) == 1, 'DRESP1 ptype=%s rtype=%s atta=%s attb=%s atti=%s' % (
                property_type, response_type, atta, attb, atti)

        elif response_type == 'CEIG':
            assert isinstance(atta, integer_types), 'DRESP1 ptype=%s rtype=%s atta=%s attb=%s atti=%s' % (
                property_type, response_type, atta, attb, atti)
            assert atta > 0, 'DRESP1 ptype=%s rtype=%s atta=%s attb=%s atti=%s' % (
                property_type, response_type, atta, attb, atti)
            if attb is None:
                attb = 'ALPHA'
            assert attb in ['ALPHA', 'OMEGA'], 'DRESP1 ptype=%s rtype=%s atta=%s attb=%s atti=%s' % (
                property_type, response_type, atta, attb, atti)
            assert len(atti) == 0, 'DRESP1 ptype=%s rtype=%s atta=%s attb=%s atti=%s' % (
                property_type, response_type, atta, attb, atti)

        elif response_type in ['EIGN', 'LAMA']: # EIGEN as well?
            assert isinstance(atta, integer_types), 'DRESP1 ptype=%s rtype=%s atta=%s attb=%s atti=%s' % (
                property_type, response_type, atta, attb, atti)
            assert atta > 0, 'DRESP1 ptype=%s rtype=%s atta=%s attb=%s atti=%s' % (
                property_type, response_type, atta, attb, atti)
            # 1 -> direct linerization (default)
            # 2 -> inverse approximation
            assert attb in [None, 1, 2], 'DRESP1 ptype=%s rtype=%s atta=%s attb=%s atti=%s' % (
                property_type, response_type, atta, attb, atti)
            assert len(atti) == 0, 'DRESP1 ptype=%s rtype=%s atta=%s attb=%s atti=%s' % (
                property_type, response_type, atta, attb, atti)
        elif response_type == 'FREQ':
            assert isinstance(atta, integer_types), 'DRESP1 ptype=%s rtype=%s atta=%s attb=%s atti=%s' % (
                property_type, response_type, atta, attb, atti)
            assert atta > 0, 'DRESP1 ptype=%s rtype=%s atta=%s attb=%s atti=%s' % (
                property_type, response_type, atta, attb, atti)
            assert attb is None, 'DRESP1 ptype=%s rtype=%s atta=%s attb=%s atti=%s' % (
                property_type, response_type, atta, attb, atti)
            assert len(atti) == 0, 'DRESP1 ptype=%s rtype=%s atta=%s attb=%s atti=%s' % (
                property_type, response_type, atta, attb, atti)

        else:
            msg = 'DRESP1 ptype=%s rtype=%s atta=%s attb=%s atti=%s' % (
                property_type, response_type, atta, attb, atti)
            raise RuntimeError(msg)

    elif response_type == 'CFAILURE':
        if property_type == 'ELEM':
            assert atta in [5, 7], 'DRESP1 ptype=%s rtype=%s atta=%s attb=%s atti=%s' % (
                property_type, response_type, atta, attb, atti)
        elif property_type == 'PCOMP':
            assert atta in [3, 5], 'DRESP1 ptype=%s rtype=%s atta=%s attb=%s atti=%s' % (
                property_type, response_type, atta, attb, atti)
        else:
            msg = 'DRESP1 ptype=%s rtype=%s atta=%s attb=%s atti=%s' % (
                property_type, response_type, atta, attb, atti)
            raise RuntimeError(msg)
        if attb is None:
            attb = 1
        assert isinstance(attb, integer_types), 'DRESP1 ptype=%s rtype=%s atta=%s attb=%s atti=%s' % (
            property_type, response_type, atta, attb, atti)
        assert attb > 0, 'DRESP1 ptype=%s rtype=%s atta=%s attb=%s atti=%s' % (
            property_type, response_type, atta, attb, atti)
        assert len(atti) > 0, 'DRESP1 ptype=%s rtype=%s atta=%s attb=%s atti=%s' % (
            property_type, response_type, atta, attb, atti)
    elif response_type in ['CSTRAIN', 'CSTRESS']:
        if attb is None:
            attb = 1
        if property_type == 'PCOMP':
            # 11 - max shear stress/strain
            assert atta in [3, 11], 'DRESP1 ptype=%s rtype=%s atta=%s attb=%s atti=%s' % (
                property_type, response_type, atta, attb, atti)
            assert len(atti) > 0, 'DRESP1 ptype=%s rtype=%s atta=%s attb=%s atti=%s' % (
                property_type, response_type, atta, attb, atti)
        else:
            msg = 'DRESP1 ptype=%s rtype=%s atta=%s attb=%s atti=%s' % (
                property_type, response_type, atta, attb, atti)
            raise RuntimeError(msg)
        assert isinstance(attb, integer_types), 'DRESP1 ptype=%s rtype=%s atta=%s attb=%s atti=%s' % (
            property_type, response_type, atta, attb, atti)
        assert attb > 0, 'DRESP1 ptype=%s rtype=%s atta=%s attb=%s atti=%s' % (
            property_type, response_type, atta, attb, atti)
        assert len(atti) > 0, 'DRESP1 ptype=%s rtype=%s atta=%s attb=%s atti=%s' % (
            property_type, response_type, atta, attb, atti)

    elif response_type == 'STRESS':
        if property_type == 'PBARL':
            assert atta in [2, 3, 4, 5, 7, 8], 'DRESP1 ptype=%s rtype=%s atta=%s attb=%s atti=%s' % (
                property_type, response_type, atta, attb, atti)
            assert attb is None, 'DRESP1 ptype=%s rtype=%s atta=%s attb=%s atti=%s' % (
                property_type, response_type, atta, attb, atti)
        elif property_type == 'PBAR':
            assert atta in [2, 6, 7, 8, 14, 15], 'DRESP1 ptype=%s rtype=%s atta=%s attb=%s atti=%s' % (
                property_type, response_type, atta, attb, atti)
            assert attb is None, 'DRESP1 ptype=%s rtype=%s atta=%s attb=%s atti=%s' % (
                property_type, response_type, atta, attb, atti)
        elif property_type == 'PBEAM':
            assert atta in [6, 8, 9, 31, 59, 108], 'DRESP1 ptype=%s rtype=%s atta=%s attb=%s atti=%s' % (
                property_type, response_type, atta, attb, atti)
            assert attb is None, 'DRESP1 ptype=%s rtype=%s atta=%s attb=%s atti=%s' % (
                property_type, response_type, atta, attb, atti)
        elif property_type == 'PROD':
            assert atta in [2, 3, 7], 'DRESP1 ptype=%s rtype=%s atta=%s attb=%s atti=%s' % (
                property_type, response_type, atta, attb, atti)
            assert attb is None, 'DRESP1 ptype=%s rtype=%s atta=%s attb=%s atti=%s' % (
                property_type, response_type, atta, attb, atti)

        elif property_type == 'PSHELL':
            assert atta in [4, 5, 6, 7, 8, 9, 15, 16, 17, 19, 26, 28, 35, 36], 'DRESP1 ptype=%s rtype=%s atta=%s attb=%s atti=%s' % (
                property_type, response_type, atta, attb, atti)
            assert attb is None, 'DRESP1 ptype=%s rtype=%s atta=%s attb=%s atti=%s' % (
                property_type, response_type, atta, attb, atti)
        elif property_type == 'PCOMP':
            # this is a SMEAR PCOMP, which is basically a PSHELL
            assert atta in [9, 11, 17], 'DRESP1 ptype=%s rtype=%s atta=%s attb=%s atti=%s' % (
                property_type, response_type, atta, attb, atti)
            assert attb is None, 'DRESP1 ptype=%s rtype=%s atta=%s attb=%s atti=%s' % (
                property_type, response_type, atta, attb, atti)
        #elif property_type in stress_types:
            #if not isinstance(atta, integer_types):
                #msg = 'DRESP1 ptype=%s rtype=%s atta=%s attb=%s atti=%s; atta should be an integer' % (
                    #property_type, response_type, atta, attb, atti)
                #raise TypeError(msg)
            #assert attb is None, 'DRESP1 ptype=%s rtype=%s atta=%s attb=%s atti=%s' % (
                #property_type, response_type, atta, attb, atti)
        else:
            msg = 'DRESP1 ptype=%s rtype=%s atta=%s attb=%s atti=%s' % (
                property_type, response_type, atta, attb, atti)
            raise RuntimeError(msg)

        assert attb is None, 'DRESP1 ptype=%s rtype=%s atta=%s attb=%s atti=%s; atta should be an integer' % (
            property_type, response_type, atta, attb, atti)
        assert len(atti) > 0, 'DRESP1 ptype=%s rtype=%s atta=%s attb=%s atti=%s' % (
            property_type, response_type, atta, attb, atti)
    elif response_type == 'FORCE':
        if property_type == 'PROD':
            assert atta in [2], 'DRESP1 ptype=%s rtype=%s atta=%s attb=%s atti=%s' % (
                property_type, response_type, atta, attb, atti)
            assert attb is None, 'DRESP1 ptype=%s rtype=%s atta=%s attb=%s atti=%s' % (
                property_type, response_type, atta, attb, atti)
        elif response_type in ['FRSTRE']:
            if property_type == 'PROD':
                assert atta in [4], 'DRESP1 ptype=%s rtype=%s atta=%s attb=%s atti=%s' % (
                    property_type, response_type, atta, attb, atti)
                assert attb is None, 'DRESP1 ptype=%s rtype=%s atta=%s attb=%s atti=%s' % (
                    property_type, response_type, atta, attb, atti)
            else:
                msg = 'DRESP1 ptype=%s rtype=%s atta=%s attb=%s atti=%s' % (
                    property_type, response_type, atta, attb, atti)
                raise RuntimeError(msg)
            assert len(atti) > 0, 'DRESP1 ptype=%s rtype=%s atta=%s attb=%s atti=%s' % (
                property_type, response_type, atta, attb, atti)
        else:
            msg = 'DRESP1 ptype=%s rtype=%s atta=%s attb=%s atti=%s' % (
                property_type, response_type, atta, attb, atti)
            raise RuntimeError(msg)
        assert len(atti) > 0, 'DRESP1 ptype=%s rtype=%s atta=%s attb=%s atti=%s' % (
            property_type, response_type, atta, attb, atti)

    else:
        msg = 'DRESP1 ptype=%s rtype=%s atta=%s attb=%s atti=%s' % (
            property_type, response_type, atta, attb, atti)
        raise RuntimeError(msg)
    return atta, attb, atti

class DRESP1(OptConstraint):
    """
    +--------+-------+---------+--------+--------+--------+-------+------+-------+
    |   1    |  2    |     3   |   4    |   5    |   6    |   7   |   8  |   9   |
    +========+=======+=========+========+========+========+=======+======+=======+
    | DRESP1 |  OID  | LABEL   | RTYPE  | PTYPE  | REGION | ATTA  | ATTB | ATTI  |
    +--------+-------+---------+--------+--------+--------+-------+------+-------+

    +--------+-------+---------+--------+--------+--------+-------+------+-------+
    | DRESP1 |  103  |  S02    | STRESS | PSHELL |        |  9    |      |   3   |
    +--------+-------+---------+--------+--------+--------+-------+------+-------+
    | DRESP1 |  1S1  | CSTRAIN | PCOMP  |        |        | 1     | 1    | 10000 |
    +--------+-------+---------+--------+--------+--------+-------+------+-------+
    """
    type = 'DRESP1'

    def __init__(self, dresp_id, label, response_type, property_type, region,
                 atta, attb, atti, comment='', validate=False):
        """
        Creates a DRESP1 card.

        A DRESP1 is used to define a "simple" output result that may be
        optimized on.  A simple result is a result like stress, strain,
        force, displacement, eigenvalue, etc. for a node/element that
        may be found in a non-optimization case.

        Parameters
        ----------
        dresp_id : int
            response id
        lable : str
            Name of the response
        response_type : str
            Response type
        property_type : str
            Element flag (PTYPE = 'ELEM'), or property entry name, or panel
            flag for ERP responses (PTYPE = 'PANEL' - See Remark 34), or
            RANDPS ID. Blank for grid point responses. 'ELEM' or property
            name used only with element type responses (stress, strain,
            force, etc.) to identify the relevant element IDs, or the property
            type and relevant property IDs.

            Must be {ELEM, PBAR, PSHELL, PCOMP, PANEL, etc.)
            PTYPE = RANDPS ID when RTYPE=PSDDISP, PSDVELO, or PSDACCL.
        region : str
            Region identifier for constraint screening
        atta : int / float / str / blank
            Response attribute
        attb : int / float / str / blank
            Response attribute
        atti : List[int / float / str]
            the response values to pull from
            List[int]:
                list of grid ids
                list of property ids
            List[str]
                'ALL'
        comment : str; default=''
            a comment for the card
        validate : bool; default=False
            should the card be validated when it's created

        Example 1
        ---------
        dresp_id = 103
        label = 'resp1'
        response_type = 'STRESS'
        property_type = 'PSHELL'
        pid = 3
        atta = 9 # von mises upper surface stress
        region = None
        attb = None
        atti = [pid]
        DRESP1(dresp_id, label, response_type, property_type, region, atta, attb, atti)

        Example 2
        ---------
        dresp_id = 104
        label = 'resp2'
        response_type = 'STRESS'
        property_type = 'PCOMP'
        pid = 3
        layer = 4
        atta = 9 # von mises upper surface stress
        region = None
        attb = layer
        atti = [pid]
        DRESP1(dresp_id, label, response_type, property_type, region, atta, attb, atti)

        Example 3 (???)
        ---------------
        dresp_id = 105
        label = 'resp3'
        response_type = 'DISP'
        #atta = ???
        #region = ???
        #attb = ???
        atti = [nid]
        DRESP1(dresp_id, label, response_type, property_type, region, atta, attb, atti)

        Example 4 (???)
        ---------------
        dresp_id = 105
        label = 'resp3'
        response_type = 'ELEM'
        #atta = ???
        #region = ???
        #attb = ???
        atti = [eid???]
        DRESP1(dresp_id, label, response_type, property_type, region, atta, attb, atti)
        """
        if comment:
            self.comment = comment
        self.dresp_id = dresp_id
        self.label = label

        # DISP, VONMISES, 14
        self.response_type = response_type

        # PSHELL, PCOMP, PBAR, etc.
        self.property_type = property_type
        self.region = region
        if isinstance(atti, integer_types):
            atti = [atti]
        assert isinstance(atti, list), 'atti=%s type=%s' % (atti, type(atti))

        if validate:
            atta, attb, atti = validate_dresp1(property_type, response_type, atta, attb, atti)

        self.atta = atta
        self.attb = attb
        self.atti = atti

    @property
    def rtype(self):
        return self.response_type

    @rtype.setter
    def rtype(self, response_type):
        self.response_type = response_type

    @property
    def ptype(self):
        return self.property_type

    @ptype.setter
    def ptype(self, property_type):
        self.property_type = property_type

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a DRESP1 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        oid = integer(card, 1, 'oid')
        label = string(card, 2, 'label')
        response_type = string(card, 3, 'rtype')

        # elem, pbar, pshell, etc. (ELEM flag or Prop Name)
        property_type = integer_string_or_blank(card, 4, 'ptype')
        region = integer_or_blank(card, 5, 'region')
        atta = integer_double_string_or_blank(card, 6, 'atta')
        attb = integer_double_string_or_blank(card, 7, 'attb')

        atti = []
        for i in range(8, len(card)):
            attii = integer_double_string_or_blank(card, i, 'atti_%i' % (i + 1))
            atti.append(attii)
        return DRESP1(oid, label, response_type, property_type, region, atta, attb, atti,
                      comment=comment, validate=False)

    def _verify(self, xref=True):
        pass

    def calculate(self, op2_model, subcase_id):
        rtype = self.rtype
        property_type = self.property_type
        if rtype == 'DISP' and property_type is None:
            msg = 'fields=%s\n' % (self.raw_fields())
            msg += 'rtype=%r ptype=%r region=%s A=%r B=%r\n' % (
                self.rtype, self.property_type, self.region,
                self.atta, self.attb)
            assert isinstance(self.atta, integer_types), self.atta
            assert self.attb is None, self.attb
            assert len(self.atti) == 1, self.atti
            #msg += str(self)
            #raise NotImplementedError(msg)
            comp = self.atta
            assert comp < 7, comp
            case = op2_model.displacements[subcase_id]
            nids = case.node_gridtype[:, 0]
            nidsi = np.array(
                [node if isinstance(node, integer_types) else node.nid
                 for node in self.atti], dtype='int32')
            inid = np.searchsorted(nids, nidsi)[0]
            itime = 0 #  TODO:  ???
            #print('atti = ', self.atti)
            atti = self.atti[0]
            out = case.data[itime, inid, comp - 1]
            #print(' DISP[itime=%s, nid=%s, comp=%i] = %s' % (itime, str(nidsi), comp, out))
        elif property_type == 'ELEM' and 0:
            pass
        else:
            msg = 'fields=%s\n' % (self.raw_fields())
            msg += 'response_type=%r property_type=%r region=%s\n' % (
                self.response_type, self.property_type, self.region)
            msg += str(self)
            raise NotImplementedError(msg)
        return out

    def OptID(self):
        return self.DRespID()

    def DRespID(self):
        return self.dresp_id

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ', which is required by DRESP1 dresp_id=%s' % (self.dresp_id)
        msg += '\n' + str(self)

        op2_results = [
            'VOLUME', 'LAMA', 'FRSPCF', 'TRIM', 'ESE', 'SPCFORCE', 'FRMASS',
            'CFAILURE', 'CSTRAT', 'STRESS', 'DIVERG', 'TOTSE', 'COMP',
            'TACCL', 'RMSACCL',
            'RMSVELO',
            'PSDDISP', 'RMSDISP',
            'TFORC', 'FRFORC',
            'TSPCF',
        ]
        if self.property_type in ['ELEM']:
            self.atti = model.Elements(self.atti, msg=msg)
            self.atti_ref = self.atti
        elif self.property_type in ['PSHELL', 'PBAR', 'PROD', 'PCOMP',
                                    'PSOLID', 'PELAS', 'PBARL', 'PBEAM',
                                    'PBEAML', 'PSHEAR', 'PTUBE',]:
            self.atti = model.Properties(self.atti, msg=msg)
            self.atti_ref = self.atti
        elif self.response_type in ['FRSTRE']:
            self.atti = model.Properties(self.atti, msg=msg)
            self.atti_ref = self.atti
        elif self.response_type in ['WEIGHT', 'STABDER', 'CEIG', 'EIGN', 'FREQ']:
            pass
        elif self.response_type == 'FLUTTER':
            # TODO: SOL-200; add check that FLFACT values exist in the FLFACT card
            #       referenced by the FLUTTER card for the given subcase
            if self.property_type == 'PKNL':
                self.atti = [
                    model.Set(self.atti[0], msg=msg),
                    model.FLFACT(self.atti[1], msg=msg),
                    model.FLFACT(self.atti[2], msg=msg),
                    model.FLFACT(self.atti[3], msg=msg),
                ]
                msgi = 'max density=%s mach=%s velocity=%s' % (
                    self.atti[1].max(), self.atti[2].max(), self.atti[3].max())
                #print(msgi)
                self.atti_ref = self.atti
            elif self.property_type is None:
                pass
            else:
                msg = 'PropertyType=%r is not supported\n' % self.property_type
                msg += str(self)
                print(msg)
        elif self.response_type in ['DISP',
                                    'TDISP',
                                    'TVELO',
                                    'FRDISP', 'FRVELO', 'FRACCL',
                                    'PSDVELO', 'PSDACCL']:
            self.atti = model.Nodes(self.atti, msg=msg)
            self.atti_ref = self.atti
        elif self.response_type in op2_results:
            pass
        elif self.response_type == 'ERP':
            assert self.property_type == 'PANEL'
        elif self.response_type == 'GPFORCP':
            assert len(self.atti) >= 1, 'atti=%r\n%s' % (self.atti, self)
            self.atta = model.Node(self.atta, msg=msg)
            self.atti = model.Nodes(self.atti, msg=msg)
        elif self.response_type == 'GPFORCE':
            #self.atta = component
            eids = [eid for eid in self.atti if eid is not None]
            self.atti = model.Elements(eids, msg=msg)
        else:
            msg = 'response_type=%r ptype=%r\n' % (self.response_type, self.property_type)
            msg += str(self)
            raise NotImplementedError(msg)

    def uncross_reference(self):
        self.atti = self.atti_values()
        if hasattr(self, 'atti_ref'):
            del self.atti_ref

    def atti_values(self):
        #return self.atti
        #if self.response_type in ['ELEM']:
            #self.atti = model.Elements(self.atti, msg=msg)
            #pass
        op2_results = [
            'VOLUME', 'CEIG', 'LAMA', 'FREQ', 'FRSPCF', 'TRIM', 'ESE', 'SPCFORCE', 'FRMASS',
            'CFAILURE', 'CSTRAT', 'STRESS', 'DIVERG', 'TOTSE', 'COMP',
            'TACCL', 'RMSACCL',
            'PSDDISP', 'RMSDISP',
            'TVELO', 'RMSVELO',
            'TFORC', 'FRFORC',
            'TSPCF',
        ]
        if self.property_type in ['ELEM']:
            data = [elem if isinstance(elem, integer_types) else elem.eid for elem in self.atti]
        elif self.property_type in ['PSHELL', 'PBAR', 'PROD', 'PCOMP',
                                    'PSOLID', 'PELAS', 'PBARL', 'PBEAM',
                                    'PBEAML', 'PSHEAR', 'PTUBE',
                                    'FRSTRE']:
            data = [prop if isinstance(prop, integer_types) else prop.pid for prop in self.atti]
            for value in data:
                assert not isinstance(value, BaseCard), value
        elif self.response_type == 'FRSTRE':
            data = [prop if isinstance(prop, integer_types) else prop.Pid() for prop in self.atti]
            for value in data:
                assert not isinstance(value, BaseCard), value
        elif self.response_type in ['WEIGHT', 'STABDER', 'EIGN', 'FREQ']:
            data = self.atti
        elif self.response_type == 'FLUTTER':
            if self.property_type == 'PKNL':
                data = [atti if isinstance(atti, integer_types) else atti.sid for atti in self.atti]
            else:
                data = self.atti
                #msg = 'PropertyType=%r is not supported\n' % self.property_type
                #msg += str(self)
                #print(msg)
                #raise NotImplementedError(msg)

        elif self.response_type in ['DISP',
                                    'TDISP', 'TVELO',
                                    'FRDISP', 'FRVELO', 'FRACCL',
                                    'PSDVELO', 'PSDACCL']:
            #self.atti = model.Nodes(self.atti, msg=msg)
            data = [node if isinstance(node, integer_types) else node.nid for node in self.atti]
        elif self.response_type in ['FRFORC', 'TFORC',
                                    'STRESS', 'ESE', 'CFAILURE', 'CSTRAIN']:
            data = [elem if isinstance(elem, integer_types) else elem.eid for elem in self.atti]
        elif self.response_type in op2_results:
            data = self.atti
            for value in data:
                assert not isinstance(value, BaseCard), 'response_type=%s value=%s' % (self.response_type, value)
        elif self.response_type in ['GPFORCP']:
            # MSC Nastran specific
            data = [node if isinstance(node, integer_types) else node.nid for node in self.atti]
        elif self.response_type in ['GPFORCE']:
            # MSC
            data = [elem if isinstance(elem, integer_types) else elem.eid for elem in self.atti
                    if elem is not None]
        elif self.response_type in ['ERP']:
            msg = 'response_type=%r property_type=%r atta=%r attb=%r atti=%r\n' % (
                self.response_type, self.property_type, self.atta, self.attb, self.atti)
            assert self.property_type == 'PANEL', msg
            # atta=None attb=3000.0 atti=[555]
            data = self.atti
        else:
            msg = 'response_type=%r property_type=%r atta=%r attb=%r atti=%r\n' % (
                self.response_type, self.property_type, self.atta, self.attb, self.atti)
            #msg += str(self)
            raise NotImplementedError(msg)
        return data

    def raw_fields(self):
        list_fields = ['DRESP1', self.dresp_id, self.label, self.response_type, self.property_type,
                       self.region, self.atta, self.attb] + self.atti_values()
        return list_fields

    def repr_fields(self):
        """
        Gets the fields in their simplified form

        Returns
        -------
        fields : List[varies]
          the fields that define the card
        """
        label = self.label.strip()
        if len(label) <= 6:
            label = ' %6s ' % label
        list_fields = ['DRESP1', self.dresp_id, self.label, self.response_type, self.property_type,
                       self.region, self.atta, self.attb] + self.atti_values()
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        if is_double:
            return self.comment + print_card_double(card)
        return self.comment + print_card_16(card)


class DRESP2(OptConstraint):
    """
    Design Sensitivity Equation Response Quantities
    Defines equation responses that are used in the design, either as
    constraints or as an objective.

    +--------+---------+--------+-----------+--------+--------+--------+--------+--------+
    |   1    |    2    |    3   |     4     |    5   |    6   |    7   |    8   |    9   |
    +========+=========+========+===========+========+========+========+========+========+
    | DRESP2 |   ID    |  LABEL | EQID/FUNC | REGION | METHOD |   C1   |   C2   |   C3   |
    +--------+---------+--------+-----------+--------+--------+--------+--------+--------+
    |        | DESVAR  | DVID1  |   DVID2   |  DVID3 |  DVID4 |  DVID5 |  DVID6 |  DVID7 |
    +--------+---------+--------+-----------+--------+--------+--------+--------+--------+
    |        |         | DVID8  |   etc.    |        |        |        |        |        |
    +--------+---------+--------+-----------+--------+--------+--------+--------+--------+
    |        | DTABLE  | LABL1  |  LABL2    |  LABL3 |  LABL4 |  LABL5 |  LABL6 |  LABL7 |
    +--------+---------+--------+-----------+--------+--------+--------+--------+--------+
    |        |         | LABL8  |  etc.     |        |        |        |        |        |
    +--------+---------+--------+-----------+--------+--------+--------+--------+--------+
    |        | DRESP1  |  NR1   |    NR2    |   NR3  |   NR4  |   NR5  |   NR6  |  NR7   |
    +--------+---------+--------+-----------+--------+--------+--------+--------+--------+
    |        |         |  NR8   |   etc.    |        |        |        |        |        |
    +--------+---------+--------+-----------+--------+--------+--------+--------+--------+
    |        | DNODE   |   G1   |    C1     |   G2   |   C2   |   G3   |   C3   |        |
    +--------+---------+--------+-----------+--------+--------+--------+--------+--------+
    |        |         |   G4   |    C4     |  etc.  |        |        |        |        |
    +--------+---------+--------+-----------+--------+--------+--------+--------+--------+
    |        | DVPREL1 | DPIP1  |   DPIP2   | DPIP3  | DPIP4  | DPIP5  | DPIP6  | DPIP7  |
    +--------+---------+--------+-----------+--------+--------+--------+--------+--------+
    |        |         | DPIP8  |   DPIP9   |  etc.  |        |        |        |        |
    +--------+---------+--------+-----------+--------+--------+--------+--------+--------+
    |        | DVCREL1 | DCIC1  |   DCIC2   | DCIC3  | DCIC4  | DCIC5  | DCIC6  | DCIC7  |
    +--------+---------+--------+-----------+--------+--------+--------+--------+--------+
    |        |         | DCIC8  |   DCIC9   |  etc.  |        |        |        |        |
    +--------+---------+--------+-----------+--------+--------+--------+--------+--------+
    |        | DVMREL1 | DMIM1  |   DMIM2   | DMIM3  | DMIM4  | DMIM5  | DMIM6  | DMIM7  |
    +--------+---------+--------+-----------+--------+--------+--------+--------+--------+
    |        |         | DMIM8  |   DMIM9   |  etc.  |        |        |        |        |
    +--------+---------+--------+-----------+--------+--------+--------+--------+--------+
    |        | DVPREL2 | DPI2P1 |   DPI2P2  | DPI2P3 | DPI2P4 | DPI2P5 | DPI2P6 | DPI2P7 |
    +--------+---------+--------+-----------+--------+--------+--------+--------+--------+
    |        |         | DPI2P8 |   DPI2P9  |  etc.  |        |        |        |        |
    +--------+---------+--------+-----------+--------+--------+--------+--------+--------+
    |        | DVCREL2 | DCI2C1 |   DCI2C2  | DCI2C3 | DCI2C4 | DCI2C5 | DCI2C6 | DCI2C7 |
    +--------+---------+--------+-----------+--------+--------+--------+--------+--------+
    |        |         | DCI2C8 |   DCI2C9  |   etc. |        |        |        |        |
    +--------+---------+--------+-----------+--------+--------+--------+--------+--------+
    |        | DVMREL2 | DMI2M1 |   DMI2M2  | DMI2M3 | DMI2M4 | DMI2M5 | DMI2M6 | DMI2M7 |
    +--------+---------+--------+-----------+--------+--------+--------+--------+--------+
    |        |         | DMI2M8 |   DMI2M9  |   etc. |        |        |        |        |
    +--------+---------+--------+-----------+--------+--------+--------+--------+--------+
    |        | DRESP2  | NRR1   |   NRR2    |  NRR3  |  NRR4  |  NRR5  |  NRR6  |  NRR7  |
    +--------+---------+--------+-----------+--------+--------+--------+--------+--------+
    |        |         | NRR8   |   etc.    |        |        |        |        |        |
    +--------+---------+--------+-----------+--------+--------+--------+--------+--------+
    |        | DVLREL1 | DLIL1  |   DLIL2   |  DLIL3 |  DLIL4 |  DLIL5 |  DLIL6 |  DLIL7 |
    +--------+---------+--------+-----------+--------+--------+--------+--------+--------+
    |        |         | DLIL8  |   etc.    |        |        |        |        |        |
    +--------+---------+--------+-----------+--------+--------+--------+--------+--------+

    C1, C2, C3 are MSC specific
    """
    type = 'DRESP2'

    def __init__(self, dresp_id, label, dequation, region, params,
                 method='MIN', c1=1., c2=0.005, c3=10., comment='',
                 validate=False):
        """
        Creates a DRESP2 card.

        A DRESP2 is used to define a "complex" output result that may be
        optimized on.  A complex result is a result that uses:
          - simple (DRESP1) results
          - complex (DRESP2) results
          - default values (DTABLE)
          - DVCRELx values
          - DVMRELx values
          - DVPRELx values
          - DESVAR values
        Then, an equation (DEQATN) is used to formulate an output response.

        Parameters
        ----------
        dresp_id : int
            response id
        lable : str
            Name of the response
        dequation : int
            DEQATN id
        region : str
            Region identifier for constraint screening
        params : dict[(index, card_type)] = values
            the storage table for the response function
            index : int
                a counter
            card_type : str
                the type of card to pull from
                DESVAR, DVPREL1, DRESP2, etc.
            values : List[int]
                the values for this response
        method : str; default=MIN
            flag used for FUNC=BETA/MATCH
            FUNC = BETA
                valid options are {MIN, MAX}
            FUNC = MATCH
                valid options are {LS, BETA}
        c1 / c2 / c3 : float; default=1. / 0.005 / 10.0
            constants for FUNC=BETA or FUNC=MATCH
        comment : str; default=''
            a comment for the card
        validate : bool; default=False
            should the card be validated when it's created

        params = {
           (0, 'DRESP1') = [10, 20],
           (1, 'DESVAR') = [30],
           (2, 'DRESP1') = [40],
        }
        """
        if comment:
            self.comment = comment
        self.func = None
        self.dequation_str = None
        self.dresp_id = dresp_id
        self.label = label
        self.dequation = dequation
        self.region = region
        self.method = method
        #if validate:
            #atta, attb, atti = validate_dresp1(
                #property_type, response_type, atta, attb, atti)
        self.c1 = c1
        self.c2 = c2
        self.c3 = c3
        self.params = params

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a DRESP2 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        dresp_id = integer(card, 1, 'dresp_id')
        label = string(card, 2, 'label')
        dequation = integer_or_string(card, 3, 'dequation_id')
        region = integer_or_blank(card, 4, 'region')
        method = string_or_blank(card, 5, 'method', 'MIN')
        c1 = double_or_blank(card, 6, 'c1', 1.)
        c2 = double_or_blank(card, 7, 'c2', 0.005)
        c3 = double_or_blank(card, 8, 'c3', 10.)

        fields = [interpret_value(field) for field in card[9:]]

        # DRESP2, dresp_id,
        #         DRESP1, 10, 20
        #         DESVAR, 30
        #         DRESP1, 40
        # params = {
        #    (0, 'DRESP1') = [10, 20],
        #    (1, 'DESVAR') = [30],
        #    (2, 'DRESP1') = [40],
        # }
        params = parse_table_fields(fields)

        #print("--DRESP2 Params--")
        #for key, value_list in sorted(iteritems(params)):
            #print("  key=%s value_list=%s" %(key, value_list))
        return DRESP2(dresp_id, label, dequation, region, params,
                      method, c1, c2, c3, comment=comment)

    def OptID(self):
        return self.DRespID()

    def DRespID(self):
        return self.dresp_id

    def _verify(self, xref=True):
        pass
        #for (j, name), value_list in sorted(iteritems(self.params)):
            #print('  DRESP2 verify - key=%s values=%s' % (name,
                #self._get_values(name, value_list)))

    def calculate(self, op2_model, subcase_id):
        argsi = []
        for key, vals in sorted(iteritems(self.params)):
            j, name = key
            if name in ['DRESP1', 'DRESP2']:
                for val in vals:
                    arg = val.calculate(op2_model, subcase_id)
                    argsi.append(arg)
            elif name in ['DVMREL1', 'DVMREL2']:
                for val in vals:
                    arg = val.calculate(op2_model, subcase_id)
                    argsi.append(arg)
            elif name in ['DVPREL1', 'DVPREL2']:
                for val in vals:
                    arg = val.calculate(op2_model, subcase_id)
                    argsi.append(arg)
            elif name == 'DTABLE':
                for val in vals:
                    arg = self.dtable[val]
                    argsi.append(arg)
            else:
                raise NotImplementedError('  TODO: xref %s' % str(key))
        #op2_model.log.info('DRESP2 args = %s' % argsi)
        out = self.func(*argsi)
        op2_model.log.info('  deqatn out = %s' % out)
        return out

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        if model.dtable is not None:
            model.log.debug(model.dtable.rstrip())
        msg = ', which is required by DRESP2 ID=%s' % (self.dresp_id)
        default_values = {}
        for key, vals in sorted(iteritems(self.params)):
            try:
                j, name = key
            except:
                raise RuntimeError(str(self))
            #print(j, name)
            if name in ['DRESP1', 'DRESP2']:
                for i, val in enumerate(vals):
                    self.params[key][i] = model.DResp(val, msg)
            elif name in ['DVMREL1', 'DVMREL2']:
                for i, val in enumerate(vals):
                    self.params[key][i] = model.DVmrel(val, msg)
            elif name in ['DVPREL1', 'DVPREL2']:
                for i, val in enumerate(vals):
                    self.params[key][i] = model.DVprel(val, msg)
            elif name == 'DESVAR':
                for i, val in enumerate(vals):
                    self.params[key][i] = model.Desvar(val, msg)
            elif name == 'DTABLE':
                #model.log.info('bdf_filename = %s' % model.bdf_filename)
                #model.log.info('\n' + model.dtable.rstrip())
                self.dtable = model.dtable
                #print('dtable =', self.dtable)
                for i, val in enumerate(vals):
                    default_values[val] = self.dtable[val]
            else:
                raise NotImplementedError('  TODO: xref %s' % str(key))
        self.params_ref = self.params

        if isinstance(self.DEquation(), integer_types):
            self.dequation = model.DEQATN(self.dequation, msg=msg)
            self.dequation_ref = self.dequation
            self.func = self.dequation_ref.func
        elif isinstance(self.dequation, str):
            self.func = fortran_to_python_short(self.dequation, default_values)
        else:
            raise NotImplementedError(self.dequation)

    def uncross_reference(self):
        if hasattr(self, 'func'):
            del self.func

        self.dequation = self.DEquation()
        if isinstance(self.dequation, integer_types):
            del self.dequation_ref

    def DEquation(self):
        if isinstance(self.dequation, (int, string_types)):
            return self.dequation
        return self.dequation_ref.equation_id

    def _get_values(self, name, values_list):
        out = []
        if name in ['DRESP1', 'DRESP2']:
            for i, val in enumerate(values_list):
                #self.params[key][i] = model.DResp(val, msg)
                if isinstance(val, integer_types):
                    out.append(val)
                else:
                    out.append(val.OptID())

        elif name in ['DVCREL1', 'DVCREL2']:
            for i, val in enumerate(values_list):
                if isinstance(val, integer_types):
                    out.append(val)
                else:
                    out.append(val.OptID())

        elif name in ['DVMREL1', 'DVMREL2']:
            for i, val in enumerate(values_list):
                if isinstance(val, integer_types):
                    out.append(val)
                else:
                    out.append(val.OptID())

        elif name in ['DVPREL1', 'DVPREL2']:
            for i, val in enumerate(values_list):
                if isinstance(val, integer_types):
                    out.append(val)
                else:
                    out.append(val.OptID())
        elif name == 'DESVAR':
            for i, val in enumerate(values_list):
                if isinstance(val, integer_types):
                    out.append(val)
                else:
                    out.append(val.desvar_id)
        elif name == 'DNODE':
            #print(values_list)
            for i in range(0, len(values_list), 2):
                val = values_list[i]
                if isinstance(val, integer_types):
                    out.append(val)
                else:
                    out.append(val.nid)
                out.append(values_list[i+1])

        elif name == 'DTABLE':
            out = values_list
        else:
            raise NotImplementedError('  TODO: _get_values %s' % str(name))
            #out = values_list
        return out

    def _pack_params(self):
        # # the amount of padding at the [beginning,end] of the 2nd line
        pack_length = {
            'DESVAR' : [1, 0],
            'DTABLE' : [1, 0],
            'DFRFNC' : [1, 0],
            'DRESP1' : [1, 0],
            'DNODE' : [1, 1],  # unique entry
            'DVPREL1' : [1, 0],
            'DVCREL1' : [1, 0],
            'DVMREL1' : [1, 0],
            'DVPREL2' : [1, 0],
            'DVCREL2' : [1, 0],
            'DVMREL2' : [1, 0],
            'DRESP2' : [1, 0],
        }

        list_fields = []
        for (j, name), value_list in sorted(iteritems(self.params)):
            values_list2 = self._get_values(name, value_list)
            fields2 = [name] + values_list2
            #try:
            (i, j) = pack_length[name]
            #except KeyError:
                #msg = 'INVALID DRESP2 name=%r fields=%s ID=%s' % (name, value_list, self.oid)
                #raise KeyError(msg)
            list_fields += build_table_lines(fields2, nstart=i, nend=j)
        return list_fields

    def raw_fields(self):
        list_fields = ['DRESP2', self.dresp_id, self.label, self.DEquation(),
                       self.region, self.method, self.c1, self.c2, self.c3]
        list_fields += self._pack_params()
        return list_fields

    def repr_fields(self):
        method = set_blank_if_default(self.method, 'MIN')
        c1 = set_blank_if_default(self.c1, 100.)
        c2 = set_blank_if_default(self.c2, 0.005)

        list_fields = ['DRESP2', self.dresp_id, self.label, self.DEquation(),
                       self.region, method, c1, c2, self.c3]
        list_fields += self._pack_params()
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        if is_double:
            return self.comment + print_card_double(card)
        return self.comment + print_card_16(card)


class DRESP3(OptConstraint):
    """
    +--------+---------+--------+--------+--------+--------+--------+--------+--------+
    |   1    |    2    |    3   |   4    |    5   |    6   |    7   |    8   |    9   |
    +========+=========+========+========+========+========+========+========+========+
    | DRESP3 |   ID    |  LABEL | GROUP  |  TYPE  | REGION |        |        |        |
    +--------+---------+--------+--------+--------+--------+--------+--------+--------+
    |        | DESVAR  | DVID1  | DVID2  |  DVID3 |  DVID4 |  DVID5 |  DVID6 |  DVID7 |
    +--------+---------+--------+--------+--------+--------+--------+--------+--------+
    |        |         | DVID8  |  etc.  |        |        |        |        |        |
    +--------+---------+--------+--------+--------+--------+--------+--------+--------+
    |        | DTABLE  | LABL1  |  LABL2 |  LABL3 |  LABL4 |  LABL5 |  LABL6 |  LABL7 |
    +--------+---------+--------+--------+--------+--------+--------+--------+--------+
    |        |         | LABL8  |  etc.  |        |        |        |        |        |
    +--------+---------+--------+--------+--------+--------+--------+--------+--------+
    |        | DRESP1  |  NR1   |  NR2   |   NR3  |   NR4  |   NR5  |   NR6  |  NR7   |
    +--------+---------+--------+--------+--------+--------+--------+--------+--------+
    |        |         |  NR8   | etc.   |        |        |        |        |        |
    +--------+---------+--------+--------+--------+--------+--------+--------+--------+
    |        | DNODE   |   G1   |   C1   |   G2   |   C2   |   G3   |   C3   |        |
    +--------+---------+--------+--------+--------+--------+--------+--------+--------+
    |        |         |   G4   |   C4   |  etc.  |        |        |        |        |
    +--------+---------+--------+--------+--------+--------+--------+--------+--------+
    |        | DVPREL1 | DPIP1  | DPIP2  | DPIP3  | DPIP4  | DPIP5  | DPIP6  | DPIP7  |
    +--------+---------+--------+--------+--------+--------+--------+--------+--------+
    |        |         | DPIP8  | DPIP9  |  etc.  |        |        |        |        |
    +--------+---------+--------+--------+--------+--------+--------+--------+--------+
    |        | DVCREL1 | DCIC1  | DCIC2  | DCIC3  | DCIC4  | DCIC5  | DCIC6  | DCIC7  |
    +--------+---------+--------+--------+--------+--------+--------+--------+--------+
    |        |         | DCIC8  | DCIC9  |  etc.  |        |        |        |        |
    +--------+---------+--------+--------+--------+--------+--------+--------+--------+
    |        | DVMREL1 | DMIM1  | DMIM2  | DMIM3  | DMIM4  | DMIM5  | DMIM6  | DMIM7  |
    +--------+---------+--------+--------+--------+--------+--------+--------+--------+
    |        |         | DMIM8  | DMIM9  |  etc.  |        |        |        |        |
    +--------+---------+--------+--------+--------+--------+--------+--------+--------+
    |        | DVPREL2 | DPI2P1 | DPI2P2 | DPI2P3 | DPI2P4 | DPI2P5 | DPI2P6 | DPI2P7 |
    +--------+---------+--------+--------+--------+--------+--------+--------+--------+
    |        |         | DPI2P8 | DPI2P9 |  etc.  |        |        |        |        |
    +--------+---------+--------+--------+--------+--------+--------+--------+--------+
    |        | DVCREL2 | DCI2C1 | DCI2C2 | DCI2C3 | DCI2C4 | DCI2C5 | DCI2C6 | DCI2C7 |
    +--------+---------+--------+--------+--------+--------+--------+--------+--------+
    |        |         | DCI2C8 | DCI2C9 |   etc. |        |        |        |        |
    +--------+---------+--------+--------+--------+--------+--------+--------+--------+
    |        | DVMREL2 | DMI2M1 | DMI2M2 | DMI2M3 | DMI2M4 | DMI2M5 | DMI2M6 | DMI2M7 |
    +--------+---------+--------+--------+--------+--------+--------+--------+--------+
    |        |         | DMI2M8 | DMI2M9 |   etc. |        |        |        |        |
    +--------+---------+--------+--------+--------+--------+--------+--------+--------+
    |        | DRESP2  | NRR1   | NRR2   |  NRR3  |  NRR4  |  NRR5  |  NRR6  |  NRR7  |
    +--------+---------+--------+--------+--------+--------+--------+--------+--------+
    |        |         | NRR8   | etc.   |        |        |        |        |        |
    +--------+---------+--------+--------+--------+--------+--------+--------+--------+
    |        | DVLREL1 | DLIL1  | DLIL2  |  DLIL3 |  DLIL4 |  DLIL5 |  DLIL6 |  DLIL7 |
    +--------+---------+--------+--------+--------+--------+--------+--------+--------+
    |        |         | DLIL8  | etc.   |        |        |        |        |        |
    +--------+---------+--------+--------+--------+--------+--------+--------+--------+
    |        | USRDATA | String                                                       |
    +--------+---------+--------------------------------------------------------------+
    |        |         | etc.                                                         |
    +--------+---------+--------------------------------------------------------------+
    """
    type = 'DRESP3'

    def __init__(self, dresp_id, label, group, Type, region, params,
                 comment=''):
        if comment:
            self.comment = comment
        self.dresp_id = dresp_id
        self.label = label
        self.group = group
        self.Type = Type
        self.region = region
        self.params = params

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a DRESP3 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        dresp_id = integer(card, 1, 'dresp_id')
        label = string(card, 2, 'label')
        group = string(card, 3, 'group')
        Type = string(card, 4, 'Type')
        region = integer_or_blank(card, 5, 'region')

        i = 0
        list_fields = [interpret_value(field) for field in card[9:]]
        key = '$NULL$'  # dummy key
        params = {key: []}
        value_list = []
        for (i, field) in enumerate(list_fields):
            if i % 8 == 0 and field is not None:
                params[key] = value_list
                key = field
                value_list = []
            elif field is not None:
                value_list.append(field)
            #else:
            #    pass
        params[key] = value_list
        del params['$NULL$']
        return DRESP3(dresp_id, label, group, Type, region, params,
                      comment=comment)

    def _get_values(self, name, values_list):
        out = []
        if name in ['DRESP1', 'DRESP2']:
            for i, val in enumerate(values_list):
                #self.params[key][i] = model.DResp(val, msg)
                if isinstance(val, integer_types):
                    out.append(val)
                else:
                    out.append(val.OptID())

        elif name in ['DVCREL1', 'DVCREL2']:
            for i, val in enumerate(values_list):
                if isinstance(val, integer_types):
                    out.append(val)
                else:
                    out.append(val.OptID())

        elif name in ['DVMREL1', 'DVMREL2']:
            for i, val in enumerate(values_list):
                if isinstance(val, integer_types):
                    out.append(val)
                else:
                    out.append(val.OptID())

        elif name in ['DVPREL1', 'DVPREL2']:
            for i, val in enumerate(values_list):
                if isinstance(val, integer_types):
                    out.append(val)
                else:
                    out.append(val.OptID())
        elif name == 'DESVAR':
            for i, val in enumerate(values_list):
                if isinstance(val, integer_types):
                    out.append(val)
                else:
                    out.append(val.desvar_id)
        elif name == 'DNODE':
            #print(values_list)
            for i in range(0, len(values_list), 2):
                val = values_list[i]
                if isinstance(val, integer_types):
                    out.append(val)
                else:
                    out.append(val.nid)
                out.append(values_list[i+1])

        elif name == 'DTABLE':
            out = values_list
        else:
            raise NotImplementedError('  TODO: _get_values %s' % str(name))
            #out = values_list
        return out

    def _pack_params(self):
        # # the amount of padding at the [beginning,end] of the 2nd line
        pack_length = {
            'DESVAR' : [1, 0],
            'DTABLE' : [1, 0],
            'DFRFNC' : [1, 0],
            'DRESP1' : [1, 0],
            'DNODE' : [1, 1],  # unique entry
            'DVPREL1' : [1, 0],
            'DVCREL1' : [1, 0],
            'DVMREL1' : [1, 0],
            'DVPREL2' : [1, 0],
            'DVCREL2' : [1, 0],
            'DVMREL2' : [1, 0],
            'DRESP2' : [1, 0],
            'USRDATA' : [1, 0],
        }
        list_fields = []
        for key, value_list in sorted(iteritems(self.params)):
            value_list2 = self._get_values(key, value_list)
            fields2 = [key] + value_list2

            try:
                (i, j) = pack_length[key]
            except KeyError:
                msg = 'INVALID DRESP3 key=%r fields=%s ID=%s' % (key, value_list, self.dresp_id)
                raise KeyError(msg)
            list_fields += build_table_lines(fields2, nstart=i, nend=j)
        return list_fields

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ', which is required by DRESP3 ID=%s' % (self.dresp_id)
        default_values = {}
        for name, vals in sorted(iteritems(self.params)):
            if name in ['DRESP1', 'DRESP2']:
                for i, val in enumerate(vals):
                    self.params[name][i] = model.DResp(val, msg)
            elif name in ['DVMREL1', 'DVMREL2']:
                for i, val in enumerate(vals):
                    self.params[name][i] = model.DVmrel(val, msg)
            elif name in ['DVPREL1', 'DVPREL2']:
                for i, val in enumerate(vals):
                    self.params[name][i] = model.DVprel(val, msg)
            elif name == 'DESVAR':
                for i, val in enumerate(vals):
                    self.params[name][i] = model.Desvar(val, msg)
            elif name == 'DTABLE':
                self.dtable = model.dtable
                for i, val in enumerate(vals):
                    default_values[val] = self.dtable[val]
            else:
                raise NotImplementedError('  TODO: xref %s' % str(name))
        self.params_ref = self.params

        #if isinstance(self.DEquation(), integer_types):
            #self.dequation = model.DEQATN(self.dequation, msg=msg)
            #self.dequation_ref = self.dequation
            #self.func = self.dequation_ref.func
        #elif isinstance(self.dequation, str):
            #self.func = fortran_to_python_short(self.dequation, default_values)
        #else:
            #raise NotImplementedError(self.dequation)

    def raw_fields(self):
        list_fields = [
            'DRESP3', self.dresp_id, self.label, self.group, self.Type, self.region,
            None, None, None]
        list_fields += self._pack_params()
        return list_fields

    def repr_fields(self):
        list_fields = [
            'DRESP3', self.dresp_id, self.label, self.group, self.Type, self.region,
            None, None, None]
        list_fields += self._pack_params()
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        if is_double:
            return self.comment + print_card_double(card)
        return self.comment + print_card_16(card)


class DCONADD(OptConstraint):
    type = 'DCONADD'

    def __init__(self, oid, dconstrs, comment=''):
        if comment:
            self.comment = comment
        self.oid = oid
        self.dconstrs = dconstrs

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a DCONADD card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        oid = integer(card, 1, 'dcid')
        dconstrs = []

        for i in range(1, len(card)):
            dconstr = integer(card, i, 'dconstr_%i' % i)
            dconstrs.append(dconstr)
        return DCONADD(oid, dconstrs, comment=comment)

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        #self.dconstrs = [model.dconstrs[oid] for oid in self.dconstr_ids]
        self.dconstrs_ref = [model.dconstrs[oid] for oid in self.dconstr_ids]

    def uncross_reference(self):
        self.dconstrs = self.dconstr_ids
        del self.dconstrs_ref

    @property
    def dconstr_ids(self):
        ids = []
        for dconstr in self.dconstrs:
            if isinstance(dconstr, list):
                for dconstri in dconstr:
                    if isinstance(dconstri, integer_types):
                        ids.append(dconstri)
                    elif isinstance(dconstri, DCONSTR):
                        ids.append(dconstri.oid)
                        break
                    else:
                        print("type=%s; dconstri=%s\n" % (type(dconstri), dconstri))
                        raise NotImplementedError(dconstri)
            elif isinstance(dconstr, integer_types):
                ids.append(dconstr)
            #elif isinstance(dconstr, DCONSTR):
                #ids.append(dconstr.oid)
            else:
                print("type=%s; dconstr=%s\n" % (type(dconstr), dconstr))
                raise NotImplementedError(dconstr)
        return ids
        #return [dconstr if isinstance(dconstr, integer_types) else dconstr.oid
                #for dconstr in self.dconstrs]

    def raw_fields(self):
        list_fields = ['DCONADD', self.oid] + self.dconstr_ids
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        if is_double:
            return self.comment + print_card_double(card)
        return self.comment + print_card_16(card)


class DSCREEN(OptConstraint):
    type = 'DSCREEN'

    def __init__(self, rType, trs=-0.5, nstr=20, comment=''):
        if comment:
            self.comment = comment

        #: Response type for which the screening criteria apply. (Character)
        self.rType = rType
        #: Truncation threshold. (Real; Default = -0.5)
        self.trs = trs
        #: Maximum number of constraints to be retained per region per load
        #: case. (Integer > 0; Default = 20)
        self.nstr = nstr

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a DSCREEN card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        rType = string(card, 1, 'rType')
        trs = double_or_blank(card, 2, 'trs', -0.5)
        nstr = integer_or_blank(card, 3, 'nstr', 20)
        assert len(card) == 4, 'len(DSCREEN card) = %i\ncard=%s' % (len(card), card)
        return DSCREEN(rType, trs=trs, nstr=nstr, comment=comment)

    def raw_fields(self):
        list_fields = ['DSCREEN', self.rType, self.trs, self.nstr]
        return list_fields

    def repr_fields(self):
        trs = set_blank_if_default(self.trs, -0.5)
        nstr = set_blank_if_default(self.nstr, 20)
        list_fields = ['DSCREEN', self.rType, trs, nstr]
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        if is_double:
            return self.comment + print_card_double(card)
        return self.comment + print_card_16(card)


class DVCREL1(OptConstraint):  # similar to DVMREL1
    type = 'DVCREL1'

    def __init__(self, oid, Type, eid, cp_name, dvids, coeffs,
                 cp_min=None, cp_max=1e20, c0=0., validate=False, comment=''):
        """
        +---------+--------+--------+--------+-----------+-------+--------+-----+---+
        |   1     |    2   |   3    |    4   |     5     |   6   |   7    |  8  | 9 |
        +=========+========+========+========+===========+=======+========+=====+===+
        | DVCREL1 |   ID   |  TYPE  |  EID   |   CPNAME  | CPMIN |  CPMAX |  C0 |   |
        +---------+--------+--------+--------+-----------+-------+--------+-----+---+
        |         | DVID1  | COEF1  | DVID2  |   COEF2   | DVID3 | -etc.- |     |   |
        +---------+--------+--------+--------+-----------+-------+--------+-----+---+

        +---------+--------+--------+--------+-------+-----+------+
        | DVCREL1 | 200000 | CQUAD4 | 1      | ZOFFS |     |  1.0 |
        +---------+--------+--------+--------+-------+-----+------+
        |         | 200000 |   1.0  |        |       |     |      |
        +---------+--------+--------+--------+-------+-----+------+
        """
        if comment:
            self.comment = comment
        if isinstance(dvids, integer_types):
            dvids = [dvids]
        if isinstance(coeffs, float):
            coeffs = [coeffs]

        self.oid = oid

        # element type (e.g. CQUAD4)
        self.Type = Type

        # element id
        self.eid = eid

        # the connectivity property (X1, X2, X3, ZOFFS, etc.)
        self.cp_name = cp_name

        # min value
        self.cp_min = cp_min

        # max value
        self.cp_max = cp_max

        # offset coefficient
        self.c0 = c0

        # DESVAR id
        self.dvids = dvids

        # scale factor for DESVAR
        self.coeffs = coeffs

        assert len(coeffs) > 0, 'len(coeffs)=%s' % len(coeffs)
        assert len(coeffs) == len(dvids), 'len(coeffs)=%s len(dvids)=%s' % (len(coeffs), len(dvids))
        validate_dvcrel(validate, Type, cp_name)

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a DVCREL1 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        oid = integer(card, 1, 'oid')
        Type = string(card, 2, 'Type')
        eid = integer(card, 3, 'eid')
        cp_name = integer_or_string(card, 4, 'cp_name')

        cp_min = double_or_blank(card, 5, 'cp_min', None)
        cp_max = double_or_blank(card, 6, 'cp_max', 1e20)
        c0 = double_or_blank(card, 7, 'c0', 0.0)

        dvids = []
        coeffs = []
        end_fields = [interpret_value(field) for field in card[9:]]

        nfields = len(end_fields) - 1
        if nfields % 2 == 1:
            end_fields.append(None)
            nfields += 1
        i = 0
        for i in range(0, nfields, 2):
            dvids.append(end_fields[i])
            coeffs.append(end_fields[i + 1])
        if nfields % 2 == 1:
            print(card)
            print("dvids = %s" % (dvids))
            print("coeffs = %s" % (coeffs))
            raise RuntimeError('invalid DVCREL1...')
        return DVCREL1(oid, Type, eid, cp_name, dvids, coeffs,
                       c0=c0, cp_min=cp_min, cp_max=cp_max, comment=comment)

    def _verify(self, xref):
        """
        Verifies all methods for this object work

        Parameters
        ----------
        xref : bool
            has this model been cross referenced
        """
        pass

    def OptID(self):
        return self.oid

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ', which is required by DVCREL1 name=%r' % self.type
        if self.Type in ['CQUAD4', 'CTRIA3', 'CBAR', 'CBEAM', 'CELAS1', 'CELAS2', 'CELAS4',
                         'CDAMP2', 'CGAP', 'CBUSH']:
            self.eid = model.Element(self.eid, msg=msg)
        elif self.Type in ['CONM1', 'CONM2', 'CMASS2', 'CMASS4']:
            self.eid = model.masses[self.eid]
        #elif Type in ['CBUSH']:
        else:
            raise NotImplementedError(self.Type)
        self.eid_ref = self.eid
        self.dvids = [model.Desvar(dvid, msg) for dvid in self.dvids]

    def uncross_reference(self):
        self.eid = self.Eid()
        del self.eid_ref

    def calculate(self, op2_model, subcase_id):
        raise NotImplementedError('\n' + str(self))

    def Eid(self):
        if isinstance(self.eid, integer_types):
            return self.eid
        return self.eid_ref.eid

    @property
    def desvar_ids(self):
        if isinstance(self.dvids[0], integer_types):
            return self.dvids
        return [desvar.desvar_id for desvar in self.dvids]

    def raw_fields(self):
        list_fields = ['DVCREL1', self.oid, self.Type, self.Eid(),
                       self.cp_name, self.cp_min, self.cp_max, self.c0, None]
        for (dvid, coeff) in zip(self.desvar_ids, self.coeffs):
            list_fields.append(dvid)
            list_fields.append(coeff)
        return list_fields

    def repr_fields(self):
        cp_max = set_blank_if_default(self.cp_max, 1e20)
        c0 = set_blank_if_default(self.c0, 0.)
        list_fields = ['DVCREL1', self.oid, self.Type, self.Eid(),
                       self.cp_name, self.cp_min, cp_max, c0, None]
        for (dvid, coeff) in zip(self.desvar_ids, self.coeffs):
            list_fields.append(dvid)
            list_fields.append(coeff)
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        else:
            return self.comment + print_card_16(card)


class DVCREL2(OptConstraint):
    type = 'DVCREL2'

    allowed_elements = [
        'CQUAD4', 'CTRIA3', 'CBAR', 'CBEAM', 'CELAS1', 'CBUSH',
        'CDAMP2',
    ]
    #allowed_masses = ['CONM2', 'CMASS2', 'CMASS4']
    #allowed_properties_mass = ['PMASS']
    def __init__(self, oid, Type, eid, cp_name, deqation, dvids, labels,
                 cp_min=None, cp_max=1e20, validate=False, comment=''):
        """
        +----------+--------+--------+-------+------------+-------+-------+-------+-------+
        |    1     |    2   |   3    |   4   |      5     |   6   |   7   |   8   |   9   |
        +==========+========+========+=======+============+=======+=======+=======+=======+
        | DVCREL2  | ID     | TYPE   | EID   | CPNAME/FID | CPMIN | CPMAX | EQID  |       |
        +----------+--------+--------+-------+------------+-------+-------+-------+-------+
        |          | DESVAR | DVID1  | DVID2 |  DVID3     | DVID4 | DVID5 | DVID6 | DVID7 |
        +----------+--------+--------+-------+------------+-------+-------+-------+-------+
        |          |        | DVID8  | etc.  |            |       |       |       |       |
        +----------+--------+--------+-------+------------+-------+-------+-------+-------+
        |          | DTABLE | LABL1  | LABL2 |  LABL3     | LABL4 | LABL5 | LABL6 | LABL7 |
        +----------+--------+--------+-------+------------+-------+-------+-------+-------+
        |          |        | LABL8  | etc.  |            |       |       |       |       |
        +----------+--------+--------+-------+------------+-------+-------+-------+-------+
        """
        if comment:
            self.comment = comment
        if dvids is None:
            dvids = []
        elif isinstance(dvids, integer_types):
            dvids = [dvids]

        if labels is None:
            labels = []
        elif isinstance(labels, str):
            labels = [labels]

        #: Unique identification number
        self.oid = oid

        #: Name of an element connectivity entry, such as CBAR, CQUAD4, etc.
        #: (Character)
        self.Type = Type

        #: Element Identification number. (Integer > 0)
        self.eid = eid

        #: Name of connectivity property, such as X1, X2, X3, ZOFFS, etc.
        #: (Character)
        self.cp_name = cp_name

        #: Minimum value allowed for this property. If CPNAME references a connectivity
        #: property that can only be positive, then the default value of CPMIN is 1.0E-15.
        #: Otherwise, it is -1.0E35. (Real)
        #: .. todo:: bad default (see DVCREL2)
        self.cp_min = cp_min

        #: Maximum value allowed for this property. (Real; Default = 1.0E20)
        self.cp_max = cp_max

        #: DEQATN entry identification number. (Integer > 0)
        self.dequation = deqation
        self.dvids = dvids
        self.labels = labels

        #assert len(coeffs) > 0, 'len(coeffs)=%s' % len(coeffs)
        #assert len(coeffs) == len(dvids), 'len(coeffs)=%s len(dvids)=%s' % (len(coeffs), len(dvids))
        validate_dvcrel(validate, Type, cp_name)

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a DVCREL2 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        oid = integer(card, 1, 'oid')
        Type = string(card, 2, 'Type')
        pid = integer(card, 3, 'pid')
        cp_name = integer_or_string(card, 4, 'cpName_FID')
        cp_min = double_or_blank(card, 5, 'cp_min')
        cp_max = double_or_blank(card, 6, 'cp_max', 1e20)
        dequation = integer_or_blank(card, 7, 'dequation') #: .. todo:: or blank?

        fields = [interpret_value(field) for field in card[9:]]
        ioffset = 9
        iend = len(fields) + ioffset

        try:
            idesvar = fields.index('DESVAR') + ioffset
        except ValueError:
            idesvar = None

        try:
            idtable = fields.index('DTABLE') + ioffset
            #iDesMax  = idtable # the index to start parsing DESVAR
            ides_stop = idtable  # the index to stop  parsing DESVAR
        except ValueError:
            idtable = None
            ides_stop = iend

        dvids = []
        if idesvar:
            n = 1
            for i in range(10, ides_stop):
                dvid_name = 'DVID' + str(n)
                dvid = integer_or_blank(card, i, dvid_name)
                #print("%s = %s" % (dvid_name, dvid))
                if dvid:
                    assert dvid is not None
                    assert dvid is not 'DESVAR'
                    dvids.append(dvid)
                    n += 1

        labels = []
        if idtable:
            n = 1
            for i in range(idtable + 1, iend):
                label_name = 'Label' + str(n)
                label = string(card, i, label_name)
                #print("%s = %s" % (label_name, label))
                if label:
                    assert label is not 'DTABLE'
                    labels.append(label)
        return DVCREL2(oid, Type, pid, cp_name, dequation, dvids, labels,
                       cp_min, cp_max, comment=comment)

    def OptID(self):
        return self.oid

    def Eid(self):
        if isinstance(self.eid, integer_types):
            return self.eid
        elif self.Type in self.allowed_elements:
            eid = self.eid_ref.eid
        #elif self.Type in self.allowed_masses:
            #pid = self.pid_ref.eid
        #elif self.Type in self.allowed_properties_mass:
            #pid = self.pid_ref.pid
        else:
            raise NotImplementedError('Type=%r is not supported' % self.Type)
        return eid

    def DEquation(self):
        if isinstance(self.dequation, integer_types):
            return self.dequation
        elif isinstance(self.dequation, string_types):
            return self.dequation
        else:
            #print('dequation=%r; type=%s' % (self.dequation, type(self.dequation)))
            return self.dequation_ref.equation_id

    def calculate(self, op2_model, subcase_id):
        """
        this should really make a call the the DEQATN;
        see the PBEAM for an example of get/set_opt_value
        """
        try:
            get = self.pid_ref.get_optimization_value(self.cp_name)
            out = self.pid_ref.set_optimization_value(self.cp_name, get)
        except:
            print('DVCREL2 calculate : %s[%r] = ???' % (self.Type, self.cp_name))
            raise

        argsi = []
        if self.dvids:
            for desvar in self.dvids: # DESVARS
                arg = desvar.calculate(op2_model, subcase_id)
                argsi.append(arg)
        if self.labels:
            for label in self.labels: # DTABLE
                arg = self.dtable[label]
                argsi.append(arg)
        #op2_model.log.info('DVPREL2; args = %s' % argsi)

        #op2_model.log.info('dvids  =', self.dvids)
        #op2_model.log.info('labels =', self.labels)
        #op2_model.log.info('%s[%r] = %s' % (self.Type, self.cp_name, out))
        out = self.func(*argsi)
        op2_model.log.info('  deqatn out = %s' % out)
        return out

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        .. todo:: add support for DEQATN cards to finish DVPREL2 xref
        """
        msg = ', which is required by DVCREL2 oid=%r' % self.oid
        #if self.Type in self.allowed_elements:
            #self.pid = model.Element(self.pid, msg=msg)
        #elif self.Type in self.allowed_masses:
            #self.pid = model.masses[self.pid]
        #elif self.Type in self.allowed_properties_mass:
            #self.pid = model.properties_mass[self.pid]
        #else:
            #raise NotImplementedError('Type=%r is not supported' % self.Type)
        self.dequation = model.DEQATN(self.dequation)

        self.eid = model.Element(self.eid, msg=msg)
        self.eid_ref = self.eid
        #self.dvids = [model.Desvar(dvid, msg) for dvid in self.dvids]

        self.eid_ref = self.eid
        self.dequation_ref = self.dequation
        assert self.eid_ref.type in self.allowed_elements, self.eid.type

    def uncross_reference(self):
        self.eid = self.Eid()
        self.dequation = self.DEquation()
        del self.eid_ref, self.dequation_ref

    def raw_fields(self):
        list_fields = ['DVCREL2', self.oid, self.Type, self.Eid(),
                       self.cp_name, self.cp_min, self.cp_max, self.DEquation(), None]
        if self.dvids:
            fields2 = ['DESVAR'] + self.dvids
            list_fields += build_table_lines(fields2, nstart=1, nend=0)
        if self.labels:
            fields2 = ['DTABLE'] + self.labels
            list_fields += build_table_lines(fields2, nstart=1, nend=0)
        return list_fields

    def repr_fields(self):
        """
        .. todo:: finish repr_fields for DVCREL2
        """
        return self.raw_fields()

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)


class DVMREL1(OptConstraint):  # similar to DVPREL1
    type = 'DVMREL1'

    def __init__(self, oid, mat_type, mid, mp_name, dvids, coeffs,
                 mp_min=None, mp_max=1e20, c0=0., validate=False, comment=''):
        """
        Design Variable to Material Relation
        Defines the relation between a material property and design variables.

        +---------+-------+-------+-------+--------+-------+-------+--------+
        |    1    |   2   |   3   |   4   |    5   |   6   |   7   |    8   |
        +=========+=======+=======+=======+========+=======+=======+========+
        | DVMREL1 |  ID   | TYPE  |  MID  | MPNAME | MPMIN | MPMAX |   C0   |
        +---------+-------+-------+-------+--------+-------+-------+--------+
        |         | DVID1 | COEF1 | DVID2 | COEF2  | DVID3 | COEF3 | -etc.- |
        +---------+-------+-------+-------+--------+-------+-------+--------+
        """
        if comment:
            self.comment = comment
        if isinstance(dvids, integer_types):
            dvids = [dvids]
        if isinstance(coeffs, float):
            coeffs = [coeffs]

        self.oid = oid
        self.mat_type = mat_type
        self.mid = mid
        self.mp_name = mp_name
        self.mp_max = mp_max
        self.mp_min = mp_min
        self.c0 = c0
        self.dvids = dvids
        self.coeffs = coeffs

        validate_dvmrel(validate, mat_type, mp_name)

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a DVMREL1 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        oid = integer(card, 1, 'oid')
        mat_type = string(card, 2, 'mat_type')
        mid = integer(card, 3, 'mid')
        mp_name = string(card, 4, 'mpName')
        #if self.mp_name in ['E', 'RHO', 'NU']:  positive values
            #self.mp_min = double_or_blank(card, 5, 'mpMin', 1e-15)
        #else: # negative
            #self.mp_min = double_or_blank(card, 5, 'mpMin', -1e-35)
        mp_min = double_or_blank(card, 5, 'mpMin')  #: .. todo:: bad default
        mp_max = double_or_blank(card, 6, 'mpMax', 1e20)
        c0 = double_or_blank(card, 7, 'c0', 0.0)

        dvids = []
        coeffs = []
        end_fields = [interpret_value(field) for field in card[9:]]
        #print("end_fields = ",end_fields)
        nfields = len(end_fields) - 1
        if nfields % 2 == 1:
            end_fields.append(None)
            nfields += 1

        i = 0
        for i in range(0, nfields, 2):
            dvids.append(end_fields[i])
            coeffs.append(end_fields[i + 1])
        if nfields % 2 == 1:
            print(card)
            print("dvids = %s" % (dvids))
            print("coeffs = %s" % (coeffs))
            raise RuntimeError('invalid DVMREL1...')
        return DVMREL1(oid, mat_type, mid, mp_name, dvids, coeffs,
                       mp_min=mp_min, mp_max=mp_max, c0=c0, comment=comment)

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        self.mid = model.Material(self.mid)
        self.mid_ref = self.mid

    def uncross_reference(self):
        self.mid = self.Mid()
        del self.mid_ref

    def _verify(self, xref):
        """
        Verifies all methods for this object work

        Parameters
        ----------
        xref : bool
            has this model been cross referenced
        """
        pass

    def OptID(self):
        return self.oid

    def Mid(self):
        if isinstance(self.mid, integer_types):
            return self.mid
        return self.mid_ref.mid

    def raw_fields(self):
        list_fields = ['DVMREL1', self.oid, self.mat_type, self.Mid(),
                       self.mp_name, self.mp_min, self.mp_max, self.c0, None]
        for (dvid, coeff) in zip(self.dvids, self.coeffs):
            list_fields.append(dvid)
            list_fields.append(coeff)
        return list_fields

    def repr_fields(self):
        mp_max = set_blank_if_default(self.mp_max, 1e20)
        c0 = set_blank_if_default(self.c0, 0.)
        list_fields = ['DVMREL1', self.oid, self.mat_type, self.Mid(),
                       self.mp_name, self.mp_min, mp_max, c0, None]
        for (dvid, coeff) in zip(self.dvids, self.coeffs):
            list_fields.append(dvid)
            list_fields.append(coeff)
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        if is_double:
            return self.comment + print_card_double(card)
        return self.comment + print_card_16(card)


class DVMREL2(OptConstraint):
    type = 'DVMREL2'

    allowed_materials = ['MAT1', 'MAT2']
    def __init__(self, oid, mat_type, mid, mp_name, deqation, dvids, labels,
                 mp_min=None, mp_max=1e20, validate=False, comment=''):
        """
        +---------+--------+--------+-------+---------+-------+-------+-------+-------+
        |    1    |    2   |   3    |   4   |     5   |   6   |   7   |   8   |   9   |
        +=========+========+========+=======+=========+=======+=======+=======+=======+
        | DVMREL2 |   ID   | TYPE   |  MID  | MPNAME  | MPMIN | MPMAX | EQID  |       |
        +---------+--------+--------+-------+---------+-------+-------+-------+-------+
        |         | DESVAR | DVID1  | DVID2 | DVID3   | DVID4 | DVID5 | DVID6 | DVID7 |
        +---------+--------+--------+-------+---------+-------+-------+-------+-------+
        |         | DVID8  | -etc.- |       |         |       |       |       |       |
        +---------+--------+--------+-------+---------+-------+-------+-------+-------+
        |         | DTABLE | LABL1  | LABL2 | LABL3   | LABL4 | LABL5 | LABL6 | LABL7 |
        +---------+--------+--------+-------+---------+-------+-------+-------+-------+
        |         | LABL8  | -etc.- |       |         |       |       |       |       |
        +---------+--------+--------+-------+---------+-------+-------+-------+-------+
        """
        if comment:
            self.comment = comment
        if dvids is None:
            dvids = []
        elif isinstance(dvids, integer_types):
            dvids = [dvids]

        if labels is None:
            labels = []
        elif isinstance(labels, string_types):
            labels = [labels]

        #: Unique identification number
        self.oid = oid

        #: Name of a material entry, such as MAT1, MAT2, etc
        self.mat_type = mat_type

        #: Property entry identification number
        self.mid = mid

        #: Property name, such as 'E', 'RHO'
        #: (Character)
        self.mp_name = mp_name

        #: Minimum value allowed for this property. If MPNAME references a material
        #: property that can only be positive, then the default value for MPMIN is 1.0E-15.
        #: Otherwise, it is -1.0E35. (Real)
        self.mp_min = mp_min

        #: Maximum value allowed for this property. (Real; Default = 1.0E20)
        self.mp_max = mp_max
        #: DEQATN entry identification number. (Integer > 0)
        self.dequation = deqation
        self.dvids = dvids
        self.labels = labels
        validate_dvmrel(validate, mat_type, mp_name)

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a DVMREL2 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        oid = integer(card, 1, 'oid')
        mat_type = string(card, 2, 'mat_type')
        mid = integer(card, 3, 'mid')
        mp_name = integer_or_string(card, 4, 'mp_name')
        mp_min = double_or_blank(card, 5, 'mp_min')
        mp_max = double_or_blank(card, 6, 'mp_max', 1e20)
        dequation = integer_or_blank(card, 7, 'dequation') #: .. todo:: or blank?

        fields = [interpret_value(field) for field in card[9:]]
        ioffset = 9
        iend = len(fields) + ioffset

        try:
            idesvar = fields.index('DESVAR') + ioffset
        except ValueError:
            idesvar = None

        try:
            idtable = fields.index('DTABLE') + ioffset
            #iDesMax  = idtable # the index to start parsing DESVAR
            ides_stop = idtable  # the index to stop  parsing DESVAR
        except ValueError:
            idtable = None
            ides_stop = iend

        dvids = []
        if idesvar:
            n = 1
            for i in range(10, ides_stop):
                dvid_name = 'DVID' + str(n)
                dvid = integer_or_blank(card, i, dvid_name)
                #print("%s = %s" % (dvid_name, dvid))
                if dvid:
                    assert dvid is not None
                    assert dvid is not 'DESVAR'
                    dvids.append(dvid)
                    n += 1

        labels = []
        if idtable:
            n = 1
            for i in range(idtable + 1, iend):
                label_name = 'Label' + str(n)
                label = string(card, i, label_name)
                #print("%s = %s" % (label_name, label))
                if label:
                    assert label is not 'DTABLE'
                    labels.append(label)
        return DVMREL2(oid, mat_type, mid, mp_name, dequation, dvids, labels,
                       mp_min=mp_min, mp_max=mp_max, comment=comment)

    def OptID(self):
        return self.oid

    def Mid(self):
        if isinstance(self.mid, integer_types):
            return self.mid
        #if self.mat_type in self.allowed_properties:
            #pid = self.pid_ref.pid
        #elif self.mat_type in self.allowed_elements:
            #pid = self.pid_ref.eid
        #elif self.mat_type in self.allowed_masses:
            #pid = self.pid_ref.eid
        #elif self.mat_type in self.allowed_properties_mass:
            #pid = self.pid_ref.pid
        elif self.mat_type in self.allowed_materials:
            mid = self.mid_ref.mid
        else:
            raise NotImplementedError('mat_type=%r is not supported' % self.mat_type)
        return mid

    def DEquation(self):
        if isinstance(self.dequation, integer_types):
            return self.dequation
        return self.dequation_ref.equation_id

    def calculate(self, op2_model, subcase_id):
        """
        this should really make a call the the DEQATN;
        see the PBEAM for an example of get/set_opt_value
        """
        try:
            get = self.mid_ref.get_optimization_value(self.mp_name)
            out = self.mid_ref.set_optimization_value(self.mp_name, get)
        except:
            print('DVMREL2 calculate : %s[%r] = ???' % (self.mat_type, self.mp_name))
            raise

        if self.dvids:
            for desvar in self.dvids: # DESVARS
                arg = desvar.calculate(op2_model, subcase_id)
                argsi.append(arg)
        if self.labels:
            for label in self.labels: # DTABLE
                arg = self.dtable[label]
                argsi.append(arg)
        #op2_model.log.info('DVMREL2; args = %s' % argsi)

        #op2_model.log.info('dvids  =', self.dvids)
        #op2_model.log.info('labels =', self.labels)
        #op2_model.log.info('%s[%r] = %s' % (self.mat_type, self.cp_name, out))
        out = self.func(*argsi)
        op2_model.log.info('  deqatn out = %s' % out)
        return out
        #raise NotImplementedError('\n' + str(self))

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        .. todo:: add support for DEQATN cards to finish DVMREL2 xref
        """
        msg = ', which is required by DVMREL2 name=%r' % self.type
        if self.mat_type in self.allowed_materials:
            self.mid = model.Material(self.mid, msg=msg)
        #elif self.mat_type in self.allowed_elements:
            #self.mid = model.Element(self.mid, msg=msg)
        #elif self.mat_type in self.allowed_masses:
            #self.mid = model.masses[self.mid]
        #elif self.mat_type in self.allowed_properties_mass:
            #self.mid = model.properties_mass[self.mid]
        else:
            raise NotImplementedError('mat_type=%r is not supported' % self.mat_type)
        self.dequation = model.DEQATN(self.dequation)

        self.mid_ref = self.mid
        self.dequation_ref = self.dequation
        #assert self.pid_ref.type not in ['PBEND', 'PBARL', 'PBEAML'], self.pid

    def uncross_reference(self):
        self.mid = self.Mid()
        self.dequation = self.DEquation()
        del self.mid_ref, self.dequation_ref

    #def OptValue(self):  #: .. todo:: not implemented
        #self.pid_ref.OptValue(self.mp_name)

    def raw_fields(self):
        list_fields = ['DVMREL2', self.oid, self.mat_type, self.Mid(),
                       self.mp_name, self.mp_min, self.mp_max, self.DEquation(), None]
        if self.dvids:
            fields2 = ['DESVAR'] + self.dvids
            list_fields += build_table_lines(fields2, nstart=1, nend=0)
        if self.labels:
            fields2 = ['DTABLE'] + self.labels
            list_fields += build_table_lines(fields2, nstart=1, nend=0)
        return list_fields

    def repr_fields(self):
        """
        .. todo:: finish repr_fields for DVMREL2
        """
        return self.raw_fields()

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)


def break_word_by_trailing_integer(pname_fid):
    nums = []
    for i, letter in enumerate(reversed(pname_fid)):
        if letter.isdigit():
            nums.append(letter)
        else:
            break
        num = ''.join(nums[::-1])
        word = pname_fid[:-i-1]
        assert len(word)+len(num) == len(pname_fid), 'word=%r num=%r pname_fid=%r' % (word, num, pname_fid)
        return word, num


class DVPREL1(OptConstraint):  # similar to DVMREL1
    """
    +---------+--------+--------+--------+-----------+-------+--------+-----+---+
    |   1     |    2   |   3    |    4   |     5     |   6   |   7    |  8  | 9 |
    +=========+========+========+========+===========+=======+========+=====+===+
    | DVPREL1 |   ID   |  TYPE  |  PID   | PNAME/FID | PMIN  |  PMAX  |  C0 |   |
    +---------+--------+--------+--------+-----------+-------+--------+-----+---+
    |         | DVID1  | COEF1  | DVID2  |   COEF2   | DVID3 | -etc.- |     |   |
    +---------+--------+--------+--------+-----------+-------+--------+-----+---+

    +---------+--------+--------+--------+-----+
    | DVPREL1 | 200000 | PCOMP  | 2000   |  T2 |
    +---------+--------+--------+--------+-----+
    |         | 200000 |   1.0  |        |     |
    +---------+--------+--------+--------+-----+
    """
    type = 'DVPREL1'

    allowed_properties = [
        'PELAS',
        'PROD', 'PTUBE',
        'PBAR', 'PBARL',
        'PBEAM', 'PBEAML',
        'PSHEAR',
        'PSHELL', 'PCOMP', 'PCOMPG',
        'PBUSH', 'PBUSH1D',
        'PGAP', 'PVISC',
        'PDAMP', 'PWELD',
        'PBMSECT',
    ]
    allowed_elements = [
        'CELAS2', 'CBAR', 'CBEAM',
        'CQUAD4',
        'CBUSH', 'CDAMP2',
    ]
    allowed_masses = ['CONM2', 'CMASS2', 'CMASS4']
    allowed_properties_mass = ['PMASS']
    def __init__(self, oid, prop_type, pid, pname_fid, dvids, coeffs,
                 p_min=None, p_max=1e20, c0=0.0, validate=False, comment=''):
        """
        Creates a DVPREL1 card

        Parameters
        ----------
        oid : int
            optimization id
        prop_type : str
            property card name (e.g., PSHELL)
        pid : int
            property id
        pname_fid : str/int
            optimization parameter as a pname (property name; T) or
            field number (fid)
        dvids : List[int]
            DESVAR ids
        coeffs : List[float]
            scale factors for DESVAR ids
        p_min : float; default=None
            minimum property value
        p_max : float; default=1e20
            maximum property value
        c0 : float; default=0.
            offset factor for the variable
        validate : bool; default=False
            should the variable be validated
        comment : str; default=''
            a comment for the card
        """
        if comment:
            self.comment = comment
        if isinstance(dvids, integer_types):
            dvids = [dvids]
        if isinstance(coeffs, float):
            coeffs = [coeffs]

        self.oid = oid

        # property type (e.g. PSHELL/PCOMP)
        self.prop_type = prop_type

        # property id
        self.pid = pid

        # the field type (e.g. 'T' on a PSHELL or the field id)
        self.pname_fid = pname_fid

        # min value for 'T'
        self.p_min = p_min

        # max value for 'T'
        self.p_max = p_max

        # offset coefficient
        self.c0 = c0

        # DESVAR id
        self.dvids = dvids

        # scale factor for DESVAR
        self.coeffs = coeffs

        if len(coeffs) == 0:
            msg = 'len(coeffs)=%s len(dvids)=%s\n' % (len(coeffs), len(dvids))
            msg += "We've added a coeff=1 and desvar_id=1 in order to look at the crashing card\n"
            self.coeffs = [1]
            self.dvids = [1]
            msg += str(self)
            raise RuntimeError(msg)
        assert len(coeffs) == len(dvids), 'len(coeffs)=%s len(dvids)=%s' % (len(coeffs), len(dvids))

        pname_fid = validate_dvprel(prop_type, pname_fid, validate)
        self.pname_fid = pname_fid

    def validate(self):
        self.pname_fid = validate_dvprel(self.prop_type, self.pname_fid, validate=True)

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a DVPREL1 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        oid = integer(card, 1, 'oid')
        prop_type = string(card, 2, 'prop_type')
        pid = integer(card, 3, 'pid')
        pname_fid = integer_or_string(card, 4, 'pName_FID')

        #: Minimum value allowed for this property.
        #: .. todo:: bad default (see DVMREL1)
        p_min = double_or_blank(card, 5, 'pMin', None)
        p_max = double_or_blank(card, 6, 'pMax', 1e20)
        c0 = double_or_blank(card, 7, 'c0', 0.0)

        dvids = []
        coeffs = []
        end_fields = [interpret_value(field) for field in card[9:]]

        nfields = len(end_fields) - 1
        if nfields % 2 == 1:
            end_fields.append(None)
            nfields += 1
        i = 0
        for i in range(0, nfields, 2):
            dvids.append(end_fields[i])
            coeffs.append(end_fields[i + 1])
        if nfields % 2 == 1:
            print(card)
            print("dvids = %s" % (dvids))
            print("coeffs = %s" % (coeffs))
            raise RuntimeError('invalid DVPREL1...')
        return DVPREL1(oid, prop_type, pid, pname_fid, dvids, coeffs,
                       p_min=p_min, p_max=p_max, c0=c0,
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

    def OptID(self):
        return self.oid

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ', which is required by DVPREL1 name=%r' % self.type
        if self.prop_type in self.allowed_properties:
            self.pid = model.Property(self.pid, msg=msg)
        #elif self.prop_type in self.allowed_elements:
            #self.pid = model.Element(self.pid, msg=msg)
        #elif self.prop_type in self.allowed_masses:
            #self.pid = model.masses[self.pid]
        elif self.prop_type in self.allowed_properties_mass:
            self.pid = model.properties_mass[self.pid]
        elif self.prop_type == 'PBUSHT':
            self.pid = model.pbusht[self.pid]
        elif self.prop_type == 'PELAST':
            self.pid = model.pelast[self.pid]
        else:
            raise NotImplementedError('prop_type=%r is not supported' % self.prop_type)
        self.pid_ref = self.pid
        self.dvids = [model.Desvar(dvid, msg) for dvid in self.dvids]

    def uncross_reference(self):
        self.pid = self.Pid()
        del self.pid_ref

    def calculate(self, op2_model, subcase_id):
        raise NotImplementedError('\n' + str(self))

    def Pid(self):
        if isinstance(self.pid, integer_types):
            return self.pid
        if self.prop_type in self.allowed_properties:
            pid = self.pid_ref.pid
        #elif self.prop_type in self.allowed_elements:
            #pid = self.pid_ref.eid
        #elif self.prop_type in self.allowed_masses:
            #pid = self.pid_ref.eid
        elif self.prop_type in self.allowed_properties_mass:
            pid = self.pid_ref.pid
        elif self.prop_type in ['PBUSHT', 'PELAST']:
            pid = self.pid_ref.pid
        else:
            raise NotImplementedError('prop_type=%r is not supported' % self.prop_type)
        return pid

    @property
    def desvar_ids(self):
        if isinstance(self.dvids[0], integer_types):
            return self.dvids
        return [desvar.desvar_id for desvar in self.dvids]

    def raw_fields(self):
        list_fields = ['DVPREL1', self.oid, self.prop_type, self.Pid(),
                       self.pname_fid, self.p_min, self.p_max, self.c0, None]
        for (dvid, coeff) in zip(self.desvar_ids, self.coeffs):
            list_fields.append(dvid)
            list_fields.append(coeff)
        return list_fields

    def repr_fields(self):
        p_max = set_blank_if_default(self.p_max, 1e20)
        c0 = set_blank_if_default(self.c0, 0.)
        list_fields = ['DVPREL1', self.oid, self.prop_type, self.Pid(),
                       self.pname_fid, self.p_min, p_max, c0, None]
        for (dvid, coeff) in zip(self.desvar_ids, self.coeffs):
            list_fields.append(dvid)
            list_fields.append(coeff)
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        else:
            return self.comment + print_card_16(card)


class DVPREL2(OptConstraint):
    """
    +----------+--------+--------+-------+-----------+-------+-------+-------+-------+
    |    1     |    2   |   3    |   4   |     5     |   6   |   7   |   8   |   9   |
    +==========+========+========+=======+===========+=======+=======+=======+=======+
    | DVPREL2  | ID     | TYPE   | PID   | PNAME/FID | PMIN  | PMAX  | EQID  |       |
    +----------+--------+--------+-------+-----------+-------+-------+-------+-------+
    |          | DESVAR | DVID1  | DVID2 |   DVID3   | DVID4 | DVID5 | DVID6 | DVID7 |
    +----------+--------+--------+-------+-----------+-------+-------+-------+-------+
    |          |        | DVID8  | etc.  |           |       |       |       |       |
    +----------+--------+--------+-------+-----------+-------+-------+-------+-------+
    |          | DTABLE | LABL1  | LABL2 |   LABL3   | LABL4 | LABL5 | LABL6 | LABL7 |
    +----------+--------+--------+-------+-----------+-------+-------+-------+-------+
    |          |        | LABL8  | etc.  |           |       |       |       |       |
    +----------+--------+--------+-------+-----------+-------+-------+-------+-------+
    """
    type = 'DVPREL2'

    allowed_properties = [
        'PELAS',
        'PROD', 'PTUBE',
        'PBAR', 'PBARL',
        'PBEAM', 'PBEAML',
        'PSHELL', 'PCOMP', 'PCOMPG',
        'PBUSH', 'PBUSH1D',
        'PGAP', 'PVISC',
        'PDAMP', 'PWELD',
    ]
    allowed_elements = [
        'CELAS2', 'CBAR', 'CBEAM',
        'CQUAD4',
        'CBUSH', 'CDAMP2',
    ]
    allowed_masses = ['CONM2', 'CMASS2', 'CMASS4']
    allowed_properties_mass = ['PMASS']
    def __init__(self, oid, Type, pid, pname_fid, deqation,
                 dvids=None, labels=None, p_min=None, p_max=1e20,
                 validate=False, comment=''):
        """
        Creates a DVPREL2 card

        Parameters
        ----------
        oid : int
            optimization id
        prop_type : str
            property card name (e.g., PSHELL)
        pid : int
            property id
        pname_fid : str/int
            optimization parameter as a pname (property name; T) or field number (fid)
        deqation : int
            DEQATN id
        dvids : List[int]; default=None
            DESVAR ids
        labels : List[str]; default=None
            DTABLE names
        p_min : float; default=None
            minimum property value
        p_max : float; default=1e20
            maximum property value
        validate : bool; default=False
            should the variable be validated
        comment : str; default=''
            a comment for the card

        .. note:: either dvids or labels is required
        """
        if comment:
            self.comment = comment
        if dvids is None:
            dvids = []
        elif isinstance(dvids, integer_types):
            dvids = [dvids]

        if labels is None:
            labels = []
        elif isinstance(labels, string_types):
            labels = [labels]

        #: Unique identification number
        self.oid = oid

        #: Name of a property entry, such as PBAR, PBEAM, etc
        self.Type = Type

        #: Property entry identification number
        self.pid = pid

        #: Property name, such as 'T', 'A', or field position of the property
        #: entry, or word position in the element property table of the
        #: analysis model. Property names that begin with an integer such as
        #: 12I/T**3 may only be referred to by field position.
        #: (Character or Integer 0)
        self.pname_fid = pname_fid

        #: Minimum value allowed for this property. If FID references a stress
        #: recovery location field, then the default value for PMIN is -1.0+35.
        #: PMIN must be explicitly set to a negative number for properties that
        #: may be less than zero (for example, field ZO on the PCOMP entry).
        #: (Real; Default = 1.E-15)
        #: .. todo:: bad default (see DVMREL1)
        self.p_min = p_min
        #: Maximum value allowed for this property. (Real; Default = 1.0E20)
        self.p_max = p_max
        #: DEQATN entry identification number. (Integer > 0)
        self.dequation = deqation
        self.dvids = dvids
        self.labels = labels

        pname_fid = validate_dvprel(Type, pname_fid, validate)
        self.pname_fid = pname_fid
        #print(self)

    def validate(self):
        assert len(self.dvids) > 0 or len(self.labels) > 0
        self.pname_fid = validate_dvprel(self.Type, self.pname_fid, validate=True)

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a DVPREL2 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        oid = integer(card, 1, 'oid')
        Type = string(card, 2, 'Type')
        pid = integer(card, 3, 'pid')
        pname_fid = integer_or_string(card, 4, 'pName_FID')
        p_min = double_or_blank(card, 5, 'pMin', None)
        p_max = double_or_blank(card, 6, 'pMax', 1e20)
        dequation = integer_or_blank(card, 7, 'dequation') #: .. todo:: or blank?

        fields = [interpret_value(field) for field in card[9:]]
        ioffset = 9
        iend = len(fields) + ioffset

        #F:\work\pyNastran\examples\femap_examples\Support\nast\tpl\d200m20.dat
        #params = parse_table_fields(fields)
        #print(params)

        try:
            idesvar = fields.index('DESVAR') + ioffset
        except ValueError:
            idesvar = None

        try:
            idtable = fields.index('DTABLE') + ioffset
            #iDesMax  = idtable # the index to start parsing DESVAR
            ides_stop = idtable  # the index to stop  parsing DESVAR
        except ValueError:
            idtable = None
            ides_stop = iend

        dvids = []
        if idesvar:
            n = 1
            for i in range(10, ides_stop):
                dvid_name = 'DVID' + str(n)
                dvid = integer_or_blank(card, i, dvid_name)
                #print("%s = %s" % (dvid_name, dvid))
                if dvid:
                    assert dvid is not None
                    assert dvid is not 'DESVAR'
                    dvids.append(dvid)
                    n += 1

        labels = []
        if idtable:
            n = 1
            for i in range(idtable + 1, iend):
                label_name = 'Label' + str(n)
                label = string(card, i, label_name)
                #print("%s = %s" % (label_name, label))
                if label:
                    assert label is not 'DTABLE'
                    labels.append(label)
        return DVPREL2(oid, Type, pid, pname_fid, dequation, dvids, labels,
                       p_min=p_min, p_max=p_max, comment=comment)

    def OptID(self):
        return self.oid

    def Pid(self):
        if isinstance(self.pid, integer_types):
            return self.pid
        if self.Type in self.allowed_properties:
            pid = self.pid_ref.pid
        elif self.Type in self.allowed_elements:
            pid = self.pid_ref.eid
        elif self.Type in self.allowed_masses:
            pid = self.pid_ref.eid
        elif self.Type in self.allowed_properties_mass:
            pid = self.pid_ref.pid
        else:
            raise NotImplementedError('Type=%r is not supported' % self.Type)
        return pid

    def DEquation(self):
        if isinstance(self.dequation, integer_types):
            return self.dequation
        return self.dequation_ref.equation_id

    def calculate(self, op2_model, subcase_id):
        """
        this should really make a call the the DEQATN;
        see the PBEAM for an example of get/set_opt_value
        """
        try:
            get = self.pid_ref.get_optimization_value(self.pname_fid)
            out = self.pid_ref.set_optimization_value(self.pname_fid, get)
        except:
            print('DVPREL2 calculate : %s[%r] = ???' % (self.Type, self.pname_fid))
            raise

        if self.dvids:
            for desvar in self.dvids: # DESVARS
                arg = desvar.calculate(op2_model, subcase_id)
                argsi.append(arg)
        if self.labels:
            for label in self.labels: # DTABLE
                arg = self.dtable[label]
                argsi.append(arg)
        #op2_model.log.info('DVPREL2; args = %s' % argsi)

        #op2_model.log.info('dvids  =', self.dvids)
        #op2_model.log.info('labels =', self.labels)
        #op2_model.log.info('%s[%r] = %s' % (self.Type, self.pname_fid, out))
        out = self.func(*argsi)
        op2_model.log.info('  deqatn out = %s' % out)
        return out

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        .. todo:: add support for DEQATN cards to finish DVPREL2 xref
        """
        msg = ', which is required by DVPREL2 name=%r' % self.type
        if self.Type in self.allowed_properties:
            self.pid = model.Property(self.pid, msg=msg)
        elif self.Type in self.allowed_elements:
            self.pid = model.Element(self.pid, msg=msg)
        elif self.Type in self.allowed_masses:
            self.pid = model.masses[self.pid]
        elif self.Type in self.allowed_properties_mass:
            self.pid = model.properties_mass[self.pid]
        else:
            raise NotImplementedError('Type=%r is not supported' % self.Type)
        self.dequation = model.DEQATN(self.dequation)

        self.pid_ref = self.pid
        self.dequation_ref = self.dequation
        assert self.pid_ref.type not in ['PBEND', 'PBARL', 'PBEAML'], self.pid

    def uncross_reference(self):
        self.pid = self.Pid()
        self.dequation = self.DEquation()
        del self.pid_ref, self.dequation_ref

    def _verify(self, xref):
        """
        Verifies all methods for this object work

        Parameters
        ----------
        xref : bool
            has this model been cross referenced
        """
        pass

    #def OptValue(self):  #: .. todo:: not implemented
        #self.pid_ref.OptValue(self.pname_fid)

    def raw_fields(self):
        list_fields = ['DVPREL2', self.oid, self.Type, self.Pid(),
                       self.pname_fid, self.p_min, self.p_max, self.DEquation(), None]
        if self.dvids:
            fields2 = ['DESVAR'] + self.dvids
            list_fields += build_table_lines(fields2, nstart=1, nend=0)
        if self.labels:
            fields2 = ['DTABLE'] + self.labels
            list_fields += build_table_lines(fields2, nstart=1, nend=0)
        return list_fields

    def repr_fields(self):
        """
        .. todo:: finish repr_fields for DVPREL2
        """
        return self.raw_fields()

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)

class DVGRID(OptConstraint):
    type = 'DVGRID'

    """
    +--------+------+-----+-----+-------+----+----+----+
    |    1   |   2  |  3  |  4  |   5   |  6 |  7 |  8 |
    +========+======+=====+=====+=======+====+====+====+
    | DVGRID | DVID | GID | CID | COEFF | N1 | N2 | N3 |
    +--------+------+-----+-----+-------+----+----+----+
    """
    def __init__(self, dvid, nid, dxyz, cid=0, coeff=1.0, comment=''):
        if comment:
            self.comment = comment
        self.dvid = dvid
        self.nid = nid
        self.cid = cid
        self.coeff = coeff
        self.dxyz = np.asarray(dxyz)

    def validate(self):
        if np.linalg.norm(self.dxyz) == 0.:
            msg = 'DVGRID Error; dvid=%s nid=%s norm(dxyz)=0.' % (self.dvid, self.nid)
            raise ValueError(msg)

    @staticmethod
    def add_card(card, comment=''):
        dvid = integer(card, 1, 'dvid')
        nid = integer(card, 2, 'nid')
        cid = integer_or_blank(card, 3, 'cid', 0)
        coeff = double_or_blank(card, 4, 'coeff', 1.0)
        dxyz = [
            double_or_blank(card, 5, 'n1', 0.),
            double_or_blank(card, 6, 'n2', 0.),
            double_or_blank(card, 7, 'n3', 0.),
        ]
        return DVGRID(dvid, nid, dxyz, cid=cid, coeff=coeff, comment=comment)

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        self.dvid_ref = model.desvars[self.dvid]
        self.nid_ref = model.Node(self.nid)
        self.cid_ref = model.Coord(self.cid)
        #self.dconstrs = [model.dconstrs[oid] for oid in self.dconstr_ids]
        #self.dconstrs_ref = [model.dconstrs[oid] for oid in self.dconstr_ids]

    def uncross_reference(self):
        return
        #self.dconstrs = self.dconstr_ids
        #del self.dconstrs_ref

    def raw_fields(self):
        list_fields = ['DVGRID', self.dvid, self.nid, self.cid, self.coeff] + list(self.dxyz)
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        if is_double:
            return self.comment + print_card_double(card)
        return self.comment + print_card_16(card)

    def _verify(self, xref=True):
        pass

def parse_table_fields(fields, allowed=None):
    """
    params = {
       (0, 'DRESP1') = [10, 20],
       (1, 'DESVAR') = [30],
       (2, 'DRESP1') = [40],
    }
    """
    params = {}
    if not fields:
        return params

    j = 0
    value_list = []
    key = None  # dummy key
    for (i, field) in enumerate(fields):
        if i % 8 == 0 and field is not None:
            if i > 0:
                assert len(value_list) > 0, 'key=%s values=%s' % (key, value_list)
                params[key] = value_list
                j += 1
            key = (j, field)
            value_list = []
            name = field
        elif field is not None:
            if name in ['DESVAR', 'DRESP1', 'DRESP2', 'DVCREL1', 'DVCREL2', 'DVMREL1', 'DVMREL2', 'DVPREL1', 'DVPREL2', 'DNODE']:
                #print('field=%s value=%r type=%r should be an integer...\ncard=%s' % (i+9, field, name, card))
                assert isinstance(field, integer_types), 'field=%i value=%r type=%s should be an integer...\ncard=%s' % (i+9, field, name, card)
            elif name == 'DTABLE':
                #print('field=%s value=%r type=%r should be an string...\ncard=%s' % (i+9, field, name, card))
                assert isinstance(field, string_types), 'field=%i value=%r type=%s should be an string...\ncard=%s' % (i+9, field, name, card)
            elif name == 'DFRFNC':
                pass
            else:
                raise NotImplementedError('name=%r\n%s' % (name, card))
            value_list.append(field)
    params[key] = value_list
    return params
