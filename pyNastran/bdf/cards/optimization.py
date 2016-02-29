# pylint: disable=C0103,R0902,R0904,R0914
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
from six import iteritems, string_types
from six.moves import zip, range
from numpy import searchsorted, array

from pyNastran.utils import integer_types
from pyNastran.bdf.field_writer_8 import set_blank_if_default
from pyNastran.bdf.cards.base_card import (BaseCard, expand_thru_by)
from pyNastran.bdf.cards.deqatn import fortran_to_python, fortran_to_python_short
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

class OptConstraint(BaseCard):
    def __init__(self):
        pass


class DCONSTR(OptConstraint):
    type = 'DCONSTR'
    def __init__(self, oid, rid, lid, uid, lowfq, highfq, comment=''):
        if comment:
            self._comment = comment
        self.oid = oid
        self.rid = rid
        self.lid = lid
        self.uid = uid
        self.lowfq = lowfq
        self.highfq = highfq

    @classmethod
    def add_card(cls, card, comment=''):
        oid = integer(card, 1, 'oid')
        rid = integer(card, 2, 'rid')
        lid = integer_double_or_blank(card, 3, 'lid', -1e20)
        uid = integer_double_or_blank(card, 4, 'uid', 1e20)
        lowfq = double_or_blank(card, 5, 'lowfq', 0.0)
        highfq = double_or_blank(card, 6, 'highfq', 1e20)
        assert len(card) <= 7, 'len(DCONSTR card) = %i' % len(card)
        return DCONSTR(oid, rid, lid, uid, lowfq, highfq, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        oid = data[0]
        rid = data[1]
        lid = data[2]
        uid = data[3]
        lowfq = data[4]
        highfq = data[5]
        return DCONSTR(oid, rid, lid, uid, lowfq, highfq, comment=comment)

    def Rid(self):
        if isinstance(self.rid, integer_types):
            return self.rid
        else:
            return self.rid_ref.oid

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
        self.rid = model.dresps[self.Rid()]
        self.rid_ref = self.rid
        if isinstance(self.lid, integer_types):
            self.lid = model.Table(self.lid)
            self.lid_ref = self.lid
        if isinstance(self.uid, integer_types):
            self.uid = model.Table(self.uid)
            self.uid_ref = self.uid

    def uncross_reference(self):
        self.rid = self.Rid()
        self.lid = self.Lid()
        self.uid = self.Uid()

        if isinstance(self.lid, integer_types):
            del self.lid
        if isinstance(self.uid, integer_types):
            del self.uid_ref
        del self.rid_ref

    def raw_fields(self):
        list_fields = ['DCONSTR', self.oid, self.Rid(), self.Lid(),
                       self.Uid(), self.lowfq, self.highfq]
        return list_fields

    def repr_fields(self):
        lid = set_blank_if_default(self.Lid(), -1e20)
        uid = set_blank_if_default(self.Uid(), 1e20)
        lowfq = set_blank_if_default(self.lowfq, 0.0)
        highfq = set_blank_if_default(self.highfq, 1e20)
        list_fields = ['DCONSTR', self.oid, self.Rid(), lid, uid, lowfq, highfq]
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
    def __init__(self, oid, label, xinit, xlb, xub, delx, ddval, comment=''):
        if comment:
            self._comment = comment
        self.oid = oid
        self.label = label
        self.xinit = xinit
        self.xlb = xlb
        self.xub = xub
        self.delx = delx
        self.ddval = ddval

    @classmethod
    def add_card(cls, card, comment=''):
        oid = integer(card, 1, 'oid')
        label = string(card, 2, 'label')
        xinit = double(card, 3, 'xinit')
        xlb = double_or_blank(card, 4, 'xlb', -1e20)
        xub = double_or_blank(card, 5, 'xub', 1e20)
        delx = double_or_blank(card, 6, 'delx', 1e20)
        ddval = integer_or_blank(card, 7, 'ddval')
        assert len(card) <= 8, 'len(DESVAR card) = %i' % len(card)
        return DESVAR(oid, label, xinit, xlb, xub, delx, ddval, comment=comment)

    def OptID(self):
        return self.oid

    def raw_fields(self):
        list_fields = ['DESVAR', self.oid, self.label, self.xinit, self.xlb,
                       self.xub, self.delx, self.ddval]
        return list_fields

    def repr_fields(self):
        xlb = set_blank_if_default(self.xlb, -1e20)
        xub = set_blank_if_default(self.xub, 1e20)
        delx = set_blank_if_default(self.delx, 1e20)
        list_fields = ['DESVAR', self.oid, self.label, self.xinit, xlb,
                       xub, delx, self.ddval]
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)


class DDVAL(OptConstraint):
    type = 'DDVAL'

    def __init__(self, oid, ddvals, comment=''):
        if comment:
            self._comment = comment
        self.oid = oid
        self.ddvals = ddvals

    @classmethod
    def add_card(cls, card, comment=''):
        oid = integer(card, 1, 'oid')
        n = 1
        ddvals = []
        for i in range(2, len(card)):
            ddval = double_string_or_blank(card, i, 'DDVAL%s' % n)
            if ddval:
                ddvals.append(ddval)
        ddvals = expand_thru_by(ddvals)
        ddvals.sort()
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
            self._comment = comment
        self.params = params

    @classmethod
    def add_card(cls, card, comment=''):
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

    def __init__(self, oid, ddvid, c0, cmult, IDv, Ci, comment=''):
        """
        Multiple Design Variable Linking
        Relates one design variable to one or more other design variables.

        +-------+------+-------+--------+-------+------+----+------+----+
        | DLINK |  ID  | DDVID |   C0   | CMULT | IDV1 | C1 | IDV2 | C2 |
        +-------+------+-------+--------+-------+------+----+------+----+
        |       | IDV3 |   C3  | -etc.- |       |      |    |      |    |
        +-------+------+-------+--------+-------+------+----+------+----+
        """
        if comment:
            self._comment = comment
        self.oid = oid
        self.ddvid = ddvid
        self.c0 = c0
        self.cmult = cmult
        self.IDv = IDv
        self.Ci = Ci

    @classmethod
    def add_card(cls, card, comment=''):
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
        return DLINK(oid, ddvid, c0, cmult, idvs, Ci, comment=comment)

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


class DRESP1(OptConstraint):
    type = 'DRESP1'

    def __init__(self, oid, label, rtype, ptype, region,
                 atta, attb, atti, comment=''):
        """
        +--------+-----------+------------+-------+------+--------+-------+-----+-------+
        | DRESP1 |       1S1 |    CSTRAIN | PCOMP |      |        | 1     | 1   | 10000 |
        +--------+-----------+------------+-------+------+--------+-------+-----+-------+
        """
        if comment:
            self._comment = comment
        self.oid = oid
        self.label = label
        self.rtype = rtype
        self.ptype = ptype
        self.region = region
        self.atta = atta
        self.attb = attb
        self.atti = atti

    @classmethod
    def add_card(cls, card, comment=''):
        oid = integer(card, 1, 'oid')
        label = string(card, 2, 'label')
        rtype = string(card, 3, 'rtype')

        # elem, pbar, pshell, etc. (ELEM flag or Prop Name)
        ptype = integer_string_or_blank(card, 4, 'ptype')
        #if 1:
            ## incomplete
            #assert self.ptype in ['ELEM', 'PSHELL', 'PBAR', 'PROD', 'PCOMP',
                                  #'PSOLID', 'PELAS', 'PBARL', 'PBEAM',
                                  #'PBEAML', 'PSHEAR', 'PTUBE',
                                  #'PKNL',
                                  #None], 'DRESP1 ptype=%s' % self.ptype
        region = integer_or_blank(card, 5, 'region')
        atta = integer_double_string_or_blank(card, 6, 'atta')
        attb = integer_double_string_or_blank(card, 7, 'attb')

        atti = []
        for i in range(8, len(card)):
            attii = integer_double_string_or_blank(card, i, 'atti_%i' % (i + 1))
            atti.append(attii)
        return DRESP1(oid, label, rtype, ptype, region, atta, attb, atti,
                      comment=comment)

    def _verify(self, xref=True):
        pass

    def calculate(self, op2_model, subcase_id):
        rtype = self.rtype
        ptype = self.ptype
        if rtype == 'DISP' and ptype is None:
            msg = 'fields=%s\n' % (self.raw_fields())
            msg += 'rtype=%r ptype=%r region=%s A=%r B=%r\n' % (self.rtype, self.ptype, self.region,
                                                                self.atta, self.attb)
            assert isinstance(self.atta, int), self.atta
            assert self.attb is None, self.attb
            assert len(self.atti) == 1, self.atti
            #msg += str(self)
            #raise NotImplementedError(msg)
            comp = self.atta
            assert comp < 7, comp
            case = op2_model.displacements[subcase_id]
            nids = case.node_gridtype[:, 0]
            nidsi = array(
                [node if isinstance(node, integer_types) else node.nid
                 for node in self.atti], dtype='int32')
            inid = searchsorted(nids, nidsi)[0]
            itime = 0 #  TODO:  ???
            #print('atti = ', self.atti)
            atti = self.atti[0]
            out = case.data[itime, inid, comp - 1]
            print(' DISP[itime=%s, nid=%s, comp=%i] = %s' % (itime, str(nidsi), comp, out))
        elif ptype == 'ELEM' and 0:
            pass
        else:
            msg = 'fields=%s\n' % (self.raw_fields())
            msg += 'rtype=%r ptype=%r region=%s\n' % (self.rtype, self.ptype, self.region)
            msg += str(self)
            raise NotImplementedError(msg)
        return out

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
        msg = ', which is required by %s oid=%s' % (self.type, self.oid)
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
        if self.ptype in ['ELEM']:
            self.atti = model.Elements(self.atti, msg=msg)
            self.atti_ref = self.atti
        elif self.ptype in ['PSHELL', 'PBAR', 'PROD', 'PCOMP',
                            'PSOLID', 'PELAS', 'PBARL', 'PBEAM',
                            'PBEAML', 'PSHEAR', 'PTUBE',]:
            self.atti = model.Properties(self.atti, msg=msg)
            self.atti_ref = self.atti
        elif self.rtype in ['FRSTRE']:
            self.atti = model.Properties(self.atti, msg=msg)
            self.atti_ref = self.atti
        elif self.rtype in ['WEIGHT', 'FLUTTER', 'STABDER', 'CEIG', 'EIGN', 'FREQ']:
            pass
        elif self.rtype in ['DISP', 'FRDISP', 'TDISP',
                            'FRVELO', 'TVELO',
                            'FRACCL', 'PSDACCL']:
            self.atti = model.Nodes(self.atti, msg=msg)
            self.atti_ref = self.atti
        elif self.rtype in op2_results:
            pass
        else:
            msg = 'rtype=%r ptype=%r\n' % (self.rtype, self.ptype)
            msg += str(self)
            raise NotImplementedError(msg)

    def uncross_reference(self):
        self.atti = self.atti_values()
        if hasattr(self, 'atti_ref'):
            del self.atti_ref

    def atti_values(self):
        #return self.atti
        #if self.rtype in ['ELEM']:
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
        if self.ptype in ['PSHELL', 'PBAR', 'PROD', 'PCOMP',
                          'PSOLID', 'PELAS', 'PBARL', 'PBEAM',
                          'PBEAML', 'PSHEAR', 'PTUBE',
                          'FRSTRE']:
            data = [prop if isinstance(prop, integer_types) else prop.pid for prop in self.atti]
            for value in data:
                assert not isinstance(value, BaseCard), value
        elif self.rtype in ['FRSTRE']:
            data = [prop if isinstance(prop, integer_types) else prop.pid for prop in self.atti]
            for value in data:
                #print('atti =', value, type(value))
                assert not isinstance(value, BaseCard), value
        elif self.rtype in ['WEIGHT', 'FLUTTER', 'STABDER', 'EIGN', 'FREQ']:
            data = self.atti
        elif self.rtype in ['DISP', 'FRDISP', 'TDISP',
                            'FRVELO',
                            'FRACCL', 'PSDACCL']:
            #self.atti = model.Nodes(self.atti, msg=msg)
            data = [node if isinstance(node, integer_types) else node.nid for node in self.atti]
        elif self.rtype in op2_results:
            data = self.atti
            for value in data:
                assert not isinstance(value, BaseCard), 'rtype=%s value=%s' % (self.rtype, value)
        elif self.rtype in ['GPFORCP']: # MSC Nastran specific
            data = self.atti
            for value in data:
                assert not isinstance(value, BaseCard), value
        else:
            msg = 'rtype=%s ptype=%s\n' % (self.rtype, self.ptype)
            #msg += str(self)
            raise NotImplementedError(msg)
        return data

    def raw_fields(self):
        list_fields = ['DRESP1', self.oid, self.label, self.rtype, self.ptype,
                       self.region, self.atta, self.attb] + self.atti_values()
        #for val in self.atti_values():
            #print(val)
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        if is_double:
            return self.comment + print_card_double(card)
        return self.comment + print_card_16(card)


class DRESP2(OptConstraint):
    type = 'DRESP2'

    def __init__(self, oid, label, dequation, region, method,
                 c1, c2, c3, params, comment=''):
        """
        Design Sensitivity Equation Response Quantities
        Defines equation responses that are used in the design, either as
        constraints or as an objective.
        """
        if comment:
            self._comment = comment
        self.func = None
        self.dequation_str = None
        self.oid = oid
        self.label = label
        self.dequation = dequation
        self.region = region
        self.method = method
        self.c1 = c1
        self.c2 = c2
        self.c3 = c3
        self.params = params

    @classmethod
    def add_card(cls, card, comment=''):
        oid = integer(card, 1, 'oid')
        label = string(card, 2, 'label')
        dequation = integer_or_string(card, 3, 'dequation_id')
        region = integer_or_blank(card, 4, 'region')
        method = string_or_blank(card, 5, 'method', 'MIN')
        c1 = double_or_blank(card, 6, 'c1', 100.)
        c2 = double_or_blank(card, 7, 'c2', 0.005)
        c3 = double_or_blank(card, 8, 'c3') #: .. todo:: or blank?

        fields = [interpret_value(field) for field in card[9:]]
        key = None  # dummy key
        params = {}
        value_list = []
        j = 0
        for (i, field) in enumerate(fields):
            if i % 8 == 0 and field is not None:
                if i > 0:
                    assert len(value_list) > 0, 'key=%s values=%s' % (key, value_list)
                    params[key] = value_list
                    j += 1
                key = (j, field)
                value_list = []
            elif field is not None:
                value_list.append(field)
        params[key] = value_list

        #print("--Params--")
        #for key, value_list in sorted(iteritems(self.params)):
            #print("  key=%s params=%s" %(key, value_list))
        return DRESP2(oid, label, dequation, region, method,
                      c1, c2, c3, params, comment=comment)

    def OptID(self):
        return self.oid

    def _verify(self, xref=True):
        pass
        #for (j, name), value_list in sorted(iteritems(self.params)):
            #print('  DRESP2 verify - key=%s values=%s' % (name,
                #self._get_values(name, value_list)))

    def calculate(self, op2_model, subcase_id):
        print(str(self))
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
        print('DRESP2 args = %s' % argsi)
        out = self.func(*argsi)
        print('  deqatn out = %s' % out)
        return out

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ', which is required by %s ID=%s' % (self.type, self.oid)
        default_values = {}
        for key, vals in sorted(iteritems(self.params)):
            j, name = key
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
                self.dtable = model.dtable
                for i, val in enumerate(vals):
                    default_values[val] = self.dtable[val]
            else:
                raise NotImplementedError('  TODO: xref %s' % str(key))
        self.params_ref = self.params

        if isinstance(self.DEquation(), int):
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
        if isinstance(self.dequation, int):
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
                if isinstance(val, int):
                    out.append(val)
                else:
                    out.append(val.OptID())

        elif name in ['DVCREL1', 'DVCREL2']:
            for i, val in enumerate(values_list):
                if isinstance(val, int):
                    out.append(val)
                else:
                    out.append(val.OptID())

        elif name in ['DVMREL1', 'DVMREL2']:
            for i, val in enumerate(values_list):
                if isinstance(val, int):
                    out.append(val)
                else:
                    out.append(val.OptID())

        elif name in ['DVPREL1', 'DVPREL2']:
            for i, val in enumerate(values_list):
                if isinstance(val, int):
                    out.append(val)
                else:
                    out.append(val.OptID())
        elif name == 'DESVAR':
            for i, val in enumerate(values_list):
                if isinstance(val, int):
                    out.append(val)
                else:
                    out.append(val.OptID())
        elif name == 'DNODE':
            #print(values_list)
            for i in range(0, len(values_list), 2):
                val = values_list[i]
                if isinstance(val, int):
                    out.append(val)
                else:
                    out.append(val.nid)
                out.append(values_list[i+1])

        elif name == 'DTABLE':
            out = values_list
        else:
            raise NotImplementedError('  TODO: _get_values %s' % str(name))
            out = values_list
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
        list_fields = ['DRESP2', self.oid, self.label, self.DEquation(),
                       self.region, self.method, self.c1, self.c2, self.c3]
        list_fields += self._pack_params()
        return list_fields

    def repr_fields(self):
        method = set_blank_if_default(self.method, 'MIN')
        c1 = set_blank_if_default(self.c1, 100.)
        c2 = set_blank_if_default(self.c2, 0.005)

        list_fields = ['DRESP2', self.oid, self.label, self.DEquation(),
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
    type = 'DRESP3'

    def __init__(self, oid, label, group, Type, region, params,
                 comment=''):
        if comment:
            self._comment = comment
        self.oid = oid
        self.label = label
        self.group = group
        self.Type = Type
        self.region = region
        self.params = params

    @classmethod
    def add_card(cls, card, comment=''):
        oid = integer(card, 1, 'ID')
        label = string(card, 2, 'label')
        group = string(card, 3, 'group')
        Type = string(card, 4, 'Type')
        region = integer(card, 5, 'Type')

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
        return DRESP3(oid, label, group, Type, region, params,
                      comment=comment)

    def _pack_params(self):
        # # the amount of padding at the [beginning,end] of the 2nd line
        packLength = {
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
            fields2 = [key] + value_list
            try:
                (i, j) = packLength[key]
            except KeyError:
                msg = 'INVALID DRESP2 key=%r fields=%s ID=%s' % (key, value_list, self.oid)
                raise KeyError(msg)
            list_fields += build_table_lines(fields2, nstart=i, nend=j)
        return list_fields

    def raw_fields(self):
        list_fields = [
            'DRESP3', self.oid, self.label, self.group, self.Type, self.region,
            None, None, None, None]
        list_fields += self._pack_params()
        return list_fields

    def repr_fields(self):
        list_fields = [
            'DRESP3', self.oid, self.label, self.group, self.Type, self.region,
            None, None, None, None]
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

    def __init__(self, dcid, dconstrs, comment=''):
        if comment:
            self._comment = comment
        self.dcid = dcid
        self.dconstrs = dconstrs

    @classmethod
    def add_card(cls, card, comment=''):
        dcid = integer(card, 1, 'dcid')
        dconstrs = []

        for i in range(1, len(card)):
            dconstr = integer(card, i, 'dconstr_%i' % i)
            dconstrs.append(dconstr)
        return DCONADD(dcid, dconstrs, comment=comment)

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        self.dconstrs = [model.dconstrs[dcid] for dcid in self.dconstr_ids]
        self.dconstrs_ref = self.dconstrs

    def uncross_reference(self):
        self.dconstrs = self.dconstr_ids
        del self.dconstrs_ref

    @property
    def dconstr_ids(self):
        return [dconstr if isinstance(dconstr, integer_types) else dconstr.oid
                for dconstr in self.dconstrs]

    def raw_fields(self):
        list_fields = ['DCONADD', self.dcid] + self.dconstr_ids
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

    def __init__(self, card=None, data=None, comment=''):
        if comment:
            self._comment = comment
        if card:
            #: Response type for which the screening criteria apply. (Character)
            self.rType = string(card, 1, 'rType')
            #: Truncation threshold. (Real; Default = -0.5)
            self.trs = double_or_blank(card, 2, 'trs', -0.5)
            #: Maximum number of constraints to be retained per region per load
            #: case. (Integer > 0; Default = 20)
            self.nstr = integer_or_blank(card, 3, 'nstr', 20)
            assert len(card) == 4, 'len(DSCREEN card) = %i' % len(card)
        else:
            raise RuntimeError(data)

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


class DVMREL1(OptConstraint):  # similar to DVPREL1
    type = 'DVMREL1'

    def __init__(self, oid, Type, mid, mpName, mpMin, mpMax, c0,
                 dvids, coeffs, comment=''):
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
            self._comment = comment
        self.oid = oid
        self.Type = Type
        self.mid = mid
        self.mpName = mpName
        self.mpMax = mpMax
        self.mpMin = mpMin
        self.c0 = c0
        self.dvids = dvids
        self.coeffs = coeffs

    @classmethod
    def add_card(cls, card, comment=''):
        oid = integer(card, 1, 'oid')
        Type = string(card, 2, 'Type')
        mid = integer(card, 3, 'mid')
        mpName = string(card, 4, 'mpName')
        #if self.mpName in ['E', 'RHO', 'NU']:  positive values
            #self.mpMin = double_or_blank(card, 5, 'mpMin', 1e-15)
        #else: # negative
            #self.mpMin = double_or_blank(card, 5, 'mpMin', -1e-35)
        mpMin = double_or_blank(card, 5, 'mpMin')  #: .. todo:: bad default
        mpMax = double_or_blank(card, 6, 'mpMax', 1e20)
        c0 = double_or_blank(card, 7, 'c0', 0.0)

        dvids = []
        coeffs = []
        end_fields = [interpret_value(field) for field in card[9:]]
        #print "end_fields = ",end_fields
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
        return DVMREL1(oid, Type, mid, mpName, mpMin, mpMax, c0,
                       dvids, coeffs, comment=comment)

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

    def OptID(self):
        return self.oid

    def Mid(self):
        if isinstance(self.mid, integer_types):
            return self.mid
        return self.mid_ref.mid

    def raw_fields(self):
        list_fields = ['DVMREL1', self.oid, self.Type, self.Mid(),
                       self.mpName, self.mpMin, self.mpMax, self.c0, None]
        for (dvid, coeff) in zip(self.dvids, self.coeffs):
            list_fields.append(dvid)
            list_fields.append(coeff)
        return list_fields

    def repr_fields(self):
        mpMax = set_blank_if_default(self.mpMax, 1e20)
        c0 = set_blank_if_default(self.c0, 0.)
        list_fields = ['DVMREL1', self.oid, self.Type, self.Mid(),
                       self.mpName, self.mpMin, mpMax, c0, None]
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


class DVPREL1(OptConstraint):  # similar to DVMREL1
    type = 'DVPREL1'

    def __init__(self, oid, Type, pid, pNameFid, pMin, pMax, c0, dvids, coeffs,
                 comment=''):
        """
        +---------+--------+--------+--------+-----+
        | DVPREL1 | 200000 | PCOMP  | 2000   |  T2 |
        +---------+--------+--------+--------+-----+
        |         | 200000 |   1.0  |        |     |
        +---------+--------+--------+--------+-----+
        """
        if comment:
            self._comment = comment
        self.oid = oid
        self.Type = Type
        self.pid = pid
        self.pNameFid = pNameFid
        self.pMin = pMin
        self.pMax = pMax
        self.c0 = c0
        self.dvids = dvids
        self.coeffs = coeffs

    @classmethod
    def add_card(cls, card, comment=''):
        oid = integer(card, 1, 'oid')
        Type = string(card, 2, 'Type')
        pid = integer(card, 3, 'pid')
        pNameFid = integer_or_string(card, 4, 'pName_FID')

        #: Minimum value allowed for this property.
        #: .. todo:: bad default (see DVMREL1)
        pMin = double_or_blank(card, 5, 'pMin')
        pMax = double_or_blank(card, 6, 'pMax', 1e20)
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
        return DVPREL1(oid, Type, pid, pNameFid, pMin, pMax, c0, dvids, coeffs,
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
        self.pid = model.Property(self.pid)
        self.pid_ref = self.pid

    def uncross_reference(self):
        self.pid = self.Pid()
        del self.pid_ref

    def calculate(self, op2_model, subcase_id):
        raise NotImplementedError('\n' + str(self))

    def Pid(self):
        if isinstance(self.pid, integer_types):
            return self.pid
        return self.pid_ref.pid

    def raw_fields(self):
        list_fields = ['DVPREL1', self.oid, self.Type, self.Pid(),
                       self.pNameFid, self.pMin, self.pMax, self.c0, None]
        for (dvid, coeff) in zip(self.dvids, self.coeffs):
            list_fields.append(dvid)
            list_fields.append(coeff)
        return list_fields

    def repr_fields(self):
        pMax = set_blank_if_default(self.pMax, 1e20)
        c0 = set_blank_if_default(self.c0, 0.)
        list_fields = ['DVPREL1', self.oid, self.Type, self.Pid(),
                       self.pNameFid, self.pMin, pMax, c0, None]
        for (dvid, coeff) in zip(self.dvids, self.coeffs):
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
    type = 'DVPREL2'

    def __init__(self, oid, Type, pid, pNameFid, pMin, pMax, deqation,
                 dvids, labels, comment=''):
        """
       +----------+--------+--------+-------+-----------+-------+-------+-------+-------+
       | DVPREL2  | ID     | TYPE   | PID   | PNAME/FID | PMIN  | PMAX  | EQID  |       |
       +----------+--------+--------+-------+-----------+-------+-------+-------+-------+
       |          | DESVAR | DVID1  | DVID2 | DVID3     | DVID4 | DVID5 | DVID6 | DVID7 |
       +----------+--------+--------+-------+-----------+-------+-------+-------+-------+
       |          | DVID8  | -etc.- |       |           |       |       |       |       |
       +----------+--------+--------+-------+-----------+-------+-------+-------+-------+
       |          | DTABLE | LABL1  | LABL2 | LABL3     | LABL4 | LABL5 | LABL6 | LABL7 |
       +----------+--------+--------+-------+-----------+-------+-------+-------+-------+
       |          | LABL8  | -etc.- |       |           |       |       |       |       |
       +----------+--------+--------+-------+-----------+-------+-------+-------+-------+
        """
        if comment:
            self._comment = comment
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
        self.pNameFid = pNameFid

        #: Minimum value allowed for this property. If FID references a stress
        #: recovery location field, then the default value for PMIN is -1.0+35.
        #: PMIN must be explicitly set to a negative number for properties that
        #: may be less than zero (for example, field ZO on the PCOMP entry).
        #: (Real; Default = 1.E-15)
        #: .. todo:: bad default (see DVMREL1)
        self.pMin = pMin
        #: Maximum value allowed for this property. (Real; Default = 1.0E20)
        self.pMax = pMax
        #: DEQATN entry identification number. (Integer > 0)
        self.dequation = deqation
        self.dvids = dvids
        self.labels = labels

    @classmethod
    def add_card(cls, card, comment=''):
        oid = integer(card, 1, 'oid')
        Type = string(card, 2, 'Type')
        pid = integer(card, 3, 'pid')
        pNameFid = integer_or_string(card, 4, 'pName_FID')
        pMin = double_or_blank(card, 5, 'pMin')
        pMax = double_or_blank(card, 6, 'pMax', 1e20)
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
        return DVPREL2(oid, Type, pid, pNameFid, pMin, pMax, dequation, dvids,
                       labels, comment=comment)

    def OptID(self):
        return self.oid

    def Pid(self):
        if isinstance(self.pid, integer_types):
            return self.pid
        return self.pid_ref.pid

    def DEquation(self):
        if isinstance(self.dequation, int):
            return self.dequation
        return self.dequation_ref.equation_id

    def calculate(self, op2_model, subcase_id):
        """
        this should really make a call the the DEQATN;
        see the PBEAM for an example of get/set_opt_value
        """
        try:
            get = self.pid_ref.get_optimization_value(self.pNameFid)
            out = self.pid_ref.set_optimization_value(self.pNameFid, get)
        except:
            print('DVPREL2 calculate : %s[%r] = ???' % (self.Type, self.pNameFid))
            raise

        if self.dvids:
            for dv in self.dvids: # DESVARS
                arg = dv.calculate(op2_model, subcase_id)
                argsi.append(arg)
        if self.labels:
            for label in self.labels: # DTABLE
                arg = self.dtable[label]
                argsi.append(arg)
        print('DVPREL2; args = %s' % argsi)

        print('dvids  =', self.dvids)
        print('labels =', self.labels)
        print('%s[%r] = %s' % (self.Type, self.pNameFid, out))
        out = self.func(*argsi)
        print('  deqatn out = %s' % out)
        return out
        #raise NotImplementedError('\n' + str(self))

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        .. todo:: add support for DEQATN cards to finish DVPREL2 xref
        """
        self.pid = model.Property(self.pid)
        self.dequation = model.DEQATN(self.dequation)

        self.pid_ref = self.pid
        self.dequation_ref = self.dequation
        assert self.pid_ref.type not in ['PBEND', 'PBARL', 'PBEAML'], self.pid

    def uncross_reference(self):
        self.pid = self.Pid()
        self.dequation = self.DEquation()
        del self.pid_ref, self.dequation_ref

    #def OptValue(self):  #: .. todo:: not implemented
        #self.pid_ref.OptValue(self.pNameFid)

    def raw_fields(self):
        list_fields = ['DVPREL2', self.oid, self.Type, self.Pid(),
                       self.pNameFid, self.pMin, self.pMax, self.DEquation(), None]
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
