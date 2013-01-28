# pylint: disable=C0103,R0902,R0904,R0914
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
from itertools import izip

from pyNastran.bdf.fieldWriter import set_blank_if_default
from pyNastran.bdf.cards.baseCard import (BaseCard, expand_thru_by,
    collapse_thru_by_float)
from pyNastran.bdf.format import (integer, integer_or_blank, integer_or_string,
                                  double, double_or_blank, string,
                                  string_or_blank, integer_double_or_blank)

class OptConstraint(BaseCard):
    def __init__(self):
        pass


class DCONSTR(OptConstraint):
    type = 'DCONSTR'
    def __init__(self, card=None, data=None, comment=''):
        if comment:
            self._comment = comment
        if card:
            self.oid = integer(card, 1, 'oid')
            self.rid = integer(card, 2, 'rid')
            self.lid = double_or_blank(card, 3, 'lid', -1e20)
            self.uid = double_or_blank(card, 4, 'uid',  1e20)
            self.lowfq = double_or_blank(card, 5, 'lowfq', 0.0)
            self.highfq = double_or_blank(card, 6, 'highfq', 1e20)
        else:
            self.oid = data[0]
            self.rid = data[1]
            self.lid = data[2]
            self.uid = data[3]
            self.lowfq = data[4]
            self.highfq = data[5]

    def rawFields(self):
        fields = ['DCONSTR', self.oid, self.rid, self.lid,
                  self.uid, self.lowfq, self.highfq]
        return fields

    def reprFields(self):
        lid = set_blank_if_default(self.lid, -1e20)
        uid = set_blank_if_default(self.uid, 1e20)
        lowfq = set_blank_if_default(self.lowfq, 0.0)
        highfq = set_blank_if_default(self.highfq, 1e20)
        fields = ['DCONSTR', self.oid, self.rid, lid, uid, lowfq, highfq]
        return fields


class DESVAR(OptConstraint):
    type = 'DESVAR'
    def __init__(self, card=None, data=None, comment=''):
        if comment:
            self._comment = comment
        self.oid = integer(card, 1, 'oid')
        self.label = string(card, 2, 'label')
        self.xinit = double(card, 3, 'xinit')
        self.xlb = double_or_blank(card, 4, 'xlb', -1e20)
        self.xub = double_or_blank(card, 5, 'xub', 1e20)
        self.delx = double_or_blank(card, 6, 'delx', 1e20)
        self.ddval = integer_or_blank(card, 7, 'ddval')

    def rawFields(self):
        fields = ['DESVAR', self.oid, self.label, self.xinit, self.xlb,
                  self.xub, self.delx, self.ddval]
        return fields

    def reprFields(self):
        xlb = set_blank_if_default(self.xlb, -1e20)
        xub = set_blank_if_default(self.xub, 1e20)
        delx = set_blank_if_default(self.delx, 1e20)
        fields = ['DESVAR', self.oid, self.label, self.xinit, xlb,
                  xub, delx, self.ddval]
        return fields


class DDVAL(OptConstraint):
    type = 'DDVAL'

    def __init__(self, card=None, data=None, comment=''):
        if comment:
            self._comment = comment
        self.oid = integer(card, 1, 'oid')
        ddvals = [ddval for ddval in card[2:] if ddval is not None]
        self.ddvals = expand_thru_by(ddvals)
        self.ddvals.sort()

    def rawFields(self):
        self.ddvals.sort()
        ddvals = collapse_thru_by_float(self.ddvals)
        fields = ['DDVAL', self.oid] + ddvals
        return fields


class DOPTPRM(OptConstraint):
    type = 'DOPTPRM'

    def __init__(self, card=None, data=None, comment=''):
        """
        Design Optimization Parameters
        Overrides default values of parameters used in design optimization

        @code
        DOPTPRM PARAM1 VAL1 PARAM2 VAL2 PARAM3 VAL3 PARAM4 VAL4
                PARAM5 VAL5 -etc.-
        @endcode
        """
        if comment:
            self._comment = comment

        nFields = len(card) - 1
        self.params = {}
        for i in xrange(0, nFields, 2):
            param = string(card, i + 1, 'param')
            val = double(card, i + 2, 'value')
            self.params[param] = val

    def rawFields(self):
        fields = ['DOPTPRM']
        for param, val in sorted(self.params.iteritems()):
            fields += [param, val]
        return fields


class DLINK(OptConstraint):
    type = 'DLINK'

    def __init__(self, card=None, data=None, comment=''):
        """
        Multiple Design Variable Linking
        Relates one design variable to one or more other design variables

        @code
        DLINK ID DDVID C0 CMULT IDV1 C1 IDV2 C2
              IDV3 C3 -etc.-
        @endcode
        """
        if comment:
            self._comment = comment
        self.oid = integer(card, 1, 'oid')
        self.ddvid = integer(card, 2, 'ddvid')
        self.c0 = double_or_blank(card, 3, 'c0', 0.)
        self.cmult = double_or_blank(card, 4, 'cmult', 1.)

        nfields = len(card) - 4
        n = nfields // 2
        self.IDv = []
        self.Ci = []

        for i in xrange(n):
            j = 2 * i + 5
            IDv = integer(card, j, 'IDv' + str(i))
            Ci = double(card, j + 1, 'Ci' + str(i))
            self.IDv.append(IDv)
            self.Ci.append(Ci)

    def rawFields(self):
        fields = ['DLINK', self.oid, self.ddvid, self.c0, self.cmult]
        for (idv, ci) in izip(self.IDv, self.Ci):
            fields += [idv, ci]
        return fields

    def reprFields(self):
        c0 = set_blank_if_default(self.c0, 0.)
        cmult = set_blank_if_default(self.cmult, 1.)
        fields = ['DLINK', self.oid, self.ddvid, c0, cmult]
        for (idv, ci) in izip(self.IDv, self.Ci):
            fields += [idv, ci]
        return fields


class DRESP1(OptConstraint):
    type = 'DRESP1'

    def __init__(self, card=None, data=None, comment=''):
        """
        @code
        DRESP1         1S1      CSTRAIN PCOMP                  1       1   10000
        @endcode
        """
        if comment:
            self._comment = comment
        self.oid = integer(card, 1, 'oid')
        self.label = string(card, 2, 'label')
        self.rtype = string(card, 3, 'rtype')
        self.ptype = string_or_blank(card, 4, 'ptype') # elem, pbar, pshell, etc. (ELEM flag or Prop Name)
        self.region = integer_or_blank(card, 5, 'region')
        self.atta = integer_double_or_blank(card, 6, 'atta')
        self.attb = integer_double_or_blank(card, 7, 'attb')
        self.atti = integer_double_or_blank(card, 8, 'atti')
        self.others = card[9:]
        #if self.others:
        #    print("self.others = %s" %(self.others))
        #    print(str(self))
        #assert len(self.others)==0

    def rawFields(self):
        fields = ['DRESP1', self.oid, self.label, self.rtype, self.ptype,
                  self.region, self.atta, self.attb, self.atti] + self.others
        return fields


class DRESP2(OptConstraint):
    type = 'DRESP2'

    def __init__(self, card=None, data=None, comment=''):
        """
        Design Sensitivity Equation Response Quantities
        Defines equation responses that are used in the design, either as
        constraints or as an objective.
        """
        if comment:
            self._comment = comment
        self.oid = integer(card, 1, 'oid')
        self.label = string(card, 2, 'label')
        self.eqidFunc = integer_or_string(card, 3, 'eqid_Func')
        self.region = integer_or_blank(card, 4, 'region')
        self.method = string_or_blank(card, 5, 'method', 'MIN')
        self.c1 = double_or_blank(card, 6, 'c1', 100.)
        self.c2 = double_or_blank(card, 7, 'c2', 0.005)
        self.c3 = double_or_blank(card, 8, 'c3') # TODO: or blank?

        i = 0
        fields = card[9:]
        key = '$NULL$'  # dummy key
        self.params = {key: []}
        valueList = []
        for (i, field) in enumerate(fields):
            if i % 8 == 0 and field is not None:
                self.params[key] = valueList
                key = field
                valueList = []
            elif field is not None:
                valueList.append(field)
            #else:
            #    pass

        self.params[key] = valueList
        del self.params['$NULL$']

        #print "--Params--"
        #for (key, valueList) in sorted(self.params.iteritems()):
        #    print("  key=%s params=%s" %(key, valueList))

        #print self

    def packParams(self):
        # # the amount of padding at the [beginning,end] of the 2nd line
        packLength = {  
             'DESVAR':  [1, 0],
             'DTABLE':  [1, 0],
             'DRESP1':  [1, 0],
             'DNODE':   [1, 1],  # unique entry
             'DVPREL1': [1, 0],
             'DVCREL1': [1, 0],
             'DVMREL1': [1, 0],
             'DVPREL2': [1, 0],
             'DVCREL2': [1, 0],
             'DVMREL2': [1, 0],
             'DRESP2':  [1, 0],
             'DESVAR':  [1, 0],
             'DESVAR':  [1, 0],
             'DESVAR':  [1, 0],
             'DESVAR':  [1, 0],
        }
        fields = []
        for (key, valueList) in sorted(self.params.iteritems()):
            fields2 = [key] + valueList
            try:
                (i, j) = packLength[key]
            except KeyError:
                msg = 'INVALID DRESP2 key=|%s| fields=%s ID=%s' % (
                    key, valueList, self.oid)
                raise KeyError(msg)
            fields += self.buildTableLines(fields2, nStart=i, nEnd=j)
        return fields

    def rawFields(self):
        fields = ['DRESP2', self.oid, self.label, self.eqidFunc,
                  self.region, self.method, self.c1, self.c2, self.c3]
        fields += self.packParams()
        return fields

    def reprFields(self):
        method = set_blank_if_default(self.method, 'MIN')
        c1 = set_blank_if_default(self.c1, 100.)
        c2 = set_blank_if_default(self.c2, 0.005)

        fields = ['DRESP2', self.oid, self.label, self.eqidFunc,
                  self.region, method, c1, c2, self.c3]
        fields += self.packParams()
        return fields


class DSCREEN(OptConstraint):
    type = 'DSCREEN'

    def __init__(self, card=None, data=None, comment=''):
        if comment:
            self._comment = comment
        ## Response type for which the screening criteria apply. (Character)
        self.rType = string(card, 1, 'rType')
        ## Truncation threshold. (Real; Default = -0.5)
        self.trs = double_or_blank(card, 2, 'trs', -0.5)
        ## Maximum number of constraints to be retained per region per load
        ## case. (Integer > 0; Default = 20)
        self.nstr = integer(card, 3, 'nstr', 20)

    def rawFields(self):
        fields = ['DSCREEN', self.rType, self.trs, self.nstr]
        return fields

    def reprFields(self):
        trs = set_blank_if_default(self.trs, -0.5)
        nstr = set_blank_if_default(self.nstr, 20)
        fields = ['DSCREEN', self.rType, trs, nstr]
        return fields


class DVMREL1(OptConstraint):  # similar to DVPREL1
    type = 'DVMREL1'

    def __init__(self, card=None, data=None, comment=''):
        """
        Design Variable to Material Relation
        Defines the relation between a material property and design variables
        @code
        DVMREL1 ID TYPE MID MPNAME MPMIN MPMAX C0
                DVID1 COEF1 DVID2 COEF2 DVID3 COEF3 -etc.-
        @endcode
        """
        if comment:
            self._comment = comment
        self.oid = integer(card, 1, 'oid')
        self.Type = string(card, 2, 'Type')
        self.mid = integer(card, 3, 'mid')
        self.mpName = string(card, 4, 'mpName')
        #if self.mpName in ['E', 'RHO', 'NU']:  positive values
            #self.mpMin = double_or_blank(card, 5, 'mpMin', 1e-15)
        #else: # negative
            #self.mpMin = double_or_blank(card, 5, 'mpMin', -1e-35)
        self.mpMin = double_or_blank(card, 5, 'mpMin')  # TODO: bad default
        self.mpMax = double_or_blank(card, 6, 'mpMax', 1e20)
        self.c0 = double_or_blank(card, 7, 'c0', 0.0)

        self.dvids = []
        self.coeffs = []
        endFields = card[9:]
        #print "endFields = ",endFields
        nFields = len(endFields) - 1
        if nFields % 2 == 1:
            endFields.append(None)
            nFields += 1

        i = 0
        for i in xrange(0, nFields, 2):
            self.dvids.append(endFields[i])
            self.coeffs.append(endFields[i + 1])
        if nFields % 2 == 1:
            print(card)
            print("dvids = %s" % (self.dvids))
            print("coeffs = %s" % (self.coeffs))
            print(str(self))
            raise RuntimeError('invalid DVMREL1...')

    def cross_reference(self, model):
        self.mid = model.Material(self.mid)

    def Mid(self):
        if isinstance(self.mid, int):
            return self.mid
        return self.mid.mid

    def rawFields(self):
        fields = ['DVMREL1', self.oid, self.Type, self.Mid(),
                  self.mpName, self.mpMin, self.mpMax, self.c0, None]
        for (dvid, coeff) in izip(self.dvids, self.coeffs):
            fields.append(dvid)
            fields.append(coeff)
        return fields

    def reprFields(self):
        mpMax = set_blank_if_default(self.mpMax, 1e20)
        c0 = set_blank_if_default(self.c0, 0.)
        fields = ['DVMREL1', self.oid, self.Type, self.Mid(),
                  self.mpName, self.mpMin, mpMax, c0, None]
        for (dvid, coeff) in izip(self.dvids, self.coeffs):
            fields.append(dvid)
            fields.append(coeff)
        return fields


class DVPREL1(OptConstraint):  # similar to DVMREL1
    type = 'DVPREL1'

    def __init__(self, card=None, data=None, comment=''):
        """
        @code
        DVPREL1   200000   PCOMP    2000      T2
                  200000     1.0
        @endcode
        """
        if comment:
            self._comment = comment
        self.oid = integer(card, 1, 'oid')
        self.Type = string(card, 2, 'Type')
        self.pid = integer(card, 3, 'pid')
        self.pNameFid = integer_or_string(card, 4, 'pName_FID')
        self.pMin = double_or_blank(card, 5, 'pMin')  # TODO bad default (see DVMREL1)
        self.pMax = double_or_blank(card, 6, 'pMax', 1e20)
        self.c0 = double_or_blank(card, 7, 'c0', 0.0)

        self.dvids = []
        self.coeffs = []
        endFields = card[9:]

        nFields = len(endFields) - 1
        if nFields % 2 == 1:
            endFields.append(None)
            nFields += 1
        i = 0
        for i in xrange(0, nFields, 2):
            self.dvids.append(endFields[i])
            self.coeffs.append(endFields[i + 1])
        if nFields % 2 == 1:
            print(card)
            print("dvids = %s" % (self.dvids))
            print("coeffs = %s" % (self.coeffs))
            print(str(self))
            raise RuntimeError('invalid DVPREL1...')

    def cross_reference(self, model):
        self.pid = model.Property(self.pid)

    def Pid(self):
        if isinstance(self.pid, int):
            return self.pid
        return self.pid.pid

    def rawFields(self):
        fields = ['DVPREL1', self.oid, self.Type, self.Pid(),
                  self.pNameFid, self.pMin, self.pMax, self.c0, None]
        for (dvid, coeff) in izip(self.dvids, self.coeffs):
            fields.append(dvid)
            fields.append(coeff)
        return fields

    def reprFields(self):
        pMax = set_blank_if_default(self.pMax, 1e20)
        c0 = set_blank_if_default(self.c0, 0.)
        fields = ['DVPREL1', self.oid, self.Type, self.Pid(),
                  self.pNameFid, self.pMin, pMax, c0, None]

        for (dvid, coeff) in izip(self.dvids, self.coeffs):
            fields.append(dvid)
            fields.append(coeff)
        return fields


class DVPREL2(OptConstraint):
    type = 'DVPREL2'

    def __init__(self, card=None, data=None, comment=''):
        """
        @code
        DVPREL2 ID TYPE PID PNAME/FID PMIN PMAX EQID
        'DESVAR' DVID1 DVID2 DVID3 DVID4 DVID5 DVID6 DVID7
                 DVID8 -etc.-
        'DTABLE' LABL1 LABL2 LABL3 LABL4 LABL5 LABL6 LABL7
                 LABL8 -etc.-
        @endcode
        """
        if comment:
            self._comment = comment
        ## Unique identification number
        self.oid = integer(card, 1, 'oid')
        ## Name of a property entry, such as PBAR, PBEAM, etc
        self.Type = string(card, 2, 'Type')
        ## Property entry identification number
        self.pid = integer(card, 3, 'pid')
        ## Property name, such as 'T', 'A', or field position of the property
        ## entry, or word position in the element property table of the
        ## analysis model. Property names that begin with an integer such as
        ## 12I/T**3 may only be referred to by field position.
        ## (Character or Integer 0)
        self.pNameFid = integer_or_string(card, 4, 'pName_FID')
        ## Minimum value allowed for this property. If FID references a stress
        ## recovery location field, then the default value for PMIN is -1.0+35.
        ## PMIN must be explicitly set to a negative number for properties that
        ## may be less than zero (for example, field ZO on the PCOMP entry).
        ## (Real; Default = 1.E-15)
        self.pMin = double_or_blank(card, 5, 'pMin')  # TODO bad default (see DVMREL1)
        ## Maximum value allowed for this property. (Real; Default = 1.0E20)
        self.pMax = double_or_blank(card, 6, 'pMax', 1e20)
        ## DEQATN entry identification number. (Integer > 0)
        self.eqID = integer_or_blank(card, 7, 'eqID') # TODO:  or blank?

        fields = card[9:]
        #print "fields = ",fields
        iOffset = 9
        iEnd = len(fields) + iOffset

        try:
            iDesvar = fields.index('DESVAR') + iOffset
        except ValueError:
            iDesvar = None

        try:
            iDTable = fields.index('DTABLE') + iOffset
            #iDesMax  = iDTable # the index to start parsing DESVAR
            iDesStop = iDTable  # the index to stop  parsing DESVAR
        except ValueError:
            iDTable = None
            iDesStop = iEnd

        self.dvids = []
        self.dtables = []
        if iDesvar:
            for i in xrange(10, iDesStop):
                dvid = integer(card, i, 'dvid')
                if dvid:
                    self.dvids.append(dvid)

        if iDTable:
            for i in xrange(iDTable + 1, iEnd):
                dtable = integer(card, i, 'dtable')
                if dtable:
                    assert dtable is not 'DTABLE'
                    self.dtables.append(dtable)

    def Pid(self):
        if isinstance(self.pid, int):
            return self.pid
        return self.pid.pid

    #def EqID(self):

    def cross_reference(self, model):
        """@todo add support for DEQATN cards to finish DVPREL2 xref"""
        self.pid = model.Property(self.pid)
        #self.eqID = model.DEquation(self.eqID)

    def OptValue(self):  # TODO not implemented
        self.pid.OptValue(self.pNameFid)

    def rawFields(self):
        fields = ['DVPREL2', self.oid, self.Type, self.Pid(),
                  self.pNameFid, self.pMin, self.pMax, self.eqID, None]

        if self.dvids:
            fields2 = ['DESVAR'] + self.dvids
            fields += self.buildTableLines(fields2, nStart=1, nEnd=0)

        if self.dtables:
            fields2 = ['DTABLE'] + self.dtables
            fields += self.buildTableLines(fields2, nStart=1, nEnd=0)
        return fields

    def reprFields(self):
        """@todo finish reprFields for DVPREL2"""
        return self.rawFields()