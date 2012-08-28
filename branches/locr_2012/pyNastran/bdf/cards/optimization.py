# pylint: disable=C0103,R0902,R0904,R0914
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
#import sys
from itertools import izip

from pyNastran.bdf.fieldWriter import set_blank_if_default
from pyNastran.bdf.cards.baseCard import (BaseCard, expand_thru_by,
    collapse_thru_by_float)


class OptConstraint(BaseCard):
    def __init__(self):
        pass


class DCONSTR(OptConstraint):
    type = 'DCONSTR'

    def __init__(self, card=None, data=None):
        if card:
            self.oid = card.field(1)
            self.rid = card.field(2)
            self.lid = card.field(3, -1e20)
            self.uid = card.field(4, 1e20)
            self.lowfq = card.field(5, 0.0)
            self.highfq = card.field(6, 1e20)
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

    def __init__(self, card=None, data=None):
        self.oid = card.field(1)
        self.label = card.field(2)
        self.xinit = card.field(3)
        self.xlb = card.field(4, -1e20)
        self.xub = card.field(5, 1e20)
        self.delx = card.field(6, 1e20)
        self.ddval = card.field(7)

    def rawFields(self):
        fields = ['DESVAR', self.oid, self.label, self.xinit, self.xlb, self.xub,
                  self.delx, self.ddval]
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

    def __init__(self, card=None, data=None):
        self.oid = card.field(1)
        ddvals = [ddval for ddval in card.fields(2) if ddval is not None]
        self.ddvals = expand_thru_by(ddvals)
        self.ddvals.sort()

    def rawFields(self):
        self.ddvals.sort()
        ddvals = collapse_thru_by_float(self.ddvals)
        fields = ['DDVAL', self.oid] + ddvals
        return fields


class DOPTPRM(OptConstraint):
    type = 'DOPTPRM'

    def __init__(self, card=None, data=None):
        """
        Design Optimization Parameters
        Overrides default values of parameters used in design optimization
        
        @code
        DOPTPRM PARAM1 VAL1 PARAM2 VAL2 PARAM3 VAL3 PARAM4 VAL4
                PARAM5 VAL5 -etc.-
        @endcode
        """
        
        fields = card.fields(1)
        nFields = len(fields)

        self.params = {}
        for i in xrange(0, nFields, 2):
            param = fields[i]
            val = fields[i + 1]
            self.params[param] = val

    def rawFields(self):
        fields = ['DOPTPRM']
        for param, val in sorted(self.params.iteritems()):
            fields += [param, val]
        return fields


class DLINK(OptConstraint):
    type = 'DLINK'

    def __init__(self, card=None, data=None):
        """
        Multiple Design Variable Linking
        Relates one design variable to one or more other design variables
        
        @code
        DLINK ID DDVID C0 CMULT IDV1 C1 IDV2 C2
              IDV3 C3 -etc.-
        @endcode
        """
        self.oid = card.field(1)
        self.ddvid = card.field(2)
        self.c0 = card.field(3, 0.)
        self.cmult = card.field(4, 1.)

        fields = card.fields(5)
        nFields = len(fields)
        self.IDv = []
        self.Ci = []

        for i in xrange(0, nFields, 2):
            self.IDv.append(fields[i])
            self.Ci.append(fields[i + 1])
        ###

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


class DSCREEN(OptConstraint):
    type = 'DSCREEN'

    def __init__(self, card=None, data=None):
        ## Response type for which the screening criteria apply. (Character)
        self.rType = card.field(1)
        ## Truncation threshold. (Real; Default = -0.5)
        self.trs = card.field(2, -0.5)
        ## Maximum number of constraints to be retained per region per load
        ## case. (Integer > 0; Default = 20)
        self.nstr = card.field(3, 20)

    def rawFields(self):
        fields = ['DSCREEN', self.rType, self.trs, self.nstr]
        return fields

    def reprFields(self):
        trs = set_blank_if_default(self.trs, -0.5)
        nstr = set_blank_if_default(self.nstr, 20)
        fields = ['DSCREEN', self.rType, trs, nstr]
        return fields


class DRESP1(OptConstraint):
    type = 'DRESP1'

    def __init__(self, card=None, data=None):
        """
        @code
        DRESP1         1S1      CSTRAIN PCOMP                  1       1   10000
        @endcode
        """
        self.oid = card.field(1)
        self.label = card.field(2)
        self.rtype = card.field(3)
        self.ptype = card.field(4)
        self.region = card.field(5)
        self.atta = card.field(6)
        self.attb = card.field(7)
        self.atti = card.field(8)
        self.others = card.fields(9)
        #if self.others:
        #    print("self.others = %s" %(self.others))
        #    print(str(self))
        #assert len(self.others)==0

    def rawFields(self):
        fields = ['DRESP1', self.oid, self.label, self.rtype, self.ptype, self.region, self.atta, self.attb, self.atti
                  ] + self.others
        return fields


class DRESP2(OptConstraint):
    type = 'DRESP2'

    def __init__(self, card=None, data=None):
        """
        Design Sensitivity Equation Response Quantities
        Defines equation responses that are used in the design, either as
        constraints or as an objective.
        """
        self.oid = card.field(1)
        self.label = card.field(2)
        self.eqidFunc = card.field(3)
        self.region = card.field(4)
        self.method = card.field(5, 'MIN')
        self.c1 = card.field(6, 100.)
        self.c2 = card.field(7, 0.005)
        self.c3 = card.field(8)

        i = 0
        fields = card.fields(9)
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
        packLength = {  # the amount of padding at the [beginning,end] of the 2nd line
                        'DESVAR': [1, 0],
                        'DTABLE': [1, 0],
                        'DRESP1': [1, 0],
                        'DNODE': [1, 1],  # unique entry
                        'DVPREL1': [1, 0],
                        'DVCREL1': [1, 0],
                        'DVMREL1': [1, 0],
                        'DVPREL2': [1, 0],
                        'DVCREL2': [1, 0],
                        'DVMREL2': [1, 0],
                        'DRESP2': [1, 0],
                        'DESVAR': [1, 0],
                        'DESVAR': [1, 0],
                        'DESVAR': [1, 0],
                        'DESVAR': [1, 0],
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
        ###
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


class DVMREL1(OptConstraint):  # similar to DVPREL1
    type = 'DVMREL1'

    def __init__(self, card=None, data=None):
        """
        Design Variable to Material Relation
        Defines the relation between a material property and design variables
        @code
        DVMREL1 ID TYPE MID MPNAME MPMIN MPMAX C0
                DVID1 COEF1 DVID2 COEF2 DVID3 COEF3 -etc.-
        @endcode
        """
        self.oid = card.field(1)
        self.Type = card.field(2)
        self.mid = card.field(3)
        self.mpName = card.field(4)
        self.mpMin = card.field(5)  ## @todo bad default
        self.mpMax = card.field(6, 1e20)
        self.c0 = card.field(7, 0.0)

        self.dvids = []
        self.coeffs = []
        endFields = card.fields(9)
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

    def __init__(self, card=None, data=None):
        """
        @code
        DVPREL1   200000   PCOMP    2000      T2
                  200000     1.0
        @endcode
        """
        self.oid = card.field(1)
        self.Type = card.field(2)
        self.pid = card.field(3)
        self.pNameFid = card.field(4)
        self.pMin = card.field(5)  ## @todo bad default
        self.pMax = card.field(6, 1e20)
        self.c0 = card.field(7, 0.0)

        self.dvids = []
        self.coeffs = []
        endFields = card.fields(9)

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

    def __init__(self, card=None, data=None):
        """
        @code
        DVPREL2 ID TYPE PID PNAME/FID PMIN PMAX EQID
        'DESVAR' DVID1 DVID2 DVID3 DVID4 DVID5 DVID6 DVID7
                 DVID8 -etc.-
        'DTABLE' LABL1 LABL2 LABL3 LABL4 LABL5 LABL6 LABL7
                 LABL8 -etc.-
        @endcode
        """
        
        ## Unique identification number
        self.oid = card.field(1)
        ## Name of a property entry, such as PBAR, PBEAM, etc
        self.Type = card.field(2)
        ## Property entry identification number
        self.pid = card.field(3)
        ## Property name, such as 'T', 'A', or field position of the property
        ## entry, or word position in the element property table of the analysis
        ## model. Property names that begin with an integer such as 12I/T**3
        ## may only be referred to by field position. (Character or Integer 0)
        self.pnameFid = card.field(4)
        ## Minimum value allowed for this property. If FID references a stress
        ## recovery location field, then the default value for PMIN is -1.0+35.
        ## PMIN must be explicitly set to a negative number for properties that
        ## may be less than zero (for example, field ZO on the PCOMP entry).
        ## (Real; Default = 1.E-15)
        self.pmin = card.field(5, 1e-15)
        ## Maximum value allowed for this property. (Real; Default = 1.0E20)
        self.pmax = card.field(6, 1e20)
        ## DEQATN entry identification number. (Integer > 0)
        self.eqID = card.field(7)

        fields = card.fields(9)
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
                dvid = card.field(i)
                if dvid:
                    self.dvids.append(dvid)

        if iDTable:
            for i in xrange(iDTable + 1, iEnd):
                dtable = card.field(i)
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

    def OptValue(self):  ## @todo not implemented
        self.pid.OptValue(self.pnameFid)

    def rawFields(self):
        fields = ['DVPREL2', self.oid, self.Type, self.Pid(),
                  self.pnameFid, self.pmin, self.pmax, self.eqID, None]

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
