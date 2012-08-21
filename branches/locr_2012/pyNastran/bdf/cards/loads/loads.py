# pylint: disable=C0103,R0902,R0904,R0914,W0231,R0201
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
#import sys
from itertools import izip

from pyNastran.bdf.fieldWriter import set_blank_if_default
from pyNastran.bdf.cards.baseCard import BaseCard


class Load(BaseCard):
    """defines the DefaultLoad class"""
    type = 'DefLoad'

    def __init__(self, card, data):
        self.cid = None
        self.nodes = None

    def Cid(self):
        if isinstance(self.cid, int):
            return self.cid
        else:
            return self.cid.cid

    def nodeIDs(self, nodes=None):
        """returns nodeIDs for repr functions"""
        if not nodes:
            nodes = self.nodes
        if isinstance(nodes[0], int):
            return [node for node in nodes]
        else:
            return [node.nid for node in nodes]


class LoadCombination(Load):  # LOAD, DLOAD
    def __init__(self, card, data):
        Load.__init__(self, card, data)

        if card:
            ## load ID
            self.sid = card.field(1)

            ## overall scale factor
            self.scale = card.field(2)

            loads = card.fields(3)  # temp list
            nLoadFields = len(loads)
            #nLoads  = nLoadFields/2
            assert nLoadFields % 2 == 0

            ## individual scale factors (corresponds to loadIDs)
            self.scaleFactors = []

            ## individual loadIDs (corresponds to scaleFactors)
            self.loadIDs = []

            # alternating of scale factor & load set ID
            for i in xrange(0, nLoadFields, 2):
                self.scaleFactors.append(loads[i])
                self.loadIDs.append(loads[i + 1])
        else:
            self.sid = data[0]
            self.scale = data[1]
            self.scaleFactors = data[2]
            self.loadIDs = data[3]
            assert len(data) == 4, '%s data=%s' % (self.type, data)
        ###

    def cross_reference(self, model):
        loadIDs2 = []
        for loadID in self.loadIDs:
            loadID2 = model.Load(loadID)
            loadIDs2.append(loadID2)
        self.loadIDs = loadIDs2

    def LoadID(self, lid):
        if isinstance(lid, int):
            return lid
        elif isinstance(lid, list):
            return lid[0].sid
        else:
            raise NotImplementedError(lid)

    #def Sid(self):
    #    try:
    #        if isinstance(self.sid, int):
    #            return self.sid
    #        elif isinstance(self.sid, list):
    #            #sys.stderr.write('type(lid[0]) = %s' %(type(self.lid[0])))
    #            #sys.stderr.write("the offending load...%s" %(self.lid[0]))
    #            return self.sid[0].sid
    #        #elif isinstance(self.lid,load)
    #        else:
    #            #sys.stderr.write("the offending load...%s" %(self.lid))
    #            return self.sid.sid
    #    except:
    #        msg = "error in loads.py - self.lid=\n %s\n" % (str(self.sid))
    #        sys.stderr.write(msg)
    #        raise

    #def LoadID(self, loadID):
        #print("load = ",loadID)
        #if isinstance(loadID, int):
        #    return loadID
        #elif isinstance(loadID, list):
        #    return loadID[0].LoadID()
        ##print("self.lid = ",load.lid)
        #asdf
        #return load.lid

    def getLoads(self):
        """@note requires a cross referenced load"""
        loads = []
        for allLoads in self.loadIDs:
            for load in allLoads:
                loads += load.getLoads()
            #loads += self.ID  # @todo:  what does this mean, was uncommented
        ###
        return loads


class LSEQ(BaseCard):  # Requires LOADSET in case control deck
    """
    Defines a sequence of static load sets
    @todo how does this work...
    """
    type = 'LSEQ'

    def __init__(self, card=None, data=None):
        if card:
            self.sid = card.field(1)
            self.exciteID = card.field(2)
            self.lid = card.field(3)
            self.tid = card.field(4)
        else:
            self.sid = data[0]
            self.exciteID = data[1]
            self.lid = data[2]
            self.tid = data[3]
            raise NotImplementedError()

    #def nodeIDs(self, nodes=None):
        #"""returns nodeIDs for repr functions"""
        #if not nodes:
        #    nodes = self.nodes
        #if isinstance(nodes[0], int):
        #    return [node for node in nodes]
        #else:
        #    return [node.nid for node in nodes]
        ###

    def cross_reference(self, model):
        self.lid = model.Load(self.lid)
        if self.tid:
            self.tid = model.Load(self.tid)

    def LoadID(self, lid):
        if isinstance(lid, int):
            return lid
        elif isinstance(lid, list):
            return self.LoadID(lid[0])
        else:
            return lid.sid
        #raise RuntimeError(lid)

    def getLoads(self):
        return self.lid

    def Lid(self):
        if isinstance(self.lid, int):
            return self.lid
        else:
            return self.LoadID(self.lid)
            #raise NotImplementedError('LSEQ ' + str(self.lid) +
            #                          '\n%s' % (type(self.lid)))

    def Tid(self):
        if self.tid is None:
            return None
        if isinstance(self.tid, int):
            return self.tid
        return self.tid.tid

    def rawFields(self):
        fields = ['LSEQ', self.sid, self.exciteID, self.Lid(), self.Tid()]
        return fields

    def reprFields(self):
        return self.rawFields()


class SLOAD(Load):
    """
    Static Scalar Load
    Defines concentrated static loads on scalar or grid points.

    @note Can be used in statics OR dynamics.

    If Si refers to a grid point, the load is applied to component T1 of the
    displacement coordinate system (see the CD field on the GRID entry).
    """
    type = 'SLOAD'

    def __init__(self, card=None, data=None):
        ## load ID
        self.sid = card.field(1)

        fields = card.fields(2)
        n = len(fields) // 2
        if len(fields) % 2 == 1:
            n += 1
            msg = 'missing last magnitude on SLOAD card=%s' % (card.fields())
            raise RuntimeError(msg)

        self.nids = []
        self.mags = []
        for i in xrange(n):
            j = 2 * i
            self.nids.append(fields[j])
            self.mags.append(fields[j + 1])
        ###

    def cross_reference(self, model):
        for (i, nid) in enumerate(self.nids):
            self.nids[i] = model.Node(nid)
        ###

    def Nid(self, node):
        if isinstance(node, int):
            return node
        return node.nid

    def rawFields(self):
        fields = ['SLOAD', self.sid]
        for (nid, mag) in izip(self.nids, self.mags):
            fields += [self.Nid(nid), mag]
        return fields

    def reprFields(self):
        return self.rawFields()


class DLOAD(LoadCombination):
    type = 'DLOAD'

    def __init__(self, card=None, data=None):
        LoadCombination.__init__(self, card, data)

    #def cross_reference(self, model):
    #    for (i, sid) in enumerate(self.sids):
    #        self.sids[i] = model.Load(sid)
    #    ###

    #def Sid(self, sid):
    #    if isinstance(sid, int):
    #        return sid
    #    return sid.lid

    def rawFields(self):
        fields = ['DLOAD', self.sid, self.scale]
        for (scaleFactor, loadID) in izip(self.scaleFactors, self.loadIDs):
            fields += [scaleFactor, self.LoadID(loadID)]
        return fields

    def reprFields(self):
        return self.rawFields()


class DAREA(BaseCard):
    """
    Defines scale (area) factors for static and dynamic loads. In dynamic
    analysis, DAREA is used in conjunction with ACSRCE, RLOADi and TLOADi
    entries.
    DAREA SID P1 C1 A1  P2 C2 A2
    DAREA 3   6   2 8.2 15 1  10.1
    """
    type = 'DAREA'

    def __init__(self, card=None, nOffset=0, data=None):
        if card:
            nOffset *= 3
            self.sid = card.field(1)
            self.p = card.field(2 + nOffset)
            self.c = card.field(3 + nOffset)
            self.scale = card.field(4 + nOffset)
        else:
            self.sid = data[0]
            self.p = data[1]
            self.c = data[2]
            self.scale = data[3]
            assert len(data) == 4, 'data = %s' % data
        ###

    def rawFields(self):
        fields = ['DAREA', self.sid, self.p, self.c, self.scale]
        return fields


class TabularLoad(BaseCard):
    def __init__(self, card, data):
        pass


class TLOAD1(TabularLoad):
    r"""
    Transient Response Dynamic Excitation, Form 1
    Defines a time-dependent dynamic load or enforced motion of the form:
    \f[ {P(t)} = {A} \cdot F(t-\tau) \f]
    for use in transient response analysis.
    """
    type = 'TLOAD1'

    def __init__(self, card=None, data=None):
        TabularLoad.__init__(self, card, data)
        ## load ID
        self.sid = card.field(1)

        ## Identification number of DAREA or SPCD entry set or a thermal load
        ## set (in heat transfer analysis) that defines {A}. (Integer > 0)
        self.exciteID = card.field(2)

        ## If it is a non-zero integer, it represents the
        ## identification number of DELAY Bulk Data entry that defines . If it
        ## is real, then it directly defines the value of that will be used for
        ## all degrees-of-freedom that are excited by this dynamic load entry.
        ## See also Remark 9. (Integer >= 0, real or blank)
        self.delay = card.field(3)

        ## Defines the type of the dynamic excitation. (LOAD,DISP, VELO, ACCE)
        self.Type = card.field(4, 'LOAD')

        ## Identification number of TABLEDi entry that gives F(t). (Integer > 0)
        self.tid = card.field(5)

        ## Factor for initial displacements of the enforced degrees-of-freedom.
        ## (Real; Default = 0.0)
        self.us0 = card.field(6, 0.0)

        ## Factor for initial velocities of the enforced degrees-of-freedom.
        ## (Real; Default = 0.0)
        self.vs0 = card.field(7, 0.0)
        if   self.Type in [0, 'L', 'LO', 'LOA', 'LOAD']:
            self.Type = 'LOAD'
        elif self.Type in [1, 'D', 'DI', 'DIS', 'DISP']:
            self.Type = 'DISP'
        elif self.Type in [2, 'V', 'VE', 'VEL', 'VELO']:
            self.Type = 'VELO'
        elif self.Type in [3, 'A', 'AC', 'ACC', 'ACCE']:
            self.Type = 'ACCE'
        else:
            msg = 'invalid TLOAD1 type  Type=|%s|' % self.Type
            raise RuntimeError(msg)

    def getLoads(self):
        return [self]

    def cross_reference(self, model):
        if self.tid:
            self.tid = model.Table(self.tid)

    def Tid(self):
        if self.tid == 0:
            return None
        elif isinstance(self.tid, int):
            return self.tid
        return self.tid.tid

    def rawFields(self):
        fields = ['TLOAD1', self.sid, self.exciteID, self.delay, self.Type,
                  self.Tid(), self.us0, self.vs0]
        return fields

    def reprFields(self):
        us0 = set_blank_if_default(self.us0, 0.0)
        vs0 = set_blank_if_default(self.vs0, 0.0)
        fields = ['TLOAD1', self.sid, self.exciteID, self.delay, self.Type,
                  self.Tid(), us0, vs0]
        return fields


class TLOAD2(TabularLoad):
    r"""
    Transient Response Dynamic Excitation, Form 1
    Defines a time-dependent dynamic load or enforced motion of the form:
    \f[ {P(t)} = {A} \cdot F(t-\tau) \f]
    for use in transient response analysis.
    """
    type = 'TLOAD2'

    def __init__(self, card=None, data=None):
        TabularLoad.__init__(self, card, data)
        ## load ID
        ## SID must be unique for all TLOAD1, TLOAD2, RLOAD1, RLOAD2, and ACSRCE entries.
        self.sid = card.field(1)

        self.exciteID = card.field(2)
        self.delay = card.field(3, 0)

        ## Defines the type of the dynamic excitation. (Integer; character
        ## or blank; Default = 0)
        self.Type = card.field(4, 0)

        ## Time constant. (Real >= 0.0)
        if self.delay == 0:
            self.T1 = card.field(5, 0.)
        else:
            self.T1 = card.field(5)
        ## Time constant. (Real; T2 > T1)
        self.T2 = card.field(6, self.T1)
        ## Frequency in cycles per unit time. (Real >= 0.0; Default = 0.0)
        self.frequency = card.field(7, 0.)
        ## Phase angle in degrees. (Real; Default = 0.0)
        self.phase = card.field(8, 0.)
        ## Exponential coefficient. (Real; Default = 0.0)
        self.c = card.field(9, 0.)
        ## Growth coefficient. (Real; Default = 0.0)
        self.b = card.field(10, 0.)
        ## Factor for initial displacements of the enforced degrees-of-freedom.
        ## (Real; Default = 0.0)
        self.us0 = card.field(11, 0.)
        ## Factor for initial velocities of the enforced degrees-of-freedom
        ## (Real; Default = 0.0)
        self.vs0 = card.field(12, 0.)

        if   self.Type in [0, 'L', 'LO', 'LOA', 'LOAD']:
            self.Type = 'LOAD'
        elif self.Type in [1, 'D', 'DI', 'DIS', 'DISP']:
            self.Type = 'DISP'
        elif self.Type in [2, 'V', 'VE', 'VEL', 'VELO']:
            self.Type = 'VELO'
        elif self.Type in [3, 'A', 'AC', 'ACC', 'ACCE']:
            self.Type = 'ACCE'
        elif self.Type in [5, 6, 7, 12, 13]:
            pass
        else:
            msg = 'invalid TLOAD1 type  Type=|%s|' % self.Type
            raise RuntimeError(msg)

    def getLoads(self):
        return [self]

    def cross_reference(self, model):
        pass
        # delay
        # exciteID

    def rawFields(self):
        fields = ['TLOAD2', self.sid, self.exciteID, self.delay, self.Type,
                  self.T1, self.T2, self.frequency, self.phase, self.c, self.b,
                  self.us0, self.vs0]
        return fields

    def reprFields(self):
        #self.Type = card.field(4,0)
        #self.T1 = card.field(5,0.)
        #self.T2 = card.field(6,self.T1)
        frequency = set_blank_if_default(self.frequency, 0.0)
        phase = set_blank_if_default(self.phase, 0.0)
        c = set_blank_if_default(self.c, 0.0)
        b = set_blank_if_default(self.b, 0.0)

        us0 = set_blank_if_default(self.us0, 0.0)
        vs0 = set_blank_if_default(self.vs0, 0.0)
        fields = ['TLOAD2', self.sid, self.exciteID, self.delay, self.Type,
                  self.T1, self.T2, frequency, phase, c, b, us0, vs0]
        return fields


class RFORCE(Load):
    type = 'RFORCE'

    def __init__(self, card=None, data=None):
        if card:
            self.sid = card.field(1)
            self.nid = card.field(2)
            self.cid = card.field(3, 0)
            self.scale = card.field(4, 1.)
            self.r1 = card.field(5, 0.)
            self.r2 = card.field(6, 0.)
            self.r3 = card.field(7, 0.)
            self.method = card.field(8, 1)
            self.racc = card.field(9, 0.)
            self.mb = card.field(10, 0)
            self.idrf = card.field(11, 0)
        else:
            self.sid = data[0]
            print("PLOADX1 = %s" % data)
            raise NotImplementedError('PLOADX1')

    def cross_reference(self, model):
        #self.nid = model.Element(self.nid)
        #self.cid = model.Coord(self.cid)
        pass

    def getLoads(self):
        return [self]

    def rawFields(self):
        fields = ['RFORCE', self.sid, self.nid, self.cid, self.scale,
                  self.r1, self.r2, self.r3, self.method, self.racc,
                  self.mb, self.idrf]
        return fields

    def reprFields(self):
        #method = set_blank_if_default(self.method,1)
        racc = set_blank_if_default(self.racc, 0.)
        mb = set_blank_if_default(self.mb, 0)
        idrf = set_blank_if_default(self.idrf, 0)
        fields = ['RFORCE', self.sid, self.nid, self.cid, self.scale,
                  self.r1, self.r2, self.r3, self.method, racc,
                  mb, idrf]
        return fields


class RLOAD1(TabularLoad):
    r"""
    Defines a frequency-dependent dynamic load of the form
    for use in frequency response problems.
    RLOAD1 SID EXCITEID DELAY DPHASE TC TD TYPE
    \f[ \large \left\{ P(f)  \right\}  = \left\{A\right\} [ C(f)+iD(f)]
        e^{  i \left\{\theta - 2 \pi f \tau \right\} } \f]
    RLOAD1 5   3                     1
    """
    type = 'RLOAD1'

    def __init__(self, card=None, data=None):
        TabularLoad.__init__(self, card, data)
        self.sid = card.field(1)  # was sid
        self.exciteID = card.field(2)
        self.delay = card.field(3)
        self.dphase = card.field(4)
        self.tc = card.field(5, 0)
        self.td = card.field(6, 0)
        self.Type = card.field(7, 'LOAD')

        if   self.Type in [0, 'L', 'LO', 'LOA', 'LOAD']:
            self.Type = 'LOAD'
        elif self.Type in [1, 'D', 'DI', 'DIS', 'DISP']:
            self.Type = 'DISP'
        elif self.Type in [2, 'V', 'VE', 'VEL', 'VELO']:
            self.Type = 'VELO'
        elif self.Type in [3, 'A', 'AC', 'ACC', 'ACCE']:
            self.Type = 'ACCE'
        else:
            msg = 'invalid RLOAD1 type  Type=|%s|' % self.Type
            raise RuntimeError(msg)

    def cross_reference(self, model):
        if self.tc:
            self.tc = model.Table(self.tc)
        if self.td:
            self.td = model.Table(self.td)

    #def LoadID(self, lid):
        #return self.Lid()

    def getLoads(self):
        return [self]

    def Tc(self):
        if self.tc == 0:
            return None
        elif isinstance(self.tc, int):
            return self.tc
        return self.tc.tid

    def Td(self):
        if self.td == 0:
            return None
        elif isinstance(self.td, int):
            return self.td
        return self.td.tid

    def rawFields(self):
        fields = ['RLOAD1', self.sid, self.exciteID, self.delay, self.dphase,
                  self.Tc(), self.Td(), self.Type]
        return fields

    def reprFields(self):
        Type = set_blank_if_default(self.Type, 'LOAD')
        fields = ['RLOAD1', self.sid, self.exciteID, self.delay, self.dphase,
                  self.Tc(), self.Td(), Type]
        return fields


class RLOAD2(TabularLoad):
    r"""
    Defines a frequency-dependent dynamic load of the form
    for use in frequency response problems.

    \f[ \large \left\{ P(f)  \right\}  = \left\{A\right\} * B(f)
        e^{  i \left\{ \phi(f) + \theta - 2 \pi f \tau \right\} } \f]
    RLOAD2 SID EXCITEID DELAY DPHASE TB TP TYPE
    RLOAD2 5   3                     1
    """
    type = 'RLOAD2'

    def __init__(self, card=None, data=None):
        TabularLoad.__init__(self, card, data)
        self.sid = card.field(1)
        self.exciteID = card.field(2)
        self.delay = card.field(3)
        self.dphase = card.field(4)
        self.tb = card.field(5, 0)
        self.tp = card.field(6, 0)
        self.Type = card.field(7, 'LOAD')

        if self.Type in [0, 'L', 'LO', 'LOA', 'LOAD']:
            self.Type = 'LOAD'
        elif self.Type in [1, 'D', 'DI', 'DIS', 'DISP']:
            self.Type = 'DISP'
        elif self.Type in [2, 'V', 'VE', 'VEL', 'VELO']:
            self.Type = 'VELO'
        elif self.Type in [3, 'A', 'AC', 'ACC', 'ACCE']:
            self.Type = 'ACCE'
        else:
            msg = 'invalid RLOAD2 type  Type=|%s|' % self.Type
            raise RuntimeError(msg)

    def cross_reference(self, model):
        if self.tb:
            self.tb = model.Table(self.tb)
        if self.tp:
            self.tp = model.Table(self.tp)

    def getLoads(self):
        return [self]

    #def Lid(self):
        #return self.lid

    def LoadID(self):
        return self.sid

    def Tb(self):
        if self.tb == 0:
            return None
        elif isinstance(self.tb, int):
            return self.tb
        return self.tb.tid

    def Tp(self):
        if self.tp == 0:
            return None
        elif isinstance(self.tp, int):
            return self.tp
        return self.tp.tid

    def rawFields(self):
        fields = ['RLOAD2', self.sid, self.exciteID, self.delay, self.dphase,
                  self.Tb(), self.Tp(), self.Type]
        return fields

    def reprFields(self):
        Type = set_blank_if_default(self.Type, 0.0)
        fields = ['RLOAD2', self.sid, self.exciteID, self.delay, self.dphase,
                  self.Tb(), self.Tp(), Type]
        return fields


class RandomLoad(BaseCard):
    def __init__(self, card, data):
        pass


class RANDPS(RandomLoad):
    r"""
    Power Spectral Density Specification
    Defines load set power spectral density factors for use in random analysis
    having the frequency dependent form:
    \f[ S_{jk}(F) = (X+iY)G(F) \f]
    """
    type = 'RANDPS'

    def __init__(self, card=None, data=None):
        if card:
            ## Random analysis set identification number. (Integer > 0)
            ## Defined by RANDOM in the Case Control Deck.
            self.sid = card.field(1)

            ## Subcase identification number of the excited load set.
            ## (Integer > 0)
            self.j = card.field(2)

            ## Subcase identification number of the applied load set.
            ## (Integer >= 0; K >= J)
            self.k = card.field(3)

            ## Components of the complex number. (Real)
            self.x = card.field(4)
            self.y = card.field(5)
            ## Identification number of a TABRNDi entry that defines G(F).
            self.tid = card.field(6, 0)

    def cross_reference(self, model):
        if self.tid:
            self.tid = model.Table(self.tid)

    def getLoads(self):
        return [self]

    def Tid(self):
        if self.tid == 0:
            return None
        elif isinstance(self.tid, int):
            return self.tid
        return self.tid.tid

    def rawFields(self):
        fields = ['RANDPS', self.sid, self.j, self.k, self.x, self.y,
                  self.Tid()]
        return fields

    def reprFields(self):
        return self.rawFields()
