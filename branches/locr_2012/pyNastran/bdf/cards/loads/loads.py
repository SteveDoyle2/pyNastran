# pylint: disable=C0103,R0902,R0904,R0914,W0231,R0201
"""
All static loads are defined in this file.  This includes:
 * LSEQ
 * DLOAD
 * DAREA
 * SLOAD
 * TLOAD1
 * TLOAD2
 * RFORCE
 * RLOAD1
 * RLOAD2
 * RANDPS
"""
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
#import sys
from itertools import izip

from pyNastran.bdf.fieldWriter import set_blank_if_default
from pyNastran.bdf.cards.baseCard import BaseCard
from pyNastran.bdf.format import (integer, integer_or_blank,
                                  double, double_or_blank,
                                  integer_string_or_blank,
                                  string_or_blank, integer_double_or_blank,
                                  components_or_blank)


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
            self.sid = integer(card, 1, 'sid')

            ## overall scale factor
            self.scale = double(card, 2, 'scale')

            ## individual scale factors (corresponds to loadIDs)
            self.scaleFactors = []

            ## individual loadIDs (corresponds to scaleFactors)
            self.loadIDs = []

            # alternating of scale factor & load set ID
            nLoads = len(card) - 3
            assert nLoads % 2 == 0
            for i in xrange(nLoads // 2):
                n = 2 * i + 3
                self.scaleFactors.append(double(card, n, 'scaleFactor'))
                self.loadIDs.append(integer(card, n + 1, 'loadID'))
        else:
            self.sid = data[0]
            self.scale = data[1]
            self.scaleFactors = data[2]
            self.loadIDs = data[3]
            assert len(data) == 4, '%s data=%s' % (self.type, data)

    def cross_reference(self, model):
        loadIDs2 = []
        msg = ' which is required by %s=%s' % (self.type, self.sid)
        for loadID in self.loadIDs:
            loadID2 = model.Load(loadID, msg=msg)
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
        return loads


class LSEQ(BaseCard):  # Requires LOADSET in case control deck
    """
    Defines a sequence of static load sets
    @todo how does this work...
    """
    type = 'LSEQ'

    def __init__(self, card=None, data=None, comment=''):
        if comment:
            self._comment = comment
        if card:
            self.sid = integer(card, 1, 'sid')
            self.exciteID = integer(card, 2, 'exciteID')
            self.lid = integer(card, 3, 'lid')
            self.tid = integer_or_blank(card, 4, 'tid')
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

    def cross_reference(self, model):
        msg = ' which is required by %s=%s' % (self.type, self.sid)
        self.lid = model.Load(self.lid, msg=msg)
        if self.tid:
            self.tid = model.Load(self.tid, msg=msg)

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
        list_fields = ['LSEQ', self.sid, self.exciteID, self.Lid(), self.Tid()]
        return list_fields

    def reprFields(self):
        return self.rawFields()


class DLOAD(LoadCombination):
    type = 'DLOAD'

    def __init__(self, card=None, data=None, comment=''):
        LoadCombination.__init__(self, card, data)
        if comment:
            self._comment = comment

    def rawFields(self):
        list_fields = ['DLOAD', self.sid, self.scale]
        for (scaleFactor, loadID) in izip(self.scaleFactors, self.loadIDs):
            list_fields += [scaleFactor, self.LoadID(loadID)]
        return list_fields

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

    def __init__(self, card=None, nOffset=0, data=None, comment=''):
        if comment:
            self._comment = comment
        if card:
            nOffset *= 3
            self.sid = integer(card, 1, 'sid')
            self.p = integer(card, 2 + nOffset, 'p')
            self.c = components_or_blank(card, 3 + nOffset, 'c', 0)
            self.scale = double(card, 4 + nOffset, 'scale')
        else:
            self.sid = data[0]
            self.p = data[1]
            self.c = data[2]
            self.scale = data[3]
            assert len(data) == 4, 'data = %s' % data

    def rawFields(self):
        list_fields = ['DAREA', self.sid, self.p, self.c, self.scale]
        return list_fields


class TabularLoad(BaseCard):
    def __init__(self, card, data):
        pass


class SLOAD(Load):
    """
    Static Scalar Load
    Defines concentrated static loads on scalar or grid points.

    @note Can be used in statics OR dynamics.

    If Si refers to a grid point, the load is applied to component T1 of the
    displacement coordinate system (see the CD field on the GRID entry).
    """
    type = 'SLOAD'

    def __init__(self, card=None, data=None, comment=''):
        if comment:
            self._comment = comment
            ## load ID
        self.sid = integer(card, 1, 'sid')

        nfields = len(card) - 2
        n = nfields // 2
        if nfields % 2 == 1:
            n += 1
            msg = 'Missing last magnitude on SLOAD card=%s' % card.fields()
            raise RuntimeError(msg)

        self.nids = []
        self.mags = []
        for i in xrange(n):
            j = 2 * i + 2
            self.nids.append(integer(card, j, 'nid' + str(i)))
            self.mags.append(double(card, j + 1, 'mag' + str(i)))

    def cross_reference(self, model):
        msg = ' which is required by %s=%s' % (self.type, self.sid)
        for (i, nid) in enumerate(self.nids):
            self.nids[i] = model.Node(nid, msg=msg)

    def Nid(self, node):
        if isinstance(node, int):
            return node
        return node.nid

    def getLoads(self):  ## @todo:  not done
        return []

    def rawFields(self):
        list_fields = ['SLOAD', self.sid]
        for (nid, mag) in izip(self.nids, self.mags):
            list_fields += [self.Nid(nid), mag]
        return list_fields

    def reprFields(self):
        return self.rawFields()


class TLOAD1(TabularLoad):
    r"""
    Transient Response Dynamic Excitation, Form 1
    Defines a time-dependent dynamic load or enforced motion of the form:
    \f[ {P(t)} = {A} \cdot F(t-\tau) \f]
    for use in transient response analysis.
    """
    type = 'TLOAD1'

    def __init__(self, card=None, data=None, comment=''):
        TabularLoad.__init__(self, card, data)
        if comment:
            self._comment = comment
        ## load ID
        self.sid = integer(card, 1, 'sid')

        ## Identification number of DAREA or SPCD entry set or a thermal load
        ## set (in heat transfer analysis) that defines {A}. (Integer > 0)
        self.exciteID = integer(card, 2, 'exciteID')

        ## If it is a non-zero integer, it represents the
        ## identification number of DELAY Bulk Data entry that defines . If it
        ## is real, then it directly defines the value of that will be used for
        ## all degrees-of-freedom that are excited by this dynamic load entry.
        ## See also Remark 9. (Integer >= 0, real or blank)
        self.delay = integer_double_or_blank(card, 3, 'delay')

        ## Defines the type of the dynamic excitation. (LOAD,DISP, VELO, ACCE)
        self.Type = integer_string_or_blank(card, 4, 'Type', 'LOAD')

        ## Identification number of TABLEDi entry that gives F(t). (Integer > 0)
        self.tid = integer(card, 5, 'tid')

        ## Factor for initial displacements of the enforced degrees-of-freedom.
        ## (Real; Default = 0.0)
        self.us0 = double_or_blank(card, 6, 'us0', 0.0)

        ## Factor for initial velocities of the enforced degrees-of-freedom.
        ## (Real; Default = 0.0)
        self.vs0 = double_or_blank(card, 7, 'vs0', 0.0)
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
            msg = ' which is required by %s=%s' % (self.type, self.sid)
            self.tid = model.Table(self.tid, msg=msg)

    def Tid(self):
        if self.tid == 0:
            return None
        elif isinstance(self.tid, int):
            return self.tid
        return self.tid.tid

    def rawFields(self):
        list_fields = ['TLOAD1', self.sid, self.exciteID, self.delay, self.Type,
                  self.Tid(), self.us0, self.vs0]
        return list_fields

    def reprFields(self):
        us0 = set_blank_if_default(self.us0, 0.0)
        vs0 = set_blank_if_default(self.vs0, 0.0)
        list_fields = ['TLOAD1', self.sid, self.exciteID, self.delay, self.Type,
                  self.Tid(), us0, vs0]
        return list_fields


class TLOAD2(TabularLoad):
    r"""
    Transient Response Dynamic Excitation, Form 1
    Defines a time-dependent dynamic load or enforced motion of the form:
    \f[ {P(t)} = {A} \cdot F(t-\tau) \f]
    for use in transient response analysis.
    """
    type = 'TLOAD2'

    def __init__(self, card=None, data=None, comment=''):
        TabularLoad.__init__(self, card, data)
        if comment:
            self._comment = comment
        ## load ID
        ## SID must be unique for all TLOAD1, TLOAD2, RLOAD1, RLOAD2, and ACSRCE entries.
        self.sid = integer(card, 1, 'sid')

        self.exciteID = integer(card, 2, 'exciteID')
        self.delay = integer_or_blank(card, 3, 'delay', 0)

        ## Defines the type of the dynamic excitation. (Integer; character
        ## or blank; Default = 0)
        self.Type = integer_string_or_blank(card, 4, 'Type', 'LOAD')

        ## Time constant. (Real >= 0.0)
        self.T1 = double_or_blank(card, 5, 'T1', 0.0)
        #if self.delay == 0:
        #self.T1 = double_or_blank(card, 5, 'T1', 0.)
        #else:
        #self.T1 = blank(card, 5, 'T1')

        ## Time constant. (Real; T2 > T1)
        self.T2 = double_or_blank(card, 6, 'T2', self.T1)
        ## Frequency in cycles per unit time. (Real >= 0.0; Default = 0.0)
        self.frequency = double_or_blank(card, 7, 'frequency', 0.)
        ## Phase angle in degrees. (Real; Default = 0.0)
        self.phase = double_or_blank(card, 8, 'phase', 0.)
        ## Exponential coefficient. (Real; Default = 0.0)
        self.c = double_or_blank(card, 9, 'c', 0.)
        ## Growth coefficient. (Real; Default = 0.0)
        self.b = double_or_blank(card, 10, 'b', 0.)
        ## Factor for initial displacements of the enforced degrees-of-freedom.
        ## (Real; Default = 0.0)
        self.us0 = double_or_blank(card, 11, 'us0', 0.)
        ## Factor for initial velocities of the enforced degrees-of-freedom
        ## (Real; Default = 0.0)
        self.vs0 = double_or_blank(card, 12, 'vs0', 0.)

        if self.Type in [0, 'L', 'LO', 'LOA', 'LOAD']:
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
        list_fields = ['TLOAD2', self.sid, self.exciteID, self.delay, self.Type,
                  self.T1, self.T2, self.frequency, self.phase, self.c, self.b,
                  self.us0, self.vs0]
        return list_fields

    def reprFields(self):
        frequency = set_blank_if_default(self.frequency, 0.0)
        phase = set_blank_if_default(self.phase, 0.0)
        c = set_blank_if_default(self.c, 0.0)
        b = set_blank_if_default(self.b, 0.0)

        us0 = set_blank_if_default(self.us0, 0.0)
        vs0 = set_blank_if_default(self.vs0, 0.0)
        list_fields = ['TLOAD2', self.sid, self.exciteID, self.delay, self.Type,
                  self.T1, self.T2, frequency, phase, c, b, us0, vs0]
        return list_fields


class RFORCE(Load):
    type = 'RFORCE'

    def __init__(self, card=None, data=None, comment=''):
        if comment:
            self._comment = comment
        if card:
            self.sid = integer(card, 1, 'sid')
            self.nid = integer(card, 2, 'nid')
            self.cid = integer_or_blank(card, 3, 'cid', 0)
            self.scale = double_or_blank(card, 4, 'scale', 1.)
            self.r1 = double_or_blank(card, 5, 'r1', 0.)
            self.r2 = double_or_blank(card, 6, 'r2', 0.)
            self.r3 = double_or_blank(card, 7, 'r3', 0.)
            self.method = integer_or_blank(card, 8, 'method', 1)
            self.racc = double_or_blank(card, 9, 'racc', 0.)
            self.mb = integer_or_blank(card, 10, 'mb', 0)
            self.idrf = integer_or_blank(card, 11, 'idrf', 0)
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
        list_fields = ['RFORCE', self.sid, self.nid, self.cid, self.scale,
                  self.r1, self.r2, self.r3, self.method, self.racc,
                  self.mb, self.idrf]
        return list_fields

    def reprFields(self):
        #method = set_blank_if_default(self.method,1)
        racc = set_blank_if_default(self.racc, 0.)
        mb = set_blank_if_default(self.mb, 0)
        idrf = set_blank_if_default(self.idrf, 0)
        list_fields = ['RFORCE', self.sid, self.nid, self.cid, self.scale,
                  self.r1, self.r2, self.r3, self.method, racc,
                  mb, idrf]
        return list_fields


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

    def __init__(self, card=None, data=None, comment=''):
        TabularLoad.__init__(self, card, data)
        if comment:
            self._comment = comment
        self.sid = integer(card, 1, 'sid')
        self.exciteID = integer(card, 2, 'exciteID')
        self.delay = integer_or_blank(card, 3, 'delay')
        self.dphase = integer_or_blank(card, 4, 'dphase')
        self.tc = integer_or_blank(card, 5, 'tc', 0)
        self.td = integer_or_blank(card, 6, 'td', 0)
        self.Type = string_or_blank(card, 7, 'Type', 'LOAD')

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
        msg = ' which is required by %s=%s' % (self.type, self.sid)
        if self.tc:
            self.tc = model.Table(self.tc, msg=msg)
        if self.td:
            self.td = model.Table(self.td, msg=msg)

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
        list_fields = ['RLOAD1', self.sid, self.exciteID, self.delay, self.dphase,
                  self.Tc(), self.Td(), self.Type]
        return list_fields

    def reprFields(self):
        Type = set_blank_if_default(self.Type, 'LOAD')
        list_fields = ['RLOAD1', self.sid, self.exciteID, self.delay, self.dphase,
                  self.Tc(), self.Td(), Type]
        return list_fields


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

    def __init__(self, card=None, data=None, comment=''):
        TabularLoad.__init__(self, card, data)
        if comment:
            self._comment = comment
        self.sid = integer(card, 1, 'sid')
        self.exciteID = integer(card, 2, 'exciteID')
        self.delay = integer_or_blank(card, 3, 'delay')
        self.dphase = integer_or_blank(card, 4, 'dphase')
        self.tb = integer_or_blank(card, 5, 'tb', 0)
        self.tp = integer_or_blank(card, 6, 'tp', 0)
        self.Type = string_or_blank(card, 7, 'Type', 'LOAD')

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
        msg = ' which is required by %s=%s' % (self.type, self.sid)
        if self.tb:
            self.tb = model.Table(self.tb, msg=msg)
        if self.tp:
            self.tp = model.Table(self.tp, msg=msg)

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
        list_fields = ['RLOAD2', self.sid, self.exciteID, self.delay, self.dphase,
                  self.Tb(), self.Tp(), self.Type]
        return list_fields

    def reprFields(self):
        Type = set_blank_if_default(self.Type, 0.0)
        list_fields = ['RLOAD2', self.sid, self.exciteID, self.delay, self.dphase,
                  self.Tb(), self.Tp(), Type]
        return list_fields


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

    def __init__(self, card=None, data=None, comment=''):
        if comment:
            self._comment = comment
        if card:
            ## Random analysis set identification number. (Integer > 0)
            ## Defined by RANDOM in the Case Control Deck.
            self.sid = integer(card, 1, 'sid')

            ## Subcase identification number of the excited load set.
            ## (Integer > 0)
            self.j = integer(card, 2, 'j')

            ## Subcase identification number of the applied load set.
            ## (Integer >= 0; K >= J)
            self.k = integer(card, 3, 'k')

            ## Components of the complex number. (Real)
            self.x = double_or_blank(card, 4, 'x', 0.0)
            self.y = double_or_blank(card, 5, 'y', 0.0)
            ## Identification number of a TABRNDi entry that defines G(F).
            self.tid = integer_or_blank(card, 6, 'tid', 0)

    def cross_reference(self, model):
        if self.tid:
            msg = ' which is required by %s=%s' % (self.type, self.sid)
            self.tid = model.Table(self.tid, msg=msg)

    def getLoads(self):
        return [self]

    def Tid(self):
        if self.tid == 0:
            return None
        elif isinstance(self.tid, int):
            return self.tid
        return self.tid.tid

    def rawFields(self):
        list_fields = ['RANDPS', self.sid, self.j, self.k, self.x, self.y,
                  self.Tid()]
        return list_fields

    def reprFields(self):
        return self.rawFields()