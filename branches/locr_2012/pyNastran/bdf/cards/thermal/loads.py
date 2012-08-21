# pylint: disable=C0103,R0902,R0904,R0914
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)

from .thermal import ThermalCard
from pyNastran.bdf.fieldWriter import set_blank_if_default
from ..baseCard import expand_thru, expand_thru_by, collapse_thru_by


class ThermalLoadDefault(ThermalCard):
    def __init__(self, card, data):
        pass


class ThermalLoad(ThermalCard):
    def __init__(self, card, data):
        pass


class QBDY1(ThermalLoad):
    """
    Defines a uniform heat flux into CHBDYj elements.
    """
    type = 'QBDY1'

    def __init__(self, card=None, data=None):
        ThermalLoad.__init__(self, card, data)

        if card:
            ## Load set identification number. (Integer > 0)
            self.sid = card.field(1)

            ## Heat flux into element (FLOAT)
            self.qFlux = card.field(2)
            eids = card.fields(3)
            ## CHBDYj element identification numbers (Integer)
            self.eids = expand_thru(eids)  ## @warning should this use expand_thru_by ???
        else:
            self.sid = data[0]
            self.qFlux = data[1]
            self.eids = data[2:]

    def cross_reference(self, model):
        self.eids = model.Elements(self.eids)

    def Eid(self):
        if isinstance(self.eid, int):
            return self.eid
        return self.eid.eid

    def nQFluxTerms(self):
        return len(self.qFlux)

    def rawFields(self):
        fields = ['QBDY1', self.sid, self.qFlux] + list(
            self.eids) + [self.qFlux]
        return fields

    def reprFields(self):
        eids = collapse_thru_by(self.eids)
        fields = ['QBDY1', self.sid, self.qFlux] + list(eids) + [self.qFlux]
        return fields


class QBDY2(ThermalLoad):  # not tested
    """
    Defines a uniform heat flux load for a boundary surface.
    """
    type = 'QBDY2'

    def __init__(self, card=None, data=None):
        ThermalLoad.__init__(self, card, data)

        if card:
            ## Load set identification number. (Integer > 0)
            self.sid = card.field(1)
            ## Identification number of an CHBDYj element. (Integer > 0)
            self.eid = card.field(2)
            ## Heat flux at the i-th grid point on the referenced CHBDYj element. (Real or blank)
            self.qFlux = self.removeTrailingNones(card.fields(3))
        else:
            self.sid = data[0]
            self.eid = data[1]
            self.qFlux = data[2]

    def cross_reference(self, model):
        self.eid = model.Element(self.eid)

    def Eid(self):
        if isinstance(self.eid, int):
            return self.eid
        return self.eid.eid

    def nQFluxTerms(self):
        return len(self.qFlux)

    def rawFields(self):
        fields = ['QBDY2', self.sid, self.Eid(), self.qFlux]
        return fields

    def reprFields(self):
        return self.rawFields()


class QBDY3(ThermalLoad):
    """
    Defines a uniform heat flux load for a boundary surface.
    """
    type = 'QBDY3'

    def __init__(self, card=None, data=None):
        ThermalLoad.__init__(self, card, data)

        if card:
            ## Load set identification number. (Integer > 0)
            self.sid = card.field(1)
            ## Heat flux into element
            self.Q0 = card.field(2)
            ## Control point for thermal flux load. (Integer > 0; Default = 0)
            self.cntrlnd = card.field(3, 0)
            ## CHBDYj element identification numbers
            self.eids = expand_thru_by(card.fields(4))
        else:
            self.sid = data[0]
            self.Q0 = data[1]
            self.cntrlnd = data[2]
            self.eids = list(data[3:])

    def cross_reference(self, model):
        for i, eid in enumerate(self.eids):
            self.eids[i] = model.Element(eid)

    def Eids(self):
        eids = []
        for eid in self.eids:
            eids.append(self.Eid(eid))
        return eids

    def Eid(self, eid):
        if isinstance(eid, int):
            return eid
        return eid.eid

    def rawFields(self):
        eids = self.Eids()
        eids.sort()
        fields = ['QBDY3', self.sid, self.Q0, self.cntrlnd
                 ] + collapse_thru_by(eids)
        return fields

    def reprFields(self):
        cntrlnd = set_blank_if_default(self.cntrlnd, 0)
        eids = self.Eids()
        eids.sort()
        fields = ['QBDY3', self.sid, self.Q0, cntrlnd] + collapse_thru_by(eids)
        return fields


class QHBDY(ThermalLoad):
    """
    Defines a uniform heat flux into a set of grid points.
    """
    type = 'QHBDY'

    def __init__(self, card=None, data=None):
        ThermalLoad.__init__(self, card, data)

        if card:
            ## Load set identification number. (Integer > 0)
            self.sid = card.field(1)

            self.flag = card.field(2)
            assert self.flag in ['POINT', 'LINE', 'REV',
                                 'AREA3', 'AREA4', 'AREA6', 'AREA8']

            ## Magnitude of thermal flux into face. Q0 is positive for heat into the surface. (Real)
            self.Q0 = card.field(3)

            ## Area factor depends on type. (Real > 0.0 or blank)
            self.af = card.field(4)
            self.grids = card.fields(5)

            ## Grid point identification of connected grid points. (Integer > 0 or blank)
            self.grids = expand_thru_by(self.grids)
        else:
            self.sid = data[0]
            self.flag = data[1]
            self.Q0 = data[2]
            self.af = data[3]
            self.grids = data[4:]

    #def cross_reference(self,model):
    #    pass

    def rawFields(self):
        fields = ['QHBDY', self.sid, self.flag, self.Q0, self.af] + self.grids
        return fields

    def reprFields(self):
        return self.rawFields()


class TEMP(ThermalLoad):
    """
    Defines temperature at grid points for determination of thermal loading,
    temperature-dependent material properties, or stress recovery.
    """
    type = 'TEMP'

    def __init__(self, card=None, data=None):
        ThermalLoad.__init__(self, card, data)

        if card:
            ## Load set identification number. (Integer > 0)
            self.sid = card.field(1)

            fields = card.fields(2)
            nFields = len(fields)
            assert nFields % 2 == 0

            ## dictionary of temperatures where the key is the grid ID (Gi) and the value is the temperature (Ti)
            self.temperatures = {}
            for i in xrange(0, nFields, 2):
                self.temperatures[fields[i]] = fields[i + 1]
        else:
            #print "TEMP data = ",data
            self.sid = data[0]
            self.temperatures = {data[1]: data[2]}

    def add(self, tempObj):
        assert self.sid == tempObj.sid
        for (gid, temp) in self.tempObj.temperatures.iteritems():
            self.temperatures[gid] = temp

    def cross_reference(self, model):
        pass

    def rawFields(self):
        """Writes the TEMP card"""
        fields = ['TEMP', self.sid]

        nTemps = len(self.temperatures) - 1
        for i, (gid, temp) in enumerate(sorted(self.temperatures.iteritems())):
            fields += [gid, temp]
            if i % 3 == 2 and nTemps > i:  # start a new TEMP card
                fields += [None, 'TEMP', self.lid]
        return fields

    def reprFields(self):
        """Writes the TEMP card"""
        return self.rawFields()

# Loads
#-------------------------------------------------------
# Default Loads


class TEMPD(ThermalLoadDefault):
    """
    Defines a temperature value for all grid points of the structural model that have not
    been given a temperature on a TEMP entry
    """
    type = 'TEMPD'

    def __init__(self, card=None, data=None):
        ThermalLoadDefault.__init__(self, card, data)
        if card:
            fields = card.fields(1)
            nFields = len(fields)
            assert nFields % 2 == 0

            ## dictionary of temperatures where the key is the set ID (SIDi) and the value is the temperature (Ti)
            self.temperatures = {}
            for i in xrange(0, nFields, 2):
                self.temperatures[fields[i]] = fields[i + 1]
        else:
            self.temperatures = {data[0]: data[1]}

    def add(self, tempdObj):
        for (lid, tempd) in self.tempdObj.temperatures.iteritems():
            self.temperatures[lid] = tempd

    def cross_reference(self, model):
        pass

    def reprFields(self):
        """Writes the TEMPD card"""
        fields = ['TEMPD']

        nTemps = len(self.temperatures) - 1
        #print "self.temperatures = ",self.temperatures
        #print "nTemps = ",nTemps
        for i, (gid, temp) in enumerate(sorted(self.temperatures.iteritems())):
            fields += [gid, temp]
            if i % 4 == 3 and nTemps > i:  # start a new TEMP card
                fields += ['TEMPD']
        return fields
