# pylint: disable=C0103,R0902,R0904,R0914
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)

import sys
#from numpy import zeros,pi

from pyNastran.bdf.fieldWriter import set_blank_if_default
from pyNastran.bdf.cards.baseCard import Property


class SpringProperty(Property):
    type = 'SpringProperty'

    def __init__(self, card, data):
        Property.__init__(self, card, data)
        pass


class PELAS(SpringProperty):
    type = 'PELAS'

    def __init__(self, card=None, nPELAS=0, data=None):
        SpringProperty.__init__(self, card, data)
        nOffset = nPELAS * 5
        if card:
            self.pid = card.field(1 + nOffset)  # 2 PELAS properties can be defined on 1 PELAS card
            self.k   = card.field(2+nOffset) # these are split into 2 separate cards
            self.ge = card.field(3 + nOffset, 0.)
            self.s = card.field(4 + nOffset, 0.)
        else:
            self.pid = data[0]
            self.k = data[1]
            self.ge = data[2]
            self.s = data[3]
        ###

    def cross_reference(self, model):
        #if self.sol in [108,129]:
            #self.pid = self.pelasts[self.pid]
        pass

    def writeCodeAster(self):
        """
        @todo check if there are 1 (DISCRET=>K_T_D_N) or 2 (DISCRET_2D=>K_T_D_L) nodes
        """
        nodes = self.nodeIDs()
        msg = ''
        msg += 'DISCRET=_F( # PELAS\n'
        if nodes[0]:
            msg += "     CARA='K_T_D_N'\n"
            msg += "     GROUP_MA=P_%s\n" % (self.Pid())
            msg += "     NOEUD=N%i,\n" % (nodes[0])

        if nodes[1]:
            msg += "     CARA='K_T_D_L'\n"
            msg += "     NOEUD=N%i,\n" % (nodes[1])
            msg += "     AMOR_HYST=%g # ge - damping\n" % (self.ge)
        msg += "     )\n"
        msg += "\n"

        if self.c1 == 1:  # @todo what is this???
            msg += "VALE=(%g,0.,0.)\n" % (self.k)
        elif self.c1 == 2:
            msg += "VALE=(0.,%g,0.)\n" % (self.k)
        elif self.c1 == 2:
            msg += "VALE=(0.,0.,%g)\n" % (self.k)
        else:
            raise ValueError('unsupported value of c1=%s' % (self.c1))
        ###
        return msg

    def rawFields(self):
        fields = ['PELAS', self.pid, self.k, self.ge, self.s]
        return fields

    def reprFields(self):
        ge = set_blank_if_default(self.ge, 0.)
        s = set_blank_if_default(self.s, 0.)
        fields = ['PELAS', self.pid, self.k, ge, s]
        return fields


class PELAST(SpringProperty):
    """
    Frequency Dependent Elastic Property
    Defines the frequency dependent properties for a PELAS Bulk Data entry.

    The PELAST entry is ignored in all solution sequences except frequency
    response (108) or nonlinear analyses (129).
    """
    type = 'PELAST'

    def __init__(self, card=None, nPELAS=0, data=None):
        SpringProperty.__init__(self, card, data)
        self.pid = card.field(1)
        ## Identification number of a TABLEDi entry that defines the force per unit
        ## displacement vs. frequency relationship. (Integer > 0; Default = 0)
        self.tkid = card.field(2, 0)
        ## Identification number of a TABLEDi entry that defines the
        ## nondimensional structural damping coefficient vs. frequency
        ## relationship. (Integer > 0; Default = 0)
        self.tgeid = card.field(3, 0)
        ## Identification number of a TABELDi entry that defines the nonlinear
        ## force vs. displacement relationship. (Integer > 0; Default = 0)
        self.tknid = card.field(4, 0)

    def cross_reference(self, model):
        self.pid = model.Property(self.pid)
        if self.tkid > 0:
            self.tkid = model.Table(self.tkid)
        if self.tgeid > 0:
            self.tgeid = model.Table(self.tgeid)
        if self.tknid > 0:
            self.tknid = model.Table(self.tknid)

    def Pid(self):
        if isinstance(self.pid, int):
            return self.pid
        return self.pid.pid

    def Tkid(self):
        if isinstance(self.tkid, int):
            return self.tkid
        return self.tkid.tid

    def Tknid(self):
        if isinstance(self.tknid, int):
            return self.tknid
        return self.tknid.tid

    def Tgeid(self):
        if isinstance(self.tgeid, int):
            return self.tgeid
        return self.tgeid.tid
