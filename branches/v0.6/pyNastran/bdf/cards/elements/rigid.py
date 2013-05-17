## GNU Lesser General Public License
## 
## Program pyNastran - a python interface to NASTRAN files
## Copyright (C) 2011-2012  Steven Doyle, Al Danial
## 
## Authors and copyright holders of pyNastran
## Steven Doyle <mesheb82@gmail.com>
## Al Danial    <al.danial@gmail.com>
## 
## This file is part of pyNastran.
## 
## pyNastran is free software: you can redistribute it and/or modify
## it under the terms of the GNU Lesser General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
## 
## pyNastran is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU Lesser General Public License
## along with pyNastran.  If not, see <http://www.gnu.org/licenses/>.
## 
# pylint: disable=C0103,R0902,R0904,R0914,C0111
"""
All rigid elements are defined in this file.  This includes:

 * RBAR
 * RBAR1
 * RBE1
 * RBE2
 * RBE3

All rigid elements are RigidElement and Element objects.
"""
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
import sys
from itertools import izip, count

from pyNastran.bdf.fieldWriter import set_blank_if_default
from pyNastran.bdf.cards.baseCard import Element
from pyNastran.bdf.bdfInterface.assign_type import (integer,
    integer_or_double, integer_double_or_blank, integer_or_blank,
    double_or_blank, components, components_or_blank, blank, fields, string)
# integer_or_double, double,

class RigidElement(Element):
    def cross_reference(self, model):
        pass


class RBAR(RigidElement):
    type = 'RBAR'

    def __init__(self, card=None, data=None, comment=''):
        """
        ::

          RBAR EID GA GB CNA    CNB CMA CMB ALPHA
          RBAR 5   1   2 123456             6.5-6
        """
        RigidElement.__init__(self, card, data)
        if comment:
            self._comment = comment
        if card:
            self.eid = integer(card, 1, 'eid')
            self.ga = integer(card, 2, 'ga')
            self.gb = integer(card, 3, 'gb')
            self.cna = components_or_blank(card, 4, 'cna')
            self.cnb = components_or_blank(card, 5, 'cnb')
            self.cma = components_or_blank(card, 6, 'cma')
            self.cmb = components_or_blank(card, 7, 'cmb')
            self.alpha = double_or_blank(card, 8, 'alpha', 0.0)
            assert len(card) <= 9, 'len(RBAR card) = %i' % len(card)
        else:
            self.eid = data[0]
            self.ga = data[1]
            self.gb = data[2]
            self.cna = data[3]
            self.cnb = data[4]
            self.cma = data[5]
            self.cmb = data[6]
            self.alpha = data[7]

    # def convertToMPC(self, mpcID):
    #     """
    #     -Ai*ui + Aj*uj = 0
    #     where ui are the base DOFs (max=6)
    #     mpc sid   g1 c1 a1  g2 c2 a2
    #     rbe2 eid  gn cm g1  g2 g3 g4
    #     """
    #     raise NotImplementedError()
    #     #i = 0
    #     nCM = len(self.cm)
    #     Ai = nCM * len(self.Gmi) / len(self.gn)  # where nGN=1
    #
    #     card = ['MPC', mpcID]
    #     for cm in self.cm:  # the minus sign is applied to the base node
    #         card += [self.gn, cm, -Ai]
    #
    #     for gm in self.Gmi:
    #         for cm in self.cm:
    #             card += [gm, cm, Ai]
    #     return card

    #def writeCodeAster(self):
        #msg = ''
        #msg += "BLOCAGE=AFFE_CHAR_MECA(  # RBAR\n"
        #msg += "        MODELE=MODELE,\n"  # rigid element
        #msg += "        \n"
        #return msg

    def rawFields(self):
        list_fields = ['RBAR', self.eid, self.ga, self.gb, self.cna,
                  self.cnb, self.cma, self.cmb, self.alpha]
        return list_fields

    def reprFields(self):
        alpha = set_blank_if_default(self.alpha, 0.0)
        list_fields = ['RBAR', self.eid, self.ga, self.gb, self.cna, self.cnb,
                  self.cma, self.cmb, alpha]
        return list_fields


class RBAR1(RigidElement):
    type = 'RBAR1'

    def __init__(self, card=None, data=None, comment=''):
        """
        ::

          RBAR1 EID GA GB CB  ALPHA
          RBAR1 5    1  2 123 6.5-6
        """
        RigidElement.__init__(self, card, data)
        if comment:
            self._comment = comment
        if card:
            self.eid = integer(card, 1, 'eid')
            self.ga = integer(card, 2, 'ga')
            self.gb = integer(card, 3, 'gb')
            self.cb = components_or_blank(card, 4, 'cb')
            self.alpha = double_or_blank(card, 5, 'alpha', 0.0)
            assert len(card) <= 6, 'len(RBAR1 card) = %i' % len(card)
        else:
            self.eid = data[0]
            self.ga = data[1]
            self.gb = data[2]
            self.cb = data[3]
            self.alpha = data[4]

    def rawFields(self):
        list_fields = ['RBAR1', self.eid, self.ga, self.gb, self.cb, self.alpha]
        return list_fields

    def reprFields(self):
        alpha = set_blank_if_default(self.alpha, 0.0)
        list_fields = ['RBAR1', self.eid, self.ga, self.gb, self.cb, alpha]
        return list_fields


class RBE1(RigidElement):  # maybe not done, needs testing
    type = 'RBE1'

    def __init__(self, card=None, data=None, comment=''):
        RigidElement.__init__(self, card, data)
        if comment:
            self._comment = comment

        self.eid = integer(card, 1, 'eid')
        self.Gni = []
        self.Cni = []

        #fields = card[2:]
        iUm = card.index('UM')
        if iUm > 0:
            assert string(card, iUm, 'UM') == 'UM'

        if isinstance(card[-1], float):
            self.alpha = card.fields[-1].pop()  # the last field is not part of fields
            #nfields = len(card) - 1
        else:
            #nfields = len(card)
            self.alpha = 0.

        # loop till UM, no field9,field10
        #print("iUm = %s" % iUm)
        
        n = 1
        i = 0
        offset = 2
        while offset + i < iUm - 1:
            #print('field(%s) = %s' % (offset + i, card.field(offset + i)))
            gni = integer_or_blank(card, offset + i, 'gn%i' % n)
            cni = components_or_blank(card, offset + i + 1, 'cn%i' % n)

            if gni:
                #print("gni=%s cni=%s" % (gni ,cni))
                self.Gni.append(gni)
                self.Cni.append(cni)
                n += 1
            else:
                assert cni is None
            i += 2
        
        #print('Gni =', self.Gni)
        #print('Cni =', self.Cni)
        self.Gmi = []
        self.Cmi = []

        # loop till alpha, no field9,field10
        n = 1
        offset = iUm + 1
        i = 0
        while offset + i < len(card):  # dont grab alpha
            gmi = integer_or_blank(card, offset + i, 'gm%i' % n)
            cmi = components_or_blank(card, offset + i + 1, 'cm%i' % n)
            if gmi:
                #print("gmi=%s cmi=%s" % (gmi ,cmi))
                self.Gmi.append(gmi)
                self.Cmi.append(cmi)
                n += 1
            else:
                assert cmi is None
            i += 2
        #print('Gmi =', self.Gmi)
        #print('Cmi =', self.Cmi)
        #print(self)
        #sys.exit()

    def rawFields(self):
        list_fields = [self.type, self.eid]

        for (i, gn, cn) in izip(count(), self.Gni, self.Cni):
            #print('i=%r gn=%r cn=%r' % (i, gn, cn))
            list_fields += [gn, cn]
            if i > 0 and i % 3 == 0:
                #print('adding blank')
                list_fields += [None]

        nSpaces = 8 - (len(list_fields) - 1) % 8  # puts UM/ALPHA onto next line
        if nSpaces < 8:
            list_fields += [None] * nSpaces

        # overly complicated loop to print the UM section
        list_fields += ['UM']
        j = 1
        for (i, gm, cm) in izip(count(), self.Gmi, self.Cmi):
            #print "j=%s gmi=%s cmi=%s" %(j,gm,cm)
            list_fields += [gm, cm]
            if i > 0 and j % 3 == 0:
                list_fields += [None, None]
                #print "---"
                j -= 3
            j += 1
        
        if self.alpha > 0.:  # handles default alpha value
            nSpaces = 8 - (len(list_fields) - 1) % 8  # puts ALPHA onto next line
            if nSpaces == 1:
                list_fields += [None, None]
            list_fields += [self.alpha]
        return list_fields

    def reprFields(self):
        return self.rawFields()


class RBE2(RigidElement):
    type = 'RBE2'

    def __init__(self, card=None, data=None, comment=''):
        """
        +-------+-----+-----+-----+------+-------+-----+-----+-----+
        |   1   |  2  |  3  |  4  |  5   |   6   |  7  |  8  |  9  |
        +-------+-----+-----+-----+------+-------+-----+-----+-----+
        |  RBE2 | EID | GN  | CM  | GM1  | GM2   | GM3 | GM4 | GM5 |
        +-------+-----+-----+-----+------+-------+-----+-----+-----+
        |       | GM6 | GM7 | GM8 | etc. | ALPHA |
        +-------+-----+-----+-----+------+-------+
        """
        RigidElement.__init__(self, card, data)
        if comment:
            self._comment = comment
        if card:
            #: Element identification number
            self.eid = integer(card, 1, 'eid')

            #: Identification number of grid point to which all six independent
            #: degrees-of-freedom for the element are assigned. (Integer > 0)
            self.gn = integer(card, 2, 'gn')

            #: Component numbers of the dependent degrees-of-freedom in the
            #: global coordinate system at grid points GMi. (Integers 1 through
            #: 6 with no embedded blanks.)
            self.cm = components_or_blank(card, 3, 'cm')

            alpha = integer_or_double(card, len(card) - 1, 'alpha')
            if isinstance(alpha, float):
                #: Grid point identification numbers at which dependent
                #: degrees-of-freedom are assigned. (Integer > 0)
                self.alpha = alpha
                
                # the last field is not part of Gmi
                n = 1
            else:
                # the last field is part of Gmi
                n = 0
                self.alpha = 0.0

            j=4
            self.Gmi = []
            for i in range(len(card)-4-n):
                gmi = integer(card, j + i, 'Gm%i' % (i + 1))
                #print('gm%i = %s' % (i + 1, gmi))
                self.Gmi.append(gmi)
        else:
            self.eid = data[0]
            self.gn = data[1]
            self.cm = data[2]
            self.Gmi = data[3]
            self.alpha = data[4]
            print("eid=%s gn=%s cm=%s Gmi=%s alpha=%s"
                  % (self.eid, self.gn, self.cm, self.Gmi, self.alpha))
            raise NotImplementedError('RBE2 data...')

        assert self.gn is not None, 'gm=%s' % self.gm
        assert self.cm is not None, 'cm=%s' % self.cm
        self.gn = str(self.gn)
        self.cm = str(self.cm)

    def convertToMPC(self, mpcID):
        """
        .. math:: -A_i u_i + A_j u_j = 0

        where :math:`u_i` are the base DOFs (max=6)::
        
          mpc sid   g1 c1 a1  g2 c2 a2
          rbe2 eid  gn cm g1  g2 g3 g4
        """
        #i = 0
        nCM = len(self.cm)
        Ai = nCM * len(self.Gmi) / len(self.gn)  # where nGN=1

        card = ['MPC', mpcID]
        for cm in self.cm:  # the minus sign is applied to the base node
            card += [self.gn, cm, -Ai]

        for gm in self.Gmi:
            for cm in self.cm:
                card += [gm, cm, Ai]
        return card

    def writeCodeAster(self):
        """
        Converts to a LIAISON SOLIDE for dofs 123456.
        For other dof combinations, general MPC equations are written
        """
        msg = ''
        msg += "BLOCAGE=AFFE_CHAR_MECA(  # RBE2 ID=%s\n" % (self.eid)
        msg += "        MODELE=MODELE,\n"  # rigid element
        if self.cm == 123456:
            msg += "        LIASON_SOLIDE=(\n"
            msg += "        _F(NOEUD=\n"
            msg += "           "
            for nid in self.Gmi:
                msg += "'N%i'," % (nid)
            msg = msg[:-1]
            msg += '\n'
        else:
            msg += "        _F(NOEUD=  # doesnt handle coordinate systems\n"
            msg += "           "
            for nid in self.Gmi:
                msg += "'N%i'," % (nid)
            msg = msg[:-1]
            msg += '\n'

            #msg += "        \n"
            #msg += "        \n"
        #msg += "        \n"
        #msg += "        \n"
        #msg += "        \n"
        return msg

    def rawFields(self):
        list_fields = ['RBE2', self.eid, self.gn, self.cm] + self.Gmi + [self.alpha]
        return list_fields

    def reprFields(self):
        alpha = set_blank_if_default(self.alpha, 0.)
        list_fields = ['RBE2', self.eid, self.gn, self.cm] + self.Gmi + [alpha]
        return list_fields


class RBE3(RigidElement):
    """
    .. todo:: not done, needs testing badly
    """
    type = 'RBE3'

    def __init__(self, card=None, data=None, comment=''):
        """
        eid
        refgrid
        refc
        WtCG_groups = [wt,ci,Gij]
        Gmi
        Cmi
        alpha
        """
        RigidElement.__init__(self, card, data)
        if comment:
            self._comment = comment
        self.eid = integer(card, 1, 'eid')
        blank(card, 2, 'blank')
        self.refgrid = integer(card, 3, 'refgrid')
        self.refc = components_or_blank(card, 4, 'refc')
        #iUM = fields.index('UM')

        fields = card[5:]
        #print "fields = ",fields
        iOffset = 5
        iWtMax = len(fields) + iOffset
        try:
            iAlpha = fields.index('ALPHA') + iOffset
            iWtMax = iAlpha  # the index to start parsing UM
            iUmStop = iAlpha  # the index to stop  parsing UM
        except ValueError:
            iAlpha = None
            iUmStop = iWtMax
        #print "iAlpha = ",iAlpha
        try:
            iUm = fields.index('UM') + iOffset
            iWtMax = iUm
        except ValueError:
            iUm = None
        #print "iAlpha=%s iUm=%s" %(iAlpha,iUm)
        #print "iAlpha=%s iWtMax=%s" %(iAlpha,iWtMax)

        #print "iUM = ",iUM
        self.WtCG_groups = []

        i = iOffset
        n = 1
        while i < iWtMax:
            Gij = []
            wtname = 'wt' + str(n)
            wt = double_or_blank(card, i, wtname)
            if wt is not None:
                cname = 'c'+str(n)
                ci = components_or_blank(card, i + 1, cname)
                
                #print("%s=%s %s=%s" % (wtname, wt, cname, ci))
                i += 2
                gij = 0
                
                j = 0
                while isinstance(gij, int) and i < iWtMax:
                    j += 1
                    gij_name = 'g%s,%s' % (n, j)
                    gij = integer_double_or_blank(card, i, gij_name)
                    if isinstance(gij, float):
                        break
                    #print("%s = %s" % (gij_name, gij))
                    if gij is not None:
                        Gij.append(gij)
                    i += 1
                wtCG_group = [wt, ci, Gij]
                self.WtCG_groups.append(wtCG_group)
                #print('----finished a group=%r----' % wtCG_group)
            else:
                i += 1

        self.Gmi = []
        self.Cmi = []
        #print ""
        if iUm:
            #print('UM =', card.field(iUm))  # UM
            i = iUm + 1
            n = 1
            #print("i=%s iUmStop=%s" % (i,iUmStop))
            for j in xrange(i, iUmStop, 2):
                
                gm_name = 'gm' + str(n)
                cm_name = 'cm' + str(n)
                gmi = integer_or_blank(card, j, gm_name)
                if gmi is not None:
                    cmi = components(card, j + 1, cm_name)
                    #print "gmi=%s cmi=%s" %(gmi,cmi)
                    self.Gmi.append(gmi)
                    self.Cmi.append(cmi)

        if iAlpha:
            self.alpha = double_or_blank(card, iAlpha + 1, 'alpha')
        else:
            #: thermal expansion coefficient
            self.alpha = 0.0
        #print self

    # def convertToMPC(self, mpcID):
    #     """
    #     -Ai*ui + Aj*uj = 0
    #     where ui are the base DOFs (max=6)
    #     mpc sid   g1 c1 a1  g2 c2 a2
    #     rbe2 eid  gn cm g1  g2 g3 g4
    #     """
    #     raise NotImplementedError('this is the code for an RBE2...not RBE3')
    #     #i = 0
    #     nCM = len(self.cm)
    #     Ai = nCM * len(self.Gmi) / len(self.gn)  # where nGN=1
    #
    #     card = ['MPC', mpcID]
    #     for cm in self.cm:  # the minus sign is applied to the base node
    #         card += [self.gn, cm, -Ai]
    #
    #     for gm in self.Gmi:
    #         for cm in self.cm:
    #             card += [gm, cm, Ai]
    #     return card

    def rawFields(self):
        list_fields = ['RBE3', self.eid, None, self.refgrid, self.refc]
        for (wt, ci, Gij) in self.WtCG_groups:
            #print 'wt=%s ci=%s Gij=%s' %(wt,ci,Gij)
            list_fields += [wt, ci] + Gij
        nSpaces = 8 - (len(list_fields) - 1) % 8  # puts UM onto next line

        if nSpaces < 8:
            list_fields += [None] * nSpaces

        if self.Gmi and 0:
            fields2 = ['UM']
            for (gmi, cmi) in izip(self.Gmi, self.Cmi):
                fields2 += [gmi, cmi]
            list_fields += self.buildTableLines(fields2, i=1, j=1)

        if self.Gmi:
            list_fields += ['UM']
        if self.Gmi:
            #print "Gmi = ",self.Gmi
            #print "Cmi = ",self.Cmi
            for (gmi, cmi) in izip(self.Gmi, self.Cmi):
                list_fields += [gmi, cmi]

        nSpaces = 8 - (len(list_fields) - 1) % 8  # puts ALPHA onto next line
        if nSpaces < 8:
            list_fields += [None] * nSpaces

        if self.alpha > 0.:  # handles the default value
            list_fields += ['ALPHA', self.alpha]
        return list_fields

    def reprFields(self):
        return self.rawFields()