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
# pylint: disable=C0103,R0902,R0904,R0914



from pyNastran.bdf.cards.baseCard import Element


class DamperElement(Element):
    def __init__(self, card, data):
        Element.__init__(self, card, data)

class LineDamper(DamperElement):
    def __init__(self, card, data):
        DamperElement.__init__(self, card, data)

    def crossReference(self,mesh):
        self.nodes = mesh.Nodes(self.nodes)
        self.pid   = mesh.Property(self.pid)

class CVISC(LineDamper):
    type = 'CVISC'
    def __init__(self,card=None,data=None):
        LineDamper.__init__(self, card, data)
        if card:
            self.eid = int(card.field(1))
            self.pid = int(card.field(2,self.eid))
            nids = card.fields(3,5)
        else:
            self.eid = data[0]
            self.pid = data[1]
            nids = data[2:4]
        ###
        self.prepareNodeIDs(nids)
        assert len(self.nodes)==2
    ###

    def B(self):
        return self.pid.ce

    def rawFields(self):
        fields = ['CVISC',self.eid,self.Pid()]+self.nodeIDs()
        return fields

    def reprFields(self):
        return self.rawFields()
###

class CDAMP1(LineDamper):
    type = 'CDAMP1'
    def __init__(self,card=None,data=None):
        LineDamper.__init__(self, card, data)
        
        if card:
            self.eid = card.field(1)
            self.pid = card.field(2)

            nids = [card.field(3,0),card.field(5,0)]

            ## component number
            self.c1 = card.field(4,0)
            self.c2 = card.field(6,0)
        else:
            self.eid = data[0]
            self.pid = data[1]
            nids     = [data[2],data[4]]
            self.c1  = data[3]
            self.c2  = data[5]
        ###
        assert self.c1 in [0,1,2,3,4,5,6],'c1=|%s| on \n%s\n is invalid validComponents=[0,1,2,3,4,5,6]' %(str(self),self.c1)
        assert self.c2 in [0,1,2,3,4,5,6],'c2=|%s| on \n%s\n is invalid validComponents=[0,1,2,3,4,5,6]' %(str(self),self.c2)
        self.prepareNodeIDs(nids,allowEmptyNodes=True)
        assert len(self.nodes)==2

    def isSameCard(self, elem, debug=False):
        if self.type!=elem.type:  return False
        fields1 = [self.eid,self.Pid()]+self.nodeIDs()+[self.c1,self.c2]
        fields2 = [elem.eid,elem.Pid()]+elem.nodeIDs()+[elem.c1,elem.c2]
        if debug:
            print("fields1=%s fields2=%s" %(fields1, fields2))
        return self.isSameFields(fields1,fields2)

    def B(self):
        return self.pid.b

    def crossReference(self,model):
        self.nodes = model.Nodes(self.nodes,allowEmptyNodes=True)
        self.pid   = model.Property(self.pid)
        
    def rawFields(self):
        nodes = self.nodeIDs(allowEmptyNodes=True)
        fields = ['CDAMP1',self.eid,self.Pid(),nodes[0],self.c1,nodes[1],self.c2]
        return fields

class CDAMP2(LineDamper):
    type = 'CDAMP2'
    def __init__(self,card=None,data=None):
        LineDamper.__init__(self, card, data)
        
        if card:
            self.eid = card.field(1)

            ## Value of the scalar damper (Real)
            self.b   = card.field(2)

            nids = [card.field(3,0),card.field(5,0)]

            ## component number
            self.c1 = card.field(4,0)
            self.c2 = card.field(6,0)
        else:
            self.eid = data[0]
            self.b   = data[1]
            nids     = [data[2],data[4]]
            self.c1  = data[3]
            self.c2  = data[5]
        ###
        assert self.c1 in [0,1,2,3,4,5,6],'c1=|%s| on \n%s\n is invalid validComponents=[0,1,2,3,4,5,6]' %(str(self),self.c1)
        assert self.c2 in [0,1,2,3,4,5,6],'c2=|%s| on \n%s\n is invalid validComponents=[0,1,2,3,4,5,6]' %(str(self),self.c2)
        self.prepareNodeIDs(nids,allowEmptyNodes=True)
        assert len(self.nodes)==2

    def B(self):
        return self.b

    def crossReference(self,model):
        self.nodes = model.Nodes(self.nodes,allowEmptyNodes=True)
        
    def rawFields(self):
        nodes = self.nodeIDs(allowEmptyNodes=True)
        fields = ['CDAMP2',self.eid,self.b,nodes[0],self.c1,nodes[1],self.c2]
        return fields

class CDAMP3(LineDamper):
    type = 'CDAMP3'
    def __init__(self,card=None,data=None):
        LineDamper.__init__(self, card, data)
        
        if card:
            self.eid = card.field(1)
            self.pid = card.field(2)
            nids = [card.field(3,0),card.field(4,0)]
        else:
            self.eid = data[0]
            self.pid = data[1]
            nids     = [data[2],data[3]]
        ###
        self.prepareNodeIDs(nids,allowEmptyNodes=True)
        assert len(self.nodes)==2

    def B(self):
        return self.pid.b

    def crossReference(self,model):
        self.nodes = model.Nodes(self.nodes,allowEmptyNodes=True)
        self.pid   = model.Property(self.pid)
        
    def rawFields(self):
        nodes = self.nodeIDs(allowEmptyNodes=True)
        fields = ['CDAMP3',self.eid,self.pid,nodes[0],nodes[1]]
        return fields

class CDAMP4(LineDamper):
    type = 'CDAMP4'
    def __init__(self,card=None,data=None):
        LineDamper.__init__(self, card, data)
        
        if card:
            self.eid = card.field(1)
            ## Value of the scalar damper (Real)
            self.b   = card.field(2)
            nids = [card.field(3,0),card.field(4,0)]
        else:
            self.eid = data[0]
            self.b   = data[1]
            nids     = [data[2],data[3]]
        ###
        self.prepareNodeIDs(nids,allowEmptyNodes=True)
        assert len(self.nodes)==2

    def B(self):
        return self.b

    def crossReference(self,model):
        self.nodes = model.Nodes(self.nodes,allowEmptyNodes=True)
        
    def rawFields(self):
        nodes = self.nodeIDs(allowEmptyNodes=True)
        fields = ['CDAMP4',self.eid,self.b,nodes[0],nodes[1]]
        return fields

class CDAMP5(LineDamper):
    type = 'CDAMP5'
    def __init__(self,card=None,data=None):
        LineDamper.__init__(self, card, data)
        
        if card:
            self.eid = card.field(1)
            ## Property ID
            self.pid = card.field(2)
            nids = [card.field(3,0),card.field(4,0)]
        else:
            self.eid = data[0]
            self.pid = data[1]
            nids     = [data[2],data[3]]
        ###
        self.prepareNodeIDs(nids,allowEmptyNodes=True)
        assert len(self.nodes)==2

    def crossReference(self,model):
        self.nodes = model.Nodes(self.nodes,allowEmptyNodes=True)
        self.pid   = model.Property(self.pid)
        
    def rawFields(self):
        nodes = self.nodeIDs(allowEmptyNodes=True)
        fields = ['CDAMP5',self.eid,self.Pid(),nodes[0],nodes[1]]
        return fields

