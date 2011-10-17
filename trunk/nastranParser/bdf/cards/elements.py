#import sys
#from numpy import array,cross,dot
#from numpy import array
#from numpy.linalg import norm

# my code
#from bdf import printCard,setBlankIfDefault,setDefaultIfBlank
from AirQuick.general.generalMath import Area,Triangle_AreaCentroidNormal,Normal
from baseCard import Element

from elementsShell import *
from elementsBars  import *
        
class CSHEAR(Element):
    type = 'CSHEAR'
    def __init__(self,card):
        Element.__init__(self,card)
        nids = card.fields(3,7)
        self.prepareNodeIDs(nids)
        assert len(self.nodes)==4

    def __repr__(self):
        fields = [self.type,self.eid,self.pid]+self.nodes
        return self.printCard(fields)

class CRAC2D(Element):
    type = 'CRAC2D'
    def __init__(self,card):
        Element.__init__(self,card)

        nids = card.fields(3,21) # caps at 18
        self.prepareNodeIDs(nids)
        assert len(self.nodes)==18

    def __repr__(self):
        fields = [self.type,self.eid,self.pid]+self.nodes
        return self.printCard(fields)

class CRAC3D(Element):
    type = 'CRAC3D'
    def __init__(self,card):
        Element.__init__(self,card)

        nids = card.fields(3,67) # cap at +3 = 67
        self.prepareNodeIDs(nids)
        assert len(self.nodes)==64

    def __repr__(self):
        fields = [self.type,self.eid,self.pid]+self.nodes
        return printCard(fields)

        
class CVISC(CROD):
    type = 'CVISC'
    def __init__(self,card):
        CROD.__init__(self,card)
    ###
###

class CONM2(Element): # not done
    type = 'CONM2'
    # 'CONM2    501274  11064          132.274'
    def __init__(self,card):
        Element.__init__(self,card)
        #self.nids  = [ card[1] ]
        #del self.nids
        self.pid = None
        self.dunno = card.field(2)
        self.blank = card.field(3)
        self.mass  = card.field(4)
        
        #print "nids       = ",self.nids
        #print 'self.dunno = ',self.dunno
        #print 'self.blank = ',self.blank
        #print "mass       = ",self.mass
        #print "card       = ",card
        #print str(self)
        #sys.exit()
    
    def __repr__(self):
        fields = [self.type,self.eid,self.dunno,self.blank,self.mass]
        #fields = [self.type,self.eid,self.blank,self.mass]
        return printCard(fields)

   
