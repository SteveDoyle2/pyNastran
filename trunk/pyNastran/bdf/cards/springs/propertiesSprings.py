import sys
from numpy import zeros,pi

# my code
from ..baseCard import Property

class SpringProperty(Property):
    type = 'SpringProperty'
    def __init__(self,card,data):
        Property.__init__(self,card,data)
        pass

class PELAS(SpringProperty):
    type = 'PELAS'
    def __init__(self,card=None,nPELAS=0,data=None):
        SpringProperty.__init__(self,card,data)
        nOffset = nPELAS*5
        if card:
            self.pid = card.field(1+nOffset) # 2 PELAS properties can be defined on 1 PELAS card
            self.k   = card.field(2+nOffset) # these are split into 2 separate cards
            self.ge  = card.field(3+nOffset)
            self.s   = card.field(4+nOffset)
        else:
            self.pid = data[0]
            self.k   = data[1]
            self.ge  = data[2]
            self.s   = data[3]
        ###
    def crossReference(self,model):
        pass

    def __repr__(self):
        fields = ['PELAS',self.pid,self.k,self.ge,self.s]
        return self.printCard(fields)
