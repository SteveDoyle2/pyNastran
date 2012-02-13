import sys
from numpy import zeros,pi

# pyNastran
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

    def writeCodeAster(self):
        """
        %PELAS - check if there are 1 (DISCRET=>K_T_D_N) or 2 (DISCRET_2D=>K_T_D_L) nodes
        """
        nodes = self.nodeIDs()
        msg = ''
        msg += 'DISCRET=_F( # PELAS\n'
        if nodes[0]:
            msg += "     CARA='K_T_D_N'\n"
            msg += "     GROUP_MA=P_%s\n" %(self.Pid())
            msg += "     NOEUD=N%i,\n" %(nodes[0])

        if nodes[1]:
            msg += "     CARA='K_T_D_L'\n"
            msg += "     NOEUD=N%i,\n" %(nodes[1])
            msg += "     AMOR_HYST=%g # ge - damping\n" %(self.ge)
        msg += "     )\n"
        msg += "\n"
        
        if self.c1==1:
            msg += "VALE=(%g,0.,0.)\n" %(self.k)
        elif self.c1==2:
            msg += "VALE=(0.,%g,0.)\n" %(self.k)
        elif self.c1==2:
            msg += "VALE=(0.,0.,%g)\n" %(self.k)
        else:
            raise Exception('unsupported value of c1=%s' %(self.c1))
        ###
        return msg

    def rawFields(self):
        fields = ['PELAS',self.pid,self.k,self.ge,self.s]
        return fields

    def reprFields(self):
        return self.rawFields()
