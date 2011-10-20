from baseCard import BaseCard

class FLFACT(BaseCard):
    """
    FLFACT SID F1 F2 F3 F4 F5 F6 F7
    F8 F9 -etc.-
    
    FLFACT 97 .3 .7 3.5
    
    FLFACT SID F1 THRU FNF NF FMID       # delta quantity approach
    FLFACT 201 .200 THRU .100 11 .133333
    """
    type = 'FLFACT'
    def __init__(self,card):
        #Material.__init__(self,card)
        self.sid     = card.field(1)
        self.factors = card.fields(2)
        if self.factors[1]=='THRU':
            raise Exception('embedded THRUs not supported yet on FLFACT card\n')
            #(a,thru,b,n,dn) = factors
            #for i in range(

    def __repr__(self):
        fields = [self.type,self.sid]+self.factors
        return self.printCard(fields)

