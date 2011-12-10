from baseCard import BaseCard

class OptConstraint(BaseCard):
    def __init__(self):
        pass

class DCONSTR(OptConstraint):
    type = 'DCONSTR'
    def __init__(self,card=None,data=None):
        if card:
            self.oid    = card.field(1)
            self.rid    = card.field(2)
            self.lid    = card.field(3,-1e20)
            self.uid    = card.field(4, 1e20)
            self.lowfq  = card.field(5,0.0)
            self.highfq = card.field(6,1e20)
        else:
            self.oid    = data[0]
            self.rid    = data[1]
            self.lid    = data[2]
            self.uid    = data[3]
            self.lowfq  = data[4]
            self.highfq = data[5]
        ###

    def __repr__(self):
        lid    = self.setBlankIfDefault(self.lid,-1e20)
        uid    = self.setBlankIfDefault(self.uid, 1e20)
        lowfq  = self.setBlankIfDefault(self.lowfq,0.0)
        highfq = self.setBlankIfDefault(self.highfq,1e20)
        fields = ['DCONSTR',self.oid,self.rid,lid,uid,lowfq,highfq]
        return self.printCard(fields)

class DESVAR(OptConstraint):
    type = 'DESVAR'
    def __init__(self,card):
        self.oid = card.field(1)
        self.label = card.field(2)
        self.xinit = card.field(3)
        self.xlb   = card.field(4,-1e20)
        self.xub   = card.field(5, 1e20)
        self.delx  = card.field(6, 1e20)
        self.ddval = card.field(7)
    
    def __repr__(self):
        #xlb = 
        #xub = 
        fields = ['DESVAR',self.oid,self.label,self.xinit,xlb,xub,
        self.delx,self.ddval]
        return self.printCard(fields)

class DDVAL(OptConstraint):
    type = 'DDVal'
    def __init__(self,card):
        self.oid = card.field(1)
        self.dval1 = card.field(2)
        self.dval2 = card.field(3)
        self.dval3 = card.field(4)
        self.dval4 = card.field(5)
        self.dval5 = card.field(6)
        self.dval6 = card.field(7)
        self.dval7 = card.field(8)
    
    def __repr__(self):
        fields = ['DDVAL',self.oid,self.dval1,self.dval2,self.dval3,
                  self.dval4,self.dval5,self.dval6,self.dval7]
        return self.printCard(fields)
