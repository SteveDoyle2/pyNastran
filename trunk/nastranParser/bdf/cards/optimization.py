from baseCard import BaseCard

class OptConstraint(BaseCard):
    def __init__(self):
        pass

class DConstr(OptConstraint):
    self.type = 'DCONSTR'
    def __init__(self,card):
        self.id = card.field(1)
        self.rid = card.field(2)
        self.lid = card.field(3) or -1e20
        self.uid = card.field(4) or 1e20
        self.lowfq = card.field(5) or 0.0
        self.highfq = card.field(6) or 1e20
    
    def __repr__(self):
        #lid = 
        #uid = 
        #lowfq = 
        #highfq = 
        fields = [self.type,self.id,self.rid,lid,uid,lowfq,highfq]
        return self.printCard(fields)

class Desvar(OptConstraint):
    self.type = 'DESVAR'
    def __init__(self,card):
        self.id = card.field(1)
        self.label = card.field(2)
        self.xinit = card.field(3)
        self.xlb   = card.field(4) or -1e20
        self.xub   = card.field(5) or 1e20
        self.delx  = card.field(6) or 1e20
        self.ddval = card.field(7)
    
    def __repr__(self):
        #xlb = 
        #xub = 
        fields = [self.type,self.id,self.label,self.xinit,xlb,xub,
        self.delx,self.ddval]
        return self.printCard(fields)

class DDVal(OptConstraint):
    self.type = 'DDVal'
    def __init__(self,card):
        self.id = card.field(1)
        self.dval1 = card.field(2)
        self.dval2 = card.field(3)
        self.dval3 = card.field(4)
        self.dval4 = card.field(5)
        self.dval5 = card.field(6)
        self.dval6 = card.field(7)
        self.dval7 = card.field(8)
    
    def __repr__(self):
        fields = [self.type,self.id,self.dval1,self.dval2,self.dval3,
                  self.dval4,self.dval5,self.dval6,self.dval7]
        return self.printCard(fields)
