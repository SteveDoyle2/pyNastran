from baseCard import Element

class RigidElement(Element):
    def __repr__(self):
        fields = ['RIGID_EL',self.eid]
        return self.printCard(fields)
    def crossReference(self,mesh):
        pass

class RBAR(RigidElement):
    type = 'RBAR'
    def __init__(self,card=None,data=None):
        """
        RBAR EID GA GB CNA    CNB CMA CMB ALPHA
        RBAR 5   1   2 123456             6.5-6
        """
        RigidElement.__init__(self,card,data)
        if card:
            self.eid   = card.field(1)
            self.ga    = card.field(2)
            self.gb    = card.field(3)
            self.cna   = card.field(4)
            self.cnb   = card.field(5)
            self.cma   = card.field(6)
            self.cmb   = card.field(7)
            self.alpha = card.field(8)
        else:
            self.eid   = data[0]
            self.ga    = data[1]
            self.gb    = data[2]
            self.cna   = data[3]
            self.cnb   = data[4]
            self.cma   = data[5]
            self.cmb   = data[6]
            self.alpha = data[7]
        ###

    def __repr__(self):
        alpha = self.setDefaultIfBlank(self.alpha,0.0)
        fields = ['RBAR',self.eid,self.ga,self.gb,self.cna,self.cnb,self.cma,self.cmb,alpha]
        return self.printCard(fields)

class RBAR1(RigidElement):
    type = 'RBAR1'
    def __init__(self,card=None,data=None):
        """
        RBAR1 EID GA GB CB  ALPHA
        RBAR1 5    1  2 123 6.5-6
        """
        RigidElement.__init__(self,card,data)
        if card:
            self.eid   = card.field(1)
            self.ga    = card.field(2)
            self.gb    = card.field(3)
            self.alpha = card.field(4)
        else:
            self.eid   = data[0]
            self.ga    = data[1]
            self.gb    = data[2]
            self.alpha = data[3]
        ###

    def __repr__(self):
        alpha = self.setDefaultIfBlank(self.alpha,0.0)
        fields = ['RBAR1',self.eid,self.ga,self.gb,alpha]
        return self.printCard(fields)

class RBE1(RigidElement):  # maybe not done, needs testing
    type = 'RBE1'
    def __init__(self,card=None,data=None):
        RigidElement.__init__(self,card,data)
        self.eid = card.field(1)
        self.Gni = []
        self.Cni = []
        fields = card.fields(4)
        iUM = fields.index('UM')
        self.alpha = fields.pop() # the last field is not part of fields

        # loop till UM, no field9,field10
        while i<iUM-4:
            self.Gni.append(card.field(i ))
            self.Cni.append(card.field(i+1))
            if i%6==0:
                i+=2
            i+=2
        ###

        self.Gni = []
        self.Cni = []
        # loop till alpha, no field9,field10
        while i <card.nFields()-1: # dont grab alpha
            self.Gmi.append(card.field(i ))
            self.Cmi.append(card.field(i+1))
            if i%6==0:
                i+=2
            i+=2
        ###
        
    def __repr__(self):
        fields = [self.type,self.eid]
        for i,(gn,cn) in enumerate(zip(self.Gni,self.Cni)):
            fields+=[gn,cn]
            if i%6==0:
                fields += None
            ###

        fields += ['UM']
        for i,(gm,cm) in enumerate(zip(self.Gmi,self.Cmi)):
            fields+=[gm,cm]
            if i%6==0:
                fields += None
            ###
        fields += [self.alpha]
        return self.printCard(fields)

class RBE2(RigidElement):
    type = 'RBE2'
    def __init__(self,card=None,data=None):
        """
        RBE2 EID GN CM GM1 GM2 GM3 GM4 GM5
        GM6 GM7 GM8 -etc.- ALPHA
        """
        RigidElement.__init__(self,card,data)
        self.eid = card.field(1)
        self.gn  = card.field(2)
        self.cm  = card.field(3)
        self.Gmi = card.fields(4) # get the rest of the fields
        self.alpha = self.Gmi.pop() # the last field is not part of self.Gmi

    def __repr__(self):
        fields = [self.type,self.eid,self.gn,self.cm]+self.Gmi+[self.alpha]
        return self.printCard(fields)

class RBE3(RigidElement):  # not done, needs testing badly
    type = 'RBE3'
    def __init__(self,card=None,data=None):
        RigidElement.__init__(self,card,data)
        self.eid     = card.field(1)
        self.refgrid = card.field(3)
        self.refc    = card.field(4)
        #fields = card.fields(5)
        #iUM = fields.index('UM')
        
        fields = card.fields(5)
        try:
            iAlpha = fields.index('ALPHA')
        except ValueError:
            iAlpha = None

        try:
            iUm = fields.index('UM')
        except ValueError:
            iUm = None
        #print "iAlpha=%s iUm=%s" %(iUm,iAlpha)

        #print "iUM = ",iUM
        self.WtCG_groups = []
        if iUm:
            for i in range(5,card.nFields()):
                Gij = []

                wt = card.field(i)
                ci = card.field(i+1)
                i+=2
                g = 0
                while isinstance(gij,int):  # does this get extra fields???
                    gij = card.field(i+1)
                    Gij.append(gij)
                    i+=1
                print "gij_stop? = ",gij
                if gij=='UM':
                    print "breaking A..."
                    break
                ###
                self.WtCG_groups.append(wt,ci,Gij)
            ###
        
        self.Gmi = []
        self.Cmi = []
        ## thermal expansion coefficient
        self.alpha = 0.0
        if iAlpha:
            for j in range(i,card.nFields()):  # does this get extra fields???
                gmi = card.field(j)
                cmi = card.field(j+1)
                nextEntry = card.field(j+2)
                self.Gmi.append(gmi)
                self.Cmi.append(cmi)
                j+=2

                print "next_stop? = ",nextEntry
                if nextEntry=='ALPHA':
                    break
                ###
            ###
            self.alpha = card.field(j)

    def __repr__(self):
        fields = [self.type,self.eid,None,self.refc]
        for (wt,ci,Gij) in self.WtCG_groups:
            fields+=[wt,ci]+Gij
        fields += ['UM']
        for (gmi,cmi) in zip(self.Gmi,self.Cmi):
            fields+=[gmi,cmi]
        fields += ['ALPHA',self.alpha]
        return self.printCard(fields)
