from pyNastran.bdf.cards.baseCard import Property

class BushingProperty(Property):
    type = 'BushingProperty'
    def __init__(self,card,data):
        Property.__init__(self,card,data)
        pass
    def crossReference(self,model):
        pass

class PBUSH(BushingProperty):
    type = 'PBUSH'
    def __init__(self,card=None,data=None):
        BushingProperty.__init__(self,card,data)
        if card:
            ## Property ID
            self.pid = card.field(1)
            
            nFields = card.nFields()
            self.vars = []
            iStart = 2
            while iStart<nFields:
                pname = card.field(iStart)
                if   pname=='K':   self.getK(card,iStart)
                elif pname=='B':   self.getB(card,iStart)
                elif pname=='GE':  self.getGE(card,iStart)
                elif pname=='RCV': self.getRCV(card,iStart)
                else:
                    break
                iStart += 8
            ###
        else:
            self.pid = data[0]
            self.b   = data[1]
            raise NotImplementedError('PBUSH data...')
        ###
        #print self
        #sys.exit()

    def getK(self,card,iStart):
        ## Flag indicating that the next 1 to 6 fields are stiffness values in the element coordinate system.
        #self.k = card.field(iStart)
        ## Nominal stiffness values in directions 1 through 6. See Remarks 2. and 3. (Real; Default = 0.0)
        self.Ki = card.fields(i=iStart+1,j=iStart+6)
        #print "Ki = ",self.Ki
        self.vars.append('K')

    def getB(self,card,iStart):
        ## Flag indicating that the next 1 to 6 fields are force-per-velocity damping.
        #self.b = card.field(iStart)
        ## Force per unit velocity (Real)
        ## Nominal damping coefficients in direction 1 through 6 in units of force per
        ## unit velocity. See Remarks 2., 3., and 9. (Real; Default = 0.0)
        self.Bi = card.fields(i=iStart+1,j=iStart+6)
        self.vars.append('B')

    def getGE(self,card,iStart):
        ## Flag indicating that the next fields, 1 through 6 are structural damping constants. See Remark 7. (Character)
        #self.ge = card.field(iStart)
        ## Nominal structural damping constant in directions 1 through 6. See
        ## Remarks 2. and 3. (Real; Default = 0.0)
        self.GEi = card.fields(i=iStart+1,j=iStart+6)
        self.vars.append('GE')

    def getRCV(self,card,iStart):
        ## Flag indicating that the next 1 to 4 fields are stress or strain coefficients. (Character)
        #self.ge = card.field(iStart)
        self.sa = card.field(iStart+1,1.)
        self.st = card.field(iStart+2,1.)
        self.ea = card.field(iStart+3,1.)
        self.et = card.field(iStart+4,1.)
        self.vars.append('RCV')

    def rawFields(self):
        fields = ['PBUSH',self.pid]
        for var in self.vars:
            if var=='K':
                fields += ['K']+self.Ki
            elif var=='B':
                fields += ['B']+self.Bi
            elif var=='GE':
                fields += ['GE']+self.GEi
            elif var=='RCV':
                fields += ['RCV',self.sa,self.st,self.ea,self.et]
            else:
                raise Exception('not supported PBUSH field...')
            nSpaces = 8-(len(fields)-1)%8
            
            #print "nSpaces = ",nSpaces
            if nSpaces<8:
                fields += [None]*(nSpaces+1)
            ###
        return fields

    def reprFields(self):
        return self.rawFields()

#class PBUSH1D
#class PBUSH2D
#class PBUSHT
