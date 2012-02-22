from baseCard import Element
import sys

class RigidElement(Element):
    #def __repr__(self):
        #fields = [self.type,self.eid]
        #return self.printCard(fields)
    def crossReference(self,model):
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
            self.alpha = card.field(8,0.0)
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

    #def writeCodeAster(self):
        #msg = ''
        #msg += "BLOCAGE=AFFE_CHAR_MECA(  # RBAR\n"
        #msg += "        MODELE=MODELE,\n"  # rigid element
        #msg += "        \n"
        #return msg
        
    def rawFields(self):
        fields = ['RBAR',self.eid,self.ga,self.gb,self.cna,self.cnb,self.cma,self.cmb,self.alpha]
        return fields

    def reprFields(self):
        alpha = self.setBlankIfDefault(self.alpha,0.0)
        fields = ['RBAR',self.eid,self.ga,self.gb,self.cna,self.cnb,self.cma,self.cmb,alpha]
        return fields

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
            self.cb    = card.field(4)
            self.alpha = card.field(5)
        else:
            self.eid   = data[0]
            self.ga    = data[1]
            self.gb    = data[2]
            self.cb    = data[3]
            self.alpha = data[4]
        ###

    def rawFields(self):
        fields = ['RBAR1',self.eid,self.ga,self.gb,self.cb,self.alpha]
        return fields

    def reprFields(self):
        alpha = self.setBlankIfDefault(self.alpha,0.0)
        fields = ['RBAR1',self.eid,self.ga,self.gb,self.cb,alpha]
        return fields

class RBE1(RigidElement):  # maybe not done, needs testing
    type = 'RBE1'
    def __init__(self,card=None,data=None):
        RigidElement.__init__(self,card,data)
        self.eid = card.field(1)
        self.Gni = []
        self.Cni = []

        fields = card.fields(2)
        iUm = fields.index('UM')+2
        if isinstance(fields[-1],float):
            self.alpha = fields.pop() # the last field is not part of fields
            nFields = card.nFields()-1
        else:
            nFields = card.nFields()
            self.alpha = 0.

        # loop till UM, no field9,field10
        i=2
        #print "iUm = %s" %(iUm)
        while i<iUm:
            gni = card.field(i  )
            #if gni:
            cni = card.field(i+1)
            self.Gni.append(gni)
            self.Cni.append(cni)
            #print "gni=%s cni=%s" %(gni,cni)
            if i%6==0:
                i+=2
            i+=2
        ###

        self.Gmi = []
        self.Cmi = []
        #print ""
        #print "i=%s iUm=%s card.field(iUm)=%s" %(i,iUm,card.field(iUm))
        # loop till alpha, no field9,field10
        while i < nFields: # dont grab alpha
            gmi = card.field(i  )
            cmi = card.field(i+1)
            if gmi:
                #print "gmi=%s cmi=%s" %(gmi,cmi)
                self.Gmi.append(gmi)
                self.Cmi.append(cmi)
            #else:
                #print "---"
            #if i%8==0:
            #    i+=2
            i+=2
        ###
        #print self
        #sys.exit()
        
    def rawFields(self):
        fields = [self.type,self.eid]
        for i,(gn,cn) in enumerate(zip(self.Gni,self.Cni)):
            fields+=[gn,cn]
            if i%6==0:
                fields += [None]
            ###
        ###
        nSpaces = 8-(len(fields)-1)%8  # puts ALPHA onto next line
        if nSpaces<8:
            fields += [None]*nSpaces

        ## overly complicated loop to print the UM section
        fields += ['UM']
        j=1
        for i,(gm,cm) in enumerate(zip(self.Gmi,self.Cmi)):
            #print "j=%s gmi=%s cmi=%s" %(j,gm,cm)
            fields+=[gm,cm]
            if i>0 and j%3==0:
                fields += [None,None]
                print "---"
                j-=3
            j+=1
            ###
        ###
        nSpaces = 8-(len(fields)-1)%8  # puts ALPHA onto next line
        if nSpaces==1:
            fields += [None,None]
        if self.alpha>0.: # handles default alpha value
            fields += [self.alpha]
        return fields

    def reprFields(self):
        return self.rawFields()

class RBE2(RigidElement):
    type = 'RBE2'
    def __init__(self,card=None,data=None):
        """
        RBE2 EID GN CM GM1 GM2 GM3 GM4 GM5
        GM6 GM7 GM8 -etc.- ALPHA
        """
        RigidElement.__init__(self,card,data)
        if card:
            ## Element identification number
            self.eid = card.field(1)
            ## Identification number of grid point to which all six independent
            ## degrees-of-freedom for the element are assigned. (Integer > 0)
            self.gn  = card.field(2)
            ## Component numbers of the dependent degrees-of-freedom in the global
            ## coordinate system at grid points GMi. (Integers 1 through 6 with no
            ## embedded blanks.)
            self.cm  = card.field(3)  # 123456 constraint or other
            ## Grid point identification numbers at which dependent degrees-offreedom
            ## are assigned. (Integer > 0)
            self.Gmi = card.fields(4) # get the rest of the fields
            if len(self.Gmi)>0 and isinstance(self.Gmi[-1],float):
                ## Thermal expansion coefficient. See Remark 11. (Real > 0.0 or blank)
                self.alpha = self.Gmi.pop() # the last field is not part of self.Gmi
            else:
                self.alpha = 0.0
            ###
        else:
            raise NotImplementedError('RBE2 data...')
        ###

    def writeCodeAster(self):
        """
        Converts to a LIAISON SOLIDE for dofs 123456.
        For other dof combinations, general MPC equations are written
        """
        msg = ''
        msg += "BLOCAGE=AFFE_CHAR_MECA(  # RBE2 ID=%s\n" %(self.eid)
        msg += "        MODELE=MODELE,\n"  # rigid element
        if self.cm==123456:
            msg += "        LIASON_SOLIDE=(\n"
            msg += "        _F(NOEUD=\n"
            msg += "           "
            for nid in self.Gmi:
                msg += "'N%i'," %(nid)
            msg = msg[:-1]
            msg += '\n'
        else:
            msg += "        _F(NOEUD=  # doesnt handle coordinate systems\n"
            msg += "           "
            for nid in self.Gmi:
                msg += "'N%i'," %(nid)
            msg = msg[:-1]
            msg += '\n'
            
            #msg += "        \n"
            #msg += "        \n"
        #msg += "        \n"
        #msg += "        \n"
        #msg += "        \n"
        return msg
        
    def rawFields(self):
        fields = ['RBE2',self.eid,self.gn,self.cm]+self.Gmi+[self.alpha]
        return fields

    def reprFields(self):
        alpha = self.setBlankIfDefault(self.alpha,0.)
        fields = ['RBE2',self.eid,self.gn,self.cm]+self.Gmi+[alpha]
        return fields

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
        #print "fields = ",fields
        iOffset = 5
        iWtMax = len(fields)+iOffset
        try:
            iAlpha = fields.index('ALPHA')+iOffset
            iWtMax  = iAlpha  # the index to start parsing UM
            iUmStop = iAlpha  # the index to stop  parsing UM
        except ValueError:
            iAlpha  = None
            iUmStop = iWtMax
        #print "iAlpha = ",iAlpha
        try:
            iUm = fields.index('UM')+iOffset
            iWtMax = iUm
        except ValueError:
            iUm = None
        #print "iAlpha=%s iUm=%s" %(iAlpha,iUm)
        #print "iAlpha=%s iWtMax=%s" %(iAlpha,iWtMax)
        #sys.stdout.flush()

        #print "iUM = ",iUM
        self.WtCG_groups = []
        
        i=iOffset
        while i<iWtMax:
            Gij = []
            wt = card.field(i)
            if wt is not None:
                ci = card.field(i+1)
                #print "wt=%s ci=%s" %(wt,ci)
                i+=2
                gij = 0
                while isinstance(gij,int) and i<iWtMax:  # does this get extra fields???
                    gij = card.field(i)
                    if gij is not None:
                        Gij.append(gij)
                    i+=1
                #print "Gij = ",Gij
                #print "gij_stop? = ",gij
                #if gij=='UM':
                    #print "breaking A..."
                    #raise Exception('error in RBE3...found UM in Weight Loop...bug...')
                    #break
                ###
                self.WtCG_groups.append([wt,ci,Gij])
            else:
                i+=1
            ###
        ###
        
        self.Gmi = []
        self.Cmi = []
        #print ""
        if iUm:
            #i+=1
            #raise NotImplementedError('parse the UM')
            i = iUm+1
            #print "i=%s iUmStop=%s" %(i,iUmStop)
            for j in range(i,iUmStop,2):  # does this get extra fields???
                gmi = card.field(j)
                if gmi is not None:
                    cmi = card.field(j+1)
                    #print "gmi=%s cmi=%s" %(gmi,cmi)
                    self.Gmi.append(gmi)
                    self.Cmi.append(cmi)
                ###
            ###
        ###
        if iAlpha:
            self.alpha = card.field(iAlpha+1)
            #print "self.alpha = ",self.alpha
        else:
            ## thermal expansion coefficient
            self.alpha = 0.0
        ###
        #print self
        #sys.exit()

    def rawFields(self):
        fields = ['RBE3',self.eid,None,self.refgrid,self.refc]
        for (wt,ci,Gij) in self.WtCG_groups:
            #print 'wt=%s ci=%s Gij=%s' %(wt,ci,Gij)
            fields+=[wt,ci]+Gij
        nSpaces = 8-(len(fields)-1)%8  # puts UM onto next line
        #print "nSpaces = ",nSpaces
        if nSpaces<8:
            fields += [None]*nSpaces
            
        if self.Gmi:
            fields += ['UM']
        if self.Gmi:
            #print "Gmi = ",self.Gmi
            #print "Cmi = ",self.Cmi
            for (gmi,cmi) in zip(self.Gmi,self.Cmi):
                fields+=[gmi,cmi]
            ###
        ###
        nSpaces = 8-(len(fields)-1)%8  # puts ALPHA onto next line
        if nSpaces<8:
            fields += [None]*nSpaces
        #print "nSpaces = ",nSpaces
        if self.alpha>0.: # handles the default value
            fields += ['ALPHA',self.alpha]
        return fields

    def reprFields(self):
        return self.rawFields()
