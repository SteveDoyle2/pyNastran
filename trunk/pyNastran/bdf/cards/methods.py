class Method(BaseCard):
    def __init__(self,card=None,data=None):
        pass

class EIGB(Method):
    """
    Defines data needed to perform buckling analysis
    """
    type = 'EIGB'
    def __init__(self,card=None,data=None):
        Method.__init__(self,card,data)
        if card:
            ## Set identification number. (Unique Integer > 0)
            self.sid = card.field(1)
            ## Method of eigenvalue extraction. (Character: “INV” for inverse power
            ## method or 'SINV' for enhanced inverse power method.)
            self.method = card.field(2,'INV')
            ## Eigenvalue range of interest. (Real, L1 < L2)
            self.L1   = card.field(3)
            self.L2   = card.field(4)
            assert self.L1 < self.L2,'L1=%s L2=%s; L1<L2 is requried' %(self.L1,self.L2)
            ## Estimate of number of roots in positive range not used for
            ## METHOD = 'SINV'. (Integer > 0)
            self.nep  = card.field(5,0)
            ## Desired number of positive and negative roots. (Integer>0; Default = 3*NEP)
            self.ndp  = card.field(6,3*self.nep)
            self.ndn  = card.field(7,3*self.nep)
            ## Method for normalizing eigenvectors. ('MAX' or 'POINT';Default='MAX')
            self.norm = card.field(9,'MAX')
            self.G    = card.field(10)
            self.C    = card.field(11)
        else:
            raise NotImplementedError('EIGB')
        ###

    def crossReference(self,model):
        pass

    def rawFields(self):
        fields = ['EIGB',self.sid,self.method,self.L1,self.L2,self.nep,self.ndn,None,
                         self.norm,self.G,self.C]
        return fields

    def reprFields(self):
        method = self.setBlankIfDefault(self.method,'INV')
        nep = self.setBlankIfDefault(self.nep,0)
        ndp = self.setBlankIfDefault(self.ndp,3*self.nep)
        ndn = self.setBlankIfDefault(self.ndn,3*self.nep)
        norm = self.setBlankIfDefault(self.norm,'MAX')
        fields = ['EIGB',self.sid,method,self.L1,self.L2,nep,ndn,None,
                         norm,self.G,self.C]
        return fields

class EIGC(Method): ## not done
    """
    Defines data needed to perform complex eigenvalue analysis
    """
    type = 'EIGC'
    def __init__(self,card=None,data=None):
        Method.__init__(self,card,data)
        if card:
            ## Set identification number. (Unique Integer > 0)
            self.sid = card.field(1)
            ## Method of complex eigenvalue extraction
            self.method = card.field(2)
            assert self.method in ['INV','HESS','CLAN']
            ## Method for normalizing eigenvectors
            self.norm   = card.field(3)

            ## Grid or scalar point identification number. Required only if NORM='POINT'. (Integer>0)
            self.G      = card.field(4)
            ## Component number. Required only if NORM='POINT' and G is a geometric grid point. (1<Integer<6)
            self.C      = card.field(5)
            ## Convergence criterion. (Real > 0.0. Default values are: 10-4 for
            ## METHOD = "INV", 10-15 for METHOD = "HESS", E is machine
            ## dependent for METHOD = "CLAN".)
            self.E      = card.field(6)
            self.ndo    = card.field(7)
            assert card.nFields()<8,'card = %s' %(card.fields(0))
        else:
            raise NotImplementedError('EIGC')
        ###

    def crossReference(self,model):
        pass

    def rawFields(self):
        fields = ['EIGC',self.sid,self.method,self.norm,self.G,self.C,self.E,self.ndo,None]
        raise Exception('EIGC not finished...')
        return fields

    def reprFields(self):
        fields = ['EIGC',self.sid,self.method,self.norm,self.G,self.C,self.E,self.ndo,None]
        raise Exception('EIGC not finished...')
        return fields

class EIGR(Method):
    """
    Defines data needed to perform real eigenvalue analysis
    """
    type = 'EIGR'
    def __init__(self,card=None,data=None):
        Method.__init__(self,card,data)
        if card:
            ## Set identification number. (Unique Integer > 0)
            self.sid = card.field(1)
            ## Method of eigenvalue extraction. (Character: “INV” for inverse power
            ## method or 'SINV' for enhanced inverse power method.)
            self.method = card.field(2,'LAN')
            ## Frequency range of interest
            self.f1   = card.field(3)
            self.f2   = card.field(4)
            ## Estimate of number of roots in range (Required for
            ## METHOD = 'INV'). Not used by 'SINV' method.
            self.ne  = card.field(5)
            ## Desired number of roots (default=600 for SINV 3*ne for INV)
            self.nd  = card.field(6)
            ## Method for normalizing eigenvectors. ('MAX' or 'POINT';Default='MAX')
            self.norm = card.field(9,'MASS')
            assert self.norm in ['POINT','MASS','MAX']
            ## Grid or scalar point identification number. Required only if NORM='POINT'. (Integer>0)
            self.G    = card.field(10)
            ## Component number. Required only if NORM='POINT' and G is a geometric grid point. (1<Integer<6)
            self.C    = card.field(11)
        else:
            raise NotImplementedError('EIGR')
        ###

    def crossReference(self,model):
        pass

    def rawFields(self):
        fields = ['EIGR',self.sid,self.method,self.f1,self.f2,self.ne,self.nd,None,None,
                         self.norm,self.G,self.C]
        return fields

    def reprFields(self):
        method = self.setBlankIfDefault(self.method,'LAN')
        norm = self.setBlankIfDefault(self.norm,'MASS')
        fields = ['EIGR',self.sid,method,self.f1,self.f2,self.ne,self.nd,None,None,
                         norm,self.G,self.C]
        return fields

class EIGP(Method):
    """
    Defines poles that are used in complex eigenvalue extraction by the Determinant method.
    """
    type = 'EIGP'
    def __init__(self,card=None,data=None):
        Method.__init__(self,card,data)
        if card:
            ## Set identification number. (Unique Integer > 0)
            self.sid = card.field(1)
            
            ## Coordinates of point in complex plane. (Real)
            self.alpha1 = card.field(2)
            ## Coordinates of point in complex plane. (Real)
            self.omega1 = card.field(3)
            ## Multiplicity of complex root at pole defined by point at ALPHAi and OMEGAi
            self.m1     = card.field(4)
            
            ## Coordinates of point in complex plane. (Real)
            self.alpha2 = card.field(5)
            ## Coordinates of point in complex plane. (Real)
            self.omega2 = card.field(6)
            ## Multiplicity of complex root at pole defined by point at ALPHAi and OMEGAi
            self.m2     = card.field(7)
        else:
            raise NotImplementedError('EIGP')
        ###

    def crossReference(self,model):
        pass

    def rawFields(self):
        fields = ['EIGP',self.alpha1,self.omega1,self.m1,self.alpha2,self.omega2,self.m2]
        return fields

    def reprFields(self):
        return self.rawFields()

class EIGRL(Method):
    """
    Defines data needed to perform real eigenvalue (vibration or buckling)
    analysis with the Lanczos method
    """
    type = 'EIGRL'
    def __init__(self,card=None,data=None):
        Method.__init__(self,card,data,sol=None)
        if card:
            ## Set identification number. (Unique Integer > 0)
            self.sid    = card.field(1)
            ## For vibration analysis: frequency range of interest. For buckling
            ## analysis: eigenvalue range of interest. See Remark 4. (Real or blank,
            ## -5 10e16 <= V1 < V2 <= 5.10e16
            self.v1     = card.field(2)
            self.v2     = card.field(3)
            ## Number of roots desired
            self.nd     = card.field(4)
            ## Diagnostic level. (0 < Integer < 4; Default = 0)
            self.msglvl = card.field(5,0)
            ## Number of vectors in block or set. Default is machine dependent
            self.maxset = card.field(6)
            ## Estimate of the first flexible mode natural frequency (Real or blank)
            self.shfscl = card.field(7)
            
            optionValues = card.fields(9)
            n = len(optionValues)
            nOptions = n//2
            assert n%2==0
            self.options = []
            self.values  = []
            for o in range(0,n,2):
                self.options.append(optionValues[o])
                self.values.append(optionValues[o+1])

            ## Method for normalizing eigenvectors
            if sol==105:
                self.norm = 'MAX'
            else:
                self.norm   = card.field(8)
            assert self.nrom in ['MASS','MAX']
            assert card.nFields()<9,'card = %s' %(card.fields(0))
        else:
            raise NotImplementedError('EIGRL')
        ###

    def crossReference(self,model):
        pass

    def rawFields(self):
        fields = ['EIGRL',self.sid,self.v1,self.v2,self.nd,self.msglvl,self.maxset,self.shfscl]
        for (option,value) in zip(self.options,self.values):
            fields += [option,value]
        return fields

    def reprFields(self):
        msglvl = self.setBlankIfDefault(self.msglvl,0)
        fields = ['EIGRL',self.sid,self.v1,self.v2,self.nd,msglvl,self.maxset,self.shfscl]
        for (option,value) in zip(self.options,self.values):
            fields += [option,value]
        return fields
