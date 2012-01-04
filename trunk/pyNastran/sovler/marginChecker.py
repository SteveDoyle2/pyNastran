import sys
from math import sqrt
from numpy import array
from numpy.linalg import norm

#pyNastran
from pyNastran.op2.op2 import OP2

class MarginChecker(object):
    def __init__(self,filenames=['fem.bdf'],subcases=[1],IDs=[None]):
        """
        Performs load case combination for:
        @param filenames  list of filenames that the subcase result will come from
        @param subcases   list of subcases to grab from from each filename
        Assumptions:
            linear static analysis
            nodes numbers are consistent across different OP2s/subcases
            only does static loading & maxDeflection
            solid elements stress ONLY
        
        User Info:
        UnitCase  Tension      Compression   Bending
        Filename  tension.op2  compBend.op2  compBend.op2
        iSubcase  1            1             2
        
        unitLoad  100.         1.            1.
        ReqLoad   100.         100.          100.
        Ratio     1.           100.          100.

        Code Inputs:
        filenames =  ['tension.op2','compBend.op2','compBend.op2']
        subcases  =  [1,            1,             2,            ]
        vmFactors = [[1.,           100.,          100.          ]] # only for VM
        IDs       =  ['tension',    'comp',       'bend'         ]
        
        @note vmFactors and IDs may have multiple levels...
            vmFactors = [[1.,100.,100.],
                         [2., 50., 75.],]
            caseNames = ['Case1','Case2']  # only for VM
        """
        self.filenames = filenames
        self.subcases  = subcases
        self.IDs = IDs
        assert len(self.filenames)==len(self.subcases)==len(self.IDs)

        print "filenames = ",self.filenames
        print "subcases  = ",self.subcases
        print "IDs       = ",self.IDs
        
        #self.loadFactors = [[1.,2.,],  # fatigue
        #                    [4.,5.,]]
        self.cases = {}

        self.solidStressResults = {}
        self.displacementResults = {}
        #self.solidStrainResults = {}

        #self.plateStressResults = {}
        #self.plateStrainResults = {}

        #self.compositePlateStressResults = {}
        #self.compositePlateStrainResults = {}
        
        for fname in self.filenames:
            self.cases[fname] = []
            #self.solidStressResults[fname] = None
            #self.plateStressResults[fname] = 
        
        #print "self.subcases2 = ",self.subcases
        for (fname,subcaseID,ID) in zip(self.filenames,self.subcases,self.IDs):
            self.cases[fname].append(subcaseID)
        
        for key in self.cases:
            print "case[%s] = %s" %(key,self.cases[key])

    def readFiles(self):
        i = 0
        for fname,subcaseList in sorted(self.cases.items()):
            #print "case[%s] = %s" %(key,self.cases[key])
            subcaseList = list(set(subcaseList))
            print "subcaseList[%s] = %s" %(fname,subcaseList)

            op2 = OP2(fname,debug=False)
            op2.setSubcases(subcaseList)
            op2.readOP2()
            
            for subcaseID in subcaseList:
                print "subcaseID = ",subcaseID
                print "i = ",i
                self.displacementResults[i] = op2.displacements[subcaseID]
                self.solidStressResults[i] = op2.solidStress[subcaseID]
                #self.solidStrainResults[i] = op2.solidStrain

                #self.plateStressResults[i] = op2.plateStress
                #self.plateStrainResults[i] = op2.plateStrain

                #self.compositePlateStressResults[i] = op2.compositePlateStress
                #self.compositePlateStrainResults[i] = op2.compositePlateStrain

                #self.cleanSolidStress(op2.solidStress[subcaseID])
                i +=1
                ###
            del op2 # makes sure that unneeded data is not stored to save memory
        ###

    def cleanSolidStress(self,stress):
        del stress.oxx
        del stress.oyy
        del stress.ozz
        del stress.txy
        del stress.tyz
        del stress.txz
        del stress.ovmShear

    def checkDeflections(self,maxDeflection):
        """
        @param maxDeflection [x,y,z,combined] abs values are used
        """
        self.maxDeflection = maxDeflection

        deflect = self.displacementResults[0] # gets the tension case
        nodesList = deflect.translations.keys()
        
        #deflections = {}
        deflectionDict = {}
        for icase,vmFactor in enumerate(self.vmFactors):
            deflectionDict[icase] = {}

            for nid in nodesList:
                deflection = array([0.,0.,0.])
                for jfact,factor in enumerate(vmFactor):
                    deflect = self.displacementResults[jfact] # gets the tension case
                    trans = deflect.translations[nid]*factor
                    deflection += trans
                #print "deflection[%s][%s][%s] = %s" %(icase,jfact,nid,str(deflection))
                deflectionDict[icase][nid] = deflection
            ###
        ###
        deflectionMargin = {}
        for nid in nodesList:
            minDeflectMargins = []
            for icase,vmFactor in enumerate(self.vmFactors):
                normT = norm(deflectionDict[icase][nid])
                
                if nid in [487,497]:
                    print "nid=%s normT=%s" %(nid,normT)
                defMargin = (self.maxDeflection-normT)/self.maxDeflection##-1
                minDeflectMargins.append(defMargin)
            ###
            minMargin = min(minDeflectMargins)
            case = minDeflectMargins.index(minMargin)

            deflectionMargin[nid] = minMargin
            print "case=%-6s minMargin[%s] = %g" %(self.caseNames[case],nid,minMargin)

    def checkVonMises(self,vmFactors=[[1.]],caseNames=['Case1'],Fty=100.):
        """
        currently only handles von mises stress for solid elements...
        @param self the object pointer
        @param vmFactors @see self.__init__
        @param caseNames @see self.__init__
        @param Fty  max allowable yeild stress (same units as FEM stress!)

        \f[ \sigma_v^2 = \tfrac{1}{2}[(\sigma_{11} - \sigma_{22})^2 + (\sigma_{22} - \sigma_{33})^2 + (\sigma_{11} - \sigma_{33})^2 + 6(\sigma_{23}^2 + \sigma_{31}^2 + \sigma_{12}^2)]  \f]
        \f[ \sigma_v   = \sqrt{\frac{(\sigma_1 - \sigma_2)^2 + (\sigma_2 - \sigma_3)^2 + (\sigma_1 - \sigma_3)^2 } {2}} \f]
        """
        self.vmFactors = vmFactors
        self.caseNames = caseNames
        self.Fty = Fty # ksi
        assert len(self.caseNames)==len(self.vmFactors)

        print "vmFactors = ",self.vmFactors
        print "caseNames = ",self.caseNames
        stressP = []
        for icase,vmFactor in enumerate(self.vmFactors):
            StressPrincipal = [] # stores cases
            for jfact,factor in enumerate(vmFactor):
                stressPrincipal = {} # stores eid
                stress = self.solidStressResults[jfact] # gets the tension case
                eidList = stress.o1.keys()
                for eid in eidList:
                    stress_principal = {}
                    nodeList = stress.o1[eid].keys()
                    for nid in sorted(nodeList):
                        #stressNode = {}
                        stress_principal[nid] = array([0.,0.,0.])
                    ###
                    stressPrincipal[eid] = stress_principal # gets all the nodes
                ###
                StressPrincipal = stressPrincipal # stores all the eids
                break
            stressP.append(StressPrincipal)

        for icase,vmFactor in enumerate(self.vmFactors):
            for jfact,factor in enumerate(vmFactor):
                stress = self.solidStressResults[jfact] # gets the tension case
                #print str(stress)[0:600]
                #print "eidList = ",eidList
                for eid in eidList:
                    nodeList = stress.o1[eid].keys()
                    for nid in sorted(nodeList):
                        o1 = stress.o1[eid][nid]
                        o2 = stress.o2[eid][nid]
                        o3 = stress.o3[eid][nid]
                        #print "s[%s][%s][%s][%s]" %(icase,jfact,eid,nid)
                        #print "vmFactor = ",vmFactor
                        #print "s[%s][%s][%s] = %s" %(icase,eid,nid,stressP[icase][eid][nid])
                        o = array([o1,o2,o3])*factor
                        stressP[icase][eid][nid] += o

                        #ov = sqrt((o1-o2)**2+(o2-o3)**2+(o1-o3)**2)
                        #ov = sqrt((o[0]-o[1])**2+(o[1]-o[2])**2+(o[0]-o[2])**2)
                    ###
                ###

        del o
        for icase,vmFactor in enumerate(self.vmFactors):
            print "icase = ",icase
            #print "vmFactor = ",vmFactor
            for eid in eidList:
                #print "eid = ",eid
                eidResults = []
                nodeList = stress.o1[eid].keys()
                #nodeList = stress[icase][eid].keys()
                for nid in sorted(nodeList):
                    #print "nid = ",nid
                    #print "o = ",o
                    #o = stressP[icase][eid]#[nid]
                    o = stressP[icase][eid]['C']
                    #print o.keys()
                    ov = sqrt((o[0]-o[1])**2+(o[1]-o[2])**2+(o[0]-o[2])**2)
                    eidResults.append(ov)
                stressP[icase][eid] = min(eidResults)
            ###
        ###
        #print stressP
        self.stressP = stressP
        
        Ftu = 150. # ksi
        minMargins = {}
        for eid in eidList:
            margins = []
            for icase,vmFactor in enumerate(self.vmFactors):
                margin = stressP[icase][eid]/Ftu-1
                margins.append(margin)
            minMargin = min(margins)
            minMargins[eid] = minMargin
            case = margins.index(minMargin)
            print "case=%-6s minMargin[%s] = %g" %(self.caseNames[case],eid,minMargin)

    def rainflow(self):
        pass


def main():
    filenames = ['fem.op2','fem.op2']
    subcases  = [1,2]
    Fty = 150. # ksi
    IDs = ['tension','compression']

    vmFactors =   [[1.,0.,],  # tension
                   [0.,1.,],  # bending
                   [0.,-1.,], # minus bending
                   [-1.,0.],  # compression
                   [-1.,1.],  # compression+bend
                   [1.,1.,],  # tension + bending
                   [-1.,-1.], # compression+mbend
                   [1.,-1.,],]# tension + mbending


    caseNames = ['tens','bend','mBend','comp','cBend','tBend','cmBend','tmBend']
    a = MarginChecker(filenames,subcases,IDs)
    a.readFiles()
    a.checkVonMises(vmFactors,caseNames,Fty)
    #a.checkDeflections(0.018)

if __name__=='__main__':
    main()
