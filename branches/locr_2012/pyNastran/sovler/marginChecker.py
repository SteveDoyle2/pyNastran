import sys
from itertools import izip
from math import sqrt, log10
from numpy import array
from numpy.linalg import norm

from pyNastran.op2.op2 import OP2


class MarginChecker(object):
    def __init__(self, filenames=['fem.bdf'], subcases=[1], IDs=[None]):
        """
        Performs load case combination for:
        @param filenames  list of filenames that the subcase result will come from
        @param subcases   list of subcases to grab from from each filename
        
        Assumptions:
            * linear static analysis
            * nodes numbers are consistent across different OP2s/subcases
            * only does static loading & maxDeflection
            * solid elements stress ONLY

        @code
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
        # only for VonMises
        vmFactors = [[1.,           100.,          100.          ]]
        
        IDs       =  ['tension',    'comp',       'bend'         ]
        @endcode

        @note vmFactors and IDs may have multiple levels...
            
        @code
             vmFactors = [[1.,100.,100.],
                          [2., 50., 75.],]
             caseNames = ['Case1','Case2']  # only for VM``
        @endcode
        """
        self.filenames = filenames
        self.subcases = subcases
        self.IDs = IDs
        assert len(self.filenames) == len(self.subcases) == len(self.IDs)

        print "filenames = ", self.filenames
        print "subcases  = ", self.subcases
        print "IDs       = ", self.IDs

        #self.loadFactors = [[1.,2.,],  # fatigue
        #                    [4.,5.,]]
        self.cases = {}

        self.displacementResults = {}
        self.solidStressResults = {}
        self.plateStressResults = {}

        #self.solidStrainResults = {}
        #self.plateStrainResults = {}

        #self.compositePlateStressResults = {}
        #self.compositePlateStrainResults = {}

        for fname in self.filenames:
            self.cases[fname] = []
            #self.solidStressResults[fname] = None
            #self.plateStressResults[fname] =

        #print "self.subcases2 = ",self.subcases
        for (fname, subcaseID, ID) in izip(self.filenames, self.subcases, self.IDs):
            self.cases[fname].append(subcaseID)

        for key in self.cases:
            print "case[%s] = %s" % (key, self.cases[key])

    def readFiles(self):
        i = 0
        for fname, subcaseList in sorted(self.cases.iteritems()):
            #print "case[%s] = %s" %(key,self.cases[key])
            subcaseList = list(set(subcaseList))
            print "subcaseList[%s] = %s" % (fname, subcaseList)

            op2 = OP2(fname, debug=False)
            op2.setSubcases(subcaseList)
            op2.readOP2()

            for subcaseID in subcaseList:
                print "subcaseID = ", subcaseID
                print "i = ", i
                self.displacementResults[i] = op2.displacements[subcaseID]

                #self.solidStrainResults[i] = op2.solidStrain
                #self.plateStrainResults[i] = op2.plateStrain

                #self.compositePlateStressResults[i] = op2.compositePlateStress
                #self.compositePlateStrainResults[i] = op2.compositePlateStrain

                self.cleanStress(i, op2, subcaseID)
                i += 1

            del op2  # makes sure that unneeded data is not stored to save memory

    def cleanStress(self, i, op2, subcaseID):
        self.cleanSolidStress(i, op2, subcaseID)
        self.cleanPlateStress(i, op2, subcaseID)

    def cleanPlateStress(self, i, op2, subcaseID):
        if subcaseID in op2.plateStress:
            self.plateStressResults[i] = op2.plateStress[subcaseID]
            stress = op2.plateStress[subcaseID]
            del stress.fiberCurvature
            #del stress.oxx
            #del stress.oyy
            #del stress.txy
            del stress.majorP
            del stress.minorP
            del stress.angle
            del stress.ovmShear

    def cleanSolidStress(self, i, op2, subcaseID):
        if subcaseID in op2.solidStress:
            self.solidStressResults[i] = op2.solidStress[subcaseID]
            stress = op2.solidStress[subcaseID]
            del stress.o1
            del stress.o2
            del stress.o3
            #del stress.oxx
            #del stress.oyy
            #del stress.ozz
            #del stress.txy
            #del stress.tyz
            #del stress.txz
            del stress.ovmShear

    def checkDeflections(self, maxDeflection):
        """
        @param maxDeflection [x,y,z,combined] abs values are used
        """
        self.maxDeflection = maxDeflection

        deflect = self.displacementResults[0]  # gets the tension case
        nodesList = deflect.translations.keys()

        deflectionDict = {}
        for icase, vmFactor in enumerate(self.vmFactors):
            deflectionDict[icase] = {}

            for nid in nodesList:
                deflection = array([0., 0., 0.])
                for jfact, factor in enumerate(vmFactor):
                    # gets the tension case
                    deflect = self.displacementResults[jfact]
                    trans = deflect.translations[nid] * factor
                    deflection += trans
                #print "deflection[%s][%s][%s] = %s" %(icase,jfact,nid,str(deflection))
                deflectionDict[icase][nid] = deflection

        deflectionMargin = {}
        for nid in nodesList:
            minDeflectMargins = []
            for icase, vmFactor in enumerate(self.vmFactors):
                normT = norm(deflectionDict[icase][nid])

                if nid in [487, 497]:
                    print "nid=%s normT=%s" % (nid, normT)
                defMargin = (self.maxDeflection - normT) / self.maxDeflection  # -1
                minDeflectMargins.append(defMargin)

            minMargin = min(minDeflectMargins)
            case = minDeflectMargins.index(minMargin)

            deflectionMargin[nid] = minMargin
            print "case=%-6s minMargin[%s] = %g" % (
                self.caseNames[case], nid, minMargin)

    def checkVonMises(self, vmFactors=[[1.]], caseNames=['Case1'], Fty=100.):
        r"""
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
        self.Fty = Fty  # ksi
        assert len(self.caseNames) == len(self.vmFactors)

        print "vmFactors = ", self.vmFactors
        print "caseNames = ", self.caseNames
        (stressP, eidList) = self.processSolidStress()
        (stressP, eidList) = self.processPlateStress()
        #print stressP
        self.stressP = stressP

        Fty = self.Fty  # ksi
        minMargins = {}
        for eid in eidList:
            margins = []
            for icase, vmFactor in enumerate(self.vmFactors):
                margin = stressP[icase][eid] / Fty - 1
                margins.append(margin)
            minMargin = min(margins)
            minMargins[eid] = minMargin
            case = margins.index(minMargin)
            print "case=%-6s minMargin[%s] = %g" % (
                self.caseNames[case], eid, minMargin)

    def processPlateStress(self):
        r"""
        \f[ \sigma_v = \sqrt{\sigma_1^2- \sigma_1\sigma_2+ \sigma_2^2 + 3\sigma_{12}^2} \f]
        ovm^2 = o1^2 - o1*o2 + o2^2 + 3*o12^2
        """
        stressP = []
        eidList = []
        for icase, vmFactor in enumerate(self.vmFactors):
            StressPrincipal = []  # stores cases
            for jfact, factor in enumerate(vmFactor):
                stressPrincipal = {}  # stores eid
                if jfact not in self.plateStressResults:
                    continue

                # gets the tension case
                stress = self.plateStressResults[jfact]
                eidList = stress.o1.keys()
                for eid in eidList:
                    stress_principal = {}
                    nodeList = stress.o1[eid].keys()
                    for nid in sorted(nodeList):
                        stress_principal[nid] = array([0., 0., 0.])

                    stressPrincipal[eid] = stress_principal  # gets all the nodes
                StressPrincipal = stressPrincipal  # stores all the eids
                break
            stressP.append(StressPrincipal)

        if eidList == []:
            return None, []
        for icase, vmFactor in enumerate(self.vmFactors):
            for jfact, factor in enumerate(vmFactor):
                # gets the tension case
                stress = self.plateStressResults[jfact]
                #print str(stress)[0:600]
                #print "eidList = ",eidList
                for eid in eidList:
                    nodeList = stress.o1[eid].keys()
                    for nid in sorted(nodeList):
                        oxx = stress.oxx[eid][nid]
                        oyy = stress.oyy[eid][nid]
                        txy = stress.txy[eid][nid]
                        #print "s[%s][%s][%s][%s]" %(icase,jfact,eid,nid)
                        #print "vmFactor = ",vmFactor
                        #print "s[%s][%s][%s] = %s" %(icase,eid,nid,stressP[icase][eid][nid])
                        o = array([oxx, oyy, txy]) * factor
                        stressP[icase][eid][nid] += o  # @todo is this legal...check with oxx, oyy, ozz, txy... ???

                        #ov = sqrt((o1-o2)**2+(o2-o3)**2+(o1-o3)**2)
                        #ov = sqrt((o[0]-o[1])**2+(o[1]-o[2])**2+(o[0]-o[2])**2)

        for icase, vmFactor in enumerate(self.vmFactors):
            print "icase = ", icase
            #print "vmFactor = ",vmFactor
            for eid in eidList:
                #print "eid = ",eid
                eidResults = []
                nodeList = stress.oxx[eid].keys()
                #nodeList = stress[icase][eid].keys()
                for nid in sorted(nodeList):
                    #print "nid = ",nid
                    #print "o = ",o
                    #o = stressP[icase][eid]#[nid]
                    (oxx, oyy, txy) = stressP[icase][eid][nid]
                    #print o.keys()
                    #ovm^2 = o1^2 - o1*o2 + o2^2 + 3*o12^2
                    ovm = sqrt(oxx ** 2 + oyy ** 2 - oxx *
                               oyy + 3 * txy ** 2)  # 2d stress
                    eidResults.append(ovm)
                stressP[icase][eid] = min(eidResults)

        return (stressP, eidList)

    def processSolidStress(self):  # stressP[icase][eid] = min(eidResults)
        stressP = []
        for icase, vmFactor in enumerate(self.vmFactors):
            StressPrincipal = []  # stores cases
            for jfact, factor in enumerate(vmFactor):
                stressPrincipal = {}  # stores eid
                # gets the tension case
                stress = self.solidStressResults[jfact]

                eidList = stress.oxx.keys()
                for eid in eidList:
                    stress_principal = {}
                    nodeList = stress.oxx[eid].keys()
                    for nid in sorted(nodeList):
                        stress_principal[nid] = array([0., 0., 0., 0., 0., 0.])
                    stressPrincipal[eid] = stress_principal  # gets all the nodes

                StressPrincipal = stressPrincipal  # stores all the eids
                break
            stressP.append(StressPrincipal)

        for icase, vmFactor in enumerate(self.vmFactors):
            for jfact, factor in enumerate(vmFactor):
                # gets the tension case
                stress = self.solidStressResults[jfact]
                
                #print str(stress)[0:600]
                #print "eidList = ",eidList
                for eid in eidList:
                    nodeList = stress.oxx[eid].keys()
                    for nid in sorted(nodeList):
                        oxx = stress.oxx[eid][nid]
                        oyy = stress.oyy[eid][nid]
                        ozz = stress.ozz[eid][nid]
                        txy = stress.txy[eid][nid]
                        tyz = stress.tyz[eid][nid]
                        txz = stress.txz[eid][nid]

                        #print "s[%s][%s][%s][%s]" %(icase,jfact,eid,nid)
                        #print "vmFactor = ",vmFactor
                        #print "s[%s][%s][%s] = %s" %(icase,eid,nid,stressP[icase][eid][nid])
                        o = array([oxx, oyy, ozz, txy, tyz, txz]) * factor
                        stressP[icase][eid][nid] += o  # @todo is this legal...check with oxx, oyy, ozz, txy... ???

                        #ov = sqrt((o1-o2)**2+(o2-o3)**2+(o1-o3)**2)
                        #ov = sqrt((o[0]-o[1])**2+(o[1]-o[2])**2+(o[0]-o[2])**2)

        del o
        for icase, vmFactor in enumerate(self.vmFactors):
            print "icase = ", icase
            #print "vmFactor = ",vmFactor
            for eid in eidList:
                #print "eid = ",eid
                eidResults = []
                nodeList = stress.oxx[eid].keys()
                #nodeList = stress[icase][eid].keys()
                for nid in sorted(nodeList):
                    #print "nid = ",nid
                    #print "o = ",o
                    #o = stressP[icase][eid]#[nid]
                    (oxx, oyy, ozz, txy, tyz, txz) = stressP[icase][eid][nid]
                    #print o.keys()
                    ovm = sqrt((oxx - oyy) ** 2 + (oyy - ozz) ** 2 + (oxx - ozz) ** 2 + 6 * (txy ** 2 + tyz ** 2 + txz ** 2))  # 3d stress
                    print "ovm = %s" % (ovm)
                    eidResults.append(ovm)
                stressP[icase][eid] = min(eidResults)

        return (stressP, eidList)

    def rainflow(self):
        """required for fatigue"""
        pass

    def damageCount(self):
        """required for fatigue"""
        # 1.  get principal stresses at different loading combinations
        # 2.  rainflow count
        # 3.  find damage caused by each load case
        # 4.  sum damage (Miner's Rule)
        # 5.  find equivalent stress
        # 6.  get damage margin per ratio of life used
        pass

    def SNcurve(stress):
        """
        @code
        stress amplitude = fatigue strength coefficient* (2 N(f))^b
        N(f) is the cycles to failure
        2N(f) is the number of load reversals to failure
        b is the fatigue strength exponent
        For an AISI Type 1015 steel
        b = -0.11
        fatigue strength coefficient where 2N(f) =1,120 Ksi


        Log Sf = Log a + b Log N = Log (1.62 Sut) + Log N-0.0851
        a = (0.9 Sut)2 / Se = 1.62 Sut
        b = -(Log (0.9 Sut / Se)) / 3 (= -0.0851 when Se = .5 Sut)
        Se = 0.5 Ftu  = Endurance Limit = Stress corresponding to 'infinite' life of 1,000,000 or more cycles.
        Sf = Stress corresponding to a fatigue life, N, of 1000 to 1,000,000 cycles inclusive.
        @endcode
        """
        Sut = 120.   # ksi
        Se = 0.5 * Sut  # ksi
        a = (0.9 * Sut) * 2 / Se
        b = -(log10(0.9 * Sut / Se)) / 3.
        #logSf = log10(a) + b*log10(N)
        #logN = (logSf-log10(a) )/b
        Ncycles = 10. ** logN

        return Ncycles


def main():
    #filenames = ['fem.op2','fem.op2']
    subcases = [1, 2]
    Fty = 150.  # ksi
    IDs = ['tension', 'compression']

    vmFactors = [  # [1.,0.,],  # tension
                   #[0.,1.,],  # bending
                   #[0.,-1.,], # minus bending
                   #[-1.,0.],  # compression
                   #[-1.,1.],  # compression+bend
                   #[1.,1.,],  # tension + bending
                   #[-1.,-1.], # compression+mbend
                   [1., -1., ], ]  # tension + mbending

    #caseNames = ['tens','bend','mBend','comp','cBend','tBend','cmBend','tmBend']
    caseNames = ['tmBend']
    a = MarginChecker(filenames, subcases, IDs)
    a.readFiles()
    a.checkVonMises(vmFactors, caseNames, Fty)
    #a.checkDeflections(0.018)

if __name__ == '__main__':
    main()
