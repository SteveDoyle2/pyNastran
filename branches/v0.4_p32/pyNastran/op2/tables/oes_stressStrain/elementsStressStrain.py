## GNU Lesser General Public License
## 
## Program pyNastran - a python interface to NASTRAN files
## Copyright (C) 2011-2012  Steven Doyle, Al Danial
## 
## Authors and copyright holders of pyNastran
## Steven Doyle <mesheb82@gmail.com>
## Al Danial    <al.danial@gmail.com>
## 
## This file is part of pyNastran.
## 
## pyNastran is free software: you can redistribute it and/or modify
## it under the terms of the GNU Lesser General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
## 
## pyNastran is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU Lesser General Public License
## along with pyNastran.  If not, see <http://www.gnu.org/licenses/>.
## 
import sys
from struct import unpack
sys.setrecursionlimit(2000)

from pyNastran.op2.op2Errors import *

#91  -> PENTANL
#2   -> BEAM
#33  -> TUBE
#92  -> CONRODNL
class ElementsStressStrain(object):

    def skipOES_Element(self):
        self.log.debug('skipping approach/table/format/sortCode=%s on %s table' %(self.atfsCode,self.tableName))
        self.skipOES_Element2()

    def skipOES_Element2(self): # works???
        self.data = b''
        self.handleResultsBuffer(self.skipOES_Element2)

    def OES_Thermal(self,debug=False): # works
        if self.makeOp2Debug:
            self.op2Debug.write('---OES_Thermal---\n')
        deviceCode = self.deviceCode
        #assert self.numWide==5,'invalid numWide...numWide=%s' %(self.numWide)
        
        while len(self.data)>=32:
            #print self.printSection(40)
            eData     = self.data[0:32]
            self.data = self.data[32: ]
            #print "len(data) = ",len(eData)

            out = unpack('iiifffff',eData)
            (eid,sideID,hbdyID,cnvCoeff,fApplied,fConv,fRad,fTotal) = out
            eid = (eid - deviceCode) // 10
            #print "eid=%s sideID=%s hbdyID=%s coeff=%s fApplied=%s fConv=%s fRad=%s fTotal=%s" %(eid,sideID,hbdyID,cnvCoeff,fApplied,fConv,fRad,fTotal)
            if self.makeOp2Debug:
                self.op2Debug.write('%s\n' %(str(out)))
            #self.obj.addNewEid(eid,axial,axialMS,torsion,torsionMS)

            #print "eid=%i axial=%i torsion=%i" %(eid,axial,torsion)
            #print "len(data) = ",len(self.data)
        ###
        self.handleResultsBuffer(self.OES_Thermal)
        #print self.rodStress[self.iSubcase]
        ###

    def basicElement(self):
        """
        genericStressReader - works on CROD_1, CELAS2_12
        stress & strain
        formatCode=1 sortCode=0 (eid,axial,axialMS,torsion,torsionMS)
        formatCode=1 sortCode=1 (eid,axial,axial,torsion,torsion)
        """
        deviceCode = self.deviceCode
        (nTotal,dataFormat) = self.obj.getLength()
        #print "nTotal=%s dataFormat=%s len(data)=%s" %(nTotal,dataFormat,len(self.data))
        
        n = 0
        nEntries = len(self.data)//nTotal
        for i in range(nEntries):
            eData = self.data[n:n+nTotal]
            out = unpack(dataFormat,eData)
            #print "out = ",out
            self.obj.addNewEid(out)
            n+=nTotal
        ###
        self.data = self.data[n: ]
        self.handleResultsBuffer(self.basicElement)

    def CBEAM_2(self): # not tested
        if self.makeOp2Debug:
            self.op2Debug.write('---BEAM_2---\n')
        deviceCode = self.deviceCode
        nNodes = 10 # 11-1

        nTotal       = self.obj.getLengthTotal()
        (n1,format1) = self.obj.getLength1()
        (n2,format2) = self.obj.getLength2()

        while len(self.data)>=nTotal:
            eData     = self.data[0:n1]
            self.data = self.data[n1: ]
            #print "len(data) = ",len(eData)

            out = unpack(format1, eData)
            #print "outA = ",out
            eid = self.obj.addNewEid(out)
            
            for iNode in range(nNodes):
                eData     = self.data[0:n2]
                self.data = self.data[n2: ]
                out = unpack(format2, eData)
                #print "outB = ",out
                self.obj.add(eid,out)

            #print "eid=%i axial=%i torsion=%i" %(eid,axial,torsion)
            #print "len(data) = ",len(self.data)
        ###
        self.handleResultsBuffer(self.CBEAM_2)
        #raise Exception('add CBEAM-2...')

    def CQUAD4_33(self): # works
        """
        GRID-ID  DISTANCE,NORMAL-X,NORMAL-Y,SHEAR-XY,ANGLE,MAJOR MINOR,VONMISES
        """
        if self.makeOp2Debug:
            self.op2Debug.write('---CQUAD4_33---\n')
        deviceCode = self.deviceCode
        nNodes = 0 # centroid + 4 corner points
        #self.printSection(20)
        #term = data[0:4] CEN/
        #data = data[4:]
        #print "*****"
        #self.printBlock(self.data)
        #print "self.numWide = ",self.numWide
        #if 0:
        
        assert self.numWide==17,'invalid numWide...numWide=%s' %(self.numWide)
        while len(self.data)>=68: # 2+17*5 = 87 -> 87*4 = 348
            #print self.printBlock(self.data[0:100])
            #(eid,) = unpack("i",self.data[0:4])
            #print "abcd=",a,b,c,d
            #self.data = self.data[8:]  # 2
            eData     = self.data[0:4*17]
            self.data = self.data[4*17: ]
            out = unpack('iffffffffffffffff',eData)  # 17
            if self.makeOp2Debug:
                self.op2Debug.write('%s\n' %(str(out)))
            (eid,fd1,sx1,sy1,txy1,angle1,major1,minor1,maxShear1,
                 fd2,sx2,sy2,txy2,angle2,major2,minor2,maxShear2) = out
            eid = (eid - deviceCode) // 10

            #print "eid=%i grid=%s fd1=%-3.1f sx1=%i sy1=%i txy1=%i angle1=%i major1=%i minor1=%i vm1=%i" %(eid,'C',fd1,sx1,sy1,txy1,angle1,major1,minor1,maxShear1)
            #print   "             fd2=%-3.1f sx2=%i sy2=%i txy2=%i angle2=%i major2=%i minor2=%i vm2=%i\n"       %(fd2,sx2,sy2,txy2,angle2,major2,minor2,maxShear2)
            #print "nNodes = ",nNodes
            self.obj.addNewEid('CQUAD4',eid,'C',fd1,sx1,sy1,txy1,angle1,major1,minor1,maxShear1)
            self.obj.add(               eid,'C',fd2,sx2,sy2,txy2,angle2,major2,minor2,maxShear2)
            
            for nodeID in range(nNodes):   #nodes pts
                eData     = self.data[0:4*17]
                self.data = self.data[4*17: ]
                out = unpack('iffffffffffffffff',eData[0:68])
                if self.makeOp2Debug:
                    self.op2Debug.write('%s\n' %(str(out)))
                (grid,fd1,sx1,sy1,txy1,angle1,major1,minor1,vm1,
                      fd2,sx2,sy2,txy2,angle2,major2,minor2,vm2,) = out

                #print "eid=%i grid=%i fd1=%i sx1=%i sy1=%i txy1=%i angle1=%i major1=%i minor1=%i vm1=%i" %(eid,grid,fd1,sx1,sy1,txy1,angle1,major1,minor1,vm1)
                #print "               fd2=%i sx2=%i sy2=%i txy2=%i angle2=%i major2=%i minor2=%i vm2=%i\n"          %(fd2,sx2,sy2,txy2,angle2,major2,minor2,vm2)
                #print "len(data) = ",len(self.data)
                #self.printBlock(self.data)
                self.obj.addNewNode(eid,grid,fd1,sx1,sy1,txy1,angle1,major1,minor1,vm1)
                self.obj.add(       eid,grid,fd2,sx2,sy2,txy2,angle2,major2,minor2,vm2)
            ###
            #print '--------------------'
            #print "len(data) = ",len(self.data)
            #print "tell = ",self.op2.tell()
            
            #self.printSection(100)
            #self.dn += 348
        ###
        self.handleResultsBuffer(self.CQUAD4_33)

    def CBAR_34(self): # ???
        if self.makeOp2Debug:
            self.op2Debug.write('---CBAR_34---\n')
        deviceCode = self.deviceCode
        #print "len(data) = ",len(self.data)
        assert self.numWide==16,'invalid numWide...numWide=%s' %(self.numWide)
        while len(self.data)>=64:
            #self.printBlock(self.data)
            eData     = self.data[0:4*16]
            self.data = self.data[4*16: ]
            #print "len(data) = ",len(eData)

            (eid,s1a,s2a,s3a,s4a,axial,smaxa,smina,MSt,
                 s1b,s2b,s3b,s4b,      smaxb,sminb,MSc)= unpack('ifffffffffffffff',eData)
            eid = (eid - deviceCode) // 10
            self.obj.addNewEid('CBAR',eid,s1a,s2a,s3a,s4a,axial,smaxa,smina,MSt,
                                          s1b,s2b,s3b,s4b,      smaxb,sminb,MSc)

            #print "eid=%i s1=%i s2=%i s3=%i s4=%i axial=%-5i smax=%i smin=%i MSt=%i MSc=%i" %(eid,s1a,s2a,s3a,s4a,axial,smaxa,smina,MSt,MSc)
            #print "         s1=%i s2=%i s3=%i s4=%i          smax=%i smin=%i" %(s1b,s2b,s3b,s4b,smaxb,sminb)
            #print "len(data) = ",len(self.data)
        ###
        self.handleResultsBuffer(self.CBAR_34)

    def CSOLID_67(self):  # works
        """
        stress is extracted at the centroid
        CTETRA_39
        CPENTA_67
        CHEXA_68
        """
        if self.makeOp2Debug:
            self.op2Debug.write('---CSOLID_67---\n')
        #print "starting solid element..."
        deviceCode = self.deviceCode
        #nNodes = 5 # 1 centroid + 4 corner points
        #self.printSection(20)
        #term      = self.data[0:4] CEN/
        #self.data = self.data[4:]
        #print "*****"
        ElementType = self.ElementType(self.elementType)
        nNodes=1  # this is a minimum, it will be reset later
        nNodesExpected = 1
        #assert self.numWide in [109,151,193],'invalid numWide...numWide=%s' %(self.numWide)
        while len(self.data)>= 16+84*nNodesExpected:
            eData     = self.data[0:16]
            self.data = self.data[16:]
            #self.printBlock(eData)

            out = unpack("iissssi",eData)
            (eid,cid,a,b,c,d,nNodes) = out
            if self.makeOp2Debug:
                self.op2Debug.write('%s\n' %(str(out)))
            #print "abcd = |%s|" %(a+b+c+d)
            #print "eid=%s cid=%s nNodes=%s nNodesExpected=%s" %(eid,cid,nNodes,nNodesExpected)
            
            assert nNodes < 21,self.printBlock(eData)
            eid = (eid - deviceCode) // 10

            if   ElementType=='TETRA':   nNodesExpected = 5
            elif ElementType=='PENTA':   nNodesExpected = 7
            elif ElementType=='HEXA':    nNodesExpected = 9
            else:
                raise AddNewElementError('not supported....EType=%s eType=%s nNodes=%s numWide=%s' %(ElementType,self.elementType,nNodes,self.numWide))

            #print "len(data) = ",len(self.data)
            for nodeID in range(nNodesExpected):   #nodes pts, +1 for centroid (???)
                #print "len(data)A = ",len(self.data)
                eData     = self.data[0:4*21]  # for the stresses
                self.data = self.data[4*21: ]
                #print "len(data)B = ",len(self.data)
                #self.printBlock(eData)

                #print "self.tableCode = ",self.tableCode
                #print "len(data) = ",len(self.data)
                
                gridDevice, = unpack('i',eData[0:4])
                #print "gridDevice = ",gridDevice
                if gridDevice==0:
                    grid = 'C'
                else:
                    #grid = (gridDevice - deviceCode) // 10
                    grid = gridDevice
                ###

                out = unpack('ffffffffffffffffffff',eData[4:4*21])
                if self.makeOp2Debug:
                    self.op2Debug.write('%s\n' %(str(out)))
                (sxx,sxy,s1,a1,a2,a3,pressure,svm,
                 syy,syz,s2,b1,b2,b3,
                 szz,sxz,s3,c1,c2,c3) = out
                
                #print "%s eid=%s cid=%s grid=%s sxx=%-6i txy=%-5i s1=%-6i a1=%i a2=%i a3=%i press=%i vm=%s" %(elementType,eid,cid,grid,sxx,sxy,s1,a1,a2,a3,pressure,svm)
                #print "%s eid=%s cid=%s grid=%s syy=%-6i tyz=%-5i s2=%-6i b1=%i b2=%i b3=%i"                %(elementType,eid,cid,grid,syy,syz,s2,b1,b2,b3)
                #print "%s eid=%s cid=%s grid=%s szz=%-6i txz=%-5i s3=%-6i c1=%i c2=%i c3=%i"                %(elementType,eid,cid,grid,szz,sxz,s3,c1,c2,c3)
                #print ""

                #smax = max(s1,s2,s3)
                #smin = min(s1,s2,s3)
                
                aCos = []
                bCos = []
                cCos = []
                if nodeID==0:
                    #print "adding new eid"
                    self.obj.addNewEid(ElementType,cid,eid,grid,sxx,syy,szz,sxy,syz,sxz,s1,s2,s3,aCos,bCos,cCos,pressure,svm)
                else:
                    self.obj.add(                      eid,grid,sxx,syy,szz,sxy,syz,sxz,s1,s2,s3,aCos,bCos,cCos,pressure,svm)
                #print "eid=%i grid=%i fd1=%i sx1=%i sy1=%i txy1=%i angle1=%i major1=%i minor1=%i vm1=%i" %(eid,grid,fd1,sx1,sy1,txy1,angle1,major1,minor1,vm1)
                #print "               fd2=%i sx2=%i sy2=%i txy2=%i angle2=%i major2=%i minor2=%i vm2=%i\n"          %(fd2,sx2,sy2,txy2,angle2,major2,minor2,vm2)
                #self.printBlock(data)
            #sys.exit('finished a CEHXA')
            #print self.solidStress[self.iSubcase]
            ###
            #print '--------------------'
            #print "len(data) = ",len(self.data)
            
            #self.printSection(100)
            #self.printBlock(self.data[0:100])
            #self.printBlock(self.data[1:100])
            #self.printBlock(self.data[2:100])
            #self.printBlock(self.data[3:100])
        ###
        self.handleResultsBuffer(self.CSOLID_67)
        #print self.solidStress[self.iSubcase]

    def CTRIA3_74(self): # works
        """
        DISTANCE,NORMAL-X,NORMAL-Y,SHEAR-XY,ANGLE,MAJOR,MINOR,VONMISES
        stress is extracted at the centroid
        """
        if self.makeOp2Debug:
            self.op2Debug.write('---CTRIA3_74---\n')
        deviceCode = self.deviceCode
        assert self.numWide==17,'invalid numWide...numWide=%s' %(self.numWide)
        while len(self.data)>=68:
            eData     = self.data[0:4*17]
            self.data = self.data[4*17: ]
            out = unpack('iffffffffffffffff',eData)

            (eid,fd1,sx1,sy1,txy1,angle1,major1,minor1,vm1,
                 fd2,sx2,sy2,txy2,angle2,major2,minor2,vm2,) = out
            eid = (eid - deviceCode) // 10
            #print "eid=%i fd1=%i sx1=%i sy1=%i txy1=%i angle1=%i major1=%i minor1=%i vm1=%i" %(eid,fd1,sx1,sy1,txy1,angle1,major1,minor1,vm1)
            #print  "      fd2=%i sx2=%i sy2=%i txy2=%i angle2=%i major2=%i minor2=%i vm2=%i\n"   %(fd2,sx2,sy2,txy2,angle2,major2,minor2,vm2)
            self.obj.addNewEid('CTRIA3',eid,'C',fd1,sx1,sy1,txy1,angle1,major1,minor1,vm1)
            self.obj.add(               eid,'C',fd2,sx2,sy2,txy2,angle2,major2,minor2,vm2)
            if self.makeOp2Debug:
                self.op2Debug.write('%s\n' %(str(out)))
        ###
        self.handleResultsBuffer(self.CTRIA3_74)

    def CSOLID_85(self):  # works
        """
        stress is extracted at the centroid
        CTETRA_85
        CPENTA_91 ???
        CHEXA_93  ???
        """
        if self.makeOp2Debug:
            self.op2Debug.write('---CSOLID_85---\n')
        #print "starting nonlinear solid element..."
        deviceCode = self.deviceCode
        #nNodes = 5 # 1 centroid + 4 corner points
        #self.printSection(20)
        #term      = self.data[0:4] CEN/
        #self.data = self.data[4:]
        #print "*****"

        nNodes=4  # this is a minimum, it will be reset later
        nNodesExpected = 1
        assert self.numWide in [82],'invalid numWide...numWide=%s' %(self.numWide)
        while len(self.data)>= 16+84*nNodesExpected:
            eData     = self.data[0:8]
            self.data = self.data[8:]
            #self.printBlock(eData)

            out = unpack("issss",eData)
            (eid,a,b,c,d) = out
            if self.makeOp2Debug:
                self.op2Debug.write('%s\n' %(str(out)))
            #print "abcd = |%s|" %(a+b+c+d)
            #print "eid=%s cid=%s nNodes=%s nNodesExpected=%s" %(eid,cid,nNodes,nNodesExpected)
            
            assert nNodes < 21,self.printBlock(eData)
            eid = (eid - deviceCode) // 10
            if(  nNodes in [4,10]):
                elementType = "CTETRA"
                nNodesExpected = 5
            #elif(nNodes in [6,15]):
            #    elementType = "CPENTA"
            #    nNodesExpected = 7
            #elif(nNodes in [8,20]):
            #    elementType = "CHEXA"
            #    nNodesExpected = 9
            #else:
            #    raise Exception('not supported....nNodes=%s' %(nNodes))

            #print "len(data) = ",len(self.data)
            for nodeID in range(nNodesExpected):   #nodes pts, +1 for centroid (???)
                #print "len(data)A = ",len(self.data)
                eData     = self.data[0:4*16]  # for the stresses
                self.data = self.data[4*16: ]
                #print "len(data)B = ",len(self.data)
                #self.printBlock(eData)

                #print "self.tableCode = ",self.tableCode
                
                #print "len(data) = ",len(self.data)
                
                #gridDevice, = unpack('i',eData[0:4])
                #print "gridDevice = ",gridDevice

                out = unpack('ifffffffffffffff',eData[:4*16])  # 18-3 = 15
                if self.makeOp2Debug:
                    self.op2Debug.write('%s\n' %(str(out)))

                (gridGauss,sxx,syy,szz,sxy,syz,sxz,se,
                   eps,ecs,exx,eyy,ezz,exy,eyz,exz) = out

                if gridGauss==0:
                    gridGauss = 'C'
                else:
                    #grid = (gridDevice - deviceCode) // 10
                    gridGauss = gridGauss
                ###

                #print "%s gridGauss=%-5s eid=%s sxx=%g syy=%g szz=%g" %(elementType,gridGauss,eid,sxx,syy,szz)
                
                #print "%s eid=%s cid=%s grid=%s sxx=%-6i txy=%-5i s1=%-6i a1=%i a2=%i a3=%i press=%i vm=%s" %(elementType,eid,cid,grid,sxx,sxy,s1,a1,a2,a3,pressure,svm)
                #print "%s eid=%s cid=%s grid=%s syy=%-6i tyz=%-5i s2=%-6i b1=%i b2=%i b3=%i"                %(elementType,eid,cid,grid,syy,syz,s2,b1,b2,b3)
                #print "%s eid=%s cid=%s grid=%s szz=%-6i txz=%-5i s3=%-6i c1=%i c2=%i c3=%i"                %(elementType,eid,cid,grid,szz,sxz,s3,c1,c2,c3)
                #print ""

                #smax = max(s1,s2,s3)
                #smin = min(s1,s2,s3)
                
                #if nodeID==0:
                #    #print "adding new eid"
                #    self.obj.addNewEid(elementType,cid,eid,grid,sxx,syy,szz,sxy,syz,sxz,aCos,bCos,cCos,pressure,svm)
                #else:
                #    self.obj.add(                      eid,grid,sxx,syy,szz,sxy,syz,sxz,aCos,bCos,cCos,pressure,svm)
                #print "eid=%i grid=%i fd1=%i sx1=%i sy1=%i txy1=%i angle1=%i major1=%i minor1=%i vm1=%i" %(eid,grid,fd1,sx1,sy1,txy1,angle1,major1,minor1,vm1)
                #print "               fd2=%i sx2=%i sy2=%i txy2=%i angle2=%i major2=%i minor2=%i vm2=%i\n"          %(fd2,sx2,sy2,txy2,angle2,major2,minor2,vm2)
                #self.printBlock(data)
            #sys.exit('finished a CTETRANL')
            #print self.solidStress[self.iSubcase]
            ###
            #print '--------------------'
            #print "len(data) = ",len(self.data)
            
            #self.printSection(100)
            #self.printBlock(self.data[0:100])
            #self.printBlock(self.data[1:100])
            #self.printBlock(self.data[2:100])
            #self.printBlock(self.data[3:100])
        ###
        self.handleResultsBuffer(self.CSOLID_85)
        #print self.solidStress[self.iSubcase]

    #def extractDt(self):
    #    pass
    #def extractInt(self):
    #    return 'i',self.scaleEid
    #def extractFloat(self):
    #    return 'f',self.scaleDt
    #def scaleEid(self,eid):
    #    return (eid-self.deviceCode)//10
    #def scaleDt(self,dt):
    #    return dt

    def OES_field1(self):
        if self.isSort1():
            #raise NotImplementedError('SORT1 is not supported')
            return 'i',self.scaleEid
        elif self.isSort2():
            if self.analysisCode in [1,2,3,4,7,8,9,11,12]: # eid
                return 'i',self.scaleEid
            elif self.analysisCode==5: # freq
                #freq
                return 'f',self.scaleDt
                #raise NotImplementedError('freq is not supported')
            elif self.analysisCode==6: # time
                #time
                return 'f',self.scaleDt
                #raise NotImplementedError('time is not supported')
            elif self.analysisCode==10: # fqts:
                #fqts # freqTime
                return 'f',self.scaleDt
                #raise NotImplementedError('freqTime is not supported')
            else:
                raise NotImplementedError('invalid SORT2 analysisCode=%s' %(self.analysisCode))
            ###
        else:
            raise NotImplementedError('invalid SORTx code')
        ###

    def CTRIAX6_53(self):
        (Format1,scaleValue) = self.OES_field1()
        Format = Format1+'ifffffff'
        while len(self.data)>=132: # (1+8*4) = 33*4 = 132
            eData     = self.data[0:4*9]
            self.data = self.data[4*9: ]
            out = unpack(Format,eData)
            (eid,loc,rs,azs,As,ss,maxp,tmax,octs) = out
            eid = scaleValue(eid)
            #print "eid=%s loc=%s rs=%s azs=%s as=%s ss=%s maxp=%s tmx=%s octs=%s" %(eid,loc,rs,azs,As,ss,maxp,tmax,octs)
            self.obj.addNewEid(eid,loc,rs,azs,As,ss,maxp,tmax,octs)

            for i in range(3):
                eData     = self.data[0:4*8]
                self.data = self.data[4*8: ]
                out = unpack('ifffffff',eData)
                (loc,rs,azs,As,ss,maxp,tmax,octs) = out
                #print "eid=%s loc=%s rs=%s azs=%s as=%s ss=%s maxp=%s tmx=%s octs=%s" %(eid,loc,rs,azs,As,ss,maxp,tmax,octs)
                self.obj.add(eid,loc,rs,azs,As,ss,maxp,tmax,octs)

            if self.makeOp2Debug:
                self.op2Debug.write('%s\n' %(str(out)))
        ###
        self.handleResultsBuffer(self.CTRIAX6_53)

    def RODNL_89_92(self):
        while len(self.data)>=28:
            eData     = self.data[0:4*7]
            self.data = self.data[4*7: ]
            out = unpack('iffffff',eData)

            (eid,axial,equivStress,totalStrain,effPlasticCreepStrain,effCreepStrain,linearTorsionalStresss) = out
            eid = (eid - self.deviceCode) // 10
            data = (eid,axial,equivStress,totalStrain,effPlasticCreepStrain,effCreepStrain,linearTorsionalStresss)
            
            #print "eid=%s axial=%s equivStress=%s totalStrain=%s effPlasticCreepStrain=%s effCreepStrain=%s linearTorsionalStresss=%s" %(eid,axial,equivStress,totalStrain,effPlasticCreepStrain,effCreepStrain,linearTorsionalStresss)
            #print "eid=%i fd1=%i sx1=%i sy1=%i txy1=%i angle1=%i major1=%i minor1=%i vm1=%i" %(eid,fd1,sx1,sy1,txy1,angle1,major1,minor1,vm1)
            #print  "      fd2=%i sx2=%i sy2=%i txy2=%i angle2=%i major2=%i minor2=%i vm2=%i\n"   %(fd2,sx2,sy2,txy2,angle2,major2,minor2,vm2)
            #self.obj.addNewEid('CTRIA3',eid,'C',fd1,sx1,sy1,txy1,angle1,major1,minor1,vm1)
            self.obj.add(self.elementType,data)
            
            if self.makeOp2Debug:
                self.op2Debug.write('%s\n' %(str(out)))
        ###
        self.handleResultsBuffer(self.RODNL_89_92)

    def CQUAD4NL_90(self):
        #print "self.numWide = ",self.numWide
        
        deviceCode = self.deviceCode
        numWide = self.numWide
        if numWide==13:
            nStep = 52  #4*13
        elif numWide==25:
            nStep = 100 #4*25
        else:
            raise Exception('invalid CQUAD4NL_90')

        while len(self.data)>=nStep:
            eData     = self.data[0:nStep]
            self.data = self.data[nStep: ]
            
            if numWide==13:
                out = unpack('iffffffffffff',eData) # numWide=13
                (eid,fd1,sx1,sy1,sz1,txy1,es1,eps1,ecs1,ex1,ey1,ez1,exy1) = out
                eid = (eid - deviceCode) // 10
                data = (eid,fd1,sx1,sy1,sz1,txy1,es1,eps1,ecs1,ex1,ey1,ez1,exy1)
                self.obj.addNewEid(self.elementType,data)
            elif numWide==25:
                out = unpack('iffffffffffffffffffffffff',eData) # numWide=25
                (eid,fd1,sx1,sy1,xxx,txy1,es1,eps1,ecs1,ex1,ey1,xxx,exy1,
                     fd2,sx2,sy2,xxx,txy2,es2,eps2,ecs2,ex2,ey2,xxx,exy2) = out
                eid = (eid - deviceCode) // 10

                data = (eid,fd1,sx1,sy1,xxx,txy1,es1,eps1,ecs1,ex1,ey1,xxx,exy1)
                self.obj.addNewEid(self.elementType,data)
                data = (eid,fd2,sx2,sy2,xxx,txy2,es2,eps2,ecs2,ex2,ey2,xxx,exy2)
                self.obj.add(data)
            
            #print "eid=%s axial=%s equivStress=%s totalStrain=%s effPlasticCreepStrain=%s effCreepStrain=%s linearTorsionalStresss=%s" %(eid,axial,equivStress,totalStrain,effPlasticCreepStrain,effCreepStrain,linearTorsionalStresss)
            
            if self.makeOp2Debug:
                self.op2Debug.write('%s\n' %(str(out)))
        ###
        self.handleResultsBuffer(self.CQUAD4NL_90)
    
    def CPENTANL_91(self):
        """
        The DMAP manual says fields 3-18 repeat 7 times. but they dont.
        They repeat 6 times.  Other DMAP cards are correct with
        their repeat statements.
        """
        #print "CHEXANL_93"
        #print "len(self.data) = ",len(self.data)

        n = 0
        while len(self.data)>=456: # 2+16*7 = 114 -> 114*4 = 456
            eData = self.data[0:8]
            self.data = self.data[8:]
            (eid,a,b,c,d) = unpack('icccc',eData)
            eid = (eid - self.deviceCode) // 10
            #out = unpack("ii",eData)
            #(eid,cType) = out
            cType = a+b+c+d

            for i in range(7):
                #print "len(self.data) = ",len(self.data)
                eData = self.data[0:64]
                self.data = self.data[64:]
                out = unpack('ifffffffffffffff',eData)
                assert len(out)==16
                (grid,sx,sy,sz,sxy,syz,sxz,se,eps,ecs,ex,ey,ez,exy,eyz,exz) = out
                #print "eid=%3s cType=%s sx=%i sy=%i sz=%i sxy=%s syz=%i szx=%i se=%s" %(eid,cType,sx,sy,sz,sxy,syz,sxz,se)
                #print "gid=%3s ecs=%.3g   ex=%.3g ey=%.3g ez=%.3g exy=%.3g eyz=%.3g ezx=%.3g"  %(grid,ecs,ex,ey,ez,exy,eyz,exz)
                #print ""
                assert a=='G'
            
            #self.data = self.data[1456:]
            #sys.exit('hexa...')
        #print "buffer time..."
        #self.firstPass = True
        self.handleResultsBuffer(self.CHEXANL_93,debug=True)

    def CHEXANL_93(self):
        """
        The DMAP manual says fields 3-18 repeat 9 times. but they dont.
        They repeat 8 times.  Other DMAP cards are correct with
        their repeat statements.
        """
        #print "CHEXANL_93"
        #print "len(self.data) = ",len(self.data)

        n = 0
        while len(self.data)>=584: # 2+16*9 = 146 -> 146*4 = 584
            eData = self.data[0:8]
            self.data = self.data[8:]
            (eid,a,b,c,d) = unpack('icccc',eData)
            eid = (eid - self.deviceCode) // 10
            #out = unpack("ii",eData)
            #(eid,cType) = out
            cType = a+b+c+d

            for i in range(9):
                #print "len(self.data) = ",len(self.data)
                eData = self.data[0:64]
                self.data = self.data[64:]
                out = unpack('ifffffffffffffff',eData)
                assert len(out)==16
                (grid,sx,sy,sz,sxy,syz,sxz,se,eps,ecs,ex,ey,ez,exy,eyz,exz) = out
                #print "eid=%3s cType=%s sx=%i sy=%i sz=%i sxy=%s syz=%i szx=%i se=%s" %(eid,cType,sx,sy,sz,sxy,syz,sxz,se)
                #print "gid=%3s ecs=%.3g   ex=%.3g ey=%.3g ez=%.3g exy=%.3g eyz=%.3g ezx=%.3g"  %(grid,ecs,ex,ey,ez,exy,eyz,exz)
                #print ""
                assert a=='G'
            
            #self.data = self.data[1456:]
            #sys.exit('hexa...')
        #print "buffer time..."
        #self.firstPass = True
        self.handleResultsBuffer(self.CHEXANL_93,debug=True)

    def CBEAM_94(self):
        if self.makeOp2Debug:
            self.op2Debug.write('---BEAM_94---\n')
        deviceCode = self.deviceCode
        nNodes = 10 # 11-1

        #nTotal       = self.obj.getLengthTotal()
        #(n1,format1) = self.obj.getLength1()
        #(n2,format2) = self.obj.getLength2()
        nTotal = 2*4 + (18-3)*9*4
        nTotal = 204

        n1 = 24
        format1 = 'ssssfffff'
        #print "len(data) = ",len(self.data)
        while len(self.data)>=nTotal:
            eData     = self.data[0:8]
            self.data = self.data[8: ]
            (eid,gridA) = unpack('ii', eData)
            #print "eid=%s gridA=%s" %(eid,gridA)

            for i in range(1):
                for j in range(4): # c,d,e,f @ A;    c,d,e,f @ B
                    eData     = self.data[0:n1]
                    self.data = self.data[n1: ]
                    out = unpack(format1, eData)
                    (loc1,loc2,loc3,loc4,nsx,nse,te,epe,ece) = out
                    loc = loc1+loc2+loc3+loc4
                    #print "loc=%s nsx=%s nse=%s te=%s epe=%s ece=%s" %(loc,nsx,nse,te,epe,ece)
                ###
                #self.obj.add(eid,out)
            ###
            #sys.exit('stoping in CBEAM_94')

        self.handleResultsBuffer(self.CBEAM_94)
        #raise Exception('add CBEAM-94...')

    def CQUAD4_95(self): # works (doesnt handle all stress/strain cases tho)
        """
        GRID-ID  DISTANCE,NORMAL-X,NORMAL-Y,SHEAR-XY,ANGLE,MAJOR MINOR,VONMISES
        """
        deviceCode = self.deviceCode
        eType = self.ElementType(self.elementType)

        #self.printSection(20)
        #term = data[0:4] CEN/
        #data = data[4:]
        assert self.numWide==11,'invalid numWide...numWide=%s' %(self.numWide)

        while len(self.data)>=44: # 2+17*5 = 87 -> 87*4 = 348
            eData     = self.data[0:4*11]
            self.data = self.data[4*11: ]
            out = unpack('iifffffffff',eData)
            if self.makeOp2Debug:
                self.op2Debug.write('%s\n' %(str(out)))
            (eid,iLayer,o1,o2,t12,t1z,t2z,angle,major,minor,ovm) = out
            #print "out =",out
            eid = (eid - deviceCode) // 10  ## @todo adjust with deviceCode...
            
            if eid!=self.eid2: # originally initialized to None, the buffer doesnt reset it, so it is the old value
                #print "1 - eid=%s iLayer=%i o1=%i o2=%i ovm=%i" %(eid,iLayer,o1,o2,ovm)
                self.obj.addNewEid(eType,eid,o1,o2,t12,t1z,t2z,angle,major,minor,ovm)
            else:
                #print "2 - eid=%s iLayer=%i o1=%i o2=%i ovm=%i" %(eid,iLayer,o1,o2,ovm)
                self.obj.add(eid,o1,o2,t12,t1z,t2z,angle,major,minor,ovm)
            ###
            self.eid2 = eid
            #self.dn += 348
        ###
        #print "3 - eid=%s iLayer=%i o1=%i o2=%i ovm=%i" %(eid,iLayer,o1,o2,ovm)
        self.printSection(100)
        self.handleResultsBuffer(self.CQUAD4_95)

    def QUAD4FD_139(self): # hyperelastic
        """
        Hyperelastic Quad
        36+4*7*4 = 148
        """
        #x = 0
        while len(self.data)>=148:
            #if x==2:
            #    sys.exit('end of hyperQuad')
            eData     = self.data[0:4*9]
            self.data = self.data[4*9: ]
            out = unpack('issssiffffff',eData)

            (eid,t1,t2,t3,t4,ID,sx,sy,sxy,angle,smj,smi) = out
            eid = (eid - self.deviceCode)//10
            Type = t1+t2+t3+t4
            self.obj.addNewEid([eid,Type,sx,sy,sxy,angle,smj,smi])
            #print "eid=%s Type=%s\n***ID=%s sx=%s sy=%s sxy=%s angle=%s major=%s minor=%s" %(eid,Type,ID,sx,sy,sxy,angle,smj,smi)
            for i in range(3):
                eData     = self.data[0:4*7]
                self.data = self.data[4*7: ]
                out = unpack('iffffff',eData)
                #(ID,sx,sy,sxy,angle,smj,smi) = out
                self.obj.add(eid,out)
                #print "***ID=%s sx=%s sy=%s sxy=%s angle=%s major=%s minor=%s" %(ID,sx,sy,sxy,angle,smj,smi)
            ###
            #self.obj.add(data)
            #x+=1
            if self.makeOp2Debug:
                self.op2Debug.write('%s\n' %(str(out)))
        ###
        #sys.exit('end of hyperQuad')
        self.handleResultsBuffer(self.QUAD4FD_139)
    
    def CQUADR_82(self): # not done...
        """
        GRID-ID  DISTANCE,NORMAL-X,NORMAL-Y,SHEAR-XY,ANGLE,MAJOR MINOR,VONMISES
        """
        if self.makeOp2Debug:
            self.op2Debug.write('---CQUADR_82---\n')
        deviceCode = self.deviceCode

        if self.elementType==82: # CQUADR
            nTotal = 348  # 2+17*5 = 87 -> 87*4 = 348
            nNodes = 4    # centroid + 4 corner points
            eType  = 'CQUADR'
        else:
            raise Exception('elementType=%s nTotal not defined...' %(self.elementType))

        while len(self.data)>=nTotal:
            (eid,_,_,_,_) = unpack("issss",self.data[0:8])
            self.data = self.data[8:]  # 2
            eid = (eid - deviceCode) // 10
            eData     = self.data[0:4*17]
            self.data = self.data[4*17: ]
            out = unpack('iffffffffffffffff',eData)  # len=17*4
            if self.makeOp2Debug:
                self.op2Debug.write('%s\n' %(str(out)))
            (grid,fd1,sx1,sy1,txy1,angle1,major1,minor1,vm1,
                  fd2,sx2,sy2,txy2,angle2,major2,minor2,vm2,) = out
            grid = 'C'
            self.obj.addNewEid(eType,eid,grid,fd1,sx1,sy1,txy1,angle1,major1,minor1,vm1)
            self.obj.add(            eid,grid,fd2,sx2,sy2,txy2,angle2,major2,minor2,vm2)

            for nodeID in range(nNodes):   #nodes pts
                eData     = self.data[0:4*17]
                self.data = self.data[4*17: ]
                out = unpack('iffffffffffffffff',eData)
                if self.makeOp2Debug:
                    self.op2Debug.write('%s\n' %(str(out)))
                (grid,fd1,sx1,sy1,txy1,angle1,major1,minor1,vm1,
                      fd2,sx2,sy2,txy2,angle2,major2,minor2,vm2,) = out

                #print "eid=%i grid=%i fd1=%i sx1=%i sy1=%i txy1=%i angle1=%i major1=%i minor1=%i vm1=%i" %(eid,grid,fd1,sx1,sy1,txy1,angle1,major1,minor1,vm1)
                #print "               fd2=%i sx2=%i sy2=%i txy2=%i angle2=%i major2=%i minor2=%i vm2=%i\n"          %(fd2,sx2,sy2,txy2,angle2,major2,minor2,vm2)
                #print "len(data) = ",len(self.data)
                #self.printBlock(self.data)
                self.obj.addNewNode(eid,grid,fd1,sx1,sy1,txy1,angle1,major1,minor1,vm1)
                self.obj.add(       eid,grid,fd2,sx2,sy2,txy2,angle2,major2,minor2,vm2)
            ###
        ###
        self.handleResultsBuffer(self.CQUAD4_82)

    def CQUAD4_144(self): # works
        """
        GRID-ID  DISTANCE,NORMAL-X,NORMAL-Y,SHEAR-XY,ANGLE,MAJOR MINOR,VONMISES
        """
        if self.makeOp2Debug:
            self.op2Debug.write('---CQUAD4_144---\n')
        deviceCode = self.deviceCode

        #self.printSection(20)
        #term = data[0:4] CEN/
        #data = data[4:]
        #print "*****"
        #self.printBlock(self.data)
        #assert self.numWide==87,'invalid numWide...numWide=%s' %(self.numWide)
        #if self.numWide==87: # 2+(17-1)*5 = 87 -> 87*4 = 348
        
        if self.elementType==144: # CQUAD4
            nTotal = 348  # 2+17*5 = 87 -> 87*4 = 348
            nNodes = 4    # centroid + 4 corner points
            eType  = 'CQUAD4'
        elif self.elementType==64: # CQUAD8
            nTotal = 348  # 2+17*5 = 87 -> 87*4 = 348
            nNodes = 4    # centroid + 4 corner points
            eType  = 'CQUAD8'
        elif self.elementType==75: # CTRIA6
            nTotal = 280  # 2+17*3 = 70 -> 70*4 = 280
            nNodes = 3    # centroid + 3 corner points
            eType  = 'CTRIA6'
        elif self.elementType==70: # CTRIAR
            nTotal = 280  # 2+17*3 = 70 -> 70*4 = 280
            nNodes = 3    # centroid + 3 corner points
            eType  = 'CTRIAR'
        else:
            raise Exception('elementType=%s nTotal not defined...' %(self.elementType))

        while len(self.data)>=nTotal:
            (eid,_,_,_,_) = unpack("issss",self.data[0:8])
            self.data = self.data[8:]  # 2
            eid = (eid - deviceCode) // 10
            eData     = self.data[0:4*17]
            self.data = self.data[4*17: ]
            out = unpack('iffffffffffffffff',eData)  # len=17*4
            if self.makeOp2Debug:
                self.op2Debug.write('%s\n' %(str(out)))
            (grid,fd1,sx1,sy1,txy1,angle1,major1,minor1,vm1,
                  fd2,sx2,sy2,txy2,angle2,major2,minor2,vm2,) = out
            grid = 'C'
            self.obj.addNewEid(eType,eid,grid,fd1,sx1,sy1,txy1,angle1,major1,minor1,vm1)
            self.obj.add(            eid,grid,fd2,sx2,sy2,txy2,angle2,major2,minor2,vm2)

            for nodeID in range(nNodes):   #nodes pts
                eData     = self.data[0:4*17]
                self.data = self.data[4*17: ]
                out = unpack('iffffffffffffffff',eData)
                if self.makeOp2Debug:
                    self.op2Debug.write('%s\n' %(str(out)))
                (grid,fd1,sx1,sy1,txy1,angle1,major1,minor1,vm1,
                      fd2,sx2,sy2,txy2,angle2,major2,minor2,vm2,) = out

                #print "eid=%i grid=%i fd1=%i sx1=%i sy1=%i txy1=%i angle1=%i major1=%i minor1=%i vm1=%i" %(eid,grid,fd1,sx1,sy1,txy1,angle1,major1,minor1,vm1)
                #print "               fd2=%i sx2=%i sy2=%i txy2=%i angle2=%i major2=%i minor2=%i vm2=%i\n"          %(fd2,sx2,sy2,txy2,angle2,major2,minor2,vm2)
                #print "len(data) = ",len(self.data)
                #self.printBlock(self.data)
                self.obj.addNewNode(eid,grid,fd1,sx1,sy1,txy1,angle1,major1,minor1,vm1)
                self.obj.add(       eid,grid,fd2,sx2,sy2,txy2,angle2,major2,minor2,vm2)
            ###
            #print '--------------------'
            #print "len(data) = ",len(self.data)
            #print "tell = ",self.op2.tell()

            #self.printSection(100)
            #sys.exit('asdf')
            #self.dn += 348
        ###
        ###
       # elif self.numWide==77:
       #     while len(self.data)>=308: # 2+15*5 = 77 -> 77*4 = 308
       #         (eid,_,_,_,_) = unpack("issss",self.data[0:8])
       #         self.data = self.data[8:]  # 2
       #         eid = (eid - deviceCode) // 10
       #         eData     = self.data[0:4*15]
       #         self.data = self.data[4*15: ]
       #
       #         out = unpack('iffffffffffffff',eData)  # 15
       #         if self.makeOp2Debug:
       #             self.op2Debug.write('%s\n' %(str(out)))
       #         (grid,fd1,sx1r,sx11,sy1r,sy11,txy1r,txy11,
       #               fd2,sx2r,sx21,sy2r,sy21,txy2r,txy21) = out
       #
       #         for nodeID in range(nNodes):   #nodes pts
       #             eData     = self.data[0:4*15]
       #             self.data = self.data[4*15: ]
       #             out = unpack('iffffffffffffff',eData)
       #             if self.makeOp2Debug:
       #                 self.op2Debug.write('%s\n' %(str(out)))
       #             (grid,fd1,sx1r,sx11,sy1r,sy11,txy1r,txy11,
       #                   fd2,sx2r,sx21,sy2r,sy21,txy2r,txy21) = out
       #         ###
       #     ###
       # ###
       # elif self.numWide==47:
       #     while len(self.data)>=188: # 2+9*5 = 47 -> 47*4 = 188
       #         (eid,_,_,_,_) = unpack("issss",self.data[0:8])
       #         self.data = self.data[8:]  # 2
       #         eid = (eid - deviceCode) // 10
       #         eData     = self.data[0:4*9]
       #         self.data = self.data[4*9: ]
       #
       #         out = unpack('iffffffff',eData)  # 9
       #         if self.makeOp2Debug:
       #             self.op2Debug.write('%s\n' %(str(out)))
       #         (grid,fd1,sx1,sy1,txy1,
       #               fd2,sx2,sy2,txy2) = out
       #
       #         for nodeID in range(nNodes):   #nodes pts
       #             eData     = self.data[0:4*9]
       #             self.data = self.data[4*9: ]
       #             out = unpack('iffffffff',eData)
       #             if self.makeOp2Debug:
       #                 self.op2Debug.write('%s\n' %(str(out)))
       #             (grid,fd1,sx1,sy1,txy1,
       #                   fd2,sx2,sy2,txy2) = out
       #             ###
       #         ###
       #     ###
       # ###
       # else:
       #     raise AddNewElementError('invalid numWide')
        self.handleResultsBuffer(self.CQUAD4_144)

