from __future__ import division
import sys
from struct import unpack
from .oef_thermalObjects import (HeatFlux_CHBDYx, HeatFlux_2D_3D, HeatFlux_1D,
                                 HeatFlux_VU, HeatFlux_VUBEAM, HeatFlux_VU_3D,
                                 HeatFlux_CONV)

class ThermalElements(object):

    def readOEF_Thermal(self):
        assert self.formatCode==1,self.codeInformation()
        #print "self.elementType = ",self.elementType
        if self.elementType in [107,108,109]: # CHBDYE, CHBDYG, CHBDYP
            assert self.numWide==8,self.codeInformation()
            self.createTransientObject(self.thermalLoad_CHBDY,HeatFlux_CHBDYx)
            self.OEF_CHBDYx()
        elif self.elementType in [33,39,67,68]: # QUAD4,TETRA,HEXA,PENTA
            assert self.numWide in [9,10],self.codeInformation()
            self.createTransientObject(self.thermalLoad_2D_3D,HeatFlux_2D_3D)
            self.OEF_2D_3D()
        elif self.elementType in [53,64,74,75]: # TRIAX6,QUAD8,TRIA3,TRIA6
            assert self.numWide==9,self.codeInformation()
            self.createTransientObject(self.thermalLoad_2D_3D,HeatFlux_2D_3D)
            self.OEF_2D_3D()
        elif self.elementType in [1,2,3,10,34,69]: # ROD,BEAM,TUBE,CONROD,BAR,BEND
            assert self.numWide==9,self.codeInformation()
            self.createTransientObject(self.thermalLoad_1D,HeatFlux_1D)
            self.OEF_1D()
        elif self.elementType in [189,190]: # VUQUAD,VUTRIA
            #assert self.numWide==27,self.codeInformation()
            self.createTransientObject(self.thermalLoad_VU,HeatFlux_VU)
            self.OEF_VU_Element()
        elif self.elementType in [191]: # VUBEAM
            #assert self.numWide==27,self.codeInformation()
            self.createTransientObject(self.thermalLoad_VUBeam,HeatFlux_VUBEAM)
            self.OEF_VUBeam_Element()
        elif self.elementType in [145,146,147]: # VUHEXA,VUPENTA,VUTETRA
            self.createTransientObject(self.thermalLoad_VU_3D,HeatFlux_VU_3D)
            self.OEF_VU_3D_Element()
        elif self.elementType in [110]:
            self.createTransientObject(self.thermalLoad_CONV,HeatFlux_CONV)
            self.OEF_CONV()
        else:
            #print self.codeInformation()
            raise NotImplementedError(self.codeInformation())
        ###

    def OEF_CHBDYx(self): # [107,108,109]  CHBDYE, CHBDYG, CHBDYP
        if self.makeOp2Debug:
            self.op2Debug.write('---OEF_CHBDYx---\n')
        #deviceCode = self.deviceCode
        
        dt = self.nonlinearFactor
        isSort1 = self.isSort1()
        if isSort1:
            #print "SORT1 - %s" %(self.ElementType(self.elementType))
            format1 = 'iccccccccfffff' # SORT1
            extract = self.extractSort1
            #dt = self.nonlinearFactor
        else:
            #print "SORT2 - %s" %(self.ElementType(self.elementType))
            format1 = 'fccccccccfffff' # SORT2
            extract = self.extractSort2
            #eid = self.nonlinearFactor

        while len(self.data)>=32: # 8*4
            eData     = self.data[0:32]
            self.data = self.data[32: ]
            #print "len(data) = ",len(eData)

            out = unpack(format1, eData)
            (eid,a,b,c,d,e,f,g,h,fApplied,freeConv,forceConv,fRad,fTotal) = out
            eid2  = extract(eid,dt)
            eType = a+b+c+d+e+f+g+h
            #print "eType=%s" %(eType)
            
            dataIn = [eid2,eType,fApplied,freeConv,forceConv,fRad,fTotal]
            #print "heatFlux %s" %(self.ElementType(self.elementType)),dataIn
            #eid = self.obj.addNewEid(out)
            self.obj.add(dt,dataIn)
            #print "len(data) = ",len(self.data)
        ###
        self.handleResultsBuffer(self.OEF_CHBDYx)
        if self.makeOp2Debug:
            print("done with OEF_CHBDYx")
        #print(self.thermalLoad_CHBDY)

    def OEF_CONV(self): # [110]  CONV
        if self.makeOp2Debug:
            self.op2Debug.write('---OEF_CONV---\n')
        #deviceCode = self.deviceCode
        
        dt = self.nonlinearFactor
        isSort1 = self.isSort1()
        if isSort1:
            #print "SORT1 - %s" %(self.ElementType(self.elementType))
            format1 = 'ifif' # SORT1
            extract = self.extractSort1
            #dt = self.nonlinearFactor
        else:
            #print "SORT2 - %s" %(self.ElementType(self.elementType))
            format1 = 'ffif' # SORT2
            extract = self.extractSort2
            #eid = self.nonlinearFactor

        while len(self.data)>=16: # 4*4
            eData     = self.data[0:16]
            self.data = self.data[16: ]
            #print "len(data) = ",len(eData)

            out = unpack(format1, eData)
            (eid,cntlNode,freeConv,freeConvK) = out
            eid2  = extract(eid,dt)
            
            dataIn = [eid2,cntlNode,freeConv,freeConvK]
            #print "heatFlux %s" %(self.ElementType(self.elementType)),dataIn
            self.obj.add(dt,dataIn)
            #print "len(data) = ",len(self.data)
        ###
        self.handleResultsBuffer(self.OEF_CONV)
        if self.makeOp2Debug:
            print("done with OEF_CONV")
        #print(self.thermalLoad_CHBDY)

    def OEF_VU_Element(self): # 189-VUQUAD 190-VUTRIA,191-VUBEAM
        dt = self.nonlinearFactor
        isSort1 = self.isSort1()
        #print "numWide = ",self.numWide

        if self.elementType in [189]:
            nNodes = 4
        elif self.elementType in [190]:
            nNodes = 3
        elif self.elementType in [191]:
            nNodes = 2
        else:
            raise NotImplementedError(self.codeInformation())

        if isSort1:
            #print "SORT1 - %s" %(self.ElementType(self.elementType))
            format1 = 'iiiccccii' # SORT1
            extract = self.extractSort1
            #dt = self.nonlinearFactor
        else:
            #print "SORT2 - %s" %(self.ElementType(self.elementType))
            format1 = 'fiiccccii' # SORT2
            extract = self.extractSort2
            #eid = self.nonlinearFactor
        ###
        formatAll = 'iffffff'

        n = 24+28*nNodes
        while len(self.data)>=n:
            eData     = self.data[0:24] # 6*4
            self.data = self.data[24: ]
            #print "len(data) = ",len(eData)

            out = unpack(format1, eData)
            (eid,parent,coord,icordA,icordB,icordC,icordD,theta,null) = out
            icord = icordA+icordB+icordC+icordD

            eid2  = extract(eid,dt)
            dataIn = [eid2,parent,coord,icord,theta]

            gradFluxes = []
            for i in range(nNodes):
                eData     = self.data[0:28] # 7*4
                self.data = self.data[28: ]
                #print "i=%s len(data)=%s" %(i,len(eData))
                out = unpack(formatAll, eData)
                gradFluxes.append(out)
            dataIn.append(gradFluxes)
            #eType = a+b+c+d+e+f+g+h
            #print "eType=%s" %(eType)
            
            #dataIn = [eid2,eType,xGrad,yGrad,zGrad,xFlux,yFlux,zFlux]
            #print "heatFlux %s" %(self.ElementType(self.elementType)),dataIn
            #eid = self.obj.addNewEid(out)            
            self.obj.add(nNodes,dt,dataIn)
            #print "len(data) = ",len(self.data)
        ###
        self.handleResultsBuffer(self.OEF_VU_Element)
        if self.makeOp2Debug:
            print("done with OEF_1D")
        #sys.exit('thermalLoad_VU')
        #print self.thermalLoad_VU

    def OEF_VUBeam_Element(self): # 191-VUBEAM
        dt = self.nonlinearFactor
        isSort1 = self.isSort1()
        #print "numWide = ",self.numWide

        if self.elementType in [191]:
            nNodes = 2
        else:
            raise NotImplementedError(self.codeInformation())

        if isSort1:
            #print "SORT1 - %s" %(self.ElementType(self.elementType))
            format1 = 'iiicccc' # SORT1
            extract = self.extractSort1
            #dt = self.nonlinearFactor
        else:
            #print "SORT2 - %s" %(self.ElementType(self.elementType))
            format1 = 'fiicccc' # SORT2
            extract = self.extractSort2
            #eid = self.nonlinearFactor
        ###
        formatAll = 'iffffff'

        n = 16+28*nNodes
        while len(self.data)>=n:
            eData     = self.data[0:16] # 4*4
            self.data = self.data[16: ]
            #print "len(data) = ",len(eData)

            out = unpack(format1, eData)
            (eid,parent,coord,icordA,icordB,icordC,icordD) = out
            icord = icordA+icordB+icordC+icordD

            eid2  = extract(eid,dt)
            dataIn = [eid2,parent,coord,icord]

            gradFluxes = []
            for i in range(nNodes):
                eData     = self.data[0:28] # 7*4
                self.data = self.data[28: ]
                #print "i=%s len(data)=%s" %(i,len(eData))
                out = unpack(formatAll, eData)
                gradFluxes.append(out)
            dataIn.append(gradFluxes)
            #eType = a+b+c+d+e+f+g+h
            #print "eType=%s" %(eType)
            
            #dataIn = [eid2,eType,xGrad,yGrad,zGrad,xFlux,yFlux,zFlux]
            #print "heatFlux %s" %(self.ElementType(self.elementType)),dataIn
            #eid = self.obj.addNewEid(out)            
            self.obj.add(nNodes,dt,dataIn)
            #print "len(data) = ",len(self.data)
        ###
        self.handleResultsBuffer(self.OEF_VUBeam_Element)
        if self.makeOp2Debug:
            print("done with OEF_1D")
        #sys.exit('thermalLoad_VUBeam')
        #print self.thermalLoad_VUBeam

    def OEF_VU_3D_Element(self): # 146-VUPENTA, 147-VUTETRA, 148-VUPENTA
        dt = self.nonlinearFactor
        isSort1 = self.isSort1()
        #print "numWide = ",self.numWide

        if self.elementType in [147]:
            nNodes = 4
        elif self.elementType in [146]:
            nNodes = 6
        elif self.elementType in [145]:
            nNodes = 8
        else:
            raise NotImplementedError(self.codeInformation())

        if isSort1:
            #print "SORT1 - %s" %(self.ElementType(self.elementType))
            format1 = 'ii' # SORT1
            extract = self.extractSort1
            #dt = self.nonlinearFactor
        else:
            #print "SORT2 - %s" %(self.ElementType(self.elementType))
            format1 = 'fi' # SORT2
            extract = self.extractSort2
            #eid = self.nonlinearFactor
        ###
        formatAll = 'iffffff'

        n = 8+28*nNodes
        while len(self.data)>=n:
            eData     = self.data[0:8] # 2*4
            self.data = self.data[8: ]
            #print "len(data) = ",len(eData)

            out = unpack(format1, eData)
            (eid,parent) = out

            eid2  = extract(eid,dt)
            dataIn = [eid2,parent]

            gradFluxes = []
            for i in range(nNodes):
                eData     = self.data[0:7*4]
                self.data = self.data[7*4: ]
                #print "i=%s len(data)=%s" %(i,len(eData))
                out = unpack(formatAll, eData)
                gradFluxes.append(out)
            dataIn.append(gradFluxes)
            
            #print "heatFlux %s" %(self.ElementType(self.elementType)),dataIn
            self.obj.add(nNodes,dt,dataIn)
            #print "len(data) = ",len(self.data)
        ###
        self.handleResultsBuffer(self.OEF_VU_3D_Element)
        if self.makeOp2Debug:
            print("done with OEF_VU_3D_Element")
        #print self.thermalLoad_VU_3D

    def OEF_1D(self): # 1-ROD, 2-BEAM, 3-TUBE, 10-CONROD, 34-BAR, 69-BEND
        dt = self.nonlinearFactor
        isSort1 = self.isSort1()
        #print "numWide = ",self.numWide
        n = 36
        if isSort1:
            #print "SORT1 - %s" %(self.ElementType(self.elementType))
            format1 = 'iccccccccffffff' # SORT1
            extract = self.extractSort1
            #dt = self.nonlinearFactor
        else:
            #print "SORT2 - %s" %(self.ElementType(self.elementType))
            format1 = 'fccccccccffffff' # SORT2
            extract = self.extractSort2
            #eid = self.nonlinearFactor
        ###

        n = 36
        while len(self.data)>=n: # 10*4
            eData     = self.data[0:n]
            self.data = self.data[n: ]
            #print "len(data) = ",len(eData)

            out = unpack(format1, eData)
            (eid,a,b,c,d,e,f,g,h,xGrad,yGrad,zGrad,xFlux,yFlux,zFlux) = out
            eid2  = extract(eid,dt)
            eType = a+b+c+d+e+f+g+h
            #print "eType=%s" %(eType)
            
            dataIn = [eid2,eType,xGrad,yGrad,zGrad,xFlux,yFlux,zFlux]
            #print "heatFlux %s" %(self.ElementType(self.elementType)),dataIn
            #eid = self.obj.addNewEid(out)            
            self.obj.add(dt,dataIn)
            #print "len(data) = ",len(self.data)
        ###
        self.handleResultsBuffer(self.OEF_1D)
        if self.makeOp2Debug:
            print("done with OEF_1D")
        #print self.thermalLoad_1D

    def OEF_2D_3D(self): # 33-QUAD4, 39-TETRA, 53-TRIAX6,64-QUAD8, 67-HEXA, 68-PENTA, 74-TRIA3, 75-TRIA6
        """numWide==10"""
        dt = self.nonlinearFactor
        isSort1 = self.isSort1()
        #print "numWide = ",self.numWide
        #print "dt = ",dt
        if self.elementType in [39,67,68]: # HEXA,PENTA
            n = 40
            if isSort1:
                #print "SORT1 - %s" %(self.ElementType(self.elementType))
                format1 = 'iccccccccffffffi' # SORT1
                extract = self.extractSort1
                #dt = self.nonlinearFactor
            else:
                #print "SORT2 - %s" %(self.ElementType(self.elementType))
                format1 = 'fccccccccffffffi' # SORT2
                extract = self.extractSort2
                #eid = self.nonlinearFactor
        elif self.elementType in [33,53,64,74,75]: # no zed on this element for some reason...
            n = 36
            if isSort1:
                #print "SORT1 - %s" %(self.ElementType(self.elementType))
                format1 = 'iccccccccffffff' # SORT1
                extract = self.extractSort1
                #dt = self.nonlinearFactor
            else:
                #print "SORT2 - %s" %(self.ElementType(self.elementType))
                format1 = 'fccccccccffffff' # SORT2
                extract = self.extractSort2
                #eid = self.nonlinearFactor
            ###
        else:
            raise NotImplementedError(self.codeInformation())
        ###

        #n = 36
        while len(self.data)>=n: # 10*4
            eData     = self.data[0:n]
            self.data = self.data[n: ]

            out = unpack(format1, eData)
            (eid,a,b,c,d,e,f,g,h,xGrad,yGrad,zGrad,xFlux,yFlux,zFlux) = out[0:15]
            eid2  = extract(eid,dt)
            eType = a+b+c+d+e+f+g+h
            #print "eType=%s" %(eType)
            
            dataIn = [eid2,eType,xGrad,yGrad,zGrad,xFlux,yFlux,zFlux]
            #print "heatFlux %s" %(self.ElementType(self.elementType)),dataIn
            #eid = self.obj.addNewEid(out)            
            self.obj.add(dt,dataIn)
            #print "len(data) = ",len(self.data)
        ###
        self.handleResultsBuffer(self.OEF_2D_3D)
        if self.makeOp2Debug:
            print("done with OEF_2D_3D")
        
        #print self.thermalLoad_2D_3D

