from __future__ import division
import sys
from struct import unpack

from oef_complexForceObjects import *

class ComplexForces(object):

    def OEF_Rod_alt(self): # 1-CROD, 3-CTUBE, 10-CONROD
        deviceCode = self.deviceCode
        
        dt = self.nonlinearFactor
        isSort1 = self.isSort1()
        if isSort1:
            #print "SORT1 - %s" %(self.ElementType(self.elementType))
            format1 = 'iffff' # SORT1
            extract = self.extractSort1
            #dt = self.nonlinearFactor
        else:
            #print "SORT2 - %s" %(self.ElementType(self.elementType))
            format1 = 'fffff' # SORT2
            extract = self.extractSort2
            #eid = self.nonlinearFactor

        self.createThermalTransientObject(self.rodForces,ComplexRodForce,isSort1)

        while len(self.data)>=20: # 5*4
            eData     = self.data[0:20]
            self.data = self.data[20: ]
            #print "len(data) = ",len(eData)

            out = unpack(format1, eData)
            (eid,axialReal,torqueReal,axialImag,torqueImag) = out
            eid2  = extract(eid,dt)
            #print "eType=%s" %(eType)
            
            dataIn = [eid2,axialReal,torqueReal,axialImag,torqueImag]
            #print "%s" %(self.ElementType(self.elementType)),dataIn
            #eid = self.obj.addNewEid(out)
            self.obj.add(dt,dataIn)
            #print "len(data) = ",len(self.data)
        ###
        self.handleResultsBuffer(self.OEF_Rod_alt)
        #print self.rodForces
        
    def OEF_Beam_alt(self): # 2-CBEAM
        deviceCode = self.deviceCode
        
        dt = self.nonlinearFactor
        isSort1 = self.isSort1()
        if isSort1:
            #print "SORT1 - %s" %(self.ElementType(self.elementType))
            format1 = 'i' # SORT1
            extract = self.extractSort1
            #dt = self.nonlinearFactor
        else:
            #print "SORT2 - %s" %(self.ElementType(self.elementType))
            format1 = 'f' # SORT2
            extract = self.extractSort2
            #eid = self.nonlinearFactor

        self.createThermalTransientObject(self.beamForces,ComplexCBEAMForce,isSort1)
        #print self.codeInformation()

        #nTotal = 16*11+1
        formatAll = 'ifffffffffffffff'
        while len(self.data)>=708: # (16*11+1)*4 = 177*4
            eData     = self.data[0:4]
            self.data = self.data[4: ]
            eidTemp, = unpack(format1, eData)
            eid2  = extract(eidTemp,dt)

            for i in range(11):
                eData     = self.data[0:64]
                self.data = self.data[64: ]
                #print "len(data) = ",len(eData)

                out = unpack(formatAll, eData)
                (nid,sd,bm1r,bm2r,ts1r,ts2r,afr,ttrqr,wtrqr,
                        bm1i,bm2i,ts1i,ts2i,afi,ttrqi,wtrqi) = out
                #print "eidTemp = ",eidTemp
                #print "nid = ",nid
                #print "sd = ",sd

                #eid = self.obj.addNewEid(out)
                if i==0: #isNewElement:
                    dataIn = [eid2,nid,sd,bm1r,bm2r,ts1r,ts2r,afr,ttrqr,wtrqr,
                                          bm1i,bm2i,ts1i,ts2i,afi,ttrqi,wtrqi]
                    print "%s cNew   " %(self.ElementType(self.elementType)),dataIn
                    self.obj.addNewElement(dt,dataIn)
                    #print
                elif sd>0.:
                    dataIn = [eid2,nid,sd,bm1r,bm2r,ts1r,ts2r,afr,ttrqr,wtrqr,
                                          bm1i,bm2i,ts1i,ts2i,afi,ttrqi,wtrqi]
                    print "%s cOld   " %(self.ElementType(self.elementType)),dataIn
                    self.obj.add(dt,dataIn)
                    #print
                #else: pass

            #print "len(data) = ",len(self.data)
        ###
        self.handleResultsBuffer(self.OEF_Beam_alt)
        #print self.beamForces

    def OEF_Shear_alt(self): # 4-CSHEAR
        deviceCode = self.deviceCode
        
        dt = self.nonlinearFactor
        isSort1 = self.isSort1()
        if isSort1:
            #print "SORT1 - %s" %(self.ElementType(self.elementType))
            format1 = 'iffffffffffffffffffffffffffffffff' # SORT1
            extract = self.extractSort1
            #dt = self.nonlinearFactor
        else:
            #print "SORT2 - %s" %(self.ElementType(self.elementType))
            format1 = 'fffffffffffffffffffffffffffffffff' # SORT2
            extract = self.extractSort2
            #eid = self.nonlinearFactor

        self.createThermalTransientObject(self.shearForces,ComplexCShearForce,isSort1)

        while len(self.data)>=132: # 33*4
            eData     = self.data[0:132]
            self.data = self.data[132: ]
            #print "len(data) = ",len(eData)

            out = unpack(format1, eData)
            (eid,f41r,f21r,f12r,f32r,f23r,f43r,f34r,f14r,kf1r,s12r,kf2r,s23r,kf3r,s34r,kf4r,s41r,
                 f41i,f21i,f12i,f32i,f23i,f43i,f34i,f14i,kf1i,s12i,kf2i,s23i,kf3i,s34i,kf4i,s41i) = out
            eid2  = extract(eid,dt)
            #print "eType=%s" %(eType)
            
            dataIn = [eid2,f41r,f21r,f12r,f32r,f23r,f43r,f34r,f14r,kf1r,s12r,kf2r,s23r,kf3r,s34r,kf4r,s41r,
                           f41i,f21i,f12i,f32i,f23i,f43i,f34i,f14i,kf1i,s12i,kf2i,s23i,kf3i,s34i,kf4i,s41i]
            print "%s" %(self.ElementType(self.elementType)),dataIn
            #eid = self.obj.addNewEid(out)
            self.obj.add(dt,dataIn)
            #print "len(data) = ",len(self.data)
        ###
        self.handleResultsBuffer(self.OEF_Shear_alt)
        #print self.shearForces
        
    def OEF_Spring_alt(self): # 11-CELAS1, 12-CELAS2, 13-CELAS3, 14-CELAS4
        deviceCode = self.deviceCode
        
        dt = self.nonlinearFactor
        isSort1 = self.isSort1()
        if isSort1:
            #print "SORT1 - %s" %(self.ElementType(self.elementType))
            format1 = 'iff' # SORT1
            extract = self.extractSort1
            #dt = self.nonlinearFactor
        else:
            #print "SORT2 - %s" %(self.ElementType(self.elementType))
            format1 = 'fff' # SORT2
            extract = self.extractSort2
            #eid = self.nonlinearFactor

        self.createThermalTransientObject(self.springForces,ComplexSpringForce,isSort1)

        while len(self.data)>=12: # 3*4
            eData     = self.data[0:12]
            self.data = self.data[12: ]
            #print "len(data) = ",len(eData)

            out = unpack(format1, eData)
            (eid,forceReal,forceImag) = out
            eid2  = extract(eid,dt)
            #print "eType=%s" %(eType)
            
            dataIn = [eid2,forceReal,forceImag]
            #print "%s" %(self.ElementType(self.elementType)),dataIn
            #eid = self.obj.addNewEid(out)
            self.obj.add(dt,dataIn)
            #print "len(data) = ",len(self.data)
        ###
        self.handleResultsBuffer(self.OEF_Spring_alt)
        #print self.springForces
        
    def OEF_CVisc_alt(self): # 24-CVISC
        deviceCode = self.deviceCode
        
        dt = self.nonlinearFactor
        isSort1 = self.isSort1()
        if isSort1:
            #print "SORT1 - %s" %(self.ElementType(self.elementType))
            format1 = 'iffff' # SORT1
            extract = self.extractSort1
            #dt = self.nonlinearFactor
        else:
            #print "SORT2 - %s" %(self.ElementType(self.elementType))
            format1 = 'fffff' # SORT2
            extract = self.extractSort2
            #eid = self.nonlinearFactor

        self.createThermalTransientObject(self.viscForces,ComplexViscForce,isSort1)

        while len(self.data)>=20: # 5*4
            eData     = self.data[0:20]
            self.data = self.data[20: ]
            #print "len(data) = ",len(eData)

            out = unpack(format1, eData)
            (eid,axialReal,torqueReal,axialImag,torqueImag) = out
            eid2  = extract(eid,dt)
            #print "eType=%s" %(eType)
            
            dataIn = [eid2,axialReal,torqueReal,axialImag,torqueImag]
            #print "%s" %(self.ElementType(self.elementType)),dataIn
            #eid = self.obj.addNewEid(out)
            self.obj.add(dt,dataIn)
            #print "len(data) = ",len(self.data)
        ###
        self.handleResultsBuffer(self.OEF_CVisc_alt)
        #print self.viscForces
        
    def OEF_CBar_alt(self): # 34-CBAR
        deviceCode = self.deviceCode
        
        dt = self.nonlinearFactor
        isSort1 = self.isSort1()
        if isSort1:
            #print "SORT1 - %s" %(self.ElementType(self.elementType))
            format1 = 'iffffffffffffffff' # SORT1
            extract = self.extractSort1
            #dt = self.nonlinearFactor
        else:
            #print "SORT2 - %s" %(self.ElementType(self.elementType))
            format1 = 'fffffffffffffffff' # SORT2
            extract = self.extractSort2
            #eid = self.nonlinearFactor

        self.createThermalTransientObject(self.barForces,ComplexCBARForce,isSort1)

        while len(self.data)>=68: # 17*4
            eData     = self.data[0:68]
            self.data = self.data[68: ]
            #print "len(data) = ",len(eData)

            out = unpack(format1, eData)
            (eid,bm1ar,bm2ar,bm1br,bm2br,ts1r,ts2r,afr,trqr,
                 bm1ai,bm2ai,bm1bi,bm2bi,ts1i,ts2i,afi,trqi) = out
            eid2  = extract(eid,dt)
            #print "eType=%s" %(eType)
            
            dataIn = [eid2,bm1ar,bm2ar,bm1br,bm2br,ts1r,ts2r,afr,trqr,
                           bm1ai,bm2ai,bm1bi,bm2bi,ts1i,ts2i,afi,trqi]
            print "%s" %(self.ElementType(self.elementType)),dataIn
            #eid = self.obj.addNewEid(out)
            self.obj.add(dt,dataIn)
            #print "len(data) = ",len(self.data)
        ###
        self.handleResultsBuffer(self.OEF_CBar_alt)
        #print self.barForces
        
    def OEF_Plate_alt(self): # 33-CQUAD4,74-CTRIA3
        deviceCode = self.deviceCode
        
        dt = self.nonlinearFactor
        isSort1 = self.isSort1()
        if isSort1:
            #print "SORT1 - %s" %(self.ElementType(self.elementType))
            format1 = 'iffffffffffffffff' # SORT1
            extract = self.extractSort1
            #dt = self.nonlinearFactor
        else:
            #print "SORT2 - %s" %(self.ElementType(self.elementType))
            format1 = 'fffffffffffffffff' # SORT2
            extract = self.extractSort2
            #eid = self.nonlinearFactor

        self.createThermalTransientObject(self.plateForces,ComplexPlateForce,isSort1)

        while len(self.data)>=68: # 17*4
            eData     = self.data[0:68]
            self.data = self.data[68: ]
            #print "len(data) = ",len(eData)

            out = unpack(format1, eData)
            (eid,mxr,myr,mxyr,bmxr,bmyr,bmxyr,txr,tyr,
                 mxi,myi,mxyi,bmxi,bmyi,bmxyi,txi,tyi) = out
            eid2  = extract(eid,dt)
            #print "eType=%s" %(eType)
            
            dataIn = [eid2,mxr,myr,mxyr,bmxr,bmyr,bmxyr,txr,tyr,
                           mxi,myi,mxyi,bmxi,bmyi,bmxyi,txi,tyi]
            print "%s" %(self.ElementType(self.elementType)),dataIn
            #eid = self.obj.addNewEid(out)
            self.obj.add(dt,dataIn)
            #print "len(data) = ",len(self.data)
        ###
        self.handleResultsBuffer(self.OEF_Plate_alt)
        #print self.plateForces

    def OEF_Plate2_alt(self): # 64-CQUAD8,70-CTRIAR,75-CTRIA6,82-CQUAD8,144-CQUAD4-bilinear
        deviceCode = self.deviceCode
        
        dt = self.nonlinearFactor
        
        if self.elementType in [70,75]: # CTRIAR,CTRIA6
            nNodes = 4
        elif self.elementType in [64,82,144]: # CQUAD8,CQUADR,CQUAD4-bilinear
            nNodes = 5
        else:
            raise NotImplementedError(self.codeInformation())
        ###
            
        isSort1 = self.isSort1()
        if isSort1:
            #print "SORT1 - %s" %(self.ElementType(self.elementType))
            format1 = 'icccc' # SORT1
            extract = self.extractSort1
            #dt = self.nonlinearFactor
        else:
            #print "SORT2 - %s" %(self.ElementType(self.elementType))
            format1 = 'fcccc' # SORT2
            extract = self.extractSort2
            #eid = self.nonlinearFactor

        self.createThermalTransientObject(self.plateForces2,ComplexPLATE2Force,isSort1)

        allFormat = 'fffffffffffffffff'
        nTotal = 8+nNodes*68
        while len(self.data)>=nTotal:
            eData     = self.data[0:76]
            self.data = self.data[76: ]
            #print self.printBlock(eData)
            #print "len(data) = ",len(eData)

            out = unpack(format1+allFormat, eData)
            (eid,a,b,c,d,nid,mxr,myr,mxyr,bmxr,bmyr,bmxyr,tx,tyr,
                             mxi,myi,mxyi,bmxi,bmyi,bmxyi,tx,tyi) = out
            term = a+b+c+d # CEN\
            #print "eType=%s" %(eType)

            eid2  = extract(eid,dt)
            
            dataIn = [term,nid,mxr,myr,mxyr,bmxr,bmyr,bmxyr,tx,tyr,
                               mxi,myi,mxyi,bmxi,bmyi,bmxyi,tx,tyi]
            print "%s" %(self.ElementType(self.elementType)),dataIn
            self.obj.addNewElement(eid2,dt,dataIn)

            for i in range(nNodes-1):
                eData     = self.data[0:68]
                self.data = self.data[68: ]
                dataIn = unpack(allFormat, eData)
                #(nid,mx,my,mxy,bmx,bmy,bmxy,tx,ty) = out
                #dataIn = [nid,mx,my,mxy,bmx,bmy,bmxy,tx,ty]
                print "***%s    " %(self.ElementType(self.elementType)),dataIn
                
                self.obj.add(eid2,dt,dataIn)
                #print "len(data) = ",len(self.data)
            ###
        ###
        #sys.exit('Plate2 stop...')
        self.handleResultsBuffer(self.OEF_Plate2_alt)
        #print self.plateForces2

    def OEF_Bend_alt(self): # 69-CBEND
        deviceCode = self.deviceCode
        
        dt = self.nonlinearFactor
        isSort1 = self.isSort1()
        if isSort1:
            #print "SORT1 - %s" %(self.ElementType(self.elementType))
            format1 = 'iifffffffffffffffffffffffff' # SORT1
            extract = self.extractSort1
            #dt = self.nonlinearFactor
        else:
            #print "SORT2 - %s" %(self.ElementType(self.elementType))
            format1 = 'fifffffffffffffffffffffffff' # SORT2
            extract = self.extractSort2
            #eid = self.nonlinearFactor

        self.createThermalTransientObject(self.bendForces,ComplexBendForce,isSort1)

        while len(self.data)>=108: # 27*4
            eData     = self.data[0:108]
            self.data = self.data[108: ]
            #print "len(data) = ",len(eData)

            out = unpack(format1, eData)
            (eid,nidA,bm1Ar,bm2Ar,ts1Ar,ts2Ar,afAr,trqAr,
                      bm1Ai,bm2Ai,ts1Ai,ts2Ai,afAi,trqAi,
                 nidB,bm1Br,bm2Br,ts1Br,ts2Br,afBr,trqBr,
                      bm1Bi,bm2Bi,ts1Bi,ts2Bi,afBi,trqBi) = out
            eid2  = extract(eid,dt)
            #print "eType=%s" %(eType)
            
            dataIn = [eid2,nidA,bm1Ar,bm2Ar,ts1Ar,ts2Ar,afAr,trqAr,
                                bm1Ai,bm2Ai,ts1Ai,ts2Ai,afAi,trqAi,
                           nidB,bm1Br,bm2Br,ts1Br,ts2Br,afBr,trqBr,
                                bm1Bi,bm2Bi,ts1Bi,ts2Bi,afBi,trqBi]
            print "%s" %(self.ElementType(self.elementType)),dataIn
            #eid = self.obj.addNewEid(out)
            self.obj.add(dt,dataIn)
            #print "len(data) = ",len(self.data)
        ###
        self.handleResultsBuffer(self.OEF_Bend_alt)
        #print self.bendForces
        
    def OEF_PentaPressure_alt(self): # 76-CHEXA_PR,77-CPENTA_PR,78-CTETRA_PR
        deviceCode = self.deviceCode
        
        dt = self.nonlinearFactor
        isSort1 = self.isSort1()
        if isSort1:
            #print "SORT1 - %s" %(self.ElementType(self.elementType))
            format1 = 'iccccccccfffffffffffff' # SORT1
            extract = self.extractSort1
            #dt = self.nonlinearFactor
        else:
            #print "SORT2 - %s" %(self.ElementType(self.elementType))
            format1 = 'fccccccccfffffffffffff' # SORT2
            extract = self.extractSort2
            #eid = self.nonlinearFactor

        self.createThermalTransientObject(self.pentaPressureForces,ComplexPentaPressureForce,isSort1)

        while len(self.data)>=64: # 16*4
            eData     = self.data[0:64]
            self.data = self.data[64: ]
            #print "len(data) = ",len(eData)

            out = unpack(format1, eData)
            (eid,a,b,c,d,e,f,g,h,axr,ayr,azr,vxr,vyr,vzr,pressure,
                                 axi,ayi,azi,vxi,vyi,vzi) = out
            eName = a+b+c+d+e+f+g+h
            eid2  = extract(eid,dt)
            #print "eType=%s" %(eType)
            
            dataIn = [eid2,eName,axr,ayr,azr,vxr,vyr,vzr,pressure,
                                 axi,ayi,azi,vxi,vyi,vzi]
            print "%s" %(self.ElementType(self.elementType)),dataIn
            #eid = self.obj.addNewEid(out)
            self.obj.add(dt,dataIn)
            #print "len(data) = ",len(self.data)
        ###
        self.handleResultsBuffer(self.OEF_PentaPressure_alt)
        #print self.bendForces

    def OEF_CBush_alt(self): # 102-CBUSH
        deviceCode = self.deviceCode
        
        dt = self.nonlinearFactor
        isSort1 = self.isSort1()
        if isSort1:
            #print "SORT1 - %s" %(self.ElementType(self.elementType))
            format1 = 'iffffffffffff' # SORT1
            extract = self.extractSort1
            #dt = self.nonlinearFactor
        else:
            #print "SORT2 - %s" %(self.ElementType(self.elementType))
            format1 = 'fffffffffffff' # SORT2
            extract = self.extractSort2
            #eid = self.nonlinearFactor

        self.createThermalTransientObject(self.bushForces,ComplexCBUSHForce,isSort1)

        while len(self.data)>=52: # 13*4
            eData     = self.data[0:52]
            self.data = self.data[52: ]
            #print "len(data) = ",len(eData)

            out = unpack(format1, eData)
            (eid,fxr,fyr,fzr,mxr,myr,mzr,
                 fxi,fyi,fzi,mxi,myi,mzi) = out
            eid2  = extract(eid,dt)
            #print "eType=%s" %(eType)
            
            dataIn = [eid2,fxr,fyr,fzr,mxr,myr,mzr,
                           fxi,fyi,fzi,mxi,myi,mzi]
            print "%s" %(self.ElementType(self.elementType)),dataIn
            #eid = self.obj.addNewEid(out)
            self.obj.add(dt,dataIn)
            #print "len(data) = ",len(self.data)
        ###
        self.handleResultsBuffer(self.OEF_CBush_alt)
        #print self.bushForces

    def OEF_Force_VU_alt(self): # 191-VUBEAM
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
        formatAll = 'ifffffffffffff'
        self.createThermalTransientObject(self.force_VU,ComplexForce_VU,isSort1)

        n = 16+56*nNodes
        while len(self.data)>=n:
            eData     = self.data[0:16] # 8*4
            self.data = self.data[16: ]
            #print "len(data) = ",len(eData)

            out = unpack(format1, eData)
            (eid,parent,coord,icordA,icordB,icordC,icordD) = out
            icord = icordA+icordB+icordC+icordD

            eid2  = extract(eid,dt)
            dataIn = [eid2,parent,coord,icord]

            forces = []
            for i in range(nNodes):
                eData     = self.data[0:56] # 14*4
                self.data = self.data[56: ]
                #print "i=%s len(data)=%s" %(i,len(eData))
                out = unpack(formatAll, eData)
                forces.append(out)
            dataIn.append(forces)
            #eType = a+b+c+d+e+f+g+h
            #print "eType=%s" %(eType)
            
            #dataIn = [vugrid,posit,forceX,shearY,shearZ,torsion,bendY,bendZ]
            print "force %s" %(self.ElementType(self.elementType)),dataIn
            #eid = self.obj.addNewEid(out)
            self.obj.add(nNodes,dt,dataIn)
            #print "len(data) = ",len(self.data)
        ###
        self.handleResultsBuffer(self.OEF_Force_VU_alt)
        if self.makeOp2Debug:
            print "done with OEF_Force_VU"
        print self.force_VU

    def OEF_Force_VUTRIA_alt(self): # 189-VUQUAD,190-VUTRIA
        dt = self.nonlinearFactor
        isSort1 = self.isSort1()
        #print "numWide = ",self.numWide

        if self.elementType in [189]: # VUQUAD
            nNodes = 4
        elif self.elementType in [190]: # VUTRIA
            nNodes = 3
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
        formatAll = 'ifffiiifffffifffiiifffffi'
        self.createThermalTransientObject(self.force_VU_2D,ComplexForce_VU_2D,isSort1)

        n = 24+100*nNodes
        while len(self.data)>=n:
            eData     = self.data[0:24] # 6*4
            self.data = self.data[24: ]
            #print "len(data) = ",len(eData)

            out = unpack(format1, eData)
            (eid,parent,coord,icordA,icordB,icordC,icordD,theta,_) = out
            icord = icordA+icordB+icordC+icordD

            eid2  = extract(eid,dt)
            dataIn = [eid2,parent,coord,icord,theta]

            forces = []
            for i in range(nNodes):
                eData     = self.data[0:100] # 13*4
                self.data = self.data[100: ]
                #print "i=%s len(data)=%s" %(i,len(eData))
                out = unpack(formatAll, eData)
                forces.append(out)
            dataIn.append(forces)
            #eType = a+b+c+d+e+f+g+h
            #print "eType=%s" %(eType)
            
            #dataIn = [vugrid,mfxr,mfyr,mfxyr,a,b,c,bmxr,bmyr,bmxyr,syzr,szxr,d,
                             #mfxi,mfyi,mfxyi,a,b,c,bmxi,bmyi,bmxyi,syzi,szxi,d]
            print "force %s" %(self.ElementType(self.elementType)),dataIn
            #eid = self.obj.addNewEid(out)
            self.obj.add(nNodes,dt,dataIn)
            #print "len(data) = ",len(self.data)
        ###
        self.handleResultsBuffer(self.OEF_Force_VUTRIA_alt)
        if self.makeOp2Debug:
            print "done with OEF_Force_VUTRIA"
        print self.force_VU_2D
