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
            print "%s" %(self.ElementType(self.elementType)),dataIn
            #eid = self.obj.addNewEid(out)
            self.obj.add(dt,dataIn)
            #print "len(data) = ",len(self.data)
        ###
        self.handleResultsBuffer(self.OEF_Rod_alt)
        #print self.rodForces
        
    def OEF_Beam_alt(self): # 2-CBEAM
        deviceCode = self.deviceCode
        raise NotImplementedError()
        dt = self.nonlinearFactor
        isSort1 = self.isSort1()
        if isSort1:
            #print "SORT1 - %s" %(self.ElementType(self.elementType))
            format1 = 'iiffffffff' # SORT1
            extract = self.extractSort1
            #dt = self.nonlinearFactor
        else:
            #print "SORT2 - %s" %(self.ElementType(self.elementType))
            format1 = 'fiffffffff' # SORT2
            extract = self.extractSort2
            #eid = self.nonlinearFactor

        self.createThermalTransientObject(self.beamForces,ComplexCBEAMForce,isSort1)

        while len(self.data)>=68: # 10*4
            eData     = self.data[0:68]
            self.data = self.data[68: ]
            #print "len(data) = ",len(eData)

            out = unpack(format1, eData)
            (eidTemp,nid,sd,bm1,bm2,ts1,ts2,af,ttrq,wtrq) = out
            #print "eidTemp = ",eidTemp
            #print "nid = ",nid
            #print "sd = ",sd
            if nid==0 or sd>0.:
                eid = self.eidOld
                isNewElement = False
            else:
                eid = eidTemp
                self.eidOld = eid
                isNewElement = True
            eid2  = extract(eid,dt)
            #print "eType=%s" %(eType)
            
            dataIn = [eid2,nid,sd,bm1,bm2,ts1,ts2,af,ttrq,wtrq]
            #print "%s" %(self.ElementType(self.elementType)),dataIn
            #eid = self.obj.addNewEid(out)
            if isNewElement:
                self.obj.addNewElement(dt,dataIn)
                #print
            elif sd>0.:
                self.obj.add(dt,dataIn)
                #print
            #else: pass

            #print "len(data) = ",len(self.data)
        ###
        self.handleResultsBuffer(self.OEF_Beam_alt)
        #print self.beamForces

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
            print "%s" %(self.ElementType(self.elementType)),dataIn
            #eid = self.obj.addNewEid(out)
            self.obj.add(dt,dataIn)
            #print "len(data) = ",len(self.data)
        ###
        self.handleResultsBuffer(self.OEF_Spring_alt)
        #print self.springForces
        
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
        print "nTotal",nTotal,nTotal/4.
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
