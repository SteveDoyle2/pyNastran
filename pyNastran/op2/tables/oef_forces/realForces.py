from __future__ import division
import sys
from struct import unpack

from oef_forceObjects import *

class RealForces(object):

    def OEF_Rod(self): # 1-CROD, 3-CTUBE, 10-CONROD
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

        self.createThermalTransientObject(self.rodForces,RealRodForce,isSort1)

        while len(self.data)>=12: # 3*4
            eData     = self.data[0:12]
            self.data = self.data[12: ]
            #print "len(data) = ",len(eData)

            out = unpack(format1, eData)
            (eid,axial,torque) = out
            eid2  = extract(eid,dt)
            #print "eType=%s" %(eType)
            
            dataIn = [eid2,axial,torque]
            #print "%s" %(self.ElementType(self.elementType)),dataIn
            #eid = self.obj.addNewEid(out)
            self.obj.add(dt,dataIn)
            #print "len(data) = ",len(self.data)
        ###
        self.handleResultsBuffer(self.OEF_Rod)
        #print self.rodForces
        
    def OEF_Beam(self): # 2-CBEAM
        deviceCode = self.deviceCode
        
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

        self.createThermalTransientObject(self.beamForces,RealCBEAMForce,isSort1)

        while len(self.data)>=40: # 10*4
            eData     = self.data[0:40]
            self.data = self.data[40: ]
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
            #print dataIn
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
        self.handleResultsBuffer(self.OEF_Beam)
        #print self.beamForces

    def OEF_Spring(self): # 11-CELAS1, 12-CELAS2, 13-CELAS3, 14-CELAS4
        deviceCode = self.deviceCode
        
        dt = self.nonlinearFactor
        isSort1 = self.isSort1()
        if isSort1:
            #print "SORT1 - %s" %(self.ElementType(self.elementType))
            format1 = 'if' # SORT1
            extract = self.extractSort1
            #dt = self.nonlinearFactor
        else:
            #print "SORT2 - %s" %(self.ElementType(self.elementType))
            format1 = 'ff' # SORT2
            extract = self.extractSort2
            #eid = self.nonlinearFactor

        self.createThermalTransientObject(self.springForces,RealSpringForce,isSort1)

        while len(self.data)>=8: # 2*4
            eData     = self.data[0:8]
            self.data = self.data[8: ]
            #print "len(data) = ",len(eData)

            out = unpack(format1, eData)
            (eid,force) = out
            eid2  = extract(eid,dt)
            #print "eType=%s" %(eType)
            
            dataIn = [eid2,force]
            #print dataIn
            #eid = self.obj.addNewEid(out)
            self.obj.add(dt,dataIn)
            #print "len(data) = ",len(self.data)
        ###
        self.handleResultsBuffer(self.OEF_Spring)
        #print self.springForces
        
    def OEF_CBar(self): # 34-CBAR
        deviceCode = self.deviceCode
        
        dt = self.nonlinearFactor
        isSort1 = self.isSort1()
        if isSort1:
            #print "SORT1 - %s" %(self.ElementType(self.elementType))
            format1 = 'iffffffff' # SORT1
            extract = self.extractSort1
            #dt = self.nonlinearFactor
        else:
            #print "SORT2 - %s" %(self.ElementType(self.elementType))
            format1 = 'fffffffff' # SORT2
            extract = self.extractSort2
            #eid = self.nonlinearFactor

        self.createThermalTransientObject(self.barForces,RealCBARForce,isSort1)

        while len(self.data)>=36: # 9*4
            eData     = self.data[0:36]
            self.data = self.data[36: ]
            #print "len(data) = ",len(eData)

            out = unpack(format1, eData)
            (eid,bm1a,bm2a,bm1b,bm2b,ts1,ts2,af,trq) = out
            eid2  = extract(eid,dt)
            #print "eType=%s" %(eType)
            
            dataIn = [eid2,bm1a,bm2a,bm1b,bm2b,ts1,ts2,af,trq]
            #print "%s" %(self.ElementType(self.elementType)),dataIn
            #eid = self.obj.addNewEid(out)
            self.obj.add(dt,dataIn)
            #print "len(data) = ",len(self.data)
        ###
        self.handleResultsBuffer(self.OEF_CBar)
        #print self.barForces
        
    def OEF_Plate(self): # 33-CQUAD4,74-CTRIA3
        deviceCode = self.deviceCode
        
        dt = self.nonlinearFactor
        isSort1 = self.isSort1()
        if isSort1:
            #print "SORT1 - %s" %(self.ElementType(self.elementType))
            format1 = 'iffffffff' # SORT1
            extract = self.extractSort1
            #dt = self.nonlinearFactor
        else:
            #print "SORT2 - %s" %(self.ElementType(self.elementType))
            format1 = 'fffffffff' # SORT2
            extract = self.extractSort2
            #eid = self.nonlinearFactor

        self.createThermalTransientObject(self.plateForces,RealPlateForce,isSort1)

        while len(self.data)>=36: # 9*4
            eData     = self.data[0:36]
            self.data = self.data[36: ]
            #print "len(data) = ",len(eData)

            out = unpack(format1, eData)
            (eid,mx,my,mxy,bmx,bmy,bmxy,tx,ty) = out
            eid2  = extract(eid,dt)
            #print "eType=%s" %(eType)
            
            dataIn = [eid2,mx,my,mxy,bmx,bmy,bmxy,tx,ty]
            print "%s" %(self.ElementType(self.elementType)),dataIn
            #eid = self.obj.addNewEid(out)
            self.obj.add(dt,dataIn)
            #print "len(data) = ",len(self.data)
        ###
        self.handleResultsBuffer(self.OEF_Plate)
        #print self.plateForces

    def OEF_Plate2(self): # 64-CQUAD8,70-CTRIAR,75-CTRIA6,82-CQUAD8,144-CQUAD4-bilinear
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

        self.createThermalTransientObject(self.plateForces2,RealPLATE2Force,isSort1)

        allFormat = 'fffffffff'
        nTotal = 44+nNodes*36
        while len(self.data)>=nTotal:
            eData     = self.data[0:44]
            self.data = self.data[44: ]
            #print self.printBlock(eData)
            #print "len(data) = ",len(eData)

            out = unpack(format1+allFormat, eData)
            (eid,a,b,c,d,nid,mx,my,mxy,bmx,bmy,bmxy,tx,ty) = out
            term = a+b+c+d # CEN\
            #print "eType=%s" %(eType)

            eid2  = extract(eid,dt)
            
            dataIn = [term,nid,mx,my,mxy,bmx,bmy,bmxy,tx,ty]
            print "%s" %(self.ElementType(self.elementType)),dataIn
            self.obj.addNewElement(eid2,dt,dataIn)

            for i in range(nNodes-1):
                eData     = self.data[0:36]
                self.data = self.data[36: ]
                dataIn = unpack(allFormat, eData)
                #(nid,mx,my,mxy,bmx,bmy,bmxy,tx,ty) = out
                #dataIn = [nid,mx,my,mxy,bmx,bmy,bmxy,tx,ty]
                print "%s    " %(self.ElementType(self.elementType)),dataIn
                
                self.obj.add(eid2,dt,dataIn)
                #print "len(data) = ",len(self.data)
            ###
        ###
        #sys.exit('Plate2 stop...')
        self.handleResultsBuffer(self.OEF_Plate2)
        #print self.plateForces2

    def OEF_CGap(self): # 38-CGAP
        deviceCode = self.deviceCode
        
        dt = self.nonlinearFactor
        isSort1 = self.isSort1()
        if isSort1:
            #print "SORT1 - %s" %(self.ElementType(self.elementType))
            format1 = 'iffffffff' # SORT1
            extract = self.extractSort1
            #dt = self.nonlinearFactor
        else:
            #print "SORT2 - %s" %(self.ElementType(self.elementType))
            format1 = 'fffffffff' # SORT2
            extract = self.extractSort2
            #eid = self.nonlinearFactor

        self.createThermalTransientObject(self.plateForces,RealCGAPForce,isSort1)

        while len(self.data)>=36: # 9*4
            eData     = self.data[0:36]
            self.data = self.data[36: ]
            #print "len(data) = ",len(eData)

            out = unpack(format1, eData)
            (eid,fx,sfy,sfz,u,v,w,sv,sw) = out
            eid2  = extract(eid,dt)
            #print "eType=%s" %(eType)
            
            dataIn = [eid2,fx,sfy,sfz,u,v,w,sv,sw]
            print "CGAP",dataIn
            #eid = self.obj.addNewEid(out)
            self.obj.add(dt,dataIn)
            #print "len(data) = ",len(self.data)
        ###
        self.handleResultsBuffer(self.OEF_CGap)
        #print self.plateForces

    def OEF_CBush(self): # 102-CBUSH
        deviceCode = self.deviceCode
        
        dt = self.nonlinearFactor
        isSort1 = self.isSort1()
        if isSort1:
            #print "SORT1 - %s" %(self.ElementType(self.elementType))
            format1 = 'iffffff' # SORT1
            extract = self.extractSort1
            #dt = self.nonlinearFactor
        else:
            #print "SORT2 - %s" %(self.ElementType(self.elementType))
            format1 = 'fffffff' # SORT2
            extract = self.extractSort2
            #eid = self.nonlinearFactor

        self.createThermalTransientObject(self.bushForces,RealCBUSHForce,isSort1)

        while len(self.data)>=28: # 7*4
            eData     = self.data[0:28]
            self.data = self.data[28: ]
            #print "len(data) = ",len(eData)

            out = unpack(format1, eData)
            (eid,fx,fy,fz,mx,my,mz) = out
            eid2  = extract(eid,dt)
            #print "eType=%s" %(eType)
            
            dataIn = [eid2,fx,fy,fz,mx,my,mz]
            #print "%s" %(self.ElementType(self.elementType)),dataIn
            #eid = self.obj.addNewEid(out)
            self.obj.add(dt,dataIn)
            #print "len(data) = ",len(self.data)
        ###
        self.handleResultsBuffer(self.OEF_CBush)
        #print self.bushForces
