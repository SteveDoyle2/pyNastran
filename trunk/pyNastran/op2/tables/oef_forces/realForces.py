from __future__ import division
import sys
from struct import unpack

from oef_forceObjects import *

class RealForces(object):

    def getOEF_FormatStart(self):
        """
        Returns an i or an f depending on if it's SORT2 or not.
        Also returns an extraction function that is called on the first argument
        """
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
        return (format1,extract)

    def OEF_Rod(self): # 1-CROD, 3-CTUBE, 10-CONROD
        deviceCode = self.deviceCode
        dt = self.nonlinearFactor

        (format1,extract) = self.getOEF_FormatStart()
        format1 += 'ff'
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
        
    def OEF_CVisc(self): # 24-CVISC
        deviceCode = self.deviceCode
        dt = self.nonlinearFactor

        (format1,extract) = self.getOEF_FormatStart()
        format1 += 'ff'
        self.createThermalTransientObject(self.viscForces,RealViscForce,isSort1)

        while len(self.data)>=12: # 3*4
            eData     = self.data[0:12]
            self.data = self.data[12: ]
            #print "len(data) = ",len(eData)

            out = unpack(format1, eData)
            (eid,axial,torque) = out
            eid2  = extract(eid,dt)
            #print "eType=%s" %(eType)
            
            dataIn = [eid2,axial,torque]
            print "%s" %(self.ElementType(self.elementType)),dataIn
            #eid = self.obj.addNewEid(out)
            self.obj.add(dt,dataIn)
            #print "len(data) = ",len(self.data)
        ###
        self.handleResultsBuffer(self.OEF_CVisc)
        #print self.viscForces
        
    def OEF_Beam(self): # 2-CBEAM
        deviceCode = self.deviceCode
        dt = self.nonlinearFactor

        (format1,extract) = self.getOEF_FormatStart()
        self.createThermalTransientObject(self.beamForces,RealCBEAMForce,isSort1)
        #print self.codeInformation()
        formatAll = 'iffffffff'
        while len(self.data)>=400: # 1+(10-1)*11=100 ->100*4 = 400
            #print "eType=%s" %(eType)

            eData     = self.data[0:4]
            self.data = self.data[4: ]
            eid, = unpack(format1, eData)
            eid2 = extract(eid,dt)

            for i in range(11):
                eData     = self.data[0:36]
                self.data = self.data[36: ]
                #print "len(data) = ",len(eData)

                out = unpack(formatAll, eData)
                (nid,sd,bm1,bm2,ts1,ts2,af,ttrq,wtrq) = out
                #print "eidTemp = ",eidTemp
                #print "nid = ",nid
                #print "sd = ",sd

                dataIn = [eid2,nid,sd,bm1,bm2,ts1,ts2,af,ttrq,wtrq]
                print "%s        " %(self.ElementType(self.elementType)),dataIn
                #eid = self.obj.addNewEid(out)
                if i==0: #isNewElement:
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

    def OEF_Shear(self): # 4-CSHEAR
        deviceCode = self.deviceCode
        dt = self.nonlinearFactor

        (format1,extract) = self.getOEF_FormatStart()
        format1 += 'ffffffffffffffff'
        self.createThermalTransientObject(self.shearForces,RealCShearForce,isSort1)

        while len(self.data)>=68: # 17*4
            eData     = self.data[0:68]
            self.data = self.data[68: ]
            #print "len(data) = ",len(eData)

            out = unpack(format1, eData)
            (eid,f41,f21,f12,f32,f23,f43,f34,f14,kf1,s12,kf2,s23,kf3,s34,kf4,s41) = out
            eid2  = extract(eid,dt)
            #print "eType=%s" %(eType)
            
            dataIn = [eid2,f41,f21,f12,f32,f23,f43,f34,f14,kf1,s12,kf2,s23,kf3,s34,kf4,s41]
            print "%s" %(self.ElementType(self.elementType)),dataIn
            #eid = self.obj.addNewEid(out)
            self.obj.add(dt,dataIn)
            #print "len(data) = ",len(self.data)
        ###
        self.handleResultsBuffer(self.OEF_Shear)
        #print self.shearForces
        
    def OEF_Spring(self): # 11-CELAS1, 12-CELAS2, 13-CELAS3, 14-CELAS4
        deviceCode = self.deviceCode
        dt = self.nonlinearFactor

        (format1,extract) = self.getOEF_FormatStart()
        format1 += 'f'
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
            #print "%s" %(self.ElementType(self.elementType)),dataIn
            #eid = self.obj.addNewEid(out)
            self.obj.add(dt,dataIn)
            #print "len(data) = ",len(self.data)
        ###
        self.handleResultsBuffer(self.OEF_Spring)
        #print self.springForces
        
    def OEF_CBar(self): # 34-CBAR
        deviceCode = self.deviceCode
        dt = self.nonlinearFactor

        (format1,extract) = self.getOEF_FormatStart()
        format1 += 'ffffffff'
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
        
    def OEF_CBar100(self): # 100-CBAR
        deviceCode = self.deviceCode
        dt = self.nonlinearFactor

        (format1,extract) = self.getOEF_FormatStart()
        format1 += 'fffffff'
        self.createThermalTransientObject(self.bar100Forces,RealCBAR100Force,isSort1)

        while len(self.data)>=36: # 9*4
            eData     = self.data[0:32]
            self.data = self.data[32: ]
            #print "len(data) = ",len(eData)

            out = unpack(format1, eData)
            (eid,sd,bm1,bm2,ts1,ts2,af,trq) = out
            eid2  = extract(eid,dt)
            #print "eType=%s" %(eType)
            
            dataIn = [eid2,sd,bm1,bm2,ts1,ts2,af,trq]
            print "%s" %(self.ElementType(self.elementType)),dataIn
            #eid = self.obj.addNewEid(out)
            self.obj.add(dt,dataIn)
            #print "len(data) = ",len(self.data)
        ###
        self.handleResultsBuffer(self.OEF_CBar)
        #print self.bar100Forces
        
    def OEF_Plate(self): # 33-CQUAD4,74-CTRIA3
        deviceCode = self.deviceCode
        dt = self.nonlinearFactor

        (format1,extract) = self.getOEF_FormatStart()
        format1 += 'ffffffff'
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

        (format1,extract) = self.getOEF_FormatStart()
        format1 += 'cccc'
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
                print "***%s    " %(self.ElementType(self.elementType)),dataIn
                
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

        (format1,extract) = self.getOEF_FormatStart()
        format1 += 'ffffffff'
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
            print "%s" %(self.ElementType(self.elementType)),dataIn
            #eid = self.obj.addNewEid(out)
            self.obj.add(dt,dataIn)
            #print "len(data) = ",len(self.data)
        ###
        self.handleResultsBuffer(self.OEF_CGap)
        #print self.plateForces

    def OEF_Bend(self): # 69-CBEND
        deviceCode = self.deviceCode
        dt = self.nonlinearFactor

        (format1,extract) = self.getOEF_FormatStart()
        format1 += 'ifffffffffffff'
        self.createThermalTransientObject(self.bendForces,RealBendForce,isSort1)

        while len(self.data)>=60: # 15*4
            eData     = self.data[0:60]
            self.data = self.data[60: ]
            #print "len(data) = ",len(eData)

            out = unpack(format1, eData)
            (eid,nidA,bm1A,bm2A,ts1A,ts2A,afA,trqA,
                 nidB,bm1B,bm2B,ts1B,ts2B,afB,trqB) = out
            eid2  = extract(eid,dt)
            #print "eType=%s" %(eType)
            
            dataIn = [eid2,nidA,bm1A,bm2A,ts1A,ts2A,afA,trqA,
                           nidB,bm1B,bm2B,ts1B,ts2B,afB,trqB]
            print "%s" %(self.ElementType(self.elementType)),dataIn
            #eid = self.obj.addNewEid(out)
            self.obj.add(dt,dataIn)
            #print "len(data) = ",len(self.data)
        ###
        self.handleResultsBuffer(self.OEF_Bend)
        #print self.bendForces
        
    def OEF_PentaPressure(self): # 77-CPENTA_PR,78-CTETRA_PR
        deviceCode = self.deviceCode
        dt = self.nonlinearFactor

        (format1,extract) = self.getOEF_FormatStart()
        format1 += 'ccccccccfffffff'
        self.createThermalTransientObject(self.pentaPressureForces,RealPentaPressureForce,isSort1)

        while len(self.data)>=40: # 10*4
            eData     = self.data[0:40]
            self.data = self.data[40: ]
            #print "len(data) = ",len(eData)

            out = unpack(format1, eData)
            (eid,a,b,c,d,e,f,g,h,ax,ay,az,vx,vy,vz,pressure) = out
            eName = a+b+c+d+e+f+g+h
            eid2  = extract(eid,dt)
            #print "eType=%s" %(eType)
            
            dataIn = [eid2,eName,ax,ay,az,vx,vy,vz,pressure]
            print "%s" %(self.ElementType(self.elementType)),dataIn
            #eid = self.obj.addNewEid(out)
            self.obj.add(dt,dataIn)
            #print "len(data) = ",len(self.data)
        ###
        self.handleResultsBuffer(self.OEF_PentaPressure)
        #print self.pentaPressureForces
        
    def OEF_CBush(self): # 102-CBUSH
        deviceCode = self.deviceCode
        
        dt = self.nonlinearFactor

        (format1,extract) = self.getOEF_FormatStart()
        format1 += 'ffffff'
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

    def OEF_Force_VU(self): # 191-VUBEAM
        dt = self.nonlinearFactor
        #print "numWide = ",self.numWide

        if self.elementType in [191]:
            nNodes = 2
        else:
            raise NotImplementedError(self.codeInformation())

        (format1,extract) = self.getOEF_FormatStart()
        format1 += 'iicccc'
        formatAll = 'ifffffff'
        self.createThermalTransientObject(self.force_VU,RealForce_VU,isSort1)

        n = 16+32*nNodes
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
                eData     = self.data[0:32] # 8*4
                self.data = self.data[32: ]
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
        self.handleResultsBuffer(self.OEF_Force_VU)
        if self.makeOp2Debug:
            print "done with OEF_Force_VU"
        print self.force_VU

    def OEF_Force_VUTRIA(self): # 189-VUQUAD,190-VUTRIA
        dt = self.nonlinearFactor
        #print "numWide = ",self.numWide

        if self.elementType in [189]: # VUQUAD
            nNodes = 4
        elif self.elementType in [190]: # VUTRIA
            nNodes = 3
        else:
            raise NotImplementedError(self.codeInformation())

        (format1,extract) = self.getOEF_FormatStart()
        format1 += 'iiccccii'
        formatAll = 'ifffiiifffffi'
        self.createThermalTransientObject(self.force_VU_2D,RealForce_VU_2D,isSort1)

        n = 24+52*nNodes
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
                eData     = self.data[0:52] # 13*4
                self.data = self.data[52: ]
                #print "i=%s len(data)=%s" %(i,len(eData))
                out = unpack(formatAll, eData)
                forces.append(out)
            dataIn.append(forces)
            #eType = a+b+c+d+e+f+g+h
            #print "eType=%s" %(eType)
            
            #dataIn = [vugrid,mfx,mfy,mfxy,a,b,c,bmx,bmy,bmxy,syz,szx,d]
            print "force %s" %(self.ElementType(self.elementType)),dataIn
            #eid = self.obj.addNewEid(out)            
            self.obj.add(nNodes,dt,dataIn)
            #print "len(data) = ",len(self.data)
        ###
        self.handleResultsBuffer(self.OEF_Force_VUTRIA)
        if self.makeOp2Debug:
            print "done with OEF_Force_VUTRIA"
        print self.force_VU_2D

