from __future__ import division
import sys
from struct import unpack

from oef_complexForceObjects import *
from pyNastran.op2.op2_helper import polarToRealImag

class ComplexForces(object):

    def OEF_Rod_alt(self): # 1-CROD, 3-CTUBE, 10-CONROD
        deviceCode = self.deviceCode
        dt = self.nonlinearFactor

        (format1,extract) = self.getOEF_FormatStart()
        format1 += 'ffff'
        self.createTransientObject(self.rodForces,ComplexRodForce)
        isMagnitudePhase = self.isMagnitudePhase()

        while len(self.data)>=20: # 5*4
            eData     = self.data[0:20]
            self.data = self.data[20: ]
            #print "len(data) = ",len(eData)

            out = unpack(format1, eData)
            (eid,axialReal,torqueReal,axialImag,torqueImag) = out

            if isMagnitudePhase:
                (axialReal,axialImag)   = polarToRealImag(axialReal,axialImag)
                (torqueReal,torqueImag) = polarToRealImag(torqueReal,torqueImag)
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

        (format1,extract) = self.getOEF_FormatStart()
        self.createTransientObject(self.beamForces,ComplexCBEAMForce)
        #print self.codeInformation()

        isMagnitudePhase = self.isMagnitudePhase()
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

                if isMagnitudePhase:
                    (bm1r,bm1i) = polarToRealImag(bm1r,bm1i)
                    (bm2r,bm2i) = polarToRealImag(bm2r,bm2i)
                    (ts1r,ts1i) = polarToRealImag(ts1r,ts1i)
                    (afr,afi)   = polarToRealImag(afr,afi)
                    (ttrqr,ttrqi) = polarToRealImag(ttrqr,ttrqi)
                    (wtrqr,wtrqi) = polarToRealImag(wtrqr,wtrqi)
                #print "eidTemp = ",eidTemp
                #print "nid = ",nid
                #print "sd = ",sd

                #eid = self.obj.addNewEid(out)
                if i==0: #isNewElement:
                    dataIn = [eid2,nid,sd,bm1r,bm2r,ts1r,ts2r,afr,ttrqr,wtrqr,
                                          bm1i,bm2i,ts1i,ts2i,afi,ttrqi,wtrqi]
                    #print "%s cNew   " %(self.ElementType(self.elementType)),dataIn
                    self.obj.addNewElement(dt,dataIn)
                    #print
                elif sd>0.:
                    dataIn = [eid2,nid,sd,bm1r,bm2r,ts1r,ts2r,afr,ttrqr,wtrqr,
                                          bm1i,bm2i,ts1i,ts2i,afi,ttrqi,wtrqi]
                    #print "%s cOld   " %(self.ElementType(self.elementType)),dataIn
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

        (format1,extract) = self.getOEF_FormatStart()
        format1 += 'ffffffffffffffffffffffffffffffff'
        self.createTransientObject(self.shearForces,ComplexCShearForce)

        isMagnitudePhase = self.isMagnitudePhase()
        while len(self.data)>=132: # 33*4
            eData     = self.data[0:132]
            self.data = self.data[132: ]
            #print "len(data) = ",len(eData)

            out = unpack(format1, eData)
            (eid,f41r,f21r,f12r,f32r,f23r,f43r,f34r,f14r,kf1r,s12r,kf2r,s23r,kf3r,s34r,kf4r,s41r,
                 f41i,f21i,f12i,f32i,f23i,f43i,f34i,f14i,kf1i,s12i,kf2i,s23i,kf3i,s34i,kf4i,s41i) = out
            if isMagnitudePhase:
                (f41r,f41i) = polarToRealImag(f41r,f41i); (kf1r,kf1i) = polarToRealImag(kf1r,kf1i)
                (f21r,f21i) = polarToRealImag(f21r,f21i); (kf2r,kf2i) = polarToRealImag(kf2r,kf2i)
                (f12r,f12i) = polarToRealImag(f12r,f12i); (kf3r,kf3i) = polarToRealImag(kf3r,kf3i)
                (f23r,f23i) = polarToRealImag(f23r,f23i); (kf4r,kf4i) = polarToRealImag(kf4r,kf4i)
                (f32r,f32i) = polarToRealImag(f32r,f32i); (s12r,s12i) = polarToRealImag(s12r,s12i)
                (f43r,f43i) = polarToRealImag(f43r,f43i); (s23r,s23i) = polarToRealImag(s23r,s23i)
                (f34r,f34i) = polarToRealImag(f34r,f34i); (s34r,s34i) = polarToRealImag(s34r,s34i)
                (f14r,f14i) = polarToRealImag(f14r,f14i); (s41r,s41i) = polarToRealImag(s41r,s41i)
                
            eid2  = extract(eid,dt)
            #print "eType=%s" %(eType)
            
            dataIn = [eid2,f41r,f21r,f12r,f32r,f23r,f43r,f34r,f14r,kf1r,s12r,kf2r,s23r,kf3r,s34r,kf4r,s41r,
                           f41i,f21i,f12i,f32i,f23i,f43i,f34i,f14i,kf1i,s12i,kf2i,s23i,kf3i,s34i,kf4i,s41i]
            #print "%s" %(self.ElementType(self.elementType)),dataIn
            #eid = self.obj.addNewEid(out)
            self.obj.add(dt,dataIn)
            #print "len(data) = ",len(self.data)
        ###
        self.handleResultsBuffer(self.OEF_Shear_alt)
        #print self.shearForces
        
    def OEF_Spring_alt(self): # 11-CELAS1, 12-CELAS2, 13-CELAS3, 14-CELAS4
        deviceCode = self.deviceCode
        dt = self.nonlinearFactor

        (format1,extract) = self.getOEF_FormatStart()
        format1 += 'ff'
        self.createTransientObject(self.springForces,ComplexSpringForce)

        isMagnitudePhase = self.isMagnitudePhase()
        while len(self.data)>=12: # 3*4
            eData     = self.data[0:12]
            self.data = self.data[12: ]
            #print "len(data) = ",len(eData)

            out = unpack(format1, eData)
            (eid,forceReal,forceImag) = out
            eid2  = extract(eid,dt)
            #print "eType=%s" %(eType)
            
            if isMagnitudePhase:
                (forceReal,forceImag) = polarToRealImag(forceReal,forceImag)

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

        (format1,extract) = self.getOEF_FormatStart()
        format1 += 'ffff'
        self.createTransientObject(self.viscForces,ComplexViscForce)

        isMagnitudePhase = self.isMagnitudePhase()
        while len(self.data)>=20: # 5*4
            eData     = self.data[0:20]
            self.data = self.data[20: ]
            #print "len(data) = ",len(eData)

            out = unpack(format1, eData)
            (eid,axialReal,torqueReal,axialImag,torqueImag) = out
            eid2  = extract(eid,dt)
            #print "eType=%s" %(eType)

            if isMagnitudePhase:
                (axialReal,axialImag)   = polarToRealImag(axialReal,axialImag)
                (torqueReal,torqueImag) = polarToRealImag(torqueReal,torqueImag)
            
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

        (format1,extract) = self.getOEF_FormatStart()
        format1 += 'ffffffffffffffff'
        self.createTransientObject(self.barForces,ComplexCBARForce)

        isMagnitudePhase = self.isMagnitudePhase()
        while len(self.data)>=68: # 17*4
            eData     = self.data[0:68]
            self.data = self.data[68: ]
            #print "len(data) = ",len(eData)

            out = unpack(format1, eData)
            (eid,bm1ar,bm2ar,bm1br,bm2br,ts1r,ts2r,afr,trqr,
                 bm1ai,bm2ai,bm1bi,bm2bi,ts1i,ts2i,afi,trqi) = out
            eid2  = extract(eid,dt)
            #print "eType=%s" %(eType)
            
            if isMagnitudePhase:
                (bm1ar,bm1ai) = polarToRealImag(bm1ar,bm1ai)
                (bm2ar,bm2ai) = polarToRealImag(bm2ar,bm2ai)
                (bm1br,bm1bi) = polarToRealImag(bm1br,bm1bi)
                (bm2br,bm2bi) = polarToRealImag(bm2br,bm2bi)
                (ts1r,ts1i)   = polarToRealImag(ts1r,ts1i)
                (ts2r,ts2i)   = polarToRealImag(ts2r,ts2i)
                (afr,afi)     = polarToRealImag(afr,afi)
                (trqr,trqi)   = polarToRealImag(trqr,trqi)

            dataIn = [eid2,bm1ar,bm2ar,bm1br,bm2br,ts1r,ts2r,afr,trqr,
                           bm1ai,bm2ai,bm1bi,bm2bi,ts1i,ts2i,afi,trqi]
            #print "%s" %(self.ElementType(self.elementType)),dataIn
            #eid = self.obj.addNewEid(out)
            self.obj.add(dt,dataIn)
            #print "len(data) = ",len(self.data)
        ###
        self.handleResultsBuffer(self.OEF_CBar_alt)
        #print self.barForces
        
    def OEF_Plate_alt(self): # 33-CQUAD4,74-CTRIA3
        deviceCode = self.deviceCode
        dt = self.nonlinearFactor

        (format1,extract) = self.getOEF_FormatStart()
        format1 += 'ffffffffffffffff'
        self.createTransientObject(self.plateForces,ComplexPlateForce)

        isMagnitudePhase = self.isMagnitudePhase()
        while len(self.data)>=68: # 17*4
            eData     = self.data[0:68]
            self.data = self.data[68: ]
            #print "len(data) = ",len(eData)

            out = unpack(format1, eData)
            (eid,mxr,myr,mxyr,bmxr,bmyr,bmxyr,txr,tyr,
                 mxi,myi,mxyi,bmxi,bmyi,bmxyi,txi,tyi) = out
            eid2  = extract(eid,dt)
            #print "eType=%s" %(eType)
            
            if isMagnitudePhase:
                (mxr,mxi)     = polarToRealImag(mxr,mxi)
                (myr,myi)     = polarToRealImag(myr,myi)
                (bmxr,bmxi)   = polarToRealImag(bmxr,bmxi)
                (bmyr,bmyi)   = polarToRealImag(bmyr,bmyi)
                (bmxyr,bmxyi) = polarToRealImag(bmxyr,bmxyi)
                (txr,txi)     = polarToRealImag(txr,txi)
                (tyr,tyi)     = polarToRealImag(tyr,tyi)

            dataIn = [eid2,mxr,myr,mxyr,bmxr,bmyr,bmxyr,txr,tyr,
                           mxi,myi,mxyi,bmxi,bmyi,bmxyi,txi,tyi]
            #print "%s" %(self.ElementType(self.elementType)),dataIn
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
            
        (format1,extract) = self.getOEF_FormatStart()
        format1 += 'cccc'
        self.createTransientObject(self.plateForces2,ComplexPLATE2Force)

        isMagnitudePhase = self.isMagnitudePhase()
        allFormat = 'fffffffffffffffff'
        nTotal = 8+nNodes*68
        while len(self.data)>=nTotal:
            eData     = self.data[0:76]
            self.data = self.data[76: ]
            #print self.printBlock(eData)
            #print "len(data) = ",len(eData)

            out = unpack(format1+allFormat, eData)
            (eid,a,b,c,d,nid,mxr,myr,mxyr,bmxr,bmyr,bmxyr,txr,tyr,
                             mxi,myi,mxyi,bmxi,bmyi,bmxyi,txi,tyi) = out
            term = a+b+c+d # CEN\
            #print "eType=%s" %(eType)

            eid2  = extract(eid,dt)
            
            if isMagnitudePhase:
                (mxr,mxi)     = polarToRealImag(mxr,mxi)
                (myr,myi)     = polarToRealImag(myr,myi)
                (mxyr,mxyi)   = polarToRealImag(mxyr,mxyi)
                (bmxr,bmxi)   = polarToRealImag(bmxr,bmxi)
                (bmyr,bmyi)   = polarToRealImag(bmyr,bmyi)
                (bmxyr,bmxyi) = polarToRealImag(bmxyr,bmxyi)
                (txr,txi)     = polarToRealImag(txr,txi)
                (tyr,tyi)     = polarToRealImag(tyr,tyi)

            dataIn = [term,nid,mxr,myr,mxyr,bmxr,bmyr,bmxyr,txr,tyr,
                               mxi,myi,mxyi,bmxi,bmyi,bmxyi,txi,tyi]
            #print "%s" %(self.ElementType(self.elementType)),dataIn
            self.obj.addNewElement(eid2,dt,dataIn)

            for i in range(nNodes-1):
                eData     = self.data[0:68]
                self.data = self.data[68: ]
                dataIn = unpack(allFormat, eData)
                #(nid,mx,my,mxy,bmx,bmy,bmxy,tx,ty) = out
                #dataIn = [nid,mx,my,mxy,bmx,bmy,bmxy,tx,ty]
                #print "***%s    " %(self.ElementType(self.elementType)),dataIn
                
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

        (format1,extract) = self.getOEF_FormatStart()
        format1 += 'ifffffffffffffffffffffffff'
        self.createTransientObject(self.bendForces,ComplexBendForce)

        isMagnitudePhase = self.isMagnitudePhase()
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
            
            if isMagnitudePhase:
                (bm1Ar,bm1Ai) = polarToRealImag(bm1Ar,bm1Ai); (bm1Br,bm1Bi) = polarToRealImag(bm1Br,bm1Bi)
                (bm2Ar,bm2Ai) = polarToRealImag(bm2Ar,bm2Ai); (bm2Br,bm2Bi) = polarToRealImag(bm2Br,bm2Bi)
                (ts1Ar,ts1Ai) = polarToRealImag(ts1Ar,ts1Ai); (ts1Br,ts1Bi) = polarToRealImag(ts1Br,ts1Bi)
                (afAr, afAi)  = polarToRealImag(afAr, afAi);  (afBr, afBi)  = polarToRealImag(afBr, afBi)
                (trqAr,trqAi) = polarToRealImag(trqAr,trqAi); (trqBr,trqBi) = polarToRealImag(trqBr,trqBi)

            dataIn = [eid2,nidA,bm1Ar,bm2Ar,ts1Ar,ts2Ar,afAr,trqAr,
                                bm1Ai,bm2Ai,ts1Ai,ts2Ai,afAi,trqAi,
                           nidB,bm1Br,bm2Br,ts1Br,ts2Br,afBr,trqBr,
                                bm1Bi,bm2Bi,ts1Bi,ts2Bi,afBi,trqBi]
            #print "%s" %(self.ElementType(self.elementType)),dataIn
            #eid = self.obj.addNewEid(out)
            self.obj.add(dt,dataIn)
            #print "len(data) = ",len(self.data)
        ###
        self.handleResultsBuffer(self.OEF_Bend_alt)
        #print self.bendForces
        
    def OEF_PentaPressure_alt(self): # 76-CHEXA_PR,77-CPENTA_PR,78-CTETRA_PR
        deviceCode = self.deviceCode
        dt = self.nonlinearFactor

        (format1,extract) = self.getOEF_FormatStart()
        format1 += 'ccccccccfffffffffffff'
        self.createTransientObject(self.pentaPressureForces,ComplexPentaPressureForce)

        isMagnitudePhase = self.isMagnitudePhase()
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
            
            if isMagnitudePhase:
                (axr,axi) = polarToRealImag(axr,axi); (vxr,vxi) = polarToRealImag(vxr,vxi)
                (ayr,ayi) = polarToRealImag(ayr,ayi); (vyr,vyi) = polarToRealImag(vyr,vyi)
                (azr,azi) = polarToRealImag(azr,azi); (vzr,vzi) = polarToRealImag(vzr,vzi)

            dataIn = [eid2,eName,axr,ayr,azr,vxr,vyr,vzr,pressure,
                                 axi,ayi,azi,vxi,vyi,vzi]
            #print "%s" %(self.ElementType(self.elementType)),dataIn
            #eid = self.obj.addNewEid(out)
            self.obj.add(dt,dataIn)
            #print "len(data) = ",len(self.data)
        ###
        self.handleResultsBuffer(self.OEF_PentaPressure_alt)
        #print self.bendForces

    def OEF_CBush_alt(self): # 102-CBUSH
        deviceCode = self.deviceCode
        dt = self.nonlinearFactor

        (format1,extract) = self.getOEF_FormatStart()
        format1 += 'ffffffffffff'
        self.createTransientObject(self.bushForces,ComplexCBUSHForce)

        isMagnitudePhase = self.isMagnitudePhase()
        while len(self.data)>=52: # 13*4
            eData     = self.data[0:52]
            self.data = self.data[52: ]
            #print "len(data) = ",len(eData)

            out = unpack(format1, eData)
            (eid,fxr,fyr,fzr,mxr,myr,mzr,
                 fxi,fyi,fzi,mxi,myi,mzi) = out
            eid2  = extract(eid,dt)
            #print "eType=%s" %(eType)
            
            if isMagnitudePhase:
                (fxr,fxi) = polarToRealImag(fxr,fxi); (mxr,mxi) = polarToRealImag(mxr,mxi)
                (fyr,fyi) = polarToRealImag(fyr,fyi); (myr,myi) = polarToRealImag(myr,myi)
                (fzr,fzi) = polarToRealImag(fzr,fzi); (mzr,mzi) = polarToRealImag(mzr,mzi)

            dataIn = [eid2,fxr,fyr,fzr,mxr,myr,mzr,
                           fxi,fyi,fzi,mxi,myi,mzi]
            #print "%s" %(self.ElementType(self.elementType)),dataIn
            #eid = self.obj.addNewEid(out)
            self.obj.add(dt,dataIn)
            #print "len(data) = ",len(self.data)
        ###
        self.handleResultsBuffer(self.OEF_CBush_alt)
        #print self.bushForces

    def OEF_Force_VU_alt(self): # 191-VUBEAM
        dt = self.nonlinearFactor
        #print "numWide = ",self.numWide

        if self.elementType in [191]:
            nNodes = 2
        else:
            raise NotImplementedError(self.codeInformation())

        (format1,extract) = self.getOEF_FormatStart()
        format1 += 'iicccc'
        formatAll = 'ifffffffffffff'
        self.createTransientObject(self.force_VU,ComplexForce_VU)

        isMagnitudePhase = self.isMagnitudePhase()
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
                [vugrid,posit,forceXr,shearYr,shearZr,torsionr,bendingYr,bendingZr,
                              forceXi,shearYi,shearZi,torsioni,bendingYi,bendingZi] = out

                if isMagnitudePhase:
                    (forceXr,forceXi)     = polarToRealImag(forceXr,forceXi)
                    (shearYr,shearYi)     = polarToRealImag(shearYr,shearYi)
                    (shearZr,shearZi)     = polarToRealImag(shearZr,shearZi)
                    (torsionr,torsioni)   = polarToRealImag(torsionr,torsioni)
                    (bendingYr,bendingYi) = polarToRealImag(bendingYr,bendingYi)
                    (bendingZr,bendingZi) = polarToRealImag(bendingZr,bendingZi)

                out2 = [vugrid,posit,forceXr,shearYr,shearZr,torsionr,bendingYr,bendingZr,
                                     forceXi,shearYi,shearZi,torsioni,bendingYi,bendingZi]
                forces.append(out)
            dataIn.append(forces)
            #eType = a+b+c+d+e+f+g+h
            #print "eType=%s" %(eType)
            
            #dataIn = [vugrid,posit,forceX,shearY,shearZ,torsion,bendY,bendZ]
            #print "force %s" %(self.ElementType(self.elementType)),dataIn
            #eid = self.obj.addNewEid(out)
            self.obj.add(nNodes,dt,dataIn)
            #print "len(data) = ",len(self.data)
        ###
        self.handleResultsBuffer(self.OEF_Force_VU_alt)
        if self.makeOp2Debug:
            print "done with OEF_Force_VU"
        #print self.force_VU

    def OEF_Force_VUTRIA_alt(self): # 189-VUQUAD,190-VUTRIA
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
        formatAll = 'ifffiiifffffifffiiifffffi'
        self.createTransientObject(self.force_VU_2D,ComplexForce_VU_2D)

        isMagnitudePhase = self.isMagnitudePhase()
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
                [vugrid,mfxr,mfyr,mfxyr,a,b,c,bmxr,bmyr,bmxyr,syzr,szxr,d,
                        mfxi,mfyi,mfxyi,a,b,c,bmxi,bmyi,bmxyi,syzi,szxi,d] = out

                if isMagnitudePhase:
                    (mfxr,mfxi)   = polarToRealImag(mfxr,mfxi)
                    (mfyr,mfyi)   = polarToRealImag(mfyr,mfyi)
                    (mfxyr,mfxyi) = polarToRealImag(mfxyr,mfxyi)
                    (bmxr,bmxi)   = polarToRealImag(bmxr,bmxi)
                    (bmyr,bmyi)   = polarToRealImag(bmyr,bmyi)
                    (bmxyr,bmxyi) = polarToRealImag(bmxyr,bmxyi)
                    (syzr,syzi)   = polarToRealImag(syzr,syzi)
                    (szxr,szxi)   = polarToRealImag(szxr,szxi)

                out2 = [vugrid,mfxr,mfyr,mfxyr,bmxr,bmyr,bmxyr,syzr,szxr,
                               mfxi,mfyi,mfxyi,bmxi,bmyi,bmxyi,syzi,szxi]
                #if isMagnitudePhase:
                    #(axr,axi) = polarToRealImag(axr,axi)
                forces.append(out2)
            dataIn.append(forces)
            #eType = a+b+c+d+e+f+g+h
            #print "eType=%s" %(eType)
            
            #dataIn = [vugrid,mfxr,mfyr,mfxyr,a,b,c,bmxr,bmyr,bmxyr,syzr,szxr,d,
                             #mfxi,mfyi,mfxyi,a,b,c,bmxi,bmyi,bmxyi,syzi,szxi,d]
            #print "force %s" %(self.ElementType(self.elementType)),dataIn
            #eid = self.obj.addNewEid(out)
            self.obj.add(nNodes,dt,dataIn)
            #print "len(data) = ",len(self.data)
        ###
        self.handleResultsBuffer(self.OEF_Force_VUTRIA_alt)
        if self.makeOp2Debug:
            print "done with OEF_Force_VUTRIA"
        #print self.force_VU_2D
