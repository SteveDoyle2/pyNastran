from __future__ import division
import sys
from struct import unpack

from pyNastran.op2.op2_helper import polarToRealImag

class ComplexForces(object):

    def OEF_Rod_alt(self): # 1-CROD, 3-CTUBE, 10-CONROD
        deviceCode = self.deviceCode
        dt = self.nonlinearFactor
        (format1,extract) = self.getOEF_FormatStart()
        format1 += 'ffff'
        isMagnitudePhase = self.isMagnitudePhase()

        while len(self.data)>=20: # 5*4
            eData     = self.data[0:20]
            self.data = self.data[20: ]
            #print "len(data) = ",len(eData)

            out = unpack(format1, eData)
            (eid,axialReal,torqueReal,axialImag,torqueImag) = out

            if isMagnitudePhase:
                (axial)   = polarToRealImag(axialReal,axialImag)
                (torque) = polarToRealImag(torqueReal,torqueImag)
            else:
                axial  = complex(axialReal,axialImag)
                torque = complex(torqueReal,torqueImag)
            eid2  = extract(eid,dt)
            #print "eType=%s" %(eType)
            
            dataIn = [eid2,axial,torque]
            #print "%s" %(self.ElementType(self.elementType)),dataIn
            #eid = self.obj.addNewEid(out)
            self.obj.add(dt,dataIn)
            #print "len(data) = ",len(self.data)
        ###
        self.handleResultsBuffer(self.OEF_Rod_alt)
        #print self.rodForces
        
    def OEF_Beam_alt(self): # 2-CBEAM
        dt = self.nonlinearFactor
        (format1,extract) = self.getOEF_FormatStart()
        isMagnitudePhase = self.isMagnitudePhase()

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

                if isMagnitudePhase:
                    bm1   = polarToRealImag(bm1r,bm1i)
                    bm2   = polarToRealImag(bm2r,bm2i)
                    ts1   = polarToRealImag(ts1r,ts1i)
                    ts2   = polarToRealImag(ts2r,ts2i)
                    af    = polarToRealImag(afr,afi)
                    ttrqr = polarToRealImag(ttrqr,ttrqi)
                    wtrq  = polarToRealImag(wtrqr,wtrqi)
                else:
                    bm1  = complex(bm1r,bm1i)
                    bm2  = complex(bm2r,bm2i)
                    ts1  = complex(ts1r,ts1i)
                    ts2  = complex(ts2r,ts2i)
                    af   = complex(afr,afi)
                    ttrq = complex(ttrqr,ttrqi)
                    wtrq = complex(wtrqr,wtrqi)
                #print "eidTemp = ",eidTemp
                #print "nid = ",nid
                #print "sd = ",sd

                #eid = self.obj.addNewEid(out)
                if i==0: #isNewElement:
                    dataIn = [eid2,nid,sd,bm1,bm2,ts1,ts2,af,ttrq,wtrq]
                    #print "%s cNew   " %(self.ElementType(self.elementType)),dataIn
                    self.obj.addNewElement(dt,dataIn)
                    #print
                elif sd>0.:
                    dataIn = [eid2,nid,sd,bm1,bm2,ts1,ts2,af,ttrq,wtrq]
                    #print "%s cOld   " %(self.ElementType(self.elementType)),dataIn
                    self.obj.add(dt,dataIn)
                    #print
                #else: pass

            #print "len(data) = ",len(self.data)
        ###
        self.handleResultsBuffer(self.OEF_Beam_alt)
        #print self.beamForces

    def OEF_Shear_alt(self): # 4-CSHEAR
        dt = self.nonlinearFactor
        (format1,extract) = self.getOEF_FormatStart()
        format1 += 'ffffffffffffffffffffffffffffffff'
        isMagnitudePhase = self.isMagnitudePhase()

        while len(self.data)>=132: # 33*4
            eData     = self.data[0:132]
            self.data = self.data[132: ]
            #print "len(data) = ",len(eData)

            out = unpack(format1, eData)
            (eid,f41r,f21r,f12r,f32r,f23r,f43r,f34r,f14r,kf1r,s12r,kf2r,s23r,kf3r,s34r,kf4r,s41r,
                 f41i,f21i,f12i,f32i,f23i,f43i,f34i,f14i,kf1i,s12i,kf2i,s23i,kf3i,s34i,kf4i,s41i) = out
            if isMagnitudePhase:
                f41r = polarToRealImag(f41r,f41i); kf1 = polarToRealImag(kf1r,kf1i)
                f21r = polarToRealImag(f21r,f21i); kf2 = polarToRealImag(kf2r,kf2i)
                f12r = polarToRealImag(f12r,f12i); kf3 = polarToRealImag(kf3r,kf3i)
                f23r = polarToRealImag(f23r,f23i); kf4 = polarToRealImag(kf4r,kf4i)
                f32r = polarToRealImag(f32r,f32i); s12 = polarToRealImag(s12r,s12i)
                f43r = polarToRealImag(f43r,f43i); s23 = polarToRealImag(s23r,s23i)
                f34r = polarToRealImag(f34r,f34i); s34 = polarToRealImag(s34r,s34i)
                f14r = polarToRealImag(f14r,f14i); s41 = polarToRealImag(s41r,s41i)
            else:
                f41r = complex(f41r,f41i); kf1 = complex(kf1r,kf1i)
                f21r = complex(f21r,f21i); kf2 = complex(kf2r,kf2i)
                f12r = complex(f12r,f12i); kf3 = complex(kf3r,kf3i)
                f23r = complex(f23r,f23i); kf4 = complex(kf4r,kf4i)
                f32r = complex(f32r,f32i); s12 = complex(s12r,s12i)
                f43r = complex(f43r,f43i); s23 = complex(s23r,s23i)
                f34r = complex(f34r,f34i); s34 = complex(s34r,s34i)
                f14r = complex(f14r,f14i); s41 = complex(s41r,s41i)
            
            eid2  = extract(eid,dt)
            #print "eType=%s" %(eType)
            
            dataIn = [eid2,f41,f21,f12,f32,f23,f43,f34,f14,kf1,s12,kf2,s23,kf3,s34,kf4,s41]
            #print "%s" %(self.ElementType(self.elementType)),dataIn
            #eid = self.obj.addNewEid(out)
            self.obj.add(dt,dataIn)
            #print "len(data) = ",len(self.data)
        ###
        self.handleResultsBuffer(self.OEF_Shear_alt)
        #print self.shearForces
        
    def OEF_Spring_alt(self): # 11-CELAS1, 12-CELAS2, 13-CELAS3, 14-CELAS4
        dt = self.nonlinearFactor
        (format1,extract) = self.getOEF_FormatStart()
        format1 += 'ff'
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
                force = polarToRealImag(forceReal,forceImag)
            else:
                force = complex(forceReal,forceImag)

            dataIn = [eid2,force]
            #print "%s" %(self.ElementType(self.elementType)),dataIn
            #eid = self.obj.addNewEid(out)
            self.obj.add(dt,dataIn)
            #print "len(data) = ",len(self.data)
        ###
        self.handleResultsBuffer(self.OEF_Spring_alt)
        #print self.springForces
        
    def OEF_CVisc_alt(self): # 24-CVISC
        dt = self.nonlinearFactor
        (format1,extract) = self.getOEF_FormatStart()
        format1 += 'ffff'
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
                axial  = polarToRealImag(axialReal,axialImag)
                torque = polarToRealImag(torqueReal,torqueImag)
            else:
                axial  = complex(axialReal,axialImag)
                torque = complex(torqueReal,torqueImag)
            
            dataIn = [eid2,axial,torque]
            #print "%s" %(self.ElementType(self.elementType)),dataIn
            #eid = self.obj.addNewEid(out)
            self.obj.add(dt,dataIn)
            #print "len(data) = ",len(self.data)
        ###
        self.handleResultsBuffer(self.OEF_CVisc_alt)
        #print self.viscForces
        
    def OEF_CBar_alt(self): # 34-CBAR
        dt = self.nonlinearFactor
        (format1,extract) = self.getOEF_FormatStart()
        format1 += 'ffffffffffffffff'
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
                bm1a = polarToRealImag(bm1ar,bm1ai)
                bm2a = polarToRealImag(bm2ar,bm2ai)
                bm1b = polarToRealImag(bm1br,bm1bi)
                bm2b = polarToRealImag(bm2br,bm2bi)
                ts1  = polarToRealImag(ts1r,ts1i)
                ts2  = polarToRealImag(ts2r,ts2i)
                af   = polarToRealImag(afr,afi)
                trq  = polarToRealImag(trqr,trqi)
            else:
                bm1a = complex(bm1ar,bm1ai)
                bm2a = complex(bm2ar,bm2ai)
                bm1b = complex(bm1br,bm1bi)
                bm2b = complex(bm2br,bm2bi)
                ts1  = complex(ts1r,ts1i)
                ts2  = complex(ts2r,ts2i)
                af   = complex(afr,afi)
                trq  = complex(trqr,trqi)

            dataIn = [eid2,bm1a,bm2a,bm1b,bm2b,ts1,ts2,af,trq]
            #print "%s" %(self.ElementType(self.elementType)),dataIn
            #eid = self.obj.addNewEid(out)
            self.obj.add(dt,dataIn)
            #print "len(data) = ",len(self.data)
        ###
        self.handleResultsBuffer(self.OEF_CBar_alt)
        #print self.barForces
        
    def OEF_Plate_alt(self): # 33-CQUAD4,74-CTRIA3
        dt = self.nonlinearFactor
        (format1,extract) = self.getOEF_FormatStart()
        format1 += 'ffffffffffffffff'
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
                mx   = polarToRealImag(mxr,mxi)
                my   = polarToRealImag(myr,myi)
                bmx  = polarToRealImag(bmxr,bmxi)
                bmy  = polarToRealImag(bmyr,bmyi)
                bmxy = polarToRealImag(bmxyr,bmxyi)
                tx   = polarToRealImag(txr,txi)
                ty   = polarToRealImag(tyr,tyi)
            else:
                mx   = complex(mxr,mxi)
                my   = complex(myr,myi)
                bmx  = complex(bmxr,bmxi)
                bmy  = complex(bmyr,bmyi)
                bmxy = complex(bmxyr,bmxyi)
                tx   = complex(txr,txi)
                ty   = complex(tyr,tyi)

            dataIn = [eid2,mx,my,mxy,bmx,bmy,bmxy,tx,ty]
            #print "%s" %(self.ElementType(self.elementType)),dataIn
            #eid = self.obj.addNewEid(out)
            self.obj.add(dt,dataIn)
            #print "len(data) = ",len(self.data)
        ###
        self.handleResultsBuffer(self.OEF_Plate_alt)
        #print self.plateForces

    def OEF_Plate2_alt(self): # 64-CQUAD8,70-CTRIAR,75-CTRIA6,82-CQUAD8,144-CQUAD4-bilinear
        dt = self.nonlinearFactor
        (format1,extract) = self.getOEF_FormatStart()
        format1 += 'cccc'
        isMagnitudePhase = self.isMagnitudePhase()
        
        if self.elementType in [70,75]: # CTRIAR,CTRIA6
            nNodes = 4
        elif self.elementType in [64,82,144]: # CQUAD8,CQUADR,CQUAD4-bilinear
            nNodes = 5
        else:
            raise NotImplementedError(self.codeInformation())
        ###
            
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
                mx   = polarToRealImag(mxr,mxi)
                my   = polarToRealImag(myr,myi)
                mxy  = polarToRealImag(mxyr,mxyi)
                bmx  = polarToRealImag(bmxr,bmxi)
                bmy  = polarToRealImag(bmyr,bmyi)
                bmxy = polarToRealImag(bmxyr,bmxyi)
                tx   = polarToRealImag(txr,txi)
                ty   = polarToRealImag(tyr,tyi)
            else:
                mx   = complex(mxr,mxi)
                my   = complex(myr,myi)
                mxy  = complex(mxyr,mxyi)
                bmx  = complex(bmxr,bmxi)
                bmy  = complex(bmyr,bmyi)
                bmxy = complex(bmxyr,bmxyi)
                tx   = complex(txr,txi)
                ty   = complex(tyr,tyi)

            dataIn = [term,nid,mx,my,mxy,bmx,bmy,bmxy,tx,ty]
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
        dt = self.nonlinearFactor
        (format1,extract) = self.getOEF_FormatStart()
        format1 += 'ifffffffffffffffffffffffff'
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
                bm1A = polarToRealImag(bm1Ar,bm1Ai); bm1B = polarToRealImag(bm1Br,bm1Bi)
                bm2A = polarToRealImag(bm2Ar,bm2Ai); bm2B = polarToRealImag(bm2Br,bm2Bi)
                ts1A = polarToRealImag(ts1Ar,ts1Ai); ts1B = polarToRealImag(ts1Br,ts1Bi)
                afA  = polarToRealImag(afAr, afAi);  afB  = polarToRealImag(afBr, afBi)
                trqA = polarToRealImag(trqAr,trqAi); trqB = polarToRealImag(trqBr,trqBi)
            else:
                bm1A = complex(bm1Ar,bm1Ai); bm1B = complex(bm1Br,bm1Bi)
                bm2A = complex(bm2Ar,bm2Ai); bm2B = complex(bm2Br,bm2Bi)
                ts1A = complex(ts1Ar,ts1Ai); ts1B = complex(ts1Br,ts1Bi)
                afA  = complex(afAr, afAi);  afB  = complex(afBr, afBi)
                trqA = complex(trqAr,trqAi); trqB = complex(trqBr,trqBi)

            dataIn = [eid2,nidA,bm1A,bm2A,ts1A,ts2A,afA,trqA,
                           nidB,bm1B,bm2B,ts1B,ts2B,afB,trqB]
            #print "%s" %(self.ElementType(self.elementType)),dataIn
            #eid = self.obj.addNewEid(out)
            self.obj.add(dt,dataIn)
            #print "len(data) = ",len(self.data)
        ###
        self.handleResultsBuffer(self.OEF_Bend_alt)
        #print self.bendForces
        
    def OEF_PentaPressure_alt(self): # 76-CHEXA_PR,77-CPENTA_PR,78-CTETRA_PR
        dt = self.nonlinearFactor
        (format1,extract) = self.getOEF_FormatStart()
        format1 += 'ccccccccfffffffffffff'
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
                ax = polarToRealImag(axr,axi); vx = polarToRealImag(vxr,vxi)
                ay = polarToRealImag(ayr,ayi); vy = polarToRealImag(vyr,vyi)
                az = polarToRealImag(azr,azi); vz = polarToRealImag(vzr,vzi)
            else:
                ax = complex(axr,axi); vx = complex(vxr,vxi)
                ay = complex(ayr,ayi); vy = complex(vyr,vyi)
                az = complex(azr,azi); vz = complex(vzr,vzi)

            dataIn = [eid2,eName,ax,ay,az,vx,vy,vz,pressure]
            #print "%s" %(self.ElementType(self.elementType)),dataIn
            #eid = self.obj.addNewEid(out)
            self.obj.add(dt,dataIn)
            #print "len(data) = ",len(self.data)
        ###
        self.handleResultsBuffer(self.OEF_PentaPressure_alt)
        #print self.bendForces

    def OEF_CBush_alt(self): # 102-CBUSH
        dt = self.nonlinearFactor
        (format1,extract) = self.getOEF_FormatStart()
        format1 += 'ffffffffffff'
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
            complex
            complex
            complex
            
            if isMagnitudePhase:
                fx = polarToRealImag(fxr,fxi); mx = polarToRealImag(mxr,mxi)
                fy = polarToRealImag(fyr,fyi); my = polarToRealImag(myr,myi)
                fz = polarToRealImag(fzr,fzi); mz = polarToRealImag(mzr,mzi)
            else:
                fx = complex(fxr,fxi); mx = complex(mxr,mxi)
                fy = complex(fyr,fyi); my = complex(myr,myi)
                fz = complex(fzr,fzi); mz = complex(mzr,mzi)

            dataIn = [eid2,fx,fy,fz,mx,my,mz]
            #print "%s" %(self.ElementType(self.elementType)),dataIn
            #eid = self.obj.addNewEid(out)
            self.obj.add(dt,dataIn)
            #print "len(data) = ",len(self.data)
        ###
        self.handleResultsBuffer(self.OEF_CBush_alt)
        #print self.bushForces

    def OEF_Force_VU_alt(self): # 191-VUBEAM
        dt = self.nonlinearFactor
        (format1,extract) = self.getOEF_FormatStart()
        format1 += 'iicccc'
        isMagnitudePhase = self.isMagnitudePhase()

        if self.elementType in [191]:
            nNodes = 2
        else:
            raise NotImplementedError(self.codeInformation())

        formatAll = 'ifffffffffffff'
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
                    forceX     = polarToRealImag(forceXr,forceXi)
                    shearY     = polarToRealImag(shearYr,shearYi)
                    shearZ     = polarToRealImag(shearZr,shearZi)
                    torsionr   = polarToRealImag(torsionr,torsioni)
                    bendingYr  = polarToRealImag(bendingYr,bendingYi)
                    bendingZr  = polarToRealImag(bendingZr,bendingZi)
                else:
                    forceX     = complex(forceXr,forceXi)
                    shearY     = complex(shearYr,shearYi)
                    shearZ     = complex(shearZr,shearZi)
                    torsionr   = complex(torsionr,torsioni)
                    bendingYr  = complex(bendingYr,bendingYi)
                    bendingZr  = complex(bendingZr,bendingZi)

                out2 = [vugrid,posit,forceX,shearY,shearZ,torsion,bendingY,bendingZ]
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
        (format1,extract) = self.getOEF_FormatStart()
        format1 += 'iiccccii'
        isMagnitudePhase = self.isMagnitudePhase()

        if self.elementType in [189]: # VUQUAD
            nNodes = 4
        elif self.elementType in [190]: # VUTRIA
            nNodes = 3
        else:
            raise NotImplementedError(self.codeInformation())

        formatAll = 'ifffiiifffffifffiiifffffi'
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
                    mfx  = polarToRealImag(mfxr,mfxi)
                    mfy  = polarToRealImag(mfyr,mfyi)
                    mfxy = polarToRealImag(mfxyr,mfxyi)
                    bmx  = polarToRealImag(bmxr,bmxi)
                    bmy  = polarToRealImag(bmyr,bmyi)
                    bmxy = polarToRealImag(bmxyr,bmxyi)
                    syz  = polarToRealImag(syzr,syzi)
                    szx  = polarToRealImag(szxr,szxi)
                else:
                    mfx  = complex(mfxr,mfxi)
                    mfy  = complex(mfyr,mfyi)
                    mfxy = complex(mfxyr,mfxyi)
                    bmx  = complex(bmxr,bmxi)
                    bmy  = complex(bmyr,bmyi)
                    bmxy = complex(bmxyr,bmxyi)
                    syz  = complex(syzr,syzi)
                    szx  = complex(szxr,szxi)

                out2 = [vugrid,mfx,mfy,mfxy,bmx,bmy,bmxy,syz,szx]
                forces.append(out2)
            ###
            dataIn.append(forces)
            #print "eType=%s" %(eType)
            
            #dataIn = [vugrid,mfxr,mfyr,mfxyr,bmxr,bmyr,bmxyr,syzr,szxr,
                             #mfxi,mfyi,mfxyi,bmxi,bmyi,bmxyi,syzi,szxi]
            #print "force %s" %(self.ElementType(self.elementType)),dataIn
            #eid = self.obj.addNewEid(out)
            self.obj.add(nNodes,dt,dataIn)
            #print "len(data) = ",len(self.data)
        ###
        self.handleResultsBuffer(self.OEF_Force_VUTRIA_alt)
        if self.makeOp2Debug:
            print "done with OEF_Force_VUTRIA"
        #print self.force_VU_2D
