import sys
from struct import unpack

class ElementsStressStrain(object):
    def CROD_1(self,stress): # works
        self.op2Debug.write('---CROD_1---\n')
        deviceCode = self.deviceCode
        while len(self.data)>=20:
            #self.printSection(40)
            eData     = self.data[0:20]
            self.data = self.data[20: ]
            #print "len(data) = ",len(eData)

            out = unpack('iffff',eData)
            (eid,axial,axialMS,torsion,torsionMS) = out
            self.op2Debug.write('%s\n' %(str(out)))
            eid = (eid - deviceCode) / 10
            stress.addNewEid(eid,axial,axialMS,torsion,torsionMS)

            print "eid=%i axial=%i torsion=%i" %(eid,axial,torsion)
            print "len(data) = ",len(self.data)

        self.skip(4) ## @todo may cause problems later...
        ###
        #print self.rodStress[self.iSubcase]
        print "done with CROD-1"

    def CBEAM_2(self,stress): # not tested; CBEAM class not written
        self.op2Debug.write('---BEAM_2---\n')
        deviceCode = self.deviceCode
        nNodes = 11
        while len(self.data)>=20:
            #self.printSection(40)
            eData     = self.data[0:44]
            self.data = self.data[44: ]
            #print "len(data) = ",len(eData)

            out = struct.unpack('iifffffffff', eData)
            (eid,nodeID,sd,sxc,sxd,sxe,sxf,smax,smin,mst,msc) = out
            self.op2Debug.write('%s\n' %(str(out)))
            eid = (eid - deviceCode) / 10
            stress.addNewEid(eid,axial,axialMS,torsion,torsionMS)
            
            for iNode in range(nNodes):
                eData     = self.data[0:40]
                self.data = self.data[40: ]
                out = struct.unpack('iifffffffff', eData)
                (nodeID,sd,sxc,sxd,sxe,sxf,smax,smin,mst,msc) = out
                self.op2Debug.write('%s\n' %(str(out)))
                stress.add(eid,axial,axialMS,torsion,torsionMS)

            print "eid=%i axial=%i torsion=%i" %(eid,axial,torsion)
            print "len(data) = ",len(self.data)

        if len(data)>4:
            print "*******there may be a problem..."
            print "lenExpected=%s len(data)=%s" %(20,len(self.data))
            self.printBlock(self.data)
            print "*******there may be a problem..."
        self.skip(4) ## @todo may cause problems later...
        ###
        print self.beamStress[self.iSubcase]
        print "done with CBEAM-2"

    def CBAR_34(self,stress):
        self.op2Debug.write('---CBAR_34---\n')
        deviceCode = self.deviceCode
        #print "len(data) = ",len(self.data)
        while len(self.data)>=64:
            #self.printBlock(self.data)
            eData     = self.data[0:4*16]
            self.data = self.data[4*16: ]
            #print "len(data) = ",len(eData)

            (eid,s1a,s2a,s3a,s4a,axial,smaxa,smina,MSt,
                 s1b,s2b,s3b,s4b,      smaxb,sminb,MSc)= unpack('ifffffffffffffff',eData)
            eid = (eid - deviceCode) / 10
            stress.addNewEid('CBAR',eid,s1a,s2a,s3a,s4a,axial,smaxa,smina,MSt,
                                        s1b,s2b,s3b,s4b,      smaxb,sminb,MSc)

            #print "eid=%i s1=%i s2=%i s3=%i s4=%i axial=%-5i smax=%i smin=%i MSt=%i MSc=%i" %(eid,s1a,s2a,s3a,s4a,axial,smaxa,smina,MSt,MSc)
            #print "         s1=%i s2=%i s3=%i s4=%i          smax=%i smin=%i" %(s1b,s2b,s3b,s4b,smaxb,sminb)
            #print "len(data) = ",len(self.data)

        if len(self.data)>4:
            print "*******there may be a problem..."
            print "lenExpected=%s len(data)=%s" %(64,len(self.data))
            self.printBlock(self.data)
            print "*******there may be a problem..."
        #sys.exit('asdf')
        self.skip(4) ## @todo may cause problems later...
        ###
        print "done with CBAR-34"

    def CTRIA3_74(self,stress): # works
        """
        DISTANCE,NORMAL-X,NORMAL-Y,SHEAR-XY,ANGLE,MAJOR,MINOR,VONMISES
        stress is extracted at the centroid
        """
        self.op2Debug.write('---CTRIA3_74---\n')
        deviceCode = self.deviceCode
        #self.printSection(20)
        #term      = self.data[0:4] CEN/
        #self.data = self.data[4:]
        print "*****"
        while self.data:
            eData     = self.data[0:4*17]
            self.data = self.data[4*17: ]
            #print "len(data) = ",len(eData)
            out = unpack('iffffffffffffffff',eData[0:68])

            (eid,fd1,sx1,sy1,txy1,angle1,major1,minor1,vm1,
                 fd2,sx2,sy2,txy2,angle2,major2,minor2,vm2,) = out
            eid = (eid - deviceCode) / 10
            #print "eid=%i fd1=%i sx1=%i sy1=%i txy1=%i angle1=%i major1=%i minor1=%i vm1=%i" %(eid,fd1,sx1,sy1,txy1,angle1,major1,minor1,vm1)
            #print  "      fd2=%i sx2=%i sy2=%i txy2=%i angle2=%i major2=%i minor2=%i vm2=%i\n"   %(fd2,sx2,sy2,txy2,angle2,major2,minor2,vm2)
            stress.addNewEid('CTRIA3',eid,'C',fd1,sx1,sy1,txy1,angle1,major1,minor1,vm1)
            stress.add(               eid,'C',fd2,sx2,sy2,txy2,angle2,major2,minor2,vm2)
            self.op2Debug.write('%s\n' %(str(out)))

            #print "len(data) = ",len(data)

            #sys.exit('asdf')
        self.skip(4) ## @todo may cause problems later...
        ###

    def CTETRA_39(self,stress):  # doesnt work..
        """
        stress is extracted at the centroid
        """
        self.op2Debug.write('---CTETRA_39---\n')
        deviceCode = self.deviceCode
        nNodes = 5 # 1 centroid + 4 corner points
        #self.printSection(20)
        #term      = self.data[0:4] CEN/
        #self.data = self.data[4:]
        print "*****"
        while self.data:
            for nodeID in range(nNodes):   #nodes pts
                if nodeID==0:
                    (eid,_,_,_,_) = unpack("issss",self.data[0:8])
                    self.data = self.data[8:]
                    eid = (eid - deviceCode) / 10


                eData     = self.data[0:4*11]
                self.data = self.data[4*11: ]

                out = unpack('iffffffffff',eData)
                (grid,sxx,txy,s1,vm,syy,txy,s2,szz,txz,s3) = out
                self.op2Debug.write('%s\n' %(str(out)))
                print "eid=%s grid=%s s1=%i s2=%i s3=%s" %(eid,grid,s1,s2,s3)
                sys.exit('CTETRA_39')
                smax = max(s1,s2,s3)
                smin = min(s1,s2,s3)
                stress.add(eid,grid,sxx,syy,txy,s1,s2,s3)  # not fully supported...

                #print "eid=%i grid=%i fd1=%i sx1=%i sy1=%i txy1=%i angle1=%i major1=%i minor1=%i vm1=%i" %(eid,grid,fd1,sx1,sy1,txy1,angle1,major1,minor1,vm1)
                #print "               fd2=%i sx2=%i sy2=%i txy2=%i angle2=%i major2=%i minor2=%i vm2=%i\n"          %(fd2,sx2,sy2,txy2,angle2,major2,minor2,vm2)
                print "len(data) = ",len(self.data)
                #self.printBlock(data)
            ###
            print '--------------------'
            print "len(data) = ",len(self.data)
            
            #self.printSection(100)
            #sys.exit('asdf')
        self.skip(4)
        ###

    def CSOLID_67(self,stress):  # kind of works...
        """
        stress is extracted at the centroid
        """
        self.op2Debug.write('---CSOLID_67---\n')
        print "starting solid element..."
        deviceCode = self.deviceCode
        #nNodes = 5 # 1 centroid + 4 corner points
        #self.printSection(20)
        #term      = self.data[0:4] CEN/
        #self.data = self.data[4:]
        print "*****"

        nNodes=1  # this is a minimum, it will be reset later
        nNodesExpected = 1
        while len(self.data)>= 16+84*nNodesExpected:
            eData     = self.data[0:16]
            self.data = self.data[16:]
            #self.printBlock(eData)

            out = unpack("iissssi",eData)
            (eid,cid,a,b,c,d,nNodes) = out
            self.op2Debug.write('%s\n' %(str(out)))
            #print "abcd = |%s|" %(a+b+c+d)
            #print "eid=%s cid=%s nNodes=%s nNodesExpected=%s" %(eid,cid,nNodes,nNodesExpected)
            
            assert nNodes < 21,self.printBlock(eData)
            eid = (eid - deviceCode) / 10

            if(  nNodes in [4,10]):
                elementType = "CTETRA"
                nNodesExpected = 5
            elif(nNodes in [6,15]):
                elementType = "CPENTA"
                nNodesExpected = 7
            elif(nNodes in [8,20]):
                elementType = "CHEXA"
                nNodesExpected = 9
            else:
                raise Exception('not supported....nNodes=%s' %(nNodes))

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
                    #grid = (gridDevice - deviceCode) / 10
                    grid = gridDevice
                ###

                out = unpack('ffffffffffffffffffff',eData[4:4*21])  ## @warning this could be a bug...
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
                    stress.addNewEid(elementType,cid,eid,grid,sxx,syy,szz,sxy,syz,sxz,aCos,bCos,cCos,pressure,svm)
                else:
                    stress.add(                      eid,grid,sxx,syy,szz,sxy,syz,sxz,aCos,bCos,cCos,pressure,svm)
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
            #sys.exit('asdf')

        if len(self.data)>20:
            print "*******there may be a problem..."
            print "offset+84*nNodesExpected=%s len(data)=%s" %(16+84*nNodesExpected,len(self.data))
            self.printBlock(self.data)
            print "*******there may be a problem..."
        
        self.skip(4)
        ###

    def CQUAD4_95(self,stress): # testing...
        """
        GRID-ID  DISTANCE,NORMAL-X,NORMAL-Y,SHEAR-XY,ANGLE,MAJOR MINOR,VONMISES
        """
        self.op2Debug.write('---CQUAD4_95---\n')
        print "getting a composite element..."
        deviceCode = self.deviceCode
        eType = self.ElementType(self.elementType)

        #self.printSection(20)
        #term = data[0:4] CEN/
        #data = data[4:]
        print "*****"
        #self.printBlock(self.data)
        while len(self.data)>=40: # 2+17*5 = 87 -> 87*4 = 348
            nodeID = 'Centroid'
            eData     = self.data[0:4*11]
            self.data = self.data[4*11: ]
            out = unpack('iifffffffff',eData)
            self.op2Debug.write('%s\n' %(str(out)))
            (eid,iLayer,o1,o2,t12,t1z,t2z,angle,major,minor,ovm) = out
            eid = (eid - deviceCode) / 10  ## @todo adjust with deviceCode...
            stress.addNewEid(eType,eid,o1,o2,t12,t1z,t2z,angle,major,minor,ovm)
            print "eid=%s iLayer=%i o1=%i o2=%i ovm=%i" %(eid,iLayer,o1,o2,ovm)
            
            nextLayer = unpack('i',self.data[0:4])
            #print "nextLayer = ",nextLayer
            #self.printBlock(self.data[:20])

            while len(self.data)>=40:   #nodes pts
                eData     = self.data[0:4*11]
                self.data = self.data[4*11: ]
                out = unpack('iifffffffff',eData)
                
                (eid2,iLayer,o1,o2,t12,t1z,t2z,angle,major,minor,ovm) = out
                self.op2Debug.write('%s\n' %(str(out)))
                eid2 = (eid2 - deviceCode) / 10  ## @todo adjust with deviceCode...
                if eid2!=eid:
                    eid = eid2
                    stress.addNewEid(eType,eid,o1,o2,t12,t1z,t2z,angle,major,minor,ovm)
                else:
                    stress.add(eid,o1,o2,t12,t1z,t2z,angle,major,minor,ovm)
                ###
                print "eid=%s iLayer=%i o1=%i o2=%i ovm=%i" %(eid,iLayer,o1,o2,ovm)

                #eid3,nextLayer = unpack('ii',self.data[0:8])
                    
                #print "nextLayer = ",nextLayer

                #print "len(data) = ",len(self.data)
                #self.printBlock(self.data)
            ###
            #sys.exit('asdf')
            #print '--------------------'
            #print "len(data) = ",len(self.data)
            #print "tell = ",self.op2.tell()
            
            #self.printSection(100)
            #sys.exit('asdf')
            #self.dn += 348
        self.skip(4)  ## @todo may be a problem later on...
        ###

    def CQUAD4_144(self,stress): # works
        """
        GRID-ID  DISTANCE,NORMAL-X,NORMAL-Y,SHEAR-XY,ANGLE,MAJOR MINOR,VONMISES
        """
        self.op2Debug.write('---CQUAD4_144---\n')
        deviceCode = self.deviceCode
        nNodes = 4 # centroid + 4 corner points
        #self.printSection(20)
        #term = data[0:4] CEN/
        #data = data[4:]
        print "*****"
        #self.printBlock(self.data)
        while len(self.data)>=348: # 2+17*5 = 87 -> 87*4 = 348
            nodeID = 'Centroid'
            (eid,_,_,_,_) = unpack("issss",self.data[0:8])
            self.data = self.data[8:]  # 2
            eid = (eid - deviceCode) / 10
            eData     = self.data[0:4*17]
            self.data = self.data[4*17: ]
            out = unpack('iffffffffffffffff',eData[0:68])  # 17
            self.op2Debug.write('%s\n' %(str(out)))
            (grid,fd1,sx1,sy1,txy1,angle1,major1,minor1,vm1,
                  fd2,sx2,sy2,txy2,angle2,major2,minor2,vm2,) = out
            grid = 'C'
            stress.addNewEid('CQUAD4',eid,grid,fd1,sx1,sy1,txy1,angle1,major1,minor1,vm1)
            stress.add(               eid,grid,fd2,sx2,sy2,txy2,angle2,major2,minor2,vm2)

            for nodeID in range(nNodes):   #nodes pts
                eData     = self.data[0:4*17]
                self.data = self.data[4*17: ]
                out = unpack('iffffffffffffffff',eData[0:68])
                self.op2Debug.write('%s\n' %(str(out)))
                (grid,fd1,sx1,sy1,txy1,angle1,major1,minor1,vm1,
                      fd2,sx2,sy2,txy2,angle2,major2,minor2,vm2,) = out

                #print "eid=%i grid=%i fd1=%i sx1=%i sy1=%i txy1=%i angle1=%i major1=%i minor1=%i vm1=%i" %(eid,grid,fd1,sx1,sy1,txy1,angle1,major1,minor1,vm1)
                #print "               fd2=%i sx2=%i sy2=%i txy2=%i angle2=%i major2=%i minor2=%i vm2=%i\n"          %(fd2,sx2,sy2,txy2,angle2,major2,minor2,vm2)
                #print "len(data) = ",len(self.data)
                #self.printBlock(self.data)
                stress.addNewNode(eid,grid,fd1,sx1,sy1,txy1,angle1,major1,minor1,vm1)
                stress.add(       eid,grid,fd2,sx2,sy2,txy2,angle2,major2,minor2,vm2)
            ###
            #print '--------------------'
            #print "len(data) = ",len(self.data)
            #print "tell = ",self.op2.tell()
            
            #self.printSection(100)
            #sys.exit('asdf')
            #self.dn += 348
        self.skip(4)  ## @todo may be a problem later on...
        ###

