import sys
from struct import unpack

class ElementsStressStrain(object):

    def handleStressBuffer(self,func,stress):
        """
        works by knowing that:
        the end of an unbuffered table has a
            - [4]
        the end of an table with a buffer has a
            - [4,4,x,4] where x is the next buffer size, which may have another buffer
        the end of the final buffer block has
            - nothing!
        
        The code knows that the large buffer is the default size and the
        only way there will be a smaller buffer is if there are no more
        buffers.  So, the op2 is shifted by 1 word (4 bytes) to account for
        this end shift.  An extra marker value is read, but no big deal.
        Beyond that it's just appending some data to the binary string
        and calling the function that's passed in
        """
        #print stress
        #print "len(data) = ",len(self.data)
        marker = self.readInts(1)
        #print "marker = ",marker
        #if marker[0]==4:
        #    print "found a 4 - end of unbuffered table"

        if len(self.data)>0:
            print "*******len(self.data)=%s...assuming a buffer block" %(len(self.data))
            markers = self.readHeader()
            data = self.readBlock()
            if len(data)<marker:
                self.goto(self.n-4) # handles last buffer not having an extra 4
            self.data += data
            func(stress)
        ###

    def CROD_1(self,stress): # works
        self.op2Debug.write('---CROD_1---\n')
        deviceCode = self.deviceCode
        assert self.numWide==5,'invalid numWide...numWide=%s' %(self.numWide)
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

            #print "eid=%i axial=%i torsion=%i" %(eid,axial,torsion)
            #print "len(data) = ",len(self.data)
        ###
        self.handleStressBuffer(self.CROD_1,stress)
        #print self.rodStress[self.iSubcase]
        print "done with CROD-1"

    def CBEAM_2(self,stress): # not tested; CBEAM class not written
        self.op2Debug.write('---BEAM_2---\n')
        deviceCode = self.deviceCode
        nNodes = 11
        assert self.numWide==12*10+1,'invalid numWide...numWide=%s' %(self.numWide)
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

            #print "eid=%i axial=%i torsion=%i" %(eid,axial,torsion)
            #print "len(data) = ",len(self.data)
        ###
        self.handleStressBuffer(self.CBEAM_2,stress)
        #print self.beamStress[self.iSubcase]
        print "done with CBEAM-2"

    def CQUAD4_33(self,stress): # works
        """
        GRID-ID  DISTANCE,NORMAL-X,NORMAL-Y,SHEAR-XY,ANGLE,MAJOR MINOR,VONMISES
        """
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
            self.op2Debug.write('%s\n' %(str(out)))
            (eid,fd1,sx1,sy1,txy1,angle1,major1,minor1,maxShear1,
                 fd2,sx2,sy2,txy2,angle2,major2,minor2,maxShear2) = out
            eid = (eid - deviceCode) / 10

            #print "eid=%i grid=%s fd1=%-3.1f sx1=%i sy1=%i txy1=%i angle1=%i major1=%i minor1=%i vm1=%i" %(eid,'C',fd1,sx1,sy1,txy1,angle1,major1,minor1,maxShear1)
            #print   "             fd2=%-3.1f sx2=%i sy2=%i txy2=%i angle2=%i major2=%i minor2=%i vm2=%i\n"       %(fd2,sx2,sy2,txy2,angle2,major2,minor2,maxShear2)
            #print "nNodes = ",nNodes
            stress.addNewEid('CQUAD4',eid,'C',fd1,sx1,sy1,txy1,angle1,major1,minor1,maxShear1)
            stress.add(               eid,'C',fd2,sx2,sy2,txy2,angle2,major2,minor2,maxShear2)
            
            #sys.exit('stopCQUAD33')
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
            #self.dn += 348
        ###
        self.handleStressBuffer(self.CQUAD4_33,stress)

    def CBAR_34(self,stress):
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
            eid = (eid - deviceCode) / 10
            stress.addNewEid('CBAR',eid,s1a,s2a,s3a,s4a,axial,smaxa,smina,MSt,
                                        s1b,s2b,s3b,s4b,      smaxb,sminb,MSc)

            #print "eid=%i s1=%i s2=%i s3=%i s4=%i axial=%-5i smax=%i smin=%i MSt=%i MSc=%i" %(eid,s1a,s2a,s3a,s4a,axial,smaxa,smina,MSt,MSc)
            #print "         s1=%i s2=%i s3=%i s4=%i          smax=%i smin=%i" %(s1b,s2b,s3b,s4b,smaxb,sminb)
            #print "len(data) = ",len(self.data)
        ###
        self.handleStressBuffer(self.CBAR_34,stress)
        print "done with CBAR-34"

    def CSOLID_67(self,stress):  # works
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
        #print "*****"

        nNodes=1  # this is a minimum, it will be reset later
        nNodesExpected = 1
        assert self.numWide in [109,151,193],'invalid numWide...numWide=%s' %(self.numWide)
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
        ###
        self.handleStressBuffer(self.CSOLID_67,stress)
        #print self.solidStress[self.iSubcase]

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
        #print "*****"
        assert self.numWide==17,'invalid numWide...numWide=%s' %(self.numWide)
        while len(self.data)>=68:
            eData     = self.data[0:4*17]
            self.data = self.data[4*17: ]
            #print "len(data) = ",len(eData)
            out = unpack('iffffffffffffffff',eData)

            (eid,fd1,sx1,sy1,txy1,angle1,major1,minor1,vm1,
                 fd2,sx2,sy2,txy2,angle2,major2,minor2,vm2,) = out
            eid = (eid - deviceCode) / 10
            #print "eid=%i fd1=%i sx1=%i sy1=%i txy1=%i angle1=%i major1=%i minor1=%i vm1=%i" %(eid,fd1,sx1,sy1,txy1,angle1,major1,minor1,vm1)
            #print  "      fd2=%i sx2=%i sy2=%i txy2=%i angle2=%i major2=%i minor2=%i vm2=%i\n"   %(fd2,sx2,sy2,txy2,angle2,major2,minor2,vm2)
            stress.addNewEid('CTRIA3',eid,'C',fd1,sx1,sy1,txy1,angle1,major1,minor1,vm1)
            stress.add(               eid,'C',fd2,sx2,sy2,txy2,angle2,major2,minor2,vm2)
            self.op2Debug.write('%s\n' %(str(out)))

            #print "len(data) = ",len(data)
        ###
        self.handleStressBuffer(self.CTRIA3_74,stress)

    def CQUAD4_95(self,stress): # works (doesnt handle all stress/strain cases tho
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
        #print "*****"
        #self.printBlock(self.data)
        assert self.numWide==11,'invalid numWide...numWide=%s' %(self.numWide)
        while len(self.data)>=40: # 2+17*5 = 87 -> 87*4 = 348
            eData     = self.data[0:4*11]
            self.data = self.data[4*11: ]
            out = unpack('iifffffffff',eData)
            self.op2Debug.write('%s\n' %(str(out)))
            (eid,iLayer,o1,o2,t12,t1z,t2z,angle,major,minor,ovm) = out
            eid = (eid - deviceCode) / 10  ## @todo adjust with deviceCode...
            stress.addNewEid(eType,eid,o1,o2,t12,t1z,t2z,angle,major,minor,ovm)
            #print "eid=%s iLayer=%i o1=%i o2=%i ovm=%i" %(eid,iLayer,o1,o2,ovm)
            
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
                #print "eid=%s iLayer=%i o1=%i o2=%i ovm=%i" %(eid,iLayer,o1,o2,ovm)

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
        ###
        self.handleStressBuffer(self.CQUAD4_95,stress)

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
        #print "*****"
        #self.printBlock(self.data)
        assert self.numWide==87,'invalid numWide...numWide=%s' %(self.numWide)
        if self.numWide==87: # 2+(17-1)*5 = 87 -> 87*4 = 348
            while len(self.data)>=348: # 2+17*5 = 87 -> 87*4 = 348
                (eid,_,_,_,_) = unpack("issss",self.data[0:8])
                self.data = self.data[8:]  # 2
                eid = (eid - deviceCode) / 10
                eData     = self.data[0:4*17]
                self.data = self.data[4*17: ]
                out = unpack('iffffffffffffffff',eData)  # len=17*4
                self.op2Debug.write('%s\n' %(str(out)))
                (grid,fd1,sx1,sy1,txy1,angle1,major1,minor1,vm1,
                      fd2,sx2,sy2,txy2,angle2,major2,minor2,vm2,) = out
                grid = 'C'
                stress.addNewEid('CQUAD4',eid,grid,fd1,sx1,sy1,txy1,angle1,major1,minor1,vm1)
                stress.add(               eid,grid,fd2,sx2,sy2,txy2,angle2,major2,minor2,vm2)

                for nodeID in range(nNodes):   #nodes pts
                    eData     = self.data[0:4*17]
                    self.data = self.data[4*17: ]
                    out = unpack('iffffffffffffffff',eData)
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
            ###
        ###
        elif self.numWide==77:
            while len(self.data)>=308: # 2+15*5 = 77 -> 77*4 = 308
                (eid,_,_,_,_) = unpack("issss",self.data[0:8])
                self.data = self.data[8:]  # 2
                eid = (eid - deviceCode) / 10
                eData     = self.data[0:4*15]
                self.data = self.data[4*15: ]

                out = unpack('iffffffffffffff',eData)  # 15
                self.op2Debug.write('%s\n' %(str(out)))
                (grid,fd1,sx1r,sx11,sy1r,sy11,txy1r,txy11,
                      fd2,sx2r,sx21,sy2r,sy21,txy2r,txy21) = out

                for nodeID in range(nNodes):   #nodes pts
                    eData     = self.data[0:4*15]
                    self.data = self.data[4*15: ]
                    out = unpack('iffffffffffffff',eData)
                    self.op2Debug.write('%s\n' %(str(out)))
                    (grid,fd1,sx1r,sx11,sy1r,sy11,txy1r,txy11,
                          fd2,sx2r,sx21,sy2r,sy21,txy2r,txy21) = out
                ###
            ###
        ###
        elif self.numWide==47:
            while len(self.data)>=188: # 2+9*5 = 47 -> 47*4 = 188
                (eid,_,_,_,_) = unpack("issss",self.data[0:8])
                self.data = self.data[8:]  # 2
                eid = (eid - deviceCode) / 10
                eData     = self.data[0:4*9]
                self.data = self.data[4*9: ]

                out = unpack('iffffffff',eData)  # 9
                self.op2Debug.write('%s\n' %(str(out)))
                (grid,fd1,sx1,sy1,txy1,
                      fd2,sx2,sy2,txy2) = out

                for nodeID in range(nNodes):   #nodes pts
                    eData     = self.data[0:4*9]
                    self.data = self.data[4*9: ]
                    out = unpack('iffffffff',eData)
                    self.op2Debug.write('%s\n' %(str(out)))
                    (grid,fd1,sx1,sy1,txy1,
                          fd2,sx2,sy2,txy2) = out
                    ###
                ###
            ###
        ###
        else:
            raise Exception('invalid numWide')
        self.handleStressBuffer(self.CQUAD4_144,stress)

