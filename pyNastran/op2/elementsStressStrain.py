import sys
from struct import unpack

class ElementsStressStrain(object):
    def CROD_1(self,stress): # works
        while len(self.data)>=20:
            self.printSection(40)
            eData     = self.data[0:20]
            self.data = self.data[20: ]
            #print "len(data) = ",len(eData)

            (eid,axial,axialMS,torsion,torsionMS) = unpack('iffff',eData)
            eid = (eid - 1) / 10
            stress.addNewEid(eid,axial,axialMS,torsion,torsionMS)

            print "eid=%i axial=%i torsion=%i" %(eid,axial,torsion)
            print "len(data) = ",len(self.data)

            #sys.exit('asdf')
        self.skip(4) ## @todo may cause problems later...
        ###
        print self.rodStress[self.iSubcase]
        print "done with CROD-1"

    def CBAR_34(self,stress):
        #print "len(data) = ",len(self.data)
        while len(self.data)>=64:
            self.printBlock(self.data)
            eData     = self.data[0:4*16]
            self.data = self.data[4*16: ]
            #print "len(data) = ",len(eData)

            (eid,s1a,s2a,s3a,s4a,axial,smaxa,smina,MSt,
                 s1b,s2b,s3b,s4b,      smaxb,sminb,MSc)= unpack('ifffffffffffffff',eData)
            eid = (eid - 1) / 10
            stress.addNewEid('CBAR',eid,s1a,s2a,s3a,s4a,axial,smaxa,smina,MSt,
                                        s1b,s2b,s3b,s4b,      smaxb,sminb,MSc)

            #print "eid=%i s1=%i s2=%i s3=%i s4=%i axial=%-5i smax=%i smin=%i MSt=%i MSc=%i" %(eid,s1a,s2a,s3a,s4a,axial,smaxa,smina,MSt,MSc)
            #print "         s1=%i s2=%i s3=%i s4=%i          smax=%i smin=%i" %(s1b,s2b,s3b,s4b,smaxb,sminb)
            #print "len(data) = ",len(self.data)

        if len(self.data)>4:
            print "there may be a problem len(self.data)=%s" %(len(self.data))
            self.printBlock(self.data)
        #sys.exit('asdf')
        self.skip(4) ## @todo may cause problems later...
        ###
        print "done with CBAR-34"

    def CTRIA3_74(self,stress): # works
        """
        DISTANCE,NORMAL-X,NORMAL-Y,SHEAR-XY,ANGLE,MAJOR,MINOR,VONMISES
        stress is extracted at the centroid
        """
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
            eid = (eid - 1) / 10
            stress.addNewEid('CTRIA3',eid,'C',fd1,sx1,sy1,txy1,angle1,major1,minor1,vm1)
            stress.add(               eid,'C',fd2,sx2,sy2,txy2,angle2,major2,minor2,vm2)

            #print "eid=%i fd1=%i sx1=%i sy1=%i txy1=%i angle1=%i major1=%i minor1=%i vm1=%i" %(eid,fd1,sx1,sy1,txy1,angle1,major1,minor1,vm1)
            #print  "      fd2=%i sx2=%i sy2=%i txy2=%i angle2=%i major2=%i minor2=%i vm2=%i\n"   %(fd2,sx2,sy2,txy2,angle2,major2,minor2,vm2)
            #print "len(data) = ",len(data)

            #sys.exit('asdf')
        self.skip(4) ## @todo may cause problems later...
        ###

    def CTETRA_39(self,stress):  # doesnt work..
        """
        stress is extracted at the centroid
        """
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
                    eid = (eid - 1) / 10


                eData     = self.data[0:4*11]
                self.data = self.data[4*11: ]

                (grid,sxx,txy,s1,vm,syy,txy,s2,szz,txz,s3) = unpack('iffffffffff',eData)
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
            
            self.printSection(100)
            #sys.exit('asdf')
        self.skip(4)
        ###

    def CHEXA_67(self,stress):  # kind of works...
        """
        stress is extracted at the centroid
        """
        print "starting solid element..."
        #nNodes = 5 # 1 centroid + 4 corner points
        #self.printSection(20)
        #term      = self.data[0:4] CEN/
        #self.data = self.data[4:]
        print "*****"

        nNodes=1  # this is a minimum, it will be reset later
        while len(self.data)>=16+84*nNodes+170:
            #self.printBlock(self.data)
           
            (eid,cid,a,b,c,d,nNodes) = unpack("iissssi",self.data[0:16])
            print "abcd = |%s|" %(a+b+c+d)
            print "eid=%s cid=%s nNodes=%s" %(eid,cid,nNodes)
            assert nNodes < 21
            assert cid >= 0
            self.data = self.data[16:]
            eid = (eid - 1) / 10
            assert eid >= 0

            if(  nNodes==4):  elementType = "CTETRA"
            elif(nNodes==8):  elementType = "CHEXA"
            elif(nNodes==6):  elementType = "CPENTA"
            else:
                raise Exception('not supported....')

            print "len(data) = ",len(self.data)
            for nodeID in range(nNodes):   #nodes pts, +1 for centroid (???)
                print "len(data)A = ",len(self.data)
                eData     = self.data[0:4*21]
                self.data = self.data[4*21: ]
                print "len(data)B = ",len(self.data)

                #print "self.tableCode = ",self.tableCode
                
                #print "len(data) = ",len(self.data)
                out = unpack('iffffffffffffffffffff',self.data[0:4*21])  ## @warning this could be a bug...
                (grid,sxx,sxy,s1,a1,a2,a3,pressure,svm,syy,syz,s2,b1,b2,b3,szz,sxz,s3,c1,c2,c3) = out

                #out = unpack('iffffffffff',self.data[0:4*11])
                #(grid,oxx,txy,o1,ovm,oyy,tyz,o2,ozz,txz,o3) = out

                print "%s eid=%s cid=%s nodef=%s grid=%s s1=%i s2=%i s3=%i svm=%i" %(elementType,eid,cid,nNodes,grid,s1,s2,s3,svm)
                smax = max(s1,s2,s3)
                smin = min(s1,s2,s3)
                
                aCos = []
                bCos = []
                cCos = []
                if nodeID==0:
                    print "adding new eid"
                    stress.addNewEid(elementType,cid,eid,grid,sxx,syy,szz,sxy,syz,sxz,aCos,bCos,cCos,pressure,svm)
                else:
                    stress.add(                      eid,grid,sxx,syy,szz,sxy,syz,sxz,aCos,bCos,cCos,pressure,svm)
                #print "eid=%i grid=%i fd1=%i sx1=%i sy1=%i txy1=%i angle1=%i major1=%i minor1=%i vm1=%i" %(eid,grid,fd1,sx1,sy1,txy1,angle1,major1,minor1,vm1)
                #print "               fd2=%i sx2=%i sy2=%i txy2=%i angle2=%i major2=%i minor2=%i vm2=%i\n"          %(fd2,sx2,sy2,txy2,angle2,major2,minor2,vm2)
                #self.printBlock(data)
            #sys.exit('finished a CEHXA')
            print self.solidStress[self.iSubcase]
            ###
            print '--------------------'
            print "len(data) = ",len(self.data)
            print "84*nNodes=",84*nNodes,nNodes,len(self.data)>=84*nNodes
            
            #self.printSection(100)
            self.printBlock(self.data[0:100])
            #self.printBlock(self.data[1:100])
            #self.printBlock(self.data[2:100])
            #self.printBlock(self.data[3:100])
            #sys.exit('asdf')

        if len(self.data)>20:
            print "*******there may be a problem..."
            self.printBlock(self.data)
        self.skip(4)
        ###

    def CQUAD4_144(self,stress): # works
        """
        GRID-ID  DISTANCE,NORMAL-X,NORMAL-Y,SHEAR-XY,ANGLE,MAJOR MINOR,VONMISES
        """
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
            eid = (eid - 1) / 10
            eData     = self.data[0:4*17]
            self.data = self.data[4*17: ]
            out = unpack('iffffffffffffffff',eData[0:68])  # 17
            (grid,fd1,sx1,sy1,txy1,angle1,major1,minor1,vm1,
                  fd2,sx2,sy2,txy2,angle2,major2,minor2,vm2,) = out
            grid = 'C'
            stress.addNewEid('CQUAD4',eid,grid,fd1,sx1,sy1,txy1,angle1,major1,minor1,vm1)
            stress.add(               eid,grid,fd2,sx2,sy2,txy2,angle2,major2,minor2,vm2)


            for nodeID in range(nNodes):   #nodes pts
                eData     = self.data[0:4*17]
                self.data = self.data[4*17: ]
                out = unpack('iffffffffffffffff',eData[0:68])
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


