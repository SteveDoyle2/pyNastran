import sys
from struct import unpack

class ElementsStressStrain(object):
    def CTRIA3_74(self,data,stress):
        """
        DISTANCE,NORMAL-X,NORMAL-Y,SHEAR-XY,ANGLE,MAJOR,MINOR,VONMISES
        stress is extracted at the centroid
        """
        #self.printSection(20)
        #term = data[0:4] CEN/
        #data = data[4:]
        print "*****"
        while data:
            eData = data[0:4*17]
            data  = data[4*17: ]
            #print "len(data) = ",len(eData)
            out = unpack('iffffffffffffffff',eData[0:68])

            (eid,fd1,sx1,sy1,txy1,angle1,major1,minor1,vm1,
                 fd2,sx2,sy2,txy2,angle2,major2,minor2,vm2,) = out
            eid = (eid - 1) / 10
            stress.addNewEid(eid,grid,sx1,sy1,txy1,major1,minor1,vm1)

            #print "eid=%i fd1=%i sx1=%i sy1=%i txy1=%i angle1=%i major1=%i minor1=%i vm1=%i" %(eid,fd1,sx1,sy1,txy1,angle1,major1,minor1,vm1)
            #print  "      fd2=%i sx2=%i sy2=%i txy2=%i angle2=%i major2=%i minor2=%i vm2=%i\n"   %(fd2,sx2,sy2,txy2,angle2,major2,minor2,vm2)
            #print "len(data) = ",len(data)
            ###
            #sys.exit('asdf')
        self.skip(4)
        ###

    def CTETRA_39(self,data,stress):  # doesnt work..
        """
        ???
        """
        nNodes = 5 # 1 centroid + 4 corner points
        #self.printSection(20)
        #term = data[0:4] CEN/
        #data = data[4:]
        print "*****"
        while data:
            for nodeID in range(nNodes):   #nodes pts
                if nodeID==0:
                    (eid,_,_,_,_) = unpack("issss",data[0:8])
                    data = data[8:]
                    eid = (eid - 1) / 10


                eData = data[0:4*11]
                data  = data[4*11: ]

                (grid,sxx,txy,s1,vm,syy,txy,s2,szz,txz,s3) = unpack('iffffffffff',eData[0:44])
                print "eid=%s grid=%s s1=%i s2=%i s3=%s" %(eid,grid,s1,s2,s3)
                sys.exit('CTETRA_39')
                smax = max(s1,s2,s3)
                smin = min(s1,s2,s3)
                stress.add(eid,grid,sxx,syy,txy,s1,s2,s3)  # not fully supported...

                #print "eid=%i grid=%i fd1=%i sx1=%i sy1=%i txy1=%i angle1=%i major1=%i minor1=%i vm1=%i" %(eid,grid,fd1,sx1,sy1,txy1,angle1,major1,minor1,vm1)
                #print "               fd2=%i sx2=%i sy2=%i txy2=%i angle2=%i major2=%i minor2=%i vm2=%i\n"          %(fd2,sx2,sy2,txy2,angle2,major2,minor2,vm2)
                print "len(data) = ",len(data)
                #self.printBlock(data)
            ###
            print '--------------------'
            print "len(data) = ",len(data)
            
            self.printSection(100)
            #sys.exit('asdf')
        self.skip(4)
        ###

    def CQUAD4_144(self,data,stress):
        """
        GRID-ID  DISTANCE,NORMAL-X,NORMAL-Y,SHEAR-XY,ANGLE,MAJOR MINOR,VONMISES
        """
        nNodes = 4 # centroid + 4 corner points
        #self.printSection(20)
        #term = data[0:4] CEN/
        #data = data[4:]
        print "*****"
        while data:
            nodeID = 'Centroid'
            (eid,_,_,_,_) = unpack("issss",data[0:8])
            data = data[8:]
            eid = (eid - 1) / 10
            eData = data[0:4*17]
            data  = data[4*17: ]
            out = unpack('iffffffffffffffff',eData[0:68])
            (grid,fd1,sx1,sy1,txy1,angle1,major1,minor1,vm1,
                  fd2,sx2,sy2,txy2,angle2,major2,minor2,vm2,) = out
            grid = 'C'
            stress.addNewEid(eid,grid,sx1,sy1,txy1,major1,minor1,vm1)


            for nodeID in range(nNodes):   #nodes pts

                eData = data[0:4*17]
                data  = data[4*17: ]
                out = unpack('iffffffffffffffff',eData[0:68])
                (grid,fd1,sx1,sy1,txy1,angle1,major1,minor1,vm1,
                      fd2,sx2,sy2,txy2,angle2,major2,minor2,vm2,) = out

                stress.add(eid,grid,sx1,sy1,txy1,major1,minor1,vm1)
                #print "eid=%i grid=%i fd1=%i sx1=%i sy1=%i txy1=%i angle1=%i major1=%i minor1=%i vm1=%i" %(eid,grid,fd1,sx1,sy1,txy1,angle1,major1,minor1,vm1)
                #print "               fd2=%i sx2=%i sy2=%i txy2=%i angle2=%i major2=%i minor2=%i vm2=%i\n"          %(fd2,sx2,sy2,txy2,angle2,major2,minor2,vm2)
                print "len(data) = ",len(data)
                #self.printBlock(data)
            ###
            print '--------------------'
            print "len(data) = ",len(data)
            
            #self.printSection(100)
            #sys.exit('asdf')
        self.skip(4)
        ###


