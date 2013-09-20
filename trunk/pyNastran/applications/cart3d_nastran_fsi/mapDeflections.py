import os
import sys
import copy
from time import time
#from cPickle import loads,dumps

# 3rd party
import numpy
from numpy import array,dot,cross,argsort
from numpy.linalg import det,solve,norm

# Steve's code
from delauney.premorph   import runPremorph
from delauneyReader      import Tet4, DelauneyReader
from mathFunctions       import ListPrint
#from f06Reader import f06Reader

from pyNastran.op2.op2 import OP2
from pyNastran.converters.cart3d.cart3d_reader import generic_cart3d_reader

from pyNastran.utils.log import get_logger
debug = True
log = get_logger(None, 'debug' if debug else 'info')

#------------------------------------------------------------------

class DeflectionReader(object):
    def __init__(self, infilename='fem.op2'):
        log.info('---starting deflectionReader.init of %s---' %(infilename))
        op2 = OP2Reader(infilename)
        terms = ['force','stress','stress_comp','strain','strain_comp','displacement','grid_point_forces']
        op2.read(terms)
        self.displacement = op2.nastranModel.displacement
        #op2.nastranModel.printDisplacement()
        self.convertDisplacements()
        log.info('---finished deflectionReader.init of %s---' %(infilename))

    def getDeflections(self,ID,n0,n1,n2,n3):
        defs = [
        self.getDeflection(n0),
        self.getDeflection(n1),
        self.getDeflection(n2),
        self.getDeflection(n3),
        ]
        #print "defs[%s]-[%s,%s,%s,%s] = %s" %(ID,n0,n1,n2,n3,defs)
        return defs

    def getDeflection(self,gridID):
        if self.deflections.has_key(gridID):
            return self.deflections[gridID]
            #return [float(gridID),]*3 # test
        else:
            return [0.,0.,0.]

    def convertDisplacements(self):
        """
        converts the deflecions from the op2reader to the format used for mapping
        """
        case = '_SUBCASE 1'
        log.info("self.displacement.keys() = %s" %(self.displacement.keys()))
        results = self.displacement[case]
        #sys.exit('stopping...')
        self.deflections = {}
        for gridID,result in results.items():
            #print "gridID = %s" %(gridID)
            layers = result.getLayers()

            for layer in layers:
                #print "layer = ",layer
                #print "  disp.mag = ",layer.mag
                self.deflections[gridID] = array(layer.mag)
        return self.deflections

    def printDisplacement(self):
        log.info("printing Displacement")
        for case,results in self.displacement.items():
            log.info("case = |%s|" %(case))
            for gid,result in results.items():
                log.info("gid = %s" %(gid))
                layers = result.getLayers()
                
                for layer in layers:
                    #print "layer = ",layer
                    log.info("  disp.mag = %s" %(layer.mag))
        log.info("")

#------------------------------------------------------------------

class DeflectionMapper(object):
    def __init__(self,aeroNodes,tets,deflectionReader):  #structuralModel
        self.aeroNodes = aeroNodes
        self.tets = tets
        self.deflectionReader = deflectionReader
        #self.structuralModel = structuralModel
        #self.setAeroInfile()
        #self.setStructuralOutfile()
        #self.structuralOutfile = structuralOutfile

    #def setAeroInfile(self,infile='cart3d.i.tri'):
    #    self.aeroInfile = infile

    #def setStructuralOutfile(self,outfile='fem.f06'):
    #    self.structuralOutfile = outfile

    def buildTetrahedralization(self):
        """runs regtet"""
        pass


    def findClosestTet(self,m,closestTet=None):
        """
        Finds the closest tet.  Has the potential to fail, but until then, skipping.
        Solution to failure:  brute force method.
        m = aeroNode
        tets = 
        """
        if closestTet==None:
            closestTet = tets[1]

        tets = self.tets
        #startingTet = tets[1]
        #closestTet = self.findClosestTet_recursion(m,startingTet,tets)
        #print "found tet = ",closestTet.ID
        #v1 = array([1.,0.,0.])
        
        #tet,tetID = self.bruteForce(m,tets)
        #return tet
        log.info("starting tet = %s" %(closestTet))
        tetID = closestTet.ID
        excluded = []
        log.info("working on point = %s" %(ListPrint(m)))
        counter = 0
        counterMax = 100
        broken = False
        
        isInternal,localVol = closestTet.isInternalNode(m)
        if isInternal:
            log.info("***already Internal")
        while isInternal==False:
            log.info("excluding ID=%s" %(tetID))
            excluded.append(tetID)
            (newTet,minValue) = self.findCloserTet(m,closestTet,tets,excluded)
            closestTet = copy.deepcopy(newTet)
            tetID = closestTet.ID
            #print "excluded = ",len(excluded)
            
            counter +=1
            if(counter==counterMax):
                break

            if tetID in excluded:
                log.info("ERROR***ID=%s was already excluded...excluded=%s" %(tetID,excluded))
                broken = True
                break
            else:
               pass
               #print "ID=%s dist=%s" %(tetID,minValue)

            (isInternal,localVol) = closestTet.isInternalNode(m)
        
        if broken:
            (closestTet,tetID) = self.bruteForce(m,tets,excluded)
        else:
            log.info("*findClosestTet worked!!!")

        #print "*tet[%s]=%s" %(tetID,closestTet)
        #print "nodes = ",closestTet.nodes

        return closestTet,closestTet.ID
        
    def bruteForce(self,m,tets,excluded=[]):
        """
        m is the point
        tets is the list of tets
        excluded is a list of tet ID's to skip
        
        Basically, loop thru the tets and put the point m in each tet.
        The tet thtat is the least bad (ideally perfect) is the one that
        will be used.
        
        If a localVol (localVolume) is negative, then the point is outside.
        However, currently localVol is actually the zeta natural coordinate, so 
        0.<zeta<1.

        If a point isnt found, the least bad point is taken (the one with the lowest
        optValue (optimization value) and a 'test' is performed to check how bad it is.
        """
        log.info("brute Forcing...")
        
        localVols = []
        counter = []
        #print "m = [%g %g %g]" %(m[0],m[1],m[2])
        foundInternalNode = False
        for i,tet in tets.items():
            #if tet.ID in excluded:
            #    pass
            #else:
            (foundInternalNode,localVol) = tet.isInternalNode(m)
            localVols.append(localVol)
            counter.append(i)
            #print "tet[%4s].internal=%s" %(i,foundInternalNode)
            if foundInternalNode:
                foundID = tet.ID
                log.info("*foundID = %s" %(foundID))
                #print self.findCloserTet(m,tets[excluded[-1]],tets)
                #print self.findCloserTet(m,tets[foundID],tets)
                
                return tet,foundID
                #raise Exception('unhandled success!')
                #break
        maxI = localVols.index(max(localVols))
        
        localI = argsort(localVols)
        localVols.sort()
        #for i,localVol in zip(localI,localVols):
        #    log.info("localVol[%s]=%s" %(i,localVol))
        log.info("guessing...closeID = %s" %(tet.ID))
        
        tetOut = tets[maxI]
        log.info('tetOut = %s' %(tetOut))
        (isInternal,optValue) = tetOut.isInternalNode(m)
        log.info("isInternalNode=%s optVal=%g" %(isInternal,optValue))
        return (tetOut,tetOut.ID)


    def findCloserTet(self,m,tet0,tets,excluded=[]):
        """
        Makes an assumption that there is a direct line from the starting tet
        to the final tet that can be minimized with nearly every subsequent tet.
        At worst, it will get the starting node one tet closer to the solution.
        
        m = starting point (mid point)
        """
        #print "findCloserTet"
        #print "working on tet = ",tet0.ID

        cent = tet0.centroid()
        dist = m-cent
        #faces = [tet0.face0,tet0.face1,tet0.face2,tet0.face3]
        dists = [9.e9]*4
        for i,neighbor in enumerate(tet0.neighbors):
            #print "i=%s neighbor=%s centroid=%s" %(i,neighbor,tets[neighbor].centroid())
            if neighbor>0:
                dists[i] = self.distance(m,tets[neighbor].centroid())
        #dists[0] = 9.e9
        #print "dists = ",dists
        
        #print dists
        minValue = min(dists)
        i = dists.index(minValue)
        neighbors = tet0.neighbors
        #print "neighbors = ",neighbors
        #print "tet0.neightbors[%s] = %s" %(i,tet0.neighbors[i])
        tetNew = tets[tet0.neighbors[i]]
        #print "tetNew = ",tetNew
        
        #closestTet = self.findClosestTet_recursion(m,tetNew,tets,excluded)
        return tetNew,minValue
        

    def distance(self,p1,p2):
        return norm(p1-p2)

    def mapDeflections(self, properTets={}):
        """
        Loops thru all the aero nodes, finds the tet it's in, interpolates
        on the deflections at the nodes and maps the deflection to the aero node
        """
        sys.stdout.flush()
        #reader = f06Reader(self.structuralOutfile)
        #d = reader.readDeflections()
        aeroNodes = self.aeroNodes
        #tets = self.tets
        d = self.deflectionReader
        tets = self.tets
        
        #tets = self.buildTetrahedralization()
        
        #aeroNodes = [array([0.0,0.,0.])]
        aeroNodes2 = []
        tet = tets[1]
        log.info("-"*80)
        #print "type(aeroNodes)=%s" %(type(aeroNodes))

        for i,aeroNode in aeroNodes.iteritems():
            if aeroNode[1]<0:
                log.info('skipping aeroNode=%s bc y<0.' % i)
                continue
            log.info("aeroNode[%s]  = %s" % (i, ListPrint(aeroNode)))
            
            #print "aeroNode  = ",aeroNode
            #continue

            if properTets.has_key(i):
                tet = tets[properTets[i]]
            else:
                (tet,ID2) = self.findClosestTet(aeroNode,tet)

            #(isInternal,localVol) = closeTet.isInternalNode(aeroNode)
            #assert isInternal==True
            #print "isInternal?  = ",isInternal

            #print "***tet = %s" %(tet)
            (n0,n1,n2,n3) = tet.nodes
            ID = tet.ID
            deflectionsTet = d.getDeflections(ID,n0,n1,n2,n3)
            aeroNode2 = tet.mapDeflections(deflectionsTet,aeroNode)
            log.info("aeroNode2 = %s" %(ListPrint(aeroNode2)))
            properTets[i] = ID
            aeroNodes2.append(aeroNode2)

            #for tet in tets:  # should select in certain order based on centroids
            #    if tet.isInternalNode(aeroNode):
            #        n0,n1,n2,n3 = tet.nodes
            #        
            #        deflectionsTet = [ d[n0],d[n1],d[n2],d[n3] ]
            #        aeroNode2 = tet.mapDeflections(deflectionsTet,aeroNode)
            #break # uncomment to run one aeroNode
            log.info("-"*80)
            sys.stdout.flush()
        #return aeroNode2
        #for key,value in properTets.items():
        #    print "pointID=%s  -> tetID=%s" %(key,value)
        sys.stdout.flush()
        return (aeroNodes2,properTets)
        
    #def writeAeroInfile(self):
        #self.aeroModel.updateNodes(nodes)
        #self.aeroModel.write(self.aeroFile)

#------------------------------------------------------------------

def loadProperTets(properTetFilename='properTets.in'):
    properTets = {}
    
    if os.path.exists(properTetFilename):
        log.info("loading tets from |%s|..." %(properTetFilename))
        infile = open(properTetFilename,'r')
        lines = infile.readlines()
        for line in lines[1:]:
            #print 
            (pointID,properTet) = line.strip().split()
            properTets[int(pointID)] = int(properTet)
    else:
        log.info("couldnt find tetFile |%s|..." %(properTetFilename))
    return properTets

def writeProperTets(workpath,properTets):
    outfilename = os.path.join(workpath,'properTets.in')
    if not(os.path.exists(outfilename)):
        log.info("writing tets...")
        msg = '#PointID    tetID\n'
        for key,value in properTets.items():
            msg += "%5s  %5s\n" %(key,value)
        outfile = open(outfilename,'wb')
        outfile.write(msg)
        outfile.close()

def test_Tet():
    b  = [10., 0., 0.]
    a  = [ 0.,10., 0.]
    c  = [ 0., 0.,10.]
    d  = [ 0., 0., 0.]
    m1 = [ 1., 1., 1.]
    m2 = [ 2., 2., 2.]
    tet = Tet4(a,b,c,d)
    #print "volume = ",tet.volume()
    log.info("isInternal = \n%s\n" %(tet.isInternalNode(m1)))
    log.info("isInternal = \n%s"   %(tet.isInternalNode(m2)))

def test_deflections():
    infilename = os.path.join('op2reader','solid_shell_bar.op2')
    deflections = {}
    op2 = deflectionReader(infilename)
    #op2.printDisplacement()
    displacements = op2.convertDisplacements()
    
    #for gridID,disp in sorted(displacements.items()):
    #    print "gridID=%s disp=%s" %(gridID,disp)

#------------------------------------------------------------------

#def makeGeometryMorphIn(basepath, bdfModel):
#    args = ['junk',bdfModel,4,4,4,-10,10,-20,20,-30,30]
#    premorphPath = os.path.join(basepath,'delauney','m4_premorph.exe')
    runPremorph(args,bdfModel,premorphPath)


def mapDeflectionsStructures_Aero(bdfModel='test_tet10.bdf',op2Filename='test_tet10.op2',tetFilename='geometry.morph.in',cart3dGeom='bJet.a.tri',cart3dOut='bJet.a.tri2',properTetFilename='properTets.in'):
    t0 = time()
    properTets = loadProperTets(properTetFilename)
    #makeGeometryMorphIn(basepath,bdfModel)


    # loading tetrahedra
    dr = DelauneyReader(tetFilename) # geometry.morph.in
    (tets,nodes,elements,volume) = dr.buildTets()
    sys.stdout.flush()

    #print "type(tets) = ",type(tets)
    #for i,tet in sorted(tets.items()):
    #    print "tet[%4s]=%s" %(i,tet)

    # loading op2 displacements
    defreader = deflectionReader(op2Filename) # test_tet10.op2
    sys.stdout.flush()
    #deflections = defreader.convertDisplacements()
    #deflections = {1:[1.,2.,3.]}
    #for key,d in deflections.items():
    #    print "d = ",d

    # loading aero nodes
    cart = Cart3DReader(cart3dGeom)  # bJet.a.tri
    (cartPoints,elements,regions,loads) = cart.read()
    #(cartPoints,elements,regions,Cp) = cart.makeHalfModel(cartPoints,elements,regions,Cp)
    sys.stdout.flush()


    #cartOutfile = os.path.join(workpath,'bJet.a.tri_test')   # test code
    #cart.writeInfile(cartOutfile,cartPoints,elements,regions)
    #for point in cartPoints:
    #    print "point = ",point


    # deflect the aero nodes
    dmap = DeflectionMapper(cartPoints,tets,defreader)
    t1 = time()
    log.info("setup time = %g sec" %(t1-t0))
    
    (aeroNodes2,properTets) = dmap.mapDeflections(properTets)
    writeProperTets(workpath,properTets)


    # write out the deflected aero nodes
    cart.writeInfile(cart3dOut,aeroNodes2,elements,regions) #bJet.a.tri_new
    log.info("done with deflection mapping!")

    #for aeroNode in aeroNodes2:
    #    print "aeroNode = ",aeroNode
    t2 = time()
    log.info("total mapDeflections.py time = %g sec" %(t2-t0))

#------------------------------------------------------------------

if __name__=='__main__':
    basepath = os.getcwd()
    configpath = os.path.join(basepath,'inputs')
    workpath   = os.path.join(basepath,'outputs')
    bdfModel   = os.path.join(configpath,'fem3.bdf')
    assert os.path.exists(bdfModel),'|%s| doesnt exist' %(bdfModel)

    os.chdir(workpath)
    log.info("basepath = %s" %(basepath))
    tetFilename = os.path.join(configpath,'geometry.morph.in')
    op2Filename = os.path.join(configpath,'fem3.op2')
    cart3dGeom  = os.path.join(configpath,'Cart3d_bwb.i.tri')
    cart3dOut   = os.path.join(workpath,  'Cart3d_bwb.i.tri2')
    properTetFilename = os.path.join(configpath,'properTets.in')  # not required to exist...

    mapDeflectionsStructures_Aero(bdfModel,op2Filename,tetFilename,cart3dGeom,cart3dOut,properTetFilename)

    #mapDeflectionsStructures_Aero()
    #test_Tet()
    #test_deflections()
    
    #model = 'test_tet10' # bdf, op2, triq, tri

    sys.exit('finished mapDeflections.py')
    #iteration = [1,2,3]
    defMapper = DeflectionMapper(aeroModel,structuralModel)
    for i in iteration:
        defMapper.setStructuralOutfile('fem.f06')
        aeroFile = 'cart3d_%s.tri' %(i)
        defMapper.setAeroInfile(aeroFile)
        defMapper.buildTetrahedralization()
        #defMapper.buildMappingMatrix()
        defMapper.mapDeflections()
        #defMapper.writeAeroInfile()
