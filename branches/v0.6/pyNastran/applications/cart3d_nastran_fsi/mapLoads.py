# standard python modules
import os
import sys
import cPickle
from time import time

from numpy import argsort, mean, array, cross

# my code
from mathFunctions import piercePlaneVector, shepardWeight, Normal, ListPrint # used to be *
#from mathFunctions import getTriangleWeights
from models import StructuralModel, AeroModel
from pyNastran.converters.cart3d.cart3d_reader import Cart3DReader
from pyNastran.bdf.bdf import BDF
from kdtree import KdTree

from pyNastran.utils.log import get_logger
debug = True
log = get_logger(None, 'debug' if debug else 'info')

#------------------------------------------------------------------

class LoadMapping(object):
    def __init__(self, aeroModel, structuralModel, configpath='', workpath=''):
        self.nCloseElements = 30

        self.configpath = configpath
        self.workpath = workpath
        self.aeroModel = aeroModel
        self.structuralModel = structuralModel

        self.mappingMatrix = {}
        self.mapfile = 'mapfile.in'
        self.bdffile = 'fem.bdf'

        self.Sref = 1582876. # inches
        self.Lref = 623.  # inch
        self.xref = 268.

    def setOutput(self,bdffile='fem.bdf', loadCase=1):
        self.bdffile = bdffile
        self.loadCase = loadCase

    def set_flight_condition(self, pInf=499.3, qInf=237.885):
        self.pInf = pInf
        self.qInf = qInf  #1.4/2.*pInf*Mach**2.

    def loadMappingMatrix(self):
        infile = open(self.mapfile,'r')
        self.mappingMatrix = cPickle.loads(infile)
        infile.close()

    def saveMappingMatrix(self):
        outString = cPickle.dumps(self.mappingMatrix)
        outfile = open(self.mapfile, 'wb')
        outfile.write(outString)
        outfile.close()

    def getPressure(self, Cp):
        """
        Caculates the pressure based on:
           Cp = (p-pInf)/qInf
        """
        p = Cp * self.qInf + self.pInf
        #p = Cp*self.qInf # TODO:  for surface with atmospheric pressure inside
        return p

    def mapLoads(self):
        """
        Loops thru the unitLoad mappingMatrix and multiplies by the
        total force F on each panel to calculate the PLOAD that will be
        printed to the output file.
        Also, performs a sum of Forces & Moments to verify the aero loads.
        """
        log.info("---starting mapLoads---")
        self.bdf = open(self.bdffile, 'wb')
        #self.buildMappingMatrix()
        log.info("self.loadCase = %s" % self.loadCase)
        self.loadCases = {self.loadCase:{}}

        #self.loadCases = {self.loadCase={}, }
        momentCenter = array([self.xref,0.,0.])
        sumMoments = array([0.,0.,0.])
        sumForces  = array([0.,0.,0.])
        sys.stdout.flush()
        for aEID, distribution in self.mappingMatrix.iteritems():
            #print "aEID = ",aEID
            #print "***distribution = ",distribution
            sumLoad = 0.
            area   = self.aeroModel.Area(aEID)
            normal = self.aeroModel.Normal(aEID)
            Cp     = self.aeroModel.Cp(aEID)
            #print "Cp = ",Cp
            #print "area[%s]=%s" %(aEID,area)

            p = self.getPressure(Cp)
            centroid = self.aeroModel.Centroid(aEID)
            r = momentCenter - centroid
            F = area * p
            Fn = F * normal
            sumMoments += cross(r, Fn)
            sumForces  += Fn
            for sNID, percentLoad in sorted(distribution.iteritems()):
                sumLoad += percentLoad

                Fxyz = Fn * percentLoad  # negative sign is to be consistent with nastran
                self.addForce(sNID, Fxyz)

                #print "Fxyz = ",Fxyz
                #print "type(structuralModel) = ", type(self.structuralModel)

                #comment = 'percentLoad=%.2f' % percentLoad
                #self.structuralModel.writeLoad(self.bdf, self.loadCase, sNID,
                #                               Fxyz[0], Fxyz[1], Fxyz[2], comment)

            #msg = '$ End of aEID=%s sumLoad=%s p=%s area=%s F=%s normal=%s\n' % (aEID, sumLoad, p, area, F, normal)
            #self.bdf.write(msg)

        self.writeLoads()  # short version of writing loads...
        self.bdf.close()

        log.info("pInf=%g [psi]; qInf= %g [psi]" % (self.pInf, self.qInf))
        log.info("sumForcesFEM  [lb]    = %s" % ListPrint(sumForces))
        log.info("sumMomentsFEM [lb-ft] = %s" % ListPrint(sumMoments/12.))  # divided by 12 to have moments in lb-ft

        Cf = sumForces /(self.Sref * self.qInf)
        log.info("Cf = %s" % ListPrint(Cf))

        Cm = sumMoments / (self.Sref * self.qInf * self.Lref)
        log.info("Cm = %s" % ListPrint(Cm*12.)) # multiply by 12 to nondimensionalize ???  maybe 144...

        #self.bdf.write('$***********\n')
        log.info("wrote loads to %r" % self.bdffile)
        log.info("---finished mapLoads---")

    def writeLoads(self):
        """writes the load in BDF format"""
        log.info("---starting writeLoads---")
        self.bdf.write('$ ***writeLoads***\n')
        self.bdf.write('$ nCloseElements=%s\n' % self.nCloseElements)
        for loadCase, loads in sorted(self.loadCases.iteritems()):
            log.info("  loadCase = %s" % loadCase)
            for (sNID, Fxyz) in sorted(loads.iteritems()):
                self.structuralModel.writeLoad(self.bdf, loadCase, sNID, *Fxyz)

        log.info("finished writeLoads---")

    def addForce(self, sNID, Fxyz):
        try:
            self.loadCases[self.loadCase][sNID]+=Fxyz
        except KeyError:
            self.loadCases[self.loadCase][sNID] =Fxyz

    def buildCentroids(self, model, eids=None):
        centroids = {}

        if eids==None:
            eids = model.ElementIDs()
        for eid in eids:
            #etype = model.etype(eid)
            #if etype=='CTRIA3' or etype=='CQUAD4':
            centroid = model.Centroid(eid)
            #print "centroid[%s] = %s" %(eid, centroid)
            if centroid is not None:
                centroids[eid] = centroid
            #elif len(centroid)!=0:
            #    print "centroid[%s] = %s" %(eid, centroid)
            #    raise Exception('error...bad centroid length')
            #else:
            #    print "skipping %s" %(etype)
        return centroids

    def buildNodalTree(self, sNodes):
        log.info("---start buildNodalTree---")
        raise Exception('DEPRECATED...buildNodalTree in mapLoads.py')
        sys.stdout.flush()
        #print "type(aCentroids)=%s type(sCentroids)=%s" %(type(aCentroids), type(sCentroids))
        self.nodalTree = KdTree('node', sNodes, nClose=self.nCloseNodes)
        log.info("---finish buildNodalTree---")
        sys.stdout.flush()

    def buildCentroidTree(self, sCentroids):
        """
        sCentroids - dict of structural centroids
        id:  element id
        """
        log.info("---start buildCentroidTree---")
        sys.stdout.flush()
        #print "type(aCentroids)=%s type(sCentroids)=%s" %(type(aCentroids), type(sCentroids))

        msg = 'Element '
        for (id,sCentroid) in sorted(sCentroids.items()):
            msg += "%s " %(id)
        log.info(msg + '\n')

        self.centroidTree = KdTree('element',sCentroids,nClose=self.nCloseElements)
        log.info("---finish buildCentroidTree---")
        sys.stdout.flush()

    def parseMapFile(self, mapFilename='mappingMatrix.new.out'):
        """
        This is used for rerunning an analysis quickly (cuts out building the mapping matrix ~1.5 hours).
        1 {8131: 0.046604568185355716, etc...}
        """
        log.info("---starting parseMapFile---")
        mappingMatrix = {}

        log.info('loading mapFile=%r' % mapFilename)
        mapFile = open(mapFilename,'r')
        lines = mapFile.readlines()
        mapFile.close()

        for (i, line) in enumerate(lines[1:]):  # dont read the first line, thats a header line
            line = line.strip()
            #print "line = %r" % line
            (aEID, dictLine) = line.split('{') # splits the dictionary from the aEID
            aEID = int(aEID)
            #assert i == int(aEID)

            # time to parse a dictionary with no leading brace
            distribution = {}
            mapSline = dictLine.strip('{}, ').split(',')
            for pair in mapSline:
                (sEID, weight) = pair.strip(' ').split(':')
                sEID = int(sEID)
                weight = float(weight)
                distribution[sEID] = weight
            mappingMatrix[aEID] = distribution
        #log.info("mappingKeys = %s" %(sorted(mappingMatrix.keys())))
        self.runMapTest(mappingMatrix)
        log.info("---finished parseMapFile---")
        return mappingMatrix

    def runMapTest(self, mappingMatrix, mapTestFilename='mapTest.out'):
        """
        Checks to see what elements loads were mapped to.
        Ideally, all elements would have a value of 1.0, given equal area.
        Instead, each element should have a value of area[i].
        """
        mapTest = {}
        for (aEID, distribution) in sorted(mappingMatrix.iteritems()):
            for (sEID,weight) in distribution.items():
                if sEID in mapTest:
                    mapTest[sEID] += weight
                else:
                    mapTest[sEID] = weight

        mapOut = open(mapTestFilename, 'wb')
        mapOut.write('# sEID  weight\n')
        for sEID, weight in sorted(mapTest.iteritems()):
            mapOut.write("%s %s\n" % (sEID, weight))
        mapOut.close()

    def buildMappingMatrix(self, debug=False):
        """
        Skips building the matrix if it already exists
        A mapping matrix translates element ID to loads on the nearby
        strucutral nodes.

        eid,distribution
        """
        if self.mappingMatrix != {}:
            return self.mappingMatrix

        log.info("---starting buildMappingMatrix---")
        #print "self.mappingMatrix = ",self.mappingMatrix
        if os.path.exists('mappingMatrix.new.out'):
            self.mappingMatrix = self.parseMapFile('mappingMatrix.new.out')
            log.info("---finished buildMappingMatrix based on mappingMatrix.new.out---")
            sys.stdout.flush()
            return self.mappingMatrix
        log.info("...couldn't find 'mappingMatrix.new.out' in %r, so going to make it..." % os.getcwd())

        # this is the else...
        log.info("creating...")
        aModel = self.aeroModel
        sModel = self.structuralModel

        #aNodes = aModel.getNodes()
        #sNodes = sModel.getNodes()
        #treeObj = Tree(nClose=5)
        #tree    = treeObj.buildTree(aNodes,sNodes) # fromNodes,toNodes

        aElementIDs  = aModel.ElementIDs() # list
        sElementIDs  = sModel.getElementIDsWithPIDs() # list
        sElementIDs2 = sModel.ElementIDs() # list

        msg = "there are no internal elements in the structural model?\n   ...len(sElementIDs)=%s len(sElementIDs2)=%s" % (len(sElementIDs), len(sElementIDs2))
        assert sElementIDs != sElementIDs2, msg
        log.info("maxAeroID=%s maxStructuralID=%s sElements=%s" % (max(aElementIDs), max(sElementIDs),len(sElementIDs2)))

        sCentroids = self.buildCentroids(sModel, sElementIDs)
        #print "sCentroids = ",sCentroids

        self.buildCentroidTree(sCentroids)
        #self.buildNodalTree(sNodes)

        aCentroids = self.buildCentroids(aModel)
        #print "aCentroids = ",aCentroids

        mapFile = open('mappingMatrix.out', 'wb')
        mapFile.write('# aEID distribution (sEID:  weight)\n')

        t0 = time()
        log.info("---start piercing---")
        nAeroElements = float(len(aElementIDs))
        log.info("nAeroElements = %s" % nAeroElements)
        for (i, aEID) in enumerate(aElementIDs):
            if i % 1000 == 0:
                log.info('  piercing %sth element' % i)
                sys.stdout.flush()

            aElement = aModel.Element(aEID)
            (aArea, aCentroid, aNormal) = aModel.get_element_properties(aEID)
            percentDone = i / nAeroElements * 100
            log.info('aEID=%s percentDone=%.2f aElement=%s aArea=%s aCentroid=%s aNormal=%s' %(aEID,percentDone,aElement,aArea,aCentroid,aNormal))
            pSource = aCentroid
            (distribution) = self.pierceElements(aCentroid, aEID, pSource, aNormal)
            #(distribution)  = self.poorMansMapping(aCentroid, aEID, pSource, aNormal)
            self.mappingMatrix[aEID] = distribution
            mapFile.write('%s %s\n' % (aEID, distribution))
            #break

            dt = (time() - t0) / 60.
            tEst = dt * nAeroElements / (i + 1)  # dtPerElement*nElements
            tLeft = tEst - dt
            log.debug("tEst=%g minutes; tLeft=%g minutes; %.3f%% done" %(tEst, tLeft, dt/tEst*100.))

        mapFile.close()
        log.info("---finish piercing---")
        self.runMapTest(self.mappingMatrix)
        #print "mappingMatrix = ", self.mappingMatrix
        log.info("---finished buildMappingMatrix---")
        sys.stdout.flush()
        return self.mappingMatrix

    def poorMansMapping(self, aCentroid, aEID, pSource, normal):
        """
        distributes load without piercing elements
        based on distance
        """
        (sElements, sDists) = self.centroidTree.getCloseElementIDs(aCentroid)
        log.info("aCentroid = %s" % aCentroid)
        log.info("sElements = %s" % sElements)
        log.info("sDists    = %s" % ListPrint(sDists))

        setNodes = set([])
        sModel = self.structuralModel
        for sEID in sElements:
            sNodes = sModel.get_element_nodes(sEID)
            setNodes.union(set(sNodes))

        nIDs = list(setNodes)
        sNodes = sModel.getNodeIDLocations(nIDs)
        weights = self.getWeights(closePoint, sNodes)
        distribution = self.createDistribution(nIDs, weights)
        return (distribution)


    def pierceElements(self, aCentroid, aEID, pSource, normal):
        """
        Pierces *1* element with a ray casted from the centroid/pSource
        in the direction of the normal vector of an aerodynamic triangle

         A  1
          \/ \
          / * \
         2---\-3
               \
                B

        *P = A+(B-A)*t
        """
        #direction = -1. # TODO: direction of normal...?
        (sElements, sDists) = self.centroidTree.getCloseElementIDs(aCentroid)
        log.info("aCentroid = %s" % aCentroid)
        log.info("sElements = %s" % sElements)
        log.info("sDists    = %s" % ListPrint(sDists))
        #(nearbySElements, nearbyDistances) = sElements
        piercedElements = []

        for sEID, sDist in zip(sElements, sDists):
            #print "aEID=%s sEID=%s" % (aEID, sEID)
            (sArea, sNormal, sCentroid) = self.structuralModel.get_element_properties(sEID)
            sNodes = self.structuralModel.get_element_nodes(sEID)
            nNodes = len(sNodes)

            pEnd = pSource + normal * 10.
            #pEnd2 = pSource - normal * 10.
            if nNodes == 3:  # TODO:  is this enough of a breakdown?
                (sA, sB, sC) = sNodes
                #pEnd = pSource+normal*10.
                tuv  = piercePlaneVector(sA, sB, sC, pSource, pEnd, piercedElements)
                #tuv2= piercePlaneVector(sA, sB, sC, pSource, pEnd2, piercedElements)
            elif nNodes == 4:
                (sA, sB, sC, sD) = sNodes
                tuv  = piercePlaneVector(sA, sB, sC, pSource, pEnd,  piercedElements)
                #tuv2= piercePlaneVector(sA, sB, sC, pSource, pEnd2, piercedElements)
                #self.pierceTriangle(sA, sB, sC, sCentroid, sNormal, piercedElements)
                #self.pierceTriangle(sA, sC, sD, sCentroid, sNormal, piercedElements)
            else:
                raise Exception('invalid element; nNodes=%s' %(nNodes))

            t1, u1, v1 = tuv
            #t2, u2, v2 = tuv2

            isInside = False
            #if self.isInside(u1, v1) or self.isInside(u2, v2):
            if self.isInside(u1, v1):
                isInside = True
                #pIntersect = pSource + (pEnd - pSource) * t1
                pIntersect = pEnd * t1 +pSource * (1 - t1)
                #P = A + (B - A) * t
                tuv  = piercePlaneVector(sA, sB, sC, pSource, pIntersect, piercedElements)
                #print "t,u,v=", tuv

                piercedElements.append([sEID, pIntersect, u1, v1, sDist])

            #t = min(t1, t2)
            #print "t1=%6.3g t2=%6.3g" % (t1, t2)
            #if isInside:
                #print "*t[%s]=%6.3g u1=%6.3g v1=%6.3g u2=%6.3g v2=%6.3g" %(sEID,t,u1,v1,u2,v2)
            #else:
                #print " t[%s]=%6.3g u1=%6.3g v1=%6.3g u2=%6.3g v2=%6.3g" %(sEID,t,u1,v1,u2,v2)

            #if isInside:
                #print "*t[%s]=%6.3g u1=%6.3g v1=%6.3g d=%g" %(sEID,t1,u1,v1,sDist)
            #else:
                #print " t[%s]=%6.3g u1=%6.3g v1=%6.3g d=%g" %(sEID,t1,u1,v1,sDist)

        log.info("avgDist = %g" % mean(sDists))
        (piercedElements, nPiercings) = self.fixPiercings(sElements, piercedElements)
        distribution = self.distributeUnitLoad(aEID, piercedElements, nPiercings)

        return (distribution)

    def fixPiercings(self, sElements, piercedElements):
        if len(piercedElements) == 0:
            piercedElements = sElements
            #for sEID in sElements:
                #piercedElements.append([sEID,None])  # TODO: why a None?
            nPiercings = 0
        else:
            dists = []
            for element in piercedElements:
                log.info("element = %s" % element)
                dist = element[-1]
                log.info("dist = %s\n" % dist)
                dists.append(dist)
            iSort = argsort(dists)
            #print "iSort = ", iSort

            piercedElements2 = []
            for iElement in iSort:
                piercedElements2.append(piercedElements[iElement])
            #piercedElements = piercedElements[iSort]

            #for element in piercedElements:
                #print "element = ",element
                #print "dist = ",dist, '\n'

            nPiercings = len(piercedElements)
        return (piercedElements, nPiercings)

    def isInside(self, u, v):
        if (0.<=u and u<=1.) and (0.<=v and v<=1.):
            return True
        return False

    def createDistribution(self, nIDs, weights):
        """
        Maps alist of structural nodes, and weights for a given aero element
        Takes the weights that are applied to a node and distributes them to the
        structural nodes
        """
        distribution = {}
        for nid,weight in zip(nIDs, weights):
            distribution[nid] = weight
        return distribution

    def getWeights(self, closePoint, nodes):
        # TODO: new weights?
        #(n1, n2, n3) = list(setNodes)
        #n = piercedPoint
        #(w1, w2, w3) = getTriangleWeights(n, n1, n2, n3)

        weights = shepardWeight(closePoint, nodes)
        return weights

    def distributeUnitLoad(self, aEID, piercedElements, nPiercings):
        """
        distribute unit loads to nearby nodes
        piercedElements is a list of piercings
        piercing = [sEID,P,u,v] or [sEID]
        where
          sEID - the structural element ID
          P    - point p was pierced
          u    - u coordinate
          v    - v coordinate
        if nPiercings==0, all the nearby nodes will recieve load
        """
        aModel = self.aeroModel
        sModel = self.structuralModel
        #print "piercedElements = ",piercedElements
        nIDs = []
        if nPiercings == 0:
            #assert len(nPiercings)==1,'fix me...'
            #print "nPiercings=0 distributing load to closest nodes...u=%g v=%g" %(-1,-1)
            log.info("nPiercings=0 distributing load to closest nodes...")
            for sEID in piercedElements:
                nIDs += sModel.get_element_node_ids(sEID)
            #print "nIDs1 = ",nIDs
            nIDs = list(set(nIDs))
            log.info("nIDs2 = %s" % nIDs)
            aCentroid = aModel.Centroid(aEID)
            nodes = sModel.getNodeIDLocations(nIDs)

            #print "nodes = ", nodes
            weights = self.getWeights(aCentroid, nodes)
            distribution = self.createDistribution(nIDs, weights)

            log.info("element aEID=%s sEID=%s weights=%s" % (aEID, sEID, ListPrint(weights)))
            #print "distribution = ", distribution
            #print "nIDs         = ", nIDs
            #print "weights      = ", weights
            #print "nodes = ", nodes
            #print "nPiercings = ", nPiercings
        else:
            log.info("mapping load to actual element...")
            nClose = 3  # number of elements to map to
            closeElements = piercedElements[:nClose]

            setCloseNodes = set([])
            for closeElement in reversed(closeElements):
                log.info("closeElement = %s" % closeElement)
                #sEID, pIntersect, u1, v1, sDist
                (sEID, P, u, v, sDist) = closeElement  # TODO:  bug here...???

                #closePoint = closeElement[1]
                #closeElement = sEID
                closePoint = P

                # get list of nearby structural nodes
                setElementNodes = set(sModel.get_element_node_ids(sEID))
                setCloseNodes = setCloseNodes.union(setElementNodes)

            # setup for weighted average
            nIDs = list(setCloseNodes)
            sNodes = sModel.getNodeIDLocations(nIDs)
            weights = self.getWeights(closePoint, sNodes)
            distribution = self.createDistribution(nIDs, weights)

            log.info("element aEID=%s sEID=%s weights=%s" %(aEID, sEID, ListPrint(weights)))
        log.info("-------------------------\n")
        sys.stdout.flush()
        return (distribution)

    def Normal(self, A, B, C):
        a = B - A
        b = C - A
        normal = Normal(a, b)
        return normal

#------------------------------------------------------------------

def run_map_loads(cart3dGeom='Components.i.triq', bdfModel='fem.bdf', bdfModelOut='fem.loads.out'):
    assert os.path.exists(bdfModel), '%r doesnt exist' % bdfModel

    t0 = time()
    mesh = Cart3DReader()
    half_model = cart3dGeom+'_half'

    result_names = ['Cp']
    (nodes, elements, regions, loads) = mesh.read_cart3d(cart3dGeom, result_names=result_names)
    #Cp = loads['Cp']
    #half_model_file =
    (nodes, elements, regions, Cp) = mesh.make_half_model(nodes, elements, regions, loads)
    #, os.path.basename(half_model)

    #(nodes, elements, regions, Cp) = mesh.renumber_mesh(nodes, elements, regions, Cp)
    mesh.write_cart3d(half_model, nodes, elements, regions)

    Mach = 0.825
    pInf = 499.3        # psf, alt=35k (per Schaufele p. 11)
    pInf = pInf / 144.  # convert to psi
    qInf = 1.4 / 2. * pInf * Mach**2.
    aeroModel = AeroModel(nodes, elements, Cp, pInf, qInf)
    log.info("elements[1] = %s" % elements[1])
    del elements, nodes, Cp


    fem = BDF(debug=True, log=log)
    fem.read_bdf(bdfModel)
    sys.stdout.flush()

    propertyRegions = [1, 1101, 1501, 1601, 1701, 1801, 1901, 2101, 2501, 2601, 2701, 2801, 2901, 10103, 10201, 10203, 10301, 10401, 10501, 10601, 10701, 10801, 10901, 20103, 20203, 20301, 20401, 20501, 20601, 20701, 20801, 20901, 701512, 801812]

    # 1 inboard
    # 1000s upper - lower inboard
    # 2000s lower - lower inboard
    # big - fin

    structuralModel = StructuralModel(fem, propertyRegions)

    mapper = LoadMapping(aeroModel, structuralModel)
    t1 = time()
    mapper.set_flight_condition(pInf, qInf)
    mapper.setOutput(bdffile=bdfModelOut, loadCase=1)
    log.info("setup time = %g sec; %g min" % (t1-t0, (t1-t0)/60.))

    mapper.buildMappingMatrix(debug=False)
    t2 = time()
    log.info("mapping matrix time = %g sec; %g min" % (t2-t1, (t2-t1)/60.))

    mapper.mapLoads()
    t3 = time()
    log.info("map loads time = %g sec" % (t3 - t2))
    log.info("total time = %g min" % ((t3 - t0) / 60.))

if __name__=='__main__':
    basepath = os.getcwd()
    configpath = os.path.join(basepath, 'inputs')
    workpath   = os.path.join(basepath, 'outputs')
    cart3dGeom = os.path.join(configpath, 'Cart3d_35000_0.825_10_0_0_0_0.i.triq')

    bdfModel   = os.path.join(configpath, 'aeroModel_mod.bdf')
    assert os.path.exists(bdfModel), '%r doesnt exist' % bdfModel
    os.chdir(workpath)
    log.info("basepath = %s" % basepath)

    bdfModelOut = os.path.join(workpath, 'fem_loads_3.bdf')
    run_map_loads(cart3dGeom, bdfModel, bdfModelOut)

