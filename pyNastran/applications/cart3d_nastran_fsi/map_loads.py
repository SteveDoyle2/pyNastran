from __future__ import print_function
from six import iteritems
from six.moves import zip
import os
import sys
import multiprocessing as mp
import cPickle
from time import time

from numpy import argsort, mean, array, cross


from pyNastran.applications.cart3d_nastran_fsi.math_functions import (
    pierce_plane_vector, shepard_weight, Normal, ListPrint)
#from pyNastran.applications.cart3d_nastran_fsi.math_functions import get_triangle_weights
from pyNastran.applications.cart3d_nastran_fsi.structural_model import StructuralModel
from pyNastran.applications.cart3d_nastran_fsi.aero_model import AeroModel
from pyNastran.applications.cart3d_nastran_fsi.kdtree import KdTree

from pyNastran.converters.cart3d.cart3d import Cart3D
from pyNastran.bdf.bdf import BDF
from pyNastran.utils.log import get_logger

debug = True
log = get_logger(None, 'debug' if debug else 'info')


def entry_exit(func):
    def new_f(*args, **kwargs):
        log.info("Entering", func.__name__)
        func(*args, **kwargs)
        log.info("Exited", func.__name__)
    return new_f

class LoadMapping(object):
    def __init__(self, aeroModel, structuralModel, configpath='', workpath=''):
        self.nCloseElements = 30

        self.configpath = configpath
        self.workpath = workpath
        self.aero_model = aeroModel
        self.structural_model = structuralModel

        self.mapping_matrix = {}
        self.mapfile = 'mapfile.in'
        self.bdffile = 'fem.bdf'

        self.centroid_tree = None
        self.nodal_tree = None
        self.pInf = None
        self.qInf = None
        self.load_case = None
        self.bdf = None
        self.load_cases = None

        self.Sref = 1.  # inches
        self.Lref = 1.  # inch
        self.xref = 0.  # inch

    #@entry_exit
    def set_output(self, bdffile='fem.bdf', load_case=1):
        self.bdffile = bdffile
        self.load_case = load_case

    #@entry_exit
    def set_flight_condition(self, pInf=499.3, qInf=237.885):
        self.pInf = pInf
        self.qInf = qInf  #1.4/2.*pInf*Mach**2.

    def set_reference_quantities(self, Sref=1, Lref=1., xref=0.):
        """
        Parameters
        ----------
        Sref : float
            the reference area (in^2)
        Lref : float
            the force/moment reference length (in)
        xref : float
            the moment reference position (in)
        """
        self.Sref = Sref
        self.Lref = Lref
        self.xref = xref

    #@entry_exit
    def load_mapping_matrix(self):
        with open(self.mapfile, 'r') as infile:
            self.mapping_matrix = cPickle.loads(infile)

    #@entry_exit
    def save_mapping_matrix(self):
        out_string = cPickle.dumps(self.mapping_matrix)
        with open(self.mapfile, 'wb') as outfile:
            outfile.write(out_string)

    #@entry_exit
    def get_pressure(self, Cp):
        """
        Caculates the pressure based on:
           Cp = (p-pInf)/qInf
        """
        p = Cp * self.qInf + self.pInf
        #p = Cp*self.qInf # TODO:  for surface with atmospheric pressure inside
        return p

    #@entry_exit
    def map_loads(self):
        """
        Loops thru the unitLoad mapping_matrix and multiplies by the
        total force F on each panel to calculate the PLOAD that will be
        printed to the output file.
        Also, performs a sum of Forces & Moments to verify the aero loads.
        """
        log.info("---starting map_loads---")
        self.bdf = open(self.bdffile, 'wb')
        #self.build_mapping_matrix()
        log.info("self.load_case = %s" % self.load_case)
        self.load_cases = {
            self.load_case : {},
        }

        #self.loadCases = {self.load_case={}, }
        moment_center = array([self.xref, 0., 0.])
        sum_forces = array([0., 0., 0.])
        sum_moments = array([0., 0., 0.])
        sys.stdout.flush()
        for aero_eid, distribution in iteritems(self.mapping_matrix):
            #print("aero_eid = ", aero_eid)
            #print("***distribution = ", distribution)
            sum_load = 0.
            area = self.aero_model.Area(aero_eid)
            normal = self.aero_model.Normal(aero_eid)
            Cp = self.aero_model.Cp(aero_eid)
            #print("Cp = ", Cp)
            #print("area[%s]=%s" % (aero_eid, area))

            p = self.get_pressure(Cp)
            centroid = self.aero_model.Centroid(aero_eid)
            r = moment_center - centroid
            F = area * p
            Fn = F * normal
            sum_moments += cross(r, Fn)
            sum_forces += Fn
            for sNID, percent_load in sorted(iteritems(distribution)):
                sum_load += percent_load

                Fxyz = Fn * percent_load  # negative sign is to be consistent with nastran
                self.add_force(sNID, Fxyz)

                #print("Fxyz = ",Fxyz)
                #print("type(structuralModel) = ", type(self.structuralModel))

                #comment = 'percent_load=%.2f' % percent_load
                #self.structuralModel.write_load(self.bdf, self.loadCase, sNID,
                #                                Fxyz[0], Fxyz[1], Fxyz[2], comment)

            #msg = '$ End of aEID=%s sumLoad=%s p=%s area=%s F=%s normal=%s\n' % (
                #aEID, sumLoad, p, area, F, normal)
            #self.bdf.write(msg)

        self.write_loads()  # short version of writing loads...
        self.bdf.close()

        log.info("pInf=%g [psi]; qInf= %g [psi]" % (self.pInf, self.qInf))
        log.info("sumForcesFEM  [lb]    = %s" % ListPrint(sum_forces))

        # divided by 12 to have moments in lb-ft
        log.info("sumMomentsFEM [lb-ft] = %s" % ListPrint(sum_moments / 12.))

        Cf = sum_forces /(self.Sref * self.qInf)
        log.info("Cf = %s" % ListPrint(Cf))

        Cm = sum_moments / (self.Sref * self.qInf * self.Lref)
        # multiply by 12 to nondimensionalize ???  maybe 144...
        log.info("Cm = %s" % ListPrint(Cm * 12.))

        #self.bdf.write('$***********\n')
        log.info("wrote loads to %r" % self.bdffile)
        log.info("---finished map_loads---")

    #@entry_exit
    def write_loads(self):
        """writes the load in BDF format"""
        log.info("---starting write_loads---")
        self.bdf.write('$ ***write_loads***\n')
        self.bdf.write('$ nCloseElements=%s\n' % self.nCloseElements)
        for load_case, loads in sorted(iteritems(self.load_cases)):
            log.info("  load_case = %s" % load_case)
            for (structural_nid, Fxyz) in sorted(iteritems(loads)):
                self.structural_model.write_load(self.bdf, load_case, structural_nid, *Fxyz)

        log.info("finished write_loads---")

    def add_force(self, structural_nid, Fxyz):
        try:
            self.load_cases[self.load_case][structural_nid] += Fxyz
        except KeyError:
            self.load_cases[self.load_case][structural_nid] = Fxyz

    #@entry_exit
    def build_centroids(self, model, eids=None):
        centroids = {}
        if eids is None:
            eids = model.ElementIDs()
        for eid in eids:
            centroid = model.Centroid(eid)
            if centroid is not None:
                centroids[eid] = centroid
        return centroids

    #@entry_exit
    def build_nodal_tree(self, sNodes):
        log.info("---start build_nodal_tree---")
        raise RuntimeError('DEPRECATED...build_nodal_tree in mapLoads.py')
        #sys.stdout.flush()
        #print("type(aCentroids)=%s type(sCentroids)=%s" %(type(aCentroids), type(sCentroids)))
        #self.nodal_tree = KdTree('node', sNodes, nclose=self.nCloseNodes)
        #log.info("---finish build_nodal_tree---")
        #sys.stdout.flush()

    #@entry_exit
    def build_centroid_tree(self, structural_centroids):
        """
        structural_centroids - dict of structural centroids
        id:  element id
        """
        log.info("---start build_centroid_tree---")
        sys.stdout.flush()
        #print("type(aCentroids)=%s type(structural_centroids)=%s" % (
            #type(aCentroids), type(structural_centroids)))

        msg = 'Element '
        for eid, structural_centroid in sorted(iteritems(structural_centroids)):
            msg += "%s " % eid
        log.info(msg + '\n')

        self.centroid_tree = KdTree('element', structural_centroids, nclose=self.nCloseElements)
        log.info("---finish build_centroid_tree---")
        sys.stdout.flush()

    #@entry_exit
    def parse_map_file(self, map_filename='mappingMatrix.new.out'):
        """
        This is used for rerunning an analysis quickly
        (cuts out building the mapping matrix ~1.5 hours).
        1 {8131: 0.046604568185355716, etc...}
        """
        log.info("---starting parse_map_file---")
        mapping_matrix = {}

        log.info('loading map_file=%r' % map_filename)
        with open(map_filename, 'r') as map_file:
            lines = map_file.readlines()

        # dont read the first line, thats a header line
        for (i, line) in enumerate(lines[1:]):
            line = line.strip()
            #print("line = %r" % line)
            (aEID, dict_line) = line.split('{') # splits the dictionary from the aEID
            aEID = int(aEID)
            #assert i == int(aEID)

            # time to parse a dictionary with no leading brace
            distribution = {}
            map_sline = dict_line.strip('{}, ').split(',')
            for pair in map_sline:
                (sEID, weight) = pair.strip(' ').split(':')
                sEID = int(sEID)
                weight = float(weight)
                distribution[sEID] = weight
            mapping_matrix[aEID] = distribution
        #log.info("mappingKeys = %s" %(sorted(mapping_matrix.keys())))
        self.run_map_test(mapping_matrix)
        log.info("---finished parse_map_file---")
        return mapping_matrix

    #@entry_exit
    def run_map_test(self, mapping_matrix, map_test_filename='map_test.out'):
        """
        Checks to see what elements loads were mapped to.
        Ideally, all elements would have a value of 1.0, given equal area.
        Instead, each element should have a value of area[i].
        """
        map_test = {}
        for (aero_eid, distribution) in sorted(iteritems(mapping_matrix)):
            for sEID, weight in distribution.items():
                if sEID in map_test:
                    map_test[sEID] += weight
                else:
                    map_test[sEID] = weight

        with open(map_test_filename, 'wb') as map_out:
            map_out.write('# sEID  weight\n')
            for sEID, weight in sorted(iteritems(map_test)):
                map_out.write("%s %s\n" % (sEID, weight))

    #@entry_exit
    def build_mapping_matrix(self, debug=False):
        """
        Skips building the matrix if it already exists
        A mapping matrix translates element ID to loads on the nearby
        strucutral nodes.

        eid,distribution
        """
        if self.mapping_matrix != {}:
            return self.mapping_matrix

        log.info("---starting build_mapping_matrix---")
        #print("self.mapping_matrix = ",self.mapping_matrix)
        mapping_matrix_filename = os.path.join(self.configpath, 'mappingMatrix.new.out')
        if os.path.exists(mapping_matrix_filename):
            self.mapping_matrix = self.parse_map_file(mapping_matrix_filename)
            log.info("---finished build_mapping_matrix based on mappingMatrix.new.out---")
            sys.stdout.flush()
            return self.mapping_matrix
        log.info("...couldn't find 'mappingMatrix.new.out' in %r"
                 ", so going to make it..." % os.getcwd())

        # this is the else...
        log.info("creating...")
        aero_model = self.aero_model
        structural_model = self.structural_model

        #aNodes = aero_model.getNodes()
        #sNodes = structural_model.getNodes()
        #treeObj = Tree(nClose=5)
        #tree    = treeObj.buildTree(aNodes,sNodes) # fromNodes,toNodes

        aero_element_ids = aero_model.ElementIDs() # list
        structural_element_ids = structural_model.getElementIDsWithPIDs() # list
        structural_element_ids2 = structural_model.ElementIDs() # list

        msg = 'there are no internal elements in the structural model?\n'
        msg += '   ...len(structural_element_ids)=%s len(structural_element_ids2)=%s' % (
            len(structural_element_ids), len(structural_element_ids2))
        assert structural_element_ids != structural_element_ids2, msg
        log.info("maxAeroID=%s maxStructuralID=%s sElements=%s" % (
            max(aero_element_ids), max(structural_element_ids), len(structural_element_ids2)))

        log.info("build_centroids - structural")
        sCentroids = self.build_centroids(structural_model, structural_element_ids)
        self.build_centroid_tree(sCentroids)
        #self.buildNodalTree(sNodes)

        log.info("build_centroids - aero")
        aero_centroids = self.build_centroids(aero_model)

        with open('mappingMatrix.out', 'wb') as map_file:
            map_file.write('# aEID distribution (sEID:  weight)\n')

            t0 = time()
            nAeroElements = float(len(aero_element_ids))
            log.info("---start piercing---")
            if debug:
                log.info("nAeroElements = %s" % nAeroElements)
            tEst = 1.
            tLeft = 1.
            percent_done = 0.

            use_multiprocessing = False
            if use_multiprocessing:
                # doesn't actually work because multiprocessing can't take self as an arguement
                # the sub-code for map_loads_mp_func needs to be updated
                num_cpus = 4
                pool = mp.Pool(num_cpus)
                result = pool.imap(map_loads_mp_func,
                                   [(aEID, aero_model) for aEID in aero_element_ids])

                for j, return_values in enumerate(result):
                    aEID, distribution = return_values
                    #self.mappingMatrix[aEID] = distribution
                    map_file.write('%s %s\n' % (aEID, distribution))
                pool.close()
                pool.join()
            else:
                for (i, aero_eid) in enumerate(aero_element_ids):
                    if i % 1000 == 0 and debug or 1:
                        #log.debug('  piercing %sth element' % i)
                        log.debug("tEst=%g minutes; tLeft=%g minutes; %.3f%% done" % (
                            tEst, tLeft, percent_done))
                        sys.stdout.flush()

                    aElement = aero_model.Element(aero_eid)
                    (aArea, aCentroid, aNormal) = aero_model.get_element_properties(aero_eid)
                    percentDone = i / nAeroElements * 100
                    if debug:
                        log.info('aEID=%s percentDone=%.2f aElement=%s aArea=%s aCentroid=%s aNormal=%s' % (
                            aero_eid, percentDone, aElement, aArea, aCentroid, aNormal))
                    pSource = aCentroid
                    (distribution) = self.pierce_elements(aCentroid, aero_eid, pSource, aNormal)
                    #(distribution)  = self.poorMansMapping(aCentroid, aero_eid, pSource, aNormal)
                    self.mapping_matrix[aero_eid] = distribution
                    map_file.write('%s %s\n' % (aero_eid, distribution))

                    dt = (time() - t0) / 60.
                    tEst = dt * nAeroElements / (i + 1)  # dtPerElement*nElements
                    tLeft = tEst - dt
                    percent_done = dt / tEst * 100.

        log.info("---finish piercing---")
        self.run_map_test(self.mapping_matrix)
        #print("mapping_matrix = ", self.mapping_matrix)
        log.info("---finished build_mapping_matrix---")
        sys.stdout.flush()
        return self.mapping_matrix

    def poor_mans_mapping(self, aCentroid, aero_eid, pSource, normal):
        """
        distributes load without piercing elements
        based on distance
        """
        (sElements, sDists) = self.centroid_tree.getCloseElementIDs(aCentroid)
        log.debug("aCentroid = %s" % aCentroid)
        log.debug("sElements = %s" % sElements)
        log.debug("sDists    = %s" % ListPrint(sDists))

        setNodes = set([])
        structural_model = self.structural_model
        for structural_eid in sElements:
            sNodes = structural_model.get_element_nodes(structural_eid)
            setNodes.union(set(sNodes))

        nIDs = list(setNodes)
        sNodes = structural_model.getNodeIDLocations(nIDs)
        weights = self.get_weights(close_point, sNodes)
        distribution = self.create_distribution(nIDs, weights)
        return distribution

    def pierce_elements(self, aCentroid, aEID, pSource, normal):
        r"""
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
        (sElements, sDists) = self.centroid_tree.getCloseElementIDs(aCentroid)
        log.info("aCentroid = %s" % aCentroid)
        log.info("sElements = %s" % sElements)
        log.info("sDists    = %s" % ListPrint(sDists))
        #(nearbySElements, nearbyDistances) = sElements
        pierced_elements = []

        for sEID, sDist in zip(sElements, sDists):
            #print("aEID=%s sEID=%s" % (aEID, sEID))
            sArea, sNormal, sCentroid = self.structural_model.get_element_properties(sEID)
            sNodes = self.structural_model.get_element_nodes(sEID)
            nNodes = len(sNodes)

            pEnd = pSource + normal * 10.
            #pEnd2 = pSource - normal * 10.
            if nNodes == 3:  # TODO:  is this enough of a breakdown?
                sA, sB, sC = sNodes
                #pEnd = pSource+normal*10.
                tuv = pierce_plane_vector(sA, sB, sC, pSource, pEnd, pierced_elements)
                #tuv2 = pierce_plane_vector(sA, sB, sC, pSource, pEnd2, pierced_elements)
            elif nNodes == 4:
                sA, sB, sC, sD = sNodes
                tuv = pierce_plane_vector(sA, sB, sC, pSource, pEnd, pierced_elements)
                #tuv2 = pierce_plane_vector(sA, sB, sC, pSource, pEnd2, pierced_elements)
                #self.pierceTriangle(sA, sB, sC, sCentroid, sNormal, pierced_elements)
                #self.pierceTriangle(sA, sC, sD, sCentroid, sNormal, pierced_elements)
            else:
                raise RuntimeError('invalid element; nNodes=%s' % nNodes)

            t1, u1, v1 = tuv
            #t2, u2, v2 = tuv2

            is_inside_bool = False
            #if is_inside(u1, v1) or is_inside(u2, v2):
            if is_inside(u1, v1):
                is_inside_bool = True
                #pIntersect = pSource + (pEnd - pSource) * t1
                pIntersect = pEnd * t1 +pSource * (1 - t1)
                #P = A + (B - A) * t
                tuv = pierce_plane_vector(sA, sB, sC, pSource, pIntersect, pierced_elements)
                #print("t,u,v=", tuv)

                pierced_elements.append([sEID, pIntersect, u1, v1, sDist])

            #t = min(t1, t2)
            #print("t1=%6.3g t2=%6.3g" % (t1, t2))
            #if is_inside_bool:
                #print("*t[%s]=%6.3g u1=%6.3g v1=%6.3g u2=%6.3g v2=%6.3g" %(sEID,t,u1,v1,u2,v2))
            #else:
                #print(" t[%s]=%6.3g u1=%6.3g v1=%6.3g u2=%6.3g v2=%6.3g" %(sEID,t,u1,v1,u2,v2))

            #if is_inside_bool:
                #print("*t[%s]=%6.3g u1=%6.3g v1=%6.3g d=%g" %(sEID,t1,u1,v1,sDist))
            #else:
                #print(" t[%s]=%6.3g u1=%6.3g v1=%6.3g d=%g" %(sEID,t1,u1,v1,sDist))

        log.info("avgDist = %g" % mean(sDists))
        (pierced_elements, npiercings) = self.fix_piercings(sElements, pierced_elements)
        distribution = self.distribute_unit_load(aEID, pierced_elements, npiercings)

        return distribution

    def fix_piercings(self, sElements, pierced_elements):
        if len(pierced_elements) == 0:
            pierced_elements = sElements
            #for sEID in sElements:
                #pierced_elements.append([sEID, None])  # TODO: why a None?
            npiercings = 0
        else:
            dists = []
            for element in pierced_elements:
                log.info("element = %s" % element)
                dist = element[-1]
                log.info("dist = %s\n" % dist)
                dists.append(dist)
            iSort = argsort(dists)
            #print("iSort = ", iSort)

            piercedElements2 = []
            for iElement in iSort:
                piercedElements2.append(pierced_elements[iElement])
            #piercedElements = piercedElements[iSort]

            #for element in pierced_elements:
                #print("element = ",element)
                #print("dist = ",dist, '\n')

            npiercings = len(pierced_elements)
        return (pierced_elements, npiercings)

    def create_distribution(self, nIDs, weights):
        """
        Maps alist of structural nodes, and weights for a given aero element
        Takes the weights that are applied to a node and distributes them to the
        structural nodes
        """
        distribution = {}
        for nid, weight in zip(nIDs, weights):
            distribution[nid] = weight
        return distribution

    def get_weights(self, closePoint, nodes):
        # TODO: new weights?
        #n1, n2, n3 = list(setNodes)
        #n = piercedPoint
        #w1, w2, w3 = getTriangleWeights(n, n1, n2, n3)

        weights = shepard_weight(closePoint, nodes)
        return weights

    def distribute_unit_load(self, aero_eid, pierced_elements, npiercings):
        """
        distribute unit loads to nearby nodes
        pierced_elements is a list of piercings
        piercing = [sEID,P,u,v] or [sEID]
        where
          sEID - the structural element ID
          P    - point p was pierced
          u    - u coordinate
          v    - v coordinate
        if npiercings==0, all the nearby nodes will recieve load
        """
        aero_model = self.aero_model
        structural_model = self.structural_model
        #print("pierced_elements = ", pierced_elements)
        nIDs = []
        if npiercings == 0:
            #assert len(npiercings)==1,'fix me...'
            #print("npiercings=0 distributing load to closest nodes...u=%g v=%g" %(-1,-1))
            log.debug("npiercings=0 distributing load to closest nodes...")
            for structural_eid in pierced_elements:
                nIDs += structural_model.get_element_node_ids(structural_eid)
            #print("nIDs1 = ", nIDs)
            nIDs = list(set(nIDs))
            log.debug("nIDs2 = %s" % nIDs)
            aCentroid = aero_model.Centroid(aero_eid)
            nodes = structural_model.getNodeIDLocations(nIDs)

            #print("nodes = ", nodes)
            weights = self.get_weights(aCentroid, nodes)
            distribution = self.create_distribution(nIDs, weights)

            log.debug("element aEID=%s sEID=%s weights=%s" % (
                aero_eid, structural_eid, ListPrint(weights)))
            #print("distribution = ", distribution)
            #print("nIDs         = ", nIDs)
            #print("weights      = ", weights)
            #print("nodes = ", nodes)
            #print("npiercings = ", npiercings)
        else:
            log.info("mapping load to actual element...")
            nclose = 3  # number of elements to map to
            close_elements = pierced_elements[:nclose]

            set_close_nodes = set([])
            for close_element in reversed(close_elements):
                log.info("close_element = %s" % close_element)
                #sEID, pIntersect, u1, v1, sdist
                structural_eid, P, u, v, sdist = close_element  # TODO:  bug here...???

                #close_point = close_element[1]
                #close_element = structural_eid
                close_point = P

                # get list of nearby structural nodes
                set_element_nodes = set(structural_model.get_element_node_ids(structural_eid))
                set_close_nodes = set_close_nodes.union(set_element_nodes)

            # setup for weighted average
            nIDs = list(set_close_nodes)
            sNodes = structural_model.getNodeIDLocations(nIDs)
            weights = self.get_weights(close_point, sNodes)
            distribution = self.create_distribution(nIDs, weights)

            log.info("element aEID=%s sEID=%s weights=%s" % (
                aero_eid, structural_eid, ListPrint(weights)))
        log.info("-------------------------\n")
        sys.stdout.flush()
        return distribution

    def Normal(self, A, B, C):
        a = B - A
        b = C - A
        normal = Normal(a, b)
        return normal

def is_inside(u, v):
    if (0. <= u <= 1.) and (0. <= v <= 1.):
        return True
    return False

#def map_loads_mp_func(aero_eid, aero_model):
    #aero_element = aero_model.Element(aero_eid)
    #(aero_area, aero_centroid, aero_normal) = aero_model.get_element_properties(aero_eid)
    ##percentDone = i / nAeroElements * 100

    #pSource = aero_centroid
    #distribution = pierce_elements(aero_centroid, aero_eid, pSource, aero_normal)
    ##mapping_matrix[aero_eid] = distribution
    #return aero_eid, distribution

#def pierce_elements(aCentroid, aEID, pSource, normal):
    #"""
    #see LoadMapping.pierce_elements
    #"""
    #(sElements, sDists) = centroid_tree.getCloseElementIDs(aCentroid)
    #pierced_elements = []

    #for sEID, sDist in zip(sElements, sDists):
        #sArea, sNormal, sCentroid = self.structural_model.get_element_properties(sEID)
        #sNodes = self.structural_model.get_element_nodes(sEID)
        #nNodes = len(sNodes)

        #pEnd = pSource + normal * 10.
        #if nNodes == 3:  # TODO:  is this enough of a breakdown?
            #sA, sB, sC = sNodes
            #tuv = pierce_plane_vector(sA, sB, sC, pSource, pEnd, pierced_elements)
        #elif nNodes == 4:
            #sA, sB, sC, sD = sNodes
            #tuv = pierce_plane_vector(sA, sB, sC, pSource, pEnd, pierced_elements)
        #else:
            #raise RuntimeError('invalid element; nNodes=%s' % nNodes)

        #t1, u1, v1 = tuv

        #is_inside_bool = False
        #if is_inside(u1, v1):
            #is_inside_bool = True
            #pIntersect = pEnd * t1 +pSource * (1 - t1)
            ##P = A + (B - A) * t
            #tuv = pierce_plane_vector(sA, sB, sC, pSource, pIntersect, pierced_elements)
            #pierced_elements.append([sEID, pIntersect, u1, v1, sDist])

    #(pierced_elements, npiercings) = fix_piercings(sElements, pierced_elements)
    #distribution = distribute_unit_load(aEID, pierced_elements, npiercings)
    #return distribution

#------------------------------------------------------------------

def run_map_loads(inputs, cart3d_geom='Components.i.triq', bdf_model='fem.bdf',
                  bdf_model_out='fem.loads.out'):
    """
    inputs : dict
        aero_format : str
            'cart3d' is the only option
        skin_property_regions : List[int, int, ...]
            list of property ids for the skin
        isubcase : int
            the SUBCASE number for the load
        pInf : float
            the static pressure (psi)
        qInf : float
            the dynamic pressure (psi)
        Sref : float
            the reference area (in^2)
        Lref : float
            the reference length (in)
        xref : float
            the reference location (in)
    cart3d_geom : str; default='Components.i.triq'
        path to the cart3d model
        only maps half model loads, so:
           '_half' on the end : use the +y side of the model
           without 'half' : cut the model and use the +y half
    bdf_model : str
        the path to the input fem
    bdf_model_out : str
        the location to save the loads
    """
    assert os.path.exists(bdf_model), '%r doesnt exist' % bdf_model

    t0 = time()
    print(inputs)
    aero_format = inputs['aero_format'].lower()

    property_regions = inputs['skin_property_regions']
    isubcase = inputs['isubcase']
    pInf = inputs['pInf']
    qInf = inputs['qInf']
    Sref = inputs['Sref']
    Lref = inputs['Lref']
    xref = inputs['xref']

    if aero_format == 'cart3d':
        cart3d_model = Cart3D(debug=True)
        half_model = cart3d_geom + '_half'
        result_names = ['Cp', 'rho', 'rhoU', 'rhoV', 'rhoW', 'E', 'Mach']
        #result_names = None

        if not os.path.exists(half_model):
            cart3d_model.read_cart3d(cart3d_geom, result_names=result_names)
            (nodes, elements, regions, loads) = cart3d_model.make_half_model(axis='y')
            cart3d_model.nodes = nodes
            cart3d_model.elements = elements
            cart3d_model.regions = regions
            cart3d_model.loads = loads
            #(nodes, elements, regions, Cp) = cart3d_model.renumber_mesh(nodes, elements, regions, Cp)
            cart3d_model.write_cart3d(half_model)
        else:
            cart3d_model.read_cart3d(half_model, result_names=result_names)
            elements = cart3d_model.elements
        Cp = cart3d_model.loads['Cp']
        mesh = cart3d_model

    else:
        raise NotImplementedError('aero_format=%r' % aero_format)

    aero_model = AeroModel(inputs, mesh.nodes, mesh.elements, Cp)
    log.info("elements[1] = %s" % elements[1])
    #del elements, nodes, Cp

    fem = BDF(debug=True, log=log)
    fem.read_bdf(bdf_model)
    sys.stdout.flush()

    structural_model = StructuralModel(fem, property_regions)

    mapper = LoadMapping(aero_model, structural_model)
    t1 = time()
    mapper.set_flight_condition(pInf, qInf)
    mapper.set_reference_quantities(Sref, Lref, xref)
    mapper.set_output(bdffile=bdf_model_out, load_case=isubcase)
    log.info("setup time = %g sec; %g min" % (t1-t0, (t1-t0)/60.))

    mapper.build_mapping_matrix(debug=False)
    t2 = time()
    log.info("mapping matrix time = %g sec; %g min" % (t2-t1, (t2-t1)/60.))

    mapper.map_loads()
    t3 = time()
    log.info("map loads time = %g sec" % (t3 - t2))
    log.info("total time = %g min" % ((t3 - t0) / 60.))


def main():
    """runs the bwb mapping problem"""
    basepath = os.getcwd()
    configpath = os.path.join(basepath, 'inputs')
    workpath = os.path.join(basepath, 'outputs')
    cart3d_geom = os.path.join(configpath, 'Cart3d_35000_0.825_10_0_0_0_0.i.triq')

    bdf_model = os.path.join(configpath, 'aeroModel_mod.bdf')
    assert os.path.exists(bdf_model), '%r doesnt exist' % bdf_model
    if not os.path.exists(workpath):
        os.makedirs(workpath)
    os.chdir(workpath)
    log.info("basepath = %s" % basepath)

    bdf_model_out = os.path.join(workpath, 'fem_loads_3.bdf')

    pInf = 499.3 / 144 # psf -> psi, alt=35k (per Schaufele p. 11)
    mach = 0.825

    # 1 inboard
    # 1000s upper - lower inboard
    # 2000s lower - lower inboard
    # big - fin
    skin_property_regions = [
        1, 1101, 1501, 1601, 1701, 1801, 1901, 2101, 2501, 2601, 2701, 2801,
        2901, 10103, 10201, 10203, 10301, 10401, 10501, 10601, 10701, 10801,
        10901, 20103, 20203, 20301, 20401, 20501, 20601, 20701, 20801, 20901,
        701512, 801812,
    ]

    inputs = {
        'aero_format' : 'cart3d',
        'isubcase' : 1,
        'Mach' : mach,
        'pInf' : pInf,  # convert to psi
        'qInf' : 1.4 / 2. * pInf * mach**2.,  # psi
        'Sref' : 1582876.,  # inch^2
        'Lref' : 623.,  # inch
        'xref' : 268.,  # inch
        'skin_property_regions' : skin_property_regions,
    }
    run_map_loads(inputs, cart3d_geom, bdf_model, bdf_model_out)


if __name__ == '__main__':
    main()
