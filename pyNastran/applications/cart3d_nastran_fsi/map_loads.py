from __future__ import print_function
import os
import sys
import multiprocessing as mp
import cPickle
from time import time

from six import iteritems
from six.moves import zip

from numpy import argsort, mean, array, cross


from pyNastran.applications.cart3d_nastran_fsi.math_functions import (
    pierce_plane_vector, shepard_weight, Normal, list_print)
#from pyNastran.applications.cart3d_nastran_fsi.math_functions import get_triangle_weights
from pyNastran.applications.cart3d_nastran_fsi.structural_model import StructuralModel
from pyNastran.applications.cart3d_nastran_fsi.aero_model import AeroModel
from pyNastran.applications.cart3d_nastran_fsi.kdtree import KdTree

from pyNastran.converters.cart3d.cart3d import Cart3D
from pyNastran.bdf.bdf import BDF
from pyNastran.utils.log import get_logger2

debug = True
log = get_logger2(None, debug=debug)


def entry_exit(func):
    def new_f(*args, **kwargs):
        log.info("Entering %s" % func.__name__)
        func(*args, **kwargs)
        log.info("Exited %s" % func.__name__)
    return new_f

class LoadMapping(object):
    def __init__(self, aero_model, structural_model, configpath, workpath):
        self.nclose_elements = 30

        self.configpath = configpath
        self.workpath = workpath
        self.aero_model = aero_model
        self.structural_model = structural_model

        self.mapping_matrix = {}
        self.map_filename = 'mapfile.in'
        self.bdf_filename = 'fem.bdf'

        self.centroid_tree = None
        self.nodal_tree = None
        self.load_case = None
        #self.bdf = None
        self.load_cases = None

        self.pInf = None
        self.qInf = None
        #self.Sref = (1., 'in^2')
        #self.Lref = (1., 'in')
        #self.xref = (0., 'in')
        self.Sref = None
        self.Lref = None
        self.xref = None
        self.force_scale = None
        self.moment_scale = None

    #@entry_exit
    def set_output(self, bdffile='fem.bdf', load_case=1):
        self.bdffile = bdffile
        self.load_case = load_case

    #@entry_exit
    def set_flight_condition(self, pInf=(499.3, 'psi'), qInf=(237.885, 'psi')):
        self.pInf = pInf
        self.qInf = qInf  #1.4/2.*pInf*Mach**2.

    def set_scale_factors(self, force_scale=(1., 'lb'), moment_scale=(12., 'ft-lb')):
        self.force_scale = force_scale
        self.moment_scale = moment_scale

    def set_reference_quantities(self, Sref=(1., 'in^2'), Lref=(1., 'in'), xref=(0., 'in')):
        """
        Parameters
        ----------
        Sref : (float, str)
            float : the reference area (in^2)
            str : unit
        Lref : (float, str)
            float : the force/moment reference length (in)
            str : unit
        xref : (float, str)
            float : the moment reference position (in)
            str : unit
        """
        self.Sref = Sref
        self.Lref = Lref
        self.xref = xref

    #@entry_exit
    def load_mapping_matrix(self):
        with open(self.map_filename, 'r') as map_file:
            self.mapping_matrix = cPickle.loads(map_file)

    #@entry_exit
    def save_mapping_matrix(self):
        out_string = cPickle.dumps(self.mapping_matrix)
        with open(self.map_filename, 'wb') as map_file:
            map_file.write(out_string)

    #@entry_exit
    def get_pressure(self, Cp):
        """
        Caculates the pressure based on:
           Cp = (p-pInf)/qInf
        """
        p = Cp * self.qInf[0] + self.pInf[0]
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

        pinf, pinf_unit = self.pInf
        qinf, qinf_unit = self.qInf
        sref, sref_unit = self.Sref
        lref, lref_unit = self.Lref
        xref, xref_unit = self.xref
        force_scale, force_scale_unit = self.force_scale
        moment_scale, moment_scale_unit = self.moment_scale

        log.info("self.load_case = %s" % self.load_case)
        self.load_cases = {
            self.load_case : {},
        }
        #self.loadCases = {self.load_case={}, }

        moment_center = array([xref, 0., 0.])
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
            for structural_nid, percent_load in sorted(iteritems(distribution)):
                sum_load += percent_load

                Fxyz = Fn * percent_load  # negative sign is to be consistent with nastran
                self.add_force(structural_nid, Fxyz)

                #print("Fxyz = ",Fxyz)
                #print("type(structuralModel) = ", type(self.structuralModel))

                #comment = 'percent_load=%.2f' % percent_load
                #self.structuralModel.write_load(bdf_file, self.loadCase, structural_nid,
                #                                Fxyz[0], Fxyz[1], Fxyz[2], comment)

            #msg = '$ End of aEID=%s sumLoad=%s p=%s area=%s F=%s normal=%s\n' % (
                #aEID, sumLoad, p, area, F, normal)
            #bdf_file.write(msg)

        with open(self.bdffile, 'wb') as bdf_file:
            #self.build_mapping_matrix()
            self.write_loads(bdf_file)  # short version of writing loads...


        log.info("pInf=%g [%s]; qInf=%g [%s]" % (
            pinf, pinf_unit,
            qinf, qinf_unit))
        log.info("sumForcesFEM  [%s]    = %s" % (
            force_scale_unit,
            list_print(sum_forces * force_scale),
        ))

        # divided by 12 to have moments in lb-ft
        log.info("sumMomentsFEM [%s] = %s" % (
            moment_scale_unit,
            list_print(sum_moments * moment_scale),
        ))

        Cf = sum_forces / (sref * qinf)
        log.info("Cf = %s" % list_print(Cf))

        Cm = sum_moments / (sref * qinf * lref)
        # multiply by 12 to nondimensionalize ???  maybe 144...
        log.info("Cm = %s" % list_print(Cm))

        #bdf_file.write('$***********\n')
        log.info("wrote loads to %r" % self.bdffile)
        log.info("---finished map_loads---")

    #@entry_exit
    def write_loads(self, bdf_file):
        """writes the load in BDF format"""
        log.info("---starting write_loads---")
        bdf_file.write('$ ***write_loads***\n')
        bdf_file.write('$ nCloseElements=%s\n' % self.nclose_elements)
        for load_case, loads in sorted(iteritems(self.load_cases)):
            log.info("  load_case = %s" % load_case)
            for (structural_nid, Fxyz) in sorted(iteritems(loads)):
                self.structural_model.write_load(bdf_file, load_case, structural_nid, *Fxyz)

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
            eids = model.element_ids
        for eid in eids:
            centroid = model.Centroid(eid)
            if centroid is not None:
                centroids[eid] = centroid
        return centroids

    #@entry_exit
    def build_nodal_tree(self, structural_nodes):
        log.info("---start build_nodal_tree---")
        raise RuntimeError('DEPRECATED...build_nodal_tree in mapLoads.py')
        #sys.stdout.flush()
        #print("type(aCentroids)=%s type(sCentroids)=%s" %(type(aCentroids), type(sCentroids)))
        #self.nodal_tree = KdTree('node', structural_nodes, nclose=self.nclose_nodes)
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

        self.centroid_tree = KdTree('element', structural_centroids, nclose=self.nclose_elements)
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
        for line in lines[1:]:
            line = line.strip()
            #print("line = %r" % line)
            (aero_eid, dict_line) = line.split('{') # splits the dictionary from the aero_eid
            aero_eid = int(aero_eid)
            #assert i == int(aero_eid)

            # time to parse a dictionary with no leading brace
            distribution = {}
            map_sline = dict_line.strip('{}, ').split(',')
            for pair in map_sline:
                (structural_eid, weight) = pair.strip(' ').split(':')
                structural_eid = int(structural_eid)
                weight = float(weight)
                distribution[structural_eid] = weight
            mapping_matrix[aero_eid] = distribution
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
            for structural_eid, weight in distribution.items():
                if structural_eid in map_test:
                    map_test[structural_eid] += weight
                else:
                    map_test[structural_eid] = weight

        with open(map_test_filename, 'wb') as map_out:
            map_out.write('# sEID  weight\n')
            for structural_eid, weight in sorted(iteritems(map_test)):
                map_out.write("%s %s\n" % (structural_eid, weight))

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
        mapping_matrix_filename = os.path.join(self.configpath, 'mappingMatrix.out')
        log.debug('self.configpath = %r' % self.configpath)
        log.debug('mapping_matrix_filename = %r' % mapping_matrix_filename)
        if os.path.exists(mapping_matrix_filename):
            self.mapping_matrix = self.parse_map_file(mapping_matrix_filename)
            log.info("---finished build_mapping_matrix based on mappingMatrix.new.out---")
            sys.stdout.flush()
            return self.mapping_matrix
        log.info("...couldn't find 'mappingMatrix.out' in %r"
                 ", so going to make it..." % self.configpath)

        # this is the else...
        log.info("creating...")
        aero_model = self.aero_model
        structural_model = self.structural_model

        #aNodes = aero_model.getNodes()
        #sNodes = structural_model.getNodes()
        #treeObj = Tree(nClose=5)
        #tree = treeObj.buildTree(aNodes,sNodes) # fromNodes, toNodes

        aero_element_ids = aero_model.element_ids # list
        structural_element_ids = structural_model.getElementIDsWithPIDs() # list
        structural_element_ids2 = structural_model.element_ids # list

        msg = 'there are no internal elements in the structural model?\n'
        msg += '   ...len(structural_element_ids)=%s len(structural_element_ids2)=%s' % (
            len(structural_element_ids), len(structural_element_ids2))
        assert structural_element_ids != structural_element_ids2, msg
        log.info("maxAeroID=%s maxStructuralID=%s sElements=%s" % (
            max(aero_element_ids), max(structural_element_ids), len(structural_element_ids2)))

        log.info("build_centroids - structural")
        structural_centroids = self.build_centroids(structural_model, structural_element_ids)
        self.build_centroid_tree(structural_centroids)
        #self.buildNodalTree(sNodes)

        #log.info("build_centroids - aero")
        #aero_centroids = self.build_centroids(aero_model)

        t0 = time()
        naero_elements = float(len(aero_element_ids))
        log.info("---start piercing---")
        if debug:
            log.info("nAeroElements = %s" % naero_elements)
        tEst = 1.
        tLeft = 1.
        percent_done = 0.
        use_multiprocessing = False

        mapping_matrix_filename_out = os.path.join(self.workpath, 'mappingMatrix.out')
        with open(mapping_matrix_filename_out, 'wb') as map_file:
            map_file.write('# aEID distribution (sEID:  weight)\n')

            if use_multiprocessing:
                # doesn't actually work because multiprocessing can't take self as an arguement
                # the sub-code for map_loads_mp_func needs to be updated
                num_cpus = 4
                pool = mp.Pool(num_cpus)
                result = pool.imap(map_loads_mp_func,
                                   [(aero_eid, aero_model) for aero_eid in aero_element_ids])

                for j, return_values in enumerate(result):
                    aero_eid, distribution = return_values
                    #self.mappingMatrix[aero_eid] = distribution
                    map_file.write('%s %s\n' % (aero_eid, distribution))
                pool.close()
                pool.join()
            else:
                for (i, aero_eid) in enumerate(aero_element_ids):
                    if i % 1000 == 0 and debug or 1:
                        #log.debug('  piercing %sth element' % i)
                        log.debug("tEst=%g minutes; tLeft=%g minutes; %.3f%% done" % (
                            tEst, tLeft, percent_done))
                        sys.stdout.flush()

                    aero_element = aero_model.Element(aero_eid)
                    (aero_area, aero_centroid, aero_normal) = aero_model.get_element_properties(
                        aero_eid)
                    percent_done = i / naero_elements * 100
                    if debug:
                        log.info('aEID=%s percentDone=%.2f aElement=%s '
                                 'aArea=%s aCentroid=%s aNormal=%s' % (
                                     aero_eid, percent_done, aero_element, aero_area,
                                     aero_centroid, aero_normal))
                    pSource = aero_centroid
                    distribution = self.pierce_elements(aero_centroid, aero_eid,
                                                        pSource, aero_normal)
                    #distribution  = self.poor_mans_mapping(aero_centroid, aero_eid,
                                                            #pSource, aero_normal)
                    self.mapping_matrix[aero_eid] = distribution
                    map_file.write('%s %s\n' % (aero_eid, distribution))

                    dt = (time() - t0) / 60.
                    tEst = dt * naero_elements / (i + 1.)  # dtPerElement*nElements
                    tLeft = tEst - dt
                    if tEst != 0.0:
                        percent_done = dt / tEst * 100.

        log.info("---finish piercing---")
        self.run_map_test(self.mapping_matrix)
        #print("mapping_matrix = ", self.mapping_matrix)
        log.info("---finished build_mapping_matrix---")
        sys.stdout.flush()
        return self.mapping_matrix

    #def poor_mans_mapping(self, aero_centroid, aero_eid, pSource, normal):
        #"""
        #distributes load without piercing elements
        #based on distance
        #"""
        #(sElements, sDists) = self.centroid_tree.get_close_element_ids(aero_centroid)
        #log.debug("aCentroid = %s" % aero_centroid)
        #log.debug("sElements = %s" % sElements)
        #log.debug("sDists    = %s" % list_print(sDists))

        #setNodes = set([])
        #structural_model = self.structural_model
        #for structural_eid in sElements:
            #sNodes = structural_model.get_element_nodes(structural_eid)
            #setNodes.union(set(sNodes))

        #nIDs = list(setNodes)
        #sNodes = structural_model.getNodeIDLocations(nIDs)
        #weights = self.get_weights(close_point, sNodes)
        #distribution = self.create_distribution(nIDs, weights)
        #return distribution

    def pierce_elements(self, aero_centroid, aero_eid, pSource, normal):
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
        (structural_elements, structural_distances) = self.centroid_tree.get_close_element_ids(
            aero_centroid)
        log.info("aCentroid = %s" % aero_centroid)
        log.info("sElements = %s" % structural_elements)
        log.info("sDists    = %s" % list_print(structural_distances))
        #(nearbySElements, nearbyDistances) = structural_elements
        pierced_elements = []

        for structural_eid, sDist in zip(structural_elements, structural_distances):
            #print("aEID=%s sEID=%s" % (aero_eid, structural_eid))
            struc = self.structural_model.get_element_properties(structural_eid)
            structural_area, structural_normal, structural_centroid = struc

            sNodes = self.structural_model.get_element_nodes(structural_eid)
            nnodes = len(sNodes)

            pEnd = pSource + normal * 10.
            #pEnd2 = pSource - normal * 10.
            if nnodes == 3:  # TODO:  is this enough of a breakdown?
                sA, sB, sC = sNodes
                #pEnd = pSource+normal*10.
                tuv = pierce_plane_vector(sA, sB, sC, pSource, pEnd, pierced_elements)
                #tuv2 = pierce_plane_vector(sA, sB, sC, pSource, pEnd2, pierced_elements)
            elif nnodes == 4:
                sA, sB, sC, sD = sNodes
                tuv = pierce_plane_vector(sA, sB, sC, pSource, pEnd, pierced_elements)
                #tuv2 = pierce_plane_vector(sA, sB, sC, pSource, pEnd2, pierced_elements)
                #self.pierceTriangle(sA, sB, sC, sCentroid, sNormal, pierced_elements)
                #self.pierceTriangle(sA, sC, sD, sCentroid, sNormal, pierced_elements)
            else:
                raise RuntimeError('invalid element; nNodes=%s' % nnodes)

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

                pierced_elements.append([structural_eid, pIntersect, u1, v1, sDist])

            #t = min(t1, t2)
            #print("t1=%6.3g t2=%6.3g" % (t1, t2))
            #if is_inside_bool:
                #print("*t[%s]=%6.3g u1=%6.3g v1=%6.3g u2=%6.3g v2=%6.3g" %(
                    #structural_eid,t,u1,v1,u2,v2))
            #else:
                #print(" t[%s]=%6.3g u1=%6.3g v1=%6.3g u2=%6.3g v2=%6.3g" %(
                    #structural_eid,t,u1,v1,u2,v2))

            #if is_inside_bool:
                #print("*t[%s]=%6.3g u1=%6.3g v1=%6.3g d=%g" %(
                    #structural_eid,t1,u1,v1,sDist))
            #else:
                #print(" t[%s]=%6.3g u1=%6.3g v1=%6.3g d=%g" %(
                    #structural_eid,t1,u1,v1,sDist))

        log.info("avgDist = %g" % mean(structural_distances))
        (pierced_elements, npiercings) = self.fix_piercings(structural_elements, pierced_elements)
        distribution = self.distribute_unit_load(aero_eid, pierced_elements, npiercings)

        return distribution

    def fix_piercings(self, structural_elements, pierced_elements):
        if len(pierced_elements) == 0:
            pierced_elements = structural_elements
            #for structural_eid in structural_elements:
                #pierced_elements.append([structural_eid, None])  # TODO: why a None?
            npiercings = 0
        else:
            dists = []
            for element in pierced_elements:
                log.info("element = %s" % element)
                dist = element[-1]
                log.info("dist = %s\n" % dist)
                dists.append(dist)
            isort = argsort(dists)
            #print("iSort = ", iSort)

            pierced_elements2 = []
            for ielement in isort:
                pierced_elements2.append(pierced_elements[ielement])
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

    def get_weights(self, close_point, nodes):
        # TODO: new weights?
        #n1, n2, n3 = list(setNodes)
        #n = piercedPoint
        #w1, w2, w3 = getTriangleWeights(n, n1, n2, n3)

        weights = shepard_weight(close_point, nodes)
        return weights

    def distribute_unit_load(self, aero_eid, pierced_elements, npiercings):
        """
        distribute unit loads to nearby nodes
        pierced_elements is a list of piercings
        piercing = [structural_eid,P,u,v] or [structural_eid]
        where
          structural_eid - the structural element ID
          P    - point p was pierced
          u    - u coordinate
          v    - v coordinate
        if npiercings==0, all the nearby nodes will recieve load
        """
        aero_model = self.aero_model
        structural_model = self.structural_model
        #print("pierced_elements = ", pierced_elements)
        node_ids = []
        if npiercings == 0:
            #assert len(npiercings)==1,'fix me...'
            #print("npiercings=0 distributing load to closest nodes...u=%g v=%g" %(-1,-1))
            log.debug("npiercings=0 distributing load to closest nodes...")
            structural_eid = None
            for structural_eid in pierced_elements:
                node_ids += structural_model.get_element_node_ids(structural_eid)
            #print("nIDs1 = ", node_ids)
            node_ids = list(set(node_ids))
            log.debug("nIDs2 = %s" % node_ids)
            aero_centroid = aero_model.Centroid(aero_eid)
            nodes = structural_model.getNodeIDLocations(node_ids)

            #print("nodes = ", nodes)
            weights = self.get_weights(aero_centroid, nodes)
            distribution = self.create_distribution(node_ids, weights)

            log.debug("element aEID=%s sEID=%s weights=%s" % (
                aero_eid, structural_eid, list_print(weights)))
            #print("distribution = ", distribution)
            #print("nIDs         = ", node_ids)
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
                #structural_eid, pIntersect, u1, v1, sdist
                structural_eid, P, u, v, sdist = close_element  # TODO:  bug here...???

                #close_point = close_element[1]
                #close_element = structural_eid
                close_point = P

                # get list of nearby structural nodes
                set_element_nodes = set(structural_model.get_element_node_ids(structural_eid))
                set_close_nodes = set_close_nodes.union(set_element_nodes)

            # setup for weighted average
            node_ids = list(set_close_nodes)
            structural_nodes = structural_model.getNodeIDLocations(node_ids)
            weights = self.get_weights(close_point, structural_nodes)
            distribution = self.create_distribution(node_ids, weights)

            log.info("element aEID=%s sEID=%s weights=%s" % (
                aero_eid, structural_eid, list_print(weights)))
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

#def pierce_elements(aCentroid, aero_eid, pSource, normal):
    #"""
    #see LoadMapping.pierce_elements
    #"""
    #(sElements, sDists) = centroid_tree.get_close_element_ids(aCentroid)
    #pierced_elements = []

    #for structural_eid, sDist in zip(sElements, sDists):
        #sArea, sNormal, sCentroid = self.structural_model.get_element_properties(structural_eid)
        #sNodes = self.structural_model.get_element_nodes(structural_eid)
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
            #pierced_elements.append([structural_eid, pIntersect, u1, v1, sDist])

    #(pierced_elements, npiercings) = fix_piercings(sElements, pierced_elements)
    #distribution = distribute_unit_load(aero_eid, pierced_elements, npiercings)
    #return distribution

#------------------------------------------------------------------

def run_map_loads(inputs, cart3d_geom='Components.i.triq', bdf_model='fem.bdf',
                  bdf_model_out='fem.loads.out', configpath='', workpath=''):
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
    force_scale = inputs['force_scale']
    moment_scale = inputs['moment_scale']

    if isinstance(bdf_model, BDF):
        pass
    else:
        assert os.path.exists(bdf_model), '%r doesnt exist' % bdf_model


    if aero_format == 'cart3d':
        cart3d_model = Cart3D(debug=True)
        result_names = ['Cp', 'rho', 'rhoU', 'rhoV', 'rhoW', 'E', 'Mach']
        #result_names = None

        if not cart3d_geom.endswith('_half'):
            half_model = cart3d_geom + '_half'
            cart3d_model.read_cart3d(cart3d_geom, result_names=result_names)
            (nodes, elements, regions, loads) = cart3d_model.make_half_model(axis='y')
            cart3d_model.nodes = nodes
            cart3d_model.elements = elements
            cart3d_model.regions = regions
            cart3d_model.loads = loads
            #(nodes, elements, regions, Cp) = cart3d_model.renumber_mesh(
                #nodes, elements, regions, Cp)
            cart3d_model.write_cart3d(half_model)
        else:
            cart3d_model.read_cart3d(cart3d_geom, result_names=result_names)
            elements = cart3d_model.elements
        Cp = cart3d_model.loads['Cp']
        mesh = cart3d_model

    else:
        raise NotImplementedError('aero_format=%r' % aero_format)

    aero_model = AeroModel(inputs, mesh.nodes, mesh.elements, Cp)
    log.info("elements[1] = %s" % elements[1])
    #del elements, nodes, Cp


    if isinstance(bdf_model, BDF):
        fem = bdf_model
    else:
        assert os.path.exists(bdf_model), '%r doesnt exist' % bdf_model
        fem = BDF(debug=True, log=log)
        fem.read_bdf(bdf_model)
        sys.stdout.flush()
    structural_model = StructuralModel(fem, property_regions)

    mapper = LoadMapping(aero_model, structural_model, configpath, workpath)
    t1 = time()
    mapper.set_flight_condition(pInf, qInf)
    mapper.set_reference_quantities(Sref, Lref, xref)
    mapper.set_scale_factors(force_scale, moment_scale)
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
    cart3d_geom = os.path.join(configpath, 'Cart3d_35000_0.825_10_0_0_0_0.i.triq_half')

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
        'pInf' : (pInf, 'psi'),
        'qInf' : (1.4 / 2. * pInf * mach**2., 'psi'),
        'Sref' : (1582876., 'in^2'),
        'Lref' : (623., 'in'),
        'xref' : (268., 'in'),
        'skin_property_regions' : skin_property_regions,
        'force_scale' : (1., 'lb'),
        'moment_scale' : (1/12., 'ft-lb'),
    }
    run_map_loads(inputs, cart3d_geom, bdf_model, bdf_model_out,
                  configpath=configpath, workpath=workpath)


if __name__ == '__main__':
    main()
