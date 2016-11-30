from __future__ import print_function
from six.moves import range
from pyNastran.applications.cart3d_nastran_fsi.model import Model
from numpy import array, cross, ndarray, isnan

from pyNastran.applications.cart3d_nastran_fsi.math_functions import (
    Triangle_AreaCentroidNormal, ListPrint)

from pyNastran.utils.log import get_logger
debug = True
log = get_logger(None, 'debug' if debug else 'info')

class AeroModel(Model):
    def __init__(self, inputs, nodes, elements, Cp):
        Model.__init__(self)

        #Bref = 623.179569341
        #Cref = 623.179569341
        #Lref = 623.179569341
        #Sref = 1582876.10613
        #Xref = 268.217444245

        self.Mach = inputs['Mach']
        self.pInf = inputs['pInf']
        self.qInf = inputs['qInf']
        self.Sref = inputs['Sref']
        self.Lref = inputs['Lref']
        self.xref = inputs['xref']

        self._nNodes = len(nodes)
        self._nElements = len(elements)
        self.nodes = nodes
        self.elements = elements
        Cp = self.prepare_Cps(Cp)  # convert nodal Cp to centroidal Cp

        self.prepare_centroid_area_normals()
        self.get_moments(Cp) # centroidal Cp

    def get_moments(self, Cp):
        moment_center = array([self.xref, 0., 0.])
        sum_forces = array([0., 0., 0.])
        sum_moments = array([0., 0., 0.])

        for (key, cp) in self.Cps.items(): # centroidal based Cp
            area = self.areas[key]
            centroid = self.centroids[key]
            normal = self.normals[key]
            p = cp * self.qInf + self.pInf
            F = area * normal * p  # negative sign is b/c the normals are flipped...
            r = moment_center - centroid
            if any(isnan(r)):
                msg = 'r=%s moment_center=%s centroid=%s' % (r, moment_center, centroid)
                raise RuntimeError(msg)
            if any(isnan(F)):
                msg = 'area=%s normal=%s p=%s' % (area, normal, p)
                raise RuntimeError(msg)

            sum_forces += F
            sum_moments += cross(r, F)
        log.info("pInf=%s [psi]; qInf= %s [psi]" % (self.pInf, self.qInf))

        log.info("sumForcesCFD  [lb]    = %s" % ListPrint(sum_forces))
        log.info("sumMomentsCFD [ft-lb] = %s" % ListPrint(sum_moments/12.))
        Cf = sum_forces  / (self.Sref * self.qInf)
        Cm = sum_moments / (self.Sref * self.qInf * self.Lref) * 12.
        log.info("Cf = %s" % ListPrint(Cf))
        log.info("Cm = %s" % ListPrint(Cm))
        return (sum_forces, sum_moments/12.)

    def prepare_centroid_area_normals(self):
        self.centroids = {}
        self.areas = {}
        self.normals = {}
        for eid, element in enumerate(self.elements):
            n1, n2, n3 = element
            n1 = self.nodes[n1 - 1]
            n2 = self.nodes[n2 - 1]
            n3 = self.nodes[n3 - 1]
            area, centroid, normal = Triangle_AreaCentroidNormal([n1, n2, n3])

            eidi = eid + 1
            if centroid is not None:
                assert len(centroid) == 3, "eid=%s centroid=%s n1=%s n2=%s n3=%s" % (eidi, centroid, n1, n2, n3)

            self.areas[eidi] = area
            self.centroids[eidi] = centroid
            self.normals[eidi] = normal

    def prepare_Cps(self, Cp):
        """
        converts Cp applied to the node -> Cp applied on the element centroid
        """
        #self.Cps = Cp
        #Cp = loads['Cp']
        assert isinstance(Cp, ndarray)

        Cp_dict = {}
        for eid, element in enumerate(self.elements):
            eidi = eid + 1
            (n1, n2, n3) = element #self.get_element_node_ids(eid)
            cp = Cp[element-1].sum() / 3.
            Cp_dict[eidi] = cp
        self.Cps = Cp_dict
        return Cp_dict

    def nNodes(self):
        return self._nNodes

    def Centroid(self, eid):
        return self.centroids[eid]
        #nodes = self.get_element_nodes(eid)
        #return Centroid(*nodes)

    def nElements(self):
        return self._nElements

    def nCps(self):
        return self._nElements

    def Node(self, nid):
        return self.nodes[nid-1]

    def Element(self, eid):
        try:
            element = self.elements[eid-1]
        except IndexError:
            print("eid=%s len(elements)=%s" % (eid, len(self.elements)))
            raise
        except TypeError:
            print("eid=%s len(elements)=%s" % (eid, len(self.elements)))
            raise
        except KeyError:
            print("eid=%s len(elements)=%s" % (eid, len(self.elements)))
            raise
        return element

    def Area(self, eid):
        return self.areas[eid]

        #nodes = self.get_element_nodes(eid)
        ##print("nodes[%s]=%s" % (eid, nodes))
        #area, normal = AreaNormal(nodes)
        #return area

    def Normal(self, eid):
        return self.normals[eid]

        #nodes = self.get_element_nodes(eid)
        #area, normal = AreaNormal(nodes)
        #return normal

    def area_normal_centroid(self, eid):
        area = self.areas[eid]
        normal = self.normals[eid]
        centroid = self.centroids[eid]

        #nodes = self.get_element_nodes(eid)
        #area, normal, centroid = Triangle_AreaNormalCentroid(nodes)
        return area, normal, centroid

    def NodeIDs(self):
        nNodes = self.nNodes()
        return range(1, nNodes + 1)

    #def getElementIDsWithPIDs(self):
        #return self.ElementIDs()

    def ElementIDs(self):
        nElements = self.nElements()
        return range(1, nElements + 1)

    def Cp(self, eid):
        cp = self.Cps[eid]
        return cp

    def get_element_node_ids(self, eid):
        e = self.Element(eid)
        nid1, nid2, nid3 = e
        return nid1, nid2, nid3

    def get_element_nodes(self, eid):
        nid1, nid2, nid3 = self.get_element_node_ids(eid)
        n1 = self.Node(nid1)
        n2 = self.Node(nid2)
        n3 = self.Node(nid3)
        #print("nids[%s] = %s %s %s" %(eid, nid1, nid2, nid3))
        #print("n[%s]=%s n[%s]=%s n[%s]=%s\n" %(nid1, n1, nid2, n2, nid3, n3))
        return n1, n2, n3

    def get_element_properties(self, eid):
        """
        Returns area, centroid, normal
        """
        nodes = self.get_element_nodes(eid)
        (area, centroid, normal) = Triangle_AreaCentroidNormal(nodes)
        return (area, centroid, normal)
