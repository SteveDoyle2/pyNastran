"""wrapper around Cart3d"""
from __future__ import print_function
from six.moves import range
from numpy import array, cross, ndarray, isnan

from pyNastran.applications.cart3d_nastran_fsi.model import Model
from pyNastran.applications.cart3d_nastran_fsi.math_functions import (
    triangle_area_centroid_normal, list_print)

from pyNastran.utils.log import get_logger2
debug = True
log = get_logger2(None, debug=debug)

class AeroModel(Model):
    """wrapper around Cart3d"""
    def __init__(self, inputs, nodes, elements, Cp):
        Model.__init__(self)
        self.centroids = {}
        self.areas = {}
        self.normals = {}

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
        self.force_scale = inputs['force_scale'] # 1 for lb -> lb
        self.moment_scale = inputs['moment_scale']  # 1/12 for in-lb to ft-lb

        self._nnodes = len(nodes)
        self._nelements = len(elements)
        self.nodes = nodes
        self.elements = elements
        Cp = self.prepare_Cps(Cp)  # convert nodal Cp to centroidal Cp

        self.prepare_centroid_area_normals()
        self.get_moments(Cp) # centroidal Cp

    def get_moments(self, Cp):
        qinf, qinf_unit = self.qInf
        pinf, pinf_unit = self.pInf
        sref, _sref_unit = self.Sref
        lref, _lref_unit = self.Lref
        xref, _xref_unit = self.xref
        force_scale, force_scale_unit = self.force_scale
        moment_scale, moment_scale_unit = self.moment_scale

        moment_center = array([xref, 0., 0.])
        sum_forces = array([0., 0., 0.])
        sum_moments = array([0., 0., 0.])

        for (key, cp) in self.Cps.items(): # centroidal based Cp
            area = self.areas[key]
            centroid = self.centroids[key]
            normal = self.normals[key]
            p = cp * qinf + pinf
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
        log.info("pInf=%s [%s]; qInf=%s [%s]" % (pinf, pinf_unit,
                                                 qinf, qinf_unit))
        log.info("force_scale=%s %s" % (force_scale, force_scale_unit))
        log.info("moment_scale=%s %s" % (moment_scale, moment_scale_unit))

        log.info("sumForcesCFD  [%s] = %s" % (
            force_scale_unit, list_print(sum_forces * force_scale)))
        log.info("sumMomentsCFD [%s] = %s" % (
            moment_scale_unit, list_print(sum_moments * moment_scale)))

        Cf = sum_forces / (sref * qinf)
        Cm = sum_moments / (sref * qinf * lref)
        log.info("Cf = %s" % list_print(Cf))
        log.info("Cm = %s" % list_print(Cm))
        return (sum_forces * force_scale, sum_moments * moment_scale)

    def prepare_centroid_area_normals(self):
        for eid, element in enumerate(self.elements):
            n1, n2, n3 = element
            n1 = self.nodes[n1 - 1]
            n2 = self.nodes[n2 - 1]
            n3 = self.nodes[n3 - 1]
            area, centroid, normal = triangle_area_centroid_normal([n1, n2, n3])

            eidi = eid + 1
            if centroid is not None:
                if len(centroid) != 3:
                    msg = "eid=%s centroid=%s n1=%s n2=%s n3=%s" % (eidi, centroid, n1, n2, n3)
                    raise RuntimeError(msg)

            self.areas[eidi] = area
            self.centroids[eidi] = centroid
            self.normals[eidi] = normal

    def prepare_Cps(self, Cp):
        """
        converts Cp applied to the node -> Cp applied on the element centroid
        """
        #self.Cps = Cp
        #Cp = loads['Cp']
        assert isinstance(Cp, ndarray), Cp

        Cp_dict = {}
        for eid, element in enumerate(self.elements):
            eidi = eid + 1
            (n1, n2, n3) = element #self.get_element_node_ids(eid)
            cp = Cp[element-1].sum() / 3.
            Cp_dict[eidi] = cp
        self.Cps = Cp_dict
        return Cp_dict

    @property
    def nnodes(self):
        return self._nnodes

    def Centroid(self, eid):
        return self.centroids[eid]
        #nodes = self.get_element_nodes(eid)
        #return Centroid(*nodes)

    @property
    def nelements(self):
        return self._nelements

    def nCps(self):
        return self._nelements

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

    #@property
    #def node_ids(self):
        #nnodes = self.nnodes
        #return range(1, nnodes + 1)

    #def getElementIDsWithPIDs(self):
        #return self.ElementIDs()

    @property
    def element_ids(self):
        nelements = self.nelements
        return range(1, nelements + 1)

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
        (area, centroid, normal) = triangle_area_centroid_normal(nodes)
        return (area, centroid, normal)
