from __future__ import print_function
from pyNastran.applications.cart3d_nastran_fsi.model import Model

from pyNastran.bdf.field_writer import print_card
from pyNastran.utils.log import get_logger2


class StructuralModel(Model):
    """
    this class gets standard mesh parameters
    """
    def __init__(self, fem, pids, debug=False):
        Model.__init__(self)

        self.log = get_logger2(None, debug=debug)

        self.debug = debug
        if self.debug:
            print("*StructuralModel.init")
        self.fem = fem
        #nodes = fem.getNodes()
        #elements = fem.getElements()
        node_ids = fem.node_ids
        element_ids = fem.element_ids

        self.nnodes = len(node_ids)
        self.nelements = len(element_ids)
        self.pids = pids
        #self.points   = points
        #self.elements = elements
        if self.debug:
            self.log.debug("***StructuralModel.init")

    @property
    def node_ids(self):
        return self.fem.node_ids

    def get_element_properties(self, eid):
        """Returns area, centroid, normal"""
        e = self.fem.Element(eid)
        area, centroid, normal = e.AreaCentroidNormal()
        return area, centroid, normal

    @property
    def element_ids(self):
        #(elements, eids) = getElementsWithPIDs(self, properties)
        eids = self.fem.element_ids
        return eids

    def get_element_node_ids(self, eid):
        e = self.fem.Element(eid)
        return e.node_ids

    def Node(self, nid):
        node = self.fem.Node(nid)
        return node.get_position()

    def Element(self, eid):
        return self.fem.Element(eid)

    def etype(self, eid):
        return self.fem.Element(eid).type

    def getElementIDsWithPIDs(self):
        self.log.info('pids =%s' % self.pids)
        return self.fem.get_element_ids_list_with_pids(self.pids)

    def get_element_nodes(self, eid):
        e = self.fem.Element(eid)
        nodes = []
        for n in e.nodes:
            nodes.append(n.get_position())
        return nodes

    def Centroid(self, eid):
        e = self.fem.Element(eid)
        nodes = self.get_element_nodes(eid)
        centroid = e.Centroid()
        return centroid

    def Centroid_Area(self, eid, nodes):
        e = self.fem.Element(eid)
        raise NotImplementedError('not implemented')
        #return centroid, area

    def Area(self, eid):
        nodes = self.fem.nodes
        e = self.fem.Element(eid)
        return e.Area(nodes)

    def Normal(self, eid):
        nodes = self.fem.nodes
        return self.fem.Normal(eid, nodes)

    def getElements(self): # dict
        return self.elements

    def Properties(self, eid):
        e = self.fem.Element(eid)
        (area, centroid, normal) = e.getAreaCentroidNormal()
        return (normal, centroid)

    def write_load(self, bdf, load_case_id, nid, Fx, Fy, Fz, comment=''):
        """
        This function takes a:
           load case
           node ID
           Force
        and writes out a properly formatted card
        """
        cid = 0
        scale_factor = 1.
        card = ['FORCE', load_case_id, nid, cid, scale_factor, Fx, Fy, Fz]
        #comment += " card=%s" % (card)
        #out = printCard(card)[:-1]+  '   $ %s\n' % comment
        out = print_card(card, size=16)
        bdf.write(out)
