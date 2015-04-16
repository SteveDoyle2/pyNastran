class BDFAttributes(object):

    def __init__(self):
        pass

    @property
    def nnodes(self):
        return len(self.nodes)

    @nnodes.setter
    def nnodes(self, value):
        raise ValueError("You cannot set node IDs like this...modify the node objects")

    @property
    def node_ids(self):
        return self.nodes.keys()

    @node_ids.setter
    def node_ids(self, value):
        raise ValueError("You cannot set node IDs like this...modify the node objects")

    #def get_nodes(self):
        #nodes = []
        #for nid, node in sorted(iteritems(self.nodes)):
            #nodes.append(node)
        #return nodes

    #--------------------
    # Elements CARDS

    @property
    def nelements(self):
        return len(self.elements)

    @nelements.setter
    def nelements(self, value):
        raise ValueError("You cannot set nnode IDs like this...modify the node objects")

    @property
    def element_ids(self):
        return self.elements.keys()

    @element_ids.setter
    def element_ids(self, value):
        raise ValueError("You cannot set node IDs like this...modify the node objects")

    #--------------------
    # Other CARDS

    @property
    def coord_ids(self):
        return self.coords.keys()

    @coord_ids.setter
    def coord_ids(self, value):
        raise ValueError("You cannot set coord IDs like this...modify the coord objects")

    @property
    def ncoord_ids(self):
        return len(self.coords)

    @ncoord_ids.setter
    def ncoord_ids(self, value):
        raise ValueError("You cannot set ncoord IDs like this...modify the coord objects")

    @property
    def ncaero_ids(self):
        return len(self.caeros)

    @ncaero_ids.setter
    def ncaero_ids(self, value):
        raise ValueError("You cannot set ncaero IDs like this...modify the caero objects")

