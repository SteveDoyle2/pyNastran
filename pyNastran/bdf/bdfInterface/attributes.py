def BDF_Attributes(object):
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

    @nnodes.setter
    def node_ids(self, value):
        raise ValueError("You cannot set node IDs like this...modify the node objects")

    def get_nodes(self):
        nodes = []
        for (nid, node) in sorted(iteritems(self.nodes)):
            nodes.append(node)
        return nodes

    #--------------------
    # Elements CARDS

    @property
    def nelements(self):
        return len(self.elements)

    @nnodes.setter
    def nelements(self, value):
        raise ValueError("You cannot set nnode IDs like this...modify the node objects")

    @property
    def element_ids(self):
        return self.elements.keys()

    @nnodes.setter
    def element_ids(self, value):
        raise ValueError("You cannot set node IDs like this...modify the node objects")

    #--------------------
    # Other CARDS

    @property
    def coord_ids(self):
        return self.coords.keys()

    @nnodes.setter
    def coord_ids(self, value):
        raise ValueError("You cannot set coord IDs like this...modify the coord objects")

    @property
    def ncoord_ids(self):
        return len(self.coords)

    @nnodes.setter
    def ncoord_ids(self, value):
        raise ValueError("You cannot set ncoord IDs like this...modify the coord objects")

    @property
    def ncaero_ids(self):
        return len(self.caeros)

    @nnodes.setter
    def ncaero_ids(self, value):
        raise ValueError("You cannot set ncaero IDs like this...modify the caero objects")

