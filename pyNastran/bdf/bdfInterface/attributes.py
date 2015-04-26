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
    # Property CARDS

    @property
    def nproperties(self):
        return len(self.properties)

    @nproperties.setter
    def nproperties(self, value):
        raise ValueError("You cannot set nproperties IDs like this...modify the properties objects")

    @property
    def property_ids(self):
        return self.properties.keys()

    @property_ids.setter
    def property_ids(self, value):
        raise ValueError("You cannot set property IDs like this...modify the properties objects")
    #--------------------
    # Material CARDS

    @property
    def material_ids(self):
        return self.materials.keys()

    @material_ids.setter
    def material_ids(self, value):
        raise ValueError("You cannot set material IDs like this...modify the materials objects")

    @property
    def nmaterials(self):
        return len(self.materials)

    @nmaterials.setter
    def nmaterials(self, value):
        raise ValueError("You cannot set nmaterials IDs like this...modify the materials objects")

    #--------------------
    # Coords CARDS

    @property
    def coord_ids(self):
        return self.coords.keys()

    @coord_ids.setter
    def coord_ids(self, value):
        raise ValueError("You cannot set coord IDs like this...modify the coord objects")

    @property
    def ncoords(self):
        return len(self.coords)

    @ncoords.setter
    def ncoords(self, value):
        raise ValueError("You cannot set ncoord IDs like this...modify the coord objects")

    #--------------------

    @property
    def ncaeros(self):
        return len(self.caeros)

    @ncaeros.setter
    def ncaeros(self, value):
        raise ValueError("You cannot set ncaero IDs like this...modify the caero objects")

