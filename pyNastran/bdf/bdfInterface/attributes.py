class BDFAttributes(object):

    def __init__(self):
        self._nastran_format = 'msc'

    def set_as_msc(self):
        self._nastran_format = 'msc'

    def set_as_nx(self):
        self._nastran_format = 'nx'

    @property
    def nastran_format(self):
        return self._nastran_format

    @nastran_format.setter
    def nastran_format(self, nastran_format):
        fmt_lower = nastran_format.lower().strip()
        if fmt_lower not in ['nx', 'msc']:
            raise RuntimeError(nastran_format)
        self._nastran_format = fmt_lower

    @property
    def is_long_ids(self):
        if self._nastran_format == 'nx':
            return True
        return False

    @property
    def caseControlDeck(self):
        self.deprecated('self.caseControlDeck', 'self.case_control_deck', '0.8')
        return self.case_control_deck

    @caseControlDeck.setter
    def caseControlDeck(self, value):
        self.deprecated('self.caseControlDeck', 'self.case_control_deck', '0.8')
        self.case_control_deck = value

    @property
    def subcases(self):
        if self.case_control_deck is None:
            return {}
        return self.case_control_deck.subcases

    @property
    def nnodes(self):
        return len(self.nodes)

    @property
    def node_ids(self):
        return self.nodes.keys()

    @property
    def point_ids(self):
        if self.spoints is not None and self.epoints is not None:
            return set(self.node_ids) | self.spoints.points | self.epoints.points
        elif self.spoints is not None:
            return set(self.node_ids) | self.spoints.points
        elif self.epoints is not None:
            return set(self.node_ids) | self.epoints.points
        return set(self.node_ids)

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

    @property
    def element_ids(self):
        return self.elements.keys()

    #--------------------
    # Property CARDS

    @property
    def nproperties(self):
        return len(self.properties)

    @property
    def property_ids(self):
        return self.properties.keys()

    #--------------------
    # Material CARDS

    @property
    def material_ids(self):
        return self.materials.keys()

    @property
    def nmaterials(self):
        return len(self.materials)

    #--------------------
    # Coords CARDS

    @property
    def coord_ids(self):
        return self.coords.keys()

    @property
    def ncoords(self):
        return len(self.coords)

    #--------------------

    @property
    def ncaeros(self):
        return len(self.caeros)

    @property
    def caero_ids(self):
        return self.caeros.keys()
