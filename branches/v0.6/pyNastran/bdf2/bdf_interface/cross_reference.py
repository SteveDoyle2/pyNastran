class XRefMesh(object):
    def __init__(self):
        pass
    
    def cross_reference(self, xref=True):
        self.grid.build()
        self.coord.build()

        #self.elements_rod.build()
        self.crod.build()
        self.conrod.build()
        #self.elements_bar.build()
        self.elements_shell.build()
        self.elements_solid.build()

        #self.properties_rod.build()
        self.prod.build()
        #self.properties_bar.build()
        self.properties_shell.build()
        self.properties_solid.build()

        self.materials.build()