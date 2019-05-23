from pyNastran.dev.bdf_vectorized.bdf_interface2.attributes import BDFAttributes
from pyNastran.dev.bdf_vectorized.cards.loads.loadcase import LoadCase


class CrossReference(BDFAttributes):
    """defines methods for writing cards"""

    def __init___(self):
        """creates methods for writing cards"""
        BDFAttributes.__init__(self)

    def cross_reference(self, xref=True):
        self.build_loadcase()

    def uncross_reference(self) -> None:
        pass

    def build_loadcase(self):
        """builds the loadcase object"""
        self.loadcase = LoadCase(self)
        self.loadcase.add_reference(self.loads.load)
        self.loadcase.add_reference(self.loads.dload)
        #self.loadcase.add_reference(self.loads.sload)
        #self.loadcase.add_reference(self.loads.lseq)

        self.loadcase.add(self.loads.force)
        #self.loadcase.add(self.loads.force1)
        #self.loadcase.add(self.loads.force2)

        self.loadcase.add(self.loads.moment)
        #self.loadcase.add(self.loads.moment1)
        #self.loadcase.add(self.loads.moment2)

        self.loadcase.add(self.loads.pload)
        self.loadcase.add(self.loads.pload1)
        self.loadcase.add(self.loads.pload2)
        #self.loadcase.add(self.loads.pload3)
        #self.loadcase.add(self.loads.pload4)

        self.loadcase.add(self.loads.ploadx1)
        self.loadcase.add(self.loads.grav)
        self.loadcase.add(self.loads.rforce)

        #self.loadcase.add(self.loads.tload1)
        #self.loadcase.add(self.loads.tload2)
        #self.loadcase.add(self.loads.rload1)
        #self.loadcase.add(self.loads.rload2)

        #self.loadcase.add(self.loads.accel1)
        #self.loadcase.add(self.loads.randps)


        #self.loadcase.resolve(2)
        #self.loadcase.resolve(1)

        # DAREA
        # RANDPS

