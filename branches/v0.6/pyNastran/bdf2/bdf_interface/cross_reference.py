from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
from ..cards.loads.loadcase import LoadCase


class XRefMesh(object):
    def __init__(self):
        pass
    
    def cross_reference(self, xref=True):
        self.grid.build()
        self.coord.build()
        
        self._build_elements_properties()
    
        self.materials.build()
        
        self._build_loads()
        self._build_constraints()
        self._build_aero()

    def _build_aero(self):
        self.caero.build()
        self.paero.build()
        self.mass.build()
        #self.trim.build()
        #self.aero.build()
        #self.aeros.build()

    def _build_elements_properties(self):
        #self.elements_rod.build()
        self.crod.build()
        self.conrod.build()
        
        self.elements_spring.build()
        #self.elements_bar.build()
        self.cbar.build()

        self.elements_shell.build()
        self.elements_solid.build()

        self.pelas.build()
        #self.properties_rod.build()
        self.prod.build()
        #self.properties_bar.build()
        self.properties_shell.build()
        self.properties_solid.build()
        
        self.mass.build()

    def _build_constraints(self):
        for t in [self.spcadd, self.spc, self.spc1, self.spcd]:
            for constraint_id, constraint in t.iteritems():
                constraint.build()

        for t in [self.mpcadd, self.mpc]:
            for constraint_id, constraint in t.iteritems():
                constraint.build()

    def _build_loads(self):
        #self.loadcase.build()
        for load_id, loads in self.load.iteritems():
            for load in loads:
                load.build()
        for load_id, loads in self.dload.iteritems():
            for load in loads:
                load.build()

        #self.loadset.build()
        self.force.build()
        #self.force1.build()
        #self.force2.build()
        self.moment.build()
        #self.moment1.build()
        #self.moment2.build()
        
        
        self.pload.build()
        self.pload1.build()
        self.pload2.build()
        #self.pload4.build()

        self.ploadx1.build()
        
        if 0:
            self.spc_object = SPCObject()
            self.spc_object.add_reference(self.spcadd)
            self.spc_object.add(self.spc)
            self.spc_object.add(self.spcd)
            self.spc_object.add(self.spc1)
            self.spc_object.add(self.spcax)

            self.mpc_object = MPCObject()
            self.mpc_object.add_reference(self.mpcadd)
            self.mpc_object.add(self.mpc)
        
        self.loadcase = LoadCase()
        self.loadcase.add_reference(self.load)
        self.loadcase.add_reference(self.dload)
        self.loadcase.add(self.force)
        #self.loadcase.add(self.force1)
        #self.loadcase.add(self.force2)

        self.loadcase.add(self.moment)
        #self.loadcase.add(self.moment1)
        #self.loadcase.add(self.moment2)
        
        self.loadcase.add(self.pload)
        self.loadcase.add(self.pload1)
        self.loadcase.add(self.pload2)
        #self.loadcase.add(self.pload4)
        
        self.loadcase.add(self.ploadx1)
        
        #self.loadcase.resolve(2)
        #self.loadcase.resolve(1)