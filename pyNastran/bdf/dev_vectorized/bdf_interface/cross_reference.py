from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
from six import iteritems
from pyNastran.bdf.dev_vectorized.cards.loads.loadcase import LoadCase


class XRefMesh(object):
    def __init__(self):
        pass

    def build(self, xref=True):
        if xref:
            self.grid.build()
            self.point.build()
            self.spoint.build()
            self.epoint.build()
            self.pointax.build()
            #======================
            self.coords.build()
            #======================
            self._build_elements_properties()
            #======================
            self.materials.build()
            #======================

            self._build_loads()
            self._build_constraints()
            self._build_aero()

    def _build_aero(self):
        self.caero.build()
        self.paero.build()
        self.spline1.build()
        #self.trim.build()
        self.aero.build()
        self.aeros.build()

    def _build_elements_properties(self):
        #self.elements_rod.build()

        if 1:
            self.elements.build()
        else:
            self.crod.build()
            self.conrod.build()
            self.ctube.build()

            self.elements_spring.build()

            #self.elements_bar.build()
            self.cbar.build()
            self.properties_bar.build()

            #self.elements_beam.build()
            self.cbeam.build()
            self.properties_beam.build()

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
            for constraint_id, constraint in iteritems(t):
                constraint.build()

        for t in [self.mpcadd, self.mpc]:
            for constraint_id, constraint in iteritems(t):
                constraint.build()

    def _build_loads(self):
        #self.loadcase.build()
        for load_id, loads in iteritems(self.loads.load):
            for load in loads:
                load.build()
        for load_id, loads in iteritems(self.loads.dload):
            for load in loads:
                load.build()

        #self.loadset.build()
        self.loads.force.build()
        #self.loads.force1.build()
        #self.loads.force2.build()
        self.loads.moment.build()
        #self.loads.moment1.build()
        #self.loads.moment2.build()


        self.loads.pload.build()
        self.loads.pload1.build()
        self.loads.pload2.build()
        #self.loads.pload4.build()

        self.loads.ploadx1.build()

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
