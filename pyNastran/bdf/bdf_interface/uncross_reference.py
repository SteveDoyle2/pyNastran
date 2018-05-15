"""
Unlinks up the various cards in the BDF.
"""
from __future__ import print_function
from typing import List, Dict, Any
from six import iteritems, itervalues
from pyNastran.bdf.bdf_interface.safe_cross_reference import SafeXrefMesh

class UnXrefMesh(SafeXrefMesh):
    """
    Unlinks up the various cards in the BDF.
    """
    def __init__(self):
        """
        The main BDF class defines all the parameters that are used.
        """
        SafeXrefMesh.__init__(self)

    def uncross_reference(self):
        """uncross references the model"""
        self.log.debug("Uncross Referencing...")
        self._uncross_reference_nodes()
        self._uncross_reference_coords()
        self._uncross_reference_elements()
        self._uncross_reference_properties()
        self._uncross_reference_materials()
        self._uncross_reference_masses()
        self._uncross_reference_aero()
        self._uncross_reference_constraints()
        self._uncross_reference_loads()
        self._uncross_reference_sets()
        self._uncross_reference_optimization()

    def _uncross_reference_nodes(self):
        # type: () -> None
        """uncross references the GRID objects"""
        for node in itervalues(self.nodes):
            node.uncross_reference()

    def _uncross_reference_coords(self):
        # type: () -> None
        """uncross references the CORDx objects"""
        for cid, coord in iteritems(self.coords):
            if cid == 0:
                continue
            coord.uncross_reference()

    def _uncross_reference_elements(self):
        # type: () -> None
        """uncross references the element objects"""
        for element in itervalues(self.elements):
            try:
                element.uncross_reference()
            except TypeError:
                raise NotImplementedError('%s.uncross_reference' % element.type)
            except AttributeError:
                print(element)
                raise
        for element in itervalues(self.rigid_elements):
            element.uncross_reference()
        for element in itervalues(self.plotels):
            element.uncross_reference()

    def _uncross_reference_properties(self):
        # type: () -> None
        """uncross references the property objects"""
        for prop in itervalues(self.properties):
            try:
                prop.uncross_reference()
            #except TypeError:
                #raise NotImplementedError('%s.uncross_reference' % prop.type)
            except AttributeError:
                print(prop)
                print('%s.uncross_reference error' % prop.type)
                raise

    def _uncross_reference_materials(self):
        # type: () -> None
        """uncross references the material objects"""
        try:
            for material in itervalues(self.materials):
                material.uncross_reference()
        except AttributeError:
            print(material)
            raise

        try:
            for material in itervalues(self.creep_materials):
                material.uncross_reference()
        except AttributeError:
            print(material)
            raise

        data = [self.MATS1, self.MATS3, self.MATS8,
                self.MATT1, self.MATT2, self.MATT3, self.MATT4, self.MATT5,
                self.MATT8, self.MATT9]
        for material_deps in data:
            for mat in itervalues(material_deps):
                try:
                    mat.uncross_reference()
                except AttributeError:
                    print(mat)
                    raise

    def _uncross_reference_masses(self):
        # type: () -> None
        """uncross references the mass objects"""
        for mass in itervalues(self.masses):
            mass.uncross_reference()
        for prop in itervalues(self.properties_mass):
            prop.uncross_reference()

    def _uncross_reference_aero(self):
        # type: () -> None
        """uncross references the aero objects"""
        for caero in itervalues(self.caeros):
            caero.uncross_reference()

        for paero in itervalues(self.paeros):
            paero.uncross_reference()

        for trim in itervalues(self.trims):
            trim.uncross_reference()

        for csschd in itervalues(self.csschds):
            csschd.uncross_reference()

        for spline in itervalues(self.splines):
            spline.uncross_reference()

        for aecomp in itervalues(self.aecomps):
            aecomp.uncross_reference()

        for aelist in itervalues(self.aelists):
            aelist.uncross_reference()

        for aeparam in itervalues(self.aeparams):
            aeparam.uncross_reference()
        for trim in itervalues(self.trims):
            trim.uncross_reference()
        for csschd in itervalues(self.csschds):
            csschd.uncross_reference()

        #for aestat in itervalues(self.aestats):
            #aestat.uncross_reference()

        for aesurf in itervalues(self.aesurf):
            aesurf.uncross_reference()

        for aesurfs in itervalues(self.aesurfs):
            aesurfs.uncross_reference()

        for flutter in itervalues(self.flutters):
            flutter.uncross_reference()

        for monitor_point in self.monitor_points:
            monitor_point.uncross_reference()

        if self.aero:
            self.aero.uncross_reference()
        if self.aeros:
            self.aeros.uncross_reference()

    def _uncross_reference_constraints(self):
        # type: () -> None
        """
        Unlinks the SPCADD, SPC, SPCAX, SPCD, MPCADD, MPC, SUPORT,
        SUPORT1, SESUPORT cards.
        """
        for spcadds in itervalues(self.spcadds):
            for spcadd in spcadds:
                spcadd.uncross_reference()
        for spc in itervalues(self.spcs):
            for spci in spc:
                spci.uncross_reference()
        for spcoffs in itervalues(self.spcoffs):
            for spcoff in spcoffs:
                spcoff.uncross_reference(self)

        for mpcadds in itervalues(self.mpcadds):
            for mpcadd in mpcadds:
                mpcadd.uncross_reference()
        for mpc in itervalues(self.mpcs):
            for mpci in mpc:
                mpci.uncross_reference()
        for suport in self.suport:
            suport.uncross_reference()
        for suport1 in itervalues(self.suport1):
            suport1.uncross_reference()
        for se_suport in self.se_suport:
            se_suport.uncross_reference()

    def _uncross_reference_loads(self):
        # type: () -> None
        """
        Unlinks the LOAD
        PLOAD1, PLOAD2, PLOAD4
        FORCE, FORCE1, FORCE2
        MOMENT, MOMENT1, MOMENT2

        DLOAD, ACSRCE, RLOAD1, RLOAD2, TLOAD1, TLOAD2
        DPHASE, DAREA

        TEMP
        """
        for (unused_lid, load_combinations) in iteritems(self.load_combinations):
            for load_combination in load_combinations:
                load_combination.uncross_reference()

        for (unused_lid, loads) in iteritems(self.loads):
            for load in loads:
                load.uncross_reference()

        for (unused_lid, dloads) in iteritems(self.dloads):
            for dload in dloads:
                dload.uncross_reference()
        for (unused_lid, dload_entries) in iteritems(self.dload_entries):
            for dload_entry in dload_entries:
                dload_entry.uncross_reference()
        for unused_key, darea in iteritems(self.dareas):
            darea.uncross_reference()
        for unused_key, dphase in iteritems(self.dphases):
            dphase.uncross_reference()

    def _uncross_reference_sets(self):
        # type: () -> None
        """uncross references the set objects"""
        for set_obj in self.asets:
            set_obj.uncross_reference()
        for set_obj in self.omits:
            set_obj.uncross_reference()
        for set_obj in self.bsets:
            set_obj.uncross_reference()
        for set_obj in self.csets:
            set_obj.uncross_reference()
        for set_obj in self.qsets:
            set_obj.uncross_reference()
        for unused_name, set_objs in iteritems(self.usets):
            for set_obj in set_objs:
                set_obj.uncross_reference()

        # superelements
        for unused_key, set_obj in iteritems(self.se_sets):
            set_obj.uncross_reference()
        for set_obj in self.se_bsets:
            set_obj.uncross_reference()
        for set_obj in self.se_csets:
            set_obj.uncross_reference()
        for set_obj in self.se_qsets:
            set_obj.uncross_reference()
        for set_obj in self.se_usets:
            set_obj.uncross_reference()

    def _uncross_reference_optimization(self):
        """uncross references the optimization objects"""
        for unused_key, deqatn in iteritems(self.dequations):
            deqatn.uncross_reference()
        for unused_key, dresp in iteritems(self.dresps):
            dresp.uncross_reference()
        for unused_key, dconstrs in iteritems(self.dconstrs):
            for dconstr in dconstrs:
                dconstr.uncross_reference()

        for unused_key, dvcrel in iteritems(self.dvcrels):
            dvcrel.uncross_reference()
        for unused_key, dvmrel in iteritems(self.dvmrels):
            dvmrel.uncross_reference()
        for unused_key, dvprel in iteritems(self.dvprels):
            dvprel.uncross_reference()
