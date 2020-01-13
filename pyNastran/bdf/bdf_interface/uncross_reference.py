"""Unlinks up the various cards in the BDF."""
from typing import List, Dict, Any
from pyNastran.bdf.bdf_interface.safe_cross_reference import SafeXrefMesh

class UnXrefMesh(SafeXrefMesh):
    """
    Unlinks up the various cards in the BDF.
    """
    def __init__(self) -> None:
        """
        The main BDF class defines all the parameters that are used.
        """
        SafeXrefMesh.__init__(self)

    def uncross_reference(self, word: str='') -> None:
        """uncross references the model"""
        self.log.debug("Uncross Referencing%s..." % word)
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
        self._uncross_reference_contact()
        self._uncross_reference_superelements()

        for super_id, superelement in sorted(self.superelement_models.items()):
            superelement.uncross_reference(word=' (Superelement %i)' % super_id)

    def _uncross_reference_nodes(self) -> None:
        """uncross references the GRID objects"""
        for node in self.nodes.values():
            node.uncross_reference()
        for point in self.points.values():
            point.uncross_reference()

    def _uncross_reference_coords(self) -> None:
        """uncross references the CORDx objects"""
        for cid, coord in self.coords.items():
            if cid == 0:
                continue
            coord.uncross_reference()

    def _uncross_reference_elements(self) -> None:
        """uncross references the element objects"""
        for element in self.elements.values():
            try:
                element.uncross_reference()
            except TypeError:
                raise NotImplementedError('%s.uncross_reference' % element.type)
            except AttributeError:
                print(element)
                raise
        for element in self.masses.values():
            element.uncross_reference()
        for element in self.rigid_elements.values():
            element.uncross_reference()
        for element in self.plotels.values():
            element.uncross_reference()

    def _uncross_reference_properties(self) -> None:
        """uncross references the property objects"""
        for prop in self.properties.values():
            try:
                prop.uncross_reference()
            #except TypeError:
                #raise NotImplementedError('%s.uncross_reference' % prop.type)
            except AttributeError:
                print(prop)
                print('%s.uncross_reference error' % prop.type)
                raise

    def _uncross_reference_materials(self) -> None:
        """uncross references the material objects"""
        try:
            for material in self.materials.values():
                material.uncross_reference()
        except AttributeError:
            print(material)
            raise

        try:
            for material in self.creep_materials.values():
                material.uncross_reference()
        except AttributeError:
            print(material)
            raise

        data = [self.MATS1, self.MATS3, self.MATS8,
                self.MATT1, self.MATT2, self.MATT3, self.MATT4, self.MATT5,
                self.MATT8, self.MATT9]
        for material_deps in data:
            for mat in material_deps.values():
                try:
                    mat.uncross_reference()
                except AttributeError:
                    print(mat)
                    raise

    def _uncross_reference_masses(self) -> None:
        """uncross references the mass objects"""
        for mass in self.masses.values():
            mass.uncross_reference()
        for prop in self.properties_mass.values():
            prop.uncross_reference()

    def _uncross_reference_aero(self) -> None:
        """uncross references the aero objects"""
        for caero in self.caeros.values():
            caero.uncross_reference()

        for paero in self.paeros.values():
            paero.uncross_reference()

        for trim in self.trims.values():
            trim.uncross_reference()

        for csschd in self.csschds.values():
            csschd.uncross_reference()

        for spline in self.splines.values():
            spline.uncross_reference()

        for aecomp in self.aecomps.values():
            aecomp.uncross_reference()

        for aelist in self.aelists.values():
            aelist.uncross_reference()

        for aeparam in self.aeparams.values():
            aeparam.uncross_reference()
        for trim in self.trims.values():
            trim.uncross_reference()
        for csschd in self.csschds.values():
            csschd.uncross_reference()

        #for aestat in self.aestats.values():
            #aestat.uncross_reference()

        for aesurf in self.aesurf.values():
            aesurf.uncross_reference()

        for aesurfs in self.aesurfs.values():
            aesurfs.uncross_reference()

        for flutter in self.flutters.values():
            flutter.uncross_reference()

        for monitor_point in self.monitor_points:
            monitor_point.uncross_reference()

        if self.aero:
            self.aero.uncross_reference()
        if self.aeros:
            self.aeros.uncross_reference()

    def _uncross_reference_constraints(self) -> None:
        """
        Unlinks the SPCADD, SPC, SPCAX, SPCD, MPCADD, MPC, SUPORT,
        SUPORT1, SESUPORT cards.
        """
        for spcadds in self.spcadds.values():
            for spcadd in spcadds:
                spcadd.uncross_reference()
        for spc in self.spcs.values():
            for spci in spc:
                spci.uncross_reference()
        for spcoffs in self.spcoffs.values():
            for spcoff in spcoffs:
                spcoff.uncross_reference()

        for mpcadds in self.mpcadds.values():
            for mpcadd in mpcadds:
                mpcadd.uncross_reference()
        for mpc in self.mpcs.values():
            for mpci in mpc:
                mpci.uncross_reference()
        for suport in self.suport:
            suport.uncross_reference()
        for suport1 in self.suport1.values():
            suport1.uncross_reference()
        for se_suport in self.se_suport:
            se_suport.uncross_reference()

    def _uncross_reference_loads(self) -> None:
        """
        Unlinks the LOAD
        PLOAD1, PLOAD2, PLOAD4
        FORCE, FORCE1, FORCE2
        MOMENT, MOMENT1, MOMENT2

        DLOAD, ACSRCE, RLOAD1, RLOAD2, TLOAD1, TLOAD2
        DPHASE, DAREA

        TEMP
        """
        for (unused_lid, load_combinations) in self.load_combinations.items():
            for load_combination in load_combinations:
                load_combination.uncross_reference()

        for (unused_lid, loads) in self.loads.items():
            for load in loads:
                load.uncross_reference()

        for (unused_lid, dloads) in self.dloads.items():
            for dload in dloads:
                dload.uncross_reference()
        for (unused_lid, dload_entries) in self.dload_entries.items():
            for dload_entry in dload_entries:
                dload_entry.uncross_reference()
        for unused_key, darea in self.dareas.items():
            darea.uncross_reference()
        for unused_key, dphase in self.dphases.items():
            dphase.uncross_reference()
        for unused_key, tic in self.tics.items():
            tic.uncross_reference()

    def _uncross_reference_sets(self) -> None:
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
        for unused_name, set_objs in self.usets.items():
            for set_obj in set_objs:
                set_obj.uncross_reference()

        # superelements
        for unused_key, set_obj in self.se_sets.items():
            set_obj.uncross_reference()
        for set_obj in self.se_bsets:
            set_obj.uncross_reference()
        for set_obj in self.se_csets:
            set_obj.uncross_reference()
        for set_obj in self.se_qsets:
            set_obj.uncross_reference()
        for set_obj in self.se_usets:
            set_obj.uncross_reference()

    def _uncross_reference_optimization(self) -> None:
        """uncross references the optimization objects"""
        for unused_key, deqatn in self.dequations.items():
            deqatn.uncross_reference()
        for unused_key, dresp in self.dresps.items():
            dresp.uncross_reference()
        for unused_key, dconstrs in self.dconstrs.items():
            for dconstr in dconstrs:
                dconstr.uncross_reference()

        for unused_key, dvcrel in self.dvcrels.items():
            dvcrel.uncross_reference()
        for unused_key, dvmrel in self.dvmrels.items():
            dvmrel.uncross_reference()
        for unused_key, dvprel in self.dvprels.items():
            dvprel.uncross_reference()
        for unused_key, desvar in self.desvars.items():
            desvar.uncross_reference()
        for unused_key, topvar in self.topvar.items():
            topvar.uncross_reference()
