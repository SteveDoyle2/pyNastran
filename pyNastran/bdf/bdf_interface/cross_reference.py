"""
Links up the various cards in the BDF.

For example, with cross referencing...

.. code-block:: python

   >>> model = BDF()
   >>> model.read_bdf(bdf_filename, xref=True)

   >>> nid1 = 1
   >>> node1 = model.nodes[nid1]
   >>> node.nid
   1

   >>> node.xyz
   [1., 2., 3.]

   >>> node.Cid()
   3

   >>> node.cid
   3

   >>> node.cid_ref
   CORD2S, 3, 1, 0., 0., 0., 0., 0., 1.,
           1., 0., 0.
   # get the position in the global frame
   >>> node.get_position()
   [4., 5., 6.]

   # get the position with respect to another frame
   >>> node.get_position_wrt(model, cid=2)
   [4., 5., 6.]


Without cross referencing...

.. code-block:: python

   >>> model = BDF()
   >>> model.read_bdf(bdf_filename, xref=True)

   >>> nid1 = 1
   >>> node1 = model.nodes[nid1]
   >>> node.nid
   1

   >>> node.xyz
   [1., 2., 3.]

   >>> node.Cid()
   3

   >>> node.cid
   3

   >>> node.cid_ref
   None

   # get the position in the global frame
   >>> node.get_position()
   Error!

Cross-referencing allows you to easily jump across cards and also helps
with calculating things like position, area, and mass.  The BDF is designed
around the idea of cross-referencing, so it's recommended that you use it.

"""
# pylint: disable=R0902,R0904,R0914
from collections import defaultdict
import traceback
from typing import List, Dict, Any

from numpy import zeros, argsort, arange, array_equal, array
from pyNastran.bdf.bdf_interface.attributes import BDFAttributes

class XrefMesh(BDFAttributes):
    """Links up the various cards in the BDF."""
    def __init__(self) -> None:
        """The main BDF class defines all the parameters that are used."""
        BDFAttributes.__init__(self)
        self._nxref_errors = 100
        self._stop_on_xref_error = True

    # def geom_check(self):
        # """
        # Performs various geometry checks
          # 1.  nodal uniqueness on elements
        # """
        # for elem in model.elements:
            # elem.check_unique_nodes()

    def cross_reference(self,
                        xref: bool=True,
                        xref_nodes: bool=True,
                        xref_elements: bool=True,
                        xref_nodes_with_elements: bool=False,
                        xref_properties: bool=True,
                        xref_masses: bool=True,
                        xref_materials: bool=True,
                        xref_loads: bool=True,
                        xref_constraints: bool=True,
                        xref_aero: bool=True,
                        xref_sets: bool=True,
                        xref_optimization: bool=True,
                        word: str='') -> None:
        """
        Links up all the cards to the cards they reference

        Parameters
        ----------
        xref : bool; default=True
           cross references the model
        xref_nodes : bool; default=True
           set cross referencing of nodes/coords
        xref_element : bool; default=True
           set cross referencing of elements
        xref_properties : bool; default=True
           set cross referencing of properties
        xref_masses : bool; default=True
           set cross referencing of CMASS/PMASS
        xref_materials : bool; default=True
           set cross referencing of materials
        xref_loads : bool; default=True
            set cross referencing of loads
        xref_constraints : bool; default=True
            set cross referencing of constraints
        xref_aero : bool; default=True
            set cross referencing of CAERO/SPLINEs
        xref_sets : bool; default=True
            set cross referencing of SETx
        word : str; default=''
            model flag

        To only cross-reference nodes:

        .. code-block:: python

           model = BDF()
           model.read_bdf(bdf_filename, xref=False)
           model.cross_reference(xref=True, xref_loads=False, xref_constraints=False,
                                            xref_materials=False, xref_properties=False,
                                            xref_aero=False, xref_masses=False,
                                            xref_sets=False)

        .. warning:: be careful if you call this method with False values
        """
        if not xref:
            return
        self.log.debug("Cross Referencing%s..." % word)
        if xref_nodes:
            self._cross_reference_nodes()
            self._cross_reference_coordinates()

        if xref_elements:
            self._cross_reference_elements()
        if xref_properties:
            self._cross_reference_properties()
        if xref_masses:
            self._cross_reference_masses()
        if xref_materials:
            self._cross_reference_materials()

        if xref_aero:
            self._cross_reference_aero()
        if xref_constraints:
            self._cross_reference_constraints()
        if xref_loads:
            self._cross_reference_loads()
        if xref_sets:
            self._cross_reference_sets()
        if xref_optimization:
            self._cross_reference_optimization()
        if xref_nodes_with_elements:
            self._cross_reference_nodes_with_elements()
        self._cross_reference_contact()
        self._cross_reference_superelements()
        #self.case_control_deck.cross_reference(self)
        self.pop_xref_errors()

        for super_id, superelement in sorted(self.superelement_models.items()):
            superelement.cross_reference(
                xref=xref, xref_nodes=xref_nodes, xref_elements=xref_elements,
                xref_nodes_with_elements=xref_nodes_with_elements,
                xref_properties=xref_properties, xref_masses=xref_masses,
                xref_materials=xref_materials, xref_loads=xref_loads,
                xref_constraints=xref_constraints, xref_aero=xref_aero,
                xref_sets=xref_sets, xref_optimization=xref_optimization,
                word=' (Superelement %i)' % super_id)

    def _cross_reference_constraints(self) -> None:
        """
        Links the SPCADD, SPC, SPCAX, SPCD, MPCADD, MPC, SUPORT,
        SUPORT1, SESUPORT cards.
        """
        for spcadds in self.spcadds.values():
            for spcadd in spcadds:
                spcadd.cross_reference(self)
        for spcs in self.spcs.values():
            for spc in spcs:
                spc.cross_reference(self)
        for spcoffs in self.spcoffs.values():
            for spcoff in spcoffs:
                spcoff.cross_reference(self)

        for mpcadds in self.mpcadds.values():
            for mpcadd in mpcadds:
                mpcadd.cross_reference(self)
        for mpcs in self.mpcs.values():
            for mpc in mpcs:
                mpc.cross_reference(self)

        for suport in self.suport:
            suport.cross_reference(self)

        for unused_suport1_id, suport1 in self.suport1.items():
            suport1.cross_reference(self)

        for se_suport in self.se_suport:
            se_suport.cross_reference(self)

    def _cross_reference_coordinates(self) -> None:
        """
        Links up all the coordinate cards to other coordinate cards and nodes
         - CORD1R, CORD1C, CORD1S
         - CORD2R, CORD2C, CORD2S
        """
        # CORD2x: links the rid to coordinate systems
        # CORD1x: links g1,g2,g3 to grid points
        for coord in self.coords.values():
            coord.cross_reference(self)

        for coord in self.coords.values():
            coord.setup()

    def _cross_reference_aero(self, check_caero_element_ids: bool=False) -> None:
        """
        Links up all the aero cards
          - CAEROx, PAEROx, SPLINEx, AECOMP, AELIST, AEPARAM, AESTAT, AESURF, AESURFS
        """
        self.zona.cross_reference()
        for caero in self.caeros.values():
            caero.cross_reference(self)

        for paero in self.paeros.values():
            paero.cross_reference(self)

        for trim in self.trims.values():
            trim.cross_reference(self)

        for csschd in self.csschds.values():
            csschd.cross_reference(self)

        for spline in self.splines.values():
            spline.cross_reference(self)

        for aecomp in self.aecomps.values():
            aecomp.cross_reference(self)

        for aelist in self.aelists.values():
            aelist.cross_reference(self)

        for aeparam in self.aeparams.values():
            aeparam.cross_reference(self)

        #for aestat in self.aestats.values(s):
            #aestat.cross_reference(self)

        for aesurf in self.aesurf.values():
            aesurf.cross_reference(self)

        for aesurfs in self.aesurfs.values():
            aesurfs.cross_reference(self)

        for flutter in self.flutters.values():
            flutter.cross_reference(self)

        for monitor_point in self.monitor_points:
            monitor_point.cross_reference(self)

        if self.aero:
            self.aero.cross_reference(self)
        if self.aeros:
            self.aeros.cross_reference(self)

        if check_caero_element_ids:  # only support CAERO1
            ncaeros = len(self.caeros)
            if ncaeros > 1:
                # we don't need to check the ncaeros=1 case
                i = 0
                min_maxs = zeros((ncaeros, 2), dtype='int32')
                for unused_eid, caero in sorted(self.caeros.items()):
                    min_maxs[i, :] = caero.min_max_eid
                    i += 1
                isort = argsort(min_maxs.ravel())
                expected = arange(ncaeros * 2, dtype='int32')
                if not array_equal(isort, expected):
                    msg = 'CAERO element ids are inconsistent\n'
                    msg += 'isort = %s' % str(isort)
                    raise RuntimeError(msg)

            #'AERO',     ## aero
            #'AEROS',    ## aeros
            #'GUST',     ## gusts
            #'FLUTTER',  ## flutters
            #'FLFACT',   ## flfacts
            #'MKAERO1', 'MKAERO2',  ## mkaeros
            #'AECOMP',   ## aecomps
            #'AEFACT',   ## aefacts
            #'AELINK',   ## aelinks
            #'AELIST',   ## aelists
            #'AEPARAM',  ## aeparams
            #'AESTAT',   ## aestats
            #'AESURF',  ## aesurfs

    def _cross_reference_nodes(self) -> None:
        """Links the nodes to coordinate systems"""
        grdset = self.grdset
        for node in self.nodes.values():
            try:
                node.cross_reference(self, grdset)
            except:
                self.log.error("Couldn't cross reference GRID.\n%s" % (str(node)))
                raise

        for point in self.points.values():
            try:
                point.cross_reference(self)
            except:
                self.log.error("Couldn't cross reference POINT.\n%s" % (str(point)))
                raise

        # SPOINTs, EPOINTs don't need xref

        # GRDPNT for mass calculations
        #if model.has_key()
        #for param_key, param in self.params:
            #if

    def _cross_reference_elements(self) -> None:
        """
        Links the elements to nodes, properties (and materials depending on
        the card).
        """
        for elem in self.elements.values():
            try:
                elem.cross_reference(self)
            except (SyntaxError, RuntimeError, AssertionError, KeyError, ValueError) as error:
                self._store_xref_error(error, elem)

        for elem in self.masses.values():
            try:
                elem.cross_reference(self)
            except (SyntaxError, RuntimeError, AssertionError, KeyError, ValueError) as error:
                self._store_xref_error(error, elem)

        for elem in self.rigid_elements.values():
            try:
                elem.cross_reference(self)
            except (SyntaxError, RuntimeError, AssertionError, KeyError, ValueError) as error:
                self._store_xref_error(error, elem)

        for elem in self.plotels.values():
            try:
                elem.cross_reference(self)
            except (SyntaxError, RuntimeError, AssertionError, KeyError, ValueError) as error:
                self._store_xref_error(error, elem)

    def _store_xref_error(self, error, card) -> None:
        self._ixref_errors += 1
        var = traceback.format_exception_only(type(error), error)
        self._stored_xref_errors.append((card, var))
        if self._ixref_errors > self._nxref_errors:
            self.pop_xref_errors()

    def _cross_reference_nodes_with_elements(self) -> None:
        """Links the nodes to all connected elements"""
        nodes = defaultdict(list)  # type: Dict[int, List[Any]]
        for element in self.elements.values():
            #if element.type in ['CONM2']:
            #    pass
            #else:
            if element.nodes is not None:
                for nid in element.node_ids:
                    if nid is None:
                        continue
                    nodes[nid].append(element)
                    #except AttributeError:
                        #print(element)
                        #print('node = %s' % str(node))
                        #raise
        for node in self.nodes.values():
            node.elements_ref = nodes[node.nid]

    def _cross_reference_masses(self) -> None:
        """
        Links the mass to nodes, properties (and materials depending on
        the card).
        """
        for mass in self.masses.values():
            try:
                mass.cross_reference(self)
            except (SyntaxError, RuntimeError, AssertionError, KeyError, ValueError) as error:
                self._store_xref_error(error, mass)

        for prop in self.properties_mass.values():
            try:
                prop.cross_reference(self)
            except (SyntaxError, RuntimeError, AssertionError, KeyError, ValueError) as error:
                self._store_xref_error(error, prop)

    def _cross_reference_properties(self) -> None:
        """Links the properties to materials"""
        for prop in self.properties.values():
            try:
                prop.cross_reference(self)
            except (SyntaxError, RuntimeError, AssertionError, KeyError, ValueError) as error:
                self._store_xref_error(error, prop)

    def _cross_reference_materials(self) -> None:
        """
        Links the materials to materials (e.g. MAT1, CREEP)
        often this is a pass statement
        """
        for mat in self.materials.values():  # MAT1
            try:
                mat.cross_reference(self)
            except (SyntaxError, RuntimeError, AssertionError, KeyError, ValueError) as error:
                self._store_xref_error(error, mat)

        for mat in self.creep_materials.values():  # CREEP
            try:
                mat.cross_reference(self)
            except (SyntaxError, RuntimeError, AssertionError, KeyError, ValueError) as error:
                self._store_xref_error(error, mat)

        # CREEP - depends on MAT1
        data = [self.MATS1, self.MATS3, self.MATS8,
                self.MATT1, self.MATT2, self.MATT3, self.MATT4, self.MATT5,
                self.MATT8, self.MATT9]
        for material_deps in data:
            for mat in material_deps.values():
                try:
                    mat.cross_reference(self)
                except (SyntaxError, RuntimeError, AssertionError, KeyError, ValueError) as error:
                    self._store_xref_error(error, mat)

    def _cross_reference_loads(self) -> None:
        """Links the loads to nodes, coordinate systems, and other loads."""
        for (unused_lid, load_combinations) in self.load_combinations.items():
            for load_combination in load_combinations:
                try:
                    load_combination.cross_reference(self)
                except (SyntaxError, RuntimeError, AssertionError, KeyError, ValueError) as error:
                    self._store_xref_error(error, load_combination)

        for (unused_lid, loads) in self.loads.items():
            for load in loads:
                try:
                    load.cross_reference(self)
                except (SyntaxError, RuntimeError, AssertionError, KeyError, ValueError) as error:
                    self._store_xref_error(error, load)

        for (unused_lid, sid) in self.dloads.items():
            for load in sid:
                #self.log.debug("  dloadi load=%s" % (load))
                try:
                    load.cross_reference(self)
                except (SyntaxError, RuntimeError, AssertionError, KeyError, ValueError) as error:
                    self._ixref_errors += 1
                    var = traceback.format_exception_only(type(error), error)
                    self._stored_xref_errors.append((load, var))
                    if self._ixref_errors > self._nxref_errors:
                        self.pop_xref_errors()

        for unused_lid, sid in self.dload_entries.items():
            for load in sid:
                #self.log.debug("  dloadi load=%s" % (load))
                try:
                    load.cross_reference(self)
                except (SyntaxError, RuntimeError, AssertionError, KeyError, ValueError) as error:
                    #raise
                    self._store_xref_error(error, load)

        for unused_key, darea in self.dareas.items():
            try:
                darea.cross_reference(self)
            except (SyntaxError, RuntimeError, AssertionError, KeyError, ValueError) as error:
                self._store_xref_error(error, darea)

        for unused_key, tic in self.tics.items():
            try:
                tic.cross_reference(self)
            except (SyntaxError, RuntimeError, AssertionError, KeyError, ValueError) as error:
                self._store_xref_error(error, tic)

        for unused_key, dphase in self.dphases.items():
            try:
                dphase.cross_reference(self)
            except (SyntaxError, RuntimeError, AssertionError, KeyError, ValueError) as error:
                self._store_xref_error(error, dphase)

    def _cross_reference_sets(self) -> None:
        """cross references the SET objects"""
        for set_obj in self.asets:
            set_obj.cross_reference(self)
        for set_obj in self.omits:
            set_obj.cross_reference(self)
        for set_obj in self.bsets:
            set_obj.cross_reference(self)
        for set_obj in self.csets:
            set_obj.cross_reference(self)
        for set_obj in self.qsets:
            set_obj.cross_reference(self)
        for unused_name, set_objs in self.usets.items():
            for set_obj in set_objs:
                set_obj.cross_reference(self)

        # superelements
        for unused_key, set_obj in self.se_sets.items():
            set_obj.cross_reference(self)
        for set_obj in self.se_bsets:
            set_obj.cross_reference(self)
        for set_obj in self.se_csets:
            set_obj.cross_reference(self)
        for set_obj in self.se_qsets:
            set_obj.cross_reference(self)
        for set_obj in self.se_usets:
            set_obj.cross_reference(self)

    def _cross_reference_optimization(self) -> None:
        """cross references the optimization objects"""
        for unused_key, deqatn in self.dequations.items():
            deqatn.cross_reference(self)
        for unused_key, dresp in self.dresps.items():
            dresp.cross_reference(self)
        for unused_key, dconstrs in self.dconstrs.items():
            for dconstr in dconstrs:
                dconstr.cross_reference(self)

        for unused_key, dvcrel in self.dvcrels.items():
            dvcrel.cross_reference(self)
        for unused_key, dvmrel in self.dvmrels.items():
            dvmrel.cross_reference(self)
        for unused_key, dvprel in self.dvprels.items():
            dvprel.cross_reference(self)
        for unused_key, desvar in self.desvars.items():
            desvar.cross_reference(self)

    def _safe_cross_reference_contact(self) -> None:
        """cross references the contact objects"""
        self._cross_reference_contact()

    def _cross_reference_contact(self) -> None:
        """cross references the contact objects"""
        for blseg in self.blseg.values():
            blseg.cross_reference(self)
        for bconp in self.bconp.values():
            bconp.cross_reference(self)

    def _uncross_reference_contact(self) -> None:
        """uncross references the contact objects"""
        for blseg in self.blseg.values():
            blseg.uncross_reference()
        for bconp in self.bconp.values():
            bconp.uncross_reference()

    def _cross_reference_superelements(self) -> None:
        """cross references the superelement objects"""
        for unused_seid, csuper in self.csuper.items():
            csuper.cross_reference(self)
        for unused_seid, csupext in self.csupext.items():
            csupext.cross_reference(self)

        for unused_seid, sebulk in self.sebulk.items():
            sebulk.cross_reference(self)
        for unused_seid, sebndry in self.sebndry.items():
            sebndry.cross_reference(self)
        for unused_seid, seconct in self.seconct.items():
            seconct.cross_reference(self)
        for unused_seid, seelt in self.seelt.items():
            seelt.cross_reference(self)
        for unused_seid, seexcld in self.seexcld.items():
            seexcld.cross_reference(self)
        for unused_seid, selabel in self.selabel.items():
            selabel.cross_reference(self)
        for unused_seid, seloc in self.seloc.items():
            seloc.cross_reference(self)
        for unused_seid, seload in self.seload.items():
            seload.cross_reference(self)
        for unused_seid, sempln in self.sempln.items():
            sempln.cross_reference(self)
        for unused_seid, setree in self.setree.items():
            setree.cross_reference(self)

        #'senqset',
        #'se_sets', 'se_usets',

    def _safe_cross_reference_superelements(
            self, create_superelement_geometry: bool=False) -> None:
        xref_errors = {}
        seloc_missing = []
        for seid, seloc in self.seloc.items():
            if seid in self.superelement_models:
                superelement = self.superelement_models[seid]
                seloc.safe_cross_reference(self, xref_errors)
                #seloc.transform(self)
            else:
                seloc_missing.append(seid)

        try:
            for unused_seid, sempln in sorted(self.sempln.items()):
                sempln.safe_cross_reference(self, xref_errors)
            for unused_seid, csuper in self.csuper.items():
                csuper.safe_cross_reference(self, xref_errors)
            for unused_seid, csupext in self.csupext.items():
                csupext.safe_cross_reference(self, xref_errors)

            if self.sebulk and create_superelement_geometry:
                #print('sebulk...')
                import os
                # we have to create the superelement in order to transform it...
                for seid, sebulk in self.sebulk.items():
                    super_filename = 'super_%i.bdf' % seid
                    if os.path.exists(super_filename):
                        os.remove(super_filename)
                    #print(sebulk)
                    rseid = sebulk.rseid
                    sebulk.safe_cross_reference(self, xref_errors)
                    mirror_model = self._create_superelement_from_sebulk(sebulk, seid, rseid)
                    if mirror_model is None:
                        continue
                    self.log.debug('made superelement %i' % seid)
                    self.superelement_models[seid] = mirror_model
                    mirror_model.write_bdf(super_filename)
            for unused_seid, sebndry in self.sebndry.items():
                sebndry.safe_cross_reference(self, xref_errors)
            for unused_seid, seconct in self.seconct.items():
                seconct.safe_cross_reference(self, xref_errors)
            for unused_seid, seelt in self.seelt.items():
                seelt.safe_cross_reference(self, xref_errors)
            for unused_seid, seexcld in self.seexcld.items():
                seexcld.safe_cross_reference(self, xref_errors)
            for unused_seid, selabel in self.selabel.items():
                selabel.safe_cross_reference(self, xref_errors)
            for seid in seloc_missing:
                seloc = self.seloc[seid]
                seloc.safe_cross_reference(self, xref_errors)
            for unused_seid, seload in self.seload.items():
                seload.safe_cross_reference(self, xref_errors)
            for unused_seid, setree in self.setree.items():
                setree.safe_cross_reference(self, xref_errors)
        except KeyError:
            if not create_superelement_geometry:
                raise
            self.write_bdf('superelement_xref.bdf')
            self.log.error('check superelement_xref.bdf')
            raise

    def _create_superelement_from_sebulk(self, sebulk, seid: int, rseid: int) -> None:
        """helper for sebulk"""
        #C:\MSC.Software\MSC.Nastran\msc20051\nast\tpl\see103q4.dat
        ref_model = self.superelement_models[rseid]
        if sebulk.superelement_type == 'MIRROR':
            from pyNastran.bdf.mesh_utils.mirror_mesh import bdf_mirror_plane
            #print('creating superelement %s from %s' % (seid, rseid))
            sempln = self.sempln[seid]
            plane = array([node.get_position() for node in sempln.nodes_ref])

            # What about seloc on the primary and sempln+seloc on the secondary?
            #  - move the primary
            #  - then apply the mirror to make the secondary
            #  - then move the secondary
            #
            # Or what about sempln+seloc on the tertiary?
            #
            # this is fine for the secondary
            if rseid in self.seloc:
                # I think this is wrong...
                seloc = self.seloc[rseid]
                plane = seloc.transform(self, plane)

            ref_model, mirror_model, unused_nid_offset, unused_eid_offset = bdf_mirror_plane(
                ref_model, plane, mirror_model=None, log=None, debug=True, use_nid_offset=False)
            mirror_model.properties = ref_model.properties
            mirror_model.materials = ref_model.materials
            new_model = mirror_model
        elif sebulk.Type in ['MANUAL', 'PRIMARY', 'COLLCTR', 'EXTERNAL']:
            self.log.info('skipping:\n%s' % sebulk)
            new_model = None
        else:  # pragma: no cover
            raise NotImplementedError(sebulk)
        return new_model

    def _uncross_reference_superelements(self) -> None:
        """cross references the superelement objects"""
        for unused_seid, csuper in self.csuper.items():
            csuper.uncross_reference()
        for unused_seid, csupext in self.csupext.items():
            csupext.uncross_reference()

        for unused_seid, sebulk in self.sebulk.items():
            sebulk.uncross_reference()
        for unused_seid, sebndry in self.sebndry.items():
            sebndry.uncross_reference()
        for unused_seid, seconct in self.seconct.items():
            seconct.uncross_reference()
        for unused_seid, seelt in self.seelt.items():
            seelt.uncross_reference()
        for unused_seid, seexcld in self.seexcld.items():
            seexcld.uncross_reference()
        for unused_seid, selabel in self.selabel.items():
            selabel.uncross_reference()
        for unused_seid, seloc in self.seloc.items():
            seloc.uncross_reference()
        for unused_seid, seload in self.seload.items():
            seload.uncross_reference()
        for unused_seid, sempln in self.sempln.items():
            sempln.uncross_reference()
        for unused_seid, setree in self.setree.items():
            setree.uncross_reference()

    def get_point_grids(self, nodes: List[Any], msg: str='') -> None:
        """gets GRID, POINT cards"""
        nodes_ref = []
        missing_nids = []
        for nid in nodes:
            if nid in self.nodes:
                node = self.nodes[nid]
            elif nid in self.points:
                node = self.points[nid]
            else:
                missing_nids.append(nid)
                continue
            nodes_ref.append(node)
        if missing_nids:
            raise KeyError('missing GRID/POINT nids=%s%s' % (missing_nids, msg))
        return nodes_ref

    def superelement_nodes(self, seid: int, nodes: List[Any], msg: str='') -> None:
        if seid == 0:
            return self.Nodes(nodes, msg=msg)
        try:
            superelement = self.superelement_models[seid]
        except KeyError:
            keys = list(self.superelement_models.keys())
            raise KeyError('cant find superelement=%i%s; seids=%s' % (seid, msg, keys))
        return superelement.Nodes(nodes, msg=msg)

    def geom_check(self, geom_check: bool, xref: bool) -> None:  # pragma: no cover
        """
        what about xref?
        """
        if geom_check:
            if xref:
                for unused_eid, element in self.elements.values():
                    #element.Mass()
                    element._verify(xref=True)
                #if 'GEOMCHECK' in self.params:  # should this be an executive control parameter?
                    #for eid, element in model.elements:
                        #element._verify()
            else:
                for unused_eid, element in self.elements.values():
                    element.verify_unique_node_ids()
                    element._verify(xref=False)

            # aspect ratio - ratio between element edges
            # warping - how planar is a face
            # taper - split a quad into 2 triangles and compare the area
            # skew - an angle, measures how skewed an element face is by drawing lines
            #        between midpoints of elements edges, finding the smallest angle
            #        between the intersecting lines and subtracting that from 90 degrees
            # Jacobian - how much does element deviate from the ideal shape by taking the
            #            determinant of the Jacobian matrix
            # quad skew <= 30.
            # quad warp >= 0.05
            # quad taper >= 0.5
            # quad iamin <= 30.
            # quad iamax >= 150.

            # tria skew <= 10.
            # tria iamax <= 160.

            # tetra ar >= 100.
            # tetra elpr <= 0.5
            # tetra detj <= 0.

            # hex ar >= 100.
            # hex elpr <= 0.5
            # hex detj <= 0.
            # hex warp <= 0.707

            # penta ar >= 100.
            # penta elpr <= 0.5
            # penta detj <= 0.
            # penta warp <= 0.707

            # pyram ar >= 100.
            # pyram elpr <= 0.5
            # pyram detj <= 0.
            # pyram warp <= 0.707
