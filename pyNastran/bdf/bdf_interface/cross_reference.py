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
  CORD2S, 3, 1, 0., 0., 0., 0., 0., 1.,
          1., 0., 0.
  # get the position in the global frame
  >>> node.get_position()
  [4., 5., 6.]

  # get the position with respect to another frame
  >>> node.PositionWRT(model, cid=2)
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

  # get the position in the global frame
  >>> node.Position()
  Error!

Cross-referencing allows you to easily jump across cards and also helps
with calculating things like position, area, and mass.  The BDF is designed
around the idea of cross-referencing, so it's recommended that you use it.
"""
# pylint: disable=R0902,R0904,R0914

from __future__ import print_function
from six import iteritems, itervalues
from collections import defaultdict
import traceback
from numpy import zeros, argsort, arange, array_equal
from pyNastran.bdf.bdf_interface.attributes import BDFAttributes

class XrefMesh(BDFAttributes):
    """
    Links up the various cards in the BDF.
    """
    def __init__(self):
        """
        The main BDF class defines all the parameters that are used.
        """
        BDFAttributes.__init__(self)
        self._ixref_errors = 0
        self._nxref_errors = 100
        self._stop_on_xref_error = True
        self._stored_xref_errors = []

    # def geom_check(self):
        # """
        # Performs various geometry checks
          # 1.  nodal uniqueness on elements
        # """
        # for elem in model.elements:
            # elem.check_unique_nodes()

    def safe_cross_reference(self, xref=True,
                             xref_elements=True,
                             xref_nodes_with_elements=True,
                             xref_properties=True,
                             xref_masses=True,
                             xref_materials=True,
                             xref_loads=True,
                             xref_constraints=True,
                             xref_aero=True,
                             xref_sets=True,
                             xref_optimization=True,
                             debug=True):
        """
        Performs cross referencing in a way that skips data gracefully.

        .. warning:: not fully implemented
        """
        self._cross_reference_nodes()
        self._cross_reference_coordinates()

        if xref_elements:
            self._safe_cross_reference_elements()
        if xref_properties:
            self._cross_reference_properties()
        if xref_masses:
            self._cross_reference_masses()
        if xref_materials:
            self._cross_reference_materials()

        if xref_sets:
            self._cross_reference_sets()
        if xref_aero:
            self._safe_cross_reference_aero()
        if xref_constraints:
            self._cross_reference_constraints()
        if xref_loads:
            self._safe_cross_reference_loads(debug=debug)
        if xref_optimization:
            self._cross_reference_optimization()
        if xref_nodes_with_elements:
            self._cross_reference_nodes_with_elements()


    def _safe_cross_reference_elements(self):
        """
        Links the elements to nodes, properties (and materials depending on
        the card).
        """
        for elem in itervalues(self.elements):
            try:
                elem.cross_reference(self)
            except (SyntaxError, RuntimeError, AssertionError, KeyError, ValueError) as e:
                self._ixref_errors += 1
                var = traceback.format_exception_only(type(e), e)
                self._stored_xref_errors.append((elem, var))
                if self._ixref_errors > self._nxref_errors:
                    self.pop_xref_errors()
                    #msg = "Couldn't cross reference Element.\n%s" % str(elem)
                    #self.log.error(msg)
                    #raise
        for elem in itervalues(self.rigid_elements):
            try:
                elem.safe_cross_reference(self)
            except AttributeError:
                elem.cross_reference(self)

    def uncross_reference(self):
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
        """cross references the GRID objects"""
        for node in itervalues(self.nodes):
            node.uncross_reference()

    def _uncross_reference_coords(self):
        """cross references the CORDx objects"""
        for cid, coord in iteritems(self.coords):
            if cid == 0:
                continue
            coord.uncross_reference()

    def _uncross_reference_elements(self):
        """cross references the element objects"""
        for element in itervalues(self.elements):
            try:
                element.uncross_reference()
            except TypeError:
                raise NotImplementedError('%s.uncross_reference' % element.type)
        for element in itervalues(self.rigid_elements):
            element.uncross_reference()
        for element in itervalues(self.plotels):
            element.uncross_reference()

    def _uncross_reference_properties(self):
        """cross references the property objects"""
        for prop in itervalues(self.properties):
            try:
                prop.uncross_reference()
            except TypeError:
                raise NotImplementedError('%s.uncross_reference' % prop.type)
            except AttributeError:
                print('%s.uncross_reference error' % prop.type)
                raise

    def _uncross_reference_materials(self):
        """cross references the material objects"""
        for material in itervalues(self.materials):
            material.uncross_reference()

    def _uncross_reference_masses(self):
        """cross references the mass objects"""
        for mass in itervalues(self.masses):
            mass.uncross_reference()
        for prop in itervalues(self.properties_mass):
            prop.uncross_reference()

    def _uncross_reference_aero(self):
        """cross references the aero objects"""
        for caero in itervalues(self.caeros):
            caero.uncross_reference()
        for paero in itervalues(self.paeros):
            paero.uncross_reference()
        for spline in itervalues(self.splines):
            spline.uncross_reference()
        for aecomp in itervalues(self.aecomps):
            aecomp.uncross_reference()
        for aelist in itervalues(self.aelists):
            aelist.uncross_reference()
        for aeparam in itervalues(self.aeparams):
            aeparam.uncross_reference()
        for aestat in itervalues(self.aestats):
            aestat.uncross_reference()
        for aesurf in itervalues(self.aesurf):
            aesurf.uncross_reference()
        for aesurfs in itervalues(self.aesurfs):
            aesurfs.uncross_reference()
        for flutter in itervalues(self.flutters):
            flutter.uncross_reference(self)

    def _uncross_reference_constraints(self):
        """
        Unlinks the SPCADD, SPC, SPCAX, SPCD, MPCADD, MPC, SUPORT,
        SUPORT1, SESUPORT cards.
        """
        for spcadd in itervalues(self.spcadds):
            spcadd.uncross_reference()
        for spc in itervalues(self.spcs):
            for spci in spc:
                spci.uncross_reference()
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
        """
        Unlinks the LOAD
        PLOAD1, PLOAD2, PLOAD4
        FORCE, FORCE1, FORCE2
        MOMENT, MOMENT1, MOMENT2

        DLOAD, ACSRCE, RLOAD1, RLOAD2, TLOAD1, TLOAD2
        DPHASE, DAREA

        TEMP
        """
        for (lid, sid) in iteritems(self.loads):
            for load in sid:
                load.uncross_reference()
        for (lid, sid) in iteritems(self.dloads):
            for load in sid:
                load.uncross_reference()
        for (lid, sid) in iteritems(self.dload_entries):
            for load in sid:
                load.uncross_reference()
        for key, darea in iteritems(self.dareas):
            darea.uncross_reference()
        for key, dphase in iteritems(self.dphases):
            dphase.uncross_reference()

    def _uncross_reference_sets(self):
        for set_obj in self.asets:
            set_obj.uncross_reference()
        for set_obj in self.bsets:
            set_obj.uncross_reference()
        for set_obj in self.csets:
            set_obj.uncross_reference()
        for set_obj in self.qsets:
            set_obj.uncross_reference()
        for name, set_objs in iteritems(self.usets):
            for set_obj in set_objs:
                set_obj.uncross_reference()

        # superelements
        for key, set_obj in iteritems(self.se_sets):
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
        for key, deqatn in iteritems(self.dequations):
            deqatn.uncross_reference()
        for key, dresp in iteritems(self.dresps):
            dresp.uncross_reference()
        for key, dconstr in iteritems(self.dconstrs):
            dconstr.uncross_reference()

        for key, dvcrel in iteritems(self.dvcrels):
            dvcrel.uncross_reference()
        for key, dvmrel in iteritems(self.dvmrels):
            dvmrel.uncross_reference()
        for key, dvprel in iteritems(self.dvprels):
            dvprel.uncross_reference()

    def cross_reference(self, xref=True,
                        xref_elements=True,
                        xref_nodes_with_elements=True,
                        xref_properties=True,
                        xref_masses=True,
                        xref_materials=True,
                        xref_loads=True,
                        xref_constraints=True,
                        xref_aero=True,
                        xref_sets=True,
                        xref_optimization=True):
        """
        Links up all the cards to the cards they reference

        Parameters
        ----------
        xref : bool
           cross references the model (default=True)
        xref_element : bool
           set cross referencing of elements (default=True)
        xref_properties : bool
           set cross referencing of properties (default=True)
        xref_masses : bool
           set cross referencing of CMASS/PMASS (default=True)
        xref_materials : bool
           set cross referencing of materials (default=True)
        xref_loads : bool
            set cross referencing of loads (default=True)
        xref_constraints : bool
            set cross referencing of constraints (default=True)
        xref_aero : bool
            set cross referencing of CAERO/SPLINEs (default=True)
        xref_sets : bool
            set cross referencing of SETx (default=True)

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
        if xref:
            self.log.debug("Cross Referencing...")
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
            #self.case_control_deck.cross_reference(self)

    def _cross_reference_constraints(self):
        """
        Links the SPCADD, SPC, SPCAX, SPCD, MPCADD, MPC, SUPORT,
        SUPORT1, SESUPORT cards.
        """
        for spcadd in itervalues(self.spcadds):
            self.spcObject.Add(spcadd)
            spcadd.cross_reference(self)

        for spcs in itervalues(self.spcs):
            for spc in spcs:
                #self.spcObject.append(spc)
                spc.cross_reference(self)

        for mpcadd in itervalues(self.mpcadds):
            #self.mpcObject.Add(mpcadd)
            mpcadd.cross_reference(self)

        for mpcs in itervalues(self.mpcs):
            for mpc in mpcs:
                #self.mpcObject.append(mpc)
                mpc.cross_reference(self)

        for suport in self.suport:
            suport.cross_reference(self)
        for suport1_id, suport1 in iteritems(self.suport1):
            suport1.cross_reference(self)
        for se_suport in self.se_suport:
            se_suport.cross_reference(self)

    def _cross_reference_coordinates(self):
        """
        Links up all the coordinate cards to other coordinate cards and nodes
         - CORD1R, CORD1C, CORD1S
         - CORD2R, CORD2C, CORD2S
        """
        # CORD2x: links the rid to coordinate systems
        # CORD1x: links g1,g2,g3 to grid points
        for coord in itervalues(self.coords):
            coord.cross_reference(self)

        for coord in itervalues(self.coords):
            coord.setup()

    def _cross_reference_aero(self):
        """
        Links up all the aero cards
          - CAEROx, PAEROx, SPLINEx, AECOMP, AELIST, AEPARAM, AESTAT, AESURF, AESURFS
        """
        for caero in itervalues(self.caeros):
            caero.cross_reference(self)
        for paero in itervalues(self.paeros):
            paero.cross_reference(self)
        for trim in itervalues(self.trims):
            trim.cross_reference(self)

        for spline in itervalues(self.splines):
            spline.cross_reference(self)
        for aecomp in itervalues(self.aecomps):
            aecomp.cross_reference(self)
        for aelist in itervalues(self.aelists):
            aelist.cross_reference(self)
        for aeparam in itervalues(self.aeparams):
            aeparam.cross_reference(self)
        for aestat in itervalues(self.aestats):
            aestat.cross_reference(self)
        for aesurf in itervalues(self.aesurf):
            aesurf.cross_reference(self)
        for aesurfs in itervalues(self.aesurfs):
            aesurfs.cross_reference(self)
        for flutter in itervalues(self.flutters):
            flutter.cross_reference(self)

        if self.aero:
            self.aero.cross_reference(self)
        if self.aeros:
            self.aeros.cross_reference(self)

        if 0:  # only support CAERO1
            ncaeros = len(self.caeros)
            if ncaeros > 1:
                # we don't need to check the ncaeros=1 case
                i = 0
                min_maxs = zeros((ncaeros, 2), dtype='int32')
                for eid, caero in sorted(iteritems(self.caeros)):
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

    def _safe_cross_reference_aero(self):
        """
        Links up all the aero cards
          - CAEROx, PAEROx, SPLINEx, AECOMP, AELIST, AEPARAM, AESTAT, AESURF, AESURFS
        """
        return self._cross_reference_aero()
        for caero in itervalues(self.caeros):
            caero.safe_cross_reference(self)
        for paero in itervalues(self.paeros):
            paero.safe_cross_reference(self)
        for trim in itervalues(self.trims):
            trim.safe_cross_reference(self)

        for spline in itervalues(self.splines):
            spline.safe_cross_reference(self)
        for aecomp in itervalues(self.aecomps):
            aecomp.safe_cross_reference(self)
        for aelist in itervalues(self.aelists):
            aelist.safe_cross_reference(self)
        for aeparam in itervalues(self.aeparams):
            aeparam.safe_cross_reference(self)
        for aestat in itervalues(self.aestats):
            aestat.safe_cross_reference(self)
        for aesurf in itervalues(self.aesurf):
            aesurf.safe_cross_reference(self)
        for aesurfs in itervalues(self.aesurfs):
            aesurfs.safe_cross_reference(self)
        for flutter in itervalues(self.flutters):
            flutter.safe_cross_reference(self)

        if 0:  # only support CAERO1
            ncaeros = len(self.caeros)
            if ncaeros > 1:
                # we don't need to check the ncaeros=1 case
                i = 0
                min_maxs = zeros((ncaeros, 2), dtype='int32')
                for eid, caero in sorted(iteritems(self.caeros)):
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


    def _cross_reference_nodes(self):
        """
        Links the nodes to coordinate systems
        """
        grid_set = self.gridSet
        for n in itervalues(self.nodes):
            try:
                n.cross_reference(self, grid_set)
            except:
                self.log.error("Couldn't cross reference GRID.\n%s" % (str(n)))
                raise

        if self.spoints:
            self.spointi = self.spoints.create_spointi()

        # GRDPNT for mass calculations
        #if model.has_key()
        #for param_key, param in self.params:
            #if

    def _cross_reference_elements(self):
        """
        Links the elements to nodes, properties (and materials depending on
        the card).
        """
        for elem in itervalues(self.elements):
            try:
                elem.cross_reference(self)
            except (SyntaxError, RuntimeError, AssertionError, KeyError, ValueError) as e:
                self._ixref_errors += 1
                var = traceback.format_exception_only(type(e), e)
                self._stored_xref_errors.append((elem, var))
                if self._ixref_errors > self._nxref_errors:
                    self.pop_xref_errors()

        for elem in itervalues(self.rigid_elements):
            try:
                elem.cross_reference(self)
            except (SyntaxError, RuntimeError, AssertionError, KeyError, ValueError) as e:
                self._ixref_errors += 1
                var = traceback.format_exception_only(type(e), e)
                self._stored_xref_errors.append((elem, var))
                if self._ixref_errors > self._nxref_errors:
                    self.pop_xref_errors()

        for elem in itervalues(self.plotels):
            try:
                elem.cross_reference(self)
            except (SyntaxError, RuntimeError, AssertionError, KeyError, ValueError) as e:
                self._ixref_errors += 1
                var = traceback.format_exception_only(type(e), e)
                self._stored_xref_errors.append((elem, var))
                if self._ixref_errors > self._nxref_errors:
                    self.pop_xref_errors()

    def _cross_reference_nodes_with_elements(self):
        """
        Links the nodes to all connected elements
        """
        nodes = defaultdict(set)
        for element in itervalues(self.elements):
            #if element.type in ['CONM2']:
            #    pass
            #else:
                if element.nodes is not None:
                    for nid in element.node_ids:
                        if nid is None:
                            continue
                        nodes[nid].add(element)
                        #except AttributeError:
                            #print(element)
                            #print('node = %s' % str(node))
                            #raise
        for node in itervalues(self.nodes):
            node.elements = nodes[node.nid]

    def _cross_reference_masses(self):
        """
        Links the mass to nodes, properties (and materials depending on
        the card).
        """
        for mass in itervalues(self.masses):
            try:
                mass.cross_reference(self)
            except (SyntaxError, RuntimeError, AssertionError, KeyError, ValueError) as e:
                self._ixref_errors += 1
                var = traceback.format_exception_only(type(e), e)
                self._stored_xref_errors.append((mass, var))
                if self._ixref_errors > self._nxref_errors:
                    self.pop_xref_errors()

        for prop in itervalues(self.properties_mass):
            try:
                prop.cross_reference(self)
            except (SyntaxError, RuntimeError, AssertionError, KeyError, ValueError) as e:
                self._ixref_errors += 1
                var = traceback.format_exception_only(type(e), e)
                self._stored_xref_errors.append((prop, var))
                if self._ixref_errors > self._nxref_errors:
                    self.pop_xref_errors()

    def _cross_reference_properties(self):
        """
        Links the properties to materials
        """
        for prop in itervalues(self.properties):
            try:
                prop.cross_reference(self)
            except (SyntaxError, RuntimeError, AssertionError, KeyError, ValueError) as e:
                self._ixref_errors += 1
                var = traceback.format_exception_only(type(e), e)
                self._stored_xref_errors.append((prop, var))
                if self._ixref_errors > self._nxref_errors:
                    self.pop_xref_errors()

    def _cross_reference_materials(self):
        """
        Links the materials to materials (e.g. MAT1, CREEP)
        often this is a pass statement
        """
        for mat in itervalues(self.materials):  # MAT1
            try:
                mat.cross_reference(self)
            except (SyntaxError, RuntimeError, AssertionError, KeyError, ValueError) as e:
                self._ixref_errors += 1
                var = traceback.format_exception_only(type(e), e)
                self._stored_xref_errors.append((mat, var))
                if self._ixref_errors > self._nxref_errors:
                    self.pop_xref_errors()

        # CREEP - depends on MAT1
        data = [self.MATS1, self.MATS3, self.MATS8,
                self.MATT1, self.MATT2, self.MATT3, self.MATT4, self.MATT5,
                self.MATT8, self.MATT9]
        for material_deps in data:
            for mat in itervalues(material_deps):
                try:
                    mat.cross_reference(self)
                except (SyntaxError, RuntimeError, AssertionError, KeyError, ValueError) as e:
                    self._ixref_errors += 1
                    var = traceback.format_exception_only(type(e), e)
                    self._stored_xref_errors.append((mat, var))
                    if self._ixref_errors > self._nxref_errors:
                        self.pop_xref_errors()

    def _cross_reference_loads(self):
        """
        Links the loads to nodes, coordinate systems, and other loads.
        """
        for (lid, sid) in iteritems(self.loads):
            #self.log.debug("load lid=%s sid=%s" %(lid, sid))
            for load in sid:
                try:
                    load.cross_reference(self)
                except (SyntaxError, RuntimeError, AssertionError, KeyError, ValueError) as e:
                    self._ixref_errors += 1
                    var = traceback.format_exception_only(type(e), e)
                    self._stored_xref_errors.append((load, var))
                    if self._ixref_errors > self._nxref_errors:
                        self.pop_xref_errors()

        for (lid, sid) in iteritems(self.dloads):
            #self.log.debug("dload lid=%s sid=%s" % (lid, sid))
            for load in sid:
                #self.log.debug("  dloadi load=%s" % (load))
                try:
                    load.cross_reference(self)
                except (SyntaxError, RuntimeError, AssertionError, KeyError, ValueError) as e:
                    self._ixref_errors += 1
                    var = traceback.format_exception_only(type(e), e)
                    self._stored_xref_errors.append((load, var))
                    if self._ixref_errors > self._nxref_errors:
                        self.pop_xref_errors()

        for (lid, sid) in iteritems(self.dload_entries):
            #self.log.debug("dload_entries lid=%s sid=%s" % (lid, sid))
            for load in sid:
                #self.log.debug("  dloadi load=%s" % (load))
                try:
                    load.cross_reference(self)
                except (SyntaxError, RuntimeError, AssertionError, KeyError, ValueError) as e:
                    #raise
                    self._ixref_errors += 1
                    var = traceback.format_exception_only(type(e), e)
                    self._stored_xref_errors.append((load, var))
                    if self._ixref_errors > self._nxref_errors:
                        self.pop_xref_errors()

        for key, darea in iteritems(self.dareas):
            try:
                darea.cross_reference(self)
            except (SyntaxError, RuntimeError, AssertionError, KeyError, ValueError) as e:
                #raise
                self._ixref_errors += 1
                var = traceback.format_exception_only(type(e), e)
                self._stored_xref_errors.append((load, var))
                if self._ixref_errors > self._nxref_errors:
                    self.pop_xref_errors()

        for key, dphase in iteritems(self.dphases):
            try:
                dphase.cross_reference(self)
            except (SyntaxError, RuntimeError, AssertionError, KeyError, ValueError) as e:
                #raise
                self._ixref_errors += 1
                var = traceback.format_exception_only(type(e), e)
                self._stored_xref_errors.append((load, var))
                if self._ixref_errors > self._nxref_errors:
                    self.pop_xref_errors()

    def _safe_cross_reference_loads(self, debug=True):
        """
        Links the loads to nodes, coordinate systems, and other loads.
        """
        for (lid, sid) in iteritems(self.loads):
            for load in sid:
                load.safe_cross_reference(self)

        for (lid, sid) in iteritems(self.dloads):
            for load in sid:
                load.safe_cross_reference(self)
        for (lid, sid) in iteritems(self.dload_entries):
            for load in sid:
                load.safe_cross_reference(self)

        for key, darea in iteritems(self.dareas):
            darea.safe_cross_reference(self)
        for key, dphase in iteritems(self.dphases):
            dphase.safe_cross_reference(self)

    def _cross_reference_sets(self):
        """cross references the SET objects"""
        for set_obj in self.asets:
            set_obj.cross_reference(self)
        for set_obj in self.bsets:
            set_obj.cross_reference(self)
        for set_obj in self.csets:
            set_obj.cross_reference(self)
        for set_obj in self.qsets:
            set_obj.cross_reference(self)
        for name, set_objs in iteritems(self.usets):
            for set_obj in set_objs:
                set_obj.cross_reference(self)

        # superelements
        for key, set_obj in iteritems(self.se_sets):
            set_obj.cross_reference(self)
        for set_obj in self.se_bsets:
            set_obj.cross_reference(self)
        for set_obj in self.se_csets:
            set_obj.cross_reference(self)
        for set_obj in self.se_qsets:
            set_obj.cross_reference(self)
        for set_obj in self.se_usets:
            set_obj.cross_reference(self)

    def _cross_reference_optimization(self):
        """cross references the optimization objects"""
        for key, deqatn in iteritems(self.dequations):
            deqatn.cross_reference(self)
        for key, dresp in iteritems(self.dresps):
            dresp.cross_reference(self)
        for key, dconstrs in iteritems(self.dconstrs):
            for dconstr in dconstrs:
                dconstr.cross_reference(self)

        for key, dvcrel in iteritems(self.dvcrels):
            dvcrel.cross_reference(self)
        for key, dvmrel in iteritems(self.dvmrels):
            dvmrel.cross_reference(self)
        for key, dvprel in iteritems(self.dvprels):
            dvprel.cross_reference(self)

    def geom_check(self, geom_check, xref):
        """
        what about xref?
        """
        if geom_check:
            if xref:
                for eid, element in iteritems(self.elements):
                    #element.Mass()
                    element._verify(xref=True)
                #if 'GEOMCHECK' in self.params:  # should this be an executive control parameter?
                    #for eid, element in model.elements:
                        #element._verify()
            else:
                for eid, element in iteritems(self.elements):
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

