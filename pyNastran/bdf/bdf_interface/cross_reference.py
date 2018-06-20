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
from collections import defaultdict
import traceback
from typing import List, Dict, Any
from six import iteritems, itervalues

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

    def cross_reference(self,
                        xref=True,
                        xref_nodes=True,
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
        # type: (bool, bool, bool, bool, bool, bool, bool, bool, bool, bool, bool, bool) -> None
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
        self.log.debug("Cross Referencing...")
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
        #self.case_control_deck.cross_reference(self)
        self.pop_xref_errors()

    def _cross_reference_constraints(self):
        # type: () -> None
        """
        Links the SPCADD, SPC, SPCAX, SPCD, MPCADD, MPC, SUPORT,
        SUPORT1, SESUPORT cards.
        """
        for spcadds in itervalues(self.spcadds):
            for spcadd in spcadds:
                spcadd.cross_reference(self)
        for spcs in itervalues(self.spcs):
            for spc in spcs:
                spc.cross_reference(self)
        for spcoffs in itervalues(self.spcoffs):
            for spcoff in spcoffs:
                spcoff.cross_reference(self)

        for mpcadds in itervalues(self.mpcadds):
            for mpcadd in mpcadds:
                mpcadd.cross_reference(self)
        for mpcs in itervalues(self.mpcs):
            for mpc in mpcs:
                mpc.cross_reference(self)

        for suport in self.suport:
            suport.cross_reference(self)

        for unused_suport1_id, suport1 in iteritems(self.suport1):
            suport1.cross_reference(self)

        for se_suport in self.se_suport:
            se_suport.cross_reference(self)

    def _cross_reference_coordinates(self):
        # type: () -> None
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

    def _cross_reference_aero(self, check_caero_element_ids=False):
        # type: () -> None
        """
        Links up all the aero cards
          - CAEROx, PAEROx, SPLINEx, AECOMP, AELIST, AEPARAM, AESTAT, AESURF, AESURFS
        """
        self.zona.cross_reference()
        for caero in itervalues(self.caeros):
            caero.cross_reference(self)

        for paero in itervalues(self.paeros):
            paero.cross_reference(self)

        for trim in itervalues(self.trims):
            trim.cross_reference(self)

        for csschd in itervalues(self.csschds):
            csschd.cross_reference(self)

        for spline in itervalues(self.splines):
            spline.cross_reference(self)

        for aecomp in itervalues(self.aecomps):
            aecomp.cross_reference(self)

        for aelist in itervalues(self.aelists):
            aelist.cross_reference(self)

        for aeparam in itervalues(self.aeparams):
            aeparam.cross_reference(self)

        #for aestat in itervalues(self.aestats):
            #aestat.cross_reference(self)

        for aesurf in itervalues(self.aesurf):
            aesurf.cross_reference(self)

        for aesurfs in itervalues(self.aesurfs):
            aesurfs.cross_reference(self)

        for flutter in itervalues(self.flutters):
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
                for unused_eid, caero in sorted(iteritems(self.caeros)):
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
        # type: () -> None
        """
        Links the nodes to coordinate systems
        """
        grdset = self.grdset
        for node in itervalues(self.nodes):
            try:
                node.cross_reference(self, grdset)
            except:
                self.log.error("Couldn't cross reference GRID.\n%s" % (str(node)))
                raise
        # SPOINTs, EPOINTs don't need xref

        # GRDPNT for mass calculations
        #if model.has_key()
        #for param_key, param in self.params:
            #if

    def _cross_reference_elements(self):
        # type: () -> None
        """
        Links the elements to nodes, properties (and materials depending on
        the card).
        """
        for elem in itervalues(self.elements):
            try:
                elem.cross_reference(self)
            except (SyntaxError, RuntimeError, AssertionError, KeyError, ValueError) as error:
                self._ixref_errors += 1
                var = traceback.format_exception_only(type(error), error)
                self._stored_xref_errors.append((elem, var))
                if self._ixref_errors > self._nxref_errors:
                    self.pop_xref_errors()

        for elem in itervalues(self.rigid_elements):
            try:
                elem.cross_reference(self)
            except (SyntaxError, RuntimeError, AssertionError, KeyError, ValueError) as error:
                self._ixref_errors += 1
                var = traceback.format_exception_only(type(error), error)
                self._stored_xref_errors.append((elem, var))
                if self._ixref_errors > self._nxref_errors:
                    self.pop_xref_errors()

        for elem in itervalues(self.plotels):
            try:
                elem.cross_reference(self)
            except (SyntaxError, RuntimeError, AssertionError, KeyError, ValueError) as error:
                self._ixref_errors += 1
                var = traceback.format_exception_only(type(error), error)
                self._stored_xref_errors.append((elem, var))
                if self._ixref_errors > self._nxref_errors:
                    self.pop_xref_errors()

    def _cross_reference_nodes_with_elements(self):
        # type: () -> None
        """
        Links the nodes to all connected elements
        """
        nodes = defaultdict(list)  # type: Dict[int, List[Any]]
        for element in itervalues(self.elements):
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
        for node in itervalues(self.nodes):
            node.elements_ref = nodes[node.nid]

    def _cross_reference_masses(self):
        # type: () -> None
        """
        Links the mass to nodes, properties (and materials depending on
        the card).
        """
        for mass in itervalues(self.masses):
            try:
                mass.cross_reference(self)
            except (SyntaxError, RuntimeError, AssertionError, KeyError, ValueError) as error:
                self._ixref_errors += 1
                var = traceback.format_exception_only(type(error), error)
                self._stored_xref_errors.append((mass, var))
                if self._ixref_errors > self._nxref_errors:
                    self.pop_xref_errors()

        for prop in itervalues(self.properties_mass):
            try:
                prop.cross_reference(self)
            except (SyntaxError, RuntimeError, AssertionError, KeyError, ValueError) as error:
                self._ixref_errors += 1
                var = traceback.format_exception_only(type(error), error)
                self._stored_xref_errors.append((prop, var))
                if self._ixref_errors > self._nxref_errors:
                    self.pop_xref_errors()

    def _cross_reference_properties(self):
        # type: () -> None
        """
        Links the properties to materials
        """
        for prop in itervalues(self.properties):
            try:
                prop.cross_reference(self)
            except (SyntaxError, RuntimeError, AssertionError, KeyError, ValueError) as error:
                self._ixref_errors += 1
                var = traceback.format_exception_only(type(error), error)
                self._stored_xref_errors.append((prop, var))
                if self._ixref_errors > self._nxref_errors:
                    self.pop_xref_errors()

    def _cross_reference_materials(self):
        # type: () -> None
        """
        Links the materials to materials (e.g. MAT1, CREEP)
        often this is a pass statement
        """
        for mat in itervalues(self.materials):  # MAT1
            try:
                mat.cross_reference(self)
            except (SyntaxError, RuntimeError, AssertionError, KeyError, ValueError) as error:
                self._ixref_errors += 1
                var = traceback.format_exception_only(type(error), error)
                self._stored_xref_errors.append((mat, var))
                if self._ixref_errors > self._nxref_errors:
                    self.pop_xref_errors()

        for mat in itervalues(self.creep_materials):  # CREEP
            try:
                mat.cross_reference(self)
            except (SyntaxError, RuntimeError, AssertionError, KeyError, ValueError) as error:
                self._ixref_errors += 1
                var = traceback.format_exception_only(type(error), error)
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
                except (SyntaxError, RuntimeError, AssertionError, KeyError, ValueError) as error:
                    self._ixref_errors += 1
                    var = traceback.format_exception_only(type(error), error)
                    self._stored_xref_errors.append((mat, var))
                    if self._ixref_errors > self._nxref_errors:
                        self.pop_xref_errors()

    def _cross_reference_loads(self):
        # type: () -> None
        """
        Links the loads to nodes, coordinate systems, and other loads.
        """
        for (unused_lid, load_combinations) in iteritems(self.load_combinations):
            for load_combination in load_combinations:
                try:
                    load_combination.cross_reference(self)
                except (SyntaxError, RuntimeError, AssertionError, KeyError, ValueError) as error:
                    self._ixref_errors += 1
                    var = traceback.format_exception_only(type(error), error)
                    self._stored_xref_errors.append((load_combination, var))
                    if self._ixref_errors > self._nxref_errors:
                        self.pop_xref_errors()

        for (unused_lid, loads) in iteritems(self.loads):
            for load in loads:
                try:
                    load.cross_reference(self)
                except (SyntaxError, RuntimeError, AssertionError, KeyError, ValueError) as error:
                    self._ixref_errors += 1
                    var = traceback.format_exception_only(type(error), error)
                    self._stored_xref_errors.append((load, var))
                    if self._ixref_errors > self._nxref_errors:
                        self.pop_xref_errors()

        for (unused_lid, sid) in iteritems(self.dloads):
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

        for (unused_lid, sid) in iteritems(self.dload_entries):
            for load in sid:
                #self.log.debug("  dloadi load=%s" % (load))
                try:
                    load.cross_reference(self)
                except (SyntaxError, RuntimeError, AssertionError, KeyError, ValueError) as error:
                    #raise
                    self._ixref_errors += 1
                    var = traceback.format_exception_only(type(error), error)
                    self._stored_xref_errors.append((load, var))
                    if self._ixref_errors > self._nxref_errors:
                        self.pop_xref_errors()

        for unused_key, darea in iteritems(self.dareas):
            try:
                darea.cross_reference(self)
            except (SyntaxError, RuntimeError, AssertionError, KeyError, ValueError) as error:
                #raise
                self._ixref_errors += 1
                var = traceback.format_exception_only(type(error), error)
                self._stored_xref_errors.append((load, var))
                if self._ixref_errors > self._nxref_errors:
                    self.pop_xref_errors()

        for unused_key, tic in iteritems(self.tics):
            try:
                tic.cross_reference(self)
            except (SyntaxError, RuntimeError, AssertionError, KeyError, ValueError) as error:
                #raise
                self._ixref_errors += 1
                var = traceback.format_exception_only(type(error), error)
                self._stored_xref_errors.append((load, var))
                if self._ixref_errors > self._nxref_errors:
                    self.pop_xref_errors()

        for unused_key, dphase in iteritems(self.dphases):
            try:
                dphase.cross_reference(self)
            except (SyntaxError, RuntimeError, AssertionError, KeyError, ValueError) as error:
                #raise
                self._ixref_errors += 1
                var = traceback.format_exception_only(type(error), error)
                self._stored_xref_errors.append((dphase, var))
                if self._ixref_errors > self._nxref_errors:
                    self.pop_xref_errors()

    def _cross_reference_sets(self):
        # type: () -> None
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
        for unused_name, set_objs in iteritems(self.usets):
            for set_obj in set_objs:
                set_obj.cross_reference(self)

        # superelements
        for unused_key, set_obj in iteritems(self.se_sets):
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
        # type: () -> None
        """cross references the optimization objects"""
        for unused_key, deqatn in iteritems(self.dequations):
            deqatn.cross_reference(self)
        for unused_key, dresp in iteritems(self.dresps):
            dresp.cross_reference(self)
        for unused_key, dconstrs in iteritems(self.dconstrs):
            for dconstr in dconstrs:
                dconstr.cross_reference(self)

        for unused_key, dvcrel in iteritems(self.dvcrels):
            dvcrel.cross_reference(self)
        for unused_key, dvmrel in iteritems(self.dvmrels):
            dvmrel.cross_reference(self)
        for unused_key, dvprel in iteritems(self.dvprels):
            dvprel.cross_reference(self)

    def geom_check(self, geom_check, xref):
        # type: (bool, bool) -> None
        """
        what about xref?
        """
        if geom_check:
            if xref:
                for unused_eid, element in iteritems(self.elements):
                    #element.Mass()
                    element._verify(xref=True)
                #if 'GEOMCHECK' in self.params:  # should this be an executive control parameter?
                    #for eid, element in model.elements:
                        #element._verify()
            else:
                for unused_eid, element in iteritems(self.elements):
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
