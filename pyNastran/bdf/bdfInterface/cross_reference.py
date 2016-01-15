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
from pyNastran.bdf.bdfInterface.attributes import BDFAttributes

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

        .. warning :: not fully implemented
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

        if xref_aero:
            self._cross_reference_aero()
        if xref_constraints:
            self._cross_reference_constraints()
        if xref_loads:
            self._safe_cross_reference_loads(debug=debug)
        if xref_sets:
            self._cross_reference_sets()
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
        #for aesurf in itervalues(self.aesurf):
            #aesurf.uncross_reference()
        for aesurfs in itervalues(self.aesurfs):
            aesurfs.uncross_reference()

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

        DLOAD, RLOAD1, RLOAD2, TLOAD1, TLOAD2
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


    def convert(self, units_to):
        """
        Parameters
        ----------
        xref : bool
           cross references the model (default=True)

        TODO: not done...
        """
        # units_start = 'in'
        # units_end = 'mm'
        # xyz_scale = in_to_mm
        # xyz_scale = 25.4
        xyz_scale, mass_scale, weight_scale = get_scale_factors(self.units, units_to)
        self._convert_nodes(xyz_scale)
        self._convert_coordinates(xyz_scale)
        self._convert_elements(xyz_scale, mass_scale)
        self._convert_properties(xyz_scale, mass_scale)
        #self._convert_masses()
        self._convert_materials(xyz_scale, mass_scale, weight_scale)

        #self._convert_aero()
        #self._convert_constraints()
        self._convert_loads(xyz_scale, weight_scale)
        #self._convert_sets()
        #self._convert_optimization()

    def _convert_nodes(self, xyz_scale):
        for node in itervalues(self.nodes):
            if node.cp_ref.type in ['CORD1R', 'CORD2R']:
                node.xyz *= xyz_scale
            else:
                # only scale R
                node.xyz[0] *= xyz_scale

    def _convert_coordinates(self, xyz_scale):
        for cid, coord in iteritems(self.coords):
            if cid == 0:
                continue
            #print(coord.object_methods())
            #print(dir(coord))
            if coord.rid_ref.type in ['CORD1R', 'CORD2R']:
                coord.origin *= xyz_scale
            else:
                # only scale R
                coord.origin[0] *= xyz_scale

    def _convert_elements(self, xyz_scale, mass_scale):
        area_scale = xyz_scale ** 2
        moi_scale = xyz_scale ** 4
        nsm_scale = mass_scale
        tri_shells = ['CTRIA3', 'CTRIAX', 'CTRIAX6']
        quad_shells = ['CQUAD4', 'CQUAD', 'CQUAD8', 'CQUADX', 'CQUADX8']
        skip_elements = ['CTETRA', 'CPENTA', 'CHEXA', 'CPYRAM', 'CROD']
        spring_elements  = ['CELAS1', 'CELAS2', 'CELAS3', 'CELAS4']
        for elem in itervalues(self.elements):
            elem_type = elem.type
            if elem_type in skip_elements:
                continue
            if elem_type in spring_elements:
                # TODO: scale k from lb/in to N/in
                continue
            elif elem_type in tri_shells:
                # thickness
                elem.T1 *= xyz_scale
                elem.T2 *= xyz_scale
                elem.T3 *= xyz_scale

                # nsm
                elem.nsm *= nsm_scale
            elif elem_type in quad_shells:
                # thickness
                elem.T1 *= xyz_scale
                elem.T2 *= xyz_scale
                elem.T3 *= xyz_scale
                elem.T4 *= xyz_scale
                # nsm
                elem.nsm *= nsm_scale
            elif elem_type == 'CONROD':
                elem.area *= area_scale
                elem.nsm *= nsm_scale
            elif elem_type == 'CBAR':
                # vector
                pass
            elif elem_type == 'CBEAM':
                # vector
                pass
            else:
                raise NotImplementedError(elem)

    def _convert_properties(self, xyz_scale, mass_scale):
        area_scale = xyz_scale ** 2
        moi_scale = xyz_scale ** 4
        nsm_scale = mass_scale

        skip_properties = ['PSOLID']
        for prop in itervalues(self.properties):
            prop_type = prop.type
            if prop_type in skip_properties:
                continue
            elif prop_type == 'PELAS':
                # TODO: stiffness
                pass
            elif prop_type == 'PROD':
                prop.area *= area_scale
                prop.moi *= moi_scale
            elif prop_type == 'PBAR':
                prop.area *= area_scale
                prop.moi *= moi_scale
            elif prop_type == 'PBEAM':
                prop.area *= area_scale
                prop.moi *= moi_scale
            elif prop_type == 'PSHELL':
                prop.t *= xyz_scale
            elif prop_type in ['PCOMP', 'PCOMPG']:
                for layer in prop.layers:
                    layer.t *= xyz_scale
            else:
                raise NotImplementedError(prop_type)

    def _convert_materials(self, xyz_scale, mass_scale, weight_scale):
        pressure_scale = weight_scale / xyz_scale ** 2
        density_scale = mass_scale / xyz_scale ** 3
        for mat in itervalues(self.materials):
            mat_type = mat.type
            if mat_type == 'MAT1':
                mat.e *= pressure_scale
                mat.g *= pressure_scale
                mat.rho *= density_scale
            else:
                raise NotImplementedError(mat)

    def _convert_loads(self, xyz_scale, weight_scale):
        force_scale = weight_scale
        moment_scale = xyz_scale * weight_scale
        pressure_scale = weight_scale / xyz_scale ** 2
        for loads in itervalues(self.loads):
            assert isinstance(loads, list), loads
            for load in loads: # list
                load_type = load.type
                if load_type in ['LOAD']:
                    pass
                elif load_type == 'FORCE':
                    load.mag *= force_scale
                elif load_type == 'MOMENT':
                    load.mag *= moment_scale
                else:
                    raise NotImplementedError(load)

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

        .. warning:: be careful if you call this method
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
                self.spcObject.append(spc)
                spc.cross_reference(self)

        for mpcadd in itervalues(self.mpcadds):
            self.mpcObject.Add(mpcadd)
            mpcadd.cross_reference(self)

        for mpcs in itervalues(self.mpcs):
            for mpc in mpcs:
                self.mpcObject.append(mpc)
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
        #for aesurf in itervalues(self.aesurf):
            #aesurf.cross_reference(self)
        for aesurfs in itervalues(self.aesurfs):
            aesurfs.cross_reference(self)

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
        gridSet = self.gridSet
        for n in itervalues(self.nodes):
            try:
                n.cross_reference(self, gridSet)
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
            if element.nodes is not None:
                for node in element.nodes:
                    if node is None:
                        continue
                    try:
                        nodes[node.nid].add(element)
                    except AttributeError:
                        print(element)
                        print('node = %s' % str(node))
                        raise
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
            #self.log.debug("lid=%s sid=%s" %(lid, sid))
            for load in sid:
                try:
                    load.cross_reference(self)
                except (SyntaxError, RuntimeError, AssertionError, KeyError, ValueError) as e:
                    raise
                    self._ixref_errors += 1
                    var = traceback.format_exception_only(type(e), e)
                    self._stored_xref_errors.append((load, var))
                    if self._ixref_errors > self._nxref_errors:
                        self.pop_xref_errors()

        for (lid, sid) in iteritems(self.dloads):
            #self.log.debug("lid=%s sid=%s" %(lid, sid))
            for load in sid:
                try:
                    load.cross_reference(self)
                except (SyntaxError, RuntimeError, AssertionError, KeyError, ValueError) as e:
                    raise
                    self._ixref_errors += 1
                    var = traceback.format_exception_only(type(e), e)
                    self._stored_xref_errors.append((load, var))
                    if self._ixref_errors > self._nxref_errors:
                        self.pop_xref_errors()

        for (lid, sid) in iteritems(self.dload_entries):
            for load in sid:
                try:
                    load.cross_reference(self)
                except (SyntaxError, RuntimeError, AssertionError, KeyError, ValueError) as e:
                    self._ixref_errors += 1
                    var = traceback.format_exception_only(type(e), e)
                    self._stored_xref_errors.append((load, var))
                    if self._ixref_errors > self._nxref_errors:
                        self.pop_xref_errors()

        for key, darea in iteritems(self.dareas):
            try:
                darea.cross_reference(self)
            except (SyntaxError, RuntimeError, AssertionError, KeyError, ValueError) as e:
                self._ixref_errors += 1
                var = traceback.format_exception_only(type(e), e)
                self._stored_xref_errors.append((load, var))
                if self._ixref_errors > self._nxref_errors:
                    self.pop_xref_errors()

        for key, dphase in iteritems(self.dphases):
            try:
                dphase.cross_reference(self)
            except (SyntaxError, RuntimeError, AssertionError, KeyError, ValueError) as e:
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
        for key, dconstr in iteritems(self.dconstrs):
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

def get_scale_factors(units_from, units_to):
    scales = []
    # in-lb-s
    # m-kg-s
    length_from = units_from[0]
    length_to = units_to[0]

    mass_from = units_from[1]
    mass_to = units_to[1]

    time_from = units_from[2]
    time_to = units_to[2]
    assert time_from == time_to, 'units_from=%s units_to=%s time_from=%r time_to=%r' % (
        units_from, units_to, time_from, time_to)
    time_scale = 1.0
    xyz_scale = convert_length(length_from, length_to)

    mass_scale = convert_mass(mass_from, mass_to)
    weight_scale = mass_scale * xyz_scale / time_scale ** 2
    return xyz_scale, mass_scale, weight_scale


def convert_length(length_from, length_to):
    xyz_scale = 1.0
    if length_from != length_to:
        if length_from == 'in':
            xyz_scale *= 0.0254
        elif length_from == 'ft':
            xyz_scale *= 0.3048
        elif length_from == 'm':
            #xyz_scale *= 1.0
            pass
        elif length_from == 'cm':
            xyz_scale *= 100.0
        elif length_from == 'mm':
            xyz_scale *= 1000.0
        else:
            raise NotImplementedError(length_from)

        if length_to == 'in':
            xyz_scale /= 0.0254
        elif length_to == 'ft':
            xyz_scale /= 0.3048
        elif length_to == 'm':
            #xyz_scale /= 1.0
            pass
        elif length_to == 'cm':
            xyz_scale /= 100.0
        elif length_to == 'mm':
            xyz_scale /= 1000.0
        else:
            raise NotImplementedError(length_to)
    return xyz_scale

def convert_mass(mass_from, mass_to):
    mass_scale = 1.0
    if mass_from != mass_to:
        if mass_from == 'kg':
            pass
        elif mass_from == 'lb':
            mass_scale *= 0.453592
        else:
            raise NotImplementedError(mass_from)

        if mass_to == 'kg':
            pass
        elif mass_to == 'lb':
            mass_scale /= 0.453592
        else:
            raise NotImplementedError(mass_to)

    return mass_scale
