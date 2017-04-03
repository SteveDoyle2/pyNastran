"""
Creates safe cross referencing

Safe cross-referencing skips failed xref's
"""
from __future__ import print_function
import traceback

from six import iteritems, itervalues
from numpy import zeros, argsort, arange, array_equal

from pyNastran.bdf.bdf_interface.cross_reference import XrefMesh

class SafeXrefMesh(XrefMesh):
    """
    Links up the various cards in the BDF.
    """
    def __init__(self):
        """
        The main BDF class defines all the parameters that are used.
        """
        XrefMesh.__init__(self)

    # def geom_check(self):
        # """
        # Performs various geometry checks
          # 1.  nodal uniqueness on elements
        # """
        # for elem in model.elements:
            # elem.check_unique_nodes()

    def safe_cross_reference(self, xref=True,
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
                             xref_optimization=True,
                             debug=True):
        """
        Performs cross referencing in a way that skips data gracefully.

        .. warning:: not fully implemented
        """
        if not xref:
            return
        self.log.debug("Safe Cross Referencing...")
        if xref_nodes:
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
            self._safe_cross_reference_constraints()
        if xref_loads:
            self._safe_cross_reference_loads(debug=debug)
        if xref_optimization:
            self._cross_reference_optimization()
        if xref_nodes_with_elements:
            self._cross_reference_nodes_with_elements()

    def _safe_cross_reference_constraints(self):
        """
        Links the SPCADD, SPC, SPCAX, SPCD, MPCADD, MPC, SUPORT,
        SUPORT1, SESUPORT cards.
        """
        for spcs in itervalues(self.spcs):
            for spc in spcs:
                spc.safe_cross_reference(self)
        for spcoffs in itervalues(self.spcoffs):
            for spcoff in spcoffs:
                spcoff.safe_cross_reference(self)

        for mpcs in itervalues(self.mpcs):
            for mpc in mpcs:
                mpc.safe_cross_reference(self)

        for suport in self.suport:
            #if hasattr(suport, 'safe_cross_reference'):
            suport.safe_cross_reference(self)
            #else:
                #suport.cross_reference(self)
                #self.log.warning('%s - add safe_cross_reference' % suport.type)

        for suport1_id, suport1 in iteritems(self.suport1):
            #if hasattr(suport1, 'safe_cross_reference'):
            suport1.safe_cross_reference(self)
            #else:
                #suport1.cross_reference(self)
                #self.log.warning('%s - add safe_cross_reference' % suport1.type)

        for se_suport in self.se_suport:
            if hasattr(se_suport, 'safe_cross_reference'):
                se_suport.safe_cross_reference(self)
            else:
                se_suport.cross_reference(self)
                self.log.warning('%s - add safe_cross_reference' % se_suport.type)

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

    def _safe_cross_reference_aero(self):
        """
        Links up all the aero cards
          - CAEROx, PAEROx, SPLINEx, AECOMP, AELIST, AEPARAM, AESTAT, AESURF, AESURFS
        """
        for caero in itervalues(self.caeros):
            caero.safe_cross_reference(self)

        for paero in itervalues(self.paeros):
            paero.safe_cross_reference(self)

        for trim in itervalues(self.trims):
            trim.safe_cross_reference(self)
        for csschd in itervalues(self.csschds):
            csschd.safe_cross_reference(self)

        for spline in itervalues(self.splines):
            spline.safe_cross_reference(self)

        for aecomp in itervalues(self.aecomps):
            if hasattr(aecomp, 'safe_cross_reference'):
                aecomp.safe_cross_reference(self)
            else:
                aecomp.cross_reference(self)
                self.log.warning('%s - add safe_cross_reference' % aecomp.type)

        for aelist in itervalues(self.aelists):
            aelist.safe_cross_reference(self)

        for aeparam in itervalues(self.aeparams):
            if hasattr(aeparam, 'safe_cross_reference'):
                aeparam.safe_cross_reference(self)
            else:
                aeparam.cross_reference(self)
                self.log.warning('%s - add safe_cross_reference' % aeparam.type)

        #for aestat in itervalues(self.aestats):
            #aestat.safe_cross_reference(self)

        for aesurf in itervalues(self.aesurf):
            aesurf.safe_cross_reference(self)

        for aesurfs in itervalues(self.aesurfs):
            if hasattr(aesurfs, 'safe_cross_reference'):
                aesurfs.safe_cross_reference(self)
            else:
                aesurfs.cross_reference(self)
                self.log.warning('%s - add safe_cross_reference' % aesurfs.type)

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

    def safe_get_elements(self, eids, msg=''):
        """safe xref version of self.Elements(eid, msg='')"""
        elements = []
        msgi = ''
        for eid in eids:
            try:
                element = self.Element(eid)
            except KeyError:
                element = eid
                msgi += msg % (eid)
            elements.append(element)
        return elements, msgi
