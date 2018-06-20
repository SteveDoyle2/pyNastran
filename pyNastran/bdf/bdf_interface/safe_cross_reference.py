"""
Creates safe cross referencing

Safe cross-referencing skips failed xref's
"""
from __future__ import print_function
from collections import defaultdict
import traceback
from typing import List, Dict, Any
from six import iteritems, itervalues

import numpy as np
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
        if xref_sets:
            self._cross_reference_sets()
        if xref_optimization:
            self._cross_reference_optimization()
        if xref_nodes_with_elements:
            self._cross_reference_nodes_with_elements()
        self.pop_xref_errors()

    def _safe_cross_reference_constraints(self):
        # type: () -> None
        """
        Links the SPCADD, SPC, SPCAX, SPCD, MPCADD, MPC, SUPORT,
        SUPORT1, SESUPORT cards.
        """
        for spcadds in itervalues(self.spcadds):
            for spcadd in spcadds:
                spcadd.safe_cross_reference(self)
        for spcs in itervalues(self.spcs):
            for spc in spcs:
                spc.safe_cross_reference(self)
        for spcoffs in itervalues(self.spcoffs):
            for spcoff in spcoffs:
                spcoff.safe_cross_reference(self)

        for mpcadds in itervalues(self.mpcadds):
            for mpcadd in mpcadds:
                mpcadd.safe_cross_reference(self)
        for mpcs in itervalues(self.mpcs):
            for mpc in mpcs:
                mpc.safe_cross_reference(self)

        for suport in self.suport:
            suport.safe_cross_reference(self)

        for suport1_id, suport1 in iteritems(self.suport1):
            suport1.safe_cross_reference(self)

        for se_suport in self.se_suport:
            se_suport.safe_cross_reference(self)

    def _safe_cross_reference_aero(self):
        # type: () -> None
        """
        Links up all the aero cards
          - CAEROx, PAEROx, SPLINEx, AECOMP, AELIST, AEPARAM, AESTAT, AESURF, AESURFS
        """
        self.zona.safe_cross_reference()
        xref_errors = defaultdict(list)
        for caero in itervalues(self.caeros):
            caero.safe_cross_reference(self, xref_errors)
        self._show_safe_xref_errors('caeros', xref_errors)

        xref_errors = defaultdict(list)
        for paero in itervalues(self.paeros):
            paero.safe_cross_reference(self, xref_errors)
        self._show_safe_xref_errors('paeros', xref_errors)

        for trim in itervalues(self.trims):
            trim.safe_cross_reference(self)

        xref_errors = defaultdict(list)
        for csschd in itervalues(self.csschds):
            csschd.safe_cross_reference(self, xref_errors)
        self._show_safe_xref_errors('csschds', xref_errors)

        xref_errors = defaultdict(list)
        for spline in itervalues(self.splines):
            spline.safe_cross_reference(self, xref_errors)
        self._show_safe_xref_errors('splines', xref_errors)

        for aecomp in itervalues(self.aecomps):
            aecomp.safe_cross_reference(self)

        for aelist in itervalues(self.aelists):
            aelist.safe_cross_reference(self)

        for aeparam in itervalues(self.aeparams):
            aeparam.safe_cross_reference(self)

        #for aestat in itervalues(self.aestats):
            #aestat.safe_cross_reference(self)

        xref_errors = defaultdict(list)
        for aesurf in itervalues(self.aesurf):
            aesurf.safe_cross_reference(self, xref_errors)
        self._show_safe_xref_errors('caeros', xref_errors)

        for aesurfs in itervalues(self.aesurfs):
            aesurfs.safe_cross_reference(self)

        for flutter in itervalues(self.flutters):
            flutter.safe_cross_reference(self)

        xref_errors = defaultdict(list)
        for monitor_point in self.monitor_points:
            monitor_point.safe_cross_reference(self, xref_errors)
        self._show_safe_xref_errors('monitor_points', xref_errors)

        if self.aero:
            xref_errors = defaultdict(list)
            self.aero.safe_cross_reference(self, xref_errors)
            self._show_safe_xref_errors('aero', xref_errors)
        if self.aeros:
            xref_errors = defaultdict(list)
            self.aeros.safe_cross_reference(self, xref_errors)
            self._show_safe_xref_errors('aeros', xref_errors)

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

    def _safe_cross_reference_elements(self):
        # type: () -> None
        """
        Links the elements to nodes, properties (and materials depending on
        the card).
        """
        xref_errors = defaultdict(list)
        missing_safe_xref = set([])
        for elem in itervalues(self.elements):
            if hasattr(elem, 'safe_cross_reference'):
                try:
                    elem.safe_cross_reference(self, xref_errors)
                except TypeError:
                    self.log.warning('element has not added xref_errors\n%s' % str(elem))
                    raise
            else:
                elem.cross_reference(self)
                missing_safe_xref.add(elem.type)

        for elem in itervalues(self.rigid_elements):
            if hasattr(elem, 'safe_cross_reference'):
                try:
                    elem.safe_cross_reference(self, xref_errors)
                except TypeError:
                    self.log.warning('element has not added xref_errors\n%s' % str(elem))
            else:
                missing_safe_xref.add(elem.type)
                elem.cross_reference(self)

        self._show_safe_xref_errors('elements', xref_errors)

        if missing_safe_xref:
            self.log.warning('These cards dont support safe_xref; %s' % str(list(missing_safe_xref)))

    def _show_safe_xref_errors(self, elements_word, xref_errors):
        """helper method to show errors"""
        if xref_errors:
            msg = 'Failed to safe xref %s\n' % elements_word
            for key, eids_pids in sorted(iteritems(xref_errors)):
                eids = [eid_pid[0] for eid_pid in eids_pids]
                eids.sort()
                pids = np.unique([eid_pid[1] for eid_pid in eids_pids]).tolist()
                msg += 'missing %r for %s = %s\n' % (key, elements_word, eids)
                msg += '%s = %s\n' % (key, pids)
            self.log.warning(msg.rstrip())

    def _safe_cross_reference_loads(self, debug=True):
        # type: (bool) -> None
        """
        Links the loads to nodes, coordinate systems, and other loads.
        """
        xref_errors = defaultdict(list)
        for (lid, load_combinations) in iteritems(self.load_combinations):
            for load_combination in load_combinations:
                try:
                    load_combination.safe_cross_reference(self, xref_errors)
                except TypeError:  # pragma: no cover
                    print(load_combination)
                    raise
        self._show_safe_xref_errors('loads', xref_errors)

        for (lid, loads) in iteritems(self.loads):
            for load in loads:
                try:
                    load.safe_cross_reference(self, xref_errors)
                except TypeError:  # pragma: no cover
                    print(load)
                    raise
        self._show_safe_xref_errors('loads', xref_errors)

        for (lid, sid) in iteritems(self.dloads):
            for load in sid:
                load.safe_cross_reference(self, xref_errors)
        for (lid, sid) in iteritems(self.dload_entries):
            for load in sid:
                try:
                    load.safe_cross_reference(self, xref_errors)
                except TypeError:  # pragma: no cover
                    print(load)
                    raise

        for key, darea in iteritems(self.dareas):
            try:
                darea.safe_cross_reference(self, xref_errors)
            except TypeError:  # pragma: no cover
                print(darea)
                raise
        for key, dphase in iteritems(self.dphases):
            try:
                dphase.safe_cross_reference(self, xref_errors)
            except TypeError:  # pragma: no cover
                print(dphase)
                raise

    def safe_empty_nodes(self, nids, msg=''):
        """safe xref version of self.Nodes(nid, msg='')"""
        nodes = []
        missing_nodes = []
        for nid in nids:
            try:
                node = self.EmptyNode(nid)
            except KeyError:
                node = nid
                missing_nodes.append(nid)
            nodes.append(node)
        if missing_nodes:
            missing_nodes.sort()
            self.log.warning('Nodes %s are missing%s' % (str(missing_nodes), msg))
        return nodes, missing_nodes

    def safe_get_nodes(self, nids, msg=''):
        """safe xref version of self.Nodes(nid, msg='')"""
        nodes = []
        error_nodes = []
        msgi = ''
        for nid in nids:
            try:
                node = self.Node(nid)
            except KeyError:
                error_nodes.append(str(nid))
                node = nid
            nodes.append(nid)
        if error_nodes:
            msgi += 'Could not find nodes %s%s\n' % (', '.join(error_nodes), msg)

        return nodes, msgi

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

    def safe_property(self, pid, ref_id, xref_errors, msg=''):
        """
        Parameters
        ----------
        ref_id : int
            the referencing value (e.g., an element references a property)
        """
        try:
            pid_ref = self.Property(pid, msg=msg)
        except KeyError:
            pid_ref = None
            #self.log.error('cant find Property=%s%s' % (mid, msg))
            xref_errors['pid'].append((ref_id, pid))
        return pid_ref

    def safe_material(self, mid, ref_id, xref_errors, msg=''):
        """
        Parameters
        ----------
        ref_id : int
            the referencing value (e.g., an property references a material)
        """
        try:
            mid_ref = self.Material(mid, msg=msg)
        except KeyError:
            mid_ref = None
            #self.log.error('cant find Material=%s%s' % (mid, msg))
            xref_errors['mid'].append((ref_id, mid))
        return mid_ref

    def safe_coord(self, cid, ref_id, xref_errors, msg=''):
        """
        Parameters
        ----------
        ref_id : int
            the referencing value (e.g., an node and element references a coord)
        """
        try:
            cid_ref = self.Coord(cid, msg=msg)
        except KeyError:
            cid_ref = None
            #self.log.error('cant find cid=%s%s' % (cid, msg))
            xref_errors['cid'].append((ref_id, cid))
        return cid_ref

    def safe_aefact(self, aefact_id, ref_id, xref_errors, msg=''):
        """
        Parameters
        ----------
        ref_id : int
            the referencing value (e.g., an CAERO eid references a AEFACT)
        """
        try:
            aefact_ref = self.AEFact(aefact_id, msg=msg)
        except KeyError:
            aefact_ref = None
            #self.log.error('cant find AFEACT=%s%s' % (aefact_id, msg))
            xref_errors['aefact'].append((ref_id, aefact_id))
        return aefact_ref
