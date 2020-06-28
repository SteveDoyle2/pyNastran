"""
Creates safe cross referencing

Safe cross-referencing skips failed xref's

"""
from collections import defaultdict
from typing import Tuple, List, Any

import numpy as np
from numpy import zeros, argsort, arange, array_equal
from pyNastran.bdf.bdf_interface.cross_reference import XrefMesh


class SafeXrefMesh(XrefMesh):
    """
    Links up the various cards in the BDF.
    """
    def __init__(self) -> None:
        """The main BDF class defines all the parameters that are used."""
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
                             xref_nodes_with_elements=False,
                             xref_properties=True,
                             xref_masses=True,
                             xref_materials=True,
                             xref_loads=True,
                             xref_constraints=True,
                             xref_aero=True,
                             xref_sets=True,
                             xref_optimization=True,
                             create_superelement_geometry=False,
                             debug=True,
                             word=''):
        """
        Performs cross referencing in a way that skips data gracefully.

        .. warning:: not fully implemented
        """
        if not xref:
            return
        self.log.debug("Safe Cross Referencing%s..." % word)
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
            self._safe_cross_reference_loads()
        if xref_sets:
            self._cross_reference_sets()
        if xref_optimization:
            self._safe_cross_reference_optimization()
        if xref_nodes_with_elements:
            self._cross_reference_nodes_with_elements()

        self._safe_cross_reference_contact()
        self._safe_cross_reference_superelements(create_superelement_geometry)

        self.pop_xref_errors()
        for super_id, superelement in sorted(self.superelement_models.items()):
            superelement.safe_cross_reference(
                xref=xref, xref_nodes=xref_nodes, xref_elements=xref_elements,
                xref_nodes_with_elements=xref_nodes_with_elements,
                xref_properties=xref_properties, xref_masses=xref_masses,
                xref_materials=xref_materials, xref_loads=xref_loads,
                xref_constraints=xref_constraints, xref_aero=xref_aero,
                xref_sets=xref_sets, xref_optimization=xref_optimization,
                word=' (Superelement %i)' % super_id)

    def _safe_cross_reference_constraints(self) -> None:
        """
        Links the SPCADD, SPC, SPCAX, SPCD, MPCADD, MPC, SUPORT,
        SUPORT1, SESUPORT cards.
        """
        for spcadds in self.spcadds.values():
            for spcadd in spcadds:
                spcadd.safe_cross_reference(self)
        for spcs in self.spcs.values():
            for spc in spcs:
                spc.safe_cross_reference(self)
        for spcoffs in self.spcoffs.values():
            for spcoff in spcoffs:
                spcoff.safe_cross_reference(self)

        for mpcadds in self.mpcadds.values():
            for mpcadd in mpcadds:
                mpcadd.safe_cross_reference(self)
        for mpcs in self.mpcs.values():
            for mpc in mpcs:
                mpc.safe_cross_reference(self)

        for suport in self.suport:
            suport.safe_cross_reference(self)

        for unused_suport1_id, suport1 in self.suport1.items():
            suport1.safe_cross_reference(self)

        for se_suport in self.se_suport:
            se_suport.safe_cross_reference(self)

    def _safe_cross_reference_aero(self) -> None:
        """
        Links up all the aero cards
          - CAEROx, PAEROx, SPLINEx, AECOMP, AELIST, AEPARAM, AESTAT, AESURF, AESURFS
        """
        self.zona.safe_cross_reference()
        xref_errors = defaultdict(list)
        for caero in self.caeros.values():
            caero.safe_cross_reference(self, xref_errors)
        self._show_safe_xref_errors('caeros', xref_errors)

        xref_errors = defaultdict(list)
        for paero in self.paeros.values():
            paero.safe_cross_reference(self, xref_errors)
        self._show_safe_xref_errors('paeros', xref_errors)

        for trim in self.trims.values():
            trim.safe_cross_reference(self)
        self._show_safe_xref_errors('trims', xref_errors)

        xref_errors = defaultdict(list)
        for csschd in self.csschds.values():
            csschd.safe_cross_reference(self, xref_errors)
        self._show_safe_xref_errors('csschds', xref_errors)

        xref_errors = defaultdict(list)
        for spline in self.splines.values():
            spline.safe_cross_reference(self, xref_errors)
        self._show_safe_xref_errors('splines', xref_errors)

        for aecomp in self.aecomps.values():
            aecomp.safe_cross_reference(self)

        for aelist in self.aelists.values():
            aelist.safe_cross_reference(self)

        for aeparam in self.aeparams.values():
            aeparam.safe_cross_reference(self)

        #for aestat in self.aestats):
            #aestat.safe_cross_reference(self)

        xref_errors = defaultdict(list)
        for aesurf in self.aesurf.values():
            aesurf.safe_cross_reference(self, xref_errors)
        self._show_safe_xref_errors('aesurf', xref_errors)

        for aesurfs in self.aesurfs.values():
            aesurfs.safe_cross_reference(self)

        for flutter in self.flutters.values():
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

    def _safe_cross_reference_elements(self) -> None:
        """
        Links the elements to nodes, properties (and materials depending on
        the card).
        """
        xref_errors = defaultdict(list)
        missing_safe_xref = set()
        for elem in self.elements.values():
            if hasattr(elem, 'safe_cross_reference'):
                elem.safe_cross_reference(self, xref_errors)
            else:
                elem.cross_reference(self)
                missing_safe_xref.add(elem.type)

        for elem in self.masses.values():
            if hasattr(elem, 'safe_cross_reference'):
                elem.safe_cross_reference(self, xref_errors)
            else:
                elem.cross_reference(self)
                missing_safe_xref.add(elem.type)

        for elem in self.rigid_elements.values():
            #if hasattr(elem, 'safe_cross_reference'):
            elem.safe_cross_reference(self, xref_errors)
            #else:
                #missing_safe_xref.add(elem.type)
                #elem.cross_reference(self)

        self._show_safe_xref_errors('elements', xref_errors)

        if missing_safe_xref:
            self.log.warning('These cards dont support safe_xref; %s' %
                             str(list(missing_safe_xref)))

    def _show_safe_xref_errors(self, elements_word, xref_errors):
        # type: (str, bool) -> None
        """helper method to show errors"""
        if xref_errors:
            msg = 'Failed to safe xref %s\n' % elements_word
            for key, eids_pids in sorted(xref_errors.items()):
                eids = [eid_pid[0] for eid_pid in eids_pids]
                eids.sort()
                pids = [eid_pid[1] for eid_pid in eids_pids]
                try:
                    upids = np.unique(pids).tolist()
                except TypeError:
                    print(msg)
                    print('key = %s' % key)
                    print(' - keys   = %s' % eids)
                    print(' - values = %s' % pids)
                    print("Make sure you don't have Nones in the values")
                    raise
                msg += 'missing %r for %s = %s\n' % (key, elements_word, eids)
                msg += '%s = %s\n' % (key, upids)
            self.log.warning(msg.rstrip())

    def _safe_cross_reference_loads(self):
        # type: (bool) -> None
        """
        Links the loads to nodes, coordinate systems, and other loads.
        """
        xref_errors = defaultdict(list)
        for unused_lid, load_combinations in self.load_combinations.items():
            for load_combination in load_combinations:
                try:
                    load_combination.safe_cross_reference(self, xref_errors)
                except TypeError:  # pragma: no cover
                    print(load_combination)
                    raise
        self._show_safe_xref_errors('loads', xref_errors)

        for unused_lid, loads in self.loads.items():
            for load in loads:
                load.safe_cross_reference(self, xref_errors)
        self._show_safe_xref_errors('loads', xref_errors)

        for unused_lid, sid in self.dloads.items():
            for load in sid:
                load.safe_cross_reference(self, xref_errors)

        for unused_lid, sid in self.dload_entries.items():
            for load in sid:
                load.safe_cross_reference(self, xref_errors)

        for unused_key, darea in self.dareas.items():
            darea.safe_cross_reference(self, xref_errors)

        for unused_key, dphase in self.dphases.items():
            dphase.safe_cross_reference(self, xref_errors)

        for unused_key, tic in self.tics.items():
            tic.safe_cross_reference(self, xref_errors)

    def _safe_cross_reference_optimization(self) -> None:
        """cross references the optimization objects"""
        #self._cross_reference_optimization()
        #return
        xref_errors = defaultdict(list)
        for unused_key, deqatn in self.dequations.items():
            deqatn.safe_cross_reference(self)

        for unused_key, dresp in self.dresps.items():
            dresp.safe_cross_reference(self, xref_errors)

        for unused_key, dconstrs in self.dconstrs.items():
            for dconstr in dconstrs:
                if hasattr(dconstr, 'safe_cross_reference'):
                    dconstr.safe_cross_reference(self)
                else:  # pragma: no cover
                    dconstr.cross_reference(self)

        for unused_key, dvcrel in self.dvcrels.items():
            if hasattr(dvcrel, 'safe_cross_reference'):
                dvcrel.safe_cross_reference(self)
            else:  # pragma: no cover
                dvcrel.cross_reference(self)

        for unused_key, dvmrel in self.dvmrels.items():
            if hasattr(dvmrel, 'safe_cross_reference'):
                dvmrel.safe_cross_reference(self)
            else:  # pragma: no cover
                dvmrel.cross_reference(self)

        for unused_key, dvprel in self.dvprels.items():
            if hasattr(dvprel, 'safe_cross_reference'):
                dvprel.safe_cross_reference(self)
            else:  # pragma: no cover
                dvprel.cross_reference(self)

        for unused_key, desvar in self.desvars.items():
            desvar.safe_cross_reference(self)

        for unused_key, topvar in self.topvar.items():
            topvar.safe_cross_reference(self)

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

    def safe_get_nodes(self, nids: List[int], msg: str='') -> Tuple[List[Any], str]:
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

    def safe_get_points(self, point_ids, msg=''):
        """safe xref version of self.Points(point_ids, msg='')"""
        points = []
        error_points = []
        msgi = ''
        for point_id in point_ids:
            try:
                point = self.Point(point_id)
            except KeyError:
                error_points.append(str(point_id))
                point = point_id
            points.append(point)
        if error_points:
            msgi += 'Could not find POINTs %s%s\n' % (', '.join(error_points), msg)
        return points, msgi

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

    def safe_element(self, eid, ref_id, xref_errors, msg=''):
        """
        Gets an element card

        Parameters
        ----------
        ref_id : int
            the referencing value (e.g., a load references an element)

        ref_id = 10 # PLOAD4
        pid = 42  # CQUAD4
        xref_errors = {'eid' : []}
        self.safe_element(eid, ref_id, xref_errors)

        """
        try:
            eid_ref = self.Element(eid, msg=msg)
        except KeyError:
            eid_ref = None
            #self.log.error('cant find Element=%s%s' % (mid, msg))
            xref_errors['eid'].append((ref_id, eid))
        return eid_ref

    def safe_elements(self, eids, ref_id, xref_errors, msg=''):
        """
        Gets an series of elements

        Doesn't get rigid (RROD, RBAR, RBE2, RBE3, RBAR, RBAR1, RSPLINE, RSSCON)
        or mass (CMASS1, CONM2)

        """
        elements = []
        bad_eids = []
        for eid in eids:
            try:
                # elements.append(self.safe_element(eid, ref_id, xref_errors, msg))
                elements.append(self.Element(eid, msg))
            except KeyError:
                bad_eids.append(eid)
                elements.append(None)
                xref_errors['eid'].append((ref_id, eid))
        #if bad_eids:
            #msg = 'eids=%s not found%s.  Allowed elements=%s' % (
                #bad_eids, msg, np.unique(list(self.elements.keys())))
            #self.log.error(msg)
            #raise KeyError(msg)
        return elements

    def safe_property(self, pid, ref_id, xref_errors, msg=''):
        """
        Parameters
        ----------
        ref_id : int
            the referencing value (e.g., an element references a property)

        ref_id = 10 # CQUAD4
        pid = 42  # PSHELL
        xref_errors = {'pid' : []}
        self.safe_property(pid, ref_id, xref_errors)
        """
        try:
            pid_ref = self.Property(pid, msg=msg)
        except KeyError:
            pid_ref = None
            #self.log.error('cant find Property=%s%s' % (mid, msg))
            xref_errors['pid'].append((ref_id, pid))
        return pid_ref

    def safe_property_mass(self, pid, ref_id, xref_errors, msg=''):
        """
        Gets a mass_property card

        Parameters
        ----------
        ref_id : int
            the referencing value (e.g., an element references a property)
        """
        try:
            pid_ref = self.PropertyMass(pid, msg=msg)
        except KeyError:
            pid_ref = None
            #self.log.error('cant find Property=%s%s' % (mid, msg))
            xref_errors['pid'].append((ref_id, pid))
        return pid_ref

    def safe_material(self, mid, ref_id, xref_errors, msg=''):
        """
        Gets a material card

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
        Gets a CORDx card

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

    def safe_paero(self, paero_id, ref_id, xref_errors, msg=''):
        """
        Gets a PAEROx card

        Parameters
        ----------
        ref_id : int
            the referencing value (e.g., a load references an element)

        ref_id = 10 # CAERO1
        pid = 42  # PAERO1
        xref_errors = {'paero' : []}
        self.safe_element(pid, ref_id, xref_errors)

        """
        try:
            paero_ref = self.PAero(paero_id, msg=msg)
        except KeyError:
            paero_ref = None
            xref_errors['paero'].append((ref_id, paero_id))
        return paero_ref

    def safe_aefact(self, aefact_id, ref_id, xref_errors, msg=''):
        """
        Gets an AEFACT card

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

    def safe_aelist(self, aelist_id, ref_id, xref_errors, msg=''):
        """
        Gets an AELIST card

        Parameters
        ----------
        ref_id : int
            the referencing value (e.g., an AESURF eid references a AELIST)

        """
        try:
            aefact_ref = self.AELIST(aelist_id, msg=msg)
        except KeyError:
            aefact_ref = None
            xref_errors['aelist'].append((ref_id, aelist_id))
        return aefact_ref

    def safe_caero(self, caero_id, ref_id, xref_errors, msg=''):
        try:
            caero_ref = self.CAero(caero_id, msg=msg)
        except KeyError:
            caero_ref = None
            xref_errors['caero'].append((ref_id, caero_id))
        return caero_ref

    def safe_tabled(self, tabled_id, ref_id, xref_errors, msg=''):
        """
        Parameters
        ----------
        ref_id : int
            the referencing value (e.g., an TLOAD1 eid references a TABLED1)
        """
        try:
            tabled_ref = self.TableD(tabled_id, msg=msg)
        except KeyError:
            tabled_ref = None
            xref_errors['tabled'].append((ref_id, tabled_id))
        return tabled_ref

    def safe_tableh(self, tableh_id, ref_id, xref_errors, msg=''):
        """
        Parameters
        ----------
        ref_id : int
            the referencing value (e.g., an MATT1 eid references a TABLEH1)
        """
        try:
            tableh_ref = self.TableH(tableh_id, msg=msg)
        except KeyError:
            tableh_ref = None
            xref_errors['tableh'].append((ref_id, tableh_id))
        return tableh_ref
