from __future__ import annotations
import os
import traceback
from io import StringIO
from collections import defaultdict
from typing import cast, TYPE_CHECKING
import numpy as np

from pyNastran.bdf.errors import CrossReferenceError
from pyNastran.bdf.cards.deqatn import DEQATN
from pyNastran.bdf.cards.optimization import DESVAR, DCONSTR
if TYPE_CHECKING:
    from pyNastran.bdf.bdf import BDF, GRID
    from pyNastran.bdf.cards.base_card import BaseCard


class CrossReference:
    def __init__(self, model: BDF) -> None:
        self.model: BDF = model
        model.log
        self.reset_errors()

    def reset_errors(self) -> None:
        """removes the errors from the model"""
        self._ixref_errors = 0
        self._stored_xref_errors = []
        self._nxref_errors = 100
        self._stop_on_xref_error = True

    def set_error_storage(self, nparse_errors: int=100,
                          stop_on_parsing_error: bool=True,
                          nxref_errors: int=100,
                          stop_on_xref_error: bool=True) -> None:
        """
        Catch parsing errors and store them up to print them out all at once
        (not all errors are caught).

        Parameters
        ----------
        nparse_errors : int
            how many parse errors should be stored
            (default=0; all=None; no storage=0)
        stop_on_parsing_error : bool
            should an error be raised if there
            are parsing errors (default=True)
        nxref_errors : int
            how many cross-reference errors
            should be stored (default=0; all=None; no storage=0)
        stop_on_xref_error : bool
            should an error be raised if there
            are cross-reference errors (default=True)

        """
        #assert isinstance(nparse_errors, int), type(nparse_errors)
        assert isinstance(nxref_errors, int), type(nxref_errors)
        #self._nparse_errors = nparse_errors
        self._nxref_errors = nxref_errors
        #self._stop_on_parsing_error = stop_on_parsing_error
        self._stop_on_xref_error = stop_on_xref_error

    def pop_xref_errors(self) -> None:
        """raises an error if there are cross-reference errors"""
        is_error = False
        #self.model
        #self.model.log
        if self._stop_on_xref_error:
            model = self.model
            if self._ixref_errors == 1 and self._nxref_errors == 0:
                raise
            if self._stored_xref_errors:
                filename_note = ''
                if model.bdf_filename and not isinstance(model.bdf_filename, StringIO):
                    filename_note = f' in {os.path.abspath(model.bdf_filename)!r}'
                msg = f'There are cross-reference errors{filename_note}.\n\n'
                for (card, lines) in self._stored_xref_errors:
                    fname = _add_file(self.model, card)
                    an_error = ''.join(lines)
                    msg += '%s%scard:\n%s\n' % (fname, an_error, card)
                    is_error = True

                if is_error and self._stop_on_xref_error:
                    raise CrossReferenceError(msg.rstrip())

    def _store_xref_error(self, error, card: BaseCard) -> None:
        """stores cross-reference errors and raises the errors when the threshold is reached"""
        self._ixref_errors += 1
        #line = traceback.format_exc()
        lines = traceback.format_exception_only(type(error), error)
        #lines2 = [line.replace(r'\\n', '\n').rstrip() for line in lines]
        #self.model.log.warning(lines)
        #self.model.log.warning(lines2)
        self._stored_xref_errors.append((card, lines))
        if self._ixref_errors > self._nxref_errors:
            self.pop_xref_errors()

    def cross_reference_constraints(self) -> None:
        """
        Links the SPCADD, SPC, SPCAX, SPCD, MPCADD, MPC, SUPORT,
        SUPORT1, SESUPORT cards.
        """
        model = self.model
        for spcadds in model.spcadds.values():
            for spcadd in spcadds:
                spcadd.cross_reference(model)
        for spcs in model.spcs.values():
            for spc in spcs:
                spc.cross_reference(model)
        for spcoffs in model.spcoffs.values():
            for spcoff in spcoffs:
                spcoff.cross_reference(model)

        for mpcadds in model.mpcadds.values():
            for mpcadd in mpcadds:
                mpcadd.cross_reference(model)
        for mpcs in model.mpcs.values():
            for mpc in mpcs:
                mpc.cross_reference(model)

        for suport in model.suport:
            suport.cross_reference(model)

        for unused_suport1_id, suport1 in model.suport1.items():
            suport1.cross_reference(model)

        for se_suport in model.se_suport:
            se_suport.cross_reference(model)

    def safe_cross_reference_constraints(self) -> None:
        """
        Links the SPCADD, SPC, SPCAX, SPCD, MPCADD, MPC, SUPORT,
        SUPORT1, SESUPORT cards.
        """
        model = self.model
        for spcadds in model.spcadds.values():
            for spcadd in spcadds:
                spcadd.safe_cross_reference(model)
        for spcs in model.spcs.values():
            for spc in spcs:
                spc.safe_cross_reference(model)
        for spcoffs in model.spcoffs.values():
            for spcoff in spcoffs:
                spcoff.safe_cross_reference(model)

        for mpcadds in model.mpcadds.values():
            for mpcadd in mpcadds:
                mpcadd.safe_cross_reference(model)
        for mpcs in model.mpcs.values():
            for mpc in mpcs:
                mpc.safe_cross_reference(model)

        for suport in model.suport:
            suport.safe_cross_reference(model)

        for unused_suport1_id, suport1 in model.suport1.items():
            suport1.safe_cross_reference(model)

        for se_suport in model.se_suport:
            se_suport.safe_cross_reference(model)

    def cross_reference_sets(self) -> None:
        """cross references the SET objects"""
        model = self.model
        set_lists = [
            model.asets, model.bsets, model.csets, model.omits,
            model.qsets,
            # superelements
            model.se_bsets, model.se_csets,
            model.se_qsets,
        ]
        set_dicts = [model.se_sets, model.se_usets,]
        names = 'abcoq1234'
        for name, set_list in zip(names, set_lists):
            assert isinstance(set_list, list), (name, set_list)
            for set_obj in set_list:
                set_obj.cross_reference(model)

        for set_dict in set_dicts:
            assert isinstance(set_dict, dict), set_dict
            for set_obj in set_dict.values():
                set_obj.cross_reference(model)

        for unused_name, set_objs in model.usets.items():
            for set_obj in set_objs:
                set_obj.cross_reference(model)

    def safe_cross_reference_sets(self) -> None:
        """cross references the SET objects"""
        model = self.model
        set_lists = [
            model.asets, model.bsets, model.csets, model.omits,
            model.qsets,
            # superelements
            model.se_bsets, model.se_csets,
            model.se_qsets,
        ]
        set_dicts = [model.se_sets, model.se_usets,]
        for set_list in set_lists:
            assert isinstance(set_list, list), set_list
            for set_obj in set_list:
                set_obj.cross_reference(model)

        for set_dict in set_dicts:
            assert isinstance(set_dict, dict), set_dict
            for set_obj in set_dict.values():
                set_obj.cross_reference(model)

        for unused_name, set_objs in model.usets.items():
            for set_obj in set_objs:
                set_obj.cross_reference(model)

    def cross_reference_aero(self, check_caero_element_ids: bool=False) -> None:
        """
        Links up all the aero cards
          - CAEROx, PAEROx, SPLINEx, AECOMP, AELIST, AEPARM, AESTAT, AESURF, AESURFS
        """
        model = self.model
        model.log
        model.zona.cross_reference()
        for caero in model.caeros.values():
            caero.cross_reference(model)

        for paero in model.paeros.values():
            paero.cross_reference(model)

        for trim in model.trims.values():
            trim.cross_reference(model)

        for csschd in model.csschds.values():
            csschd.cross_reference(model)

        for spline in model.splines.values():
            spline.cross_reference(model)

        for aecomp in model.aecomps.values():
            aecomp.cross_reference(model)

        for aelist in model.aelists.values():
            aelist.cross_reference(model)

        for aeparam in model.aeparams.values():
            aeparam.cross_reference(model)

        #for aestat in self.aestats.values():
            #aestat.cross_reference(self)

        for aelinks in model.aelinks.values():
            for aelink in aelinks:
                aelink.cross_reference(model)

        for aesurf in model.aesurf.values():
            aesurf.cross_reference(model)

        for aesurfs in model.aesurfs.values():
            aesurfs.cross_reference(model)

        for flutter in model.flutters.values():
            flutter.cross_reference(model)

        for monitor_point in model.monitor_points:
            monitor_point.cross_reference(model)

        if model.aero:
            model.aero.cross_reference(model)
        if model.aeros:
            model.aeros.cross_reference(model)

        if len(model.caeros) and 'WKK' in model.dmi:
            check_caero_element_ids = True

        if check_caero_element_ids:
            naeroboxes = _check_caero_box_overlap(model)

            if 'W2GJ' in model.dmi:
                w2gj = model.dmi['W2GJ']
                assert w2gj.shape == (naeroboxes, 1), f'naeroboxes={naeroboxes}; w2gj.shape={str(w2gj.shape)}'
            if 'FA2J' in model.dmi:
                fa2j = model.dmi['FA2J']
                assert fa2j.shape == (naeroboxes, 1), f'naeroboxes={naeroboxes}; fa2j.shape={str(fa2j.shape)}'
            if 'WKK' in model.dmi:
                wkk = model.dmi['WKK']
                assert wkk.shape in [(naeroboxes*2, 1), (naeroboxes*2, naeroboxes*2)], f'naeroboxes*2={naeroboxes*2}; wkk.shape={str(wkk.shape)}'

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
            #'AEPARM',  ## aeparams
            #'AESTAT',   ## aestats
            #'AESURF',  ## aesurfs

    def safe_cross_reference_aero(self) -> None:
        """
        Links up all the aero cards
          - CAEROx, PAEROx, SPLINEx, AECOMP, AELIST, AEPARM, AESTAT, AESURF, AESURFS
        """
        model = self.model
        model.zona.safe_cross_reference()
        xref_errors = defaultdict(list)
        for caero in model.caeros.values():
            caero.safe_cross_reference(model, xref_errors)
        self._show_safe_xref_errors('caeros', xref_errors)

        xref_errors = defaultdict(list)
        for paero in model.paeros.values():
            paero.safe_cross_reference(model, xref_errors)
        self._show_safe_xref_errors('paeros', xref_errors)

        for trim in model.trims.values():
            trim.safe_cross_reference(model)
        self._show_safe_xref_errors('trims', xref_errors)

        xref_errors = defaultdict(list)
        for csschd in model.csschds.values():
            csschd.safe_cross_reference(model, xref_errors)
        self._show_safe_xref_errors('csschds', xref_errors)

        xref_errors = defaultdict(list)
        for spline in model.splines.values():
            spline.safe_cross_reference(model, xref_errors)
        self._show_safe_xref_errors('splines', xref_errors)

        for aecomp in model.aecomps.values():
            aecomp.safe_cross_reference(model)

        for aelist in model.aelists.values():
            aelist.safe_cross_reference(model)

        for aeparam in model.aeparams.values():
            aeparam.safe_cross_reference(model)

        #for aestat in model.aestats):
            #aestat.safe_cross_reference(self)

        for aelinks in model.aelinks.values():
            for aelink in aelinks:
                aelink.cross_reference(model)

        xref_errors = defaultdict(list)
        for aesurf in model.aesurf.values():
            aesurf.safe_cross_reference(model, xref_errors)
        self._show_safe_xref_errors('aesurf', xref_errors)

        for aesurfs in model.aesurfs.values():
            aesurfs.safe_cross_reference(model)

        for flutter in model.flutters.values():
            flutter.safe_cross_reference(model)

        xref_errors = defaultdict(list)
        for monitor_point in model.monitor_points:
            monitor_point.safe_cross_reference(model, xref_errors)
        self._show_safe_xref_errors('monitor_points', xref_errors)

        if model.aero:
            xref_errors = defaultdict(list)
            model.aero.safe_cross_reference(model, xref_errors)
            self._show_safe_xref_errors('aero', xref_errors)
        if model.aeros:
            xref_errors = defaultdict(list)
            model.aeros.safe_cross_reference(model, xref_errors)
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
            #'AEPARM',  ## aeparams
            #'AESTAT',   ## aestats
            #'AESURF',  ## aesurfs

    def cross_reference_coordinates(self) -> None:
        """
        Links up all the coordinate cards to other coordinate cards and nodes
         - CORD1R, CORD1C, CORD1S
         - CORD2R, CORD2C, CORD2S
        """
        # CORD2x: links the rid to coordinate systems
        # CORD1x: links g1,g2,g3 to grid points
        model = self.model
        for coord in model.coords.values():
            coord.cross_reference(model)

        for coord in model.coords.values():
            coord.setup()

    def safe_cross_reference_coordinates(self) -> None:
        """
        Links up all the coordinate cards to other coordinate cards and nodes
         - CORD1R, CORD1C, CORD1S
         - CORD2R, CORD2C, CORD2S
        """
        # CORD2x: links the rid to coordinate systems
        # CORD1x: links g1,g2,g3 to grid points
        model = self.model
        xref_errors = {}
        for coord in model.coords.values():
            coord.safe_cross_reference(model, xref_errors)

        for coord in model.coords.values():
            coord.setup()

    def cross_reference_nodes(self) -> None:
        """Links the nodes to coordinate systems"""
        model = self.model
        model.log
        grdset = model.grdset
        for node in model.nodes.values():
            try:
                node.cross_reference(model, grdset)
            except Exception:
                model.log.error("Couldn't cross reference GRID.\n%s" % (str(node)))
                raise

        for point in model.points.values():
            try:
                point.cross_reference(model)
            except Exception:
                model.log.error("Couldn't cross reference POINT.\n%s" % (str(point)))
                raise

        # SPOINTs, EPOINTs don't need xref

        # GRDPNT for mass calculations
        #if model.has_key()
        #for param_key, param in self.params:
            #if

    def safe_cross_reference_nodes(self) -> None:
        """Links the nodes to coordinate systems"""
        xref_errors = defaultdict(list)
        model = self.model
        grdset = model.grdset
        for node in model.nodes.values():
            node.safe_cross_reference(model, xref_errors, grdset)

        for point in model.points.values():
            try:
                point.cross_reference(model)
            except Exception:
                model.log.error("Couldn't cross reference POINT.\n%s" % (str(point)))
                raise

        # SPOINTs, EPOINTs don't need xref

        # GRDPNT for mass calculations
        #if model.has_key()
        #for param_key, param in self.params:
            #if
        self._show_safe_xref_errors('nodes', xref_errors)


    def cross_reference_masses(self) -> None:
        """
        Links the mass to nodes, properties (and materials depending on
        the card).
        """
        model = self.model
        for mass in model.masses.values():
            try:
                mass.cross_reference(model)
            except (SyntaxError, RuntimeError, AssertionError, KeyError, ValueError) as error:  # pragma: no cover
                self._store_xref_error(error, mass)

        for prop in model.properties_mass.values():
            try:
                prop.cross_reference(model)
            except (SyntaxError, RuntimeError, AssertionError, KeyError, ValueError) as error:  # pragma: no cover
                self._store_xref_error(error, prop)

        for nsms in model.nsms.values():
            for nsm in nsms:
                nsm.cross_reference(model)

        for nsmadds in model.nsmadds.values():
            for nsmadd in nsmadds:
                nsmadd.cross_reference(model)

    def safe_cross_reference_masses(self) -> None:
        """
        Links the mass to nodes, properties (and materials depending on
        the card).
        """
        xref_errors = defaultdict(list)
        model = self.model
        for mass in model.masses.values():
            try:
                mass.safe_cross_reference(model, xref_errors)
            except (SyntaxError, RuntimeError, AssertionError, KeyError, ValueError) as error:  # pragma: no cover
                self._store_xref_error(error, mass)

        for prop in model.properties_mass.values():
            try:
                prop.safe_cross_reference(model, xref_errors)
            except (SyntaxError, RuntimeError, AssertionError, KeyError, ValueError) as error:  # pragma: no cover
                self._store_xref_error(error, prop)

    def cross_reference_elements(self) -> None:
        """
        Links the elements to nodes, properties (and materials depending on
        the card).
        """
        model = self.model
        for elem in model.elements.values():
            try:
                elem.cross_reference(model)
            except (SyntaxError, RuntimeError, AssertionError, KeyError, ValueError) as error:  # pragma: no cover
                self._store_xref_error(error, elem)

        for elem in model.masses.values():
            try:
                elem.cross_reference(model)
            except (SyntaxError, RuntimeError, AssertionError, KeyError, ValueError) as error:  # pragma: no cover
                self._store_xref_error(error, elem)

    def safe_cross_reference_elements(self) -> None:
        """
        Links the elements to nodes, properties (and materials depending on
        the card).
        """
        model = self.model
        xref_errors = defaultdict(list)
        missing_safe_xref = set()
        for elem in model.elements.values():
            if hasattr(elem, 'safe_cross_reference'):
                elem.safe_cross_reference(model, xref_errors)
            else:
                elem.cross_reference(model)
                missing_safe_xref.add(elem.type)

        for elem in model.masses.values():
            if hasattr(elem, 'safe_cross_reference'):
                elem.safe_cross_reference(model, xref_errors)
            else:
                elem.cross_reference(model)
                missing_safe_xref.add(elem.type)

        for elem in model.rigid_elements.values():
            elem.safe_cross_reference(model, xref_errors)

        self._show_safe_xref_errors('elements', xref_errors)
        if missing_safe_xref:
            model.log.warning('These cards dont support safe_xref; %s' %
                             str(list(missing_safe_xref)))

    def cross_reference_rigid_elements(self) -> None:
        model = self.model
        for elem in model.rigid_elements.values():
            try:
                elem.cross_reference(model)
            except (SyntaxError, RuntimeError, AssertionError, KeyError, ValueError) as error:  # pragma: no cover
                self._store_xref_error(error, elem)

        for elem in model.plotels.values():
            try:
                elem.cross_reference(model)
            except (SyntaxError, RuntimeError, AssertionError, KeyError, ValueError) as error:  # pragma: no cover
                self._store_xref_error(error, elem)

    def cross_reference_properties(self) -> None:
        """Links the properties to materials"""
        model = self.model
        for prop in model.properties.values():
            try:
                prop.cross_reference(model)
            except (SyntaxError, RuntimeError, AssertionError, KeyError, ValueError) as error:  # pragma: no cover
                self._store_xref_error(error, prop)

    def safe_cross_reference_properties(self) -> None:
        """Links the properties to materials"""
        xref_errors = {}
        model = self.model
        for prop in model.properties.values():
            if hasattr(prop, 'safe_cross_reference'):
                try:
                    prop.safe_cross_reference(model, xref_errors)
                except (SyntaxError, RuntimeError, AssertionError, KeyError, ValueError) as error:  # pragma: no cover
                    self._store_xref_error(error, prop)
            else:
                try:
                    prop.cross_reference(model)
                except (SyntaxError, RuntimeError, AssertionError, KeyError, ValueError) as error:  # pragma: no cover
                    self._store_xref_error(error, prop)

    def cross_reference_materials(self) -> None:
        """
        Links the materials to materials (e.g. MAT1, CREEP)
        often this is a pass statement
        """
        model = self.model
        for mat in model.materials.values():  # MAT1
            try:
                mat.cross_reference(model)
            except (SyntaxError, RuntimeError, AssertionError, KeyError, ValueError) as error:  # pragma: no cover
                self._store_xref_error(error, mat)

        for mat in model.creep_materials.values():  # CREEP
            try:
                mat.cross_reference(model)
            except (SyntaxError, RuntimeError, AssertionError, KeyError, ValueError) as error:  # pragma: no cover
                self._store_xref_error(error, mat)

        # CREEP - depends on MAT1
        data = [model.MATS1, model.MATS3, model.MATS8,
                model.MATT1, model.MATT2, model.MATT3, model.MATT4, model.MATT5,
                model.MATT8, model.MATT9, model.MATT11]
        for material_deps in data:
            for mat in material_deps.values():
                try:
                    mat.cross_reference(model)
                except (SyntaxError, RuntimeError, AssertionError, KeyError, ValueError) as error:  # pragma: no cover
                    self._store_xref_error(error, mat)

    def safe_cross_reference_materials(self) -> None:
        """
        Links the materials to materials (e.g. MAT1, CREEP)
        often this is a pass statement
        """
        model = self.model
        xref_errors = defaultdict(list)
        missing_safe_xref = set()
        for mat in model.materials.values():  # MAT1
            #if hasattr(mat, 'safe_cross_reference'):
            try:
                mat.safe_cross_reference(model, xref_errors)
            except (SyntaxError, RuntimeError, AssertionError, KeyError, ValueError) as error:  # pragma: no cover
                self._store_xref_error(error, mat)
            # else:
            #     missing_safe_xref.add(mat.type)
            #     try:
            #         mat.cross_reference(model)
            #     except (SyntaxError, RuntimeError, AssertionError, KeyError, ValueError) as error:  # pragma: no cover
            #         self._store_xref_error(error, mat)

        for mat in model.creep_materials.values():  # CREEP
            try:
                mat.cross_reference(model)
            except (SyntaxError, RuntimeError, AssertionError, KeyError, ValueError) as error:  # pragma: no cover
                self._store_xref_error(error, mat)

        # CREEP - depends on MAT1
        data = [model.MATS1, model.MATS3, model.MATS8,
                model.MATT1, model.MATT2, model.MATT3, model.MATT4, model.MATT5,
                model.MATT8, model.MATT9, model.MATT11]
        for material_deps in data:
            for mat in material_deps.values():
                #if hasattr(mat, 'safe_cross_reference'):
                try:
                    mat.safe_cross_reference(model, xref_errors)
                except (SyntaxError, RuntimeError, AssertionError, KeyError, ValueError) as error:  # pragma: no cover
                    self._store_xref_error(error, mat)
                # else:
                #     missing_safe_xref.add(mat.type)
                #     try:
                #         mat.cross_reference(model)
                #     except (SyntaxError, RuntimeError, AssertionError, KeyError, ValueError) as error:  # pragma: no cover
                #         self._store_xref_error(error, mat)

        self._show_safe_xref_errors('materials', xref_errors)
        if missing_safe_xref:
            model.log.warning('These cards dont support safe_xref; %s' %
                             str(list(missing_safe_xref)))

    def _show_safe_xref_errors(self, elements_word: str, xref_errors: bool) -> None:
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
            self.model.log.warning(msg.rstrip())

    def cross_reference_loads(self) -> None:
        """Links the loads to nodes, coordinate systems, and other loads."""
        model = self.model
        for (unused_lid, load_combinations) in model.load_combinations.items():
            for load_combination in load_combinations:
                try:
                    load_combination.cross_reference(model)
                except (SyntaxError, RuntimeError, AssertionError, KeyError, ValueError) as error:  # pragma: no cover
                    self._store_xref_error(error, load_combination)

        for (unused_lid, loads) in model.loads.items():
            for load in loads:
                try:
                    load.cross_reference(model)
                except (SyntaxError, RuntimeError, AssertionError, KeyError, ValueError) as error:  # pragma: no cover
                    self._store_xref_error(error, load)

        for (unused_lid, sid) in model.dloads.items():
            for load in sid:
                #self.log.debug("  dloadi load=%s" % (load))
                try:
                    load.cross_reference(model)
                except (SyntaxError, RuntimeError, AssertionError, KeyError, ValueError) as error:  # pragma: no cover
                    self._ixref_errors += 1
                    var = traceback.format_exception_only(type(error), error)
                    self._stored_xref_errors.append((load, var))
                    if self._ixref_errors > self._nxref_errors:
                        self.pop_xref_errors()

        for unused_lid, sid in model.dload_entries.items():
            for load in sid:
                #self.log.debug("  dloadi load=%s" % (load))
                try:
                    load.cross_reference(model)
                except (SyntaxError, RuntimeError, AssertionError, KeyError, ValueError) as error:  # pragma: no cover
                    #raise
                    self._store_xref_error(error, load)

        for unused_key, darea in model.dareas.items():
            try:
                darea.cross_reference(model)
            except (SyntaxError, RuntimeError, AssertionError, KeyError, ValueError) as error:  # pragma: no cover
                self._store_xref_error(error, darea)

        for unused_key, tic in model.tics.items():
            try:
                tic.cross_reference(model)
            except (SyntaxError, RuntimeError, AssertionError, KeyError, ValueError) as error:  # pragma: no cover
                self._store_xref_error(error, tic)

        for unused_key, dphase in model.dphases.items():
            try:
                dphase.cross_reference(model)
            except (SyntaxError, RuntimeError, AssertionError, KeyError, ValueError) as error:  # pragma: no cover
                self._store_xref_error(error, dphase)

    def safe_cross_reference_loads(self) -> None:
        """
        Links the loads to nodes, coordinate systems, and other loads.
        """
        model = self.model
        xref_errors = defaultdict(list)
        for unused_lid, load_combinations in model.load_combinations.items():
            for load_combination in load_combinations:
                try:
                    load_combination.safe_cross_reference(model, xref_errors)
                except TypeError:  # pragma: no cover
                    print(load_combination)
                    raise
        self._show_safe_xref_errors('loads', xref_errors)

        for unused_lid, loads in model.loads.items():
            for load in loads:
                load.safe_cross_reference(model, xref_errors)
        self._show_safe_xref_errors('loads', xref_errors)

        for unused_lid, sid in model.dloads.items():
            for load in sid:
                load.safe_cross_reference(model, xref_errors)

        for unused_lid, sid in model.dload_entries.items():
            for load in sid:
                load.safe_cross_reference(model, xref_errors)

        for unused_key, darea in model.dareas.items():
            darea.safe_cross_reference(model, xref_errors)

        for unused_key, dphase in model.dphases.items():
            dphase.safe_cross_reference(model, xref_errors)

        for unused_key, tic in model.tics.items():
            tic.safe_cross_reference(model, xref_errors)

    def cross_reference_optimization(self) -> None:
        """cross references the optimization objects"""
        model = self.model
        remove_missing_optimization = True
        dconstrs_to_remove = []

        for unused_key, deqatn in model.dequations.items():
            deqatn.cross_reference(model)
        for unused_key, dresp in model.dresps.items():
            dresp.cross_reference(model)

        for key, dconstrs in model.dconstrs.items():
            for i, dconstr in enumerate(dconstrs):
                try:
                    dconstr.cross_reference(model)
                except:
                    if not remove_missing_optimization:
                        raise
                    dconstrs_to_remove.append((key, i))

        for unused_key, dvcrel in model.dvcrels.items():
            dvcrel.cross_reference(model)
        for unused_key, dvmrel in model.dvmrels.items():
            dvmrel.cross_reference(model)
        for unused_key, dvprel in model.dvprels.items():
            dvprel.cross_reference(model)
        for unused_key, desvar in model.desvars.items():
            desvar.cross_reference(model)
        for unused_key, topvar in model.topvar.items():
            topvar.cross_reference(model)

        for key, i in dconstrs_to_remove:
            del model.dconstrs[key][i]

    def safe_cross_reference_optimization(self) -> None:
        """cross references the optimization objects"""
        #self._cross_reference_optimization()
        #return
        model = self.model
        xref_errors = defaultdict(list)
        for unused_key, deqatn in model.dequations.items():
            deqatn = cast(DEQATN, deqatn)
            deqatn.safe_cross_reference(model, xref_errors)

        for unused_key, dresp in model.dresps.items():
            dresp.safe_cross_reference(model, xref_errors)

        for unused_key, dconstrs in model.dconstrs.items():
            for dconstr in dconstrs:
                dconstr = cast(DCONSTR, dconstr)
                if hasattr(dconstr, 'safe_cross_reference'):
                    dconstr.safe_cross_reference(model, xref_errors)
                else:  # pragma: no cover
                    dconstr.cross_reference(model)

        for unused_key, dvcrel in model.dvcrels.items():
            if hasattr(dvcrel, 'safe_cross_reference'):
                dvcrel.safe_cross_reference(model, xref_errors)
            else:  # pragma: no cover
                dvcrel.cross_reference(model)

        for unused_key, dvmrel in model.dvmrels.items():
            if hasattr(dvmrel, 'safe_cross_reference'):
                dvmrel.safe_cross_reference(model, xref_errors)
            else:  # pragma: no cover
                dvmrel.cross_reference(model)

        for unused_key, dvprel in model.dvprels.items():
            if hasattr(dvprel, 'safe_cross_reference'):
                dvprel.safe_cross_reference(model, xref_errors)
            else:  # pragma: no cover
                dvprel.cross_reference(model)

        for unused_key, desvar in model.desvars.items():
            desvar = cast(DESVAR, desvar)
            desvar.safe_cross_reference(model, xref_errors)

        for unused_key, topvar in model.topvar.items():
            topvar.safe_cross_reference(model, xref_errors)

    def cross_reference_bolts(self) -> None:
        model = self.model
        for bolt_dict in (model.bolt, model.boltfor, model.boltseq,
                          model.boltld, model.boltfrc):
            for bolt in bolt_dict.values():
                bolt.cross_reference(model)

    def cross_reference_contact(self) -> None:
        """cross references the contact objects"""
        model = self.model
        for blseg in model.blseg.values():
            blseg.cross_reference(model)
        for bconp in model.bconp.values():
            bconp.cross_reference(model)

        # bgset
        # bctset
        #for bgadd in self.bgadds.values():
            #bgadd.cross_reference(self)
        #for bctadd in self.bctadds.values():
            #bctadd.cross_reference(self)

    def safe_cross_reference_contact(self) -> None:
        """cross references the contact objects"""
        self.cross_reference_contact()

    def cross_reference_superelements(self) -> None:
        """cross references the superelement objects"""
        model = self.model
        for unused_seid, csuper in model.csuper.items():
            csuper.cross_reference(model)
        for unused_seid, csupext in model.csupext.items():
            csupext.cross_reference(model)

        for unused_seid, sebulk in model.sebulk.items():
            sebulk.cross_reference(model)
        for unused_seid, sebndry in model.sebndry.items():
            sebndry.cross_reference(model)
        for unused_seid, seconct in model.seconct.items():
            seconct.cross_reference(model)
        for unused_seid, seelt in model.seelt.items():
            seelt.cross_reference(model)
        for unused_seid, seexcld in model.seexcld.items():
            seexcld.cross_reference(model)
        for unused_seid, selabel in model.selabel.items():
            selabel.cross_reference(model)
        for unused_seid, seloc in model.seloc.items():
            seloc.cross_reference(model)
        for unused_seid, seload in model.seload.items():
            seload.cross_reference(model)
        for unused_seid, sempln in model.sempln.items():
            sempln.cross_reference(model)
        for unused_seid, setree in model.setree.items():
            setree.cross_reference(model)

        #'senqset',
        #'se_sets', 'se_usets',

    def safe_cross_reference_superelements(
            self, create_superelement_geometry: bool=False) -> None:
        model = self.model
        xref_errors = {}
        seloc_missing = []
        for seid, seloc in model.seloc.items():
            super_key = ('SUPER', seid, '')
            if super_key in model.superelement_models:
                superelement = model.superelement_models[super_key]
                seloc.safe_cross_reference(model, xref_errors)
                #seloc.transform(self)
            else:
                seloc_missing.append(seid)

        try:
            for unused_seid, sempln in sorted(model.sempln.items()):
                sempln.safe_cross_reference(model, xref_errors)
            for unused_seid, csuper in model.csuper.items():
                csuper.safe_cross_reference(model, xref_errors)
            for unused_seid, csupext in model.csupext.items():
                csupext.safe_cross_reference(model, xref_errors)

            if model.sebulk and create_superelement_geometry:
                #print('sebulk...')
                import os
                # we have to create the superelement in order to transform it...
                for seid, sebulk in model.sebulk.items():
                    super_key = ('SUPER', seid, '')
                    super_filename = f'super_{seid:d}.bdf'
                    if os.path.exists(super_filename):
                        os.remove(super_filename)
                    #print(sebulk)
                    rseid = sebulk.rseid
                    sebulk.safe_cross_reference(model, xref_errors)
                    mirror_model = model._create_superelement_from_sebulk(sebulk, seid, rseid)
                    if mirror_model is None:
                        continue
                    model.log.debug(f'made superelement {seid:d}')
                    model.superelement_models[super_key] = mirror_model
                    mirror_model.write_bdf(super_filename)
            for unused_seid, sebndry in model.sebndry.items():
                sebndry.safe_cross_reference(model, xref_errors)
            for unused_seid, seconct in model.seconct.items():
                seconct.safe_cross_reference(model, xref_errors)
            for unused_seid, seelt in model.seelt.items():
                seelt.safe_cross_reference(model, xref_errors)
            for unused_seid, seexcld in model.seexcld.items():
                seexcld.safe_cross_reference(model, xref_errors)
            for unused_seid, selabel in model.selabel.items():
                selabel.safe_cross_reference(model, xref_errors)
            for seid in seloc_missing:
                seloc = model.seloc[seid]
                seloc.safe_cross_reference(model, xref_errors)
            for unused_seid, seload in model.seload.items():
                seload.safe_cross_reference(model, xref_errors)
            for unused_seid, setree in model.setree.items():
                setree.safe_cross_reference(model, xref_errors)
        except KeyError:
            if not create_superelement_geometry:
                raise
            model.write_bdf('superelement_xref.bdf')
            model.log.error('check superelement_xref.bdf')
            raise

    def cross_reference_nodes_with_elements(self) -> None:
        """Links the nodes to all connected elements"""
        nodes: dict[int, list[GRID]] = defaultdict(list)
        model = self.model
        for element in model.elements.values():
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
        for node in model.nodes.values():
            node.elements_ref = nodes[node.nid]


def _check_caero_box_overlap(model: BDF) -> int:
    """
    only support CAERO1
    check that aerobox ids don't overlap
    """
    ncaeros = len(model.caeros)
    if ncaeros == 0:
        return 0
    # we don't need to check the ncaeros=1 case
    i = 0
    naeroboxes = 0
    #print('eid, naeroboxesi')
    min_maxs = np.zeros((ncaeros, 2), dtype='int32')
    for eid, caero in sorted(model.caeros.items()):
        min_maxs[i, :] = caero.min_max_eid
        npointsi, naeroboxesi = caero.get_panel_npoints_nelements()
        naeroboxes += naeroboxesi
        #print(caero)
        #print(eid, naeroboxesi)
        i += 1
    naeroboxes_per_caeros = naeroboxes
    mins = min_maxs[:, 0]
    maxs = min_maxs[:, 1]
    delta = maxs - mins
    #print('delta:')
    #print(f'delta = {delta}')
    assert delta.min() > 0, delta
    naeroboxes = delta.sum()
    assert naeroboxes == naeroboxes_per_caeros, f'naeroboxes={naeroboxes} != naeroboxes_per_caeros={naeroboxes_per_caeros}'

    isort = np.argsort(min_maxs.ravel())
    expected = np.arange(ncaeros * 2, dtype='int32')
    if not np.array_equal(isort, expected):
        msg = 'CAERO element ids are inconsistent\n'
        msg += 'isort = %s' % str(isort)
        raise RuntimeError(msg)
    return naeroboxes


def _add_file(model: BDF, obj: BaseCard) -> str:
    """adds the tag of the active file"""
    msg = ''
    if obj is not None and hasattr(obj, 'ifile'):
        filename = model.active_filenames[obj.ifile]
        # tag = f' ({obj.type}={obj.eid:d})' if hasattr(obj, 'eid') else ''
        # msg = f'file[{obj.ifile:d}]={filename}{tag}\n'
        #msg += f'active: {model.active_filenames}\n'
        msg = f'file={filename}\n'
    return msg
