"""
This file defines:
  - WriteMesh

"""
import sys
from io import StringIO, IOBase

from numpy import array, unique, concatenate, intersect1d, where

from pyNastran.bdf.bdf_interface.utils import print_filename
from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.field_writer_16 import print_card_16
from pyNastran.dev.bdf_vectorized.bdf_interface2.attributes import BDFAttributes


class WriteMesh(BDFAttributes):
    """
    Defines methods for writing cards

    Major methods:
      - model.write_bdf(...)
    """
    def __init__(self):
        """creates methods for writing cards"""
        BDFAttributes.__init__(self)

    def get_encoding(self, encoding=None):
        if encoding is not None:
            pass
        else:
            encoding = self._encoding
            if encoding is None:
                encoding = sys.getdefaultencoding()
        return encoding

    def _output_helper(self, out_filename, interspersed, size, is_double):
        """
        Performs type checking on the write_bdf inputs
        """
        if out_filename is None:
            from pyNastran.utils.gui_io import save_file_dialog
            wildcard_wx = "Nastran BDF (*.bdf; *.dat; *.nas; *.pch)|" \
                "*.bdf;*.dat;*.nas;*.pch|" \
                "All files (*.*)|*.*"
            wildcard_qt = "Nastran BDF (*.bdf *.dat *.nas *.pch);;All files (*)"
            title = 'Save BDF/DAT/PCH'
            out_filename = save_file_dialog(title, wildcard_wx, wildcard_qt)
            assert out_filename is not None, out_filename

        if not(hasattr(out_filename, 'read') and hasattr(out_filename, 'write')) or isinstance(out_filename, IOBase):
            return out_filename
        elif not isinstance(out_filename, str):
            msg = 'out_filename=%r must be a string; type=%s' % (
                out_filename, type(out_filename))
            raise TypeError(msg)

        if size == 8:
            assert is_double is False, 'is_double=%r' % is_double
        elif size == 16:
            assert is_double in [True, False], 'is_double=%r' % is_double
        else:
            assert size in [8, 16], size

        assert isinstance(interspersed, bool)
        fname = print_filename(out_filename, self._relpath)
        self.log.debug("***writing %s" % fname)
        return out_filename

    def write_bdf(self, out_filename=None, encoding=None,
                  size=8, is_double=False,
                  interspersed=False, enddata=None, close=True):
        """
        Writes the BDF.

        Parameters
        ----------
        out_filename : varies; default=None
            str        - the name to call the output bdf
            file       - a file object
            StringIO() - a StringIO object
            None       - pops a dialog
        encoding : str; default=None -> system specified encoding
            the unicode encoding
            latin1, and utf8 are generally good options
        size : int; {8, 16}
            the field size
        is_double : bool; default=False
            False : small field
            True : large field
        interspersed : bool; default=False
            Writes a bdf with properties & elements
            interspersed like how Patran writes the bdf.  This takes
            slightly longer than if interspersed=False, but makes it
            much easier to compare to a Patran-formatted bdf and is
            more clear.
        enddata : bool; default=None
            bool - enable/disable writing ENDDATA
            None - depends on input BDF
        close : bool; default=True
            should the output file be closed
        """
        out_filename = self._output_helper(out_filename,
                                           interspersed, size, is_double)
        self.log.debug('---starting BDF.write_bdf of %s---' % out_filename)
        encoding = self.get_encoding(encoding)
        #assert encoding.lower() in ['ascii', 'latin1', 'utf8'], encoding

        if hasattr(out_filename, 'read') and hasattr(out_filename, 'write'):
            bdf_file = out_filename
        else:
            bdf_file = open(out_filename, 'w', encoding=encoding)
        self._write_header(bdf_file, encoding)
        self._write_params(bdf_file, size, is_double)
        self._write_nodes(bdf_file, size, is_double)

        self.write_elements_properties(bdf_file, size, is_double, interspersed)
        self._write_materials(bdf_file, size, is_double)
        self._write_rigid_elements(bdf_file, size, is_double) # split out for write_bdf_symmetric
        self._write_aero(bdf_file, size, is_double) # split out for write_bdf_symmetric
        self._write_common(bdf_file, size, is_double)

        if (enddata is None and 'ENDDATA' in self.card_count) or enddata:
            bdf_file.write('ENDDATA\n')
        if close:
            bdf_file.close()

    def write_bdf_symmetric(self, out_filename=None, encoding=None,
                            size=8, is_double=False,
                            enddata=None, close=True, plane='xz'):
        """
        Writes the BDF.

        Parameters
        ----------
        out_filename : varies; default=None
            str        - the name to call the output bdf
            file       - a file object
            StringIO() - a StringIO object
            None       - pops a dialog
        encoding : str; default=None -> system specified encoding
            the unicode encoding
            latin1, and utf8 are generally good options
        size : int; {8, 16}
            the field size
        is_double : bool; default=False
            False : small field
            True : large field
        enddata : bool; default=None
            bool - enable/disable writing ENDDATA
            None - depends on input BDF
        close : bool; default=True
            should the output file be closed
        plane : str; {'xy', 'yz', 'xz'}; default='xz'
            the plane to mirror about
        """
        raise NotImplementedError()

    def _write_header(self, bdf_file, encoding):
        """
        Writes the executive and case control decks.
        """
        if self.punch is None:
            # writing a mesh without using read_bdf
            if self.executive_control_lines or self.case_control_deck:
                self.punch = False
            else:
                self.punch = True

        if self.nastran_format:
            bdf_file.write('$pyNastran: version=%s\n' % self.nastran_format)
            bdf_file.write('$pyNastran: punch=%s\n' % self.punch)
            bdf_file.write('$pyNastran: encoding=%s\n' % encoding)
            bdf_file.write('$pyNastran: nnodes=%s\n' % self.grid.n)
            bdf_file.write('$pyNastran: nelements=%s\n' % len(self.elements))

        if not self.punch:
            self._write_executive_control_deck(bdf_file)
            self._write_case_control_deck(bdf_file)

    def _write_executive_control_deck(self, bdf_file):
        """
        Writes the executive control deck.
        """
        if self.executive_control_lines:
            msg = '$EXECUTIVE CONTROL DECK\n'
            if self.sol == 600:
                new_sol = 'SOL 600,%s' % self.sol_method
            else:
                new_sol = 'SOL %s' % self.sol

            if self.sol_iline is not None:
                self.executive_control_lines[self.sol_iline] = new_sol

            for line in self.executive_control_lines:
                msg += line + '\n'
            bdf_file.write(msg)

    def _write_case_control_deck(self, bdf_file):
        """
        Writes the Case Control Deck.
        """
        if self.case_control_deck:
            msg = '$CASE CONTROL DECK\n'
            msg += str(self.case_control_deck)
            assert 'BEGIN BULK' in msg, msg
            bdf_file.write(''.join(msg))

    def _write_elements_properties(self, bdf_file, size):
        """
        Writes the elements and properties in and interspersed order
        """
        #return self._write_elements_properties2(f, size)
        #msg = []
        #missing_properties = []
        if self.properties:
            bdf_file.write('$ELEMENTS_WITH_PROPERTIES\n')

        #eids_written = []
        #pids = sorted(self.properties.keys())

        ptypes = [
            self.properties_shell.pshell,
            self.properties_shell.pcomp,
            self.pshear,
            self.prod,
            self.properties_solid.psolid,

            #self.properties_bar.pbar,
            #self.properties_bar.pbarl,
            #self.properties_beam.pbeam,
            #self.properties_beam.pbeaml,
        ]

        n = 0
        pids_all = None  # the actual properties
        for t in ptypes:
            if t.n and n == 0:
                pids_all = t.property_id
                n = 1
            elif t.n:
                self.log.debug(pids_all)
                self.log.debug(t.property_id)
                try:
                    pids_all = concatenate(pids_all, t.property_id)
                except ValueError:
                    pids_all = array(list(pids_all) + list(t.property_id))

        etypes = (self.elements_shell._get_types() +
                  self.elements_solid._get_types() +
                  [self.crod, self.cshear])

        #pids_set = None
        if pids_all is None:
            bdf_file.write('$MISSING_ELEMENTS because there are no properties\n')
            for t in etypes:
                #print("t.type =", t.type)
                t.write_card(bdf_file, size=size)
            return

        # there are properties
        pids_set = set(list(pids_all))

        n = 0
        pids = None
        for t in etypes:
            #print("t.type =", t.type)
            if t.n and n == 0:
                eids = t.element_id
                pids = t.property_id
                n = 1
            elif t.n:
                try:
                    eids = concatenate(eids, t.element_id)
                #except AttributeError:
                    #eids = array(list(eids) + list(t.element_id))
                except TypeError:
                    #print eids
                    #print t.element_id
                    eids = array(list(eids) + list(t.element_id))
                except ValueError:
                    #print eids
                    #print t.element_id
                    eids = array(list(eids) + list(t.element_id))

                try:
                    pids = concatenate(pids, t.property_id)
                except AttributeError:
                    pids = array(list(pids) + list(t.property_id))
                except TypeError:
                    pids = array(list(pids) + list(t.property_id))
                except ValueError:
                    pids = array(list(pids) + list(t.property_id))
            #else:
                #print t.type

        elements_by_pid = {}
        if pids is not None:
            pids_unique = unique(pids)
            self.log.debug("pids_unique = %s" % pids_unique)
            pids_unique.sort()
            if len(pids_unique) > 0:
                bdf_file.write('$ELEMENTS_WITH_PROPERTIES\n')

            for pid in pids_all:
                i = where(pid == pids)[0]
                eids2 = eids[i]

                for t in ptypes:
                    if t.n and pid in t.property_id:
                        self.log.debug("prop.type = %s" % t.type)
                        t.write_card(bdf_file, size=size, property_ids=[pid])
                        pids_set.remove(pid)
                n = 0
                for t in etypes:
                    if not t.n:
                        continue
                    eids3 = intersect1d(t.element_id, eids2, assume_unique=False)
                    #print("eids3[pid=%s]" %(pid), eids3)
                    if n == 0 and len(eids3):
                        elements_by_pid[pid] = eids3
                        n = 1
                    elif len(eids3):
                        try:
                            c = concatenate(elements_by_pid[pid], eids3)
                        except TypeError:
                            c = array(list(elements_by_pid[pid]) + list(eids3))
                        except ValueError:
                            c = array(list(elements_by_pid[pid]) + list(eids3))
                        elements_by_pid[pid] = c
                    else:
                        continue
                    try:
                        t.write_card(bdf_file, size=size, element_ids=eids3)
                    except TypeError:
                        print("t.type = %s" % t.type)
                        raise
                    del eids3
            #for pid, elements in elements_by_pid.items():
                #print("pid=%s n=%s" % (pid, len(elements)))
            #print elements_by_pid

        # missing properties
        if pids_set:
            pids_list = list(pids_set)
            bdf_file.write('$UNASSOCIATED_PROPERTIES\n')
            for pid in pids_list:
                for prop in ptypes:
                    if prop.n and pid in prop.property_id:
                        prop.write_card(bdf_file, size=size, property_ids=[pid])

        #.. todo:: finish...
        bdf_file.write('$UNASSOCIATED_ELEMENTS\n')
        # missing elements...

    def write_elements_properties(self, bdf_file, size, is_double, interspersed):
        self.elements.write_card(bdf_file, size=size, is_double=is_double,
                                 include_properties=True, interspersed=interspersed)

    def _write_aero(self, bdf_file, size=8, is_double=False):
        """Writes the aero cards"""
        if self.caeros or self.paeros or self.monitor_points or self.splines:
            msg = ['$AERO\n']
            for (unused_id, caero) in sorted(self.caeros.items()):
                msg.append(caero.write_card(size, is_double))
            for (unused_id, paero) in sorted(self.paeros.items()):
                msg.append(paero.write_card(size, is_double))
            for (unused_id, spline) in sorted(self.splines.items()):
                msg.append(spline.write_card(size, is_double))
            for monitor_point in self.monitor_points:
                msg.append(monitor_point.write_card(size, is_double))
            bdf_file.write(''.join(msg))

    def _write_aero_control(self, bdf_file, size=8, is_double=False):
        """Writes the aero control surface cards"""
        if(self.aecomps or self.aefacts or self.aeparams or self.aelinks or
           self.aelists or self.aestats or self.aesurf or self.aesurfs):
            msg = ['$AERO CONTROL SURFACES\n']
            for unused_id, aelinks in sorted(self.aelinks.items()):
                for aelink in aelinks:
                    msg.append(aelink.write_card(size, is_double))

            for unused_id, aecomp in sorted(self.aecomps.items()):
                msg.append(aecomp.write_card(size, is_double))
            for unused_id, aeparam in sorted(self.aeparams.items()):
                msg.append(aeparam.write_card(size, is_double))
            for unused_id, aestat in sorted(self.aestats.items()):
                msg.append(aestat.write_card(size, is_double))

            for unused_id, aelist in sorted(self.aelists.items()):
                msg.append(aelist.write_card(size, is_double))
            for unused_id, aesurf in sorted(self.aesurf.items()):
                msg.append(aesurf.write_card(size, is_double))
            for unused_id, aesurfs in sorted(self.aesurfs.items()):
                msg.append(aesurfs.write_card(size, is_double))
            for unused_id, aefact in sorted(self.aefacts.items()):
                msg.append(aefact.write_card(size, is_double))
            bdf_file.write(''.join(msg))

    def _write_static_aero(self, bdf_file, size=8, is_double=False):
        """Writes the static aero cards"""
        if self.aeros or self.trims or self.divergs:
            msg = ['$STATIC AERO\n']
            # static aero
            if self.aeros:
                msg.append(self.aeros.write_card(size, is_double))
            for (unused_id, trim) in sorted(self.trims.items()):
                msg.append(trim.write_card(size, is_double))
            for (unused_id, diverg) in sorted(self.divergs.items()):
                msg.append(diverg.write_card(size, is_double))
            bdf_file.write(''.join(msg))

    def _find_aero_location(self):
        """Determines where the AERO card should be written"""
        write_aero_in_flutter = False
        write_aero_in_gust = False
        if self.aero:
            if self.flfacts or self.flutters or self.mkaeros:
                write_aero_in_flutter = True
            elif self.gusts:
                write_aero_in_gust = True
            else:
                # an AERO card exists, but no FLUTTER, FLFACT, MKAEROx or GUST card
                write_aero_in_flutter = True
        return write_aero_in_flutter, write_aero_in_gust

    def _write_flutter(self, bdf_file, size=8, is_double=False, write_aero_in_flutter=True):
        """Writes the flutter cards"""
        if (write_aero_in_flutter and self.aero) or self.flfacts or self.flutters or self.mkaeros:
            msg = ['$FLUTTER\n']
            if write_aero_in_flutter:
                msg.append(self.aero.write_card(size, is_double))
            for (unused_id, flutter) in sorted(self.flutters.items()):
                msg.append(flutter.write_card(size, is_double))
            for (unused_id, flfact) in sorted(self.flfacts.items()):
                msg.append(flfact.write_card(size, is_double))
            for mkaero in self.mkaeros:
                msg.append(mkaero.write_card(size, is_double))
            bdf_file.write(''.join(msg))

    def _write_gust(self, bdf_file, size=8, is_double=False, write_aero_in_gust=True):
        """Writes the gust cards"""
        if (write_aero_in_gust and self.aero) or self.gusts:
            msg = ['$GUST\n']
            if write_aero_in_gust:
                for (unused_id, aero) in sorted(self.aero.items()):
                    msg.append(aero.write_card(size, is_double))
            for (unused_id, gust) in sorted(self.gusts.items()):
                msg.append(gust.write_card(size, is_double))
            bdf_file.write(''.join(msg))

    def _write_common(self, bdf_file, size=8, is_double=False):
        """
        Write the common outputs so none get missed...

        Parameters
        ----------
        bdf_file : file
            the file object
        size : int (default=8)
            the field width
        is_double : bool (default=False)
            is this double precision

        Returns
        -------
        msg : str
            part of the bdf
        """
        self._write_rigid_elements(bdf_file, size, is_double)
        self._write_dmigs(bdf_file, size, is_double)
        self._write_loads(bdf_file, size, is_double)

        self._write_dynamic(bdf_file, size, is_double)
        self._write_aero_control(bdf_file, size, is_double)
        self._write_static_aero(bdf_file, size, is_double)

        write_aero_in_flutter, write_aero_in_gust = self._find_aero_location()
        self._write_flutter(bdf_file, size, is_double, write_aero_in_flutter)
        self._write_gust(bdf_file, size, is_double, write_aero_in_gust)

        #self._write_thermal(bdf_file, size, is_double)
        #self._write_thermal_materials(bdf_file, size, is_double)

        self._write_constraints(bdf_file, size, is_double)
        self._write_optimization(bdf_file, size, is_double)
        #self._write_tables(bdf_file, size, is_double)
        self._write_sets(bdf_file, size, is_double)
        self._write_superelements(bdf_file, size, is_double)
        self._write_contact(bdf_file, size, is_double)
        self._write_rejects(bdf_file, size, is_double)
        self._write_coords(bdf_file, size, is_double)

    def _write_constraints(self, bdf_file, size=8, is_double=False):
        """Writes the constraint cards sorted by ID"""
        spcs = [self.spcadd, self.spc, self.spcd, self.spc1]
        mpcs = [self.mpcadd, self.mpc]
        self._write_constraints_spc_mpc(bdf_file, size, spcs, 'SPC')
        self._write_constraints_spc_mpc(bdf_file, size, mpcs, 'MPC')

    def _write_constraints_spc_mpc(self, bdf_file, size, constraint_dicts, constraint_type):
        interspersed = False
        if interspersed:
            raise NotImplementedError()
        else:
            ids = []
            for constraint_dict in constraint_dicts:
                ids += constraint_dict.keys()
            #self.log.debug(ids)
            ids = unique(ids)
            ids.sort()
            self.log.debug('%s ids = %s' % (constraint_type, ids))
            if len(ids) > 0:
                bdf_file.write('$CONSTRAINTS\n')
                for idi in ids:
                    for constraint_dict in constraint_dicts:
                        for constraint_id, constraint in sorted(constraint_dict.items()):
                            if idi == constraint_id:
                                #self.log.debug('writing %s' % constraint.type)
                                constraint.write_card(bdf_file, size=size)

    def _write_contact(self, bdf_file, size=8, is_double=False):
        """Writes the contact cards sorted by ID"""
        is_contact = (self.bcrparas or self.bctadds or self.bctparas
                      or self.bctsets or self.bsurf or self.bsurfs)
        if is_contact:
            msg = ['$CONTACT\n']
            for (unused_id, bcrpara) in sorted(self.bcrparas.items()):
                msg.append(bcrpara.write_card(size, is_double))
            for (unused_id, bctadds) in sorted(self.bctadds.items()):
                msg.append(bctadds.write_card(size, is_double))
            for (unused_id, bctpara) in sorted(self.bctparas.items()):
                msg.append(bctpara.write_card(size, is_double))

            for (unused_id, bctset) in sorted(self.bctsets.items()):
                msg.append(bctset.write_card(size, is_double))
            for (unused_id, bsurfi) in sorted(self.bsurf.items()):
                msg.append(bsurfi.write_card(size, is_double))
            for (unused_id, bsurfsi) in sorted(self.bsurfs.items()):
                msg.append(bsurfsi.write_card(size, is_double))
            bdf_file.write(''.join(msg))

    def _write_coords(self, bdf_file, size=8, is_double=False):
        """Writes the coordinate cards in a sorted order"""
        self.coords.write_card(bdf_file, size, is_double)

    def _write_dmigs(self, bdf_file, size=8, is_double=False):
        """
        Writes the DMIG cards

        Parameters
        ----------
        size : int
            large field (16) or small field (8)

        Returns
        -------
        msg : str
            string representation of the DMIGs
        """
        msg = []
        for unused_name, dmig in sorted(self.dmigs.items()):
            msg.append(dmig.write_card(size, is_double))
        for unused_name, dmi in sorted(self.dmis.items()):
            msg.append(dmi.write_card(size, is_double))
        for unused_name, dmij in sorted(self.dmijs.items()):
            msg.append(dmij.write_card(size, is_double))
        for unused_name, dmiji in sorted(self.dmijis.items()):
            msg.append(dmiji.write_card(size, is_double))
        for unused_name, dmik in sorted(self.dmiks.items()):
            msg.append(dmik.write_card(size, is_double))
        bdf_file.write(''.join(msg))

    def _write_dynamic(self, bdf_file, size=8, is_double=False):
        """Writes the dynamic cards sorted by ID"""
        is_dynamic = (self.nlparms or self.frequencies or self.methods or
                      self.cMethods or self.tsteps or self.tstepnls)
        if is_dynamic:
            msg = ['$DYNAMIC\n']
            for unused_id, method in sorted(self.methods.items()):
                msg.append(method.write_card(size, is_double))
            for unused_id, cmethod in sorted(self.cMethods.items()):
                msg.append(cmethod.write_card(size, is_double))
            #for (unused_id, darea) in sorted(self.dareas.items()):
                #msg.append(darea.write_card(size, is_double))
            for unused_id, nlparm in sorted(self.nlparms.items()):
                msg.append(nlparm.write_card(size, is_double))
            for unused_id, nlpci in sorted(self.nlpcis.items()):
                msg.append(nlpci.write_card(size, is_double))
            for unused_id, tstep in sorted(self.tsteps.items()):
                msg.append(tstep.write_card(size, is_double))
            for unused_id, tstepnl in sorted(self.tstepnls.items()):
                msg.append(tstepnl.write_card(size, is_double))
            for unused_id, freq in sorted(self.frequencies.items()):
                msg.append(freq.write_card(size, is_double))
            bdf_file.write(''.join(msg))

    def _write_nonlinear(self, bdf_file, size):
        for unused_key, card in sorted(self.nlparms.items()):
            card.write_card(bdf_file, size)
        for unused_key, card in sorted(self.nlpcis.items()):
            card.write_card(bdf_file, size)
        #self.tables1.write_card(bdf_file, size)

    def _write_loads(self, bdf_file, size=8, is_double=False):
        """Writes the load cards sorted by ID"""
        self.loads.write_card(bdf_file, size, is_double, sort_data=False)
        self.temps.write_card(bdf_file, size, is_double, sort_data=False)

    def _write_masses(self, bdf_file, size, is_double):
        """Writes the mass cards sorted by ID"""
        pass

    def _write_materials(self, bdf_file, size=8, is_double=False):
        """Writes the materials in a sorted order"""
        self.materials.write_card(bdf_file, size, is_double)

    def _write_nodes(self, bdf_file, size=8, is_double=False):
        """
        Writes the NODE-type cards
        """
        self.grdset.write_card(bdf_file, size=size, is_double=is_double)
        self.grid.write_card(bdf_file, size=size, is_double=is_double)
        self.point.write_card(bdf_file, size=size, is_double=is_double)
        self.epoint.write_card(bdf_file, size=size, is_double=is_double)
        self.spoint.write_card(bdf_file, size=size, is_double=is_double)
        self.pointax.write_card(bdf_file, size=size, is_double=is_double)
    def _write_nodes_symmetric(self, bdf_file, size=8, is_double=False, plane='xz'):
        """
        Writes the NODE-type cards

        .. warning:: doesn't consider coordinate systems;
                      it could, but you'd need 20 new coordinate systems
        .. warning:: doesn't mirror SPOINTs, EPOINTs
        """
        raise NotImplementedError()


    def _write_optimization(self, bdf_file, size=8, is_double=False):
        """Writes the optimization cards sorted by ID"""
        is_optimization = (self.dconstrs or self.desvars or self.ddvals or
                           self.dresps or
                           self.dvprels or self.dvmrels or self.dvcrels or self.doptprm or
                           self.dlinks or
                           self.dvgrids)
        if is_optimization:
            msg = ['$OPTIMIZATION\n']
            for (unused_id, dconstrs) in sorted(self.dconstrs.items()):
                for dconstr in dconstrs:
                    msg.append(dconstr.write_card(size, is_double))
            for (unused_id, desvar) in sorted(self.desvars.items()):
                msg.append(desvar.write_card(size, is_double))
            for (unused_id, ddval) in sorted(self.ddvals.items()):
                msg.append(ddval.write_card(size, is_double))
            for (unused_id, dlink) in sorted(self.dlinks.items()):
                msg.append(dlink.write_card(size, is_double))
            for (unused_id, dresp) in sorted(self.dresps.items()):
                msg.append(dresp.write_card(size, is_double))

            for (unused_id, dvcrel) in sorted(self.dvcrels.items()):
                msg.append(dvcrel.write_card(size, is_double))
            for (unused_id, dvmrel) in sorted(self.dvmrels.items()):
                msg.append(dvmrel.write_card(size, is_double))
            for (unused_id, dvprel) in sorted(self.dvprels.items()):
                msg.append(dvprel.write_card(size, is_double))
            for (unused_id, dvgrid) in sorted(self.dvgrids.items()):
                msg.append(dvgrid.write_card(size, is_double))
            for (unused_id, equation) in sorted(self.dequations.items()):
                msg.append(str(equation))
            if self.doptprm is not None:
                msg.append(self.doptprm.write_card(size, is_double))
            bdf_file.write(''.join(msg))

    def _write_params(self, bdf_file, size=8, is_double=False):
        """
        Writes the PARAM cards
        """
        if self.params:
            msg = ['$PARAMS\n']
            for (unused_key, param) in sorted(self.params.items()):
                msg.append(param.write_card(size))
            bdf_file.write(''.join(msg))

    def _write_rejects(self, bdf_file, size=8, is_double=False):
        """
        Writes the rejected (processed) cards and the rejected unprocessed
        cardlines
        """
        if size == 8:
            print_func = print_card_8
        else:
            print_func = print_card_16
        msg = []
        if self.reject_cards:
            msg.append('$REJECTS\n')
            for reject_card in self.reject_cards:
                try:
                    msg.append(print_func(reject_card))
                except RuntimeError:
                    for field in reject_card:
                        if field is not None and '=' in field:
                            raise SyntaxError('cannot reject equal signed '
                                              'cards\ncard=%s\n' % reject_card)
                    raise

        if self.rejects:
            msg.append('$REJECT_LINES\n')
        for reject_lines in self.reject_lines:
            if isinstance(reject_lines, (list, tuple)):
                for reject in reject_lines:
                    reject2 = reject.rstrip()
                    if reject2:
                        msg.append('%s\n' % reject2)
            elif isinstance(reject_lines, str):
                reject2 = reject_lines.rstrip()
                if reject2:
                    msg.append('%s\n' % reject2)
            else:
                raise TypeError(reject_lines)
        bdf_file.write(''.join(msg))

    def _write_rigid_elements(self, bdf_file, size=8, is_double=False):
        """Writes the rigid elements in a sorted order"""
        #self.rbe2.write_card(bdf_file, size, is_double)
        #self.rbe3.write_card(bdf_file, size, is_double)
        if self.rigid_elements:
            bdf_file.write('$RIGID ELEMENTS\n')
            if self.is_long_ids:
                for (eid, element) in sorted(self.rigid_elements.items()):
                    try:
                        bdf_file.write(element.write_card_16(is_double))
                    except:
                        print('failed printing element...'
                              'type=%s eid=%s' % (element.type, eid))
                        raise
            else:
                for (eid, element) in sorted(self.rigid_elements.items()):
                    try:
                        bdf_file.write(element.write_card(size, is_double))
                    except:
                        print('failed printing element...'
                              'type=%s eid=%s' % (element.type, eid))
                        raise
        if self.plotels:
            bdf_file.write('$PLOT ELEMENTS\n')
            for (eid, element) in sorted(self.plotels.items()):
                bdf_file.write(element.write_card(size, is_double))


    def _write_sets(self, bdf_file, size=8, is_double=False):
        """Writes the SETx cards sorted by ID"""
        is_sets = bool(self.sets or self.asets or self.omits or self.bsets or
                       self.csets or self.qsets or self.usets)
        if is_sets:
            msg = ['$SETS\n']
            for (unused_id, set_obj) in sorted(self.sets.items()):  # dict
                msg.append(set_obj.write_card(size, is_double))
            for set_obj in self.asets:  # list
                msg.append(set_obj.write_card(size, is_double))
            for set_obj in self.omits:  # list
                msg.append(set_obj.write_card(size, is_double))
            for set_obj in self.bsets:  # list
                msg.append(set_obj.write_card(size, is_double))
            for set_obj in self.csets:  # list
                msg.append(set_obj.write_card(size, is_double))
            for set_obj in self.qsets:  # list
                msg.append(set_obj.write_card(size, is_double))
            for unused_name, usets in sorted(self.usets.items()):  # dict
                for set_obj in usets:  # list
                    msg.append(set_obj.write_card(size, is_double))
            bdf_file.write(''.join(msg))

    def _write_superelements(self, bdf_file, size=8, is_double=False):
        """Writes the SETx cards sorted by ID"""
        is_sets = bool(self.se_sets or self.se_bsets or self.se_csets or
                       self.se_qsets or self.se_usets)
        if is_sets:
            msg = ['$SUPERELEMENTS\n']
            for set_obj in self.se_bsets:  # list
                msg.append(set_obj.write_card(size, is_double))
            for set_obj in self.se_csets:  # list
                msg.append(set_obj.write_card(size, is_double))
            for set_obj in self.se_qsets:  # list
                msg.append(set_obj.write_card(size, is_double))
            for unused_set_id, set_obj in sorted(self.se_sets.items()):  # dict
                msg.append(set_obj.write_card(size, is_double))
            for unused_name, usets in sorted(self.se_usets.items()):  # dict
                for set_obj in usets:  # list
                    msg.append(set_obj.write_card(size, is_double))
            for suport in self.se_suport:  # list
                msg.append(suport.write_card(size, is_double))
            bdf_file.write(''.join(msg))

    def _write_tables(self, bdf_file, size=8, is_double=False):
        """Writes the TABLEx cards sorted by ID"""
        if self.tables:
            msg = ['$TABLES\n']
            for (unused_id, table) in sorted(self.tables.items()):
                msg.append(table.write_card(size, is_double))
            bdf_file.write(''.join(msg))

        if self.random_tables:
            msg = ['$RANDOM TABLES\n']
            for (unused_id, table) in sorted(self.random_tables.items()):
                msg.append(table.write_card(size, is_double))
            bdf_file.write(''.join(msg))

    def _write_thermal(self, bdf_file, size=8, is_double=False):
        """Writes the thermal cards"""
        # PHBDY
        if self.phbdys or self.convection_properties or self.bcs:
            # self.thermalProperties or
            msg = ['$THERMAL\n']

            for unused_key, phbdy in sorted(self.phbdys.items()):
                msg.append(phbdy.write_card(size, is_double))

            #for unused_key, prop in sorted(self.thermalProperties.items()):
            #    msg.append(str(prop))
            for unused_key, prop in sorted(self.convection_properties.items()):
                msg.append(prop.write_card(size, is_double))

            # BCs
            for unused_key, bcs in sorted(self.bcs.items()):
                for bc in bcs:  # list
                    msg.append(bc.write_card(size, is_double))
            bdf_file.write(''.join(msg))

    def _write_thermal_materials(self, bdf_file, size=8, is_double=False):
        """Writes the thermal materials in a sorted order"""
        if self.thermal_materials:
            msg = ['$THERMAL MATERIALS\n']
            for unused_mid, material in sorted(self.thermal_materials.items()):
                msg.append(material.write_card(size, is_double))
            bdf_file.write(''.join(msg))
