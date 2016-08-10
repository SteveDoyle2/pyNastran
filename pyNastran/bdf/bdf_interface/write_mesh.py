# coding: utf-8
"""
This file defines:
  - WriteMesh
"""
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
import sys
import io
from codecs import open
from six import string_types, iteritems, itervalues, PY2, StringIO

#from pyNastran.utils import is_file_obj
from pyNastran.bdf.utils import print_filename
from pyNastran.utils.gui_io import save_file_dialog
from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.field_writer_16 import print_card_16
from pyNastran.bdf.field_writer_double import print_card_double
from pyNastran.bdf.bdf_interface.attributes import BDFAttributes


class WriteMesh(BDFAttributes):
    """
    Major methods:
      - model.write_bdf(...)
      - model.echo_bdf(...)
      - model.auto_reject_bdf(...)
    """
    def __init__(self):
        BDFAttributes.__init__(self)
        self._auto_reject = True
        self.cards_to_read = set([])

    def echo_bdf(self, infile_name):
        """
        This method removes all comment lines from the bdf
        A write method is stil required.

        .. todo:: maybe add the write method

        .. code-block:: python

          model = BDF()
          model.echo_bdf(bdf_filename)
        """
        self.deprecated('self.echo_bdf()', 'removed...', '0.8')
        self.cards_to_read = set([])
        return self.read_bdf(infile_name)

    def auto_reject_bdf(self, infile_name):
        """
        This method parses supported cards, but does not group them into
        nodes, elements, properties, etc.

        .. todo:: maybe add the write method
        """
        self.deprecated('self.auto_reject_bdf()', 'removed...', '0.8')
        self._auto_reject = True
        return self.read_bdf(infile_name)

    def _output_helper(self, out_filename, interspersed, size, is_double):
        """
        Performs type checking on the write_bdf inputs
        """
        if out_filename is None:
            wildcard_wx = "Nastran BDF (*.bdf; *.dat; *.nas; *.pch)|" \
                "*.bdf;*.dat;*.nas;*.pch|" \
                "All files (*.*)|*.*"
            wildcard_qt = "Nastran BDF (*.bdf *.dat *.nas *.pch);;All files (*)"
            title = 'Save BDF/DAT/PCH'
            out_filename = save_file_dialog(title, wildcard_wx, wildcard_qt)
            assert out_filename is not None, out_filename

        if PY2:
            #if not is_file_obj(out_filename):
            if not (hasattr(out_filename, 'read') and hasattr(out_filename, 'write')) or isinstance(out_filename, (file, StringIO)):
                return out_filename
            elif not isinstance(out_filename, string_types):
                raise TypeError('out_filename=%r must be a string; type=%s' % (out_filename, type(out_filename)))
        else:
            if not (hasattr(out_filename, 'read') and hasattr(out_filename, 'write')) or isinstance(out_filename, io.IOBase):
                return out_filename
            elif not isinstance(out_filename, string_types):
                raise TypeError('out_filename=%r must be a string; type=%s' % (out_filename, type(out_filename)))

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

    def write_caero_model(self, caero_bdf_filename='caero.bdf'):
        """write the CAERO cards as CQUAD4s that can be visualized"""
        bdf_file = open(caero_bdf_filename, 'w')
        bdf_file.write('CEND\n')
        bdf_file.write('BEGIN BULK\n')
        bdf_file.write('$ punch=True\n')
        i = 1

        mid = 1
        bdf_file.write('MAT1,%s,3.0E7,,0.3\n' % mid)
        for aesurf_id, aesurf in iteritems(self.aesurf):
            cid = aesurf.cid1
            bdf_file.write('PSHELL,%s,%s,0.1\n' % (aesurf_id, aesurf_id))
            #print(cid)
            #ax, ay, az = cid.i
            #bx, by, bz = cid.j
            #cx, cy, cz = cid.k
            #bdf_file.write('CORD2R,%s,,%s,%s,%s,%s,%s,%s\n' % (cid, ax, ay, az, bx, by, bz))
            #bdf_file.write(',%s,%s,%s\n' % (cx, cy, cz))
            #print(cid)
            bdf_file.write(str(cid))
            #aesurf.elements
        for eid, caero in sorted(iteritems(self.caeros)):
            assert eid != 1, 'CAERO eid=1 is reserved for non-flaps'
            scaero = str(caero).rstrip().split('\n')
            bdf_file.write('$ ' + '\n$ '.join(scaero) + '\n')
            points, elements = caero.panel_points_elements()
            npoints = points.shape[0]
            #nelements = elements.shape[0]
            for ipoint, point in enumerate(points):
                x, y, z = point
                bdf_file.write('GRID,%s,,%s,%s,%s\n' % (i + ipoint, x, y, z))

            pid = eid
            mid = eid
            #if 0:
                #bdf_file.write('PSHELL,%s,%s,0.1\n' % (pid, mid))
                #bdf_file.write('MAT1,%s,3.0E7,,0.3\n' % mid)
            #else:
            bdf_file.write('PSHELL,%s,%s,0.1\n' % (1, 1))
            bdf_file.write('MAT1,%s,3.0E7,,0.3\n' % 1)

            j = 0
            for elem in elements + i:
                p1, p2, p3, p4 = elem
                eid2 = j + eid
                pidi = None
                for aesurf_id, aesurf in iteritems(self.aesurf):
                    aelist_id = aesurf.AELIST_id1()
                    aelist = self.aelists[aelist_id]
                    if eid2 in aelist.elements:
                        pidi = aesurf_id
                        break
                if pidi is None:
                    #pidi = pid
                    pidi = 1
                bdf_file.write('CQUAD4,%s,%s,%s,%s,%s,%s\n' % (j + eid, pidi, p1, p2, p3, p4))
                j += 1
            i += npoints
            #break
            #j += nelements
        bdf_file.write('ENDDATA\n')

    def get_encoding(self, encoding=None):
        if encoding is not None:
            pass
        else:
            encoding = self._encoding
            if encoding is None:
                encoding = sys.getdefaultencoding()
        return encoding

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
        interspersed : bool; default=True
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
        #self.write_caero_model()
        out_filename = self._output_helper(out_filename,
                                           interspersed, size, is_double)
        encoding = self.get_encoding(encoding)
        #assert encoding.lower() in ['ascii', 'latin1', 'utf8'], encoding

        if hasattr(out_filename, 'read') and hasattr(out_filename, 'write'):
            outfile = out_filename
        else:
            if PY2:
                outfile = open(out_filename, 'wb', encoding=encoding)
            else:
                outfile = open(out_filename, 'w', encoding=encoding)
        self._write_header(outfile, encoding)
        self._write_params(outfile, size, is_double)
        self._write_nodes(outfile, size, is_double)

        if interspersed:
            self._write_elements_properties(outfile, size, is_double)
        else:
            self._write_elements(outfile, size, is_double)
            self._write_properties(outfile, size, is_double)
        self._write_materials(outfile, size, is_double)

        self._write_masses(outfile, size, is_double)
        self._write_common(outfile, size, is_double)
        if (enddata is None and 'ENDDATA' in self.card_count) or enddata:
            outfile.write('ENDDATA\n')
        if close:
            outfile.close()

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
        interspersed = False
        #self.write_caero_model()
        out_filename = self._output_helper(out_filename,
                                           interspersed, size, is_double)
        if encoding is not None:
            pass
        else:
            encoding = self._encoding
            if encoding is None:
                encoding = sys.getdefaultencoding()
        #assert encoding.lower() in ['ascii', 'latin1', 'utf8'], encoding

        if hasattr(out_filename, 'read') and hasattr(out_filename, 'write'):
            outfile = out_filename
        else:
            if PY2:
                outfile = open(out_filename, 'wb', encoding=encoding)
            else:
                outfile = open(out_filename, 'w', encoding=encoding)
        self._write_header(outfile, encoding)
        self._write_params(outfile, size, is_double)
        self._write_nodes_symmetric(outfile, size, is_double, plane=plane)

        if interspersed:
            raise RuntimeError(interspersed)
            #self._write_elements_properties(outfile, size, is_double)
        else:
            self._write_elements_symmetric(outfile, size, is_double)
            self._write_properties(outfile, size, is_double)
        self._write_materials(outfile, size, is_double)

        self._write_masses(outfile, size, is_double)
        self._write_common(outfile, size, is_double)
        if (enddata is None and 'ENDDATA' in self.card_count) or enddata:
            outfile.write('ENDDATA\n')
        if close:
            outfile.close()

    def _write_header(self, outfile, encoding):
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
            outfile.write('$pyNastran: version=%s\n' % self.nastran_format)
            outfile.write('$pyNastran: punch=%s\n' % self.punch)
            outfile.write('$pyNastran: encoding=%s\n' % encoding)
            outfile.write('$pyNastran: nnodes=%s\n' % len(self.nodes))
            outfile.write('$pyNastran: nelements=%s\n' % len(self.elements))

        if not self.punch:
            self._write_executive_control_deck(outfile)
            self._write_case_control_deck(outfile)

    def _write_executive_control_deck(self, outfile):
        """
        Writes the executive control deck.
        """
        if self.executive_control_lines:
            msg = '$EXECUTIVE CONTROL DECK\n'
            if self.sol == 600:
                new_sol = 'SOL 600,%s' % self.solMethod
            else:
                new_sol = 'SOL %s' % self.sol

            if self.iSolLine is not None:
                self.executive_control_lines[self.iSolLine] = new_sol

            for line in self.executive_control_lines:
                msg += line + '\n'
            outfile.write(msg)

    def _write_case_control_deck(self, outfile):
        """
        Writes the Case Control Deck.
        """
        if self.case_control_deck:
            msg = '$CASE CONTROL DECK\n'
            msg += str(self.case_control_deck)
            assert 'BEGIN BULK' in msg, msg
            outfile.write(''.join(msg))

    def _write_elements(self, outfile, size=8, is_double=False):
        """
        Writes the elements in a sorted order
        """
        if self.elements:
            outfile.write('$ELEMENTS\n')
            if self.is_long_ids:
                for (eid, element) in sorted(iteritems(self.elements)):
                    outfile.write(element.write_card_16(is_double))
            else:
                for (eid, element) in sorted(iteritems(self.elements)):
                    try:
                        outfile.write(element.write_card(size, is_double))
                    except:
                        print('failed printing element...'
                              'type=%s eid=%s' % (element.type, eid))
                        raise

    def _write_elements_symmetric(self, outfile, size=8, is_double=False):
        """
        Writes the elements in a sorted order
        """
        nid_offset = max(self.nodes.keys())
        eid_offset = max(self.elements.keys())
        if self.elements:
            outfile.write('$ELEMENTS\n')
            if self.is_long_ids:
                for (eid, element) in sorted(iteritems(self.elements)):
                    nodes = element.node_ids
                    outfile.write(element.write_card_16(is_double))
                    element.eid += eid_offset
                    nodes = [node_id + nid_offset for node_id in nodes]
                    element.nodes = nodes
                    outfile.write(element.write_card_16(is_double))
            else:
                for (eid, element) in sorted(iteritems(self.elements)):
                    nodes = element.node_ids
                    outfile.write(element.write_card(size, is_double))
                    nodes = [node_id + nid_offset for node_id in nodes]
                    element.eid += eid_offset
                    element.nodes = nodes
                    outfile.write(element.write_card(size, is_double))

    def _write_elements_properties(self, outfile, size=8, is_double=False):
        """
        Writes the elements and properties in and interspersed order
        """
        missing_properties = []
        if self.properties:
            outfile.write('$ELEMENTS_WITH_PROPERTIES\n')

        eids_written = []
        pids = sorted(self.properties.keys())
        pid_eids = self.get_element_ids_dict_with_pids(pids)

        msg = []
        #failed_element_types = set([])
        for (pid, eids) in sorted(iteritems(pid_eids)):
            prop = self.properties[pid]
            if eids:
                msg.append(prop.write_card(size, is_double))
                eids.sort()
                for eid in eids:
                    element = self.Element(eid)
                    try:
                        msg.append(element.write_card(size, is_double))
                    except:
                        print('failed printing element...' 'type=%r eid=%s'
                              % (element.type, eid))
                        raise
                eids_written += eids
            else:
                missing_properties.append(prop.write_card(size, is_double))
        outfile.write(''.join(msg))

        eids_missing = set(self.elements.keys()).difference(set(eids_written))
        if eids_missing:
            msg = ['$ELEMENTS_WITH_NO_PROPERTIES '
                   '(PID=0 and unanalyzed properties)\n']
            for eid in sorted(eids_missing):
                element = self.Element(eid, msg='')
                try:
                    msg.append(element.write_card(size, is_double))
                except:
                    print('failed printing element...'
                          'type=%s eid=%s' % (element.type, eid))
                    raise
            outfile.write(''.join(msg))

        if missing_properties or self.pdampt or self.pbusht or self.pelast:
            msg = ['$UNASSOCIATED_PROPERTIES\n']
            for card in sorted(itervalues(self.pbusht)):
                msg.append(card.write_card(size, is_double))
            for card in sorted(itervalues(self.pdampt)):
                msg.append(card.write_card(size, is_double))
            for card in sorted(itervalues(self.pelast)):
                msg.append(card.write_card(size, is_double))
            for card in missing_properties:
                # this is a string...
                #print("missing_property = ", card
                msg.append(card)
            outfile.write(''.join(msg))

    def _write_aero(self, outfile, size=8, is_double=False):
        """Writes the aero cards"""
        if self.caeros or self.paeros or self.monitor_points or self.splines:
            msg = ['$AERO\n']
            for (unused_id, caero) in sorted(iteritems(self.caeros)):
                msg.append(caero.write_card(size, is_double))
            for (unused_id, paero) in sorted(iteritems(self.paeros)):
                msg.append(paero.write_card(size, is_double))
            for (unused_id, spline) in sorted(iteritems(self.splines)):
                msg.append(spline.write_card(size, is_double))
            for monitor_point in self.monitor_points:
                msg.append(monitor_point.write_card(size, is_double))
            outfile.write(''.join(msg))

    def _write_aero_control(self, outfile, size=8, is_double=False):
        """Writes the aero control surface cards"""
        if(self.aecomps or self.aefacts or self.aeparams or self.aelinks or
           self.aelists or self.aestats or self.aesurf or self.aesurfs):
            msg = ['$AERO CONTROL SURFACES\n']
            for (unused_id, aelinks) in sorted(iteritems(self.aelinks)):
                for aelink in aelinks:
                    msg.append(aelink.write_card(size, is_double))

            for (unused_id, aecomp) in sorted(iteritems(self.aecomps)):
                msg.append(aecomp.write_card(size, is_double))
            for (unused_id, aeparam) in sorted(iteritems(self.aeparams)):
                msg.append(aeparam.write_card(size, is_double))
            for (unused_id, aestat) in sorted(iteritems(self.aestats)):
                msg.append(aestat.write_card(size, is_double))

            for (unused_id, aelist) in sorted(iteritems(self.aelists)):
                msg.append(aelist.write_card(size, is_double))
            for (unused_id, aesurf) in sorted(iteritems(self.aesurf)):
                msg.append(aesurf.write_card(size, is_double))
            for (unused_id, aesurfs) in sorted(iteritems(self.aesurfs)):
                msg.append(aesurfs.write_card(size, is_double))
            for (unused_id, aefact) in sorted(iteritems(self.aefacts)):
                msg.append(aefact.write_card(size, is_double))
            outfile.write(''.join(msg))

    def _write_static_aero(self, outfile, size=8, is_double=False):
        """Writes the static aero cards"""
        if self.aeros or self.trims or self.divergs:
            msg = ['$STATIC AERO\n']
            # static aero
            if self.aeros:
                msg.append(self.aeros.write_card(size, is_double))
            for (unused_id, trim) in sorted(iteritems(self.trims)):
                msg.append(trim.write_card(size, is_double))
            for (unused_id, diverg) in sorted(iteritems(self.divergs)):
                msg.append(diverg.write_card(size, is_double))
            outfile.write(''.join(msg))

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


    def _write_flutter(self, outfile, size=8, is_double=False, write_aero_in_flutter=True):
        """Writes the flutter cards"""
        if (write_aero_in_flutter and self.aero) or self.flfacts or self.flutters or self.mkaeros:
            msg = ['$FLUTTER\n']
            if write_aero_in_flutter:
                msg.append(self.aero.write_card(size, is_double))
            for (unused_id, flutter) in sorted(iteritems(self.flutters)):
                msg.append(flutter.write_card(size, is_double))
            for (unused_id, flfact) in sorted(iteritems(self.flfacts)):
                msg.append(flfact.write_card(size, is_double))
            for mkaero in self.mkaeros:
                msg.append(mkaero.write_card(size, is_double))
            outfile.write(''.join(msg))

    def _write_gust(self, outfile, size=8, is_double=False, write_aero_in_gust=True):
        """Writes the gust cards"""
        if (write_aero_in_gust and self.aero) or self.gusts:
            msg = ['$GUST\n']
            if write_aero_in_gust:
                for (unused_id, aero) in sorted(iteritems(self.aero)):
                    msg.append(aero.write_card(size, is_double))
            for (unused_id, gust) in sorted(iteritems(self.gusts)):
                msg.append(gust.write_card(size, is_double))
            outfile.write(''.join(msg))

    def _write_common(self, outfile, size=8, is_double=False):
        """
        Write the common outputs so none get missed...

        Parameters
        ----------
        outfile : file
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
        self._write_rigid_elements(outfile, size, is_double)
        self._write_dmigs(outfile, size, is_double)
        self._write_loads(outfile, size, is_double)
        self._write_dynamic(outfile, size, is_double)
        self._write_aero(outfile, size, is_double)
        self._write_aero_control(outfile, size, is_double)
        self._write_static_aero(outfile, size, is_double)

        write_aero_in_flutter, write_aero_in_gust = self._find_aero_location()
        self._write_flutter(outfile, size, is_double, write_aero_in_flutter)
        self._write_gust(outfile, size, is_double, write_aero_in_gust)

        self._write_thermal(outfile, size, is_double)
        self._write_thermal_materials(outfile, size, is_double)

        self._write_constraints(outfile, size, is_double)
        self._write_optimization(outfile, size, is_double)
        self._write_tables(outfile, size, is_double)
        self._write_sets(outfile, size, is_double)
        self._write_superelements(outfile, size, is_double)
        self._write_contact(outfile, size, is_double)
        self._write_rejects(outfile, size, is_double)
        self._write_coords(outfile, size, is_double)

    def _write_constraints(self, outfile, size=8, is_double=False):
        """Writes the constraint cards sorted by ID"""
        if self.suport or self.suport1:
            msg = ['$CONSTRAINTS\n']
            for suport in self.suport:
                msg.append(suport.write_card(size, is_double))
            for suport_id, suport in sorted(iteritems(self.suport1)):
                msg.append(suport.write_card(size, is_double))
            outfile.write(''.join(msg))

        if self.spcs or self.spcadds:
            #msg = ['$SPCs\n']
            #str_spc = str(self.spcObject) # old
            #if str_spc:
                #msg.append(str_spc)
            #else:
            msg = ['$SPCs\n']
            for (unused_id, spcadd) in sorted(iteritems(self.spcadds)):
                msg.append(str(spcadd))
            for (unused_id, spcs) in sorted(iteritems(self.spcs)):
                for spc in spcs:
                    msg.append(str(spc))
            outfile.write(''.join(msg))

        if self.mpcs or self.mpcadds:
            msg = ['$MPCs\n']
            for (unused_id, mpcadd) in sorted(iteritems(self.mpcadds)):
                msg.append(str(mpcadd))
            for (unused_id, mpcs) in sorted(iteritems(self.mpcs)):
                for mpc in mpcs:
                    msg.append(mpc.write_card(size, is_double))
            outfile.write(''.join(msg))


    def _write_contact(self, outfile, size=8, is_double=False):
        """Writes the contact cards sorted by ID"""
        if(self.bcrparas or self.bctadds or self.bctparas or self.bctsets
           or self.bsurf or self.bsurfs):
            msg = ['$CONTACT\n']
            for (unused_id, bcrpara) in sorted(iteritems(self.bcrparas)):
                msg.append(bcrpara.write_card(size, is_double))
            for (unused_id, bctadds) in sorted(iteritems(self.bctadds)):
                msg.append(bctadds.write_card(size, is_double))
            for (unused_id, bctpara) in sorted(iteritems(self.bctparas)):
                msg.append(bctpara.write_card(size, is_double))

            for (unused_id, bctset) in sorted(iteritems(self.bctsets)):
                msg.append(bctset.write_card(size, is_double))
            for (unused_id, bsurfi) in sorted(iteritems(self.bsurf)):
                msg.append(bsurfi.write_card(size, is_double))
            for (unused_id, bsurfsi) in sorted(iteritems(self.bsurfs)):
                msg.append(bsurfsi.write_card(size, is_double))
            outfile.write(''.join(msg))

    def _write_coords(self, outfile, size=8, is_double=False):
        """Writes the coordinate cards in a sorted order"""
        msg = []
        if len(self.coords) > 1:
            msg.append('$COORDS\n')
        for (unused_id, coord) in sorted(iteritems(self.coords)):
            if unused_id != 0:
                msg.append(coord.write_card(size, is_double))
        outfile.write(''.join(msg))

    def _write_dmigs(self, outfile, size=8, is_double=False):
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
        for (unused_name, dmig) in sorted(iteritems(self.dmigs)):
            msg.append(dmig.write_card(size, is_double))
        for (unused_name, dmi) in sorted(iteritems(self.dmis)):
            msg.append(dmi.write_card(size, is_double))
        for (unused_name, dmij) in sorted(iteritems(self.dmijs)):
            msg.append(dmij.write_card(size, is_double))
        for (unused_name, dmiji) in sorted(iteritems(self.dmijis)):
            msg.append(dmiji.write_card(size, is_double))
        for (unused_name, dmik) in sorted(iteritems(self.dmiks)):
            msg.append(dmik.write_card(size, is_double))
        outfile.write(''.join(msg))

    def _write_dynamic(self, outfile, size=8, is_double=False):
        """Writes the dynamic cards sorted by ID"""
        if(self.dareas or self.dphases or self.nlparms or self.frequencies or self.methods or
           self.cMethods or self.tsteps or self.tstepnls or self.transfer_functions or
           self.delays):
            msg = ['$DYNAMIC\n']
            for (unused_id, method) in sorted(iteritems(self.methods)):
                msg.append(method.write_card(size, is_double))
            for (unused_id, cmethod) in sorted(iteritems(self.cMethods)):
                msg.append(cmethod.write_card(size, is_double))
            for (unused_id, darea) in sorted(iteritems(self.dareas)):
                msg.append(darea.write_card(size, is_double))
            for (unused_id, dphase) in sorted(iteritems(self.dphases)):
                msg.append(dphase.write_card(size, is_double))
            for (unused_id, nlparm) in sorted(iteritems(self.nlparms)):
                msg.append(nlparm.write_card(size, is_double))
            for (unused_id, nlpci) in sorted(iteritems(self.nlpcis)):
                msg.append(nlpci.write_card(size, is_double))
            for (unused_id, tstep) in sorted(iteritems(self.tsteps)):
                msg.append(tstep.write_card(size, is_double))
            for (unused_id, tstepnl) in sorted(iteritems(self.tstepnls)):
                msg.append(tstepnl.write_card(size, is_double))
            for (unused_id, freq) in sorted(iteritems(self.frequencies)):
                msg.append(freq.write_card(size, is_double))
            for (unused_id, delay) in sorted(iteritems(self.delays)):
                msg.append(delay.write_card(size, is_double))

            for (unused_id, tfs) in sorted(iteritems(self.transfer_functions)):
                for tf in tfs:
                    msg.append(tf.write_card(size, is_double))
            outfile.write(''.join(msg))

    def _write_loads(self, outfile, size=8, is_double=False):
        """Writes the load cards sorted by ID"""
        if self.loads or self.tempds:
            msg = ['$LOADS\n']
            for (key, loadcase) in sorted(iteritems(self.loads)):
                for load in loadcase:
                    try:
                        msg.append(load.write_card(size, is_double))
                    except:
                        print('failed printing load...type=%s key=%r'
                              % (load.type, key))
                        raise
            for key, tempd in sorted(iteritems(self.tempds)):
                msg.append(tempd.write_card(size, is_double))
            outfile.write(''.join(msg))

        if self.dloads or self.dload_entries:
            msg = ['$DLOADS\n']
            for (key, loadcase) in sorted(iteritems(self.dloads)):
                for load in loadcase:
                    try:
                        msg.append(load.write_card(size, is_double))
                    except:
                        print('failed printing load...type=%s key=%r'
                              % (load.type, key))
                        raise

            for (key, loadcase) in sorted(iteritems(self.dload_entries)):
                for load in loadcase:
                    try:
                        msg.append(load.write_card(size, is_double))
                    except:
                        print('failed printing load...type=%s key=%r'
                              % (load.type, key))
                        raise
            outfile.write(''.join(msg))


    def _write_masses(self, outfile, size=8, is_double=False):
        """Writes the mass cards sorted by ID"""
        if self.properties_mass:
            outfile.write('$PROPERTIES_MASS\n')
            for (pid, mass) in sorted(iteritems(self.properties_mass)):
                try:
                    outfile.write(mass.write_card(size, is_double))
                except:
                    print('failed printing mass property...'
                          'type=%s eid=%s' % (mass.type, pid))
                    raise

        if self.masses:
            outfile.write('$MASSES\n')
            for (eid, mass) in sorted(iteritems(self.masses)):
                try:
                    outfile.write(mass.write_card(size, is_double))
                except:
                    print('failed printing masses...'
                          'type=%s eid=%s' % (mass.type, eid))
                    raise

    def _write_materials(self, outfile, size=8, is_double=False):
        """Writes the materials in a sorted order"""
        if(self.materials or self.hyperelastic_materials or self.creep_materials or
           self.MATS1 or self.MATS3 or self.MATS8 or self.MATT1 or
           self.MATT2 or self.MATT3 or self.MATT4 or self.MATT5 or
           self.MATT8 or self.MATT9):
            msg = ['$MATERIALS\n']
            for (unused_mid, material) in sorted(iteritems(self.materials)):
                msg.append(material.write_card(size, is_double))
            for (unused_mid, material) in sorted(iteritems(self.hyperelastic_materials)):
                msg.append(material.write_card(size, is_double))
            for (unused_mid, material) in sorted(iteritems(self.creep_materials)):
                msg.append(material.write_card(size, is_double))

            for (unused_mid, material) in sorted(iteritems(self.MATS1)):
                msg.append(material.write_card(size, is_double))
            for (unused_mid, material) in sorted(iteritems(self.MATS3)):
                msg.append(material.write_card(size, is_double))
            for (unused_mid, material) in sorted(iteritems(self.MATS8)):
                msg.append(material.write_card(size, is_double))

            for (unused_mid, material) in sorted(iteritems(self.MATT1)):
                msg.append(material.write_card(size, is_double))
            for (unused_mid, material) in sorted(iteritems(self.MATT2)):
                msg.append(material.write_card(size, is_double))
            for (unused_mid, material) in sorted(iteritems(self.MATT3)):
                msg.append(material.write_card(size, is_double))
            for (unused_mid, material) in sorted(iteritems(self.MATT4)):
                msg.append(material.write_card(size, is_double))
            for (unused_mid, material) in sorted(iteritems(self.MATT5)):
                msg.append(material.write_card(size, is_double))
            for (unused_mid, material) in sorted(iteritems(self.MATT8)):
                msg.append(material.write_card(size, is_double))
            for (unused_mid, material) in sorted(iteritems(self.MATT9)):
                msg.append(material.write_card(size, is_double))
            outfile.write(''.join(msg))

    def _write_nodes(self, outfile, size=8, is_double=False):
        """
        Writes the NODE-type cards
        """
        if self.spoints:
            msg = []
            msg.append('$SPOINTS\n')
            msg.append(self.spoints.write_card(size, is_double))
            outfile.write(''.join(msg))
        if self.epoints:
            msg = []
            msg.append('$EPOINTS\n')
            msg.append(self.epoints.write_card(size, is_double))
            outfile.write(''.join(msg))

        if self.nodes:
            msg = []
            msg.append('$NODES\n')
            if self.gridSet:
                msg.append(self.gridSet.print_card(size))

            if self.is_long_ids:
                for (unused_nid, node) in sorted(iteritems(self.nodes)):
                    msg.append(node.write_card_16(is_double))
            else:
                for (unused_nid, node) in sorted(iteritems(self.nodes)):
                    msg.append(node.write_card(size, is_double))
            outfile.write(''.join(msg))
        #if 0:  # not finished
            #self._write_nodes_associated(outfile, size, is_double)

    #def _write_nodes_associated(self, outfile, size=8, is_double=False):
        #"""
        #Writes the NODE-type in associated and unassociated groups.

        #.. warning:: Sometimes crashes, probably on invalid BDFs.
        #"""
        #msg = []
        #associated_nodes = set([])
        #for (eid, element) in iteritems(self.elements):
            #associated_nodes = associated_nodes.union(set(element.node_ids))

        #all_nodes = set(self.nodes.keys())
        #unassociated_nodes = list(all_nodes.difference(associated_nodes))
        ##missing_nodes = all_nodes.difference(

        ## TODO: this really shouldn't be a list...???
        #associated_nodes = list(associated_nodes)

        #if associated_nodes:
            #msg += ['$ASSOCIATED NODES\n']
            #if self.gridSet:
                #msg.append(self.gridSet.write_card(size, is_double))
            ## TODO: this really shouldn't be a dictionary...???
            #for key, node in sorted(iteritems(associated_nodes)):
                #msg.append(node.write_card(size, is_double))

        #if unassociated_nodes:
            #msg.append('$UNASSOCIATED NODES\n')
            #if self.gridSet and not associated_nodes:
                #msg.append(self.gridSet.write_card(size, is_double))
            #for key, node in sorted(iteritems(unassociated_nodes)):
                #if key in self.nodes:
                    #msg.append(node.write_card(size, is_double))
                #else:
                    #msg.append('$ Missing NodeID=%s' % key)
        #outfile.write(''.join(msg))

    def _write_nodes_symmetric(self, outfile, size=8, is_double=False, plane='xz'):
        """
        Writes the NODE-type cards

        .. warning:: doesn't consider coordinate systems;
                      it could, but you'd need 20 new coordinate systems
        .. warning:: doesn't mirror SPOINTs, EPOINTs
        """
        if self.spoints:
            msg = []
            msg.append('$SPOINTS\n')
            msg.append(self.spoints.write_card(size, is_double))
            outfile.write(''.join(msg))
        if self.epoints:
            msg = []
            msg.append('$EPOINTS\n')
            msg.append(self.epoints.write_card(size, is_double))
            outfile.write(''.join(msg))

        plane = plane.strip().lower()
        if plane == 'xz':
            iy = 4
        elif plane == 'xy':
            iy = 5
        elif plane == 'yz':
            iy = 3
        else:
            raise NotImplementedError(plane)
        if self.nodes:
            msg = []
            msg.append('$NODES\n')
            if self.gridSet:
                msg.append(self.gridSet.print_card(size))

            nid_offset = max(self.nodes.keys())
            if self.is_long_ids:
                print_card_long = print_card_double if is_double else print_card_16
                for (unused_nid, node) in sorted(iteritems(self.nodes)):
                    repr_fields = node.repr_fields()
                    msg.append(print_card_long(repr_fields))
                    repr_fields = node.repr_fields()
                    # grid, nid, cp, x, y, z
                    repr_fields[1] += nid_offset
                    repr_fields[iy] *= -1.0
                    msg.append(print_card_long(repr_fields))
            else:
                for (unused_nid, node) in sorted(iteritems(self.nodes)):
                    repr_fields = node.repr_fields()
                    msg.append(print_card_8(repr_fields))
                    repr_fields = node.repr_fields()
                    # grid, nid, cp, x, y, z
                    repr_fields[1] += nid_offset
                    repr_fields[iy] *= -1.0
                    msg.append(print_card_8(repr_fields))
            outfile.write(''.join(msg))
        #if 0:  # not finished
            #self._write_nodes_associated(outfile, size, is_double)

    def _write_optimization(self, outfile, size=8, is_double=False):
        """Writes the optimization cards sorted by ID"""
        if(self.dconadds or self.dconstrs or self.desvars or self.ddvals or
           self.dresps or self.ddvals or
           self.dvprels or self.dvmrels or self.dvcrels or self.doptprm or
           self.dlinks or self.dequations or self.dtable is not None):
            msg = ['$OPTIMIZATION\n']
            for (unused_id, dconadd) in sorted(iteritems(self.dconadds)):
                msg.append(dconadd.write_card(size, is_double))
            for (unused_id, dconstrs) in sorted(iteritems(self.dconstrs)):
                for dconstr in dconstrs:
                    msg.append(dconstr.write_card(size, is_double))
            for (unused_id, desvar) in sorted(iteritems(self.desvars)):
                msg.append(desvar.write_card(size, is_double))
            for (unused_id, ddval) in sorted(iteritems(self.ddvals)):
                msg.append(ddval.write_card(size, is_double))
            for (unused_id, dlink) in sorted(iteritems(self.dlinks)):
                msg.append(dlink.write_card(size, is_double))
            for (unused_id, dresp) in sorted(iteritems(self.dresps)):
                msg.append(dresp.write_card(size, is_double))

            for (unused_id, dvcrel) in sorted(iteritems(self.dvcrels)):
                msg.append(dvcrel.write_card(size, is_double))
            for (unused_id, dvmrel) in sorted(iteritems(self.dvmrels)):
                msg.append(dvmrel.write_card(size, is_double))
            for (unused_id, dvprel) in sorted(iteritems(self.dvprels)):
                msg.append(dvprel.write_card(size, is_double))
            for (unused_id, equation) in sorted(iteritems(self.dequations)):
                msg.append(str(equation))

            if self.dtable is not None:
                msg.append(self.dtable.write_card(size, is_double))
            if self.doptprm is not None:
                msg.append(self.doptprm.write_card(size, is_double))
            outfile.write(''.join(msg))

    def _write_params(self, outfile, size=8, is_double=False):
        """
        Writes the PARAM cards
        """
        if self.params:
            msg = ['$PARAMS\n']
            if self.is_long_ids:
                for (unused_key, param) in sorted(iteritems(self.params)):
                    msg.append(param.write_card(16, is_double))
            else:
                for (unused_key, param) in sorted(iteritems(self.params)):
                    msg.append(param.write_card(size, is_double))
            outfile.write(''.join(msg))

    def _write_properties(self, outfile, size=8, is_double=False):
        """Writes the properties in a sorted order"""
        if self.properties:
            msg = ['$PROPERTIES\n']
            prop_groups = (self.properties, self.pelast, self.pdampt, self.pbusht)
            if self.is_long_ids:
                for prop_group in prop_groups:
                    for unused_pid, prop in sorted(iteritems(prop_group)):
                        msg.append(prop.write_card_16(is_double))
                #except:
                    #print('failed printing property type=%s' % prop.type)
                    #raise
            else:
                for prop_group in prop_groups:
                    for unused_pid, prop in sorted(iteritems(prop_group)):
                        msg.append(prop.write_card(size, is_double))
            outfile.write(''.join(msg))

    def _write_rejects(self, outfile, size=8, is_double=False):
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
            elif isinstance(reject_lines, string_types):
                reject2 = reject_lines.rstrip()
                if reject2:
                    msg.append('%s\n' % reject2)
            else:
                raise TypeError(reject_lines)
        outfile.write(''.join(msg))

    def _write_rigid_elements(self, outfile, size=8, is_double=False):
        """Writes the rigid elements in a sorted order"""
        if self.rigid_elements:
            outfile.write('$RIGID ELEMENTS\n')
            if self.is_long_ids:
                for (eid, element) in sorted(iteritems(self.rigid_elements)):
                    try:
                        outfile.write(element.write_card_16(is_double))
                    except:
                        print('failed printing element...'
                              'type=%s eid=%s' % (element.type, eid))
                        raise
            else:
                for (eid, element) in sorted(iteritems(self.rigid_elements)):
                    try:
                        outfile.write(element.write_card(size, is_double))
                    except:
                        print('failed printing element...'
                              'type=%s eid=%s' % (element.type, eid))
                        raise
        if self.plotels:
            outfile.write('$PLOT ELEMENTS\n')
            for (eid, element) in sorted(iteritems(self.plotels)):
                outfile.write(element.write_card(size, is_double))

    def _write_sets(self, outfile, size=8, is_double=False):
        """Writes the SETx cards sorted by ID"""
        if(self.sets or self.asets or self.bsets or self.csets or self.qsets
           or self.usets):
            msg = ['$SETS\n']
            for (unused_id, set_obj) in sorted(iteritems(self.sets)):  # dict
                msg.append(set_obj.write_card(size, is_double))
            for set_obj in self.asets:  # list
                msg.append(set_obj.write_card(size, is_double))
            for set_obj in self.bsets:  # list
                msg.append(set_obj.write_card(size, is_double))
            for set_obj in self.csets:  # list
                msg.append(set_obj.write_card(size, is_double))
            for set_obj in self.qsets:  # list
                msg.append(set_obj.write_card(size, is_double))
            for name, usets in sorted(iteritems(self.usets)):  # dict
                for set_obj in usets:  # list
                    msg.append(set_obj.write_card(size, is_double))
            outfile.write(''.join(msg))

    def _write_superelements(self, outfile, size=8, is_double=False):
        """Writes the SETx cards sorted by ID"""
        if(self.se_sets or self.se_bsets or self.se_csets or self.se_qsets
           or self.se_usets):
            msg = ['$SUPERELEMENTS\n']
            for set_obj in self.se_bsets:  # list
                msg.append(set_obj.write_card(size, is_double))
            for set_obj in self.se_csets:  # list
                msg.append(set_obj.write_card(size, is_double))
            for set_obj in self.se_qsets:  # list
                msg.append(set_obj.write_card(size, is_double))
            for (set_id, set_obj) in sorted(iteritems(self.se_sets)):  # dict
                msg.append(set_obj.write_card(size, is_double))
            for name, usets in sorted(iteritems(self.se_usets)):  # dict
                for set_obj in usets:  # list
                    msg.append(set_obj.write_card(size, is_double))
            for suport in self.se_suport:  # list
                msg.append(suport.write_card(size, is_double))
            outfile.write(''.join(msg))

    def _write_tables(self, outfile, size=8, is_double=False):
        """Writes the TABLEx cards sorted by ID"""
        if self.tables or self.tables_sdamping:
            msg = ['$TABLES\n']
            for (unused_id, table) in sorted(iteritems(self.tables)):
                msg.append(table.write_card(size, is_double))
            for (unused_id, table) in sorted(iteritems(self.tables_sdamping)):
                msg.append(table.write_card(size, is_double))
            outfile.write(''.join(msg))

        if self.random_tables:
            msg = ['$RANDOM TABLES\n']
            for (unused_id, table) in sorted(iteritems(self.random_tables)):
                msg.append(table.write_card(size, is_double))
            outfile.write(''.join(msg))

    def _write_thermal(self, outfile, size=8, is_double=False):
        """Writes the thermal cards"""
        # PHBDY
        if self.phbdys or self.convection_properties or self.bcs:
            # self.thermalProperties or
            msg = ['$THERMAL\n']

            for (unused_key, phbdy) in sorted(iteritems(self.phbdys)):
                msg.append(phbdy.write_card(size, is_double))

            #for unused_key, prop in sorted(iteritems(self.thermalProperties)):
            #    msg.append(str(prop))
            for (unused_key, prop) in sorted(iteritems(self.convection_properties)):
                msg.append(prop.write_card(size, is_double))

            # BCs
            for (unused_key, bcs) in sorted(iteritems(self.bcs)):
                for boundary_condition in bcs:  # list
                    msg.append(boundary_condition.write_card(size, is_double))
            outfile.write(''.join(msg))

    def _write_thermal_materials(self, outfile, size=8, is_double=False):
        """Writes the thermal materials in a sorted order"""
        if self.thermal_materials:
            msg = ['$THERMAL MATERIALS\n']
            for (unused_mid, material) in sorted(iteritems(self.thermal_materials)):
                msg.append(material.write_card(size, is_double))
            outfile.write(''.join(msg))
