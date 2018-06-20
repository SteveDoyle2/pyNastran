# coding: utf-8
"""
This file defines:
  - WriteMesh
"""
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
import sys
import io
from typing import List, Dict, Union, Optional, Tuple, Any, cast
from codecs import open
from six import string_types, iteritems, itervalues, PY2, StringIO

from pyNastran.bdf.utils import print_filename
from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.field_writer_16 import print_card_16
from pyNastran.bdf.bdf_interface.attributes import BDFAttributes
from pyNastran.bdf.cards.nodes import write_xpoints


class WriteMesh(BDFAttributes):
    """
    Defines methods for writing cards

    Major methods:
      - model.write_bdf(...)
      - model.echo_bdf(...)
      - model.auto_reject_bdf(...)
    """
    def __init__(self):
        """creates methods for writing cards"""
        BDFAttributes.__init__(self)
        self._auto_reject = True
        self.cards_to_read = set([])

    def get_encoding(self, encoding=None):
        # type: (Optional[str]) -> str
        """gets the file encoding"""
        if encoding is not None:
            pass
        else:
            encoding = self._encoding
            if encoding is None:
                encoding = sys.getdefaultencoding()
        encoding = cast(str, encoding)
        return encoding

    def _output_helper(self, out_filename, interspersed, size, is_double):
        # type: (Optional[str], bool, int, bool) -> str
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

        if PY2:
            if not (hasattr(out_filename, 'read') and hasattr(out_filename, 'write')
                   ) or isinstance(out_filename, (file, StringIO)):
                return out_filename
            elif not isinstance(out_filename, string_types):
                msg = 'out_filename=%r must be a string; type=%s' % (
                    out_filename, type(out_filename))
                raise TypeError(msg)
        else:
            if not(hasattr(out_filename, 'read') and hasattr(out_filename, 'write')
                  ) or isinstance(out_filename, io.IOBase):
                return out_filename
            elif not isinstance(out_filename, string_types):
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
        fname = print_filename(out_filename)
        self.log.debug("***writing %s" % fname)
        return out_filename

    def write_caero_model(self, caero_bdf_filename='caero.bdf'):
        # type: (str) -> None
        """write the CAERO cards as CQUAD4s that can be visualized"""
        self.log.debug('---starting BDF.write_caero_model of %s---' % caero_bdf_filename)
        with open(caero_bdf_filename, 'w') as bdf_file:
            #bdf_file.write('$ pyNastran: punch=True\n')
            bdf_file.write('CEND\n')
            bdf_file.write('BEGIN BULK\n')
            i = 1
            mid = 1
            bdf_file.write('MAT1,%s,3.0E7,,0.3\n' % mid)
            for aesurf_id, aesurf in iteritems(self.aesurf):
                #cid = aesurf.cid1
                bdf_file.write('PSHELL,%s,%s,0.1\n' % (aesurf_id, aesurf_id))
                #print(cid)
                #ax, ay, az = cid.i
                #bx, by, bz = cid.j
                #cx, cy, cz = cid.k
                #bdf_file.write('CORD2R,%s,,%s,%s,%s,%s,%s,%s\n' % (cid, ax, ay, az, bx, by, bz))
                #bdf_file.write(',%s,%s,%s\n' % (cx, cy, cz))
                #print(cid)
                #aesurf.elements

            for eid, caero in sorted(iteritems(self.caeros)):
                #assert eid != 1, 'CAERO eid=1 is reserved for non-flaps'
                scaero = str(caero).rstrip().split('\n')
                bdf_file.write('$ ' + '\n$ '.join(scaero) + '\n')
                points, elements = caero.panel_points_elements()
                npoints = points.shape[0]
                #nelements = elements.shape[0]
                for ipoint, point in enumerate(points):
                    x, y, z = point
                    bdf_file.write(print_card_8(['GRID', i+ipoint, None, x, y, z]))

                #pid = eid
                #mid = eid
                bdf_file.write('PSHELL,%s,%s,0.1\n' % (1, 1))
                bdf_file.write('MAT1,%s,3.0E7,,0.3\n' % 1)

                j = 0
                for elem in elements + i:
                    p1, p2, p3, p4 = elem
                    eid2 = j + eid
                    pidi = None
                    for aesurf_id, aesurf in iteritems(self.aesurf):
                        aelist_id = aesurf.aelist_id1()
                        aelist = self.aelists[aelist_id]
                        if eid2 in aelist.elements:
                            pidi = aesurf_id
                            break
                    if pidi is None:
                        #pidi = pid
                        pidi = 1
                    bdf_file.write(print_card_8(['CQUAD4', j + eid, pidi, p1, p2, p3, p4]))
                    j += 1
                i += npoints
            bdf_file.write('ENDDATA\n')

    def write_bdf(self, out_filename=None, encoding=None,
                  size=8, is_double=False,
                  interspersed=False, enddata=None, close=True):
        # type: (Optional[Union[str, StringIO]], Optional[str], int, bool, bool, Optional[bool], bool) -> None
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
        is_long_ids = False

        # required for MasterModelTaxi
        if self.is_bdf_vectorized:
            pass
        else:
            is_long_ids = (
                self.nodes and max(self.nodes) > 100000000 or
                self.coords and max(self.coords) > 100000000 or
                self.elements and max(self.elements) > 100000000 or
                self.properties and max(self.properties) > 100000000 or
                self.materials and max(self.materials) > 100000000 or
                self.thermal_materials and max(self.thermal_materials) > 100000000 or
                self.nsms and max(self.nsms) > 100000000 or
                self.nsmadds and max(self.nsmadds) > 100000000)
            if is_long_ids:
                size = 16

        #self.write_caero_model()
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
        self._write_params(bdf_file, size, is_double, is_long_ids=is_long_ids)
        self._write_nodes(bdf_file, size, is_double, is_long_ids=is_long_ids)

        if interspersed:
            self._write_elements_interspersed(bdf_file, size, is_double, is_long_ids=is_long_ids)
        else:
            self._write_elements(bdf_file, size, is_double, is_long_ids=is_long_ids)
            self._write_properties(bdf_file, size, is_double, is_long_ids=is_long_ids)
        self._write_materials(bdf_file, size, is_double, is_long_ids=is_long_ids)

        self._write_masses(bdf_file, size, is_double, is_long_ids=is_long_ids)
        self._write_rigid_elements(bdf_file, size, is_double, is_long_ids=is_long_ids) # split out for write_bdf_symmetric
        self._write_aero(bdf_file, size, is_double, is_long_ids=is_long_ids)  # split out for write_bdf_symmetric

        self._write_common(bdf_file, size, is_double, is_long_ids=is_long_ids)
        if (enddata is None and 'ENDDATA' in self.card_count) or enddata:
            bdf_file.write('ENDDATA\n')
        if close:
            bdf_file.close()

    def _write_header(self, bdf_file, encoding):
        # type: (Any, bool) -> None
        """
        Writes the executive and case control decks.
        """
        if self.punch is None:
            # writing a mesh without using read_bdf
            if self.system_command_lines or self.executive_control_lines or self.case_control_deck:
                self.punch = False
            else:
                self.punch = True

        if self.nastran_format:
            bdf_file.write('$pyNastran: version=%s\n' % self.nastran_format)
            bdf_file.write('$pyNastran: punch=%s\n' % self.punch)
            bdf_file.write('$pyNastran: encoding=%s\n' % encoding)
            bdf_file.write('$pyNastran: nnodes=%s\n' % len(self.nodes))
            bdf_file.write('$pyNastran: nelements=%s\n' % len(self.elements))

        if not self.punch:
            self._write_executive_control_deck(bdf_file)
            self._write_case_control_deck(bdf_file)

    def _write_executive_control_deck(self, bdf_file):
        # type: (Any) -> None
        """
        Writes the executive control deck.
        """
        msg = ''
        for line in self.system_command_lines:
            msg += line + '\n'

        if self.executive_control_lines:
            msg += '$EXECUTIVE CONTROL DECK\n'
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
        # type: (Any) -> None
        """
        Writes the Case Control Deck.
        """
        if self.case_control_deck:
            msg = '$CASE CONTROL DECK\n'
            msg += str(self.case_control_deck)
            assert 'BEGIN BULK' in msg, msg
            bdf_file.write(''.join(msg))

    def _write_elements(self, bdf_file, size=8, is_double=False, is_long_ids=None):
        # type: (Any, int, bool) -> None
        """
        Writes the elements in a sorted order
        """
        size, is_long_ids = self._write_mesh_long_ids_size(size, is_long_ids)
        if self.elements:
            bdf_file.write('$ELEMENTS\n')
            if is_long_ids:
                for (eid, element) in sorted(iteritems(self.elements)):
                    bdf_file.write(element.write_card_16(is_double))
            else:
                for (eid, element) in sorted(iteritems(self.elements)):
                    try:
                        bdf_file.write(element.write_card(size, is_double))
                    except:
                        print('failed printing element...'
                              'type=%s eid=%s' % (element.type, eid))
                        raise
        if self.ao_element_flags:
            for (eid, element) in sorted(iteritems(self.ao_element_flags)):
                bdf_file.write(element.write_card(size, is_double))
        if self.normals:
            for (nid, snorm) in sorted(iteritems(self.normals)):
                bdf_file.write(snorm.write_card(size, is_double))
        self._write_nsm(bdf_file, size, is_double)

    def _write_nsm(self, bdf_file, size=8, is_double=False, is_long_ids=None):
        # type: (Any, int, bool) -> None
        """
        Writes the nsm in a sorted order
        """
        if self.nsms or self.nsmadds:
            bdf_file.write('$NSM\n')
            for (unused_id, nsmadds) in sorted(iteritems(self.nsmadds)):
                for nsmadd in nsmadds:
                    bdf_file.write(str(nsmadd))
            for (key, nsms) in sorted(iteritems(self.nsms)):
                for nsm in nsms:
                    try:
                        bdf_file.write(nsm.write_card(size, is_double))
                    except:
                        print('failed printing nsm...type=%s key=%r'
                              % (nsm.type, key))
                        raise

    def _write_elements_interspersed(self, bdf_file, size=8, is_double=False, is_long_ids=None):
        # type: (Any, int, bool) -> None
        """
        Writes the elements and properties in and interspersed order
        """
        missing_properties = []
        if self.properties:
            bdf_file.write('$ELEMENTS_WITH_PROPERTIES\n')

        eids_written = []  # type: List[int]
        pids = sorted(self.properties.keys())
        pid_eids = self.get_element_ids_dict_with_pids(pids, stop_if_no_eids=False)

        #failed_element_types = set([])
        for (pid, eids) in sorted(iteritems(pid_eids)):
            prop = self.properties[pid]
            if eids:
                bdf_file.write(prop.write_card(size, is_double))
                eids.sort()
                for eid in eids:
                    element = self.elements[eid]
                    try:
                        bdf_file.write(element.write_card(size, is_double))
                    except:
                        print('failed printing element...' 'type=%r eid=%s'
                              % (element.type, eid))
                        raise
                eids_written += eids
            else:
                missing_properties.append(prop.write_card(size, is_double))

        eids_missing = set(self.elements.keys()).difference(set(eids_written))
        if eids_missing:
            bdf_file.write('$ELEMENTS_WITH_NO_PROPERTIES '
                           '(PID=0 and unanalyzed properties)\n')
            for eid in sorted(eids_missing):
                element = self.elements[eid]
                try:
                    bdf_file.write(element.write_card(size, is_double))
                except:
                    print('failed printing element...'
                          'type=%s eid=%s' % (element.type, eid))
                    raise

        if missing_properties or self.pdampt or self.pbusht or self.pelast:
            bdf_file.write('$UNASSOCIATED_PROPERTIES\n')
            for card in sorted(itervalues(self.pbusht)):
                bdf_file.write(card.write_card(size, is_double))
            for card in sorted(itervalues(self.pdampt)):
                bdf_file.write(card.write_card(size, is_double))
            for card in sorted(itervalues(self.pelast)):
                bdf_file.write(card.write_card(size, is_double))
            for card in missing_properties:
                # this is a string...
                #print("missing_property = ", card
                bdf_file.write(card)

        if self.ao_element_flags:
            for (eid, element) in sorted(iteritems(self.ao_element_flags)):
                bdf_file.write(element.write_card(size, is_double))
        if self.normals:
            for (nid, snorm) in sorted(iteritems(self.normals)):
                bdf_file.write(snorm.write_card(size, is_double))
        self._write_nsm(bdf_file, size, is_double)

    def _write_aero(self, bdf_file, size=8, is_double=False, is_long_ids=None):
        # type: (Any, int, bool) -> None
        """Writes the aero cards"""
        if self.caeros or self.paeros or self.monitor_points or self.splines:
            bdf_file.write('$AERO\n')
            for (unused_id, caero) in sorted(iteritems(self.caeros)):
                bdf_file.write(caero.write_card(size, is_double))
            for (unused_id, paero) in sorted(iteritems(self.paeros)):
                bdf_file.write(paero.write_card(size, is_double))
            for (unused_id, spline) in sorted(iteritems(self.splines)):
                bdf_file.write(spline.write_card(size, is_double))
            for monitor_point in self.monitor_points:
                bdf_file.write(monitor_point.write_card(size, is_double))
        self.zona.write_bdf(bdf_file, size=8, is_double=False)

    def _write_aero_control(self, bdf_file, size=8, is_double=False, is_long_ids=None):
        # type: (Any, int, bool) -> None
        """Writes the aero control surface cards"""
        if(self.aecomps or self.aefacts or self.aeparams or self.aelinks or
           self.aelists or self.aestats or self.aesurf or self.aesurfs):
            bdf_file.write('$AERO CONTROL SURFACES\n')
            for (unused_id, aelinks) in sorted(iteritems(self.aelinks)):
                for aelink in aelinks:
                    bdf_file.write(aelink.write_card(size, is_double))

            for (unused_id, aecomp) in sorted(iteritems(self.aecomps)):
                bdf_file.write(aecomp.write_card(size, is_double))
            for (unused_id, aeparam) in sorted(iteritems(self.aeparams)):
                bdf_file.write(aeparam.write_card(size, is_double))
            for (unused_id, aestat) in sorted(iteritems(self.aestats)):
                bdf_file.write(aestat.write_card(size, is_double))

            for (unused_id, aelist) in sorted(iteritems(self.aelists)):
                bdf_file.write(aelist.write_card(size, is_double))
            for (unused_id, aesurf) in sorted(iteritems(self.aesurf)):
                bdf_file.write(aesurf.write_card(size, is_double))
            for (unused_id, aesurfs) in sorted(iteritems(self.aesurfs)):
                bdf_file.write(aesurfs.write_card(size, is_double))
            for (unused_id, aefact) in sorted(iteritems(self.aefacts)):
                bdf_file.write(aefact.write_card(size, is_double))

    def _write_static_aero(self, bdf_file, size=8, is_double=False, is_long_ids=None):
        # type: (Any, int, bool) -> None
        """Writes the static aero cards"""
        if self.aeros or self.trims or self.divergs:
            bdf_file.write('$STATIC AERO\n')
            # static aero
            if self.aeros:
                bdf_file.write(self.aeros.write_card(size, is_double))
            for (unused_id, trim) in sorted(iteritems(self.trims)):
                bdf_file.write(trim.write_card(size, is_double))
            for (unused_id, diverg) in sorted(iteritems(self.divergs)):
                bdf_file.write(diverg.write_card(size, is_double))

    def _find_aero_location(self):
        # type: () -> Tuple[bool, bool]
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

    def _write_flutter(self, bdf_file, size=8, is_double=False, write_aero_in_flutter=True, is_long_ids=None):
        # type: (Any, int, bool, bool) -> None
        """Writes the flutter cards"""
        if (write_aero_in_flutter and self.aero) or self.flfacts or self.flutters or self.mkaeros:
            bdf_file.write('$FLUTTER\n')
            if write_aero_in_flutter:
                bdf_file.write(self.aero.write_card(size, is_double))
            for (unused_id, flutter) in sorted(iteritems(self.flutters)):
                bdf_file.write(flutter.write_card(size, is_double))
            for (unused_id, flfact) in sorted(iteritems(self.flfacts)):
                bdf_file.write(flfact.write_card(size, is_double))
            for mkaero in self.mkaeros:
                bdf_file.write(mkaero.write_card(size, is_double))

    def _write_gust(self, bdf_file, size=8, is_double=False, write_aero_in_gust=True, is_long_ids=None):
        # type: (Any, int, bool, bool) -> None
        """Writes the gust cards"""
        if (write_aero_in_gust and self.aero) or self.gusts:
            bdf_file.write('$GUST\n')
            if write_aero_in_gust:
                for (unused_id, aero) in sorted(iteritems(self.aero)):
                    bdf_file.write(aero.write_card(size, is_double))
            for (unused_id, gust) in sorted(iteritems(self.gusts)):
                bdf_file.write(gust.write_card(size, is_double))

    def _write_common(self, bdf_file, size=8, is_double=False, is_long_ids=None):
        # type: (Any, int, bool) -> None
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

        """
        self._write_dmigs(bdf_file, size, is_double, is_long_ids=is_long_ids)
        self._write_loads(bdf_file, size, is_double, is_long_ids=is_long_ids)
        self._write_dynamic(bdf_file, size, is_double, is_long_ids=is_long_ids)
        self._write_aero_control(bdf_file, size, is_double, is_long_ids=is_long_ids)
        self._write_static_aero(bdf_file, size, is_double, is_long_ids=is_long_ids)

        write_aero_in_flutter, write_aero_in_gust = self._find_aero_location()
        self._write_flutter(bdf_file, size, is_double, write_aero_in_flutter, is_long_ids=is_long_ids)
        self._write_gust(bdf_file, size, is_double, write_aero_in_gust, is_long_ids=is_long_ids)

        self._write_thermal(bdf_file, size, is_double, is_long_ids=is_long_ids)
        self._write_thermal_materials(bdf_file, size, is_double, is_long_ids=is_long_ids)

        self._write_constraints(bdf_file, size, is_double, is_long_ids=is_long_ids)
        self._write_optimization(bdf_file, size, is_double, is_long_ids=is_long_ids)
        self._write_tables(bdf_file, size, is_double, is_long_ids=is_long_ids)
        self._write_sets(bdf_file, size, is_double, is_long_ids=is_long_ids)
        self._write_superelements(bdf_file, size, is_double, is_long_ids=is_long_ids)
        self._write_contact(bdf_file, size, is_double, is_long_ids=is_long_ids)
        self._write_rejects(bdf_file, size, is_double, is_long_ids=is_long_ids)
        self._write_coords(bdf_file, size, is_double, is_long_ids=is_long_ids)

    def _write_constraints(self, bdf_file, size=8, is_double=False, is_long_ids=None):
        # type: (Any, int, bool) -> None
        """Writes the constraint cards sorted by ID"""
        size, is_long_ids = self._write_mesh_long_ids_size(size, is_long_ids)
        if self.suport or self.suport1:
            bdf_file.write('$CONSTRAINTS\n')
            for suport in self.suport:
                bdf_file.write(suport.write_card(size, is_double))
            for unused_suport_id, suport in sorted(iteritems(self.suport1)):
                bdf_file.write(suport.write_card(size, is_double))

        if self.spcs or self.spcadds or self.spcoffs:
            #bdf_file.write('$SPCs\n')
            #str_spc = str(self.spcObject) # old
            #if str_spc:
                #bdf_file.write(str_spc)
            #else:
            bdf_file.write('$SPCs\n')
            for (unused_id, spcadds) in sorted(iteritems(self.spcadds)):
                for spcadd in spcadds:
                    bdf_file.write(str(spcadd))
            for (unused_id, spcs) in sorted(iteritems(self.spcs)):
                for spc in spcs:
                    bdf_file.write(str(spc))
            for (unused_id, spcoffs) in sorted(iteritems(self.spcoffs)):
                for spc in spcoffs:
                    bdf_file.write(str(spc))

        if self.mpcs or self.mpcadds:
            bdf_file.write('$MPCs\n')
            for (unused_id, mpcadds) in sorted(iteritems(self.mpcadds)):
                for mpcadd in mpcadds:
                    bdf_file.write(str(mpcadd))
            for (unused_id, mpcs) in sorted(iteritems(self.mpcs)):
                for mpc in mpcs:
                    bdf_file.write(mpc.write_card(size, is_double))

    def _write_contact(self, bdf_file, size=8, is_double=False, is_long_ids=None):
        # type: (Any, int, bool) -> None
        """Writes the contact cards sorted by ID"""
        is_contact = (self.bcrparas or self.bctadds or self.bctparas
                      or self.bctsets or self.bsurf or self.bsurfs)
        if is_contact:
            bdf_file.write('$CONTACT\n')
            for (unused_id, bcrpara) in sorted(iteritems(self.bcrparas)):
                bdf_file.write(bcrpara.write_card(size, is_double))
            for (unused_id, bctadds) in sorted(iteritems(self.bctadds)):
                bdf_file.write(bctadds.write_card(size, is_double))
            for (unused_id, bctpara) in sorted(iteritems(self.bctparas)):
                bdf_file.write(bctpara.write_card(size, is_double))

            for (unused_id, bctset) in sorted(iteritems(self.bctsets)):
                bdf_file.write(bctset.write_card(size, is_double))
            for (unused_id, bsurfi) in sorted(iteritems(self.bsurf)):
                bdf_file.write(bsurfi.write_card(size, is_double))
            for (unused_id, bsurfsi) in sorted(iteritems(self.bsurfs)):
                bdf_file.write(bsurfsi.write_card(size, is_double))

    def _write_coords(self, bdf_file, size=8, is_double=False, is_long_ids=None):
        # type: (Any, int, bool) -> None
        """Writes the coordinate cards in a sorted order"""
        size, is_long_ids = self._write_mesh_long_ids_size(size, is_long_ids)
        if len(self.coords) > 1:
            bdf_file.write('$COORDS\n')
        for (unused_id, coord) in sorted(iteritems(self.coords)):
            if unused_id != 0:
                try:
                    bdf_file.write(coord.write_card(size, is_double))
                except RuntimeError:
                    bdf_file.write(coord.write_card(16, is_double))

    def _write_dmigs(self, bdf_file, size=8, is_double=False, is_long_ids=None):
        # type: (Any, int, bool) -> None
        """
        Writes the DMIG cards

        Parameters
        ----------
        size : int
            large field (16) or small field (8)

        """
        for (unused_name, dmig) in sorted(iteritems(self.dmigs)):
            bdf_file.write(dmig.write_card(size, is_double))
        for (unused_name, dmi) in sorted(iteritems(self.dmis)):
            bdf_file.write(dmi.write_card(size, is_double))
        for (unused_name, dmij) in sorted(iteritems(self.dmijs)):
            bdf_file.write(dmij.write_card(size, is_double))
        for (unused_name, dmiji) in sorted(iteritems(self.dmijis)):
            bdf_file.write(dmiji.write_card(size, is_double))
        for (unused_name, dmik) in sorted(iteritems(self.dmiks)):
            bdf_file.write(dmik.write_card(size, is_double))

    def _write_dynamic(self, bdf_file, size=8, is_double=False, is_long_ids=None):
        # type: (Any, int, bool) -> None
        """Writes the dynamic cards sorted by ID"""
        is_dynamic = (self.dareas or self.dphases or self.nlparms or self.frequencies or
                      self.methods or self.cMethods or self.tsteps or self.tstepnls or
                      self.transfer_functions or self.delays or self.rotors or self.tics or
                      self.nlpcis)
        if is_dynamic:
            bdf_file.write('$DYNAMIC\n')
            for (unused_id, method) in sorted(iteritems(self.methods)):
                bdf_file.write(method.write_card(size, is_double))
            for (unused_id, cmethod) in sorted(iteritems(self.cMethods)):
                bdf_file.write(cmethod.write_card(size, is_double))
            for (unused_id, darea) in sorted(iteritems(self.dareas)):
                bdf_file.write(darea.write_card(size, is_double))
            for (unused_id, dphase) in sorted(iteritems(self.dphases)):
                bdf_file.write(dphase.write_card(size, is_double))
            for (unused_id, nlparm) in sorted(iteritems(self.nlparms)):
                bdf_file.write(nlparm.write_card(size, is_double))
            for (unused_id, nlpci) in sorted(iteritems(self.nlpcis)):
                bdf_file.write(nlpci.write_card(size, is_double))
            for (unused_id, tstep) in sorted(iteritems(self.tsteps)):
                bdf_file.write(tstep.write_card(size, is_double))
            for (unused_id, tstepnl) in sorted(iteritems(self.tstepnls)):
                bdf_file.write(tstepnl.write_card(size, is_double))
            for (unused_id, freqs) in sorted(iteritems(self.frequencies)):
                for freq in freqs:
                    bdf_file.write(freq.write_card(size, is_double))
            for (unused_id, delay) in sorted(iteritems(self.delays)):
                bdf_file.write(delay.write_card(size, is_double))
            for (unused_id, rotor) in sorted(iteritems(self.rotors)):
                bdf_file.write(rotor.write_card(size, is_double))
            for (unused_id, tic) in sorted(iteritems(self.tics)):
                bdf_file.write(tic.write_card(size, is_double))

            for (unused_id, tfs) in sorted(iteritems(self.transfer_functions)):
                for transfer_function in tfs:
                    bdf_file.write(transfer_function.write_card(size, is_double))

    def _write_mesh_long_ids_size(self, size, is_long_ids):
        """helper method"""
        if is_long_ids and size == 16 or is_long_ids is False:
            return size, is_long_ids

        if size == 16 and is_long_ids is None or self.is_long_ids:
            size = 16
            is_long_ids = True
        else:
            is_long_ids = False
        return size, is_long_ids

    def _write_loads(self, bdf_file, size=8, is_double=False, is_long_ids=None):
        # type: (Any, int, bool) -> None
        """Writes the load cards sorted by ID"""
        size, is_long_ids = self._write_mesh_long_ids_size(size, is_long_ids)
        if self.load_combinations or self.loads or self.tempds:
            bdf_file.write('$LOADS\n')
            for (key, load_combinations) in sorted(iteritems(self.load_combinations)):
                for load_combination in load_combinations:
                        try:
                            bdf_file.write(load_combination.write_card(size, is_double))
                        except:
                            print('failed printing load...type=%s key=%r'
                                  % (load_combination.type, key))
                            raise
            for (key, loadcase) in sorted(iteritems(self.loads)):
                for load in loadcase:
                    try:
                        bdf_file.write(load.write_card(size, is_double))
                    except:
                        print('failed printing load...type=%s key=%r'
                              % (load.type, key))
                        raise
            for unused_key, tempd in sorted(iteritems(self.tempds)):
                bdf_file.write(tempd.write_card(size, is_double))
        self._write_dloads(bdf_file, size=size, is_double=is_double, is_long_ids=is_long_ids)

    def _write_dloads(self, bdf_file, size=8, is_double=False, is_long_ids=None):
    # type: (Any, int, bool) -> None
        """Writes the dload cards sorted by ID"""
        size, is_long_ids = self._write_mesh_long_ids_size(size, is_long_ids)
        if self.dloads or self.dload_entries:
            bdf_file.write('$DLOADS\n')
            for (key, loadcase) in sorted(iteritems(self.dloads)):
                for load in loadcase:
                    try:
                        bdf_file.write(load.write_card(size, is_double))
                    except:
                        print('failed printing load...type=%s key=%r'
                              % (load.type, key))
                        raise

            for (key, loadcase) in sorted(iteritems(self.dload_entries)):
                for load in loadcase:
                    try:
                        bdf_file.write(load.write_card(size, is_double))
                    except:
                        print('failed printing load...type=%s key=%r'
                              % (load.type, key))
                        raise


    def _write_masses(self, bdf_file, size=8, is_double=False, is_long_ids=None):
        # type: (Any, int, bool) -> None
        """Writes the mass cards sorted by ID"""
        size, is_long_ids = self._write_mesh_long_ids_size(size, is_long_ids)
        if self.properties_mass:
            bdf_file.write('$PROPERTIES_MASS\n')
            for (pid, mass) in sorted(iteritems(self.properties_mass)):
                try:
                    bdf_file.write(mass.write_card(size, is_double))
                except:
                    print('failed printing mass property...'
                          'type=%s eid=%s' % (mass.type, pid))
                    raise

        if self.masses:
            bdf_file.write('$MASSES\n')
            for (eid, mass) in sorted(iteritems(self.masses)):
                try:
                    bdf_file.write(mass.write_card(size, is_double))
                except:
                    print('failed printing masses...'
                          'type=%s eid=%s' % (mass.type, eid))
                    raise

    def _write_materials(self, bdf_file, size=8, is_double=False, is_long_ids=None):
        # type: (Any, int, bool) -> None
        """Writes the materials in a sorted order"""
        size, is_long_ids = self._write_mesh_long_ids_size(size, is_long_ids)
        is_materials = (self.materials or self.hyperelastic_materials or self.creep_materials or
                        self.MATS1 or self.MATS3 or self.MATS8 or self.MATT1 or
                        self.MATT2 or self.MATT3 or self.MATT4 or self.MATT5 or
                        self.MATT8 or self.MATT9 or self.nxstrats)
        if is_materials:
            bdf_file.write('$MATERIALS\n')
            for (unused_mid, material) in sorted(iteritems(self.materials)):
                bdf_file.write(material.write_card(size, is_double))
            for (unused_mid, material) in sorted(iteritems(self.hyperelastic_materials)):
                bdf_file.write(material.write_card(size, is_double))
            for (unused_mid, material) in sorted(iteritems(self.creep_materials)):
                bdf_file.write(material.write_card(size, is_double))

            for (unused_mid, material) in sorted(iteritems(self.MATS1)):
                bdf_file.write(material.write_card(size, is_double))
            for (unused_mid, material) in sorted(iteritems(self.MATS3)):
                bdf_file.write(material.write_card(size, is_double))
            for (unused_mid, material) in sorted(iteritems(self.MATS8)):
                bdf_file.write(material.write_card(size, is_double))

            for (unused_mid, material) in sorted(iteritems(self.MATT1)):
                bdf_file.write(material.write_card(size, is_double))
            for (unused_mid, material) in sorted(iteritems(self.MATT2)):
                bdf_file.write(material.write_card(size, is_double))
            for (unused_mid, material) in sorted(iteritems(self.MATT3)):
                bdf_file.write(material.write_card(size, is_double))
            for (unused_mid, material) in sorted(iteritems(self.MATT4)):
                bdf_file.write(material.write_card(size, is_double))
            for (unused_mid, material) in sorted(iteritems(self.MATT5)):
                bdf_file.write(material.write_card(size, is_double))
            for (unused_mid, material) in sorted(iteritems(self.MATT8)):
                bdf_file.write(material.write_card(size, is_double))
            for (unused_mid, material) in sorted(iteritems(self.MATT9)):
                bdf_file.write(material.write_card(size, is_double))
            for (unused_sid, nxstrat) in sorted(iteritems(self.nxstrats)):
                bdf_file.write(nxstrat.write_card(size, is_double))

    def _write_nodes(self, bdf_file, size=8, is_double=False, is_long_ids=None):
        # type: (Any, int, bool) -> None
        """Writes the NODE-type cards"""
        if self.spoints:
            bdf_file.write('$SPOINTS\n')
            bdf_file.write(write_xpoints('SPOINT', self.spoints))
        if self.epoints:
            bdf_file.write('$EPOINTS\n')
            bdf_file.write(write_xpoints('EPOINT', self.epoints))
        if self.points:
            bdf_file.write('$POINTS\n')
            for unused_point_id, point in sorted(iteritems(self.points)):
                bdf_file.write(point.write_card(size, is_double))
        if self.axic:
            bdf_file.write(self.axic.write_card(size, is_double))
            for unused_nid, ringax_pointax in iteritems(self.ringaxs):
                bdf_file.write(ringax_pointax.write_card(size, is_double))

        self._write_grids(bdf_file, size=size, is_double=is_double)
        if self.seqgp:
            bdf_file.write(self.seqgp.write_card(size, is_double))

        #if 0:  # not finished
            #self._write_nodes_associated(bdf_file, size, is_double)

    def _write_grids(self, bdf_file, size=8, is_double=False, is_long_ids=None):
        # type: (Any, int, bool) -> None
        """Writes the GRID-type cards"""
        size, is_long_ids = self._write_mesh_long_ids_size(size, is_long_ids)
        if self.nodes:
            bdf_file.write('$NODES\n')
            if self.grdset:
                bdf_file.write(self.grdset.write_card(size))
            if is_long_ids:
                for (unused_nid, node) in sorted(iteritems(self.nodes)):
                    bdf_file.write(node.write_card_16(is_double))
            else:
                for (unused_nid, node) in sorted(iteritems(self.nodes)):
                    bdf_file.write(node.write_card(size, is_double))

    #def _write_nodes_associated(self, bdf_file, size=8, is_double=False):
        #"""
        #Writes the NODE-type in associated and unassociated groups.

        #.. warning:: Sometimes crashes, probably on invalid BDFs.
        #"""
        #associated_nodes = set([])
        #for (eid, element) in iteritems(self.elements):
            #associated_nodes = associated_nodes.union(set(element.node_ids))

        #all_nodes = set(self.nodes.keys())
        #unassociated_nodes = list(all_nodes.difference(associated_nodes))
        ##missing_nodes = all_nodes.difference(

        ## TODO: this really shouldn't be a list...???
        #associated_nodes = list(associated_nodes)

        #if associated_nodes:
            #bdf_file.write('$ASSOCIATED NODES\n')
            #if self.grdset:
                #bdf_file.write(self.grdset.write_card(size, is_double))
            ## TODO: this really shouldn't be a dictionary...???
            #for key, node in sorted(iteritems(associated_nodes)):
                #bdf_file.write(node.write_card(size, is_double))

        #if unassociated_nodes:
            #bdf_file.write('$UNASSOCIATED NODES\n')
            #if self.grdset and not associated_nodes:
                #v(self.grdset.write_card(size, is_double))
            #for key, node in sorted(iteritems(unassociated_nodes)):
                #if key in self.nodes:
                    #bdf_file.write(node.write_card(size, is_double))
                #else:
                    #bdf_file.write('$ Missing NodeID=%s' % key)

    def _write_optimization(self, bdf_file, size=8, is_double=False, is_long_ids=None):
        # type: (Any, int, bool) -> None
        """Writes the optimization cards sorted by ID"""
        is_optimization = (self.dconadds or self.dconstrs or self.desvars or self.ddvals or
                           self.dresps or
                           self.dvprels or self.dvmrels or self.dvcrels or self.doptprm or
                           self.dlinks or self.dequations or self.dtable is not None or
                           self.dvgrids or self.dscreen)
        if is_optimization:
            bdf_file.write('$OPTIMIZATION\n')
            for (unused_id, dconadd) in sorted(iteritems(self.dconadds)):
                bdf_file.write(dconadd.write_card(size, is_double))
            for (unused_id, dconstrs) in sorted(iteritems(self.dconstrs)):
                for dconstr in dconstrs:
                    bdf_file.write(dconstr.write_card(size, is_double))
            for (unused_id, desvar) in sorted(iteritems(self.desvars)):
                bdf_file.write(desvar.write_card(size, is_double))
            for (unused_id, ddval) in sorted(iteritems(self.ddvals)):
                bdf_file.write(ddval.write_card(size, is_double))
            for (unused_id, dlink) in sorted(iteritems(self.dlinks)):
                bdf_file.write(dlink.write_card(size, is_double))
            for (unused_id, dresp) in sorted(iteritems(self.dresps)):
                bdf_file.write(dresp.write_card(size, is_double))

            for (unused_id, dvcrel) in sorted(iteritems(self.dvcrels)):
                bdf_file.write(dvcrel.write_card(size, is_double))
            for (unused_id, dvmrel) in sorted(iteritems(self.dvmrels)):
                bdf_file.write(dvmrel.write_card(size, is_double))
            for (unused_id, dvprel) in sorted(iteritems(self.dvprels)):
                bdf_file.write(dvprel.write_card(size, is_double))
            for (unused_id, dvgrids) in sorted(iteritems(self.dvgrids)):
                for dvgrid in dvgrids:
                    bdf_file.write(dvgrid.write_card(size, is_double))
            for (unused_id, dscreen) in sorted(iteritems(self.dscreen)):
                bdf_file.write(str(dscreen))

            for (unused_id, equation) in sorted(iteritems(self.dequations)):
                bdf_file.write(str(equation))

            if self.dtable is not None:
                bdf_file.write(self.dtable.write_card(size, is_double))
            if self.doptprm is not None:
                bdf_file.write(self.doptprm.write_card(size, is_double))

    def _write_params(self, bdf_file, size=8, is_double=False, is_long_ids=None):
        # type: (Any, int, bool) -> None
        """
        Writes the PARAM cards
        """
        size, is_long_ids = self._write_mesh_long_ids_size(size, is_long_ids)
        if self.params or self.dti:
            bdf_file.write('$PARAMS\n')  # type: List[str]
            for unused_name, dti in sorted(iteritems(self.dti)):
                bdf_file.write(dti.write_card(size=size, is_double=is_double))

            for (unused_key, param) in sorted(iteritems(self.params)):
                bdf_file.write(param.write_card(size, is_double))

    def _write_properties(self, bdf_file, size=8, is_double=False, is_long_ids=None):
        # type: (Any, int, bool) -> None
        """Writes the properties in a sorted order"""
        size, is_long_ids = self._write_mesh_long_ids_size(size, is_long_ids)
        if self.properties:
            bdf_file.write('$PROPERTIES\n')
            prop_groups = (self.properties, self.pelast, self.pdampt, self.pbusht)
            if is_long_ids:
                for prop_group in prop_groups:
                    for unused_pid, prop in sorted(iteritems(prop_group)):
                        bdf_file.write(prop.write_card_16(is_double))
                #except:
                    #print('failed printing property type=%s' % prop.type)
                    #raise
            else:
                for prop_group in prop_groups:
                    for unused_pid, prop in sorted(iteritems(prop_group)):
                        bdf_file.write(prop.write_card(size, is_double))

    def _write_rejects(self, bdf_file, size=8, is_double=False, is_long_ids=None):
        # type: (Any, int, bool) -> None
        """
        Writes the rejected (processed) cards and the rejected unprocessed
        cardlines
        """
        if size == 8:
            print_func = print_card_8
        else:
            print_func = print_card_16

        if self.reject_cards:
            bdf_file.write('$REJECT_CARDS\n')
            for reject_card in self.reject_cards:
                try:
                    bdf_file.write(print_func(reject_card))
                except RuntimeError:
                    for field in reject_card:
                        if field is not None and '=' in field:
                            raise SyntaxError('cannot reject equal signed '
                                              'cards\ncard=%s\n' % reject_card)
                    raise

        if self.reject_lines:
            bdf_file.write('$REJECT_LINES\n')

            for reject_lines in self.reject_lines:
                if isinstance(reject_lines, (list, tuple)):
                    for reject in reject_lines:
                        reject2 = reject.rstrip()
                        if reject2:
                            bdf_file.write('%s\n' % reject2)
                elif isinstance(reject_lines, string_types):
                    reject2 = reject_lines.rstrip()
                    if reject2:
                        bdf_file.write('%s\n' % reject2)
                else:
                    raise TypeError(reject_lines)

    def _write_rigid_elements(self, bdf_file, size=8, is_double=False, is_long_ids=None):
        # type: (Any, int, bool) -> None
        """Writes the rigid elements in a sorted order"""
        size, is_long_ids = self._write_mesh_long_ids_size(size, is_long_ids)
        if self.rigid_elements:
            bdf_file.write('$RIGID ELEMENTS\n')
            if is_long_ids:
                for (eid, element) in sorted(iteritems(self.rigid_elements)):
                    try:
                        bdf_file.write(element.write_card_16(is_double))
                    except:
                        print('failed printing element...'
                              'type=%s eid=%s' % (element.type, eid))
                        raise
            else:
                for (eid, element) in sorted(iteritems(self.rigid_elements)):
                    try:
                        bdf_file.write(element.write_card(size, is_double))
                    except:
                        print('failed printing element...'
                              'type=%s eid=%s' % (element.type, eid))
                        raise
        if self.plotels:
            bdf_file.write('$PLOT ELEMENTS\n')
            for (eid, element) in sorted(iteritems(self.plotels)):
                bdf_file.write(element.write_card(size, is_double))

    def _write_sets(self, bdf_file, size=8, is_double=False, is_long_ids=None):
        # type: (Any, int, bool) -> None
        """Writes the SETx cards sorted by ID"""
        is_sets = (self.sets or self.asets or self.omits or self.bsets or self.csets or self.qsets
                   or self.usets)
        if is_sets:
            bdf_file.write('$SETS\n')  # type: List[str]
            for (unused_id, set_obj) in sorted(iteritems(self.sets)):  # dict
                bdf_file.write(set_obj.write_card(size, is_double))
            for set_obj in self.asets:  # list
                bdf_file.write(set_obj.write_card(size, is_double))
            for set_obj in self.omits:  # list
                bdf_file.write(set_obj.write_card(size, is_double))
            for set_obj in self.bsets:  # list
                bdf_file.write(set_obj.write_card(size, is_double))
            for set_obj in self.csets:  # list
                bdf_file.write(set_obj.write_card(size, is_double))
            for set_obj in self.qsets:  # list
                bdf_file.write(set_obj.write_card(size, is_double))
            for unused_name, usets in sorted(iteritems(self.usets)):  # dict
                for set_obj in usets:  # list
                    bdf_file.write(set_obj.write_card(size, is_double))

    def _write_superelements(self, bdf_file, size=8, is_double=False, is_long_ids=None):
        # type: (Any, int, bool) -> None
        """Writes the SETx cards sorted by ID"""
        is_sets = (self.se_sets or self.se_bsets or self.se_csets or self.se_qsets
                   or self.se_usets)
        if is_sets:
            bdf_file.write('$SUPERELEMENTS\n')  # type: List[str]
            for set_obj in self.se_bsets:  # list
                bdf_file.write(set_obj.write_card(size, is_double))
            for set_obj in self.se_csets:  # list
                bdf_file.write(set_obj.write_card(size, is_double))
            for set_obj in self.se_qsets:  # list
                bdf_file.write(set_obj.write_card(size, is_double))
            for (unused_set_id, set_obj) in sorted(iteritems(self.se_sets)):  # dict
                bdf_file.write(set_obj.write_card(size, is_double))
            for unused_name, usets in sorted(iteritems(self.se_usets)):  # dict
                for set_obj in usets:  # list
                    bdf_file.write(set_obj.write_card(size, is_double))
            for suport in self.se_suport:  # list
                bdf_file.write(suport.write_card(size, is_double))

    def _write_tables(self, bdf_file, size=8, is_double=False, is_long_ids=None):
        # type: (Any, int, bool) -> None
        """Writes the TABLEx cards sorted by ID"""
        if self.tables or self.tables_d or self.tables_m or self.tables_sdamping:
            bdf_file.write('$TABLES\n')  # type: List[str]
            for (unused_id, table) in sorted(iteritems(self.tables)):
                bdf_file.write(table.write_card(size, is_double))
            for (unused_id, table) in sorted(iteritems(self.tables_d)):
                bdf_file.write(table.write_card(size, is_double))
            for (unused_id, table) in sorted(iteritems(self.tables_m)):
                bdf_file.write(table.write_card(size, is_double))
            for (unused_id, table) in sorted(iteritems(self.tables_sdamping)):
                bdf_file.write(table.write_card(size, is_double))

        if self.random_tables:
            bdf_file.write('$RANDOM TABLES\n')
            for (unused_id, table) in sorted(iteritems(self.random_tables)):
                bdf_file.write(table.write_card(size, is_double))

    def _write_thermal(self, bdf_file, size=8, is_double=False, is_long_ids=None):
        # type: (Any, int, bool) -> None
        """Writes the thermal cards"""
        # PHBDY
        is_thermal = (self.phbdys or self.convection_properties or self.bcs or
                      self.views or self.view3ds or self.radset or self.radcavs)
        if is_thermal:
            bdf_file.write('$THERMAL\n')

            for (unused_key, phbdy) in sorted(iteritems(self.phbdys)):
                bdf_file.write(phbdy.write_card(size, is_double))

            #for unused_key, prop in sorted(iteritems(self.thermal_properties)):
            #    bdf_file.write(str(prop))
            for (unused_key, prop) in sorted(iteritems(self.convection_properties)):
                bdf_file.write(prop.write_card(size, is_double))

            # BCs
            for (unused_key, bcs) in sorted(iteritems(self.bcs)):
                for boundary_condition in bcs:  # list
                    bdf_file.write(boundary_condition.write_card(size, is_double))

            for (unused_key, view) in sorted(iteritems(self.views)):
                bdf_file.write(view.write_card(size, is_double))
            for (unused_key, view3d) in sorted(iteritems(self.view3ds)):
                bdf_file.write(view3d.write_card(size, is_double))
            if self.radset:
                bdf_file.write(self.radset.write_card(size, is_double))
            for unused_icavity, radcav in iteritems(self.radcavs):
                bdf_file.write(radcav.write_card(size, is_double))


    def _write_thermal_materials(self, bdf_file, size=8, is_double=False, is_long_ids=None):
        # type: (Any, int, bool) -> None
        """Writes the thermal materials in a sorted order"""
        if self.thermal_materials:
            bdf_file.write('$THERMAL MATERIALS\n')
            for (unused_mid, material) in sorted(iteritems(self.thermal_materials)):
                bdf_file.write(material.write_card(size, is_double))
