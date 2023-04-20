# coding: utf-8
"""
This file defines:
  - WriteMesh

"""
from __future__ import annotations
import sys
from io import IOBase
from pathlib import PurePath
from typing import Union, Optional, Any, cast, TYPE_CHECKING

from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.field_writer_16 import print_card_16
from pyNastran.bdf.bdf_interface.attributes import BDFAttributes
from pyNastran.bdf.bdf_interface.write_mesh_utils import (
    find_aero_location, write_dict, get_properties_by_element_type)
from pyNastran.bdf.cards.nodes import write_xpoints

try:
    from natsort import natsorted
except ImportError:
    natsorted = sorted

if TYPE_CHECKING:  # pragma: no cover
    from io import StringIO
    from cpylog import SimpleLogger
    from pyNastran.bdf.bdf import BDF, DESVAR


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
        self.cards_to_read = set()

    def get_encoding(self, encoding: Optional[str]=None) -> str:
        """gets the file encoding"""
        if encoding is not None:
            pass
        else:
            encoding = self._encoding
            if encoding is None:
                encoding = sys.getdefaultencoding()
            elif isinstance(encoding, bytes):
                # needed for hdf5 loader for some reason...
                encoding = self._encoding.decode('latin1')
                self._encoding = encoding
        encoding = cast(str, encoding)  # just for typing
        assert isinstance(encoding, str), encoding
        return encoding

    def _get_long_ids(self, size) -> tuple[bool, int]:
        """calculate is_long_ids"""
        # required for MasterModelTaxi
        is_long_ids = (
            any((len(dicti) and max(dicti) > 100000000 for dicti in (
                self.nodes, self.coords, self.elements, self.properties,
                self.materials, self.thermal_materials,
                #self.loads,
                self.load_combinations, self.masses,
                self.nsms, self.nsmadds,
            ))))
        if not is_long_ids:
            for loads in self.loads.values():
                for load in loads:
                    if load.type == 'TEMP':
                        if max(load.temperatures) > 100000000:
                            is_long_ids = True
                            break
                        #print(load.get_stats())
        if is_long_ids:
            size = 16
        return is_long_ids, size

    def write_bdf(self, out_filename: Optional[Union[str, StringIO]]=None,
                  encoding: Optional[str]=None,
                  size: int=8,
                  nodes_size: Optional[int]=None,
                  elements_size: Optional[int]=None,
                  loads_size: Optional[int]=None,
                  is_double: bool=False,
                  interspersed: bool=False, enddata: Optional[bool]=None,
                  write_header: bool=True, close: bool=True) -> None:
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
        write_header : bool; default=True
            flag for writing the pyNastran header
        close : bool; default=True
            should the output file be closed

        """
        if self.is_bdf_vectorized:
            is_long_ids = False
        else:
            is_long_ids, size = self._get_long_ids(size)

        out_filename, size = _output_helper(
            out_filename, interspersed, size, is_double, self.log)
        encoding = self.get_encoding(encoding)
        #assert encoding.lower() in ['ascii', 'latin1', 'utf8'], encoding

        has_read_write = hasattr(out_filename, 'read') and hasattr(out_filename, 'write')
        if has_read_write:
            bdf_file = out_filename
        else:
            self.log.debug(f'---starting BDF.write_bdf of {out_filename}---')
            assert isinstance(encoding, str), encoding
            bdf_file = open(out_filename, 'w', encoding=encoding)
        self._write_header(bdf_file, encoding, write_header=write_header)


        if self.superelement_models:
            bdf_file.write('$' + '*'*80+'\n')
            for superelement_id, superelement in sorted(self.superelement_models.items()):
                bdf_file.write(f'BEGIN SUPER={superelement_id}\n')
                superelement.write_bdf(out_filename=bdf_file, encoding=encoding,
                                       size=size, is_double=is_double,
                                       interspersed=interspersed, enddata=False,
                                       write_header=False, close=False)
                bdf_file.write('$' + '*'*80+'\n')
            bdf_file.write('BEGIN BULK\n')
        self.write_bulk_data(bdf_file, size=size, is_double=is_double,
                             interspersed=interspersed,
                             enddata=enddata, close=close,
                             nodes_size=nodes_size, elements_size=elements_size, loads_size=loads_size,
                             is_long_ids=is_long_ids)

    def write_bulk_data(self, bdf_file,
                        size: int=8, is_double: bool=False,
                        interspersed: bool=False,
                        enddata: Optional[bool]=None, close: bool=True,
                        nodes_size: Optional[int]=None,
                        elements_size: Optional[int]=None,
                        loads_size: Optional[int]=None,
                        is_long_ids: bool=False) -> None:
        """
        Writes the BDF.

        Parameters
        ----------
        bdf_file : varies
            file       - a file object
            StringIO() - a StringIO object
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

        .. note:: is_long_ids is only needed if you have ids longer
                  than 8 characters. It's an internal parameter, but if
                  you're calling the new sub-function, you might need
                  it.  Chances are you won't.
        """
        size, nodes_size, elements_size, loads_size = _fix_sizes(
            size, nodes_size, elements_size, loads_size)

        self._write_params(bdf_file, size, is_double, is_long_ids=is_long_ids)
        self._write_nodes(bdf_file, nodes_size, is_double, is_long_ids=is_long_ids)

        if interspersed:
            self._write_elements_interspersed(bdf_file, elements_size, is_double, is_long_ids=is_long_ids)
        else:
            self._write_elements(bdf_file, elements_size, is_double, is_long_ids=is_long_ids)
            self._write_properties(bdf_file, size, is_double, is_long_ids=is_long_ids)
            #self._write_properties_by_element_type(bdf_file, size, is_double, is_long_ids)

        for cards in (self.bolt, self.boltseq, self.boltfor, self.boltfrc, self.boltld):
            for key, card in cards.items():
                bdf_file.write(card.write_card(size, is_double))

        self._write_materials(bdf_file, size, is_double, is_long_ids=is_long_ids)

        self._write_masses(bdf_file, size, is_double, is_long_ids=is_long_ids)

        # split out for write_bdf_symmetric
        self._write_rigid_elements(bdf_file, size, is_double, is_long_ids=is_long_ids)
        self._write_aero(bdf_file, size, is_double, is_long_ids=is_long_ids)

        self._write_common(bdf_file, loads_size, is_double, is_long_ids=is_long_ids)
        if (enddata is None and 'ENDDATA' in self.card_count) or enddata:
            bdf_file.write('ENDDATA\n')
        if close:
            bdf_file.close()

    def _write_header(self, bdf_file: Any, encoding: str, write_header: bool=True) -> None:
        """Writes the executive and case control decks."""
        self._set_punch()

        if self.nastran_format and write_header:
            bdf_file.write(f'$pyNastran: version={self.nastran_format}\n')
            bdf_file.write(f'$pyNastran: punch={self.punch}\n')
            bdf_file.write(f'$pyNastran: encoding={encoding}\n')
            bdf_file.write(f'$pyNastran: nnodes={len(self.nodes):d}\n')
            bdf_file.write(f'$pyNastran: nelements={len(self.elements):d}\n')

        if not self.punch:
            self._write_executive_control_deck(bdf_file)
            self._write_case_control_deck(bdf_file)

    def _write_executive_control_deck(self, bdf_file: Any) -> None:
        """Writes the executive control deck."""
        msg = ''
        for line in self.system_command_lines:
            msg += line + '\n'

        if self.executive_control_lines:
            msg += '$EXECUTIVE CONTROL DECK\n'

            if self.sol_iline is not None:
                if self.sol == 600 and self.sol_method:
                    new_sol = f'SOL 600,{self.sol_method}'
                else:
                    new_sol = f'SOL {self.sol}'
                self.executive_control_lines[self.sol_iline] = new_sol

            for line in self.executive_control_lines:
                msg += line + '\n'
            bdf_file.write(msg)

    def _write_case_control_deck(self, bdf_file: Any) -> None:
        """Writes the Case Control Deck."""
        if self.case_control_deck:
            msg = '$CASE CONTROL DECK\n'
            if self.superelement_models:
                msg += self.case_control_deck.write(write_begin_bulk=False)
            else:
                msg += str(self.case_control_deck)
                assert 'BEGIN BULK' in msg, msg
            bdf_file.write(''.join(msg))

    def _write_elements(self, bdf_file: Any, size: int=8, is_double: bool=False,
                        is_long_ids: Optional[bool]=None) -> None:
        """Writes the elements in a sorted order"""
        size, is_long_ids = self._write_mesh_long_ids_size(size, is_long_ids)
        if self.elements:
            bdf_file.write('$ELEMENTS\n')
            if is_long_ids:
                for (eid, element) in sorted(self.elements.items()):
                    bdf_file.write(element.write_card_16(is_double))
            else:
                for (eid, element) in sorted(self.elements.items()):
                    try:
                        bdf_file.write(element.write_card(size, is_double))
                    except Exception:
                        print(f'failed printing element...type={element.type} eid={eid}')
                        raise
        if self.ao_element_flags:
            for (eid, element) in sorted(self.ao_element_flags.items()):
                bdf_file.write(element.write_card(size, is_double))
        if self.normals:
            for (unused_nid, snorm) in sorted(self.normals.items()):
                bdf_file.write(snorm.write_card(size, is_double))
        self._write_nsm(bdf_file, size, is_double)

    def _write_nsm(self, bdf_file: Any, size: int=8, is_double: bool=False,
                   is_long_ids: Optional[bool]=None) -> None:
        """Writes the nsm in a sorted order"""
        if self.nsms or self.nsmadds:
            bdf_file.write('$NSM\n')
            for (unused_id, nsmadds) in sorted(self.nsmadds.items()):
                for nsmadd in nsmadds:
                    bdf_file.write(str(nsmadd))
            for (key, nsms) in sorted(self.nsms.items()):
                for nsm in nsms:
                    try:
                        bdf_file.write(nsm.write_card(size, is_double))
                    except Exception:
                        print(f'failed printing nsm...type={nsm.type} key={key!r}')
                        raise

    def _write_elements_interspersed(self, bdf_file: Any, size: int=8, is_double: bool=False,
                                     is_long_ids: Optional[bool]=None) -> None:
        """Writes the elements and properties in and interspersed order"""
        missing_properties = []
        if self.properties:
            bdf_file.write('$ELEMENTS_WITH_PROPERTIES\n')

        eids_written = []  # type: list[int]
        pids = sorted(self.properties.keys())
        pid_eids = self.get_element_ids_dict_with_pids(pids, stop_if_no_eids=False)

        #failed_element_types = set()
        for (pid, eids) in sorted(pid_eids.items()):
            prop = self.properties[pid]
            if eids:
                bdf_file.write(prop.write_card(size, is_double))
                eids.sort()
                for eid in eids:
                    element = self.elements[eid]
                    try:
                        bdf_file.write(element.write_card(size, is_double))
                    except Exception:
                        print(f'failed printing element...type={element.type!r} eid={eid}')
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
                except Exception:
                    print(f'failed printing element...type={element.type} eid={eid}')
                    raise

        if missing_properties or self.pdampt or self.pbusht or self.pelast:
            bdf_file.write('$UNASSOCIATED_PROPERTIES\n')
            for card in sorted(self.pbusht.values()):
                bdf_file.write(card.write_card(size, is_double))
            for card in sorted(self.pdampt.values()):
                bdf_file.write(card.write_card(size, is_double))
            for card in sorted(self.pelast.values()):
                bdf_file.write(card.write_card(size, is_double))
            for card in missing_properties:
                # this is a string...
                #print("missing_property = ", card
                bdf_file.write(card)

        if self.ao_element_flags:
            for (eid, element) in sorted(self.ao_element_flags.items()):
                bdf_file.write(element.write_card(size, is_double))
        if self.normals:
            for (unused_nid, snorm) in sorted(self.normals.items()):
                bdf_file.write(snorm.write_card(size, is_double))
        self._write_nsm(bdf_file, size, is_double)

    def _write_aero(self, bdf_file: Any, size: int=8, is_double: bool=False,
                    is_long_ids: Optional[bool]=None) -> None:
        """Writes the aero cards"""
        if self.caeros or self.paeros or self.monitor_points or self.splines:
            bdf_file.write('$AERO\n')
            for (unused_id, caero) in sorted(self.caeros.items()):
                bdf_file.write(caero.write_card(size, is_double))
            for (unused_id, paero) in sorted(self.paeros.items()):
                bdf_file.write(paero.write_card(size, is_double))
            for (unused_id, spline) in sorted(self.splines.items()):
                bdf_file.write(spline.write_card(size, is_double))
            for monitor_point in self.monitor_points:
                bdf_file.write(monitor_point.write_card(size, is_double))
        self.zona.write_bdf(bdf_file, size=8, is_double=False)

    def _write_aero_control(self, bdf_file: Any, size: int=8, is_double: bool=False,
                            is_long_ids: Optional[bool]=None) -> None:
        """Writes the aero control surface cards"""
        if(self.aecomps or self.aefacts or self.aeparams or self.aelinks or
           self.aelists or self.aestats or self.aesurf or self.aesurfs):
            bdf_file.write('$AERO CONTROL SURFACES\n')
            for (unused_id, aelinks) in sorted(self.aelinks.items()):
                for aelink in aelinks:
                    bdf_file.write(aelink.write_card(size, is_double))

            for (unused_id, aecomp) in sorted(self.aecomps.items()):
                bdf_file.write(aecomp.write_card(size, is_double))
            for (unused_id, aeparam) in sorted(self.aeparams.items()):
                bdf_file.write(aeparam.write_card(size, is_double))
            for (unused_id, aestat) in sorted(self.aestats.items()):
                bdf_file.write(aestat.write_card(size, is_double))

            for (unused_id, aelist) in sorted(self.aelists.items()):
                bdf_file.write(aelist.write_card(size, is_double))
            for (unused_id, aesurf) in sorted(self.aesurf.items()):
                bdf_file.write(aesurf.write_card(size, is_double))
            for (unused_id, aesurfs) in sorted(self.aesurfs.items()):
                bdf_file.write(aesurfs.write_card(size, is_double))
            for (unused_id, aefact) in sorted(self.aefacts.items()):
                bdf_file.write(aefact.write_card(size, is_double))

    def _write_static_aero(self, bdf_file: Any, size: int=8, is_double: bool=False,
                           is_long_ids: Optional[bool]=None) -> None:
        """Writes the static aero cards"""
        if self.aeros or self.trims or self.divergs:
            bdf_file.write('$STATIC AERO\n')
            # static aero
            if self.aeros:
                bdf_file.write(self.aeros.write_card(size, is_double))
            for (unused_id, trim) in sorted(self.trims.items()):
                bdf_file.write(trim.write_card(size, is_double))
            for (unused_id, diverg) in sorted(self.divergs.items()):
                bdf_file.write(diverg.write_card(size, is_double))


    def _write_flutter(self, bdf_file: Any, size: int=8, is_double: bool=False,
                       write_aero_in_flutter: bool=True,
                       is_long_ids: Optional[bool]=None) -> None:
        """Writes the flutter cards"""
        if (write_aero_in_flutter and self.aero) or self.flfacts or self.flutters or self.mkaeros:
            bdf_file.write('$FLUTTER\n')
            if write_aero_in_flutter:
                bdf_file.write(self.aero.write_card(size, is_double))
            for (unused_id, flutter) in sorted(self.flutters.items()):
                bdf_file.write(flutter.write_card(size, is_double))
            for (unused_id, flfact) in sorted(self.flfacts.items()):
                bdf_file.write(flfact.write_card(size, is_double))
            for mkaero in self.mkaeros:
                bdf_file.write(mkaero.write_card(size, is_double))

    def _write_gust(self, bdf_file: Any, size: int=8, is_double: bool=False,
                    write_aero_in_gust: bool=True,
                    is_long_ids: Optional[bool]=None) -> None:
        """Writes the gust cards"""
        if (write_aero_in_gust and self.aero) or self.gusts:
            bdf_file.write('$GUST\n')
            if write_aero_in_gust:
                for (unused_id, aero) in sorted(self.aero.items()):
                    bdf_file.write(aero.write_card(size, is_double))
            for (unused_id, gust) in sorted(self.gusts.items()):
                bdf_file.write(gust.write_card(size, is_double))

    def _write_common(self, bdf_file: Any, size: int=8, is_double: bool=False,
                      is_long_ids: Optional[bool]=None) -> None:
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

        write_aero_in_flutter, write_aero_in_gust = find_aero_location(self)
        self._write_flutter(bdf_file, size, is_double, write_aero_in_flutter,
                            is_long_ids=is_long_ids)
        self._write_gust(bdf_file, size, is_double, write_aero_in_gust, is_long_ids=is_long_ids)

        self._write_thermal(bdf_file, size, is_double, is_long_ids=is_long_ids)
        self._write_thermal_materials(bdf_file, size, is_double, is_long_ids=is_long_ids)
        self._write_constraints(bdf_file, size, is_double, is_long_ids=is_long_ids)
        self._write_optimization(bdf_file, size, is_double, is_long_ids=is_long_ids)
        self._write_tables(bdf_file, size, is_double, is_long_ids=is_long_ids)
        self._write_sets(bdf_file, size, is_double, is_long_ids=is_long_ids)
        self._write_superelements(bdf_file, size, is_double, is_long_ids=is_long_ids)
        self._write_contact(bdf_file, size, is_double, is_long_ids=is_long_ids)
        self._write_parametric(bdf_file, size, is_double, is_long_ids=is_long_ids)
        self._write_rejects(bdf_file, size, is_double, is_long_ids=is_long_ids)
        self._write_coords(bdf_file, size, is_double, is_long_ids=is_long_ids)
        self._write_matcids(bdf_file, size, is_double, is_long_ids=is_long_ids)

        if self.acmodl:
            bdf_file.write(self.acmodl.write_card(size, is_double))


    def _write_constraints(self, bdf_file: Any, size: int=8, is_double: bool=False,
                           is_long_ids: Optional[bool]=None) -> None:
        """Writes the constraint cards sorted by ID"""
        size, is_long_ids = self._write_mesh_long_ids_size(size, is_long_ids)
        if self.suport or self.suport1:
            bdf_file.write('$CONSTRAINTS\n')
            for suport in self.suport:
                bdf_file.write(suport.write_card(size, is_double))
            for unused_suport_id, suport in sorted(self.suport1.items()):
                bdf_file.write(suport.write_card(size, is_double))

        if self.spcs or self.spcadds or self.spcoffs:
            #bdf_file.write('$SPCs\n')
            #str_spc = str(self.spcObject) # old
            #if str_spc:
                #bdf_file.write(str_spc)
            #else:
            bdf_file.write('$SPCs\n')
            for (unused_id, spcadds) in sorted(self.spcadds.items()):
                for spcadd in spcadds:
                    bdf_file.write(str(spcadd))
            for (unused_id, spcs) in sorted(self.spcs.items()):
                for spc in spcs:
                    bdf_file.write(str(spc))
            for (unused_id, spcoffs) in sorted(self.spcoffs.items()):
                for spc in spcoffs:
                    bdf_file.write(str(spc))

        if self.mpcs or self.mpcadds:
            bdf_file.write('$MPCs\n')
            for (unused_id, mpcadds) in sorted(self.mpcadds.items()):
                for mpcadd in mpcadds:
                    bdf_file.write(str(mpcadd))
            for (unused_id, mpcs) in sorted(self.mpcs.items()):
                for mpc in mpcs:
                    bdf_file.write(mpc.write_card(size, is_double))

    def _write_contact(self, bdf_file: Any, size: int=8, is_double: bool=False,
                       is_long_ids: Optional[bool]=None) -> None:
        """Writes the contact cards sorted by ID"""
        is_contact = (self.bcrparas or self.bctadds or self.bctparas
                      or self.bctsets or self.bsurf or self.bsurfs
                      or self.bconp or self.blseg or self.bfric
                      or self.bgadds or self.bgsets or self.bctparms or self.bcbodys or self.bcparas)
        if is_contact:
            bdf_file.write('$CONTACT\n')
            for (unused_id, bcbody) in sorted(self.bcbodys.items()):
                bdf_file.write(bcbody.write_card(size, is_double))
            for (unused_id, bcpara) in sorted(self.bcparas.items()):
                bdf_file.write(bcpara.write_card(size, is_double))

            for (unused_id, bcrpara) in sorted(self.bcrparas.items()):
                bdf_file.write(bcrpara.write_card(size, is_double))
            for (unused_id, bctparam) in sorted(self.bctparms.items()):
                bdf_file.write(bctparam.write_card(size, is_double))
            for (unused_id, bctadds) in sorted(self.bctadds.items()):
                bdf_file.write(bctadds.write_card(size, is_double))
            for (unused_id, bctpara) in sorted(self.bctparas.items()):
                bdf_file.write(bctpara.write_card(size, is_double))

            for (unused_id, bctset) in sorted(self.bctsets.items()):
                bdf_file.write(bctset.write_card(size, is_double))
            for (unused_id, bsurfi) in sorted(self.bsurf.items()):
                bdf_file.write(bsurfi.write_card(size, is_double))
            for (unused_id, bsurfsi) in sorted(self.bsurfs.items()):
                bdf_file.write(bsurfsi.write_card(size, is_double))
            for (unused_id, bconp) in sorted(self.bconp.items()):
                bdf_file.write(bconp.write_card(size, is_double))
            for (unused_id, blseg) in sorted(self.blseg.items()):
                bdf_file.write(blseg.write_card(size, is_double))
            for (unused_id, bfric) in sorted(self.bfric.items()):
                bdf_file.write(bfric.write_card(size, is_double))
            for (unused_id, bgadd) in sorted(self.bgadds.items()):
                bdf_file.write(bgadd.write_card(size, is_double))
            for (unused_id, bgset) in sorted(self.bgsets.items()):
                bdf_file.write(bgset.write_card(size, is_double))

    def _write_coords(self, bdf_file: Any, size: int=8, is_double: bool=False,
                      is_long_ids: Optional[bool]=None) -> None:
        """Writes the coordinate cards in a sorted order"""
        size, is_long_ids = self._write_mesh_long_ids_size(size, is_long_ids)
        if len(self.coords) > 1:
            bdf_file.write('$COORDS\n')
        for (coord_id, coord) in sorted(self.coords.items()):
            if coord_id == 0:
                continue
            try:
                bdf_file.write(coord.write_card(size, is_double))
            except RuntimeError:
                bdf_file.write(coord.write_card_16(is_double))

    def _write_matcids(self, bdf_file: Any, size: int=8, is_double: bool=False,
                      is_long_ids: Optional[bool]=None) -> None:
        """Writes the MATCID cards in a sorted order"""
        size, is_long_ids = self._write_mesh_long_ids_size(size, is_long_ids)

        if len(self.MATCID):
            bdf_file.write('$MATCID\n')
        for (cid, matcids) in sorted(self.MATCID.items()):
            for matcid in matcids:
                try:
                    bdf_file.write(matcid.write_card(size, is_double))
                except RuntimeError:
                    bdf_file.write(matcid.write_card_16(is_double))

    def _write_dmigs(self, bdf_file: Any, size: int=8, is_double: bool=False,
                     is_long_ids: Optional[bool]=None) -> None:
        """
        Writes the DMIG cards

        Parameters
        ----------
        size : int
            large field (16) or small field (8)

        """
        for (unused_name, dmig) in natsorted(self.dmig.items()):
            bdf_file.write(dmig.write_card(size, is_double))
        for (unused_name, dmi) in natsorted(self.dmi.items()):
            bdf_file.write(dmi.write_card(size, is_double))
        for (unused_name, dmij) in natsorted(self.dmij.items()):
            bdf_file.write(dmij.write_card(size, is_double))
        for (unused_name, dmiji) in natsorted(self.dmiji.items()):
            bdf_file.write(dmiji.write_card(size, is_double))
        for (unused_name, dmik) in natsorted(self.dmik.items()):
            bdf_file.write(dmik.write_card(size, is_double))
        for (unused_name, dmiax) in natsorted(self.dmiax.items()):
            bdf_file.write(dmiax.write_card(size, is_double))

    def _write_dynamic(self, bdf_file: Any, size: int=8, is_double: bool=False,
                       is_long_ids: Optional[bool]=None) -> None:
        """Writes the dynamic cards sorted by ID"""
        is_dynamic = (self.dareas or self.dphases or self.nlparms or self.frequencies or
                      self.methods or self.cMethods or self.tsteps or self.tstepnls or
                      self.transfer_functions or self.delays or self.rotors or self.tics or
                      self.nlpcis)
        if is_dynamic:
            bdf_file.write('$DYNAMIC\n')
            for (unused_id, method) in sorted(self.methods.items()):
                bdf_file.write(method.write_card(size, is_double))
            for (unused_id, cmethod) in sorted(self.cMethods.items()):
                bdf_file.write(cmethod.write_card(size, is_double))
            for (unused_id, darea) in sorted(self.dareas.items()):
                bdf_file.write(darea.write_card(size, is_double))
            for (unused_id, dphase) in sorted(self.dphases.items()):
                bdf_file.write(dphase.write_card(size, is_double))
            for (unused_id, nlparm) in sorted(self.nlparms.items()):
                bdf_file.write(nlparm.write_card(size, is_double))
            for (unused_id, nlpci) in sorted(self.nlpcis.items()):
                bdf_file.write(nlpci.write_card(size, is_double))
            for (unused_id, tstep) in sorted(self.tsteps.items()):
                bdf_file.write(tstep.write_card(size, is_double))
            for (unused_id, tstepnl) in sorted(self.tstepnls.items()):
                bdf_file.write(tstepnl.write_card(size, is_double))
            for (unused_id, freqs) in sorted(self.frequencies.items()):
                for freq in freqs:
                    bdf_file.write(freq.write_card(size, is_double))
            for (unused_id, delay) in sorted(self.delays.items()):
                bdf_file.write(delay.write_card(size, is_double))
            for (unused_id, rotor) in sorted(self.rotors.items()):
                bdf_file.write(rotor.write_card(size, is_double))
            for (unused_id, tic) in sorted(self.tics.items()):
                bdf_file.write(tic.write_card(size, is_double))

            for (unused_id, tfs) in sorted(self.transfer_functions.items()):
                for transfer_function in tfs:
                    bdf_file.write(transfer_function.write_card(size, is_double))

    def _write_mesh_long_ids_size(self, size: bool, is_long_ids: bool) -> tuple[int, bool]:
        """helper method"""
        if is_long_ids and size == 16 or is_long_ids is False:
            return size, is_long_ids

        if size == 16 and is_long_ids is None or self.is_long_ids:
            size = 16
            is_long_ids = True
        else:
            is_long_ids = False
        return size, is_long_ids

    def _write_loads(self, bdf_file: Any, size: int=8, is_double: bool=False,
                     is_long_ids: Optional[bool]=None) -> None:
        """Writes the load cards sorted by ID"""
        size, is_long_ids = self._write_mesh_long_ids_size(size, is_long_ids)
        if self.load_combinations or self.loads or self.tempds or self.cyjoin:
            bdf_file.write('$LOADS\n')
            for (key, load_combinations) in sorted(self.load_combinations.items()):
                for load_combination in load_combinations:
                    try:
                        bdf_file.write(load_combination.write_card(size, is_double))
                    except Exception:
                        print(f'failed printing load...type={load_combination.type} key={key!r}')
                        raise

            if is_long_ids:
                for (key, loadcase) in sorted(self.loads.items()):
                    for load in loadcase:
                        try:
                            bdf_file.write(load.write_card_16(is_double))
                        except Exception:
                            print(f'failed printing load...type={load.type} key={key!r}')
                            raise
            else:
                for (key, loadcase) in sorted(self.loads.items()):
                    for load in loadcase:
#                        if load.type == 'PLOAD2':
#                            try:
#                                load.raw_fields()
#                            except Exception:
#                                bdf_file.write(load.write_card_separate(self, size, is_double))
#                                continue
                        try:
                            bdf_file.write(load.write_card(size, is_double))
                        except Exception:
                            print(f'failed printing load...type={load.type} key={key!r}')
                            raise

            for unused_key, tempd in sorted(self.tempds.items()):
                bdf_file.write(tempd.write_card(size, is_double))
            for unused_key, cyjoin in sorted(self.cyjoin.items()):
                bdf_file.write(cyjoin.write_card(size, is_double))
        self._write_dloads(bdf_file, size=size, is_double=is_double, is_long_ids=is_long_ids)

    def _write_dloads(self, bdf_file: Any, size: int=8, is_double: bool=False,
                      is_long_ids: Optional[bool]=None) -> None:
        """Writes the dload cards sorted by ID"""
        size, is_long_ids = self._write_mesh_long_ids_size(size, is_long_ids)
        if self.dloads or self.dload_entries:
            bdf_file.write('$DLOADS\n')
            for (key, loadcase) in sorted(self.dloads.items()):
                for load in loadcase:
                    try:
                        bdf_file.write(load.write_card(size, is_double))
                    except Exception:
                        print(f'failed printing load...type={load.type} key={key!r}')
                        raise

            for (key, loadcase) in sorted(self.dload_entries.items()):
                for load in loadcase:
                    try:
                        bdf_file.write(load.write_card(size, is_double))
                    except Exception:
                        print(f'failed printing load...type={load.type} key={key!r}')
                        raise


    def _write_masses(self, bdf_file: Any, size: int=8, is_double: bool=False,
                      is_long_ids: Optional[bool]=None) -> None:
        """Writes the mass cards sorted by ID"""
        size, is_long_ids = self._write_mesh_long_ids_size(size, is_long_ids)
        if self.properties_mass:
            bdf_file.write('$PROPERTIES_MASS\n')
            for (pid, mass) in sorted(self.properties_mass.items()):
                try:
                    bdf_file.write(mass.write_card(size, is_double))
                except Exception:
                    print(f'failed printing mass property...type={mass.type} pid={pid}')
                    raise

        if self.masses:
            bdf_file.write('$MASSES\n')
            for (eid, mass) in sorted(self.masses.items()):
                try:
                    bdf_file.write(mass.write_card(size, is_double))
                except Exception:
                    print(f'failed printing masses...type={mass.type} eid={eid}')
                    raise

    def _write_materials(self, bdf_file: Any, size: int=8, is_double: bool=False,
                         is_long_ids: Optional[bool]=None) -> None:
        """Writes the materials in a sorted order"""
        size, is_long_ids = self._write_mesh_long_ids_size(size, is_long_ids)
        is_big_materials = hasattr(self, 'big_materials') and self.big_materials
        is_materials = (self.materials or self.hyperelastic_materials or self.creep_materials or
                        self.MATS1 or self.MATS3 or self.MATS8 or self.MATT1 or
                        self.MATT2 or self.MATT3 or self.MATT4 or self.MATT5 or
                        self.MATT8 or self.MATT9 or self.nxstrats or is_big_materials)
        if is_materials:
            bdf_file.write('$MATERIALS\n')
            for (unused_mid, material) in sorted(self.materials.items()):
                bdf_file.write(material.write_card(size, is_double))
            for (unused_mid, material) in sorted(self.hyperelastic_materials.items()):
                bdf_file.write(material.write_card(size, is_double))
            for (unused_mid, material) in sorted(self.creep_materials.items()):
                bdf_file.write(material.write_card(size, is_double))

            for (unused_mid, material) in sorted(self.MATS1.items()):
                bdf_file.write(material.write_card(size, is_double))
            for (unused_mid, material) in sorted(self.MATS3.items()):
                bdf_file.write(material.write_card(size, is_double))
            for (unused_mid, material) in sorted(self.MATS8.items()):
                bdf_file.write(material.write_card(size, is_double))

            for (unused_mid, material) in sorted(self.MATT1.items()):
                bdf_file.write(material.write_card(size, is_double))
            for (unused_mid, material) in sorted(self.MATT2.items()):
                bdf_file.write(material.write_card(size, is_double))
            for (unused_mid, material) in sorted(self.MATT3.items()):
                bdf_file.write(material.write_card(size, is_double))
            for (unused_mid, material) in sorted(self.MATT4.items()):
                bdf_file.write(material.write_card(size, is_double))
            for (unused_mid, material) in sorted(self.MATT5.items()):
                bdf_file.write(material.write_card(size, is_double))
            for (unused_mid, material) in sorted(self.MATT8.items()):
                bdf_file.write(material.write_card(size, is_double))
            for (unused_mid, material) in sorted(self.MATT9.items()):
                bdf_file.write(material.write_card(size, is_double))
            for (unused_mid, material) in sorted(self.MATDMG.items()):
                bdf_file.write(material.write_card(size, is_double))
            for (unused_sid, nxstrat) in sorted(self.nxstrats.items()):
                bdf_file.write(nxstrat.write_card(size, is_double))

            if is_big_materials:
                for unused_mid, mat in sorted(self.big_materials.items()):
                    bdf_file.write(mat.write_card_16(is_double))

    def _write_nodes(self, bdf_file: Any, size: int=8, is_double: bool=False,
                     is_long_ids: Optional[bool]=None) -> None:
        """Writes the NODE-type cards"""
        if self.spoints:
            bdf_file.write('$SPOINTS\n')
            bdf_file.write(write_xpoints('SPOINT', self.spoints))
        if self.epoints:
            bdf_file.write('$EPOINTS\n')
            bdf_file.write(write_xpoints('EPOINT', self.epoints))
        if self.points:
            bdf_file.write('$POINTS\n')
            for unused_point_id, point in sorted(self.points.items()):
                bdf_file.write(point.write_card(size, is_double))

        if self._is_axis_symmetric:
            if self.axic:
                bdf_file.write(self.axic.write_card(size, is_double))
            if self.axif:
                bdf_file.write(self.axif.write_card(size, is_double))
            for unused_nid, ringax_pointax in sorted(self.ringaxs.items()):
                bdf_file.write(ringax_pointax.write_card(size, is_double))
            for unused_ringfl, ringfl in sorted(self.ringfl.items()):
                bdf_file.write(ringfl.write_card(size, is_double))
            for unused_nid, gridb in sorted(self.gridb.items()):
                bdf_file.write(gridb.write_card(size, is_double))
        if self.cyax:
            bdf_file.write(self.cyax.write_card(size, is_double))

        self._write_grids(bdf_file, size=size, is_double=is_double)
        if self.seqgp:
            bdf_file.write(self.seqgp.write_card(size, is_double))

        #if 0:  # not finished
            #self._write_nodes_associated(bdf_file, size, is_double)

    def _write_grids(self, bdf_file: Any, size: int=8, is_double: bool=False,
                     is_long_ids: Optional[bool]=None) -> None:
        """Writes the GRID-type cards"""
        size, is_long_ids = self._write_mesh_long_ids_size(size, is_long_ids)
        if self.nodes:
            bdf_file.write('$NODES\n')
            if self.grdset:
                bdf_file.write(self.grdset.write_card(size))
            write_dict(bdf_file, self.nodes, size, is_double, is_long_ids)

    #def _write_nodes_associated(self, bdf_file, size=8, is_double=False):
        #"""
        #Writes the NODE-type in associated and unassociated groups.

        #.. warning:: Sometimes crashes, probably on invalid BDFs.
        #"""
        #associated_nodes = set()
        #for (eid, element) in self.elements.items():
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
            #for key, node in sorted(associated_nodes.items()):
                #bdf_file.write(node.write_card(size, is_double))

        #if unassociated_nodes:
            #bdf_file.write('$UNASSOCIATED NODES\n')
            #if self.grdset and not associated_nodes:
                #v(self.grdset.write_card(size, is_double))
            #for key, node in sorted(unassociated_nodes.items()):
                #if key in self.nodes:
                    #bdf_file.write(node.write_card(size, is_double))
                #else:
                    #bdf_file.write('$ Missing NodeID=%s' % key)

    def _write_optimization(self, bdf_file: Any, size: int=8, is_double: bool=False,
                            is_long_ids: Optional[bool]=None) -> None:
        """Writes the optimization cards sorted by ID"""
        is_optimization = (
            self.dconadds or self.dconstrs or self.desvars or self.ddvals or
            self.dresps or
            self.dvprels or self.dvmrels or self.dvcrels or self.doptprm or
            self.dlinks or self.dequations or self.dtable is not None or
            self.dvgrids or self.dscreen or self.topvar or self.modtrak or
            # nx optimization
            self.dvtrels or self.group or self.dmncon
        )
        if is_optimization:
            bdf_file.write('$OPTIMIZATION\n')
            for (unused_id, dconadd) in sorted(self.dconadds.items()):
                bdf_file.write(dconadd.write_card(size, is_double))
            for (unused_id, dconstrs) in sorted(self.dconstrs.items()):
                for dconstr in dconstrs:
                    bdf_file.write(dconstr.write_card(size, is_double))
            for (unused_id, desvar) in sorted(self.desvars.items()):
                bdf_file.write(desvar.write_card(size, is_double))
            for (unused_id, topvar) in sorted(self.topvar.items()):
                bdf_file.write(topvar.write_card(size, is_double))
            for (unused_id, ddval) in sorted(self.ddvals.items()):
                bdf_file.write(ddval.write_card(size, is_double))
            for (unused_id, dlink) in sorted(self.dlinks.items()):
                bdf_file.write(dlink.write_card(size, is_double))
            for (unused_id, dresp) in sorted(self.dresps.items()):
                bdf_file.write(dresp.write_card(size, is_double))

            for (unused_id, dvcrel) in sorted(self.dvcrels.items()):
                bdf_file.write(dvcrel.write_card(size, is_double))
            for (unused_id, dvmrel) in sorted(self.dvmrels.items()):
                bdf_file.write(dvmrel.write_card(size, is_double))
            for (unused_id, dvprel) in sorted(self.dvprels.items()):
                bdf_file.write(dvprel.write_card(size, is_double))
            for (unused_id, dvgrids) in sorted(self.dvgrids.items()):
                for dvgrid in dvgrids:
                    bdf_file.write(dvgrid.write_card(size, is_double))
            for (unused_id, dscreen) in sorted(self.dscreen.items()):
                bdf_file.write(str(dscreen))

            for (unused_id, equation) in sorted(self.dequations.items()):
                bdf_file.write(str(equation))

            if self.dtable is not None:
                bdf_file.write(self.dtable.write_card(size, is_double))
            if self.doptprm is not None:
                bdf_file.write(self.doptprm.write_card(size, is_double))
            if self.modtrak is not None:
                bdf_file.write(self.modtrak.write_card(size, is_double))

            # nx optimization
            for (unused_id, dvtrel) in sorted(self.dvtrels.items()):
                bdf_file.write(dvtrel.write_card(size, is_double))
            for (unused_id, group) in sorted(self.group.items()):
                bdf_file.write(group.write_card(size, is_double))
            for (unused_id, dmncon) in sorted(self.dmncon.items()):
                bdf_file.write(dmncon.write_card(size, is_double))

    def _write_parametric(self, bdf_file: Any, size: int=8, is_double: bool=False,
                          is_long_ids: Optional[bool]=None) -> None:
        """Writes the optimization cards sorted by ID"""
        is_parametric = self.pset or self.pval or self.gmcurv or self.feedge or self.feface
        if is_parametric:
            for (unused_id, pset) in sorted(self.pset.items()):
                bdf_file.write(pset.write_card(size, is_double))

            for (unused_adapt_id, pvals) in sorted(self.pval.items()):
                for pval in pvals:
                    bdf_file.write(pval.write_card(size, is_double))

            for (unused_id, gmcurv) in sorted(self.gmcurv.items()):
                bdf_file.write(gmcurv.write_card(size, is_double))

            for (unused_id, gmsurf) in sorted(self.gmsurf.items()):
                bdf_file.write(gmsurf.write_card(size, is_double))

            for (unused_id, feedge) in sorted(self.feedge.items()):
                bdf_file.write(feedge.write_card(size, is_double))

            for (unused_id, feface) in sorted(self.feface.items()):
                bdf_file.write(feface.write_card(size, is_double))

    def _write_params(self, bdf_file: Any, size: int=8, is_double: bool=False,
                      is_long_ids: Optional[bool]=None) -> None:
        """Writes the PARAM cards"""
        size, is_long_ids = self._write_mesh_long_ids_size(size, is_long_ids)
        if self.params or self.dti or self.mdlprm:
            bdf_file.write('$PARAMS\n')
            for unused_name, dti in sorted(self.dti.items()):
                bdf_file.write(dti.write_card(size=size, is_double=is_double))

            for (unused_key, param) in sorted(self.params.items()):
                bdf_file.write(param.write_card(size, is_double))
            if self.mdlprm:
                bdf_file.write(self.mdlprm.write_card(size, is_double))

    def _write_properties(self, bdf_file: Any, size: int=8, is_double: bool=False,
                          is_long_ids: Optional[bool]=None) -> None:
        """Writes the properties in a sorted order"""
        size, is_long_ids = self._write_mesh_long_ids_size(size, is_long_ids)
        is_big_properties = hasattr(self, 'big_properties') and self.big_properties
        is_properties = (self.properties or self.pelast or
                         self.pdampt or self.pbusht or is_big_properties)
        if is_properties:
            bdf_file.write('$PROPERTIES\n')
            prop_groups = (self.properties, self.pelast, self.pdampt, self.pbusht)
            if is_long_ids:
                for prop_group in prop_groups:
                    for unused_pid, prop in sorted(prop_group.items()):
                        bdf_file.write(prop.write_card_16(is_double))
                #except Exception:
                    #print('failed printing property type=%s' % prop.type)
                    #raise
            else:
                for prop_group in prop_groups:
                    for unused_pid, prop in sorted(prop_group.items()):
                        bdf_file.write(prop.write_card(size, is_double))

            if is_big_properties:
                for unused_pid, prop in sorted(self.big_properties.items()):
                    bdf_file.write(prop.write_card_16(is_double))

    def _write_properties_by_element_type(self, bdf_file: Any, size: int=8, is_double: bool=False,
                                          is_long_ids: Optional[bool]=None) -> None:
        """
        Writes the properties in a sorted order by property type grouping

        TODO: Missing some property types.
        """
        size, is_long_ids = self._write_mesh_long_ids_size(size, is_long_ids)
        is_properties = self.properties or self.pelast or self.pdampt or self.pbusht
        if not is_properties:
            return

        out = get_properties_by_element_type(self)
        propertys_class_to_property_types, property_type_to_property_class, properties_by_class = out

        bdf_file.write('$PROPERTIES\n')
        for prop_class, prop_types in propertys_class_to_property_types.items():
            #print(prop_class, prop_types)
            #for prop_type in prop_types:
                #if prop_type not in properties_by_class:
                    #continue
                #print('  ', prop_type)
            props = properties_by_class[prop_class]
            if not props:
                continue
            bdf_file.write('$' + '-' * 80 + '\n')
            bdf_file.write('$ %s\n' % prop_class)

            for prop in props:
                bdf_file.write(prop.write_card(size, is_double))
        bdf_file.write('$' + '-' * 80 + '\n')

    def _write_rejects(self, bdf_file: Any, size: int=8, is_double: bool=False,
                       is_long_ids: Optional[bool]=None) -> None:
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
                    if len(reject_card) > 0:
                        line0 = reject_card[0].upper()
                        if line0.startswith('ADAPT'):
                            for line in reject_card:
                                assert isinstance(line, str), line
                                bdf_file.write(line+'\n')
                            continue
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
                elif isinstance(reject_lines, str):
                    reject2 = reject_lines.rstrip()
                    if reject2:
                        bdf_file.write('%s\n' % reject2)
                else:
                    raise TypeError(reject_lines)

    def _write_rigid_elements(self, bdf_file: Any, size: int=8, is_double: bool=False,
                              is_long_ids: Optional[bool]=None) -> None:
        """Writes the rigid elements in a sorted order"""
        size, is_long_ids = self._write_mesh_long_ids_size(size, is_long_ids)
        if self.rigid_elements:
            bdf_file.write('$RIGID ELEMENTS\n')
            if is_long_ids:
                for (eid, element) in sorted(self.rigid_elements.items()):
                    try:
                        bdf_file.write(element.write_card_16(is_double))
                    except Exception:
                        print(f'failed printing element...type={element.type} eid={eid}')
                        raise
            else:
                for (eid, element) in sorted(self.rigid_elements.items()):
                    try:
                        bdf_file.write(element.write_card(size, is_double))
                    except Exception:
                        print(f'failed printing element...type={element.type} eid={eid}')
                        raise
        if self.plotels:
            bdf_file.write('$PLOT ELEMENTS\n')
            write_dict(bdf_file, self.plotels, size, is_double, is_long_ids)

    def _write_sets(self, bdf_file: Any, size: int=8, is_double: bool=False,
                    is_long_ids: Optional[bool]=None) -> None:
        """Writes the SETx cards sorted by ID"""
        is_sets = (self.sets or self.asets or self.omits or self.bsets or self.csets or self.qsets
                   or self.usets)
        if is_sets:
            bdf_file.write('$SETS\n')
            for (unused_id, set_obj) in sorted(self.sets.items()):  # dict
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
            for unused_name, usets in sorted(self.usets.items()):  # dict
                for set_obj in usets:  # list
                    bdf_file.write(set_obj.write_card(size, is_double))

    def _write_superelements(self, bdf_file: Any, size: int=8, is_double: bool=False,
                             is_long_ids: Optional[bool]=None) -> None:
        """
        Writes the Superelement cards

        Parameters
        ----------
        size : int
            large field (16) or small field (8)

        """
        is_sets = (self.se_sets or self.se_bsets or self.se_csets or self.se_qsets
                   or self.se_usets)
        if is_sets:
            bdf_file.write('$SUPERELEMENTS\n')
            for set_obj in self.se_bsets:  # list
                bdf_file.write(set_obj.write_card(size, is_double))
            for set_obj in self.se_csets:  # list
                bdf_file.write(set_obj.write_card(size, is_double))
            for set_obj in self.se_qsets:  # list
                bdf_file.write(set_obj.write_card(size, is_double))
            for (unused_set_id, set_obj) in sorted(self.se_sets.items()):  # dict
                bdf_file.write(set_obj.write_card(size, is_double))
            for unused_name, usets in sorted(self.se_usets.items()):  # dict
                for set_obj in usets:  # list
                    bdf_file.write(set_obj.write_card(size, is_double))
            for suport in self.se_suport:  # list
                bdf_file.write(suport.write_card(size, is_double))

        for unused_seid, csuper in sorted(self.csuper.items()):
            bdf_file.write(csuper.write_card(size, is_double))
        for unused_seid, csupext in sorted(self.csupext.items()):
            bdf_file.write(csupext.write_card(size, is_double))

        for unused_seid, sebulk in sorted(self.sebulk.items()):
            bdf_file.write(sebulk.write_card(size, is_double))
        for unused_seid, seconct in sorted(self.seconct.items()):
            bdf_file.write(seconct.write_card(size, is_double))

        for unused_seid, sebndry in sorted(self.sebndry.items()):
            bdf_file.write(sebndry.write_card(size, is_double))
        for unused_seid, seelt in sorted(self.seelt.items()):
            bdf_file.write(seelt.write_card(size, is_double))
        for unused_seid, seexcld in sorted(self.seexcld.items()):
            bdf_file.write(seexcld.write_card(size, is_double))

        for unused_seid, selabel in sorted(self.selabel.items()):
            bdf_file.write(selabel.write_card(size, is_double))
        for unused_seid, seloc in sorted(self.seloc.items()):
            bdf_file.write(seloc.write_card(size, is_double))
        for unused_seid, seload in sorted(self.seload.items()):
            bdf_file.write(seload.write_card(size, is_double))
        for unused_seid, sempln in sorted(self.sempln.items()):
            bdf_file.write(sempln.write_card(size, is_double))
        for unused_setid, senqset in sorted(self.senqset.items()):
            bdf_file.write(senqset.write_card(size, is_double))
        for unused_seid, setree in sorted(self.setree.items()):
            bdf_file.write(setree.write_card(size, is_double))
        for unused_seid, release in sorted(self.release.items()):
            bdf_file.write(release.write_card(size, is_double))


    def _write_tables(self, bdf_file: Any, size: int=8, is_double: bool=False,
                      is_long_ids: Optional[bool]=None) -> None:
        """Writes the TABLEx cards sorted by ID"""
        if self.tables or self.tables_d or self.tables_m or self.tables_sdamping:
            bdf_file.write('$TABLES\n')
            for (unused_id, table) in sorted(self.tables.items()):
                bdf_file.write(table.write_card(size, is_double))
            for (unused_id, table) in sorted(self.tables_d.items()):
                bdf_file.write(table.write_card(size, is_double))
            for (unused_id, table) in sorted(self.tables_m.items()):
                bdf_file.write(table.write_card(size, is_double))
            for (unused_id, table) in sorted(self.tables_sdamping.items()):
                bdf_file.write(table.write_card(size, is_double))

        if self.random_tables:
            bdf_file.write('$RANDOM TABLES\n')
            for (unused_id, table) in sorted(self.random_tables.items()):
                bdf_file.write(table.write_card(size, is_double))

    def _write_thermal(self, bdf_file: Any, size: int=8, is_double: bool=False,
                       is_long_ids: Optional[bool]=None) -> None:
        """Writes the thermal cards"""
        # PHBDY
        is_thermal = (self.phbdys or self.convection_properties or self.bcs or
                      self.views or self.view3ds or self.radset or self.radcavs)
        if is_thermal:
            bdf_file.write('$THERMAL\n')

            for (unused_key, phbdy) in sorted(self.phbdys.items()):
                bdf_file.write(phbdy.write_card(size, is_double))

            #for unused_key, prop in sorted(self.thermal_properties.items()):
            #    bdf_file.write(str(prop))
            for (unused_key, prop) in sorted(self.convection_properties.items()):
                bdf_file.write(prop.write_card(size, is_double))

            # BCs
            for (unused_key, bcs) in sorted(self.bcs.items()):
                for boundary_condition in bcs:  # list
                    bdf_file.write(boundary_condition.write_card(size, is_double))

            for (unused_key, view) in sorted(self.views.items()):
                bdf_file.write(view.write_card(size, is_double))
            for (unused_key, view3d) in sorted(self.view3ds.items()):
                bdf_file.write(view3d.write_card(size, is_double))
            if self.radset:
                bdf_file.write(self.radset.write_card(size, is_double))
            for unused_icavity, radcav in self.radcavs.items():
                bdf_file.write(radcav.write_card(size, is_double))


    def _write_thermal_materials(self, bdf_file: Any, size: int=8, is_double: bool=False,
                                 is_long_ids: Optional[bool]=None) -> None:
        """Writes the thermal materials in a sorted order"""
        if self.thermal_materials:
            bdf_file.write('$THERMAL MATERIALS\n')
            for (unused_mid, material) in sorted(self.thermal_materials.items()):
                bdf_file.write(material.write_card(size, is_double))

def _fix_sizes(size: int,
               nodes_size: Optional[int],
               elements_size: Optional[int],
               loads_size: Optional[int]) -> tuple[int, int]:
    if nodes_size is None:
        nodes_size = size
    if elements_size is None:
        elements_size = size
    if loads_size is None:
        loads_size = size
    return size, nodes_size, elements_size, loads_size

def _output_helper(out_filename: Optional[str], interspersed: bool,
                   size: int, is_double: bool, log: SimpleLogger) -> tuple[str, int]:
    """Performs type checking on the write_bdf inputs"""
    if out_filename is None:
        from pyNastran.utils.gui_io import save_file_dialog
        wildcard_wx = "Nastran BDF (*.bdf; *.dat; *.nas; *.pch)|" \
            "*.bdf;*.dat;*.nas;*.pch|" \
            "All files (*.*)|*.*"
        wildcard_qt = "Nastran BDF (*.bdf *.dat *.nas *.pch);;All files (*)"
        title = 'Save BDF/DAT/PCH'
        out_filename = save_file_dialog(title, wildcard_wx, wildcard_qt)
        assert out_filename is not None, out_filename

    has_read_write = hasattr(out_filename, 'read') and hasattr(out_filename, 'write')
    if has_read_write or isinstance(out_filename, IOBase):
        return out_filename, size
    if not isinstance(out_filename, (str, PurePath)):
        msg = f'out_filename={out_filename!r} must be a string; type={type(out_filename)}'
        raise TypeError(msg)

    assert size in {8, 16}, f'size={size!r}'
    assert is_double in {True, False}, f'is_double={is_double!r}'
    if size == 8:
        if is_double is True:
            log.warning('is_double=True...changing size from 8 to 16...')
            size = 16
    else:
        assert is_double in {True, False}, f'is_double={is_double!r}'

    assert isinstance(interspersed, bool)
    #fname = print_filename(out_filename)
    #self.log.debug("***writing %s" % fname)
    return out_filename, size

def get_optimization_include(model: BDF) -> tuple[list[int], list[int]]:
    """gets the properties and materials refereced by DVPRELx/DVMRELx"""
    property_types = {'PELAS', 'PDAMP', 'PGAP', 'PBUSH', 'PBUSH1D', 'PVISC', 'PWELD',
                  'PROD', 'PTUBE', 'PBAR', 'PBARL', 'PBEAM', 'PBEAML', 'PBMSECT',
                  'PSHEAR', 'PSHELL', 'PCOMP', 'PCOMPG', }
    material_types = {'MAT1', 'MAT8', 'MAT9'}

    pids_to_remove = set([])
    for dvprel_id, dvprel in model.dvprels.items():
        if dvprel.prop_type in property_types:
            pid = dvprel.pid
            pids_to_remove.add(pid)

    mids_to_remove = set([])
    for dvmrel_id, dvmrel in model.dvmrels.items():
        if dvmrel.mat_type in material_types:
            mid = dvprel.mid
            mids_to_remove.add(mid)
    return list(pids_to_remove), list(mids_to_remove)


def delete_optimization_data(model: BDF) -> tuple[list[int], list[int]]:
    """removes optimization referenced data (that will be in the PCH file)"""
    pids_to_remove, mids_to_remove = get_optimization_include(model)
    properties_to_write, materials_to_write, desvars = _delete_optimization_data(
        model, pids_to_remove, mids_to_remove)
    return properties_to_write, materials_to_write, desvars

def _delete_optimization_data(model: BDF,
                              pids_to_remove: list[int],
                              mids_to_remove: list[int]) -> tuple[list[Any], list[Any], list[DESVAR]]:
    """heper method"""
    properties_to_write = []
    materials_to_write = []

    #--------------------------------------------------------------
    prop_types = set([])
    for pid in sorted(pids_to_remove):
        try:
            prop = model.properties[pid]
        except KeyError:
            continue
        prop_types.add(prop.type)

    pids_to_delete = []
    for pid, prop in model.properties.items():
        if prop.type in prop_types:
            properties_to_write.append(prop)
            pids_to_delete.append(pid)
    for pid in pids_to_delete:
        del model.properties[pid]

    # ----------------------------------------------------------------
    #mids_to_remove = list(mids_to_remove)
    mat_types = set([])
    for mid in sorted(mids_to_remove):
        try:
            mat = model.materials[mid]
        except KeyError:
            continue
        mat_types.add(mat.type)

    mids_to_delete = []
    for mid, mat in model.materials.items():
        if mat.type in mat_types:
            materials_to_write.append(mat)
            del model.materials[pid]
            mids_to_delete.append(mid)
    for mid in mids_to_delete:
        del model.mids_to_delete[mid]

    # ----------------------------------------------------------------
    desvars = []
    for desvar_id, desvar in sorted(model.desvars.items()):
        desvars.append(desvar)
    model.desvars = {}
    return properties_to_write, materials_to_write, desvars


def write_optimization_include(model: BDF, pch_include_filename: str, size: int=8) -> None:
    """writes an optimization deck"""
    properties_to_write, materials_to_write, desvars = delete_optimization_data(
        model)
    encoding = model.get_encoding()
    with open(pch_include_filename, 'w', encoding=encoding) as pch_file:
        for prop in properties_to_write:
            pch_file.write(prop.write_card(size=size))
        for mat in materials_to_write:
            pch_file.write(mat.write_card(size=size))

        for desvar in desvars:
            pch_file.write(desvar.write_card(size=size))
    return

