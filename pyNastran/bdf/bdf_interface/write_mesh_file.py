# coding: utf-8
"""
This file defines:
  - WriteMesh
"""
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
#import sys
#import io
from typing import List, Dict, Union, Optional, Tuple, Any, cast
from codecs import open
from six import string_types #, iteritems, PY2, StringIO

#from pyNastran.bdf.bdf_interface.utils import print_filename
from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.field_writer_16 import print_card_16
from pyNastran.bdf.bdf_interface.write_mesh import WriteMesh
#from pyNastran.bdf.cards.nodes import write_xpoints


class WriteMeshs(WriteMesh):
    """
    Defines methods for writing cards

    Major methods:
      - model.write_bdf(...)
      - model.echo_bdf(...)
      - model.auto_reject_bdf(...)
    """
    def __init__(self):
        """creates methods for writing cards"""
        WriteMesh.__init__(self)

    def write_bdfs(self, out_filenames, encoding=None,
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

        if self.is_bdf_vectorized:  # pragma: no cover
            raise NotImplementedError()
        else:
            # required for MasterModelTaxi
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
        out_filename = self._output_helper(out_filenames[0],
                                           interspersed, size, is_double)
        self.log.debug('---starting BDF.write_bdf of %s---' % out_filename)
        encoding = self.get_encoding(encoding)

        class DevNull(object):
            def write(self, *_):
                pass

        devnull = DevNull()
        bdf_files = {i : devnull for i in range(len(self.active_filenames))}
        for ifile, out_filename in out_filenames.items():
            if hasattr(out_filename, 'read') and hasattr(out_filename, 'write'):
                bdf_file = out_filename
            else:
                bdf_file = open(out_filename, 'w', encoding=encoding)
            bdf_files[ifile] = bdf_file

        self._write_header(bdf_files[0], encoding)
        self._write_params_file(bdf_files, size, is_double, is_long_ids=is_long_ids)
        self._write_nodes_file(bdf_files, size, is_double, is_long_ids=is_long_ids)

        self._write_elements_file(bdf_files, size, is_double, is_long_ids=is_long_ids)
        self._write_properties_file(bdf_files, size, is_double, is_long_ids=is_long_ids)
        self._write_materials_file(bdf_files, size, is_double, is_long_ids=is_long_ids)

        self._write_masses_file(bdf_files, size, is_double, is_long_ids=is_long_ids)
        self._write_rigid_elements_file(bdf_files, size, is_double, is_long_ids=is_long_ids) # split out for write_bdf_symmetric
        self._write_aero_file(bdf_files, size, is_double, is_long_ids=is_long_ids)  # split out for write_bdf_symmetric

        self._write_common_file(bdf_files, size, is_double, is_long_ids=is_long_ids)
        if (enddata is None and 'ENDDATA' in self.card_count) or enddata:
            bdf_file.write('ENDDATA\n')
        if close:
            bdf_file.close()


    def _write_elements_file(self, bdf_files, size=8, is_double=False, is_long_ids=None):
        # type: (Any, int, bool) -> None
        """
        Writes the elements in a sorted order
        """
        size, is_long_ids = self._write_mesh_long_ids_size(size, is_long_ids)
        if self.elements:
            if is_long_ids:
                for (eid, element) in sorted(self.elements.items()):
                    bdf_files[element.ifile].write(element.write_card_16(is_double))
            else:
                for (eid, element) in sorted(self.elements.items()):
                    try:
                        bdf_files[element.ifile].write(element.write_card(size, is_double))
                    except:
                        print('failed printing element...'
                              'type=%s eid=%s' % (element.type, eid))
                        raise
        if self.ao_element_flags:
            for (eid, element) in sorted(self.ao_element_flags.items()):
                bdf_files[element.ifile].write(element.write_card(size, is_double))
        if self.normals:
            for (unused_nid, snorm) in sorted(self.normals.items()):
                bdf_files[snorm.ifile].write(snorm.write_card(size, is_double))
        self._write_nsm_file(bdf_files, size, is_double)

    def _write_nsm_file(self, bdf_files, size=8, is_double=False, is_long_ids=None):
        # type: (Any, int, bool) -> None
        """
        Writes the nsm in a sorted order
        """
        if self.nsms or self.nsmadds:
            for (unused_id, nsmadds) in sorted(self.nsmadds.items()):
                for nsmadd in nsmadds:
                    bdf_files[nsmadd.ifile].write(str(nsmadd))
            for (key, nsms) in sorted(self.nsms.items()):
                for nsm in nsms:
                    try:
                        bdf_files[nsm.ifile].write(nsm.write_card(size, is_double))
                    except:
                        print('failed printing nsm...type=%s key=%r'
                              % (nsm.type, key))
                        raise


    def _write_aero_file(self, bdf_files, size=8, is_double=False, is_long_ids=None):
        # type: (Any, int, bool) -> None
        """Writes the aero cards"""
        if self.caeros or self.paeros or self.monitor_points or self.splines:
            for (unused_id, caero) in sorted(self.caeros.items()):
                bdf_files[caero.ifile].write(caero.write_card(size, is_double))
            for (unused_id, paero) in sorted(self.paeros.items()):
                bdf_files[paero.ifile].write(paero.write_card(size, is_double))
            for (unused_id, spline) in sorted(self.splines.items()):
                bdf_files[spline.ifile].write(spline.write_card(size, is_double))
            for monitor_point in self.monitor_points:
                bdf_files[monitor_point.ifile].write(monitor_point.write_card(size, is_double))
        self.zona.write_bdf(bdf_files[0], size=8, is_double=False)

    def _write_aero_control_file(self, bdf_files, size=8, is_double=False, is_long_ids=None):
        # type: (Any, int, bool) -> None
        """Writes the aero control surface cards"""
        if(self.aecomps or self.aefacts or self.aeparams or self.aelinks or
           self.aelists or self.aestats or self.aesurf or self.aesurfs):
            for (unused_id, aelinks) in sorted(self.aelinks.items()):
                for aelink in aelinks:
                    bdf_files[aelink.ifile].write(aelink.write_card(size, is_double))

            for (unused_id, aecomp) in sorted(self.aecomps.items()):
                bdf_files[aecomp.ifile].write(aecomp.write_card(size, is_double))
            for (unused_id, aeparam) in sorted(self.aeparams.items()):
                bdf_files[aeparam.ifile].write(aeparam.write_card(size, is_double))
            for (unused_id, aestat) in sorted(self.aestats.items()):
                bdf_files[aestat.ifile].write(aestat.write_card(size, is_double))

            for (unused_id, aelist) in sorted(self.aelists.items()):
                bdf_files[aelist.ifile].write(aelist.write_card(size, is_double))
            for (unused_id, aesurf) in sorted(self.aesurf.items()):
                bdf_files[aesurf.ifile].write(aesurf.write_card(size, is_double))
            for (unused_id, aesurfs) in sorted(self.aesurfs.items()):
                bdf_files[aesurfs.ifile].write(aesurfs.write_card(size, is_double))
            for (unused_id, aefact) in sorted(self.aefacts.items()):
                bdf_files[aefact.ifile].write(aefact.write_card(size, is_double))

    def _write_static_aero_file(self, bdf_files, size=8, is_double=False, is_long_ids=None):
        # type: (Any, int, bool) -> None
        """Writes the static aero cards"""
        if self.aeros or self.trims or self.divergs:
            # static aero
            if self.aeros:
                bdf_files[self.aeros.ifile].write(self.aeros.write_card(size, is_double))
            for (unused_id, trim) in sorted(self.trims.items()):
                bdf_files[trim.ifile].write(trim.write_card(size, is_double))
            for (unused_id, diverg) in sorted(self.divergs.items()):
                bdf_files[diverg.ifile].write(diverg.write_card(size, is_double))


    def _write_flutter_file(self, bdf_files, size=8, is_double=False, write_aero_in_flutter=True,
                            is_long_ids=None):
        # type: (Any, int, bool, bool) -> None
        """Writes the flutter cards"""
        if (write_aero_in_flutter and self.aero) or self.flfacts or self.flutters or self.mkaeros:
            if write_aero_in_flutter:
                bdf_files[self.aero.ifile].write(self.aero.write_card(size, is_double))
            for (unused_id, flutter) in sorted(self.flutters.items()):
                bdf_files[flutter.ifile].write(flutter.write_card(size, is_double))
            for (unused_id, flfact) in sorted(self.flfacts.items()):
                bdf_files[flfact.ifile].write(flfact.write_card(size, is_double))
            for mkaero in self.mkaeros:
                bdf_files[mkaero.ifile].write(mkaero.write_card(size, is_double))

    def _write_gust_file(self, bdf_files, size=8, is_double=False, write_aero_in_gust=True, is_long_ids=None):
        # type: (Any, int, bool, bool) -> None
        """Writes the gust cards"""
        if (write_aero_in_gust and self.aero) or self.gusts:
            if write_aero_in_gust:
                for (unused_id, aero) in sorted(self.aero.items()):
                    bdf_files[aero.ifile].write(aero.write_card(size, is_double))
            for (unused_id, gust) in sorted(self.gusts.items()):
                bdf_files[gust.ifile].write(gust.write_card(size, is_double))

    def _write_common_file(self, bdf_files, size=8, is_double=False, is_long_ids=None):
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
        self._write_dmigs_file(bdf_files, size, is_double, is_long_ids=is_long_ids)
        self._write_loads_file(bdf_files, size, is_double, is_long_ids=is_long_ids)
        self._write_dynamic_file(bdf_files, size, is_double, is_long_ids=is_long_ids)
        self._write_aero_control_file(bdf_files, size, is_double, is_long_ids=is_long_ids)
        self._write_static_aero_file(bdf_files, size, is_double, is_long_ids=is_long_ids)

        write_aero_in_flutter, write_aero_in_gust = self._find_aero_location()
        self._write_flutter_file(bdf_files, size, is_double, write_aero_in_flutter, is_long_ids=is_long_ids)
        self._write_gust_file(bdf_files, size, is_double, write_aero_in_gust, is_long_ids=is_long_ids)

        self._write_thermal_file(bdf_files, size, is_double, is_long_ids=is_long_ids)
        self._write_thermal_materials_file(bdf_files, size, is_double, is_long_ids=is_long_ids)
        self._write_constraints_file(bdf_files, size, is_double, is_long_ids=is_long_ids)
        self._write_optimization_file(bdf_files, size, is_double, is_long_ids=is_long_ids)
        self._write_tables_file(bdf_files, size, is_double, is_long_ids=is_long_ids)
        self._write_sets_file(bdf_files, size, is_double, is_long_ids=is_long_ids)
        self._write_superelements_file(bdf_files, size, is_double, is_long_ids=is_long_ids)
        self._write_contact_file(bdf_files, size, is_double, is_long_ids=is_long_ids)
        self._write_rejects_file(bdf_files, size, is_double, is_long_ids=is_long_ids)
        self._write_coords_file(bdf_files, size, is_double, is_long_ids=is_long_ids)

    def _write_constraints_file(self, bdf_files, size=8, is_double=False, is_long_ids=None):
        # type: (Any, int, bool) -> None
        """Writes the constraint cards sorted by ID"""
        size, is_long_ids = self._write_mesh_long_ids_size(size, is_long_ids)
        if self.suport or self.suport1:
            for suport in self.suport:
                bdf_files[suport.ifile].write(suport.write_card(size, is_double))
            for unused_suport_id, suport in sorted(self.suport1.items()):
                bdf_files[suport.ifile].write(suport.write_card(size, is_double))

        if self.spcs or self.spcadds or self.spcoffs:
            #bdf_file.write('$SPCs\n')
            #str_spc = str(self.spcObject) # old
            #if str_spc:
                #bdf_file.write(str_spc)
            #else:
            for (unused_id, spcadds) in sorted(self.spcadds.items()):
                for spcadd in spcadds:
                    bdf_files[spcadd.ifile].write(str(spcadd))
            for (unused_id, spcs) in sorted(self.spcs.items()):
                for spc in spcs:
                    bdf_files[spc.ifile].write(str(spc))
            for (unused_id, spcoffs) in sorted(self.spcoffs.items()):
                for spc in spcoffs:
                    bdf_files[spc.ifile].write(str(spc))

        if self.mpcs or self.mpcadds:
            for (unused_id, mpcadds) in sorted(self.mpcadds.items()):
                for mpcadd in mpcadds:
                    bdf_files[mpcadd.ifile].write(str(mpcadd))
            for (unused_id, mpcs) in sorted(self.mpcs.items()):
                for mpc in mpcs:
                    bdf_files[mpc.ifile].write(mpc.write_card(size, is_double))

    def _write_contact_file(self, bdf_files, size=8, is_double=False, is_long_ids=None):
        # type: (Any, int, bool) -> None
        """Writes the contact cards sorted by ID"""
        is_contact = (self.bcrparas or self.bctadds or self.bctparas
                      or self.bctsets or self.bsurf or self.bsurfs)
        if is_contact:
            for (unused_id, bcrpara) in sorted(self.bcrparas.items()):
                bdf_files[bcrpara.ifile].write(bcrpara.write_card(size, is_double))
            for (unused_id, bctadds) in sorted(self.bctadds.items()):
                bdf_files[bctadds.ifile].write(bctadds.write_card(size, is_double))
            for (unused_id, bctpara) in sorted(self.bctparas.items()):
                bdf_files[bctpara.ifile].write(bctpara.write_card(size, is_double))

            for (unused_id, bctset) in sorted(self.bctsets.items()):
                bdf_files[bctset.ifile].write(bctset.write_card(size, is_double))
            for (unused_id, bsurfi) in sorted(self.bsurf.items()):
                bdf_files[bsurfi.ifile].write(bsurfi.write_card(size, is_double))
            for (unused_id, bsurfsi) in sorted(self.bsurfs.items()):
                bdf_files[bsurfsi.ifile].write(bsurfsi.write_card(size, is_double))

    def _write_coords_file(self, bdf_files, size=8, is_double=False, is_long_ids=None):
        # type: (Any, int, bool) -> None
        """Writes the coordinate cards in a sorted order"""
        size, is_long_ids = self._write_mesh_long_ids_size(size, is_long_ids)

        for (unused_id, coord) in sorted(self.coords.items()):
            if unused_id != 0:
                bdf_file = bdf_files[coord.ifile]
                try:
                    bdf_file.write(coord.write_card(size, is_double))
                except RuntimeError:
                    bdf_file.write(coord.write_card(16, is_double))

    def _write_dmigs_file(self, bdf_files, size=8, is_double=False, is_long_ids=None):
        # type: (Any, int, bool) -> None
        """
        Writes the DMIG cards

        Parameters
        ----------
        size : int
            large field (16) or small field (8)

        """
        for (unused_name, dmig) in sorted(self.dmigs.items()):
            bdf_files[dmig.ifile].write(dmig.write_card(size, is_double))
        for (unused_name, dmi) in sorted(self.dmis.items()):
            bdf_files[dmi.ifile].write(dmi.write_card(size, is_double))
        for (unused_name, dmij) in sorted(self.dmijs.items()):
            bdf_files[dmij.ifile].write(dmij.write_card(size, is_double))
        for (unused_name, dmiji) in sorted(self.dmijis.items()):
            bdf_files[dmiji.ifile].write(dmiji.write_card(size, is_double))
        for (unused_name, dmik) in sorted(self.dmiks.items()):
            bdf_files[dmik.ifile].write(dmik.write_card(size, is_double))

    def _write_dynamic_file(self, bdf_files, size=8, is_double=False, is_long_ids=None):
        # type: (Any, int, bool) -> None
        """Writes the dynamic cards sorted by ID"""
        is_dynamic = (self.dareas or self.dphases or self.nlparms or self.frequencies or
                      self.methods or self.cMethods or self.tsteps or self.tstepnls or
                      self.transfer_functions or self.delays or self.rotors or self.tics or
                      self.nlpcis)
        if is_dynamic:
            for (unused_id, method) in sorted(self.methods.items()):
                bdf_files[method.ifile].write(method.write_card(size, is_double))
            for (unused_id, cmethod) in sorted(self.cMethods.items()):
                bdf_files[cmethod.ifile].write(cmethod.write_card(size, is_double))
            for (unused_id, darea) in sorted(self.dareas.items()):
                bdf_files[darea.ifile].write(darea.write_card(size, is_double))
            for (unused_id, dphase) in sorted(self.dphases.items()):
                bdf_files[dphase.ifile].write(dphase.write_card(size, is_double))
            for (unused_id, nlparm) in sorted(self.nlparms.items()):
                bdf_files[nlparm.ifile].write(nlparm.write_card(size, is_double))
            for (unused_id, nlpci) in sorted(self.nlpcis.items()):
                bdf_files[nlpci.ifile].write(nlpci.write_card(size, is_double))
            for (unused_id, tstep) in sorted(self.tsteps.items()):
                bdf_files[tstep.ifile].write(tstep.write_card(size, is_double))
            for (unused_id, tstepnl) in sorted(self.tstepnls.items()):
                bdf_files[tstepnl.ifile].write(tstepnl.write_card(size, is_double))
            for (unused_id, freqs) in sorted(self.frequencies.items()):
                for freq in freqs:
                    bdf_files[freq.ifile].write(freq.write_card(size, is_double))
            for (unused_id, delay) in sorted(self.delays.items()):
                bdf_files[delay.ifile].write(delay.write_card(size, is_double))
            for (unused_id, rotor) in sorted(self.rotors.items()):
                bdf_files[rotor.ifile].write(rotor.write_card(size, is_double))
            for (unused_id, tic) in sorted(self.tics.items()):
                bdf_files[tic.ifile].write(tic.write_card(size, is_double))

            for (unused_id, tfs) in sorted(self.transfer_functions.items()):
                for transfer_function in tfs:
                    bdf_files[transfer_function.ifile].write(transfer_function.write_card(size, is_double))


    def _write_loads_file(self, bdf_files, size=8, is_double=False, is_long_ids=None):
        # type: (Any, int, bool) -> None
        """Writes the load cards sorted by ID"""
        size, is_long_ids = self._write_mesh_long_ids_size(size, is_long_ids)
        if self.load_combinations or self.loads or self.tempds:
            for (key, load_combinations) in sorted(self.load_combinations.items()):
                for load_combination in load_combinations:
                    try:
                        bdf_files[load_combination.ifile].write(load_combination.write_card(size, is_double))
                    except:
                        print('failed printing load...type=%s key=%r'
                              % (load_combination.type, key))
                        raise
            for (key, loadcase) in sorted(self.loads.items()):
                for load in loadcase:
                    try:
                        bdf_files[load.ifile].write(load.write_card(size, is_double))
                    except:
                        print('failed printing load...type=%s key=%r'
                              % (load.type, key))
                        raise
            for unused_key, tempd in sorted(self.tempds.items()):
                bdf_files[tempd.ifile].write(tempd.write_card(size, is_double))
        self._write_dloads_file(bdf_files, size=size, is_double=is_double, is_long_ids=is_long_ids)

    def _write_dloads_file(self, bdf_files, size=8, is_double=False, is_long_ids=None):
    # type: (Any, int, bool) -> None
        """Writes the dload cards sorted by ID"""
        size, is_long_ids = self._write_mesh_long_ids_size(size, is_long_ids)
        if self.dloads or self.dload_entries:
            for (key, loadcase) in sorted(self.dloads.items()):
                for load in loadcase:
                    try:
                        bdf_files[load.ifile].write(load.write_card(size, is_double))
                    except:
                        print('failed printing load...type=%s key=%r'
                              % (load.type, key))
                        raise

            for (key, loadcase) in sorted(self.dload_entries.items()):
                for load in loadcase:
                    try:
                        bdf_files[load.ifile].write(load.write_card(size, is_double))
                    except:
                        print('failed printing load...type=%s key=%r'
                              % (load.type, key))
                        raise


    def _write_masses_file(self, bdf_files, size=8, is_double=False, is_long_ids=None):
        # type: (Any, int, bool) -> None
        """Writes the mass cards sorted by ID"""
        size, is_long_ids = self._write_mesh_long_ids_size(size, is_long_ids)
        if self.properties_mass:
            for (pid, mass) in sorted(self.properties_mass.items()):
                try:
                    bdf_files[mass.ifile].write(mass.write_card(size, is_double))
                except:
                    print('failed printing mass property...'
                          'type=%s eid=%s' % (mass.type, pid))
                    raise

        if self.masses:
            for (eid, mass) in sorted(self.masses.items()):
                try:
                    bdf_files[mass.ifile].write(mass.write_card(size, is_double))
                except:
                    print('failed printing masses...'
                          'type=%s eid=%s' % (mass.type, eid))
                    raise

    def _write_materials_file(self, bdf_files, size=8, is_double=False, is_long_ids=None):
        # type: (Any, int, bool) -> None
        """Writes the materials in a sorted order"""
        size, is_long_ids = self._write_mesh_long_ids_size(size, is_long_ids)
        is_materials = (self.materials or self.hyperelastic_materials or self.creep_materials or
                        self.MATS1 or self.MATS3 or self.MATS8 or self.MATT1 or
                        self.MATT2 or self.MATT3 or self.MATT4 or self.MATT5 or
                        self.MATT8 or self.MATT9 or self.nxstrats)
        if is_materials:
            for (unused_mid, material) in sorted(self.materials.items()):
                bdf_files[material.ifile].write(material.write_card(size, is_double))
            for (unused_mid, material) in sorted(self.hyperelastic_materials.items()):
                bdf_files[material.ifile].write(material.write_card(size, is_double))
            for (unused_mid, material) in sorted(self.creep_materials.items()):
                bdf_files[material.ifile].write(material.write_card(size, is_double))

            for (unused_mid, material) in sorted(self.MATS1.items()):
                bdf_files[material.ifile].write(material.write_card(size, is_double))
            for (unused_mid, material) in sorted(self.MATS3.items()):
                bdf_files[material.ifile].write(material.write_card(size, is_double))
            for (unused_mid, material) in sorted(self.MATS8.items()):
                bdf_files[material.ifile].write(material.write_card(size, is_double))

            for (unused_mid, material) in sorted(self.MATT1.items()):
                bdf_files[material.ifile].write(material.write_card(size, is_double))
            for (unused_mid, material) in sorted(self.MATT2.items()):
                bdf_files[material.ifile].write(material.write_card(size, is_double))
            for (unused_mid, material) in sorted(self.MATT3.items()):
                bdf_files[material.ifile].write(material.write_card(size, is_double))
            for (unused_mid, material) in sorted(self.MATT4.items()):
                bdf_files[material.ifile].write(material.write_card(size, is_double))
            for (unused_mid, material) in sorted(self.MATT5.items()):
                bdf_files[material.ifile].write(material.write_card(size, is_double))
            for (unused_mid, material) in sorted(self.MATT8.items()):
                bdf_files[material.ifile].write(material.write_card(size, is_double))
            for (unused_mid, material) in sorted(self.MATT9.items()):
                bdf_files[material.ifile].write(material.write_card(size, is_double))
            for (unused_sid, nxstrat) in sorted(self.nxstrats.items()):
                bdf_files[material.ifile].write(nxstrat.write_card(size, is_double))

    def _write_nodes_file(self, bdf_files, size=8, is_double=False, is_long_ids=None):
        # type: (Any, int, bool) -> None
        """Writes the NODE-type cards"""
        if self.spoints:
            write_xpoints_file(bdf_files, 'SPOINT', self.spoints)
        if self.epoints:
            write_xpoints_file(bdf_files, 'EPOINT', self.epoints)
        if self.points:
            for unused_point_id, point in sorted(self.points.items()):
                bdf_files[point.ifile].write(point.write_card(size, is_double))

        if self._is_axis_symmetric:
            if self.axic:
                bdf_files[self.axic.ifile].write(self.axic.write_card(size, is_double))
            if self.axif:
                bdf_files[self.axif.ifile].write(self.axif.write_card(size, is_double))
            for unused_nid, ringax_pointax in sorted(self.ringaxs.items()):
                bdf_files[ringax_pointax.ifile].write(ringax_pointax.write_card(size, is_double))
            for unused_ringfl, ringfl in sorted(self.ringfl.items()):
                bdf_files[ringfl.ifile].write(ringfl.write_card(size, is_double))
            for unused_nid, gridb in sorted(self.gridb.items()):
                bdf_files[gridb.ifile].write(gridb.write_card(size, is_double))

        self._write_grids_file(bdf_files, size=size, is_double=is_double)
        if self.seqgp:
            bdf_files[self.seqgp.ifile].write(self.seqgp.write_card(size, is_double))

        #if 0:  # not finished
            #self._write_nodes_associated(bdf_file, size, is_double)

    def _write_grids_file(self, bdf_files, size=8, is_double=False, is_long_ids=None):
        # type: (Any, int, bool) -> None
        """Writes the GRID-type cards"""
        size, is_long_ids = self._write_mesh_long_ids_size(size, is_long_ids)
        if self.nodes:
            if self.grdset:
                bdf_files[self.grdset.ifile].write(self.grdset.write_card(size))
            if is_long_ids:
                for (unused_nid, node) in sorted(self.nodes.items()):
                    bdf_files[node.ifile].write(node.write_card_16(is_double))
            else:
                for (unused_nid, node) in sorted(self.nodes.items()):
                    bdf_files[node.ifile].write(node.write_card(size, is_double))

    def _write_optimization_file(self, bdf_files, size=8, is_double=False, is_long_ids=None):
        # type: (Any, int, bool) -> None
        """Writes the optimization cards sorted by ID"""
        is_optimization = (self.dconadds or self.dconstrs or self.desvars or self.ddvals or
                           self.dresps or
                           self.dvprels or self.dvmrels or self.dvcrels or self.doptprm or
                           self.dlinks or self.dequations or self.dtable is not None or
                           self.dvgrids or self.dscreen)
        if is_optimization:
            for (unused_id, dconadd) in sorted(self.dconadds.items()):
                bdf_files[dconadd.ifile].write(dconadd.write_card(size, is_double))
            for (unused_id, dconstrs) in sorted(self.dconstrs.items()):
                for dconstr in dconstrs:
                    bdf_files[dconstr.ifile].write(dconstr.write_card(size, is_double))
            for (unused_id, desvar) in sorted(self.desvars.items()):
                bdf_files[desvar.ifile].write(desvar.write_card(size, is_double))
            for (unused_id, ddval) in sorted(self.ddvals.items()):
                bdf_files[ddval.ifile].write(ddval.write_card(size, is_double))
            for (unused_id, dlink) in sorted(self.dlinks.items()):
                bdf_files[dlink.ifile].write(dlink.write_card(size, is_double))
            for (unused_id, dresp) in sorted(self.dresps.items()):
                bdf_files[dresp.ifile].write(dresp.write_card(size, is_double))

            for (unused_id, dvcrel) in sorted(self.dvcrels.items()):
                bdf_files[dvcrel.ifile].write(dvcrel.write_card(size, is_double))
            for (unused_id, dvmrel) in sorted(self.dvmrels.items()):
                bdf_files[dvmrel.ifile].write(dvmrel.write_card(size, is_double))
            for (unused_id, dvprel) in sorted(self.dvprels.items()):
                bdf_files[dvprel.ifile].write(dvprel.write_card(size, is_double))
            for (unused_id, dvgrids) in sorted(self.dvgrids.items()):
                for dvgrid in dvgrids:
                    bdf_files[dvgrid.ifile].write(dvgrid.write_card(size, is_double))
            for (unused_id, dscreen) in sorted(self.dscreen.items()):
                bdf_files[dscreen.ifile].write(str(dscreen))

            for (unused_id, equation) in sorted(self.dequations.items()):
                bdf_files[equation.ifile].write(str(equation))

            if self.dtable is not None:
                bdf_files[self.dtable.ifile].write(self.dtable.write_card(size, is_double))
            if self.doptprm is not None:
                bdf_files[self.doptprm.ifile].write(self.doptprm.write_card(size, is_double))

    def _write_params_file(self, bdf_files, size=8, is_double=False, is_long_ids=None):
        # type: (Any, int, bool) -> None
        """
        Writes the PARAM cards
        """
        size, is_long_ids = self._write_mesh_long_ids_size(size, is_long_ids)
        if self.params or self.dti:
            for unused_name, dti in sorted(self.dti.items()):
                bdf_files[dti.ifile].write(dti.write_card(size=size, is_double=is_double))

            for (unused_key, param) in sorted(self.params.items()):
                bdf_files[param.ifile].write(param.write_card(size, is_double))

    def _write_properties_file(self, bdf_files, size=8, is_double=False, is_long_ids=None):
        # type: (Any, int, bool) -> None
        """Writes the properties in a sorted order"""
        size, is_long_ids = self._write_mesh_long_ids_size(size, is_long_ids)
        if self.properties:
            prop_groups = (self.properties, self.pelast, self.pdampt, self.pbusht)
            if is_long_ids:
                for prop_group in prop_groups:
                    for unused_pid, prop in sorted(prop_group.items()):
                        bdf_files[prop.ifile].write(prop.write_card_16(is_double))
            else:
                for prop_group in prop_groups:
                    for unused_pid, prop in sorted(prop_group.items()):
                        bdf_files[prop.ifile].write(prop.write_card(size, is_double))

    def _write_rejects_file(self, bdf_files, size=8, is_double=False, is_long_ids=None):
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
            for reject_card in self.reject_cards:
                try:
                    bdf_files[0].write(print_func(reject_card))
                except RuntimeError:
                    for field in reject_card:
                        if field is not None and '=' in field:
                            raise SyntaxError('cannot reject equal signed '
                                              'cards\ncard=%s\n' % reject_card)
                    raise

        if self.reject_lines:
            for reject_lines in self.reject_lines:
                if isinstance(reject_lines, (list, tuple)):
                    for reject in reject_lines:
                        reject2 = reject.rstrip()
                        if reject2:
                            bdf_files[0].write('%s\n' % reject2)
                elif isinstance(reject_lines, string_types):
                    reject2 = reject_lines.rstrip()
                    if reject2:
                        bdf_files[0].write('%s\n' % reject2)
                else:
                    raise TypeError(reject_lines)

    def _write_rigid_elements_file(self, bdf_files, size=8, is_double=False, is_long_ids=None):
        # type: (Any, int, bool) -> None
        """Writes the rigid elements in a sorted order"""
        size, is_long_ids = self._write_mesh_long_ids_size(size, is_long_ids)
        if self.rigid_elements:
            if is_long_ids:
                for (eid, element) in sorted(self.rigid_elements.items()):
                    try:
                        bdf_files[element.ifile].write(element.write_card_16(is_double))
                    except:
                        print('failed printing element...'
                              'type=%s eid=%s' % (element.type, eid))
                        raise
            else:
                for (eid, element) in sorted(self.rigid_elements.items()):
                    try:
                        bdf_files[element.ifile].write(element.write_card(size, is_double))
                    except:
                        print('failed printing element...'
                              'type=%s eid=%s' % (element.type, eid))
                        raise
        if self.plotels:
            for (eid, element) in sorted(self.plotels.items()):
                bdf_files[element.ifile].write(element.write_card(size, is_double))

    def _write_sets_file(self, bdf_files, size=8, is_double=False, is_long_ids=None):
        # type: (Any, int, bool) -> None
        """Writes the SETx cards sorted by ID"""
        is_sets = (self.sets or self.asets or self.omits or self.bsets or self.csets or self.qsets
                   or self.usets)
        if is_sets:
            for (unused_id, set_obj) in sorted(self.sets.items()):  # dict
                bdf_files[set_obj.ifile].write(set_obj.write_card(size, is_double))
            for set_obj in self.asets:  # list
                bdf_files[set_obj.ifile].write(set_obj.write_card(size, is_double))
            for set_obj in self.omits:  # list
                bdf_files[set_obj.ifile].write(set_obj.write_card(size, is_double))
            for set_obj in self.bsets:  # list
                bdf_files[set_obj.ifile].write(set_obj.write_card(size, is_double))
            for set_obj in self.csets:  # list
                bdf_files[set_obj.ifile].write(set_obj.write_card(size, is_double))
            for set_obj in self.qsets:  # list
                bdf_files[set_obj.ifile].write(set_obj.write_card(size, is_double))
            for unused_name, usets in sorted(self.usets.items()):  # dict
                for set_obj in usets:  # list
                    bdf_files[set_obj.ifile].write(set_obj.write_card(size, is_double))

    def _write_superelements_file(self, bdf_files, size=8, is_double=False, is_long_ids=None):
        # type: (Any, int, bool) -> None
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
            for set_obj in self.se_bsets:  # list
                bdf_files[set_obj.ifile].write(set_obj.write_card(size, is_double))
            for set_obj in self.se_csets:  # list
                bdf_files[set_obj.ifile].write(set_obj.write_card(size, is_double))
            for set_obj in self.se_qsets:  # list
                bdf_files[set_obj.ifile].write(set_obj.write_card(size, is_double))
            for (unused_set_id, set_obj) in sorted(self.se_sets.items()):  # dict
                bdf_files[set_obj.ifile].write(set_obj.write_card(size, is_double))
            for unused_name, usets in sorted(self.se_usets.items()):  # dict
                for set_obj in usets:  # list
                    bdf_files[set_obj.ifile].write(set_obj.write_card(size, is_double))
            for se_suport in self.se_suport:  # list
                bdf_files[se_suport.ifile].write(se_suport.write_card(size, is_double))

        for unused_seid, csuper in sorted(self.csuper.items()):
            bdf_files[csuper.ifile].write(csuper.write_card(size, is_double))
        for unused_seid, csupext in sorted(self.csupext.items()):
            bdf_files[csupext.ifile].write(csupext.write_card(size, is_double))

        for unused_seid, sebndry in sorted(self.sebndry.items()):
            bdf_files[sebndry.ifile].write(sebndry.write_card(size, is_double))
        for unused_seid, seelt in sorted(self.seelt.items()):
            bdf_files[seelt.ifile].write(seelt.write_card(size, is_double))
        for unused_seid, seexcld in sorted(self.seexcld.items()):
            bdf_files[seexcld.ifile].write(seexcld.write_card(size, is_double))

        for unused_seid, seloc in sorted(self.seloc.items()):
            bdf_files[seloc.ifile].write(seloc.write_card(size, is_double))
        for unused_seid, seload in sorted(self.seload.items()):
            bdf_files[seload.ifile].write(seload.write_card(size, is_double))
        for unused_seid, sempln in sorted(self.sempln.items()):
            bdf_files[sempln.ifile].write(sempln.write_card(size, is_double))
        for unused_setid, senqset in sorted(self.senqset.items()):
            bdf_files[senqset.ifile].write(senqset.write_card(size, is_double))
        for unused_seid, setree in sorted(self.setree.items()):
            bdf_files[setree.ifile].write(setree.write_card(size, is_double))


    def _write_tables_file(self, bdf_files, size=8, is_double=False, is_long_ids=None):
        # type: (Any, int, bool) -> None
        """Writes the TABLEx cards sorted by ID"""
        if self.tables or self.tables_d or self.tables_m or self.tables_sdamping:
            for (unused_id, table) in sorted(self.tables.items()):
                bdf_files[table.ifile].write(table.write_card(size, is_double))
            for (unused_id, table) in sorted(self.tables_d.items()):
                bdf_files[table.ifile].write(table.write_card(size, is_double))
            for (unused_id, table) in sorted(self.tables_m.items()):
                bdf_files[table.ifile].write(table.write_card(size, is_double))
            for (unused_id, table) in sorted(self.tables_sdamping.items()):
                bdf_files[table.ifile].write(table.write_card(size, is_double))

        if self.random_tables:
            for (unused_id, table) in sorted(self.random_tables.items()):
                bdf_files[table.ifile].write(table.write_card(size, is_double))

    def _write_thermal_file(self, bdf_files, size=8, is_double=False, is_long_ids=None):
        # type: (Any, int, bool) -> None
        """Writes the thermal cards"""
        # PHBDY
        is_thermal = (self.phbdys or self.convection_properties or self.bcs or
                      self.views or self.view3ds or self.radset or self.radcavs)
        if is_thermal:
            for (unused_key, phbdy) in sorted(self.phbdys.items()):
                bdf_files[phbdy.ifile].write(phbdy.write_card(size, is_double))

            #for unused_key, prop in sorted(self.thermal_properties.items()):
            #    bdf_file.write(str(prop))
            for (unused_key, prop) in sorted(self.convection_properties.items()):
                bdf_files[prop.ifile].write(prop.write_card(size, is_double))

            # BCs
            for (unused_key, bcs) in sorted(self.bcs.items()):
                for boundary_condition in bcs:  # list
                    bdf_files[boundary_condition.ifile].write(boundary_condition.write_card(size, is_double))

            for (unused_key, view) in sorted(self.views.items()):
                bdf_files[view.ifile].write(view.write_card(size, is_double))
            for (unused_key, view3d) in sorted(self.view3ds.items()):
                bdf_files[view3d.ifile].write(view3d.write_card(size, is_double))
            if self.radset:
                bdf_files[self.radset.ifile].write(self.radset.write_card(size, is_double))
            for unused_icavity, radcav in self.radcavs.items():
                bdf_files[radcav.ifile].write(radcav.write_card(size, is_double))


    def _write_thermal_materials_file(self, bdf_files, size=8, is_double=False, is_long_ids=None):
        # type: (Any, int, bool) -> None
        """Writes the thermal materials in a sorted order"""
        if self.thermal_materials:
            for (unused_mid, material) in sorted(self.thermal_materials.items()):
                bdf_files[material.ifile].write(material.write_card(size, is_double))

def write_xpoints_file(bdf_files, cardtype, points, comment=''):
    """writes SPOINTs/EPOINTs"""
    assert isinstance(points, dict), points
    for point in points:
        bdf_files[point.ifile].write(point.write_card())
