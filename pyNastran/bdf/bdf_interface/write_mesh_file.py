# coding: utf-8
"""
This file defines:
  - WriteMesh

"""
from __future__ import annotations
import os
from typing import Any, Union, Optional, Any
from collections import defaultdict
from typing import Optional, TYPE_CHECKING
if TYPE_CHECKING:  # pragma: no cover
    from io import StringIO

import numpy as np
from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.field_writer_16 import print_card_16
from pyNastran.bdf.bdf_interface.write_mesh import WriteMesh
from pyNastran.bdf.write_path import write_include


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

    def write_bdfs(self, out_filenames: Optional[Union[str, StringIO]],
                   relative_dirname: Optional[str]=None, encoding: Optional[str]=None,
                   size: int=8, is_double: bool=False,
                   enddata: Optional[bool]=None, close: bool=True,
                   is_windows: Optional[bool]=None) -> None:
        """
        Writes the BDF.

        Parameters
        ----------
        out_filename : varies; default=None
            str        - the name to call the output bdf
            file       - a file object
            StringIO() - a StringIO object
            None       - pops a dialog
        relative_dirname : str; default=None -> os.curdir
            A relative path to reference INCLUDEs.
            ''   : relative to the main bdf
            None : use the current directory
            path : absolute path
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
        is_windows : bool; default=None
            True/False : Windows has a special format for writing INCLUDE
                files, so the format for a BDF that will run on Linux and
                Windows is different.
            None : Check the platform
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

        ifile_out_filenames = _map_filenames_to_ifile_filname_dict(
            out_filenames, self.active_filenames)
        ifile0 = list(sorted(ifile_out_filenames))[0]
        #print('ifile_out_filenames =', ifile_out_filenames)

        out_filename0 = ifile_out_filenames[ifile0]
        #print("out_filename0 =", out_filename0)

        interspersed = False
        out_filename = self._output_helper(out_filename0,
                                           interspersed, size, is_double)
        self.log.debug('---starting BDF.write_bdf of %s---' % out_filename)
        encoding = self.get_encoding(encoding)

        #class DevNull:
            #def write(self, *_):
                #pass

        #devnull = DevNull()

        bdf_files, bdf_file0 = _open_bdf_files(ifile_out_filenames, self.active_filenames, encoding)

        if bdf_file0 is not None:
            self._write_header(bdf_file0, encoding)

        self._write_bdf_includes(out_filenames, bdf_files, relative_dirname=relative_dirname,
                                 is_windows=is_windows)

        self._write_params_file(bdf_files, size, is_double, is_long_ids=is_long_ids)
        self._write_nodes_file(bdf_files, size, is_double, is_long_ids=is_long_ids)

        self._write_elements_file(bdf_files, size, is_double, is_long_ids=is_long_ids)
        self._write_properties_file(bdf_files, size, is_double, is_long_ids=is_long_ids)
        self._write_materials_file(bdf_files, size, is_double, is_long_ids=is_long_ids)

        self._write_masses_file(bdf_files, size, is_double, is_long_ids=is_long_ids)
        self._write_rigid_elements_file(bdf_files, size, is_double, is_long_ids=is_long_ids)
        self._write_aero_file(bdf_files, size, is_double, is_long_ids=is_long_ids)

        self._write_common_file(bdf_files, size, is_double, is_long_ids=is_long_ids)
        if (enddata is None and 'ENDDATA' in self.card_count) or enddata:
            if bdf_file0:
                bdf_file0.write('ENDDATA\n')
        if close:
            for bdf_file in bdf_files.values():
                if bdf_file is not None:
                    bdf_file.close()
        del bdf_files

    def _write_bdf_includes(self, out_filenames, bdf_files, relative_dirname=None, is_windows=True):
        """
        Writes the INCLUDE files

        Parameters
        ----------
        out_filenames : dict[fname] : fname2
            fname_in - the nominal bdf that was read
            fname_out - the bdf that will be written
        relative_dirname : str; default=None -> os.curdir
            A relative path to reference INCLUDEs.
            ''   : relative to the main bdf
            None : use the current directory
            path : absolute path
        is_windows : bool; default=None
            True/False : Windows has a special format for writing INCLUDE
                files, so the format for a BDF that will run on Linux and
                Windows is different.
            None : Check the platform
        """
        if relative_dirname is None:
            relative_dirname = os.curdir
        elif relative_dirname == '':
            out_filename0 = list(out_filenames.keys())[0]
            relative_dirname = os.path.dirname(os.path.abspath(out_filename0))
            self.log.debug('relative_dirname = %s' % relative_dirname)

        self.log.debug('include_filenames:')
        for ifile, include_filenames in self.include_filenames.items():
            self.log.debug('ifile=%i %s' % (ifile, include_filenames))
            assert len(include_filenames) > 0, include_filenames
            bdf_file = bdf_files[ifile]
            if bdf_file is None:
                continue
            #self.log.info('ifile=%s include_files=%s' % (ifile, include_filenames))
            for include_filename in include_filenames:
                assert len(include_filename) > 0, include_filename
                #print('***', include_filename, '***')

                mapped_include_filename = include_filename
                if include_filename in out_filenames:
                    mapped_include_filename = out_filenames[include_filename]

                if relative_dirname == '':
                    # absolute path
                    rel_include_filename = mapped_include_filename
                else:
                    rel_include_filename = os.path.relpath(mapped_include_filename, relative_dirname)
                bdf_file.write(write_include(rel_include_filename, is_windows=is_windows))
                #print('* %r *' % (include_filename))
                #print('** %r **' % (relative_dirname))
                #print('***', rel_include_filename, '***', '')
            #msg = '\n        '.join(include_lines) + '\n'
            #print(msg)
            #bdf_file.write(msg)

    def _write_elements_file(self, bdf_files: Any, size: int=8, is_double: bool=False,
                             is_long_ids: Optional[bool]=None) -> None:
        """
        Writes the elements in a sorted order
        """
        size, is_long_ids = self._write_mesh_long_ids_size(size, is_long_ids)

        if self.elements:
            write_bdfs_dict(bdf_files, self.elements, size, is_double, is_long_ids)

        if self.ao_element_flags:
            write_bdfs_dict(bdf_files, self.ao_element_flags, size, is_double, is_long_ids)
        if self.normals:
            write_bdfs_dict(bdf_files, self.normals, size, is_double, is_long_ids)
        self._write_nsm_file(bdf_files, size, is_double)

    def _write_nsm_file(self, bdf_files: Any, size: int=8, is_double: bool=False,
                        is_long_ids: Optional[bool]=None) -> None:
        """
        Writes the nsm in a sorted order
        """
        if self.nsms or self.nsmadds:
            write_bdfs_dict(bdf_files, self.nsmadds, size, is_double, is_long_ids)
            write_bdfs_dict(bdf_files, self.nsms, size, is_double, is_long_ids)

    def _write_aero_file(self, bdf_files: Any, size: int=8, is_double: bool=False,
                         is_long_ids: Optional[bool]=None) -> None:
        """Writes the aero cards"""
        if self.caeros or self.paeros or self.monitor_points or self.splines:
            write_bdfs_dict(bdf_files, self.caeros, size, is_double, is_long_ids)
            write_bdfs_dict(bdf_files, self.splines, size, is_double, is_long_ids)
            write_bdfs_dict(bdf_files, self.nsmadds, size, is_double, is_long_ids)
            write_bdfs_dict(bdf_files, self.nsmadds, size, is_double, is_long_ids)
            for monitor_point in self.monitor_points:
                bdf_files[monitor_point.ifile].write(monitor_point.write_card(size, is_double))
        self.zona.write_bdf(bdf_files[0], size=8, is_double=False)

    def _write_aero_control_file(self, bdf_files: Any, size: int=8, is_double: bool=False,
                                 is_long_ids: Optional[bool]=None) -> None:
        """Writes the aero control surface cards"""
        if(self.aecomps or self.aefacts or self.aeparams or self.aelinks or
           self.aelists or self.aestats or self.aesurf or self.aesurfs):
            write_bdfs_dict_list(bdf_files, self.aecomps, size, is_double, is_long_ids)

            write_bdfs_dict(bdf_files, self.aecomps, size, is_double, is_long_ids)
            write_bdfs_dict(bdf_files, self.aeparams, size, is_double, is_long_ids)
            write_bdfs_dict(bdf_files, self.aestats, size, is_double, is_long_ids)
            write_bdfs_dict(bdf_files, self.aelists, size, is_double, is_long_ids)
            write_bdfs_dict(bdf_files, self.aesurf, size, is_double, is_long_ids)
            write_bdfs_dict(bdf_files, self.aesurfs, size, is_double, is_long_ids)
            write_bdfs_dict(bdf_files, self.aefacts, size, is_double, is_long_ids)

    def _write_static_aero_file(self, bdf_files: Any, size: int=8, is_double: bool=False,
                                is_long_ids: Optional[bool]=None) -> None:
        """Writes the static aero cards"""
        if self.aeros or self.trims or self.divergs:
            # static aero
            if self.aeros:
                bdf_files[self.aeros.ifile].write(self.aeros.write_card(size, is_double))

            write_bdfs_dict(bdf_files, self.trims, size, is_double, is_long_ids)
            write_bdfs_dict(bdf_files, self.divergs, size, is_double, is_long_ids)


    def _write_flutter_file(self, bdf_files: Any, size: int=8, is_double: bool=False,
                            write_aero_in_flutter: bool=True,
                            is_long_ids: Optional[bool]=None) -> None:
        """Writes the flutter cards"""
        if (write_aero_in_flutter and self.aero) or self.flfacts or self.flutters or self.mkaeros:
            if write_aero_in_flutter:
                bdf_files[self.aero.ifile].write(self.aero.write_card(size, is_double))
            write_bdfs_dict(bdf_files, self.flutters, size, is_double, is_long_ids)
            write_bdfs_dict(bdf_files, self.flfacts, size, is_double, is_long_ids)
            write_bdfs_list(bdf_files, self.mkaeros, size, is_double, is_long_ids)

    def _write_gust_file(self, bdf_files: Any, size: int=8, is_double: bool=False,
                         write_aero_in_gust: bool=True, is_long_ids: Optional[bool]=None) -> None:
        """Writes the gust cards"""
        if (write_aero_in_gust and self.aero) or self.gusts:
            if write_aero_in_gust:
                for (unused_id, aero) in sorted(self.aero.items()):
                    bdf_files[aero.ifile].write(aero.write_card(size, is_double))
            write_bdfs_dict(bdf_files, self.gusts, size, is_double, is_long_ids)

    def _write_common_file(self, bdf_files: Any, size: int=8, is_double: bool=False,
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
        self._write_dmigs_file(bdf_files, size, is_double, is_long_ids=is_long_ids)
        self._write_loads_file(bdf_files, size, is_double, is_long_ids=is_long_ids)
        self._write_dynamic_file(bdf_files, size, is_double, is_long_ids=is_long_ids)
        self._write_aero_control_file(bdf_files, size, is_double, is_long_ids=is_long_ids)
        self._write_static_aero_file(bdf_files, size, is_double, is_long_ids=is_long_ids)

        write_aero_in_flutter, write_aero_in_gust = self._find_aero_location()
        self._write_flutter_file(bdf_files, size, is_double, write_aero_in_flutter,
                                 is_long_ids=is_long_ids)
        self._write_gust_file(bdf_files, size, is_double, write_aero_in_gust,
                              is_long_ids=is_long_ids)

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

    def _write_constraints_file(self, bdf_files: Any, size: int=8, is_double: bool=False,
                                is_long_ids: Optional[bool]=None) -> None:
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
            write_bdfs_dict_list(bdf_files, self.spcadds, size, is_double, is_long_ids)
            write_bdfs_dict_list(bdf_files, self.spcs, size, is_double, is_long_ids)
            write_bdfs_dict_list(bdf_files, self.spcoffs, size, is_double, is_long_ids)


        if self.mpcs or self.mpcadds:
            write_bdfs_dict_list(bdf_files, self.mpcadds, size, is_double, is_long_ids)
            write_bdfs_dict_list(bdf_files, self.mpcs, size, is_double, is_long_ids)

    def _write_contact_file(self, bdf_files: Any, size: int=8, is_double: bool=False,
                            is_long_ids: Optional[bool]=None) -> None:
        """Writes the contact cards sorted by ID"""
        is_contact = (self.bcrparas or self.bctadds or self.bctparas
                      or self.bctsets or self.bsurf or self.bsurfs
                      or self.bconp or self.blseg or self.bfric)
        if is_contact:
            write_bdfs_dict(bdf_files, self.bcrparas, size, is_double, is_long_ids)
            write_bdfs_dict(bdf_files, self.bctadds, size, is_double, is_long_ids)
            write_bdfs_dict(bdf_files, self.bctparas, size, is_double, is_long_ids)
            write_bdfs_dict(bdf_files, self.bctsets, size, is_double, is_long_ids)
            write_bdfs_dict(bdf_files, self.bsurf, size, is_double, is_long_ids)
            write_bdfs_dict(bdf_files, self.bsurfs, size, is_double, is_long_ids)

            write_bdfs_dict(bdf_files, self.bconp, size, is_double, is_long_ids)
            write_bdfs_dict(bdf_files, self.blseg, size, is_double, is_long_ids)
            write_bdfs_dict(bdf_files, self.bfric, size, is_double, is_long_ids)

    def _write_coords_file(self, bdf_files: Any, size: int=8, is_double: bool=False,
                           is_long_ids: Optional[bool]=None) -> None:
        """Writes the coordinate cards in a sorted order"""
        size, is_long_ids = self._write_mesh_long_ids_size(size, is_long_ids)

        for (unused_id, coord) in sorted(self.coords.items()):
            if unused_id != 0:
                bdf_file = bdf_files[coord.ifile]
                try:
                    bdf_file.write(coord.write_card(size, is_double))
                except RuntimeError:
                    bdf_file.write(coord.write_card(16, is_double))

    def _write_dmigs_file(self, bdf_files: Any, size: int=8, is_double: bool=False,
                          is_long_ids: Optional[bool]=None) -> None:
        """
        Writes the DMIG cards

        Parameters
        ----------
        size : int
            large field (16) or small field (8)

        """
        write_bdfs_dict(bdf_files, self.dmig, size, is_double, is_long_ids)
        write_bdfs_dict(bdf_files, self.dmi, size, is_double, is_long_ids)
        write_bdfs_dict(bdf_files, self.dmij, size, is_double, is_long_ids)
        write_bdfs_dict(bdf_files, self.dmiji, size, is_double, is_long_ids)
        write_bdfs_dict(bdf_files, self.dmik, size, is_double, is_long_ids)
        write_bdfs_dict(bdf_files, self.dmiax, size, is_double, is_long_ids)

    def _write_dynamic_file(self, bdf_files: Any, size: int=8, is_double: bool=False,
                            is_long_ids: Optional[bool]=None) -> None:
        """Writes the dynamic cards sorted by ID"""
        is_dynamic = (self.dareas or self.dphases or self.nlparms or self.frequencies or
                      self.methods or self.cMethods or self.tsteps or self.tstepnls or
                      self.transfer_functions or self.delays or self.rotors or self.tics or
                      self.nlpcis)
        if is_dynamic:
            write_bdfs_dict(bdf_files, self.methods, size, is_double, is_long_ids)
            write_bdfs_dict(bdf_files, self.cMethods, size, is_double, is_long_ids)
            write_bdfs_dict(bdf_files, self.dareas, size, is_double, is_long_ids)
            write_bdfs_dict(bdf_files, self.dphases, size, is_double, is_long_ids)
            write_bdfs_dict(bdf_files, self.nlparms, size, is_double, is_long_ids)
            write_bdfs_dict(bdf_files, self.nlpcis, size, is_double, is_long_ids)
            write_bdfs_dict(bdf_files, self.tsteps, size, is_double, is_long_ids)
            write_bdfs_dict(bdf_files, self.tstepnls, size, is_double, is_long_ids)

            write_bdfs_dict_list(bdf_files, self.frequencies, size, is_double, is_long_ids)

            write_bdfs_dict(bdf_files, self.delays, size, is_double, is_long_ids)
            write_bdfs_dict(bdf_files, self.rotors, size, is_double, is_long_ids)
            write_bdfs_dict(bdf_files, self.tics, size, is_double, is_long_ids)

            write_bdfs_dict_list(bdf_files, self.transfer_functions, size, is_double, is_long_ids)


    def _write_loads_file(self, bdf_files: Any, size: int=8, is_double: bool=False,
                          is_long_ids: Optional[bool]=None) -> None:
        """Writes the load cards sorted by ID"""
        size, is_long_ids = self._write_mesh_long_ids_size(size, is_long_ids)
        if self.load_combinations or self.loads or self.tempds:
            write_bdfs_dict_list(bdf_files, self.load_combinations, size, is_double, is_long_ids)
            write_bdfs_dict_list(bdf_files, self.loads, size, is_double, is_long_ids)
            write_bdfs_dict(bdf_files, self.tempds, size, is_double, is_long_ids)
        self._write_dloads_file(bdf_files, size=size, is_double=is_double, is_long_ids=is_long_ids)

    def _write_dloads_file(self, bdf_files: Any, size: int=8, is_double: bool=False,
                           is_long_ids: Optional[bool]=None) -> None:
        """Writes the dload cards sorted by ID"""
        size, is_long_ids = self._write_mesh_long_ids_size(size, is_long_ids)
        if self.dloads or self.dload_entries:
            write_bdfs_dict_list(bdf_files, self.dloads, size, is_double, is_long_ids)
            write_bdfs_dict_list(bdf_files, self.dload_entries, size, is_double, is_long_ids)

    def _write_masses_file(self, bdf_files: Any, size: int=8, is_double: bool=False,
                           is_long_ids: Optional[bool]=None) -> None:
        """Writes the mass cards sorted by ID"""
        size, is_long_ids = self._write_mesh_long_ids_size(size, is_long_ids)
        if self.properties_mass:
            write_bdfs_dict(bdf_files, self.properties_mass, size, is_double, is_long_ids)
        if self.masses:
            write_bdfs_dict(bdf_files, self.masses, size, is_double, is_long_ids)

    def _write_materials_file(self, bdf_files: Any, size: int=8, is_double: bool=False,
                              is_long_ids: Optional[bool]=None) -> None:
        """Writes the materials in a sorted order"""
        size, is_long_ids = self._write_mesh_long_ids_size(size, is_long_ids)
        is_materials = (self.materials or self.hyperelastic_materials or self.creep_materials or
                        self.MATS1 or self.MATS3 or self.MATS8 or self.MATT1 or
                        self.MATT2 or self.MATT3 or self.MATT4 or self.MATT5 or
                        self.MATT8 or self.MATT9 or self.nxstrats)
        if is_materials:
            write_bdfs_dict(bdf_files, self.materials, size, is_double, is_long_ids)
            write_bdfs_dict(bdf_files, self.hyperelastic_materials, size, is_double, is_long_ids)
            write_bdfs_dict(bdf_files, self.creep_materials, size, is_double, is_long_ids)
            write_bdfs_dict(bdf_files, self.MATS1, size, is_double, is_long_ids)
            write_bdfs_dict(bdf_files, self.MATS3, size, is_double, is_long_ids)
            write_bdfs_dict(bdf_files, self.MATS8, size, is_double, is_long_ids)
            write_bdfs_dict(bdf_files, self.MATT1, size, is_double, is_long_ids)
            write_bdfs_dict(bdf_files, self.MATT2, size, is_double, is_long_ids)
            write_bdfs_dict(bdf_files, self.MATT3, size, is_double, is_long_ids)
            write_bdfs_dict(bdf_files, self.MATT4, size, is_double, is_long_ids)
            write_bdfs_dict(bdf_files, self.MATT5, size, is_double, is_long_ids)
            write_bdfs_dict(bdf_files, self.MATT8, size, is_double, is_long_ids)
            write_bdfs_dict(bdf_files, self.MATT9, size, is_double, is_long_ids)
            write_bdfs_dict(bdf_files, self.nxstrats, size, is_double, is_long_ids)

    def _write_nodes_file(self, bdf_files: Any, size: int=8, is_double: bool=False,
                          is_long_ids: Optional[bool]=None) -> None:
        """Writes the NODE-type cards"""
        if self.spoints:
            write_xpoints_file(bdf_files, 'SPOINT', self.spoints)
        if self.epoints:
            write_xpoints_file(bdf_files, 'EPOINT', self.epoints)
        if self.points:
            write_bdfs_dict(bdf_files, self.points, size, is_double, is_long_ids)

        if self._is_axis_symmetric:
            if self.axic:
                bdf_files[self.axic.ifile].write(self.axic.write_card(size, is_double))
            if self.axif:
                bdf_files[self.axif.ifile].write(self.axif.write_card(size, is_double))
            write_bdfs_dict(bdf_files, self.ringaxs, size, is_double, is_long_ids)
            write_bdfs_dict(bdf_files, self.ringfl, size, is_double, is_long_ids)
            write_bdfs_dict(bdf_files, self.gridb, size, is_double, is_long_ids)

        self._write_grids_file(bdf_files, size=size, is_double=is_double)
        if self.seqgp:
            bdf_files[self.seqgp.ifile].write(self.seqgp.write_card(size, is_double))

        #if 0:  # not finished
            #self._write_nodes_associated(bdf_file, size, is_double)

    def _write_grids_file(self, bdf_files: Any, size: int=8, is_double: bool=False,
                          is_long_ids: Optional[bool]=None) -> None:
        """Writes the GRID-type cards"""
        size, is_long_ids = self._write_mesh_long_ids_size(size, is_long_ids)
        if self.nodes:
            if self.grdset:
                bdf_files[self.grdset.ifile].write(self.grdset.write_card(size))
            write_bdfs_dict(bdf_files, self.nodes, size, is_double, is_long_ids)

    def _write_optimization_file(self, bdf_files: Any, size: int=8, is_double: bool=False,
                                 is_long_ids: Optional[bool]=None) -> None:
        """Writes the optimization cards sorted by ID"""
        is_optimization = (self.dconadds or self.dconstrs or self.desvars or self.ddvals or
                           self.dresps or
                           self.dvprels or self.dvmrels or self.dvcrels or self.doptprm or
                           self.dlinks or self.dequations or self.dtable is not None or
                           self.dvgrids or self.dscreen or self.topvar)
        if is_optimization:
            write_bdfs_dict(bdf_files, self.dconadds, size, is_double, is_long_ids)
            write_bdfs_dict_list(bdf_files, self.dconadds, size, is_double, is_long_ids)

            write_bdfs_dict(bdf_files, self.desvars, size, is_double, is_long_ids)
            write_bdfs_dict(bdf_files, self.topvar, size, is_double, is_long_ids)
            write_bdfs_dict(bdf_files, self.ddvals, size, is_double, is_long_ids)
            write_bdfs_dict(bdf_files, self.dlinks, size, is_double, is_long_ids)
            write_bdfs_dict(bdf_files, self.dresps, size, is_double, is_long_ids)
            write_bdfs_dict(bdf_files, self.dvcrels, size, is_double, is_long_ids)
            write_bdfs_dict(bdf_files, self.dvmrels, size, is_double, is_long_ids)
            write_bdfs_dict(bdf_files, self.dvprels, size, is_double, is_long_ids)

            write_bdfs_dict_list(bdf_files, self.dvgrids, size, is_double, is_long_ids)

            for (unused_id, dscreen) in sorted(self.dscreen.items()):
                bdf_files[dscreen.ifile].write(str(dscreen))

            for (unused_id, equation) in sorted(self.dequations.items()):
                bdf_files[equation.ifile].write(str(equation))

            if self.dtable is not None:
                bdf_files[self.dtable.ifile].write(self.dtable.write_card(size, is_double))
            if self.doptprm is not None:
                bdf_files[self.doptprm.ifile].write(self.doptprm.write_card(size, is_double))
            if self.modtrak is not None:
                bdf_files[self.modtrak.ifile].write(self.modtrak.write_card(size, is_double))

    def _write_params_file(self, bdf_files: Any, size: int=8, is_double: bool=False,
                           is_long_ids: Optional[bool]=None) -> None:
        """
        Writes the PARAM cards
        """
        size, is_long_ids = self._write_mesh_long_ids_size(size, is_long_ids)
        if self.params or self.dti:
            write_bdfs_dict(bdf_files, self.params, size, is_double, is_long_ids)
            write_bdfs_dict(bdf_files, self.dti, size, is_double, is_long_ids)

    def _write_properties_file(self, bdf_files: Any, size: int=8, is_double: bool=False,
                               is_long_ids: Optional[bool]=None) -> None:
        """Writes the properties in a sorted order"""
        size, is_long_ids = self._write_mesh_long_ids_size(size, is_long_ids)
        is_properties = self.properties or self.pelast or self.pdampt or self.pbusht
        if is_properties:
            write_bdfs_dict(bdf_files, self.properties, size, is_double, is_long_ids)
            write_bdfs_dict(bdf_files, self.pelast, size, is_double, is_long_ids)
            write_bdfs_dict(bdf_files, self.pdampt, size, is_double, is_long_ids)
            write_bdfs_dict(bdf_files, self.pbusht, size, is_double, is_long_ids)

    def _write_rejects_file(self, bdf_files: Any, size: int=8, is_double: bool=False,
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
            #print(self.reject_lines)
            for reject_lines in self.reject_lines:
                if isinstance(reject_lines, (list, tuple)):
                    for reject in reject_lines:
                        reject2 = reject.rstrip()
                        if reject2:
                            bdf_files[0].write('%s\n' % reject2)
                elif isinstance(reject_lines, str):
                    reject2 = reject_lines.rstrip()
                    if reject2:
                        bdf_files[0].write('%s\n' % reject2)
                else:
                    raise TypeError(reject_lines)

    def _write_rigid_elements_file(self, bdf_files: Any, size: int=8, is_double: bool=False,
                           is_long_ids: Optional[bool]=None) -> None:
        """Writes the rigid elements in a sorted order"""
        size, is_long_ids = self._write_mesh_long_ids_size(size, is_long_ids)
        if self.rigid_elements:
            write_bdfs_dict(bdf_files, self.rigid_elements, size, is_double, is_long_ids)

        if self.plotels:
            write_bdfs_dict(bdf_files, self.plotels, size, is_double, is_long_ids)

    def _write_sets_file(self, bdf_files: Any, size: int=8, is_double: bool=False,
                           is_long_ids: Optional[bool]=None) -> None:
        """Writes the SETx cards sorted by ID"""
        is_sets = (self.sets or self.asets or self.omits or self.bsets or self.csets or self.qsets
                   or self.usets)
        if is_sets:
            write_bdfs_dict(bdf_files, self.sets, size, is_double, is_long_ids)
            write_bdfs_list(bdf_files, self.asets, size, is_double, is_long_ids)
            write_bdfs_list(bdf_files, self.omits, size, is_double, is_long_ids)
            write_bdfs_list(bdf_files, self.bsets, size, is_double, is_long_ids)
            write_bdfs_list(bdf_files, self.csets, size, is_double, is_long_ids)
            write_bdfs_list(bdf_files, self.qsets, size, is_double, is_long_ids)

            write_bdfs_dict_list(bdf_files, self.usets, size, is_double, is_long_ids)

    def _write_superelements_file(self, bdf_files: Any, size: int=8, is_double: bool=False,
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
            write_bdfs_list(bdf_files, self.se_bsets, size, is_double, is_long_ids)
            write_bdfs_list(bdf_files, self.se_csets, size, is_double, is_long_ids)
            write_bdfs_list(bdf_files, self.se_qsets, size, is_double, is_long_ids)

            write_bdfs_dict(bdf_files, self.se_sets, size, is_double, is_long_ids)
            write_bdfs_dict_list(bdf_files, self.se_usets, size, is_double, is_long_ids)
            write_bdfs_list(bdf_files, self.se_suport, size, is_double, is_long_ids)

        write_bdfs_dict(bdf_files, self.csuper, size, is_double, is_long_ids)
        write_bdfs_dict(bdf_files, self.csupext, size, is_double, is_long_ids)

        write_bdfs_dict(bdf_files, self.sebndry, size, is_double, is_long_ids)
        write_bdfs_dict(bdf_files, self.sebulk, size, is_double, is_long_ids)
        write_bdfs_dict(bdf_files, self.seconct, size, is_double, is_long_ids)
        write_bdfs_dict(bdf_files, self.seelt, size, is_double, is_long_ids)
        write_bdfs_dict(bdf_files, self.seexcld, size, is_double, is_long_ids)
        write_bdfs_dict(bdf_files, self.seloc, size, is_double, is_long_ids)
        write_bdfs_dict(bdf_files, self.seload, size, is_double, is_long_ids)
        write_bdfs_dict(bdf_files, self.sempln, size, is_double, is_long_ids)
        write_bdfs_dict(bdf_files, self.senqset, size, is_double, is_long_ids)
        write_bdfs_dict(bdf_files, self.setree, size, is_double, is_long_ids)


    def _write_tables_file(self, bdf_files: Any, size: int=8, is_double: bool=False,
                           is_long_ids: Optional[bool]=None) -> None:
        """Writes the TABLEx cards sorted by ID"""
        if self.tables or self.tables_d or self.tables_m or self.tables_sdamping:
            write_bdfs_dict(bdf_files, self.tables, size, is_double, is_long_ids)
            write_bdfs_dict(bdf_files, self.tables_d, size, is_double, is_long_ids)
            write_bdfs_dict(bdf_files, self.tables_m, size, is_double, is_long_ids)
            write_bdfs_dict(bdf_files, self.tables_sdamping, size, is_double, is_long_ids)

        if self.random_tables:
            write_bdfs_dict(bdf_files, self.random_tables, size, is_double, is_long_ids)

    def _write_thermal_file(self, bdf_files: Any, size: int=8, is_double: bool=False,
                            is_long_ids: Optional[bool]=None) -> None:
        """Writes the thermal cards"""
        # PHBDY
        is_thermal = (self.phbdys or self.convection_properties or self.bcs or
                      self.views or self.view3ds or self.radset or self.radcavs)
        if is_thermal:
            write_bdfs_dict(bdf_files, self.phbdys, size, is_double, is_long_ids)

            #for unused_key, prop in sorted(self.thermal_properties.items()):
            #    bdf_file.write(str(prop))
            write_bdfs_dict(bdf_files, self.convection_properties, size, is_double, is_long_ids)

            # BCs
            write_bdfs_dict_list(bdf_files, self.bcs, size, is_double, is_long_ids)

            write_bdfs_dict(bdf_files, self.views, size, is_double, is_long_ids)
            write_bdfs_dict(bdf_files, self.view3ds, size, is_double, is_long_ids)
            if self.radset:
                bdf_files[self.radset.ifile].write(self.radset.write_card(size, is_double))
            write_bdfs_dict(bdf_files, self.radcavs, size, is_double, is_long_ids)


    def _write_thermal_materials_file(self, bdf_files: Any, size: int=8, is_double: bool=False,
                                      is_long_ids: Optional[bool]=None) -> None:
        """Writes the thermal materials in a sorted order"""
        if self.thermal_materials:
            write_bdfs_dict(bdf_files, self.thermal_materials, size, is_double, is_long_ids)


def _get_ifiles_dict(cards_dict):
    """gets the ids for a dictionary by file number"""
    assert isinstance(cards_dict, dict), cards_dict
    ifiles_dict = defaultdict(list)
    for unused_id, card in sorted(cards_dict.items()):
        ifiles_dict[card.ifile].append(card)
    return ifiles_dict

def write_bdf_dict_ids(bdf_file, cards, ids, size, is_double, is_long_ids):
    """writes a dictionary by ifile"""
    assert isinstance(cards, dict), cards
    assert isinstance(cards, (list, tuple, np.ndarray)), ids
    if bdf_file is None:
        return
    if is_long_ids:
        for idi in ids:
            bdf_file.write(cards[idi].write_card_16(is_double))
    else:
        for idi in ids:
            bdf_file.write(cards[idi].write_card(size, is_double))

def write_bdfs_dict(bdf_files, cards, size, is_double, is_long_ids):
    """writes a dictionary by ifile"""
    assert isinstance(cards, dict), cards
    ifiles_dict = _get_ifiles_dict(cards)
    for file_id, file_cards in ifiles_dict.items():
        bdf_file = bdf_files[file_id]
        _write_bdf_dict_cards(bdf_file, file_cards, size, is_double, is_long_ids)

def _write_bdf_dict_cards(bdf_file, cards, size, is_double, is_long_ids):
    """writes a dictionary"""
    if bdf_file is None:
        return
    if is_long_ids:
        for card in cards:
            bdf_file.write(card.write_card_16(is_double))
    else:
        for card in cards:
            bdf_file.write(card.write_card(size, is_double))

def _get_ifiles_dict_list(cards):
    """gets the ids for a dictionary of lists by file number"""
    assert isinstance(cards, dict), cards
    ifiles_dict_list = defaultdict(list)
    for (unused_id, cardsi) in sorted(cards.items()):
        assert isinstance(cardsi, list), cardsi
        for card in cardsi:
            ifiles_dict_list[card.ifile].append(card)
    return ifiles_dict_list

def write_bdfs_dict_list(bdf_files, cards, size, is_double, is_long_ids):
    """writes a dictionary of lists by ifile"""
    ifiles_dict_list = _get_ifiles_dict_list(cards)
    for file_id, file_cards in ifiles_dict_list.items():
        bdf_file = bdf_files[file_id]
        _write_bdf_dict_cards(bdf_file, file_cards, size, is_double, is_long_ids)

def write_bdfs_list(bdf_files, cards, size, is_double, is_long_ids):
    """writes a list by ifile"""
    assert isinstance(cards, list), cards
    if is_long_ids:
        for card in cards:
            bdf_files[card.ifile].write_card_16(size, is_double)
    else:
        for card in cards:
            bdf_files[card.ifile].write_card(size, is_double)

def _map_filenames_to_ifile_filname_dict(out_filenames, active_filenames):
    """
    Converts a old_filename->new_filename dict to a
    ifile->new_filename dict.
    """
    #print('active_filenames = %s' % active_filenames)
    active_filenames_abspath = [os.path.abspath(path) for path in active_filenames]
    ifile_out_filenames = {}
    unused_out_filename0 = None
    for filename, new_filename in out_filenames.items():
        assert isinstance(filename, str), 'filename=%r' % filename
        #print('filename = %r' % filename)
        abs_filename = os.path.abspath(filename)
        #print('abs_filename = %r' % abs_filename)
        ifile = active_filenames_abspath.index(abs_filename)
        #print('ifile = %r' % ifile)
        #print('new_filename = %r' % new_filename)
        ifile_out_filenames[ifile] = new_filename
    return ifile_out_filenames

def _open_bdf_files(ifile_out_filenames, active_filenames, encoding):
    """opens N bdf files"""
    bdf_files = {i : None for i in range(len(active_filenames))}
    for ifile, out_filename in ifile_out_filenames.items():
        if hasattr(out_filename, 'read') and hasattr(out_filename, 'write'):
            bdf_file = out_filename
        else:
            bdf_file = open(out_filename, 'w', encoding=encoding)
        bdf_files[ifile] = bdf_file
    bdf_file0 = bdf_files[0]
    return bdf_files, bdf_file0

def write_xpoints_file(bdf_files, cardtype, points, comment=''):
    """writes SPOINTs/EPOINTs"""
    assert isinstance(points, dict), points
    for point_id, point in points.items():
        bdf_files[point.ifile].write(point.write_card())
