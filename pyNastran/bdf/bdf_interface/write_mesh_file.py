# coding: utf-8
"""
This file defines:
  - WriteMesh

"""
from __future__ import annotations
import os
import sys
import warnings
from pathlib import PurePath
from io import IOBase
from collections import defaultdict
from typing import TextIO, Optional, Any, cast, TYPE_CHECKING

import numpy as np

from pyNastran.utils import PathLike
from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.field_writer_16 import print_card_16
from pyNastran.bdf.bdf_interface.attributes import BDFAttributes
#from pyNastran.bdf.bdf_interface.write_mesh_utils import (
#    find_aero_location, write_dict, write_list, get_properties_by_element_type)
#from pyNastran.bdf.bdf_interface.utils import sorteddict
from pyNastran.bdf.write_path import write_include

#try:
#    from natsort import natsorted
#except ModuleNotFoundError:
#    natsorted = sorted

if TYPE_CHECKING:  # pragma: no cover
    from io import StringIO
    from cpylog import SimpleLogger
    from pyNastran.utils import PathLike
    from pyNastran.bdf.bdf import BDF, DESVAR
    TextFile = StringIO | TextIO


class TrashWriter:
    def __init__(self, log, filename: PathLike=''):
        self.log = log
        self.filename = str(filename)

    def write(self, msg):
        if self.filename:
            self.log.warning(f'trying to write to {self.filename}\n{str(msg)}')

    def close(self):
        pass


class WriteMesh(BDFAttributes):
    """
    Defines methods for writing cards

    Major methods:
      - model.write_bdf(...)
      - model.write_bdfs(...)
      - model.echo_bdf(...)
      - model.auto_reject_bdf(...)

    """
    def __init__(self):
        """creates methods for writing cards"""
        BDFAttributes.__init__(self)
        self._auto_reject = True
        self.cards_to_read = set()

    def get_encoding(self, encoding: Optional[str]=None) -> str:
        """gets the TextFile encoding"""
        return self.writer.get_encoding(encoding)

    def _reorganize_sets(self) -> None:
        """sorts the sets because ..."""
        for set_id, set_obj in self.sets.items():
            assert isinstance(set_obj.ids, list), set_obj.get_stats()
            set_obj.ids.sort()

    def apply_wtmass(self) -> None:
        #self._reorganize_sets()
        if 'WTMASS' not in self.params:
            return
        wtmass = self.params['WTMASS'].values[0]
        assert isinstance(wtmass, float), wtmass
        if wtmass == 1.0:
            return
        self.params['WTMASS'].values[0] = 1.0

        is_error = False
        for sid, nsms in self.nsms.items():
            for nsm in nsms:
                nsm.value *= wtmass

        for mid, mat in self.materials.items():
            try:
                mat.rho *= wtmass
            except:
                print(mat.get_stats())
                is_error = True

        for eid, elem in self.masses.items():
            try:
                if elem.type == 'CONM2':
                    elem.mass *= wtmass
                    elem.I *= wtmass
                else:
                    raise NotImplementedError(elem.get_stats())
            except:
                print(elem.get_stats())

        skip_properties = {'PSOLID'}
        for pid, prop in self.properties.items():
            prop_type = prop.type
            if prop_type in skip_properties:
                continue

            if prop_type == 'PBUSH':
                if prop.mass is not None:
                    prop.mass *= wtmass
            else:
                try:
                    prop.nsm *= wtmass
                except:
                    print(prop.get_stats())
                    is_error = True
        if is_error:
            raise RuntimeError('stopping...')

    def write_bdf(self, out_filename: Optional[PathLike | StringIO]=None,
                  encoding: Optional[str]=None,
                  size: int=8,
                  nodes_size: Optional[int]=None,
                  elements_size: Optional[int]=None,
                  loads_size: Optional[int]=None,
                  #table_size: Optional[int]=None,
                  flfact_size: int=0,
                  is_double: bool=False,
                  is_csv: bool=False,
                  sort_cards: bool=True,
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
        is_csv : bool; default=False
            False : write in standard format
            True ; write in CSV format (currently very limited)
        sort_cards : bool; default=True
            sort the nodes, elements, ... to make finding things easier
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

        .. note:: If you want to drop the executive & case control
                  decks, set model.punch = False

        """
        is_long_ids, size = self.writer.get_long_ids(size)

        out_filename, size = _output_helper(
            out_filename, interspersed, size, is_double, self.log)
        encoding = self.writer.get_encoding(encoding)
        #assert encoding.lower() in ['ascii', 'latin1', 'utf8'], encoding

        has_read_write = hasattr(out_filename, 'read') and hasattr(out_filename, 'write')
        if has_read_write:
            bdf_file = out_filename
        else:
            self.log.debug(f'---starting BDF.write_bdf of {out_filename}---')
            assert isinstance(encoding, str), encoding
            bdf_file = open(out_filename, 'w', encoding=encoding)

        writer = self.writer
        writer.write_header(
            bdf_file, encoding, write_header=write_header)
        #self.apply_wtmass()

        if self.superelement_models:
            bdf_file.write('$' + '*'*80+'\n')
            for superelement_tuple, superelement in self.superelement_models.items():
                if isinstance(superelement_tuple, int):
                    superelement_id = superelement_tuple
                    bdf_file.write(f'BEGIN SUPER={superelement_id}\n')
                else:
                    word, value, label = superelement_tuple
                    if label:
                        bdf_file.write(f'BEGIN {word}={value:d} LABEL={label}\n')
                    else:
                        bdf_file.write(f'BEGIN {word}={value:d}\n')
                superelement.write_bdf(out_filename=bdf_file, encoding=encoding,
                                       size=size, is_double=is_double,
                                       interspersed=interspersed, enddata=False,
                                       sort_cards=sort_cards,
                                       write_header=False, close=False)
                bdf_file.write('$' + '*'*80+'\n')
            bdf_file.write('BEGIN BULK\n')
        table_size = None
        writer.write_bulk_data(
            bdf_file, size=size, is_double=is_double,
            interspersed=interspersed,
            enddata=enddata, close=close,
            nodes_size=nodes_size, elements_size=elements_size, loads_size=loads_size,
            table_size=table_size,
            flfact_size=flfact_size,
            sort_cards=sort_cards,
            is_long_ids=is_long_ids,
            is_csv=is_csv)

    def write_bdfs(self, out_files_map: dict[str, str],
                   relative_dirname: PathLike='',
                   encoding: Optional[str]=None,
                   size: int=8, is_double: bool=False,
                   enddata: Optional[bool]=None, close: bool=True,
                   is_windows: Optional[bool]=None) -> None:
        """
        Writes the BDF.

        Parameters
        ----------
        out_files_map : dict[source_bdf, out_bdf]
            source_bdf : str
                the name of the original bdf
            out_bdf : str
                the name of the output bdf
        relative_dirname : str; default=''
            A relative path to reference INCLUDEs.
            ''   : relative to the first bdf in out_files_map
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

        out_files_map = {}
        out_files_map[fem.active_filenames[0]] = bdf_filename[:-4] + "_NEW" + bdf_filename[-4:]
        for ifile, include_filenames in model.include_filenames.items():
            for include_filename in include_filenames:
                base, ext = os.path.splitext(include_filename)
                new_filename = base + "_NEW" + ext
                out_files_map[include_filename] = new_filename

        """
        log = self.log
        assert isinstance(out_files_map, dict), out_files_map
        for key, value in out_files_map.items():
            assert isinstance(key, PathLike), type(key)
            assert isinstance(value, PathLike), type(value)
            # assert isinstance(key, PathLike), key
            # assert isinstance(value, PathLike), value
        #is_long_ids = False

        is_long_ids, size = self.writer.get_long_ids(size)

        ifile_out_filenames = _map_filenames_to_ifile_filname_dict(
            out_files_map, self.active_filenames)
        if len(ifile_out_filenames) == 0:
            log.warning(f'active_filenames = {self.active_filenames}')
            log.warning(f'ifile_out_filenames = {ifile_out_filenames}')
            raise RuntimeError(ifile_out_filenames)
            # bdf_filename = 'temp.bdf'
            # self.write_bdf(bdf_filename, size=size, is_double=is_double, encoding=encoding,
            #                close=close)
            return
        ifiles = list(sorted(ifile_out_filenames))
        #print(f'ifiles = {ifiles}')
        ifile0 = ifiles[0]
        #print('ifile_out_filenames =', ifile_out_filenames)

        out_filename0 = ifile_out_filenames[ifile0]
        #print("out_filename0 =", out_filename0)

        interspersed = False
        out_filename, size = _output_helper(
            out_filename0, interspersed, size, is_double, self.log)
        self.log.debug(f'---starting BDF.write_bdf of {out_filename}---')
        encoding = self.get_encoding(encoding)

        # class DevNull:
        #     def write(self, *_):
        #         pass
        # devnull = DevNull()

        bdf_files, bdf_file0 = _open_bdf_files(
            ifile_out_filenames, self.active_filenames, encoding, self.log)
        bdf_files[-1] = bdf_file0 # TrashWriter()

        writer = self.writer
        if bdf_file0 is not None:
            writer.write_header(bdf_file0, encoding)

        self.write_bdf_includes(out_files_map, bdf_files, relative_dirname=relative_dirname,
                                is_windows=is_windows)

        writer.write_params_file(bdf_files, size, is_double, is_long_ids=is_long_ids)
        writer.write_nodes_file(bdf_files, size, is_double, is_long_ids=is_long_ids)

        writer.write_elements_file(bdf_files, size, is_double, is_long_ids=is_long_ids)
        writer.write_properties_file(bdf_files, size, is_double, is_long_ids=is_long_ids)
        writer.write_materials_file(bdf_files, size, is_double, is_long_ids=is_long_ids)

        writer.write_masses_file(bdf_files, size, is_double, is_long_ids=is_long_ids)
        writer.write_rigid_elements_file(bdf_files, size, is_double, is_long_ids=is_long_ids)
        writer.write_aero_file(bdf_files, size, is_double, is_long_ids=is_long_ids)

        writer.write_common_file(bdf_files, size, is_double, is_long_ids=is_long_ids)
        if (enddata is None and 'ENDDATA' in self.card_count) or enddata:
            if bdf_file0:
                bdf_file0.write('ENDDATA\n')
        if close:
            for bdf_file in bdf_files.values():
                if bdf_file is not None:
                    bdf_file.close()
        del bdf_files

    def write_bdf_includes(self, out_files_map: dict[str, str],
                           bdf_files,
                           relative_dirname: str='',
                           is_windows: bool=True):
        """
        Writes the INCLUDE files

        Parameters
        ----------
        out_files_map : dict[fname] : fname2
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
            out_filename0 = list(out_files_map.keys())[0]
            relative_dirname = os.path.dirname(os.path.abspath(out_filename0))
            self.log.debug(f'relative_dirname = {relative_dirname}')

        self.log.debug('include_filenames:')
        for ifile, include_filenames in self.include_filenames.items():
            self.log.debug(f'ifile={ifile:d} {include_filenames}')
            assert len(include_filenames) > 0, include_filenames
            bdf_file = bdf_files[ifile]
            if bdf_file is None:
                continue
            #self.log.info('ifile=%s include_files=%s' % (ifile, include_filenames))
            for include_filename in include_filenames:
                assert len(include_filename) > 0, include_filename
                #print('***', include_filename, '***')

                mapped_include_filename = include_filename
                if include_filename in out_files_map:
                    mapped_include_filename = out_files_map[include_filename]

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


def _map_filenames_to_ifile_filname_dict(out_filenames: dict[str, str],
                                         active_filenames: list[str]) -> dict[int, str]:
    """
    Converts a old_filename->new_filename dict to a
    ifile->new_filename dict.
    """
    #print('active_filenames = %s' % active_filenames)
    active_filenames_abspath = [os.path.abspath(path) for path in active_filenames]
    ifile_out_filenames = {}
    unused_out_filename0 = None
    for filename, new_filename in out_filenames.items():
        assert isinstance(filename, PathLike), 'filename=%r' % filename
        assert isinstance(new_filename, PathLike), 'new_filename=%r' % new_filename
        #print('filename = %r' % filename)
        abs_filename = os.path.abspath(filename)
        #print('abs_filename = %r' % abs_filename)
        if abs_filename not in active_filenames_abspath:
            continue
        ifile = active_filenames_abspath.index(abs_filename)
        #print('ifile = %r' % ifile)
        #print('new_filename = %r' % new_filename)
        ifile_out_filenames[ifile] = new_filename
    return ifile_out_filenames


def _open_bdf_files(ifile_out_filenames, active_filenames, encoding, log):
    """opens N bdf files"""
    bdf_files = {i: TrashWriter(log, fname) for i, fname in enumerate(active_filenames)}
    for ifile, out_filename in ifile_out_filenames.items():
        log.debug(f'opening {str(out_filename)}')
        if hasattr(out_filename, 'read') and hasattr(out_filename, 'write'):
            bdf_file = out_filename
        else:
            bdf_file = open(out_filename, 'w', encoding=encoding)
        bdf_files[ifile] = bdf_file
    bdf_file0 = bdf_files[0]
    return bdf_files, bdf_file0


def _ifile(card) -> int:
    try:
        ifile = card.ifile
    except AttributeError:
        warnings.warn(f'cant find ifile in\n{str(card)}')
        ifile = -1
    return ifile

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
