# coding: utf-8
"""
Main BDF class.  Defines:
  - BDFInputPy

"""
import os
from collections import defaultdict
from itertools import count
from typing import List, Tuple, Optional, Union, Any, cast
from io import StringIO

import numpy as np
from cpylog import get_logger2
from pyNastran.nptyping import NDArrayN2int
from pyNastran.utils import print_bad_path, _filename

from pyNastran.bdf import BULK_DATA_CARDS, FLAGGED_CARDS
from pyNastran.bdf.errors import MissingDeckSections
from pyNastran.bdf.bdf_interface.utils import _parse_pynastran_header
from pyNastran.bdf.bdf_interface.include_file import get_include_filename

FILE_MANAGEMENT = (
    'ACQUIRE ', 'ASSIGN ', 'CONNECT ', 'DBCLEAN ', 'DBDICT ', 'DBDIR ',
    'DBFIX ', 'DBLOAD ', 'DBLOCATE ', 'DBSETDEL ', 'DBUNLOAD ',
    'DBUPDATE ', 'ENDJOB ', 'EXPAND ', 'INCLUDE ', 'INIT ', 'NASTRAN ',
    'PROJ ',
)
EXECUTIVE_CASE_SPACES = tuple(list(FILE_MANAGEMENT) + ['SOL ', 'SET ', 'SUBCASE '])


class BDFInputPy:
    """BDF reader class that only handles lines and not building cards or parsing cards"""
    def __init__(self, read_includes: bool, dumplines: bool,
                 encoding: str, nastran_format: str='msc',
                 consider_superelements: bool=True,
                 log: Any=None, debug: bool=False):
        """
        BDF reader class that only handles lines and not building cards or parsing cards

        Parameters
        ----------
        read_includes : bool
            should include files be read
        dumplines : bool
            Writes 'pyNastran_dump.bdf' up to some failed line index
        encoding : str
            the character encoding (e.g., utf8, latin1, cp1252)
        nastran_format : str; default='msc'
            'zona' has a special read method
            {msc, nx, zona}
        consider_superelements : bool; default=True
            parse 'begin super=2'
        log : logger(); default=None
            a logger for printing INCLUDE files that are loadaed
        debug : bool; default=False
            used when testing; for the logger

        """
        self.dumplines = dumplines
        self.encoding = encoding
        self.nastran_format = nastran_format
        self.include_dir = ''

        self.include_lines = defaultdict(list)
        self.read_includes = read_includes
        self.active_filenames = []
        self.active_filename = None

        self.consider_superelements = consider_superelements
        self.debug = debug
        self.log = get_logger2(log, debug)

    def get_lines(self, bdf_filename: Union[str, StringIO],
                  punch: Optional[bool]=False,
                  make_ilines: bool=True) -> List[str]:
        """
        Opens the bdf and extracts the lines by group

        Parameters
        ----------
        bdf_filename : str
            the main bdf_filename
        punch : bool / None; default=False
            is this a punch file
            None : guess
            True : no executive/case control decks
            False : executive/case control decks exist
        make_ilines : bool; default=True
            flag for bulk_data_ilines

        Returns
        -------
        system_lines : List[str]
            the system control lines (typically empty; used for alters)
        executive_control_lines : List[str]
            the executive control lines (stores SOL 101)
        case_control_lines : List[str]
            the case control lines (stores subcases)
        bulk_data_lines : List[str]
            the bulk data lines (stores geometry, boundary conditions, loads, etc.)
        bulk_data_ilines : None / (nlines, 2) int ndarray
            if make_ilines = True:
                the [ifile, iline] pair for each line in the file
            if make_ilines = False:
                 ilines = None

        """
        main_lines = self.get_main_lines(bdf_filename)
        all_lines, ilines = self.lines_to_deck_lines(main_lines, make_ilines=make_ilines)

        out = _lines_to_decks(all_lines, ilines, punch, self.log,
                              keep_enddata=True,
                              consider_superelements=self.consider_superelements)
        (
            system_lines, executive_control_lines, case_control_lines,
            bulk_data_lines, bulk_data_ilines,
            superelement_lines, superelement_ilines) = out
        if self.nastran_format in ['msc', 'nx', 'nasa95']:
            pass
        elif self.nastran_format == 'zona':
            bulk_data_lines, bulk_data_ilines, system_lines = self._get_lines_zona(
                system_lines, bulk_data_lines, bulk_data_ilines, punch)
        else:
            msg = 'nastran_format=%r and must be msc, nx, or zona' % self.nastran_format
            raise NotImplementedError(msg)
        return (system_lines, executive_control_lines, case_control_lines,
                bulk_data_lines, bulk_data_ilines,
                superelement_lines, superelement_ilines)

    def _get_lines_zona(self, system_lines: List[str], bulk_data_lines: List[str],
                        bulk_data_ilines: NDArrayN2int,
                        punch: bool) -> Tuple[List[str], NDArrayN2int, List[str]]:
        """load and update the lines for ZONA"""
        system_lines2 = []
        for system_line in system_lines:
            if system_line.upper().startswith('ASSIGN'):
                split_system = system_line.split(',')
                header = split_system[0]
                header_upper = header.upper()
                if header_upper.startswith('ASSIGN FEM'):
                    unused_fem, filename = header.split('=')
                    filename = filename.strip('"\'')
                    self.log.debug('reading %s' % filename)
                    if filename.lower().endswith('.f06'):
                        filename = os.path.splitext(filename)[0] + '.bdf'
                    if not filename.endswith('.bdf'):
                        raise RuntimeError('filename must end in bdf; %s' % filename)


                    _main_lines = self.get_main_lines(filename)
                    make_ilines = bulk_data_ilines is not None
                    _all_lines, _ilines = self.lines_to_deck_lines(
                        _main_lines, make_ilines=make_ilines)
                    _out = _lines_to_decks(_all_lines, _ilines, punch, self.log,
                                           keep_enddata=False,
                                           consider_superelements=self.consider_superelements)
                    (
                        _system_lines, _executive_control_lines, _case_control_lines,
                        bulk_data_lines2, bulk_data_ilines2,
                        _superelement_lines, _superelement_ilines,
                    ) = _out
                    bulk_data_lines = bulk_data_lines2 + bulk_data_lines
                    #print("bulk_data_ilines2 =", bulk_data_ilines2, bulk_data_ilines2.shape)
                    #print("bulk_data_ilines =", bulk_data_ilines, bulk_data_ilines.shape)
                    bulk_data_ilines = np.vstack([bulk_data_ilines2, bulk_data_ilines])
                    continue
                elif header_upper.startswith('ASSIGN MATRIX'):
                    pass
                elif header_upper.startswith('ASSIGN OUTPUT4'):
                    pass
                else:  # pragma: no cover
                    raise NotImplementedError(system_line)
            system_lines2.append(system_line)
        system_lines = system_lines
        return bulk_data_lines, bulk_data_ilines, system_lines

    def get_main_lines(self, bdf_filename: Union[str, StringIO]) -> List[str]:
        """
        Opens the bdf and extracts the lines

        Parameters
        ----------
        bdf_filename : str / StringIO
            the main bdf_filename

        Returns
        -------
        lines : List[str]
            all the lines packed into a single line stream

        """
        #print('bdf_filename_main =', bdf_filename)
        if hasattr(bdf_filename, 'read') and hasattr(bdf_filename, 'write'):
            bdf_filename = cast(StringIO, bdf_filename)
            lines = bdf_filename.readlines()
            if len(lines) == 0:
                raise RuntimeError('lines in %s is empty' % bdf_filename)
            return lines

        bdf_filename = cast(str, bdf_filename)
        self.bdf_filename = bdf_filename
        # the directory of the 1st BDF (include BDFs are relative to this one)
        self.include_dir = os.path.dirname(os.path.abspath(bdf_filename))

        with self._open_file(bdf_filename, basename=True) as bdf_file:
            try:
                lines = bdf_file.readlines()
            except UnicodeDecodeError:
                _show_bad_file(self, bdf_filename, encoding=self.encoding)
        return lines

    def lines_to_deck_lines(self, lines: List[str], make_ilines: bool=True) -> Tuple[List[str], int]:
        """
        Merges the includes into the main deck.

        Parameters
        ----------
        lines : List[str]
            the lines from the main BDF
        make_ilines : bool; default=True
            flag for ilines

        Returns
        -------
        active_lines : List[str]
            all the active lines in the deck
        ilines : (nlines, 2) int ndarray
            if make_ilines = True:
                the [ifile, iline] pair for each line in the file
            if make_ilines = False:
                 ilines = None

        """
        nlines = len(lines)
        #bdf_filenames = [self.bdf_filename]

        ilines = None
        if make_ilines:
            ilines = _make_ilines(nlines, ifile=0)

        i = 0
        ifile = 1
        while i < nlines:
            try:
                line = lines[i].rstrip('\r\n\t')
            except IndexError:
                break
            uline = line.upper()
            if uline.startswith('INCLUDE'):
                j, include_lines = self._get_include_lines(lines, line, i, nlines)
                bdf_filename2 = get_include_filename(include_lines, include_dir=self.include_dir)
                #bdf_filenames.append(bdf_filename2)
                jfile = ilines[i, 0]
                # these are the lines associated with the 1st/2nd include file found
                self.include_lines[jfile].append((include_lines, bdf_filename2))

                if self.read_includes:
                    lines, nlines, ilines = self._update_include(
                        lines, nlines, ilines,
                        include_lines, bdf_filename2, i, j, ifile, make_ilines=make_ilines)
                    ifile += 1
                else:
                    lines = lines[:i] + lines[j:]
                    if make_ilines:
                        ilines = np.vstack([
                            ilines[:i, :],
                            ilines[j:, :],
                        ])
                        #assert len(ilines[:, 1]) == len(np.unique(ilines[:, 1]))
                    #self.include_lines[ifile].append(include_lines)
                    #self.reject_lines.append(write_include(bdf_filename2))
            i += 1

        if self.dumplines:
            self._dump_file('pyNastran_dump.bdf', lines, i)

        #print(bdf_filenames)
        #if make_ilines:
            #nilines = ilines.shape[0]
            #assert nlines == ilines.shape[0], 'nlines=%s nilines=%s' % (nlines, nilines)
        return lines, ilines

    def _update_include(self, lines: List[str], nlines: int, ilines,
                        include_lines: List[str], bdf_filename2: str, i: int, j: int, ifile: int,
                        make_ilines: bool=False):
        """incorporates an include file into the lines"""
        try:
            self._open_file_checks(bdf_filename2)
        except IOError:
            crash_name = 'pyNastran_crash.bdf'
            self._dump_file(crash_name, lines, j)
            msg = 'There was an invalid filename found while parsing.\n'
            msg += 'Check the end of %r\n' % crash_name
            msg += 'bdf_filename2 = %r\n' % bdf_filename2
            msg += 'abs_filename2 = %r\n' % os.path.abspath(bdf_filename2)
            msg += 'include_lines = %s' % include_lines
            #msg += 'len(bdf_filename2) = %s' % len(bdf_filename2)
            print(msg)
            raise
            #raise IOError(msg)

        read_again = False
        with self._open_file(bdf_filename2, basename=False) as bdf_file:
            #print('bdf_file.name = %s' % bdf_file.name)
            try:
                lines2 = bdf_file.readlines()
            except UnicodeDecodeError:
                #try:
                bdf_file.seek(0)
                try:
                    encoding2 = _check_pynastran_encoding(bdf_filename2, encoding=self.encoding)
                except UnicodeDecodeError:
                    encoding2 = self.encoding

                #print('***encoding=%s encoding2=%s' % (self.encoding, encoding2))
                if self.encoding != encoding2:
                    read_again = True
                else:
                    msg = (
                        'Invalid Encoding: encoding=%r.  Fix it by:\n'
                        '  1.  try a different encoding (e.g., latin1, cp1252, utf8)\n'
                        "  2.  call read_bdf(...) with `encoding`'\n"
                        "  3.  Add '$ pyNastran : encoding=latin1"
                        ' (or other encoding) to the top of the main/INCLUDE file\n' % (
                            self.encoding))
                    raise RuntimeError(msg)

        if read_again:
            self.active_filenames.pop()
            with self._open_file(bdf_filename2, basename=False, encoding=encoding2) as bdf_file:
                #print('bdf_file.name = %s' % bdf_file.name)
                try:
                    lines2 = bdf_file.readlines()
                except UnicodeDecodeError:
                    msg = (
                        'Incorrect Encoding: encoding=%r.  Fix it by:\n'
                        '  1.  try a different encoding (e.g., latin1, cp1252, utf8)\n'
                        "  2.  call read_bdf(...) with `encoding`'\n"
                        "  3.  Add '$ pyNastran : encoding=latin1"
                        ' (or other encoding) to the top of the main/INCLUDE file\n' % encoding2)
                    raise RuntimeError(msg)

        #print('lines2 = %s' % lines2)

        #line2 = lines[j].split('$')
        #if not line2[0].isalpha():
            #print('** %s' % line2)

        include_comment = '\n$ INCLUDE processed:  %s\n' % bdf_filename2
        #for line in lines2:
            #print("  ?%s" % line.rstrip())

            #for ii, line in enumerate(lines):
                #print('  %i %r' % (ii, line))

        #for ii in range(i):
            #print('i=%i %r' % (ii, lines[ii]))
        #print('---------')
        #for jj in range(j, len(lines)):
            #print('j=%i %r' % (jj, lines[jj]))
        #print('include_comment = %r' % include_comment)

        nlines2 = len(lines2)
        if make_ilines:
            ilines2 = _make_ilines(nlines2, ifile)
            #n_ilines = ilines.shape[0]
            #print(ilines[j:, :])
            #assert len(lines[:i]) == ilines[:i+1, :].shape[0] - ifile, 'A: nlines=%s nilines=%s' % (len(lines[:i]), ilines[:i+1, :].shape[0])
            #assert len(lines[j:]) == ilines[j:, :].shape[0],           'B: nlines=%s nilines=%s' % (len(lines[j:]), ilines[j:, :].shape[0])
            #assert len(lines2) == ilines2.shape[0],                    'C: nlines=%s nilines=%s' % (len(lines2),    ilines2.shape[0])
            #assert nlines == ilines.shape[0], 'B: nlines=%s nilines=%s' % (nlines, nilines)
            #assert nlines == ilines.shape[0], 'C: nlines=%s nilines=%s' % (nlines, nilines)

            #print(nlines-ifile+1, nlines2)
            #print(
                #len(lines[:i]),
                #len([include_comment]),
                #len(lines2),
                #len(lines[j:]),
            #)
            #print(
                #ilines[:i+1, :].shape[0],
                #ilines2.shape[0],
                #ilines[j:, :].shape[0],
            #)
            #print('ilines:\n%s' % ilines)
            #print('ilines2:\n%s' % ilines2)
            #print('lines2:\n%s' % ''.join(lines2))
            ilines = np.vstack([
                ilines[:i+1, :],
                ilines2,
                ilines[j:, :],
            ])
            #dij = j - i

        nlines += nlines2
        lines = lines[:i] + [include_comment] + lines2 + lines[j:]
        #print('*lines:\n%s' % ''.join(lines))
        #for ifile_iline, line in zip(ilines, lines):
            #print(ifile_iline, line.rstrip())
        #if make_ilines:
            #n_ilines = ilines.shape[0]
            #ncompare = n_ilines - dij
            #assert n_ilines == n_ilines, 'nlines=%s dij=%s n_ilines=%s' % (n_ilines, dij, n_ilines)
        #for line in lines:
            #print("  *%s" % line.rstrip())
        return lines, nlines, ilines

    def _get_include_lines(self, lines: List[str], line: str,
                           i: int, nlines: int) -> Tuple[int, List[str]]:
        """
        gets the lines for the include file

        INCLUDE 'Satellite_V02_INCLUDE:Satellite_V02_Panneau_Externe.dat'
        INCLUDE '../../BULK/COORDS/satellite_V02_Coord.blk'
        """
        j = i + 1
        line_base = line.split('$')[0]
        include_lines = [line_base.strip()]
        if "'" not in line_base:
            pass
        else:
            #print('----------------------')

            line_base = line_base[8:].strip()
            if line_base.startswith("'") and line_base.endswith("'"):
                pass
            else:
                while not line.split('$')[0].endswith("'") and j < nlines:
                    #print('j=%s nlines=%s less?=%s'  % (j, nlines, j < nlines))
                    try:
                        line = lines[j].split('$')[0].strip()
                    except IndexError:
                        #print('bdf_filename=%r' % bdf_filename)
                        crash_name = 'pyNastran_crash.bdf'
                        self._dump_file(crash_name, lines, i+1)
                        msg = 'There was an invalid filename found while parsing (index).\n'
                        msg += 'Check the end of %r\n' % crash_name
                        #msg += 'bdf_filename2 = %r\n' % bdf_filename
                        msg += 'include_lines = %s' % include_lines
                        raise IndexError(msg)
                    #print('endswith_quote=%s; %r' % (
                        #line.split('$')[0].strip().endswith(""), line.strip()))
                    include_lines.append(line.strip())
                    j += 1
                #print('j=%s nlines=%s less?=%s'  % (j, nlines, j < nlines))

                #print('*** %s' % line)
                #bdf_filename2 = line[7:].strip(" '")
                #include_lines = [line] + lines[i+1:j]
        #print(include_lines)
        return j, include_lines

    def _dump_file(self, bdf_dump_filename: str,
                   lines: List[str],
                   i: int) -> None:
        """
        Writes a BDF up to some failed line index

        Parameters
        ----------
        bdf_dump_filename : str
            the bdf filename to dump
        lines : List[str]
            the entire list of lines
        i : int
            the last index to write

        """
        with open(_filename(bdf_dump_filename),
                  'w', encoding=self.encoding) as crash_file:
            for line in lines[:i]:
                crash_file.write(line)

    def _open_file_checks(self, bdf_filename: str, basename: bool=False) -> None:
        """
        Verifies that the BDF about to be opened:
           1.  Exists
           2.  Is Unique
           3.  Isn't an OP2
           4.  Is a File

        Parameters
        ----------
        bdf_filename : str
            the bdf filename to open
        basename : bool; default=False
            only take the basename of the bdf

        """
        if basename:
            bdf_filename_inc = os.path.join(self.include_dir, os.path.basename(bdf_filename))
        else:
            bdf_filename_inc = os.path.join(self.include_dir, bdf_filename)

        if not os.path.exists(_filename(bdf_filename_inc)):
            msg = 'No such bdf_filename: %r\n' % bdf_filename_inc
            msg += 'cwd: %r\n' % os.getcwd()
            msg += 'include_dir: %r\n' % self.include_dir
            msg += print_bad_path(bdf_filename_inc)
            self.log.error(msg)
            raise IOError(msg)
        elif bdf_filename_inc.endswith('.op2'):
            msg = 'Invalid filetype: bdf_filename=%r' % bdf_filename_inc
            self.log.error(msg)
            raise IOError(msg)
        bdf_filename = bdf_filename_inc

        if bdf_filename in self.active_filenames:
            msg = 'bdf_filename=%s is already active.\nactive_filenames=%s' \
                % (bdf_filename, self.active_filenames)
            self.log.error(msg)
            raise RuntimeError(msg)
        elif os.path.isdir(_filename(bdf_filename)):
            current_filename = self.active_filename if len(self.active_filenames) > 0 else 'None'
            msg = 'Found a directory: bdf_filename=%r\ncurrent_file=%s' % (
                bdf_filename_inc, current_filename)
            self.log.error(msg)
            raise IOError(msg)
        elif not os.path.isfile(_filename(bdf_filename)):
            msg = 'Not a file: bdf_filename=%r' % bdf_filename
            self.log.error(msg)
            raise IOError(msg)

    def _open_file(self, bdf_filename: Union[str, StringIO],
                   basename: bool=False, check: bool=True, encoding: Optional[str]=None) -> Any:
        """
        Opens a new bdf_filename with the proper encoding and include directory

        Parameters
        ----------
        bdf_filename : str
            the filename to open
        basename : bool (default=False)
            should the basename of bdf_filename be appended to the include directory
        check : bool; default=True
            you can disable the checks

        Returns
        -------
        bdf_file : file
            a file object

        """
        if encoding is None:
            encoding = self.encoding
        if basename:
            bdf_filename_inc = os.path.join(self.include_dir, os.path.basename(bdf_filename))
        else:
            bdf_filename_inc = os.path.join(self.include_dir, bdf_filename)

        self._validate_open_file(bdf_filename, bdf_filename_inc, check)


        self.log.debug('opening %r' % bdf_filename_inc)
        self.active_filenames.append(bdf_filename_inc)

        #print('ENCODING - _open_file=%r' % self.encoding)
        #self._check_pynastran_header(lines)
        bdf_file = open(_filename(bdf_filename_inc), 'r', encoding=encoding)
        return bdf_file

    def _validate_open_file(self, bdf_filename: Union[str, StringIO],
                            bdf_filename_inc: str,
                            check: bool=True) -> None:
        """
        checks that the file doesn't have obvious errors
         - hasn't been used
         - not a directory
         - is a file

        Parameters
        ----------
        bdf_filename : str
           the current bdf filename
        bdf_filename_inc : str
           the next bdf filename
        check : bool; default=True
            you can disable the checks

        Raises
        ------
        RuntimeError : file is active
        IOError : Invalid file type

        """
        if check:
            if not os.path.exists(_filename(bdf_filename_inc)):
                msg = 'No such bdf_filename: %r\n' % bdf_filename_inc
                msg += 'cwd: %r\n' % os.getcwd()
                msg += 'include_dir: %r\n' % self.include_dir
                msg += print_bad_path(bdf_filename_inc)
                raise IOError(msg)
            elif bdf_filename_inc.endswith('.op2'):
                raise IOError('Invalid filetype: bdf_filename=%r' % bdf_filename_inc)

            bdf_filename = bdf_filename_inc
            if bdf_filename in self.active_filenames:
                msg = 'bdf_filename=%s is already active.\nactive_filenames=%s' \
                    % (bdf_filename, self.active_filenames)
                raise RuntimeError(msg)
            elif os.path.isdir(_filename(bdf_filename)):
                current_fname = self.active_filename if len(self.active_filenames) > 0 else 'None'
                raise IOError('Found a directory: bdf_filename=%r\ncurrent_file=%s' % (
                    bdf_filename_inc, current_fname))
            elif not os.path.isfile(_filename(bdf_filename)):
                raise IOError('Not a file: bdf_filename=%r' % bdf_filename)


def _is_bulk_data_line(text: str) -> bool:
    """
    Returns True if there is a Bulk Data Deck

    Parameters
    ----------
    text : str
        a line in the deck

    Returns
    -------
    is_bulk_line : bool
        is this a bulk data line

    """
    #
    # Stripping the data isn't ideal as a GRID card cannot have a leading space.
    #
    # We strip the data because we need to support:
    # '    SUPORT1 = 10'
    # if you use text[0:8].strip('* ').upper()
    # we get:
    # '    SUPO', which is not a case control card called 'SUPORT1'
    #
    text2 = text.split('$')[0].rstrip()
    #card_name = text2.strip().replace(' ', '')[0:8].split(',')[0].upper().rstrip('*')
    card_name = text2.strip()[0:8].split(',')[0].upper().rstrip('*')
    #card_name2 = text2.split(',')[0].upper()
    #print('card_name =%r' % card_name)
    #print('card_name2=%r' % card_name2)

    # bulk data cards
    if card_name in BULK_DATA_CARDS:
        # case control + bulk data cards
        if '=' in text2 and card_name in FLAGGED_CARDS or text2.startswith(' '):
            return False
        elif card_name == 'PARAM':
            # The PARAM card can have a comma or tab, but no equals sign.
            # If there is a PARAM card, we have to assume we're not in the
            # case control.
            return False
        return True
    return False


def _check_pynastran_encoding(bdf_filename: Union[str, StringIO], encoding: str) -> str:
    """updates the $pyNastran: key=value variables"""
    line = '$pyNastran: punch=False'
    #line_temp = u'é à è ê'.encode('utf8').decode('ascii')

    skip_keys = [
        'version', 'punch', 'nnodes', 'nelements', 'dumplines',
        'is_superelements', 'skip_cards', 'units']

    with open(bdf_filename, 'rb') as bdf_file:
        line = bdf_file.readline()
        line_str = line.decode('ascii')
        while '$' in line_str:
            #if not line.startswith('$'):
                #break

            key, value = _parse_pynastran_header(line_str)
            if not key:
                break

            # key/value are lowercase
            if key == 'encoding':
                encoding = value
                break
            elif key in skip_keys or 'skip ' in key:
                pass
            else:
                raise NotImplementedError(key)
            line = bdf_file.readline()
            line_str = line.decode('ascii')
    return encoding


IGNORE_COMMENTS = (
    '$EXECUTIVE CONTROL DECK',
    '$CASE CONTROL DECK',
    'NODES', 'SPOINTS', 'EPOINTS', 'ELEMENTS',
    'PARAMS', 'PROPERTIES', 'ELEMENTS_WITH_PROPERTIES',
    'ELEMENTS_WITH_NO_PROPERTIES (PID=0 and unanalyzed properties)',
    'UNASSOCIATED_PROPERTIES',
    'MATERIALS', 'THERMAL MATERIALS',
    'CONSTRAINTS', 'SPCs', 'MPCs', 'RIGID ELEMENTS',
    'LOADS', 'AERO', 'STATIC AERO', 'AERO CONTROL SURFACES',
    'FLUTTER', 'GUST', 'DYNAMIC', 'OPTIMIZATION',
    'COORDS', 'THERMAL', 'TABLES', 'RANDOM TABLES',
    'SETS', 'CONTACT', 'REJECTS', 'REJECT_CARDS', 'REJECT_LINES',
    'PROPERTIES_MASS', 'MASSES')


def _clean_comment(comment: str) -> Optional[str]:
    """
    Removes specific pyNastran comment lines so duplicate lines aren't
    created.

    Parameters
    ----------
    comment : str
        the comment to possibly remove

    Returns
    -------
    updated_comment : str
        the comment

    """
    if comment == '':
        pass
    elif comment in IGNORE_COMMENTS:
        comment = None
    elif 'pynastran' in comment.lower():
        csline = comment.lower().split('pynastran', 1)
        if csline[1].strip()[0] == ':':
            comment = None

    #if comment:
        #print(comment)
    return comment


def _lines_to_decks(lines: List[str],
                    ilines: NDArrayN2int,
                    punch: Optional[bool],
                    log: Any,
                    keep_enddata: bool=True,
                    consider_superelements: bool=False) -> Tuple[
                        List[str], List[str], List[str], List[str], NDArrayN2int,
                        List[str], List[str], List[str]]:
    """
    Splits the BDF lines into:
     - system lines
     - executive control deck
     - case control deck
     - bulk data deck

    Parameters
    ----------
    lines : List[str]
        all the active lines in the deck
    ilines : None / (nlines, 2) int ndarray
        None : the old behavior
        narray : the [iline, ifile] pair for each line in the file
    punch : bool / None
        None : guess
        True : starts from the bulk data deck
        False : read the entire deck
    keep_enddata : bool; default=True
        True : don't throw away the enddata card
        False : throw away the enddata card

    Returns
    -------
    system_lines : List[str]
        the system control lines (typically empty; used for alters)
    executive_control_lines : List[str]
        the executive control lines (stores SOL 101)
    case_control_lines : List[str]
        the case control lines (stores subcases)
    bulk_data_lines : List[str]
        the bulk data lines (stores geometry, boundary conditions, loads, etc.)
    bulk_data_ilines : None / (nlines, 2) int ndarray
        None : the old behavior
        narray : the [ifile, iline] pair for each line in the file
    superelement_lines : List[str]
        ???
    superelement_ilines : List[str]
        ???
    auxmodel_lines : List[str]
        ???

    """
    if punch: # True
        system_lines = []
        executive_control_lines = []
        case_control_lines = []
        bulk_data_lines = lines
        bulk_data_ilines = ilines
        superelement_lines = {}
        superelement_ilines = {}
        #auxmodel_lines = {}
        return (
            system_lines, executive_control_lines, case_control_lines,
            bulk_data_lines, bulk_data_ilines,
            superelement_lines, superelement_ilines)

    # typical deck
    out = _lines_to_decks_main(lines, ilines, log, punch=punch,
                               keep_enddata=keep_enddata,
                               consider_superelements=consider_superelements)
    (executive_control_lines, case_control_lines,
     bulk_data_lines, bulk_data_ilines,
     superelement_lines, superelement_ilines,
     auxmodel_lines, afpm_lines) = out


    # break out system commands
    system_lines, executive_control_lines = _break_system_lines(executive_control_lines)

    for super_id, _lines in superelement_lines.items():
        # cqrsee101b2.bdf
        if len(_lines) == 0:
            raise RuntimeError('lines in superelement %i is empty' % super_id)

        #assert len(_lines) == len(superelement_ilines[super_id]), 'superelement %i ilines is the wrong length' % (super_id)

    for auxmodel_id, _lines in auxmodel_lines.items():
        # C:\MSC.Software\MSC.Nastran2005r3\msc20055\nast\tpl\motion21.dat
        # C:\MSC.Software\MSC.Nastran2005r3\msc20055\nast\tpl\d200am1.dat
        # C:\MSC.Software\MSC.Nastran2005r3\msc20055\nast\tpl\d200am2.dat
        log.warning('skipping auxmodel=%i' % auxmodel_id)
        raise RuntimeError('lines in auxmodel %i is empty' % auxmodel_id)

    for afpm_id, _lines in afpm_lines.items():
        log.warning('skipping AFPM=%i' % afpm_id)
        raise RuntimeError('lines in AFPM %i is empty' % afpm_id)

    # clean comments
    system_lines = [_clean_comment(line) for line in system_lines
                    if _clean_comment(line) is not None]
    executive_control_lines = [_clean_comment(line) for line in executive_control_lines
                               if _clean_comment(line) is not None]
    case_control_lines = [_clean_comment(line) for line in case_control_lines
                          if _clean_comment(line) is not None]
    return (
        system_lines, executive_control_lines, case_control_lines,
        bulk_data_lines, bulk_data_ilines,
        superelement_lines, superelement_ilines)

def _lines_to_decks_main(lines: List[str],
                         ilines: Any, log: Any,
                         punch: Optional[bool]=False,
                         keep_enddata: bool=True,
                         consider_superelements: bool=False) -> Tuple[
                        List[str], List[str], List[str], List[str], NDArrayN2int,
                        List[str], List[str], List[str]]:
    """
    Splits the BDF lines into:
     - system lines
     - executive control deck
     - case control deck
     - bulk data deck

    Parameters
    ----------
    lines : List[str]
        all the active lines in the deck
    ilines : None / (nlines, 2) int ndarray
        None : the old behavior
        narray : the [iline, ifile] pair for each line in the file
    punch : bool / None; default=False
        None : guess
        True : punch file (skipped previously, so this can't be True)
        False : not a punch file
    keep_enddata : bool; default=True
        True : don't throw away the enddata card
        False : throw away the enddata card
    consider_superelements : bool; default=True
        parse 'begin super=2'

    Returns
    -------
    system_executive_control_lines : List[str]
        the system control lines (typically empty; used for alters)
        and the executive control lines (stores SOL 101)
    case_control_lines : List[str]
        the case control lines (stores subcases)
    bulk_data_lines : List[str]
        the bulk data lines (stores geometry, boundary conditions, loads, etc.)
    bulk_data_ilines : None / (nlines, 2) int ndarray
        None : the old behavior
        narray : the [ifile, iline] pair for each line in the file
    superelement_lines : List[str]
        ???
    superelement_ilines : List[str]
        ???
    auxmodel_lines : List[str]
        ???

    """
    make_ilines = ilines is not None
    guess_deck_sections = punch is None

    executive_control_lines = []
    case_control_lines = []
    bulk_data_lines = []
    superelement_lines = defaultdict(list)
    superelement_ilines = defaultdict(list)
    auxmodel_lines = defaultdict(list)
    afpm_lines = defaultdict(list)
    auxmodels_found = set()
    afpms_found = set()
    auxmodels_to_find = []
    afpms_to_find = []
    is_auxmodel = False
    is_afpm = False
    is_superelement = False
    is_auxmodel_active = False
    is_afpm_active = False
    auxmodel_id = None
    afpm_id = None
    #---------------------------------------------
    current_lines = executive_control_lines

    #flag_word = 'executive'
    flag = 1
    old_flags = []
    bulk_data_ilines = []
    if ilines is None:
        ilines = count()

    #print('guess_deck_sections =', guess_deck_sections, punch)
    for i, ifile_iline, line in zip(count(), ilines, lines):
        #print('%s %-8s %s' % (ifile_iline, flag_word, line.rstrip()))
        #print('%s %i %s' % (ifile_iline, flag, line.rstrip()))
        uline = line.split('$')[0].upper().strip()

        if guess_deck_sections and flag == 1 and uline.startswith('BEGIN'):
            section_name_map = {
                1 : 'executive control',
                2 : 'case control',
            }
            section_name = section_name_map[flag]

            if _is_begin_bulk(uline):
                #old_flags.append(flag)
                log.warning('currently in %s deck and skipping directly '
                            'to bulk data section' % section_name)
                flag = 3
                current_ilines = bulk_data_ilines
                current_lines = bulk_data_lines
                bulk_data_ilines = _bulk_data_lines_extract(
                    lines, ilines, bulk_data_lines, i,
                    make_ilines=make_ilines, keep_enddata=keep_enddata)
            else:
                raise RuntimeError('currently in %s deck and unexpectedly found the following '
                                   'line:\n%s' % (section_name, line))
            break

        if guess_deck_sections and flag in [1, 2] and _is_bulk_data_line(line):
            section_name_map = {
                1 : 'executive control',
                2 : 'case control',
            }
            section_name = section_name_map[flag]
            log.warning('currently in %s deck and skipping directly '
                        'to bulk data section\n%s' % (section_name, line))
            log.warning(line)
            flag = 3
            current_ilines = bulk_data_ilines
            current_lines = bulk_data_lines
            bulk_data_ilines = _bulk_data_lines_extract(
                lines, ilines, bulk_data_lines, i-1,
                make_ilines=make_ilines, keep_enddata=keep_enddata)
            break

        elif flag == 1:
            # I don't think we need to handle the comment because
            # this uses a startswith
            if line.upper().startswith('CEND'):
                # case control
                old_flags.append(flag)
                if flag != 1:
                    raise RuntimeError('expected a flag of 1 (executive control deck) '
                                       'when going to the case control deck')

                flag = 2
                #flag_word = 'case'
                current_lines = case_control_lines
            #print('executive: ', line.rstrip())
            executive_control_lines.append(line.rstrip())

        elif flag == 2 or flag < 0:
            # we're in the case control deck right now and looking
            # for one of the following:
            #  - 'BEGIN BULK'
            #  - 'BEGIN SUPER=1'
            #  - 'BEGIN BULK AUXMODEL=200'
            #  - 'BEGIN BULK AFPM=300'
            #
            # There's a special case for 'BEGIN AUXMODEL=1', so we flag
            # AUXCASE/AUXMODEL, and do some extra parsing in flag=3.
            #
            # flag=2 (BEGIN BULK)
            # flag=-1 (BEGIN SUPER=1)
            # flag=-2 (BEGIN SUPER=2)
            # ...
            #
            # We haven't yet tried to handle the AFPM special case

            # we have to handle the comment because we could incorrectly
            # flag the model as flipping to the BULK data section if we
            # have BEGIN BULK in a comment
            if '$' in line:
                line, comment = line.split('$', 1)
                current_lines.append('$' + comment.rstrip())
                #print('%s: %s' % (flag_word, '$' + comment.rstrip()))

            uline = line.upper().strip()
            if uline.startswith('BEGIN'):
                if _is_begin_bulk(uline):
                    old_flags.append(flag)
                    #assert flag == 2, flag

                    # we're about to break because we found begin bulk
                    flag = 3
                    current_ilines = bulk_data_ilines

                    #or not keep_enddata
                    is_extra_bulk = (is_auxmodel or is_afpm or
                                     is_superelement or consider_superelements)

                    if not is_extra_bulk:
                        #print('breaking begin bulk...')
                        bulk_data_ilines = _bulk_data_lines_extract(
                            lines, ilines, bulk_data_lines, i,
                            make_ilines=make_ilines, keep_enddata=keep_enddata)
                        break
                    #print('setting lines to bulk---')
                    current_lines = bulk_data_lines
                    #flag_word = 'bulk'
                    #print('case: %s' % (line.rstrip()))
                    case_control_lines.append(line.rstrip())
                    continue

                elif 'SUPER' in uline and '=' in uline:
                    super_id = _get_super_id(line, uline)
                    old_flags.append(flag)
                    flag = -super_id
                    #flag_word = 'SUPER=%s' % super_id
                    current_lines = superelement_lines[super_id]
                    current_ilines = superelement_ilines[super_id]

                elif ('AUXMODEL' in uline or 'AFPM' in uline) and '=' in uline:
                    out = _read_bulk_for_auxmodel(
                        ifile_iline, line, flag, bulk_data_lines,
                        current_lines, current_ilines,
                        old_flags,
                        is_auxmodel, auxmodel_lines, auxmodels_to_find, auxmodels_found,
                        is_afpm, afpm_lines, afpms_to_find, afpms_found,
                        superelement_lines, superelement_ilines,
                        is_auxmodel_active, auxmodel_id,
                        is_afpm_active, afpm_id,
                        bulk_data_ilines)
                    (is_broken,
                     auxmodel_id, is_auxmodel_active,
                     afpm_id, is_afpm_active,
                     flag, current_lines) = out
                    if is_broken:
                        break
                else:
                    msg = 'expected "BEGIN BULK" or "BEGIN SUPER=1"\nline = %s' % line
                    raise RuntimeError(msg)

                #print('%s: %s' % (flag_word, line.rstrip()))
                current_lines.append(line.rstrip())
            elif uline.startswith('SUPER'):
                # case control line
                # SUPER = ALL
                #auxmodel_idi = int(uline.split('=')[1])
                #auxmodels_to_find.append(auxmodel_idi)
                if flag != 2:
                    raise RuntimeError('expected a flag of 2 (case control deck) '
                                       'when going to an SUPER model')

                is_superelement = True
            elif uline.startswith('AUXMODEL'):
                # case control line
                # AUXMODEL = 10
                auxmodel_idi = int(uline.split('=')[1])
                auxmodels_to_find.append(auxmodel_idi)
                if flag != 2:
                    raise RuntimeError('expected a flag of 2 (case control deck) '
                                       'when going to an AUXMODEL')
                is_auxmodel = True
            elif uline.startswith('AFPM'):
                # case control line
                # AFPM = 10
                afpm_idi = int(uline.split('=')[1])
                afpms_to_find.append(afpm_idi)
                if flag != 2:
                    raise RuntimeError('expected a flag of 2 (case control deck) '
                                       'when going to an AFPM model')
                is_afpm = True

            #print('%s: %s' % (flag_word, line.rstrip()))
            current_lines.append(line.rstrip())
        elif flag == 3:
            if not(is_auxmodel is True or is_superelement is True or consider_superelements):
                raise RuntimeError(f'one must be True: is_auxmodel={is_auxmodel}; '
                                   'is_superelement={is_superelement}; '
                                   'consider_superelements={consider_superelements}')
            #assert is_auxmodel is True or is_superelement is True or consider_superelements

            # we have to handle the comment because we could incorrectly
            # flag the model as flipping to the BULK data section if we
            # have BEGIN BULK in a comment
            if '$' in line:
                line, comment = line.split('$', 1)
                current_lines.append('$' + comment.rstrip())
                if bulk_data_ilines != current_ilines:
                    raise RuntimeError('bulk_data_ilines != current_ilines')

                current_ilines.append(ifile_iline)
                #bulk_data_ilines.append(ifile_iline)

            out = _read_bulk_for_auxmodel(
                ifile_iline, line, flag, bulk_data_lines,
                current_lines, current_ilines,
                old_flags,
                is_auxmodel, auxmodel_lines, auxmodels_to_find, auxmodels_found,
                is_afpm, afpm_lines, afpms_to_find, afpms_found,
                superelement_lines, superelement_ilines,
                is_auxmodel_active, auxmodel_id,
                is_afpm_active, afpm_id,
                bulk_data_ilines)
            (is_broken,
             auxmodel_id, is_auxmodel_active,
             afpm_id, is_afpm_active,
             flag, current_lines) = out
            if is_broken:
                #print('breaking...')
                break
        else:
            raise RuntimeError(line)

    _check_valid_deck(flag, old_flags)

    if len(bulk_data_lines) == 0:
        raise RuntimeError('no bulk data lines were found')
    #print('nbulk=%s nilines=%s' % (len(bulk_data_lines),
                                   #len(bulk_data_ilines)), bulk_data_ilines.shape)

    #if bulk_data_ilines is not None and len(bulk_data_lines) != len(bulk_data_ilines):
        #raise RuntimeError('nbulk=%s nilines=%s' % (len(bulk_data_lines), len(bulk_data_ilines)))
        #print('nbulk=%s nilines=%s' % (
            #len(bulk_data_lines), len(bulk_data_ilines)), bulk_data_ilines.shape)

    bulk_data_ilines = np.asarray(bulk_data_ilines)

    out = (
        executive_control_lines, case_control_lines,
        bulk_data_lines, bulk_data_ilines,
        superelement_lines, superelement_ilines,
        auxmodel_lines, afpm_lines,
    )
    return out

def _bulk_data_lines_extract(lines: List[str],
                             ilines: Any,
                             bulk_data_lines: List[str],
                             i: int,
                             make_ilines: bool=True,
                             keep_enddata: bool=True) -> NDArrayN2int:
    """grabs the bulk data lines and ilines when we're breaking"""
    if keep_enddata:
        for line in lines[i+1:]:
            bulk_data_lines.append(line.rstrip())
        if make_ilines:
            bulk_data_ilines = ilines[i+1:, :]
    else:
        bulk_data_ilines = None
        j = 0
        for j, line in enumerate(lines[i+1:]):
            rline = line.rstrip()
            if rline.upper().startswith('ENDDATA'):
                break
            bulk_data_lines.append(rline)
        if make_ilines:
            bulk_data_ilines = ilines[i+1:i+j+1, :]

    #if not len(bulk_data_lines) == len(bulk_data_ilines):
        #msg = 'len(bulk_data_lines)=%s len(bulk_data_ilines)=%s' % (
            #len(bulk_data_lines), len(bulk_data_ilines))
        #raise RuntimeError(msg)
    return bulk_data_ilines

def _is_begin_bulk(uline: str) -> bool:
    """
    is this a:
      'BEGIN BULK'
    but not:
      'BEGIN BULK SUPER=2'
      'BEGIN BULK AUXMODEL=2'
      'BEGIN BULK AFPM=2'

    """
    is_begin_bulk = 'BULK' in uline and (
        'AUXMODEL' not in uline and
        'AFPM' not in uline and
        'SUPER' not in uline)
    return is_begin_bulk

def _read_bulk_for_auxmodel(ifile_iline, line, flag: int, bulk_data_lines: List[str],
                            current_lines, current_ilines,
                            old_flags,
                            unused_is_auxmodel, auxmodel_lines, auxmodels_to_find, auxmodels_found,
                            unused_is_afpm, afpm_lines, afpm_to_find, afpm_found,
                            superelement_lines, superelement_ilines,
                            is_auxmodel_active: bool, auxmodel_id: int,
                            is_afpm_active: bool, afpm_id: int,
                            bulk_data_ilines):
    """
    Reads a BEGIN BULK section searching for 'BEGIN AUXMODEL=1' and BEGIN SUPER=1'
    """
    # we're in the bulk data deck right now and looking
    # for a 'BEGIN AUXMODEL=1' or ???.
    #
    #print(len(bulk_data_lines), len(bulk_data_ilines))
    if len(bulk_data_lines) != len(bulk_data_ilines):
        raise RuntimeError('len(bulk_data_lines)=%s len(bulk_data_ilines)=%s are not equal' % (
            len(bulk_data_lines), len(bulk_data_ilines)))
    is_broken = False

    #if not is_auxmodel:
        #print('breaking B', flag)
        #is_broken = True
        #return is_broken, auxmodel_id, is_auxmodel_active, flag, current_lines

    uline = line.upper().strip()
    if uline.startswith('BEGIN'):
        if 'AUXMODEL' in uline:
            is_auxmodel_active = True
            auxmodel_id = _get_auxmodel_id(line, uline)
            old_flags.append(flag)
            flag = -auxmodel_id
            current_lines = auxmodel_lines[auxmodel_id]
            current_ilines = []
            auxmodels_found.add(auxmodel_id)
            if len(auxmodels_found) == len(auxmodels_to_find) and len(auxmodels_found):
                #print('broken...final', len(bulk_data_lines), len(bulk_data_ilines))
                is_broken = True
                out = (is_broken,
                       auxmodel_id, is_auxmodel_active,
                       afpm_id, is_afpm_active,
                       flag, current_lines)
                return out
        elif 'SUPER' in uline:
            super_id = _get_super_id(line, uline)
            old_flags.append(flag)
            flag = -super_id
            current_lines = superelement_lines[super_id]
            current_ilines = superelement_ilines[super_id]
        elif 'AFPM' in uline:
            is_afpm_active = True
            afpm_id = _get_afpm_id(line, uline)
            old_flags.append(flag)
            flag = -afpm_id
            current_lines = afpm_lines[afpm_id]
            current_ilines = []
            afpm_found.add(afpm_id)
            if len(afpm_found) == len(afpm_to_find) and len(afpm_found):
                #print('broken...final', len(bulk_data_lines), len(bulk_data_ilines))
                is_broken = True
                return is_broken, auxmodel_id, is_auxmodel_active, flag, current_lines
        else:
            msg = 'expected "BEGIN SUPER=1", "BEGIN AUXMODEL=1" or "BEGIN AFPM=1"\nline = %s' % line
            raise RuntimeError(msg)
    rline = line.rstrip()
    if rline:
        #if flag == 3:
            #bulk_data_ilines.append(ifile_iline)
        current_lines.append(rline)
        current_ilines.append(ifile_iline)

    out = (
        is_broken,
        auxmodel_id, is_auxmodel_active,
        afpm_id, is_afpm_active,
        flag, current_lines)
    return out

def _break_system_lines(executive_control_lines: List[str]) -> Tuple[List[str], List[str]]:
    """
    Extracts the Nastran system lines.
    System lines may be interspersed with executive lines.

    Per NX Nastran 10:

    Header    Description
    ======    ===========
    ACQUIRE   Selects NDDL schema and NX Nastran Delivery Database.
    ASSIGN    Assigns physical files to DBset members or special FORTRAN files.
    CONNECT   Groups geometry data by evaluator and database.
    DBCLEAN   Deletes selected database version(s) and/or projects.
    DBDICT    Prints the database directory in user-defined format.
    DBDIR     Prints the database directory.
    DBFIX     Identifies and optionally corrects errors found in the database.
    DBLOAD    Loads a database previously unloaded by DBUNLOAD.
    DBLOCATE  Obtains data blocks and parameters from databases.
    DBSETDEL  Deletes DBsets.
    DBUNLOAD  Unloads a database for compression, transfer, or archival storage.
    DBUPDATE  Specifies the time between updates of the database directory.
    ENDJOB    Terminates a job upon completion of FMS statements.
    EXPAND    Concatenates additional DBset members to an existing DBset.
    INCLUDE   Inserts an external file in the input file.
    INIT      Creates a temporary or permanent DBset.
    NASTRAN   Specifies values for system cells.
    PROJ      Defines the current or default project identifier.

    F:\\Program Files\\Siemens\\NXNastran\\nxn10p1\\nxn10p1\\nast\\tpl\\mdb01.dat

    """
    j = None
    sol_line = None
    isol_line = None
    system_lines = []
    executive_control_lines2 = []

    # add all the lines before and including the file management section
    # to the system lines
    #
    # add the other lines (and the SOL 101) to the executive control lines
    for i, line in enumerate(executive_control_lines):
        line_upper = line.strip().upper()
        if line_upper.startswith('SOL '):
            isol_line = i+1
            sol_line = line
        if line_upper.startswith(FILE_MANAGEMENT):
            system_lines += executive_control_lines[j:i+1]
            j = i+1

    # remove SOL 101 from the system lines if it's there
    system_lines2 = [
        line for line in system_lines
        if not line.upper().strip().startswith('SOL ') and
        not line.upper().strip().startswith('CEND')
    ]

    if j is None:
        # no system lines
        executive_control_lines2 = executive_control_lines
    else:
        # append SOL 101 to the executive lines if we found it
        # inside the system section
        append_sol_line = isol_line is not None and j is not None and isol_line < j
        if append_sol_line:
            executive_control_lines2.append(sol_line)

        # add the rest of the executive control cards
        for iline in range(j, len(executive_control_lines)):
            executive_control_lines2.append(executive_control_lines[iline])

    #for line in executive_control_lines:
        #print('eline = %r' % line)
    #for line in system_lines2:
        #print('sline2 = %r' % line)
    #for line in executive_control_lines2:
        #print('eline2 = %r' % line)
    return system_lines2, executive_control_lines2


def _check_valid_deck(flag: int, old_flags: List[int]) -> None:
    """Crashes if the flag is set wrong"""
    if flag != 3:
        if flag == 1:
            found = ' - Executive Control Deck\n'
            missing = ' - Case Control Deck\n'
            missing += ' - Bulk Data Deck\n'
        elif flag == 2:
            found = ' - Executive Control Deck\n'
            found += ' - Case Control Deck\n'
            missing = ' - Bulk Data Deck\n'
        elif flag < 0:
            # superelement/auxmodel
            found = str('old_flags=%s' % old_flags)
            missing = '???'
            return
        else:
            raise RuntimeError('flag=%r is not [1, 2, 3]' % flag)

        msg = 'This is not a valid BDF (a BDF capable of running Nastran).\n\n'
        msg += 'The following sections were found:\n%s\n' % found
        msg += 'The following sections are missing:\n%s\n' % missing
        msg += 'If you do not have an Executive Control Deck or a Case Control Deck:\n'
        msg += '  1.  call read_bdf(...) with `punch=True`\n'
        msg += "  2.  Add '$ pyNastran : punch=True' to the top of the main file\n"
        msg += '  3.  Name your file *.pch\n\n'
        msg += 'You cannot read a deck that has an Executive Control Deck, but\n'
        msg += 'not a Case Control Deck (or vice versa), even if you have a Bulk Data Deck.\n'
        raise MissingDeckSections(msg)
    return

def _show_bad_file(self: Any, bdf_filename: Union[str, StringIO],
                   encoding: str,
                   nlines_previous: int=10) -> None:
    """
    Prints the 10 lines before the UnicodeDecodeError occurred.

    Parameters
    ----------
    bdf_filename : str
        the filename to print the lines of
    encoding : str
        the file encoding
    nlines_previous : int; default=10
        the number of lines to show

    """
    lines = []  # type: List[str]
    print('ENCODING - show_bad_file=%r' % encoding)

    with open(_filename(bdf_filename), 'r', encoding=encoding) as bdf_file:
        iline = 0
        nblank = 0
        while 1:
            try:
                line = bdf_file.readline().rstrip()
            except UnicodeDecodeError:
                iline0 = max([iline - nlines_previous, 0])
                self.log.error('filename=%s' % bdf_filename)
                for iline1, line in enumerate(lines[iline0:iline]):
                    self.log.error('lines[%i]=%r' % (iline0 + iline1, line))
                msg = "\n%s encoding error on line=%s of %s; not '%s'" % (
                    encoding, iline, bdf_filename, encoding)
                raise RuntimeError(msg)

            if line:
                nblank = 0
            else:
                nblank += 1
            if nblank == 20:
                raise RuntimeError('20 blank lines')
            iline += 1
            lines.append(line)

def _get_auxmodel_id(line: str, uline: str) -> int:
    """
    parses the superelement header::

        BEGIN AUXMODEL=2
        BEGIN BULK AUXMODEL=2
        BEGIN BULK AUXMODEL = 2

    """
    #if '=' in uline:
    sline = uline.split('=')
    #else:
        #sline = uline.split()
    try:
        auxmodel_id = int(sline[1])
    except (IndexError, ValueError):
        msg = 'expected "BEGIN AUXMODEL=1"\nline = %s' % line
        raise SyntaxError(msg)

    if auxmodel_id < 0:
        raise SyntaxError('auxmodel_id=%i must be greater than 0; line=%s' % (
            auxmodel_id, line))
    return auxmodel_id

def _get_afpm_id(line: str, uline: str) -> int:
    """
    parses the superelement header::

        BEGIN AFPM=2
        BEGIN BULK AFPM=2
        BEGIN BULK AFPM = 2

    """
    sline = uline.split('=')
    try:
        afpm_id = int(sline[1])
    except (IndexError, ValueError):
        msg = 'expected "BEGIN AFPM=1"\nline = %s' % line
        raise SyntaxError(msg)

    if afpm_id < 0:
        raise SyntaxError('afpm_id=%i must be greater than 0; line=%s' % (
            afpm_id, line))
    return afpm_id

def _get_super_id(line: str, uline: str) -> int:
    """
    parses the superelement header::

        BEGIN SUPER=2
        BEGIN BULK SUPER=2
        BEGIN BULK SUPER = 2
        BEGIN BULK SUPER 2

    """
    if '=' in uline:
        sline = uline.split('=')
        super_id_str = sline[1]
        if len(sline) != 2:
            msg = 'expected "BEGIN SUPER=1"\nline = %s' % line
            raise SyntaxError(msg)
    else:
        sline = uline.split()
        if len(sline) not in [3, 4]:
            msg = 'expected "BEGIN SUPER=1"\nline = %s' % line
            raise SyntaxError(msg)
        if len(sline) == 3:
            # BEGIN SUPER 2
            super_id_str = sline[2]
        elif len(sline) == 4:
            super_id_str = sline[3]

    try:
        super_id = int(super_id_str)
    except (IndexError, ValueError):
        msg = 'expected "BEGIN SUPER=1"\nline = %s' % line
        raise SyntaxError(msg)

    if super_id < 0:
        raise SyntaxError('super_id=%i must be greater than 0; line=%s' % (
            super_id, line))
    return super_id

def _clean_comment_bulk(comment: str) -> str:
    """
    Removes specific pyNastran comment lines so duplicate lines aren't
    created.

    Parameters
    ----------
    comment : str
        the comment to possibly remove

    Returns
    -------
    updated_comment : str
        the comment

    """
    if comment == '':
        pass
    elif comment in IGNORE_COMMENTS:
        comment = ''
    elif 'pynastran' in comment.lower():
        csline = comment.lower().split('pynastran', 1)
        if csline[1].strip() == ':':
            comment = ''

    #if comment:
        #print(comment)
    return comment

def _make_ilines(nlines: int, ifile: int) -> NDArrayN2int:
    """helper method"""
    ilines = np.empty((nlines, 2), dtype='int32')
    ilines[:, 0] = ifile
    ilines[:, 1] = np.arange(nlines) # 0 to N-1
    return ilines
