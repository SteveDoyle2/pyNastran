# coding: utf-8
"""
Main BDF class.  Defines:
  - BDFInputPy
"""
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
import os
from io import open
from typing import List, Dict, Optional, Union, Set, Any, cast
from six import StringIO

import numpy as np
from pyNastran.utils import print_bad_path, _filename
from pyNastran.utils.log import get_logger2
from pyNastran.bdf.bdf_interface.include_file import get_include_filename
from pyNastran.bdf.errors import MissingDeckSections

FILE_MANAGEMENT = (
    'ACQUIRE ', 'ASSIGN ', 'CONNECT ', 'DBCLEAN ', 'DBDICT ', 'DBDIR ',
    'DBFIX ', 'DBLOAD ', 'DBLOCATE ', 'DBSETDEL ', 'DBUNLOAD ',
    'DBUPDATE ', 'ENDJOB ', 'EXPAND ', 'INCLUDE ', 'INIT ', 'NASTRAN ',
    'PROJ ',
)
EXECUTIVE_CASE_SPACES = tuple(list(FILE_MANAGEMENT) + ['SOL ', 'SET ', 'SUBCASE '])


class BDFInputPy(object):
    def __init__(self, read_includes, dumplines, encoding, nastran_format='msc',
                 log=None, debug=False):
        """
        Parameters
        ----------
        read_includes : bool
            should include files be read
        dumplines : bool
            ???
        encoding : str
            ???
        nastran_format : str; default='msc'
            'zona' has a special read method
            {msc, nx, zona}
        log : logger(); default=None
            ???
        debug : bool; default=False
            ???
        """
        self.dumplines = dumplines
        self.encoding = encoding
        self.nastran_format = nastran_format
        self.include_dir = ''

        self.reject_lines = []
        self.read_includes = read_includes
        self.active_filenames = []
        self.active_filename = None

        self.debug = debug
        self.log = get_logger2(log, debug)

    def _get_lines(self, bdf_filename, punch=False, make_ilines=True):
        # type: (Union[str, StringIO], bool) -> List[str]
        """
        Opens the bdf and extracts the lines by group

        Parameters
        ----------
        bdf_filename : str
            the main bdf_filename
        punch : bool; default=False
            is this a punch file
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
        main_lines = self._get_main_lines(bdf_filename)
        all_lines, ilines = self._lines_to_deck_lines(main_lines, make_ilines=make_ilines)

        out = _lines_to_decks(all_lines, ilines, punch, keep_enddata=True)
        system_lines, executive_control_lines, case_control_lines, bulk_data_lines, bulk_data_ilines = out
        if self.nastran_format in ['msc', 'nx']:
            pass
        elif self.nastran_format == 'zona':
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
                        assert filename.endswith('.bdf'), filename

                        _main_lines = self._get_main_lines(filename)
                        _all_lines, _ilines = self._lines_to_deck_lines(_main_lines, make_ilines=False)
                        _out = _lines_to_decks(_all_lines, _ilines, punch, keep_enddata=False)
                        _system_lines, _executive_control_lines, _case_control_lines, bulk_data_lines2, _bulk_data_ilines2 = _out
                        bulk_data_lines = bulk_data_lines2 + bulk_data_lines
                        continue
                    elif header_upper.startswith('ASSIGN MATRIX'):
                        pass
                    elif header_upper.startswith('ASSIGN OUTPUT4'):
                        pass
                    else:  # pragma: no cover
                        raise NotImplementedError(system_line)
                system_lines2.append(system_line)
            system_lines = system_lines
        else:
            msg = 'nastran_format=%r and must be msc, nx, or zona' % self.nastran_format
            raise NotImplementedError(msg)
        return system_lines, executive_control_lines, case_control_lines, bulk_data_lines, bulk_data_ilines

    def _get_main_lines(self, bdf_filename):
        # type: (Union[str, StringIO]) -> List[str]
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
            assert len(lines) > 0, lines
            return lines

        bdf_filename = cast(str, bdf_filename)

        # the directory of the 1st BDF (include BDFs are relative to this one)
        self.include_dir = os.path.dirname(os.path.abspath(bdf_filename))

        with self._open_file(bdf_filename, basename=True) as bdf_file:
            try:
                lines = bdf_file.readlines()
            except:
                _show_bad_file(self, bdf_filename, encoding=self.encoding)
        return lines

    def _lines_to_deck_lines(self, lines, make_ilines=True):
        # type: (List[str]) -> List[str], int
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
                if self.read_includes:
                    lines, nlines, ilines = self._update_include(
                        lines, nlines, ilines,
                        bdf_filename2, i, j, ifile, make_ilines=make_ilines)
                    ifile += 1
                else:
                    lines = lines[:i] + lines[j:]
                    if make_ilines:
                        ilines = np.vstack([
                            ilines[:i, :],
                            ilines[j:, :],
                        ])
                    self.reject_lines.append(include_lines)
                    #self.reject_lines.append(write_include(bdf_filename2))
            i += 1

        if self.dumplines:
            self._dump_file('pyNastran_dump.bdf', lines, i)

        #if make_ilines:
            #nilines = ilines.shape[0]
            #assert nlines == ilines.shape[0], 'nlines=%s nilines=%s' % (nlines, nilines)
        return lines, ilines

    def _update_include(self, lines, nlines, ilines,
                        bdf_filename2, i, j, ifile, make_ilines=False):
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
            #msg += 'len(bdf_filename2) = %s' % len(bdf_filename2)
            print(msg)
            raise
            #raise IOError(msg)

        with self._open_file(bdf_filename2, basename=False) as bdf_file:
            #print('bdf_file.name = %s' % bdf_file.name)
            try:
                lines2 = bdf_file.readlines()
            except UnicodeDecodeError:
                msg = (
                    'Invalid Encoding: encoding=%r.  Fix it by:\n'
                    '  1.  try a different encoding (e.g., latin1, cp1252, utf8)\n'
                    "  2.  call read_bdf(...) with `encoding`'\n"
                    "  3.  Add '$ pyNastran : encoding=latin1"
                    ' (or other encoding) to the top of the main file\n' % self.encoding)
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
            ilines = np.vstack([
                ilines[:i+1, :],
                ilines2,
                ilines[j:, :],
            ])
            #dij = j - i

        nlines += nlines2
        lines = lines[:i] + [include_comment] + lines2 + lines[j:]
        #if make_ilines:
            #n_ilines = ilines.shape[0]
            #ncompare = n_ilines - dij
            #assert n_ilines == n_ilines, 'nlines=%s dij=%s n_ilines=%s' % (n_ilines, dij, n_ilines)
        #for line in lines:
            #print("  *%s" % line.rstrip())
        return lines, nlines, ilines

    def _get_include_lines(self, lines, line, i, nlines):
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

    def _dump_file(self, bdf_dump_filename, lines, i):
        # type: (str, List[str], int) -> None
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

    def _open_file_checks(self, bdf_filename, basename=False):
        # type: (str, bool) -> None
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
            print(msg)
            raise IOError(msg)
        elif bdf_filename_inc.endswith('.op2'):
            msg = 'Invalid filetype: bdf_filename=%r' % bdf_filename_inc
            print(msg)
            raise IOError(msg)
        bdf_filename = bdf_filename_inc

        if bdf_filename in self.active_filenames:
            msg = 'bdf_filename=%s is already active.\nactive_filenames=%s' \
                % (bdf_filename, self.active_filenames)
            print(msg)
            raise RuntimeError(msg)
        elif os.path.isdir(_filename(bdf_filename)):
            current_filename = self.active_filename if len(self.active_filenames) > 0 else 'None'
            msg = 'Found a directory: bdf_filename=%r\ncurrent_file=%s' % (
                bdf_filename_inc, current_filename)
            print(msg)
            raise IOError(msg)
        elif not os.path.isfile(_filename(bdf_filename)):
            msg = 'Not a file: bdf_filename=%r' % bdf_filename
            print(msg)
            raise IOError(msg)

    def _open_file(self, bdf_filename, basename=False, check=True):
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
        """
        if basename:
            bdf_filename_inc = os.path.join(self.include_dir, os.path.basename(bdf_filename))
        else:
            bdf_filename_inc = os.path.join(self.include_dir, bdf_filename)

        self._validate_open_file(bdf_filename, bdf_filename_inc, check)


        self.log.debug('opening %r' % bdf_filename_inc)
        self.active_filenames.append(bdf_filename_inc)

        #print('ENCODING - _open_file=%r' % self.encoding)
        bdf_file = open(_filename(bdf_filename_inc), 'r', encoding=self.encoding)
        return bdf_file

    def _validate_open_file(self, bdf_filename, bdf_filename_inc, check=True):
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


def _clean_comment(comment):
    # type: (str) -> Optional[str]
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


def _lines_to_decks(lines, ilines, punch, keep_enddata=True):
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
    punch : bool
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
    """
    make_ilines = ilines is not None
    executive_control_lines = []
    case_control_lines = []
    bulk_data_lines = []
    bulk_data_ilines = None
    if punch:
        bulk_data_lines = lines
        bulk_data_ilines = ilines
    else:
        flag = 1
        iline = 0
        for iline, line in enumerate(lines):
            #print(flag, line.rstrip())
            if flag == 1:
                # I don't think we need to handle the comment because
                # this uses a startswith
                if line.upper().startswith('CEND'):
                    assert flag == 1
                    flag = 2
                executive_control_lines.append(line.rstrip())

            elif flag == 2:
                # we have to handle the comment because we could incorrectly
                # flag the model as flipping to the BULK data section if we
                # have BEGIN BULK in a comment
                if '$' in line:
                    line, comment = line.split('$', 1)
                    case_control_lines.append('$' + comment.rstrip())

                uline = line.upper().strip()
                if uline.startswith('BEGIN') and ('BULK' in uline or 'SUPER' in uline):
                    assert flag == 2
                    flag = 3
                case_control_lines.append(line.rstrip())
            else:
                break

        if keep_enddata:
            for line in lines[iline:]:
                bulk_data_lines.append(line.rstrip())
            #bulk_data_lines = [line.rstrip() for line in lines[iline:]]
        else:
            for line in lines[iline:]:
                if line.upper().startswith('ENDDATA'):
                    continue
                bulk_data_lines.append(line.rstrip())

        _check_valid_deck(flag)
        if make_ilines:
            bulk_data_ilines = ilines[iline:, :]

    del lines
    #for line in bulk_data_lines:
        #print(line)

    # break out system commands
    system_lines, executive_control_lines = _break_system_lines(executive_control_lines)

    # clean comments
    system_lines = [_clean_comment(line) for line in system_lines
                    if _clean_comment(line) is not None]
    executive_control_lines = [_clean_comment(line) for line in executive_control_lines
                               if _clean_comment(line) is not None]
    case_control_lines = [_clean_comment(line) for line in case_control_lines
                          if _clean_comment(line) is not None]
    return system_lines, executive_control_lines, case_control_lines, bulk_data_lines, bulk_data_ilines

def _break_system_lines(executive_control_lines):
    """
    Extracts the Nastran system lines

    Per NX Nastran 10:

    ACQUIRE Selects NDDL schema and NX Nastran Delivery Database.
    ASSIGN Assigns physical files to DBset members or special FORTRAN
    files.
    CONNECT Groups geometry data by evaluator and database.
    DBCLEAN Deletes selected database version(s) and/or projects.
    DBDICT Prints the database directory in user-defined format.
    DBDIR Prints the database directory.
    DBFIX Identifies and optionally corrects errors found in the database.
    DBLOAD Loads a database previously unloaded by DBUNLOAD.
    DBLOCATE Obtains data blocks and parameters from databases.
    DBSETDEL Deletes DBsets.
    DBUNLOAD Unloads a database for compression, transfer, or archival
    storage.
    DBUPDATE Specifies the time between updates of the database directory.
    ENDJOB Terminates a job upon completion of FMS statements.
    EXPAND Concatenates additional DBset members to an existing DBset.
    INCLUDE Inserts an external file in the input file.
    INIT Creates a temporary or permanent DBset.
    NASTRAN Specifies values for system cells.
    PROJ Defines the current or default project identifier.

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


def _check_valid_deck(flag):
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


def _show_bad_file(self, bdf_filename, encoding, nlines_previous=10):
    # type: (Union[str, StringIO]) -> None
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


def _clean_comment_bulk(comment):
    # type: (str) -> str
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

def _make_ilines(nlines, ifile):
    """helper method"""
    ilines = np.empty((nlines, 2), dtype='int32')
    ilines[:, 0] = ifile
    ilines[:, 1] = np.arange(nlines) # 0 to N-1
    return ilines
