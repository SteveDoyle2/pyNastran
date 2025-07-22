# coding: utf-8
"""
Main BDF class.  Defines:
  - BDFInputPy

BEGIN [BULK] = [
    AFPM = afpmid
    ARBMODEL = arbmid
    AUXMODEL = auxmind
    MASSID = massid[LABEL = masslabel]
    MODULE= moduleid[APPEND][LABEL = modlabel]
    FLXBDY = flexbody
    SUPER = seid
    TRMC = trimid
    UDS
]

BEGIN BULK
BEGIN AUXMODEL=22
BEGIN BULK TRMC=101
BEGIN TRMC=102
"""
import os
import shlex
import warnings
from io import StringIO
from itertools import count
from collections import defaultdict
from typing import Optional, Any, cast

import numpy as np
# from cpylog import get_logger, SimpleLogger
from pyNastran.nptyping_interface import NDArrayN2int
from pyNastran.utils import print_bad_path, PathLike

from pyNastran.bdf import BULK_DATA_CARDS, CASE_BULK_CARDS
from pyNastran.bdf.errors import AuxModelError, MissingDeckSections, SuperelementFlagError
from pyNastran.bdf.bdf_interface.utils import _parse_pynastran_header
from pyNastran.bdf.bdf_interface.include_file import get_include_filename, parse_include_lines
from cpylog import SimpleLogger, __version__ as CPYLOG_VERSION
if CPYLOG_VERSION > '1.6.0':
    from cpylog import get_logger
else:  # pragma: no cover
    from cpylog import get_logger2 as get_logger


# these allow spaces
FILE_MANAGEMENT = (
    'ACQUIRE ', 'ASSIGN ', 'CONNECT ', 'DBCLEAN ', 'DBDICT ', 'DBDIR ',
    'DBFIX ', 'DBLOAD ', 'DBLOCATE ', 'DBSETDEL ', 'DBUNLOAD ',
    'DBUPDATE ', 'ENDJOB ', 'EXPAND ', 'INCLUDE ', 'INIT ', 'NASTRAN ',
    'PROJ ',
)
EXECUTIVE_SPACES = ('ALTER ', 'APP ', 'COMPILE ', 'COMPILER ', 'DIAG ',
                    'GEOMCHECK ', 'ID ', 'LINK ', 'MALTER ', 'SOL ', 'TIME ')

#CASE_BULK_CARDS = {
    ## of the form 'LOAD = 5', so 'PARAM,POST,-1' doesn't count
    #'LOAD', 'SPC', 'FREQ', 'MPC',  # case control + bulk data cards
    #'FORCE', 'TRIM', 'DESVAR', 'TSTEP', 'TSTEPNL', 'NSM', 'CLOAD', 'SUPORT1',
    #'CSSCHD', 'SDAMPING', 'DLOAD', 'TRIM',
    #'SUPORT', # short for SUPORT1
    #'ACCEL',  # short for ACCELERATION
    ## 'PARAM', # equals sign is problematic
#}

CASE_CARDS_NO_BULK = (
    # 'LOAD'
    # NX 2019.2
    'A2GG', 'ACCE',  # ACCELERATION
    'ACINTENSITY', 'ACORDCHK', 'ACPOWER', 'ACVELOCITY', 'ADAMSMNF',
    'ADAPTERR', 'ADMRECVR', 'AECONFIG', 'AEROF', 'AESYMXY', 'AESYMXZ', 'ALOAD',
    'ANALYSIS', 'APRESSURE', 'ATVOUT', 'AUXCASE', 'B2GG', 'B2PP', 'BC', 'BCRESULTS',
    'BCSET', 'BGRESULTS', 'BGSET', 'BOLTID', 'BOLTRESULTS', 'BOUTPUT', 'CKGAP',
    #CLOAD - in bulk
    'CMETHOD', 'CRSTRN', 'CSMSET',
    #CSSCHD - in bulk
    'CYCFORCES', 'CYCSET', 'CZRESULTS',
    'DEFORM', 'DESGLB', 'DESOBJ', 'DESSUB',
    'DISP',  # DISPLACEMENT
    # DIVERG - in bulk
    # DLOAD - in bulk
    'DMTRCOEF', 'DMTRLOSS', 'DRSPAN', 'DSAPRT', 'DSYM', 'DTEMP', 'EBDSET',
    # 'ECHO'
    'EDE', 'EFLOAD', 'EKE', 'ELAR', 'ELAROUT', 'ELSDCON', 'ELSTRN', 'ELSUM',
    'ENTHALPY', 'ERP', 'ESE', 'EXTSEOUT', 'FLSFSEL', 'FLSPOUT', 'FLSTCNT',
    'FLUX', 'FLXSLI', 'FLXRESULTS', 'FMETHOD',
    #FORCE - in bulk
    'FREQU',  # 'FREQUENCY',
    'FRFIN',
    'GCRSTRN', 'GELSTRN', 'GKRESULTS', 'GPFORCE', 'GPKE', 'GPLSTRN', 'GPRSORT',
    'GPSDCON', 'GPSTRAIN', 'GPSTRESS', 'GRDCON', 'GROUNDCHECK', 'GSTRAIN', 'GSTRESS',
    'GTHSTRN', 'GUST', 'HARMONICS', 'HDOT', 'HOUTPUT', 'IC', 'IMPERF',
    #INCLUDE - already handled
    'INITS', 'INPOWER', 'JCONSET', 'JRESULTS', 'JINTEG', 'K2GG', 'K2PP', 'LABEL',
    'LINE',
    #LOAD - in bulk
    'LOADSET', 'M2GG', 'M2PP', 'MASTER', 'MAXLINES', 'MAXMIN', 'MBDEXPORT', 'MBDRECVR',
    'MEFFMASS', 'METHOD', 'MFLUID', 'MODALE', 'MODCON', 'MODES', 'MODSEL',
    'MODTRAK', 'MONITOR', 'MONVAR',
    #MPC - in bulk
    'MPCF',  # MPCFORCES
    'MPRES', 'NLARCL', 'NLCNTL', 'NLLOAD', 'NLPARM', 'NLSTRESS', 'NONLINEAR',
    'NOUTPUT',
    #NSM - in bulk
    'OFREQUENCY', 'OLOAD', 'OMODES', 'OPRESS', 'OSTNINI', 'OTEMP', 'OTIME',
    'OTMFORC', 'OUTPUT', 'P2G', 'PANCON',
    #PARAM - in bulk
    'PARTN', 'PEAKOUT', 'PFRESULTS', 'PLOTID', 'PLSTRN', 'PRESSURE',
    'RANDOM', 'RCROSS', 'REPCASE', 'RESVEC', 'RIGID', 'RMAXMIN', 'RMETHOD',
    'RSMETHOD',
    'SACCEL'  # SACCELERATION
    'SDAMPING',
    'SDISP',  # SDISPLACEMENT
    'SEALL', 'SEDR', 'SEDV',
    'SEEXCLUDE', 'SEFINAL',
    'SEKREDUCE',
    'SELGENERATE', 'SELREDUCE',
    'SEMGENERATE', 'SEMREDUCE',
    'SEQDEP', 'SERESP',
    #SET - in bulk
    'SETMC', 'SETMCNAME', 'SETS DEFINITION', 'SHELLTHK', 'SKIP',
    'SMETHOD',
    # SPC - in bulk
    'SPCF',  # SPCFORCES
    'STATSUB', 'STATVAR', 'STRAIN', 'STRESS', 'STRFIELD', 'SUBCASE', 'SUBCOM',
    'SUBSEQ',
    'SUBT',  # SUBTITLE
    #'SUPER',
    # SUPORT - in bulk
    # SUPORT1 - in bulk
    'SVECTOR',
    'SVELO',  # SVELOCITY
    'SYM', 'SYMCOM', 'SYMSEQ', 'TEMPERATURE', 'TFL', 'THERMAL', 'TITLE',
    # TRIM - in bulk
    'TRLOSS', 'TRPOWER',
    # TSTEP - in bulk
    # TSTEPNL - in bulk
    'TSTRU', 'VATVOUT',
    'VELO',  # VELOCITY
    'VOLUME', 'WEIGHTCHECK',
)

EXECUTIVE_CASE_SPACES = tuple(
    list(FILE_MANAGEMENT) +
    list(EXECUTIVE_SPACES) +
    ['SET ', 'SUBCASE '],
)


class BDFInputPy:
    """BDF reader class that only handles lines and not building cards or parsing cards"""
    def __init__(self, read_includes: bool, dumplines: bool,
                 encoding: str, nastran_format: str='msc',
                 replace_includes: Optional[dict[str, str]]=None,
                 consider_superelements: bool=True,
                 log: SimpleLogger=None, debug: bool=False):
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
        if replace_includes is None:
            replace_includes = {}

        # not used
        self.replace_includes: dict[str, str] = replace_includes

        self.dumplines: bool = dumplines
        self.encoding: str = encoding
        self.nastran_format: str = nastran_format
        self.include_dir: str = ''

        self.include_lines = defaultdict(list)
        self.read_includes: bool = read_includes
        # the list of files in the order they were opened
        self.loaded_filenames: list[str] = []
        self.active_filenames: list[str] = []
        self.active_filename = None

        self.consider_superelements: bool = consider_superelements
        self.debug: bool = debug
        self.log = get_logger(log, debug)
        self.use_new_parser: bool = False

    def get_lines(self, bdf_filename: PathLike | StringIO,
                  punch: Optional[bool]=False,
                  make_ilines: bool=True) -> tuple[list[str], list[str], list[str],
                                                   list[str], Optional[np.ndarray],
                                                   dict[tuple[str, str], list[str]],
                                                    #list[str], Optional[np.ndarray],
                                                   ]:
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
        system_lines : list[str]
            the system control lines (typically empty; used for alters)
        executive_control_lines : list[str]
            the executive control lines (stores SOL 101)
        case_control_lines : list[str]
            the case control lines (stores subcases)
        bulk_data_lines : list[str]
            the bulk data lines (stores geometry, boundary conditions, loads, etc.)
        bulk_data_ilines : None / (nlines, 2) int ndarray
            if make_ilines = True:
                the [ifile, iline] pair for each line in the file
            if make_ilines = False:
                 ilines = None

        superelement_lines : dict[int, list[str]]
            the superelement data lines (stores geometry, boundary conditions, loads, etc.)
        superelement_ilines : dict[int, np.ndarray]
            data: None / (nlines, 2) int ndarray  ???
            if make_ilines = True:
                the [ifile, iline] pair for each line in the file
            if make_ilines = False:
                 ilines = None

        """
        self.include_dir = (
            '' if isinstance(bdf_filename, StringIO) else
            os.path.dirname(bdf_filename))

        main_lines = self.get_main_lines(bdf_filename)
        all_lines, ilines = self.lines_to_deck_lines(main_lines)
        if not make_ilines:
            ilines = None

        nastran_format = self.nastran_format
        #superelement_lines = {}
        #superelement_ilines = None

        log = self.log
        if self.use_new_parser and nastran_format != 'optistruct':  # make_ilines is False
            out = lines_to_decks2(all_lines, ilines, punch, log,
                                  nastran_format=nastran_format)
            (system_lines, executive_control_lines, case_control_lines,
             bulk_data_lines, bulk_data_ilines,
             additional_deck_lines, additional_deck_ilines) = out
            if len(bulk_data_lines) != len(bulk_data_ilines):
                msg = 'len(bulk_data_lines)=%s len(bulk_data_ilines)=%s' % (
                    len(bulk_data_lines), len(bulk_data_ilines))
                self.log.warning(msg)
        else:
            #print(self.use_new_parser, nastran_format != 'optistruct')
            self.log.debug('using old deck splitter')
            out = _lines_to_decks(all_lines, ilines, punch, log,
                                  keep_enddata=True,
                                  consider_superelements=self.consider_superelements,
                                  nastran_format=nastran_format)
            (
                system_lines, executive_control_lines, case_control_lines,
                bulk_data_lines, bulk_data_ilines,
                superelement_lines, superelement_ilines) = out
            if len(bulk_data_lines) != len(bulk_data_ilines):
                msg = 'len(bulk_data_lines)=%s len(bulk_data_ilines)=%s' % (
                    len(bulk_data_lines), len(bulk_data_ilines))
                self.log.warning(msg)

            additional_deck_lines = {}
            additional_deck_ilines = {}
            for superelement_id, super_lines in superelement_lines.items():
                assert isinstance(super_lines, list), super_lines
                super_tuple = ('SUPER', superelement_id, '')
                additional_deck_lines[super_tuple] = super_lines
            del superelement_lines, superelement_ilines

        if nastran_format in {'msc', 'nx', 'optistruct', 'mystran'}:
            pass
        elif nastran_format == 'zona':
            bulk_data_lines, bulk_data_ilines, system_lines = self._get_lines_zona(
                system_lines, bulk_data_lines, bulk_data_ilines, punch)
        else:
            msg = f'nastran_format={nastran_format!r} and must be msc, nx, optistruct, mystran, or zona'
            raise NotImplementedError(msg)

        #for line in bulk_data_lines:
            #print(line.rstrip())
        return (system_lines, executive_control_lines, case_control_lines,
                bulk_data_lines, bulk_data_ilines,
                additional_deck_lines, additional_deck_ilines)

    def _get_lines_zona(self, system_lines: list[str],
                        bulk_data_lines: list[str],
                        bulk_data_ilines: NDArrayN2int,
                        punch: bool) -> tuple[list[str], NDArrayN2int, list[str]]:
        """load and update the lines for ZONA"""
        system_lines2 = []
        log = self.log
        for system_line in system_lines:
            if system_line.upper().startswith('ASSIGN'):
                split_system = system_line.split(',')
                header = split_system[0]
                header_upper = header.upper()
                if header_upper.startswith('ASSIGN FEM'):
                    unused_fem, filename = header.split('=')
                    filename = filename.strip('"\'')
                    log.debug(f'reading {filename}')
                    filename_lower = filename.lower()
                    if filename_lower.endswith('.f06'):
                        if os.path.exists(filename):
                            log.debug(f'reading geom from f06: {filename}')
                            from pyNastran.f06.parse_geom import parse_f06_geom
                            out = parse_f06_geom(filename)
                            (_system_lines, _executive_control_lines,
                             _case_control_lines, bulk_data_lines2) = out
                            is_bdf = False
                            bulk_data_ilines2 = np.arange(len(bulk_data_ilines2))
                            bulk_data_lines = bulk_data_lines2 + bulk_data_lines
                            bulk_data_ilines = np.vstack([bulk_data_ilines2, bulk_data_ilines])
                        else:
                            log.debug(f'reading geom from bdf: {filename}')
                            is_bdf = True
                            filename = os.path.splitext(filename)[0] + '.bdf'
                    else:
                        if filename_lower.endswith(('.fre', '.mod')):
                            continue
                        is_bdf = True
                        if not filename_lower.endswith('.bdf'):
                            raise RuntimeError(f'filename must end in bdf; {filename}')

                    if is_bdf:
                        if not os.path.exists(filename):
                            log.error(f'skipping {filename}; missing')
                            continue
                        _main_lines = self.get_main_lines(filename)
                        make_ilines = bulk_data_ilines is not None
                        all_lines, ilines = self.lines_to_deck_lines(_main_lines)
                        if not make_ilines:
                            ilines = None
                        _out = _lines_to_decks(all_lines, ilines, punch, self.log,
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

    def get_main_lines(self, bdf_filename: PathLike | StringIO) -> list[str]:
        """
        Opens the bdf and extracts the lines

        Parameters
        ----------
        bdf_filename : str / StringIO
            the main bdf_filename

        Returns
        -------
        lines : list[str]
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

    def lines_to_deck_lines(self, lines: list[str]) -> tuple[list[str], np.ndarray]:
        """
        Merges the includes into the main deck.

        Parameters
        ----------
        lines : list[str]
            the lines from the main BDF

        Returns
        -------
        active_lines : list[str]
            all the active lines in the deck
        ilines : (nlines, 2) int ndarray
            the [ifile, iline] pair for each line in the file

        """
        log = self.log
        replace_includes = self.replace_includes
        nlines = len(lines)
        #bdf_filenames = [self.bdf_filename]

        make_ilines = True
        ilines = _make_ilines(nlines, ifile=0)

        i = 0
        ifile = 1
        while i < nlines:
            try:
                line = lines[i].rstrip('\r\n\t')
            except IndexError:
                break
            line_upper = line.upper()
            if line_upper.startswith('INCLUDE'):
                j, include_lines = self._get_include_lines(lines, line, i, nlines)

                jfile = ilines[i, 0]
                jline = ilines[i, 1]
                source_filename = os.path.abspath(self.loaded_filenames[jfile])
                include_dir = os.path.abspath(os.path.dirname(source_filename))
                if self.debug and 0:  # pragma: no cover
                    print(f'jfile={jfile} jline={jline}')
                    print(f'-> include_line = {line_upper}')
                    print(f'-> source_file  = {source_filename}')
                    print(f'-> include_dir   = {include_dir}')
                    print(f'-? include_diri  = {self.include_dir}')
                include_dirs = [include_dir]
                if self.include_dir not in include_dirs:
                    include_dirs.append(self.include_dir)
                if self.read_includes:
                    bdf_filename2 = get_include_filename(
                        log, include_lines,
                        include_dirs=include_dirs,
                        replace_includes=replace_includes,
                        source_filename=source_filename,
                        write_env_on_error=False,
                    )
                    if bdf_filename2 == '':
                        # replace filename
                        log.warning(f'dropping {include_lines}')
                        lines[i] = f'$ removed {include_lines}\n'
                        i += 1
                        continue

                    # these are the lines associated with the 1st/2nd include file found
                    self.include_lines[jfile].append((include_lines, bdf_filename2))

                    lines, nlines, ilines = self._update_include(
                        lines, nlines, ilines,
                        include_lines, bdf_filename2, i, j, ifile)
                    ifile += 1
                else:
                    bdf_filename2 = parse_include_lines(include_lines)

                    # these are the lines associated with the 1st/2nd include file found
                    self.include_lines[jfile].append((include_lines, bdf_filename2))

                    # remove the include lines
                    lines = lines[:i] + lines[j:]
                    if make_ilines:
                        ilines = np.vstack([
                            ilines[:i, :],
                            ilines[j:, :],
                        ])
                        #assert len(ilines[:, 1]) == len(np.unique(ilines[:, 1]))
                    # we have to reprocess the line just in case there is
                    # an include on the current line
                    i -= 1
                    #self.include_lines[ifile].append(include_lines)
                    #self.reject_lines.append(write_include(bdf_filename2))
            i += 1

        bdf_dump_filename = os.path.join(self.include_dir, 'pyNastran_dump.bdf')
        if self.dumplines:
            _dump_file(bdf_dump_filename, lines, i, self.encoding)

        #print(bdf_filenames)
        #if make_ilines:
            #nilines = ilines.shape[0]
            #assert nlines == ilines.shape[0], 'nlines=%s nilines=%s' % (nlines, nilines)
        return lines, ilines

    def _update_include(self, lines: list[str], nlines: int, ilines: np.ndarray,
                        include_lines: list[str], bdf_filename2: str,
                        i: int, j: int, ifile: int,
                        ) -> tuple[list[str], int, np.ndarray]:
        """incorporates an include file into the lines

        Returns
        -------
        lines : list[str]
           ???
        nlines : int
           ???
        ilines : np.ndarray
           ???
        """
        try:
            self._open_file_checks(bdf_filename2)
        except IOError:
            crash_name = os.path.join(self.include_dir, 'pyNastran_include_error.bdf')
            _dump_file(crash_name, lines, j, self.encoding)
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
        if 1:  # make_ilines:
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

    def _get_include_lines(self, lines: list[str], line: str,
                           i: int, nlines: int) -> tuple[int, list[str]]:
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
                        msg += f'Check the end of {crash_name!r}\n'
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

    def _open_file_checks(self, bdf_filename: str,
                          basename: bool=False) -> None:
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

        log = self.log
        if not os.path.exists(bdf_filename_inc):
            msg = f'No such bdf_filename: {bdf_filename_inc!r}\n'
            msg += f'cwd: {os.getcwd()!r}\n'
            msg += f'include_dir: {self.include_dir!r}\n'
            msg += print_bad_path(bdf_filename_inc)
            log.error(msg)
            raise IOError(msg)
        elif bdf_filename_inc.endswith('.op2'):
            msg = f'Invalid filetype: bdf_filename={bdf_filename_inc!r}'
            log.error(msg)
            raise IOError(msg)
        bdf_filename = bdf_filename_inc

        current_filename_str = f'active_file={self.active_filename!r}' if len(self.active_filenames) > 0 else ''
        if bdf_filename in self.active_filenames:
            active_filenames_str = '\n - '.join(self.active_filenames)
            msg = (f'bdf_filename={bdf_filename} is already active.\n'
                   f'{current_filename_str}\n'
                   #f'referenced_file={referenced_file}\n'
                   f'active_filenames:\n - {active_filenames_str}')
            log.error(msg)
            raise RuntimeError(msg)
        elif os.path.isdir(bdf_filename):
            msg = f'Found a directory: bdf_filename={bdf_filename_inc!r}\n{current_filename_str}'
            log.error(msg)
            raise IOError(msg)
        elif not os.path.isfile(bdf_filename):
            msg = 'Not a file: bdf_filename=%r' % bdf_filename
            log.error(msg)
            raise IOError(msg)

    def _open_file(self, bdf_filename: str | StringIO,
                   basename: bool=False, check: bool=True,
                   encoding: Optional[str]=None) -> Any:
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

        # bdf_filename
        self._validate_open_file(bdf_filename_inc, check)

        self.log.debug('opening %r' % bdf_filename_inc)
        self.active_filenames.append(bdf_filename_inc)
        self.loaded_filenames.append(bdf_filename_inc)

        #print('ENCODING - _open_file=%r' % self.encoding)
        #self._check_pynastran_header(lines)
        bdf_file = open(bdf_filename_inc, 'r', encoding=encoding)
        return bdf_file

    def _validate_open_file(self,
                            bdf_filename_inc: str,
                            check: bool=True) -> None:
        """
        checks that the file doesn't have obvious errors
         - hasn't been used
         - not a directory
         - is a file

        Parameters
        ----------
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
            if not os.path.exists(bdf_filename_inc):
                msg = 'No such bdf_filename: %r\n' % bdf_filename_inc
                msg += 'cwd: %r\n' % os.getcwd()
                msg += 'include_dir: %r\n' % self.include_dir
                msg += print_bad_path(bdf_filename_inc)
                raise IOError(msg)
            elif bdf_filename_inc.endswith('.op2'):
                raise IOError('Invalid filetype: bdf_filename=%r' % bdf_filename_inc)

            bdf_filename = bdf_filename_inc
            if bdf_filename in self.active_filenames:
                active_filenames_str = '\n - '.join(self.active_filenames)
                msg = 'bdf_filename=%s is already active.\nactive_filenames:\n - %s' \
                    % (bdf_filename, active_filenames_str)
                raise RuntimeError(msg)
            elif os.path.isdir(bdf_filename):
                current_fname = self.active_filename if len(self.active_filenames) > 0 else 'None'
                raise IOError('Found a directory: bdf_filename=%r\ncurrent_file=%s' % (
                    bdf_filename_inc, current_fname))
            elif not os.path.isfile(bdf_filename):
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
    card_name = text2.strip()[0:8].split(',')[0].upper().rstrip(' *')
    #card_name2 = text2.split(',')[0].upper()
    #print('card_name =%r' % card_name)
    #print('card_name2=%r' % card_name2)

    # bulk data cards
    if card_name in BULK_DATA_CARDS:
        # case control + bulk data cards
        if '=' in text2 and card_name in CASE_BULK_CARDS or text2.startswith(' '):
            return False
        elif card_name == 'PARAM':
            # The PARAM card can have a comma or tab, but no equals sign.
            # If there is a PARAM card, we have to assume we're not in the
            # case control.
            return False
        return True
    return False


def _is_case_control_line(text: str) -> bool:
    """
    Returns True if there is a Case Control Deck

    Parameters
    ----------
    text : str
        a line in the deck

    Returns
    -------
    is_case_line : bool
        is this a case control line

    """
    line_upper = text.split('$')[0].strip().upper()
    if line_upper.startswith(CASE_CARDS_NO_BULK):
        return True
    return False


def _check_pynastran_encoding(bdf_filename: PathLike | StringIO,
                              encoding: str) -> str:
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
        end = csline[1].strip()
        if len(end) == 0:
            return None

        if end[0] == ':':
            comment = None

    #if comment:
        #print(comment)
    return comment


def _lines_to_decks(lines: list[str],
                    ilines: NDArrayN2int,
                    punch: Optional[bool],
                    log: SimpleLogger,
                    keep_enddata: bool=True,
                    consider_superelements: bool=False,
                    nastran_format: str='msc') -> tuple[
                        list[str], list[str], list[str], list[str], NDArrayN2int,
                        dict[int, list[str]], dict[int, np.ndarray],
                        dict[int, list[str]],
                    ]:
    """
    Splits the BDF lines into:
     - system lines
     - executive control deck
     - case control deck
     - bulk data deck

    Parameters
    ----------
    lines : list[str]
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
    system_lines : list[str]
        the system control lines (typically empty; used for alters)
    executive_control_lines : list[str]
        the executive control lines (stores SOL 101)
    case_control_lines : list[str]
        the case control lines (stores subcases)
    bulk_data_lines : list[str]
        the bulk data lines (stores geometry, boundary conditions, loads, etc.)
    bulk_data_ilines : None / (nlines, 2) int ndarray
        None : the old behavior
        narray : the [ifile, iline] pair for each line in the file
    superelement_lines : dict[int, list[str]]
        ???
    superelement_ilines : dict[int, np.ndarray]
        ???
    auxmodel_lines : list[str]
        ???

    """
    start_flag = 0  # default
    if punch and not consider_superelements:  # True
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
    elif punch and consider_superelements:  # True
        # start_flag=2 forces the model into 'begin bulk' mode
        start_flag = 2

    # typical deck
    out = _lines_to_decks_main(lines, ilines, log, punch=punch,
                               keep_enddata=keep_enddata,
                               consider_superelements=consider_superelements,
                               start_flag=start_flag,
                               nastran_format=nastran_format)
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
        raise AuxModelError('lines in auxmodel %i is empty' % auxmodel_id)

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
    assert isinstance(superelement_lines, dict), superelement_lines
    assert isinstance(superelement_ilines, dict), superelement_ilines
    return (
        system_lines, executive_control_lines, case_control_lines,
        bulk_data_lines, bulk_data_ilines,
        superelement_lines, superelement_ilines)


def _lines_to_decks_main(lines: list[str],
                         ilines: np.ndarray,
                         log: SimpleLogger,
                         punch: Optional[bool]=False,
                         keep_enddata: bool=True,
                         consider_superelements: bool=False,
                         start_flag: int=0,
                         nastran_format: str='msc') -> tuple[
                        list[str], list[str], list[str], list[str], NDArrayN2int,
                        list[str], list[str], list[str]]:
    """
    Splits the BDF lines into:
     - system lines
     - executive control deck
     - case control deck
     - bulk data deck

    Parameters
    ----------
    lines : list[str]
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
    system_executive_control_lines : list[str]
        the system control lines (typically empty; used for alters)
        and the executive control lines (stores SOL 101)
    case_control_lines : list[str]
        the case control lines (stores subcases)
    bulk_data_lines : list[str]
        the bulk data lines (stores geometry, boundary conditions, loads, etc.)
    bulk_data_ilines : None / (nlines, 2) int ndarray
        None : the old behavior
        narray : the [ifile, iline] pair for each line in the file
    superelement_lines : list[str]
        ???
    superelement_ilines : list[str]
        ???
    auxmodel_lines : list[str]
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
    is_module = False
    is_auxmodel_active = False
    is_afpm_active = False
    auxmodel_id = None
    afpm_id = None
    #---------------------------------------------
    current_lines = executive_control_lines

    if nastran_format in {'msc', 'nx', 'mystran', 'zona'}:
        flag_word = 'executive'
        flag = 1  # start from executive control deck
    elif nastran_format == 'optistruct':
        flag_word = 'case control'
        flag = 2  # case from control deck
    else:  # pragma: no cover
        raise RuntimeError(nastran_format)
    if start_flag != 0:
        flag = start_flag

    #flag = 1
    old_flags: list[int] = []
    bulk_data_ilines = []
    if ilines is None:
        ilines = count()

    #guess_deck_sections = True
    #print('guess_deck_sections =', guess_deck_sections, punch)
    for i, ifile_iline, line in zip(count(), ilines, lines):
        #print('%s %-8s %s' % (ifile_iline, flag_word, line.rstrip()))
        #print('%s %i %s' % (ifile_iline, flag, line.rstrip()))
        line_upper = line.split('$')[0].upper().strip()

        if guess_deck_sections and flag == 1 and line_upper.startswith('BEGIN'):
            # we're in the executive deck and found the bulk data deck
            section_name_map = {
                1: 'executive control',
                2: 'case control',
            }
            section_name = section_name_map[flag]

            if _is_begin_bulk(line_upper):
                #old_flags.append(flag)
                log.warning(f'currently in {section_name} deck and skipping '
                            'directly to bulk data section')
                flag = 3
                current_ilines = bulk_data_ilines
                current_lines = bulk_data_lines
                bulk_data_ilines = _bulk_data_lines_extract(
                    lines, ilines, bulk_data_lines, i,
                    make_ilines=make_ilines, keep_enddata=keep_enddata)
            else:
                raise RuntimeError(f'currently in {section_name} deck and unexpectedly '
                                   f'found the following line:\n{line}')
            break

        if guess_deck_sections and flag in [1, 2] and _is_bulk_data_line(line):
            # we found the case control deck successfully from the executive deck
            # then we found the bulk data deck unexpectedly
            section_name_map = {
                1: 'executive control',
                2: 'case control',
            }
            section_name = section_name_map[flag]
            log.warning(f'currently in {section_name} deck and skipping directly '
                        f'to bulk data section\n{line}')
            flag = 3
            current_ilines = bulk_data_ilines
            current_lines = bulk_data_lines
            bulk_data_ilines = _bulk_data_lines_extract(
                lines, ilines, bulk_data_lines, i-1,
                make_ilines=make_ilines, keep_enddata=keep_enddata)
            break

        elif flag == 1:
            # handles ' CEND'
            if line_upper.startswith('CEND'):
                # case control
                old_flags.append(flag)
                if flag != 1:
                    raise RuntimeError('expected a flag of 1 (executive control deck) '
                                       'when going to the case control deck')

                flag = 2
                flag_word = 'case'
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

            # just reuse the existing one
            #line_upper = line.upper().strip()
            if line_upper.startswith('BEGIN'):
                if _is_begin_bulk(line_upper):
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
                    flag_word = 'bulk'
                    #print('case: %s' % (line.rstrip()))
                    case_control_lines.append(line.rstrip())
                    continue

                elif 'SUPER' in line_upper:  # and '=' in line_upper:
                    super_id = _get_super_id(line, line_upper)
                    log.info(f'super_id={super_id}')
                    old_flags.append(flag)
                    flag = -super_id
                    flag_word = 'SUPER=%s' % super_id
                    current_lines = superelement_lines[super_id]
                    current_ilines = superelement_ilines[super_id]
                #elif 'SUPER' in line_upper:
                    #super_id = -1
                    #flag = -1
                    #flag_word = 'SUPER'
                    #current_lines = superelement_lines[super_id]
                    #current_ilines = superelement_ilines[super_id]

                elif ('AUXMODEL' in line_upper or 'AFPM' in line_upper) and '=' in line_upper:
                    out = _read_bulk_for_model(
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
                else:  # pragma: no cover
                    msg = f'expected "BEGIN BULK" or "BEGIN SUPER=1"\nline = {line}'
                    raise RuntimeError(msg)

                #print('%s: %s' % (flag_word, line.rstrip()))
                current_lines.append(line.rstrip())
            elif line_upper.startswith('SUPER'):
                # case control line
                # SUPER = ALL
                #auxmodel_idi = int(line_upper.split('=')[1])
                #auxmodels_to_find.append(auxmodel_idi)
                if flag != 2:
                    raise RuntimeError('expected a flag of 2 (case control deck) '
                                       'when going to an SUPER model')

                is_superelement = True
            elif line_upper.startswith('AUXMODEL'):
                # case control line
                # AUXMODEL = 10
                auxmodel_idi = int(line_upper.split('=')[1])
                auxmodels_to_find.append(auxmodel_idi)
                if flag != 2:
                    raise RuntimeError('expected a flag of 2 (case control deck) '
                                       'when going to an AUXMODEL')
                is_auxmodel = True
            elif line_upper.startswith('AFPM'):
                # case control line
                # AFPM = 10
                afpm_idi = int(line_upper.split('=')[1])
                afpms_to_find.append(afpm_idi)
                if flag != 2:
                    raise RuntimeError('expected a flag of 2 (case control deck) '
                                       'when going to an AFPM model')
                is_afpm = True

            #print('%s: %s' % (flag_word, line.rstrip()))
            current_lines.append(line.rstrip())
        elif flag == 3:
            is_special_flag = (
                is_module is True or
                is_auxmodel is True or
                is_superelement is True or
                consider_superelements)
            if not is_special_flag:
                raise RuntimeError(f'one must be True: is_auxmodel={is_auxmodel}; '
                                   f'is_superelement={is_superelement}; '
                                   f'consider_superelements={consider_superelements}')
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

            out = _read_bulk_for_model(
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
        else:  # pragma: no cover
            raise RuntimeError(line)

    _check_valid_deck(flag, old_flags, nastran_format)

    if start_flag == 0 and len(bulk_data_lines) == 0:  # and flag != -1:
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
        dict(superelement_lines), dict(superelement_ilines),
        auxmodel_lines, afpm_lines,
    )
    return out


def _bulk_data_lines_extract(lines: list[str],
                             ilines: Any,
                             bulk_data_lines: list[str],
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


def _is_begin_bulk(line_upper: str) -> bool:
    """
    is this a:
      'BEGIN BULK'
    but not:
      'BEGIN BULK SUPER=2'
      'BEGIN BULK AUXMODEL=2'
      'BEGIN BULK AFPM=2'

    """
    is_begin_bulk = 'BULK' in line_upper and (
        'AUXMODEL' not in line_upper and
        'AFPM' not in line_upper and
        'SUPER' not in line_upper)
    return is_begin_bulk


def _read_bulk_for_model(ifile_iline, line: str, flag: int, bulk_data_lines: list[str],
                         current_lines: list[str], current_ilines,
                         old_flags: list[int],
                         unused_is_auxmodel: bool, auxmodel_lines: list[str], auxmodels_to_find, auxmodels_found,
                         unused_is_afpm: bool, afpm_lines: list[str], afpm_to_find, afpm_found,
                         superelement_lines: list[str], superelement_ilines,
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

    is_module_active = False
    module_lines: dict[int, int] = {}
    modules_found: set[int] = set()

    line_upper = line.upper().strip()
    if line_upper.startswith('BEGIN'):
        parse_begin(line_upper)
        if 'MODULE' in line_upper:
            is_module_active = True
            module_id, label = _get_module_id(line, line_upper)
            old_flags.append(flag)
            flag = module_id
            #current_lines = module_lines[(module_id, label)]
            #current_ilines = []
            #modules_found.add(module_id)
            raise NotImplementedError(line_upper)
            #if len(modules_found) == len(modules_to_find) and len(modules_found):
            #    #print('broken...final', len(bulk_data_lines), len(bulk_data_ilines))
            #    is_broken = True
            #    out = (is_broken,
            #           auxmodel_id, is_auxmodel_active,
            #           afpm_id, is_afpm_active,
            #           flag, current_lines)
            #    return out
        elif 'AUXMODEL' in line_upper:
            is_auxmodel_active = True
            auxmodel_id = _get_auxmodel_id(line, line_upper)
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
        elif 'SUPER' in line_upper:
            super_id = _get_super_id(line, line_upper)
            old_flags.append(flag)
            flag = -super_id
            current_lines = superelement_lines[super_id]
            current_ilines = superelement_ilines[super_id]
        elif 'AFPM' in line_upper:
            is_afpm_active = True
            afpm_id = _get_afpm_id(line, line_upper)
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


def _break_system_lines(executive_control_lines: list[str]) -> tuple[list[str],
                                                                     list[str]]:
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
    ID        Flag to name the run.
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


def _check_valid_deck(flag: int, old_flags: list[int],
                      nastran_format: str) -> None:
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

        msg = f'This is not a valid {nastran_format} BDF (a BDF capable of running Nastran).\n\n'
        msg += f'The following sections were found:\n{found}\n'
        msg += f'The following sections are missing:\n{missing}\n'
        msg += 'If you do not have an Executive Control Deck or a Case Control Deck:\n'
        msg += '  1.  call read_bdf(...) with `punch=True`\n'
        msg += "  2.  Add '$ pyNastran : punch=True' to the top of the main file\n"
        msg += '  3.  Name your file *.pch\n\n'
        msg += 'You cannot read a deck that has an Executive Control Deck, but\n'
        msg += 'not a Case Control Deck (or vice versa), even if you have a Bulk Data Deck.\n'
        raise MissingDeckSections(msg)
    return


def _show_bad_file(self: Any, bdf_filename: str | StringIO,
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
    lines: list[str] = []
    print('ENCODING - show_bad_file=%r' % encoding)

    with open(bdf_filename, 'r', encoding=encoding) as bdf_file:
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


def _get_module_id(line: str, line_upper: str) -> tuple[int, str]:
    """
    parses the module header::

        BEGIN MODULE=2  LABEL='MODULE2'

    """
    #if '=' in line_upper:
    sline = split_quoted_string(line_upper)
    assert len(sline) == 3, sline
    begin, module_num, label_name = sline
    assert begin == 'BEGIN', begin
    assert module_num.startswith('MODULE='), f'module_num={module_num!r}; line_upper={line_upper!r}'
    assert label_name.startswith("LABEL="), f'label_name={label_name!r}; line_upper={line_upper!r}'
    module_id_str = module_num.split('=')[1]
    module_id = int(module_id_str)

    label = label_name.split('=')[1].strip("'")
    print(f'module_id={module_id} label={label!r}')

    if module_id < 0:
        raise SyntaxError(f'module_id={module_id:d} must be greater than 0; line={line!r}')
    return module_id, label


def split_quoted_string(line: str) -> list[str]:
    """
    this is "a test"
    ['this','is','a test']
    """
    sline = shlex.split(line)
    return sline


def parse_begin(line_upper: str) -> None:
    """
    BEGIN MODULE=2  LABEL='MODULE2'

    BEGIN AUXMODEL=2
    BEGIN BULK AUXMODEL=2
    BEGIN BULK AUXMODEL = 2

    BEGIN AFPM=2
    BEGIN BULK AFPM=2
    BEGIN BULK AFPM = 2

    BEGIN SUPER=2
    BEGIN BULK SUPER=2
    BEGIN BULK SUPER = 2
    BEGIN BULK SUPER 2

    BEGIN BULK TRMC=101
    BEGIN TRMC=102

    """
    sline = split_quoted_string(line_upper.replace('=', ' '))

    #['BEGIN', 'MODULE=1', 'LABEL=MODULE1']
    #['BEGIN', 'MODULE', '1', 'LABEL', 'MODULE1']

    #['BEGIN', 'SUPER', '1']
    if sline[1] == 'BULK':
        sline.pop(1)
    if len(sline) == 1:
        # BEGIN BULK -> popped to BEGIN
        return

    word = sline[1]
    if word == 'MODULE':
        assert len(sline) == 5, sline
    elif word in ['SUPER', 'AFPM', 'AUXMODEL', 'TRMC']:
        assert len(sline) == 3, sline
    else:
        raise RuntimeError(sline)

    #line_upper = replace_multiple_spaces_with_single(line_upper)
    return


def _get_auxmodel_id(line: str, line_upper: str) -> int:
    """
    parses the superelement header::

        BEGIN AUXMODEL=2
        BEGIN BULK AUXMODEL=2
        BEGIN BULK AUXMODEL = 2

    """
    #if '=' in line_upper:
    sline = line_upper.split('=')
    #else:
        #sline = line_upper.split()
    try:
        auxmodel_id = int(sline[1])
    except (IndexError, ValueError):
        msg = 'expected "BEGIN AUXMODEL=1"\nline = %s' % line
        raise SyntaxError(msg)

    if auxmodel_id < 0:
        raise SyntaxError('auxmodel_id=%i must be greater than 0; line=%s' % (
            auxmodel_id, line))
    return auxmodel_id


def _get_afpm_id(line: str, line_upper: str) -> int:
    """
    parses the superelement header::

        BEGIN AFPM=2
        BEGIN BULK AFPM=2
        BEGIN BULK AFPM = 2

    """
    sline = line_upper.split('=')
    try:
        afpm_id = int(sline[1])
    except (IndexError, ValueError):
        msg = 'expected "BEGIN AFPM=1"\nline = %s' % line
        raise SyntaxError(msg)

    if afpm_id < 0:
        raise SyntaxError('afpm_id=%i must be greater than 0; line=%s' % (
            afpm_id, line))
    return afpm_id


def _get_super_id(line: str, line_upper: str) -> int:
    """
    parses the superelement header::

        BEGIN SUPER=2
        BEGIN BULK SUPER=2
        BEGIN BULK SUPER = 2
        BEGIN BULK SUPER 2

    """
    if '=' in line_upper:
        sline = line_upper.split('=')
        super_id_str = sline[1]
        if len(sline) != 2:
            msg = 'expected "BEGIN SUPER=1"\nline = %s' % line
            raise SyntaxError(msg)
    else:
        sline = line_upper.split()
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
    ilines[:, 1] = np.arange(nlines)  # 0 to N-1
    return ilines


def _check_for_spaces(card_name: str, card_lines: list[str], comment: str,
                      log: SimpleLogger) -> None:
    if ' ' in card_name:
        if card_name.startswith(EXECUTIVE_CASE_SPACES):  # TODO verify upper
            msg = (
                f'No spaces allowed in card name {card_name!r}.\n'
                'Did you mean to call read_bdf(punch=False) instead of '
                f'read_bdf(punch=True)?\n{card_lines}')
            raise RuntimeError(msg)
        elif card_name.startswith('BEGIN '):
            uline = card_lines[0].upper()
            if 'SUPER' in uline:
                msg = (
                    'Misindentified Superelement section.  Use:\n'
                    '$ pyNastran: is_superelements=True\n')
            else:
                msg = (
                    'Is there a second BEGIN BULK in your deck?\n'
                    'Another possibility is that punch=True and there is a '
                    'BEGIN BULK in your deck.\n')
            msg += '%s\n' % card_lines
            log.error(msg)
            raise SuperelementFlagError(msg)
        else:
            msg = (
                'No spaces allowed in card name %r.\n'
                'Should this be a comment?\n%s%s' % (
                    card_name, comment, card_lines))
        raise RuntimeError(msg)

    if card_name in ['SUBCASE ', 'CEND']:
        raise RuntimeError('No executive/case control deck was defined.')

#-------------------------------------------------------------------------------


DECK_TAGS = ('AFPM', 'ARBMODEL', 'AUXMODEL', 'MASSID', 'MODULE', 'FLXBDY', 'SUPER', 'TRMC', 'UDS')
ALLOW_LABEL = ('MASSID', 'MODULE')


def lines_to_decks2(lines: list[str],
                    ilines: NDArrayN2int,
                    punch: Optional[bool],
                    log: SimpleLogger,
                    nastran_format: str='msc') -> tuple[list[str],
                                                        list[str], list[str], list[str],
                                                        dict[tuple[str, str], list[str]],
                                                        dict[tuple[str, str], np.ndarray],]:
    """
    Splits the BDF lines into:
     - system lines
     - executive control deck
     - case control deck
     - bulk data deck

    Parameters
    ----------
    lines : list[str]
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
    system_lines : list[str]
        the system control lines (typically empty; used for alters)
    executive_control_lines : list[str]
        the executive control lines (stores SOL 101)
    case_control_lines : list[str]
        the case control lines (stores subcases)
    bulk_data_lines : list[str]
        the bulk data lines (stores geometry, boundary conditions, loads, etc.)
    bulk_data_ilines : None / (nlines, 2) int ndarray
        None : the old behavior
        narray : the [ifile, iline] pair for each line in the file
    superelement_lines : list[str]
        ???
    superelement_ilines : list[str]
        ???
    auxmodel_lines : list[str]
        ???
    auxmodel_ilines : dict
        ???

    """
    fake_ilines = ilines is None
    if fake_ilines:
        #log.warning('lines_to_decks2: fake_ilines')
        nlines = len(lines)
        ilines = np.zeros((nlines, 2), dtype='int32')
        ilines[:, 1] = np.arange(nlines)
    assert ilines is not None, ilines

    system_lines = []
    executive_control_lines = []   # flag = 1
    case_control_lines = []        # flag = 2
    bulk_data_lines = []           # flag = 3
    additional_deck_lines = {}
    additional_deck_ilines = {}
    if punch:  # True
        bulk_data_lines = lines
        bulk_data_ilines = ilines
        return (
            system_lines, executive_control_lines, case_control_lines,
            bulk_data_lines, bulk_data_ilines,
            additional_deck_lines, additional_deck_ilines,
        )

    # flag = -1 (dunno)
    flag = 'N/A'

    # AFPM = afpmid
    # ARBMODEL = arbmid
    # AUXMODEL = auxmind
    # MASSID = massidLABEL = masslabel
    # MODULE= moduleidAPPENDLABEL = modlabel
    # FLXBDY = flexbody
    # SUPER = seid
    # TRMC = trimid
    # UDS

    #AFPM = afpmid
    #ARBMODEL = arbmid
    #AUXMODEL = auxmind
    #MASSID = massid LABEL = masslabel
    #MODULE= moduleid APPENDLABEL = modlabel
    #FLXBDY = flexbody
    #SUPER = seid
    #TRMC = trimid
    #UDS

    #isol_line = -1
    #sol_line = ''
    #is_sol = False
    guess_deck_sections = punch is None
    flags = []
    executive_control_ilines_list = []
    case_control_ilines_list = []
    bulk_data_ilines_list = []
    active_ilines = executive_control_ilines_list
    active_lines = executive_control_lines
    ibulk0 = -1
    save_comment_flag = False
    for i, line in enumerate(lines):
        # handle empty comments
        file_iline = ilines[i]
        line = line.rstrip()
        if '$' in line:
            idollar = line.index('$')
            line_upper = line[:idollar].rstrip().upper()
            comment = line[idollar:].rstrip()
            #log.debug(f'comment={comment!r}')

            if len(line_upper) == 0:
                # don't drop the comment line if it goes in the bulk section
                # if it goes in the executive/case control, we need to drop it
                if comment and save_comment_flag:  # bulk or super
                    active_lines.append(comment)
                    active_ilines.append(file_iline)
                    #log.info(f'{i}: skipping {flag}: {line.rstrip()!r}')
                elif comment:
                    # log.debug(f'{i}: skipping {flag}: {line.rstrip()!r}')
                    # don't forget to reset it
                    comment = ''
                continue
        else:
            line_upper = line.rstrip().upper()
            if len(line_upper) == 0:  # no comment
                continue
        #----------------------------------------------------------------------
        #log.info(f'i={i} flag={flag!r} guess={guess_deck_sections} line={line_upper!r}')
        #assert len(bulk_data_lines) == len(bulk_data_ilines_list)
        #if line_upper.startswith('SOL '):
            #isol_line = i+1
            #sol_line = line
            #assert flag == -1, flag
            #executive_control_lines.append(line)
        #if line_upper.startswith(FILE_MANAGEMENT):
            #assert flag == -1, flag
            #executive_control_lines.append(line)
        if line_upper.startswith('CEND'):
            #log.debug('CEND')
            # start of case control deck
            #executive_control_lines = lines[:i]
            #log.warning(f'{i}: found CEND -> case_control')
            #active_lines.append(line)
            if flag != 'N/A':
                warnings.warn('No Executive Control Deck')
            flag = 'case_control'
            #active_lines.append(line)
            active_lines = case_control_lines
            active_ilines = case_control_ilines_list
            save_comment_flag = False

        elif line_upper.startswith('BEGIN'):
            #log.debug(f'***BEGIN; {i}')
            #active_lines.append(line)
            begin_tags = _get_begin_flag(line_upper)
            assert len(begin_tags) == 1, begin_tags
            begin_tag = begin_tags[0]
            begin_flag, idi, label = begin_tag
            assert isinstance(idi, int), (begin_flag, idi, label)

            begin_flag_old = f'{begin_flag}={idi:d}'

            if ibulk0 != -1:  # save last set of lines
                #if flag == 'bulk':
                # log.debug(f'flag={flag!r}; begin_tag={begin_tag!r}')
                #log.info(f'flag={flag!r}; begin_tag={begin_tag!r}')
                additional_deck_ilines[flag] = ilines[ibulk0:i]
                #log.info(f'additional_deck_ilines.keys = {list(additional_deck_ilines.keys())}')
            save_comment_flag = True
            ibulk0 = i

            #log.warning(f'{i}: ibulk0={ibulk0}; found begin_flag_old={begin_flag_old!r} -> begin_flag={begin_flag!r}')
            if begin_flag == 'BULK':
                assert flag in {'N/A', 'case_control'} or flag.startswith(DECK_TAGS), flag
                flag = 'bulk'
                active_lines = bulk_data_lines
                active_ilines = bulk_data_ilines_list
            elif begin_flag.startswith(DECK_TAGS):
                log.info(f'{i}: found deck_tag -> {begin_tag}')
                #tag = (begin_flag, idi, label)
                active_lines = []
                active_ilines = []
                additional_deck_lines[begin_tag] = active_lines
                additional_deck_ilines[begin_tag] = active_ilines
                flag = begin_flag_old
            else:  # pragma: no cover
                raise RuntimeError(begin_flag)
            flags.append(' '.join(str(value) for value in begin_tag))
            #log.info(f'begin_tag={begin_tag!r}')
            continue

        elif guess_deck_sections and _is_bulk_data_line(line) and flag in {'N/A', 'case_control'}:  # and flag in [1, 2]
            #log.debug('bulk')
            section_name = flag
            log.warning(f'currently in {section_name} deck and skipping directly '
                        f'to bulk data section\n{line}')

            flag = 'bulk'
            save_comment_flag = True
            active_lines = bulk_data_lines
            active_ilines = bulk_data_ilines_list
            active_ilines.append(file_iline)
            active_lines.append(line)
        else:
            #log.debug('else')
            #is_flag = flag in {'N/A', 'case_control'}
            #log.debug(f'  line={line!r}')
            #log.debug(f'  guess_deck_sections={guess_deck_sections}')
            #log.debug(f'  _is_bulk_data_line(line)={_is_bulk_data_line(line)}')
            #log.debug(f'  flag={flag!r} in N/A, case_control={is_flag}')
            #mytag = line_upper.startswith('BEGIN')
            #log.warning(f'{i}: {flag} -> {line_upper!r}; {mytag}')
            active_ilines.append(file_iline)
            active_lines.append(line)

        assert isinstance(flag, str), flag
        #print(flag)
        flags.append(flag)
        #log.debug('end; unhandled?')

    bulk_data_ilines = np.zeros((0, 2), dtype='int32')
    if len(bulk_data_ilines_list):
        bulk_data_ilines = np.vstack(bulk_data_ilines_list)
        assert len(bulk_data_lines) == len(bulk_data_ilines)
        del bulk_data_ilines_list

    if ibulk0 != -1:
        #if flag == 'bulk':
            #log.info(f'flag={flag!r}; begin_tag={begin_tag!r}')
        additional_deck_ilines[flag] = ilines[ibulk0:]

    #bulk_data_ilines = None
    bulk_data_flag = ('BULK', 0, '')
    #log.info(f'additional_deck_ilines.keys = {list(additional_deck_ilines.keys())}')

    if 'bulk' in additional_deck_ilines:
        bulk_data_flag = 'bulk'
    if bulk_data_flag in additional_deck_ilines:
        #bulk_data_ilines = additional_deck_ilines[bulk_data_flag]
        assert len(bulk_data_lines) == len(bulk_data_ilines)
        assert bulk_data_ilines is not None, bulk_data_ilines
        del additional_deck_ilines[bulk_data_flag]

    assert len(bulk_data_lines) == len(bulk_data_ilines)
    if len(bulk_data_lines) != len(bulk_data_ilines):
        msg = 'len(bulk_data_lines)=%s len(bulk_data_ilines)=%s' % (
            len(bulk_data_lines), len(bulk_data_ilines))
        log.warning(msg)
    system_lines, executive_control_lines = _break_system_lines(executive_control_lines)

    # clean comments
    system_lines = [_clean_comment(line) for line in system_lines
                    if _clean_comment(line) is not None]
    executive_control_lines = [_clean_comment(line) for line in executive_control_lines
                               if _clean_comment(line) is not None]
    case_control_lines = [_clean_comment(line) for line in case_control_lines
                          if _clean_comment(line) is not None]

    #for line in bulk_data_lines:
        #print(line.rstrip())

    if fake_ilines:
        bulk_data_ilines = {}
        additional_deck_ilines = {}
    out = (
        #isol_line, sol_line,
        system_lines, executive_control_lines, case_control_lines,
        bulk_data_lines, bulk_data_ilines,
        additional_deck_lines, additional_deck_ilines,
    )
    return out


def split_words_by_spaces(line: str) -> list[str]:
    """
    A B C='def'
    ['A', 'B', 'C=', 'def']
    """
    words = []
    word = ''
    label = ''
    is_free = True
    for char in line:
        if is_free:
            if char == ' ':
                words.append(word)
                word = ''
            elif char == "'":
                assert len(label) == 0, label
                is_free = False
                label += char
                if word:
                    words.append(word)
                    word = ''
            elif char == '"':
                raise RuntimeError('expected single quotes')
            else:
                word += char
        else:
            if char == "'":
                is_free = True
                label += char
                words.append(label)
            elif char == '"':
                raise RuntimeError('expected single quotes')
            else:
                label += char
    if word:
        words.append(word)

    words2 = []
    iword = 0
    while iword < len(words):
        word = words[iword]
        if word.endswith('='):
            iword += 1
            word += words[iword]
        else:
            words2.append(word)
            word = ''
        iword += 1
    if word:
        words2.append(word)

    assert len(''.join(words)) == len(''.join(words2)), (words, words2)
    return words2


def _get_begin_flag(line_upper: str) -> list[tuple[str, int, str]]:
    """
    begin massid=1 label='cat dog'
    'BEGIN   SUPER=7 MASSID=300'
    """
    if line_upper in {'BEGIN', 'BEGINBULK', 'BEGIN BULKS'}:
        line_upper = 'BEGIN BULK'
    if line_upper == 'BEGIN BULK':
        flag = [('BULK', 0, '')]
        return flag

    #sline0 = shlex.split(line_upper)
    sline0 = split_words_by_spaces(line_upper)
    words1 = _remove_bulk_words(sline0)

    # split the words by =
    words2 = []
    for word in words1:
        if '=' in word:
            sline = word.split('=', 1)
            words2.extend(sline)
        else:
            words2.append(word)
    words3 = _remove_bulk_words(words2)

    i = 0
    word = ''
    active_key = ''
    out_words = []
    temp_words = []
    while i < len(words3):
        word = words3[i].strip()
        if word == '':
            i += 1
            continue
        if active_key and word.isdigit():
            active_key = ''
            temp_words[-2] = int(word)
            out_words.append(tuple(temp_words))
            temp_words = []
            i += 1
            continue
        if active_key == 'LABEL':
            active_key = ''
            strip_label = word.strip('" ').strip("'")
            out_wordi = out_words.pop()

            out_wordi = list(out_wordi)
            out_wordi[-1] = strip_label
            out_words.append(tuple(out_wordi))
            #out_words.extend(temp_words)
            i += 1
            continue

        if word in DECK_TAGS:
            assert len(temp_words) == 0, temp_words
            active_key = word
            temp_words = [word, 0, '']
        elif word == 'LABEL':
            assert len(temp_words) == 0, temp_words
            active_key = word
        elif word == 'APPEND':
            out_words.append(('APPEND', -1, ''))
        else:
            assert len(temp_words) == 0, temp_words
            raise RuntimeError(f'i={i}; word={word!r}; words={words3}')
        i += 1

    if active_key:
        #active_key = ''
        out_words.append(tuple(temp_words))
        #i += 1
        #continue

    if len(out_words) == 0:
        return [('BULK', 0, '')]
    write_tag(out_words)
    return out_words


def _remove_bulk_words(words, *args) -> list[str]:
    words2 = [word.strip(' =') for word in words
              if 'BEGIN' != word and 'BULK' != word and 'BULKDATA' != word and '=' != word]
    return words2

#def _is_begin_bulk2(line_upper: str) -> bool:
    #"""
    #is this a:
      #'BEGIN'
      #'BEGIN BULK'
    #but not:
      #'BEGIN SUPER'
      #'BEGIN BULK SUPER=2'
      #'BEGIN BULK AUXMODEL=2'
      #'BEGIN BULK AFPM=2'

    #"""
    #words = line_upper.split()
    #if words in ('BEGIN', ('BEGIN', 'BULK')):
        #return True
    #return False


def add_superelements_from_deck_lines(self,
                                      BDF,
                                      additional_deck_lines: dict[tuple[str, str], list[str]]) -> dict:
    for superelement_key, superelement_lines in sorted(additional_deck_lines.items()):
        #if
        #superelement_id_str = superelement_key[0].split('=')[1]
        #superelement_id = int(superelement_id_str)
        assert len(superelement_key) == 3, superelement_key
        #superelement_id = superelement_key[1]
        assert isinstance(superelement_lines, list), superelement_lines

        # hack to get rid of extra 'BEGIN SUPER=2' lines
        iminus = 0
        for line in superelement_lines:
            uline = line.upper()
            if not uline.startswith('BEGIN '):
                continue
            iminus += 1

        nlines = len(superelement_lines) - iminus
        model = BDF()
        if hasattr(self, 'is_strict_card_parser'):
            model.is_strict_card_parser = self.is_strict_card_parser
        model.active_filenames = self.active_filenames
        model.log = self.log
        model.punch = True

        superelement_ilines = np.zeros((nlines, 2), dtype='int32')  ## TODO: calculate this
        model._parse_all_cards(superelement_lines[iminus:], superelement_ilines)
        self.superelement_models[superelement_key] = model
        self.initial_superelement_models.append(superelement_key)


def write_tag(begin_tags: list[tuple[str, int, str]]) -> str:
    word = 'BEGIN'
    for begin_tag in begin_tags:
        name, value, label = begin_tag
        if name in DECK_TAGS:
            word += f' {name}={value:d}'
        elif name == 'APPEND':
            word += ' APPEND'
        else:
            raise RuntimeError(name)

        if label:
            word += f" LABEL='{label}'"
    return word


def _dump_file(bdf_filename: str,
               lines: list[str],
               i: int, encoding: Optional[str]) -> None:
    """
    Writes a BDF up to some failed line index

    Parameters
    ----------
    bdf_filename : str
        the bdf filename to dump
    lines : list[str]
        the entire list of lines
    i : int
        the last index to write

    """
    with open(bdf_filename, 'w', encoding=encoding) as crash_file:
        for line in lines[:i]:
            crash_file.write(line)
