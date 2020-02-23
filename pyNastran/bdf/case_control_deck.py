# pylint: disable=R0904,R0902,C0103
"""
CaseControlDeck parsing and extraction class

CaseControlDeck:
----------------
    get_subcase_parameter(self, isubcase, param_name)
    has_subcase(self, isubcase)
    create_new_subcase(self, isubcase)
    delete_subcase(self, isubcase)
    copy_subcase(self, i_from_subcase, i_to_subcase, overwrite_subcase=True)
    get_subcase_list(self)
    get_local_subcase_list(self)
    update_solution(self, isubcase, sol)
    add_parameter_to_global_subcase(self, param)
    add_parameter_to_local_subcase(self, isubcase, param)
    finish_subcases(self)
    convert_to_sol_200(self, model)

"""
from __future__ import annotations
import re
import sys
import copy
from typing import List, Tuple, Dict, Any, Optional, TYPE_CHECKING

from cpylog import get_logger

#from pyNastran.bdf import subcase
from pyNastran.bdf.subcase import Subcase, update_param_name
from pyNastran.bdf.bdf_interface.subcase_cards import (
    EXTSEOUT, WEIGHTCHECK, GROUNDCHECK,
    MODCON, SET, SETMC, #AXISYMMETRIC,
    #INT_CARD_DICT, INT_CARD_NAMES,
    #INTSTR_CARD_DICT, INTSTR_CARD_NAMES,
    #STR_CARD_DICT, STR_CARD_NAMES,
    #CHECK_CARD_DICT, CHECK_CARD_NAMES,
    split_by_mixed_commas_parentheses,
)
from pyNastran.utils import object_attributes
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf import BDF


class CaseControlDeck:
    """CaseControlDeck parsing and extraction class"""
    type = 'CaseControlDeck'

    def __getstate__(self):
        # Copy the object's state from self.__dict__ which contains
        # all our instance attributes. Always use the dict.copy()
        # method to avoid modifying the original state.
        state = self.__dict__.copy()
        # Remove the unpicklable entries.
        del state['log']
        return state

    def __init__(self, lines: List[str], log: Optional[Any]=None) -> None:
        """
        Creates the CaseControlDeck from a set of lines

        Parameters
        ----------
        lines : List[str]
            list of lines that represent the case control deck
            ending with BEGIN BULK
        log : log()
            a :mod: `logging` object

        """
        # pulls the logger from the BDF object
        self.log = get_logger(log, level="debug")
        self.debug = False

        self.sol_200_map = {
            #101 - Linear Static
            #103 - Modal
            #105 - Buckling
            #106 - Non-Linear Static
            #107 - Direct Complex Eigenvalue
            #108 - Direct Frequency Response
            #109 - Direct Transient Response
            #110 - Modal Complex Eigenvalue
            #111 - Modal Frequency Response
            #112 - Modal Transient Response
            #129 - Nonlinear Transient
            #144 - Static Aeroelastic Analysis
            #145 - Flutter / Aeroservoelastic analysis
            #146 - Dynamic Aeroelastic Analysis
            #153 - Non-Linear static coupled with heat transfer
            #159 - Nonlinear Transient coupled with Heat transfer
            #187 - Dynamic Design Analysis Method
            #190 - DBTRANS Database Transfer
            #200 - Design Optimization and Sensitivity analysis
            #400 - Non-Linear Static and Dynamic (implicit) (MSC.NASTRAN native, supersedes 106, 129, 153 and 159 - part of MSC.NASTRAN)
            #401 - Non-Linear Static (SAMCEF based for NX.NASTRAN)
            #402 - Non-Linear Static and Dynamic (implicit) (SAMCEF based for NX.NASTRAN)
            #600 - Non-Linear Static and Dynamic (implicit) (front end to MSC.Marc - part of MSC.NASTRAN)
            #601 - Implicit Non-Linear (ADINA for NX Nastran, will no longer be available in NX NASTRAN after 2020)
            #700 - Explicit Non-Linear (LS Dyna plus MSC.Dytran - part of MSC.NASTRAN)
            #701 - Explicit Non-Linear (ADINA for NX Nastran, will no longer be available in NX NASTRAN after 2020)
            'STATICS' : 101,
            'STATIC' : 101,

            'MODES' : 103,
            'MODE' : 103,

            'BUCK' : 105,
            'BUCKLING' : 105,

            'DFREQ' : 108,
            'MFREQ' : 111,
            'SAERO' : 144,

            'FLUTTER' : 145,
            'FLUT' : 145,

            'DIVERGE' : 144,
            'DIVERG' : 145,

            # 'HEAT' : ,
            # 'STRUCTURE' : ,
            'NLSTATICS' : 400,
            'LNSTATICS' : 400,
            'MTRAN' : 112,
            'DCEIG' : 107,
        }
        # 'HEAT', 'ANALYSIS', 'MFREQ', 'STATICS', 'MODES', 'DFREQ',
        # 'MTRAN', 'BUCK', 'MCEIG', 'DCEIG', 'SAERO', 'NLSTATIC', 'NLSTAT',
        # 'STATIC', 'MTRANS', 'MODE', 'FLUTTER', 'DIVERG', 'NLTRAN', 'FLUT',
        #self.debug = True

        #: stores a single copy of 'BEGIN BULK' or 'BEGIN SUPER'
        self.reject_lines = []  # type: List[str]
        self.begin_bulk = ['BEGIN', 'BULK']

        # allows BEGIN BULK to be turned off
        self.write_begin_bulk = True
        self._begin_count = 0

        self.lines = lines
        self.subcases = {0: Subcase(id=0)}  # type: Dict[int, Subcase]
        try:
            self._read(self.lines)
        except:
            self.log.error('Invalid Case Control Deck:\n' + '\n'.join(self.lines))
            raise

    def load_hdf5_file(self, hdf5_file, encoding: str) -> None:
        """loads the case control deck section from a hdf5 file"""
        from pyNastran.utils.dict_to_h5py import _cast

        keys = list(hdf5_file.keys())
        for key in keys:
            if key in ['_begin_count', 'debug', 'write_begin_bulk', 'use_card_dict']: # scalars
                value = _cast(hdf5_file[key])
                setattr(self, key, value)
            elif key in ['reject_lines', 'begin_bulk', 'lines', 'output_lines']: # lists of strings
                value_bytes = _cast(hdf5_file[key]).tolist()
                unused_value_str = [line.decode(encoding) for line in value_bytes]
            elif key == 'subcases':
                subcase_group = hdf5_file[key]
                keys = list(subcase_group.keys())
                keys.remove('keys')
                subcases = {}
                for key2 in keys:
                    sub_group = subcase_group[key2]
                    ikey2 = int(key2)
                    subcase = Subcase(id=ikey2)
                    subcase.log = self.log
                    subcase.load_hdf5_file(sub_group, encoding)
                    subcases[ikey2] = subcase
                    str(subcase)
                self.subcases = subcases
                #print(value_bytes)
            else:  # pragma: no cover
                self.log.warning('skipping CaseControlDeck/%s' % key)
                raise RuntimeError('error loading hdf5 CaseControlDeck/%s' % key)

    def export_to_hdf5(self, hdf5_file, model: BDF, encoding: str) -> None:
        """exports the case control deck section to an hdf5 file"""
        keys_to_skip = ['type', 'log', 'sol_200_map', 'rsolmap_to_str', 'solmap_to_value']

        # scalars----
        #_begin_count
        #debug
        #write_begin_bulk

        # lines----
        #begin_bulk
        #lines
        #output_lines
        #reject_lines

        # subcases-----
        #subcases

        h5attrs = object_attributes(self, mode='both', keys_to_skip=keys_to_skip)
        for h5attr in h5attrs:
            value = getattr(self, h5attr)
            if h5attr in ['_begin_count', 'debug', 'write_begin_bulk']: # scalars
                # simple export
                hdf5_file.create_dataset(h5attr, data=value)
            elif h5attr in ['reject_lines', 'begin_bulk', 'lines', 'output_lines']:
                # lists of strings
                if len(value) == 0:
                    continue
                value_bytes = [line.encode(encoding) for line in value]
                #print(value_bytes)
                hdf5_file.create_dataset(h5attr, data=value_bytes)
            elif h5attr == 'subcases':
                keys = list(self.subcases.keys())
                subcase_group = hdf5_file.create_group('subcases')
                subcase_group.create_dataset('keys', data=keys)
                for key, subcase in self.subcases.items():
                    #print('***key =', key)
                    sub_group = subcase_group.create_group(str(key))
                    subcase.export_to_hdf5(sub_group, encoding)
            else:
                self.log.warning('skipping CaseControlDeck/%s' % h5attr)
                raise RuntimeError('error exporting hdf5 CaseControlDeck/%s' % h5attr)

    def suppress_output(self) -> None:
        """
        Replaces F06 printing with OP2 printing

        Converts:
            STRESS(PRINT,SORT1,REAL)
            FORCE(PRINT,PLOT,SORT1,REAL)

        to:
            STRESS(PLOT,SORT1,REAL)
            FORCE(PLOT,SORT1,REAL)

        .. warning:: most case control types are not supported

        """
        for subcase in self.subcases.values():
            # if isubcase == 0:
                # continue
            subcase.suppress_output()

    def has_parameter(self, isubcase: int, *param_names: List[str]) -> bool:
        """
        Checks to see if a parameter (e.g. STRESS) is defined in a certain
        subcase ID.

        Parameters
        ----------
        isubcase : int
            the subcase ID to check
        param_names : List[str]
            the parameter name to look for

        Returns
        -------
        has_parameter : bool
            does any subcase have a parameter
        """
        if self.has_subcase(isubcase):
            return any(self.subcases[isubcase].has_parameter(*param_names))
        return False

    def get_subcase_parameter(self, isubcase: int, param_name: str, obj: bool=False) -> Any:
        """
        Get the [value, options] of a subcase's parameter.  For example, for
        STRESS(PLOT,POST)=ALL:
            param_name=STRESS
            value=ALL
            options=['PLOT', 'POST']

        Parameters
        ----------
        isubcase : int
            the subcase ID to check
        param_name : str
            the parameter name to look for
        obj : bool; default=False
            should the object be returned

        """
        if self.has_subcase(isubcase):
            return self.subcases[isubcase].get_parameter(param_name.upper(), obj=obj)
        msg = ('isubcase=%r does not exist...subcases=%s'
               % (isubcase, str(sorted(self.subcases.keys()))))
        raise RuntimeError(msg)

    def has_subcase(self, isubcase: int) -> bool:
        """
        Checks to see if a subcase exists.

        Parameters
        ----------
        isubcase : int
            the subcase ID

        Returns
        -------
        val : bool
            does_subcase_exist (type = bool)

        """
        if isubcase in self.subcases:
            return True
        return False

    def create_new_subcase(self, isubcase: int) -> Subcase:
        """
        Method create_new_subcase:

        Parameters
        ----------
        isubcase : int
            the subcase ID

        Returns
        -------
        subcase : Subcase()
            the new subcase

        .. warning ::  be careful you dont add data to the global subcase
                       after running this...is this True???

        """
        #print("creating subcase=%s" % isubcase)
        if self.has_subcase(isubcase):
            sys.stderr.write('subcase=%s already exists...skipping\n' %
                             isubcase)
            return self.subcases[isubcase]
        subcase = self.copy_subcase(i_from_subcase=0, i_to_subcase=isubcase,
                                    overwrite_subcase=True)
        #self.subcases[isubcase] = Subcase(id=isubcase)
        return subcase

    def delete_subcase(self, isubcase: int) -> None:
        """
        Deletes a subcase.

        Parameters
        ----------
        isubcase : int
            the Subcase to delete

        """
        if not self.has_subcase(isubcase):
            sys.stderr.write('subcase %s doesnt exist...skipping\n' % isubcase)
        del self.subcases[isubcase]

    def copy_subcase(self, i_from_subcase: int, i_to_subcase: int,
                     overwrite_subcase: bool=True) -> Subcase:
        """
        Overwrites the parameters from one subcase to another.

        Parameters
        ----------
        i_from_subcase : int
            the Subcase to pull the data from
        i_to_subcase : int
            the Subcase to map the data to
        overwrite_subcase : bool; default=True
            NULLs i_to_subcase before copying i_from_subcase

        Returns
        -------
        subcase : Subcase()
            the new subcase

        """
        #print("copying subcase from=%s to=%s overwrite=%s" % (
            #i_from_subcase, i_to_subcase, overwrite_subcase))
        if not self.has_subcase(i_from_subcase):
            msg = 'i_from_subcase=%r does not exist...subcases=%s' % (
                i_from_subcase, str(sorted(self.subcases.keys())))
            raise RuntimeError(msg)
        if overwrite_subcase:
            subcase_from = self.subcases[i_from_subcase]
            subcase_to = copy.deepcopy(subcase_from)
            subcase_to.id = i_to_subcase
            #for key, param in sorted(subcase_from.params.items()):
                #print("going to copy key=%s param=%s" % (key, param))
            self.subcases[i_to_subcase] = subcase_to
        else:
            if not self.has_subcase(i_to_subcase):
                msg = ('i_from_subcase=%r does not exist...subcases=%s'
                       % (i_to_subcase, str(sorted(self.subcases.keys()))))
                raise RuntimeError(msg)
            subcase_to = self.subcases[i_to_subcase]

            for key, param in sorted(subcase_to.items()):
                #print('copying key=%s param=%s' % (key, param))
                if key == 'BEGIN':
                    pass
                subcase_to[key] = copy.deepcopy(param)
        return subcase_to

    def get_subcase_list(self) -> List[int]:
        """
        Gets the list of subcases including the global subcase ID (0)
        """
        return sorted(self.subcases.keys())

    def get_local_subcase_list(self) -> List[int]:
        """
        Gets the list of subcases that aren't the global subcase ID
        """
        id_list = [idi for idi in self.subcases if idi != 0]  # skip the global
        return sorted(id_list)

    def update_solution(self, isubcase: int, sol: str) -> None:
        """
        sol = STATICS, FLUTTER, MODAL, etc.

        Parameters
        ----------
        isubcase : int
            the subcase ID to update
        sol : str
            the solution type to change the solution to

        >>> bdf.case_control
        SUBCASE 1
            DISP = ALL

        >>> bdf.case_control.update_solution(1, 'FLUTTER')
        >>> bdf.case_control
        SUBCASE 1
            ANALYSIS FLUTTER
            DISP = ALL
        >>>

        """
        self.add_parameter_to_local_subcase(isubcase, 'ANALYSIS %s' % sol)

    def add_parameter_to_global_subcase(self, param):
        """
        Takes in a single-lined string and adds it to the global subcase.

        Parameters
        ----------
        param : str
            the variable to add

        Notes
        -----
        dont worry about overbounding the line

        Examples
        --------
        >>> bdf = BDF()
        >>> bdf.read_bdf(bdf_filename)
        >>> bdf.case_control.add_parameter_to_global_subcase('DISP=ALL')
        >>> bdf.case_control
        TITLE = DUMMY LINE
        DISP = ALL

        """
        (unused_j, key, value, options, param_type) = self._parse_data_from_user(param)
        subcase_list = self.get_subcase_list()
        for isubcase in subcase_list:
            self._add_parameter_to_subcase(key, value, options, param_type,
                                           isubcase)

    def add_parameter_to_local_subcase(self, isubcase: int, param: List[str]) -> None:
        """
        Takes in a single-lined string and adds it to a single Subcase.

        Parameters
        ----------
        isubcase : int
            the subcase ID to add
        param_name : List[str]
            the parameter name to add

        Notes
        -----
        dont worry about overbounding the line

        Examples
        --------
        >>> bdf = BDF()
        >>> bdf.read_bdf(bdf_filename)
        >>> bdf.case_control.add_parameter_to_local_subcase(1, 'DISP=ALL')
        >>> print(bdf.case_control)
        TITLE = DUMMY LINE
        SUBCASE 1
            DISP = ALL
        >>>

        """
        (unused_j, key, value, options, param_type) = self._parse_data_from_user(param)
        self._add_parameter_to_subcase(key, value, options, param_type,
                                       isubcase)

    def _parse_data_from_user(self, param):
        """
        Parses a case control line

        Parameters
        ----------
        param : str
            the variable to add

        """
        if '\n' in param or '\r' in param or '\t' in param:
            msg = 'doesnt support embedded endline/tab characters\n'
            msg += '  param = %r' % param
            raise SyntaxError(msg)
        #self.read([param])
        lines = _clean_lines([param])
        (j, key, value, options, param_type) = self._parse_entry(lines)
        return (j, key, value, options, param_type)

    def _read(self, lines: List[str]) -> None:
        """
        Reads the case control deck

        Notes
        -----
        supports comment lines

        .. warning:: doesnt check for 72 character width lines, but will
                     follow that when it's written out

        """
        isubcase = 0
        lines = _clean_lines(lines)
        self.output_lines = []
        i = 0
        #is_output_lines = False
        while i < len(lines):
            line = lines[i] #[:72]
            #comment = lines[i][72:]

            lines2 = [line]
            while ',' in lines[i][-1]:
                i += 1
                lines2.append(lines[i])
                #comment = lines[i][72:]
            (j, key, value, options, param_type) = self._parse_entry(lines2)
            i += 1

            line_upper = line.upper()
            if key == 'BEGIN':
                #if ('NLSTEP' in line_upper or 'SUBSTEP' in line_upper or 'SUBSEQ' in line_upper or
                    #'SUBCOM' in line_upper or 'REPCASE' in line_upper or 'NONLINEAR' in line_upper):
                    #asdfadsf
                if 'BULK' not in line_upper and 'SUPER' not in line_upper:
                    raise NotImplementedError('line=%r' % line)
                if self._begin_count == 1:
                    raise NotImplementedError('multiple BEGIN lines are defined...')
                self.begin_bulk = [key, value]
                self._begin_count += 1
                continue
            elif line_upper.startswith('OUTPUT'):
                #is_output_lines = True
                #output_line = '%s(%s) = %s\n' % (key, options, value)
                key = 'OUTPUT'

                # OUTPUT(POST) -> POST
                post = line_upper.split('OUTPUT')[1].strip('( )')
                options = [post]
                value = None
                param_type = 'STRESS-type'

                isubcase = self._add_parameter_to_subcase(
                    key, value, options, param_type, isubcase)
                self.output_lines.append(line)
                continue
            #print("key=%-12r icase=%i value=%r options=%r param_type=%r" % (
                #key, isubcase, value, options, param_type))
            isubcase = self._add_parameter_to_subcase(
                key, value, options, param_type, isubcase)

        #print(str(self))
        self.finish_subcases()

    def _parse_entry(self, lines):
        r"""
        Internal method for parsing a card of the case control deck

        Parses a single case control deck card into 4 sections

        1.  param_name - obvious
        2.  Value      - still kind of obvious
        3.  options    - rarely used data
        4.  param_type - STRESS-type, SUBCASE-type, PARAM-type, SET-type, BEGIN_BULK-type

        It's easier with examples:

        param_type = SUBCASE-type
          SUBCASE 1              ->   paramName=SUBCASE  value=1            options=[]
        param_type = STRESS-type
          STRESS       = ALL     ->   paramName=STRESS    value=ALL         options=[]
          STRAIN(PLOT) = 5       ->   paramName=STRAIN    value=5           options=[PLOT]
          TITLE        = stuff   ->   paramName=TITLE     value=stuff       options=[]
        param_type = SET-type
          SET 1 = 10,20,30       ->   paramName=SET       value=[10,20,30]  options = 1
        param_type = BEGIN_BULK-type
          BEGIN BULK             ->   paramName=BEGIN     value=BULK        options = []
        param_type = CSV-type
          PARAM,FIXEDB,-1        ->   paramName=PARAM     value=FIXEDB      options = [-1]

        The param_type is the "macro" form of the data (similar to integer, float, string).
        The value is generally whats on the RHS of the equals sign (assuming it's there).
        Options are modifiers on the data.  Form things like the PARAM card or the SET card
        they arent as clear, but the param_type lets the program know how to format it
        when writing it out.

        Parameters
        ----------
        lines : List[str, str, ...]
            list of lines

        Returns
        -------
        paramName : str
            see brief
        value : List[...]
            see brief
        options : List[str/int/float]
            see brief
        param_type : str/int/float/List
            see brief

        """
        i = 0
        options = []
        value = None
        key = None
        param_type = None

        line = lines[i]
        #print(line)
        #print("*****lines = ", lines)

        equals_count = 0
        for letter in line:
            if letter == '=':
                equals_count += 1
        line_upper = line.upper().expandtabs()

        #print("line_upper = %r" % line)
        #print('  equals_count = %s' % equals_count)
        if line_upper.startswith('SUBCASE'):
            param_type = split_equal_space(line_upper, 'SUBCASE', 'SUBCASE = 5')
            if ' ' in param_type:
                # SUBCASE 1 STATIC
                param_type = param_type.split()[0]
                self.log.debug(f"key={key!r} param_type={param_type!r} line_upper={line_upper!r}")

            key = 'SUBCASE'
            value = integer(param_type, line_upper)
            #self.isubcase = int(isubcase)
            param_type = 'SUBCASE-type'
            assert key.upper() == key, key

        elif line_upper.startswith(('LABEL', 'SUBT', 'TITL')):  # SUBTITLE/TITLE
            try:
                eindex = line.index('=')
            except ValueError:
                msg = "cannot find an = sign in LABEL/SUBTITLE/TITLE line\n"
                msg += "line = %r" % line_upper.strip()
                raise RuntimeError(msg)

            key = line_upper[0:eindex].strip()
            value = line[eindex + 1:].strip()
            options = []
            param_type = 'STRING-type'
        elif line_upper.startswith('SET ') and equals_count == 1:
            obj = SET.add_from_case_control(line_upper, lines, i)
            #if 0:
                #key = obj.key
                #options = None
                #value = obj
                #param_type = 'OBJ-type'
            #else:
            key = obj.key
            options = obj.set_id
            value = obj.value
            param_type = 'SET-type'
        elif line_upper.startswith('SETMC ') and equals_count == 1:
            obj = SETMC.add_from_case_control(line_upper, lines, i)
            #if 0:
                #key = obj.key
                #options = None
                #value = obj
                #param_type = 'OBJ-type'
            #else:
            key = value.key  # type: str
            options = obj.set_id  # type: List[int]
            value = obj.value  # type: int
            param_type = 'SET-type'

        #elif line_upper.startswith(CHECK_CARD_NAMES) and self.use_card_dict:
            #if '(' in line:
                #key = line_upper.strip().split('(', 1)[0].strip()
            #elif '=' in line:
                #key = line_upper.strip().split('=', 1)[0].strip()
            #else:
                #msg = 'expected item of form "name = value"   line=%r' % line.strip()
                #raise RuntimeError(msg)

            #key = update_param_name(key)
            #obj = CHECK_CARD_DICT[key].add_from_case_control(line, line_upper, lines, i)
            #value = obj.value
            #options = obj.options
            #param_type = 'STRESS-type'
            #key = obj.type

        #elif line_upper.startswith(INT_CARD_NAMES) and self.use_card_dict:
            #if '=' in line:
                #(name, value) = line_upper.strip().split('=')
            #else:
                #msg = 'expected item of form "name = value"   line=%r' % line.strip()
                #raise RuntimeError(msg)
            #name = update_param_name(name)
            #obj = INT_CARD_DICT[name].add_from_case_control(line, line_upper, lines, i)
            #key = obj.type
            ##if 0:
                ##value = obj.value
                ##options = []
                ##param_type = 'STRESS-type'
            ##else:
            #value = obj
            #options = None
            #param_type = 'OBJ-type'
            #key = obj.type

        #elif line_upper.startswith(INTSTR_CARD_NAMES) and self.use_card_dict:
            #if '=' in line:
                #(name, value) = line_upper.strip().split('=')
            #else:
                #msg = 'expected item of form "name = value"   line=%r' % line.strip()
                #raise RuntimeError(msg)
            #name = name.strip()
            #obj = INTSTR_CARD_DICT[name].add_from_case_control(line, line_upper, lines, i)
            #key = obj.type
            ##if 0:
                ##value = obj.value
                ##options = []
                ##param_type = 'STRESS-type'
            ##else:
            #value = obj
            #options = None
            #param_type = 'OBJ-type'
            #key = obj.type

        elif line_upper.startswith('EXTSEOUT'):
            options = None
            param_type = 'OBJ-type'
            obj = EXTSEOUT.add_from_case_control(line_upper.strip())
            value = obj
            key = obj.type
        elif line_upper.startswith('WEIGHTCHECK'):
            options = None
            param_type = 'OBJ-type'
            obj = WEIGHTCHECK.add_from_case_control(line, line_upper, lines, i)
            value = obj
            key = obj.type
        elif line_upper.startswith('GROUNDCHECK'):
            options = None
            param_type = 'OBJ-type'
            obj = GROUNDCHECK.add_from_case_control(line, line_upper, lines, i)
            value = obj
            key = obj.type
        elif line_upper.startswith('MODCON'):
            options = None
            param_type = 'OBJ-type'
            obj = MODCON.add_from_case_control(line, line_upper, lines, i)
            value = obj
            key = obj.type
        #elif line_upper.startswith('AUXMODEL'):
            #options = None
            #param_type = 'OBJ-type'
            #value = AUXMODEL.add_from_case_control(line, line_upper, lines, i)
            #key = value.type

        #elif line_upper.startswith(STR_CARD_NAMES):
            #if '=' in line:
                #(name, value) = line_upper.strip().split('=')
            #else:
                #msg = 'expected item of form "name = value"   line=%r' % line.strip()
                #raise RuntimeError(msg)
            #name = name.strip()
            #obj = STR_CARD_DICT[name].add_from_case_control(line, line_upper, lines, i)
            #value = obj
            #options = None
            #param_type = 'OBJ-type'
            #key = obj.type

        elif line_upper.startswith('TEMP'):
            if '=' in line:
                (key, value) = line_upper.strip().split('=')
            else:
                msg = 'expected item of form "name = value"   line=%r' % line.strip()
                raise RuntimeError(msg)
            assert equals_count == 1, line_upper
            if '(' in line_upper:
                options = None
                param_type = 'STRESS-type'
                key = key.strip().upper()
                value = value.strip()
                if self.debug:
                    self.log.debug("key=%r value=%r" % (key, value))
                param_type = 'STRESS-type'
                assert key.upper() == key, key

                #param_type = 'STRESS-type'
                sline = key.strip(')').split('(')
                key = sline[0]
                options = sline[1].split(',')

                assert len(options) == 1, line_upper
                # handle TEMPERATURE(INITIAL) and TEMPERATURE(LOAD) cards
                key = 'TEMPERATURE(%s)' % options[0]
                value = value.strip()
                options = []
            else:
                key = 'TEMPERATURE(BOTH)'
                options = []
                param_type = 'STRESS-type'
            value = int(value)

        elif line_upper.startswith('RIGID'):
            if '=' in line:
                (key, value) = line_upper.strip().split('=')
                key = key.strip()
                value = value.strip()
            else:
                msg = 'expected item of form "name = value"   line=%r' % line.strip()
                raise RuntimeError(msg)
            #print('line_upper=%r' % line_upper)
            assert key == 'RIGID', 'key=%r value=%r line=%r'  % (key, value, line)
            param_type = 'STRESS-type'
            options = []
            #RIGID = LAGR, LGELIM, LAGRANGE, STIFF, LINEAR
            if value in ['LAGR', 'LAGRAN']:
                value = 'LAGRANGE'
            elif value in ['LGELIM', 'LAGRANGE', 'STIFF', 'LINEAR', 'AUTO']:
                pass
            else:
                raise NotImplementedError('key=%r value=%r line=%r'  % (key, value, line))

        elif equals_count == 1:  # STRESS
            if '=' in line:
                (key, value) = line_upper.strip().split('=')
            else:
                msg = 'expected item of form "name = value"   line=%r' % line.strip()
                raise RuntimeError(msg)

            key = key.strip().upper()
            value = value.strip()
            if self.debug:
                self.log.debug("key=%r value=%r" % (key, value))
            param_type = 'STRESS-type'
            assert key.upper() == key, key

            if '(' in key:  # comma may be in line - STRESS-type
                #param_type = 'STRESS-type'
                sline = key.strip(')').split('(')
                key = sline[0]
                options = sline[1].split(',')

            elif ',' in value:  # STRESS-type; special TITLE = stuffA,stuffB
                #print('A ??? line = ',line)
                #raise RuntimeError(line)
                pass
            else:  # STRESS-type; TITLE = stuff
                #print('B ??? line = ',line)
                pass

            key = update_param_name(key)
            verify_card(key, value, options, line)
            assert key.upper() == key, key
        elif equals_count > 2 and '(' in line and 'FLSPOUT' not in line:
            #GROUNDCHECK(PRINT,SET=(G,N,N+AUTOSPC,F,A),DATAREC=NO)=YES
            #print('****', lines)
            assert len(lines) == 1, lines
            line = lines[0]
            try:
                key, value_options = line.split('(', 1)
            except ValueError:
                msg = 'Expected a "(", but did not find one.\n'
                msg += 'Looking for something of the form:\n'
                msg += '   GROUNDCHECK(PRINT,SET=(G,N,N+AUTOSPC,F,A),DATAREC=NO)=YES\n'
                msg += '%r' % line
                raise ValueError(msg)

            try:
                options_paren, value = value_options.rsplit('=', 1)
            except ValueError:
                msg = 'Expected a "=", but did not find one.\n'
                msg += 'Looking for something of the form:\n'
                msg += '   GROUNDCHECK(PRINT,SET=(G,N,N+AUTOSPC,F,A),DATAREC=NO)=YES\n'
                msg += 'value_options=%r\n' % value_options
                msg += '%r' % line
                raise ValueError(msg)
            options_paren = options_paren.strip()

            value = value.strip()
            if value.isdigit():
                value = int(value)
            if not options_paren.endswith(')'):
                raise RuntimeError(line)
            str_options = options_paren[:-1]

            if '(' in str_options:
                options = split_by_mixed_commas_parentheses(str_options)
            else:
                options = str_options.split(',')
            param_type = 'STRESS-type'
            key = key.upper()

        elif line_upper.startswith('BEGIN'):  # begin bulk
            try:
                (key, value) = line_upper.split(' ')
            except ValueError:
                msg = 'excepted "BEGIN BULK" found=%r' % line
                raise RuntimeError(msg)
            key = key.upper()
            param_type = 'BEGIN_BULK-type'
            assert key.upper() == key, key
        elif 'PARAM' in line_upper:  # param
            key, value, options, param_type = _split_param(line, line_upper)
        elif ' ' not in line:
            key = line.strip().upper()
            value = line.strip()
            options = None
            param_type = 'KEY-type'
            assert key.upper() == key, key
        else:
            msg = 'generic catch all...line=%r' % line
            key = ''
            value = line
            options = None
            param_type = 'KEY-type'
            assert key.upper() == key, key
        i += 1
        assert key.upper() == key, 'key=%s param_type=%s' % (key, param_type)

        return (i, key, value, options, param_type)

    def finish_subcases(self):
        """
        Removes any unwanted data in the subcase...specifically the SUBCASE
        data member.  Otherwise it will print out after a key like stress.

        """
        for subcase in self.subcases.values():
            subcase.finish_subcase()

    def convert_to_sol_200(self, model: BDF) -> None:
        """
        Takes a case control deck and changes it from a SOL xxx to a SOL 200

        Parameters
        ----------
        model : BDF()
            the BDF object

        .. todo:: not done...

        """
        analysis = model.rsolmap_to_str[model.sol]
        model.sol = 200

        subcase0 = self.subcases[0]
        subcase0.add_parameter_to_global_subcase('ANALYSIS', analysis)
        #subcase.add_parameter_to_global_subcase('DESSUB', dessub)

    def _add_parameter_to_subcase(self, key, value, options, param_type, isubcase):
        """Internal method"""
        if self.debug:
            a = 'key=%r' % key
            b = 'value=%r' % value
            c = 'options=%r' % options
            d = 'param_type=%r' % param_type
            msg = "_adding isubcase=%s %-12s %-12s %-12s %-12s" % (isubcase, a,
                                                                   b, c, d)
            self.log.debug(msg)

        if key == 'SUBCASE':
            assert value not in self.subcases, 'key=%s value=%s already exists' % (key, value)
            assert isinstance(value, int)
            isubcase = value
            self.copy_subcase(i_from_subcase=0, i_to_subcase=isubcase,
                              overwrite_subcase=True)
            if self.debug:
                msg = "copied subcase i_from_subcase=%r to i_to_subcase=%r" % (0, isubcase)
                self.log.debug(msg)
        elif isubcase not in self.subcases:  # initialize new subcase
            #self.isubcase += 1 # is handled in the read code
            msg = 'isubcase=%r is not a valid subcase...subcases=%s' % (
                isubcase, str(sorted(self.subcases.keys())))
            raise RuntimeError(msg)

        subcase = self.subcases[isubcase]
        subcase._add_data(key, value, options, param_type)

        #print("\n%s\n" % (self.subcases[isubcase]))
        return isubcase

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        for unused_isubcase, subcase in sorted(self.subcases.items()):
            subcase.cross_reference(model)

    def get_op2_data(self) -> Dict[int, Any]:
        """
        Gets the relevant op2 parameters required for a given subcase

        .. todo:: not done...
        """
        cases = {}
        for isubcase, subcase in sorted(self.subcases.items()):
            if isubcase:
                cases[isubcase] = subcase.get_op2_data(self.sol, subcase.solmap_to_value)
        return cases

    def __repr__(self) -> str:
        return self.write()

    def write(self, write_begin_bulk: Optional[bool]=None) -> str:
        """
        Writes the case control deck.  Has an option to not write the begin bulk line

        Parameters
        ----------
        write_begin_bulk : bool; default=None
            None: use the value in the original deck

        Returns
        -------
        msg : str
            the deck as a string
        """
        if write_begin_bulk is None:
            write_begin_bulk = self.write_begin_bulk
        msg = ''
        subcase0 = self.subcases[0]
        for unused_subcase_id, subcase in sorted(self.subcases.items()):
            msg += subcase.write_subcase(subcase0)
        #if len(self.subcases) == 1:
            #msg += 'BEGIN BULK\n'

        if self.output_lines:
            msg += '\n'.join(self.output_lines) + '\n'
        msg += '\n'.join(self.reject_lines)
        if write_begin_bulk:
            msg += ' '.join(self.begin_bulk) + '\n'
        return msg

def _split_param(line: str, line_upper: str) -> Tuple[str, str, str]:
    """parses a PARAM card"""
    tabbed_line_upper = line_upper.expandtabs().rstrip()
    if ',' in tabbed_line_upper:
        sline = tabbed_line_upper.split(',')
    elif '\t' in tabbed_line_upper:
        sline = tabbed_line_upper.expandtabs().split()
    elif ' ' in tabbed_line_upper:
        sline = tabbed_line_upper.split()
    else:
        raise SyntaxError("trying to parse %r..." % line)

    if len(sline) == 2:
        if ' ' in tabbed_line_upper:
            sline = tabbed_line_upper.replace(',', ' ').split()

    if len(sline) != 3:
        raise SyntaxError("trying to parse %r..." % line)
    (key, value, options) = sline
    param_type = 'CSV-type'
    assert key.upper() == key, key
    return key, value, options, param_type

def verify_card(key: int, value: Any, options: Any, line: str) -> None:
    """Make sure there are no obvious errors"""
    if key in ['AUXMODEL', 'BC', 'BCHANGE', 'BCMOVE', 'CAMPBELL', 'CLOAD',
               'CMETHOD', 'CSSCHD', 'DEACTEL', 'DEFORM', 'DESGLB', 'DESSUB',
               'DIVERG', 'DLOAD', 'DRSPAN', 'FMETHOD', 'FREQUENCY', 'GUST',
               'HADAPART', 'LINE', 'LOAD', 'LOADSET', 'MAXLINES', 'MCHSTAT',
               'MFLUID', 'MODES', 'MODTRAK', 'MPC', 'NLHARM',]:
        value2 = integer(value, line)
        assert value2 > 0, 'line=%r is invalid; value=%r must be greater than 0.' % (line, value2)

def verify_card2(key, value, options, line):
    """Make sure there are no obvious errors"""
    # this is purposely made overly strict to catch all the cases
    int_cards = [
        'SPC', 'MPC', 'TRIM', 'FMETHOD', 'METHOD', 'LOAD',
        'SUPORT', 'SUPORT1', 'TEMPERATURE(INITIAL)', 'TEMPERATURE(LOAD)',
        'DLOAD', 'MFLUID', 'CLOAD', 'NLPARM', 'CMETHOD',
        'FREQUENCY', 'TSTEP', 'TSTEPNL', 'SDAMPING', 'DESOBJ',
        'TEMPERATURE(INIT)', 'RANDOM', 'DESSUB', 'ADAPT', 'MAXLINES',
        'TFL', 'DESGLB', 'SMETHOD', 'DYNRED', 'GUST', 'TEMPERATURE(MATE)',
        'OTIME', 'NONLINEAR', 'AUXM', 'IC', 'BC', 'OUTRCV', 'DIVERG',
        'DATAREC', 'TEMPERATURE(BOTH)', 'DEFORM', 'MODES', 'CASE',
        'SEDR', 'SELG', 'SEFINAL', 'SEKR', 'TEMPERATURE(ESTIMATE)',
        'GPSDCON', 'AUXMODEL',
        'MODTRAK', 'OFREQ', 'DRSPAN', 'OMODES', 'ADACT', 'SERESP', 'STATSUB',
        'CURVESYM', 'ELSDCON', 'CSSCHD', 'NSM', 'TSTRU', 'RANDVAR',
        'RGYRO', 'SELR', 'TEMPERATURE(ESTI)', 'RCROSS', 'SERE', 'SEMR',
    ]

    # these may only be integers
    #print("key =", key)

    pass_headers = [
        'SUBTITLE', 'TITLE',
        'A2GG', 'M2GG', 'K2GG',
        'K2PP', 'M2PP',
        'K42GG',

        'XMIN', 'XMAX', 'XTITLE', 'XPAPE', 'XPAPER', 'XAXIS', 'XGRID', 'XGRID LINES', 'XLOG',
        'YMIN', 'YMAX', 'YTITLE', 'YPAPE', 'YPAPER', 'YAXIS', 'YGRID', 'YGRID LINES', 'YLOG',
        'XTMIN', 'XTMAX', 'XTGRID', 'XTTITLE', 'XTAXIS', 'XTGRID LINES', 'XTLOG',
        'YTMIN', 'YTMAX', 'YTGRID', 'YTTITLE', 'YTAXIS', 'YTGRID LINES', 'YTLOG',
        'XBMIN', 'XBMAX', 'XBGRID', 'XBAXIS', 'XBGRID LINES', 'XBTITLE', 'XBLOG',
        'YBMIN', 'YBMAX', 'YBGRID', 'YBAXIS', 'YBGRID LINES', 'YBTITLE', 'YBLOG',

        'RIGHT TICS', 'UPPER TICS',
        'TRIGHT TICS',
        'BRIGHT TICS',

        'PLOTTER', 'XYPLOT',

        'PTITLE',
        'HOUTPUT', 'PLOTID',
        'AXISYMMETRIC', 'CURVELINESYMBOL', 'CURVELINESYMB', 'AECONFIG',
        'B2GG', 'B2PP', 'AESYMXZ', 'TEMP', 'DSAPRT', 'MEFFMASS',
        'MAXMIN', 'RESVEC', 'MODESELECT', 'RIGID', 'TCURVE',
        'SUPER', 'MAXI DEFO', 'P2G',
        'EXTSEOUT', 'FLSTCNT PREFDB', 'AESYMXY',
        'DSYM',
    ]
    all_none_cards = [
        'STRESS', 'STRAIN', 'SPCFORCES', 'DISPLACEMENT', 'MPCFORCES', 'SVECTOR',
        'VELOCITY', 'ACCELERATION', 'FORCE', 'ESE', 'OLOAD', 'SEALL', 'GPFORCE',
        'GPSTRESS', 'GPSTRAIN', 'FLUX', 'AEROF', 'THERMAL', 'STRFIELD',
        'NOUTPUT', 'SEDV', 'APRES', 'HTFLOW', 'NLSTRESS', 'GPKE',
        'SACCELERATION', 'SDISPLACEMENT', 'SEMG', 'HARMONICS', 'PRESSURE', 'VUGRID',
        'ELSUM', 'SVELOCITY', 'STRFIELD REAL', 'SENSITY', 'MONITOR',
        'NLLOAD', 'GPSDCON', 'BOUTPUT',
    ]

    if key in ['BCONTACT', 'CURVELINESYMBOL']:
        value2 = integer(value, line)

    # these may only be integers greater than 0
    elif key in int_cards:
        value2 = integer(value, line)
        assert value2 > 0, 'line=%r is invalid; value=%r must be greater than 0.' % (line, value2)

    # these may have a value of all/none/integer, nothing else
    # except commas are allowed
    # 'DISP=ALL', 'DISP=NONE', 'DISP=1', 'DISP=1,2'
    elif key in all_none_cards:
        if value not in ['ALL', 'NONE']:
            if ',' in value:
                sline = value.split(',')
                for spot in sline:
                    value2 = integer(spot, line)
            else:
                value2 = integer(value, line)
                if value2 <= 0:
                    msg = 'line=%r is invalid; value=%r must be greater than 0.' % (line, value2)
                    raise ValueError(msg)
    elif key in ['ECHO']:
        #assert value in ['NONE','BOTH','UNSORT','SORT', 'NOSORT', 'PUNCH',
                         #''], 'line=%r is invalid; value=%r.' % (line, value)
        pass
    elif key in ['CSCALE', 'SUBSEQ', 'SYMSEQ', 'DEFORMATION SCALE', '']:
        # floats
        pass
    elif 'SET' in key:
        pass

    # weird cards
    elif key in pass_headers:
        pass
    elif key == 'ANALYSIS':
        assert value in ['HEAT', 'ANALYSIS', 'MFREQ', 'STATICS', 'MODES', 'DFREQ',
                         'MTRAN', 'BUCK', 'MCEIG', 'DCEIG', 'SAERO', 'NLSTATIC', 'NLSTAT',
                         'STATIC', 'MTRANS', 'MODE', 'FLUTTER', 'DIVERG', 'NLTRAN',
                         ], 'line=%r is invalid; value=%r' % (line, value)
    elif key == 'AUTOSPC':
        assert value in ['YES'], 'line=%r is invalid; value=%r' % (line, value)
    else:
        raise NotImplementedError('key=%r line=%r' % (key, line))


def _clean_lines(lines: List[str]) -> List[str]:
    """
    Removes comment characters defined by a *$*.

    Parameters
    ----------
    lines : List[str, ...]
        the lines to clean.

    """
    lines2 = []  # type: List[str]
    for line in lines:
        line = line.strip(' \n\r').split('$')[0].rstrip()
        if line:
            lines2.append(line)

    lines3 = []  # TODO: line, comment
    lines_pack = []
    for line in lines2:
        #print(line)
        if len(lines_pack) == 0:
            #print('0--', line)
            lines_pack.append(line)
            if not line.endswith(','):
                #print('next...')
                lines3.append(lines_pack)
                lines_pack = []
        elif line.endswith(','):
            #print('C--', line)
            lines_pack.append(line)
        else:
            if lines_pack[-1][-1] == ',':  # continued
                #print('xx--', line)
                lines_pack.append(line)
                lines3.append(lines_pack)
                #print('pack =', lines_pack)
                lines_pack = []
            else:  # new card
                #print('new--', line)
                lines3.append(lines_pack)
                lines_pack = [line]
    return [''.join(pack) for pack in lines3]

def split_equal_space(line: str, word: str, example: str) -> str:
    """
    Splits a case insensative line by an

    reads:
     - 'SUBCASE = 5'
     - 'SUBCASE 5'
    """
    if ' ' not in line and '=' not in line:
        raise SyntaxError("expected data of the form '%s', not %r" % (example, line))
    out = re.split(r'\s*%s\s*=?\s*' % word, line, maxsplit=1, flags=re.IGNORECASE)
    return out[1]

def integer(str_value: str, line: str) -> int:
    """casts the value as an integer"""
    try:
        value = int(str_value)
    except ValueError:
        raise ValueError('%r is not an integer; line:\n%r' % (str_value, line))
    return value
