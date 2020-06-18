"""
defines the F06Writer class and:
 - write_f06(...)
"""
#pylint: disable=W0201,C0301,C0111
import os
import sys
import copy
import getpass
from datetime import date
from collections import defaultdict
from traceback import print_exc
from typing import Optional, List, Dict, Union

import pyNastran
from pyNastran.op2.op2_interface.op2_f06_common import OP2_F06_Common
from pyNastran.op2.op2_interface.result_set import ResultSet

def make_stamp(title, today=None, build=None):
    if title is None:
        title = ''

    #lenghts = [7, 8, 5, 5, 3, 4, 4, 6, 9, 7, 8, 8]
    months = [' January', 'February', 'March', 'April', 'May', 'June',
              'July', 'August', 'September', 'October', 'November', 'December']
    if today is None:
        today = date.today()
        str_month = months[today.month - 1].upper()
        str_today = '%-9s %2s, %4s' % (str_month, today.day, today.year)
    else:
        (month, day, year) = today
        str_month = months[month - 1].upper()
        str_today = '%-9s %2s, %4s' % (str_month, day, year)
    #str_today = str_today  #.strip()

    #release_date = '02/08/12'  # pyNastran.__releaseDate__
    if build is None:
        release_date = ''
        build = 'pyNastran v%s %s' % (pyNastran.__version__, release_date)
    out = '1    %-67s   %-19s %-22s PAGE %%5i\n' % (title.strip(), str_today, build)
    return out


def make_f06_header() -> str:
    spaces = ''
    lines1 = [
        spaces + '/* -------------------------------------------------------------------  */\n',
        spaces + '/*                              PYNASTRAN                               */\n',
        spaces + '/*                      - NASTRAN FILE INTERFACE -                      */\n',
        spaces + '/*                                                                      */\n',
        spaces + '/*              A Python reader/editor/writer for the various           */\n',
        spaces + '/*                        NASTRAN file formats.                         */\n',
        spaces + '/*                       Copyright (C) 2011-2020                        */\n',
        spaces + '/*                             Steven Doyle                             */\n',
        spaces + '/*                                                                      */\n',
        spaces + '/*    This program is free software; you can redistribute it and/or     */\n',
        spaces + '/*    modify it under the terms of the GNU Lesser General Public        */\n',
        spaces + '/*    License as published by the Free Software Foundation;             */\n',
        spaces + '/*    version 3 of the License.                                         */\n',
        spaces + '/*                                                                      */\n',
        spaces + '/*    This program is distributed in the hope that it will be useful,   */\n',
        spaces + '/*    but WITHOUT ANY WARRANTY; without even the implied warranty of    */\n',
        spaces + '/*    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the      */\n',
        spaces + '/*    GNU Lesser General Public License for more details.               */\n',
        spaces + '/*                                                                      */\n',
        spaces + '/*    You should have received a copy of the GNU Lesser General Public  */\n',
        spaces + '/*    License along with this program; if not, write to the             */\n',
        spaces + '/*    Free Software Foundation, Inc.,                                   */\n',
        spaces + '/*    675 Mass Ave, Cambridge, MA 02139, USA.                           */\n',
        spaces + '/* -------------------------------------------------------------------  */\n',
        '\n']

    spaces = 46 * ' '
    version = 'Version %8s' % pyNastran.__version__
    lines2 = [
        spaces + '* * * * * * * * * * * * * * * * * * * *\n',
        spaces + '* * * * * * * * * * * * * * * * * * * *\n',
        spaces + '* *                                 * *\n',
        spaces + '* *                                 * *\n',
        spaces + '* *                                 * *\n',
        spaces + '* *                                 * *\n',
        spaces + '* *            pyNastran            * *\n',
        spaces + '* *                                 * *\n',
        spaces + '* *                                 * *\n',
        spaces + '* *                                 * *\n',
        spaces + '* *%s* *\n' % version.center(33),
        spaces + '* *                                 * *\n',
        spaces + '* *                                 * *\n',
        spaces + '* *          %15s        * *\n' % pyNastran.__releaseDate2__,
        spaces + '* *                                 * *\n',
        spaces + '* *            Questions            * *\n',
        spaces + '* *        mesheb82@gmail.com       * *\n',
        spaces + '* *                                 * *\n',
        spaces + '* *                                 * *\n',
        spaces + '* *                                 * *\n',
        spaces + '* * * * * * * * * * * * * * * * * * * *\n',
        spaces + '* * * * * * * * * * * * * * * * * * * *\n\n\n']
    return ''.join(lines1 + lines2)


def sorted_bulk_data_header() -> str:
    """creates the bulk data echo header"""
    msg = '0                                                 S O R T E D   B U L K   D A T A   E C H O                                         \n'
    msg += '                 ENTRY                                                                                                              \n'
    msg += '                 COUNT        .   1  ..   2  ..   3  ..   4  ..   5  ..   6  ..   7  ..   8  ..   9  ..  10  .                      \n'
    return msg


def make_end(end_flag: bool=False,
             options: Optional[Dict[str, str]]=None) -> str:
    """creates the F06 footer"""
    lines = []
    lines2 = []
    if options is None:
        options = {}
    if end_flag:
        lines = [
            '', '',
            '0                                   * * * *  A N A L Y S I S  S U M M A R Y  T A B L E  * * * *',
            '0 SEID  PEID PROJ VERS APRCH      SEMG SEMR SEKR SELG SELR MODES DYNRED SOLLIN PVALID SOLNL LOOPID DESIGN CYCLE SENSITIVITY',
            ' --------------------------------------------------------------------------------------------------------------------------']
        #0     0    1    1 '        '    T    T    T    T    T     F      F      T      0     F     -1            0           F

        seid = 0
        peid = 0
        proj = 1
        vers = 1
        approach = "'        '"

        SELG = 'T'
        SEMG = 'T'
        SEMR = 'F'
        if 'SEMR' in options:
            SEMR = 'T' # modal
        SEKR = 'T'
        SELR = 'T'
        MODES = 'F'
        DYNRED = 'F'

        SOLLIN = 'F'
        if 'SEMR' in options:
            SOLLIN = 'T' # p-elements
        PVALID = 0
        SOLNL = 'F'
        LOOPID = -1
        CYCLE = 0
        SENSITIVITY = 'F'

        msg = '     %s     %s    %s    %s %8s    %s    %s    %s    %s    %s     %s      %s      %s      %s     %s     %s            %s           %s' % (
            seid, peid, proj, vers, approach, SEMG, SEMR, SEKR, SELG, SELR, MODES, DYNRED, SOLLIN, PVALID, SOLNL,
            LOOPID, CYCLE, SENSITIVITY)
        lines.append(msg)

        lines2 = [
            '0SEID = SUPERELEMENT ID.',
            ' PEID = PRIMARY SUPERELEMENT ID OF IMAGE SUPERELEMENT.',
            ' PROJ = PROJECT ID NUMBER.',
            ' VERS = VERSION ID.',
            ' APRCH = BLANK FOR STRUCTURAL ANALYSIS.  HEAT FOR HEAT TRANSFER ANALYSIS.',
            ' SEMG = STIFFNESS AND MASS MATRIX GENERATION STEP.',
            ' SEMR = MASS MATRIX REDUCTION STEP (INCLUDES EIGENVALUE SOLUTION FOR MODES).',
            ' SEKR = STIFFNESS MATRIX REDUCTION STEP.',
            ' SELG = LOAD MATRIX GENERATION STEP.',
            ' SELR = LOAD MATRIX REDUCTION STEP. ',
            ' MODES = T (TRUE) IF NORMAL MODES OR BUCKLING MODES CALCULATED.',
            ' DYNRED = T (TRUE) MEANS GENERALIZED DYNAMIC AND/OR COMPONENT MODE REDUCTION PERFORMED.',
            ' SOLLIN = T (TRUE) IF LINEAR SOLUTION EXISTS IN DATABASE.',
            ' PVALID = P-DISTRIBUTION ID OF P-VALUE FOR P-ELEMENTS',
            ' LOOPID = THE LAST LOOPID VALUE USED IN THE NONLINEAR ANALYSIS.  USEFUL FOR RESTARTS.',
            ' SOLNL = T (TRUE) IF NONLINEAR SOLUTION EXISTS IN DATABASE.',
            ' DESIGN CYCLE = THE LAST DESIGN CYCLE (ONLY VALID IN OPTIMIZATION).',
            ' SENSITIVITY = SENSITIVITY MATRIX GENERATION FLAG.',
            ' ',
            ' No PARAM values were set in the Control File.'
        ]

    lines3 = [
        ' ',
        '1                                        * * * END OF JOB * * *',
        ' ',
        ' '
    ]
    return '\n'.join(lines + lines2 + lines3)


class F06Writer(OP2_F06_Common):
    def __init__(self):
        OP2_F06_Common.__init__(self)
        self.card_count = {}
        self.additional_matrices = {}
        self.subcase_key = defaultdict(list)
        self.end_options = {}

        self._results = ResultSet(
            self.get_all_results(),
            self.op2_results._get_sum_objects_map(),
            self.log)

    def get_all_results(self) -> List[str]:
        all_results = ['stress', 'strain', 'element_forces', 'constraint_forces', 'thermal_load'] + self.get_table_types()
        return all_results

    def clear_results(self) -> None:
        self._results.clear()

    def add_results(self, results: Union[str, List[str]]) -> None:
        if isinstance(results, str):
            results = [results]
        all_results = self.get_all_results()
        for result in results:
            result = str(result)
            if result not in all_results:
                raise RuntimeError('%r is not a valid result to remove; all_results=%s' % (result, all_results))
            if result == 'stress':
                stress_results = []
                for result in all_results:
                    if 'stress' in result.lower():
                        stress_results.append(result)
                #stress_results = [result if 'stress' in result.lower() for result in all_results]
                self._results.update(stress_results)
            elif result == 'strain':
                strain_results = []
                for result in all_results:
                    if 'strain' in result.lower():
                        strain_results.append(result)
                #strain_results = [result if 'strain' in result.lower() for result in all_results]
                self._results.update(strain_results)
            elif 'stress' in result.lower():
                self._results.add('stress')
            elif 'strain' in result.lower():
                self._results.add('strain')
            elif result in ('spc_forces', 'mpc_forces', 'constraint_forces'):
                self._results.add('constraint_forces')
            elif 'force' in result.lower(): # could use more validation...
                self._results.add('element_forces')
            # thermalLoad_VU_3D, thermalLoad_1D, conv_thermal_load, thermalLoad_2D_3D
            self._results.add(result)

    def set_results(self, results: Union[str, List[str]]) -> None:
        self.clear_results()
        self.add_results(results)

    def remove_results(self, results: Union[str, List[str]]) -> None:
        self._results.remove(results)

    def make_f06_header(self):
        """If this class is inherited, the F06 Header may be overwritten"""
        return make_f06_header()

    def make_stamp(self, title, today, build=None):
        """If this class is inherited, the PAGE stamp may be overwritten"""
        return make_stamp(title, today, build=None)

    def make_grid_point_singularity_table(self, failed):
        """
        creates a grid point singularity table

        Parameters
        ----------
        failed : List[(nid, component), ...]
            defines failed degrees of freedom
        """
        msg = ''
        if failed:
            msg += '0                                         G R I D   P O I N T   S I N G U L A R I T Y   T A B L E\n'
            msg += '0                             POINT    TYPE   FAILED      STIFFNESS       OLD USET           NEW USET\n'
            msg += '                               ID            DIRECTION      RATIO     EXCLUSIVE  UNION   EXCLUSIVE  UNION\n'
            for (nid, dof) in failed:
                msg += '                         %8s        G      %s         0.00E+00          B        F         SB       SB   *\n' % (nid, dof)
        else:
            msg += 'No constraints have been applied...\n'

        page_stamp = self.make_stamp(self.title, self.date)
        msg += page_stamp % self.page_num
        self.page_num += 1
        return msg

    def _write_summary(self, f06_file, card_count=None):
        """writes the F06 card summary table"""
        summary_header = '                                        M O D E L   S U M M A R Y\n\n'
        summary = ''

        self.cards_to_read = {

            # rigid elements
            'RBAR', 'RBAR1', 'RBE1', 'RBE2', 'RBE3',

            # spc/mpc constraints
            'SPC', 'SPCADD', 'SPC1', 'SPCD', 'SPCAX',
            'MPC', 'MPCADD',
            'SUPORT', 'SUPORT1',

            # aero cards
            'CAERO1', 'CAERO2', 'CAERO3', 'CAERO4', 'CAERO5',

            # temperature cards
            'CHBDYE', 'CHBDYG', 'CHBDYP',
            'CONV',
        }


        blocks = [
            ['POINTS', ['GRID', 'GRDSET', ]],
            ['ENTRIES', ['SPOINT']],

            ['ELEMENTS',
             [
                 # these are sorted
                 # elements
                 'CONM1', 'CONM2', 'CMASS1', 'CMASS2', 'CMASS3', 'CMASS4',

                 # springs
                 'CELAS1', 'CELAS2', 'CELAS3', 'CELAS4', 'CELAS5',

                 # bushings
                 'CBUSH', 'CBUSH1D', 'CBUSH2D',

                 # dampers
                 'CDAMP1', 'CDAMP2', 'CDAMP3', 'CDAMP4', 'CDAMP5',

                 # bar flags
                 'BAROR', 'CBARAO',
                 # bars
                 'CBAR', 'CROD', 'CTUBE', 'BEAMOR', 'CBEAM', 'CBEAM3', 'CONROD', 'CBEND',

                 # shells
                 'CTRIA3', 'CTRIA6', 'CTRIAR', 'CTRIAX', 'CTRIAX6',
                 'CQUAD4', 'CQUAD8', 'CQUADR', 'CQUADX', 'CQUAD',

                 # solids
                 'CTETRA', 'CPENTA', 'CHEXA',

                 # other
                 'CSHEAR', 'CVISC', 'CRAC2D', 'CRAC3D',
                 'CGAP', 'CFAST', 'RBE2', 'RBE3',

                 # thermal
                 'CHBDYP', 'CHBDYG', 'CONV',
             ]],
        ]
        #print("self.card_count", self.card_count)
        if card_count is None:
            card_count = self.card_count

        for block in blocks:
            block_name, keys = block
            key_count = 0
            for key in sorted(keys):
                try:
                    value = card_count[key]
                    summary += '                                   NUMBER OF %-8s %-8s = %8s\n' % (key, block_name, value)
                    key_count += 1
                except KeyError:
                    pass
            if key_count:
                summary += ' \n'

        if summary:
            f06_file.write(summary_header)
            f06_file.write(summary)

            page_stamp = self.make_stamp(self.title, self.date)
            f06_file.write(page_stamp % self.page_num)
            self.page_num += 1

    def write_f06(self, f06_outname: str, matrix_filename: Optional[str]=None,
                  is_mag_phase: bool=False, is_sort1: bool=True,
                  delete_objects: bool=True, end_flag: bool=False,
                  quiet: bool=True, repr_check: bool=False,
                  close: bool=True) -> None:
        """
        Writes an F06 file based on the data we have stored in the object

        Parameters
        ----------
        f06_outname : str
            the name of the F06 file to write
        matrix_filename : str; default=None
            str : the name of the .mat file to write
            None : based on f06_outname
        is_mag_phase : bool; default=False
            should complex data be written using Magnitude/Phase
            instead of Real/Imaginary
            Real objects don't use this parameter
        is_sort1 : bool; default=True
            writes output in SORT1 format if the output is transient;
            ignored for static analyses
        delete_objects : bool; default=True
            should objects be deleted after they're written to reduce memory
        end_flag : bool; default=False
            should a dummy Nastran "END" table be made
        quiet : bool; default=False
            suppress print messages
        repr_check: bool; defualt=False
            calls the object repr as a validation test (prints nothing)
        close : bool; default=True
            close the f06 file
        """
        if not quiet:
            print("F06:")

        if isinstance(f06_outname, str):
            if matrix_filename is None:
                matrix_filename = os.path.splitext(f06_outname)[0] + '.mat'
            #print("matrix_filename =", matrix_filename)
            #mat = open(matrix_filename, 'wb')

            f06 = open(f06_outname, 'w')
            self._write_summary(f06)
        elif hasattr(f06_outname, 'read') and hasattr(f06_outname, 'write'):
            #f06 = f06_outname
        #else:
            #print('type(f06_outname) =', type(f06_outname))
            #assert isinstance(f06_outname, file), 'type(f06_outname)= %s' % f06_outname
            f06 = f06_outname
            f06_outname = f06.name
            if matrix_filename is None:
                matrix_filename = os.path.splitext(f06_outname)[0] + '.mat'
            if not quiet:
                print('f06_outname =', f06_outname)

        page_stamp = self.make_stamp(self.title, self.date)
        if self.grid_point_weight:
            if not quiet:
                print(" grid_point_weight")
            for key, weight in self.grid_point_weight.items():
                self.page_num = weight.write_f06(f06, page_stamp, self.page_num)

            if repr_check:
                str(self.grid_point_weight)
            assert isinstance(self.page_num, int), self.grid_point_weight.__class__.__name__

        if self.oload_resultant is not None:
            self.page_num = self.oload_resultant.write_f06(f06, page_stamp, self.page_num)
            if repr_check:
                str(self.oload_resultant)
            assert isinstance(self.page_num, int), self.oload_resultant.__class__.__name__

        # writes all results for
        self._write_f06_subcase_based(f06, page_stamp, delete_objects=delete_objects,
                                      is_mag_phase=is_mag_phase, is_sort1=is_sort1,
                                      quiet=quiet, repr_check=repr_check)

        self.op2_results.psds.write_f06(f06)

        #self._write_f06_time_based(f06, page_stamp)
        self.write_matrices(f06, matrix_filename, page_stamp, self.page_num, quiet=quiet)
        f06.write(make_end(end_flag, self.end_options))
        if close:
            f06.close()

    def write_matrices(self, f06, matrix_filename, page_stamp, page_num, quiet=True):
        """writes the f06 matrices"""
        if len(self.matrices):
            if hasattr(self, 'monitor1'):
                page_num = self.monitor1.write(f06, page_stamp=page_stamp, page_num=page_num)
                self.log.debug('MONPNT1 from [PMRF, PERF, PFRF, AGRF]')

            with open(matrix_filename, 'wb') as mat:
                for name, matrix in self.matrices.items():
                    if name == 'MP3F':
                        page_num = self.monitor3.write(f06, page_stamp=page_stamp, page_num=page_num)
                        self.log.debug('MONPNT3 from MP3F')
                    elif name in ['PMRF', 'PERF', 'PFRF', 'AGRF']:
                        pass
                    else:
                        if not quiet:
                            print(matrix)
                        matrix.write(mat)

    def _write_f06_subcase_based(self, f06, page_stamp, delete_objects=True,
                                 is_mag_phase=False, is_sort1=True, quiet=False,
                                 repr_check=False):
        """
        Helper function for ``write_f06`` that does the real work

        Parameters
        ----------
        f06 : file
            the opened file object
        page_stamp : str
            the format string stamp is the ending to every F06 page that
            contains the version, date, and page number
            (e.g., 'pyNastran 0.8   1/1/2016  PAGE %%i')
        is_mag_phase : bool; default=False
            should complex data be written using Magnitude/Phase
            instead of Real/Imaginary
            Real objects don't use this parameter
        is_sort1 : bool; default=True
            writes output in SORT1 format if the output is transient;
            ignored for static analyses
        delete_objects : bool; default=True
            should objects be deleted after they're written to reduce memory
        quiet : bool; default=False
            suppress print messages
        repr_check: bool; defualt=False
            calls the object repr as a validation test (prints nothing)
        """
        is_failed = False
        header = ['     DEFAULT                                                                                                                        \n',
                  '\n', '']

        # eigenvalues are written first
        f06.write(page_stamp % self.page_num)
        self.page_num += 1
        for ikey, result in sorted(self.eigenvalues.items()):
            if not quiet:
                print('%-18s case=%r' % (result.__class__.__name__, ikey))
            self.page_num = result.write_f06(f06, header, page_stamp,
                                             page_num=self.page_num)
            if repr_check:
                str(result)
            assert isinstance(self.page_num, int), 'pageNum=%r' % str(self.page_num)
            if delete_objects:
                del result
            self.page_num += 1

        # then eigenvectors
        # has a special header
        # isubcases = sorted(self.isubcase_name_map.keys())

        # TODO: superelement version...need the nominal...
        res_keys_subcase = self.subcase_key
        if len(res_keys_subcase) == 0:
            self.log.warning('no cases to write...self.subcase_key=%r' % self.subcase_key)
            return
        for isubcase, res_keys in sorted(res_keys_subcase.items()):
            for res_key in res_keys:
                if isinstance(res_key, tuple):
                    pass
                    #is_compressed = False
                else:
                    # int
                    #is_compressed = True
                    isubcase = res_key

                if res_key not in self.eigenvectors:
                    continue
                result = self.eigenvectors[res_key]
                if repr_check:
                    str(result)
                subtitle = result.subtitle
                header[0] = '     %s\n' % subtitle
                header[1] = '0                                                                                                            SUBCASE %i\n' % isubcase
                #header[2] = complex/nonlinear

                res_length = 18
                res_format = '*%%-%is SUBCASE=%%i' % res_length
                res_format_vectorized = ' %%-%is SUBCASE=%%i SUBTITLE=%%s' % res_length
                class_name = result.__class__.__name__
                if hasattr(result, 'data'):
                    if not quiet:
                        print(res_format_vectorized % (class_name, isubcase, subtitle))
                else:
                    print(res_format % (class_name, isubcase))

                self.page_num = result.write_f06(f06, header, page_stamp,
                                                 self.page_num, is_mag_phase=is_mag_phase, is_sort1=True)
                assert isinstance(self.page_num, int), 'pageNum=%r' % str(self.page_num)
                if delete_objects:
                    del result
                self.page_num += 1

        # finally, we writte all the other tables
        # nastran puts the tables in order of the Case Control deck,
        # but we're lazy so we just hardcode the order

        # subcase name, subcase ID, transient word & value
        header_old = ['     DEFAULT                                                                                                                        \n',
                      '\n', ' \n']
        header = copy.deepcopy(header_old)
        unallowed_results = ['eigenvectors', 'eigenvalues', 'params', 'gpdt', 'bgpdt', 'eqexin',
                             'grid_point_weight', 'psds']
        res_types = list(self.get_result(table_type) for table_type in sorted(self.get_table_types())
                         if table_type not in unallowed_results and not table_type.startswith('responses.'))

        for isubcase, res_keys in sorted(res_keys_subcase.items()):
            for res_key in res_keys:
                #if isinstance(res_key, tuple):
                    #is_compressed = False
                #else:
                    #is_compressed = True

                res_length = self._get_result_length(res_types, res_key)
                if res_length == 0:
                    # skipped subcase; no saved results
                    continue

                res_format = '*%%-%is SUBCASE=%%i%%s' % res_length
                res_format_vectorized = ' %%-%is SUBCASE=%%i SUBTITLE=%%s %%s' % res_length

                for res_type in res_types:
                    if res_key not in res_type:
                        continue

                    result = res_type[res_key]
                    if repr_check:
                        str(result)
                    subtitle = result.subtitle
                    label = result.label

                    header = ['', '']
                    header[0] = '      %-126s\n' % subtitle
                    header[1] = '0     %-32s                                                                       SUBCASE %-15i\n \n' % (label, isubcase)

                    if result.nonlinear_factor is not None:
                        header.append('')
                    try:
                        class_name = result.__class__.__name__

                        element_name = ''
                        if hasattr(result, 'element_name'):
                            element_name = ' - ' + result.element_name
                            is_ignored = 'StrainEnergy' not in class_name and 'GridPointForces' not in class_name
                            has_nnodes = not hasattr(result, 'nnodes_per_element')
                            if has_nnodes and is_ignored and getpass.getuser() == 'sdoyle':
                                self.log.error(f'{class_name} is missing nnodes_per_element')

                        if hasattr(result, 'data'):
                            if not quiet:
                                print(res_format_vectorized % (
                                    class_name, isubcase, subtitle, element_name))
                        else:
                            print(res_format % (class_name, isubcase, element_name))

                        result.is_complex
                        result.is_real

                        try:
                            self.page_num = result.write_f06(
                                f06, header, page_stamp, page_num=self.page_num,
                                is_mag_phase=is_mag_phase, is_sort1=is_sort1)
                        except Exception as error:
                            print_exc(file=sys.stdout)
                            print(''.join(result.get_stats()))
                            raise

                        #assert 'table_name=' in ''.join(result.get_stats())
                        assert isinstance(self.page_num, int), 'result=%s pageNum=%r' % (result, str(self.page_num))
                    except:
                        #print("result name = %r" % result.name())
                        raise
                    if delete_objects:
                        del result
                    self.page_num += 1
