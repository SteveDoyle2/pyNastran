"""
defines the CSVWriter class and:
 - write_csv(...)
"""
#pylint: disable=W0201,C0301,C0111
from __future__ import annotations
import os
import getpass
from typing import Optional, cast, TextIO, TYPE_CHECKING

import numpy as np

#import pyNastran
from pyNastran.op2.tables.oee_energy.oee_objects import RealStrainEnergyArray
from pyNastran.op2.tables.ogf_gridPointForces.ogf_objects import RealGridPointForcesArray
#from pyNastran.op2.tables.onmd import NormalizedMassDensity
#from pyNastran.op2.op2_interface.op2_f06_common import OP2_F06_Common
#from pyNastran.op2.op2_interface.result_set import ResultSet
from pyNastran.op2.result_objects.matrix import Matrix #, MatrixDict
if TYPE_CHECKING:
    from pyNastran.op2.op2 import OP2

def make_csv_header() -> str:
    #spaces = ''
    lines = """csv_version = 0.1

Key Table
---------
FLAG, ID, SubcaseID, Type, BLANK, <BLANK>, <BLANK>, <BLANK>, <BLANK>, <BLANK>, <BLANK>
0,   105,    1,       101,     0,       0,       0,       0,       0,       0,       0

Displacement Table
------------------
Flag, SubcaseID,  iTime, NID,       dx,      dy,       dz,      rx,       ry,      rz,  cd,  PointType
1,            1,  0,     101, 0.014159, 0.03448, 0.019135, 0.00637, 0.008042, 0.00762,   0,  1
# TODO: uses cd=-1 for unknown cd...requires geometry

Grid Point Forces Table
-----------------------
Flag, NID, SubcaseID, iTime, EID, TYPE,    Fx,      Fy,      Fz,      Mx,      My,      Mz
13,   101,         1,     1, 0,   APPLIED, 30.9864, 19.7278, 70.2515, 53.3872, 80.9687, 77.4302
13,   101,         1,     1, 301, RBE3,    41.9012, 53.6651, 0.09483, 76.041,  67.506,  98.0225
13,   101,         1,     1, 0,   SPC,     71.6306, 97.0527, 89.8733, 5.89262, 61.0523, 48.9043
13,   101,         1,     1, 301, CTRIA3,  84.5273, 69.36,   92.3295, 52.7074, 77.9904, 68.905
13,   101,         1,     1, 302, CQUAD4,  97.7843, 11.7545, 99.3901, 44.9476, 70.818,  7.47876
# TODO: add cd=-1 for unknown cd...requires geometry

Stress Table - PSHELL
---------------------
2,  stress,  ,  ,  ,  ,  ,  ,  ,  ,
,  E_type,  CQUAD4,  (PSHELL),  ,  ,  ,  ,  ,  ,
Flag, SubcaseID, iTime,  EID,  NID,      FD,      Sxx,        Syy,  Szz,       Sxy,  Syz,  Szx
2,            1,     0,  301,    0,   0.125,  265.173,   1535.666,    0,   169.811,    0,    0
2,            1,     0,  301,  101,   0.125,  62.7342,   1021.736,    0,  759.3948,    0,    0
2,            1,     0,  301,  102,   0.125,  12671707,  1078.741,    0,  1352.053,    0,    0
2,            1,     0,  301,  103,   0.125,  1797.972,  1944.449,    0,  719.5833,    0,    0
2,            1,     0,  301,  104,   0.125,  1109.873,  1651.793,    0,  1187.893,    0,    0
2,            1,     0,  301,    0,  -0.125,  389.2484,  1939.577,    0,  1270.192,    0,    0
2,            1,     0,  301,  101,  -0.125,  500.0936,  832.8021,    0,  1562.421,    0,    0
2,            1,     0,  301,  102,  -0.125,  1118429,   441.4652,    0,  303.9695,    0,    0
2,            1,     0,  301,  103,  -0.125,  530.9528,  116.0281,    0,  522.3726,    0,    0
2,            1,     0,  301,  104,  -0.125,  295.108,   1303.783,    0,  1465.683,    0,    0
# TODO: element or material coordinate system?

,  E_type,  CTRIA3,  (PSHELL),  ,  ,  ,  ,  ,  ,
Flag, SubcaseID, iTime, EID,  NID,      FD,       Sxx,       Syy,  Szz,       Sxy,  Syz,  Szx
2,            1,     0, 302,    0,  -0.125,  1289.590,  640.5084,    0,  1822.057,    0,    0
2,            1,     0, 302,    0,   0.125,  1851.625,  1957.094,    0,  1276.033,    0,    0
# TODO: element or material coordinate system?

Stress Table - PCOMP
--------------------
Flag,  SubcaseID, iTime, EID, Layer,      FD,      Sxx,      Syy, Szz,      Sxy, Syz, Szx
10,            1,     1, 301,     1, -0.1250, 1208.084, 290.0204,   0, 1594.263,   0,   0
10,            1,     1, 301,     2, -0.0625, 291.0046, 991.1379,   0, 916.3578,   0,   0
10,            1,     1, 301,     3,       0, 809.2431, 413.8515,   0, 966.4908,   0,   0
10,            1,     1, 301,     4,  0.0625, 443.0045, 1707.213,   0, 1897.417,   0,   0
10,            1,     1, 301,     5,   0.125, 370.4785, 1253.329,   0, 1221.529,   0,   0
# TODO: FD is just 0.0 for now...requires geometry
# TODO: element or material coordinate system?

,  E_type,  CTRIA3,  (PCOMP),  ,  ,  ,  ,  ,  ,
Flag,  SubcaseID, iTime, EID, Layer,      FD,      Sxx,      Syy, Szz,      Sxy, Syz, Szx
2,             1,     0, 302,     1,  -0.125,  189.859, 640.5084,   0, 1822.057,   0,   0
2,             1,     0, 302,     2, -0.0625, 125.6374, 404.5608,   0, 266.9504,   0,   0
2,             1,     0, 302,     3,       0, 16.93197, 1278.958,   0, 1397.945,   0,   0
2,             1,     0, 302,     4,  0.0625, 135.7806, 1048.905,   0, 1802.966,   0,   0
2,             1,     0, 302,     5,   0.125, 181.2625, 1957.094,   0, 1276.033,   0,   0
# TODO: FD is just 0.0 for now...requires geometry
# TODO: element or material coordinate system?

""".split('\n')

    lines2 = """
,  E_type,  CBAR,  ,  ,  ,  ,  ,  ,  ,
,  EID,  SubcaseID,  NID,  End,  S_axial,  S1,  S2,  S3,  S4,  <BLANK>
2,  303,  1,  101,  0,  1237.016753,  465.1125253,  1984.450331,  530.3698657,  1473.697261,  0
2,  303,  1,  102,  1,  1237.016753,  1862.546613,  358.6294203,  464.2576025,  1360.340506,  0
,  E_type,  CBUSH,  ,  ,  ,  ,  ,  ,  ,
,  EID,  SubcaseID,  <BLANK>,  <BLANK>,  Sxx,  Syy,  Szz,  Sxy,  Syz,  Szx
2,  304,  1,  0,  0,  724.0590054,  1012.171456,  982.1300666,  1830.004284,  1490.710065,  1694.931041
,  E_type,  CELAS1,  ,  ,  ,  ,  ,  ,  ,
,  EID,  SubcaseID,  <BLANK>,  <BLANK>,  Sxx,  <BLANK>,  <BLANK>,  <BLANK>,  <BLANK>,  <BLANK>
2,  305,  1,  0,  0,  724.0590054,  0,  0,  0,  0,  0
,  E_type,  CELAS2,  ,  ,  ,  ,  ,  ,  ,
,  EID,  SubcaseID,  <BLANK>,  <BLANK>,  Sxx,  <BLANK>,  <BLANK>,  <BLANK>,  <BLANK>,  <BLANK>
2,  306,  1,  0,  0,  724.0590054,  0,  0,  0,  0,  0
,  E_type,  CELAS3,  ,  ,  ,  ,  ,  ,  ,
,  EID,  SubcaseID,  <BLANK>,  <BLANK>,  Sxx,  <BLANK>,  <BLANK>,  <BLANK>,  <BLANK>,  <BLANK>
2,  307,  1,  0,  0,  724.0590054,  0,  0,  0,  0,  0
,  E_type,  CELAS4,  ,  ,  ,  ,  ,  ,  ,
,  EID,  SubcaseID,  <BLANK>,  <BLANK>,  Sxx,  <BLANK>,  <BLANK>,  <BLANK>,  <BLANK>,  <BLANK>
2,  308,  1,  0,  0,  724.0590054,  0,  0,  0,  0,  0
,  E_type,  CHEXA,  ,  ,  ,  ,  ,  ,  ,
,  EID,  SubcaseID,  NID,  <BLANK>,  Sxx,  Syy,  Szz,  Sxy,  Syz,  Szx
2,  309,  1,  0,  0,  1642.50393,  831.2196577,  1589.597775,  1676.545171,  1138.796188,  1213.684753
2,  309,  1,  101,  0,  524.3349312,  321.418164,  1553.325622,  1354.905702,  87.05448736,  520.19246
2,  309,  1,  102,  0,  1513.504576,  280.8425876,  253.6939871,  929.7662109,  1033.050714,  819.1099685
2,  309,  1,  103,  0,  1276.470632,  1721.405323,  417.1731873,  1301.056147,  374.7267389,  1372.351951
2,  309,  1,  104,  0,  271.1603737,  281.7424088,  506.6380522,  921.3316392,  783.4949079,  1796.671493
2,  309,  1,  105,  0,  1914.81801,  1833.386833,  1271.391108,  565.4908999,  1773.047636,  974.3552834
2,  309,  1,  106,  0,  930.7110598,  284.4549776,  1783.545275,  1573.502801,  1623.116685,  808.3474934
2,  309,  1,  107,  0,  416.5564825,  497.2385557,  1540.910353,  1777.870271,  1824.504737,  613.6814374
2,  309,  1,  108,  0,  877.779611,  1923.949045,  1211.281142,  64.96084846,  1387.493442,  1945.862149
,  E_type,  CTETRA,  ,  ,  ,  ,  ,  ,  ,
,  EID,  SubcaseID,  NID,  <BLANK>,  Sxx,  Syy,  Szz,  Sxy,  Syz,  Szx
2,  310,  1,  0,  0,  1642.50393,  831.2196577,  1589.597775,  1676.545171,  1138.796188,  1213.684753
2,  310,  1,  101,  0,  524.3349312,  321.418164,  1553.325622,  1354.905702,  87.05448736,  520.19246
2,  310,  1,  102,  0,  1513.504576,  280.8425876,  253.6939871,  929.7662109,  1033.050714,  819.1099685
2,  310,  1,  103,  0,  1276.470632,  1721.405323,  417.1731873,  1301.056147,  374.7267389,  1372.351951
2,  310,  1,  104,  0,  271.1603737,  281.7424088,  506.6380522,  921.3316392,  783.4949079,  1796.671493
,  E_type,  CPENTA,  ,  ,  ,  ,  ,  ,  ,
,  EID,  SubcaseID,  NID,  <BLANK>,  Sxx,  Syy,  Szz,  Sxy,  Syz,  Szx
2,  311,  1,  0,  0,  1642.50393,  831.2196577,  1589.597775,  1676.545171,  1138.796188,  1213.684753
2,  311,  1,  101,  0,  524.3349312,  321.418164,  1553.325622,  1354.905702,  87.05448736,  520.19246
2,  311,  1,  102,  0,  1513.504576,  280.8425876,  253.6939871,  929.7662109,  1033.050714,  819.1099685
2,  311,  1,  103,  0,  1276.470632,  1721.405323,  417.1731873,  1301.056147,  374.7267389,  1372.351951
2,  311,  1,  104,  0,  271.1603737,  281.7424088,  506.6380522,  921.3316392,  783.4949079,  1796.671493
2,  311,  1,  105,  0,  271.1603737,  281.7424088,  506.6380522,  921.3316392,  783.4949079,  1796.671493
2,  311,  1,  106,  0,  271.1603737,  281.7424088,  506.6380522,  921.3316392,  783.4949079,  1796.671493
,  E_type,  CROD,  ,  ,  ,  ,  ,  ,  ,
,  EID,  SubcaseID,  <BLANK>,  <BLANK>,  Sxx,  Syy,  Szz,  Sxy,  Syz,  Szx
2,  312,  1,  0,  0,  1642.50393,  0,  0,  1676.545171,  0,  0
,  E_type,  CONROD,  ,  ,  ,  ,  ,  ,  ,
,  EID,  SubcaseID,  <BLANK>,  <BLANK>,  Sxx,  Syy,  Szz,  Sxy,  Syz,  Szx
2,  313,  1,  0,  0,  1642.50393,  0,  0,  1676.545171,  0,  0
,  E_type,  CSHEAR,  ,  ,  ,  ,  ,  ,  ,
,  EID,  SubcaseID,  <BLANK>,  <BLANK>,  Sxx,  Syy,  Szz,  Sxy,  Syz,  Szx
2,  314,  1,  0,  0,  0,  0,  0,  1676.545171,  0,  0"""
    out = '#' + '\n#'.join(lines)
    return out

class CSVWriter:
    def __init__(self, op2: OP2):
        self.op2 = op2

    def write(self, csv_filename: str,
              matrix_filename: Optional[str]=None,
              is_exponent_format: bool=False,
              is_mag_phase: bool=False, is_sort1: bool=True,
              quiet: bool=True, repr_check: bool=False,
              close: bool=True) -> None:
        """
        Writes an F06 file based on the data we have stored in the object

        Parameters
        ----------
        csv_filename : str
            the name of the F06 file to write
        matrix_filename : str; default=None
            str : the name of the .mat file to write
            None : based on csv_filename
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
        repr_check: bool; default=False
            calls the object repr as a validation test (prints nothing)
        close : bool; default=True
            close the f06 file

        """
        model = self.op2
        log = model.log
        if not quiet:
            print("CSV:")

        csv, csv_filename, matrix_filename = _get_file_obj(
            csv_filename, matrix_filename, quiet=quiet)
        csv.write(make_csv_header())

        if model.grid_point_weight:
            if not quiet:
                print(" grid_point_weight")
            for key, weight in model.grid_point_weight.items():
                weight.write_csv(csv)

            if repr_check:
                str(model.grid_point_weight)

        #if self.oload_resultant is not None:
            #self.oload_resultant.write_csv(csv)
            #if repr_check:
                #str(self.oload_resultant)

        # writes all results for
        self._write_csv_subcase_based(csv,
                                      is_exponent_format=is_exponent_format,
                                      is_mag_phase=is_mag_phase, is_sort1=is_sort1,
                                      quiet=quiet, repr_check=repr_check)

        #model.op2_results.psds.write_csv(csv)
        model._write_normalized_mass_density(csv)

        #model._write_csv_time_based(csv)
        self.write_matrices(csv, matrix_filename, quiet=quiet)
        if close:
            csv.close()
        log.info(f'finished writing {csv.name}')
        return

    def write_matrices(self, csv: TextIO, matrix_filename: str,
                       quiet: bool=True):  # pragma: no cover
        """writes the f06 matrices"""
        model = self.op2
        if len(model.matrices) == 0:
            return
        log = model.log
        log.warning('write_matrices is not supported')
        return
        #raise NotImplementedError('write_matrices')
        results = model.op2_results
        #if results.monitor1 is not None:
            #results.monitor1.write(csv)
            #log.debug('MONPNT1 from [PMRF, PERF, PFRF, AGRF]')

        with open(matrix_filename, 'wb') as mat_file:
            for name, matrix in model.matrices.items():  #type: Matrix
                matrix = cast(Matrix, matrix)
                if name == 'MP3F':
                    results.monitor3.write(csv)
                    log.debug('MONPNT3 from MP3F')
                elif name in ['PMRF', 'PERF', 'PFRF', 'AGRF']:
                    pass
                else:
                    if not quiet:
                        print(matrix)
                    matrix.write(mat_file)

            responses = results.responses
            desvars = responses.desvars
            dscmcol = responses.dscmcol
            print(responses)
            if 'DSCM2' in model.matrices and desvars is not None and dscmcol is not None:
                data = matrix.data.todense()
                #print('dscmcol =', responses.dscmcol)
                row_data1 = ', '.join(str(val) for val in dscmcol.external_ids)
                row_data2 = ', '.join(str(val) for val in dscmcol.names)
                col_data = ', '.join(desvars.label)
                header = (
                    f'rows (external DRESPx ID): {row_data1}\n'
                    f'rows (external names): {row_data2}\n'
                    f'columns: {col_data}\n'
                )
                np.savetxt(mat_file, data, header=header, delimiter=',')

    def _write_csv_subcase_based(self, csv: TextIO,
                                 is_exponent_format: bool=False,
                                 is_mag_phase: bool=False,
                                 is_sort1: bool=True,
                                 quiet: bool=False,
                                 repr_check: bool=False):
        """
        Helper function for ``write_f06`` that does the real work

        Parameters
        ----------
        csv : file
            the opened file object
        is_mag_phase : bool; default=False
            should complex data be written using Magnitude/Phase
            instead of Real/Imaginary
            Real objects don't use this parameter
        is_sort1 : bool; default=True
            writes output in SORT1 format if the output is transient;
            ignored for static analyses
        quiet : bool; default=False
            suppress print messages
        repr_check: bool; default=False
            calls the object repr as a validation test (prints nothing)

        """
        model = self.op2
        log = model.log
        #is_failed = False
        header = ['     DEFAULT                                                                                                                        \n',
                  '\n', '']

        # eigenvalues are written first
        if 0:  # pragma: no cover
            for ikey, result in sorted(model.eigenvalues.items()):
                if not quiet:
                    print('%-18s case=%r' % (result.__class__.__name__, ikey))
                result.write_csv(csv, is_exponent_format=is_exponent_format)
                if repr_check:
                    str(result)

        # then eigenvectors
        # has a special header
        # isubcases = sorted(self.isubcase_name_map.keys())

        # TODO: superelement version...need the nominal...
        res_keys_subcase = model.subcase_key
        if len(res_keys_subcase) == 0:
            log.warning('no cases to write...subcase_key=%r' % model.subcase_key)
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

                if res_key not in model.eigenvectors:
                    continue
                result = model.eigenvectors[res_key]
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

                result.write_csv(csv, is_exponent_format=is_exponent_format,
                                 is_mag_phase=is_mag_phase, is_sort1=True)
                #check_element_node(result)

        # finally, we writte all the other tables
        # nastran puts the tables in order of the Case Control deck,
        # but we're lazy so we just hardcode the order

        # subcase name, subcase ID, transient word & value
        unallowed_results = ['eigenvectors', 'eigenvalues', 'params', 'gpdt', 'bgpdt', 'eqexin',
                             'grid_point_weight', 'psds', 'monitor1', 'monitor3']
        res_types = list(model.get_result(table_type) for table_type in sorted(model.get_table_types())
                         if table_type not in unallowed_results and not table_type.startswith('responses.'))

        for isubcase, res_keys in sorted(res_keys_subcase.items()):
            for res_key in res_keys:
                #if isinstance(res_key, tuple):
                    #is_compressed = False
                #else:
                    #is_compressed = True

                res_length = model._get_result_length(res_types, res_key)
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
                                log.error(f'{class_name} is missing nnodes_per_element')

                        if hasattr(result, 'data'):
                            if not quiet:
                                print(res_format_vectorized % (
                                    class_name, isubcase, subtitle, element_name))
                        else:
                            print(res_format % (class_name, isubcase, element_name))

                        result.is_complex
                        result.is_real

                        #check_element_node(result)
                        if not hasattr(result, 'write_csv'):
                            log.warning(f'missing write_csv for {result.class_name}')
                            continue
                        try:
                            result.write_csv(csv, is_exponent_format=is_exponent_format,
                                             is_mag_phase=is_mag_phase, is_sort1=is_sort1)
                        except Exception as error:
                            #print_exc(file=sys.stdout)
                            print(''.join(result.get_stats()))
                            #log.warning('')
                            #continue
                            raise

                        #assert 'table_name=' in ''.join(result.get_stats())
                    except Exception:
                        #print("result name = %r" % result.name())
                        raise

    def _write_normalized_mass_density(self, csv: TextIO):  # pragma: no cover
        model = self.op2
        normalized_mass_density = model.op2_results.responses.normalized_mass_density
        if normalized_mass_density is None:
            return
        normalized_mass_density0 = normalized_mass_density[0]  # type: NormalizedMassDensity

        csv.write('NORMALIZED MASS DENSITY HISTORY\n')
        for mass in normalized_mass_density0:
            csv.write(f'0    DESIGN_CYCLE={mass.dcycle:d} OBJ={mass.robj:g} RCON={mass.rcon:g}\n')

        for mass in normalized_mass_density0:
            csv.write('\nNORMALIZED MASS DENSITY\n')
            csv.write(f'0    DESIGN_CYCLE={mass.dcycle:d} OBJ={mass.robj:g} RCON={mass.rcon:g}\n')
            csv.write(' EID DENSITY\n')
            for eid, density in zip(mass.eids, mass.data):
                csv.write(f' {eid:-8d} {density:.8f}\n')


def _get_file_obj(csv_filename: str,
                  matrix_filename: Optional[str],
                  quiet: bool=True) -> tuple[TextIO, str, str]:
    if isinstance(csv_filename, str):
        if matrix_filename is None:
            matrix_filename = os.path.splitext(csv_filename)[0] + '.mat'
        #print("matrix_filename =", matrix_filename)
        #mat = open(matrix_filename, 'wb')

        csv = open(csv_filename, 'w')
    elif hasattr(csv_filename, 'read') and hasattr(csv_filename, 'write'):
        #f06 = f06_outname
    #else:
        #print('type(f06_outname) =', type(f06_outname))
        #assert isinstance(f06_outname, file), 'type(f06_outname)= %s' % f06_outname
        csv = csv_filename
        csv_filename = csv.name
        if matrix_filename is None:
            matrix_filename = os.path.splitext(csv_filename)[0] + '.mat'
        if not quiet:
            print('csv_filename =', csv_filename)
    return csv, csv_filename, matrix_filename

def write_csv(model: OP2, csv_filename: str, is_exponent_format: bool=True):
    csv = CSVWriter(model)
    csv.write(csv_filename, matrix_filename=None,
              is_exponent_format=is_exponent_format,
              is_mag_phase=False,
              is_sort1=True, quiet=True, repr_check=False, close=True)
