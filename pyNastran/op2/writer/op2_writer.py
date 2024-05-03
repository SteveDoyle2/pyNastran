#pylint: disable=W0201,C0301,C0111
from __future__ import annotations
import datetime
from collections import defaultdict
from struct import pack, Struct
from typing import Optional, TYPE_CHECKING
from cpylog import get_logger2

#import pyNastran
from pyNastran.op2.op2_interface.op2_f06_common import OP2_F06_Common
from pyNastran.op2.op2_interface.write_utils import _write_markers
#from pyNastran.op2.errors import FatalError
from .case_writer import write_casecc
from .geom1_writer import write_geom1
from .geom2_writer import write_geom2
from .geom3_writer import write_geom3
from .geom4_writer import write_geom4
from .ept_writer import write_ept
from .mpt_writer import write_mpt
from .edt_writer import write_edt
from .edom_writer import write_edom
from .dit_writer import write_dit
#from .dynamic_writer import write_dynamic
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.op2.op2 import OP2


class TrashWriter:
    def __init__(self, *args, **kwargs):
        pass
    def write(self, data):
        pass
    def close(self):
        pass


class OP2Writer(OP2_F06_Common):
    def __init__(self, log=None, debug: bool=False):
        self.log = get_logger2(log, debug)
        OP2_F06_Common.__init__(self)
        self.card_count = {}

    def write_op2(self, op2_out_filename: str,
                  post: int=-1,
                  endian: bytes=b'<',
                  includes: Optional[list[str]]=None,
                  skips: Optional[list[str]]=None,
                  nastran_format: Optional[str]=None) -> int:
        """
        Writes an OP2 file based on the data we have stored in the object

        Parameters
        ----------
        op2_out_filename : str
            the name of the F06 file to write
        post : int; default=-1
            the PARAM,POST flag
        endian : bytes; default='<'
            little endian is strongly recommended
        includes : list[str]; default=None
            list of results to include; exclusive with skips; default for both includes/skips=None -> all included
        skips : list[str]; default=None
            list of results to skip; exclusive with includes
        nastran_format : str; default=None -> 'msc'
            supported formats: ['msc', 'nx', 'optistruct']
        #is_mag_phase : bool; default=False
            #should complex data be written using Magnitude/Phase
            #instead of Real/Imaginary (default=False; Real/Imag)
            #Real objects don't use this parameter.

        """
        if nastran_format is None:
            nastran_format = self._nastran_format
        assert nastran_format in {'msc', 'nx', 'optistruct'}, nastran_format
        skips = _set_skips(self, includes, skips)

        #print('writing %s' % op2_outname)

        if isinstance(op2_out_filename, str):
            op2_file = open(op2_out_filename, 'wb')
            #fop2_ascii = open(op2_outname + '.txt', 'w')
            fop2_ascii = TrashWriter()
            #print('op2 out = %r' % op2_outname)
            close = True
        else:
            assert isinstance(op2_out_filename, file), f'type(op2_out_filename) = {op2_out_filename}'
            op2_file = op2_out_filename
            op2_outname = op2_out_filename.name
            close = False
            #print('op2_outname =', op2_outname)

        try:
            total_case_count, table_names = _write_op2(
                op2_file, fop2_ascii, self,
                skips,
                post=post, endian=endian,
                nastran_format=nastran_format)
        except Exception:  # NotImplementedError
            if close:
                op2_file.close()
                fop2_ascii.close()
            raise
        return total_case_count

def _set_skips(model: OP2Writer, includes: Optional[list[str]], skips: Optional[list[str]]):
    """
    Helper method for ``write_op2``

    Parameters
    ----------
    includes : list[str]; default=None
        list of results to include; exclusive with skips; default for both includes/skips=None -> all included
    skips : list[str]; default=None
        list of results to skip; exclusive with includes

    """
    if skips is not None and includes is not None:
        raise RuntimeError('Neither skips nor includes is None')
    elif skips is None and includes is None:
        skips = set()
    elif skips is not None:
        skips = set(skips)
    elif includes is not None:
        includes = set(includes)

        all_tables = {
            'GEOM1', 'GEOM2', 'GEOM3', 'GEOM4', 'EPT', 'MPT', 'EDT', 'EDOM', 'DIT', 'DYNAMIC',
            'params', 'grid_point_weight',
            'eigenvalues', 'eigenvalues_fluid',

            'displacements', 'mpc_forces', 'spc_forces', 'load_vectors',
            'cbeam_force', 'cbeam_stress',
            'grid_point_forces',
        }
        bad_tables = includes - all_tables
        if bad_tables:
            model.log.error(f'The following tables are not part of "all_tables"\n{bad_tables}')
        skips = all_tables - includes
        #model.log.warning("includes is poorly supported!")
    else:
        print(f'skips={skips}')
        print(f'includes={includes}')
        raise RuntimeError('Can this happen?')

    #if skips is None and includes is None:
        #skips = set([])
    #if skips is None:
        #skips = set([])
    #else:
        #skips = set(skips)
    skips.add('params')
    return skips



def _write_op2(op2_file, fop2_ascii, obj: OP2,
               skips: set[str],
               post: int=-1, endian: bytes=b'<',
               nastran_format: str='nx') -> tuple[int, list[str]]:
    """actually writes the op2"""
    date = obj.date
    #op2_ascii.write('writing [3, 7, 0] header\n')

    struct_3i = Struct(endian + b'3i')
    write_op2_header(obj, op2_file, fop2_ascii, struct_3i, post=post, endian=endian)
    #if 'CASECC' not in skips:
        #write_casecc(op2_file, fop2_ascii, obj, endian=endian, nastran_format=nastran_format)
    obj.log.debug(f'nastran_format={nastran_format}')
    if 'GEOM1' not in skips:  # nodes
        write_geom1(op2_file, fop2_ascii, obj, endian=endian)
    if 'GEOM2' not in skips:  # elements
        write_geom2(op2_file, fop2_ascii, obj, endian=endian)
    if 'GEOM3' not in skips:  # constraints
        write_geom3(op2_file, fop2_ascii, obj, endian=endian, nastran_format=nastran_format)
    if 'GEOM4' not in skips:  # loads
        write_geom4(op2_file, fop2_ascii, obj, endian=endian, nastran_format=nastran_format)
    if 'EPT' not in skips:    # properties
        write_ept(op2_file, fop2_ascii, obj, endian=endian, nastran_format=nastran_format)
    if 'MPT' not in skips:    # materials
        write_mpt(op2_file, fop2_ascii, obj, endian=endian, nastran_format=nastran_format)

    if 'EDT' not in skips:  # aero
        write_edt(op2_file, fop2_ascii, obj, endian=endian, nastran_format=nastran_format)
    if 'EDOM' not in skips:  # optimization
        write_edom(op2_file, fop2_ascii, obj, endian=endian, nastran_format=nastran_format)
    if 'DIT' not in skips:  # tables
        write_dit(op2_file, fop2_ascii, obj, endian=endian, nastran_format=nastran_format)
    #if 'DYNAMIC' not in skips:
        #write_dynamic(op2_file, fop2_ascii, obj)
    if 'grid_point_weight' not in skips:
        for key, weight in obj.grid_point_weight.items():
            weight.write_op2(op2_file, fop2_ascii, date, endian=endian)

    #is_mag_phase = False
    # we write all the other tables
    # nastran puts the tables in order of the Case Control deck,
    # but we're lazy so we just hardcode the order

    case_count, table_names = _write_result_tables(obj, op2_file, fop2_ascii, struct_3i, endian, skips)
    return case_count, table_names

def _write_result_tables(obj: OP2, op2_file, fop2_ascii,
                         struct_3i,
                         endian, skips: set[str]) -> tuple[int, list[str]]:
    """writes the op2 result tables"""
    table_names_found = []
    date = obj.date
    log = obj.log
    res_categories2 = defaultdict(list)
    table_order = [
        'OUGV1', 'BOUGV1',
        'OUGF1', 'OUG1F', 'BOUGF1',
        'TOUGV1', 'OTEMP1',

        'OAG1', 'OVG1', 'OUG1', 'OUGV1PAT',

        # eigenvectors, basic eigenvectors, basic fluid eigenvectors
        'OPHIG', 'BOPHIG', 'BOPHIGF',
        # eigenvectors in SOL 401 - NX
        'OUXY1', 'OUXY2', 'OPHSA',

        # spc/mpc forces
        'OQG1',
        'OQGV1',
        'OQP1',
        'OQMG1',
        # contact/glue forces
        'OQGCF1', 'OQGGF1',
        'OBC1',
        # load vectors
        'OPG1', 'BOPG1', 'OPGV1', 'OPNL1',

        # ---------------
        # random displacement-style tables
        'OUGATO1', 'OUGCRM1', 'OUGNO1', 'OUGPSD1', 'OUGRMS1', # disp/vel/acc/eigenvector
        'OUGATO2', 'OUGCRM2', 'OUGNO2', 'OUGPSD2', 'OUGRMS2',

        'OVGATO1', 'OVGCRM1', 'OVGNO1', 'OVGPSD1', 'OVGRMS1', # velocity
        'OVGATO2', 'OVGCRM2', 'OVGNO2', 'OVGPSD2', 'OVGRMS2',

        'OAGATO1', 'OAGCRM1', 'OAGNO1', 'OAGPSD1', 'OAGRMS1', # acceleration
        'OAGATO2', 'OAGCRM2', 'OAGNO2', 'OAGPSD2', 'OAGRMS2',

        'OQGATO1', 'OQGCRM1', 'OQGNO1', 'OQGPSD1', 'OQGRMS1', # spc/mpc forces
        'OQGATO2', 'OQGCRM2', 'OQGNO2', 'OQGPSD2', 'OQGRMS2',

        'OQMATO1', 'OQMCRM1', 'OQMNO1', 'OQMPSD1', 'OQMRMS1', # mpc forces
        'OQMATO2', 'OQMCRM2', 'OQMNO2', 'OQMPSD2', 'OQMRMS2',

        'OPGATO1', 'OPGCRM1', 'OPGNO1', 'OPGPSD1', 'OPGRMS1', # load vector
        'OPGATO2', 'OPGCRM2', 'OPGNO2', 'OPGPSD2', 'OPGRMS2',

        'OPRATO1', 'OPRCRM1', 'OPRNO1', 'OPRPSD1', # pressure

        'RADCONS', 'RADEATC', 'RADEFFM',

        # ---------------
        # force/heat flux
        'HOEF1',
        'OEF1', 'OEF1X', 'OEF2',
        'OEFATO1', 'OEFCRM1', 'OEFNO1', 'OEFPSD1', 'OEFRMS1',
        'OEFATO2', 'OEFCRM2', 'OEFNO2', 'OEFPSD2', 'OEFRMS2',

        # ---------------
        # stress
        'OESNLXD', 'OESNLXR', 'OESNL1X',
        'OES1', 'OES1X', 'OES1X1', 'OES1C', 'OESVM1', 'OESVM1C',
        'OES1A',
        'OES2', 'OESVM2',
        'OESCP',
        'OESNL1',

        'OCRPG', 'OCRUG',
        'OESATO1', 'OESCRM1', 'OESNO1', 'OESPSD1', 'OESRMS1', 'OESXNO1', 'OESXRMS1',
        'OESATO2', 'OESCRM2', 'OESNO2', 'OESPSD2', 'OESRMS2',


        # ---------------
        #strain
        'OSTR1', 'OSTR1X', 'OSTR1C', 'OSTRVM1', 'OSTRVM1C',
        'OSTR2',
        'OESTRCP',

        'OSTRATO1', 'OSTRCRM1', 'OSTRNO1', 'OSTRPSD1', 'OSTRRMS1',
        'OSTRATO2', 'OSTRCRM2', 'OSTRNO2', 'OSTRPSD2', 'OSTRRMS2',
        'OSTRVM2',
        # ---------------
        # shock response spectra
        'OUPV1', # displacment
        'DOEF1', # force
        'DOES1', # stress

        # ---------------
        'OGPFB1',
        'ONRGY', 'ONRGY1',
        'OGS1', 'OGSTR1',
        'OEFIT', 'OESRT',
        'OSPDS1', 'OSPDSI1',
        'OSLIDE1', 'OSLIDEG1',
        'OBOLT1',
        'RAFCONS', 'RAFEATC',
        'RAGCONS', 'RAGEATC',
        'RANCONS', 'RANEATC',

    ]
    skip_results = {'gpdt', 'bgpdt', 'eqexin', 'psds', 'monitor1', 'monitor3', 'cstm'}
    for table_type in obj.get_table_types():
        if table_type in skip_results or table_type in skips or table_type.startswith('responses.'):
            continue

        res_dict = obj.get_result(table_type)
        if not isinstance(res_dict, dict):
            raise TypeError(table_type)

        for key, res in res_dict.items():
            #print(key)
            #print(res)
            _fix_subcase_id(key, res)

            if hasattr(res, 'table_name_str') and res.table_name_str not in {'LAMA'}: # params
                #print(table_type)
                res_categories2[res.table_name_str].append(res)
            else:
                # grid_point
                #print(table_type, unused_key)
                #print(res)
                class_name = res.__class__.__name__
                log.debug(class_name)
                if class_name not in ['GridPointWeight', 'PARAM',
                                     'RealEigenvalues', 'ComplexEigenvalues']:
                    log.warning(class_name)

    for table_name, results in sorted(res_categories2.items()):
        assert table_name in table_order, table_name

    total_case_count = 0
    pretables = ['LAMA', 'BLAMA', ] # 'CLAMA'
    if 'eigenvalues' not in skips:
        for unused_title, eigenvalue in obj.eigenvalues.items():
            res_categories2[eigenvalue.table_name].append(eigenvalue)
    if 'eigenvalues_fluid' not in skips:
        for unused_title, eigenvalue in obj.eigenvalues_fluid.items():
            res_categories2[eigenvalue.table_name].append(eigenvalue)
    # res_categories2 is a dictionary:
    #  - key: table_name
    #  - value:  displacment/stress/etc.
    #assert len(res_categories2['OUGV1']) == 13, (len(res_categories2['OUGV1']), res_categories2['OUGV1'])

    skip_tables = [
        'LAMA', 'BLAMA',
        #'OGPFB1',
        #'OESVM1', 'OESVM1C',
        #'OSTRVM1',
        #'OEF1X',
        #'OPG1',
        #'OQG1',
        #'OUGV1',
    ]

    footer = [4, 0, 4]
    footer_bytes = struct_3i.pack(*footer)

    struct_9i = Struct(endian + b'9i')
    for table_name in pretables + table_order:
        if table_name not in res_categories2:
            # no LAMA table in a static run
            continue
        if table_name in skip_tables:
            log.warning(f'skipping table={table_name}')
            continue
        results = res_categories2[table_name]
        itable = -1
        case_count = 0
        for result in results:
            element_name = ''
            new_result = True
            if hasattr(result, 'element_name'):
                element_name = ' - ' + result.element_name
                #print(element_name)

            #print(result.class_name)
            if hasattr(result, 'write_op2'):
                if result.table_name not in table_names_found:
                    table_names_found.append(result.table_name)

                fop2_ascii.write('-' * 60 + '\n')
                #if hasattr(result, 'is_bilinear') and result.is_bilinear():
                    #obj.log.warning("  *op2 - %s (%s) not written" % (
                        #result.__class__.__name__, result.element_name))
                    #continue
                isubcase = ''
                if hasattr(result, 'isubcase'): # no for eigenvalues
                    isubcase = result.isubcase
                    #print(f' {result.__class__.__name__} - isubcase={isubcase}')

                try:
                    #print(' %-6s - %s - isubcase=%s%s; itable=%s %s' % (
                        #table_name, result.__class__.__name__,
                        #isubcase, element_name, itable, new_result))
                    itable = result.write_op2(op2_file, fop2_ascii, itable, new_result,
                                              date, is_mag_phase=False, endian=endian)
                except Exception:
                    print(f' {result.__class__.__name__} - isubcase={isubcase}{element_name}')
                    raise
            elif hasattr(result, 'element_name'):
                if result.element_name in ['CBAR', 'CBEND']:
                    log.warning('skipping:\n%s' % result)
                    continue
                #log.warning('????:\n%s' % result)
            else:
                raise NotImplementedError(f'  *op2 - {result.__class__.__name__} not written')
                log.warning(f'  *op2 - {result.__class__.__name__} not written')
                #continue

            case_count += 1
            header = [
                4, itable, 4,
                4, 1, 4,
                4, 0, 4,
            ]
            #print('writing itable=%s' % itable)
            assert itable is not None, '%s itable is None' % result.__class__.__name__
            op2_file.write(struct_9i.pack(*header))
            fop2_ascii.write('footer2 = %s\n' % header)
            new_result = False

        #assert case_count > 0, case_count
        if case_count:
            #print(result.table_name, itable)
            #print('res_category_name=%s case_count=%s'  % (res_category_name, case_count))
            # close off the result - [4, 0, 4]
            op2_file.write(footer_bytes)
            fop2_ascii.write('close_a = %s\n' % footer)
            fop2_ascii.write('---------------\n')
            total_case_count += case_count
        total_case_count += case_count

    #if total_case_count == 0:
        #raise FatalError('total_case_count = 0')

    # close off the op2 - [4, 0, 4]
    op2_file.write(footer_bytes)
    fop2_ascii.write('close_b = %s\n' % footer)
    op2_file.close()
    fop2_ascii.close()
    return total_case_count, table_names_found

def _fix_subcase_id(key: Union[int, tuple[Any]], res: Any) -> None:
    """
    The subcase id may be inconsistent between the op2 result object
    and the key. The key takes priority to make doing a load case
    combination easier.
    """
    subcase_id = key
    if isinstance(key, tuple):
        subcase_id = key[0]
    res.isubcase = subcase_id

def write_op2_header(model: OP2, op2_file, fop2_ascii,
                     struct_3i: Struct,
                     post: int=-1, endian: bytes=b'<'):
    """writes the op2 header"""
    is_nx = model.is_nx
    is_msc = model.is_msc
    #is_nasa95 = model.is_nasa95
    is_optistruct = model.is_optistruct
    if model.date == (1, 1, 2000):  # (7, 24, 2020)
        today = datetime.datetime.today()
        model.date = (today.month, today.day, today.year)

    if post == -1:
    #_write_markers(op2_file, op2_ascii, [3, 0, 7])
        op2_file.write(struct_3i.pack(*[4, 3, 4,]))
        tape_code = b'NASTRAN FORT TAPE ID CODE - '
        if is_nx:
            op2_file.write(pack(endian + b'7i 28s i', *[4, 1, 4,
                                                        4, 7, 4,
                                                        28, tape_code, 28]))
            nastran_version = b'NX8.5   '
        elif is_msc or is_optistruct:
            day, month, year = model.date
            op2_file.write(pack(endian + b'9i 28s i', *[12,
                                                        day, month, year - 2000,
                                                        12, 4, 7, 4,
                                                        28, tape_code, 28]))
            nastran_version = b'XXXXXXXX'
        else:
            raise NotImplementedError(model._nastran_format)

        op2_file.write(pack(endian + b'4i 8s i', *[4, 2, 4,
                                                   #4, 2, 4,
                                                   #4, 1, 4,
                                                   #4, 8, 4,
                                                   8, nastran_version, 8]))
        op2_file.write(pack(endian + b'6i', *[4, -1, 4,
                                              4, 0, 4,]))
    elif post == -2:
        _write_markers(op2_file, fop2_ascii, [2, 4])
    else:
        raise RuntimeError(f'post = {post:d}; use -1 or -2')
