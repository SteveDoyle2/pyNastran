#pylint: disable=W0201,C0301,C0111
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
#import sys
#import copy
from datetime import date
from collections import defaultdict
from struct import pack, Struct

import pyNastran
from pyNastran.op2.op2_interface.op2_f06_common import OP2_F06_Common
from pyNastran.op2.op2_interface.write_utils import _write_markers
from pyNastran.op2.errors import FatalError
from .writer.geom1 import write_geom1
from .writer.geom2 import write_geom2
from .writer.geom3 import write_geom3
from .writer.geom4 import write_geom4
from .writer.ept import write_ept
from .writer.mpt import write_mpt


def make_stamp(title, today=None):
    if 'Title' is None:
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
    str_today = str_today  #.strip()

    release_date = '02/08/12'  # pyNastran.__releaseDate__
    release_date = ''
    build = 'pyNastran v%s %s' % (pyNastran.__version__, release_date)
    if title is None:
        title = ''
    out = '1    %-67s   %-19s %-22s PAGE %%5i\n' % (title.strip(), str_today, build)
    return out


class OP2Writer(OP2_F06_Common):
    def __init__(self, op2):
        self.log = op2.log
        OP2_F06_Common.__init__(self)
        self.card_count = {}

    #def make_f06_header(self):
        #"""If this class is inherited, the F06 Header may be overwritten"""
        #return make_f06_header()

    def make_stamp(self, title, today):
        """If this class is inherited, the PAGE stamp may be overwritten"""
        return make_stamp(title, today)

    def write_op2(self, op2_outname, obj=None, is_mag_phase=False,
                  post=-1, endian=b'<'):
        """
        Writes an OP2 file based on the data we have stored in the object

        Parameters
        ----------
        op2_outname : str
            the name of the F06 file to write
        obj : OP2(); default=None -> self
            the OP2 object if you didn't inherit the class
        is_mag_phase : bool; default=False
            should complex data be written using Magnitude/Phase
            instead of Real/Imaginary (default=False; Real/Imag)
            Real objects don't use this parameter.
        """
        #print('writing %s' % op2_outname)
        struct_3i = Struct(endian + b'3i')

        if obj is None:
            obj = self
        if isinstance(op2_outname, str):
            fop2 = open(op2_outname, 'wb')
            fop2_ascii = open(op2_outname + '.txt', 'w')
            #print('op2 out = %r' % op2_outname)
        else:
            assert isinstance(op2_outname, file), 'type(op2_outname)= %s' % op2_outname
            fop2 = op2_outname
            op2_outname = op2_outname.name
            print('op2_outname =', op2_outname)

        #op2_ascii.write('writing [3, 7, 0] header\n')
        #if markers == [3,]:  # PARAM, POST, -1
            #self.op2_reader.read_markers([3])
            #data = self.read_block()

            #self.op2_reader.read_markers([7])
            #data = self.read_block()
            #data = self._read_record()
            #self.op2_reader.read_markers([-1, 0])
        #elif markers == [2,]:  # PARAM, POST, -2
        if post == -1:
        #_write_markers(op2, op2_ascii, [3, 0, 7])
            fop2.write(struct_3i.pack(*[4, 3, 4,]))
            tape_code = b'NASTRAN FORT TAPE ID CODE - '
            fop2.write(pack('7i 28s i', *[4, 1, 4,
                                          4, 7, 4,
                                          28, tape_code, 28]))

            nastran_version = b'NX8.5   ' if obj.is_nx else b'XXXXXXXX'
            fop2.write(pack(b'4i 8s i', *[4, 2, 4,
                                          #4, 2, 4,
                                          #4, 1, 4,
                                          #4, 8, 4,
                                          8, nastran_version, 8]))
            fop2.write(pack(b'6i', *[4, -1, 4,
                                     4, 0, 4,]))
        elif post == -2:
            _write_markers(fop2, fop2_ascii, [2, 4])
        else:
            raise RuntimeError('post = %r; use -1 or -2' % post)

        write_geom1(fop2, fop2_ascii, obj)
        write_geom2(fop2, fop2_ascii, obj)
        write_geom3(fop2, fop2_ascii, obj)
        write_geom4(fop2, fop2_ascii, obj)
        write_ept(fop2, fop2_ascii, obj)
        write_mpt(fop2, fop2_ascii, obj)
        #write_dit(fop2, fop2_ascii, obj)
        #write_dynamic(fop2, fop2_ascii, obj)
        if obj.grid_point_weight.reference_point is not None:
            if hasattr(obj.grid_point_weight, 'write_op2'):
                print("grid_point_weight")
                obj.grid_point_weight.write_op2(fop2, endian=endian)
            else:
                print("*op2 - grid_point_weight not written")


        #is_mag_phase = False

        # we writte all the other tables
        # nastran puts the tables in order of the Case Control deck,
        # but we're lazy so we just hardcode the order


        #if markers == [3,]:  # PARAM, POST, -2
            #self.op2_reader.read_markers([3])
            #data = self.read_block()
            #self.op2_reader.read_markers([7])
            #data = self.read_block()
            ##self.show(100)
            #data = self._read_record()

        self._write(obj, fop2, fop2_ascii, struct_3i, endian)

    def _write(self, obj, fop2, fop2_ascii, struct_3i, endian):
        res_categories2 = defaultdict(list)
        table_order = [
            'OUGV1',
            'BOUGV1', 'BOPHIG',
            'OUPV1',
            'TOUGV1',
            'OAGATO1', 'OAGCRM1', 'OAGNO1', 'OAGPSD1', 'OAGRMS1',
            'OAGATO2', 'OAGCRM2', 'OAGNO2', 'OAGPSD2', 'OAGRMS2',

            'OQG1',
            'OQGV1',
            'OQP1',
            'OQMG1',
            'OPGV1', 'OPG1', 'OPNL1',
            'OUGATO1', 'OUGCRM1', 'OUGNO1', 'OUGPSD1', 'OUGRMS1',
            'OUGATO2', 'OUGCRM2', 'OUGNO2', 'OUGPSD2', 'OUGRMS2',

            'OVGATO1', 'OVGCRM1', 'OVGNO1', 'OVGPSD1', 'OVGRMS1',
            'OVGATO2', 'OVGCRM2', 'OVGNO2', 'OVGPSD2', 'OVGRMS2',

            'OQGATO1', 'OQGCRM1', 'OQGNO1', 'OQGPSD1', 'OQGRMS1',
            'OQGATO2', 'OQGCRM2', 'OQGNO2', 'OQGPSD2', 'OQGRMS2',

            'OQMATO1', 'OQMCRM1', 'OQMNO1', 'OQMPSD1', 'OQMRMS1',
            'OQMATO2', 'OQMCRM2', 'OQMNO2', 'OQMPSD2', 'OQMRMS2',

            'OPGATO1', 'OPGCRM1', 'OPGNO1', 'OPGPSD1', 'OPGRMS1',
            'OPGATO2', 'OPGCRM2', 'OPGNO2', 'OPGPSD2', 'OPGRMS2',

            'DOEF1', 'HOEF1',
            'OEF1', 'OEF1X',
            'OEFATO1', 'OEFCRM1', 'OEFNO1', 'OEFPSD1', 'OEFRMS1',
            'OEFATO2', 'OEFCRM2', 'OEFNO2', 'OEFPSD2', 'OEFRMS2',

            'OESNLXD', 'OESNLXR', 'OESNL1X',
            'OES1', 'OES1X', 'OES1X1',
            'OES1C',
            'OESCP',
            'OCRPG', 'OCRUG',
            'OESATO1', 'OESCRM1', 'OESNO1', 'OESPSD1', 'OESRMS1',
            'OESATO2', 'OESCRM2', 'OESNO2', 'OESPSD2', 'OESRMS2',

            'OSTRATO1', 'OSTRCRM1', 'OSTRNO1', 'OSTRPSD1', 'OSTRRMS1',
            'OSTRATO2', 'OSTRCRM2', 'OSTRNO2', 'OSTRPSD2', 'OSTRRMS2',

            'OESTRCP',
            'OSTR1C',
            'OSTR1X',

            'OGPFB1',
            'ONRGY1',
            'OGS1',
        ]
        for table_type in obj.get_table_types():
            if table_type in ['gpdt', 'eqexin']:
                continue
            res_dict = obj.get_result(table_type)
            for unused_key, res in res_dict.items():
                if hasattr(res, 'table_name_str'): # params
                    res_categories2[res.table_name_str].append(res)

        for table_name, results in sorted(res_categories2.items()):
            assert table_name in table_order, table_name

        total_case_count = 0
        pretables = ['LAMA', 'BLAMA']
        for unused_title, eigenvalue in obj.eigenvalues.items():
            res_categories2[eigenvalue.table_name].append(eigenvalue)

        for table_name in pretables + table_order:
            if table_name not in res_categories2:
                continue
            results = res_categories2[table_name]
            itable = -1
            case_count = 0
            for result in results:
                element_name = ''
                new_result = True
                if hasattr(result, 'element_name'):
                    element_name = ' - ' + result.element_name

                if hasattr(result, 'write_op2'):
                    #if hasattr(result, 'is_bilinear') and result.is_bilinear():
                        #obj.log.warning("  *op2 - %s (%s) not written" % (
                            #result.__class__.__name__, result.element_name))
                        #continue
                    isubcase = ''
                    if hasattr(result, 'isubcase'): # no for eigenvalues
                        isubcase = result.isubcase
                    try:
                        #print(' %-6s - %s - isubcase=%s%s; itable=%s %s' % (
                            #table_name, result.__class__.__name__,
                            #isubcase, element_name, itable, new_result))
                        itable = result.write_op2(fop2, fop2_ascii, itable, new_result,
                                                  obj.date, is_mag_phase=False, endian=endian)
                    except:
                        print(' %s - isubcase=%s%s' % (result.__class__.__name__,
                                                       isubcase, element_name))
                        raise
                else:
                    raise NotImplementedError("  *op2 - %s not written" % result.__class__.__name__)
                    #obj.log.warning("  *op2 - %s not written" % result.__class__.__name__)
                    #continue

                case_count += 1
                header = [
                    4, itable, 4,
                    4, 1, 4,
                    4, 0, 4,
                ]
                #print('writing itable=%s' % itable)
                assert itable is not None, '%s itable is None' % result.__class__.__name__
                fop2.write(pack(b'9i', *header))
                fop2_ascii.write('footer2 = %s\n' % header)
                new_result = False

            assert case_count > 0, case_count
            if case_count:
                #print(result.table_name, itable)
                #print('res_category_name=%s case_count=%s'  % (res_category_name, case_count))
                # close off the result
                footer = [4, 0, 4]
                fop2.write(struct_3i.pack(*footer))
                fop2_ascii.write('close_a = %s\n' % footer)
                fop2_ascii.write('---------------\n')
                total_case_count += case_count
            total_case_count += case_count

        if total_case_count == 0:
            raise FatalError('total_case_count = 0')
        # close off the op2
        footer = [4, 0, 4]
        fop2.write(struct_3i.pack(*footer))
        fop2_ascii.write('close_b = %s\n' % footer)
        fop2.close()
        fop2_ascii.close()
