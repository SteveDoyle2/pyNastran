from __future__ import annotations
from struct import Struct
from typing import Any, TYPE_CHECKING
import numpy as np

from pyNastran.op2.op2_interface.utils import reshape_bytes_block_strip
# from pyNastran.op2.tables.utils import get_eid_dt_from_eid_device
from pyNastran.op2.tables.utils import get_is_slot_saved, get_eid_dt_from_eid_device

from pyNastran.op2.tables.oef_forces.oef_force_objects import (
    FailureIndicesArray,
)
# from pyNastran.op2.tables.oef_forces.oef_complex_force_objects import (
#     ComplexCShearForceArray,
# )

if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.op2.op2 import OP2


def oef_shells_composite(op2: OP2, data, ndata, dt, unused_is_magnitude_phase,
                         result_type, prefix, postfix):
    """
    95 - CQUAD4
    96 - CQUAD8
    97 - CTRIA3
    98 - CTRIA6 (composite)
    232 - CQUADR
    233 - CTRIAR

    """
    if op2.element_type == 95:
        result_name = prefix + 'cquad4_composite_force' + postfix
    elif op2.element_type == 96:
        result_name = prefix + 'cquad8_composite_force' + postfix
    elif op2.element_type == 97:
        result_name = prefix + 'ctria3_composite_force' + postfix
    elif op2.element_type == 98:
        result_name = prefix + 'ctria6_composite_force' + postfix
    elif op2.element_type == 232:
        result_name = prefix + 'cquadr_composite_force' + postfix
    elif op2.element_type == 233:
        result_name = prefix + 'ctriar_composite_force' + postfix
    else:  # pragma: no cover
        raise NotImplementedError(op2.code_information())

    is_saved, slot = get_is_slot_saved(op2, result_name)
    if not is_saved:
        return ndata, None, None

    n = 0
    factor = op2.factor
    if op2.format_code == 1 and op2.num_wide == 9:  # real
        ntotal = 36 * factor # 9 * 4
        nelements = ndata // ntotal


        auto_return, is_vectorized = op2._create_oes_object4(
            nelements, result_name, slot, FailureIndicesArray)
        #print('read_mode ', op2.read_mode, auto_return, is_vectorized)
        if auto_return:
            #op2._data_factor = nnodes_all
            return nelements * ntotal, None, None

        obj = op2.obj
        nelements = ndata // ntotal

        ## TODO: add
        #if op2.is_debug_file:
            #op2.binary_debug.write('  [cap, element1, element2, ..., cap]\n')
            #op2.binary_debug.write('  cap = %i  # assume 1 cap when there could have been multiple\n' % ndata)
            ##op2.binary_debug.write('  #centeri = [eid_device, j, grid, fd1, sx1, sy1, txy1, angle1, major1, minor1, vm1,\n')
            ##op2.binary_debug.write('  #                                fd2, sx2, sy2, txy2, angle2, major2, minor2, vm2,)]\n')
            ##op2.binary_debug.write('  #nodeji = [eid, ilayer, o1, o2, t12, t1z, t2z, angle, major, minor, ovm)]\n')
            #op2.binary_debug.write('  nelements=%i; nnodes=1 # centroid\n' % nelements)
        n = oef_shells_composite_real_9(op2, data, obj, nelements, ntotal, dt)
    elif op2.format_code in {2, 3} and op2.num_wide == 9:  # complex
        #  device_code   = 1   Print
        #  analysis_code = 5   Frequency
        #  table_code    = 25  OEFIT-OEF - Composite failure indices
        #  format_code   = 2   Real/Imaginary
        #  result_type   = 1   Complex
        #  sort_method   = 1
        #  sort_code     = 1
        #      sort_bits   = (1, 0, 1)
        #      data_format = 1   Real/Imaginary
        #      sort_type   = 0   Sort1
        #      is_random   = 1   Random Responses
        #  random_code   = 0
        #  element_type  = 95  QUAD4LC-composite
        #  s_code        = None ???
        #  thermal       = 0   isHeatTransfer = False
        #  thermal_bits  = [0, 0, 0, 0, 0]
        #  num_wide      = 9
        #  isubcase      = 1
        #  MSC Nastran
        msg = op2.code_information()
        msg = (f'etype={op2.element_name} ({op2.element_type}) '
               f'{op2.table_name_str}-COMP-random-numwide={op2.num_wide} '
               f'numwide_real=11 numwide_imag=9 result_type={result_type}')
        if data is None:
            return op2._not_implemented_or_skip(data, ndata, msg), None, None
        '      ELEMENT-ID =      11'
        '          F A I L U R E   I N D I C E S   F O R   L A Y E R E D   C O M P O S I T E   E L E M E N T S   ( Q U A D 4 )'
        '    PLY     FAILURE              FP=FAILURE INDEX FOR PLY    FB=FAILURE INDEX FOR BONDING   FAILURE INDEX FOR ELEMENT      FLAG'
        '     ID      THEORY      FREQ  (DIRECT STRESSES/STRAINS)     (INTER-LAMINAR STRESSES)      MAX OF FP,FB FOR ALL PLIES  '
        '         1   HFABRIC  2.0000E+01        0.1217    1 '
        '                                                                           0.0000                                               '
        '                      4.0000E+01        0.5472   -1 '
        '                                                                           0.0001                                               '
        '                      6.0000E+01        0.2454    1 '
        '                                                                           0.0000                                               '
        '                      8.0000E+01        0.3217    1 '
        '                                                                           0.0000                                               '
        '                      2.0000E+02        0.7284    1 '
        '                                                                           0.0000                        0.7284                 '
        assert op2.sort_method == 1, op2.code_information()
        ntotal = 36 * factor
        nelements = ndata // ntotal
        sf = Struct(op2._endian + b'i8s if i ff 4s')  # if if i
        si = Struct(op2._endian + b'i8s if i fi 4s')  # if if f
        sf2 = Struct(op2._endian + b'f')
        for unused_i in range(nelements):
            edata = data[n:n + ntotal]
            out = si.unpack(edata)
            (ply_id, failure_theory_bytes,
             c, d,
             e, f,
             g,
            #ply_fp, failure_index_for_ply,
            #ply_fb, failure_index_for_bonding,
            #failure_index_for_element,
            end) = out
            failure_theory = failure_theory_bytes.decode(op2._encoding).rstrip()
            #print(f'ply_id={ply_id} failure_theory={failure_theory!r} ply_fp={ply_fp} failure_index_for_ply={failure_index_for_ply:.3g} '
                  #f'ply_fb={ply_fb} fi_bonding={failure_index_for_bonding:.3g} fi_for_element={failure_index_for_element:g}')
            #if failure_index_for_bonding != -1:
                #failure_index_for_bonding = -1000.
                #failure_index_for_element = f.unpack(edata[-8:-4])

                #*junk, failure_index_for_bonding, failure_index_for_element, end = sf.unpack(edata)
                #print(f'ply_id={ply_id} failure_theory={failure_theory!r} ply_fp={ply_fp} fi_ply={failure_index_for_ply:.3g} '
                      #f'ply_fb={ply_fb} fi_bonding={failure_index_for_bonding:.3g} fi_for_element={failure_index_for_element:g}')
            #else:
            if g != -1:
                g, = sf2.unpack(edata[-8:-4])
                #op2.show_data(edata)
                out = (ply_id, failure_theory_bytes, c, d, e, f, g)
            #print(out)
            assert end in {b'    ', b'*** '}, end
            n += ntotal
        #aaa
        return op2._not_implemented_or_skip(data, ndata, msg), None, None
    else:  # pragma: no cover
        raise RuntimeError(op2.code_information())
        #msg = op2.code_information()
        #print(msg)
        #return op2._not_implemented_or_skip(data, ndata, msg), None, None
    return n, nelements, ntotal


def oef_shells_composite_real_9(op2: OP2, data: bytes,
                                obj: FailureIndicesArray,
                                nelements: int, ntotal: int,
                                dt: Any) -> int:
    """
    2 THEORY(2) CHA/R4
    4 PLY       I
    5 DIRECT    RS
    6 INDEX     CHA/R4
    7 LAMIN     RS (I)
    8 MAX       RS (I)
    9 FLAG      CHA
    """
    n = 0
    size = op2.size
    if size == 4:
        #                                5 6  7 8-i/f 9
        s1 = Struct(op2._endian + b'i8sif  i f i     4s')
        s2 = Struct(op2._endian + b'i8sif 4s f f     4s')
    elif size == 8:
        s1 = Struct(op2._endian + b'q16sqd  q d q     8s')
        s2 = Struct(op2._endian + b'q16sqd 8s d d     8s')
    else:  # pragma: no cover
        raise RuntimeError(size)

    eid_old = None
    #print(op2.element_type)
    add_sort_x = getattr(obj, 'add_sort' + str(op2.sort_method))
    for unused_i in range(nelements):
        #2 THEORY(2) CHAR4 Theory
        #4 LAMID     I Lamina number

        #5 FP       RS Failure index for direct stresses
        #6 FM       RS Failure mode for maximum strain theory
        #7 FB       RS Failure index for interlaminar shear stress or -1
        #8 FMAX     RS Maximum of FP and FB or -1.
        #9 FFLAG CHAR4 Failure flag
        edata = data[n:n+ntotal]  # 4*9
        #print(self.show_data(edata[4+8+4:-4-4], types='ifs'))
        out1 = s1.unpack(edata)
        out2 = s2.unpack(edata)

        # failure_stress_for_ply = failure_strain_for_ply = failure_index_for_ply???
        # i    8s               i      f
        (eid, failure_theoryb, ply_id, failure_stress_for_ply,
         # 4s   7-f/i                8-f/i      9-4s
         flagi, interlaminar_stress, max_value, failure_flagb,
         #failure_index_for_bonding,
         #failure_index_for_element,
         #flag,
         #direct_stress_or_strain,
         #interlaminar_stress,
         #max_of_fb_fp_for_all_plies
        ) = out1

        (_eid, _failure_theoryb, _ply_id, _failure_stress_for_ply,
         # 4s   7-f/i                8-f/i      9-4s
         flagb, _interlaminar_stress, max_value_float, _failure_flagb,
         #failure_index_for_bonding,
         #failure_index_for_element,
         #flag,
         #direct_stress_or_strain,
         #interlaminar_stress,
         #max_of_fb_fp_for_all_plies
        ) = out2
        #print(interlaminar_stress, _interlaminar_stress)

        #print('failure_flagb = %r' % failure_flagb)
        if flagi == 0:
            flag = ''
        else:
            #print('flagb = %r' % flagb)
            flag = reshape_bytes_block_strip(flagb, size=size)

        failure_theory = reshape_bytes_block_strip(failure_theoryb, size=size)
        failure_flag = reshape_bytes_block_strip(failure_flagb, size=size)
        #print('flag = %r' % flag)
        #print('failure_flag = %r' % failure_flag)

        if max_value == -1:
            max_value = np.nan
        else:
            max_value = max_value_float # out2[6]

        if eid == -1:
            #print(f'  ply_id={ply_id} failure_stress_for_ply={failure_stress_for_ply:g} '
                  #f'flag={flag!r} interlaminar_stress={interlaminar_stress} '
                  #f'max_value={max_value:g} failure_flag={failure_flag!r}')
            eid = eid_old
        else:
            #print(f"eid={eid} ft={failure_theory!r}\n"
                  #f'  ply_id={ply_id} failure_stress_for_ply={failure_stress_for_ply:g} '
                  #f'flag={flag!r} interlaminar_stress={interlaminar_stress} '
                  #f'max_value={max_value} failure_flag={failure_flag!r}')
            eid_old = eid
        assert flag in ['', '-1', '-2', '-12', 'IN'], f'flag={flag!r} flagb={flagb!r}'

        # 'HILL' for the Hill theory.
        # 'HOFF' for the Hoffman theory.
        # 'TSAI' for the Tsai-Wu theory.
        # 'STRN' for the Maximum Strain theory.
        # 'HFAIL' for the Hashin failure criterion
        # 'HTAPE' for the Hashin tape criterion
        # 'HFABR' for the Hashin fabric criterion
        assert failure_theory in ['TSAI-WU', 'STRAIN', 'HILL', 'HOFFMAN', 'HFAIL', 'HFABRIC', 'HTAPE', ''], f'failure_theory={failure_theory!r}'
        assert failure_flag in ['', '***'], 'failure_flag=%r' % failure_flag
        add_sort_x(dt, eid, failure_theory, ply_id, failure_stress_for_ply, flag,
                   interlaminar_stress, max_value, failure_flag)
        n += ntotal

    #add_sort_x = getattr(obj, 'add_sort' + str(op2.sort_method))
    #s = Struct(op2._endian + b'i8si4f4s')
    #for i in range(nelements):
        #if i % 10000 == 0:
            #print 'i = ', i
        #edata = data[n:n+ntotal]  # 4*9
        #out = s.unpack(edata)
        #(eid_device, theory, lamid, failure_index_direct_stress, failure_mode_max_shear,
                 #failure_index_interlaminar_shear, fmax, failure_flag) = out
        #eid, dt = get_eid_dt_from_eid_device(
            #eid_device, op2.nonlinear_factor, op2.sort_method)
        #if op2.is_debug_file:
            #if eid > 0:
                #op2.binary_debug.write('  eid=%i; C=[%s]\n' % (', '.join(['%r' % di for di in out]) ))
            #else:
                #op2.binary_debug.write('      %s  C=[%s]\n' % (' ' * len(str(eid)), ', '.join(['%r' % di for di in out]) ))

        #if eid > 0:
            #obj.add_new_eid_sort1(eType, dt, eid, o1, o2, t12, t1z, t2z, angle, major, minor, ovm)
        #else:
            #add_sort_x(dt, eid, o1, o2, t12, t1z, t2z, angle, major, minor, ovm)
        #n += ntotal
    return n
