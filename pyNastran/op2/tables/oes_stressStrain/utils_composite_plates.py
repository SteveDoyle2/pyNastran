from __future__ import annotations
from struct import Struct
from typing import Any, TYPE_CHECKING
import numpy as np

from pyNastran.op2.op2_interface.op2_reader import mapfmt
# from pyNastran.op2.op2_helper import polar_to_real_imag

from pyNastran.op2.tables.utils import get_is_slot_saved, get_eid_dt_from_eid_device
from pyNastran.op2.tables.oes_stressStrain.complex.oes_composite_plates import (
    ComplexLayeredCompositeStrainArray, ComplexLayeredCompositeStressArray,
    ComplexLayeredCompositesArray)

from pyNastran.op2.tables.oes_stressStrain.real.oes_composite_plates import RealCompositePlateStressArray, RealCompositePlateStrainArray
from pyNastran.op2.tables.oes_stressStrain.real.oes_composite_plates_strength_ratio import RealCompositePlateStressStrengthRatioArray # , RealCompositePlateStrainStrengthRatioArray
#RealCompositePlateStressStrengthRatioArray, RealCompositePlateStrainStrengthRatioArray = None, None
#RealCompositePlateStrainStrengthRatioArray = None
from pyNastran.op2.tables.oes_stressStrain.random.oes_composite_plates import (
    RandomCompositePlateStressArray, RandomCompositePlateStrainArray,
    RandomCompositePlateStressVMArray, RandomCompositePlateStrainVMArray)

if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.op2.op2 import OP2


def oes_shells_composite(op2: OP2, data, ndata: int, dt, is_magnitude_phase: bool,
                         result_type: int, prefix: str, postfix: str) -> tuple[int, Any, Any]:
    """
    reads stress/strain for element type:
     - 95 : CQUAD4
     - 96 : CQUAD8
     - 97 : CTRIA3
     - 98 : CTRIA6 (composite)
     - 232 : QUADRLC (CQUADR-composite)
     - 233 : TRIARLC (CTRIAR-composite)

    """
    factor = op2.factor
    table_name_bytes = op2.table_name
    assert isinstance(table_name_bytes, bytes), table_name_bytes
    n = 0

    composite_element_name_map = {
        95: 'cquad4',
        96: 'cquad8',
        97: 'ctria3',
        98: 'ctria6',
        232: 'cquadr',
        233: 'ctriar',
    }
    try:
        element_name = composite_element_name_map[op2.element_type]
    except KeyError:  # pragma: no cover
        raise KeyError(op2.code_information())

    if op2.is_optistruct:
        is_stress = ('stress' in prefix)
    else:
        is_stress = op2.is_stress

    if is_stress:
        stress_strain = 'stress'
        obj_vector_real = RealCompositePlateStressArray
        obj_vector_strength = RealCompositePlateStressStrengthRatioArray
        # obj_vector_complex = ComplexCompositePlateStressArray
        layered_cls = ComplexLayeredCompositeStressArray
    else:
        stress_strain = 'strain'
        obj_vector_real = RealCompositePlateStrainArray
        obj_vector_strength = None  # RealCompositePlateStrainStrengthRatioArray
        # obj_vector_complex = ComplexCompositePlateStrainArray
        # obj_vector_random = RandomCompositePlateStrainArray
        layered_cls = ComplexLayeredCompositeStrainArray

    if op2._results.is_not_saved(prefix.rstrip('.')):
        return ndata, None, None

    result_name = prefix + f'{element_name}_composite_{stress_strain}' + postfix
    is_saved, slot = get_is_slot_saved(op2, result_name)
    if not is_saved:
        return ndata, None, None

    etype = op2.element_name
    sort_method = op2.sort_method

    num_wide = op2.num_wide
    # print('table_name =', table_name_bytes)
    if result_type == 0 and num_wide == 11:  # real
        #                    S T R E S S E S   I N   L A Y E R E D   C O M P O S I T E   E L E M E N T S   ( T R I A R )
        #   ELEMENT      PLY STRESSES IN FIBER AND MATRIX DIRECTIONS   INTER-LAMINAR  STRESSES  PRINCIPAL STRESSES (ZERO SHEAR)      MAX
        #     ID          ID   NORMAL-1     NORMAL-2     SHEAR-12    SHEAR XZ-MAT  SHEAR YZ-MAT  ANGLE    MAJOR        MINOR        SHEAR
        #      7070        1  7.50000E-01  3.00000E+00  9.86167E-08  -6.58903E-08  3.00000E+00   90.00  3.00000E+00  7.50000E-01  1.12500E+00
        #      7070        2 -7.50000E-01 -3.00000E+00 -9.86167E-08   0.0          0.0           -0.00 -7.50000E-01 -3.00000E+00  1.12500E+00
        ntotal = 44 * factor
        nelements = ndata // ntotal
        auto_return, is_vectorized = op2._create_oes_object4(
            nelements, result_name, slot, obj_vector_real)
        if auto_return:
            return nelements * ntotal, None, None

        obj = op2.obj
        if op2.is_debug_file:
            op2.binary_debug.write('  [cap, element1, element2, ..., cap]\n')
            op2.binary_debug.write(f'  cap = {ndata:d}  # assume 1 cap when there could have been multiple\n')
            op2.binary_debug.write(
                '  element1 = [eid_device, layer, o1, o2, t12, t1z, t2z, angle, major, minor, ovm)]\n')
            op2.binary_debug.write(f'  nelements={nelements:d}; nnodes=1 # centroid\n')

        if op2.use_vector and is_vectorized and sort_method == 1:
            n = nelements * op2.num_wide * 4

            istart = obj.itotal
            iend = istart + nelements
            obj._times[obj.itime] = dt

            if obj.itime == 0:
                ints = np.frombuffer(data, dtype=op2.idtype8).reshape(nelements, 11).copy()
                eids = ints[:, 0] // 10
                nids = ints[:, 1]
                obj.element_layer[istart:iend, 0] = eids
                obj.element_layer[istart:iend, 1] = nids

            floats = np.frombuffer(data, dtype=op2.fdtype8).reshape(nelements, 11)
            # [o1, o2, t12, t1z, t2z, angle, major, minor, ovm]
            obj.data[obj.itime, istart:iend, :] = floats[:, 2:].copy()
        else:
            if is_vectorized and op2.use_vector:  # pragma: no cover
                op2.log.debug(f'vectorize COMP_SHELL real SORT{sort_method}')

            n = oes_comp_shell_real_11(op2, data, ndata, obj,
                                       ntotal, nelements, etype, dt)

    # elif result_type == 1 and op2.num_wide == 9:  # TODO: imag? - not done...
    # TODO: vectorize
    # raise NotImplementedError('imaginary composite stress?')
    # msg = op2.code_information()
    # nelements = ndata // ntotal
    # obj_vector_complex = None
    # auto_return, is_vectorized = op2._create_oes_object4(
    # nelements, result_name, slot, obj_vector_complex)
    # if auto_return:
    # assert ntotal == op2.num_wide * 4
    # return nelements * ntotal, None, None

    ## TODO: this is an OEF result???
    ##    furthermore the actual table is calle dout as
    ##    'i8si4f4s', not 'i8si3fi4s'
    # ntotal = 36
    # nelements = ndata // ntotal
    # s = self.struct_i
    # s2 = Struct(op2._endian + b'8si3fi4s')
    # s3 = Struct(op2._endian + b'8si4f4s')
    # for i in range(nelements):
    ##out = s.unpack(data[n:n + ntotal])
    # eid_device, = s.unpack(data[n:n+4])
    ##t, = s.unpack(data[n:n+4])

    # if eid_device > 0:
    # out = s2.unpack(data[n+4:n+ntotal])
    # else:
    # unused_out1 = s2.unpack(data[n+4:n+ntotal])
    # out = s3.unpack(data[n+4:n+ntotal])
    # (theory, lamid, fp, fm, fb, fmax, fflag) = out

    # if op2.is_debug_file:
    # op2.binary_debug.write('%s-%s - (%s) + %s\n' % (op2.element_name, op2.element_type, eid_device, str(out)))
    # obj.add_new_eid_sort1(dt, eid, theory, lamid, fp, fm, fb, fmax, fflag)
    # n += ntotal
    # raise NotImplementedError('this is a really weird case...')
    elif result_type == 1 and num_wide == 11 and table_name_bytes in [b'OESCP', b'OESTRCP']:  # complex
        # OESCP - STRAINS IN LAYERED COMPOSITE ELEMENTS (QUAD4)
        ntotal = 44 * factor
        nelements = ndata // ntotal

        op2.log.warning(f'skipping complex {op2.table_name_str}-PCOMP')
        complex_obj = ComplexLayeredCompositeStressArray if op2.is_stress else ComplexLayeredCompositeStrainArray
        return nelements * ntotal, None, None

        auto_return, is_vectorized = op2._create_oes_object4(
            nelements, result_name, slot, complex_obj)
        if auto_return:
            return nelements * ntotal, None, None

        obj = op2.obj
        n = oes_shell_composite_complex_11(op2, data, obj,
                                           ntotal, nelements, sort_method,
                                           dt, is_magnitude_phase)
        return n, nelements, ntotal

    elif table_name_bytes in {b'OESRT', b'OESRTN'}:
        n, nelements, ntotal = oes_shells_composite_oesrt(
            op2, result_name, slot,
            result_type, sort_method, obj_vector_strength,
            data, ndata, dt)

    elif result_type == 1 and num_wide == 13 and table_name_bytes in [b'OESVM1C', b'OSTRVM1C']:  # complex
        # op2.log.warning(f'skipping complex {op2.table_name_str}-PCOMP (numwide=13)')
        ntotal = 52 * factor
        nelements = ndata // ntotal
        # return nelements * ntotal, None, None
        op2.table_name = table_name_bytes
        auto_return, is_vectorized = op2._create_oes_object4(
            nelements, result_name, slot, layered_cls)
        if auto_return:
            return nelements * ntotal, None, None

        obj = op2.obj
        if is_vectorized and op2.use_vector:
            # op2.log.warning(f'vectorize COMP_SHELL complex SORT{sort_method} (numwide=13)')
            n = len(data)

            istart = obj.itotal
            iend = istart + nelements
            obj._times[obj.itime] = dt

            if obj.itime == 0:
                ints = np.frombuffer(data, dtype=op2.idtype8).reshape(nelements, 13).copy()
                eids = ints[:, 0] // 10
                nids = ints[:, 1]
                obj.element_layer[istart:iend, 0] = eids
                obj.element_layer[istart:iend, 1] = nids

            floats = np.frombuffer(data, dtype=op2.fdtype8).reshape(nelements, 13)
            # [o1a, o2a, t12a, o1za, o2za, o1b, o2b, t12b, o1zb, e2zb, ovm]
            obj.data[obj.itime, istart:iend, :] = floats[:, 2:].copy()
            # struct1 = Struct(op2._endian + op2._analysis_code_fmt + b'i9f ff')
            # add_sort_x = getattr(obj, 'add_sort' + str(op2.sort_method))
            # for unused_i in range(nelements):
            #     edata = data[n:n + ntotal]
            #     out = struct1.unpack(edata)
            #
            #     (eid_device, ply_id,
            #      o1a, o2a, t12a, o1za, o2za,
            #      o1b, o2b, t12b, o1zb, e2zb, ovm,) = out
        else:
            n = oes_shell_composite_complex_13(op2, data, obj,
                                               ntotal, nelements, sort_method,
                                               dt, is_magnitude_phase)
        # return nelements * ntotal, None, None

    elif num_wide == 7 and (result_type == 2 or table_name_bytes == b'OSTRMS1C'):  # random (no VM)
        # TCODE,7 =0 Real
        # 2 PLY I Lamina Number
        # 3 EX1 RS Normal-1
        # 4 EY1 RS Normal-2
        # 5 ET1 RS Shear-12
        # 6 EL1 RS Shear-1Z
        # 7 EL2 RS Shear-2Z
        # 8 A1 RS Shear angle
        # 9 EMJRP1 RS Major Principal
        # 10 EMNRP1 RS Minor Principal
        # 11 ETMAX1 RS von Mises or Maximum shear
        #
        # TCODE,7 =1 Real/imaginary
        # 2 PLY I Lamina Number
        # 3 EX1 RS Normal-1
        # 4 EY1 RS Normal-2
        # 5 ET1 RS Shear-12
        # 6 EL1 RS Shear-1Z
        # 7 EL2 RS Shear-2Z
        # 8 EX1I RS Normal-1
        # 9 EY1I RS Normal-2
        # 10 ET1I RS Shear-12
        # 11 EL1I RS Shear-1Z
        # 12 EL2I RS Shear-2Z
        #
        # TCODE,7 =2 Random Response
        # 2 PLY I Lamina Number
        # 3 EX1 RS Normal-1
        # 4 EY1 RS Normal-2
        # 5 ET1 RS Shear-12
        # 6 EL1 RS Shear-1Z
        # 7 EL2 RS Shear-2Z
        ntotal = 28 * factor
        nelements = ndata // ntotal

        obj_vector_random = RandomCompositePlateStressArray if op2.is_stress else RandomCompositePlateStrainArray

        auto_return, is_vectorized = op2._create_oes_object4(
            nelements, result_name, slot, obj_vector_random)
        if auto_return:
            assert ntotal == op2.num_wide * 4
            return nelements * ntotal, None, None

        obj = op2.obj
        n = oes_composite_shells_nx_random_7(op2, data, obj, nelements, ntotal)

    elif result_type == 2 and num_wide == 8:  # random VM
        # analysis_code = 5   Frequency
        # table_code    = 805 OESXRM1C-OESXRMS1 - element RMS stresses for random analysis that includes von Mises stress output.
        # format_code   = 1   Real
        # result_type   = 2   Random
        # sort_method   = 1
        # sort_code     = 4
        #     sort_bits   = (0, 0, 1)
        #     data_format = 0   Real
        #     sort_type   = 0   Sort1
        #     is_random   = 1   Random Responses
        # random_code   = 0
        # element_type  = 95  QUAD4-nonlinear
        # num_wide      = 8
        # freq          = 0.0
        # NX Nastran
        if op2.is_stress:
            obj_vector_random = RandomCompositePlateStressVMArray
        else:
            obj_vector_random = RandomCompositePlateStrainVMArray
        msg = op2.code_information()
        msg = (f'etype={op2.element_name} ({op2.element_type}) '
               f'{op2.table_name_str}-COMP-random-numwide={op2.num_wide} '
               f'numwide_real=11 numwide_imag=9 result_type={result_type}')
        return op2._not_implemented_or_skip(data, ndata, msg), None, None

    elif result_type == 1 and num_wide == 12 and op2.is_msc:
        # analysis_code = 5   Frequency
        # table_code    = 5   OES1C-OES - Element Stress
        # format_code   = 2   Real/Imaginary
        # result_type   = 1   Complex
        # sort_method   = 1
        # sort_code     = 1
        #     sort_bits   = (1, 0, 1)
        #     data_format = 1   Real/Imaginary
        #     sort_type   = 0   Sort1
        #     is_random   = 1   Random Responses
        # random_code   = 0
        # element_type  = 95  QUAD4LC-composite
        # s_code        = 0   Coordinate Element - Stress Max Shear (Octahedral)
        # thermal       = 0   isHeatTransfer = False
        # thermal_bits  = [0, 0, 0, 0, 0]
        # num_wide      = 12
        # freq          = 990.0
        # MSC Nastran
        msg = (f'etype={op2.element_name} ({op2.element_type}) '
               f'{op2.table_name_str}-COMP-complex-numwide={num_wide} '
               # f'numwide_real=11 numwide_imag=9 result_type={result_type}'
               )
        # op2.log.warning(f'skipping complex {op2.table_name_str}-PCOMP-12')
        return op2._not_implemented_or_skip(data, ndata, msg), None, None

    elif result_type == 1 and num_wide == 11:
        # analysis_code = 9   Complex eigenvalues
        # table_code    = 5   OESCP-OES - Element Stress
        # format_code   = 2   Real/Imaginary
        # result_type   = 1   Complex
        # sort_method   = 1
        # sort_code     = 0
        #     sort_bits   = (1, 0, 0)
        #     data_format = 1   Real/Imaginary
        #     sort_type   = 0   Sort1
        #     is_random   = 0   Sorted Responses
        # random_code   = 0
        # element_type  = 95  QUAD4-nonlinear
        # num_wide      = 11
        # mode          = 0
        # eigr          = 0.0
        # eigi          = 0.0
        # NX Nastran
        msg = op2.code_information()
        msg = (f'etype={op2.element_name} ({op2.element_type}) '
               f'{op2.table_name_str}-COMP-random-numwide={num_wide} '
               f'numwide_real=11 numwide_imag=9 result_type={result_type}')
        return op2._not_implemented_or_skip(data, ndata, msg), None, None
    else:
        raise RuntimeError(op2.code_information())
        # msg = op2.code_information()
        # msg = (f'etype={op2.element_name} ({op2.element_type}) '
        # f'{op2.table_name_str}-COMP-random-numwide={op2.num_wide} '
        # f'numwide_real=11 numwide_imag=9 result_type={result_type}')
        # return op2._not_implemented_or_skip(data, ndata, msg), None, None
    return n, nelements, ntotal


def oes_shells_composite_oesrt(op2: OP2, result_name: str, slot: dict[Any, Any],
                               result_type: int, sort_method: int,
                               obj_vector_strength,
                               data: bytes, ndata: int, dt: Any):
    """Table of composite element ply strength ratio. Output by SDRCOMP"""
    n = 0
    factor = op2.factor
    if result_type == 0 and op2.num_wide == 9:  # real
        # strength_ratio.cquad4_composite_stress
        ntotal = 36 * factor
        nelements = ndata // ntotal

        auto_return, is_vectorized = op2._create_oes_object4(
            nelements, result_name, slot, obj_vector_strength)
        if auto_return:
            return nelements * ntotal, None, None

        obj = op2.obj
        if op2.is_debug_file:
            op2.binary_debug.write('  [cap, element1, element2, ..., cap]\n')
            op2.binary_debug.write('  cap = %i  # assume 1 cap when there could have been multiple\n' % ndata)
            op2.binary_debug.write(
                '  element1 = [eid_device, failure_theory, ply_id, strength_ratio_ply, failure_index_bonding, strength_ratio_bonding, flag, flag2)]\n')
            op2.binary_debug.write('  nelements=%i; nnodes=1 # centroid\n' % nelements)

        element_type = op2.element_type
        if op2.use_vector and is_vectorized and sort_method == 1 and 0:
            n = nelements * op2.num_wide * 4 * factor
        else:
            op2.log.warning(f'need to vectorize oes_shell_composite; {op2.element_name}-{op2.element_type} '
                            f'(numwide={op2.num_wide}) {op2.table_name_str}')
            n = oesrt_comp_shell_real_9(op2, data, ndata, obj,
                                        ntotal, nelements, element_type, dt)
    elif result_type == 1 and op2.num_wide == 9:  # complex
        # '          S T R E N G T H   R A T I O S   F O R   L A Y E R E D   C O M P O S I T E   E L E M E N T S   ( Q U A D 4 )'
        # '   ELEMENT  FAILURE        PLY  SRP-STRENGTH RATIO FOR PLY  SRB-STRENGTH RATIO FOR BONDING  STRENGTH RATIO FOR ELEMENT     FLAG'
        # '     ID      THEORY         ID  (DIRECT STRESSES/STRAINS)     (INTER-LAMINAR STRESSES)      MIN OF SRP,SRB FOR ALL PLIES'
        # '         3   HILL            1      4.510528E-02      '
        # '                                                                     4.160837E+01                                               '
        # '                             2      1.414114E-01      '
        # '                                                                     4.807957E+01                                               '
        # '                             3      2.013647E-01      '
        # '                                                                                                   4.510528E-02             *** '
        ntotal = 36 * factor
        nelements = ndata // ntotal
        if data is None:
            return ndata, None, None  # eid hill ply_id ? i f i blank

            auto_return, is_vectorized = op2._create_oes_object4(
                nelements, result_name, slot, obj_vector_strength)
            if auto_return:
                return nelements * ntotal, None, None

        is_vectorized = False
        op2.log.warning(f'OESRT: faking result; {op2.element_name}-{op2.element_type}')
        if op2.use_vector and is_vectorized and sort_method == 1 and 0:
            n = nelements * op2.num_wide * 4
            asdf
        else:
            op2.log.warning(f'need to vectorize oes_shell_composite; {op2.element_name}-{op2.element_type} '
                            f'(numwide={op2.num_wide}) {op2.table_name_str}')

            structi = Struct(op2._endian + b'i   8s   i      f i f i 4s')
            structf = Struct(op2._endian + b'i   8s   i      f i f f 4s')
            for unused_i in range(nelements):
                edata = data[n:n + ntotal]  # 4*9
                outs = structi.unpack(edata)
                eid, hill_bytes, ply_id, f1, i1, f2, minus_1, blank_bytes = outs
                hill = hill_bytes.decode('latin1').strip()
                blank = blank_bytes.decode('latin1').strip()
                if minus_1 != -1:
                    outs = structf.unpack(edata)
                    eid, hill_bytes, ply_id, f1, i1, f2, minus_1, blank_bytes = outs
                    # op2.log.error(f'minus1={minus_1}')
                assert hill in ['HILL', ''], hill
                assert blank in ['', '***'], blank
                # op2.show_data(edata)
                if eid != -1:
                    print(f'eid={eid} hill={hill!r} ply={ply_id} i1={i1} f2={f2:.4e} minus1={minus_1} blank={blank!r}')
                else:
                    print(f'    hill={hill!r} ply={ply_id} i1={i1} f2={f2:.4e} minus1={minus_1:.4e} blank={blank!r}')
                n += ntotal
    else:  # pragma: no cover
        raise RuntimeError(op2.code_information())
    return n, nelements, ntotal


def oes_comp_shell_real_11(op2: OP2, data: bytes, ndata: int,
                           obj: RealCompositePlateStressArray | RealCompositePlateStrainArray,
                           ntotal: int, nelements: int, etype: str, dt: Any) -> int:
    n = 0
    eid_old = 0
    fmt = mapfmt(op2._analysis_code_fmt + b'i9f', op2.size)  # 11
    struct1 = Struct(op2._endian + fmt)
    sort_method = op2.sort_method
    add_eid_sort_x = getattr(obj, 'add_eid_sort' + str(op2.sort_method))
    add_sort_x = getattr(obj, 'add_sort' + str(op2.sort_method))

    for unused_i in range(nelements):
        edata = data[n:n+ntotal]  # 4*11
        out = struct1.unpack(edata)
        (eid_device, layer, o1, o2, t12, t1z, t2z, angle, major, minor, ovm) = out
        eid, dt = get_eid_dt_from_eid_device(
            eid_device, op2.nonlinear_factor, sort_method)

        if op2.is_debug_file:
            op2.binary_debug.write('  eid=%i; layer=%i; C=[%s]\n' % (eid, layer, ', '.join(['%r' % di for di in out])))

        if eid != eid_old:
            # originally initialized to None, the buffer doesnt reset it, so it is the old value
            add_eid_sort_x(etype, dt, eid, layer, o1, o2,
                           t12, t1z, t2z, angle, major, minor, ovm)
        else:
            add_sort_x(dt, eid, layer, o1, o2,
                       t12, t1z, t2z, angle, major, minor, ovm)
        eid_old = eid
        n += ntotal
    return n


def oes_shell_composite_complex_11(op2: OP2,
                                   data: bytes,
                                   obj: ComplexLayeredCompositeStressArray | ComplexLayeredCompositeStrainArray,
                                   ntotal: int, nelements: int, sort_method: int,
                                   dt: Any, is_magnitude_phase: bool) -> int:
    """OESCP, OESTRCP"""
    n = 0
    fmt = mapfmt(op2._analysis_code_fmt + b'i9f', op2.size)
    struct1 = Struct(op2._endian + fmt)
    add_sort_x = getattr(obj, 'add_sort' + str(op2.sort_method))
    for unused_i in range(nelements):
        edata = data[n:n+ntotal]
        out = struct1.unpack(edata)

        (eid_device, ply_id, oxx, oyy, txy, txz, tyz, angle, omax, omin, max_shear) = out
        eid, dt = get_eid_dt_from_eid_device(
            eid_device, op2.nonlinear_factor, sort_method)
        #print(eid, out)

        #if op2.is_debug_file:
            #op2.binary_debug.write('%s-%s - (%s) + %s\n' % (op2.element_name, op2.element_type, eid_device, str(out)))
        #add_sort_x(dt, eid, theory, lamid, fp, fm, fb, fmax, fflag)
        add_sort_x(dt, eid, ply_id, oxx, oyy, txy, txz, tyz, angle, omax, omin, max_shear)
        n += ntotal
    return n


def oes_shell_composite_complex_13(op2: OP2,
                                   data: bytes,
                                   obj: ComplexLayeredCompositesArray,
                                   ntotal: int, nelements: int, sort_method: int,
                                   dt: Any, is_magnitude_phase: bool) -> int:
    """OESVM1C, OSTRVM1C"""
    # OESCP - STRAINS IN LAYERED COMPOSITE ELEMENTS (QUAD4)
    n = 0
    fmt = mapfmt(op2._analysis_code_fmt + b'i9f ff', op2.size)
    struct1 = Struct(op2._endian + fmt)
    add_sort_x = getattr(obj, 'add_sort' + str(op2.sort_method))
    for unused_i in range(nelements):
        edata = data[n:n+ntotal]
        out = struct1.unpack(edata)

        (eid_device, ply_id,
         o1a, o2a, t12a, o1za, o2za,
         o1b, o2b, t12b, o1zb, e2zb, ovm,) = out
        eid, dt = get_eid_dt_from_eid_device(
            eid_device, op2.nonlinear_factor, sort_method)
        #print(eid, out)

        #print('%s-%s - (%s) + %s\n' % (op2.element_name, op2.element_type, eid_device, str(out)))
        #if op2.is_debug_file:
            #op2.binary_debug.write('%s-%s - (%s) + %s\n' % (op2.element_name, op2.element_type, eid_device, str(out)))
        add_sort_x(dt, eid, ply_id,
                   o1a, o2a, t12a, o1za, o2za,
                   o1b, o2b, t12b, o1zb, e2zb, ovm)
        n += ntotal
    return n


def oes_composite_shells_nx_random_7(op2: OP2, data: bytes,
                                     obj: RandomCompositePlateStressArray | RandomCompositePlateStrainArray,
                                     nelements: int, ntotal: int) -> int:
    n = 0
    size = op2.size
    assert size == 4, size
    fmt = mapfmt(op2._analysis_code_fmt + b'i5f', op2.size)
    struct1 = Struct(op2._endian + fmt)
    add_sort_x_7words = getattr(obj, f'add_sort{op2.sort_method}_7words')
    #obj.add_sort1_7words
    #print(f'random nelements={nelements}')
    for unused_i in range(nelements):
        edata = data[n:n+ntotal]
        out = struct1.unpack(edata)

        (eid_device, ply_id, oxx, oyy, txy, txz, tyz) = out
        eid, dt = get_eid_dt_from_eid_device(
            eid_device, op2.nonlinear_factor, op2.sort_method)
        #print(eid, out)

        #if op2.is_debug_file:
            #op2.binary_debug.write('%s-%s - (%s) + %s\n' % (op2.element_name, op2.element_type, eid_device, str(out)))
        #print(obj)
        add_sort_x_7words(dt, eid, ply_id, oxx, oyy, txy, txz, tyz)
        n += ntotal
    #print(f'all_times n={len(obj.all_times)}')
    #print(f'all_elements n={len(obj.all_elements)}')
    return n


def oesrt_comp_shell_real_9(op2: OP2, data: bytes, ndata: int,
                            obj: RealCompositePlateStressStrengthRatioArray,
                            ntotal: int, nelements: int, etype: str, dt: Any) -> int:
    """
    # (1, b'TSAI-WU ', 1, 9.118781089782715, 1.3563156426940112e-19, 310.244140625, nan, 1.3563156426940112e-19)
    # (-1, b'        ', 2, 11.119629859924316, 1.3563156426940112e-19, 202.90777587890625, nan, 1.3563156426940112e-19)
    # (-1, b'        ', 3, 9.412755012512207, 1.3563156426940112e-19, 124.20950317382812, nan, 1.3563156426940112e-19)
    # (-1, b'        ', 4, 9.56696605682373, 1.3563156426940112e-19, 89.90276336669922, nan, 1.3563156426940112e-19)
    # (-1, b'        ', 5, 11.674513816833496, 1.3563156426940112e-19, 78.44961547851562, nan, 1.3563156426940112e-19)
    # (-1, b'        ', 6, 9.891059875488281, 1.3563156426940112e-19, 63.61910629272461, nan, 1.3563156426940112e-19)
    # (-1, b'        ', 7, 8.669649959069946e+22, 1.3563156426940112e-19, 63.619102478027344, nan, 1.3563156426940112e-19)
    # (-1, b'        ', 8, 1.0372562759530309e+23, 1.3563156426940112e-19, 63.61910629272461, nan, 1.3563156426940112e-19)
    # (-1, b'        ', 9, 26.741111755371094, 1.3563156426940112e-19, 78.44961547851562, nan, 1.3563156426940112e-19)
    # (-1, b'        ', 10, 48.35708236694336, 1.3563156426940112e-19, 89.90276336669922, nan, 1.3563156426940112e-19)
    # (-1, b'        ', 11, 25.865455627441406, 1.3563156426940112e-19, 124.20950317382812, nan, 1.3563156426940112e-19)
    # (-1, b'        ', 12, 25.44878578186035, 1.3563156426940112e-19, 202.90777587890625, nan, 1.3563156426940112e-19)
    # (-1, b'        ', 13, 46.05017852783203, 1.3563156426940112e-19, 310.244140625, nan, 1.3563156426940112e-19)
    # (-1, b'        ', 14, 24.654462814331055, 1.3563156426940112e-19, nan, 9.118781089782715, 1.3563156426940112e-19)
    # (2, b'TSAI-WU ', 1, 12.055706024169922, 1.3563156426940112e-19, 310.244140625, nan, 1.3563156426940112e-19)
    """
    n = 0
    eid_old = 0
    if op2.size == 4:
                                                                    # ft ply  sr,fi,sr
        structs = Struct(op2._endian + op2._analysis_code_fmt + b' 8s  i    f  4s  2f 4s')
        structf = Struct(op2._endian + op2._analysis_code_fmt + b' 8s  i    f  f   2f 4s')
    else:
        structs = Struct(op2._endian + op2._analysis_code_fmt + b' 16s q    d  8s  2d 4s')
        structf = Struct(op2._endian + op2._analysis_code_fmt + b' 16s q    d  d   2d 4s')
    sort_method = op2.sort_method
    assert sort_method == 1, 'oesrt_comp_shell_real_9'
    add_eid_sort_x = getattr(obj, 'add_eid_sort' + str(op2.sort_method))
    add_sort_x = getattr(obj, 'add_sort' + str(op2.sort_method))

    for unused_i in range(nelements):
        edata = data[n:n+ntotal]  # 4*9
        outs = structs.unpack(edata)

        # strings = (b'\x01\x00\x00\x00HFAIL   \x01\x00\x00\x00\x00\x00\xc8B\x00\x00\x00\x00\xff\xff\xff\xff\xff\xff\xff\xff    ',)
        #sr = 100.0; fi = 0.; sr = nan;
        #           eid ft         ply sr    fi     sr min  flag
        #ints    = (1, 'HFAIL   ', 1, 100.0, 0.0,   -1, -1, '    ')
        #floats  = (1, 'HFAIL   ', 1, 100.0, 0.0, nan, nan, '    ')
        (eid, failure_theory_bytes, ply_id, strength_ratio_ply, failure_index_bonding_bytes, strength_ratio_bonding,
         min_sr_bonding_fi_bonding, flag_bytes) = outs

        out = outs
        if failure_index_bonding_bytes != b'    ':
            outf = structf.unpack(edata)
            (eid, failure_theory_bytes, ply_id, strength_ratio_ply, failure_index_bonding, strength_ratio_bonding,
             min_sr_bonding_fi_bonding, flag_bytes) = outf
            out = outf
        else:
            failure_index_bonding = np.nan
        failure_theory = failure_theory_bytes.rstrip().decode('latin1')
        flag = flag_bytes.rstrip().decode('latin1')
        assert len(flag) <= 3, 'flag=%r' % flag

        #print(f'eid={eid}; ft={failure_theory!r}; sr_ply={strength_ratio_ply:g}; fi_bonding={failure_index_bonding!r} '
              #f'sr_bonding={strength_ratio_bonding}; min(sr_bonding, fi_bonding)={min_sr_bonding_fi_bonding:g} flag={flag!r}')
        #self.show_data(edata, types='ifs')

        #if op2.is_debug_file:
            #op2.binary_debug.write('  eid=%i; layer=%i; C=[%s]\n' % (eid, layer, ', '.join(['%r' % di for di in out])))
        if op2.is_debug_file:
            op2.binary_debug.write('%s-%s - (%s) + %s\n' % (op2.element_name, op2.element_type, eid, str(out)))
        #print(f'{op2.element_name}-{op2.element_type} - eid={eid}; ply_id={ply_id}; out={out}\n')

        #if ply_id == 1:  # not sufficient; ply_id may start at 101

        if eid == -1:
            add_sort_x(dt, eid_old, failure_theory, ply_id, strength_ratio_ply,
                       failure_index_bonding, strength_ratio_bonding, min_sr_bonding_fi_bonding, flag)
        #if eid != -1:
            # originally initialized to None, the buffer doesnt reset it, so it is the old value
            #add_eid_sort_x(etype, dt, eid, layer, o1, o2,
                           #t12, t1z, t2z, angle, major, minor, ovm)
        else:
            add_eid_sort_x(etype, dt, eid, failure_theory, ply_id, strength_ratio_ply,
                           failure_index_bonding, strength_ratio_bonding, min_sr_bonding_fi_bonding, flag)
            eid_old = eid
        n += ntotal
    return n




def oes_composite_solid_nx_real_center(op2: OP2, data: bytes,
                                       obj: RandomCompositePlateStressArray | RandomCompositePlateStrainArray,
                                       nelements: int, ntotal: int) -> int:
    n = 0
    size = op2.size

    fmt = mapfmt(op2._analysis_code_fmt + b'i 4s i 7f', size) # 11
    struct11 = Struct(op2._endian + fmt)
    #sort_method = op2.sort_method
    #add_eid_sort_x = getattr(obj, 'add_eid_sort' + str(op2.sort_method))
    #add_sort_x = getattr(obj, 'add_sort' + str(op2.sort_method))

    for unused_i in range(nelements):
        edata = data[n:n+ntotal]  # 4*11
        out = struct11.unpack(edata)
        #print(out)
        (eid_device, layer, location_bytes, grid, o11, o22, o33, t12, t23, t13, ovm) = out
        eid, dt = get_eid_dt_from_eid_device(
            eid_device, op2.nonlinear_factor, op2.sort_method)

        location = location_bytes.strip().decode('latin1')
        assert location == 'MID', out
        #print(f'eid,layer=({eid},{layer}) location={location!r} grid={grid} o11={o11:g} o22={o22:g} o33={o33:g} t12={t12:g} t1z={t13:g} t2z={t23:g} ovm={ovm:g}')
        obj.add_sort1(dt, eid, layer, location, grid, o11, o22, o33, t12, t23, t13, ovm)
        n += ntotal
    return n
