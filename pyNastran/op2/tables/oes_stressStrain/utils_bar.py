from __future__ import annotations
from struct import Struct
from typing import Any, TYPE_CHECKING
import numpy as np

from pyNastran.op2.op2_interface.op2_reader import mapfmt
from pyNastran.op2.op2_helper import polar_to_real_imag

from pyNastran.op2.tables.utils import get_eid_dt_from_eid_device

if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.op2.tables.oes_stressStrain.real.oes_bars import RealBarStressArray, RealBarStrainArray
    from pyNastran.op2.tables.oes_stressStrain.complex.oes_bars import ComplexBarStressArray, ComplexBarStrainArray
    from pyNastran.op2.tables.oes_stressStrain.random.oes_bars import RandomBarStressArray, RandomBarStrainArray
    from pyNastran.op2.op2 import OP2


def oes_cbar_real_16(op2: OP2, data: bytes,
                     obj: RealBarStressArray | RealBarStrainArray,
                     nelements: int, ntotal: int, dt: Any) -> int:
    n = 0
    fmt = mapfmt(op2._endian + op2._analysis_code_fmt + b'15f', op2.size)
    struct1 = Struct(fmt)
    add_sort_x = getattr(obj, 'add_new_eid_sort' + str(op2.sort_method))
    for i in range(nelements):
        edata = data[n:n+ntotal]
        out = struct1.unpack(edata)
        (eid_device,
         s1a, s2a, s3a, s4a, axial, smaxa, smina, margin_tension,
         s1b, s2b, s3b, s4b, smaxb, sminb, margin_compression) = out
        eid, dt = get_eid_dt_from_eid_device(
            eid_device, op2.nonlinear_factor, op2.sort_method)
        if op2.is_debug_file:
            op2.binary_debug.write('  eid=%i; C%i=[%s]\n' % (
                eid, i, ', '.join(['%r' % di for di in out])))
        n += ntotal
        add_sort_x(
            dt, eid,
            s1a, s2a, s3a, s4a, axial, smaxa, smina, margin_tension,
            s1b, s2b, s3b, s4b, smaxb, sminb, margin_compression)

    if op2.is_sort2:
        #print(add_sort_x)
        #print(''.join(obj.get_stats()))
        #print(f'{self.table_name} sort_method={op2.sort_method}', obj._times)
        assert len(np.unique(obj._times)) == len(obj._times), obj._times.tolist()
    return n


def oes_cbar_complex_19(op2: OP2,
                        data: bytes,
                        obj: ComplexBarStressArray | ComplexBarStrainArray,
                        nelements: int, ntotal: int,
                        is_magnitude_phase: bool) -> int:
    n = 0
    struct1 = Struct(mapfmt(op2._endian + op2._analysis_code_fmt + b'18f', op2.size))
    for i in range(nelements):
        edata = data[n:n+ntotal]
        n += ntotal
        out = struct1.unpack(edata)
        (eid_device,
         s1ar, s2ar, s3ar, s4ar, axialr,
         s1ai, s2ai, s3ai, s4ai, axiali,
         s1br, s2br, s3br, s4br,
         s1bi, s2bi, s3bi, s4bi) = out

        eid, dt = get_eid_dt_from_eid_device(
            eid_device, op2.nonlinear_factor, op2.sort_method)
        if op2.is_debug_file:
            op2.binary_debug.write('  eid=%i; C%i=[%s]\n' % (
                eid, i, ', '.join(['%r' % di for di in out])))
        if is_magnitude_phase:
            s1a = polar_to_real_imag(s1ar, s1ai)
            s1b = polar_to_real_imag(s1br, s1bi)
            s2a = polar_to_real_imag(s2ar, s2ai)
            s2b = polar_to_real_imag(s2br, s2bi)
            s3a = polar_to_real_imag(s3ar, s3ai)
            s3b = polar_to_real_imag(s3br, s3bi)
            s4a = polar_to_real_imag(s4ar, s4ai)
            s4b = polar_to_real_imag(s4br, s4bi)
            axial = polar_to_real_imag(axialr, axiali)
        else:
            s1a = complex(s1ar, s1ai)
            s1b = complex(s1br, s1bi)
            s2a = complex(s2ar, s2ai)
            s2b = complex(s2br, s2bi)
            s3a = complex(s3ar, s3ai)
            s3b = complex(s3br, s3bi)
            s4a = complex(s4ar, s4ai)
            s4b = complex(s4br, s4bi)
            axial = complex(axialr, axiali)

        obj.add_new_eid_sort1(dt, eid,
                              s1a, s2a, s3a, s4a, axial,
                              s1b, s2b, s3b, s4b)
    return n


def oes_cbar_random_10(op2: OP2, data: bytes,
                       obj: RandomBarStressArray | RandomBarStrainArray,
                       nelements: int, ntotal: int) -> int:
    n = 0
    #print(op2.code_information())
    #print('op2._analysis_code_fmt =', op2._analysis_code_fmt)
    struct1 = Struct(op2._endian + op2._analysis_code_fmt + b'9f')
    add_sort_x = getattr(obj, 'add_sort' + str(op2.sort_method))

    #self.log.info('op2.nonlinear_factor = %s' % op2.nonlinear_factor)
    #assert op2.sort_method == 2, op2.code_information()
    #if sort_method == 2:
        #obj.node_id = 42
    nonlinear_factor = op2.nonlinear_factor
    #print(f'CBAR: nelements={nelements}')
    for i in range(nelements):
        edata = data[n:n+ntotal]
        n += ntotal

        out = struct1.unpack(edata)
        (eid_device,
         s1a, s2a, s3a, s4a, axial,
         s1b, s2b, s3b, s4b) = out
        eid, dt = get_eid_dt_from_eid_device(
            eid_device, nonlinear_factor, op2.sort_method)
        #print(f'  eid_device={eid_device} eid={eid} dt={nonlinear_factor} nf={nonlinear_factor} -> {obj.data.shape}')
        #continue
        #print('  eid=%i; C%i=[%s]\n' % (eid, i, ', '.join(['%r' % di for di in out])))
        if op2.table_name_str == 'OESXRMS1':
            #assert sort_method == 2
            assert op2.sort_method == 1, op2.code_information()

        if op2.is_debug_file:
            op2.binary_debug.write('  eid=%i; C%i=[%s]\n' % (
                eid, i, ', '.join(['%r' % di for di in out])))

        assert eid > 0, "dt=%s eid=%s" % (dt, eid)
        add_sort_x(dt, eid,
                   s1a, s2a, s3a, s4a, axial,
                   s1b, s2b, s3b, s4b)
    return n


def oes_cbar100_real_10(op2: OP2, data: bytes, obj, nelements: int, ntotal: int, dt: Any) -> int:
    n = 0
    struct1 = Struct(op2._endian + mapfmt(op2._analysis_code_fmt + b'9f', op2.size))
    for i in range(nelements):
        edata = data[n:n+ntotal]
        out = struct1.unpack(edata)
        (eid_device, sd, sxc, sxd, sxe, sxf, axial, smax, smin, MS) = out
        eid, dt = get_eid_dt_from_eid_device(
            eid_device, op2.nonlinear_factor, op2.sort_method)
        if op2.is_debug_file:
            op2.binary_debug.write('  eid=%i; C%i=[%s]\n' % (
                eid, i, ', '.join(['%r' % di for di in out])))
        n += ntotal
        obj.add_new_eid_sort1(op2.element_name, dt, eid,
                              sd, sxc, sxd, sxe, sxf, axial, smax, smin, MS)
    return n
