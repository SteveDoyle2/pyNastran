from __future__ import annotations
from struct import Struct
from typing import Any, TYPE_CHECKING
import numpy as np

from pyNastran.op2.op2_interface.op2_reader import mapfmt
from pyNastran.op2.op2_helper import polar_to_real_imag
# from pyNastran.op2.op2_interface.utils import reshape_bytes_block_strip
from pyNastran.op2.tables.utils import get_eid_dt_from_eid_device

from pyNastran.op2.tables.oef_forces.oef_force_objects import (
    RealSolidPressureForceArray,
)
from pyNastran.op2.tables.oef_forces.oef_complex_force_objects import (
    ComplexSolidPressureForceArray,
)

if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.op2.op2 import OP2

def oef_csolid_pressure(op2: OP2, data, ndata, dt, is_magnitude_phase,
                        result_type, prefix, postfix):
    """
    76-HEXPR
    77-PENPR
    78-TETPR
    """
    if op2.element_type == 76:
        result_name = prefix + 'chexa_pressure_force' + postfix
    elif op2.element_type == 77:
        result_name = prefix + 'cpenta_pressure_force' + postfix
    elif op2.element_type == 78:
        result_name = prefix + 'ctetra_pressure_force' + postfix
    elif op2.element_type == 79:
        result_name = prefix + 'cpyram_pressure_force' + postfix
    else:
        msg = op2.code_information()
        return op2._not_implemented_or_skip(data, ndata, msg), None, None
    slot = op2.get_result(result_name)

    op2._results._found_result(result_name)
    factor = op2.factor
    if op2.format_code == 1 and op2.num_wide == 10:  # real
        ntotal = 40 * factor
        nelements = ndata // ntotal
        #nelements = ndata // ntotal

        obj_real = RealSolidPressureForceArray
        auto_return, is_vectorized = op2._create_oes_object4(
            nelements, result_name, slot, obj_real)
        if auto_return:
            return nelements * ntotal, None, None

        obj = op2.obj
        #if op2.is_debug_file:
            #op2.binary_debug.write('  [cap, element1, element2, ..., cap]\n')
            #op2.binary_debug.write('  cap = %i  # assume 1 cap when there could have been multiple\n' % ndata)
            #op2.binary_debug.write('  #elementi = [eid_device, axial, torque]\n')
            #op2.binary_debug.write('  nelements=%i; nnodes=1 # centroid\n' % nelements)
        if op2.use_vector and is_vectorized and op2.sort_method == 1:
            # self.itime = 0
            # self.ielement = 0
            # self.itotal = 0
            #self.ntimes = 0
            #self.nelements = 0
            n = nelements * ntotal
            itotal = obj.ielement
            ielement2 = obj.itotal + nelements
            itotal2 = ielement2

            floats = np.frombuffer(data, dtype=op2.fdtype8).reshape(nelements, 10)
            obj._times[obj.itime] = dt
            if obj.itime == 0:
                ints = np.frombuffer(data, dtype=op2.idtype8).reshape(nelements, 10)
                eids = ints[:, 0] // 10
                assert eids.min() > 0, eids.min()
                obj.element[itotal:itotal2] = eids

            #[axial_force, torque]
            obj.data[obj.itime, itotal:itotal2, :] = floats[:, 3:].copy()
            obj.itotal = itotal2
            obj.ielement = ielement2
        else:
            n = oef_csolid_pressure_10(op2, data, obj,
                                       nelements, ntotal, dt)

    elif op2.format_code in [2, 3] and op2.num_wide == 16:  # imag
        ntotal = 64 * factor
        nelements = ndata // ntotal
        auto_return, is_vectorized = op2._create_oes_object4(
            nelements, result_name, slot, ComplexSolidPressureForceArray)
        if auto_return:
            return nelements * ntotal, None, None

        obj = op2.obj
        #if op2.is_debug_file:
            #op2.binary_debug.write('  [cap, element1, element2, ..., cap]\n')
            #op2.binary_debug.write('  cap = %i  # assume 1 cap when there could have been multiple\n' % ndata)
            #op2.binary_debug.write('  #elementi = [eid_device, axial, torque]\n')
            #op2.binary_debug.write('  nelements=%i; nnodes=1 # centroid\n' % nelements)

        if op2.use_vector and is_vectorized and op2.sort_method == 1:
            n = nelements * ntotal
            itotal = obj.ielement
            ielement2 = obj.itotal + nelements
            itotal2 = ielement2

            floats = np.frombuffer(data, dtype=op2.fdtype8).reshape(nelements, 16).copy()
            obj._times[obj.itime] = dt
            if obj.itime == 0:
                ints = np.frombuffer(data, dtype=op2.idtype8).reshape(nelements, 16)
                eids = ints[:, 0] // 10
                assert eids.min() > 0, eids.min()
                obj.element[itotal:itotal2] = eids

            #[xaccr, yaccr, zaccr, xvelr, yvelr, zvelr, pressure,
            # xacci, yacci, zacci, xveli, yveli, zveli]
            if is_magnitude_phase:
                mag = floats[:, [3, 4, 5, 6, 7, 8, 9]]
                phase = np.hstack([
                    floats[:, [10, 11, 12, 13, 14, 15]],
                    np.zeros((len(floats), 1), dtype='float32')
                ])
                rtheta = np.radians(phase)
                real_imag = mag * (np.cos(rtheta) + 1.j * np.sin(rtheta))
            else:
                real = floats[:, [3, 4, 5, 6, 7, 8, 9]]
                imag = np.hstack([
                    floats[:, [10, 11, 12, 13, 14, 15]],
                    np.zeros((len(floats), 1), dtype='float32')
                ])
                real_imag = real + 1.j * imag
            obj.data[obj.itime, itotal:itotal2, :] = real_imag
            obj.itotal = itotal2
            obj.ielement = ielement2
        else:
            n = oef_csolid_imag_16(op2, data, obj,
                                   nelements, ntotal,
                                   is_magnitude_phase)
    else:  # pragma: no cover
        raise RuntimeError(op2.code_information())
        #msg = op2.code_information()
        #return op2._not_implemented_or_skip(data, ndata, msg), None, None
    return n, nelements, ntotal


def oef_csolid_pressure_10(self, data: bytes,
                           obj: RealSolidPressureForceArray,
                           nelements: int, ntotal: int, dt: Any) -> int:
    op2 = self
    n = 0
    if self.size == 4:
        fmt = op2._endian + op2._analysis_code_fmt + b'8s7f'
    else:
        fmt = op2._endian + mapfmt(op2._analysis_code_fmt, self.size) + b'16s7d'

    s = Struct(fmt)
    add_sort_x = getattr(obj, 'add_sort' + str(op2.sort_method))
    for unused_i in range(nelements):
        edata = data[n : n + ntotal]
        n += ntotal
        out = s.unpack(edata)
        if op2.is_debug_file:
            op2.binary_debug.write('OEF_PentaPressure-%s %s\n' % (op2.element_type, str(out)))
        (eid_device, ename, ax, ay, az, vx, vy, vz, pressure) = out
        eid, dt = get_eid_dt_from_eid_device(
            eid_device, op2.nonlinear_factor, op2.sort_method)
        add_sort_x(dt, eid, ename, ax, ay, az, vx, vy, vz, pressure)
    return n

def oef_csolid_imag_16(self, data: bytes,
                       obj: ComplexSolidPressureForceArray,
                       nelements: int, ntotal: int,
                       is_magnitude_phase: bool) -> int:
    op2 = self
    n = 0
    if self.size == 4:
        s = Struct(op2._endian + op2._analysis_code_fmt + b'8s 13f')
    else:
        s = Struct(mapfmt(op2._endian + op2._analysis_code_fmt, self.size) + b'16s 13d')

    add_sort_x = getattr(obj, 'add_sort' + str(op2.sort_method))
    for unused_i in range(nelements):
        edata = data[n:n+ntotal]
        n += ntotal

        #print(len(edata))
        out = s.unpack(edata)
        if op2.is_debug_file:
            op2.binary_debug.write('_oef_csolid_pressure-%s %s\n' % (op2.element_type, str(out)))
        (eid_device, ename,
         axr, ayr, azr, vxr, vyr, vzr, pressure,
         axi, ayi, azi, vxi, vyi, vzi) = out
        eid, dt = get_eid_dt_from_eid_device(
            eid_device, op2.nonlinear_factor, op2.sort_method)
        ename = ename.decode('utf-8').strip()

        if is_magnitude_phase:
            ax = polar_to_real_imag(axr, axi)
            vx = polar_to_real_imag(vxr, vxi)
            ay = polar_to_real_imag(ayr, ayi)
            vy = polar_to_real_imag(vyr, vyi)
            az = polar_to_real_imag(azr, azi)
            vz = polar_to_real_imag(vzr, vzi)
        else:
            ax = complex(axr, axi)
            vx = complex(vxr, vxi)
            ay = complex(ayr, ayi)
            vy = complex(vyr, vyi)
            az = complex(azr, azi)
            vz = complex(vzr, vzi)
        cpressure = complex(pressure, 0.)
        add_sort_x(dt, eid, ename, ax, ay, az, vx, vy, vz, cpressure)
    return n
