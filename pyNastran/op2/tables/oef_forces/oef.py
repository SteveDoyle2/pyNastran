#pylint: disable=C0301,W0212,C0103,W0201
"""
Defines the Real/Complex Forces created by:
    FORCE = ALL

NX Case Control  Block         Description
===============  ==========    ===========
FORCE            OEF1          Element forces or heat flux (linear elements only)
FORCE            OEF1X         Element forces (nonlinear elements only)
????             HOEF1         ???
FORCE            DOEF1         Scaled Response Spectra
MODCON           OEFMC         Modal contributions
FORCE            OEF1X         Element forces with intermediate (CBAR and CBEAM)
                               station forces and forces on nonlinear elements
FLUX             HOEF1         Element heat flux

"""
from __future__ import annotations
from struct import Struct
from typing import TYPE_CHECKING

import numpy as np

from pyNastran.op2.op2_interface.function_codes import func1, func7
from pyNastran.op2.op2_interface.op2_reader import mapfmt
from pyNastran.op2.tables.utils import get_is_slot_saved, get_eid_dt_from_eid_device
from pyNastran.op2.op2_interface.msc_tables import MSC_OEF_REAL_MAPPER, MSC_OEF_IMAG_MAPPER
from pyNastran.op2.op2_interface.nx_tables import NX_OEF_REAL_MAPPER, NX_OEF_IMAG_MAPPER
from pyNastran.op2.op2_interface.op2_codes import SORT1_TABLES_BYTES, TABLES_BYTES


from pyNastran.op2.tables.oef_forces.utils_celas_cdamp import oef_celas_cdamp
from pyNastran.op2.tables.oef_forces.utils_cbar import oef_cbar_34, oef_cbar_100
from pyNastran.op2.tables.oef_forces.utils_crod import oef_crod
from pyNastran.op2.tables.oef_forces.utils_cvisc import oef_cvisc
from pyNastran.op2.tables.oef_forces.utils_cgap import oef_cgap
from pyNastran.op2.tables.oef_forces.utils_cbend import oef_cbend
from pyNastran.op2.tables.oef_forces.utils_cbush import oef_cbush
from pyNastran.op2.tables.oef_forces.utils_cshear import oef_cshear
from pyNastran.op2.tables.oef_forces.utils_cbeam import oef_cbeam
from pyNastran.op2.tables.oef_forces.utils_shells import oef_shells_centroidal
from pyNastran.op2.tables.oef_forces.utils_shells_nodal import oef_shells_nodal
from pyNastran.op2.tables.oef_forces.utils_composite_plates import oef_shells_composite
from pyNastran.op2.tables.oef_forces.utils_solid import oef_csolid_pressure
from pyNastran.op2.tables.oef_forces.utils_cconeax import oef_cconeax

from pyNastran.op2.tables.oef_forces.oef_thermal_objects import (
    Real1DHeatFluxArray,
    RealHeatFlux_2D_3DArray,
    RealChbdyHeatFluxArray,
    RealConvHeatFluxArray,
)
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.op2.op2 import OP2


class OEF:
    """Defines OEFx table reading for element forces/heat flux"""
    def __init__(self, op2: OP2):
        self.op2 = op2

    def _oef_force_code(self):
        """
        Gets the numwide codes for the element to determine if
        the real or complex result should be found.
        The format and sort codes do not always give the right answer...
        """
        op2 = self.op2
        if op2.is_nx:
            real_mapper = NX_OEF_REAL_MAPPER
            imag_mapper = NX_OEF_IMAG_MAPPER
        else:
            real_mapper = MSC_OEF_REAL_MAPPER
            imag_mapper = MSC_OEF_IMAG_MAPPER


        try:
            real = real_mapper[op2.element_type]
        except KeyError:
            real = None

        try:
            imag = imag_mapper[op2.element_type]
        except KeyError:
            imag = None
        return real, imag

    def _read_oef1_3(self, data: bytes, ndata: int):
        """Table 3 parser for OEF1 table"""
        op2 = self.op2
        op2._analysis_code_fmt = b'i'
        op2._data_factor = 1
        op2.words = [
            'aCode', 'tCode', 'element_type', 'isubcase',
            '???', '???', '???', '???',
            'format_code', 'num_wide', 'o_code', '???',
            '???', '???', '???', '???',
            '???', '???', '???', '???',
            '???', '???', '???', '???',
            '???', 'Title', 'subtitle', 'label']

        op2.parse_approach_code(data)

        #: element type
        op2.element_type = op2.add_data_parameter(data, 'element_type', b'i', 3, False)

        # dynamic load set ID/random code
        #self.dLoadID = op2.add_data_parameter(data, 'dLoadID', b'i', 8, False)

        #: format code
        op2.format_code = op2.add_data_parameter(data, 'format_code', b'i', 9, False)

        #: number of words per entry in record
        #: .. note: is this needed for this table ???
        op2.num_wide = op2.add_data_parameter(data, 'num_wide', b'i', 10, False)

        #: undefined in DMAP...
        op2.o_code = op2.add_data_parameter(data, 'o_code', b'i', 11, False)

        #: thermal flag; 1 for heat ransfer, 0 otherwise
        op2.thermal = op2.add_data_parameter(data, 'thermal', b'i', 23, False)

        ## assuming tCode=1
        if op2.analysis_code == 1:   # statics
            op2.loadID = op2.add_data_parameter(data, 'loadID', b'i', 5, False)  # load set ID number
            op2.data_names = op2.apply_data_code_value('data_names', ['loadID'])
            op2.setNullNonlinearFactor()
        elif op2.analysis_code == 2:  # normal modes/buckling (real eigenvalues)
            #: mode number
            op2.mode = op2.add_data_parameter(data, 'mode', b'i', 5)
            #: eigenvalue
            op2.eign = op2.add_data_parameter(data, 'eign', b'f', 6, False)
            op2.cycle = 0.
            op2._op2_readers.reader_oug.update_mode_cycle('cycle')
            op2.data_names = op2.apply_data_code_value('data_names', ['mode', 'eign', 'cycle'])
            # TODO: mode_cycle is not defined?
            #op2.data_names = op2.apply_data_code_value('data_names', ['mode', 'eign', 'mode_cycle'])
        elif op2.analysis_code == 3:  # differential stiffness 0
            #: load set ID number
            op2.loadID = op2.add_data_parameter(data, 'loadID', b'i', 5)
            op2.data_names = op2.apply_data_code_value('data_names', ['loadID'])
        elif op2.analysis_code == 4:  # differential stiffness 1
            #: load set ID number
            op2.loadID = op2.add_data_parameter(data, 'loadID', b'i', 5)
            op2.data_names = op2.apply_data_code_value('data_names', ['loadID'])
        elif op2.analysis_code == 5:   # frequency
            self.freq = op2.add_data_parameter(data, 'freq', b'f', 5)  # frequency
            op2.data_names = op2.apply_data_code_value('data_names', ['freq'])
        elif op2.analysis_code == 6:  # transient
            self.time = op2.add_data_parameter(data, 'time', b'f', 5)  # time step
            op2.data_names = op2.apply_data_code_value('data_names', ['time'])
        elif op2.analysis_code == 7:  # pre-buckling
            #: load set ID number
            op2.loadID = op2.add_data_parameter(data, 'loadID', b'i', 5)
            #op2.apply_data_code_value('data_names',['lsdvmn'])
            op2.data_names = op2.apply_data_code_value('data_names', ['loadID'])
        elif op2.analysis_code == 8:  # post-buckling
            #: load set ID number
            op2.loadID = op2.add_data_parameter(data, 'loadID', b'i', 5)
            #: real eigenvalue
            op2.eigr = op2.add_data_parameter(data, 'eigr', b'f', 6, False)
            op2.data_names = op2.apply_data_code_value('data_names', ['loadID', 'eigr'])
        elif op2.analysis_code == 9:  # complex eigenvalues
            #: mode number
            op2.mode = op2.add_data_parameter(data, 'mode', b'i', 5)
            #: real eigenvalue
            op2.eigr = op2.add_data_parameter(data, 'eigr', b'f', 6, False)
            #: imaginary eigenvalue
            op2.eigi = op2.add_data_parameter(data, 'eigi', b'f', 7, False)
            op2.data_names = op2.apply_data_code_value('data_names', ['mode', 'eigr', 'eigi'])
        elif op2.analysis_code == 10:  # nonlinear statics
            #: load step
            self.load_step = op2.add_data_parameter(data, 'load_step', b'f', 5)
            op2.data_names = op2.apply_data_code_value('data_names', ['load_step'])
        elif op2.analysis_code == 11:  # geometric nonlinear statics
            #: load set ID number
            op2.loadID = op2.add_data_parameter(data, 'loadID', b'i', 5)
            op2.data_names = op2.apply_data_code_value('data_names', ['loadID'])
        else:
            msg = 'invalid analysis_code...analysis_code=%s' % str(op2.analysis_code)
            raise RuntimeError(msg)

        op2.fix_format_code()
        op2._parse_thermal_code()
        self._set_force_stress_element_name()

        if op2.is_debug_file:
            op2.binary_debug.write('  %-14s = %r\n' % ('element_name', op2.element_name))
            op2.binary_debug.write('  %-14s = %r %s\n' % ('approach_code', op2.approach_code,
                                                           op2.approach_code_str(op2.approach_code)))
            op2.binary_debug.write('  %-14s = %r\n' % ('tCode', op2.tCode))
            op2.binary_debug.write('  %-14s = %r\n' % ('isubcase', op2.isubcase))

        op2._read_title(data)
        if op2.element_type not in op2.element_mapper:
            msg = 'element_type = %s' % op2.element_type
            return op2._not_implemented_or_skip(data, ndata, msg)
        op2._write_debug_bits()
        assert op2.num_wide != 146, op2.code_information()
        #print('OEF-%s' % op2.element_name)
        #self._check_result_type()

    def _set_force_stress_element_name(self):
        """
        Not all cards can have OES/OEF output, so if they do,
        we have in the wrong solver, specifically:
        - RBAR

        """
        op2 = self.op2
        try:
            op2.element_name = op2.element_mapper[op2.element_type]
        except KeyError:
            op2.log.error(op2.code_information())
            raise

        if op2.element_type == 227 and op2.element_name == 'RBAR' and op2.is_msc:
            op2.to_nx(' because element_type=227 was found')
            op2.element_name = op2.element_mapper[op2.element_type]
        assert op2.element_name != '', op2.code_information()

        op2.data_code['element_name'] = op2.element_name

    def _read_oef2_3(self, data, unused_ndata):
        """Table 3 parser for OEF2 table"""
        op2 = self.op2
        op2._data_factor = 1
        op2.words = [
            'aCode', 'tCode', 'element_type', 'isubcase',
            '???', '???', '???', '???',
            'format_code', 'num_wide', 'o_code', '???',
            '???', '???', '???', '???',
            '???', '???', '???', '???',
            '???', '???', '???', '???',
            '???', 'Title', 'subtitle', 'label']

        op2.parse_approach_code(data)  # 3
        op2.sort_method = 2

        #: element type
        op2.element_type = op2.add_data_parameter(data, 'element_type', b'i', 3, False)

        # dynamic load set ID/random code
        #self.dLoadID = op2.add_data_parameter(data, 'dLoadID', b'i', 8, False)

        #: format code
        op2.format_code = op2.add_data_parameter(data, 'format_code', b'i', 9, False)

        #: number of words per entry in record
        #: .. note: is this needed for this table ???
        op2.num_wide = op2.add_data_parameter(data, 'num_wide', b'i', 10, False)

        #: undefined in DMAP...
        op2.o_code = op2.add_data_parameter(data, 'o_code', b'i', 11, False)

        #: thermal flag; 1 for heat ransfer, 0 otherwise
        op2.thermal = op2.add_data_parameter(data, 'thermal', b'i', 23, False)

        op2.element_id = op2.add_data_parameter(data, 'element_id', b'i', 5, fix_device_code=True)
        op2._element_id = op2.add_data_parameter(data, '_element_id', b'f', 5, apply_nonlinear_factor=False, add_to_dict=True)

        if op2.analysis_code == 1:  # static...because reasons.
            op2._analysis_code_fmt = b'f'
            op2.data_names = op2.apply_data_code_value('data_names', ['element_id'])
            #op2.apply_data_code_value('analysis_method', 'static')
        elif op2.analysis_code == 2:  # real eigenvalues
            op2._analysis_code_fmt = b'i'
            op2.eign = op2.add_data_parameter(data, 'eign', b'f', 6, False)
            op2.mode_cycle = op2.add_data_parameter(data, 'mode_cycle', b'i', 7, False)  # mode or cycle .. todo:: confused on the type - F1???
            op2.data_names = op2.apply_data_code_value('data_names', ['element_id', 'eign', 'mode_cycle'])
        elif op2.analysis_code == 5:   # frequency
            op2._analysis_code_fmt = b'f'
            op2.data_names = op2.apply_data_code_value('data_names', ['element_id'])
            op2.apply_data_code_value('analysis_method', 'freq')
        elif op2.analysis_code == 6:  # transient
            op2._analysis_code_fmt = b'f'
            op2.data_names = op2.apply_data_code_value('data_names', ['element_id'])
            op2.apply_data_code_value('analysis_method', 'time')
        elif op2.analysis_code == 7:  # pre-buckling
            op2._analysis_code_fmt = b'i'
            op2.data_names = op2.apply_data_code_value('data_names', ['element_id'])
            op2.apply_data_code_value('analysis_method', 'lsdvmn')
        elif op2.analysis_code == 8:  # post-buckling
            op2._analysis_code_fmt = b'f'
            op2.eigr = op2.add_data_parameter(data, 'eigr', b'f', 6, False)
            op2.data_names = op2.apply_data_code_value('data_names', ['element_id', 'eigr'])
            op2.apply_data_code_value('analysis_method', 'eigr')
        elif op2.analysis_code == 9:  # complex eigenvalues
            # mode number
            op2._analysis_code_fmt = b'i'
            op2.eigr = op2.add_data_parameter(data, 'eigr', b'f', 6, False)
            op2.eigi = op2.add_data_parameter(data, 'eigi', b'f', 7, False)
            op2.data_names = op2.apply_data_code_value('data_names', ['element_id', 'eigr', 'eigi'])
            op2.apply_data_code_value('analysis_method', 'mode')
        elif op2.analysis_code == 10:  # nonlinear statics
            # load step
            op2._analysis_code_fmt = b'f'
            op2.data_names = op2.apply_data_code_value('data_names', ['element_id'])
            op2.apply_data_code_value('analysis_method', 'lftsfq')
        elif op2.analysis_code == 11:  # old geometric nonlinear statics
            # load set number
            op2.data_names = op2.apply_data_code_value('data_names', ['element_id'])
        elif op2.analysis_code == 12:  # contran ? (may appear as aCode=6)  --> straight from DMAP...grrr...
            op2.data_names = op2.apply_data_code_value('data_names', ['element_id'])
        else:
            msg = 'invalid analysis_code...analysis_code=%s' % op2.analysis_code
            raise RuntimeError(msg)

        op2.fix_format_code()
        op2._op2_readers.reader_oes._fix_oes_sort2(data)
        self._set_force_stress_element_name()
        #assert isinstance(op2.nonlinear_factor, int), op2.nonlinear_factor
        #self._check_result_type()

    def _read_oef1_4(self, data: bytes, ndata: int):
        """Table 4 parser for OEF1 table"""
        op2 = self.op2
        if op2.thermal == 0:
            op2._setup_op2_subcase('FORCE')
            n = self._read_oef1_loads(data, ndata)
        elif op2.thermal == 1:
            n = self._read_oef1_thermal(data, ndata)
        elif op2.thermal in [2, 4, 8]: # 2=ABS, 4=SRSS, 8=NRL
            #C:\NASA\m4\formats\git\examples\move_tpl\ms103.op2 # SRSS
            n = self._read_oef1_loads(data, ndata)
        else:
            n = op2._not_implemented_or_skip(data, ndata, 'thermal=%s' % op2.thermal)
        return n

    def _read_oef2_4(self, data: bytes, ndata: int):
        op2 = self.op2
        if op2.thermal == 0: # and op2.element_type not in [77]:
            op2._setup_op2_subcase('FORCE')
            n = self._read_oef1_loads(data, ndata)
        else:
            n = op2._not_implemented_or_skip(data, ndata, 'thermal=%s' % op2.thermal)
        assert n is not None, op2.code_information()
        return n

    def _read_oef1_thermal(self, data: bytes, ndata: int):
        """Table 4 parser for OEF1 thermal table"""
        op2 = self.op2
        if op2._results.is_not_saved('force'):
            # self.log.debug('skipping OEF due to force')
            return ndata
        prefix, postfix = get_oef_prefix_postfix(op2)

        n = 0
        #thermal
        #is_magnitude_phase = op2.is_magnitude_phase()
        dt = op2.nonlinear_factor

        #flag = 'element_id'
        if op2.element_type in [1, 2, 3, 10, 34, 69]:  # ROD,BEAM,TUBE,CONROD,BAR,BEND
            n, nelements, ntotal = self._thermal_1d_heat_flux(data, ndata, dt, prefix, postfix)
            if nelements is None:
                return n

        elif op2.element_type in [33, 53, 64, 74, 75,  # CQUAD4, CTRIAX6, CQUAD8, CTRIA3, CTRIA6
                                   39, 67, 68]:  # TETRA, HEXA, PENTA
            n, nelements, ntotal = self._thermal_2d_3d_heat_flux(data, ndata, dt, prefix, postfix)

        elif op2.element_type in [107, 108, 109]:  # CHBDYE, CHBDYG, CHBDYP
            n, nelements, ntotal = self._thermal_chbdy(data, ndata, dt, prefix, postfix)

        elif op2.element_type == 110:  # CONV
            n, nelements, ntotal = self._thermal_conv(data, ndata, dt, prefix, postfix)

        elif op2.element_type in [145, 146, 147,  # VUHEXA,VUPENTA,VUTETRA
                                  189, 190,  # VUQUAD,VUTRIA
                                  191]:  # VUBEAM
            # removed by msc/nx
            msg = f'{op2.table_name_str} {op2.element_name}-{op2.element_type} has been removed'
            op2.log.warning(msg)
            return ndata
            # return op2._not_implemented_or_skip(data, ndata, msg)
        else:
            msg = 'OEF sort1 thermal Type=%s num=%s' % (op2.element_name, op2.element_type)
            return op2._not_implemented_or_skip(data, ndata, msg)
        if nelements is None:
            return n


        assert op2.thermal == 1, op2.thermal
        assert ndata > 0, ndata
        try:
            assert nelements > 0, 'nelements=%r element_type=%s element_name=%r\n%s' % (nelements, op2.element_type, op2.element_name, op2.code_information())
        except UnboundLocalError:
            raise UnboundLocalError('element_name=%r' % op2.element_name)
        #assert ndata % ntotal == 0, '%s n=%s nwide=%s len=%s ntotal=%s' % (op2.element_name, ndata % ntotal, ndata % op2.num_wide, ndata, ntotal)
        assert op2.num_wide * 4 == ntotal, 'numwide*4=%s ntotal=%s' % (op2.num_wide*4, ntotal)
        assert n > 0, 'n=%s element_type=%s element_name=%s numwide=%s' % (
            n, op2.element_type, op2.element_name, op2.num_wide)
        return n

    def _thermal_1d_heat_flux(self, data, ndata, dt, prefix, postfix):
        """
        1-CROD
        2-CBEAM
        3-CTUBE
        10-CONROD
        34-CBAR
        69-CBEND
        """
        op2 = self.op2
        n = 0
        obj_vector_real = Real1DHeatFluxArray
        #if op2.element_type == 1: # CROD
        element_type = op2.element_type
        if element_type == 1:
            result_name = prefix + 'crod_thermal_load' + postfix
        elif element_type == 2:
            result_name = prefix + 'cbeam_thermal_load' + postfix
        elif element_type == 3:
            result_name = prefix + 'ctube_thermal_load' + postfix
        elif element_type == 10:
            result_name = prefix + 'conrod_thermal_load' + postfix
        elif element_type == 34:
            result_name = prefix + 'cbar_thermal_load' + postfix
        elif element_type == 69:
            result_name = prefix + 'cbend_thermal_load' + postfix
        else:
            raise NotImplementedError('element_type=%s element_name=%s' % (
                element_type, op2.element_name))

        is_saved, slot = get_is_slot_saved(op2, result_name)
        if not is_saved:
            return ndata, None, None

        factor = op2.factor
        if op2.format_code == 1 and op2.num_wide == 9:  # real
            ntotal = 36 * factor
            nelements = ndata // ntotal
            auto_return, is_vectorized = op2._create_oes_object4(
                nelements, result_name, slot, obj_vector_real)
            if auto_return:
                return nelements * op2.num_wide * 4, None, None
            obj = op2.obj
            #if op2.is_debug_file:
                #op2.binary_debug.write('  [cap, element1, element2, ..., cap]\n')
                #op2.binary_debug.write('  cap = %i  # assume 1 cap when there could have been multiple\n' % ndata)
                #op2.binary_debug.write('  #elementi = [eid_device, axial, torque]\n')
                #op2.binary_debug.write('  nelements=%i; nnodes=1 # centroid\n' % nelements)

            if op2.use_vector and is_vectorized and op2.sort_method == 1:
                n = nelements * 4 * op2.num_wide
                itotal = obj.ielement
                ielement2 = obj.itotal + nelements
                itotal2 = ielement2

                floats = np.frombuffer(data, dtype=op2.fdtype8).reshape(nelements, 9)
                obj._times[obj.itime] = dt

                strings = np.frombuffer(data, dtype=op2._uendian + 'S4').reshape(nelements, 9)
                s = np.array([s1+s2 for s1, s2 in zip(strings[:, 1], strings[:, 2])])
                #print(s)
                #print('itime = ', obj.itime)
                #print('---------')
                if obj.itime == 0:
                    ints = np.frombuffer(data, dtype=op2.idtype8).reshape(nelements, 9)
                    eids = ints[:, 0] // 10
                    assert eids.min() > 0, eids.min()
                    obj.element[itotal:itotal2] = eids
                    obj.element_data_type[itotal:itotal2] = s
                    #obj.element_type[obj.itime, itotal:itotal2, :] = strings[:, 3:]

                #[etype, xgrad, ygrad, zgrad, xflux, yflux, zflux]
                obj.data[obj.itime, itotal:itotal2, :] = floats[:, 3:].copy()
                obj.itotal = itotal2
                obj.ielement = ielement2
            else:
                if op2.size == 4:
                    s = Struct(op2._endian + op2._analysis_code_fmt + b'8s6f')
                else:
                    s = Struct(op2._endian + op2._analysis_code_fmt + b'16s6d')
                add_sort_x = getattr(obj, 'add_sort' + str(op2.sort_method))
                for unused_i in range(nelements):
                    edata = data[n:n+ntotal]
                    out = s.unpack(edata)
                    (eid_device, etype, xgrad, ygrad, zgrad, xflux, yflux, zflux) = out
                    eid, dt = get_eid_dt_from_eid_device(
                        eid_device, op2.nonlinear_factor, op2.sort_method)
                    add_sort_x(dt, eid, etype, xgrad, ygrad, zgrad, xflux, yflux, zflux)
                    n += ntotal
        else:  # pragma: no cover
            msg = op2.code_information()
            return op2._not_implemented_or_skip(data, ndata, msg), None, None
        return n, nelements, ntotal

    def _thermal_2d_3d_heat_flux(self, data, ndata, dt, prefix, postfix):
        """
        33-QUAD4-centroidal
        53-TRIAX6
        64-QUAD8
        74-TRIA3
        75-TRIA6

        39-TETRA
        67-HEXA
        68-PENTA
        """
        op2 = self.op2
        n = 0
        if op2.element_type == 33:
            result_name = prefix + 'cquad4_thermal_load' + postfix
        elif op2.element_type == 53:
            result_name = prefix + 'ctriax6_thermal_load' + postfix
        elif op2.element_type == 64:
            result_name = prefix + 'cquad8_thermal_load' + postfix
        elif op2.element_type == 74:
            result_name = prefix + 'ctria3_thermal_load' + postfix
        elif op2.element_type == 75:
            result_name = prefix + 'ctria6_thermal_load' + postfix
        elif op2.element_type == 39:
            result_name = prefix + 'ctetra_thermal_load' + postfix
        elif op2.element_type == 67:
            result_name = prefix + 'chexa_thermal_load' + postfix
        elif op2.element_type == 68:
            result_name = prefix + 'cpenta_thermal_load' + postfix
        else:
            raise NotImplementedError('element_type=%s element_name=%s' % (
                op2.element_type, op2.element_name))

        if op2._results.is_not_saved(result_name):
            return ndata, None, None

        obj_vector_real = RealHeatFlux_2D_3DArray
        #if op2.element_type == 1: # CROD
        #result_name = 'thermalLoad_2D_3D'

        is_saved, slot = get_is_slot_saved(op2, result_name)
        if not is_saved:
            return ndata, None, None

        factor = op2.factor
        size = op2.size
        if op2.format_code == 1 and op2.num_wide == 9:  # real - 2D
            # [33, 53, 64, 74, 75]
            ntotal = 36 * factor
            nelements = ndata // ntotal
            auto_return, is_vectorized = op2._create_oes_object4(
                nelements, result_name, slot, obj_vector_real)
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

                floats = np.frombuffer(data, dtype=op2.fdtype8).reshape(nelements, 9)
                obj._times[obj.itime] = dt
                #if obj.itime == 0:
                ints = np.frombuffer(data, dtype=op2.idtype8).reshape(nelements, 9)
                eids = ints[:, 0] // 10
                assert eids.min() > 0, eids.min()
                obj.element[itotal:itotal2] = eids
                strings = np.frombuffer(data, dtype=op2._uendian + 'S4').reshape(nelements, 9)
                obj.element_data_type[itotal:itotal2] = np.array([s1+s2 for s1, s2 in zip(strings[:, 1], strings[:, 2])])

                #[etype, xgrad, ygrad, zgrad, xflux, yflux, zflux]
                obj.data[obj.itime, itotal:itotal2, :] = floats[:, 3:].copy()
                obj.itotal = itotal2
                obj.ielement = ielement2
            else:
                # no zed on this element for some reason...
                if size == 4:
                    fmt = op2._endian + op2._analysis_code_fmt + b'8s 6f'
                else:
                    fmt = op2._endian + mapfmt(op2._analysis_code_fmt, size) + b'16s 6d'
                s = Struct(fmt)
                add_sort_x = getattr(obj, 'add_sort' + str(op2.sort_method))
                for unused_i in range(nelements):
                    edata = data[n:n+ntotal]
                    n += ntotal
                    out = s.unpack(edata)
                    (eid_device, etype, xgrad, ygrad, zgrad, xflux, yflux, zflux) = out
                    eid, dt = get_eid_dt_from_eid_device(
                        eid_device, op2.nonlinear_factor, op2.sort_method)
                    add_sort_x(dt, eid, etype, xgrad, ygrad, zgrad, xflux, yflux, zflux)

        elif op2.format_code == 1 and op2.num_wide == 10:  # real - 3D
            # [39, 67, 68]:  # HEXA,PENTA
            ntotal = 40 * factor
            nelements = ndata // ntotal
            auto_return, is_vectorized = op2._create_oes_object4(
                nelements, result_name, slot, obj_vector_real)
            if auto_return:
                return nelements * ntotal, None, None
            obj = op2.obj
            assert nelements > 0, 'ndata=%s ntotal=%s' % (ndata, ntotal)
            if op2.use_vector and is_vectorized and op2.sort_method == 1:
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
                    strings = np.frombuffer(data, dtype=op2._uendian + 'S4').reshape(nelements, 10)
                    obj.element_data_type[itotal:itotal2] = np.array([s1+s2 for s1, s2 in zip(strings[:, 1], strings[:, 2])])

                #[etype, xgrad, ygrad, zgrad, xflux, yflux, zflux, zed]
                obj.data[obj.itime, itotal:itotal2, :] = floats[:, 3:-1].copy()
                obj.itotal = itotal2
                obj.ielement = ielement2
            else:
                s = Struct(op2._endian + op2._analysis_code_fmt + b'8s6fi')
                add_sort_x = getattr(obj, 'add_sort' + str(op2.sort_method))
                for unused_i in range(nelements):
                    edata = data[n:n+ntotal]
                    n += ntotal
                    out = s.unpack(edata)
                    (eid_device, etype, xgrad, ygrad, zgrad, xflux, yflux, zflux, unused_zed) = out
                    eid, dt = get_eid_dt_from_eid_device(
                        eid_device, op2.nonlinear_factor, op2.sort_method)
                    add_sort_x(dt, eid, etype, xgrad, ygrad, zgrad, xflux, yflux, zflux)
        else:
            raise RuntimeError(op2.code_information())
        return n, nelements, ntotal

    def _thermal_chbdy(self, data, ndata, dt, prefix, postfix):
        """
        107-CHBDYE
        108-CHBDYG
        109-CHBDYP

        """
        op2 = self.op2
        n = 0
        #if op2.table_name in ['OEF1X']:
        element_type = op2.element_type
        if element_type == 107:
            result_name = prefix + 'chbdye_thermal_load' + postfix
        elif element_type == 108:
            result_name = prefix + 'chbdyg_thermal_load' + postfix
        elif element_type == 109:
            result_name = prefix + 'chbdyp_thermal_load' + postfix
        else:
            raise NotImplementedError('element_type=%s element_name=%s' % (
                element_type, op2.element_name))

        if op2._results.is_not_saved(result_name):
            return ndata, None, None

        #elif op2.table_name in ['HOEF1']:
            #if op2.element_type == 107:
                #result_name = 'chbdye_thermal_flux'
            #elif op2.element_type == 108:
                #result_name = 'chbdyg_thermal_flux'
            #elif op2.element_type == 109:
                #result_name = 'chbdyp_thermal_flux'
            #else:
                #raise NotImplementedError('element_type=%s element_name=%s' % (
                    #op2.element_type, op2.element_name))
        #else:
            #raise NotImplementedError(msg)

        factor = op2.factor
        size = op2.size
        if op2.format_code == 1 and op2.num_wide == 8:  # real
            #result_name = 'thermalLoad_CHBDY'
            is_saved, slot = get_is_slot_saved(op2, result_name)
            if not is_saved:
                return ndata, None, None

            if op2.format_code == 1 and op2.num_wide == 8:  # real
                obj_vector_real = RealChbdyHeatFluxArray
                ntotal = 32 * factor
                nelements = ndata // ntotal
                auto_return, is_vectorized = op2._create_oes_object4(
                    nelements, result_name, slot, obj_vector_real)
                if auto_return:
                    return nelements * op2.num_wide * 4, None, None
                obj = op2.obj
                #if op2.is_debug_file:
                    #op2.binary_debug.write('  [cap, element1, element2, ..., cap]\n')
                    #op2.binary_debug.write('  cap = %i  # assume 1 cap when there could have been multiple\n' % ndata)
                    #op2.binary_debug.write('  #elementi = [eid_device, etype, fapplied, free_conv, force_conv, frad, ftotal]\n')
                    #op2.binary_debug.write('  nelements=%i; nnodes=1 # centroid\n' % nelements)

                if op2.use_vector and is_vectorized and op2.sort_method == 1:
                    n = nelements * ntotal
                    itotal = obj.ielement
                    ielement2 = obj.itotal + nelements
                    itotal2 = ielement2

                    floats = np.frombuffer(data, dtype=op2.fdtype8).reshape(nelements, 8)
                    obj._times[obj.itime] = dt
                    if obj.itime == 0:
                        ints = np.frombuffer(data, dtype=op2.idtype8).reshape(nelements, 8)
                        eids = ints[:, 0] // 10
                        assert eids.min() > 0, eids.min()
                        obj.element[itotal:itotal2] = eids
                        #obj.element_type[obj.itime, itotal:itotal2, :] = strings[:, 3:]

                    #[fapplied, free_conv, force_conv, frad, ftotal]
                    obj.data[obj.itime, itotal:itotal2, :] = floats[:, 3:].copy()
                    obj.itotal = itotal2
                    obj.ielement = ielement2
                else:
                    if size == 4:
                        s1 = Struct(op2._endian + op2._analysis_code_fmt + b'8s5f')
                    else:
                        s1 = Struct(op2._endian + mapfmt(op2._analysis_code_fmt, 8) + b'16s5d')
                    add_sort_x = getattr(obj, 'add_sort' + str(op2.sort_method))
                    for unused_i in range(nelements):
                        edata = data[n:n+ntotal]
                        n += ntotal
                        out = s1.unpack(edata)
                        (eid_device, etype, fapplied, free_conv, force_conv, frad, ftotal) = out
                        eid, dt = get_eid_dt_from_eid_device(
                            eid_device, op2.nonlinear_factor, op2.sort_method)

                        if op2.is_debug_file:
                            op2.binary_debug.write('  %s -> [%s, %s, %s, %s, %s, %s, %s]\n'
                                                    % (eid, eid_device, etype, fapplied, free_conv, force_conv, frad, ftotal))
                        add_sort_x(dt, eid, etype, fapplied, free_conv, force_conv, frad, ftotal)
        else:  # pragma: no cover
            msg = op2.code_information()
            return op2._not_implemented_or_skip(data, ndata, msg), None, None
        return n, nelements, ntotal

    def _thermal_conv(self, data, ndata, dt, prefix, postfix):
        op2 = self.op2
        n = 0
        # 110-CONV
        result_name = prefix + 'conv_thermal_load' + postfix
        is_saved, slot = get_is_slot_saved(op2, result_name)
        if not is_saved:
            return ndata, None, None

        factor = op2.factor
        if op2.format_code == 1 and op2.num_wide == 4:
            ntotal = 16 * factor
            nelements = ndata // ntotal

            auto_return, is_vectorized = op2._create_oes_object4(
                nelements, result_name, slot, RealConvHeatFluxArray)
            if auto_return:
                return nelements * ntotal, None, None
            obj = op2.obj
            #if op2.is_debug_file:
                #op2.binary_debug.write('  [cap, element1, element2, ..., cap]\n')
                #op2.binary_debug.write('  cap = %i  # assume 1 cap when there could have been multiple\n' % ndata)
                #op2.binary_debug.write('  #elementi = [eid_device, etype, fapplied, free_conv, force_conv, frad, ftotal]\n')
                #op2.binary_debug.write('  nelements=%i; nnodes=1 # centroid\n' % nelements)

            if op2.use_vector and is_vectorized and op2.sort_method == 1:
                n = nelements * ntotal
                ielement = obj.ielement
                ielement2 = ielement + nelements

                floats = np.frombuffer(data, dtype=op2.fdtype8).reshape(nelements, 4).copy()
                obj._times[obj.itime] = dt
                if obj.itime == 0:
                    ints = np.frombuffer(data, dtype=op2.idtype8).reshape(nelements, 4).copy()
                    eids = ints[:, 0] // 10
                    nids = ints[:, 2]
                    assert eids.min() > 0, eids.min()
                    assert nids.min() >= 0, nids.min()
                    obj.element_node[ielement:ielement2, 0] = eids
                    obj.element_node[ielement:ielement2, 1] = nids

                #[eid, free_conv, cntl_node, free_conv_k]
                obj.data[obj.itime, ielement:ielement2, :] = floats[:, [1, 3]]
                obj.itotal = ielement2
                obj.ielement = ielement2
            else:
                fmt = mapfmt(op2._analysis_code_fmt + b'fif', op2.size)
                s1 = Struct(op2._endian + fmt)
                add_sort_x = getattr(obj, 'add_sort' + str(op2.sort_method))
                for unused_i in range(nelements):
                    edata = data[n:n+ntotal]
                    n += ntotal
                    out = s1.unpack(edata)
                    (eid_device, free_conv, cntl_node, free_conv_k) = out
                    eid, dt = get_eid_dt_from_eid_device(
                        eid_device, op2.nonlinear_factor, op2.sort_method)
                    assert cntl_node >= 0, cntl_node
                    add_sort_x(dt, eid, cntl_node, free_conv, free_conv_k)
        else:  # pragma: no cover
            msg = op2.code_information()
            return op2._not_implemented_or_skip(data, ndata, msg), None, None
        return n, nelements, ntotal

    def _print_obj_name_on_crash(func):
        """
        Decorator debugging function to print the object name and an needed parameters
        """
        def new_func(self, data):
            """
            The actual function exec'd by the decorated function.
            """
            try:
                n = func(self, data)
            except Exception:
                raise
                #print("----------")
                #try:
                    #print(op2.obj)
                #except Exception:
                    #print("error printing %r" % op2.obj.__class__.__name__)
                #print(op2.data_code)
                #if op2.obj is not None:
                    ##from pyNastran.utils import object_attributes
                    ##print object_attributes(op2.obj)
                    #print(op2.obj.data_code)
                #print("----------")
                #raise
            return n
        return new_func

    # @_print_obj_name_on_crash
    def _read_oef1_loads(self, data: bytes, ndata: int) -> int:
        """Reads the OEF1 table; stores the element forces/heat flux."""
        op2 = self.op2
        #self._apply_oef_ato_crm_psd_rms_no('') # TODO: just testing
        if op2._results.is_not_saved('force'):
            return ndata

        prefix, postfix = get_oef_prefix_postfix(op2)
        if prefix and op2._results.is_not_saved(prefix.strip('.')):
            # op2.log.info(f'skipping {op2.table_name} due to prefix={prefix}')
            return ndata

        _sort_method = func1(op2.tCode)
        result_type = op2.result_type # func7(op2.tCode)

        #print('prefix=%r postfix=%s element_name=%s' % (prefix, postfix, op2.element_name))
        unused_flag = 'element_id'
        (num_wide_real, num_wide_imag) = self._oef_force_code()
        if op2.is_debug_file:
            op2.binary_debug.write(f'  num_wide_real = {num_wide_real!r}\n')
            op2.binary_debug.write(f'  num_wide_imag = {num_wide_imag!r}\n')

        n = 0
        is_magnitude_phase = op2.is_magnitude_phase()
        dt = op2.nonlinear_factor

        element_type = op2.element_type
        if element_type in [1, 3, 10]:  # rods
            n, nelements, ntotal = oef_crod(self.op2, data, ndata, dt, is_magnitude_phase,
                                            result_type, prefix, postfix)

        elif element_type == 2:  # cbeam
            #2-CBEAM
            n, nelements, ntotal = oef_cbeam(self.op2, data, ndata, dt, is_magnitude_phase,
                                             result_type, prefix, postfix)

        elif element_type in [11, 12, 13, 14,   # springs
                                  20, 21, 22, 23]:  # dampers
            # 11-CELAS1
            # 12-CELAS2
            # 13-CELAS3
            # 14-CELAS4

            # 20-CDAMP1
            # 21-CDAMP2
            # 22-CDAMP3
            # 23-CDAMP4
            n, nelements, ntotal = oef_celas_cdamp(self.op2, data, ndata, dt, is_magnitude_phase,
                                                   result_type, prefix, postfix)
            # if op2.table_name_str.startswith('OEFCR'):
            #     raise RuntimeError(op2.table_name)

        elif element_type == 24:  # CVISC
            n, nelements, ntotal = oef_cvisc(self.op2, data, ndata, dt, is_magnitude_phase,
                                             result_type, prefix, postfix)

        elif element_type == 34:  # cbar-34
            n, nelements, ntotal = oef_cbar_34(self.op2, data, ndata, dt, is_magnitude_phase,
                                               result_type, prefix, postfix)

        elif element_type == 100:  # cbar-100
            n, nelements, ntotal = oef_cbar_100(self.op2, data, ndata, dt, is_magnitude_phase,
                                                result_type, prefix, postfix)

        elif element_type in [33, 74]:  # centroidal shells
            # 33-CQUAD4
            # 74-CTRIA3
            n, nelements, ntotal = oef_shells_centroidal(self.op2, data, ndata, dt, is_magnitude_phase,
                                                         result_type, prefix, postfix)
        elif op2.is_nx and element_type in [227, 228]:  # centroidal shells
            # 227-CTRIAR? (C:\MSC.Software\simcenter_nastran_2019.2\tpl_post1\cqrdbx102.op2)
            # 228-CQUADR
            n, nelements, ntotal = oef_shells_centroidal(self.op2, data, ndata, dt, is_magnitude_phase,
                                                         result_type, prefix, postfix)

        elif element_type in [64, 70, 75, 82, 144]: # bilinear shells
            # 64-CQUAD8
            # 70-CTRIAR
            # 75-CTRIA6
            # 82-CQUADR
            # 144-CQUAD4-bilinear
            n, nelements, ntotal = oef_shells_nodal(self.op2, data, ndata, dt, is_magnitude_phase,
                                                    result_type, prefix, postfix)

        elif element_type in [95, 96, 97, 98]: # composites
            # 95 - CQUAD4
            # 96 - CQUAD8
            # 97 - CTRIA3
            # 98 - CTRIA6 (composite)
            n, nelements, ntotal = oef_shells_composite(self.op2, data, ndata, dt, is_magnitude_phase,
                                                        result_type, prefix, postfix)
        elif op2.is_nx and element_type in [232, 233]: # composites
            # 232 - CQUADR
            # 233 - CTRIAR
            n, nelements, ntotal = oef_shells_composite(self.op2, data, ndata, dt, is_magnitude_phase,
                                                        result_type, prefix, postfix)

        elif element_type in [39, 67, 68]: # solids
            # 39-CTETRA
            # 67-CHEXA
            # 68-CPENTA
            if op2.read_mode == 1:
                return ndata
            #op2._results._found_result('solid_forces')
            raise RuntimeError(op2.code_information())
            #if op2.format_code == 1 and op2.num_wide == 0:  # real
                ##self.create_transient_object(result_name, self.solidForces, RealCSolidForce)
                #raise RuntimeError(op2.code_information())
            #else:
                #msg = op2.code_information()
                #return op2._not_implemented_or_skip(data, ndata, msg)

        elif element_type == 53:  # ctriax6
            # 53-CTRIAX6
            op2._results._found_result('ctriax_force')
            #if op2.format_code == 1 and op2.num_wide == 0:  # real
                #pass
                #self.create_transient_object(self.ctriax_force, RealCTriaxForce)  # undefined
            #else:  # pragma: no cover
            raise NotImplementedError(op2.code_information())
                #msg = op2.code_information()
                #return op2._not_implemented_or_skip(data, ndata, msg)
            #return ndata

        elif element_type == 4:  # cshear
            n, nelements, ntotal = oef_cshear(self.op2, data, ndata, dt, is_magnitude_phase,
                                              result_type, prefix, postfix)

        elif element_type == 35:  # coneax
            n, nelements, ntotal = oef_cconeax(self.op2, data, ndata, dt, is_magnitude_phase,
                                               result_type, prefix, postfix)

        elif element_type == 38:  # cgap
            n, nelements, ntotal = oef_cgap(self.op2, data, ndata, dt, is_magnitude_phase,
                                            result_type, prefix, postfix)

        elif element_type == 69:  # cbend
            n, nelements, ntotal = oef_cbend(self.op2, data, ndata, dt, is_magnitude_phase,
                                             result_type, prefix, postfix)

        elif element_type in [76, 77, 78, 79]:
            # 76-HEXPR
            # 77-PENPR
            # 78-TETPR
            # 79-CPYRAM
            n, nelements, ntotal = oef_csolid_pressure(self.op2, data, ndata, dt, is_magnitude_phase,
                                                       result_type, prefix, postfix)

        elif element_type in [102, 280]:
            # 102: cbush
            # 280: cbear
            n, nelements, ntotal = oef_cbush(self.op2, data, ndata, dt, is_magnitude_phase,
                                             result_type, prefix, postfix)

        elif element_type == 126 and op2.is_msc:
            # 119-CFAST-MSC
            n, nelements, ntotal = oef_cbush(self.op2, data, ndata, dt, is_magnitude_phase,
                                             result_type, prefix, postfix)
        elif element_type == 119 and op2.is_nx:
            # 119-CFAST-NX
            n, nelements, ntotal = oef_cbar_34(self.op2, data, ndata, dt, is_magnitude_phase,
                                               result_type, prefix, postfix)
        elif element_type in [117, 200]:
            # 117-CWELDC
            # 200-CWELD
            n, nelements, ntotal = oef_cbar_34(self.op2, data, ndata, dt, is_magnitude_phase,
                                               result_type, prefix, postfix)
        #elif element_type == 119 and op2.is_msc:
            #raise NotImplementedError(op2.code_information())
        elif element_type == 235:
            # 235-CQUADR
            return op2._not_implemented_or_skip(data, ndata, op2.code_information())

        elif element_type in [145, 146, 147,
                              189, 190, 191]:
            # removed from msc/nx
            # 145-VUHEXA
            # 146-VUPENTA
            # 147-VUTETRA
            # 189-VUQUAD
            # 190-VTRIA
            # 191-VUBEAM
            # n, nelements, ntotal = self._oef_vu_shell(data, ndata, dt, is_magnitude_phase,
            #                                           result_type, prefix, postfix)
            # n, nelements, ntotal = self._oef_vu_beam(data, ndata, dt, is_magnitude_phase,
            #                                          result_type, prefix, postfix)
            if op2.read_mode == 1:
                msg = f'{op2.table_name_str} {op2.element_name}-{element_type} has been removed'
                op2.log.warning(msg)
            return ndata
            # return op2._not_implemented_or_skip(data, ndata, msg)
        elif op2.is_nx:
            if element_type in [118, 343, 344, 345, 346, 347, 348, 349, 350, 351, 352,
                                356, 357, 363]:
                # 118 CWELDP
                # 343 CTRIA6 SOL 401
                #344 CQUAD8 SOL 401
                #345 CTRIAR SOL 401
                # 347 CBAR SOL 401
                # 348 CBEAM SOL 401
                #346 CQUADR SOL 401
                # 349 CBUSH1D SOL 401
                # 350 CELAS1 SOL 401
                # 351 CELAS2 SOL 401
                # 352 CBUSH SOL 401
                #356 Composite quadrilateral shell element (CQUAD8); SOL 402?
                #357 Composite triangular shell element (CTRIAR); SOL 402?
                # 363 CROD SOL 402
                return op2._not_implemented_or_skip(data, ndata, op2.code_information())
            #print(op2.code_information())
            #nx_missing
            return op2._not_implemented_or_skip(data, ndata, op2.code_information())
        else:
            if element_type == 118:
                # 118-WELDP
                n, nelements, ntotal = oef_cbar_34(
                    self.op2, data, ndata, dt, is_magnitude_phase,
                    result_type, prefix, postfix)
            elif element_type in [184]:
                # 184-CBEAM3
                if op2.read_mode == 1:
                    msg = f'{op2.table_name_str} {op2.element_name}-{element_type} is not supported'
                    op2.log.warning(msg)
                return ndata
            else:
                #print(op2.code_information())
                #msc_missing
                return op2._not_implemented_or_skip(data, ndata, op2.code_information())

        if nelements is None:
            return n
        if element_type != 2:
            # CBEAM-2: has a finalize step
            op2._op2_readers.reader_oes.check_element_ids()

        #assert op2.thermal == 0, op2.thermal
        assert ndata > 0, ndata
        assert nelements > 0, 'nelements=%r element_type=%s element_name=%r num_wide=%s' % (
            nelements, op2.element_type, op2.element_name, op2.num_wide)
        #assert ndata % ntotal == 0, '%s n=%s nwide=%s len=%s ntotal=%s' % (op2.element_name, ndata % ntotal, ndata % op2.num_wide, ndata, ntotal)
        assert op2.num_wide * 4 * op2.factor == ntotal, f'numwide*4={op2.num_wide*4} ntotal={ntotal}'
        assert n is not None and n > 0, op2.code_information()
        return n


def shock_response_prefix(thermal: int) -> str:
    prefix = ''
    if thermal == 0:
        pass
    elif thermal == 2:
        prefix = 'abs.'  # Scaled response spectra ABS
    elif thermal == 4:
        #D:\NASA\git\examples\move_tpl\ms103.op2
        #C:\NASA\m4\formats\git\examples\move_tpl\ms103.op2
        prefix = 'srss.' #   # Scaled response spectra SRSS
    elif thermal == 8:
        prefix = 'nrl.'  # Scaled response spectra NRL
    else:
        assert thermal in [2, 4, 8], thermal # , op2.code_information() # abs
    return prefix


def get_oef_prefix_postfix(op2: OP2) -> tuple[str, str]:
    """
    NX Case Control  Block         Description
    ===============  ==========    ===========
    NLSTRESS         OESNLXR       Nonlinear static stresses
    BOUTPUT          OESNLBR       Slideline stresses
    STRESS           OESNLXD       Nonlinear Transient Stresses
    STRESS           OES1C/OSTR1C  Ply stresses/strains
    STRESS           OES1X         Element stresses with intermediate (CBAR and CBEAM)
                                   station stresses and stresses on nonlinear elements
    STRESS           OES/OESVM     Element stresses (linear elements only)
    STRAIN           OSTR1         Element strains
    STRESS/STRAIN    DOES1/DOSTR1  Scaled Response Spectra
    MODCON           OSTRMC        Modal contributions
    """
    prefix = ''
    postfix = ''
    table_name_bytes = op2.table_name
    assert isinstance(table_name_bytes, bytes), table_name_bytes
    is_sort1 = table_name_bytes in SORT1_TABLES_BYTES
    assert table_name_bytes in TABLES_BYTES, table_name_bytes

    if table_name_bytes in [b'OEF1X', b'OEF1', b'OEF2']:
        if op2.thermal == 0:
            prefix = 'force.'
        elif op2.thermal == 1:
            prefix = 'thermal_load.'
        else:
            raise NotImplementedError(op2.code_information())
    elif table_name_bytes in [b'HOEF1']:
        postfix = '_flux'
    #elif op2.table_name in ['OESNLXR']:
        #prefix = 'sideline_'
    #elif op2.table_name in ['OESNLXD', 'OESNL1X', 'OESNLXR']:
        #prefix = 'nonlinear_'
    #elif op2.table_name == 'OESNLBR':
        #prefix = 'sideline_'
    #elif op2.table_name == 'OESRT':
        #prefix = 'strength_ratio.'
    #elif op2.table_name in ['OESCP', 'OESTRCP']:
        #pass # TODO: update
    elif table_name_bytes in [b'OEFCRM1', b'OEFCRM2']:
        assert op2.table_code in [4, 504], op2.code_information()
        prefix = 'crm.'
        op2._op2_readers.reader_oes._set_as_random()
    elif table_name_bytes in [b'OEFPSD1', b'OEFPSD2']:
        assert op2.table_code in [4, 604], op2.code_information()
        op2._op2_readers.reader_oes._set_as_random()
        prefix = 'psd.'
    elif table_name_bytes in [b'OEFRMS1', b'OEFRMS2', b'OEFPK1']:
        assert op2.table_code in [4, 404, 804], op2.code_information()
        op2._op2_readers.reader_oes._set_as_random()
        is_sort1 = True
        op2._analysis_code_fmt = b'i'
        prefix = 'rms.'
    elif table_name_bytes in [b'OEFNO1', b'OEFNO2']:
        assert op2.table_code in [4, 904], op2.code_information()
        op2._op2_readers.reader_oes._set_as_random()
        op2.sort_method = 1
        op2.data_code['nonlinear_factor'] = None
        op2._analysis_code_fmt = b'i'
        assert op2.sort_method == 1, op2.code_information()
        prefix = 'no.'
    elif table_name_bytes in [b'OEFATO1', b'OEFATO2']:
        assert op2.table_code in [4], op2.code_information()
        prefix = 'ato.'

    elif table_name_bytes in [b'RAFCONS']:
        prefix = 'RAFCONS.'
    elif table_name_bytes in [b'RAFEATC']:
        prefix = 'RAFEATC.'
    elif table_name_bytes in [b'DOEF1']:
        assert op2.table_code in [4], op2.code_information()
        prefix = shock_response_prefix(op2.thermal)
    elif table_name_bytes in [b'OEFIT']:
        assert op2.table_code in [25], op2.code_information()
        prefix = 'failure_indices.'
        #raise NotImplementedError(op2.code_information())
    elif table_name_bytes in [b'OEFITSTN']: # composite failure indicies
        assert op2.table_code in [25], op2.code_information()
        prefix = 'failure_indices.'
    else:
        raise NotImplementedError('%r' % op2.table_name)

    op2.sort_bits.is_sort1 = is_sort1 # sort2
    op2.sort_method = 1 if is_sort1 else 2
    return prefix, postfix
