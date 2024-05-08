"""
defines:
  - read_descyc(op2_reader)
  - read_dbcopt(op2_reader)
  - read_destab(op2_reader)
  - read_dscmcol(op2_reader)
  - read_r1tabrg(op2_reader)
  - read_hisadd(op2_reader)

"""
from __future__ import annotations
from struct import unpack, Struct # , error as struct_error
from typing import TYPE_CHECKING

import numpy as np

#from cpylog import SimpleLogger
#from pyNastran.utils.numpy_utils import integer_types
#from pyNastran.f06.errors import FatalError
#from pyNastran.op2.errors import FortranMarkerError, SortCodeError, EmptyRecordError

from pyNastran.op2.op2_interface.utils import (
    mapfmt, reshape_bytes_block,
    #reshape_bytes_block_size,
)

from pyNastran.op2.result_objects.design_response import (
    WeightResponse, DisplacementResponse, StressResponse, StrainResponse, ForceResponse,
    FlutterResponse, FractionalMassResponse, Convergence, Desvars, DSCMCOL)
if TYPE_CHECKING:
    from pyNastran.op2.op2_interface.op2_reader import OP2Reader
    from pyNastran.op2.op2 import OP2


def read_dbcopt(op2_reader: OP2Reader) -> None:
    """reads the DBCOPT table, which is a design variable history table"""
    #C:\MSC.Software\simcenter_nastran_2019.2\tpl_post1\cc577.op2
    # ints    = (101, 11, 10, 3, 4, 0, 0)
    # objective_function = [1175749.5, 711181.875, 369194.03125, 112453.1406, 229.4625, 50.9286, 50.8863, 50.8316, 50.7017, 49.7571, 49.3475]
    # approx = [0.0, 651307.3125, 388093.0625, 150222.453125, 5894.2626, 50.9092, 50.8834, 50.8252, 49.4544, 49.6455, 49.2533]
    # max_value_of_constraint = [nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan]
    # desvar_ids = [1, 2, 3]
    # cycle_1_values = [1.0, 1.0, 1.0]
    # cycle_n_values = [1.4, 0.6, 1.4,
    #                   1.96, 0.36, 1.96,
    #                   2.744, 0.216, 2.744,
    #                   3.8416, 0.1296, 3.8416, ...]

    # 1 NFEA I Number of finite element analyses
    # 2 NAOP I Number of optimization cycles w.r.t. approximate model
    # 3 NDV I Number of design variables
    # 4 NCC I Convergence criterion
    # 5 UNDEF(2 ) None
    op2: OP2 = op2_reader.op2
    endian = op2_reader._endian
    size = op2_reader.size

    op2.table_name = op2_reader._read_table_name(rewind=False)
    op2_reader.read_markers([-1])
    data = op2_reader._read_record()

    fmt = mapfmt(endian + b'7i', size)
    num, nopt, napprox, nvars, one, zeroa, zerob = Struct(fmt).unpack(data)
    assert num == 101, num
    assert zeroa == 0, zeroa
    assert zerob == 0, zerob

    # (101, 11, 10, 3, 4, 0, 0)
    #op2_reader.show_data(data)
    op2_reader.read_3_markers([-2, 1, 0])
    data = op2_reader._read_record()

    if size == 4:
        name, = Struct(endian + b'8s').unpack(data)
    else:
        name, = Struct(endian + b'16s').unpack(data)
        name = reshape_bytes_block(name)
    assert name == b'DBCOPT  ', name

    op2_reader.read_3_markers([-3, 1, 0])
    data = op2_reader._read_record()
    #ndata = len(data) // 4
    if size == 4:
        fdtype = 'float32'
        idtype = 'int32'
    else:
        fdtype = 'float64'
        idtype = 'int64'
    objective_function = np.frombuffer(data, dtype=fdtype) # .tolist()
    #print(f'  objective_function = {objective_function}; n={len(objective_function)}')
    assert len(objective_function) == nopt, f'len(objective_function)={len(objective_function)} nopt={nopt}'

    op2_reader.read_3_markers([-4, 1, 0])
    data = op2_reader._read_record()
    approx = np.frombuffer(data, dtype=fdtype).copy()# .tolist()
    napprox_actual = len(approx)
    if approx[0] == 0.0:
        approx[0] = np.nan
        napprox_actual -= 1
    #print(approx.tolist())
    #assert napprox_actual == napprox, f'napprox_actual={napprox_actual} napprox={napprox}'
    #print(f'  approx = {approx}; n={len(approx)}')

    op2_reader.read_3_markers([-5, 1, 0])
    data = op2_reader._read_record()
    max_value_of_constraint = np.frombuffer(data, dtype=fdtype).tolist()
    #print(f'  max_value_of_constraint = {max_va/lue_of_constraint}; n={len(max_value_of_constraint)}')


    op2_reader.read_3_markers([-6, 1, 0])
    data = op2_reader._read_record()
    desvar_ids = np.frombuffer(data, dtype=idtype).tolist()
    assert len(desvar_ids) == nvars, f'len(desvars)={len(desvars)} nvars={nvars}'

    op2_reader.read_3_markers([-7, 1, 0])
    data = op2_reader._read_record()
    cycle_1_values = np.frombuffer(data, dtype=fdtype).tolist()

    op2_reader.read_3_markers([-8, 1, 0])
    marker0 = op2_reader.get_marker1(rewind=True)
    if marker0 == 0:
        op2_reader.read_markers([0])
        return
    data = op2_reader._read_record()
    cycle_n_values = np.frombuffer(data, dtype=fdtype)
    cycle_n_values2 = cycle_n_values.reshape(len(cycle_n_values) // nvars, nvars)
    cycle_n_values3 = np.vstack([cycle_1_values, cycle_n_values2])

    approx_obj_constraint = np.vstack([approx, objective_function, max_value_of_constraint]).T
    #print(f'  desvar_ids = {desvar_ids}')
    #print(f'  approx_obj_constraint; {approx_obj_constraint.shape}:\n{approx_obj_constraint}')
    #print(f'  cycle_n_values {cycle_n_values3.shape}:\n{cycle_n_values3}')

    op2_reader.read_markers([-9, 1, 0, 0])
    #data = op2_reader._read_record()
    #op2_reader.show_data(data)
    #op2_reader.show_ndata(100)

def read_descyc(op2_reader: OP2Reader) -> None:
    """reads the DESCYC table"""
    op2: OP2 = op2_reader.op2
    endian = op2_reader._endian
    size = op2_reader.size

    # TODO: I think this is used to handle optimization
    op2._count += 1

    #op2.log.debug("table_name = %r" % op2.table_name)
    op2.table_name = op2_reader._read_table_name(rewind=False)
    op2_reader.read_markers([-1])
    data = op2_reader._read_record()
    fmt = mapfmt(op2_reader._endian + b'7i', size)
    unused_ints = Struct(fmt).unpack(data)
    op2_reader.read_3_markers([-2, 1, 0])
    data = op2_reader._read_record()
    if size == 4:
        name, = Struct(op2_reader._endian + b'8s').unpack(data)
    else:
        name, = Struct(endian + b'16s').unpack(data)
        name = reshape_bytes_block(name)
    assert name == b'DESCYC  ', name

    op2_reader.read_3_markers([-3, 1, 0])
    data = op2_reader._read_record()
    if size == 4:
        design_cycle, design_cycle_type_bytes = Struct(endian + b'i8s').unpack(data)
    else:
        design_cycle, design_cycle_type_bytes = Struct(endian + b'q16s').unpack(data)
        design_cycle_type_bytes = reshape_bytes_block(design_cycle_type_bytes)

    #Design cycle type; 'D' for discretized design cycle
    #Blank for continuous design cycle
    if design_cycle_type_bytes == b'        ':
        design_cycle_type_str = 'continuous'
    elif design_cycle_type_bytes == b'D       ':
        design_cycle_type_str = 'discrete'
    else:  # pragma: no cover
        raise NotImplementedError(design_cycle_type_bytes)

    op2_reader.read_markers([-4, 1, 0, 0])
    #print('descyc')
    #print('  ints =', ints)
    #print('  name =', name)
    #print(f'  design_cycle={design_cycle} type={design_cycle_type_str!r} count={op2._count}')

def read_destab(self) -> None:
    """reads the DESTAB table"""
    op2: OP2 = self.op2

    result_name = 'responses.desvars'
    if op2._results.is_not_saved(result_name):
        data = self._skip_table('DESTAB', warn=False)
        return

    #if self.read_mode == 1:
        #return ndata
    op2.table_name = self._read_table_name(rewind=False)
    #self.log.debug('table_name = %r' % op2.table_name)
    if self.is_debug_file:
        self.binary_debug.write('_read_destab - %s\n' % op2.table_name)

    self.read_markers([-1])
    data = self._read_record()

    # (101, 3, 3, 0, 3, 0, 0)

    itable = -2
    markers = self.read_3_markers([itable, 1, 0])
    data = self._read_record()
    if self.size == 4:
        destab = op2.struct_8s.unpack(data)[0].rstrip()
        structi = Struct('2i 8s 4f')
    else:
        destab = op2.struct_16s.unpack(data)[0]
        destab = reshape_bytes_block(destab).rstrip()
        structi = Struct('2q 16s 4d')
    assert destab == b'DESTAB', destab

    itable -= 1
    markers = self.read_3_markers([itable, 1, 0])

    desvars = []
    while 1:
        markers = self.get_nmarkers(1, rewind=True)
        if markers == [0]:
            break

        data = self._read_record()

        #self.show(100, types='ifs', endian=None)
        #self.show_data(data[:8], types='ifs', endian=None)

        #1 IDVID I Internal design variable identification number
        #2 DVID I External design variable identification number
        #3 LABEL1 CHAR4 First part of design Variable
        #4 LABEL2 CHAR4 Second part of design Variable
        #5 VMIN RS Lower bound
        #6 VMAX RS Upper bound
        #7 DELX RS Move limit for a design cycle
        # 8 ???

        #C:\NASA\m4\formats\git\examples\move_tpl\accopt3.op2
        #---------------------------------------------------------------------------------------------------------
            #INTERNAL       DESVAR                         LOWER                               UPPER
                #ID            ID          LABEL            BOUND             VALUE             BOUND
        #---------------------------------------------------------------------------------------------------------
            #1             1      THICK           2.0000E-02        5.0000E-02        8.0000E-02
            #2             2      SPRING          1.0000E-02        6.2500E-02        7.5000E-02
            #3             3      SPRING          5.0000E-02        1.2500E-01        1.5000E-01
        # internal, desvar, label, lower,   upper,  ???,  ???
        # (1, 1, b'THICK   ', 0.019999, 0.07999999, -0.5, 0.0)
        # (2, 2, b'SPRING  ', 0.009999, 0.07500000, -0.5, 0.0)
        # (3, 3, b'SPRING  ', 0.050000, 0.15000000, -0.5, 0.0)

        ##        id       label   xinit   xlb   xub delxv
        #DESVAR  1       THICK   .05     0.02    .08
        #DESVAR  2       SPRING  .0625   .01     .075
        #DESVAR  3       SPRING  .125    .05     .15

        # C:\NASA\m4\formats\git\examples\move_tpl\betaadj.op2
        #---------------------------------------------------------------------------------------------------------
            #INTERNAL       DESVAR                         LOWER                               UPPER
                #ID            ID          LABEL            BOUND             VALUE             BOUND
        #---------------------------------------------------------------------------------------------------------
            #11            20      BETA            1.0000E-03        8.0000E-01        1.0000E+20
        # int, des, label,     lower,     upper,     ???,  ???
        #(11, 20, b'BETA    ', 0.0010000, 1.000e+20, 0.200, 0.0)
        #        id       label   xinit   xlb   xub delxv
        #desvar  20       beta    0.8     0.001	xub	0.20
        desvar = structi.unpack(data)
        #internal_id, desvar_id, label, lower, upper, delxv, dunno = desvar
        #print(desvar)
        #assert np.allclose(desvar[5], -0.5), desvar  # -0.5 is the default
        assert np.allclose(desvar[6], 0.0), desvar
        desvars.append(desvar)
        itable -= 1
        markers = self.read_3_markers([itable, 1, 0])

    self.op2.op2_results.responses.desvars = Desvars(desvars)
    #if self.read_mode == 2:
        #self.log.warning('DESTAB results were read, but not saved')
    markers = self.read_markers([0])


def read_dscmcol(op2_reader: OP2Reader) -> None:
    """reads the DSCMCOL table, which defines the columns? for the DSCM2 table"""
    op2: OP2 = op2_reader.op2
    log = op2.log
    endian = op2_reader._endian
    size = op2_reader.size

    log.debug("table_name = %r" % op2.table_name)
    op2.table_name = op2_reader._read_table_name(rewind=False)
    op2_reader.read_markers([-1])
    data = op2_reader._read_record()
    #fmt = mapfmt(endian + b'7i', size)
    #num, ndesvars, one_zero, zeroa, zerob, zeroc, zerod = Struct(fmt).unpack(data)
    #print(num, ndesvars, one_zero, zeroa, zerob, zeroc, zerod)
    # (101, 3, 1, 0, 0, 0, 0)
    #op2_reader.show_data(data)
    op2_reader.read_3_markers([-2, 1, 0])
    data = op2_reader._read_record()
    if size == 4:
        name, = Struct(endian + b'8s').unpack(data)
    else:
        name, = Struct(endian + b'16s').unpack(data)
        name = reshape_bytes_block(name)
    assert name == b'DSCMCOL ', name

    op2_reader.read_3_markers([-3, 1, 0])
    data = op2_reader._read_record()

    responses = {}
    if op2.read_mode == 2:
        ints = np.frombuffer(data, dtype=op2.idtype8)
        floats = np.frombuffer(data, dtype=op2.fdtype8)
        nresponses_dresp1 = len(ints) // 9
        dscmcol_dresp1(responses, nresponses_dresp1, ints, floats)

    #op2_reader.show_data(data[4*idata:])
    op2_reader.read_3_markers([-4, 1, 0])
    nfields = op2_reader.get_marker1(rewind=True)
    if nfields == 0:
        _save_dscmcol_response(op2_reader, responses)
        op2_reader.read_markers([0])
        return

    data = op2_reader._read_record()
    if op2.read_mode == 2:
        # read the DRESP2 columns
        ints = np.frombuffer(data, dtype=op2.idtype8)
        floats = np.frombuffer(data, dtype=op2.fdtype8)
        nresponses_dresp2 = len(ints) // 6
        dscmcol_dresp2(responses, nresponses_dresp2, ints, floats)

    _save_dscmcol_response(op2_reader, responses)
    op2_reader.read_markers([-5, 1, 0, 0])

def _save_dscmcol_response(op2_reader: OP2Reader, responses):
    """saves the DSCMCOL dictionary"""
    op2: OP2 = op2_reader.op2
    if op2.read_mode == 2:
        assert len(responses) > 0
    if responses:
        if op2.op2_results.responses.dscmcol is not None:
            op2_reader.log.warning('overwriting DSCMCOL')
        respi = DSCMCOL(responses)
        str(respi)
        op2.op2_results.responses.dscmcol = respi

def dscmcol_dresp1(responses: dict[int, dict[str, Any]],
                   nresponses_dresp1: int,
                   ints: np.ndarray,
                   floats: np.ndarray) -> None:
    """helper for DSCMCOL"""
    idata = 0
    if nresponses_dresp1 == 0:
        return

    for iresp in range(nresponses_dresp1):
        # Record – Type 1 Responses
        #   1 IRID  I Internal response identification number
        #   2 RID   I External response identification number
        #   3 RTYPE I Response Type
        # RTYPE=8 Static force
        #   4 EID I Element identification number
        #   5 COMP I Force component number
        #   6 SUBCASE I Subcase identification number
        #   7 UNDEF None
        #   8 VIEWID I View element identification number
        #   9 SEID I Superelement identification number
        # RTYPE=10 Composite stress
        #   4 EID I Element identification number
        #   5 COMP I Stress component number
        #   6 SUBCASE I Subcase identification number
        #   7 UNDEF None
        #   8 PLY I Ply number
        #   9 SEID I Superelement identification number
        # RTYPE=13 Static SPC force
        #   4 GRID I Grid identification number
        #   5 COMP I SPC force component number
        #   6 SUBCASE I Subcase identification number
        #   7 UNDEF(2) None
        #   9 SEID I Superelement identification number
        # RTYPE=14 Element static strain energy
        #   4 EID I Element identification number
        #   5 COMP I Strain energy component number
        #   6 SUBCASE I Subcase identification number
        #   7 UNDEF(2) None
        #   9 SEID I Superelement identification number
        # RTYPE=17 Compliance
        #   4 UNDEF(2) None
        #   6 SUBCASE I Subcase identification number
        #   7 UNDEF(3) None
        # RTYPE=21 Frequency response velocity
        #   4 GRID I Grid identification number
        #   5 COMP I Velocity component number
        #   6 SUBCASE I Subcase identification number
        #   7 FREQ RS Frequency
        #   8 UNDEF None
        #   9 SEID I Superelement identification number
        # RTYPE=22 Frequency response acceleration
        #   4 GRID I Grid identification number
        #   5 COMP I Acceleration component number
        #   6 SUBCASE I Subcase identification number
        #   7 FREQ RS Frequency
        #   8 UNDEF None
        #   9 SEID I Superelement identification number
        # RTYPE=23 Frequency response SPC Force
        #   4 GRID I Grid identification number
        #   5 COMP I SPC Force component number
        #   6 SUBCASE I Subcase identification number
        #   7 FREQ RS Frequency
        #   8 UNDEF None
        #   9 SEID I Superelement identification number
        # RTYPE=24 Frequency response stress
        #   4 EID I Element identification number
        #   5 COMP I Stress component number
        #   6 SUBCASE I Subcase identification number
        #   7 FREQ RS Frequency
        #   8 UNDEF None
        #   9 SEID I Superelement identification number
        # RTYPE=25 Frequency response force
        #   4 EID I Element identification number
        #   5 COMP I Force component number
        #   6 SUBCASE I Subcase identification number
        #   7 FREQ RS Frequency
        #   8 UNDEF None
        #   9 SEID I Superelement identification number
        # RTYPE=26 RMS displacement
        #   4 GRID I Grid identification number
        #   5 COMP I RMS displacement component number
        #   6 SUBCASE I Subcase identification number
        #   7 UNDEF None
        #   8 RANDPS I RANDPS ID
        #   9 SEID I Superelement identification number
        # RTYPE=27 RMS velocity
        #   4 GRID I Grid identification number
        #   5 COMP I RMS velocity component number
        #   6 SUBCASE I Subcase identification number
        #   7 UNDEF None
        #   8 RANDPS I RANDPS ID
        #   9 SEID I Superelement identification number
        # RTYPE=28 RMS acceleration
        #   4 GRID I Grid identification number
        #   5 COMP I RMS acceleration component number
        #   6 SUBCASE I Subcase identification number
        #   7 UNDEF None
        #   8 RANDPS I RANDPS ID
        #   9 SEID I Superelement identification number
        # RTYPE=30 PSD velocity
        # 4 GRID I Grid identification number
        # 5 COMP I PSD velocity component number
        # 6 SUBCASE I Subcase identification number
        # 7 FREQ RS Frequency
        # 8 RANDPS I RANDPS ID
        # 9 SEID I Superelement identification number
        # RTYPE=60 Transient response displacement
        # 4 GRID I Grid identification number
        # 5 COMP I Displacement component number
        # 6 SUBCASE I Subcase identification number
        # 7 TIME RS Time
        # 8 UNDEF None
        # 9 SEID I Superelement identification number
        # RTYPE=61 Transient response velocity
        # 4 GRID I Grid identification number
        # 5 COMP I Velocity component number
        # 6 SUBCASE I Subcase identification number
        # 7 TIME RS Time
        # 8 UNDEF None
        # 9 SEID I Superelement identification number
        # RTYPE=62 Transient response acceleration
        # 4 GRID I Grid identification number
        # 5 COMP I Acceleration component number
        # 6 SUBCASE I Subcase identification number
        # 7 TIME RS Time
        # 8 UNDEF None
        # 9 SEID I Superelement identification number
        # RTYPE=63 Transient response SPC Force
        # 4 GRID I Grid identification number
        # 5 COMP I SPC force component number
        # 6 SUBCASE I Subcase identification number
        # 7 TIME RS Time
        # 8 UNDEF None
        # 9 SEID I Superelement identification number
        # RTYPE=64 Transient response stress
        # 4 EID I Element identification number
        # 5 COMP I Stress component number
        # 6 SUBCASE I Subcase identification number
        # 7 TIME RS Time
        # 8 UNDEF None
        # 9 SEID I Superelement identification number
        # RTYPE=65 Transient response force
        #   4 EID I Element identification number
        #   5 COMP I Force component number
        #   6 SUBCASE I Subcase identification number
        #   7 TIME RS Time
        #   8 UNDEF None
        #   9 SEID I Superelement identification number
        # RTYPE=81 Aeroelastic divergence
        #   4 SUBCASE I Subcase identification number
        #   5 UNDEF None
        #   6 ROOT I Root
        #   7 MACH RS Mach number
        #   8 UNDEF None
        #   9 SEID I Superelement identification number
        # RTYPE=82 Aeroelastic trim
        #   4 SUBCASE I Subcase identification number
        #   5 UNDEF None
        #   6 XID I XID
        #   7 UNDEF(2) None
        #   9 SEID I Superelement identification number
        # RTYPE=83 Aeroelastic stability derivative
        #   4 SUBCASE I Subcase identification number
        #   5 RU I R/U
        #   6 COMP I Component number
        #   7 UNDEF None
        #   8 XID I XID
        #   9 SEID I Superelement identification number
        # RTYPE=84 Aeroelastic flutter damping
        #   4 SUBCASE I Subcase identification number
        #   5 MODE I Mode number
        #   6 DENSITY RS Density
        #   7 MACH RS Mach number
        #   8 VEL RS Velocity
        #   9 SEID I Superelement identification number
        # End RTYPE

        internal_response_id = ints[idata]
        external_response_id = ints[idata+1]
        response_type = ints[idata+2]
        #print(f'internal_response_id={internal_response_id} '
              #f'external_response_id={external_response_id} response_type={response_type}')

        if response_type == 1:
            # RTYPE=1 Weight
            #   4 UNDEF(5) None
            #   9 SEID I Superelement identification number
            seid = ints[idata+8]
            response = {'name': 'weight', 'seid': seid}
            #print(f'  seid={seid} (weight)')
        elif response_type == 2:
            # RTYPE=2 Volume
            #   4 UNDEF(5) None
            #   9 SEID I Superelement identification number
            seid = ints[idata+8]
            response = {'name': 'volume', 'seid': seid}
            #print(f'  seid={seid} (volume)')
        elif response_type == 3:
            # RTYPE=3 Buckling
            #   4 MODE    I Mode number
            #   5 UNDEF None
            #   6 SUBCASE I Subcase identification number
            #   7 UNDEF(2) None
            #   9 SEID    I Superelement identification number
            mode = ints[idata+3]
            subcase = ints[idata+5]
            seid = ints[idata+8]
            response = {'name': 'buckling', 'mode_num': mode, 'subcase': subcase, 'seid': seid}
            #print(f'  mode={mode} subcase={subcase} seid={seid} (buckling)')

        elif response_type == 4:
            # RTYPE=4 Normal modes
            #   4 MODE    I Mode number
            #   5 UNDEF None
            #   6 SUBCASE I Subcase identification number
            #   7 UNDEF(2) None
            #   9 SEID    I Superelement identification number
            mode_num = ints[idata+3]
            # blank
            subcase = ints[idata+5]
            seid = ints[idata+8]
            response = {'name': 'normal modes', 'mode_num': mode_num, 'subcase': subcase, 'seid': seid}
            #print(f'  mode_num={mode_num} subcase={subcase} seid={seid} (normal modes)')
        elif response_type == 5:
            # RTYPE=5 Static displacement
            #   4 GRID    I Grid identification number
            #   5 COMP    I Displacement component number
            #   6 SUBCASE I Subcase identification number
            #   7 UNDEF(2) None
            #   9 SEID    I Superelement identification number
            grid = ints[idata+3]
            comp = ints[idata+4]
            subcase = ints[idata+5]
            seid = ints[idata+8]
            response = {'name': 'static displacement', 'grid': grid, 'component': comp,
                        'subcase': subcase, 'seid': seid}
            #print(f'  grid={grid} comp={comp} subcase={subcase} seid={seid} (static displacement)')
        elif response_type == 6:
            # RTYPE=6 Static stress
            # 4 EID     I Element identification number
            # 5 COMP    I Stress component number
            # 6 SUBCASE I Subcase identification number
            # 7 UNDEF(2) None
            # 9 SEID    I Superelement identification number
            eid = ints[idata+3]
            comp = ints[idata+4]
            subcase = ints[idata+5]
            seid = ints[idata+8]
            response = {'name': 'static stress', 'eid': eid, 'component': comp,
                        'subcase': subcase, 'seid': seid}
            #print(f'  eid={eid} comp={comp} subcase={subcase} seid={seid} (static stress)')
        elif response_type == 7:
            # RTYPE=7 Static strain
            # 4 EID     I Element identification number
            # 5 COMP    I Strain component number
            # 6 SUBCASE I Subcase identification number
            # 7 UNDEF None
            # 8 VIEWID  I View element identification number
            # 9 SEID    I Superelement identification number
            eid = ints[idata+3]
            comp = ints[idata+4]
            subcase = ints[idata+5]
            view_id = ints[idata+7]
            seid = ints[idata+8]
            response = {'name': 'static strain', 'eid': eid, 'component': comp,
                        'subcase': subcase, 'view id': view_id, 'seid': seid}
            #print(f'  eid={eid} comp={comp} subcase={subcase} view_id={view_id} seid={seid} (static strain)')
        elif response_type == 9:
            # RTYPE=9 Composite failure
            #   4 EID     I Element identification number
            #   5 COMP    I Failure component number
            #   6 SUBCASE I Subcase identification number
            #   7 UNDEF None
            #   8 PLY     I Ply number
            #   9 SEID    I Superelement identification number
            eid = ints[idata+3]
            comp = ints[idata+4]
            subcase = ints[idata+5]
            ply = ints[idata+7]
            seid = ints[idata+8]
            response = {'name': 'composite failure', 'eid': eid, 'component': comp,
                        'subcase': subcase, 'ply': ply, 'seid': seid}
            #print(f'  eid={eid} comp={comp} subcase={subcase} ply={ply} seid={seid} (composite failure)')
        elif response_type == 11:
            # RTYPE=11 Composite strain
            #   4 EID     I Element identification number
            #   5 COMP    I Strain component number
            #   6 SUBCASE I Subcase identification number
            #   7 UNDEF None
            #   8 PLY     I Ply number
            #   9 SEID    I Superelement identification number
            eid = ints[idata+3]
            comp = ints[idata+4]
            subcase = ints[idata+5]
            ply = ints[idata+7]
            seid = ints[idata+8]
            response = {'name': 'composite strain', 'eid': eid, 'component': comp,
                        'subcase': subcase, 'ply': ply, 'seid': seid}
            #print(f'  eid={eid} comp={comp} subcase={subcase} ply={ply} seid={seid} (composite strain)')
        elif response_type == 15:
            # CEIG - complex eigenvalues
            mode_num = ints[idata+3]
            subcase = ints[idata+5]
            seid = ints[idata+8]
            #print(ints[idata+4:idata+9])
            #print(floats[idata+4:idata+9])
            #print(f'internal_response_id={internal_response_id} '
                  #f'external_response_id={external_response_id} response_type={response_type}')
            #print(f'  mode_num={mode_num} subcase={subcase} seid={seid} (CEIG)')
            response = {'name': 'ceig', 'mode_num': mode_num, 'subcase': subcase, 'seid': seid}

        elif response_type == 19:
            # RTYPE=19 Equivalent radiated power
            # 4 PANEL      I Element identification number
            # 5 FLAG       I A subcase ID based code. +: magnitude; –:density
            # 6 SUBCASE    I Subcase identification number
            # 7 FREQUENCY RS Frequency
            # 8 UNDEF None
            # 9 SEID       I Superelement identification number
            panel = ints[idata+3]
            flag = ints[idata+4]
            subcase = ints[idata+5]
            freq = floats[idata+6]
            seid = ints[idata+8]
            response = {'name': 'equivalent radiated power', 'panel': panel, 'flag': flag,
                        'subcase': subcase, 'freq': freq, 'seid': seid}
            #print(f'  panel={panel} flag={flag} subcase={subcase} freq={freq} seid={seid} '
                  #'(equivalent radiated power)')
        elif response_type == 20:
            # RTYPE=20 Frequency response displacement
            # 4 GRID    I Grid identification number
            # 5 COMP    I Displacement component number
            # 6 SUBCASE I Subcase identification number
            # 7 FREQ    RS Frequency
            # 8 UNDEF None
            # 9 SEID    I Superelement identification number
            grid = ints[idata+3]
            comp = ints[idata+4]
            subcase = ints[idata+5]
            freq = floats[idata+6]
            ply = ints[idata+7]
            unused_freq8 = floats[idata+7] # also freq despite DMAP saying this is unused...
            seid = ints[idata+8]
            response = {'name': 'frequency response displacement', 'grid': grid, 'component': comp,
                        'subcase': subcase, 'ply': ply, 'seid': seid}
            #print(f'  grid={grid} comp={comp} subcase={subcase} freq={freq} '
                  #f'ply={ply} seid={seid} (frequency response displacement)')
        elif response_type == 21:
            # RTYPE=21 frequency response stress???
            # 4 GRID    I Grid identification number
            # 5 COMP    I Displacement component number
            # 6 SUBCASE I Subcase identification number
            # 7 FREQ    RS Frequency
            # 8 UNDEF None
            # 9 SEID    I Superelement identification number
            grid = ints[idata+3]
            comp = ints[idata+4]
            subcase = ints[idata+5]
            freq = floats[idata+6]
            #sixI = ints[idata+6]
            ply = ints[idata+7]
            seid = ints[idata+8]
            response = {'name': 'frequency response stress?', 'grid': grid, 'component': comp,
                        'subcase': subcase, 'ply': ply, 'seid': seid}
            #print(f'  grid={grid} comp={comp} subcase={subcase} freq={freq} '
                  #f'ply={ply} seid={seid} (frequency response stress?)')
        elif response_type == 29:
            # RTYPE=29 PSD displacement
            # 4 GRID    I Grid identification number
            # 5 COMP    I PSD displacement component number
            # 6 SUBCASE I Subcase identification number
            # 7 FREQ    RS Frequency
            # 8 RANDPS  I RANDPS ID
            # 9 SEID I   Superelement identification number
            grid = ints[idata+3]
            comp = ints[idata+4]
            subcase = ints[idata+5]
            freq = floats[idata+6]
            randps = ints[idata+7]
            seid = ints[idata+8]
            response = {'name': 'psd displacement', 'grid': grid, 'component': comp,
                        'subcase': subcase, 'freq': freq, 'randps': randps, 'seid': seid}
            #print(f'  grid={grid} comp={comp} subcase={subcase} freq={freq} '
                  #f'randps={randps} seid={seid} (PSD displacement)')
        elif response_type == 31:
            # RTYPE=31 PSD acceleration
            # 4 GRID    I Grid identification number
            # 5 COMP    I PSD acceleration component number
            # 6 SUBCASE I Subcase identification number
            # 7 FREQ   RS Frequency
            # 8 RANDPS  I RANDPS ID
            # 9 SEID    I Superelement identification number
            grid = ints[idata+3]
            comp = ints[idata+4]
            subcase = ints[idata+5]
            freq = floats[idata+6]
            randps = ints[idata+7]
            seid = ints[idata+8]
            response = {'name': 'psd acceleration', 'grid': grid, 'component': comp,
                        'subcase': subcase, 'freq': freq, 'randps': randps, 'seid': seid}
            #print(f'  grid={grid} comp={comp} subcase={subcase} freq={freq} '
                  #f'randps={randps} seid={seid} (PSD acceleration)')
        elif response_type == 84:
            # RTYPE=84 Aeroelastic flutter damping
            #   4 SUBCASE I Subcase identification number
            #   5 MODE I Mode number
            #   6 DENSITY RS Density
            #   7 MACH RS Mach number
            #   8 VEL RS Velocity
            #   9 SEID I Superelement identification number
            subcase = ints[idata+3]
            mode = ints[idata+4]
            density = floats[idata+5]
            mach = floats[idata+6]
            vel = floats[idata+7]
            seid = ints[idata+8]
            response= {'name': 'aeroelastic flutter damping', 'subcase': subcase, 'mode': mode,
                       'density': density, 'mach': mach, 'velocity': vel, 'seid': seid}

        else:  # pragma: no cover
            print(f'internal_response_id={internal_response_id} '
                  f'external_response_id={external_response_id} response_type={response_type}')
            raise NotImplementedError(response_type)
        response['internal_response_id'] = internal_response_id
        response['external_response_id'] = external_response_id
        response['response_type'] = response_type
        response['iresponse'] = iresp
        response['response_number'] = 1
        responses[iresp] = response
        idata += 9
    return

def dscmcol_dresp2(responses: dict[int, dict[str, Any]],
                   nresponses_dresp2: int,
                   ints: np.ndarray,
                   floats: np.ndarray) -> None:
    """helper for DSCMCOL"""
    if nresponses_dresp2 == 0:
        return
    nresponses = len(responses)

    idata = 0
    for iresp in range(nresponses_dresp2):
        # Word Name Type Description
        # 1 IRID      I Internal response identification number
        # 2 RID       I External response identification number
        # 3 SUBCASE   I Subcase identification number
        # 4 DFLAG     I Dynamic response flag (See Note)
        # 5 FREQTIME RS Frequency or time step
        # 6 SEID      I Superelement identification number
        internal_response_id = ints[idata]
        external_response_id = ints[idata+1]
        subcase = ints[idata+2]
        dflag = ints[idata+3]
        freqtime = floats[idata+4]
        seid = ints[idata+5]
        iresp2 = iresp + nresponses
        response = {
            'iresponse': iresp2,
            'response_number': 2,
            'internal_response_id': internal_response_id,
            'external_response_id': external_response_id,
            'subcase': subcase, 'dflag': dflag,
            'freq': freqtime, 'seid': seid}
        responses[iresp2] = response
        #print(f'internal_response_id={internal_response_id} '
              #f'external_response_id={external_response_id} '
              #f'subcase={subcase} dflag={dflag} freq/time={freqtime} seid={seid}')
        idata += 6
    return


def read_hisadd(op2_reader: OP2Reader) -> None:
    """optimization history (SOL200) table"""
    op2: OP2 = op2_reader.op2
    log = op2.log
    endian = op2_reader._endian
    size = op2_reader.size
    factor: int = op2_reader.factor

    result_name = 'responses.convergence_data'
    is_saved_result = op2._results.is_saved(result_name)
    #read_mode = self.read_mode
    if is_saved_result:
        op2._results._found_result(result_name)
    else:
        op2_reader._skip_table('HISADD', warn=False)
        return
    is_not_saved_result = not is_saved_result

    #slot = op2.get_result(result_name)

    responses = op2.op2_results.responses
    op2.table_name = op2_reader._read_table_name(rewind=False)

    if op2.read_mode == 1 or is_not_saved_result:
        op2_reader.read_markers([-1])
        op2_reader._skip_record()
        op2_reader.read_3_markers([-2, 1, 0])
        op2_reader._skip_record()
        op2_reader.read_3_markers([-3, 1, 0])

        if is_saved_result and responses.convergence_data is None:
            data = op2_reader._read_record()
            ndvs = len(data) // 4 - 7
            responses.convergence_data = Convergence(ndvs)
        #elif is_not_saved_result:
            #self._skip_record()
        else:
            op2_reader._skip_record()
            responses.convergence_data.n += 1

        op2_reader.read_markers([-4, 1, 0, 0])
        return

    if op2_reader.is_debug_file:
        op2_reader.binary_debug.write('_read_geom_table - %s\n' % op2.table_name)
    #self.log.info('----marker1----')
    op2_reader.read_markers([-1])
    if op2_reader.is_debug_file:
        op2_reader.binary_debug.write('---markers = [-1]---\n')
    data = op2_reader._read_record()  # ()102, 303, 0, 0, 0, 0, 0) date???
    #print('hisadd data1')
    #op2_reader.show_data(data)

    #op2_reader.log.info('----marker2----')
    markers = op2_reader.get_nmarkers(1, rewind=True)
    if op2_reader.is_debug_file:
        op2_reader.binary_debug.write('---marker0 = %s---\n' % markers)
    op2_reader.read_3_markers([-2, 1, 0])
    data = op2_reader._read_record()  # ('HISADD', )
    #print('hisadd data2')
    #op2_reader.show_data(data)

    #log.info('----marker3----')
    op2_reader.read_3_markers([-3, 1, 0])
    data = op2_reader._read_record()

    fmt = mapfmt(endian + b'3i3fi', size)
    (design_iter, iconvergence, conv_result, obj_intial, obj_final,
     constraint_max, row_constraint_max) = unpack(fmt, data[:28 * factor])

    if iconvergence == 1:
        iconvergence = 'soft'
    elif iconvergence == 2:
        iconvergence = 'hard'
    elif iconvergence == 6:
        log.warning('HISADD iconverge=6')
        iconvergence = '???'
    else:  # pragma: no cover
        msg = 'iconvergence=%s\n' % iconvergence
        op2_reader.show_data(data, types='ifs', endian=None)
        raise NotImplementedError(msg)

    if conv_result == 0:
        conv_result = 'no'
    elif conv_result == 1:
        conv_result = 'soft'
    elif conv_result == 2:
        conv_result = 'hard'
    elif conv_result in [3, 4]:
        #self.log.warning('HISADD conv_result=%s' % conv_result)
        # not sure why this happens, but the field is wrong
        # it seems to apply to one step before this one
        conv_result = 'best_design'
    else:
        log.debug('design_iter=%s iconvergence=%s conv_result=%s obj_intial=%s '
                  'obj_final=%s constraint_max=%s row_constraint_max=%s' % (
                      design_iter, iconvergence, conv_result, obj_intial,
                      obj_final, constraint_max, row_constraint_max))
        raise NotImplementedError('conv_result=%s' % conv_result)
    #log.debug('design_iter=%s iconvergence=%s conv_result=%s obj_intial=%s '
              #'obj_final=%s constraint_max=%s row_constraint_max=%s' % (
                  #design_iter, iconvergence, conv_result, obj_intial,
                  #obj_final, constraint_max, row_constraint_max))

    ndvs = len(data) // 4 - 7
    desvar_values = unpack('%sf' % ndvs, data[28:])

    responses.convergence_data.append(
        design_iter, iconvergence, conv_result, obj_intial,
        obj_final, constraint_max, row_constraint_max, desvar_values)
    op2_reader.read_markers([-4, 1, 0, 0])

def read_r1tabrg(op2_reader: OP2Reader) -> None:
    """Reads the R1TABRG design response optimization table"""
    op2: OP2 = op2_reader.op2
    endian = op2_reader._endian
    size = op2_reader.size
    factor: int = op2_reader.factor

    # TODO: I think this is used to handle optimization
    op2._count += 1

    op2.table_name = op2_reader._read_table_name(rewind=False)
    op2_reader.read_markers([-1])

    read_record_ndata = op2_reader.get_skip_read_record_ndata()

    # (101, 221355, 0, 0, 0, 0, 0)
    # (???, nvalues,?, ?, ?, ?, ?)
    data = op2_reader._read_record()
    values = unpack(mapfmt(endian + b'7i', size), data)
    unused_nvalues = values[1]
    #print(values)
    #self.show_data(data, types='i', endian=None)

    #'R1TAB   '
    op2_reader.read_3_markers([-2, 1, 0])
    data, ndata = read_record_ndata()
    assert ndata == 8 * factor, ndata

    itable = -3
    while 1:
        op2_reader.read_3_markers([itable, 1, 0])
        stop_marker = op2_reader.get_marker1(rewind=True)
        if stop_marker == 0:
            break

        data, ndata = read_record_ndata()
        _read_r1tabrg(op2_reader, data, ndata)
        itable -= 1
    stop_marker = op2_reader.get_marker1(rewind=False)

def _read_r1tabrg(op2_reader: OP2Reader, data: bytes, ndata: int) -> int:
    """
    Design Responses:
      - Weight
      - Flutter Speed
      - Stress
      - Strain
      - Displacement

    """
    op2: OP2 = op2_reader.op2
    endian = op2_reader._endian
    size = op2_reader.size
    log = op2.log

    result_name = 'responses'
    if op2._results.is_not_saved(result_name):
        return ndata

    #op2._results._found_result(result_name)
    responses = op2.op2_results.responses

    # create the result object
    if op2.read_mode == 1:
        assert data is not None, data
        assert len(data) > 12, len(data)
        if size == 4:
            response_type, = op2.struct_i.unpack(data[8:12])
        else:
            response_type, = op2.struct_q.unpack(data[16:24])

        #assert response_type in [1, 6, 10, 84], response_type
        if response_type == 1:
            if responses.weight_response is None:
                responses.weight_response = WeightResponse()
            else:
                responses.weight_response.n += 1
        elif response_type == 4:
            #TYPE =4 EIGN or FREQ
            #8 MODE I Mode number
            #9 APRX I Approximation code
            pass
        elif response_type == 5:
            #TYPE =5 DISP
            #8 COMP I Displacement component
            #9 UNDEF None
            #10 GRID I Grid identification number
            if responses.displacement_response is None:
                responses.displacement_response = DisplacementResponse()
            else:
                responses.displacement_response.n += 1
        elif response_type == 6:
            if responses.stress_response is None:
                responses.stress_response = StressResponse()
            else:
                responses.stress_response.n += 1
        elif response_type == 7:
            if responses.strain_response is None:
                responses.strain_response = StrainResponse()
            else:
                responses.strain_response.n += 1

        elif response_type == 8:
            if responses.force_response is None:
                responses.force_response = ForceResponse()
            else:
                responses.force_response.n += 1
        elif response_type == 15:
            # CEIG
            #8 MODE I Mode number
            #9 ICODE I 1: Real component or 2: Imaginary component
            pass
        elif response_type == 23:
            if responses.flutter_response is None:
                responses.fractional_mass_response = FractionalMassResponse()
            else:
                responses.fractional_mass_response.n += 1
        elif response_type == 84:
            if responses.flutter_response is None:
                responses.flutter_response = FlutterResponse()
            else:
                responses.flutter_response.n += 1
        else:
            log.warning('skipping response_type=%d' % response_type)
        return ndata
        #else: # response not added...
            #pass

    # fill the result object
    read_r1tabrg = True
    if read_r1tabrg:
        #self.show_data(data, types='ifs', endian=None)
        fmt = endian + b'3i 8s 4i i 5i' if size == 4 else b'3q 16s 4q q 5q'
        out = unpack(fmt, data)
        # per the R1TAB DMAP page:
        #   all indicies are downshift by 1
        #   indices above out[3] are off by +2 because of the 2 field response_label
        internal_id = out[0]
        dresp_id = out[1]
        response_type = out[2]
        response_label = out[3].strip()
        # -1 for 2 field wide response_label
        region = out[4]
        subcase = out[5]
        type_flag = out[12]  # no meaning per MSC DMAP 2005
        seid = out[13]

        response_label = response_label.decode(op2_reader._encoding)

        if response_type == 1:
            responses.weight_response.add_from_op2(out, log)
        elif response_type == 5:  # DISP
            # out = (1, 101, 5, 'DISP1   ', 101, 1, 3, 0, 1, 0, 0, 0, 0, 0)
            comp = out[6]
            nid = out[8]
            #msg = f'DISP - label={response_label!r} region={region} subcase={subcase} nid={nid}'
            #print(msg)
            #print(out[6:])
            # (3,   0,  1,    0,   0,   0,   0,   0)
            # (???, NA, comp, ???, ???, ???, ???, ???)
            responses.displacement_response.append(
                internal_id, dresp_id, response_label, region,
                subcase, type_flag, seid,
                nid, comp)
        elif response_type == 6:  # STRESS
            #                              -----   STRESS RESPONSES   -----
            #  -------------------------------------------------------------------------------------------
            #   INTERNAL  DRESP1  RESPONSE  ELEMENT   VIEW   COMPONENT  LOWER   INPUT    OUTPUT    UPPER
            #      ID       ID     LABEL       ID    ELM ID     NO.     BOUND   VALUE     VALUE    BOUND
            #  -------------------------------------------------------------------------------------------
            #         21      209  S09L      144747             17       N/A   4.85E+04  5.00E+04  5.00E+04
            # (21, 209, 6, 'S09L    ', 30, 1011, 17, 0, 144747, 0, 0, 0, 0, 0)
            stress_code = out[6]
            pid = out[8]
            #msg = ('STRESS - response_type=%r label=%r region=%s subcase=%s '
                   #'stress_code=%s pid=%s' % (
                       #response_type, response_label, region, subcase,
                       #stress_code, pid))
            responses.stress_response.append(
                internal_id, dresp_id, response_label, region,
                subcase, type_flag, seid,
                stress_code, pid)

        elif response_type == 7:  # STRAIN
            strain_code = out[6]
            pid = out[8]
            responses.strain_response.append(
                internal_id, dresp_id, response_label, region,
                subcase, type_flag, seid,
                strain_code, pid)
        elif response_type == 8:  # FORCE
            #print('internal_id=%s dresp_id=%s response_type=%s response_label=%s'
                  #' region=%s subcase=%s type_flag=%s seid=%s' % (
                #internal_id, dresp_id, response_type, response_label,
                #region, subcase, type_flag, seid
            #))
            force_code = out[6]
            pid = out[8]
            #msg = 'FORCE - label=%r region=%s subcase=%s force_code=%s pid=%s' % (
                #response_label, region, subcase, force_code, pid)
            #print(msg)
            #print(out)
            responses.force_response.append(
                internal_id, dresp_id, response_label, region,
                subcase, type_flag, seid,
                force_code, pid)

        elif response_type == 10:  # CSTRESS
            stress_code = out[6]
            #ply = out[7]
            #pid = out[8]  # is this element id?
            #msg = 'CSTRESS - label=%r region=%s subcase=%s stress_code=%s ply=%s pid=%s' % (
                #response_label, region, subcase, stress_code, ply, pid)
            #print(msg)
        #elif response_type == 10:  # CSTRAIN
            #pass
        elif response_type == 23:  # fractional mass response
            #msg = f'??? - label={response_label!r} region={region} subcase={subcase}'
            #print(msg)
            #print(out[6:])
            responses.fractional_mass_response.append(
                internal_id, dresp_id, response_label, region,
                subcase, type_flag, seid)

        elif response_type == 24:  # FRSTRE
            #8 ICODE I Stress item code
            #9 UNDEF None
            #10 ELID I Element identification number
            #11 FREQ RS Frequency
            #12 IFLAG I Integrated response flag. See Remark 20 of DRESP1.
            #Value is -1 to -6, for SUM, AVG, SSQ,
            pass
        elif response_type == 28:  # RMSACCL
            #8 COMP I RMS Acceleration component
            #9 RANDPS I RANDPS entry identification number
            #10 GRID I Grid identification number
            #11 DMFREQ RS Dummy frequency for internal use
            pass
        elif response_type == 84:
            # FLUTTER  (iii, label, mode, (Ma, V, rho), flutter_id, fff)
            out = unpack(endian + b'iii 8s iii fff i fff', data)
            mode = out[6]
            mach = out[7]
            velocity = out[8]
            density = out[9]
            flutter_id = out[10]
            #msg = ('FLUTTER - _count=%s label=%r region=%s subcase=%s mode=%s '
                   #'mach=%s velocity=%s density=%s flutter_id=%s' % (
                       #op2._count, response_label, region, subcase, mode,
                       #mach, velocity, density, flutter_id))
            responses.flutter_response.append(
                internal_id, dresp_id, response_label, region,
                subcase, type_flag, seid,
                mode, mach, velocity, density, flutter_id)
            #print(msg)
            #self.log.debug(msg)
        #else:
            #self.log.debug('R1TABRG response response_type=%s not supported' % response_type)
            #raise NotImplementedError(response_type)
        assert len(out) == 14, len(out)
    #self.response1_table[self._count] = out
    return ndata
