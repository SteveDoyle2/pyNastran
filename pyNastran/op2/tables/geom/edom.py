"""
defines readers for BDF objects in the OP2 EDOM/EDOMS table
"""
from __future__ import annotations
from struct import Struct
from itertools import count
from typing import Union, TYPE_CHECKING
import numpy as np

from pyNastran.op2.tables.geom.geom_common import GeomCommon
from pyNastran.op2.op2_interface.op2_reader import mapfmt, reshape_bytes_block, reshape_bytes_block_size
from pyNastran.op2.tables.geom.geom4 import _read_spcadd_mpcadd
from .utils import get_minus1_start_end

#if TYPE_CHECKING:  # pragma: no cover
from pyNastran.bdf.cards.optimization import DVPREL1, DVPREL2, DVMREL2, DCONSTR
from pyNastran.op2.errors import DoubleCardError
DSCREEN_INT_TO_RTYPE = {
    1: 'WEIGHT',  # goland_final_test.op2
    3 : 'LAMA',
    4 : 'EIGN',
    5 : 'DISP',
    6 : 'STRESS',
    #8: '???a',  force?
    #9: '???b',
    12: 'FREQ',  # goland_final_test.op2
}

FLAG_TO_RESP_MSC = {
    1 : 'WEIGHT',
    2 : 'VOLUME',
    3 : 'LAMA',
    4 : 'EIGN',
    5 : 'DISP',
    6 : 'STRESS',
    7 : 'STRAIN',
    8 : 'FORCE',
}
FLAG_TO_RESP_NX = {
    1 : 'WEIGHT',
    2 : 'VOLUME',
    3 : 'LAMA',
    4 : 'EIGN',
    5 : 'DISP',
    6 : 'STRESS',
    7 : 'STRAIN',
    8 : 'FORCE',
    9 : 'CFAILURE',
    10 : 'CSTRESS',
    11 : 'CSTRAIN',
    12 : 'FREQ',
    13 : 'SPCFORCE',
    14 : 'ESE',
    15 : 'CEIG',
    17 : 'Compliance',
    19 : 'ERP',
    20: 'FRDISP',
    21: 'FRVELO',
    22: 'FRACCL',
    23: 'FRSPCF',
    24: 'FRSTRE',
    25: 'FRFORC',
    26: 'RMSDISP',
    27: 'RMSVELO',
    28: 'RMSACCL',
    29: 'PSDDISP',
    30: 'PSDVELO',
    31: 'PSDACCL',

    60 : 'TDISP',
    61 : 'TVELO',
    62 : 'TACCL',

    # nx
    #84 : 'FLUTTER',
}

DSCREEN_RTYPE_TO_INT = {value: key for key, value in DSCREEN_INT_TO_RTYPE.items()}

if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.op2.op2_geom import OP2Geom

class EDOM(GeomCommon):
    """defines methods for reading op2 properties"""

    def read_edom4_4(self, data: bytes, ndata: int):
        """
        reads the EDOM table
        SOL 200 design optimization and sensitivity analysis bulk entries.

        """
        return self.op2._read_geom_4(self.edom_map, data, ndata)

    @property
    def size(self) -> int:
        return self.op2.size
    @property
    def factor(self) -> int:
        return self.op2.factor

    def _read_fake(self, data: bytes, n: int) -> int:
        return self.op2._read_fake(data, n)

    def __init__(self, op2: OP2Geom):
        self.op2 = op2

        # F:\work\pyNastran\pyNastran\master2\pyNastran\bdf\test\nx_spike\out_altmdtku4.op2
        # F:\work\pyNastran\pyNastran\master2\pyNastran\bdf\test\nx_spike\out_altd200x7.op2
        # F:\work\pyNastran\pyNastran\master2\pyNastran\bdf\test\nx_spike\out_mdtku1.op2
        # F:\work\pyNastran\pyNastran\master2\pyNastran\bdf\test\nx_spike\out_mcso14.op2
        # F:\work\pyNastran\pyNastran\master2\pyNastran\bdf\test\nx_spike\out_ds105.op2
        # F:\work\pyNastran\pyNastran\master2\pyNastran\bdf\test\nx_spike\out_altcc574.op2
        # F:\work\pyNastran\pyNastran\master2\pyNastran\bdf\test\nx_spike\out_adjoint.op2
        # F:\work\pyNastran\pyNastran\master2\pyNastran\bdf\test\nx_spike\out_mcso18.op2
        # F:\work\pyNastran\pyNastran\master2\pyNastran\bdf\test\nx_spike\out_cqr4optstdis.op2
        # F:\work\pyNastran\pyNastran\master2\pyNastran\bdf\test\nx_spike\out_d200ce12.op2
        #: Optimization Table (I think this is NX-specifc)
        self.edom_map = {
            # are these 3 really EDOM?
            #MAT1DOM(103,1,9944)
            #MAT10DOM(2801,28,9945)
            #MODTRAK(6006,60,477)
            (103, 1, 9944) : ['MAT1DOM', self._read_mat1dom],
            (304, 3, 276) : ['DSCONS', self._read_dscons],
            (404, 4, 277) : ['DVAR', self._read_dvar],
            (504, 5, 246) : ['DVSET', self._read_dvset],

            (4106, 41, 362) : ['DCONSTR', self._read_dconstr],
            #DDVAL(7000,70,563)
            #DRESP3(6700,67,433)
            #(504, 5, 246) : ['???', self._read_fake],

            (3106, 31, 352) : ['DESVAR', self._read_desvar],
            (3206, 32, 353) : ['DLINK', self._read_dlink],
            (3306, 33, 354) : ['DVPREL1', self._read_dvprel1],
            (3406, 34, 355) : ['DVPREL2', self._read_dvprel2],
            #DOPTPRM(4306,43,364)
            (3706, 37, 358) : ['DTABLE', self._read_dtable],
            #(3806, 38, 359) : ['DRESP1', self._read_fake],
            (3806, 38, 359) : ['DRESP1', self._read_dresp1],
            (3906, 39, 360) : ['DRESP2', self._read_fake],
            (4206, 42, 363) : ['DSCREEN', self._read_dscreen],
            (4306, 43, 364) : ['DOPTPRM', self._read_doptprm],
            (4406, 44, 372) : ['DVGRID', self._read_dvgrid],
            #DVSHAP(5006,50,470)
            (5106, 51, 471) : ['DCONADD', self._read_dconadd],
            (5806, 58, 474) : ['DVBSHAP', self._read_fake],
            #DVGEOM(5906,59,356)
            (6006, 60, 477) : ['MODTRAK', self._read_fake],
            #DRESP3(6700,67,433)
            (6100, 61, 429) : ['DVCREL1', self._read_dvcrel1],
            (6200, 62, 430) : ['DVCREL2', self._read_dvcrel2],
            (6300, 63, 431) : ['DVMREL1', self._read_dvmrel1],
            (6400, 64, 432) : ['DVMREL2', self._read_dvmrel2],
            (6006, 60, 477) : ['???', self._read_fake],
            (7000, 70, 563) : ['DCONSTR/DDVAL?', self._read_fake],

            # C:\MSC.Software\simcenter_nastran_2019.2\tpl_post2\s200tpgchbc1.op2
            (6903, 69, 637) : ['DMNCON', self._read_dmncon], # nx
            (7102, 71, 645) : ['DMRLAW', self._read_fake], # nx
            (6803, 68, 636) : ['DVTREL1', self._read_dvtrel1], # nx


            (2801, 28, 9945) : ['MAT10DOM', self._read_fake],
            (5706, 57, 634) : ['DVEREL1', self._read_fake],
            #(6903, 69, 637) : ['???', self._read_fake],
            #(6903, 69, 637) : ['???', self._read_fake],
            #(6903, 69, 637) : ['???', self._read_fake],
            #(6903, 69, 637) : ['???', self._read_fake],
        }

    def _read_dconadd(self, data: bytes, n: int) -> int:
        op2 = self.op2
        datai = np.frombuffer(data[n:], op2.idtype8).copy()
        _read_spcadd_mpcadd(op2, 'DCONADD', datai)
        return len(data)

    def _read_dmncon(self, data: bytes, n: int) -> int:
        """
        Record – DMNCON(6903,69,637)
        NX

        Word Name Type Description
        1 ID      I Unique entry identifier
        2 GID     I GROUP identification number
        3 IOPTION I Number of element layers to consider or for CYCLIC_SYMMETRY
          REPREF=1 indicates repeated sector
          2 indicates reflected sector
        4 FLAG I
        FLAG = 1 Type of constraint = PLANE_SYMMETRY
          5   X RS X axis of point on plane
          6   Y RS Y axis of point on plane
          7   Z RS Z axis of point on plane
          8  N1 RS X component of vector normal to plane
          9  N2 RS Y component of vector normal to plane
          10 N3 RS Z component of vector normal to plane
          11 UNDEF(5) None
        FLAG = 2 Type of constraint = CYCLIC_SYMMETRY
          5   X RS X component of point on axis
          6   Y RS Y component of point on axis
          7   Z RS Z component of point on axis
          8  N1 RS X component of vector defining the rotational axis
          9  N2 RS Y component of vector defining the rotational axis
          10 N3 RS Z component of vector defining the rotational axis
          11 M1 RS X component of vector defining symmetry plane
          12 M2 RS Y component of vector defining symmetry plane
          13 M3 RS Z component of vector defining symmetry plane
          14 NSECT I Number of sectors
          15 UNDEF(5) None
        FLAG = 3 Type of constraint = EXTRUSION
          5 N1 RS X component of vector to define extrusion
          6 N2 RS Y component of vector to define extrusion
          7 N3 RS Z component of vector to define extrusion
          8 UNDEF(5) None
        FLAG = 4 Type of constraint = CASTING
          5 X RS X component of point on casting plane
          6 Y RS Y component of point on casting plane
          7 Z RS Z component of point on casting plane
          8 N1 RS X component of a vector normal to the casting plane
          9 N2 RS Y component of a vector normal to the casting plane
          10 N3 RS Z component of a vector normal to the casting plane
          11 D11 RS X component of a vector which defines the mold removal direction 1 of casting
          12 D12 RS Y component of a vector which defines the mold removal direction 1 of casting
          13 D13 RS Z component of a vector which defines the mold removal direction 1 of casting
          14 D21 RS X component of a vector which defines the mold removal direction 2 of casting
          15 D22 RS Y component of a vector which defines the mold removal direction 2 of casting
          16 D23 RS Z component of a vector which defines the mold removal direction 2 of casting
          17 UNDEF(5) None
        FLAG = 5 Type of constraint = MAX_SIZE
          5 MSIZE RS Maximum size
          6 UNDEF(5) None
        FLAG = 6 Type of constraint = MIN_SIZE
          5 MSIZE RS Minimum size
          6 UNDEF(5) None
        FLAG = 7 Type of constraint = ADDITIVE
          5 ANGLE RS Maximum angle measured from the vector N
          6 MIND RS Minimum allowed dimension
          7 X RS X coordinate of point on base plate
          8 Y RS Y coordinate of point on base plate
          9 Z RS Z coordinate of point on base plate
          10 N1 RS X component of a vector normal to the casting plate in the direction of material addition
          11 N2 RS Y component of a vector normal to the casting plate in the direction of material addition
          12 N3 RS Z component of a vector normal to the casting plate in the direction of material addition
          13 UNDEF(5) None
        FLAG = 8 Type of constraint = CHECKER-BOARD CONTROL (CHBC)
          5 OFF-FLAG RS Negative real number indicates CHBC is off. CHBC is on by default if no negative real OFF-FLAG or if CHBC record segment does not exist.
          6 UNDEF(5) None

        ints    = (6903, 69, 637,
        1, 0, 0, 1,
          252.0, 0,   0,   0,   1.0, 0,   0,   0,   0,   0,   0,
        2, 0, 0, 6,
          0.40, 0, 0, 0, 0, 0)
        floats  = (6903, 69, 637,
        1, 0, 0, 1,
          252.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        2, 0, 0, 6,
          0.40, 0.0, 0.0, 0.0, 0.0, 0.0)
        """
        op2 = self.op2
        op2.to_nx('; DMNCON found')
        ntotal0 = 16 * self.factor # 4 * 4
        ntotal6 = 24 * self.factor # 6 * 4
        ntotal8 = 32 * self.factor # 8 * 4
        ntotal11 = 44 * self.factor # 11 * 4
        ntotal13 = 52 * self.factor # 13 * 4
        ntotal15 = 60 * self.factor # 15 * 4
        ntotal17 = 68 * self.factor # 17 * 4

        struct_9f_6i = Struct(mapfmt(op2._endian + b'9f 6i', self.size))
        struct_4i = Struct(mapfmt(op2._endian + b'4i', self.size))
        struct_3f_5i = Struct(mapfmt(op2._endian + b'3f 5i', self.size))
        struct_6f_5i = Struct(mapfmt(op2._endian + b'6f 5i', self.size))
        struct_f_5i = Struct(mapfmt(op2._endian + b'f 5i', self.size))
        struct_8f_5i = Struct(mapfmt(op2._endian + b'8f 5i', self.size))
        struct_12f_5i = Struct(mapfmt(op2._endian + b'12f 5i', self.size))
        ncards = 0
        #self.show_data(data, types='ifs')
        while n < len(data):
            edata1 = data[n:n+ntotal0]
            constraint_id, group_id, ioption, flag = struct_4i.unpack(edata1)
            op2.log.debug(f'   constraint_id={constraint_id} group_id={group_id} ioption={ioption} flag={flag}')
            n += ntotal0
            if flag == 1:
                # plane symmetry
                constraint_type = 'SYMP'
                # 5   X RS X axis of point on plane
                # 6   Y RS Y axis of point on plane
                # 7   Z RS Z axis of point on plane
                # 8  N1 RS X component of vector normal to plane
                # 9  N2 RS Y component of vector normal to plane
                # 10 N3 RS Z component of vector normal to plane
                edata = data[n:n+ntotal11]
                out = struct_6f_5i.unpack(edata)
                x, y, z, nx, ny, nz, *zeros = out
                op2.log.debug(f'    xyz=[{x:g}, {y:g}, {z:g}]; nxyz=[{nx:g}, {ny:g}, {nz:g}]; zeros={zeros}')

                xyz = np.array([x, y, z])
                normal = np.array([nx, ny, nz])
                dmncon = op2.add_dmncon(constraint_id, constraint_type, xyz=xyz, normal=normal)
                n += ntotal11
            elif flag == 2: # CYCLIC_SYMMETRY
                constraint_type = 'SYMC'
                # 5   X RS X component of point on axis
                # 6   Y RS Y component of point on axis
                # 7   Z RS Z component of point on axis
                # 8  N1 RS X component of vector defining the rotational axis
                # 9  N2 RS Y component of vector defining the rotational axis
                # 10 N3 RS Z component of vector defining the rotational axis
                # 11 M1 RS X component of vector defining symmetry plane
                # 12 M2 RS Y component of vector defining symmetry plane
                # 13 M3 RS Z component of vector defining symmetry plane
                # 14 NSECT I Number of sectors
                # 15 UNDEF(5) None
                edata = data[n:n+ntotal15]
                out = struct_9f_6i.unpack(edata)
                x, y, z, n1, n2, n3, m1, m2, m3, nsections, *zeros = out
                xyz = np.array([x, y, z])
                normal = np.array([n1, n2, n3])
                m = np.array([m1, m2, m3])
                op2.log.debug(f'    xyz=[{x:g}, {y:g}, {z:g}]; nxyz={normal}; m={m} zeros={zeros}')
                dmncon = op2.add_dmncon(constraint_id, constraint_type, xyz=xyz,
                                        normal=normal, m=m, nsections=nsections)

                n += ntotal15

            elif flag == 3: # EXTRUSION
                edata = data[n:n+ntotal8]
                constraint_type = 'EXTC'
                out = struct_3f_5i.unpack(edata)
                n1, n2, n3, *null = out
                normal = np.array([n1, n2, n3])
                op2.log.debug(f'    normal=[{n1:g}, {n2:g}, {n3:g}]; null={null}')
                dmncon = op2.add_dmncon(constraint_id, constraint_type, normal=normal)
                n += ntotal8
                # 5 N1 RS X component of vector to define extrusion
                # 6 N2 RS Y component of vector to define extrusion
                # 7 N3 RS Z component of vector to define extrusion
                # 8 UNDEF(5) None

            elif flag == 4: # CASTING
                # 5 X RS X component of point on casting plane
                # 6 Y RS Y component of point on casting plane
                # 7 Z RS Z component of point on casting plane
                # 8 N1 RS X component of a vector normal to the casting plane
                # 9 N2 RS Y component of a vector normal to the casting plane
                # 10 N3 RS Z component of a vector normal to the casting plane
                # 11 D11 RS X component of a vector which defines the mold removal direction 1 of casting
                # 12 D12 RS Y component of a vector which defines the mold removal direction 1 of casting
                # 13 D13 RS Z component of a vector which defines the mold removal direction 1 of casting
                # 14 D21 RS X component of a vector which defines the mold removal direction 2 of casting
                # 15 D22 RS Y component of a vector which defines the mold removal direction 2 of casting
                # 16 D23 RS Z component of a vector which defines the mold removal direction 2 of casting
                # 17 UNDEF(5) None
                edata = data[n:n+ntotal17]
                n += ntotal17
                constraint_type = 'CDID'
                x, y, z, n1, n2, n3, d11, d12, d13, d21, d22, d23, *null = struct_12f_5i.unpack(edata)
                xyz = np.array([x, y, z])
                normal = np.array([n1, n2, n3])
                d = [d11, d12, d13, d21, d22, d23]
                op2.log.warning(f' DMNCON: CDID/CASTING: xyz=[{x}, {y}, {z}]; n=[{n1}, {n2}, {n3}]; d={d}; null={null}')
                dmncon = op2.add_dmncon(constraint_id, constraint_type, xyz=xyz, normal=normal, d=d)

            elif flag in {5, 6}:
                # FLAG = 5 Type of constraint = MAX_SIZE
                # FLAG = 6 Type of constraint = MIN_SIZE
                #  5 MSIZE RS Minimum size
                #  6 UNDEF(5) None
                if flag == 5:
                    constraint_type = 'MAXS'
                elif flag == 6:
                    constraint_type = 'MINS'
                else:  # pragma: no cover
                    raise RuntimeError(flag)

                edata = data[n:n+ntotal6]
                out = struct_f_5i.unpack(edata)
                size, *zeros = out
                op2.log.debug(f'    size={size:g}; zeros={zeros}')
                dmncon = op2.add_dmncon(constraint_id, constraint_type, size=size)
                n += ntotal6
            elif flag == 7:  # additive / ADDM
                # 5 ANGLE RS Maximum angle measured from the vector N
                # 6 MIND RS Minimum allowed dimension
                # 7 X RS X coordinate of point on base plate
                # 8 Y RS Y coordinate of point on base plate
                # 9 Z RS Z coordinate of point on base plate
                # 10 N1 RS X component of a vector normal to the casting plate in the direction of material addition
                # 11 N2 RS Y component of a vector normal to the casting plate in the direction of material addition
                # 12 N3 RS Z component of a vector normal to the casting plate in the direction of material addition
                # 13 UNDEF(5) None
                edata13 = data[n:n+ntotal13]
                constraint_type = 'ADDM'
                angle, mind, x, y, z, n1, n2, n3, *null = struct_8f_5i.unpack(edata13)
                op2.log.debug(f'    angle={angle:g}; mind={mind}; xyz=[{x}, {y}, {z}]; n=[{n1}, {n2}, {n3}]; null={null}')
                xyz = np.array([x, y, z])
                normal = np.array([n1, n2, n3])
                dmncon = op2.add_dmncon(constraint_id, constraint_type, angle=angle, mind=mind, xyz=xyz, normal=normal)
                n += ntotal13
            elif flag == 8:
                # FLAG = 8 Type of constraint = CHECKER-BOARD CONTROL (CHBC)
                #   5 OFF-FLAG RS Negative real number indicates CHBC is off. CHBC is on by default if no negative real OFF-FLAG or if CHBC record segment does not exist.
                #   6 UNDEF(5) None
                constraint_type = 'CHBC'
                edata = data[n:n+ntotal6]
                off_flag, *zeros = struct_f_5i.unpack(edata)
                op2.log.debug(f'    off_flag={off_flag:g}; zeros={zeros}')
                dmncon = op2.add_dmncon(constraint_id, constraint_type, off_flag=off_flag)
                n += ntotal6
            else:
                raise RuntimeError(flag)
            str(dmncon)
            ncards += 1
            del flag
        op2.card_count['DMNCON'] = ncards
        assert n == len(data), f'n={n}; ndata={len(data)}'
        return n

    def _read_dvtrel1(self, data: bytes, n: int) -> int:
        """
        Record – DVTREL1(6803,68,636)
        NX

        Word Name Type Description
        1 ID       I     Identification number
        2 LABEL(2) CHAR4 Label
        4 GID      I     Referenced GROUP bulk entry identification number
        5 STATE    I     0=ACTIVE, 1=FROZEN
        6 DSVFLG   I     Flag normally to be added by way of the punch file bulk output
        7 UNDEF(3) None
        10 DVID   I Blank or identification number for the 1st auto-generated DESVAR entry
        11 COEF  RS Coefficient in the expression P=COEF*DV (currently inactive)
        12 UNDEF(6) None
        """
        op2 = self.op2
        op2.to_nx('; DVTREL1 found')
        ntotal = 68 * self.factor # 17 * 4
        if self.size == 4:
            struct1 = Struct(op2._endian + b'i 8s 7i f 6i')
        else:
            struct1 = Struct(op2._endian + b'q 16s 7q d 6q')
        ndatai = len(data) - n
        ncards = ndatai // ntotal

        for unused_icard in range(ncards):
            edata = data[n:n+ntotal]
            #self.show_data(edata)
            out = struct1.unpack(edata)
            dvtrel_id, label_bytes, group_id, state_int, dsv_flag, undef1, undef2, undef3, dvid, coeff, *undef = out
            state = 'ACTIVE' if state_int == 0 else 'FROZEN'
            label_bytes = label_bytes.replace(b'\x00', b' ')
            label = label_bytes.decode('latin1').strip()
            #print("label = %r" % label)
            op2.add_dvtrel1(dvtrel_id, label, group_id,
                            state=state, dsv_flag=dsv_flag, dvid1=dvid)
            #op2.log.warning(f'DVTREL1 dvtrel_id={dvtrel_id} label_bytes={label!r} group_id={group_id} state={state!r} dsv_flag={dsv_flag} undef123=[{undef1}, {undef2}, {undef3}], dvid={dvid} coeff={coeff} undef={undef}')
            n += ntotal
            ncards += 1
        op2.card_count['DVTREL1'] = ncards
        return n

    def _read_dconstr(self, data: bytes, n: int) -> int:
        op2 = self.op2
        card_name = 'DCONSTR'
        card_obj = DCONSTR
        methods = {
            28 : self._read_dconstr_28,  # msc
            32 : self._read_dconstr_32,  # nx
        }
        try:
            n = op2.reader_geom2._read_double_card(
                card_name, card_obj,
                op2._add_methods._add_dconstr_object,
                methods, data, n)
        except DoubleCardError:
            raise
            #self.op2.log.warning(f'try-except {card_name}')
            #n = self._read_split_card(data, n,
                                      #self._read_cquad8_current, self._read_cquad8_v2001,
                                      #card_name, op2.add_op2_element)
        return n

    def _read_dconstr_28(self, card_obj: DCONSTR, data: bytes, n: int) ->  tuple[int, list[DCONSTR]]:
        r"""
        Record – DCONSTR(4106,41,362) - MSC

        Word Name Type Description
        1 DCID   I
        2 RID    I
        3 LALLOW RS
        4 UALLOW RS
        5 LOWFQ  RS Freq Range/Low End
        6 HIGHFQ RS Freq Range/High End
        7 DTYPE  I  Data type for allowables

        # bugs\msc_dscmcol>test_op2 goland_final_test.op2
        ints    = (70002, 4, -1.0e+20, 1.0e+35, 0,   1.0e+20, 0)
        floats  = (70002, 3, -1.0e+20, 1.0e+35, 0.0, 1.0e+20, 0.0)
        """
        op2 = self.op2

        ntotal = 28 * self.factor # 7 * 4
        struct1 = Struct(mapfmt(op2._endian + b'ii 4f i', self.size))
        ndatai = len(data) - n
        ncards = ndatai // ntotal
        assert ndatai % ntotal == 0

        dconstrs = []
        for unused_icard in range(ncards):
            edata = data[n:n+ntotal]
            out = struct1.unpack(edata)
            oid, dresp_id, lallow, uallow, lowfq, highfq, dtype = out
            #print(oid, dresp_id, lallow, uallow, lowfq, highfq, ltid, utid)
            assert oid > 0
            assert dtype == 0, dtype
            #lid = ltid if ltid != 0 else lallow
            #uid = utid if utid != 0 else uallow

            dconstr = DCONSTR(oid, dresp_id, # lid=lid, uid=uid,
                              lowfq=lowfq, highfq=highfq)
            #dconstr = op2.add_dconstr(oid, dresp_id, # lid=lid, # uid=uid,
                                       #lowfq=lowfq, highfq=highfq)
            dconstr.validate()
            str(dconstr)
            #print(dconstr)
            n += ntotal
            dconstrs.append(dconstr)
        assert n == len(data), f'n={n} ndata={len(data)}'
        op2.to_msc('; DCONSTR-28 found')
        return n, dconstrs

    def _read_dconstr_32(self, card_obj: DCONSTR, data: bytes, n: int) -> tuple[int, list[DCONSTR]]:
        """
        Record – DCONSTR(4106,41,362) - NX
        Design constraints.

        Word Name Type Description
        1 DCID    I Design constraint set identification number
        2 RID     I DRESPi entry identification number
        3 LALLOW RS Lower bound on the response quantity. Undefined if
                    LTID is nonzero
        4 UALLOW RS Upper bound on the response quantity. Undefined if
                    UTID is nonzero
        5 LOWFQ  RS Low end of frequency range in Hz
        6 HIGHFQ RS High end of frequency range in Hz
        7 LTID    I Identification number of TABLEDi entry giving lower
                    bound on the response quantity as a function of
                    frequency or 0 if not specified
        8 UTID    I Identification number of TABLEDi entry giving upper
                    bound on the response quantity as a function of
                    frequency or 0 if not specified

        data  = (50, 2, 0.0016, 0.0018, 0.0, 1.0e+20, 0, 0)

        """
        op2 = self.op2
        ntotal = 32 * self.factor # 8 * 4
        struct1 = Struct(mapfmt(op2._endian + b'ii 4f ii', self.size))
        ndatai = len(data) - n
        ncards = ndatai // ntotal
        assert ndatai % ntotal == 0

        dconstrs = []
        for unused_icard in range(ncards):
            edata = data[n:n+ntotal]
            out = struct1.unpack(edata)
            oid, dresp_id, lallow, uallow, lowfq, highfq, ltid, utid = out
            #print(oid, dresp_id, lallow, uallow, lowfq, highfq, ltid, utid)
            assert oid > 0
            lid = ltid if ltid != 0 else lallow
            uid = utid if utid != 0 else uallow
            dconstr = DCONSTR(oid, dresp_id, lid=lid, uid=uid,
                              lowfq=lowfq, highfq=highfq)
            dconstr.validate()
            str(dconstr)
            #print(dconstr)
            n += ntotal
            dconstrs.append(dconstr)
        assert n == len(data), f'n={n} ndata={len(data)}'
        op2.to_msc('; DCONSTR-32 found')
        return n, dconstrs

    def _read_dscreen(self, data: bytes, n: int) -> int:
        """
        DSCREEN(4206, 42, 363)
        Design constraint screening data.

        Word Name Type Description
        1 RTYPE I Response type for which the screening criteria apply
        2 TRS  RS Truncation threshold
        3 NSTR  I Maximum number of constraints to be retained per region per load case

        data = (5, -0.70, 10)
        """
        op2 = self.op2
        ntotal = 12 * self.factor # 3*4
        struct1 = Struct(mapfmt(op2._endian + b'ifi', self.size))
        ndatai = len(data) - n
        ncards = ndatai // ntotal

        msg = ''
        for unused_icard in range(ncards):
            edata = data[n:n+ntotal]
            out = struct1.unpack(edata)
            rtype_int, trs, nstr = out
            n += ntotal

            if rtype_int in DSCREEN_INT_TO_RTYPE:
                rtype = DSCREEN_INT_TO_RTYPE[rtype_int]
            elif rtype_int == 7:  # STRAIN/FORCE/EQUA?
                # C:\MSC.Software\simcenter_nastran_2019.2\tpl_post1\mereiglc.op2
                #DSCREEN,DISP,-1000.0,20
                #DSCREEN,STRESS,-1000.0,20
                #DSCREEN,STRAIN,-1000.0,20
                #DSCREEN,FORCE,-1000.0,20
                #DSCREEN,EQUA,-1000.0,20
                #DSCREEN,EIGN,-1000.0
                #DSCREEN,LAMA,-1000.0

                rtype = 'STRAIN?'
                msg += f'rtype_int={rtype_int}? trs={trs} nstr={nstr}\n'
                continue
            elif rtype_int == 8:  # STRAIN/FORCE/EQUA?
                rtype = 'FORCE?'
                msg += f'rtype_int={rtype_int}? trs={trs} nstr={nstr}\n'
                continue
            elif rtype_int == 9:
                msg += f'rtype_int={rtype_int}? trs={trs} nstr={nstr}\n'
                continue
            elif rtype_int == 91:  # STRAIN/FORCE/EQUA?
                rtype = 'EQUA?'
                #DSCREEN,FORCE,-1000.0,20
                #DSCREEN,EQUA,-1000.0,20
                msg += f'rtype_int={rtype_int}? trs={trs} nstr={nstr}\n'
                continue
            else:
                rtype = "?"
                msg += f'rtype_int={rtype_int}? trs={trs} nstr={nstr}\n'
                continue
                #raise NotImplementedError(f'rtype_int={rtype_int}? trs={trs} nstr={nstr}')
            #op2.log.info(f'rtype_int={rtype_int} trs={trs} nstr={nstr}')
            dscreen = op2.add_dscreen(rtype, trs=trs, nstr=nstr)
            dscreen.validate()
            str(dscreen)
            #print(dscreen.rstrip())
        assert n == len(data), f'n={n} ndata={len(data)}'
        if msg:
            msg2 = 'Error reading DSCREEN\n' + msg
            op2.log.error(msg2)
            raise RuntimeError(msg2)
        return n

    def _read_doptprm(self, data: bytes, n: int) -> int:
        """
        Record – DOPTPRM(4306,43,364)
        Design optimization parameters.

        Word Name Type Description
        1 APRCOD   I Approximation method
        2 IPRINT   I Print control during approximate optimization phase with DOT
        3 DESMAX   I Maximum number of design cycles
        4 METHOD   I DOT optimization method
        5 DELP    RS Fractional change allowed in each property during any
                     optimization design cycle
        6 DPMIN   RS Minimum move limit imposed
        7 PTOL    RS Maximum tolerance on differences allowed between the
                     property values on property entries and the property
                     values calculated from the design variable values on
                     the DESVAR entry
        8 CONV1   RS Relative objective function convergence criterion
        9 CONV2   RS Absolute objective function convergence criterion
        10 GMAX   RS Maximum constraint violation allowed at the
                     converged optimum
        11 DELX   RS Fractional change allowed in each design variable
                     during any optimization cycle
        12 DXMIN  RS Minimum absolute limit on design variable move
        13 DELB   RS Relative finite difference move parameter
        14 GSCAL  RS Constraint normalization factor
        15 CONVDV RS Relative convergence criterion on design variables
        16 CONVPR RS Relative convergence criterion on properties

        17 P1      I Design cycles in which output is printed
        18 P2      I Items to be printed at the design cycles defined
                     by P1
        19 CT     RS Constraint tolerance
        20 CTMIN  RS Constraint violation threshold
        21 DABOBJ RS DOT absolute objective function convergence criterion
        22 DELOBJ RS DOT relative objective function convergence criterion
        23 DOBJ1  RS DOT 1–D search absolute objective limit
        24 DOBJ2  RS DOT 1–D search relative objective limit
        25 DX1    RS DOT 1–D search absolute DV limit
        26 DX2    RS DOT 1–D search relative DV limit
        27 ISCAL   I Design variables are rescaled every ISCAL iterations
        28 ITMAX   I Maximum DOT MFD iterations per cycle
        29 ITRMOP  I Maximum consecutive DOT MFD iterations at convergence
        30 IWRITE  I File number for DOT optimizer printout
        31 IGMAX   I Active constraint counter

        32 JTMAX   I Maximum DOT SLP iterations per cycle
        33 ITRMST  I Maximum consecutive DOT SLP iterations at convergence
        34 JPRINT  I SLP subproblem print within DOT
        35 IPRNT1  I Print scaling factors for design variable vector within DOT
        36 IPRNT2  I DOT 1–D search or miscellaneous information print
        37 JWRITE  I File number on which iteration history is written within DOT
        38 STPSCL RS Scale factor for shape finite difference step sizes
                     applied to all shape design variables
        39 FSDMAX  I Number of FSD cycles to be performed
        40 FSDALP RS Relaxation parameter applied in FSD
        41 DISCOD  I Discrete processing method code
        42 DISBEG  I Design cycle ID for discrete variable processing initiation
        43 PLVIOL  I Flag for handling property limit violation
        44 P2RSET  I ID of a SET1 entry listing constrained responses to
                     be printed if retained
        45 EDVOUT RS Fraction of DVEREL1 DESVARs to be output in f06 file
        46 MXCRTRSP I Flag to handle CSV output

        """
        op2 = self.op2
        #if self.size == 4:
        ntotal = 184 * self.factor # 46 * 4
        struct1 = Struct(mapfmt(op2._endian + b'4i 12f 2i 8f 11i f i f 4i f i', self.size))
        ndatai = len(data) - n
        ncards = ndatai // ntotal

        for unused_icard in range(ncards):
            edata = data[n:n+ntotal]
            out = struct1.unpack(edata)
            (aprcod, iprint, desmax, method, # ints
             delp, dpmin, ptol, conv1, conv2, gmax, delx, dxmin, delb, gscal, convdv, convpr, # floats
             p1, p2, # ints
             ct, ctmin, dabobj, delobj, dobj1, dobj2, dx1, dx2, # floats
             iscal, itmax, itrmop, iwrite, igmax, jtmax, itrmst, jprint, iprnt1, iprnt2, jwrite, # ints
             stpscl, # float
             fsdmax, # int
             fsdalp, # float
             discod, disbeg, plviol, p2rset, # ints
             edvout, # float
             mxcrtrsp) = out # int
            params = {
                # ints
                'APRCOD' : aprcod,
                'IPRINT' : iprint,
                'DESMAX' : desmax,
                'METHOD' : method,
                # floats
                'DELP' : delp,
                'DPMIN' : dpmin,
                'PTOL' : ptol,
                'CONV1' : conv1,
                'CONV2' : conv2,
                'GMAX' : gmax,
                'DELX' : delx,
                'DXMIN' : dxmin,
                'DELB' : delb,
                'GSCAL' : gscal,
                'CONVDV' : convdv,
                'CONVPR' : convpr,
                #  ints
                'P1' : p1,
                'P2' : p2,
                # floats
                'CT' : ct,
                'CTMIN' : ctmin,
                'DABOBJ' : dabobj,
                'DELOBJ' : delobj,
                'DOBJ1' : dobj1,
                'DOBJ2' : dobj2,
                'DX1' : dx1,
                'DX2' : dx2,
                # ints
                'ISCAL' : iscal,
                'ITMAX' : itmax,
                'ITRMOP' : itrmop,
                'IWRITE' : iwrite,
                'IGMAX' : igmax,
                'JTMAX' : jtmax,
                'ITRMST' : itrmst,
                'JPRINT' : jprint,
                'IPRNT1' : iprnt1,
                'IPRNT2' : iprnt2,
                'JWRITE' : jwrite,
                'STPSCL' : stpscl, # float
                'FSDMAX' : fsdmax,
                'FSDALP' : fsdalp, # float
                'DISCOD' : discod,
                'DISBEG' : disbeg,
                'PLVIOL' : plviol,
                'P2RSET' : p2rset,
                'EDVOUT' : edvout, # float
                'MXCRTRSP' : mxcrtrsp,
            }
            doptprm = op2.add_doptprm(params)
            for key, default_value in doptprm.defaults.items():
                if default_value is None:
                    continue
                if key not in params:
                    continue
                value_actual = params[key]
                assert isinstance(default_value, type(value_actual)), f'key={key!r} value={default_value!r} value_actual={value_actual!r}'
                if isinstance(value_actual, int) and value_actual == default_value:
                    del doptprm.params[key]
                elif isinstance(default_value, float) and np.allclose(value_actual, default_value):
                    del doptprm.params[key]
            str(doptprm)
            n += ntotal
        op2.card_count['DOPTPRM'] = ncards
        return n

    def _read_dtable(self, data: bytes, n: int) -> int:
        """
        Record – DTABLE(3706,37,358)
        Table constants.
        Word Name Type Description
        1 LABLi(2) CHAR4 Label for the constant
        3 VALUi       RS Value of the constant
        Words 1 thru 3 repeat until -1 occurs

        """
        op2 = self.op2
        if self.size == 4:
            struct1 = Struct(op2._endian + b'8s f')
        else:
            aaa

        ints = np.frombuffer(data[n:], op2.idtype8).copy()
        #floats = np.frombuffer(data[n:], self.fdtype8).copy()
        istart, iend = get_minus1_start_end(ints)

        ncards = 0
        size = self.size
        ntotal = 12 * self.factor # 3*4
        for (i0, i1) in zip(istart, iend):
            assert ints[i1] == -1, ints[i1]

            default_values = {}
            nfields = (i1 - i0) // 3
            for unused_i in range(nfields):
                edata = data[n:n+ntotal]
                key_bytes, value = struct1.unpack(edata)
                if size == 4:
                    key = key_bytes.decode('latin1').rstrip()
                default_values[key] = value
                n += ntotal
                assert n <= len(data), n
            dtable = op2.add_dtable(default_values)
            str(dtable)
            n += size
            ncards += 1
        op2.card_count['DTABLE'] = ncards
        return n

    def _read_mat1dom(self, data: bytes, n: int) -> int:
        """
        If one or more properties from a MAT1 entry are used as design
        variables, the MAT1DOM record is written to the EDOM data block.
        This is different than the MAT1 record in the MPT data block.

        Word Name Type Description
        1 MID  I MAT1 identification number
        2 FTE  I Format code for Young’s modulus
        3 FTG  I Format code for shear modulus
        4 FTNU I Format code for Poisson’s ratio

        (2, 2, 0, 2)
        """
        assert len(data) == 28, len(data)
        return len(data)

    def _read_dvgrid(self, data: bytes, n: int) -> int:
        """
        Design variable to grid point relation.
        Word Name Type Description
        1 DVID   I DESVAR entry identification number
        2 GID    I Grid point or geometric point identification number
        3 CID    I Coordinate system identification number
        4 COEFF RS Multiplier of the vector defined by N(3)
        5 N1    RS Component of the vector measured in the coordinate system defined by CID
        6 N2    RS Component of the vector measured in the coordinate system defined by CID
        7 N3    RS Component of the vector measured in the coordinate system defined by CID

        """
        op2 = self.op2
        ntotal = 28 * self.factor # 7*4
        struct1 = Struct(mapfmt(op2._endian + b'3i 4f', self.size))

        ncards = (len(data) - n) // ntotal
        for unused_i in range(ncards):
            edata = data[n:n + ntotal]
            dvgrid_id, nid, cid, coeff, *dxyz = struct1.unpack(edata)
            assert len(dxyz) == 3, dxyz

            dvgrid = op2.add_dvgrid(dvgrid_id, nid, dxyz,
                                    cid=cid, coeff=coeff)
            dvgrid.write_card_16()
            n += ntotal
        return n

    def _read_dvprel1(self, data: bytes, n: int) -> int:
        """
        Word Name Type Description
        1 ID          I Unique identification number
        2 TYPE(2) CHAR4 Name of a property entry
        4 PID         I Property entry identification number
        5 FID         I FID number input. Otherwise, either 0 if property
                        name is input, or frequency (RS) if entry is for
                        frequency dependent property. (See Words 9 and 10)
        6 PMIN       RS Minimum value allowed for this property
        7 PMAX       RS Maximum value allowed for this property
        8 C0         RS Constant term of relation
        9 PNAME1  CHAR4 First word of property name, if any, or blanks if
                        FID number is nonzero in Word 5
        10 PNAME2 CHAR4 Second word of property name, if any. Otherwise,
                        either blanks if FID number is nonzero in Word 5,
                        or frequency (RS) if entry is for frequency
                        dependent property. (See Word 5)
        11 DVIDi I DESVAR entry identification number
        12 COEFi RS Coefficient of linear relation
        Words 11 and 12 repeat until -1 occurs
        """
        op2 = self.op2
        ints = np.frombuffer(data[n:], op2.idtype8).copy()
        floats = np.frombuffer(data[n:], op2.fdtype8).copy()
        iminus1 = np.where(ints == -1)[0]

        ncards = 0
        istart = [0] + list(iminus1[:-1] + 1)
        iend = iminus1
        size = self.size
        for (i0, i1) in zip(istart, iend):
            #print(self.show_data(data[n+i0*4:n+i1*4]))
            assert ints[i1] == -1, ints[i1]
            dvprel_id = ints[i0]
            type_bytes = data[n+size:n+3*size]
            prop_type = reshape_bytes_block_size(type_bytes, size=size)

            pid, fid = ints[i0+3:i0+5]
            pmin, pmax, c0 = floats[i0+5:i0+8]
            property_name_bytes = data[n+8*size:n+10*size]
            if fid == 0:
                fid = reshape_bytes_block_size(property_name_bytes, size=size)

            # fid = fidi
            #print(f'dvprel_id={dvprel_id} prop_type={prop_type} pid={pid} fid={fid} (pmin, pmax, c0)=({pmin,pmax,c0})')
            desvar_ids = ints[i0+10:i1:2]
            coeffs = floats[i0+11:i1:2]
            # 2 TYPE(2) CHAR4 Name of a property entry
            # 4 PID         I Property entry identification number
            # 5 FID         I FID number input. Otherwise, either 0 if property
            #                 name is input, or frequency (RS) if entry is for
            #                 frequency dependent property. (See Words 9 and 10)
            # 6 PMIN       RS Minimum value allowed for this property
            # 7 PMAX       RS Maximum value allowed for this property
            # 8 C0         RS Constant term of relation
            dvprel = DVPREL1(dvprel_id, prop_type, pid, fid,
                             desvar_ids, coeffs,
                             p_min=pmin, p_max=pmax, c0=c0,
                             validate=True)
            if dvprel_id in op2.dvprels:
                dvprel_old = op2.dvprels[dvprel_id]
                if dvprel == dvprel_old:
                    pass
                else:
                    op2._add_methods._add_dvprel_object(dvprel)
                    ncards += 1
            dvprel.write_card_16()
            n += (i1 - i0 + 1) * size
        op2.card_count['DVPREL1'] = ncards
        return n

    def _read_dvcrel1(self, data: bytes, n: int) -> int:
        """
        Record – DVCREL1(6100,61,429)
        Design variable to connectivity property relation.
        Word Name   Type Description
        1 ID        I     Unique identification number
        2 TYPE(2)   CHAR4 Name of an element connectivity entry
        4 EID       I     Element identification number
        5 FID       I     Entry is 0
        6 CPMIN     RS    Minimum value allowed for this property
        7 CPMAX     RS    Maximum value allowed for this property
        8 C0        RS    Constant term of relation
        9 CPNAME(2) CHAR4 Name of connectivity property
        11 DVIDi    I     DESVAR entry identification number
        12 COEFi    RS    Coefficient of linear relation
        Words 11 and 12 repeat until -1 occurs
        """
        op2 = self.op2
        ints = np.frombuffer(data[n:], op2.idtype8).copy()
        floats = np.frombuffer(data[n:], op2.fdtype8).copy()
        istart, iend = get_minus1_start_end(ints)
        size = self.size
        for (i0, i1) in zip(istart, iend):
            #self.show_data(data[n+i0*size:n+i1*size], types='ifs')
            assert ints[i1] == -1, ints[i1]
            dvcrel_id = ints[i0]
            elem_type_bytes = data[n+size:n+3*size]
            eid, fid = ints[i0+3:i0+5]
            cp_min, cp_max, c0 = floats[i0+5:i0+8]

            cp_name_bytes = data[n+8*size:n+10*size]
            elem_type = reshape_bytes_block_size(elem_type_bytes, size=size)
            cp_name = reshape_bytes_block_size(cp_name_bytes, size=size)
            assert fid == 0, (dvcrel_id, mid, mat_type_bytes, mp_name, pid, fid)
            #print(dvmrel_id, mid, mat_type_bytes, mp_name, pid, fid)

            desvar_ids = ints[i0+10:i1:2]
            coeffs = floats[i0+11:i1:2]
            dvcrel = op2.add_dvcrel1(
                dvcrel_id, elem_type, dvcrel_id, cp_name, desvar_ids, coeffs,
                cp_min=cp_min, cp_max=cp_max, c0=c0, validate=True)
            dvcrel.write_card_16()
            n += (i1 - i0 + 1) * size
        return n

    def _read_dvcrel2(self, data: bytes, n: int) -> int:
        read_dvcrel2

    def _read_dvmrel1(self, data: bytes, n: int) -> int:
        """
        Design variable to material relation.
        Word Name Type Description
        1 ID            I Unique identification number
        2 TYPE(2)   CHAR4 Name of a material property entry
        4 MID           I Material identification number
        5 FID           I Entry is 0
        6 MPMIN        RS Minimum value allowed for this property
        7 MPMAX        RS Maximum value allowed for this property
        8 C0           RS Constant term of relation
        9 MPNAME(2) CHAR4 Name of material property
        11 DVIDi        I DESVAR entry identification number
        12 COEFi       RS Coefficient of linear relation
        Words 11 and 12 repeat until -1 occurs
        """
        op2 = self.op2
        ints = np.frombuffer(data[n:], op2.idtype8).copy()
        floats = np.frombuffer(data[n:], op2.fdtype8).copy()
        istart, iend = get_minus1_start_end(ints)
        size = self.size
        for (i0, i1) in zip(istart, iend):
            #self.show_data(data[n+i0*size:n+i1*size], types='ifs')
            assert ints[i1] == -1, ints[i1]
            dvmrel_id = ints[i0]
            mat_type_bytes = data[n+size:n+3*size]
            mid = ints[i0+3]
            mp_name_bytes = data[n+8*size:n+10*size]
            mat_type = reshape_bytes_block_size(mat_type_bytes, size=size)
            mp_name = reshape_bytes_block_size(mp_name_bytes, size=size)
            pid, fid = ints[i0+3:i0+5]
            mp_min, mp_max, c0 = floats[i0+5:i0+8]
            assert fid == 0, (dvmrel_id, mid, mat_type_bytes, mp_name, pid, fid)
            #print(dvmrel_id, mid, mat_type_bytes, mp_name, pid, fid)

            desvar_ids = ints[i0+10:i1:2]
            coeffs = floats[i0+11:i1:2]
            dvmrel = op2.add_dvmrel1(dvmrel_id, mat_type, mid, mp_name,
                                     desvar_ids, coeffs,
                                     mp_min=mp_min, mp_max=mp_max, c0=c0,
                                     validate=True)
            dvmrel.write_card_16()
            n += (i1 - i0 + 1) * size
        return n

    def _read_dvprel2(self, data: bytes, n: int) -> int:
        """
        Record – DVPREL2(3406,34,355)
        Design variable to property relation based on a user-supplied equation.

        Word Name Type Description
        1 ID          I Unique identification number
        2 TYPE(2) CHAR4 Name of a property entry
        4 PID         I Property entry identification number
        5 FID         I FID number input. Otherwise, either 0 if property
                        name is input, or frequency (RS) if entry is for
                        frequency dependent property. (See Words 9 and 10)
        6 PMIN       RS Minimum value allowed for this property
        7 PMAX       RS Maximum value allowed for this property
        8 EQID        I DEQATN entry identification number
        9 PNAME1  CHAR4 First word of property name, if any, or blank if
                        FID number is nonzero (Word 5)
        10 PNAME2 CHAR4 Second word of property name, if any. Otherwise,
                        either blanks if FID number is nonzero (See Word 5),
                        or frequency (RS) if entry is for frequency
                        dependent property. (See Word 5)
        11 FLAG I DESVAR/DTABLE
        FLAG = 1000 DESVAR
          12 DVIDi I A DESVAR entry identification number
          Word 12 repeats until -1000
        FLAG = 2000 DTABLE
          12 LABLi(2) CHAR4 Label for a constant on the DTABLE entry
          Words 12 and 13 repeat until -2000
        End flag when -1 occurs

        data = (2, PROD, 101, 4, -1.0e+35, 1.0e+20, 2, '', 1000, 2, -1000,
                                                           2000, L1, -2000)
        """
        op2 = self.op2
        #return  self._read_dvxrel2(data, n, DVPREL2)

        n0 = n
        ints = np.frombuffer(data[n:], op2.idtype8).copy()
        floats = np.frombuffer(data[n:], op2.fdtype8).copy()
        istart, iend = get_minus1_start_end(ints)
        size = self.size
        for (i0, i1) in zip(istart, iend):
            #self.show_data(data[n+i0*size:n+i1*size], types='ifs')
            assert ints[i1] == -1, ints[i1]
            dvprel_id = ints[i0]
            prop_type_bytes = data[n0+(i0+1)*size:n0+(i0+3)*size]
            pid, fid = ints[i0+3:i0+5]
            p_min, p_max = floats[i0+5:i0+7]
            deqation = ints[i0+7]

            #data[n0+iflag*size:n0+(iflag+1)*size])
            prop_name_bytes = data[n0+(i0+8)*size:n0+(i0+10)*size]

            if size == 4:
                prop_type = prop_type_bytes.decode('latin1').rstrip()
                prop_name = prop_name_bytes.decode('latin1').rstrip()
            else:
                asdf
            if prop_name_bytes == b'        ':
                assert fid != 0
                pname_fid = fid
            else:
                assert fid == 0, f'fid={fid} prop_name_bytes={prop_name_bytes}'
                pname_fid = prop_name

            #print(dvprel_id, prop_type, pid, pname_fid, deqation)
            iend, dvids, labels = _read_dvxrel2_flag(data, n0, i0, i1, size, ints)

            #print(dvids, labels)
            dvprel = op2.add_dvprel2(dvprel_id, prop_type, pid,
                                      pname_fid, deqation,
                                      dvids=dvids,
                                      labels=labels,
                                      p_min=p_min, p_max=p_max,
                                      validate=True)
            dvprel.validate()
            #print(dvprel)
            #print('--------------------')

            dvprel.write_card_16()
            n += (i1 - i0 + 1) * size
        return n

    def _read_dvmrel2(self, data: bytes, n: int) -> int:
        """
        Record – DVMREL2(6400,64,432)
        Design variable to material relation based on a user-supplied equation.
        Word Name Type Description
        1 ID            I Unique identification number
        2 TYPE(2)   CHAR4 Name of a material property entry
        4 MID           I Material identification number
        5 FID           I Entry is 0
        6 MPMIN        RS Minimum value allowed for this property
        7 MPMAX        RS Maximum value allowed for this property
        8 EQID          I DEQATN entry identification number
        9 MPNAME(2) CHAR4 Name of material property
        11 FLAG         I DESVAR/DTABLE
        FLAG = 1000 DESVAR
          12 DVIDi I A DESVAR entry identification number
          Word 12 repeats until -1000
        FLAG = 2000 DTABLE
          12 LABLi(2) CHAR4 Label for a constant on the DTABLE entry
          Words 12 and 13 repeat until -2000
        End flag when -1 occurs

        """
        op2 = self.op2
        cls = DVMREL2
        #return  self._read_dvxrel2(data, n, DVMREL2)

    #def _read_dvxrel2(self, data: bytes, n: int, cls) -> int:
        n0 = n
        ints = np.frombuffer(data[n:], op2.idtype8).copy()
        floats = np.frombuffer(data[n:], op2.fdtype8).copy()
        iminus1 = np.where(ints == -1)[0]

        istart = [0] + list(iminus1[:-1] + 1)
        iend = iminus1
        size = self.size
        for (i0, i1) in zip(istart, iend):
            #self.show_data(data[n+i0*size:n+i1*size], types='ifs')
            assert ints[i1] == -1, ints[i1]
            dvmrel_id = ints[i0]
            mat_type_bytes = data[n0+(i0+1)*size:n0+(i0+3)*size]
            mid, fid = ints[i0+3:i0+5]
            mp_min, mp_max = floats[i0+5:i0+7]
            deqation = ints[i0+7]

            #data[n0+iflag*size:n0+(iflag+1)*size])
            mp_name_bytes = data[n0+(i0+8)*size:n0+(i0+10)*size]

            if size == 4:
                mat_type = mat_type_bytes.decode('latin1').rstrip()
                mp_name = mp_name_bytes.decode('latin1').rstrip()
            else:
                asdf
            if mp_name_bytes == b'        ':
                assert fid != 0
                mpname_fid = fid
            else:
                assert fid == 0, f'fid={fid} mp_name_bytes={mp_name_bytes}'

            #print(dvmrel_id, mat_type, (mid, fid), (mp_min, mp_max), deqation, mp_name, flag)
            iend, dvids, labels = _read_dvxrel2_flag(data, n0, i0, i1, size, ints)

            #labels = labels.
            #print(dvids, labels)
            card_name = cls.type
            if card_name == 'DVPREL2':
                pid = mid
                dvprel_id = dvmrel_id
                prop_type = mat_type
                pname_fid = mp_name
                dvxrel = op2.add_dvprel2(dvprel_id, prop_type, pid,
                                         pname_fid, deqation,
                                         dvids=dvids,
                                         labels=labels,
                                         p_min=mp_min, p_max=mp_max,
                                         validate=True)

            elif card_name == 'DVMREL2':
                dvxrel = op2.add_dvmrel2(dvmrel_id, mat_type, mid, mp_name,
                                         deqation,
                                         dvids=dvids,
                                         labels=labels,
                                         mp_min=mp_min, mp_max=mp_max,
                                         validate=True)
            dvxrel.validate()
            #print(dvxrel)
            #print('--------------------')

            dvxrel.write_card_16()
            n += (i1 - i0 + 1) * size
        return n

    def _read_dresp1(self, data: bytes, n: int) -> int:
        """
        Word Name Type Description
        1 ID           I Unique entry identifier
        2 LABEL(2) CHAR4 User-defined label
        4 FLAG         I Flag indicating response type
        FLAG = 1 WEIGHT
          5 UNDEF(2) None
          7 REGION I Region identifier for constraint screening
          8 ATTA   I Response attribute (-10 for DWEIGHT which is the topology optimization design weight
          9 ATTB   I Response attribute
          10 MONE  I Entry is -1
        FLAG = 2 VOLUME
          5 UNDEF(2) None
          7 REGION I Region identifier for constraint screening
          8 ATTA   I Response attribute
          9 ATTB   I Response attribute
          10 MONE  I Entry is -1
        FLAG = 3 LAMA
          5 UNDEF(2) None
          7 REGION I Region identifier for constraint screening
          8 ATTA   I Response attribute
          9 ATTB   I Response attribute
          10 MONE  I Entry is -1
        FLAG = 4 EIGN
          5 UNDEF(2) None
          7 REGION I Region identifier for constraint screening
          8 ATTA   I Response attribute
          9 ATTB   I Response attribute
          10 MONE  I Entry is -1
        FLAG = 5 DISP
          5 UNDEF(2) None
          7 REGION I Region identifier for constraint screening
          8 ATTA   I Response attribute
          9 ATTB   I Response attribute
          10 ATTi  I Grid point IDs
          Word 10 repeats until -1 occurs
        FLAG = 6 STRESS
          5 PTYPE(2) CHAR4 Element flag (ELEM) or property entry name
          7 REGION I Region identifier for constraint screening
          8 ATTA   I Response attribute
          9 ATTB   I Response attribute
          10 ATTi  I Element numbers (if Word 5 is ELEM) or property IDs
          Word 10 repeats until -1 occurs
        FLAG = 7 STRAIN
          5 PTYPE(2) CHAR4 Element flag (ELEM) or property entry name
          7 REGION I Region identifier for constraint screening
          8 ATTA   I Response attribute
          9 ATTB   I Response attribute
          10 ATTi  I Element numbers (if Word 5 is ELEM) or property IDs
          Word 10 repeats until -1 occurs
        FLAG = 8 FORCE
          5 PTYPE(2) CHAR4 Element flag (ELEM) or property entry name
          7 REGION I Region identifier for constraint screening
          8 ATTA   I Response attribute
          9 ATTB   I Response attribute
          10 ATTi  I Element numbers (if Word 5 is ELEM) or property IDs
          Word 10 repeats until -1 occurs
        FLAG = 9 CFAILURE
          5 PTYPE(2) CHAR4 Element flag (ELEM) or composite property entry name
          7 REGION I Region identifier for constraint screening
          8 ATTA   I Response attribute
          9 ATTB   I Response attribute
          10 ATTi I Element numbers (if Word 5 is ELEM) or composite property IDs
          Word 10 repeats until -1 occurs
        FLAG = 10 CSTRESS
          5 PTYPE(2) CHAR4 Element flag (ELEM) or composite property entry name
          7 REGION I Region identifier for constraint screening
          8 ATTA I Response attribute
          9 ATTB I Response attribute
          10 ATTi I Element numbers (if Word 5 is ELEM) or composite property IDs
          Word 10 repeats until -1 occurs
        FLAG = 11 CSTRAIN
          5 PTYPE(2) CHAR4 Element flag (ELEM) or composite property entry
          name
          7 REGION I Region identifier for constraint screening
          8 ATTA   I Response attribute
          9 ATTB   I Response attribute
          10 ATTi  I Element numbers (if Word 5 is ELEM) or composite property IDs
          Word 10 repeats until -1 occurs
        FLAG = 12 FREQ
          5 UNDEF(2) None
          7 REGION I Region identifier for constraint screening
          8 ATTA   I Response attribute
          9 ATTB   I Response attribute
          10 MONE  I Entry is -1
        FLAG = 13 SPCFORCE
          5 UNDEF(2) None
          7 REGION I Region identifier for constraint screening
          8 ATTA   I Response attribute
          9 ATTB   I Response attribute
          10 ATTi  I Grid point IDs
          Word 10 repeats until -1 occurs
        FLAG = 14 ESE
          5 PTYPE(2) CHAR4 Element flag (ELEM) or property entry name
          7 REGION I Region identifier for constraint screening
          8 ATTA   I Response attribute
          9 ATTB   I Response attribute
          10 ATTi  I Element numbers (if Word 5 is ELEM) or property IDs
          Word 10 repeats until -1 occurs
        FLAG = 15 CEIG
          5 UNDEF(2) None
          7 REGION I Region identifier for constraint screening
          8 ATTA I Response attribute
          9 ATTB I Response attribute
          10 MONE I Entry is -1
        FLAG = 17 Compliance
          5 UNDEF(2) None
          7 UNDEF I Reserved for SEID for compliance DRESP1
          8 UNDEF(2) None
          10 MONE I Entry is -1
        FLAG = 19 ERP
          5 UNDEF(2) None
          7 REGION I Region identifier
          8 ATTA   I Response attribute
          9 ATTB   I Frequency or real code for character input, or -1=spawn)
          10 ATTi  I Panel SET3 IDs
          Word 10 repeats until -1 occurs
        FLAG = 20 FRDISP
          5 UNDEF(2) None
          7 REGION I Region identifier for constraint screening
          8 ATTA   I Response attribute
          9 ATTB  RS Frequency value; -1 (integer) spawn for all
          frequencies in set; -1.10000E+08 for SUM;
          -1.20000E+08 for AVG; -1.30000E+08 for SSQ;
          -1.40000E+08 for RSS; -1.50000E+08 for MAX;
          -1.60000E+08 for MIN
          10 ATTi I Grid point IDs
          Word 10 repeats until -1 occurs
        FLAG = 21 FRVELO
          5 UNDEF(2) None
          7 REGION I Region identifier for constraint screening
          8 ATTA   I Response attribute
          9 ATTB  RS Frequency value; -1 (integer) spawn for all
          frequencies in set; -1.10000E+08 for SUM;
          -1.20000E+08 for AVG; -1.30000E+08 for SSQ;
          -1.40000E+08 for RSS; -1.50000E+08 for MAX;
          -1.60000E+08 for MIN
          10 ATTi I Grid point IDs
          Word 10 repeats until -1 occurs
        FLAG = 22 FRACCL
          5 UNDEF(2) None
          7 REGION I Region identifier for constraint screening
          8 ATTA   I Response attribute
          9 ATTB  RS Frequency value; -1 (integer) spawn for all
          frequencies in set; -1.10000E+08 for SUM;
          -1.20000E+08 for AVG; -1.30000E+08 for SSQ;
          -1.40000E+08 for RSS; -1.50000E+08 for MAX;
          -1.60000E+08 for MIN
          10 ATTi I Grid point IDs
          Word 10 repeats until -1 occurs
        FLAG = 23 FRSPCF
          5 UNDEF(2) None
          7 REGION I Region identifier for constraint screening
          8 ATTA   I Response attribute
          9 ATTB  RS Frequency value; -1 (integer) spawn for all
          frequencies in set; -1.10000E+08 for SUM;
          -1.20000E+08 for AVG; -1.30000E+08 for SSQ;
          -1.40000E+08 for RSS; -1.50000E+08 for MAX;
          -1.60000E+08 for MIN
          10 ATTi I Grid point IDs
          Word 10 repeats until -1 occurs
        FLAG = 24 FRSTRE
          5 PTYPE(2) CHAR4 Element flag (ELEM) or property entry name
          7 REGION I Region identifier for constraint screening
          8 ATTA   I Response attribute
          9 ATTB  RS Frequency value; -1 (integer) spawn for all
          frequencies in set; -1.10000E+08 for SUM;
          -1.20000E+08 for AVG; -1.30000E+08 for SSQ;
          -1.40000E+08 for RSS; -1.50000E+08 for MAX;
          -1.60000E+08 for MIN
          10 ATTi I Element numbers (if Word 5 is ELEM) or property IDs
          Word 10 repeats until -1 occurs
        FLAG = 25 FRFORC
          5 PTYPE(2) CHAR4 Element flag (ELEM) or property entry name
          7 REGION I Region identifier for constraint screening
          8 ATTA I Response attribute
          9 ATTB RS Frequency value; -1 (integer) spawn for all
          frequencies in set; -1.10000E+08 for SUM;
          -1.20000E+08 for AVG; -1.30000E+08 for SSQ;
          -1.40000E+08 for RSS; -1.50000E+08 for MAX;
          -1.60000E+08 for MIN
          10 ATTi I Element numbers (if Word 5 is ELEM) or property IDs
          Word 10 repeats until -1 occurs
        FLAG = 26 RMSDISP
          5 UNDEF(2) None
          7 REGION I Region identifier for constraint screening
          8 ATTA I Response attribute
          9 ATTB I Random ID
          10 ATTi I Grid point IDs
          Word 10 repeats until -1 occurs
        FLAG = 27 RMSVELO
          5 UNDEF(2) None
          7 REGION I Region identifier for constraint screening
          8 ATTA I Response attribute
          9 ATTB I Random ID
          10 ATTi I Grid point IDs
          Word 10 repeats until -1 occurs
        FLAG = 28 RMSACCL
          5 UNDEF(2) None
          7 REGION I Region identifier for constraint screening
          8 ATTA I Response attribute
          9 ATTB I Random ID
          10 ATTi I Grid point IDs
          Word 10 repeats until -1 occurs
        FLAG = 29 PSDDISP
          5 UNDEF None
          6 PTYPE I Random ID
          7 REGION I Region identifier for constraint screening
          8 ATTA I Response attribute
          9 ATTB RS Frequency value; -1 (integer) spawn for all
          frequencies in set; -1.10000E+08 for SUM;
          -1.20000E+08 for AVG; -1.30000E+08 for SSQ;
          -1.40000E+08 for RSS; -1.50000E+08 for MAX;
          -1.60000E+08 for MIN
          10 ATTi I Grid point IDs
          Word 10 repeats until -1 occurs
        FLAG = 30 PSDVELO
          5 UNDEF None
          6 PTYPE I Random ID
          7 REGION I Region identifier for constraint screening
          8 ATTA I Response attribute
          9 ATTB RS Frequency value; -1 (integer) spawn for all
          frequencies in set; -1.10000E+08 for SUM;
          -1.20000E+08 for AVG; -1.30000E+08 for SSQ;
          -1.40000E+08 for RSS; -1.50000E+08 for MAX;
          -1.60000E+08 for MIN
          10 ATTi I Grid point IDs
          Word 10 repeats until -1 occurs
        FLAG = 60 TDISP
          5 UNDEF(2) None
          7 REGION I Region identifier for constraint screening
          8 ATTA I Response attribute
          9 ATTB RS Time value; -1 (integer) spawn for all time steps
          in set; -1.10000E+08 for SUM; -1.20000E+08 for
          AVG; -1.30000E+08 for SSQ; -1.40000E+08 for
          RSS; -1.50000E+08 for MAX; -1.60000E+08 for MIN
          10 ATTi I Grid point IDs
          Word 10 repeats until -1 occurs
        FLAG = 61 TVELO
          5 UNDEF(2) None
          7 REGION I Region identifier for constraint screening
          8 ATTA I Response attribute
          9 ATTB RS Time value; -1 (integer) spawn for all time steps
          in set; -1.10000E+08 for SUM; -1.20000E+08 for
          AVG; -1.30000E+08 for SSQ; -1.40000E+08 for
          RSS; -1.50000E+08 for MAX; -1.60000E+08 for MIN
          10 ATTi I Grid point IDs
          Word 10 repeats until -1 occurs
        FLAG = 62 TACCL

        [31, 538981700, 538976288, 5, 538976288, 538976288, 0, 3, 0, 5, -1,
         32, 538981444, 538976288, 60, 538976288, 538976288, 0, 3, -1, 4, -1,
         33, 1195984215, 538989640, 1, 538976288, 538976288, 0, 33, -9999, -1]
        """
        op2 = self.op2
        is_nx = False
        is_msc = False
        if op2.is_msc:
            is_msc = True
            flag_to_resp = FLAG_TO_RESP_NX
        else:
            # NX
            is_nx = True
            flag_to_resp = FLAG_TO_RESP_NX

        #self.show_data(data[n:], types='qds')
        ints = np.frombuffer(data[n:], op2.idtype8).copy()
        floats = np.frombuffer(data[n:], op2.fdtype8).copy()
        #print(ints.tolist())
        istart, iend = get_minus1_start_end(ints)
        #if self.size == 4:
            #struct1 = Struct(op2._endian + b'i 8s i')
            #strs = np.frombuffer(data[n:], dtype='|S4')
        #else:
            #struct1 = Struct(op2._endian + b'q 16s q')
            #strs = np.frombuffer(data[n:], dtype='|S8')
        #6i
        #ntotal1 = 16 * self.factor # 4*4

        size = self.size

        def _pick_attbi_attbf(attbi: int, attbf: float) -> Union[float, str]:
            """
            9 ATTB  RS Frequency value; -1 (integer) spawn for all
            frequencies in set; -1.10000E+08 for SUM;
            -1.20000E+08 for AVG; -1.30000E+08 for SSQ;
            -1.40000E+08 for RSS; -1.50000E+08 for MAX;
            -1.60000E+08 for MIN

            """
            if attbi == -1:
                attb = 'ALL'
            else:
                attb = attbf
                #assert attb > -1.0e+8, '%g' % attb
                #print(attbf)
                #ddd
            return attb

        size = self.size
        idresps_to_skip = set()
        for (idresp, i0, i1) in zip(count(), istart, iend):
            assert ints[i1] == -1, ints[i1]
            if idresp in idresps_to_skip:
                #print(f'skipping idresp={idresp}')
                n += (i1 - i0 + 1) * self.size
                continue
            #print(i0, i1)
            #print('ints: ', ints[i0:i1])
            dresp_id = ints[i0]
            label_bytes = data[n+size:n+3*size]
            label = reshape_bytes_block_size(label_bytes, size=size)
            flag = ints[i0+3]
            #d0 = i0 * size
            d1 = i1 * size + size
            op2.show_data(data[n:n+d1], types='ifs')
            try:
                response_type = flag_to_resp[flag]
            except KeyError:
                op2.show_data(data[n:], types='ifs')
                raise RuntimeError(f'dresp_id={dresp_id} label={label!r}')
            print(f'dresp_id={dresp_id} flag={flag}->response_type={response_type!r} label={label!r}')
            if flag == 1:
                # WEIGHT
                # 5 UNDEF(2) None
                # 7 REGION I Region identifier for constraint screening
                # 8 ATTA   I Response attribute (-10 for DWEIGHT which is the topology optimization design weight
                # 9 ATTB   I Response attribute
                # 10 MONE  I Entry is -1
                region, atta, attb = ints[i0+6:i0+9]
                property_type = None
                #response_type = 'WEIGHT'
                assert atta == 33, atta
                assert attb == -9999, attb
                atta = None
                attb = None
                atti = None
            elif flag == 2:
                # FLAG = 2 VOLUME
                #   5 UNDEF(2) None
                #   7 REGION I Region identifier for constraint screening
                #   8 ATTA   I Response attribute
                #   9 ATTB   I Response attribute
                #   10 MONE  I Entry is -1
                property_type = None
                region, atta, attb = ints[i0+6:i0+9]
                atti = None
                attb = None
            elif flag == 3:
                # FLAG = 3 LAMA
                #   5 UNDEF(2) None
                #   7 REGION I Region identifier for constraint screening
                #   8 ATTA   I Response attribute
                #   9 ATTB   I Response attribute
                #   10 MONE  I Entry is -1
                #print(response_type, ints[i0+6:i1], floats[i0+6:i1])
                region, atta, attb = ints[i0+6:i0+9]
                property_type = None
                #response_type = 'EIGN'
                assert atta == 1, atta
                assert attb == 0, attb
                atti = None

            elif flag == 4:
                # FLAG = 4 EIGN
                #   5 UNDEF(2) None
                #   7 REGION I Region identifier for constraint screening
                #   8 ATTA   I Response attribute
                #   9 ATTB   I Response attribute
                #   10 MONE  I Entry is -1
                region, atta, attb = ints[i0+6:i0+9]
                property_type = None
                #response_type = 'EIGN'
                #assert atta == 1, atta
                assert attb == 0, attb
                atti = None
                #atta = None
                #attb = None
            elif flag == 5:
                # DISP
                # 5 UNDEF(2) None
                # 7 REGION I Region identifier for constraint screening
                # 8 ATTA   I Response attribute
                # 9 ATTB   I Response attribute
                # 10 ATTi  I Grid point IDs
                # Word 10 repeats until -1 occurs
                property_type = None
                response_type = 'DISP'
                region, atta, attb = ints[i0+6:i0+9]
                atti = ints[i0+9:i1].tolist()

            elif flag in [6, 7, 9, 11]:
                # FLAG = 6 STRESS
                #FLAG = 9 CFAILURE
                #  5 PTYPE(2) CHAR4 Element flag (ELEM) or property entry name
                #  7 REGION I Region identifier for constraint screening
                #  8 ATTA   I Response attribute
                #  9 ATTB   I Response attribute
                #  10 ATTi  I Element numbers (if Word 5 is ELEM) or property IDs
                #  Word 10 repeats until -1 occurs
                property_type_bytes = data[n+4*size:n+6*size]
                property_type = reshape_bytes_block_size(property_type_bytes, size=size)
                region, atta, attb = ints[i0+6:i0+9]
                atti = ints[i0+9:i1].tolist()
                #print('ptype =', property_type)
                #print('region =', region)
                #print(atta, attb, atti)
                #response_type = 'STRESS'
            #elif flag == 7:
                # FLAG = 7 STRAIN
                #   5 PTYPE(2) CHAR4 Element flag (ELEM) or property entry name
                #   7 REGION I Region identifier for constraint screening
                #   8 ATTA   I Response attribute
                #   9 ATTB   I Response attribute
                #   10 ATTi  I Element numbers (if Word 5 is ELEM) or property IDs
                #   Word 10 repeats until -1 occurs
            #elif flag == 11:
                # FLAG = 11 CSTRAIN
                #   5 PTYPE(2) CHAR4 Element flag (ELEM) or composite property entry name
                #   7 REGION I Region identifier for constraint screening
                #   8 ATTA   I Response attribute
                #   9 ATTB   I Response attribute
                #   10 ATTi  I Element numbers (if Word 5 is ELEM) or composite property IDs
                #   Word 10 repeats until -1 occurs

            elif flag == 12: # FREQ; no is_nx
                # FLAG = 12 FREQ
                #   5 UNDEF(2) None
                #   7 REGION I Region identifier for constraint screening
                #   8 ATTA   I Response attribute
                #   9 ATTB   I Response attribute
                #   10 MONE  I Entry is -1
                property_type = None
                region, atta, attb = ints[i0+6:i0+9]
                atti = None
            elif flag == 15 and is_nx: # CEIG
                # FLAG = 15 CEIG
                #   5 UNDEF(2) None
                #   7 REGION I Region identifier for constraint screening
                #   8 ATTA I Response attribute
                #   9 ATTB I Response attribute
                #   10 MONE I Entry is -1
                #print(ints[i0+6:i1])
                #print(floats[i0+6:i1])
                property_type = None
                region, atta, attb = ints[i0+6:i0+9]
                atti = None
            elif flag == 17 and is_nx: # Compliance
                ## TODO: is this right?
                # FLAG = 17 Compliance
                #   5 UNDEF(2) None
                #   7 UNDEF I Reserved for SEID for compliance DRESP1
                #   8 UNDEF(2) None
                #   10 MONE I Entry is -1
                property_type = None
                region, atta, attb = ints[i0+6:i0+9]
                atti = None
                #print(17, region, atta, attb)
            elif flag == 19 and is_nx: # ERP
                # FLAG = 19 ERP
                #   5 UNDEF(2) None
                #   7 REGION I Region identifier
                #   8 ATTA   I Response attribute
                #   9 ATTB   I Frequency or real code for character input, or -1=spawn)
                #   10 ATTi  I Panel SET3 IDs
                #   Word 10 repeats until -1 occurs
                property_type = None
                region, atta, attb = ints[i0+6:i0+9]
                atti = ints[i0+9:i1].tolist()

            elif flag == 20 and is_nx: # FRDISP
                property_type = None
                # FLAG = 20 FRDISP
                #   5 UNDEF(2) None
                #   7 REGION I Region identifier for constraint screening
                #   8 ATTA   I Response attribute
                #   9 ATTB  RS Frequency value; -1 (integer) spawn for all
                #   frequencies in set; -1.10000E+08 for SUM;
                #   -1.20000E+08 for AVG; -1.30000E+08 for SSQ;
                #   -1.40000E+08 for RSS; -1.50000E+08 for MAX;
                #   -1.60000E+08 for MIN
                #   10 ATTi I Grid point IDs
                #   Word 10 repeats until -1 occurs
                #print(ints[i0+5:i1])
                #print(floats[i0+5:i1])
                #print(ints[i0+6:i1])
                #print(floats[i0+6:i1])
                region, atta, attbi = ints[i0+6:i0+9]
                attbf = floats[i0+8]
                atti = ints[i0+9:i1].tolist()

                attb = _pick_attbi_attbf(attbi, attbf)
                #print(region, atta, attb, atti)
            elif flag == 22 and is_nx: # FRACCL
                # FLAG = 22 FRACCL
                #   5 UNDEF(2) None
                #   7 REGION I Region identifier for constraint screening
                #   8 ATTA   I Response attribute
                #   9 ATTB  RS Frequency value; -1 (integer) spawn for all
                #   frequencies in set; -1.10000E+08 for SUM;
                #   -1.20000E+08 for AVG; -1.30000E+08 for SSQ;
                #   -1.40000E+08 for RSS; -1.50000E+08 for MAX;
                #   -1.60000E+08 for MIN
                #   10 ATTi I Grid point IDs
                #   Word 10 repeats until -1 occurs
                property_type = None
                region, atta, attbi = ints[i0+6:i0+9]
                attbf = floats[i0+8]
                attb = _pick_attbi_attbf(attbi, attbf)
                atti = ints[i0+9:i1].tolist()

            elif flag in [24, 25] and is_nx: # FRSTRE, FRFORC
                # FLAG = 24 FRSTRE
                #   5 PTYPE(2) CHAR4 Element flag (ELEM) or property entry name
                #   7 REGION I Region identifier for constraint screening
                #   8 ATTA   I Response attribute
                #   9 ATTB  RS Frequency value; -1 (integer) spawn for all
                #   frequencies in set; -1.10000E+08 for SUM;
                #   -1.20000E+08 for AVG; -1.30000E+08 for SSQ;
                #   -1.40000E+08 for RSS; -1.50000E+08 for MAX;
                #   -1.60000E+08 for MIN
                #   10 ATTi I Element numbers (if Word 5 is ELEM) or property IDs
                #   Word 10 repeats until -1 occurs
                #
                # FRDISP
                #  5 PTYPE(2) CHAR4 Element flag (ELEM) or property entry name
                #  7 REGION I Region identifier for constraint screening
                #  8 ATTA   I Response attribute
                #  9 ATTB  RS Frequency value; -1 (integer) spawn for all
                #  frequencies in set; -1.10000E+08 for SUM;
                #  -1.20000E+08 for AVG; -1.30000E+08 for SSQ;
                #  -1.40000E+08 for RSS; -1.50000E+08 for MAX;
                #  -1.60000E+08 for MIN
                #  10 ATTi I Element numbers (if Word 5 is ELEM) or property IDs
                #  Word 10 repeats until -1 occurs
                #
                property_type_bytes = data[n+4*size:n+6*size]
                property_type = reshape_bytes_block_size(property_type_bytes, size=size)

                region, atta, attbi = ints[i0+6:i0+9]
                attbf = floats[i0+8]
                attb = _pick_attbi_attbf(attbi, attbf)
                atti = ints[i0+9:i1].tolist()
                print(property_type, region, atta, attb, atti)
            elif flag in {20} and is_msc: # PSDDISP
                # DRESP1       ID   LABEL   RTYPE   PTYPE  REGION    ATTA    ATTB    ATT1
                # DRESP1        11      L1 PSDDISP      91               3   60.00       3
                #5 NTUSED CHAR4
                #6 RPSID I RANDPS ID
                #7 REGION I
                #8 ATTA I
                #9 ATTB RS
                #10 ATTI I
                #Word 10 repeats until End of Record
                property_type_bytes = data[n+4*size:n+6*size]
                property_type = reshape_bytes_block_size(property_type_bytes, size=size)
                print(data[n:+10*size])
                print(ints[i0:i1+1])
                print(floats[i0:i1+1])
                op2.show_data(data[n:])
                randps_id, region, atta = ints[i0+6:i0+9]
                attb = floats[i0+10]
                atti = ints[i0+9:i1].tolist()
                print(f'property_type={property_type!r} randps_id={randps_id} '
                      f'region={region} atta={atta} attb={attb} atti={atti}')
                raise RuntimeError('not done...')
            elif flag in {29} and is_msc: # PSDDISP
                # DRESP1       ID   LABEL   RTYPE   PTYPE  REGION    ATTA    ATTB    ATT1
                # DRESP1        11      L1 PSDDISP      91               3   60.00       3
                property_type, region, atta, attbi = ints[i0+5:i0+9]
                #print(ints[i0+4:i1+5])
                #print(floats[i0+4:i1+5])
                attbf = floats[i0+8]
                attb = _pick_attbi_attbf(attbi, attbf)
                atti = ints[i0+9:i1].tolist()
                print(property_type, region, atta, attb, atti)

                asdf
            elif flag in {29} and is_nx: # PSDDISP
                #FLAG = 29 PSDDISP
                #  5 UNDEF None
                #  6 PTYPE  I Random ID
                #  7 REGION I Region identifier for constraint screening
                #  8 ATTA   I Response attribute
                #  9 ATTB  RS Frequency value; -1 (integer) spawn for all
                #  frequencies in set; -1.10000E+08 for SUM;
                #  -1.20000E+08 for AVG; -1.30000E+08 for SSQ;
                #  -1.40000E+08 for RSS; -1.50000E+08 for MAX;
                #  -1.60000E+08 for MIN
                #  10 ATTi I Grid point IDs
                #  Word 10 repeats until -1 occurs
                property_type, region, atta, attbi = ints[i0+5:i0+9]
                #print(ints[i0+4:i1+5])
                #print(floats[i0+4:i1+5])
                attbf = floats[i0+8]
                attb = _pick_attbi_attbf(attbi, attbf)
                atti = ints[i0+9:i1].tolist()
            elif flag == 31 and is_nx:
                #FLAG = 31 PSDACCL
                #  5 UNDEF None
                #  6 PTYPE I Random ID
                #  7 REGION I Region identifier for constraint screening
                #  8 ATTA I Response attribute
                #  9 ATTB RS Frequency value; -1 (integer) spawn for all
                #  frequencies in set; -1.10000E+08 for SUM;
                #  -1.20000E+08 for AVG; -1.30000E+08 for SSQ;
                #  -1.40000E+08 for RSS; -1.50000E+08 for MAX;
                #  -1.60000E+08 for MIN
                #  10 ATTi I Grid point IDs
                #  Word 10 repeats until -1 occurs
                property_type, region, atta, attbi = ints[i0+5:i0+9]
                #print(ints[i0+4:i1+5])
                #print(floats[i0+4:i1+5])
                attbf = floats[i0+8]
                attb = _pick_attbi_attbf(attbi, attbf)
                atti = ints[i0+9:i1].tolist()


            elif flag == 60 and is_nx:
                #FLAG = 60 TDISP
                #  5 UNDEF(2) None
                #  7 REGION I Region identifier for constraint screening
                #  8 ATTA I Response attribute
                #  9 ATTB RS Time value; -1 (integer) spawn for all time steps
                #  in set; -1.10000E+08 for SUM; -1.20000E+08 for
                #  AVG; -1.30000E+08 for SSQ; -1.40000E+08 for
                #  RSS; -1.50000E+08 for MAX; -1.60000E+08 for MIN
                #  10 ATTi I Grid point IDs
                #  Word 10 repeats until -1 occurs
                property_type = None
                region, atta, attb_int = ints[i0+6:i0+9]
                #print('ints: region, atta, attb=', region, atta, attb_int)

                #region, atta, attb_float = floats[i0+6:i0+9]
                #print('floats: region, atta, attb=', region, atta, attb_float)
                attb = attb_int
                if attb_int == -1:
                    attb = None
                else:
                    attb = floats[i0+8]
                    #print('attb =', attb)
                    raise RuntimeError(('attb', attb))

                #print('ints =', ints[i0+6:i1].tolist())
                #print('floats =', floats[i0+6:i1].tolist())

                #print('---')
                # the grids are on the next idresp, so we'll just skip it on the next round
                i0b = istart[idresp+1]
                i1b = iend[idresp+1]
                atti = grids = ints[i0b:i1b].tolist()
                #print('grids=', grids)
                del grids
                idresps_to_skip.add(idresp+1)
            elif flag == 84 and is_nx:
                # nx flutter
                print('ints =', ints)
                print('floats =', floats)
                continue
            else:
                raise NotImplementedError(flag)

            #print(response_type)
            if property_type == '':
                property_type = None
            if atta == 0:
                atta = None
            if attb == 0:
                attb = None
            if atta is not None:
                atta = int(atta)

            #print(dresp_id, label,
                  #response_type, property_type, region,
                  #atta, attb, atti)
            dresp1 = op2.add_dresp1(dresp_id, label,
                                    response_type, property_type, region,
                                    atta, attb, atti, validate=True)
            print(dresp1)
            dresp1.write_card_16()
            n += (i1 - i0 + 1) * self.size
            del dresp_id, label, response_type, property_type, region, atta, attb, atti

        #for i in range(10):
            #edata = data[n:n+ntotal1]
            #dresp1_id, label_bytes, flag = struct1.unpack(edata)
            ##, undef1, undef2, atta, attb, mone
            #label = reshape_bytes_block(label_bytes).decode('latin1').rstrip()
            #print(dresp1_id, label, flag)
            #self.show_data
        #ddd
        return n

    def _read_dvset(self, data: bytes, n: int) -> int:
        """
        DVSET   13013   PSHELL  4       .02     1.0     13013
        DVSET   13016   PSHELL  4       .02     1.0     13016

        (11013,  902,  9, 4, 2, 1.0,  1.0, 11013, -1)
        (11014,  902,  9, 4, 2, 1.0,  1.0, 11014, -1)
        (13013, 2302, 23, 4, 2, 0.02, 1.0, 13013, -1)
        (13016, 2302, 23, 4, 2, 0.02, 1.0, 13016, -1)

        MSC 2018.2
        Word Name Type Description
        1 VID     I
        2 TYPE(2) I
        4 FIELD   I
        5 I
        =1
          6 PREF I
          7 ALPHA I
        =2
          6 PREF RS
          7 ALPHA RS
        End
        8 PID I
        Word 8 repeats until End of Record

        data = (
            41, 902,    9, 4, 2, 1.0, 1.0, 21, -1,
            42, 302,    3, 5, 2, 1.0, 1.0, 22, -1,
            43, 902,    9, 4, 2, 1.0, 1.0, 23, -1,
            44, 902,    9, 4, 2, 1.0, 1.0, 24, -1,
            45, 302,    3, 5, 2, 1.0, 1.0, 25, -1,
            46, 902,    9, 4, 2, 1.0, 1.0, 26, -1,
            47, 52,    20, 4, 2, 1.0, 1.0, 27, -1,
            48, 5402,  54, -7, 2, 1.0, 1.0, 28, -1,
            48, 5402,  54, -167, 2, 1.0, 1.0, 28, -1,
            49, 5402,  54, -7, 2, 1.0, 1.0, 29, -1,
            49, 5402,  54, -167, 2, 1.0, 1.0, 29, -1,
            50, 52,    20, 4, 2, 1.0, 1.0, 30, -1,
            99, 52,    20, 3, 1, 91, 0/0.0, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, -1,
            99, 5402,  54, 3, 1, 91, 0/0.0, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, -1,
            99, 902,    9, 3, 1, 91, 0/0.0, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, -1,
            100, 402,   4, 3, 3, 538976288,              1.0, 100, -1,   ??? PMASS,None,None
            100, 402,   4, 3, 3, 1.3563156426940112e-19, 1.0, 100, -1,  ???
            410, 902,   9, 4, 2, 1.0, -1.0, 21, -1,
            430, 902,   9, 4, 2, 1.0, -1.0, 23, -1,
            440, 902,   9, 4, 2, 1.0, -1.0, 24, -1,
            460, 902,   9, 4, 2, 1.0, -1.0, 26, -1,
            470, 52,   20, 4, 2, 1.0, -1.0, 27, -1,
            480, 5402, 54, -7, 2, 1.0, -1.0, 28, -1,
            480, 5402, 54, -167, 2, 1.0, -1.0, 28, -1,
            490, 5402, 54, -7, 2, 1.0, -1.0, 29, -1,
            490, 5402, 54, -167, 2, 1.0, -1.0, 29, -1,
            500, 52,   20, 4, 2, 1.0, -1.0, 30, -1,
            999, 52,   20, 3, 1, 91, 0/0.0, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, -1,

            )

        999, 302, 3, 3, 3, 538976288, 1.0, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, -1,  ???
        999, 302, 3, 3, 3, 1.3563156426940112e-19, 1.0, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, -1,  ???

        999, 5402, 54, 3, 1, 91, 0/0.0, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, -1,
        999, 902, 9, 3, 1, 91, 0/0.0, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, -1)
        """
        op2 = self.op2
        #self.show_data(data[12:50], types='ifs')
        #n0 = n
        size = self.size
        #structi = Struct(op2._endian + b'iii ii ff ii')

        ints = np.frombuffer(data[n:], op2.idtype8).copy()
        floats = np.frombuffer(data[n:], op2.fdtype8).copy()
        istart, iend = get_minus1_start_end(ints)
        for (i0, i1) in zip(istart, iend):
            assert ints[i1] == -1, ints[i1]
            #edata = data[n:n + ntotal]
            #out = structi.unpack(edata)
            #print(out)

            dvset_id, dvset_ptype1, dvset_ptype2, field, flag = ints[i0:i0+5]
            if flag == 1:
                pref, alpha = ints[i0+5:i0+7]
            elif flag == 2:
                pref, alpha = floats[i0+5:i0+7]
            elif flag == 3:
                #print(dvset_id, dvset_ptype1, dvset_ptype2, field, flag)
                #print('  ? =', ints[i0+5:i0+7], floats[i0+5:i0+7], data[n0+(i0+4)*size:n0+(i0+7)*size])
                #pref, alpha = '???', '???'
                pref = None
                alpha = None
                #flag3
            else:
                print(dvset_id, dvset_ptype1, dvset_ptype2, field, flag)
                raise NotImplementedError(flag)
            pids = ints[i0+7:i1].tolist()
            #assert field in [3, 4], field

            dvset_ptype = (dvset_ptype1, dvset_ptype2)
            #if dvset_ptype == (902, 9):
                #ptype = 'PSHELL'
            if dvset_ptype == (902, 9):
                ptype = 'PROD'
            elif dvset_ptype == (5402, 54):
                ptype = 'PBEAM'
            elif dvset_ptype == (302, 3):
                ptype = 'PELAS'
            elif dvset_ptype == (52, 20):
                ptype = 'PBAR'
            elif dvset_ptype == (402, 4):
                ptype = 'PMASS'
            elif dvset_ptype == (2302, 23):
                ptype = 'PSHELL'
            elif dvset_ptype == (1002, 10):
                ptype = 'PSHEAR'
            #elif dvset_ptype == (402, 4):
                #ptype = 'PMASS'
            #elif dvset_ptype == (402, 4):
                #ptype = 'PMASS'

            #elif dvset_ptype == (2302, 23):
                #ptype = 'PROD'
            #elif dvset_ptype == (1002, 10):
                #ptype = 'PSHEAR'
            else:
                raise NotImplementedError(f'DVSET={dvset_id} dvset_ptype={dvset_ptype}')
            #print(dvset_id, (ptype, field), flag, (pref, alpha), pids)
            op2.add_dvset(dvset_id, ptype, field, pref, pids, alpha=alpha)
            n += (i1 - i0 + 1) * size
        #op2.log.info(f'geom skipping {self.card_name} in {self.table_name}; ndata={len(data)-12}')
        return n

    def _read_dvar(self, data: bytes, n: int) -> int:
        """
        DVAR    13013   SPARPNL .01     13013

          data = (404, 4, 277,
          11013, 'SPARPNL ', 0.01, 11013, -1)

        """
        op2 = self.op2
        ntotal = 24
        ndatai = len(data) - n
        ncards = ndatai // ntotal
        assert ndatai % ntotal == 0
        structi = Struct(op2._endian + b'i 8s fii')
        for unused_i in range(ncards):
            edata = data[n:n + ntotal]
            #self.show_data(edata, types='ifs')
            #(11013, b'SPRCAPS ', 0.01, 11013, -1)
            #(11014, b'SPRCAPS ', 0.01, 11014, -1)
            #(11015, b'SPRCAPS ', 0.01, 11015, -1)
            #(11016, b'SPRCAPS ', 0.01, 11016, -1)
            out = structi.unpack(edata)
            bid, label_bytes, deltab, vid, minus1 = out
            assert minus1 == -1, out

            assert isinstance(deltab, float), deltab
            label = label_bytes.decode('latin1').rstrip()
            vids = [vid]
            op2.add_dvar(bid, label, vids, deltab=deltab)
            n += ntotal
        return n

    def _read_dscons(self, data: bytes, n: int) -> int:
        """DSCONS

        DSCONS  110131  SPRCAPS STRESS  11013   2       25000.  MAX
        """
        op2 = self.op2
        ndatai = len(data) - n
        # !12
        ntotal = 32
        ncards = ndatai // ntotal
        assert ndatai % ntotal == 0
        structi = Struct(op2._endian + b'i 8s i 2i fi')
        constraint_map = {
            1 : 'DISP',
            2 : 'STRESS',
            3 : 'FORCE',
            4 : 'LAMA',
            5 : 'FREQ',
        }

        for unused_i in range(ncards):
            edata = data[n:n + ntotal]
            #self.show_data(edata, types='ifs')
            #(110131, b'SPRCAPS ', 2, 11013, 2, 25000.0, 0)
            #(110132, b'SPRCAPS ', 2, 11013, 2, -25000.0, 1)
            #(110141, b'SPRCAPS ', 2, 11014, 2, 25000.0, 0)
            #(110142, b'SPRCAPS ', 2, 11014, 2, -25000.0, 1)
            #(110151, b'SPRCAPS ', 2, 11015, 2, 25000.0, 0)
            #(110152, b'SPRCAPS ', 2, 11015, 2, -25000.0, 1)
            #(110161, b'SPRCAPS ', 2, 11016, 2, 25000.0, 0)
            out = structi.unpack(edata)
            dscid, label_bytes, constraint_int, nid_eid, comp, limit, min_max = out
            label = label_bytes.decode('latin1').rstrip()
            try:
                constraint_type = constraint_map[constraint_int]
            except KeyError:
                raise NotImplementedError(f'disp_stress_force={disp_stress_force} out={out}')
            assert min_max in [0, 1], min_max
            out = list(out)

            #print(dscid, label, constraint_type, nid_eid, comp, limit, min_max)
            if min_max == 0:
                opt = 'MAX'
            elif min_max == 1:
                opt = 'MIN'

            layer_id = 1
            op2.add_dscons(dscid, label, constraint_type, nid_eid, comp,
                           limit=limit, opt=opt, layer_id=layer_id)

            n += ntotal
        return n

    def _read_dlink(self, data: bytes, n: int) -> int:
        """
        DLINK(3206,32,353)

        Word Name Type Description
        1 ID     I
        2 DVID   I
        3 C0    RS
        4 CMULT RS
        5 INDV   I
        6 C     RS
        Words 5 through 6 repeat until End of Record

        ints    = (1, 2, 0,   1.0, 1, 1.0, -1)
        floats  = (1, 2, 0.0, 1.0, 1, 1.0, nan)
        """
        op2 = self.op2
        ints = np.frombuffer(data[n:], op2.idtype8).copy()
        floats = np.frombuffer(data[n:], op2.fdtype8).copy()
        istart, iend = get_minus1_start_end(ints)
        for (i0, i1) in zip(istart, iend):
            assert ints[i1] == -1, ints[i1]
            dlink_id, dependent_desvar = ints[i0:i0+2]
            c0, cmult = floats[i0+2:i0+4]
            independent_desvars = ints[i0+4:i1:2]
            coeffs = floats[i0+5:i1:2]
            #print(dlink_id, dependent_desvar, c0, cmult)
            #print(independent_desvars, coeffs)
            assert len(independent_desvars) == len(coeffs)
            assert len(independent_desvars) > 0, independent_desvars
            dlink = op2.add_dlink(dlink_id, dependent_desvar,
                                  independent_desvars,
                                  coeffs,
                                  c0=c0, cmult=cmult)
            #print(dlink)
            str(dlink)
            n += (i1 - i0 + 1) * self.size
        return n

    def _read_desvar(self, data: bytes, n: int) -> int:
        """
        (3106, 31, 352)
        Word Name  Type  Description
        1 ID       I     Unique design variable identification number
        2 LABEL(2) CHAR4 User-supplied name for printing purposes
        4 XINIT    RS    Initial value
        5 XLB      RS    Lower bound
        6 XUB      RS    Upper bound
        7 DELXV    RS    Fractional change allowed for the design variable
                         during approximate optimization
        8 DDVAL    I     ID of a DDVAL entry that provides a set of allowable
                         discrete values
        """
        op2 = self.op2
        if self.size == 4:
            ntotal = 32  # 8*4
            structi = Struct(op2._endian + b'i8s ffff i')
        else:
            ntotal = 64
            structi = Struct(op2._endian + b'q16s dddd q')

        ncards = (len(data) - n) // ntotal
        for unused_i in range(ncards):
            edata = data[n:n + ntotal]
            desvar_id, blabel, xinit, xlb, xub, delx, ddval = structi.unpack(edata)
            label = blabel.decode('ascii')
            if delx == 0:
                delx = None
            if ddval == 0:
                ddval = None
            if desvar_id not in op2.desvars:
                unused_desvar = op2.add_desvar(desvar_id, label, xinit, xlb=xlb, xub=xub,
                                               delx=delx, ddval=ddval, comment='')
            else:
                # duplicate DESVAR
                desvar_temp = op2.add_desvar(1.0, label, xinit, xlb=xlb, xub=xub,
                                             delx=delx, ddval=ddval, comment='')
                del op2.desvars[1.0]
                desvar_temp.desvar_id = desvar_id
                assert desvar_temp == op2.desvars[desvar_id]
            n += ntotal
            #print(desvar)
        op2.card_count['DESVAR'] = ncards
        return n

def _read_dvxrel2_flag(data: bytes, n0: int,
                       i0: int, i1: int,
                       size: int,
                       ints: np.ndarray) -> tuple[list[int], list[str]]:
    """reads the DVxREL2 flag table"""
    flag = ints[i0+10]
    #print(ints[i0+11:])
    #print(floats[i0+11:])
    assert flag in [1000, 2000], flag

    iflag = i0 + 10
    dvids = []
    labels = []
    flags_found = []
    while flag != -1:
        flags_found.append(flag)
        #print(f'i0={i0} iflag={iflag} i1={i1}')
        flag2 = ints[iflag]
        assert flag == flag2
        flag_test, = Struct(b'i').unpack(data[n0+iflag*size:n0+(iflag+1)*size])
        assert flag == flag_test, f'flag={flag} flag_test={flag_test}; n={n}'
        if flag == 1000:
            assert ints[iflag] == 1000, ints[iflag]
            #print('  ', ints[iflag:i1])
            iend = np.where(ints[iflag+1:i1] == -1000)[0][0] + (iflag+1)
            dvids = ints[iflag+1:iend].tolist()
            assert ints[iend] == -1000, ints[iflag+1:i1]
        elif flag == 2000:
            #print('  ', ints[iflag:i1])
            iend = np.where(ints[iflag+1:i1] == -2000)[0][0] + (iflag+1)
            assert ints[iflag] == 2000, ints[iflag]
            assert ints[iend] == -2000, ints[iflag+1:i1]

            labels_bytes = data[n0+(iflag+1)*size:n0+iend*size]
            labels_bytes2 = data[n0+(iflag+1)*size:n0+(iend+1)*size]
            #print('labels_bytes =', labels_bytes)
            nbytes = len(labels_bytes)
            nlabels = nbytes // 8
            assert nbytes % 8 == 0
            assert nlabels > 0, nlabels
            for ilabel in range(nlabels):
                #print(ilabel*size*2, (ilabel+2)*size)
                labels_bytesi = labels_bytes[ilabel*size*2:(ilabel+2)*size]
                label = labels_bytesi.decode('latin1').rstrip()
                assert 1 <= len(str(label)) <= 8, f'label={label}; labels_bytesi={labels_bytesi} labels_bytes={labels_bytes2}'
                labels.append(label)
            #print(labels)
        else:
            raise RuntimeError(flag)
        iflag = iend + 1
        flag = ints[iflag]
        #print(f'\nflag={flag}')
    assert len(flags_found) in [1, 2], flags_found
    return iend, dvids, labels
