"""
defines readers for BDF objects in the OP2 GEOM3/GEOM3S table
"""
#pylint: disable=C0103,C0111,C0301,W0612,R0914,C0326
from __future__ import annotations
from struct import Struct
from typing import TYPE_CHECKING
import numpy as np

from pyNastran.bdf.field_writer_16 import print_card_16
#from pyNastran.bdf.cards.loads.static_loads import (
    #FORCE, FORCE1, FORCE2, GRAV,
    #MOMENT, MOMENT1, MOMENT2,
    #LOAD, PLOAD, PLOAD1, PLOAD2,  #PLOAD3,
    #PLOAD4, )  # PLOAD3,
#from pyNastran.bdf.cards.axisymmetric.loads import PLOADX1 # , PRESAX, TEMPAX, FORCEAX
#from pyNastran.bdf.cards.loads.loads import LSEQ, SLOAD, RFORCE #, DAREA, RANDPS, RFORCE1, LOADCYN
#from pyNastran.bdf.cards.thermal.loads import (
    #QBDY1, QBDY2, QBDY3, TEMP, TEMPD, TEMPP1, QVOL, QHBDY, QVECT)
from pyNastran.op2.op2_interface.op2_reader import mapfmt, reshape_bytes_block
from .utils import get_ints_floats, get_ints, get_ints_floats_strings
from .geom2 import NID_CAP, EID_CAP # , filter_large_element_ids
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.dev.op2_vectorized3.op2_geom import OP2Geom


class GEOM3:
    """defines methods for reading op2 loads"""

    def _add_op2_rigid_element(self, elem):
        """helper method for op2"""
        op2 = self.op2
        ntables = op2.table_names.count(b'GEOM4') + op2.table_names.count(b'GEOM4S')
        eid = elem.eid
        allow_overwrites = (
            ntables > 1 and
            eid in op2.rigid_elements and
            op2.rigid_elements[eid].type == elem.type)
        op2._add_methods.add_rigid_element_object(elem, allow_overwrites=allow_overwrites)

    @property
    def size(self) -> int:
        return self.op2.size
    @property
    def factor(self) -> int:
        return self.op2.factor

    def _read_fake(self, data: bytes, n: int) -> int:
        return self.op2._read_fake(data, n)

    def read_geom3_4(self, data: bytes, ndata: int):
        return self.op2._read_geom_4(self.geom3_map, data, ndata)

    def __init__(self, op2: OP2Geom):
        self.op2 = op2
        self.geom3_map = {
            (11302, 113, 600): ['ACCEL', self._read_accel],    # record 2 - not done
            (11402, 114, 601): ['ACCEL1', self._read_accel1],  # record 3 - not done
            (4201, 42,  18): ['FORCE', self._read_force],      # record 3
            (4001, 40,  20): ['FORCE1', self._read_force1],    # record 4
            (4101, 41,  22): ['FORCE2', self._read_force2],    # record 5
            (4401, 44,  26): ['GRAV', self._read_grav],        # record 7 - buggy
            (4551, 61,  84): ['LOAD', self._read_load],        # record 8
            (3709, 37, 331): ['LOADCYH', self._read_loadcyh],  # record 9 - not done
            (3609, 36, 188): ['LSEQ', self._read_lseq],        # record 12 - not done
            (4801, 48,  19): ['MOMENT', self._read_moment],    # record 13 - not tested
            (4601, 46,  21): ['MOMENT1', self._read_moment1],  # record 14 - not tested
            (4701, 47,  23): ['MOMENT2', self._read_moment2],  # record 15 - not tested
            (5101, 51,  24): ['PLOAD', self._read_pload],      # record 16 - not done
            (6909, 69, 198): ['PLOAD1', self._read_pload1],    # record 17 - buggy
            (6802, 68, 199): ['PLOAD2', self._read_pload2],    # record 18 - buggy
            (7109, 81, 255): ['PLOAD3', self._read_pload3],    # record 19 - not done
            (7209, 72, 299): ['PLOAD4', self._read_pload4],    # record 20 - buggy - g1/g3/g4
            (7001, 70, 278) : ['PLOADX', self._read_ploadx],
            (7309, 73, 351): ['PLOADX1', self._read_ploadx1],  # record 22
            (4509, 45, 239): ['QBDY1', self._read_qbdy1],      # record 24
            (4909, 49, 240): ['QBDY2', self._read_qbdy2],      # record 25
            (2109, 21, 414): ['QBDY3', self._read_qbdy3],      # record 26
            (5509, 55, 190): ['RFORCE', self._read_rforce],    # record 30 - not done
            (5401, 54,  25): ['SLOAD', self._read_sload],      # record 31 - not done
            (5701, 57,  27): ['TEMP', self._read_temp],        # record 32
            (5641, 65,  98): ['TEMPD', self._read_tempd],      # record 33
            (8409, 84, 204): ['TEMPRB', self._read_temprb],    # record 40 - not done
            (8109, 81, 201): ['TEMPP1', self._read_tempp1],    # record 37 - not done
            (8209, 82, 202): ['TEMPP2', self._read_tempp2],    # record 38 - not done
            (8309, 83, 203): ['TEMPP3', self._read_tempp3],    # record 39 - not done
            (8409, 84, 204): ['TEMP4', self._read_tempp4],     # record 40 - not done
            (2309, 23, 416): ['QVOL', self._read_qvol],
            (4309, 43, 233): ['QHBDY', self._read_qhbdy],
            (6609, 66, 9031): ['PEDGE', self._read_pedge],
            (8100, 81, 381): ['CHACAB', self._read_fake],
            (2209, 22, 241): ['QVECT', self._read_qvect],
            (6409, 64, 9032): ['PFACE', self._read_pface],
            (3809, 38, 332): ['LOADCYN', self._read_loadcyn],    # record
            (6209, 62, 390): ['TEMPF', self._read_tempf],    # record
            (10901, 109, 427): ['', self._read_fake],  # record
            (10801, 108, 428): ['GMQVOL', self._read_fake],  # record
            (11329, 113, 9602): ['', self._read_fake],  # record
            (11429, 114, 9603): ['', self._read_fake],  # record
            (11529, 115, 9604): ['', self._read_fake],  # record
            (7002, 70, 254) : ['BOLTFOR', self._read_boltfor],  # record
            (7601, 76, 608) : ['BOLTLD', self._read_boltld],  # record

            # ???
            (6701,67,978): ['PLOADE1', self._read_fake],  # record
            (7002, 85, 254): ['BOLTFOR', self._read_fake],  # record
            (7701, 77, 619): ['DTEMP', self._read_fake],  # record
            (5215, 52, 154): ['PRESAX', self._read_fake],  # record
            (7401, 74, 601): ['ACCEL', self._read_fake],  # record
            (7501, 75, 602): ['ACCEL1', self._read_fake],  # record
            (17600, 176, 627): ['RFORCE2', self._read_fake],  # record
            #(7002, 85, 254): ['BOLTFOR', self._read_fake],  # record
            #(7002, 85, 254): ['BOLTFOR', self._read_fake],  # record

            # nx-specific
            (3909, 39, 333): ['LOADCYT', self._read_fake],  # record
            (9709, 97, 635): ['BOLTFRC', self._read_fake],


            (17300, 173, 615): ['RFORCE1', self._read_fake],
            (7801, 78, 968): ['CRAKTP', self._read_fake],
            (5001, 50, 646): ['FORCDST', self._read_fake],
            (1101, 11, 626): ['INITADD', self._read_fake],
            (8701, 87, 625): ['INITS', self._read_fake],
            (11601, 116, 625): ['PLOADB3', self._read_ploadb3],
            (7901, 79, 967): ['VCEV', self._read_fake],
            (2901, 29, 638): ['INITSO', self._read_fake],
            (9801, 98, 695): ['DRIVER', self._read_fake],
            (11501, 115, 624): ['TEMPBC', self._read_fake],
            #(9709, 97, 635): ['???', self._read_fake],
            #(9709, 97, 635): ['???', self._read_fake],
            #(9709, 97, 635): ['???', self._read_fake],
            #(9709, 97, 635): ['???', self._read_fake],
            #(9709, 97, 635): ['???', self._read_fake],
            #(9709, 97, 635): ['???', self._read_fake],
            #(9709, 97, 635): ['???', self._read_fake],
            #(9709, 97, 635): ['???', self._read_fake],
            #(9709, 97, 635): ['???', self._read_fake],
            #(9709, 97, 635): ['???', self._read_fake],


        }

    def _read_ploadb3(self, data: bytes, n: int) -> int:
        """
        PLOADB3 SID     EID     CID     N1      N2      N3      TYPE    SCALE
                P(A)    P(B)    P(C)
        PLOADB3 10      1       LOCAL   0.0     1.0     0.0     FORCE   1.0     +
        +       100.    100.    100.0
        PLOADB3 10      2       LOCAL   0.0     1.0     0.0     FORCE   1.0     +
        +       100.    100.    100.0

        3i         3i          3f         i  4f
        ints    = (10, 1,  -1, 0, 1.0, 0, 1, 1.0, 100.0, 100.0, 100.0,
                   10, 2,  -1, 0, 1.0, 0, 1, 1.0, 100.0, 100.0, 100.0,
                   10, 3,  -1, 0, 1.0, 0, 1, 1.0, 100.0, 100.0, 100.0,
                   10, 4,  -1, 0, 1.0, 0, 1, 1.0, 100.0, 100.0, 100.0,
                   10, 5,  -1, 0, 1.0, 0, 1, 1.0, 100.0, 100.0, 100.0,
                   10, 6,  -1, 0, 1.0, 0, 1, 1.0, 100.0, 100.0, 100.0,
                   10, 7,  -1, 0, 1.0, 0, 1, 1.0, 100.0, 100.0, 100.0,
                   10, 8,  -1, 0, 1.0, 0, 1, 1.0, 100.0, 100.0, 100.0,
                   10, 9,  -1, 0, 1.0, 0, 1, 1.0, 100.0, 100.0, 100.0,
                   10, 10, -1, 0, 1.0, 0, 1, 1.0, 100.0, 100.0, 100.0,
                   10, 11, -1, 0, 1.0, 0, 1, 1.0, 100.0, 100.0, 100.0,
                   10, 12, -1, 0, 1.0, 0, 1, 1.0, 100.0, 100.0, 100.0)
        floats  = (10, 1, nan, 0.0, 1.0, 0.0, 1, 1.0, 100.0, 100.0, 100.0,
                   10, 2, -1, 0.0, 1.0, 0.0, 1, 1.0, 100.0, 100.0, 100.0,
                   10, 3, nan, 0.0, 1.0, 0.0, 1, 1.0, 100.0, 100.0, 100.0,
                   10, 4, nan, 0.0, 1.0, 0.0, 1, 1.0, 100.0, 100.0, 100.0,
                   10, 5, nan, 0.0, 1.0, 0.0, 1, 1.0, 100.0, 100.0, 100.0,
                   10, 6, nan, 0.0, 1.0, 0.0, 1, 1.0, 100.0, 100.0, 100.0,
                   10, 7, nan, 0.0, 1.0, 0.0, 1, 1.0, 100.0, 100.0, 100.0,
                   10, 8, nan, 0.0, 1.0, 0.0, 1, 1.0, 100.0, 100.0, 100.0,
                   10, 9, nan, 0.0, 1.0, 0.0, 1, 1.0, 100.0, 100.0, 100.0,
                   10, 10, nan, 0.0, 1.0, 0.0, 1, 1.0, 100.0, 100.0, 100.0,
                   10, 11, nan, 0.0, 1.0, 0.0, 1, 1.0, 100.0, 100.0, 100.0,
                   10, 12, nan, 0.0, 1.0, 0.0, 1, 1.0, 100.0, 100.0, 100.0)
        """
        op2 = self.op2
        ntotal = 44 * self.factor  # 11*4
        nelements = (len(data) - n) // ntotal
        structi = Struct(mapfmt(op2._endian + b'3i 3f i 4f', self.size))

        for i in range(nelements):
            datai = data[n:n+ntotal]
            out = structi.unpack(datai)
            (sid, eid, cid, n1, n2, n3, load_type, scale, pa, pb, pc) = out
            #print(out)
            if cid == -1:
                cid = 'LOCAL'
            elif cid == 0:
                cid = 'BASIC'
            elif cid == -2:
                cid = 'ELEMENT'
            elif cid > 1:
                pass
            else:
                raise NotImplementedError(f'PLOADB3 cid={cid}')

            if load_type == 1:
                load_type = 'FORCE'
            elif load_type == 2:
                load_type = 'MOMENT'
            else:
                raise NotImplementedError(f'PLOADB3 load_type={load_type}')
            #self.add_cbeam3(eid, pid, nids, x, g0, wa, wb, wc, tw, s)
            n += ntotal
        #self.show_data(data[n:])
        op2.log.warning('geom skipping PLOADB3')
        op2.card_count['PLOADB3'] = nelements
        return n

    def _read_tempb3(self, data: bytes, n: int) -> int:
        r"""
        $TEMPB3 SID     EID     TA      TB      TC      TPYA    TPZA    TPYB
        $       TPZB    TPYC    TPZC    TCA     TDA     TEA     TFA     TCB
        $       TDB     TEB     TFB     TCC     TDC     TEC     TFC
        $       List of element IDs
        $-------2-------3-------4-------5-------6-------7-------8-------9-------0-------
        TEMPB3  10      1201    10.0    10.0    10.0                            +
        +                               10.0    10.0    10.0    10.0    10.0    +
        +       10.0    10.0    10.0    10.0    10.0    10.0    10.0            +
        +       1202

        C:\MSC.Software\msc_nastran_runs\b3temp01.op2
        floats  = (10, 1201, 10.0, 10.0, 10.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0,
                   10, 1202, 10.0, 10.0, 10.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0,
                   100, 1201, 110.0, 110.0, 110.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 110.0, 110.0, 110.0, 110.0, 110.0, 110.0, 110.0, 110.0, 110.0, 110.0, 110.0, 110.0,
                   100, 1202, 110.0, 110.0, 110.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 110.0, 110.0, 110.0, 110.0, 110.0, 110.0, 110.0, 110.0, 110.0, 110.0, 110.0, 110.0,
                   200, 1201, 100.0, 100.0, 100.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 90.0, 100.0, 110.0, 100.0, 90.0, 100.0, 110.0, 100.0, 90.0, 100.0, 110.0, 100.0,
                   200, 1202, 100.0, 100.0, 100.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 90.0, 100.0, 110.0, 100.0, 90.0, 100.0, 110.0, 100.0, 90.0, 100.0, 110.0, 100.0,
                   300, 1201, 68.0, 23.0, 45.5, 0.0, 0.0, 28.0, 2.5, 14.0, 1.25, 68.0, 91.0, 45.0, 68.0, 48.0, 80.0, 20.0, 23.0, 58.0, 85.5, 32.5, 0.0,
                   300, 1202, 68.0, 23.0, 45.5, 0.0, 0.0, 28.0, 2.5, 14.0, 1.25, 68.0, 91.0, 45.0, 68.0, 48.0, 80.0, 20.0, 23.0, 58.0, 85.5, 32.5, 0.0,
                   300, 1203, 10.0, 50.0, 30.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

        """
        op2 = self.op2
        op2.show_data(data[n:], types='ifs')
        sss

    def _read_accel(self, data: bytes, n: int) -> int:
        """
        Record - ACCEL(7401,74,601)

        Word Name Type Description
        1 SID     I Load set identification number
        2 CID     I Coordinate system identification number
        3 N(3)   RS Components of a vector coordinate system defined by CID
        6 DIR CHAR1 Component direction of acceleration variation
        7 LOCi   RS Location along direction DIR in coordinate system
        8 VALi   RS The load scale factor associated with location LOCi
        Words 7 through 8 repeat until (-1,-1) occurs.
        """
        op2 = self.op2
        op2.show_data(data)
        op2.log.info('geom skipping ACCEL in GEOM3')
        return len(data)

    def _read_accel1(self, data: bytes, n: int) -> int:
        """
        ACCEL1(7501,75,602)

        1 SID    I Load set identification number
        2 CID    I Coordinate system identification number
        3 A     RS Acceleration vector scale factor
        4 N(3)  RS Components of a vector coordinate system defined by CID
        7 GRIDID I Grid ID or THRU or BY code
        Words 7 repeats until (-1) occurs.
        NX/MSC
        """
        op2 = self.op2
        #ntotal = 28  # 7*4
        ints = np.frombuffer(data[n:], dtype='int32').copy()
        floats = np.frombuffer(data[n:], dtype='float32').copy()
        strings = np.frombuffer(data[n:], dtype='|S4')
        i_minus_1s = np.where(ints == -1)[0]

        i0 = 0
        #self.show_data(data[n:])
        for i_minus_1 in i_minus_1s:
            sid = ints[i0]
            cid = ints[i0 + 1]
            scale = floats[i0 + 2]
            n1 = floats[i0 + 3]
            n2 = floats[i0 + 4]
            n3 = floats[i0 + 5]
            nids = [
                strings[i].decode('utf8').strip() if strings[i] in [b'THRU', b'BY  '] else ints[i]
                for i in range(i0+6, i_minus_1)]
            assert nids[-1] > 0
            if op2.is_debug_file:
                op2.binary_debug.write('  ACCEL1=%s\n' % str([sid, cid, scale, n1, n2, n3, nids]))
            accel = op2.add_accel1(sid, scale, [n1, n2, n3], nids, cid=cid)
            accel.validate()
            i0 = i_minus_1 + 1
        return len(data)

    def _read_force(self, data: bytes, n: int) -> int:
        """
        FORCE(4201,42,18) - the marker for Record 3
        """
        op2 = self.op2
        ntotal = 28 * op2.factor # 7*4
        nloads = (len(data) - n) // ntotal
        n, ints, floats = get_ints_floats(data, n, nloads, 7,
                                          size=op2.size, endian=op2._endian)
        load_id = ints[:, 0]
        node_id = ints[:, 1]
        coord_id = ints[:, 2]
        magnitude = floats[:, 3]
        xyz = floats[:, [4, 5, 6]]
        load = op2.force
        load._save(load_id, node_id, coord_id, magnitude, xyz)
        filter_large_node_ids(load)
        load.write()
        #s = Struct(op2._endian + b'iiiffff')
        #for unused_i in range(nloads):
            #out = s.unpack(data[n:n + 28])
            #(sid, g, cid, f, n1, n2, n3) = out
            #if op2.is_debug_file:
                #op2.binary_debug.write('  FORCE=%s\n' % str(out))
            #force = FORCE(sid, g, f, cid=cid, xyz=np.array([n1, n2, n3]))
            #op2._add_methods.add_load_object(force)
            #n += 28
        op2.card_count['FORCE'] = nloads
        return n

    def _read_force1(self, data: bytes, n: int) -> int:
        """
        FORCE1(4001,40,20) - the marker for Record 4
        """
        op2 = self.op2
        ntotal = 20 * self.factor  # 5*4
        nloads = (len(data) - n) // ntotal
        n, ints, floats = get_ints_floats(data, n, nloads, 5,
                                          size=op2.size, endian=op2._endian)
        load_id = ints[:, 0]
        node_id = ints[:, 1]
        magnitude = floats[:, 2]
        nodes = ints[:, [3, 4]]
        load = op2.force1
        load._save(load_id, node_id, magnitude, nodes)
        filter_large_node_ids(load)
        load.write()
        #s = Struct(op2._endian + b'iifii')
        #for unused_i in range(nloads):
            #edata = data[n:n + 20]
            #out = s.unpack(edata)
            #(sid, g, f, n1, n2) = out
            #if op2.is_debug_file:
                #op2.binary_debug.write('  FORCE1=%s\n' % str(out))
            #load = FORCE1.add_op2_data([sid, g, f, n1, n2])
            #op2._add_methods.add_load_object(load)
            #n += 20
        op2.card_count['FORCE1'] = nloads
        return n

    def _read_force2(self, data: bytes, n: int) -> int:
        """
        FORCE2(4101,41,22) - the marker for Record 5
        """
        op2 = self.op2
        ntotal = 28 * self.factor  # 7*4
        nloads = (len(data) - n) // ntotal
        n, ints, floats = get_ints_floats(data, n, nloads, 7,
                                          size=op2.size, endian=op2._endian)
        load_id = ints[:, 0]
        node_id = ints[:, 1]
        magnitude = floats[:, 2]
        nodes = ints[:, [3, 4, 5, 6]]
        nloads = (len(data) - n) // ntotal
        load = op2.force2
        load._save(load_id, node_id, magnitude, nodes)
        filter_large_node_ids(load)
        load.write()

        #s = Struct(op2._endian + b'iif4i')
        #for unused_i in range(nloads):
            #out = s.unpack(data[n:n + 28])
            #(sid, g, f, n1, n2, n3, n4) = out
            #if op2.is_debug_file:
                #op2.binary_debug.write('  FORCE2=%s\n' % str(out))
            #load = FORCE2.add_op2_data([sid, g, f, n1, n2, n3, n4])
            #op2._add_methods.add_load_object(load)
            #n += 28
        op2.card_count['FORCE2'] = nloads
        return n

    def _read_gmload(self, data: bytes, n: int) -> int:
        """GMLOAD"""
        self.op2.log.info('geom skipping GMLOAD in GEOM3')
        return len(data)

    def _read_grav(self, data: bytes, n: int) -> int:
        """
        GRAV(4401,44,26) - the marker for Record 7

        Word Name Type Description
        1 SID I Load set identification number
        2 CID I Coordinate system identification number
        3 A RS Acceleration vector scale factor
        4 N(3) RS Components of a vector coordinate system defined by CID
        7 MB I Bulk Data Section with CID definition: -1=main, 0=partitioned
        """
        op2 = self.op2
        ntotal = 28 * self.factor  # 7*4
        nloads = (len(data) - n) // ntotal

        n, ints, floats = get_ints_floats(data, n, nloads, 7,
                                          size=op2.size, endian=op2._endian)
        load_id = ints[:, 0]
        coord_id = ints[:, 1]
        scale = floats[:, 2]
        N = floats[:, [3, 4, 5]]
        main_bulk = ints[:, 6]
        load = op2.grav
        load._save(load_id, coord_id, scale, main_bulk, N)
        load.write()

        #s = Struct(mapfmt(op2._endian + b'ii4fi', self.size))
        #for unused_i in range(nloads):
            #edata = data[n:n + ntotal]
            #out = s.unpack(edata)
            #if op2.is_debug_file:
                #op2.binary_debug.write('  GRAV=%s\n' % str(out))
            ##(sid, cid, a, n1, n2, n3, mb) = out
            #grav = GRAV.add_op2_data(out)
            #op2._add_methods.add_load_object(grav)
            #n += ntotal
        op2.card_count['GRAV'] = nloads
        return n

    def _read_load(self, data: bytes, n: int) -> int:
        """
        (4551, 61, 84) - the marker for Record 8
        .. todo:: add object
        """
        op2 = self.op2
        ntotal = 16 * self.factor # 4*4
        ntotal2 = 8 * self.factor
        #nloads = (len(data) - n) // ntotal
        #count = 0
        struct_i2fi = Struct(mapfmt(op2._endian + b'iffi', self.size))
        struct_fi = Struct(mapfmt(op2._endian + b'fi', self.size))

        nloads_actual = 0
        while (len(data) - n) >= ntotal:
            edata = data[n:n+ntotal]
            n += ntotal
            out = struct_i2fi.unpack(edata)
            (sid, s, si, l1) = out
            if op2.is_debug_file:
                op2.binary_debug.write('  LOAD=%s\n' % str(out))
            Si = [si]
            L1 = [l1]
            while 1:
                edata = data[n:n+ntotal2]
                n += ntotal2
                (si, l1) = struct_fi.unpack(edata)
                si_test, = op2.struct_i.unpack(edata[0:4])

                if [si_test, l1] == [-1, -1]:
                    break
                Si.append(si)
                L1.append(l1)
                if op2.is_debug_file:
                    op2.binary_debug.write('       [%s,%s]\n' % (si, l1))

            load = op2.add_load(sid, s, Si, L1)
            #op2._add_methods.add_load_combination_object(load)
            nloads_actual += 1
            #count += 1
            #if count > 1000:
                #raise RuntimeError('Iteration limit...probably have a bug.')
        op2.card_count['LOAD'] = nloads_actual
        return n

    def _read_loadcyh(self, data: bytes, n: int) -> int:
        """LOADCYH"""
        self.op2.log.info('geom skipping LOADCYH in GEOM3')
        return len(data)

    def _read_loadcyn(self, data: bytes, n: int) -> int:
        """LOADCYN"""
        self.op2.log.info('geom skipping LOADCYN in GEOM3')
        return len(data)

    def _read_loadcyt(self, data: bytes, n: int) -> int:
        """LOADCYT"""
        self.op2.log.info('geom skipping LOADCYT in GEOM3')
        return len(data)

    def _read_lseq(self, data: bytes, n: int) -> int:
        return len(data)
        op2 = self.op2
        ntotal = 20 * self.factor  # 5*4
        nloads = (len(data) - n) // ntotal

        #op2.show_data(data[n:], types='iq')
        n, ints = get_ints(data, n, nloads, 5,
                           size=op2.size, endian=op2._endian)
        lseq_id = ints[:, 0]
        darea_id = ints[:, 1]
        load_id = ints[:, 2]
        temperature_id = ints[:, 3]
        load = op2.lseq
        load._save(lseq_id, darea_id, load_id, temperature_id)
        #struct_5i = Struct(op2._endian + b'5i')
        #for unused_i in range(nloads):
            #out = struct_5i.unpack(data[n:n + ntotal])
            ##(sid, darea, load_id, temperature_id, undef) = out
            #if op2.is_debug_file:
                #op2.binary_debug.write('  LSEQ=%s\n' % str(out))
            #load = LSEQ.add_op2_data(out)
            #op2._add_methods.add_lseq_object(load)
            #n += ntotal
        op2.card_count['LSEQ'] = nloads
        return n

    def _read_moment(self, data: bytes, n: int) -> int:
        """
        MOMENT(4801,48,19) - the marker for Record 13
        """
        op2 = self.op2
        ntotal = 28 * self.factor # 7*4
        nloads = (len(data) - n) // ntotal
        n, ints, floats = get_ints_floats(data, n, nloads, 7,
                                          size=op2.size, endian=op2._endian)
        load_id = ints[:, 0]
        node_id = ints[:, 1]
        coord_id = ints[:, 2]
        magnitude = floats[:, 3]
        xyz = floats[:, [4, 5, 6]]
        load = op2.moment
        load._save(load_id, node_id, coord_id, magnitude, xyz)
        filter_large_node_ids(load)
        load.write()
        #s = Struct(op2._endian + b'3i4f')
        #for unused_i in range(nloads):
            #edata = data[n:n + 28]
            #out = s.unpack(edata)
            #if op2.is_debug_file:
                #op2.binary_debug.write('  MOMENT=%s\n' % str(out))
            ##(sid, g, cid, m, n1, n2, n3) = out
            #load = MOMENT.add_op2_data(out)
            #op2._add_methods.add_load_object(load)
            #n += 28
        op2.card_count['MOMENT'] = nloads
        return n

    def _read_moment1(self, data: bytes, n: int) -> int:
        """
        MOMENT1(4601,46,21) - the marker for Record 14
        """
        op2 = self.op2
        ntotal = 20 * self.factor  # 5*4
        nloads = (len(data) - n) // ntotal
        n, ints, floats = get_ints_floats(data, n, nloads, 5,
                                          size=op2.size, endian=op2._endian)
        load_id = ints[:, 0]
        node_id = ints[:, 1]
        magnitude = floats[:, 2]
        nodes = ints[:, [3, 4]]
        nloads = (len(data) - n) // ntotal
        load = op2.moment1
        load._save(load_id, node_id, magnitude, nodes)
        filter_large_node_ids(load)
        load.write()

        #s = Struct(op2._endian + b'iifii')
        #for unused_i in range(nloads):
            #edata = data[n:n + 20]
            #out = s.unpack(edata)
            #if op2.is_debug_file:
                #op2.binary_debug.write('  MOMENT1=%s\n' % str(out))
            ##(sid, g, m, n1, n2) = out
            #load = MOMENT1.add_op2_data(out)
            #op2._add_methods.add_load_object(load)
            #n += 20
        op2.card_count['MOMENT1'] = nloads
        return n

    def _read_moment2(self, data: bytes, n: int) -> int:
        """
        MOMENT2(4701,47,23) - the marker for Record 15
        """
        op2 = self.op2
        ntotal = 28 * self.factor  # 7*4
        nloads = (len(data) - n) // ntotal

        n, ints, floats = get_ints_floats(data, n, nloads, 7,
                                          size=op2.size, endian=op2._endian)
        load_id = ints[:, 0]
        node_id = ints[:, 1]
        magnitude = floats[:, 2]
        nodes = ints[:, [3, 4, 5, 6]]
        nloads = (len(data) - n) // ntotal
        load = op2.moment2
        load._save(load_id, node_id, magnitude, nodes)
        filter_large_node_ids(load)
        load.write()

        #s = Struct(op2._endian + b'iif4i')
        #for unused_i in range(nloads):
            #edata = data[n:n + 28]
            #out = s.unpack(edata)
            #if op2.is_debug_file:
                #op2.binary_debug.write('  MOMENT2=%s\n' % str(out))
            ##(sid, g, m, n1, n2, n3, n4) = out
            #load = MOMENT2.add_op2_data(out)
            #op2._add_methods.add_load_object(load)
            #n += 28
        op2.card_count['MOMENT2'] = nloads
        return n

    def _read_pload(self, data: bytes, n: int) -> int:
        """
        PLOAD(5101,51,24)
        """
        op2 = self.op2
        ntotal = 24 * self.factor  # 6*4

        nloads = (len(data) - n) // ntotal
        n, ints, floats = get_ints_floats(data, n, nloads, 6,
                                          size=op2.size, endian=op2._endian)
        load_id = ints[:, 0]
        pressure = floats[:, 1]
        nodes = ints[:, [2, 3, 4, 5]]
        nloads = (len(data) - n) // ntotal
        load = op2.pload
        load._save(load_id, pressure, nodes)
        load.write()
        if 0:  # pragma: no cover
            headers = ['$', 'sid', 'press', 'n1', 'n2',
                       'n3', 'n4']
            print('%-8s'*len(headers) % tuple(headers))
            print(load.write())

        #s = Struct(op2._endian + b'i f 4i')
        #for unused_i in range(nloads):
            #edata = data[n:n + 24]
            #out = s.unpack(edata)
            #if op2.is_debug_file:
                #op2.binary_debug.write('  PLOAD=%s\n' % str(out))
            ##(sid, pressure, n1, n2, n3, n4) = out
            #load = PLOAD.add_op2_data(out)
            #op2._add_methods.add_load_object(load)
            #n += 24
        op2.card_count['PLOAD1'] = nloads
        return n

    def _read_pload1(self, data: bytes, n: int) -> int:
        """
        PLOAD1(6909, 69, 198) - the marker for Record 17
        """
        op2 = self.op2
        ntotal = 32 * self.factor  # 8*4
        nloads = (len(data) - n) // ntotal

        n, ints, floats = get_ints_floats(data, n, nloads, 8,
                                          size=op2.size, endian=op2._endian)

        valid_load_types = ['FX', 'FY', 'FZ', 'FXE', 'FYE', 'FZE',
                       'MX', 'MY', 'MZ', 'MXE', 'MYE', 'MZE']

        # LE: length-based; FR: fractional; PR:projected
        valid_scales = ['LE', 'FR', 'LEPR', 'FRPR']

        load_id = ints[:, 0]
        element_id = ints[:, 1]
        load_type = ints[:, 2]
        scale = ints[:, 3]

        load_type_str = np.zeros(nloads, dtype='|U8')
        scale_str = np.zeros(nloads, dtype='|U8')

        for load_typei in np.unique(load_type):
            iload_type = np.where(load_typei == load_type)[0]
            load_type_stri = valid_load_types[load_typei-1]
            load_type_str[iload_type] = load_type_stri

        for scalei in np.unique(scale):
            iscale = np.where(scalei == scale)[0]
            scale_stri = valid_scales[scalei-1]
            scale_str[iscale] = scale_stri

        x = floats[:, [4, 6]]
        pressure = floats[:, [5, 7]]
        nloads = (len(data) - n) // ntotal
        load = op2.pload1
        load._save(load_id, element_id, load_type_str, scale_str, x, pressure)
        if 0:  # pragma: no cover
            headers = ['$', 'sid', 'eid', 'load_type', 'scale',
                       'x1', 'p1', 'x2', 'p2']
            print('%-8s'*len(headers) % tuple(headers))
            print(load.write())
        load.write()

        #s = Struct(mapfmt(op2._endian + b'4i4f', self.size))
        #for unused_i in range(nloads):
            #edata = data[n:n + ntotal]
            #out = s.unpack(edata)
            #if op2.is_debug_file:
                #op2.binary_debug.write('  PLOAD1=%s\n' % str(out))
            ##(sid, eid, load_type, scale, x1, p1, x2, p2) = out
            #load = PLOAD1.add_op2_data(out)
            #op2._add_methods.add_load_object(load)
            #n += ntotal
        op2.card_count['PLOAD1'] = nloads
        return n

    def _read_pload2(self, data: bytes, n: int) -> int:
        """
        PLOAD2(6802,68,199) - the marker for Record 18
        """
        op2 = self.op2
        ntotal = 12 * self.factor  # 3*4
        nloads = (len(data) - n) // ntotal

        n, ints, floats = get_ints_floats(data, n, nloads, 3,
                                          size=op2.size, endian=op2._endian)
        load_id = ints[:, 0]
        pressure = floats[:, 1]
        element_id = ints[:, 2]
        load = op2.pload2
        nelement = np.ones(nloads, dtype='int32')
        load._save(load_id, pressure, element_id, nelement)
        load.write()

        #struct_ifi = Struct(op2._endian + b'ifi')
        #for unused_i in range(nloads):
            #edata = data[n:n + 12]
            #out = struct_ifi.unpack(edata)
            #if op2.is_debug_file:
                #op2.binary_debug.write('  PLOAD2=%s\n' % str(out))
            ##(sid, p, eid) = out
            #load = PLOAD2.add_op2_data(out)
            #op2._add_methods.add_load_object(load)
            #n += 12
        op2.card_count['PLOAD2'] = nloads
        return n

    def _read_pload3(self, data: bytes, n: int) -> int:
        """PLOAD3(7109,71,255) - the marker for Record 19"""
        op2 = self.op2
        ntotal = 20  # 5*4
        nloads = (len(data) - n) // ntotal
        s = Struct(op2._endian + b'if3i')
        for unused_i in range(nloads):
            edata = data[n:n + 20]
            out = s.unpack(edata)
            if op2.is_debug_file:
                op2.binary_debug.write('  PLOAD3=%s\n' % str(out))
            #(sid, p, eid, n1, n2) = out
            load = PLOAD3.add_op2_data(out)  # undefined
            op2._add_methods.add_load_object(load)
            n += 20
        op2.card_count['PLOAD3'] = nloads
        return n

    def _read_pload4(self, data: bytes, n: int) -> int:
        """PLOAD4(7209,72,299) - the marker for Record 20"""
        n = self.op2.reader_geom2._read_dual_card(
            data, n,
            self._read_pload4_nx, self._read_pload4_msc,
            'PLOAD4', None)
        return n

    def _read_pload4_msc(self, data, n):  ## inconsistent with DMAP
        """
        PLOAD4(7209,72,299) - the marker for Record 20

        Word Name Type Description
        1 SID          I Load set identification number
        2 EID          I Element identification number
        3 P(4)        RS Pressures
        7 G1           I Grid point identification number at a corner of the face
        8 G34          I Grid point identification number at a diagonal from G1 or CTETRA corner
        9  CID         I Coordinate system identification number
        10 N(3)       RS Components of a vector coordinate system defined by CID
        13 SDRL(2) CHAR4 Load set on element SURF or LINE
        15 LDIR(2) CHAR4 Load direction
        """
        op2 = self.op2
        ntotal = 64 * self.factor # 16*4
        nloads = (len(data) - n) // ntotal
        assert (len(data) - n) % ntotal == 0
        #print('nloads =', nloads)
        loads = []
        #if self.size == 4:
            #s = Struct(op2._endian + b'2i 4f 3i 3f 8s 8s')
        #else:
            #s = Struct(op2._endian + b'2q 4d 3q 3d 16s 16s')

        op2.show_data(data[n:], types='ifs')
        n, ints, floats, strings = get_ints_floats_strings(
            data, n, nloads, 16,
            size=op2.size, endian=op2._endian)
        load_id = ints[:, 0]
        element_id = ints[:, 1]
        pressure = floats[:, [2, 3, 4, 5]]
        g1_34 = ints[:, [6, 7]]
        coord_id = ints[:, 8]
        nvector = floats[:, [9, 10, 11]]
        surf_or_line = np.core.defchararray.decode(
            np.core.defchararray.add(strings[:, 12], strings[:, 13]),
            'latin1')
        line_load_dir = strings[:, [12, 13]]

        assert element_id.min() > 0, element_id
        assert g1_34.min() >= 0, g1_34
        assert coord_id.min() >= 0, coord_id

        load = op2.pload4
        nelement = np.ones(nloads, dtype='int32')
        load._save(load_id, element_id, coord_id, pressure, g1_34,
                   surf_or_line=surf_or_line,
                   line_load_dir=line_load_dir,
                   nvector=nvector, nelement=nelement)
        load.write()

        #for unused_i in range(nloads):
            #edata = data[n:n + ntotal]
            #out = s.unpack(edata)
            #if op2.is_debug_file:
                #op2.binary_debug.write('  PLOAD4=%s\n' % str(out))
            #(sid, eid, p1, p2, p3, p4, g1, g34, cid, n1, n2, n3, surf_or_line, line_load_dir) = out

            #surf_or_line = surf_or_line.rstrip().decode('latin1')
            #line_load_dir = line_load_dir.rstrip().decode('latin1')
            #if line_load_dir == '':
                ## TODO: not 100%
                #line_load_dir = 'NORM'

            ## forces NX pload4 function to get called if it should be
            #assert surf_or_line in ['SURF', 'LINE']
            #assert line_load_dir in ['LINE', 'X', 'Y', 'Z', 'TANG', 'NORM'], 'line_load_dir=%r' % line_load_dir

            #load = PLOAD4.add_op2_data(
                #[sid, eid, [p1, p2, p3, p4], g1, g34,
                 #cid, [n1, n2, n3], surf_or_line, line_load_dir])
            #load.validate()
            #loads.append(load)
            #n += ntotal
        op2.card_count['PLOAD4'] = nloads
        return n, loads

    def _read_pload4_nx(self, data, n):  ## inconsistent with DMAP
        """
        PLOAD4(7209,72,299) - the marker for Record 20

        Word Name Type Description
        1 SID          I Load set identification number
        2 EID          I Element identification number
        3 P(4)        RS Pressures
        7 G1           I Grid point identification number at a corner of the face
        8 G34          I Grid point identification number at a diagonal from G1 or CTETRA corner
        9  CID         I Coordinate system identification number
        10 N(3)       RS Components of a vector coordinate system defined by CID
        """
        op2 = self.op2
        ntotal = 48 * self.factor  # 12*4
        nloads = (len(data) - n) // ntotal
        assert (len(data) - n) % ntotal == 0
        loads = []

        n, ints, floats = get_ints_floats(
            data, n, nloads, 12,
            size=op2.size, endian=op2._endian)
        load_id = ints[:, 0]
        element_id = ints[:, 1]
        pressure = floats[:, [2, 3, 4, 5]]
        g1_34 = ints[:, [6, 7]]
        coord_id = ints[:, 8]
        nvector = floats[:, [9, 10, 11]]
        load = op2.pload4
        assert element_id.min() > 0, element_id
        assert g1_34.min() >= 0, g1_34
        assert coord_id.min() >= 0, coord_id

        nelement = np.ones(nloads, dtype='int32')
        load._save(load_id, element_id, coord_id, pressure, g1_34,
                   line_load_dir=None, surf_or_line=None,
                   nvector=nvector, nelement=nelement)
        filter_large_element_ids(load, word='element_ids')
        #load.nelement = np.ones(load.n, dtype='int32')
        load.write()

        #s = Struct(mapfmt(op2._endian + b'2i 4f 3i 3f', self.size))
        #for unused_i in range(nloads):
            #edata = data[n:n + ntotal]
            #out = s.unpack(edata)
            #if op2.is_debug_file:
                #op2.binary_debug.write('  PLOAD4=%s\n' % str(out))
            #(sid, eid, p1, p2, p3, p4, g1, g34, cid, n1, n2, n3) = out

            #surf_or_line = None
            #line_load_dir = None
            #load = PLOAD4.add_op2_data(
                #[sid, eid, [p1, p2, p3, p4], g1, g34,
                 #cid, [n1, n2, n3], surf_or_line, line_load_dir])
            #load.validate()
            #loads.append(load)
            #n += ntotal
        #op2.card_count['PLOAD4'] = nloads
        return n, loads

    def _read_ploadx(self, data: bytes, n: int) -> int:
        """
        Record - PLOADX(7001,70,278)
        This record is obsolete

        Word Name Type Description
        1 SID   I Load set identification number
        2 P(2) RS Pressure
        4 G(3)  I Grid point identification numbers
        """
        op2 = self.op2
        ntotal = 24 * self.factor  # 6*4
        nloads = (len(data) - n) // ntotal
        assert (len(data) - n) % ntotal == 0
        struc = Struct(mapfmt(op2._endian + b'iff3i', self.size))
        for unused_i in range(nloads):
            edata = data[n:n + ntotal]
            out = struc.unpack(edata)
            if op2.is_debug_file:
                op2.binary_debug.write('  PLOADX=%s\n' % str(out))
            fields = ['PLOADX'] + list(out)
            op2.reject_lines.append(print_card_16(fields))
            #load = PLOADX.add_op2_data(out)
            #op2._add_methods.add_load_object(load)

            #n += ntotal + 28 * self.factor
            n += ntotal
        op2.card_count['PLOADX'] = nloads
        return n

    def _read_ploadx1(self, data: bytes, n: int) -> int:
        """
        Record - PLOADX1(7309,73,351)

        Word Name Type Description
        1 SID    I Load set identification number
        2 EID    I Element identification number
        3 PA    RS Surface traction at grid point GA
        4 PB    RS Surface traction at grid point GB
        5 G(2)   I Corner grid point identification numbers
        7 THETA RS Angle between surface traction and inward normal
        """
        op2 = self.op2
        if 'PLOADX1' not in op2.cards_to_read:
            return len(data)
        ntotal = 28 * self.factor  # 7*4
        nloads = (len(data) - n) // ntotal

        n, ints, floats = get_ints_floats(
            data, n, nloads, 7,
            size=op2.size, endian=op2._endian)
        load_id = ints[:, 0]
        element_id = ints[:, 1]
        surface_traction = floats[:, [2, 3]]
        nodes = ints[:, [4, 5]]
        theta = floats[:, 6]
        load = op2.ploadx1
        load._save(load_id, element_id, nodes, surface_traction, theta)
        filter_large_element_ids(load, word='element_id')
        load.write()

        #struc = Struct(op2._endian + b'2i2f iif')
        #for unused_i in range(nloads):
            #edata = data[n:n + 28]
            #out = struc.unpack(edata)
            #if op2.is_debug_file:
                #op2.binary_debug.write('  PLOADX1=%s\n' % str(out))
            #load = PLOADX1.add_op2_data(out)
            #op2._add_methods.add_load_object(load)
            #n += 28
        op2.card_count['PLOADX1'] = nloads
        return n

# PRESAX

    def _read_qbdy1(self, data: bytes, n: int) -> int:
        """
        QBDY1(4509,45,239) - the marker for Record 24
        """
        op2 = self.op2
        ntotal = 12 * self.factor  # 3*4
        nloads = (len(data) - n) // ntotal

        n, ints, floats = get_ints_floats(
            data, n, nloads, 3,
            size=op2.size, endian=op2._endian)
        load_id = ints[:, 0]
        q0 = floats[:, 1]
        element_id = ints[:, 2]
        nelement = np.ones(nloads, dtype='int32')
        load = op2.qbdy1
        load._save(load_id, q0, element_id, nelement)
        load.write()

        #struct_ifi = Struct(op2._endian + b'ifi')
        #for unused_i in range(nloads):
            #edata = data[n:n + 12]
            #out = struct_ifi.unpack(edata)
            #if op2.is_debug_file:
                #op2.binary_debug.write('  QBDY1=%s\n' % str(out))
            ##(sid, q0, eid) = out
            #load = QBDY1.add_op2_data(out)
            #op2._add_methods.add_thermal_load_object(load)
            #n += 12
        op2.card_count['QBDY1'] = nloads
        return n

    def _read_qbdy2(self, data: bytes, n: int) -> int:
        """
        QBDY2(4909,49,240) - the marker for Record 25
        """
        op2 = self.op2
        ntotal = 40 * self.factor  # 10*4
        nloads = (len(data) - n) // ntotal

        n, ints, floats = get_ints_floats(
            data, n, nloads, 10,
            size=op2.size, endian=op2._endian)
        load_id = ints[:, 0]
        element_id = ints[:, 1]
        q0 = floats[:, [2, 3, 4, 5,
                        6, 7, 8, 9]].flatten()
        nflux = np.full(nloads, 8, dtype='int32')
        load = op2.qbdy2
        load._save(load_id, element_id, q0, nflux)
        load.write()

        #struct_2i8f = Struct(op2._endian + b'ii8f')
        #for unused_i in range(nloads):
            #edata = data[n:n + 40]
            #out = struct_2i8f.unpack(edata)
            #if op2.is_debug_file:
                #op2.binary_debug.write('  QBDY2=%s\n' % str(out))
            ##(sid, eid, q1, q2, q3, q4, q5, q6, q7, q8) = out
            #load = QBDY2.add_op2_data(out)
            #op2._add_methods.add_thermal_load_object(load)
            #n += 40
        op2.card_count['QBDY2'] = nloads
        return n

    def _read_qbdy3(self, data: bytes, n: int) -> int:
        """
        QBDY3(2109,21,414) - the marker for Record 26
        """
        op2 = self.op2
        ntotal = 16 * self.factor  # 4*4
        nloads = (len(data) - n) // ntotal

        n, ints, floats = get_ints_floats(
            data, n, nloads, 4,
            size=op2.size, endian=op2._endian)
        load_id = ints[:, 0]
        q0 = floats[:, 1]
        control_node = ints[:, 2]
        element_id = ints[:, 3]
        load = op2.qbdy3
        nelement = np.ones(nloads, dtype='int32')
        load._save(load_id, q0, control_node, element_id, nelement)
        load.write()

        #struct_if2i = Struct(op2._endian + b'ifii')
        #for unused_i in range(nloads):
            #edata = data[n:n + 16]
            #out = struct_if2i.unpack(edata)
            ##(sid, q0, control_node, eid) = out
            #load = QBDY3.add_op2_data(out)
            #op2._add_methods.add_thermal_load_object(load)
            #n += 16
        op2.card_count['QBDY3'] = nloads
        return n

    def _read_temp(self, data: bytes, n: int) -> int:
        """
        TEMP(5701,57,27) - the marker for Record 32
        .. warning:: buggy
        """
        op2 = self.op2
        ntotal = 12 * self.factor # 3*4
        nloads = (len(data) - n) // ntotal

        n, ints, floats = get_ints_floats(
            data, n, nloads, 3,
            size=op2.size, endian=op2._endian)
        load_id = ints[:, 0]
        node_id = ints[:, 1]
        temperature = floats[:, 2]
        load = op2.temp
        nnode = np.arange(nloads, dtype='int32')
        load._save(load_id, node_id, temperature, nnode)
        load.write()

        #struct_2if = Struct(mapfmt(op2._endian + b'iif', self.size))
        #for unused_i in range(nloads):
            #edata = data[n:n + ntotal]
            #out = struct_2if.unpack(edata)
            #if op2.is_debug_file:
                #op2.binary_debug.write('  TEMP=%s\n' % str(out))
            #(sid, g, T) = out
            #if g < 10000000:
                #load = TEMP.add_op2_data(out)
                #op2._add_methods.add_thermal_load_object(load)
            #else:
                #op2.log.debug('TEMP = %s' % str(out))
            #n += ntotal
        op2.card_count['TEMP'] = nloads
        return n

    def _read_tempd(self, data: bytes, n: int) -> int:
        """
        TEMPD(5641,65,98) - the marker for Record 33
        .. todo:: add object
        """
        return len(data)
        op2 = self.op2
        ntotal = 8 * self.factor  # 2*4
        nloads = (len(data) - n) // ntotal

        n, ints, floats = get_ints_floats(
            data, n, nloads, 2,
            size=op2.size, endian=op2._endian)
        load_id = ints[:, 0]
        temperature = floats[:, 1]
        load = op2.tempd
        load._save(load_id, temperature)
        load.write()

        #struct_if = Struct(mapfmt(op2._endian + b'if', self.size))
        #for unused_i in range(nloads):
            #edata = data[n:n + ntotal]
            #out = struct_if.unpack(edata)
            #if op2.is_debug_file:
                #op2.binary_debug.write('  TEMPD=%s\n' % str(out))
            ##(sid, T) = out
            #load = TEMPD.add_op2_data(out)
            #op2._add_methods.add_tempd_object(load)
            #n += ntotal
        op2.card_count['TEMPD'] = nloads
        return n

    def _read_qhbdy(self, data: bytes, n: int) -> int:
        """
        (4309,43,233)

        Word Name Type Description
        1 SID   I Load set identification number
        2 FLAG  I Face type
        3 Q0   RS Magnitude of thermal flux into face
        4 AF   RS Area factor
        5 G(8)  I Grid point identification numbers

        """
        op2 = self.op2
        ntotal = 48 * self.factor  # 12*4
        nloads = (len(data) - n) // ntotal

        n, ints, floats = get_ints_floats(
            data, n, nloads, 12,
            size=op2.size, endian=op2._endian)
        load_id = ints[:, 0]
        flag = ints[:, 1]
        q0 = floats[:, 2]
        area_factor = floats[:, 3]
        nodes = ints[:, [4, 5, 6, 7,
                         8, 9, 10, 11]].flatten()
        nnodes = np.ones(nloads, dtype='int32') * 8

        flag_int_to_str = {
            1 : 'POINT',
            2: 'LINE',
            3: 'REV',
            5: 'AREA4',
            9: 'AREA8',
        }
        #if flag_int == 1:
            #flag = 'POINT'
            #nnodes = 1
        #elif flag_int == 2:
            #flag = 'LINE'
            #nnodes = 2
        #elif flag_int == 3:
            #flag = 'REV'
            #nnodes = 2
        #elif flag_int == 5:
            #flag = 'AREA4'
            #nnodes = 4
        #elif flag_int == 9:
            #flag = 'AREA8'
            #nnodes = 8
        #else:  # pragma: no cover
            #raise NotImplementedError(f'QHBDY sid={sid} flag_int={flag_int} data={data[2:]}')

        flag_str = np.zeros(nloads, dtype='|U8')
        for flag_inti in np.unique(flag):
            i = np.where(flag_inti == flag)[0]
            flag_stri = flag_int_to_str[flag_inti]
            flag_str[i] = flag_stri
        load = op2.qhbdy
        load._save(load_id, q0, flag_str, area_factor, nodes, nnodes)
        load.write()

        #struct_if = Struct(op2._endian + b'2i2f 8i')
        #for unused_i in range(nloads):
            #edata = data[n:n + ntotal]
            #out = struct_if.unpack(edata)
            #if op2.is_debug_file:
                #op2.binary_debug.write('  QHBDY=%s\n' % str(out))
            ##(sid, flag, q0, af, n1-n8) = out
            #load = QHBDY.add_op2_data(out)
            ##self.add_thermal_load(load)
            #op2._add_methods.add_load_object(load)
            #n += ntotal
        op2.card_count['QHBDY'] = nloads
        return n

    def _read_qvect(self, data: bytes, n: int) -> int:
        """
        MSC
        Record 29 -- QVECT(2209,22,241)

        Word Name Type Description
        1 SID   I  Load set identification number
        2 Q0    RS Magnitude of thermal flux vector into face
        3 TSOUR RS Temperature of the radiant source
        4 CE    I  Coordinate system identification number for thermal vector flux
        5 FLAG  I
        6 E     RS Vector component of flux in coordinate system CE
        Words 5 through 6 repeat 3 times
        5b, 6b
        5c, 6c
        7 CONTROL_NODE I Control point
        8 EID     I Element identification number
        ints    = (200, 442.0, 10400.0, 0,   0,   0,   0,   0,   0,   -1.0, 0,   10)
        floats  = (200, 442.0, 10400.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.0, 10)
        """
        op2 = self.op2
        ntotal = 48 * self.factor  # 12*4
        ndatai = len(data) - n
        nloads = ndatai // ntotal
        assert nloads > 0, nloads
        assert ndatai % ntotal == 0, ndatai

        n, ints, floats = get_ints_floats(
            data, n, nloads, 12,
            size=op2.size, endian=op2._endian)
        load_id = ints[:, 0]
        q0 = floats[:, 1]
        t_source = floats[:, 2]

        ce = ints[:, 3]

        #e1 = floats[:, 5]
        #e2 = floats[:, 7]
        #e3 = floats[:, 9]
        #flag1 = ints[:, 4]
        #flag2 = ints[:, 6]
        #flag3 = ints[:, 8]

        control_node = ints[:, 10]
        element_id = ints[:, 11]
        nelement = np.ones(nloads, dtype='int32')
        load = op2.qvect
        #load_id, q0, t_source, control_node, ce, vector_tableds,
                      #element, nelement

        vector = floats[:, [5, 7, 9]]
        tableds = ints[:, [4, 6, 8]]
        #print(vector)
        #print(tableds)
        load._save(load_id, q0, t_source, control_node, ce,
                   vector, tableds,
                   element_id, nelement)
        load.write()
                                             # ce
        #struct_if = Struct(op2._endian + b'i 2f i if if if 2i')
        #for unused_i in range(nloads):
            #edata = data[n:n + ntotal]
            #out = struct_if.unpack(edata)
            #if op2.is_debug_file:
                #op2.binary_debug.write('  QHBDY=%s\n' % str(out))
            #(sid, q0, t_source, ce, flag1, e1, flag2, e2, flag3, e3, control_node, eid) = out

            #load = QVECT.add_op2_data(out)
            #str(load)
            ##self.add_thermal_load(load)
            #op2._add_methods.add_load_object(load)
            #n += ntotal
        op2.card_count['QHBDY'] = nloads
        return n

    def _read_qvol(self, data: bytes, n: int) -> int:
        """
        Record 30 -- QVOL(2309,23,416)

        1 SID      I Load set identification number
        2 QVOL    RS Power input per unit volume produced by a
                     conduction element
        3 CONTROL_NODE  I Control point used for controlling heat generation
        4 EID      I Element identification number

        """
        op2 = self.op2
        ntotal = 16 * self.factor  # 4*4
        nloads = (len(data) - n) // ntotal

        n, ints, floats = get_ints_floats(
            data, n, nloads, 4,
            size=op2.size, endian=op2._endian)
        load_id = ints[:, 0]
        qvol = floats[:, 1]
        control_node = ints[:, 2]
        element_id = ints[:, 2]
        nelement = np.ones(nloads, dtype='int32')
        load = op2.qvol
        load._save(load_id, qvol, control_node, element_id, nelement)
        load.write()

        #struc = Struct(op2._endian + b'if2i')
        #for unused_i in range(nloads):
            #edata = data[n:n + ntotal]
            #out = struc.unpack(edata)
            #if op2.is_debug_file:
                #op2.binary_debug.write('  QVOL=%s\n' % str(out))
            ##(sid, qvol, control_node, eid) = out
            #load = QVOL.add_op2_data(out)
            #op2._add_methods.add_load_object(load)
            #n += ntotal
        op2.card_count['QVOL'] = nloads
        return n

    def _read_rforce(self, data: bytes, n: int) -> int:
        op2 = self.op2
        ntotal = 40 * self.factor  # 10*4
        nloads = (len(data) - n) // ntotal

        n, ints, floats = get_ints_floats(
            data, n, nloads, 10,
            size=op2.size, endian=op2._endian)
        load_id = ints[:, 0]
        node_id = ints[:, 1]
        coord_id = ints[:, 2]
        scale = floats[:, 3]
        r123 = floats[:, [4, 5, 6]]
        method = ints[:, 7]
        racc = floats[:, 8]
        main_bulk = ints[:, 9]
        idrf = np.zeros(nloads, dtype='int32')
        load = op2.rforce
        load._save(load_id, node_id, coord_id, scale, r123, method, racc, main_bulk, idrf)
        load.write()

        #struc = Struct(mapfmt(op2._endian + b'3i 4f ifi', self.size))
        #for unused_i in range(nloads):
            #edata = data[n:n + ntotal]
            #out = struc.unpack(edata)
            #if op2.is_debug_file:
                #op2.binary_debug.write('  RFORCE=%s\n' % str(out))
            ##(sid, nid, cid, a, r1, r2, r3, method, racc, mb) = out
            #load = RFORCE.add_op2_data(out)
            #op2._add_methods.add_load_object(load)
            #n += ntotal
        op2.card_count['RFORCE'] = nloads
        return n

    def _read_sload(self, data: bytes, n: int) -> int:
        """SLOAD(5401, 54, 25)"""
        op2 = self.op2
        ntotal = 12 * self.factor  # 3*4
        nloads = (len(data) - n) // ntotal

        n, ints, floats = get_ints_floats(
            data, n, nloads, 3,
            size=op2.size, endian=op2._endian)
        load_id = ints[:, 0]
        node_id = ints[:, 1]
        scale_factor = floats[:, 2]
        load = op2.sload
        load._save(load_id, node_id, scale_factor)
        filter_large_node_ids(load, word='nodes')
        load.write()

        #struc = Struct(op2._endian + b'2i f')
        #for unused_i in range(nloads):
            #edata = data[n:n + ntotal]
            #out = struc.unpack(edata)
            #if op2.is_debug_file:
                #op2.binary_debug.write('  SLOAD=%s\n' % str(out))
            ##(sid, nid, scale_factor) = out
            #load = SLOAD.add_op2_data(out)
            #op2._add_methods.add_load_object(load)
            #n += ntotal
        op2.card_count['SLOAD'] = nloads
        return n

# TEMP(5701,57,27) # 32
# TEMPD(5641,65,98) # 33
# TEMPEST
    def _read_tempf(self, data: bytes, n: int) -> int:
        self.op2.log.info('geom skipping TEMPF in GEOM3')
        return len(data)
# TEMP1C

    def _read_tempp1(self, data: bytes, n: int) -> int:
        op2 = self.op2
        op2.log.warning('geom skipping TEMPP2 in GEOM3')
        ntotal = 24 * self.factor  # 6*4
        nloads = (len(data) - n) // ntotal
        n, ints, floats = get_ints_floats(
            data, n, nloads, 6,
            size=op2.size, endian=op2._endian)
        return len(data)
        struc = Struct(op2._endian + b'2i 4f')
        for unused_i in range(nloads):
            edata = data[n:n + ntotal]
            out = struc.unpack(edata)
            if op2.is_debug_file:
                op2.binary_debug.write('  TEMPP1=%s\n' % str(out))
            #sid, eid, t, tprime, ts1, ts2 = data
            load = TEMPP1.add_op2_data(out)
            op2._add_methods.add_load_object(load)
            n += ntotal
        op2.card_count['TEMPP1'] = nloads
        return n

    def _read_tempp2(self, data: bytes, n: int) -> int:
        op2 = self.op2
        op2.log.info('geom skipping TEMPP2 in GEOM3')
        if op2.is_debug_file:
            op2.binary_debug.write('geom skipping TEMPP2 in GEOM3\n')
        return len(data)

    def _read_tempp3(self, data: bytes, n: int) -> int:
        op2 = self.op2
        op2.log.info('geom skipping TEMPP3 in GEOM3')
        if op2.is_debug_file:
            op2.binary_debug.write('geom skipping TEMPP3 in GEOM3\n')
        return len(data)

    def _read_tempp4(self, data: bytes, n: int) -> int:
        """
        TEMPP4(4201,42,18) - the marker for Record 40
        """
        self.op2.log.info('geom skipping TEMPP4 in GEOM3')
        return len(data)

    def _read_temprb(self, data: bytes, n: int) -> int:
        op2 = self.op2
        op2.log.info('geom skipping TEMPRB in GEOM3')
        if op2.is_debug_file:
            op2.binary_debug.write('geom skipping TEMPRB in GEOM3\n')
        return len(data)

    def _read_pface(self, data: bytes, n: int) -> int:
        self.op2.log.info('geom skipping PFACE in GEOM3')
        return len(data)

    def _read_pedge(self, data: bytes, n: int) -> int:
        self.op2.log.info('geom skipping PEDGE in GEOM3')
        return len(data)

    def _read_boltfor(self, data: bytes, n: int) -> int:
        self.op2.log.info('geom skipping BOLTFOR in GEOM3')
        return len(data)

    def _read_boltld(self, data: bytes, n: int) -> int:
        self.op2.log.info('geom skipping BOLTLD in GEOM3')
        return len(data)

def filter_large_node_ids(load, word='node_id'):
    node_id = getattr(load, word)
    if node_id.max() < NID_CAP:
        return
    ibig = np.where(node_id >= NID_CAP)[0]
    i = np.where(node_id < EID_CAP)[0]
    removed_nodes = node_id[ibig]
    load.model.log.warning(f'{load.type}: removing node ids > {NID_CAP}\n'
                           f'{word}={removed_nodes}\n')
    load.__apply_slice__(load, i)

def filter_large_element_ids(load, word='element_id'):
    element_id = getattr(load, word)
    if element_id.max() < EID_CAP:
        return
    ibig = np.where(element_id >= EID_CAP)[0]
    i = np.where(element_id < EID_CAP)[0]
    removed_elements = element_id[ibig]
    load.model.log.warning(f'{load.type}: removing element ids > {EID_CAP}\n'
                           f'element_ids={removed_elements}\n')
    load.__apply_slice__(load, i)
