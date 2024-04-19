"""
defines readers for BDF objects in the OP2 CONTACT/CONTACTS table
"""
from __future__ import annotations
from struct import Struct
from typing import TYPE_CHECKING

import numpy as np
#from pyNastran.op2.op2_interface.op2_reader import mapfmt, reshape_bytes_block

if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.op2.op2_geom import OP2Geom


class CONTACT:
    """defines methods for reading contact geometry"""

    @property
    def size(self) -> int:
        return self.op2.size
    @property
    def factor(self) -> int:
        return self.op2.factor

    def _read_fake(self, data: bytes, n: int) -> int:
        return self.op2._read_fake(data, n)

    def read_stop(self, data: bytes, n: int) -> int:
        return self.op2.reader_geom1.read_stop(data, n)

    def read_contact_4(self, data: bytes, ndata: int):
        """
        reads the CONTACT/CONTACTS table
        Table of Bulk Data entry related to surface contact

        """
        return self.op2._read_geom_4(self.contact_map, data, ndata)

    def __init__(self, op2: OP2Geom):
        self.op2 = op2

        # F:\Program Files\Siemens\NXNastran\nxn10p1\nxn10p1\nast\tpl\fsw_eng.op2
        # F:\work\pyNastran\pyNastran\master2\pyNastran\bdf\test\nx_spike\out_boltld04i.op2
        # F:\work\pyNastran\pyNastran\master2\pyNastran\bdf\test\nx_spike\out_eliter17.op2
        # F:\work\pyNastran\pyNastran\master2\pyNastran\bdf\test\nx_spike\out_weld01i.op2
        # F:\work\pyNastran\examples\Dropbox\move_tpl\ac10901a_new.op2

        # F:\work\pyNastran\pyNastran\master2\pyNastran\bdf\test\nx_spike\out_boltsold11b.op2
        # F:\work\pyNastran\pyNastran\master2\pyNastran\bdf\test\nx_spike\out_conedg01b.op2
        # F:\work\pyNastran\pyNastran\master2\pyNastran\bdf\test\nx_spike\out_conprop06.op2
        # F:\work\pyNastran\pyNastran\master2\pyNastran\bdf\test\nx_spike\out_glueac103a.op2
        # F:\work\pyNastran\pyNastran\master2\pyNastran\bdf\test\nx_spike\out_conedg01s.op2
        # F:\work\pyNastran\pyNastran\master2\pyNastran\bdf\test\nx_spike\out_sline5.op2
        self.contact_map = {
            # C:\MSC.Software\simcenter_nastran_2019.2\tpl_post1\femao8rand.op2
            (7110, 71, 588) : ['BSURFS', self._read_bsurfs],
            (724, 7, 441) : ['BSURF', self._read_fake],
            (224, 2, 436) : ['BLSEG', self._read_fake],
            (1224, 12, 446) : ['BGSET', self._read_fake],
            (7210, 72, 589) : ['BCPROP', self._read_fake],
            (7410, 74, 591) : ['BCTSET', self._read_fake],
            (7510, 75, 592) : ['BCTADD', self._read_fake],
            (8810, 88, 603) : ['BGADD', self._read_fake],
            (8920, 89, 614) : ['BEDGE', self._read_fake],
            (124, 1, 435) : ['BCONP', self._read_fake],
            (7710, 77, 594) : ['BCRPARA', self._read_fake],
            (8110, 81, 598) : ['BCTPARM', self._read_bctparm],
            (8301, 83, 605) : ['BCPROPS', self._read_fake],

            #(124, 1, 435) : ['???', self._read_fake],

            # Record – ACTRAD(5907,60,654)
            # Record – AMLREG(811,8,628)  .
            # Record – ATVFS(6571,65,657)  .
            # Record – BCMATL (7310,73,590)
            # Record – BFLUID(9001,90,964)
            # Record – BCTPAR2(6621,66,662)
            # Record – BGPARA (7610,76,593)
            # Record – CSMADD(6700,67,670)
            # Record – CSMSET(6590,62,659)
            # Record – EBDADD (8610,86,448)
            # Record – EBDSET (8510,85,447)
            # Record – FLXADD (9201,92,694)
            # Record – FLXSLI(9101,91,693)
            # Record – IPLANE(911,9,629)  .
            # Record – NXSTRAT (7810,78,595)
            # Record – PACTRAD(6581,61,658)
            # Record – TMCPARA (7910,79,989)
            # Record – VATVFS(6801,68,680)
            (8710, 87, 449) : ['???', self._read_fake],
            (424, 4, 438) : ['???', self._read_fake],

            (624, 6, 440) : ['???', self._read_fake],
            (1124, 11, 445) : ['BCPARA', self._read_fake],
            (20029, 29, 493) : ['BCPROP', self._read_fake],
            (1024, 10, 444) : ['BCBODY', self._read_fake],
            (811, 8, 628) : ['AMLREG', self._read_amlreg],

            (6621, 66, 662) : ['BCTPAR2', self._read_fake],
            (7610, 76, 593) : ['BGPARA', self._read_bgpara],
            (9101, 91, 693) : ['FLXSLI', self._read_fake],
            (5524, 55, 897) : ['BCONPRG', self._read_fake],
            (924, 9, 443) : ['???', self._read_fake],

            (6571, 65, 657) : ['???', self._read_fake],
            (4624, 46, 888) : ['BCONPRP', self._read_fake],
            (4524, 45, 887) : ['BCONECT', self._read_fake],
            (4424, 44, 886) : ['BCTABL1', self._read_fake],
            (7124, 71, 992) : ['BCAUTOP', self._read_fake],
            (4724, 47, 889) : ['BCBODY1', self._read_fake],
            (4824, 48, 890) : ['BCBDPRP', self._read_fake],
            (6724, 67, 948) : ['BCSCAP', self._read_fake],
            (4924, 49, 891) : ['???', self._read_fake],
            (5024, 50, 892) : ['BCRGSRF', self._read_fake],
            (5424, 54, 896) : ['BCNURBS', self._read_fake],
            (5724, 57, 899) : ['BCTRIM', self._read_fake],
            (7324, 73, 996) : ['BCSURF', self._read_fake],
            (7424, 74, 997) : ['BSURF', self._read_fake],
            (3824, 38, 742) : ['BCBMRAD', self._read_fake],
            (7224, 72, 995) : ['BCGRID', self._read_fake],
            (20032, 32, 496) : ['GMNURB', self._read_fake],

            (5224, 52, 894) : ['BCBZIER', self._read_fake],
            (5324, 53, 895) : ['BCNURB2', self._read_fake],
            (5124, 51, 893) : ['BCPATCH', self._read_fake],
        }

    def _read_amlreg(self, data: bytes, n: int) -> int:
        """
        Record – AMLREG(811,8,628)

        Word Name Type Description
        1 RID          I AML region ID
        2 SID          I BSURFS ID for surface definition
        3 DESC(12) CHAR4 Description - 48 character maximum
        15 NLAY        I Number of layers
        16 RADTYPE     I Radiation surface type:
                           0=None
                           1=AML
                           2=Physical boundary
        17 INFID1      I Infinite plane-1 ID
        18 INFID2      I Infinite plane-2 ID
        19 INFID3      I Infinite plane-3 ID
        20 UNDEF(3) None

        """
        op2: OP2Geom = self.op2
        size = self.size

        ntotal = 22 * size
        nentries = (len(data) - n) // ntotal
        assert (len(data) - n) % ntotal == 0, 'AMLREG'

        #loads = []
        if size == 4:
            struc = Struct(op2._endian + b'2i 48s 5i 3i')
        else:
            raise RuntimeError('size=8 AMLREG')
        for unused_i in range(nentries):
            edata = data[n:n+ntotal]
            out = struc.unpack(edata)
            (rid, sid, name_bytes, nlayers, radtype,
             infid1, infid2, infid3,
             dummy1, dummy2, dummy3) = out
            assert rid > 0, rid
            assert sid > 0, sid

            if radtype == 0:
                radsurf = 'NONE'
            elif radtype == 1:
                radsurf = 'AML'
            elif radtype == 2:
                radsurf = 'PHYB'
            else:
                raise NotImplementedError(radtype)

            name = name_bytes.decode('utf8').rstrip()
            assert (dummy1, dummy2, dummy3) == (0, 0, 0), (dummy1, dummy2, dummy3)

            infid = [infid1, infid2, infid3]
            # AMLREG RID SID     Name/Descriptor
            #        NL  RADSURF INFID1 INFID2 INFID3
            op2.add_amlreg(rid, sid, name,
                           infid,
                           nlayers=nlayers, radsurf=radsurf)
            #acsrce = ACSRCE(sid, excite_id, rho, b, delay=delay, dphase=dphase, power=0)
            n += ntotal
            #loads.append(acsrce)
        return n # , loads


    #def _read_bcrpara(self, data: bytes, n: int) -> int:
        #return self._read_bc_param(data, n, 'BCRPARA')

    def _read_bgpara(self, data: bytes, n: int) -> int:
        return self._read_bc_param(data, n, 'BGPARA')

    def _read_bctparm(self, data: bytes, n: int) -> int:
        return self._read_bc_param(data, n, 'BCTPARM')

    def _read_bc_param(self, data: bytes, n: int, name: str) -> int:
        """
        Record – BCTPARM (8110,81,598)

        Word Name Type Description
        1 CSID           I Contact set ID
        2-3 Param(i) CHAR4 Parameter name
        4 TYPE           I Parameter data type
        5 Value(i)    I/RS Parameter value (See Note 1 below)
        6 -1             I Delimiter

        Note: entry 2-5 repeats for each parameter

        Note
        ----
        1. The parameter NCHG can now be entered as a real or an integer value, but is always written to
        the OP2 file as a real value. For example, the integer value NCHG=1 is written to the OP2 file
        as NCHG=1.0.

        ndata = 24:
        strings = (b'e\x00\x00\x00INIPENE \x01\x00\x00\x00\x00\x00\x00\x00\xff\xff\xff\xff',)
        ints    = (101, 'INIPENE', 1, 0,   -1)
        floats  = (101, 'INIPENE', 1, 0.0, -1)
        """
        op2 = self.op2
        size = self.size
        #assert size == 4, f'{name} size={size} is not supported'

        #self.show_data(data[n:], types='ifs')
        ints = np.frombuffer(data[n:], op2.idtype8)
        floats = np.frombuffer(data[n:], op2.fdtype8)
        iminus1 = np.where(ints == -1)[0]

        istart = [0] + list(iminus1[:-1] + 1)
        iend = iminus1
        #istart = [0] + list(iend + 1)

        #print(istart)
        #print(iend)
        for i0, i1 in zip(istart, iend):
            strings = data[n+i0*size:n+i1*size]
            intsi = ints[i0:i1]
            floatsi = floats[i0:i1]
            contact_id = intsi[0]
            #print(contact_id, strings)
            j = 4
            nj = (len(intsi) - 1) // 4
            jn = 1
            j = size
            # ints    = (1,
            #            PENN, 2, 0.01,
            #            PENT, 2, 1008981770,
            #            CTOL, 2, 970045207,
            #            MAXS, 1, 50,
            #            MAXF, 1, 30,
            #            PENSCAL, 1, 1,
            #            1229342034, 538985806, 1, 2)
            # floats  = (PENN, 2, 0.01,
            #            PENT, 2, 0.01,
            #            CTOL, 2, 0.0004,
            #            MAXS, 1, 50,
            #            MAXF, 1, 30,
            #            PENSCAL, 1, 1,
            #            REFINE,  1, 2)
            #self.show_data(strings)
            #PENN 2 0.01
            #PENT 2 0.01
            #CTOL 2 0.0004
            #MAXS 1 7e-44
            #MAXF 1 4.2e-44
            #PENSCAL 1 1e-45
            #REFINE 1 3e-45
            #print(strings)
            #print(floatsi.tolist())
            #print(contact_id)
            params = {}
            for ji in range(nj):
                #assert len(intsi) == 5, f'nints={len(intsi)} ints={intsi}'
                word = strings[j:j+size*2].decode('latin1').rstrip()
                value_type = intsi[jn + 2]
                if value_type == 1:
                    value = intsi[jn + 3]
                elif value_type == 2:
                    value = floatsi[jn + 3]
                else:
                    raise NotImplementedError(f'value_type = {value_type}')
                #print('  ', word, value_type, value)
                params[word] = value
                j += size * 4
                jn += 4

            #1 CSID           I Contact set ID
            # 2-3 Param(i) CHAR4 Parameter name
            # 4 TYPE           I Parameter data type
            # 5 Value(i)    I/RS Parameter value (See Note 1 below)
            #self.add_bctparm(contact_id, params)
            #print(contact_id, params)
        op2.log.warning(f'skipping {name}; id={contact_id} params={params}')
        return len(data)

    def _read_bsurfs(self, data: bytes, n: int) -> int:
        """
        BSURFS

        (7110, 71, 588,

        1,
        24, 190, 198, 189,
        44, 188, 197, 190,
        64, 106, 189, 196,
        84, 195, 188, 106, -1)
                BSURFS         1                              24     190     198     189+
        $           EID2      G1      G2      G3    EID3      G1      G2      G3
        +             44     188     197     190      64     106     189     196+
        $           EID4      G1      G2      G3
        +             84     195     188     106

          ints    = (7110, 71, 588,

          1,
          1, 97, 1716, 1706,
          36, 358, 1713, 359,
          38, 97, 5, 1715,
          43, 357, 1717, 1713, 93, 1753, 1763, 1754, 110, 1696, 1697, 1705, 115, 346, 1705, 347, 117, 6, 1774, 1772, 118, 346, 347, 1768, 122, 333, 1774, 6, 124, 6, 1772, 130, 132, 1722, 1729, 1728, 212, 342, 343, 1766, 215, 92, 1686, 1682, 218, 342, 1703, 343, 220, 1731, 1737, 1736, 230, 1678, 1684, 1685, 336, 333, 1700, 1702, 389, 359, 1713, 1714, 436, 336, 1702, 1701, 444, 333, 1702, 334, 447, 362, 1715, 363, 454, 334, 1702, 335, 467, 1696, 1708, 1697, 475, 1765, 1770, 1766, 486, 1679, 1681, 1680, 688, 1748, 1749, 1755, 776, 1708, 1711, 1712, 777, 1757, 1758, 1759, 830, 130, 1772, 1771, 837, 1749, 1750, 1756, 838, 1750, 1771, 1765, 874, 359, 360, 1758, 879, 361, 362, 1760, 882, 361, 1715, 362, 920, 360, 361, 1760, 921, 1692, 1694, 1693, 922, 5, 1760, 363, 923, 1674, 1681, 1675, 928, 151, 1681, 152, 929, 151, 1675, 1681, 944, 335, 1773, 1774, 947, 1713, 1716, 1714, 953, 81, 1719, 1725, 954, 80, 81, 1725, 957, 81, 82, 1719, 962, 94, 95, 1690, 981, 149, 1671, 1677, 985, 1668, 1669, 1673, 986, 1670, 1671, 1672, 1006, 137, 1734, 1737, 1021, 95, 96, 1692, 1053, 347, 1705, 1697, 1073, 158, 1691, 1693, 1076, 157, 1689, 1691, 1077, 1688, 1691, 1689, 1080, 88, 89, 1680, 1081, 1674, 1680, 1681, 1095, 1769, 1772, 1773, 1124, 350, 351, 1762, 1125, 349, 350, 1762, 1126, 1753, 1762, 1761, 1129, 1735, 1738, 1739, 1130, 74, 1736, 1735, 1131, 73, 1735, 1739, 1133, 1695, 1709, 1710, 1137, 144, 1724, 1718, 1154, 1723, 1725, 1724, 1186, 136, 1734, 137, 1188, 137, 1737, 138, 1192, 352, 353, 1761, 1193, 351, 1761, 1762, 1194, 131, 1748, 1746, 1195, 1747, 1748, 1755, 1197, 73, 74, 1735, 1198, 1734, 1738, 1735, 1200, 74, 75, 1736, 1201, 342, 1704, 1703, 1202, 1772, 1774, 1773, 1229, 1746, 1748, 1747, 1233, 69, 1745, 1747, 1234, 349, 1762, 1752, 1235, 350, 1712, 351, 1237, 1698, 1703, 1704, 1238, 340, 1699, 1704, 1239, 95, 1692, 1690, 1240, 1688, 1690, 1691, 1242, 1686, 1688, 1689, 1244, 93, 94, 1688, 1247, 1696, 1710, 1708, 1250, 1697, 1708, 1712, 1260, 138, 1737, 1730, 1261, 1734, 1736, 1737, 1269, 144, 1718, 145, 1270, 1718, 1719, 1721, 1271, 143, 1724, 144, 1314, 145, 1718, 1721, 1320, 130, 1748, 131, 1374, 351, 1712, 352, 1377, 82, 1720, 1719, 1378, 1718, 1725, 1719, 1379, 80, 1725, 1723, 1382, 349, 1712, 350, 1383, 1696, 1705, 1698, 1391, 1698, 1705, 1703, 1402, 345, 1703, 1705, 1408, 1682, 1687, 1683, 1427, 1752, 1762, 1753, 1429, 335, 1702, 336, 1431, 150, 1675, 151, 1435, 150, 1677, 1675, 1439, 1743, 1744, 1745, 1485, 360, 1714, 361, 1487, 5, 363, 1715, 1488, 362, 363, 1760, 1491, 132, 1746, 1744, 1496, 353, 354, 1763, 1534, 131, 1746, 132, 1535, 1745, 1746, 1747, 1536, 70, 1743, 1745, 1538, 145, 1721, 146,
          1539, 1668, 1720, 1669,
          1541, 1668, 1721, 1720,
          1547, 156, 1689, 157)

        """
        op2 = self.op2
        #self.show_data(data[n:], types='i')
        ints = np.frombuffer(data[n:], op2.idtype8)
        iminus1 = np.where(ints == -1)[0]

        istart = [0] + list(iminus1[:-1] + 1)
        iend = iminus1
        #print(ints.tolist())
        #istart = [0] + list(iend + 1)

        #print(istart)
        #print(iend)
        for i0, i1 in zip(istart, iend):
            intsi = ints[i0:i1]
            bsurfs_id = intsi[0]

            nints = len(intsi) - 1
            eids_grids = intsi[1:].reshape(nints // 4, 4)
            eids = eids_grids[:, 0]
            g1s = eids_grids[:, 1]
            g2s = eids_grids[:, 2]
            g3s = eids_grids[:, 3]
            bsurfs = op2.add_bsurfs(bsurfs_id, eids, g1s, g2s, g3s)
            str(bsurfs)
        return len(data)
