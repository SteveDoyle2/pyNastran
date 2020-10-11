"""
defines readers for BDF objects in the OP2 CONTACT/CONTACTS table
"""
from struct import Struct
from typing import Tuple, Any

import numpy as np
from pyNastran.op2.tables.geom.geom_common import GeomCommon
from pyNastran.op2.op2_interface.op2_reader import mapfmt, reshape_bytes_block


class CONTACT(GeomCommon):
    """defines methods for reading contact geometry"""

    def _read_contact_4(self, data, ndata):
        """
        reads the CONTACT/CONTACTS table
        Table of Bulk Data entry related to surface contact

        """
        return self._read_geom_4(self._contact_map, data, ndata)

    def __init__(self):
        GeomCommon.__init__(self)

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
        self._contact_map = {
            # C:\MSC.Software\simcenter_nastran_2019.2\tpl_post1\femao8rand.op2
            (7110, 71, 588) : ['BSURFS', self._read_fake],
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
            (8110, 81, 598) : ['BCTPARM', self._read_fake],
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
            (811, 8, 628) : ['AMLREG', self._read_fake],

            (6621, 66, 662) : ['BCTPAR2', self._read_fake],
            (7610, 76, 593) : ['BGPARA', self._read_fake],
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
            #(4824, 48, 890) : ['???', self._read_fake],
            #(4824, 48, 890) : ['???', self._read_fake],
            #(4824, 48, 890) : ['???', self._read_fake],
            #(4824, 48, 890) : ['???', self._read_fake],
            #(4824, 48, 890) : ['???', self._read_fake],
        }

    #def _read_924(self, data: bytes, n: int) -> int:
        #self.show_data(data[n:])
        #sds

    def _read_bsurfs(self, data: bytes, n: int) -> int:
        """
        BSURFS

        (7110, 71, 588, 1, 24, 190, 198, 189, 44, 188, 197, 190, 64, 106, 189, 196, 84, 195, 188, 106, -1)
                BSURFS         1                              24     190     198     189+
        $           EID2      G1      G2      G3    EID3      G1      G2      G3
        +             44     188     197     190      64     106     189     196+
        $           EID4      G1      G2      G3
        +             84     195     188     106

        """
        bsurfs
