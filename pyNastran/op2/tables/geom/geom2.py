"""
defines readers for BDF objects in the OP2 GEOM2/GEOM2S table
"""
# pylint: disable=C0103
from __future__ import annotations
from struct import Struct
from functools import partial
from typing import Union, Any, TYPE_CHECKING
import numpy as np

from pyNastran.bdf.errors import UnsupportedCard
from pyNastran.bdf.cards.elements.elements import CGAP, PLOTEL
from pyNastran.bdf.cards.elements.damper import (CDAMP1, CDAMP2, CDAMP3,
                                                 CDAMP4, CDAMP5, CVISC)
from pyNastran.bdf.cards.elements.springs import CELAS1, CELAS2, CELAS3, CELAS4
from pyNastran.bdf.cards.elements.axisymmetric_shells import (
    CQUADX, CTRIAX6, CTRAX3, CTRAX6, CQUADX8, CTRIAX)
from pyNastran.bdf.cards.elements.shell import (CTRIA3, CQUAD4, CTRIA6,
                                                CQUADR, CTRIAR,
                                                CQUAD8, CQUAD,
                                                CSHEAR)
from pyNastran.bdf.cards.elements.rods import CROD, CTUBE, CONROD
from pyNastran.bdf.cards.elements.bars import CBAR, CBEND, CBEAM3
from pyNastran.bdf.cards.elements.beam import CBEAM
from pyNastran.bdf.cards.elements.mass import (CONM1, CONM2, CMASS1, CMASS2,
                                               CMASS3, CMASS4)
from pyNastran.bdf.cards.elements.solid import (CTETRA4, CPYRAM5, CPENTA6, CHEXA8,
                                                CTETRA10, CPYRAM13, CPENTA15, CHEXA20,
                                                CPENTCZ, CHEXCZ)
from pyNastran.bdf.cards.thermal.thermal import CHBDYG, CONV, CHBDYP, CHBDYE, CONVM
from pyNastran.bdf.cards.thermal.radiation import RADBC # , RADM, RADCAV, RADLST, RADMTX, VIEW, VIEW3D
from pyNastran.bdf.cards.nodes import SPOINTs
from pyNastran.bdf.cards.elements.bush import CBUSH
from pyNastran.bdf.cards.parametric.geometry import FEEDGE
from pyNastran.bdf.cards.elements.acoustic import CHACAB, CHACBR, CAABSF
from pyNastran.op2.errors import MixedVersionCard
from pyNastran.op2.op2_interface.op2_reader import mapfmt # , reshape_bytes_block
from pyNastran.op2.tables.geom.geom4 import RBE3

from pyNastran.op2.errors import DoubleCardError, EmptyCardError
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.op2.op2_geom import OP2Geom

def _map_offt(num: int) -> str:
    offt = ['G', 'G', 'G']
    # 1G->1B: 4^3; 1+64=65
    if num > 64:
        num -= 64 # 4^3
        offt[0] = 'B'

    # 3G->3O: 4^2; 1+16=17
    # 3G->3B: 2*4^2; 1+32=33
    if num > 32:
        num -= 32 # 4^2
        offt[2] = 'O'
    elif num > 16:
        num -= 16 # 4^2
        offt[2] = 'B'

    # 2G->2O: 4*1; 1+4=5
    # 2G->2B: 4*2; 1+8=9
    if num > 8:
        num -= 8 # 4*2
        offt[1] = 'O'
    elif num > 4:
        num -= 4 # 4
        offt[1] = 'B'

    offt_str = ''.join(offt)
    return offt_str

BAR_FE_MAP = {
    #B/G   2
    #B/G/O 3
    #GGG = 1
    #1: 'GGG',
    #5: 'GOG', # 2G->2O: 4*1; 1+4=5
    #9: 'GBG', # 2G->2B: 4*2; 1+8=9
    #17: 'GGO', # 3G->3O: 4^2; 1+16=17
    #21: 'GOO', # 2G->2O: 4*1; 17+4=21
    #33: 'GGB', # 3G->3B: 2*4^2; 1+32=33
    #41: 'GBB', # 2G->2B: 4*2; 33+8=41

    #65: 'BGG', # 1G->1B: 4^3; 1+64=65
    #69: 'BOG', # 2G->2O: 4*1; 65+4=69
    #73: 'BBG', # 2G->2B: 4*2; 65+8=73
    #81: 'BGO', # 3G->3O: 4^2; 65+16=81
    #85: 'BOO', # 2G->2O: 4*1; 81+4=85
    #97: 'BGB', # 3G->3B: 2*4^2; 65+32=97
    #105: 'BBB', # 2G->2B: 4*2; 97+8=105
    #
    # 1G->1B: 4^3; 1+64=65
    #
    # 2G->2O: 4*1; 1+4=5
    # 2G->2B: 4*2; 1+8=9
    #
    # 3G->3O: 4^2; 1+16=17
    # 3G->3B: 2*4^2; 1+32=33

    #fprime = fe-f
    # f=1    f=2
    0: 'GGG', # NX?
    1: 'GGG', 2: 'GGG',
    5: 'GOG', 6: 'GOG',
    9: 'GBG', 10: 'GBG',
    17: 'GGO', 18: 'GGO',
    21: 'GOO', 22: 'GOO',
    33: 'GGB', 34: 'GGB',
    41: 'GBB', 42: 'GBB',

    65: 'BGG',
    69: 'BOG',
    73: 'BBG',
    81: 'BGO',
    85: 'BOO',
    97: 'BGB',
    105: 'BBB',


    #https://docs.plm.automation.siemens.com/data_services/resources/scnastran/2020_1/help/tdoc/en_US/pdf/release_guide.pdf
    #F = 1 XYZ option - global or basic coordinate system
    #F = 2 Grid option
    #F = 5 XYZ option - global or basic coordinate system
    #F = 6 XYZ option - global or basic coordinate system
    #F = 17 XYZ option - global or basic coordinate system
    #F = 18 XYZ option - global or basic coordinate system
    #F = 21 XYZ option - global or basic coordinate system
    #F = 22 XYZ option - global or basic coordinate system
    #F = 65 XYZ option - global or basic coordinate system
    #F = 69 XYZ option - global or basic coordinate system
    #F = 81 XYZ option - global or basic coordinate system
    #F = 85 XYZ option - global or basic coordinate system

}

class GEOM2:
    """defines methods for reading op2 elements"""
    #@property
    #def struct_i(self) -> Struct:
        #return self.op2.struct_i
    @property
    def struct_q(self) -> Struct:
        return self.op2.struct_q

    #@property
    #def idtype(self) -> str:
        #return self.op2.idtype
    #@property
    #def idtype8(self) -> str:
        #return self.op2.idtype8
    #@property
    #def fdtype8(self) -> str:
        #return self.op2.fdtype8

    #@property
    #def log(self) -> Any:
        #return self.op2.log
    #@property
    #def is_debug_file(self) -> bool:
        #return self.op2.is_debug_file
    #@property
    #def card_count(self) -> dict[str, int]:
        #return self.op2.card_count
    #@property
    #def binary_debug(self) -> Any:
        #return self.op2.binary_debug

    #@property
    #def _endian(self) -> bytes:
        #return self.op2._endian
    @property
    def size(self) -> int:
        return self.op2.size
    @property
    def factor(self) -> int:
        return self.op2.factor

    def read_fake(self, data: bytes, n: int) -> int:
        return self.op2._read_fake(data, n)

    def read_stop(self, data: bytes, n: int) -> int:
        return self.op2.reader_geom1.read_stop(data, n)

    def read_geom2_4(self, data: bytes, ndata: int):
        return self.op2._read_geom_4(self.geom2_map, data, ndata)

    def __init__(self, op2: OP2Geom):
        self.op2 = op2
        self.geom2_map = {
            # per dmap-nx-10.pdf or nx12.pdf
            (15200, 152, 9912): ['ACFACE3', self.read_fake],
            (15500, 155, 9913): ['ACFACE4', self.read_fake],
            (15600, 156, 9914): ['ACFACE6', self.read_fake],
            (15700, 157, 9915): ['ACFACE8', self.read_fake],
            (2601, 26, 0): ['BEAMAERO', self.read_fake],
            (2708, 27, 59): ['CAABSF', self.read_caabsf],
            (2108, 21, 224): ['CAXIF2', self.read_fake],
            (2208, 22, 225): ['CAXIF3', self.read_fake],
            (2308, 23, 226): ['CAXIF4', self.read_fake],

            (2408, 24, 180): ['CBAR', self.read_cbar],         # record 8
            (4001, 40, 275): ['CBARAO', self.read_cbarao],     # record 9  - not done
            (5408, 54, 261): ['CBEAM', self.read_cbeam],       # record 10
            (11401, 114, 9016): ['CBEAMP', self.read_cbeamp],  # record 11 - not done
            (4601, 46, 298): ['CBEND', self.read_cbend],     # record 12 - not done
            (2608, 26, 60): ['CBUSH', self.read_cbush],      # record 13
            (5608, 56, 218): ['CBUSH1D', self.read_cbush1d], # record 14 - not done
            (5609, 60, 9899): ['CBUSH1DNL', self.read_fake],
            (14801, 148, 956): ['CCHOCK3', self.read_fake],
            (14901, 149, 957): ['CCHOCK4', self.read_fake],
            (15001, 150, 958): ['CCHOCK6', self.read_fake],
            (15101, 151, 959): ['CCHOCK8', self.read_fake],
            (2315, 23, 146): ['CCONE-10', self.read_ccone],  # nx10
            (2315, 23, 0): ['CCONE-12', self.read_fake],     # nx12
            (201, 2, 69): ['CDAMP1', self.read_cdamp1],
            (301, 3, 70): ['CDAMP2', self.read_cdamp2],
            (401, 4, 71): ['CDAMP3', self.read_cdamp3],
            (501, 5, 72): ['CDAMP4', self.read_cdamp4],
            (10608, 106, 404): ['CDAMPS', self.read_cdamp5],
            (6208, 62, 108): ['CDUM2', self.read_fake],
            (6308, 63, 109): ['CDUM3', self.read_fake],
            (6408, 64, 110): ['CDUM4', self.read_fake],
            (6508, 65, 111): ['CDUM5', self.read_fake],
            (6608, 66, 112): ['CDUM6', self.read_fake],
            (6708, 67, 113): ['CDUM7', self.read_fake],
            (6808, 68, 114): ['CDUM8', self.read_cdum8],
            (6908, 69, 115): ['CDUM9', self.read_cdum9],
            (601, 6, 73): ['CELAS1', self.read_celas1],
            (6010, 53, 9900): ['CELAS1NL', self.read_fake],
            (701, 7, 74): ['CELAS2', self.read_celas2],
            (7010, 20, 9898): ['CELAS2NL', self.read_fake],
            (801, 8, 75): ['CELAS3', self.read_celas3],
            (901, 9, 76): ['CELAS4', self.read_celas4],

            (9801, 98, 506): ['CFAST-nx10', self._read_cfast_msc_nx10],  # nx10
            (13801, 138, 566): ['CFAST-nx12', self.read_cfast], # nx12

            (8515, 85, 209): ['CFLUID2-nx10', self.read_cfluid2],  # nx10
            (8615, 86, 210): ['CFLUID3-nx10', self.read_cfluid3],  # nx10
            (8715, 87, 211): ['CFLUID4-nx10', self.read_cfluid4],  # nx10
            (8515, 85, 0): ['CFLUID2-nx12', self.read_cfluid2], # nx12
            (8615, 86, 0): ['CFLUID3-nx12', self.read_cfluid3], # nx12
            (8715, 87, 0): ['CFLUID4-nx12', self.read_cfluid4], # nx12

            (1908, 19, 104): ['CGAP', self.read_cgap],
            #(13101, 131, 9901): ['CGPLSTN3', self.read_fake],
            #(13401, 134, 9904): ['CGPLSTN8', self.read_fake],
            #(13301, 133, 9903): ['CGPLSTN6', self.read_fake],
            #(13201, 132, 9902): ['CGPLSTN4', self.read_fake],
            (8308, 83, 405): ['CHBDYE', self.read_chbdye],
            (10808, 108, 406): ['CHBDYG', self.read_chbdyg],
            (10908, 109, 407): ['CHBDYP', self.read_chbdyp],
            #(7308, 73, 253): ['CHEXA', self.read_fake],
            #(16300, 163, 9989): ['CHEXA20F', self.read_fake],
            #(14100, 141, 9990): ['CHEXAF', self.read_fake],
            #(14000, 140, 9990): ['CHEXAFD', self.read_fake],
            #(7708, 77, 9944): ['CHEXAL', self.read_fake],
            (11801, 118, 907): ['CHEXCZ', self.read_chexa_cz],
            #(12001, 120, 9011): ['CHEXP', self.read_fake],
            #(7409, 74, 9991): ['CHEXPR', self.read_fake],
            (1001, 10, 65): ['CMASS1', self.read_cmass1],    # record 52
            (1101, 11, 66): ['CMASS2', self.read_cmass2],    # record 53
            (1201, 12, 67): ['CMASS3', self.read_cmass3],    # record 54
            (1301, 13, 68): ['CMASS4', self.read_cmass4],    # record 55
            (2508, 25, 0): ['CMFREE', self.read_cmfree],     # record 56 - not done
            (1401, 14, 63): ['CONM1', self.read_conm1],      # record 57
            (1501, 15, 64): ['CONM2', self.read_conm2],      # record 58
            (1601, 16, 47): ['CONROD', self.read_conrod],    # record 59
            (12701, 127, 408): ['CONV', self.read_conv],
            (8908, 89, 422): ['CONVM', self.read_convm],
            (12101, 121, 9012): ['CPENP', self.read_fake],
            (4108, 41, 280): ['CPENTA', self.read_cpenta],
            (14200, 142, 9906): ['CPENTAF', self.read_cpenta],
            #(7108, 71, 9943): ['CPENTAL', self.read_fake],
            (7509, 75, 9992): ['CPENPR', self.read_fake],
            #(16500, 165, 9987): ['CPENT15F', self.read_fake],
            #(16000, 160, 9988): ['CPENT6FD', self.read_fake],
            (11901, 119, 908): ['CPENTCZ', self.read_cpenta_cz],
            (1701, 17, 980): ['CPLSTN3', self.read_cplstn3],
            (5701, 57, 981): ['CPLSTN4', self.read_cplstn4],
            (5801, 58, 982): ['CPLSTN6', self.read_cplstn6],
            (7201, 72, 983): ['CPLSTN8', self.read_cplstn8],
            (8801, 88, 984): ['CPLSTS3', self.read_cplsts3],
            (8401, 84, 985): ['CPLSTS4', self.read_cplsts4],
            (1801, 18, 986): ['CPLSTS6', self.read_cplsts6],
            (3601, 36, 987): ['CPLSTS8', self.read_cplsts8],
            (17200, 172, 1000) : ['CPYRAM', self.read_cpyram], # nx-specific
            (14400, 144, 9908): ['CPYRAMF', self.read_cpyram], # nx-specific
            (25700, 257, 9948) : ['CPYRA5FD', self.read_cpyram], # nx-specific
            (25800, 258, 9947) : ['CPYRA13F', self.read_cpyram], # nx-specific
            (7909, 79, 9946) : ['CPYRAMPR', self.read_cpyram], # nx-specific
            (17000, 170, 9980): ['CQDX4FD', self.read_cquad],
            (17100, 171, 9979): ['CQDX9FD', self.read_cquadx],
            #(25110, 170, 9951): ['CQDX4FDN', self.read_fake],
            #(25310, 171, 9949): ['CQDX8FDN', self.read_fake],
            (9108, 91, 507): ['CQUAD', self.read_cquad],
            (2958, 51, 177): ['CQUAD4', self.read_cquad4],
            #(14600, 146, 9910): ['CQUADF', self.read_fake],
            (13900, 139, 9984): ['CQUAD4FD', self.read_cquad],
            #(4701, 47, 326): ['CQUAD8', self.read_fake],
            #(3302, 33, 1694): ['CQUAD8L', self.read_fake],
            #(15901, 159, 9956): ['CQUAD8N', self.read_fake],
            (16400, 164, 9983): ['CQUAD9FD', self.read_cquad],
            #(11101, 111, 9014): ['CQUADP', self.read_fake],
            #(8009, 80, 367): ['CQUADR', self.read_fake],
            #(13002, 130, 1692): ['CQUADRL', self.read_fake],
            #(15401, 154, 9954): ['CQUADRN', self.read_fake],
            #(9008, 90, 508): ['CQUADX', self.read_fake],
            (6112, 61, 997): ['CQUADX4', self.read_fake],
            (6114, 61, 999): ['CQUADX8', self.read_cquadx8],
            (3001, 30, 48): ['CROD', self.read_crod],         # record 81
            (14500, 145, 9909): ['CRODF', self.read_crod],
            (3501, 35, 1): ['CSBOLT', self.read_fake],
            (3101, 31, 61): ['CSHEAR', self.read_cshear],     # record 84
            (4408, 44, 227): ['CSLOT3', self.read_fake],
            (4508, 45, 228): ['CSLOT4', self.read_fake],
            #(12201, 122, 9013): ['CTETP', self.read_fake],
            #(5508, 55, 217): ['CTETRA', self.read_fake],
            (14300, 143, 9907): ['CTETRAF', self.read_ctetra],
            #(7609, 76, 9993): ['CTETPR', self.read_fake],
            #(16600, 166, 9985): ['CTETR10F', self.read_fake],
            #(16100, 161, 9986): ['CTETR4FD', self.read_fake],
            #(25010, 168, 9952): ['CTRAX3FDN', self.read_fake],
            #(25210, 169, 9950): ['CTRAX6FDN', self.read_fake],
            #(5959, 59, 282): ['CTRIA3', self.read_fake],
            #(14700, 147, 9911): ['CTRIAF', self.read_fake],
            #(16200, 162, 9982): ['CTRIA3FD', self.read_fake],
            #(4801, 48, 327): ['CTRIA6', self.read_fake],
            #(16700, 167, 9981): ['CTRIA6FD', self.read_fake],
            #(3202, 32, 1693): ['CTRIA6L', self.read_fake],
            (15801, 158, 9955): ['CTRIA6N', self.read_ctria6],
            #(11301, 113, 9015): ['CTRIAP', self.read_fake],
            #(9200, 92, 385): ['CTRIAR', self.read_fake],
            #(12902, 129, 1691): ['CTRIARL', self.read_fake],
            #(15301, 153, 9953): ['CTRIARN', self.read_fake],
            #(10108, 101, 512): ['CTRIAX', self.read_fake],
            #(6108, 61, 107): ['CTRIAX6', self.read_fake],
            (6111, 61, 996): ['CTRAX3', self.read_ctrax3],
            (6113, 61, 998): ['CTRAX6', self.read_ctrax6],
            #(16800, 168, 9978): ['CTRIX3FD', self.read_fake],
            #(16900, 169, 9977): ['CTRIX6FD', self.read_fake],
            (3701, 37, 49): ['CTUBE', self.read_ctube], # record 104
            (3901, 39, 50): ['CVISC', self.read_cvisc], # record 105
            #(11701, 117, 559): ['CWELD', self.read_fake],
            #(13501, 135, 564): ['CWELDC', self.read_fake],
            (13701, 137, 565): ['CWELDP', self.read_fake],
            #(13601, 136, 562): ['CWELDG', self.read_fake], # This record is no longer used
            (4301, 43, 28): ['GENEL', self.read_genel],
            #(3201, 32, 478): ['GMBNDC', self.read_fake],
            (12901, 129, 482): ['GMBNDS', self.read_gmbnds],
            (3301, 33, 479): ['GMINTC', self.read_fake],
            (13001, 130, 483): ['GMINTS', self.read_fake],
            (2801, 28, 630): ['MICPNT', self.read_micpnt],
            (5201, 52, 11): ['PLOTEL', self.read_plotel],
            #(5202, 52, 669): ['PLOTEL3', self.read_fake],
            (5203, 52, 670): ['PLOTEL4', self.read_fake],
            #(5204, 52, 671): ['PLOTEL6', self.read_fake],
            #(5205, 52, 672): ['PLOTEL8', self.read_fake],
            #(5206, 52, 673): ['PLOTHEX', self.read_fake],
            #(5208, 52, 675): ['PLOTPEN', self.read_fake],
            #(5209, 52, 676): ['PLOTPYR', self.read_fake],
            #(5207, 52, 674): ['PLOTTET', self.read_fake],
            #(3002, 46, 0): ['Q4AERO', self.read_fake],
            (12801, 128, 417): ['RADBC', self.read_radbc],
            (7801, 78, 8883): ['SINT', self.read_sint],
            (5551, 49, 105): ['SPOINT', self.read_spoint], # record 119
            #(2701, 27, 0): ['T3AERO', self.read_fake],
            #(11601, 116, 9942): ['VUBEAM', self.read_fake],
            #(12301, 123, 145): ['VUHEXA', self.read_fake],
            #(12401, 124, 146): ['VUPENTA', self.read_fake],
            #(11201, 112, 9940): ['VUQUAD4', self.read_fake],
            #(12501, 125, 147): ['VUTETRA', self.read_fake],
            #(11501, 115, 9941): ['VUTRIA3', self.read_fake],
            #(65535, 65535, 65535): ['EODB', self.read_fake],


            #------------------------------------------------------------------------------
            # MSC - DMAP 2016.1
            (2002, 20, 0): ['AEROQ4', self.read_fake], # 2
            (1801, 18, 0): ['AEROT3', self.read_fake], # 3
            (1701, 17, 0): ['BEAMAERO', self.read_fake], # 4
            (2708, 27, 59): ['CAABSF', self.read_caabsf], # 5
            (2108, 21, 224): ['CAXIF2', self.read_fake], # 6
            (2208, 22, 225): ['CAXIF3', self.read_fake], # 7
            (2308, 23, 226): ['CAXIF4', self.read_fake], # 8

            #(4001, 40, 275): ['CBARAO', self.read_fake],    # record 10
            #(5408, 54, 261): ['CBEAM', self.read_fake],     # record 11
            #(11401, 114, 9016): ['CBEAMP', self.read_fake], # record 12
            #(4601, 46, 298): ['CBEND', self.read_fake],   # record 13
            #(2608, 26, 60): ['CBUSH', self.read_fake],    # record 14
            #(5608, 56, 218): ['CBUSH1D', self.read_fake], # record 15
            (2315, 23, 0): ['CCONE-msc', self.read_fake],  # record 16
            #(201, 2, 69): ['CDAMP1', self._read_cdamp1], # record 17
            #(301, 3, 70): ['CDAMP2', self._read_cdamp2], # record 18
            #(401, 4, 71): ['CDAMP3', self._read_cdamp3], # record 19
            #(501, 5, 72): ['CDAMP4', self._read_cdamp4],      # record 20
            #(10608, 106, 404): ['CDAMP5', self._read_cdamp5], # record 21
            #(6208, 62, 108): ['CDUM2', self.read_fake], # 22
            #(6308, 63, 109): ['CDUM3', self.read_fake], # 23
            #(6408, 64, 110): ['CDUM4', self.read_fake], # 24
            #(6508, 65, 111): ['CDUM5', self.read_fake], # 25
            #(6608, 66, 112): ['CDUM6', self.read_fake], # 26
            #(6708, 67, 113): ['CDUM7', self.read_fake], # 27
            #(6808, 68, 114): ['CDUM8', self._read_cdum8], # 28
            #(6908, 69, 115): ['CDUM9', self._read_cdum9], # 29

            #(601, 6, 73): ['CELAS1', self.read_fake], # record 30
            #(701, 7, 74): ['CELAS2', self.read_fake], # record 31
            #(801, 8, 75): ['CELAS3', self.read_fake], # record 32
            #(901, 9, 76): ['CELAS4', self.read_fake], # record 33
            #(9801, 98, 506): ['CFAST', self.read_fake], # record 34
            (9301, 93, 607): ['CFASTP', self.read_cfastp],    # 35
            #(8515, 85, 0): ['CFLUID2', self._read_cfluid2],  # record 36 - not done
            #(8615, 86, 0): ['CFLUID3', self._read_cfluid3],  # record 37 - not done
            #(8715, 87, 0): ['CFLUID4', self._read_cfluid4],  # record 38 - not done
            (7701, 77, 8881): ['CINT', self.read_cint],      # record 39 - not done

            #(1908, 19, 104): ['CGAP', self._read_cgap],      # record 40 - buggy
            (8100, 81, 381): ['CHACAB', self.read_chacab],    # 41
            (8200, 82, 383): ['CHACBR', self.read_chacbr],    # 42 - not done
            #(8308, 83, 405): ['CHBDYE', self._read_chbdye],  # record 43
            #(10808, 108, 406): ['CHBDYG', self.read_fake],  # 44
            #(10908, 109, 407): ['CHBDYP', self.read_fake],  # 45
            (7308, 73, 253): ['CHEXA', self.read_chexa],    # record 46
            (16300, 163, 9999): ['CHEXA20F', self.read_fake], # 47
            (14000, 140, 9990): ['CHEXAFD', self.read_chexa], # record 48
            (7908, 79, 369): ['CHEXAL', self.read_fake],     # 49

            (12001, 120, 9011): ['CHEXP', self.read_fake],  # record 50
            (7409, 74, 9991): ['CHEXPR', self.read_chexpr],  # record 51
            #(1001, 10, 65): ['CMASS1', self._read_cmass1],   # record 52
            #(1101, 11, 66): ['CMASS2', self._read_cmass2],   # record 53
            #(1201, 12, 67): ['CMASS3', self._read_cmass3],   # record 54
            #(1301, 13, 68): ['CMASS4', self._read_cmass4],   # record 55
            #(2508, 25, 0): ['CMFREE', self._read_cmfree],    # record 56 - not done
            #(1401, 14, 63): ['CONM1', self._read_conm1],     # record 57
            #(1501, 15, 64): ['CONM2', self._read_conm2],     # record 58
            #(1601, 16, 47): ['CONROD', self._read_conrod],   # record 59

            #(12701, 127, 408): ['CONV', self._read_conv],     # record 60 - not tested
            #(8908, 89, 422): ['CONVM', self._read_convm],     # record 61 - not tested
            #(12101, 121, 9012): ['CPENP', self.read_fake],   # 62
            #(4108, 41, 280): ['CPENTA', self._read_cpenta],   # record 63
            #(7509, 75, 9992): ['CPENPR', self.read_fake],    # 64
            (16500, 165, 9999): ['CPENT15F', self.read_fake], # 65
            (16000, 160, 9999): ['CPENT6FD', self.read_fake], # 66
            (17000, 170, 9999): ['CQDX4FD', self.read_fake],  # 67
            (17100, 171, 9999) : ['CQDX9FD', self.read_cquadx],  # record 68
            #(9108, 91, 507): ['CQUAD', self._read_cquad],        # record 69 - not tested

            #(2958, 51, 177): ['CQUAD4', self._read_cquad4],     # record 70
            (13900, 139, 9989): ['CQUAD4FD', self.read_cquad4], # record 71
            (4701, 47, 326): ['CQUAD8', self.read_cquad8], # record 72
            (16400, 164, 9999): ['CQUAD9FD', self.read_fake], # 73
            (11101, 111, 9014): ['CQUADP', self.read_fake], # 74
            (8009, 80, 367): ['CQUADR', self.read_cquadr], # record 75
            (9008, 90, 508): ['CQUADX', self.read_cquadx], # record 76
            (14700, 147, 6662): ['CRBAR', self.read_crbar], # 77
            (17300, 173, 6664): ['CRBE1', self.read_crbe1], # 78
            (17200, 172, 6663): ['CRBE3', self.read_crbe3], # 79

            (11000, 110, 6667): ['CRJOINT', self.read_crjoint],  # 80
            #(3001, 30, 48): ['CROD', self._read_crod],         # record 81
            (12600, 126, 6661): ['CRROD', self.read_crrod],    # 82
            (13801, 138, 570): ['CSEAM', self.read_fake],     # 83
            #(3101, 31, 61): ['CSHEAR', self._read_cshear],     # record 84
            #(4408, 44, 227): ['CSLOT3', self.read_fake],      # 85
            #(4508, 45, 228): ['CSLOT4', self.read_fake],      # 86
            (12201, 122, 9013): ['CTETP', self.read_ctetrap], # record 87
            (5508, 55, 217): ['CTETRA', self.read_ctetra],    # record 88
            (7609, 76, 9993): ['CTETPR', self.read_ctetra],   # record 89

            (16600, 166, 9999): ['CTETR10F', self.read_ctetra], # record 90
            (16100, 161, 9999): ['CTETR4FD', self.read_ctetra], # record 91
            (14801, 148, 643): ['CTQUAD', self.read_fake], # 92
            (14901, 149, 644): ['CTTRIA', self.read_fake], # 93
            (5959, 59, 282): ['CTRIA3', self.read_ctria3], # record 94
            (16200, 162, 9999): ['CTRIA3FD', self.read_fake], # 95
            (4801, 48, 327): ['CTRIA6', self.read_ctria6], # record 96 - buggy
            (16700, 167, 9999): ['CTRIA6FD', self.read_ctria6], # 97
            (11301, 113, 9015): ['CTRIAP', self.read_fake], # 98
            (9200, 92, 385): ['CTRIAR', self.read_ctriar], # record 99

            (6108, 61, 107): ['CTRIAX6', self.read_ctriax6], # 101
            (16800, 168, 9978): ['CTRIX3FD', self.read_fake], # 102
            (16900, 169, 9977): ['CTRIX6FD', self.read_fake], # 103
            #(3701, 37, 49): ['CTUBE', self.read_ctube], # record 104
            #(3901, 39, 50): ['CVISC', self.read_cvisc], # record 105
            (11701, 117, 559): ['CWELD', self.read_fake], # 106; same as cfast
            (13501, 135, 564): ['CWELDC', self.read_fake], # 107
            (13601, 136, 562): ['CWELDG', self.read_fake], # 108
            (14600, 146, 630): ['CWSEAM', self.read_fake], # 109

            #(4301, 43, 28): ['GENEL', self.read_genel],     # 110
            (3201, 32, 478): ['GMBNDC', self.read_gmbndc],   # 111
            #(12901, 129, 482): ['GMBNDS', self.read_gmbnds], # 112
            #(3301, 33, 479): ['GMINTC', self.read_fake],   # 113
            #(13001, 130, 483): ['GMINTS', self.read_fake], # 114
            #(5201, 52, 11): ['PLOTEL', self.read_plotel],  # record 115
            #(12801, 128, 417): ['RADBC', self.read_radbc], # record 116
            (15501, 155, 634): ['RADINT', self.read_fake], # 117
            #(7801, 78, 8883): ['SINT', self.read_fake],    # 118
            #(5551, 49, 105): ['SPOINT', self.read_spoint], # record 119

            (11601, 116, 9942): ['VUBEAM', self.read_vubeam], # record 120
            (12301, 123, 145): ['VUHEXA', self.read_fake],    # 121
            (11201, 112, 9940): ['VUQUAD4', self.read_vuquad4],  # 122
            (12401, 124, 146): ['VUPENTA', self.read_fake],   # 123
            (12501, 125, 147): ['VUTETRA', self.read_fake],   # 124
            (11501, 115, 9941): ['VUTRIA3', self.read_vutria3],  # 125
            (13701, 137, 569): ['WELDP', self.read_fake],     # 126; same as CFASTP

            #----------------------------------------------

            # unorganized
            (6113, 61, 998): ['CTRAX6', self.read_ctrax6],
            (10108, 101, 512) : ['CTRIAX', self.read_ctriax],
            (2108, 21, 224): ['CAXIF2', self.read_fake],
            (5601, 56, 296): ['SESET', self.read_fake],
            (7509, 75, 9992): ['CPENPR', self.read_cpenta],
            (16000, 160, 9988): ['CPENTA6FD', self.read_cpenta],
            (16100, 161, 9986): ['CTETRAFD', self.read_ctetra],
            (16300, 163, 9989): ['CHEXA20F', self.read_chexa],
            (16700, 167, 9981): ['CTRI6FD', self.read_ctria6fd],
            (16800, 168, 9978): ['CTRIAX3FD', self.read_ctriax3fd],  # same as ctria6fd
            (16500, 165, 9987): ['CPENT15F', self.read_cpenta],
            (5008, 50, 258): ['CNGRET', self.read_cngret],
            (12301, 123, 9921): ['ADAPT card', self.read_adapt],
            (12401, 124, 9922): ['FEFACE/PVAL?', self.read_feface_pval],
            (7309, 73, 0): ['CaseControl SET?', self.read_fake],
            (12501, 125, 9923): ['ADAPT card 2', self.read_adapt],    # record
            (3401, 34, 9600): ['GMCONV?', self.read_fake],    # record
            (2901, 29, 9601): ['FEEDGE', self.read_feedge2],  # record
            (16600, 166, 9985) : ['CTETRA?', self.read_ctetra],  # record
            (16200, 162, 9982) : ['CTRIA3', self.read_ctria3fd],  # record
            (16900, 169, 9977) : ['CTRIAX', self.read_ctriax],  # record
            (23500, 235, 6662) : ['', self.read_fake],  # record
            (23800, 238, 6665) : ['', self.read_fake],  # record
            (23900, 239, 6666) : ['', self.read_fake],  # record

            (1976, 1, 1996) : ['', self.read_fake],  # record
            (6120, 1, 60434) : ['', self.read_fake],  # record
            (2024, 1001, 2024) : ['', self.read_fake],  # record
            (801, 1, 572) : ['', self.read_fake],  # record

            (1001, 100, 10000) : ['', self.read_fake],  # record
            (1118, 1, 1874) : ['', self.read_fake],  # record

            (7708, 77, 9944): ['CHEXAL', self.read_fake],  # record
            (7108, 71, 9943): ['CPENTAL', self.read_fake],  # record
            (11001, 110, 8881): ['???', self.read_fake],
            (15301, 153, 9953): ['CTRIARN', self.read_ctria3],
            (15401, 154, 9954): ['CQUADRN', self.read_cquad4],

            (9508, 95, 9801): ['CQUADX', self.read_cquadx_9508],

            (15418, 154, 610): ['CBEAM3', self.read_cbeam3],
            (15901, 159, 9956): ['CQUAD8N', self.read_cquad8],
            (14600, 146, 9910): ['CQUAD4F', self.read_cquad4],
            (7908, 79, 9702): ['CSEAM?', self.read_cseam_maybe],

            (14100, 141, 9905): ['CHEXA', self.read_chexa20],
            (14700, 147, 9911): ['CTRIAF', self.read_ctria3],
            (9301, 93, 690): ['CJOINT', self.read_fake],
            #(14200, 142, 9906): ['???', self.read_fake],
            #(15801, 158, 9955): ['???', self.read_fake],
            #(15801, 158, 9955): ['???', self.read_fake],
            (17608, 176, 993): ['CPYRAM', self.read_cpyram15],
            (8408, 84, 994): ['CPYRAM', self.read_cpyram15],

            #C:\MSC.Software\msc_nastran_runs\sl_1blk.op2
           (11908, 119, 469): ['CSEAM-a', self.read_fake],
           (15701, 157, 639): ['CSEAM-b', self.read_fake],
           (15601, 156, 635): ['CBUSH2D', self.read_cbush2d],
        }

    def read_cpyram15(self, data: bytes, n: int) -> int:
        """15 fields"""
        op2: OP2Geom = self.op2
        ntotal = 15 * self.size
        struct_15i = Struct(mapfmt(op2._endian + b'15i', self.size))
        nelements = (len(data) - n) // ntotal
        for unused_i in range(nelements):
            edata = data[n:n + ntotal]
            out = struct_15i.unpack(edata)
            if op2.is_debug_file:
                op2.binary_debug.write('  CPENTA=%s\n' % str(out))
            (eid, pid, g1, g2, g3, g4, g5, g6, g7, g8, g9, g10,
             g11, g12, g13) = out

            data_in = [eid, pid, g1, g2, g3, g4, g5]
            big_nodes = [g6, g7, g8, g9, g10, g11, g12, g13]
            if sum(big_nodes) > 0:
                elem = CPYRAM13.add_op2_data(data_in + big_nodes)
            else:
                elem = CPYRAM5.add_op2_data(data_in)
            self.add_op2_element(elem)
            n += ntotal
        op2.card_count['CPYRAM'] = nelements
        return n

    def read_chexa20(self, data: bytes, n: int) -> int:
        """
        ints    = (14100, 141, 9905,
                   8013, 4, 6035, 6034, 6036, 6037, 6039, 6038, 6040, 6041,
                   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
        ]
        """
        op2: OP2Geom = self.op2
        ints = np.frombuffer(data[n:], dtype=op2.idtype8)
        nfields = len(ints)
        nelements = nfields // 22
        ints = ints.reshape(nelements, 22)
        for eid, pid, g1, g2, g3, g4, g5, g6, g7, g8, *big_nodes in ints:
            data_in = [eid, pid, g1, g2, g3, g4, g5, g6, g7, g8, ]
            #big_nodes = [g9, g10, g11, g12, g13, g14, g15, g16,
                         #g17, g18, g19, g20]
            if sum(big_nodes) > 0:
                elem = CHEXA20.add_op2_data(data_in + big_nodes)
            else:
                elem = CHEXA8.add_op2_data(data_in)
            self.add_op2_element(elem)
        op2.card_count['CHEXA'] = nelements
        return len(data)

    def read_cquadx_9508(self, data: bytes, n: int) -> int:
        r"""
        ints    = (1, 1, [1, 2,  8,  7], [0, 0, 0, 0, 0, 0, -1],
                   2, 1, [2, 3,  9,  8], [0, 0, 0, 0, 0, 0, -1],
                   3, 1, [3, 4, 10,  9], [0, 0, 0, 0, 0, 0, -1],
                   4, 1, [4, 5, 11, 10], [0, 0, 0, 0, 0, 0, -1],
                   5, 1, [5, 6, 12, 11], [0, 0, 0, 0, 0, 0, -1])
        C:\MSC.Software\msc_nastran_runs\axh101a2.op2
        """
        op2: OP2Geom = self.op2
        ntotal = 52 * self.factor # 16*4
        nelements = (len(data) - n) // ntotal
        s = Struct(mapfmt(op2._endian + b'2i 4i 7i', self.size))

        for unused_i in range(nelements):
            edata = data[n:n + ntotal]
            out = s.unpack(edata)
            #print(out)
            (eid, pid,
             n1, n2, n3, n4,
             n5, n6, n7, n8, n9, f0, gm1) = out
            assert (n5, n6, n7, n8, n9, f0, gm1) == (0, 0, 0, 0, 0, 0, -1)
            nids = [n1, n2, n3, n4,
                    n5, n6, n7, n8, n9]
            elem = CQUADX(eid, pid, nids)
            op2._add_methods._add_element_object(elem)
            n += ntotal
        op2.card_count['CQUADX'] = nelements
        return n

    #def _show_geom2_fake(self, data: bytes, n: int):
        #"""
        #ints    = (1, 2, 1, 2, 2, 2, 1, 2, 11, 12, 16, 21, 25, 2, 3, 28, 29, 34, 41, 45, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
        #floats  = (1, 2.802596928649634e-45, 1, 2.802596928649634e-45, 2.802596928649634e-45, 2.802596928649634e-45, 1, 2.802596928649634e-45, 1.5414283107572988e-44, 1.6815581571897805e-44, 2.2420775429197073e-44, 2.942726775082116e-44, 3.5032461608120427e-44, 2.802596928649634e-45, 3, 3.923635700109488e-44, 4.0637655465419695e-44, 4.764414778704378e-44, 5.74532370373175e-44, 6.305843089461677e-44, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
        #"""
        #self.show_data(data[n:])

    def add_op2_element(self, elem):
        """checks that eids are positive and that -1 node ids become None"""
        op2: OP2Geom = self.op2
        if elem.eid <= 0:
            op2.log.debug(str(elem))
            raise ValueError(elem)
            #return

        #if elem.eid > 100000000:
            #raise RuntimeError('bad parsing...elem:\n%s' % elem)

        if elem.type in ['CTRIA6', 'CQUAD8']:
            for nid in elem.nodes:
                if nid == -1:
                    nid = None
        else:
            for nid in elem.nodes:
                if nid == -1:
                    assert nid > 0, elem
        op2._add_methods._add_element_object(elem, allow_overwrites=False)
        #print(str(elem)[:-1])

# 1-AEROQ4 (???)
# AEROT3   (???)
# 1-BEAMAERO (1701,17,0)
    def read_caabsf(self, data: bytes, n: int) -> int:
        """2-CAABSF (2708,27,59)

        Word Name Type Description
        1 EID I Element identification number
        2 PID I Property identification number
        3 G1  I Grid point 1 identification number
        4 G2  I Grid point 2 identification number
        5 G3  I Grid point 3 identification number
        6 G4  I Grid point 4 identification number
        """
        return self._run_4nodes(CAABSF, data, n)

# 3-CAXIF2 (2108,21,224)
# 4-CAXIF3 (2208,22,225)
# 5-CAXIF4 (2308,23,226)

    def read_cbar(self, data: bytes, n: int) -> int:
        r"""
        CBAR(2408,24,180) - the marker for Record 8

        MSC/NX
        Word Name Type Description
        1 EID    I  Element identification number
        2 PID    I  Property identification number
        3 GA     I  Grid point identification number at end A
        4 GB     I  Grid point identification number at end B

        F=0* XYZ option -- basic coordinate system
           5 X1 RS  T1 component of orientation vector from GA
           6 X2 RS  T2 component of orientation vector from GA
           7 X3 RS  T3 component of orientation vector from GA
           8 FE  I  Orientation vector flag (encoded)
        F=1* XYZ option -- global coordinate system
           5 X1 RS  T1 component of orientation vector from GA
           6 X2 RS  T2 component of orientation vector from GA
           7 X3 RS  T3 component of orientation vector from GA
           8 FE  I   Orientation vector flag (encoded)
        F=2* Grid option
           5 GO I Grid point identification number at end of orientation vector
           6 UNDEF(2) none Not used
           8 FE I Orientation vector flag (encoded)
        *F = FE bit-wise AND with 3
        End F

        9  PA   I Pin flags for end A
        10 PB   I Pin flags for end B
        11 W1A RS T1 component of offset vector from GA
        12 W2A RS T2 component of offset vector from GA
        13 W3A RS T3 component of offset vector from GA
        14 W1B RS T1 component of offset vector from GB
        15 W2B RS T2 component of offset vector from GB
        16 W3B RS T3 component of offset vector from GB
        F:\work\pyNastran\pyNastran\master2\pyNastran\bdf\test\nx_spike\out_sebload1.op2
        """
        op2: OP2Geom = self.op2
        ntotal = 64 * self.factor # 16*4
        fe1 = 28 * self.factor
        fe2 = 32 * self.factor
        nelements = (len(data) - n) // ntotal
        struct_i = op2.struct_i if self.size == 4 else self.struct_q
        s1 = Struct(mapfmt(op2._endian + b'4i3f3i6f', self.size))
        #s2 = Struct(op2._endian + b'4i3f3i6f')
        s2 = s1
        s3 = Struct(mapfmt(op2._endian + b'7ii2i6f', self.size))
        for unused_i in range(nelements):
            edata = data[n:n + ntotal]
            # we need this flag before we can figure out how to read f
            # per DMAP: F = FE bit-wise AND with 3
            fe, = struct_i.unpack(edata[fe1:fe2])
            f = fe & 3
            #if f not in [0, 1, 2]: f = 0

            # CBAR    EID     PID     GA      GB      X1       X2     X3        OFFT
            #          PA      PB     W1A     W2A     W3A      W1B    W2B       W3B
            # CBAR    401     3       2217    81      .2769987-.931498-.235759  B
            # 'fe = 65; fe&3=1'
            # ints    = (401, 3, 2217, 81, 0.277, -0.931498, -0.235759, 65, 0, 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
            #print(f'eid={eid} fe={fe} f={f}')
            if f == 0:
                # nodes defined in cid=0
                # XYZ option -- basic coordinate system
                out = s1.unpack(edata)
                (eid, pid, ga, gb, x1, x2, x3, _f, pa, pb,
                 w1a, w2a, w3a, w1b, w2b, w3b) = out
                data_in = [[eid, pid, ga, gb, pa, pb, w1a, w2a, w3a, w1b, w2b, w3b],
                           [f, x1, x2, x3]]
                #print(f'eid={eid} fe,f={fe},{f} x=[{x1},{x2},{x3}]')
            elif f == 1:
                # nodes defined in a coordinate system
                #XYZ option -- global coordinate system
                out = s2.unpack(edata)
                (eid, pid, ga, gb, x1, x2, x3, _f, pa, pb,
                 w1a, w2a, w3a, w1b, w2b, w3b) = out
                data_in = [[eid, pid, ga, gb, pa, pb, w1a, w2a, w3a, w1b, w2b, w3b],
                           [f, x1, x2, x3]]
                #print(f'eid={eid} fe,f={fe},{f} x=[{x1},{x2},{x3}]')
            elif f == 2:
                # Grid option
                out = s3.unpack(edata)
                (eid, pid, ga, gb, g0, unused_junk1, unused_junk2, _f, pa,
                 pb, w1a, w2a, w3a, w1b, w2b, w3b) = out
                data_in = [[eid, pid, ga, gb, pa, pb, w1a,
                            w2a, w3a, w1b, w2b, w3b], [f, g0]]
                #print(f'eid={eid} fe,f={fe},{f} g0={g0}')
            else:
                raise RuntimeError('invalid f value...f=%s' % (f))
            elem = CBAR.add_op2_data(data_in)
            try:
                offt = BAR_FE_MAP[fe]
            except Exception:
                print(elem)
                raise
            elem.offt = offt
            #print(f'eid={eid} f={f} fe={fe} offt={offt}')
            #assert f == fe, 'f=%s type(f)=%s fe=%s\n%s' % (f, type(f), fe, elem)

            self.add_op2_element(elem)
            n += ntotal
        op2.card_count['CBAR'] = nelements
        return n

    def read_cbarao(self, data: bytes, n: int) -> int:
        """
        CBARAO(4001,40,275) - the marker for Record 9

        1 EID   I Element identification number
        2 SCALE I Scale of Xi values
        3 X1 RS 1st intermediate station for data recovery
        4 X2 RS 2nd intermediate station for data recovery
        5 X3 RS 3rd intermediate station for data recovery
        6 X4 RS 4th intermediate station for data recovery
        7 X5 RS 5th intermediate station for data recovery
        8 X6 RS 6th intermediate station for data recovery
        9 UNDEF none Not used
        """
        op2: OP2Geom = self.op2
        nelements = (len(data) - n) // 36
        s = Struct(op2._endian + b'2i7f')
        for unused_i in range(nelements):
            edata = data[n:n + 36]  # 9*4
            out = s.unpack(edata)
            if op2.is_debug_file:
                op2.binary_debug.write('  CBARAO=%s\n' % str(out))
            (eid, scale, x1, x2, x3, x4, x5, x6, unused_null) = out
            if scale == 2:
                scale = 'FR'
            else:
                NotImplementedError('CBARAO scale=%r; 2=FR' % scale)
            x = [x1, x2, x3, x4, x5, x6]
            op2.add_cbarao(eid, scale, x, comment='')
            n += 36
        op2.card_count['CBARAO'] = nelements
        return n

    def read_cbeam3(self, data: bytes, n: int) -> int:
        """Common method for reading CBEAM3s"""
        card_name = 'CBEAM3'
        card_obj = CBEAM3
        methods = {
            104 : self._read_cbeam3_104,
            108 : self._read_cbeam3_108,
        }
        try:
            n = self._read_double_card(
                card_name, card_obj, self.add_op2_element,
                methods, data, n)
        except DoubleCardError:
            raise
        return n

    def _read_cbeam3_104(self, card_obj, data: bytes, n: int) -> int:
        """
        CBEAM3(15418, 154, 610)

        $       eid     pid     ga      gb      gc      g0
        CBEAM3  1       1       1       2       21      100

                  eid pid ga gb   gc  5  6  7  g0   x2 x3 ?
        ints    = (1,  1, 1, 2,   21, 0, 0, 0, 100, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                   2,  1, 2, 3,   22, 0, 0, 0, 100, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                   3,  1, 3, 4,   23, 0, 0, 0, 100, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                   4,  1, 4, 5,   24, 0, 0, 0, 100, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                   5,  1, 5, 6,   25, 0, 0, 0, 100, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                   6,  1, 6, 7,   26, 0, 0, 0, 100, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                   7,  1, 1, 8,   31, 0, 0, 0, 100, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                   8,  1, 8, 9,   32, 0, 0, 0, 100, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                   9,  1, 9, 10,  33, 0, 0, 0, 100, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                   10, 1, 10, 11, 34, 0, 0, 0, 100, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                   11, 1, 11, 12, 35, 0, 0, 0, 100, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                   12, 1, 12, 13, 36, 0, 0, 0, 100, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)

               eid, pid, ga, gb  gc,  x1  x2  x3
        CBEAM3, 11,   2,  2,  4, ,   1.0,0.0,0.0
                    eid pid ga gb, gc  5    6    7    8/x1 9/x2 10/x3
          ints    = (11, 2, 2, 4, 0,   0,   0,   0,   1.0, 0,   0,   1, 0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0)
          floats  = (11, 2, 2, 4, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

          CBEAM3  11      2201    1       3       2       1001                    +
          +                                                                       +
          +                                       101     103     102
          sa, sb, sc = (101, 103, 102)

                    eid  pid  ga gb, gc  sa   sb   sc    8/x1 9/x2 10/x3
          ints    = (11, 2201, 1, 3, 2, 101, 103, 102, 1001, 0,   0,   2, 0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0)
          floats  = (11, 2201, 1, 3, 2, 101, 103, 102, 1001, 0.0, 0.0, 2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

          CBEAM3  13      11      3       4       7        10
          $       W1A     W2A     W3A     W1B     W2B     W3B     W1C     W2C
          +       0.0     0.0     1.0     0.0     0.0     1.0     0.0     0.0     +
          $       W3C
          +       1.0
        """
        op2: OP2Geom = self.op2
        #self.show_data(data, types='if')
        ntotal = 104 * self.factor  # 26*4
        nelements = (len(data) - n) // ntotal
        structi = Struct(mapfmt(op2._endian + b'8i 3i i 9f 5i', self.size))
        structf = Struct(mapfmt(op2._endian + b'8i 3f i 9f 5i', self.size))

        elements = []
        for unused_i in range(nelements):
            datai = data[n:n+ntotal]
            out = structi.unpack(datai)
            #self.show_data(datai, types='if')
            (eid, pid, ga, gb, gc, sa, sb, sc,
             g0_x1, x2, x3, flag,
             w1a, w2a, w3a, w1b, w2b, w3b, w1c, w2c, w3c,
             *other) = out
            #print([eid, pid], [ga, gb, gc], [sa, sb, sc], [g0_x1, x2, x3, flag],
                  #[w1a, w2a, w3a], [w1b, w2b, w3b], [w1c, w2c, w3c],
                  #*other)
            if flag == 1:
                (eid, pid, ga, gb, gc, sa, sb, sc,
                 x1, x2, x3, flag,
                 w1a, w2a, w3a, w1b, w2b, w3b, w1c, w2c, w3c,
                 *other) = structf.unpack(datai)
                x = [x1, x2, x3]
                g0 = None
            elif flag == 2:
                g0 = g0_x1
                x = None
            else:
                raise NotImplementedError(flag)

            assert ga > 0, (ga, gb, gc)
            assert gb > 0, (ga, gb, gc)
            assert gc >= 0, (ga, gb, gc)
            if gc == 0:
                gc = None
            #print(eid, pid, ga, gb, gc, sa, sb, sc, g0, sum(other))
            #print(other)
            assert sum(other) == 0, other # self.show_data(datai, types='if')

            nids = [ga, gb, gc]
            wa = [w1a, w2a, w3a]
            wb = [w1b, w2b, w3b]
            wc = [w1c, w2c, w3c]
            tw = None # TWA TWB TWC
            s = [sa, sb, sc]
            elem = CBEAM3(eid, pid, nids, x, g0, wa, wb, wc, tw, s)
            assert eid > 0, elem.get_stats()
            assert pid > 0, elem.get_stats()
            elements.append(elem)
            n += ntotal
        #self.show_data(data[n:])
        return n, elements

    def _read_cbeam3_108(self, card_obj, data: bytes, n: int) -> int:
        """
        CBEAM3(15418, 154, 610)

        $       eid     pid     ga      gb      gc      g0
        CBEAM3  1       1       1       2       21      100

                  eid pid ga gb   gc  5  6  7  g0   x2 x3 ?
        ints    = (1,  1, 1, 2,   21, 0, 0, 0, 100, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                   2,  1, 2, 3,   22, 0, 0, 0, 100, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                   3,  1, 3, 4,   23, 0, 0, 0, 100, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                   4,  1, 4, 5,   24, 0, 0, 0, 100, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                   5,  1, 5, 6,   25, 0, 0, 0, 100, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                   6,  1, 6, 7,   26, 0, 0, 0, 100, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                   7,  1, 1, 8,   31, 0, 0, 0, 100, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                   8,  1, 8, 9,   32, 0, 0, 0, 100, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                   9,  1, 9, 10,  33, 0, 0, 0, 100, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                   10, 1, 10, 11, 34, 0, 0, 0, 100, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                   11, 1, 11, 12, 35, 0, 0, 0, 100, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                   12, 1, 12, 13, 36, 0, 0, 0, 100, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)

               eid, pid, ga, gb  gc,  x1  x2  x3
        CBEAM3, 11,   2,  2,  4, ,   1.0,0.0,0.0
                    eid pid ga gb, gc  5    6    7    8/x1 9/x2 10/x3
          ints    = (11, 2, 2, 4, 0,   0,   0,   0,   1.0, 0,   0,   1, 0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0)
          floats  = (11, 2, 2, 4, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

          CBEAM3  11      2201    1       3       2       1001                    +
          +                                                                       +
          +                                       101     103     102
          sa, sb, sc = (101, 103, 102)

                    eid  pid  ga gb, gc  sa   sb   sc    8/x1 9/x2 10/x3
          ints    = (11, 2201, 1, 3, 2, 101, 103, 102, 1001, 0,   0,   2, 0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0)
          floats  = (11, 2201, 1, 3, 2, 101, 103, 102, 1001, 0.0, 0.0, 2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

          CBEAM3  13      11      3       4       7        10
          $       W1A     W2A     W3A     W1B     W2B     W3B     W1C     W2C
          +       0.0     0.0     1.0     0.0     0.0     1.0     0.0     0.0     +
          $       W3C
          +       1.0
                                                                         wa             wb             wc
                    eid  pid ga gb,gc sa   sb   sc   x1  x2  x3    flag [?    ?    ?]  [?    ?    ?]  [?    ?    ?]
          ints    = (13, 11, 3, 4, 7, 0,   0,   0,   10, 0,   0,   2,   0,   0,   1.0, 0,   0,   1.0, 0,   0,   1.0, 0, 0, 0, 0, 0)
          floats  = (13, 11, 3, 4, 7, 0.0, 0.0, 0.0, 10, 0.0, 0.0, 2,   0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0)

                    eid   pid    ga    gb,     gc   sa   sb   sc    x1     x2 x3 flag [?    ?    ?]  [?    ?    ?]  [?    ?    ?]
        ints    = (2901,  2901,  2901,  2902,  2903, 0,   0,   0,   10.0,    0, 0, 1, 0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
                   3901,  2901,  3901,  3903,  3904, 0,   0,   0,   10.0,    0, 0, 1, 0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
                   3902,  2901,  3903,  3902,  3905, 0,   0,   0,   10.0,    0, 0, 1, 0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
                   4901,  2901,  4901,  4904,  4903, 0,   0,   0,   10.0,    0, 0, 1, 0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
                   4902,  2901,  4904,  4906,  4905, 0,   0,   0,   10.0,    0, 0, 1, 0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
                   4903,  2901,  4906,  4902,  4907, 0,   0,   0,   10.0,    0, 0, 1, 0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
                   12901, 2901, 12901, 12902, 12903, 0,   0,   0,      0, 10.0, 0, 1, 0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0)
        floats  = (2901, 2901, 2901, 2902, 2903, 0.0, 0.0, 0.0, 10.0, 0.0, 0.0,    1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                   3901, 2901, 3901, 3903, 3904, 0.0, 0.0, 0.0, 10.0, 0.0, 0.0,    1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 5.467866607795436e-42, 2901, 3903, 5.467866607795436e-42, 5.4720705031884107e-42, 0.0, 0.0, 0.0, 10.0, 0.0, 0.0, 1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 6.867763773655928e-42, 2901, 6.867763773655928e-42, 6.871967669048903e-42, 6.870566370584578e-42, 0.0, 0.0, 0.0, 10.0, 0.0, 0.0, 1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 6.869165072120253e-42, 2901, 6.871967669048903e-42, 6.874770265977553e-42, 6.873368967513228e-42, 0.0, 0.0, 0.0, 10.0, 0.0, 0.0, 1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 6.870566370584578e-42, 2901, 6.874770265977553e-42, 6.869165072120253e-42, 6.876171564441877e-42, 0.0, 0.0, 0.0, 10.0, 0.0, 0.0, 1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.8078151488254465e-41, 2901, 1.8078151488254465e-41, 1.807955278671879e-41, 1.8080954085183115e-41, 0.0, 0.0, 0.0, 0.0, 10.0, 0.0, 1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
        """
        op2: OP2Geom = self.op2
        #self.show_data(data, types='if')
        ntotal = 108 * self.factor  # 26*4
        nelements = (len(data) - n) // ntotal
        structi = Struct(mapfmt(op2._endian + b'8i 3i i 9f 5i i', self.size))
        structf = Struct(mapfmt(op2._endian + b'8i 3f i 9f 5i i', self.size))

        elements = []
        for unused_i in range(nelements):
            datai = data[n:n+ntotal]
            out = structi.unpack(datai)
            #self.show_data(datai, types='if')
            (eid, pid, ga, gb, gc, sa, sb, sc,
             g0_x1, x2, x3, flag,
             w1a, w2a, w3a, w1b, w2b, w3b, w1c, w2c, w3c,
             *other) = out
            #print([eid, pid], [ga, gb, gc], [sa, sb, sc], [g0_x1, x2, x3, flag],
                  #[w1a, w2a, w3a], [w1b, w2b, w3b], [w1c, w2c, w3c],
                  #*other)
            if flag == 1:
                (eid, pid, ga, gb, gc, sa, sb, sc,
                 x1, x2, x3, flag,
                 w1a, w2a, w3a, w1b, w2b, w3b, w1c, w2c, w3c,
                 *other) = structf.unpack(datai)
                x = [x1, x2, x3]
                g0 = None
            elif flag == 2:
                g0 = g0_x1
                x = None
            else:
                raise NotImplementedError(flag)
            #print(eid, pid, ga, gb, gc, sa, sb, sc, g0, sum(other))
            #print(other)
            assert sum(other) == 0, other # self.show_data(datai, types='if')

            nids = [ga, gb, gc]
            wa = [w1a, w2a, w3a]
            wb = [w1b, w2b, w3b]
            wc = [w1c, w2c, w3c]
            tw = None # TWA TWB TWC
            s = [sa, sb, sc]
            elem = CBEAM3(eid, pid, nids, x, g0, wa, wb, wc, tw, s)
            assert eid > 0, elem.get_stats()
            assert pid > 0, elem.get_stats()
            elements.append(elem)
            n += ntotal
        #self.show_data(data[n:])
        return n, elements

    def read_cbeam(self, data: bytes, n: int) -> int:
        """CBEAM(5408,54,261) - the marker for Record 10"""
        op2: OP2Geom = self.op2
        size = self.size
        ntotal = 18 * size  # 18*4
        fe1 = 10 * size
        fe2 = 11 * size
        nelements = (len(data) - n) // ntotal
        struct_i = op2.struct_i if size == 4 else self.struct_q
        #print(mapfmt(op2._endian + b'6i3f3i6f', self.size))
        s1 = Struct(mapfmt(op2._endian + b'6i3f3i6f', size))
        s3 = Struct(mapfmt(op2._endian + b'12i6f', size))

        list_warnings = []
        for unused_i in range(nelements):
            edata = data[n:n + ntotal]
            # 23
            fe, = struct_i.unpack(edata[fe1:fe2])
            # per DMAP: F = FE bit-wise AND with 3
            f = fe & 3
            #f = fe
            #while f > 3:
                #f &= 3
            if f == 0:  # basic cid
                out = s1.unpack(edata)
                (eid, pid, ga, gb, sa, sb, x1, x2, x3, fe,
                 pa, pb, w1a, w2a, w3a, w1b, w2b, w3b) = out
                #op2.log.info('CBEAM: eid=%s fe=%s f=%s; basic cid' % (eid, fe, f))

                data_in = [[eid, pid, ga, gb, sa, sb, pa, pb, w1a, w2a, w3a, w1b, w2b, w3b],
                           [f, x1, x2, x3]]
            elif f == 1:  # global cid
                # CBEAM    89616   5       384720  384521  0.      0.     -1.
                out = s1.unpack(edata)
                (eid, pid, ga, gb, sa, sb, x1, x2, x3, fe,
                 pa, pb, w1a, w2a, w3a, w1b, w2b, w3b) = out
                #op2.log.info('CBEAM: eid=%s fe=%s f=%s; global cid' % (eid, fe, f))
                data_in = [[eid, pid, ga, gb, sa, sb, pa, pb, w1a, w2a, w3a, w1b, w2b, w3b],
                           [f, x1, x2, x3]]
            elif f == 2:  # grid option
                out = s3.unpack(edata)
                (eid, pid, ga, gb, sa, sb, g0, xxa, xxb, fe,
                 pa, pb, w1a, w2a, w3a, w1b, w2b, w3b) = out
                #op2.log.info('CBEAM: eid=%s fe=%s f=%s; grid option '
                              #'(g0=%s xxa=%s xxb=%s)' % (eid, fe, f, g0, xxa, xxb))
                if g0 <= 0 or g0 >= 100000000 or xxa != 0 or xxb != 0:
                    # Nastran set this wrong...MasterModelTaxi
                    #CBEAM    621614  2672    900380  900379 .197266 -.978394.0600586
                    #    6


                    f = 1
                    out = s1.unpack(edata)
                    (eid, pid, ga, gb, sa, sb, x1, x2, x3, fe, pa,
                     pb, w1a, w2a, w3a, w1b, w2b, w3b) = out
                    #op2.log.info('CBEAM: eid=%s fe=%s f=%s; global cid' % (eid, fe, f))
                    data_in = [[eid, pid, ga, gb, sa, sb, pa, pb, w1a, w2a, w3a, w1b, w2b, w3b],
                               [f, x1, x2, x3]]
                    #op2.log.info('   (x1=%s x2=%s x3=%s)' % (x1, x2, x3))
                else:
                    data_in = [[eid, pid, ga, gb, sa, sb, pa, pb, w1a, w2a, w3a, w1b, w2b, w3b],
                               [f, g0]]
            else:
                out = s1.unpack(edata)
                (eid, pid, ga, gb, sa, sb, x1, x2, x3, fe,
                 pa, pb, w1a, w2a, w3a, w1b, w2b, w3b) = out

                out = s3.unpack(edata)
                (eid, pid, ga, gb, sa, sb, g0, xxa, xxb, fe,
                 pa, pb, w1a, w2a, w3a, w1b, w2b, w3b) = out
                #op2.show_data(data, types='ifsdq')
                #if g0 == 0:

                list_warnings.append(f'CBEAM: invalid f={f} value (fe={fe}) -> f=1...eid={eid} pid={pid} '
                                     f'nodes=[{ga},{gb}] x=[{x1},{x2},{x3}]; g0={g0} xx=[{xxa},{xxb}]')
                f = 1
            #print(f'eid={eid} fe={fe} f={f}')
            if op2.is_debug_file:
                op2.binary_debug.write('  CBEAM eid=%s f=%s fe=%s %s\n' % (
                    eid, f, fe, str(data_in)))

            elem = CBEAM.add_op2_data(data_in, f)
            try:
                offt = BAR_FE_MAP[fe]
            except Exception:
                print(elem)
                raise
            elem.offt = offt

            elem.validate()
            self.add_op2_element(elem)
            n += ntotal
        if len(list_warnings):
            op2.log.warning('\n'.join(list_warnings))
        op2.card_count['CBEAM'] = nelements
        return n

    def read_cbeamp(self, data: bytes, n: int) -> int:
        """
        CBEAMP(11401,114,9016) - the marker for Record 11
        """
        self.op2.log.info('geom skipping CBEAMP in GEOM2')
        return len(data)

    def read_cbend(self, data: bytes, n: int) -> int:
        """
        CBEND(4601,46,298) - the marker for Record 12

        1 EID I Element identification number
        2 PID I Property identification number
        3 GA  I Grid point End A identification number
        4 GB  I Grid point End B identification number

        F = 0 Z
          5 X1 RS T1 component of orientation vector from GA
          6 X2 RS T2 component of orientation vector from GA
          7 X3 RS T3 component of orientation vector from GA
        8 F    I     Orientation vector flag = 0
        F = 1 XYZ option - global coordinate system
          5 X1 RS T1 component of orientation vector from GA
          6 X2 RS T2 component of orientation vector from GA
          7 X3 RS T3 component of orientation vector from GA
          8 F   I    Orientation vector flag = 1
        F = 2 Grid option
          5 GO       I Grid point ID at end of orientation vector
          6 UNDEF(2)    None
          8 F        I Orientation vector flag = 2
        End F
        9 UNDEF(4) None
        13 GEOM I Element geometry option
        """
        op2: OP2Geom = self.op2
        ntotal = 52 * self.factor # 4*13
        nentries = (len(data) - n) // ntotal
        fstruc = Struct(op2._endian + b'4i 3f 6i')
        istruc = Struct(op2._endian + b'4i 3i 6i')

        for unused_i in range(nentries):
            edata = data[n:n + 52]  # 13*4
            fe, = op2.struct_i.unpack(edata[28:32])
            # per DMAP: F = FE bit-wise AND with 3
            f = fe & 3
            if f == 0:
                out = fstruc.unpack(edata)
                (eid, pid, ga, gb, x1, x2, x3, fe,
                 unused_dunnoa, unused_dunnob, unused_dunnoc, unused_dunnod, geom) = out
                data_in = [[eid, pid, ga, gb, geom],
                           [f, x1, x2, x3]]
            elif f == 1:
                out = fstruc.unpack(edata)
                (eid, pid, ga, gb, x1, x2, x3, fe,
                 unused_dunnoa, unused_dunnob, unused_dunnoc, unused_dunnod, geom) = out
                data_in = [[eid, pid, ga, gb, geom],
                           [f, x1, x2, x3]]
            elif f == 2:
                out = istruc.unpack(edata)
                (eid, pid, ga, gb, g0, unused_junk1, unused_junk2, fe,
                 unused_dunnoa, unused_dunnob, unused_dunnoc, unused_dunnod, geom) = out
                data_in = [[eid, pid, ga, gb, geom],
                           [f, g0]]
            else:
                raise RuntimeError('invalid f value...f=%s' % (f))
            elem = CBEND.add_op2_data(data_in)
            elem.validate()
            assert f == fe, 'f=%s type(f)=%s fe=%s\n%s' % (f, type(f), fe, elem)

            self.add_op2_element(elem)
            n += 52
        self.op2.increase_card_count('CBEND', nentries)
        return n

    def read_cbush(self, data: bytes, n: int) -> int:
        """
        CBUSH(2608,26,60) - the marker for Record 13
        """
        op2: OP2Geom = self.op2
        ntotal = 56 * self.factor
        nelements = (len(data) - n) // ntotal
        struct_obj1 = Struct(mapfmt(op2._endian + b'4i iii i ifi3f', self.size))
        struct_obj2 = Struct(mapfmt(op2._endian + b'4i fff i ifi3f', self.size))
        for unused_i in range(nelements):
            edata = data[n:n + ntotal]  # 14*4
            out = struct_obj1.unpack(edata)
            eid, pid, ga, gb, five, unused_sixi, unused_seven, f, cid, s, ocid, s1, s2, s3 = out
            si = [s1, s2, s3]
            if f == -1: # Use Element CID below for orientation
                x = [None, None, None]
                g0 = None
            elif f in [0, 1]:
                # 5, 6, 7, f
                #0:4 4:8 8:12 12:16 16:20 20:24 24:28 28:32 32:36
                #0   1   2    3     4     5     6     7     8
                #x1, x2, x3, f2 = unpack('3f i', edata[20:36])
                out = struct_obj2.unpack(edata)
                eid, pid, ga, gb, x1, x2, x3, f2, cid, s, ocid, s1, s2, s3 = out

                assert f == f2, 'f=%s f2=%s' % (f, f2)
                x = [x1, x2, x3]
                g0 = None
            elif f == 2:
                x = [None, None, None]
                g0 = five
            else:
                raise RuntimeError('invalid f value...f=%r' % f)
            if cid == -1:
                cid = None
            data_in = [[eid, pid, ga, gb, cid, s, ocid, si], x, g0]

            elem = CBUSH.add_op2_data(data_in, f)
            self.add_op2_element(elem)
            n += ntotal
        op2.card_count['CBUSH'] = nelements
        return n

    def read_cbush1d(self, data: bytes, n: int) -> int:
        """
        CBUSH1D(5608,56,218) - the marker for Record 14

        1 EID  I Element identification number
        2 PID  I Property identification number
        3 G(2) I Grid point identification numbers
        5 CID  I Coordinate system identification number
        6 UNDEF(3) none

        """
        op2: OP2Geom = self.op2
        ntotal = 32 * self.factor # 4*8
        nelements = (len(data) - n) // ntotal
        struct_6i = Struct(mapfmt(op2._endian + b'8i', self.size))
        for unused_i in range(nelements):
            edata = data[n:n + ntotal]
            out = struct_6i.unpack(edata)
            if op2.is_debug_file:
                op2.binary_debug.write('  CBUSH1D=%s\n' % str(out))
            (eid, pid, g1, g2, cid, unused_a, unused_b, unused_c) = out
            if cid == -1:
                cid = None
            op2.add_cbush1d(eid, pid, [g1, g2], cid=cid)
            n += ntotal
        op2.card_count['CBUSH1D'] = nelements
        return n

    def read_cbush2d(self, data: bytes, n: int) -> int:
        """
        CBUSH2D

        1 EID  I Element identification number
        2 PID  I Property identification number
        3 G(2) I Grid point identification numbers
        5 CID  I Coordinate system identification number
        6 UNDEF(3) none

        CBUSH2D  EID  PID  GA   GB  CID PLANE
        CBUSH2D  201  202  101  100  0  XY
        (201, 202, 101, 100, 0, 1, 0, 0)
        """
        op2: OP2Geom = self.op2
        ntotal = 8 * self.size # 4*8
        nelements = (len(data) - n) // ntotal
        struct_6i = Struct(mapfmt(op2._endian + b'8i', self.size))
        for unused_i in range(nelements):
            edata = data[n:n + ntotal]
            out = struct_6i.unpack(edata)
            if op2.is_debug_file:
                op2.binary_debug.write('  CBUSH2D=%s\n' % str(out))
            (eid, pid, n1, n2, cid, plane_int, undef1, undef2) = out
            assert (undef1, undef2) == (0, 0), (undef1, undef2)
            assert plane_int == 1, plane_int
            plane = 'XY'
            assert cid >= 0, cid
            sptid = None
            nids = [n1, n2]
            op2.add_cbush2d(eid, pid, nids, cid=cid, plane=plane, sptid=sptid, comment='')
            n += ntotal
        op2.card_count['CBUSH2D'] = nelements
        return n

    def read_ccone(self, data: bytes, n: int) -> int:
        """
        CCONE(2315,23,0) - the marker for Record 15
        """
        op2: OP2Geom = self.op2
        op2.log.info('geom skipping CCONE in GEOM2')
        if op2.is_debug_file:
            op2.binary_debug.write('geom skipping CCONE in GEOM2\n')
        return len(data)

    def read_cdamp1(self, data: bytes, n: int) -> int:
        """
        CDAMP1(201,2,69) - the marker for Record 16
        """
        op2: OP2Geom = self.op2
        ntotal = 24 * self.factor # 6*4
        nelements = (len(data) - n) // ntotal
        struct_6i = Struct(mapfmt(op2._endian + b'6i', self.size))
        for unused_i in range(nelements):
            edata = data[n:n + ntotal]
            out = struct_6i.unpack(edata)
            if op2.is_debug_file:
                op2.binary_debug.write('  CDAMP1=%s\n' % str(out))
            #(eid, pid, g1, g2, c1, c2) = out
            elem = CDAMP1.add_op2_data(out)
            self.add_op2_element(elem)
            n += ntotal
        op2.card_count['CDAMP1'] = nelements
        return n

    def read_cdamp2(self, data: bytes, n: int) -> int:
        """
        CDAMP2(301,3,70) - the marker for Record 17
        """
        op2: OP2Geom = self.op2
        ntotal = 24 * self.factor # 6*4
        nelements = (len(data) - n) // ntotal
        s = Struct(mapfmt(op2._endian + b'if4i', self.size))
        for unused_i in range(nelements):
            edata = data[n:n + ntotal]
            out = s.unpack(edata)
            if op2.is_debug_file:
                op2.binary_debug.write('  CDAMP2=%s\n' % str(out))
            #(eid, bdamp, g1, g2, c1, c2) = out
            elem = CDAMP2.add_op2_data(out)
            self.add_op2_element(elem)
            n += ntotal
        op2.card_count['CDAMP2'] = nelements
        return n

    def read_cdamp3(self, data: bytes, n: int) -> int:
        """
        CDAMP3(401,4,71) - the marker for Record 18
        """
        op2: OP2Geom = self.op2
        struct_4i = Struct(op2._endian + b'4i')
        nelements = (len(data) - n) // 16
        for unused_i in range(nelements):
            edata = data[n:n + 16]  # 4*4
            out = struct_4i.unpack(edata)
            if op2.is_debug_file:
                op2.binary_debug.write('  CDAMP3=%s\n' % str(out))
            #(eid, pid, s1, s2) = out
            elem = CDAMP3.add_op2_data(out)
            self.add_op2_element(elem)
            n += 16
        op2.card_count['CDAMP3'] = nelements
        return n

    def read_cdamp4(self, data: bytes, n: int) -> int:
        """
        CDAMP4(501,5,72) - the marker for Record 19
        """
        op2: OP2Geom = self.op2
        s = Struct(op2._endian + b'ifii')
        nelements = (len(data) - n) // 16
        for unused_i in range(nelements):
            edata = data[n:n + 16]  # 4*4
            out = s.unpack(edata)
            if op2.is_debug_file:
                op2.binary_debug.write('  CDAMP4=%s\n' % str(out))
            #(eid, bdamp, s1, s2) = out
            elem = CDAMP4.add_op2_data(out)
            self.add_op2_element(elem)
            n += 16
        op2.card_count['CDAMP4'] = nelements
        return n

    def read_cdamp5(self, data: bytes, n: int) -> int:
        """
        CDAMP5(10608,106,404) - the marker for Record 20
        """
        op2: OP2Geom = self.op2
        s = Struct(op2._endian + b'4i')
        nelements = (len(data) - n) // 16
        for unused_i in range(nelements):
            edata = data[n:n + 16]  # 4*4
            out = s.unpack(edata)
            if op2.is_debug_file:
                op2.binary_debug.write('  CDAMP5=%s\n' % str(out))
            #(eid, pid, s1, s2) = out
            elem = CDAMP5.add_op2_data(out)
            self.add_op2_element(elem)
            n += 16
        op2.card_count['CDAMP5'] = nelements
        return n

# CDUM2
# CDUM3
# CDUM4
# CDUM5
# CDUM6
# CDUM7

    def read_cdum8(self, data: bytes, n: int) -> int:
        self.op2.log.info('geom skipping CDUM9 in GEOM2')
        #ints = np.frombuffer(data[n:], dtype='int32').copy()
        #print('CDUM8', ints)
        return n

    def read_cdum9(self, data: bytes, n: int) -> int:
        self.op2.log.info('geom skipping CDUM9 in GEOM2')
        #ints = np.frombuffer(data[n:], dtype='int32').copy()
        #print('CDUM9', ints)
        return n

    def read_celas1(self, data: bytes, n: int) -> int:
        """
        CELAS1(601,6,73) - the marker for Record 29
        """
        op2: OP2Geom = self.op2
        ntotal = 24 * self.factor  # 6*4
        struct_4i = Struct(mapfmt(op2._endian + b'6i', self.size))
        nelements = (len(data) - n) // ntotal
        for unused_i in range(nelements):
            edata = data[n:n+ntotal]
            out = struct_4i.unpack(edata)
            if op2.is_debug_file:
                op2.binary_debug.write('  CELAS1=%s\n' % str(out))
            #(eid, pid, g1, g2, c1, c2) = out
            elem = CELAS1.add_op2_data(out)
            self.add_op2_element(elem)
            n += ntotal
        op2.card_count['CELAS1'] = nelements
        return n

    def read_celas2(self, data: bytes, n: int) -> int:
        """
        CELAS2(701,7,74) - the marker for Record 30
        """
        op2: OP2Geom = self.op2
        s1 = Struct(mapfmt(op2._endian + b'if4iff', self.size))
        ntotal = 32 * self.factor # 8*4
        nelements = (len(data) - n) // ntotal
        for unused_i in range(nelements):
            edata = data[n:n+ntotal]
            out = s1.unpack(edata)
            if op2.is_debug_file:
                op2.binary_debug.write('  CELAS2=%s\n' % str(out))
            #(eid, k, g1, g2, c1, c2, ge, s) = out
            elem = CELAS2.add_op2_data(out)
            self.add_op2_element(elem)
            n += ntotal
        op2.card_count['CELAS2'] = nelements
        return n

    def read_celas3(self, data: bytes, n: int) -> int:
        """
        CELAS3(801,8,75) - the marker for Record 31
        """
        op2: OP2Geom = self.op2
        ntotal = 16 * self.factor  # 4*4
        struct_4i = Struct(mapfmt(op2._endian + b'4i', self.size))
        ndatai = len(data) - n
        nelements = ndatai // ntotal
        for unused_i in range(nelements):
            edata = data[n:n+ntotal]
            out = struct_4i.unpack(edata)
            if op2.is_debug_file:
                op2.binary_debug.write('  CELAS3=%s\n' % str(out))
            #(eid, pid, s1, s2) = out
            elem = CELAS3.add_op2_data(out)
            self.add_op2_element(elem)
            n += ntotal
        op2.card_count['CELAS3'] = nelements
        return n

    def read_celas4(self, data: bytes, n: int) -> int:
        """
        CELAS4(901,9,76) - the marker for Record 32
        """
        op2: OP2Geom = self.op2
        ntotal = 16 * self.factor  # 4*4
        s = Struct(mapfmt(op2._endian + b'ifii', self.size))
        ndatai = len(data) - n
        nelements = ndatai // ntotal
        for unused_i in range(nelements):
            edata = data[n:n + ntotal]
            out = s.unpack(edata)
            if op2.is_debug_file:
                op2.binary_debug.write('  CELAS4=%s\n' % str(out))
            #(eid, k, s1, s2) = out
            elem = CELAS4.add_op2_data(out)
            self.add_op2_element(elem)
            n += ntotal
        op2.card_count['CELAS4'] = nelements
        return n

    def read_cfastp(self, data: bytes, n: int) -> int:
        """MSC 2020"""
        op2: OP2Geom = self.op2
        ntotal = 328 * self.factor  # 82*4
        s = Struct(mapfmt(op2._endian + b'40i 8f 8i 26i', self.size))
        ndatai = len(data) - n
        nelements = ndatai // ntotal
        for unused_i in range(nelements):
            edata = data[n:n + ntotal]
            #op2.show_data(edata, types='if')
            out = s.unpack(edata)
            if op2.is_debug_file:
                op2.binary_debug.write('  CFASTP=%s\n' % str(out))
            #(eid, k, s1, s2) = out
            eid, pid, elem_prop_flag, *other = out
            assert eid > 0, eid
            #elem = CELAS4.add_op2_data(out)
            #self.add_op2_element(elem)
            #print(out)
            n += ntotal
        #op2.card_count['CFAST'] = nelements
        #raise RuntimeError('CFASTP')
        return n

    def read_cfast(self, data: bytes, n: int) -> int:
        r"""
        RECORD – CFAST(13801,138,566) - NX

        Word Name Type Description
        1 EID       I Element identification number
        2 PID       I Property identification number
        3 GS        I Spot weld master node identification number GS
        4 FORMAT(C) I Connection format (9=elpat, 10=partpat)
        5 GA        I Identification number of GA
        6 GB        I Identification number of GB
        7–8 UNDEF(2)
        9  GUPPER(8)  I Grid identification numbers of the upper shell
        17 GLOWER(8)  I Grid identification numbers of the lower shell
        25 GUACT(32)  I Unique set of grid IDs of the active shells in upper patch
        57 GLACT(32)  I Unique set of grid IDs of the active shells in lower patch
        89 NUG        I Number of active grids in upper patch
        90 NLG        I Number of active grids in lower patch
        91  GUELE(32) I Grid IDs of the active shells in upper patch
        123 GLELE(32) I Grid IDs of the active shells in lower patch
        155 GHA(12)   RS Coordinates of 4 GHA points
        167 GHB(12)   RS Coordinates of 4 GHB points
        179 TAVG      RS Average shell thickness
        FORMAT=9 ELPAT
        180 EIDUP  I Element ID of upper shell
        181 EIDLOW I Element ID of lower shell
        FORMAT=PARTPAT PARTPAT
        180 PIDUP  I Property ID of upper shell
        181 PIDLOW I Property ID of lower shell
        END FORMAT
        182 TMIN RS Minimum shell thickness
        183 XS   RS X coordinate of spot weld location
        184 YS   RS Y coordinate of spot weld location
        185 ZS   RS Z coordinate of spot weld location
        186 XGA  RS X coordinate of point ga
        187 YGA  RS Y coordinate of point ga
        188 ZGA  RS Z coordinate of point ga
        189 XGB  RS X coordinate of point gb
        190 YGB  RS Y coordinate of point gb
        191 ZGB  RS Z coordinate of point gb
        ints = (
            101, 3, 100, 9, 0,   0,   44,  0,
            9, 14, 13, 8,   0,   0,   0,   0,
            29, 34, 33, 28, 0,   0,   0,   0,
            9, 14, 13, 8,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            29, 34, 33, 28, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            4, 4, 9, 14, 13, 8, 0, 0, 0, 0, 9, 14, 13, 8, 0, 0, 0, 0, 9, 14, 13, 8, 0, 0, 0, 0, 9, 14, 13, 8, 0, 0, 0, 0, 29, 34, 33, 28, 0, 0, 0, 0, 29, 34, 33, 28, 0, 0, 0, 0, 29, 34, 33, 28, 0, 0, 0, 0, 29, 34, 33, 28, 0, 0, 0, 0,
            1077697529, 1065830415, 0, 1077697529, 1073264625, 0,
            1073980423, 1073264625, 0,
            1073980423, 1065830415, 0, 1077697529,
            1065830415, 1036831949, 1077697529, 1073264625, 1036831949, 1073980423, 1073264625, 1036831949, 1073980423, 1065830415, 1036831949, -1.0,
            7, 19,
            0.001, 2.5, 1.5, 0.1, 2.5, 1.5, 0.0, 2.5, 1.5,
            0.1)

        CFAST EID PID TYPE IDA IDB GS GA GB
              XS  YS  ZS
                     eid     pid    type       ida    idb    gs    ga    gb
        CFAST        101       3    ELEM       7      19     100
        floats = (
            101, 3, 100, 9, 0.0, 0.0, 44, 0.0,
            9, 14, 13, 8,   0.0, 0.0, 0.0, 0.0,
            29, 34, 33, 28, 0.0, 0.0, 0.0, 0.0,
            9, 14, 13, 8,   0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            29, 34, 33, 28, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            5.64.605193857299268e-45, 1.2611686178923354e-44, 1.961817850054744e-44, 1.8216880036222622e-44, 1.1210387714598537e-44, 0.0, 0.0, 0.0, 0.0, 1.2611686178923354e-44, 1.961817850054744e-44, 1.8216880036222622e-44, 1.1210387714598537e-44, 0.0, 0.0, 0.0, 0.0, 1.2611686178923354e-44, 1.961817850054744e-44, 1.8216880036222622e-44, 1.1210387714598537e-44, 0.0, 0.0, 0.0, 0.0, 1.2611686178923354e-44, 1.961817850054744e-44, 1.8216880036222622e-44, 1.1210387714598537e-44, 0.0, 0.0, 0.0, 0.0, 4.0637655465419695e-44, 4.764414778704378e-44, 4.624284932271896e-44, 3.923635700109488e-44, 0.0, 0.0, 0.0, 0.0, 4.0637655465419695e-44, 4.764414778704378e-44, 4.624284932271896e-44, 3.923635700109488e-44, 0.0, 0.0, 0.0, 0.0, 4.0637655465419695e-44, 4.764414778704378e-44, 4.624284932271896e-44, 3.923635700109488e-44, 0.0, 0.0, 0.0, 0.0, 4.0637655465419695e-44, 4.764414778704378e-44, 4.624284932271896e-44, 3.923635700109488e-44, 0.0, 0.0, 0.0, 0.0,
            2.9431135654449463, 1.056, 0.0, 2.9431135654449463, 1.943, 0.0,
            2.056, 1.943, 0.0, 2.056, 1.056, 0.0, 2.9431135654449463,
            1.056, 0.1, 2.9431135654449463, 1.943, 0.1, 2.056, 1.943, 0.1, 2.056, 1.056, 0.1, -1.0,
            7, 19,
            0.001, 2.5, 1.5, 0.1, 2.5, 1.5, 0.0, 2.5, 1.5,
            0.1)

        C:\MSC.Software\simcenter_nastran_2019.2\tpl_post2\cfast04.op2
        $            eid      pid   type       ida     idb   gs      ga      gb
        CFAST        101       3    PROP       1       2             100     101
        CFAST        102       3    PROP       1       2     200

        C:\MSC.Software\simcenter_nastran_2019.2\tpl_post2\cfast01.op2
        ints    = (eid=101, pid3, gs=100, format=9, ga=0, gb=0,
                   44, 0,
                   gupper=9, 14,  13, 8,  0, 0, 0, 0,
                   glower=29, 34, 33, 28, 0, 0, 0, 0,
                   # guact
                   5,  10,  9,  4, 15, 20, 19, 14,
                   18, 17, 12, 13,  8,  7, 2,   3,
                   0,   0,  0,  0,  0,  0, 0,   0,
                   0,   0,   0, 0,  0,  0, 0,   0,
                   # glact
                   25, 30, 29, 24, 35, 40, 39, 34,
                   38, 37, 32, 33, 28, 27, 22, 23,
                   0,   0, 0,   0,  0,  0,  0,   0,
                   0,   0, 0,   0,  0,  0,  0,   0,
                   # nug, nlg
                   16, 16,
                   #guele
                   5,  10,  9,  4, 0, 0, 0, 0,
                   15, 20, 19, 14, 0, 0, 0, 0,
                   18, 17, 12, 13, 0, 0, 0, 0,
                   8,   7,  2,  3, 0, 0, 0, 0,
                   # glele
                   25, 30, 29, 24, 0, 0, 0, 0,
                   35, 40, 39, 34, 0, 0, 0, 0,
                   38, 37, 32, 33, 0, 0, 0, 0,
                   28, 27, 22, 23, 0, 0, 0, 0,
                   # gba/ghb
                   3.03178, 0.9682, 0.0, 3.031, 2.031, 0.0,
                   1.968, 2.031, 0.0, 1.968, 0.968, 0.0,
                   1078069239, 1064820772, 1036831949, 1078069239, 1073874935, 1036831949, 1073475602, 1073874935, 1036831949, 1073475602, 1064820772, 1036831949, -1082130432,
                   upper=7, lower=19,
                   981668463, 1075838976, 1069547520, 1036831949, 1075838976, 1069547520, 0, 1075838976, 1069547520, 1036831949)
        floats  = (101, 3, 100, 9, 0, 0,
                   44, 0, 9, 14, 13, 8, 0, 0, 0, 0,
                   29, 34, 33, 28, 0, 0, 0, 0,
                   5, 10, 9, 4, 15, 20, 19, 14, 18, 17, 12, 13, 8, 7, 2,
                   3, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                   25, 30, 29, 24, 35, 40, 39, 34, 38,
                   37, 32, 33, 28, 27, 22, 23,
                   0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                   16, 16, 5, 10, 9, 4, 0.0, 0.0, 0.0, 0.0,
                   15, 20, 19, 14, 0.0, 0.0, 0.0, 0.0,
                   18, 17, 12, 13, 0.0, 0.0, 0.0, 0.0,
                   8, 7, 2, 3, 0.0, 0.0, 0.0, 0.0,
                   25, 30, 29, 24, 0.0, 0.0, 0.0, 0.0,
                   35, 40, 39, 34, 0.0, 0.0, 0.0, 0.0,
                   38, 37, 32, 33, 0.0, 0.0, 0.0, 0.0,
                   28, 27, 22, 23, 0.0, 0.0, 0.0, 0.0,
                   3.03178, 0.9682, 0.0, 3.031, 2.031, 0.0,
                   1.968, 2.031, 0.0, 1.968, 0.968, 0.0, 3.031, 0.968, 0.10, 3.031, 2.031, 0.10, 1.968263864517212, 2.031736135482788, 0.10, 1.968263864517212, 0.9682638645172119, 0.10, -1.0,
                   upper=7, lower=19,
                   tmin=0.10,
                   xs=2.5, 1.5, 0.10,
                   xga=2.5, 1.5, 0.0,
                   xgb=2.5, 1.5, 0.10)
        CFAST        101       3    ELEM       7      19     100

        """
        op2: OP2Geom = self.op2
        size = self.size

        #op2.show_data(data[n:], 'if')
        ntotal = 191 * size  # 191*4
        s = Struct(mapfmt(op2._endian + b'8i 2i 144i 13f 12f 2i 10f', size))
        ndatai = len(data) - n
        nelements = ndatai // ntotal
        assert ndatai % ntotal == 0, f'ndatai={ndatai}'
        for unused_i in range(nelements):
            edata = data[n:n + ntotal]
            #op2.show_data(edata, 'if')

            out = s.unpack(edata)
            if op2.is_debug_file:
                op2.binary_debug.write('  CFAST=%s\n' % str(out))
            eid, pid, gs, elem_grid_flag, ga, gb, *other = out
            gupper = out[8:16]
            glower = out[16:24]
            guact = out[24:56]
            glact = out[56:88]
            nug, nlg = out[88:90]

            guele = out[90:122]
            glele = out[122:154]
            gha = out[154:166]
            ghb = out[166:178]
            tavg, = out[178:179]
            upper, lower = out[179:181]
            tmin, = out[181:182]
            xs, ys, zs = out[182:185]
            xga, yga, zga = out[185:188]
            xgb, ygb, zgb = out[188:191]

            if elem_grid_flag == 9:
                elem_grid_flag = 'ELEM'
            elif elem_grid_flag == 10:
                elem_grid_flag = 'PROP'
                #ida, idb
            else:  # pragma: no cover
                raise NotImplementedError(elem_grid_flag)
            ida = upper
            idb = lower
            #print(out)
            id_gs = None if gs == 0 else gs
    #def add_cfast(self, eid: int, pid: int, Type: str, ida: int, idb: int,
                  #gs=None, ga=None, gb=None,
                  #xs=None, ys=None, zs=None, comment: str='') -> CFAST:
            assert eid > 0, eid
            element = op2.add_cfast(
                eid, pid, elem_grid_flag, ida, idb, gs=id_gs,
                xs=xs, ys=ys, zs=zs)
            #print(element)
            element.validate()
            #elem = CFAST.add_op2_data(out)
            #self.add_op2_element(elem)
            n += ntotal
        op2.card_count['CFAST'] = nelements
        return n
        #raise RuntimeError('CFAST')

    def _read_cfast_msc_nx10(self, data: bytes, n: int) -> int:
        """
        Record 34 -- CFAST(9801,98,506)

        MSC 2005r2 -> MSC 2016

        CFAST(9801,98,506) - the marker for Record 34
        1 EID       I Element identification number
        2 PID       I Property identification number
        3 GS        I Spot weld master node identification numberGS
        4 FORMAT(C) I Connection format (0=gridid)
        5 GA        I ID of GA
        6 GB        I ID of GB
        7 TYPE      I Types of upper and lower elements for
        FORM="GRIDID"
        8 CID       I C
        9 GUPPER(8) I Grid IDs of the upper shell
        FORMAT =0 GRIDID of GBI
          17 GLOWER(8) I
        FORMAT =1 ALIGN (not used)
          17 GLOWER(8) I
        FORMAT =2 ELEMID (not used)
          17 GLOWER(8) I
        FORMAT =9 ELPAT for xyz
          17 XYZ(3) RS
          20 UNDEF(5) none Not used
        FORMAT =10 PARTPAT for xyz
          17 XYZ(3) RS
          20 UNDEF(5) none Not used
        End FORMAT
        25 TAVG RS Average shell thickness
        26 UNDEF(2) none Not used
        28 TMIN RS Minimum shell thickness

        """
        raise RuntimeError('CFAST-28 fields')
        self.op2.log.info('geom skipping CFAST in GEOM2')
        return len(data)

# CFASTP

    def read_cfluid2(self, data: bytes, n: int) -> int:
        """
        CFLUID2(8515,85,209) - the marker for Record 35

        1 EID       I Element identification number
        2 IDF1      I RINGFL point 1 identification number
        3 IDF2      I RINGFL point 2 identification number
        4 RHO      RS Mass density
        5 B        RS Bulk modulus
        6 HARMINDX  I Harmonic index
        """
        op2: OP2Geom = self.op2
        s = Struct(op2._endian + b'3i2fi')
        ntotal = 24 * self.factor
        ndatai = len(data) - n
        nelements = ndatai // ntotal
        for unused_i in range(nelements):
            edata = data[n:n + 24]  # 6*4
            out = s.unpack(edata)
            if op2.is_debug_file:
                op2.binary_debug.write('  CFLUID2=%s\n' % str(out))
            eid, idf1, idf2, rho, bi, harmonic = out
            op2.add_cfluid2(eid, [idf1, idf2], rho, bi, harmonic)
            n += 24
        op2.card_count['CFLUID2'] = nelements
        return n

    def read_cfluid3(self, data: bytes, n: int) -> int:
        """
        CFLUID3(8615,86,210) - the marker for Record 36

        1 EID       I Element identification number
        2 IDF1      I RINGFL point 1 identification number
        3 IDF2      I RINGFL point 2 identification number
        4 IDF3      I RINGFL point 3 identification number
        5 RHO      RS Mass density
        6 B        RS Bulk modulus
        7 HARMINDX  I Harmonic index
        """
        op2: OP2Geom = self.op2
        s = Struct(op2._endian + b'4i2fi')
        ntotal = 28 * self.factor
        ndatai = len(data) - n
        nelements = ndatai // ntotal
        for unused_i in range(nelements):
            edata = data[n:n + 28]  # 7*4
            out = s.unpack(edata)
            if op2.is_debug_file:
                op2.binary_debug.write('  CFLUID3=%s\n' % str(out))
            eid, idf1, idf2, idf3, rho, b, harmonic = out
            op2.add_cfluid3(eid, [idf1, idf2, idf3], rho, b, harmonic)
            n += 28
        op2.card_count['CFLUID3'] = nelements
        return n

    def read_cfluid4(self, data: bytes, n: int) -> int:
        """
        CFLUID4(8715,87,211) - the marker for Record 37

        1 EID       I Element identification number
        2 IDF1      I RINGFL point 1 identification number
        3 IDF2      I RINGFL point 2 identification number
        4 IDF3      I RINGFL point 3 identification number
        5 IDF4      I RINGFL point 4 identification number
        6 RHO      RS Mass density
        7 B        RS Bulk modulus
        8 HARMINDX  I Harmonic index
        """
        op2: OP2Geom = self.op2
        s = Struct(op2._endian + b'5i2fi')
        ntotal = 32 * self.factor
        ndatai = len(data) - n
        nelements = ndatai // ntotal
        for unused_i in range(nelements):
            edata = data[n:n + 32]  # 8*4
            out = s.unpack(edata)
            if op2.is_debug_file:
                op2.binary_debug.write('  CFLUID4=%s\n' % str(out))
            eid, idf1, idf2, idf3, idf4, rho, bi, harmonic = out
            op2.add_cfluid4(eid, [idf1, idf2, idf3, idf4], rho, bi, harmonic)
            n += 32
        op2.card_count['CFLUID4'] = nelements
        return n

    def read_cint(self, data: bytes, n: int) -> int:  # pragma: no cover
        """
        Word Name Type Description
        1 EID    I Element identification number
        2 PID    I Property identification number
        3 PTELC  I Pointer to element identification number
        4 NSEG   I Number of segments
        5 PTSGR  I Pointer to segment displacements
        6 NBOUND I Number of boundaries
        7 BID    I Boundary identification number
        8 NEDGE  I Number of edges
        9 PTBND  I Pointer to boundary identification number
        10 PTBGR I Pointer to boundary grid displacements
        11 PTBED I Pointer to boundary edge displacements
        12 PTBGL I Pointer to boundary grid Lagrange Multipliers
        13 PTBEL I Pointer to boundary edge Lagrange Multipliers
        Words 7 through 13 repeat 6 times
        14 UNDEF(2 ) none
        """
        raise UnsupportedCard('CINT')
        self.op2.log.info('geom skipping CINT in GEOM2')
        # C:\NASA\m4\formats\git\examples\move_tpl\ifcq12p.op2
        # doesn't seem to be a card, more of a general info on the geometry...
        #ints = np.frombuffer(data[n:], dtype=op2.idtype).copy()
        return len(data)

    def read_cgap(self, data: bytes, n: int) -> int:
        """
        CGAP(1908,19,104) - the marker for Record 39
        """
        op2: OP2Geom = self.op2
        ntotal = 36 * self.factor  # 9*4
        s1 = Struct(mapfmt(op2._endian + b'4i3fii', self.size))
        struct_i = op2.struct_i if self.size == 4 else self.struct_q
        f2a = 28 * self.factor
        f2b = 32 * self.factor
        g0a = 16 * self.factor
        g0b = 20 * self.factor

        ndatai = len(data) - n
        nelements = ndatai // ntotal
        for unused_i in range(nelements):
            edata = data[n:n + ntotal]
            out = s1.unpack(edata)
            (eid, pid, ga, gb, x1, x2, x3, f, cid) = out  # f=0,1
            g0 = None
            f2, = struct_i.unpack(edata[f2a:f2b])
            assert f == f2, 'f=%s f2=%s' % (f, f2)
            if f == 2:
                g0, = struct_i.unpack(edata[g0a:g0b])
                x1 = None
                x2 = None
                x3 = None
            else:
                assert f == 1, 'CGAP - f=%r f2=%r' % (f, f2)
                assert f2 == 1, 'CGAP - f=%r f2=%r' % (f, f2)
                if op2.is_debug_file:
                    op2.binary_debug.write('  CGAP eid=%s pid=%s gab0=[%s,%s,%s] x123=[%s,%s,%s] '
                                            'cid=%s f=%s\n' % (eid, pid, ga, gb,
                                                               g0, x1, x2, x3, cid, f))
                #raise NotImplementedError('CGAP - f=%r f2=%r' % (f, f2))
            #print('  CGAP eid=%s pid=%s gab0=[%s,%s,%s] x123=[%s,%s,%s] '
                  #'cid=%s f=%s\n' % (eid, pid, ga, gb,
                                     #g0, x1, x2, x3, cid, f))

            data_in = [eid, pid, ga, gb, g0, x1, x2, x3, cid]
            elem = CGAP.add_op2_data(data_in)
            self.add_op2_element(elem)
            n += ntotal
        op2.card_count['CGAP'] = nelements
        return n

    def read_chacab(self, data: bytes, n: int) -> int:
        """CHACAB(8100,81,381)"""
        return self._run_20nodes(CHACAB, data, n)

    def read_chacbr(self, data: bytes, n: int) -> int:
        """CHACAB(8100,81,381)"""
        return self._run_20nodes(CHACBR, data, n)

    def _run_20nodes(self, element: CHEXA20, data: bytes, n: int) -> int:
        """
        common method for
        Word Name Type Description
        1 EID I Element identification number
        2 PID I Property identification number
        3 G(20) I Grid point identification numbers of connection points
        """
        op2: OP2Geom = self.op2
        nelements = (len(data) - n) // 88
        s = Struct(op2._endian + b'22i')
        #if op2.is_debug_file:
            #op2.binary_debug.write('ndata=%s\n' % (nelements * 44))

        #if op2.is_debug_file:
            #op2.binary_debug.write(f'  {element.type}=(eid, pid, [n1, n2]')

        for unused_i in range(nelements):
            edata = data[n:n + 88]  # 22*4
            out = s.unpack(edata)
            (eid, pid, *nodes) = out
            nodes = list(nodes)
            if op2.is_debug_file:
                op2.binary_debug.write(f'  {element.type}=({eid}, {pid}, {nodes}')

            elem = element(eid, pid, nodes)
            self.add_op2_element(elem)
            n += 88
        #if stop:
            #raise RuntimeError('theta is too large...make the quad wrong')
        op2.card_count[elem.type] = nelements
        return n

# CHACBR

    def read_chbdye(self, data: bytes, n: int) -> int:
        """
        CHBDYE(8308,83,405) - the marker for Record ???
        """
        op2: OP2Geom = self.op2
        ntotal = 28 * self.factor  # 7*4
        ndatai = len(data) - n
        nelements = ndatai // ntotal
        s = Struct(op2._endian + b'7i')
        for unused_i in range(nelements):
            edata = data[n:n+28]
            out = s.unpack(edata)
            (eid, eid2, side, iviewf, iviewb, radmidf, radmidb) = out
            if op2.is_debug_file:
                op2.binary_debug.write('  CHBDYE=%s\n' % str(out))
            #op2.log.debug('  CHBDYE=%s' % str(out))
            data_in = [eid, eid2, side, iviewf, iviewb, radmidf, radmidb]
            elem = CHBDYE.add_op2_data(data_in)
            op2._add_methods._add_thermal_element_object(elem)
            n += ntotal
        op2.card_count['CHBDYE'] = nelements
        return n

    def read_chbdyg(self, data: bytes, n: int) -> int:
        """
        CHBDYG(10808,108,406) - the marker for Record 43
        """
        op2: OP2Geom = self.op2
        ntotal = 64 * self.factor  # 16*4
        ndatai = len(data) - n
        nelements = ndatai // ntotal
        s = Struct(op2._endian + b'16i')
        for unused_i in range(nelements):
            edata = data[n:n+64]
            out = s.unpack(edata)
            (eid, unused_blank, Type, iviewf, iviewb, radmidf, radmidb, unused_blank2,
             g1, g2, g3, g4, g5, g6, g7, g8) = out
            if op2.is_debug_file:
                op2.binary_debug.write('  CHBDYG=%s\n' % str(out))
            #op2.log.debug('  CHBDYG=%s' % str(out))
            data_in = [eid, Type, iviewf, iviewb, radmidf, radmidb,
                       g1, g2, g3, g4, g5, g6, g7, g8]
            elem = CHBDYG.add_op2_data(data_in)
            op2._add_methods._add_thermal_element_object(elem)
            n += ntotal
        op2.card_count['CHBDYG'] = nelements
        return n

    def read_chbdyp(self, data: bytes, n: int) -> int:
        """
        CHBDYP(10908,109,407)
        """
        op2: OP2Geom = self.op2
        ntotal = 60  # 16*4
        s = Struct(op2._endian + b'12i 3f')
        ndatai = len(data) - n
        nelements = ndatai // ntotal
        for unused_i in range(nelements):
            edata = data[n:n+60]
            out = s.unpack(edata)
            (eid, pid, Type, iviewf, iviewb, g1, g2, g0, radmidf, radmidb,
             dislin, ce, e1, e2, e3) = out
            if op2.is_debug_file:
                op2.binary_debug.write('  CHBDYP=%s\n' % str(out))
            #op2.log.debug('  CHBDYP=%s' % str(out))
            data_in = [eid, pid, Type, iviewf, iviewb, g1, g2, g0, radmidf, radmidb,
                       dislin, ce, e1, e2, e3]
            elem = CHBDYP.add_op2_data(data_in)
            self._add_thermal_element_object_safe(elem)
            n += ntotal
        op2.card_count['CHBDYP'] = nelements
        return n

    def _add_thermal_element_object_safe(self, obj):
        op2: OP2Geom = self.op2
        if obj.eid in op2.elements:
            op2.reject_lines.append(obj.write_card(size=16))
        else:
            op2._add_methods._add_element_object(obj)
        #raise RuntimeError('this should be overwritten by the BDF class')

    def read_chexa(self, data: bytes, n: int) -> int:
        """
        CHEXA(7308,73,253) - the marker for Record 45
        """
        op2: OP2Geom = self.op2
        s = Struct(mapfmt(op2._endian + b'22i', self.size))
        ntotal = 88 * self.factor  # 22*4
        nelements = (len(data) - n) // ntotal
        for unused_i in range(nelements):
            edata = data[n:n+ntotal]
            out = s.unpack(edata)
            if op2.is_debug_file:
                op2.binary_debug.write('  CHEXA=%s\n' % str(out))
            (eid, pid, g1, g2, g3, g4, g5, g6, g7, g8, g9, g10,
             g11, g12, g13, g14, g15, g16, g17, g18, g19, g20) = out

            data_in = [eid, pid, g1, g2, g3, g4, g5, g6, g7, g8, ]
            big_nodes = [g9, g10, g11, g12, g13, g14, g15, g16,
                         g17, g18, g19, g20]
            if sum(big_nodes) > 0:
                elem = CHEXA20.add_op2_data(data_in + big_nodes)
            else:
                elem = CHEXA8.add_op2_data(data_in)
            self.add_op2_element(elem)
            n += ntotal
        op2.card_count['CHEXA'] = nelements
        return n

    def read_chexa_cz(self, data: bytes, n: int) -> int:
        """
        CHEXCZ(11801,118,907)
        """
        op2: OP2Geom = self.op2
        s = Struct(mapfmt(op2._endian + b'22i 2i', self.size))
        ntotal = 96 * self.factor  # 24*4
        nelements = (len(data) - n) // ntotal
        for unused_i in range(nelements):
            edata = data[n:n+ntotal]
            out = s.unpack(edata)
            if op2.is_debug_file:
                op2.binary_debug.write('  CHEXCZ=%s\n' % str(out))
            (eid, pid, g1, g2, g3, g4, g5, g6, g7, g8, g9, g10,
             g11, g12, g13, g14, g15, g16, g17, g18, g19, g20,
             dummy1, dummy2) = out
            dummy = (dummy1, dummy2)
            assert dummy == (0, 0), dummy
            big_nodes = [g9, g10, g11, g12, g13, g14, g15, g16,
                         g17, g18, g19, g20]
            #print(eid, pid, big_nodes)

            data_in = [eid, pid, g1, g2, g3, g4, g5, g6, g7, g8, ] + big_nodes
            elem = CHEXCZ.add_op2_data(data_in)
            self.add_op2_element(elem)
            n += ntotal
        op2.card_count['CHEXCZ'] = nelements
        return n

# CHEXA20F
# CHEXAFD
# CHEXAL
# CHEXP
    #def _read_chexp(self, data: bytes, n: int) -> int:
        #"""
        #CHEXP(12001,120,9011) - the marker for Record 50
        #"""
        #op2.log.info('geom skipping CHEXP in GEOM2')
        #if op2.is_debug_file:
            #op2.binary_debug.write('geom skipping CHEXP in GEOM2\n')
        #return len(data)

    def read_chexpr(self, data: bytes, n: int) -> int:
        """
        CHEXPR(7409,74,9991) - the marker for Record 48
        """
        return self.read_chexa(data, n)

    def read_cmass1(self, data: bytes, n: int) -> int:
        """
        CMASS1(1001,10,65) - the marker for Record 51
        """
        op2: OP2Geom = self.op2
        ntotal = 24 * self.factor  # 6*4
        nelements = (len(data) - n) // ntotal
        struct_6i = Struct(mapfmt(op2._endian + b'6i', self.size))
        for unused_i in range(nelements):
            edata = data[n:n + ntotal]
            out = struct_6i.unpack(edata)
            if op2.is_debug_file:
                op2.binary_debug.write('  CMASS1=%s\n' % str(out))
            #(eid, pid, g1, g2, c1, c2) = out
            elem = CMASS1.add_op2_data(out)
            op2._add_methods._add_mass_object(elem)
            n += ntotal
        op2.card_count['CMASS1'] = nelements
        return n

    def read_cmass2(self, data: bytes, n: int) -> int:
        """
        CMASS2(1101,11,66) - the marker for Record 52
        """
        op2: OP2Geom = self.op2
        ntotal = 24 * self.factor  # 6*4
        nelements = (len(data) - n) // ntotal
        s = Struct(mapfmt(op2._endian + b'if4i', self.size))
        for unused_i in range(nelements):
            edata = data[n:n + ntotal]
            out = s.unpack(edata)
            if op2.is_debug_file:
                op2.binary_debug.write('  CMASS2=%s\n' % str(out))
            #(eid, m, g1, g2, c1, c2) = out
            elem = CMASS2.add_op2_data(out)
            op2._add_methods._add_mass_object(elem)
            n += ntotal
        op2.card_count['CMASS2'] = nelements
        return n

    def read_cmass3(self, data: bytes, n: int) -> int:
        """
        CMASS3(1201,12,67) - the marker for Record 53
        """
        op2: OP2Geom = self.op2
        struct_4i = Struct(op2._endian + b'4i')
        nelements = (len(data) - n) // 16
        for unused_i in range(nelements):
            edata = data[n:n + 16]  # 4*4
            out = struct_4i.unpack(edata)
            if op2.is_debug_file:
                op2.binary_debug.write('  CMASS3=%s\n' % str(out))
            #(eid, pid, s1, s2) = out
            elem = CMASS3.add_op2_data(out)
            op2._add_methods._add_mass_object(elem)
            n += 16
        op2.card_count['CMASS3'] = nelements
        return n

    def read_cmass4(self, data: bytes, n: int) -> int:
        """
        CMASS4(1301,13,68) - the marker for Record 54
        """
        op2: OP2Geom = self.op2
        nelements = (len(data) - n) // 16
        struct_if2i = Struct(op2._endian + b'ifii')
        for unused_i in range(nelements):
            edata = data[n:n + 16]  # 4*4
            out = struct_if2i.unpack(edata)
            #(eid, m,s 1, s2) = out
            elem = CMASS4.add_op2_data(out)
            op2._add_methods._add_mass_object(elem)
            n += 16
        op2.card_count['CMASS4'] = nelements
        return n

    def read_cmfree(self, data: bytes, n: int) -> int:
        """
        CMFREE(2508,25,0) - the marker for Record 55

        1 EID  I Element identification number
        2   S  I
        3  S2  I
        4   Y RS
        5   N  I
        """
        op2: OP2Geom = self.op2
        assert n == 12, n
        nelements = (len(data) - n) // 20
        assert (len(data) - n) % 20 == 0
        struct_3ifi = Struct(op2._endian + b'3ifi')
        for unused_i in range(nelements):
            edata = data[n:n + 20]  # 5*4
            out = struct_3ifi.unpack(edata)
            eid, s, s2, y, ncm = out
            op2.add_cmfree(eid, s, s2, y, ncm)
            n += 20
        op2.card_count['CMFREE'] = nelements
        return n

    def read_conm1(self, data: bytes, n: int) -> int:
        """
        CONM1(1401,14,63) - the marker for Record 56
        """
        op2: OP2Geom = self.op2
        ntotal = 96 * self.factor   # 24*4
        s = Struct(mapfmt(op2._endian + b'3i21f', self.size))
        nelements = (len(data) - n) // ntotal
        for unused_i in range(nelements):
            edata = data[n:n + ntotal]
            out = s.unpack(edata)
            if op2.is_debug_file:
                op2.binary_debug.write('  CONM1=%s\n' % str(out))
            #(eid, g, cid, m1, m2a, m2b, m3a, m3b, m3c, m4a, m4b, m4c, m4d,
             #m5a, m5b, m5c, m5d, m5e, m6a, m6b, m6c, m6d, m6e, m6f) = out
            elem = CONM1.add_op2_data(out)
            op2._add_methods._add_mass_object(elem)
            n += ntotal
        op2.card_count['CONM1'] = nelements
        return n

    def read_conm2(self, data: bytes, n: int) -> int:
        """
        CONM2(1501,15,64) - the marker for Record 57
        """
        op2: OP2Geom = self.op2
        ntotal = 52 * self.factor  # 13*4
        nelements = (len(data) - n) // ntotal
        s = Struct(mapfmt(op2._endian + b'3i10f', self.size))
        for unused_i in range(nelements):
            edata = data[n:n+ntotal]
            out = s.unpack(edata)
            if op2.is_debug_file:
                op2.binary_debug.write('  CONM2=%s\n' % str(out))
            #(eid, g, cid, m, x1, x2, x3, i1, i2a, i2b, i3a, i3b, i3c) = out
            elem = CONM2.add_op2_data(out)
            op2._add_methods._add_mass_object(elem)
            n += ntotal
        op2.card_count['CONM2'] = nelements
        return n

    def read_conrod(self, data: bytes, n: int) -> int:
        """
        CONROD(1601,16,47) - the marker for Record 58
        """
        op2: OP2Geom = self.op2
        ntotal = 32 *self.factor  # 8*4
        s = Struct(mapfmt(op2._endian + b'4i4f', self.size))
        nelements = (len(data) - n) // ntotal
        for unused_i in range(nelements):
            edata = data[n:n+ntotal]
            out = s.unpack(edata)
            if op2.is_debug_file:
                op2.binary_debug.write('  CONROD=%s\n' % str(out))
            #(eid, n1, n2, mid, a, j, c, nsm) = out
            elem = CONROD.add_op2_data(out)
            self.add_op2_element(elem)
            n += ntotal
        op2.card_count['CONROD'] = nelements
        return n

    def read_conv(self, data: bytes, n: int) -> int:
        """
        The CONV card is different between MSC and NX Nastran.
        The MSC version is 8 fields longer.
        """
        n0 = n
        assert self.factor == 1, self.factor
        op2: OP2Geom = self.op2
        if op2.is_nx:
            try:
                n, elements = self._read_conv_nx(data, n)
            except (AssertionError, MixedVersionCard):
                n, elements = self._read_conv_msc(data, n0)
        else:
            try:
                n, elements = self._read_conv_msc(data, n)
            except (AssertionError, MixedVersionCard):
                n, elements = self._read_conv_nx(data, n0)

        nelements = len(elements)
        for elem in elements:
            op2._add_methods._add_thermal_bc_object(elem, elem.eid)
        op2.card_count['CONV'] = nelements
        assert n == len(data), f'ndata={len(data)} n={n}'
        return n

    def _read_split_card(self, data, n, read1, read2, card_name, card_obj, add_method):
        """
        generalization of multi read methods for different
        versions of MSC Nastran
        """
        op2: OP2Geom = self.op2
        n0 = n
        try:
            n, elements = read1(card_obj, data, n)
        except AssertionError:
            op2.log.info(f'AssertionError...try again reading {card_name!r}')
            n, elements = read2(card_obj, data, n0)

        nelements = len(elements)
        for elem in elements:
            add_method(elem)
        op2.card_count[card_name] = nelements
        return n

    def _read_dual_card(self, data, n, nx_read, msc_read, card_name, add_method) -> int:
        """
        generalization of multi read methods (MSC, NX)
        """
        n, elements = self._read_dual_card_load(
            data, n, nx_read, msc_read)
        nelements = len(elements)
        for elem in elements:
            add_method(elem)
        self.op2.card_count[card_name] = nelements
        return n

    def _read_dual_card_load(self, data, n,
                             nx_read, msc_read) -> tuple[int, list[Any]]:
        n0 = n
        op2: OP2Geom = self.op2
        if op2.is_nx:
            try:
                n, elements = nx_read(data, n0)
            except (AssertionError, MixedVersionCard):
                #raise
                n, elements = msc_read(data, n0)
        else:
            try:
                n, elements = msc_read(data, n0)
            except (AssertionError, MixedVersionCard):
                #raise
                n, elements = nx_read(data, n0)
        assert n is not None
        return n, elements

    def _read_conv_nx(self, data: bytes, n: int) -> int:
        """
        CONV(12701,127,408) - the marker for Record 59
        """
        op2: OP2Geom = self.op2
        ntotal = 48 * self.factor  # 12*4
        s = Struct(op2._endian + b'4i 8i')
        nelements = (len(data) - n) // ntotal
        assert (len(data) - n) % ntotal == 0
        elements = []
        for unused_i in range(nelements):
            edata = data[n:n+ntotal]
            out = s.unpack(edata)
            if op2.is_debug_file:
                op2.binary_debug.write('  CONV=%s; len=%s\n' % (str(out), len(out)))
            (eid, pcon_id, flmnd, cntrlnd,
             ta1, ta2, ta3, ta4, ta5, ta6, ta7, ta8) = out
            assert eid > 0, out
            #assert eid > 0, out

            #ta = [ta1, ta2, ta3, ta5, ta6, ta7, ta8]
            weights = [None] * 8
            data_in = [eid, pcon_id, flmnd, cntrlnd,
                       [ta1, ta2, ta3, ta4, ta5, ta6, ta7, ta8],
                       weights]

            elem = CONV.add_op2_data(data_in)
            elements.append(elem)
            n += ntotal
        return n, elements

    def _read_conv_msc(self, data: bytes, n: int) -> int:
        """
        CONV(12701,127,408) - the marker for Record 60
        """
        op2: OP2Geom = self.op2
        ntotal = 80 * self.factor  # 20*4
        s = Struct(op2._endian + b'12i 8f')
        nelements = (len(data) - n) // ntotal
        assert (len(data) - n) % ntotal == 0
        elements = []
        for unused_i in range(nelements):
            edata = data[n:n+80]
            out = s.unpack(edata)
            if op2.is_debug_file:
                op2.binary_debug.write('  CONV=%s; len=%s\n' % (str(out), len(out)))
            (eid, pcon_id, flmnd, cntrlnd,
             # TODO: why is ta4 and wt4 unused?
             ta1, ta2, ta3, unused_ta4, ta5, ta6, ta7, ta8,
             wt1, wt2, wt3, unused_wt4, wt5, wt6, wt7, wt8) = out
            assert eid > 0, out
            data_in = [eid, pcon_id, flmnd, cntrlnd,
                       [ta1, ta2, ta3, ta5, ta6, ta7, ta8],
                       [wt1, wt2, wt3, wt5, wt6, wt7, wt8]]
            elem = CONV.add_op2_data(data_in)
            elements.append(elem)
            n += ntotal
        return n, elements

    def read_convm(self, data: bytes, n: int) -> int:
        """
        CONVM(8908,89,422) - the marker for Record 60

        MSC
        1 EID I Element identification number
        2 PCONID I Convection property identification number
        3 FLMND I Point for film convection fluid property temperature
        4 CNTMDOT I Control point used for controlling mass flow.
        5 TA I Ambient points used for convection
        Word 5 repeats 2 times

        NX
        1 EID I Element identification number
        2 PCONID I Convection property identification number
        3 FLMND I Point for film convection fluid property temperature
        4 CNTMDOT I Control point used for controlling mass flow.
        5 TA I Ambient points used for convection
        Word 5 repeats 2 times

        [110, 200, 0, 50000, 99999, 99999, 1.0,
        111, 200, 0, 50000, 99999, 99999, 1.0,
        112, 200, 0, 50000, 99999, 99999, 1.0,
        113, 200, 0, 50000, 99999, 99999, 1.0,
        114, 200, 0, 50000, 99999, 99999, 1.0,
        115, 200, 0, 50000, 99999, 99999, 1.0,
        116, 200, 0, 50000, 99999, 99999, 1.0,
        117, 200, 0, 50000, 99999, 99999, 1.0,
        118, 200, 0, 50000, 99999, 99999, 1.0,
        119, 200, 0, 50000, 99999, 99999, 1.0,
        130, 200, 0, 50000, 99999, 99999, 1.0,
        131, 200, 0, 50000, 99999, 99999, 1.0,
        132, 200, 0, 50000, 99999, 99999, 1.0,
        133, 200, 0, 50000, 99999, 99999, 1.0,
        134, 200, 0, 50000, 99999, 99999, 1.0,
        135, 200, 0, 50000, 99999, 99999, 1.0,
        136, 200, 0, 50000, 99999, 99999, 1.0,
        137, 200, 0, 50000, 99999, 99999, 1.0,
        138, 200, 0, 50000, 99999, 99999, 1.0,
        139, 200, 0, 50000, 99999, 99999, 1.0,
        150, 200, 0, 50000, 99999, 99999, 1.0,
        151, 200, 0, 50000, 99999, 99999, 1.0,
        152, 200, 0, 50000, 99999, 99999, 1.0,
        153, 200, 0, 50000, 99999, 99999, 1.0,
        154, 200, 0, 50000, 99999, 99999, 1.0,
        155, 200, 0, 50000, 99999, 99999, 1.0,
        156, 200, 0, 50000, 99999, 99999, 1.0,
        157, 200, 0, 50000, 99999, 99999, 1.0,
        158, 200, 0, 50000, 99999, 99999, 1.0,
        159, 200, 0, 50000, 99999, 99999, 1.0,
        170, 200, 0, 50000, 99999, 99999, 1.0,
        171, 200, 0, 50000, 99999, 99999, 1.0,
        172, 200, 0, 50000, 99999, 99999, 1.0,
        173, 200, 0, 50000, 99999, 99999, 1.0,
        174, 200, 0, 50000, 99999, 99999, 1.0,
        175, 200, 0, 50000, 99999, 99999, 1.0,
        176, 200, 0, 50000, 99999, 99999, 1.0,
        177, 200, 0, 50000, 99999, 99999, 1.0,
        178, 200, 0, 50000, 99999, 99999, 1.0,
        179, 200, 0, 50000, 99999, 99999, 1.0,
        190, 200, 0, 50000, 99999, 99999, 1.0,
        191, 200, 0, 50000, 99999, 99999, 1.0,
        192, 200, 0, 50000, 99999, 99999, 1.0,
        193, 200, 0, 50000, 99999, 99999, 1.0,
        194, 200, 0, 50000, 99999, 99999, 1.0,
        195, 200, 0, 50000, 99999, 99999, 1.0,
        196, 200, 0, 50000, 99999, 99999, 1.0,
        197, 200, 0, 50000, 99999, 99999, 1.0,
        198, 200, 0, 50000, 99999, 99999, 1.0,
        199, 200, 0, 50000, 99999, 99999, 1.0]

        [8908, 89, 422,
        101, 101, 0, 50000, 99999, 99999, 1.0,
        102, 102, 0, 50001, 99999, 99999, 1.0,
        103, 102, 0, 50001, 99999, 99999, 1.0,
        104, 104, 0, 50002, 99999, 99999, 1.0,
        105, 105, 0, 50003, 99999, 99999, 1.0,
        CONVM EID PCONID FLMND CNTMDOT TA1 TA2
        106, 105, 0,     50003, 99999, 99999, 1.0)
        """
        op2: OP2Geom = self.op2
        #C:\Users\sdoyle\Dropbox\move_tpl\ht15330.op2
        ntotal6 = 24  # 6*4
        ntotal7 = 28  # 7*4
        ndata = len(data)
        nelements6 = (ndata - n) // ntotal6
        nelements7 = (ndata - n) // ntotal7

        is_six = (ndata - n) % ntotal6 == 0
        is_seven = (ndata - n) % ntotal7 == 0
        assert self.factor == 1, self.factor

        if is_six and is_seven:
            try:
                op2.log.warning('CONVM: assuming 6 fields')
                elements, n = self.read_convm6(data, nelements6, n)
                op2.log.debug('CONVM: 6 fields')
            except RuntimeError:  # eid < 0
                op2.log.warning('CONVM: assuming 7 fields')
                elements, n = self.read_convm7(data, nelements7, n)
                op2.log.debug('CONVM: 7 fields')
        elif is_six:
            elements, n = self.read_convm6(data, nelements6, n)
        elif is_seven:
            elements, n = self.read_convm7(data, nelements7, n)
        else:  # pragma: no cover
            raise RuntimeError('CONVM is_six=%s is_seven=%s' % (is_six, is_seven))
        nelements = len(elements)
        for element in elements:
            op2._add_methods._add_thermal_bc_object(element, element.eid)
        op2.card_count['CONVM'] = nelements
        return n

    def read_convm6(self, data, nelements, n):
        op2: OP2Geom = self.op2
        structi = Struct(op2._endian + b'6i')
        elements = []
        for unused_i in range(nelements):
            edata = data[n:n+24]
            out = structi.unpack(edata)
            if op2.is_debug_file:
                op2.binary_debug.write('  CONVM=%s\n' % str(out))
            (eid, pcon_id, flmnd, cntrlnd, ta1, ta2) = out
            #if eid <= 0:
            if eid <= 0 or pcon_id <= 0 or flmnd < 0 or cntrlnd <= 0 or ta1 <= 0 or ta2 <= 0:
                #self.show_data(data, 'if')
                # TODO: I'm not sure that this really has 6 fields...
                raise RuntimeError(f'eid={eid} pconid={pcon_id} flmnd={flmnd} cntrlnd={cntrlnd} ta1={ta1} ta2={ta2} < 0')
            mdot = 0.
            data_in = [eid, pcon_id, flmnd, cntrlnd, ta1, ta2, mdot]
            elem = CONVM.add_op2_data(data_in)
            n += 24
            elements.append(elem)
        return elements, n

    def read_convm7(self, data, nelements, n):
        op2: OP2Geom = self.op2
        structi = Struct(op2._endian + b'6if')
        elements = []
        for unused_i in range(nelements):
            edata = data[n:n+28]
            out = structi.unpack(edata)
            if op2.is_debug_file:
                op2.binary_debug.write('  CONVM=%s\n' % str(out))
            (eid, pcon_id, flmnd, cntrlnd, ta1, ta2, mdot) = out

            if eid <= 0 or pcon_id <= 0 or flmnd < 0 or cntrlnd <= 0 or ta1 <= 0 or ta2 <= 0:
                op2.show_data(data, 'if')
                # TODO: I'm not sure that this really has 7 fields...
                raise RuntimeError(f'eid={eid} pconid={pcon_id} flmnd={flmnd} '
                                   f'cntrlnd={cntrlnd} ta1={ta1} ta2={ta2} < 0')
            data_in = [eid, pcon_id, flmnd, cntrlnd, ta1, ta2, mdot]
            elem = CONVM.add_op2_data(data_in)
            n += 28
            elements.append(elem)
        return elements, n

    def read_cplsts3(self, data: bytes, n: int) -> int:
        """
        RECORD - CPLSTS3(8801,88,984)
        Word Name Type Description
        1 EID  I Element identification number
        2 PID  I Property identification number
        3 G(3) I Grid point identification numbers of connection points
        6 UNDEF None
        7 THETA RS Material property orientation angle or coordinate system ID
        8 UNDEF(4) None
        12 TFLAG I Flag signifying meaning of T(3) values
        13 T(3) RS Membrane thickness of element at grid points
        16 UNDEF None
        """
        op2: OP2Geom = self.op2
        #self.show_data(data[n:], types='if')
        op2.to_nx(' because CPLSTS3 was found')
        ntotal = 64 * self.factor  # 16*4
        struct_16i = Struct(mapfmt(op2._endian + b'6i f 4i i3f i', self.size))
        ndatai = len(data)
        nelements = (ndatai - n) // ntotal
        leftover = (ndatai - n) % ntotal
        assert leftover == 0, leftover
        for unused_i in range(nelements):
            edata = data[n:n + ntotal]
            out = struct_16i.unpack(edata)
            if op2.is_debug_file:
                op2.binary_debug.write('  CPLSTS3=%s\n' % str(out))
            (eid, pid, n1, n2, n3, undef6, theta, undef8, undef9, undef10, undef11,
             tflag, t1, t2, t3, undef16) = out
            #print(eid, pid, (n1, n2, n3), theta,
                  #tflag, t1, t2, t3)
            nids = [n1, n2, n3]
            undefs = (undef6, undef8, undef9, undef10, undef11, undef16)
            assert min(undefs) == 0
            assert max(undefs) == 0
            cplsts3 = op2.add_cplsts3(eid, pid, nids, theta=theta,
                                       tflag=tflag, T1=t1, T2=t2, T3=t3)
            str(cplsts3)
            n += ntotal
        return n

    def read_cplsts4(self, data: bytes, n: int) -> int:
        """
        RECORD - CPLSTS4(8401,84,985)

        Word Name Type Description
        1 EID  I Element identification number
        2 PID  I Property identification number
        3 G(4) I Grid point identification numbers of connection points
        7 THETA RS Material property orientation angle or coordinate system ID
        8 UNDEF(4) None
        12 TFLAG I Flag signifying meaning of T(4) values
        13 T(4) RS Membrane thickness of element at grid points
        """
        op2: OP2Geom = self.op2
        op2.to_nx(' because CPLSTS4-NX was found')
        ntotal = 64 * self.factor  # 16*4
        struct_16i = Struct(mapfmt(op2._endian + b'6i f 4i i4f', self.size))
        ndatai = len(data)
        nelements = (ndatai - n) // ntotal
        leftover = (ndatai - n) % ntotal
        assert leftover == 0, leftover
        for unused_i in range(nelements):
            edata = data[n:n + ntotal]
            out = struct_16i.unpack(edata)
            if op2.is_debug_file:
                op2.binary_debug.write('  CPLSTS4=%s\n' % str(out))
            (eid, pid, n1, n2, n3, n4, theta, undef8, undef9, undef10, undef11,
             tflag, t1, t2, t3, t4) = out
            #print(eid, pid, (n1, n2, n3), theta,
                  #tflag, t1, t2, t3)
            nids = [n1, n2, n3, n4]
            undefs = (undef8, undef9, undef10, undef11)
            assert min(undefs) == 0, undefs
            assert max(undefs) == 0, undefs
            cplsts4 = op2.add_cplsts4(eid, pid, nids, theta=theta,
                                       tflag=tflag, T1=t1, T2=t2, T3=t3, T4=t4)
            str(cplsts4)
            n += ntotal
        return n

    def read_cplsts6(self, data: bytes, n: int) -> int:
        """
        RECORD - CPLSTS6(1801,18,986)
        Word Name Type Description
        1 EID  I Element identification number
        2 PID  I Property identification number
        3 G(6) I Grid point identification numbers of connection points
        9 UNDEF(2) None
        11 THETA RS Material property orientation angle or coordinate system ID
        12 TFLAG I Flag signifying meaning of T(3) values
        13 TC(3) RS Membrane thickness of element at corner grid points
        16 UNDEF(5) None
        21 TM(3) RS Membrane thickness of element at mid-side grid points 24 UNDEF None
        """
        #1728 / 4 = 432
        #432 = 16 * 27
        op2: OP2Geom = self.op2
        op2.to_nx(' because CPLSTS6-NX was found')
        struct_16i = Struct(mapfmt(op2._endian + b'2i 8i fi 4f 4i 4f', self.size))
        ntotal = 96 * self.factor  # 24*4
        #ntotal = 128 * self.factor  # 16*4
        #struct_16i = Struct(mapfmt(op2._endian + b'2i 8i f i4f 4i 4f', self.size))
        ndatai = len(data) - n
        nelements = ndatai // ntotal
        leftover = ndatai % ntotal
        assert leftover == 0, leftover
        for unused_i in range(nelements):
            edata = data[n:n + ntotal]
            #self.show_data(edata, types='if')
            out = struct_16i.unpack(edata)

            if op2.is_debug_file:
                op2.binary_debug.write('  CPLSTS6=%s\n' % str(out))
            (eid, pid, n1, n2, n3, n4, n5, n6, undef7, undef8, theta,
             tflag, t1, t2, t3, zero1, zero2, zero3, zero4, zero5, t4, t5, t6, undef8b) = out
            nids = [n1, n2, n3, n4, n5, n6]
            undef = (zero1, zero2, zero3, zero4, zero5, undef7, undef8, undef8b)
            #print(eid, pid, nids, theta,
                  #(tflag, t1, t2, t3, t4, t5, t6), undef)
            assert min(nids) >= 0, nids
            assert min(t1, t2, t3, t4, t5, t6) == -1.0, (t1, t2, t3, t4, t5, t6)
            assert max(t1, t2, t3, t4, t5, t6) == -1.0, (t1, t2, t3, t4, t5, t6)
            cplsts6 = op2.add_cplsts6(eid, pid, nids, theta=theta,
                                       tflag=tflag,
                                       T1=t1, T2=t2, T3=t3, T4=t4,
                                       T5=t5, T6=t6)
            #print(cplsts6)
            str(cplsts6)
            n += ntotal
        return n

    def read_cplsts8(self, data: bytes, n: int) -> int:
        """
        RECORD - CPLSTS8(3601,36,987)

        Word Name Type Description
        1 EID     I Element identification number
        2 PID     I Property identification number
        3 G(8)    I Grid point identification numbers of connection points
        11 THETA RS Material property orientation angle or coordinate system ID
        12 TFLAG  I Flag signifying meaning of T(4) values
        13 TC(4) RS Membrane thickness of element at corner grid points
        17 UNDEF(4) None
        21 TM(4) RS Membrane

        64:
          ints    = (39, 4, 43, 41, 114, 115, 54, 55, 116, 56, 0, 0,     -1.0, -1.0, -1.0, -1.0)
          floats  = (39, 4, 43, 41, 114, 115, 54, 55, 116, 56, 0.0, 0.0, -1.0, -1.0, -1.0, -1.0)
        """
        op2: OP2Geom = self.op2
        op2.to_nx(' because CPLSTS8 was found')
        struct_16i = Struct(mapfmt(op2._endian + b'2i 8i fi 4f 4i 4f', self.size))
        ntotal = 96 * self.factor  # 24*4
        ndatai = len(data) - n
        nelements = ndatai // ntotal
        leftover = ndatai % ntotal
        assert leftover == 0, leftover
        for unused_i in range(nelements):
            edata = data[n:n + ntotal]
            out = struct_16i.unpack(edata)

            if op2.is_debug_file:
                op2.binary_debug.write('  CPLSTS8=%s\n' % str(out))
            (eid, pid, n1, n2, n3, n4, n5, n6, n7, n8, theta,
             tflag, t1, t2, t3, t4, zero1, zero2, zero3, zero4, t5, t6, t7, t8) = out
            nids = [n1, n2, n3, n4, n5, n6, n7, n8]
            #print(eid, pid, nids, theta,
                  #(tflag, t1, t2, t3, t4, t5, t6, t7, t8))
            assert min(nids) >= 0, nids
            assert min(t5, t6, t7, t8) == -1.0, (t5, t6, t7, t8)
            assert max(t5, t6, t7, t8) == -1.0, (t5, t6, t7, t8)
            cplsts8 = op2.add_cplsts8(eid, pid, nids, theta=theta,
                                       tflag=tflag,
                                       T1=t1, T2=t2, T3=t3, T4=t4,
                                       T5=t5, T6=t6, T7=t7, T8=t8)
            #print(cplsts8)
            str(cplsts8)
            n += ntotal
        return n

    def read_cplstn3(self, data: bytes, n: int) -> int:
        """
        RECORD - CPLSTN3(1701,17,980)

        Word Name Type Description
        1 EID    I Element identification number
        2 PID    I Property identification number
        3 G(3)   I Grid point identification numbers of connection points
        6 THETA RS Material property orientation angle or coordinate system ID
        7 UNDEF(10) None

        """
        op2: OP2Geom = self.op2
        op2.to_nx(' because CPLSTN3 was found')
        struct_16i = Struct(mapfmt(op2._endian + b'2i 3i f 10i', self.size))
        ntotal = 64 * self.factor  # 16*4
        ndatai = len(data) - n
        nelements = ndatai // ntotal
        leftover = ndatai % ntotal
        assert leftover == 0, leftover
        for unused_i in range(nelements):
            edata = data[n:n + ntotal]
            out = struct_16i.unpack(edata)

            if op2.is_debug_file:
                op2.binary_debug.write('  CPLSTN3=%s\n' % str(out))
            (eid, pid, n1, n2, n3, theta, *undef) = out
            nids = [n1, n2, n3]
            assert min(nids) > 0, nids
            assert min(undef) == 0, undef
            assert max(undef) == 0, undef
            cplstn3 = op2.add_cplstn3(eid, pid, nids, theta=theta)
            #print(cplstn3)
            str(cplstn3)
            n += ntotal
        return n

    def read_cplstn4(self, data: bytes, n: int) -> int:
        """
        RECORD - CPLSTN4(5701,57,981)
        Word Name Type Description
        1 EID    I Element identification number
        2 PID    I Property identification number
        3 G(4)   I Grid point identification numbers of connection points
        7 THETA RS Material property orientation angle or coordinate system ID
        8 UNDEF(9) None

        """
        op2: OP2Geom = self.op2
        op2.to_nx(' because CPLSTN4 was found')
        struct_16i = Struct(mapfmt(op2._endian + b'2i 4i f 9i', self.size))
        ntotal = 64 * self.factor  # 16*4
        ndatai = len(data) - n
        nelements = ndatai // ntotal
        leftover = ndatai % ntotal
        assert leftover == 0, leftover
        for unused_i in range(nelements):
            edata = data[n:n + ntotal]
            out = struct_16i.unpack(edata)

            if op2.is_debug_file:
                op2.binary_debug.write('  CPLSTN4=%s\n' % str(out))
            (eid, pid, n1, n2, n3, n4, theta, *undef) = out
            nids = [n1, n2, n3, n4]
            assert min(nids) > 0, nids
            assert min(undef) == 0, undef
            assert max(undef) == 0, undef
            cplstn4 = op2.add_cplstn4(eid, pid, nids, theta=theta)
            #print(cplstn4)
            str(cplstn4)
            n += ntotal
        return n

    def read_cplstn6(self, data: bytes, n: int) -> int:
        """
        RECORD - CPLSTN6(5801,58,982)

        Word Name Type Description
        1 EID  I Element identification number
        2 PID  I Property identification number
        3 G(6) I Grid point identification numbers of connection points
        9 THETA RS Material property orientation angle or coordinate system ID
        10 UNDEF(7) None

        """
        op2: OP2Geom = self.op2
        op2.to_nx(' because CPLSTN6 was found')
        struct_16i = Struct(mapfmt(op2._endian + b'2i 6i f 7i', self.size))
        ntotal = 64 * self.factor  # 16*4
        ndatai = len(data) - n
        nelements = ndatai // ntotal
        leftover = ndatai % ntotal
        assert leftover == 0, leftover
        for unused_i in range(nelements):
            edata = data[n:n + ntotal]
            out = struct_16i.unpack(edata)

            if op2.is_debug_file:
                op2.binary_debug.write('  CPLSTN6=%s\n' % str(out))
            (eid, pid, n1, n2, n3, n4, n5, n6, theta, *undef) = out
            nids = [n1, n2, n3, n4, n5, n6]
            assert min(nids) > 0, nids
            assert min(undef) == 0, undef
            assert max(undef) == 0, undef
            cplstn6 = op2.add_cplstn6(eid, pid, nids, theta=theta)
            #print(cplstn6)
            str(cplstn6)
            n += ntotal
        return n

    def read_cplstn8(self, data: bytes, n: int) -> int:
        """
        RECORD - CPLSTN8(7201,72,983)

        Word Name Type Description
        1 EID     I Element identification number
        2 PID     I Property identification number
        3 G(8)    I Grid point identification numbers of connection points
        11 THETA RS Material property orientation angle or coordinate system ID
        12 UNDEF(5) None

        """
        op2: OP2Geom = self.op2
        op2.to_nx(' because CPLSTN8 was found')
        struct_16i = Struct(mapfmt(op2._endian + b'2i 8i f 5i', self.size))
        ntotal = 64 * self.factor  # 16*4
        ndatai = len(data) - n
        nelements = ndatai // ntotal
        leftover = ndatai % ntotal
        assert leftover == 0, leftover
        for unused_i in range(nelements):
            edata = data[n:n + ntotal]
            out = struct_16i.unpack(edata)

            if op2.is_debug_file:
                op2.binary_debug.write('  CPLSTN8=%s\n' % str(out))
            (eid, pid, n1, n2, n3, n4, n5, n6, n7, n8, theta, *undef) = out
            nids = [n1, n2, n3, n4, n5, n6, n7, n8]
            assert min(nids) > 0, nids
            assert min(undef) == 0, undef
            assert max(undef) == 0, undef
            cplstn8 = op2.add_cplstn8(eid, pid, nids, theta=theta)
            str(cplstn8)
            n += ntotal
        return n

    def read_cpyram(self, data: bytes, n: int) -> int:
        """
        CPYRAM(17200,172,1000) - the marker for Record ???

        Specific to NX Nastran
        """
        op2: OP2Geom = self.op2
        ntotal = 64 * self.factor  # 16*4
        struct_16i = Struct(mapfmt(op2._endian + b'16i', self.size))
        nelements = (len(data) - n) // ntotal
        for unused_i in range(nelements):
            edata = data[n:n + ntotal]
            out = struct_16i.unpack(edata)
            if op2.is_debug_file:
                op2.binary_debug.write('  CPENTA=%s\n' % str(out))
            (eid, pid, g1, g2, g3, g4, g5, g6, g7, g8, g9, g10,
             g11, g12, g13, _g14) = out

            data_in = [eid, pid, g1, g2, g3, g4, g5]
            big_nodes = [g6, g7, g8, g9, g10, g11, g12, g13]
            if sum(big_nodes) > 0:
                elem = CPYRAM13.add_op2_data(data_in + big_nodes)
            else:
                elem = CPYRAM5.add_op2_data(data_in)
            self.add_op2_element(elem)
            n += ntotal
        op2.card_count['CPYRAM'] = nelements
        return n

# CPENP

    def read_cpenta(self, data: bytes, n: int) -> int:
        """
        CPENTA(4108,41,280)      - the marker for Record 63
        CPENPR(7509,75,9992)     - the marker for Record 64
        CPENT15F(16500,165,9999) - the marker for Record 65
        CPENT6FD(16000,160,9999) - the marker for Record 66
        """
        op2: OP2Geom = self.op2
        ntotal = 68 * self.factor  # 17*4
        s = Struct(mapfmt(op2._endian + b'17i', self.size))
        nelements = (len(data) - n) // ntotal
        for unused_i in range(nelements):
            edata = data[n:n + ntotal]
            out = s.unpack(edata)
            if op2.is_debug_file:
                op2.binary_debug.write('  CPENTA=%s\n' % str(out))
            (eid, pid, g1, g2, g3, g4, g5, g6, g7, g8, g9, g10,
             g11, g12, g13, g14, g15) = out

            data_in = [eid, pid, g1, g2, g3, g4, g5, g6]
            big_nodes = [g7, g8, g9, g10, g11, g12, g13, g14, g15]
            if sum(big_nodes) > 0:
                elem = CPENTA15.add_op2_data(data_in + big_nodes)
            else:
                elem = CPENTA6.add_op2_data(data_in)
            self.add_op2_element(elem)
            n += ntotal
        op2.card_count['CPENTA'] = nelements
        return n

    def read_cpenta_cz(self, data: bytes, n: int) -> int:
        """
        CPENTCZ
        """
        op2: OP2Geom = self.op2
        ntotal = 96 * self.factor  # 19*4
        s = Struct(mapfmt(op2._endian + b'17i 7i', self.size))
        nelements = (len(data) - n) // ntotal
        for unused_i in range(nelements):
            edata = data[n:n + ntotal]

            out = s.unpack(edata)
            if op2.is_debug_file:
                op2.binary_debug.write('  CPENTCZ=%s\n' % str(out))
            (eid, pid, g1, g2, g3, g4, g5, g6, g7, g8, g9, g10,
             g11, g12, g13, g14, g15,
             dummy1, dummy2, dummy3, dummy4, dummy5, dummy6, dummy7) = out
            dummy = (dummy1, dummy2, dummy3, dummy4, dummy5, dummy6, dummy7)
            assert dummy == (0, 0, 0, 0, 0, 0, 0), dummy
            #print(eid, pid, dummy)

            data_in = [eid, pid, g1, g2, g3, g4, g5, g6]
            big_nodes = [g7, g8, g9, g10, g11, g12, g13, g14, g15]
            #if sum(big_nodes) > 0:
            elem = CPENTCZ.add_op2_data(data_in + big_nodes)
            #else:
                #elem = CPENTA6.add_op2_data(data_in)
            self.add_op2_element(elem)
            n += ntotal
        op2.card_count['CPENTCZ'] = nelements
        return n

# CQDX4FD
# CQDX9FD - same as CQDX4FD

    def read_cquad(self, data: bytes, n: int) -> int:
        """CQUAD(9108,91,507)  - the marker for Record 69"""
        return self._run_cquad(CQUAD, data, n)

    def _read_cquad4_nasa95(self, data: bytes, n: int) -> int:
        """
        CQUAD4(5408,54,261) - the marker for Record 10

        CQUAD4       350       5       5     241     329      13
        CQUAD4       351       5      13     329     330      14
        CQUAD4       352       5      14     330     331      15
        ints    = (5408, 54, 261,
                   eid pid n1  n2    n3  n4
                   350, 5, 5, 241, 329,  13, 0, 0, 0, 0, 0, 0, 0,
                   351, 5, 13, 329, 330, 14, 0, 0, 0, 0, 0, 0, 0,
                   352, 5, 14, 330, 331, 15, 0, 0, 0, 0, 0, 0, 0)
        """
        op2: OP2Geom = self.op2
        elements = []
        ntotal = 52 * self.factor  # 13*4
        nelements = (len(data) - n) // ntotal
        leftover = (len(data) - n) % ntotal
        assert leftover == 0, leftover

        #  TODO: not sure...
        #   6i is correct
        #   3f-i zeros as float/int???
        #   4f correct
        #   i correct
        s = Struct(mapfmt(op2._endian + b'6i 7i', self.size))
        #if op2.is_debug_file:
            #op2.binary_debug.write('ndata=%s\n' % (nelements * 44))

        if op2.is_debug_file:
            op2.binary_debug.write(f'  {element.type}=(eid, pid, [n1, n2, n3, n4], theta, zoffs, '
                                    'unused_blank, [tflag, t1, t2, t3, t4]); theta_mcid\n')

        for unused_i in range(nelements):
            edata = data[n:n + ntotal]
            out = s.unpack(edata)
            # theta, zoffs, blank, tflag, t1, t2, t3, t4
            (eid, pid, n1, n2, n3, n4,
             theta, zoffs, tflag,
             t1, t2, t3, t4) = out

            msg = ('eid=%s pid=%s nodes=%s '
                   'theta=%s zoffs=%s '
                   'tflag=%s '
                   't1-t3=%s' % (
                       eid, pid, (n1, n2, n3, n4),
                       theta, zoffs,
                       tflag,
                       (t1, t2, t3, t4)))

            #assert theta == 0, msg
            assert zoffs == 0, msg
            assert tflag == 0, msg
            assert t1 == 0, msg
            assert t2 == 0, msg
            assert t3 == 0, msg
            assert t4 == 0, msg

            #theta_mcid = convert_theta_to_mcid(theta)
            #if op2.is_debug_file:
                #op2.binary_debug.write(
                    #f'  {element.type}=({eid}, {pid}, [{n1}, {n2}, {n3}, {n4}], '
                    #f'{theta}, {zoffs}, [{tflag}, {t1}, {t2}, {t3})]; {theta_mcid}\n')

            data_init = [
                eid, pid, n1, n2, n3, n4, theta, zoffs,
                tflag, t1, t2, t3, t4]
            elem = CQUAD4.add_op2_data(data_init)
            #elements.append(elem)
            self.add_op2_element(elem)
            n += ntotal
        op2.card_count['CQUAD4'] = nelements
        assert n == len(data), f'n={n} ndata={len(data)}'
        return n

    def read_cquad4(self, data: bytes, n: int) -> int:
        """
        CQUAD4(2958,51,177)    - the marker for Record 70
        CQUAD4(13900,139,9989) - the marker for Record 71
        """
        card_name = 'CQUAD4'
        card_obj = CQUAD4
        methods = {
            56 : self._run_cquad4_nx_56,
            60 : self._run_cquad4_msc_60,
        }
        try:
            n = self._read_double_card(card_name, card_obj, self.add_op2_element,
                                       methods, data, n)
        except DoubleCardError:
            nx_method = partial(self._run_cquad4_nx_56, card_obj)
            msc_method = partial(self._run_cquad4_msc_60, card_obj)
            n = self._read_dual_card(
                data, n,
                nx_method, msc_method,
                card_name, self.add_op2_element)
        #nentries = len(elements)
        #for elem in elements:
            #op2.add_op2_element(elem)
        #op2.card_count[card_name] = nentries

        return n

    def _read_double_card(self, card_name: str, card_obj, add_method,
                          methods, data: bytes, n: int) -> int:
        op2: OP2Geom = self.op2
        n, elements = self._read_double_card_load(
            card_name, card_obj,
            methods, data, n)

        nentries = len(elements)
        for elem in elements:
            add_method(elem)
        op2.card_count[card_name] = nentries
        return n

    def _read_double_card_load(self, card_name: str, card_obj,
                               methods, data: bytes, n: int) -> int:
        op2: OP2Geom = self.op2
        assert isinstance(data, bytes), type(data)
        ndatai = (len(data) - n) // self.factor
        keys = np.array(list(methods.keys()))

        errors = ndatai % keys
        izero = np.where(errors == 0)[0]
        if len(izero) == 0:
            op2.show_data(data, types='ifs')
            print(f'ndatai={ndatai} keys={keys} -> errors={errors}')
            print('izero = ', izero)
            raise EmptyCardError()

        elif len(izero) == 1:
            key0 = keys[izero[0]]
            method = methods[key0]
            #print(method)
            n, elements = method(card_obj, data, n)
        else:
            #print(keys, errors)
            methods2 = {key : methods[key] for i, key in enumerate(keys)
                        if i in izero}
            elements = None
            for key, method in methods2.items():
                try:
                    n, elements = method(card_obj, data, n)
                    #print('breaking')
                    #print(methods2)
                    break
                except Exception as e:
                    #print(str(e))
                    #print('error')
                    pass
                else:
                    raise DoubleCardError(f'No {card_name} cards found')
        #else:
        if elements is None:
            op2.show_data(data, types='ifs')
            print(f'ndatai={ndatai} keys={keys} -> errors={errors}')
            print('izero = ', izero)
            raise EmptyCardError()
        return n, elements

    def read_vutria3(self, data: bytes, n: int) -> int:
        return self._run_3nodes(CTRIA3, data, n)

    def read_vuquad4(self, data: bytes, n: int) -> int:
        """VUQUAD4(11201,112,9940)"""
        return self._run_4nodes(CQUAD4, data, n)

    def _run_cquad(self, element: Union[CQUAD, CQUADX], data: bytes, n: int) -> int:
        """common method for CQUAD, CQUADX"""
        op2: OP2Geom = self.op2
        ntotal = 44 * self.factor  # 11*4
        s = Struct(mapfmt(op2._endian + b'11i', self.size))
        nelements = (len(data) - n) // ntotal
        if op2.is_debug_file:
            op2.binary_debug.write('ndata=%s\n' % (nelements * ntotal))
        for unused_i in range(nelements):
            edata = data[n:n + ntotal]
            out = s.unpack(edata)
            (eid, pid, n1, n2, n3, n4, n5, n6, n7, n8, n9) = out
            if op2.is_debug_file:
                op2.binary_debug.write('  %s=%s\n' % (element.type, str(out)))
            #print('CQUAD eid=%s pid=%s n1=%s n2=%s n3=%s n4=%s n5=%s n6=%s n7=%s n8=%s' % (
                #eid, pid, n1, n2, n3, n4, n5, n6, n7, n8))
            datai = [eid, pid, n1, n2, n3, n4, n5, n6, n7, n8, n9]
            elem = element.add_op2_data(datai)
            self.add_op2_element(elem)
            n += ntotal
        op2.card_count[element.type] = nelements
        return n

    def _run_2nodes(self, element, data: bytes, n: int) -> int:
        """common method for VUBEAM"""
        op2: OP2Geom = self.op2
        nelements = (len(data) - n) // 16
        s = Struct(op2._endian + b'4i')
        #if op2.is_debug_file:
            #op2.binary_debug.write('ndata=%s\n' % (nelements * 44))

        if op2.is_debug_file:
            op2.binary_debug.write(f'  {element.type}=(eid, pid, [n1, n2]')

        for unused_i in range(nelements):
            edata = data[n:n + 16]  # 4*4
            out = s.unpack(edata)
            (eid, pid, n1, n2) = out
            if op2.is_debug_file:
                op2.binary_debug.write(
                    f'  {element.type}=({eid}, {pid}, [{n1}, {n2}]')

            nids = [n1, n2]
            elem = element(eid, pid, nids)
            self.add_op2_element(elem)
            n += 16
        #if stop:
            #raise RuntimeError('theta is too large...make the quad wrong')
        op2.card_count[elem.type] = nelements
        return n

    def _run_3nodes(self, element: CTRIA3, data: bytes, n: int) -> int:
        """common method for CTRIA3, VUTRIA3"""
        op2: OP2Geom = self.op2
        nelements = (len(data) - n) // 20
        s = Struct(op2._endian + b'5i')
        #if op2.is_debug_file:
            #op2.binary_debug.write('ndata=%s\n' % (nelements * 44))

        if op2.is_debug_file:
            op2.binary_debug.write(f'  {element.type}=(eid, pid, [n1, n2, n3]')

        for unused_i in range(nelements):
            edata = data[n:n + 20]  # 5*4
            out = s.unpack(edata)
            (eid, pid, n1, n2, n3) = out
            if op2.is_debug_file:
                op2.binary_debug.write(
                    f'  {element.type}=({eid}, {pid}, [{n1}, {n2}, {n3}]')

            nids = [n1, n2, n3]
            elem = element(eid, pid, nids)
            self.add_op2_element(elem)
            n += 20
        #if stop:
            #raise RuntimeError('theta is too large...make the quad wrong')
        op2.card_count[element.type] = nelements
        return n

    def _run_4nodes(self, element: Union[CQUAD4, CAABSF], data: bytes, n: int) -> int:
        """common method for CQUAD4, CQUADR"""
        op2: OP2Geom = self.op2
        nelements = (len(data) - n) // 24
        s = Struct(op2._endian + b'6i')
        #if op2.is_debug_file:
            #op2.binary_debug.write('ndata=%s\n' % (nelements * 44))

        if op2.is_debug_file:
            op2.binary_debug.write(f'  {element.type}=(eid, pid, [n1, n2, n3, n4]')

        for unused_i in range(nelements):
            edata = data[n:n + 24]  # 6*4
            out = s.unpack(edata)
            (eid, pid, n1, n2, n3, n4) = out
            if op2.is_debug_file:
                op2.binary_debug.write(
                    f'  {element.type}=({eid}, {pid}, [{n1}, {n2}, {n3}, {n4}]')

            nids = [n1, n2, n3, n4]
            elem = element(eid, pid, nids)
            self.add_op2_element(elem)
            n += 24
        #if stop:
            #raise RuntimeError('theta is too large...make the quad wrong')
        op2.card_count[element.type] = nelements
        return n

    def _run_cquad4_msc_60(self, element: CQUAD4, data: bytes, n: int) -> tuple[int, Any]:
        r"""
        buggy MSC 2018.2 version

        #TODO is this right?  What's with the intermediate zeros?
        data = (
            2958, 51, 177,
            1, 1, 1, 2, 73, 72, 0, 0, 0, 0, -1.0, -1.0, -1.0, -1.0, -1,
            2, 1, 2, 3, 74, 73, 0, 0, 0, 0, -1.0, -1.0, -1.0, -1.0, -1,
            3, 1, 3, 4, 75, 74, 0, 0, 0, 0, -1.0, -1.0, -1.0, -1.0, -1,
        )

        C:\MSC.Software\msc_nastran_runs\cc179h.op2

        data = (
            eid   pid n1     n2     n3     n4     theta zoffs  blank, tflag   t1    t2    t3    t4   -1/mcid?
            101,  30, 30000, 30001, 30101, 30100, 19.2, 0.0,   0,     0,      0.2,  0.2,  0.2,  0.2, -1,
            102,  30, 30001, 30002, 30102, 30101, 19.2, 0.0,   0,     0,      0.2,  0.2,  0.2,  0.2, -1,
            103,  30, 30002, 30003, 30103, 30102, 19.2, 0.0,   0,     0,      0.2,  0.2,  0.2,  0.2, -1,
            104,  30, 30100, 30101, 30201, 30200, 19.2, 0.0,   0,     0,      0.2,  0.2,  0.2,  0.2, -1,
            1001, 20, 20000, 20001, 20101, 20100, 0,    0.0,   0,     0,     -1.0, -1.0, -1.0, -1.0, -1,
            1002, 20, 20001, 20002, 20102, 20101, 0,    0.0,   0,     0,     -1.0, -1.0, -1.0, -1.0, -1,
            1003, 20, 20002, 20003, 20103, 20102, 0,    0.0,   0,     0,     -1.0, -1.0, -1.0, -1.0, 42,
            1004, 20, 20100, 20101, 20201, 20200, 0,    0.0,   0,     0,     -1.0, -1.0, -1.0, -1.0, -1)
        """
        op2: OP2Geom = self.op2
        elements = []
        ntotal = 60 * self.factor  # 15*4
        nelements = (len(data) - n) // ntotal
        leftover = (len(data) - n) % ntotal
        assert leftover == 0, leftover


        #  TODO: not sure...
        #   6i is correct
        #   3f-i zeros as float/int???
        #   4f correct
        #   i correct
        structi = Struct(mapfmt(op2._endian + b'6i 3fi 4f i', self.size))
        #if op2.is_debug_file:
            #op2.binary_debug.write('ndata=%s\n' % (nelements * 44))

        if op2.is_debug_file:
            op2.binary_debug.write(f'  {element.type}=(eid, pid, [n1, n2, n3, n4], theta, zoffs, '
                                    'unused_blank, [tflag, t1, t2, t3, t4]); theta_mcid\n')

        cids = set(list(self.op2.coords.keys()))
        check_cids = max(cids) > 0 if len(cids) else False
        for unused_i in range(nelements):
            edata = data[n:n + ntotal]
            #self.show_data(edata)
            out = structi.unpack(edata)
            # theta, zoffs, blank, tflag, t1, t2, t3, t4
            (eid, pid, n1, n2, n3, n4,
             theta, zoffs, blank, tflag,
             t1, t2, t3, t4,
             minus1_mcid) = out
            assert blank == 0.0
            assert isinstance(tflag, int), tflag
            #minus1 = out[-1]
            #eid, pid, n1, n2, n3, n4, theta, a, b, c, t1, t2, t3, t4, minus1 = out
            #print(eid, pid)
            theta_mcid = theta
            if minus1_mcid != -1:
                theta_mcid = minus1_mcid
                assert theta == 0., theta
                if check_cids:
                    assert minus1_mcid in cids, minus1_mcid
            msg = (
                f'eid={eid} pid={pid} nodes=({n1}, {n2}, {n3}, {n4}) '
                f'theta={theta} zoffs={zoffs} blank={blank} '
                f'tflag={tflag} t1-t4=({t1}, {t2}, {t3}, {t4}) '
                f'minus1_mcid={minus1_mcid} -> theta_mcid={theta_mcid}')
            #print(msg)
            #assert theta == 0, msg
            assert zoffs == 0, msg
            assert blank == 0, msg
            assert tflag == 0, msg

            if op2.is_debug_file:
                op2.binary_debug.write(
                    f'  {element.type}=({eid}, {pid}, [{n1}, {n2}, {n3}, {n4}], '
                    f'{theta}, {zoffs}, {blank}, [{tflag}, {t1}, {t2}, {t3}, {t4})]; {theta_mcid}\n')

            data_init = [
                eid, pid, n1, n2, n3, n4, theta_mcid, zoffs,
                tflag, t1, t2, t3, t4]
            elem = element.add_op2_data(data_init)
            elements.append(elem)
            n += ntotal
        #if stop:
            #raise RuntimeError('theta is too large...make the quad wrong')
        #op2.card_count[element.type] = nelements
        return n, elements

    def _run_cquad4_nx_56(self, element: Union[CQUAD4, CQUADR],
                          data: bytes, n: int) -> tuple[int, Any]:
        """
        common method for CQUAD4, CQUADR

        data = (
            2958, 51, 177,
            1, 1, 1, 2, 8, 7, 0, 0, 0, 0, -1.0, -1.0, -1.0, -1.0
            2, 1, 2, 3, 9, 8, 0, 0, 0, 0, -1.0, -1.0, -1.0, -1.0
            3, 1, 3, 4, 10, 9, 0, 0, 0, 0, -1.0, -1.0, -1.0, -1.0
        )
        """
        op2: OP2Geom = self.op2
        elements = []
        ntotal = 56 * self.factor  # 14*4
        nelements = (len(data) - n) // ntotal
        leftover = (len(data) - n) % ntotal
        assert leftover == 0, leftover
        s = Struct(mapfmt(op2._endian + b'6i ff ii 4f', self.size))
        if op2.is_debug_file:
            op2.binary_debug.write('ndata=%s\n' % (nelements * 44))

        if op2.is_debug_file:
            op2.binary_debug.write(f'  {element.type}=(eid, pid, [n1, n2, n3, n4], theta, zoffs, '
                                    'unused_blank, [tflag, t1, t2, t3, t4]); theta_mcid\n')

        for unused_i in range(nelements):
            edata = data[n:n + ntotal]
            out = s.unpack(edata)
            (eid, pid, n1, n2, n3, n4, theta, zoffs, unused_blank, tflag,
             t1, t2, t3, t4) = out
            theta_mcid = convert_theta_to_mcid(theta)
            if op2.is_debug_file:
                op2.binary_debug.write(
                    f'  {element.type}=({eid}, {pid}, [{n1}, {n2}, {n3}, {n4}], '
                    f'{theta}, {zoffs}, {unused_blank}, [{tflag}, {t1}, {t2}, {t3}, {t4})]; {theta_mcid}\n')

            #print('eid=%s pid=%s n1=%s n2=%s n3=%s n4=%s theta=%s zoffs=%s '
                  #'blank=%s tflag=%s t1=%s t2=%s t3=%s t4=%s' % (
                      #eid, pid, n1, n2, n3, n4, theta, zoffs,
                      #blank, tflag, t1, t2, t3, t4))

            data_init = [
                eid, pid, n1, n2, n3, n4, theta_mcid, zoffs,
                tflag, t1, t2, t3, t4]
            elem = element.add_op2_data(data_init)
            #op2.add_op2_element(elem)
            elements.append(elem)
            n += ntotal
        #if stop:
            #raise RuntimeError('theta is too large...make the quad wrong')
        #op2.card_count[element.type] = nelements
        #self.to_nx()  # not really an nx specific thing...
        return n, elements

# CQUAD4FD

    def read_cquad8(self, data: bytes, n: int) -> int:
        """common method for reading CQUAD8s"""
        card_name = 'CQUAD8'
        card_obj = CQUAD8
        methods = {
            68 : self._read_cquad8_current,
            64 : self._read_cquad8_v2001,
            72 : self._read_cquad8_72,
        }
        try:
            n = self._read_double_card(card_name, card_obj, self.add_op2_element,
                                       methods, data, n)
        except DoubleCardError:
            raise
            #self.op2.log.warning(f'try-except {card_name}')
            #n = self._read_split_card(data, n,
                                      #self._read_cquad8_current, self._read_cquad8_v2001,
                                      #card_name, op2.add_op2_element)
        #nelements = op2.card_count['CQUAD8']
        #op2.log.debug(f'nCQUAD8 = {nelements}')
        return n

    def _read_cquad8_current(self, card_obj, data: bytes, n: int) -> int:
        """
        CQUAD8(4701,47,326)  - the marker for Record 72
        .. warning:: inconsistent with dmap manual

        1 EID     I Element identification number
        2 PID     I Property identification number
        3 G(8)    I Grid point identification numbers of
                    connection points
        11 T(4)  RS Membrane thickness of element at grid
                    points
        15 THETA RS Material property orientation angle or
                    coordinate system identification number
        16 ZOFFS RS Offset from the surface of grid points
                    reference plane
        17 TFLAG  I Relative thickness flag
        """
        op2: OP2Geom = self.op2
        ntotal = 68 * self.factor # 17*4
        ndatai = len(data) - n
        nelements = ndatai // ntotal
        assert ndatai % ntotal == 0
        s = Struct(mapfmt(op2._endian + b'10i 6f i', self.size))
        elements = []
        for unused_i in range(nelements):
            edata = data[n:n + ntotal]
            out = s.unpack(edata)
            if op2.is_debug_file:
                op2.binary_debug.write('  CQUAD8=%s\n' % str(out))
            #(eid, pid, n1, n2, n3, n4, n5, n6, n7, n8, t1, t2,
             #t3, t4, theta, zoffs, tflag) = out
            tflag = out[-1]
            #op2.log.info('cquad8 tflag = %s' % tflag)
            assert isinstance(tflag, int), tflag
            assert tflag in [-1, 0, 1], tflag
            #print('eid=%s pid=%s n1=%s n2=%s n3=%s n4=%s theta=%s zoffs=%s '
                  #'tflag=%s t1=%s t2=%s t3=%s t4=%s' % (
                      #eid, pid, n1, n2, n3, n4, theta, zoffs,
                      #tflag, t1, t2, t3, t4))
            #data_init = [eid,pid,n1,n2,n3,n4,theta,zoffs,tflag,t1,t2,t3,t4]
            elem = CQUAD8.add_op2_data(out)
            elements.append(elem)
            n += ntotal
        return n, elements

    def _read_cquad8_v2001(self, card_obj, data: bytes, n: int) -> int:
        """
        CQUAD8(4701,47,326)  - the marker for Record 72

        1 EID     I Element identification number
        2 PID     I Property identification number
        3 G(8)    I Grid point identification numbers of
                    connection points
        11 T(4)  RS Membrane thickness of element at grid
                    points
        15 THETA RS Material property orientation angle or
                    coordinate system identification number
        16 ZOFFS RS Offset from the surface of grid points
                    reference plane
        """
        op2: OP2Geom = self.op2
        elements = []
        #self.show_data(data, types='if')

        ntotal = 64 * self.factor  # 16*4
        ndatai = (len(data) - n)
        nelements = ndatai // ntotal
        assert ndatai % ntotal == 0
        sf = Struct(mapfmt(op2._endian + b'10i 6f', self.size))
        for unused_i in range(nelements):
            edata = data[n:n + ntotal]
            out = sf.unpack(edata)
            if op2.is_debug_file:
                op2.binary_debug.write('  CQUAD8=%s\n' % str(out))
            (eid, pid, n1, n2, n3, n4, n5, n6, n7, n8, t1, t2,
             t3, t4, theta, zoffs) = out
            assert eid > 0
            assert pid > 0
            tflag = None
            out = (eid, pid, n1, n2, n3, n4, n5, n6, n7, n8, t1, t2,
                   t3, t4, theta, zoffs, tflag)
            #print('eid=%s pid=%s n1=%s n2=%s n3=%s n4=%s theta=%g zoffs=%s '
                  #'tflag=%s t1=%g t2=%g t3=%g t4=%g' % (
                      #eid, pid, n1, n2, n3, n4, theta, zoffs, tflag, t1, t2, t3, t4))
            #data_init = [eid,pid,n1,n2,n3,n4,theta,zoffs,tflag,t1,t2,t3,t4]
            elem = CQUAD8.add_op2_data(out)
            elements.append(elem)
            self.add_op2_element(elem)
            n += ntotal
        return n, elements


    def _read_cquad8_72(self, card_obj, data: bytes, n: int) -> int:
        """
        data = (
            eid   pid n1     n2     n3     n4     n5     n6     n7     n8     t1   t2   t3   t4   theta ? ?
            301,  30, 40000, 40002, 40202, 40200, 40001, 40102, 40201, 40100, 0.2, 0.2, 0.2, 0.2, 19.2, 0, 0, -1,
            302,  30, 40002, 40004, 40204, 40202, 40003, 40104, 40203, 40102, 0.2, 0.2, 0.2, 0.2, 19.2, 0, 0, -1,
            303,  30, 40200, 40202, 40402, 40400, 40201, 40302, 40401, 40300, 0.2, 0.2, 0.2, 0.2, 19.2, 0, 0, -1,
            304,  30, 40202, 40204, 40404, 40402, 40203, 40304, 40403, 40302, 0.2, 0.2, 0.2, 0.2, 19.2, 0, 0, -1,
            1301, 20, 10000, 10002, 10202, 10200, 10001, 10102, 10201, 10100, -1.0, -1.0, -1.0, -1.0, 0, 0, 0, -1,
            1302, 20, 10002, 10004, 10204, 10202, 10003, 10104, 10203, 10102, -1.0, -1.0, -1.0, -1.0, 0, 0, 0, -1,
            1303, 20, 10200, 10202, 10402, 10400, 10201, 10302, 10401, 10300, -1.0, -1.0, -1.0, -1.0, 0, 0, 0, -1,
            1304, 20, 10202, 10204, 10404, 10402, 10203, 10304, 10403, 10302, -1.0, -1.0, -1.0, -1.0, 0, 0, 0, -1)
        """
        op2: OP2Geom = self.op2
        elements = []
        #self.show_data(data, types='if')
        #ss
        ntotal = 72 * self.factor  # 16*4
        ndatai = (len(data) - n)
        nelements = ndatai // ntotal
        assert ndatai % ntotal == 0
        s = Struct(mapfmt(op2._endian + b'10i 5f 3i', self.size))

        #sf = Struct(mapfmt(op2._endian + b'10i 6f', self.size))
        edata0 = data[n:n + ntotal]
        flag = s.unpack(edata0)[-1]
        if flag == -1:
            for unused_i in range(nelements):
                edata = data[n:n + ntotal]
                out = s.unpack(edata)
                if op2.is_debug_file:
                    op2.binary_debug.write('  CQUAD8=%s\n' % str(out))
                (eid, pid, n1, n2, n3, n4, n5, n6, n7, n8, t1, t2,
                 t3, t4, theta, tflag, zoffs, flag) = out
                assert eid > 0
                assert pid > 0
                assert tflag == 0
                assert zoffs == 0
                assert flag == -1, flag
                tflag = None
                out = (eid, pid, n1, n2, n3, n4, n5, n6, n7, n8, t1, t2,
                       t3, t4, theta, zoffs, tflag)
                #print('eid=%s pid=%s n1=%s n2=%s n3=%s n4=%s theta=%g zoffs=%s '
                      #'tflag=%s t1=%g t2=%g t3=%g t4=%g' % (
                          #eid, pid, n1, n2, n3, n4, theta, zoffs, tflag, t1, t2, t3, t4))
                #data_init = [eid,pid,n1,n2,n3,n4,theta,zoffs,tflag,t1,t2,t3,t4]
                elem = CQUAD8.add_op2_data(out)
                elements.append(elem)
                self.add_op2_element(elem)
                n += ntotal
        else:
            for unused_i in range(nelements):
                edata = data[n:n + ntotal]
                out = s.unpack(edata)
                if op2.is_debug_file:
                    op2.binary_debug.write('  CQUAD8=%s\n' % str(out))
                (eid, pid, n1, n2, n3, n4, n5, n6, n7, n8, t1, t2,
                 t3, t4, theta, zoffs) = out
                assert eid > 0
                assert pid > 0
                tflag = None
                out = (eid, pid, n1, n2, n3, n4, n5, n6, n7, n8, t1, t2,
                       t3, t4, theta, zoffs, tflag)
                #print('eid=%s pid=%s n1=%s n2=%s n3=%s n4=%s theta=%g zoffs=%s '
                      #'tflag=%s t1=%g t2=%g t3=%g t4=%g' % (
                          #eid, pid, n1, n2, n3, n4, theta, zoffs, tflag, t1, t2, t3, t4))
                #data_init = [eid,pid,n1,n2,n3,n4,theta,zoffs,tflag,t1,t2,t3,t4]
                elem = CQUAD8.add_op2_data(out)
                elements.append(elem)
                self.add_op2_element(elem)
                n += ntotal
        return n, elements

# CQUAD9FD
# CQUADP

    def read_cquadr(self, data: bytes, n: int) -> int:
        """CQUADR(8009,80,367)  - the marker for Record 75

        ints    = (8009, 80, 367,
            pshell go from 1-242
            nids   go from 1-29382
            1, 1,   3600, 3310, 3797, 731, 0, 0, 0, 0, -1.0, -1.0, -1.0, -1.0, 1,
            2, 1, 731, 3797, 3798, 732, 0, 0, 0, 0, -1.0, -1.0, -1.0, -1.0, 1,
            3, 1, 732, 3798, 3799, 733, 0, 0, 0, 0, -1.0, -1.0, -1.0, -1.0, 1,
            ...
        )
        """
        #return len(data)
        #self.op2.show_data(data, types='if')
        card_name = 'CQUADR'
        card_obj = CQUADR
        methods = {
            56 : self._run_cquad4_nx_56,
            60 : self._run_cquad4_msc_60,  # 15*4
        }
        try:
            n = self._read_double_card(card_name, card_obj, self.add_op2_element,
                                       methods, data, n)
        except DoubleCardError:
            nx_method = partial(self._run_cquad4_nx_56, card_obj)
            msc_method = partial(self._run_cquad4_msc_60, card_obj)
            n = self._read_dual_card(
                data, n,
                nx_method, msc_method,
                card_name, self.add_op2_element)
        return n

    def read_cquadx(self, data: bytes, n: int) -> int:
        """CQUADX(9008,90,508)  - the marker for Record 76"""
        return self._run_cquad(CQUADX, data, n)

    def read_crbar(self, data: bytes, n: int) -> int:
        """
        Word Name Type Description
        1 EID     I Element identification number
        2 GA      I Grid point A identification number
        3 GB      I Grid point B identification number
        4 LMID1   I Lagrange multiplier identification number
        5 NDOFS   I Number of DOF for Lagrange multiplier
        6 CNA     I Component numbers of independent DOFs at end A
        7 CNB     I Component numbers of independent DOFs at end B
        8 CMA     I Component numbers of dependent DOFs at end A
        9 CMB     I Component numbers of dependent DOFs at end B
        10 ALPHA RS Thermal expansion coefficient
        ints    = (3, 3, 41, 101000041, 6, 123456, 0, 0,     123456, 2.e-6,
                   4, 3, 4, 101000004, 6, 123456, 0, 0, 654321, 0)
        floats  = (3, 3, 41, 101000041, 6, 123456, 0.0, 0.0, 123456, 2.e-6,
                   4, 3, 4,  101000004, 6, 123456, 0.0, 0.0, 654321, 0.0)

        """
        # C:\NASA\m4\formats\git\examples\move_tpl\nrgd20c.op2
        op2: OP2Geom = self.op2
        idata = np.frombuffer(data[n:], op2.idtype8).copy()
        fdata = np.frombuffer(data[n:], op2.fdtype8).copy()
        #op2.show_data(data[n:], types='ifs', endian=None, force=False)
        nrows = len(idata) // 10
        ncols = 10
        idata = idata.reshape(nrows, ncols)
        fdata = fdata.reshape(nrows, ncols)
        eids = idata[:, 0]
        nids = idata[:, [1, 2]]
        #lagrange_id = idata[:, 3]
        #ndof_lagrange = idata[:, 4]
        cnas = idata[:, 5].astype('U6').tolist()
        cnbs = idata[:, 6].astype('U6').tolist()
        cmas = idata[:, 7].astype('U6').tolist()
        cmbs = idata[:, 8].astype('U6').tolist()
        alphas = fdata[:, 9]
        for eid, nidsi, cna, cnb, cma, cmb, alpha in zip(eids, nids, cnas, cnbs, cmas, cmbs, alphas):
            cna = ''.join(sorted(cna))
            cnb = ''.join(sorted(cnb))
            cma = ''.join(sorted(cma))
            cmb = ''.join(sorted(cmb))
            if cna == '0':
                cna = ''
            if cnb == '0':
                cnb = ''
            if cma == '0':
                cma = ''
            if cna == '0':
                cmb = ''
            elem = op2.add_rbar(eid, nidsi, cna, cnb, cma, cmb, alpha=alpha, comment='')
            str(elem)
        return len(data)

    def read_crbe1(self, data: bytes, n: int) -> int:
        """
        Word Name Type Description
        1 EID    I Element identification number
        2 NWE    I Number of words for the element
        3 ELTYPE I Element type: 1-RBE1 2-RBE2 3-RTRPLT 4-RTRPLT1
        4 GN     I Grid point identification number for independent degrees-of-freedom
        5 CN     I Component numbers of independent degrees-of-freedom
        Words 4 through 5 repeat until (-2,-2) occurs
        6 GM     I Grid point identification number for dependent degrees-of-freedom
        7 CM     I Component numbers of dependent degrees of-freedom
        Words 6 through 7 repeat until (-3,-3) occurs
        8 ALPHA RS  Thermal expansion coefficient
        9 UNDEF none Not used
        Words 8 through 9 repeat until (-4,-4) occurs
        10 LMID I Lagrange multiplier identification number
        11 NDOF I Number of DOF for each Lagrange multiplier identification number
        Words 10 through 11 repeat until (-1,-1) occurs

        ints   = [35, 22, 2,
                     3,  123456, -2, -2,
                     41, 123456, 51, 123456, -3, -3,
                     990057071, 0, -4, -4,
                  101000041, 6, 101000051, 6, -1, -1,
                  46, 22, 2, 3, 123456, -2, -2,
                  4, 123456, 5, 123456, -3, -3,
                  906377149, 0, -4, -4,
                  101000004, 6, 101000005, 6, -1, -1]
        floats = [4.90454462513686e-44, 3.0828566215145976e-44, 2.802596928649634e-45, 4.203895392974451e-45, 1.7299870321168462e-40, nan, nan, 5.74532370373175e-44, 1.7299870321168462e-40, 7.146622168056567e-44, 1.7299870321168462e-40, nan, nan, 0.0020000000949949026, 0.0, nan, nan, 2.5040535059644284e-35, 8.407790785948902e-45, 2.5040563758236833e-35, 8.407790785948902e-45, nan, nan, 6.445972935894159e-44, 3.0828566215145976e-44, 2.802596928649634e-45, 4.203895392974451e-45, 1.7299870321168462e-40, nan, nan, 5.605193857299268e-45, 1.7299870321168462e-40, 7.006492321624085e-45, 1.7299870321168462e-40, nan, nan, 1.9999999949504854e-06, 0.0, nan, nan, 2.504042887485185e-35, 8.407790785948902e-45, 2.5040431744711106e-35, 8.407790785948902e-45, nan, nan]

        """
        # C:\NASA\m4\formats\git\examples\move_tpl\nrgd406a.op2
        op2: OP2Geom = self.op2
        idata = np.frombuffer(data[n:], op2.idtype8).copy()
        fdata = np.frombuffer(data[n:], op2.fdtype8).copy()
        i = 0
        nints = len(idata)
        while i < nints:
            #print('i = ', i)
            eid = idata[i]
            nwe = idata[i+1]
            eltype = idata[i+2]
            assert eltype in {1, 2, 3, 4}, eltype
            i += 3

            Gni = []
            Cni = []
            while i < nints:
                gn, cn = idata[i:i+2]
                i += 2
                if (gn, cn) == (-2, -2):
                    break
                Gni.append(gn)
                Cni.append(cn)
            assert len(Gni), Gni

            Gmi = []
            Cmi = []
            while i < nints:
                gm, cm = idata[i:i+2]
                i += 2
                if (gm, cm) == (-3, -3):
                    break
                Gmi.append(gm)
                Cmi.append(cm)
            assert len(Gmi), Gmi

            alphas = []
            while i < nints:
                i1, i2 = idata[i:i+2]
                alpha = fdata[i]
                # undef
                i += 2
                if (i1, i2) == (-4, -4):
                    break
                alphas.append(alpha)
            assert len(alphas) == 1, alphas


            #10 LMID I Lagrange multiplier identification number
            #11 NDOF I Number of DOF for each Lagrange multiplier identification number
            lmid = []
            ndof = []
            while i < nints:
                lmidi, ndofi = idata[i:i+2]
                i += 2
                if (lmidi, ndofi) == (-1, -1):
                    break
                lmid.append(lmidi)
                ndof.append(ndofi)
            assert len(ndof), ndof

            #cna = ''.join(sorted(cna))
            #cnb = ''.join(sorted(cnb))
            #cma = ''.join(sorted(cma))
            #cmb = ''.join(sorted(cmb))
            #if cna == '0':
                #cna = ''
            #if cnb == '0':
                #cnb = ''
            #if cma == '0':
                #cma = ''
            #if cna == '0':
                #cmb = ''
            elem = op2.add_rbe1(eid, Gni, Cni, Gmi, Cmi, alpha=alpha, comment='')
            str(elem)
        return len(data)

    def read_crbe3(self, data: bytes, n: int) -> int:
        """
        This card is an internal RBE3 that's used for Langrage elements.
        It's not what the user entered.

        Word Name Type Description
        1 EID  I  Element identification number
        2 NWE  I  Number of words for the element
        3 REFG I  Reference grid point identification number
        4 REFC I  Component numbers at the reference grid point
        5 WT1  RS Weighting factor for components of motion at G
        6 C    I  Component numbers
        7 G    I  Grid point identification number
        Word 7 repeats until End of Record (-1)
        Words 5 through 7 repeat until End of Record (-2)

        8 GM I Grid point identification number for dependent DOFs
        9 CM I Component numbers of dependent DOFs
        Words 8 through 9 repeat until End of Record (-4?)

        10 ALPHA RS Thermal expansion coefficient
        Word 10 repeats until End of Record (-5?)

        11 LMID1 I Lagrange multiplier identification number
        12 NDOFS I Number of DOF for Lagrange multiplier
        Words 11 through 12 repeat until End of Record (-3?)

        data = (3, 14, 41, 123456, 1.0, 123456, 3, -1,
                -2,
                -4,
                0.002, -5,
                101000041, 6, -3,
                4, 14, 4,  123456, 1.0, 123456, 3, -1, -2, -4, 2.0e-6, -5, 101000004, 6, -3)
          """
        op2: OP2Geom = self.op2
        # C:\NASA\m4\formats\git\examples\move_tpl\ngd720a.op2
        idata = np.frombuffer(data[n:], op2.idtype8).copy()
        fdata = np.frombuffer(data[n:], op2.fdtype8).copy()

        iminus3 = np.where(idata == -3)[0]
        if idata[-1] == -3:
            is_alpha = False
            i = np.hstack([[0], iminus3[:-1]+1])
        else:
            is_alpha = True
            i = np.hstack([[0], iminus3[:-1]+2])
        j = np.hstack([iminus3[:-1], len(idata)])
        #print('i = ', i)
        #print('j = ', j)
        #print('is_alpha =', is_alpha)

        assert len(i) == len(j)
        for ii, jj in zip(i, j):

            idatai = idata[ii:jj]
            #print(idatai)
            #print(fdata[ii:ii+5])
            #1 EID  I  Element identification number
            #2 NWE  I  Number of words for the element
            #3 REFG I  Reference grid point identification number
            #4 REFC I  Component numbers at the reference grid point
            #5 WT1  RS Weighting factor for components of motion at G
            #6 C    I  Component numbers
            #7 G    I  Grid point identification number
            eid, nwords, refg, refc, unused_weight, comp, grid = idatai[:7]
            #nwords += 1
            #print(f'nwords={nwords}; len(datai)={len(idatai)}', len(idatai))
            weight = fdata[ii+4]
            weights = [weight]
            comps = [comp]
            gijs = [grid]

            iji = ii
            values = []
            values_dict = {}
            value = idata[iji]
            while value != -3:
                if value < 0:
                    if len(values):
                        values_dict[value] = values
                    values = []
                else:
                    values.append(iji)
                iji += 1
                value = idata[iji]
            if len(values):
                values_dict[value] = values

            #print('values_dict =', values_dict)
            del values_dict[-1]
            #{-1: [3, 14, 41, 123456, 1065353216, 123456, 3], -5: [990057071], -3: [101000041, 6]}
            if -5 in values_dict:
                ialpha_list = values_dict[-5]
                assert len(ialpha_list) == 1, ialpha_list
                alpha = fdata[ialpha_list[0]]
                #print("alpha =", alpha)
                del values_dict[-5]
            if -3 in values_dict:
                ilagrange_ndof = values_dict[-3]
                assert len(ilagrange_ndof) >= 2, ilagrange_ndof
                assert len(ilagrange_ndof) % 2 == 0, ilagrange_ndof
                lagrange_ndof = idata[ilagrange_ndof]
                #print("lagrange_ndof =", lagrange_ndof)
                del values_dict[-3]
            if values_dict:
                print('values_dict =', values_dict)
                assert len(values_dict) == 0, values_dict

            nrows = len(lagrange_ndof) // 2
            lagrange_ndof = lagrange_ndof.reshape(nrows, 2)
            gmi = lagrange_ndof[:, 0]
            cmi = lagrange_ndof[:, 1]
            in_data = [eid, refg, refc, weights, comps, gijs,
                       gmi, cmi, alpha]
            if op2.is_debug_file:
                op2.binary_debug.write('  RBE3=%s\n' % str(in_data))
            #print('rbe3 =', in_data)
            rbe3 = RBE3.add_op2_data(in_data)
            str(rbe3)
            continue
            #self._add_op2_rigid_element(rbe3)
        return len(data)

    def read_crjoint(self, data: bytes, n: int) -> int:
        """
        Word Name Type Description
        1 EID   I Element identification number
        2 GA    I Grid point A identification number
        3 GB    I Grid point B identification number
        4 LMID1 I Lagrange multiplier identification number
        5 NDOFS I Number of DOF for Lagrange multiplier
        6 CB    I Component numbers of dependent DOFs at end B
        """
        # no rjoint...
        # C:\NASA\m4\formats\git\examples\move_tpl\ngd720a.op2
        self.op2.log.info('geom skipping RJOINT in GEOM2')
        return len(data)

    def read_crod(self, data: bytes, n: int) -> int:
        """
        CROD(3001,30,48)    - the marker for Record 81
        """
        op2: OP2Geom = self.op2
        ntotal = 16 * self.factor # 4*4
        struct_4i = Struct(mapfmt(op2._endian + b'4i', self.size))
        nelements = (len(data) - n) // ntotal
        #is_long_ids = False
        for unused_i in range(nelements):
            edata = data[n:n + ntotal]
            out = struct_4i.unpack(edata)
            if op2.is_debug_file:
                op2.binary_debug.write('  CROD=%s\n' % str(out))
            #(eid, pid, n1, n2) = out
            #if n1 > 100000000 or n2 > 100000000:
                #is_long_ids = True
            elem = CROD.add_op2_data(out)
            self.add_op2_element(elem)
            n += ntotal
        #self._is_long_ids = is_long_ids
        op2.card_count['CROD'] = nelements
        return n

    def read_crrod(self, data: bytes, n: int) -> int:
        """
        Word Name Type Description
        1 EID    I Element identification number
        2 GA     I Grid point A identification number
        3 GB     I Grid point B identification number
        4 LMID1  I Lagrange multiplier identification number
        5 CMA    I Component numbers of dependent DOFs at end A
        6 CMB    I Component numbers of dependent DOFs at end B
        7 ALPHA RS Thermal expansion cofficient
        """
        op2: OP2Geom = self.op2
        structi = Struct(op2._endian + b'6if')
        nelements = (len(data) - n) // 28  # 7*4
        #is_long_ids = False
        for unused_i in range(nelements):
            edata = data[n:n + 28]
            out = structi.unpack(edata)
            if op2.is_debug_file:
                op2.binary_debug.write('  RROD=%s\n' % str(out))

            #print(out)
            (eid, n1, n2, unused_lagrange_multiplier_id, cma, cmb, alpha) = out
            assert cma == 0, (eid, n1, n2, unused_lagrange_multiplier_id, cma, cmb, alpha)
            assert cmb == 0, (eid, n1, n2, unused_lagrange_multiplier_id, cma, cmb, alpha)
            #if n1 > 100000000 or n2 > 100000000:
                #is_long_ids = True

            #eid, nids, cma='', cmb='', alpha=0.0
            elem = op2.add_rrod(eid, [n1, n2], cma=cma, cmb=cmb, alpha=alpha)
            self.add_op2_element(elem)
            n += 28
        #self._is_long_ids = is_long_ids
        op2.card_count['RROD'] = nelements
        return n

# CSEAM

    def read_cshear(self, data: bytes, n: int) -> int:
        """
        CSHEAR(3101,31,61)    - the marker for Record 84
        """
        op2: OP2Geom = self.op2
        ntotal = 24 * self.factor  # 6*4
        struct_6i = Struct(mapfmt(op2._endian + b'6i', self.size))
        nelements = (len(data) - n) // ntotal
        for unused_i in range(nelements):
            edata = data[n:n + ntotal]
            out = struct_6i.unpack(edata)
            if op2.is_debug_file:
                op2.binary_debug.write('  CSHEAR=%s\n' % str(out))
            #(eid, pid, n1, n2, n3, n4) = out
            elem = CSHEAR.add_op2_data(out)
            self.add_op2_element(elem)
            n += ntotal
        op2.card_count['CSHEAR'] = nelements
        return n

# CSLOT3
# CSLOT4

    def read_ctetrap(self, data: bytes, n: int) -> int:
        """
        CTETP(12201,122,9013)    - the marker for Record 87
        .. todo:: needs work
        """
        op2: OP2Geom = self.op2
        op2.log.info('poor reading of CTETRAP in GEOM2')
        nelements = (len(data) - n) // 108  # 27*4
        struct_27i = Struct(op2._endian + b'27i')
        for unused_i in range(nelements):
            edata = data[n:n+108]
            out = struct_27i.unpack(edata)
            if op2.is_debug_file:
                op2.binary_debug.write('  CTETP=%s\n' % str(out))

            #print(out)
            eid, pid, n1, n2, n3, n4 = out[:6]
            #(eid, pid, n1, n2, n3, n4,
             #e1, e2, e3, e4, e5, e6, e7, e8, e9, e10, e11, e12,
             #f1, f2, f3, f4, b1, ee1, ee2, ee3, ee4) = out
            #print("out = ",out)
            #e = [e1, e2, e3, e4, e5, e6, e7, e8, e9, e10, e11, e12]
            #f = [f1, f2, f3, f4]
            #ee = [ee1, ee2, ee3, ee4]

            #print("e  = ",e)
            #print("f  = ",f)
            #print("b1  = ",b1)
            #print("ee = ",ee)
            data_in = [eid, pid, n1, n2, n3, n4]
            elem = CTETRA4.add_op2_data(data_in)
            self.add_op2_element(elem)
            n += 108
        return n

    def read_ctetra(self, data: bytes, n: int) -> int:
        """
        CTETRA(5508,55,217)      - the marker for Record 88
        CTETPR(7609,76,9993)     - the marker for Record 89
        CTETR10F(16600,166,9999) - the marker for Record 90
        CTETR4FD(16100,161,9999) - the marker for Record 91
        """
        op2: OP2Geom = self.op2
        ntotal = 48 * self.factor  # 12*4
        s = Struct(mapfmt(op2._endian + b'12i', self.size))
        nelements = (len(data) - n) // ntotal
        for unused_i in range(nelements):
            edata = data[n:n + ntotal]
            out = s.unpack(edata)
            if op2.is_debug_file:
                op2.binary_debug.write('  CTETRA=%s\n' % str(out))
            (eid, pid, n1, n2, n3, n4, n5, n6, n7, n8, n9, n10) = out

            data_in = [eid, pid, n1, n2, n3, n4]
            big_nodes = [n5, n6, n7, n8, n9, n10]
            if sum(big_nodes) > 0:
                elem = CTETRA10.add_op2_data(data_in + big_nodes)
            else:
                elem = CTETRA4.add_op2_data(data_in)
            try:
                elem.validate()
                str(elem)
            except Exception:
                print(data_in, big_nodes)
                raise
            self.add_op2_element(elem)
            n += ntotal
        op2.card_count['CTETRA'] = nelements
        return n

# CTQUAD - 92
# CTTRIA - 93

    def read_ctria3(self, data: bytes, n: int) -> int:
        """Common method for reading CTRIA3s"""
        card_name = 'CTRIA3'
        card_obj = CTRIA3
        methods = {
            52 : self._read_ctria3_52,
            56 : self._read_ctria3_56,
        }
        try:
            n = self._read_double_card(card_name, card_obj, self.add_op2_element,
                                       methods, data, n)
        except DoubleCardError:
            raise
        return n

    def _read_ctria3_52(self, card_obj, data: bytes, n: int) -> int:
        """
        CTRIA3(5959,59,282) - the marker for Record 94

        """
        op2: OP2Geom = self.op2
        ntotal = 52 * self.factor  # 13*4
        s = Struct(mapfmt(op2._endian + b'5iff3i3f', self.size))
        nelements = (len(data) - n)// ntotal
        elements = []
        for unused_i in range(nelements):
            edata = data[n:n+ntotal]
            out = s.unpack(edata)
            #print('eid=%s pid=%s n1=%s n2=%s n3=%s theta=%s zoffs=%s '
                  #'blank1=%s blank2=%s tflag=%s t1=%s t2=%s t3=%s' % (
                      #eid, pid, n1, n2, n3, theta, zoffs,
                      #blank1, blank2, tflag, t1, t2, t3))
            (eid, pid, n1, n2, n3, theta, zoffs, unused_blank1,
             unused_blank2, tflag, t1, t2, t3) = out
            if op2.is_debug_file:
                op2.binary_debug.write('  CTRIA3=%s\n' % str(out))

            theta_mcid = convert_theta_to_mcid(theta)
            data_in = [eid, pid, n1, n2, n3, theta_mcid, zoffs, tflag, t1, t2, t3]
            elem = CTRIA3.add_op2_data(data_in)
            elements.append(elem)
            n += ntotal
        return n, elements

    def _read_ctria3_56(self, card_obj, data: bytes, n: int) -> int:
        r"""
        CTRIA3(5959,59,282) - the marker for Record 94
            eid  pid n1     n2     n3     theta?      ?  ?  ?  ?
            201, 35, 50000, 50001, 50101, 19.2, 0, 0, 0, 0, 0.2, 0.2, 0.2, -1,

        C:\MSC.Software\msc_nastran_runs\cc179h.op2
        ints = (
            5959, 59, 282,
            5i                            f     4i          3f             i
            eid  pid n1     n2     n3     theta ?  ?  ?  ?  t1   t2   t3
            201, 35, 50000, 50001, 50101, 19.2, 0, 0, 0, 0, 0.2, 0.2, 0.2, -1,
            202, 35, 50101, 50100, 50000, 19.2, 0, 0, 0, 0, 0.2, 0.2, 0.2, -1,
            203, 35, 50001, 50002, 50102, 19.2, 0, 0, 0, 0, 0.2, 0.2, 0.2, -1,
            204, 35, 50102, 50101, 50001, 19.2, 0, 0, 0, 0, 0.2, 0.2, 0.2, -1,
            205, 35, 50002, 50003, 50103, 19.2, 0, 0, 0, 0, 0.2, 0.2, 0.2, -1, ...)

        """
        op2: OP2Geom = self.op2
        ntotal = 56 * self.factor  # 13*4
        s = Struct(mapfmt(op2._endian + b'5i f 4i 3f i', self.size))
        nelements = (len(data) - n)// ntotal
        cids = set(list(self.op2.coords.keys()))
        check_cids = max(cids) > 0 if len(cids) else False
        elements = []
        for unused_i in range(nelements):
            edata = data[n:n+ntotal]
            out = s.unpack(edata)
            #print('eid=%s pid=%s n1=%s n2=%s n3=%s theta=%s zoffs=%s '
                  #'blank1=%s blank2=%s tflag=%s t1=%s t2=%s t3=%s' % (
                      #eid, pid, n1, n2, n3, theta, zoffs,
                      #blank1, blank2, tflag, t1, t2, t3))
            (eid, pid, n1, n2, n3, theta, a, b, c, d,
             t1, t2, t3, minus1_mcid) = out
            abcd = (a, b, c, d)
            assert abcd == (0, 0, 0, 0), abcd
            if op2.is_debug_file:
                op2.binary_debug.write('  CTRIA3=%s\n' % str(out))

            zoffs = 0.0
            tflag = 0
            theta_mcid = theta
            if minus1_mcid != -1:
                theta_mcid = minus1_mcid
                assert theta == 0., theta
                if check_cids:
                    assert minus1_mcid in cids, minus1_mcid

            data_in = [eid, pid, n1, n2, n3, theta_mcid, zoffs, tflag, t1, t2, t3]
            elem = CTRIA3.add_op2_data(data_in)
            elements.append(elem)
            n += ntotal
        return n, elements


# CTRIAFD - 95

    def read_ctria6(self, data: bytes, n: int) -> int:
        """
        common method for reading CTRIA6

        CTRIA6(4801,48,327) # MSC 2005 - GEOM201
        Word Name Type Description
        1  EID    I Element identification number
        2  PID    I Property identification number
        3  G(6)   I Grid point identification numbers of connection points
        9 THETA  RS Material property orientation angle or coordinate system identification number
        10 ZOFFS RS Offset from the surface of grid points reference plane
        11 T(3)  RS Membrane thickness of element at grid points

        Record 90 -- CTRIA6(4801,48,327) # MSC 2005 - GEOM2
        CTRIA6(4801,48,327)
        Word Name Type Description
        1 EID     I Element identification number
        2 PID     I Property identification number
        3 G(6)    I Grid point identification numbers of connection points
        9 THETA  RS Material property orientation angle or coordinate system identification number
        10 ZOFFS RS Offset from the surface of grid points reference plane
        11 T(3)  RS Membrane thickness of element at grid points
        14 TFLAG  I Relative thickness flag
        """
        card_name = 'CTRIA6'
        card_obj = CTRIA6
        methods = {
            52 : self._read_ctria6_v2001_52,
            56 : self._read_ctria6_current_56,
            60 : self._read_ctria6_60,
        }
        try:
            n = self._read_double_card(card_name, card_obj, self.add_op2_element,
                                       methods, data, n)
        except DoubleCardError:
            raise

        #n = self._read_split_card(data, n,
                                  #self._read_ctria6_current, self._read_ctria6_v2001,
                                  #'CTRIA6', CTRIA6, self.add_op2_element)
        return n

    def read_ctria3fd(self, data: bytes, n: int) -> int:
        """
        data= (
            16200, 16201, 16201, 16202, 16203, 0, 0, 0, 0, -1,
            16201, 16201, 16201, 16203, 16204, 0, 0, 0, 0, -1) - 10

        """
        """
        Common method for reading CTRIA3s

        """
        card_name = 'CTRIA3'
        card_obj = CTRIA3
        methods = {
            32 : self._read_ctria3fd_32,
            40 : self._read_ctria3fd_40,
        }
        try:
            n = self._read_double_card(card_name, card_obj, self.add_op2_element,
                                       methods, data, n)
        except DoubleCardError:
            raise
        return n

    def _read_ctria3fd_40(self, card_obj, data: bytes, n: int) -> int:
        op2: OP2Geom = self.op2
        ntotal = 40 * self.factor
        nelements = (len(data) - n) // ntotal  # 10*4
        assert (len(data) - n) % ntotal == 0
        elements = []
        s = Struct(op2._endian + b'10i')
        for unused_i in range(nelements):
            edata = data[n:n + ntotal]
            out = s.unpack(edata)
            if op2.is_debug_file:
                op2.binary_debug.write('  CTRIA3=%s\n' % str(out))
            (eid, pid, n1, n2, n3, n4, n5, n6, a, minus1) = out
            assert (a, minus1) == (0, -1), (a, minus1)
            assert minus1 == -1
            assert n4 == 0, out
            assert n5 == 0, out
            assert n6 == 0, out
            nids = [n1, n2, n3]
            elem = CTRIA3(eid, pid, nids,
                          theta_mcid=0., zoffset=0., tflag=0,
                          T1=None, T2=None, T3=None, comment='')
            elements.append(elem)
            n += ntotal
        return n, elements

    def _read_ctria3fd_32(self, card_obj, data: bytes, n: int) -> int:
        """
        Word Name Type Description
        1 EID  I Element identification number
        2 PID  I Property identification number
        3 G(6) I Grid point identification numbers of connection points
        """
        op2: OP2Geom = self.op2
        ntotal = 32 * self.factor
        nelements = (len(data) - n) // ntotal  # 8*4
        assert (len(data) - n) % ntotal == 0
        s = Struct(op2._endian + b'8i')
        elements = []
        for unused_i in range(nelements):
            edata = data[n:n + ntotal]
            out = s.unpack(edata)
            if op2.is_debug_file:
                op2.binary_debug.write('  CTRIA3=%s\n' % str(out))
            (eid, pid, n1, n2, n3, n4, n5, n6) = out
            assert n4 == 0, out
            assert n5 == 0, out
            assert n6 == 0, out
            nids = [n1, n2, n3]
            elem = CTRIA3(eid, pid, nids,
                          theta_mcid=0., zoffset=0., tflag=0,
                          T1=None, T2=None, T3=None, comment='')
            elements.append(elem)
            n += ntotal
        return n, elements

    def read_ctriax3fd(self, data: bytes, n: int) -> int:
        """
        Common method for reading CTRIAX3
        """
        card_name = 'CTRIAX'
        card_obj = CTRIAX
        methods = {
            32 : self._read_ctriax3fd_32,
            40 : self._read_ctriax3fd_40,
        }
        try:
            n = self._read_double_card(card_name, card_obj, self.add_op2_element,
                                       methods, data, n)
        except DoubleCardError:
            raise
        return n

    def _read_ctriax3fd_32(self, card_obj, data: bytes, n: int) -> int:
        """
        Word Name Type Description
        1 EID  I Element identification number
        2 PID  I Property identification number
        3 G(6) I Grid point identification numbers of connection points

        (331, 111, 331, 332, 333, 0, 0, 0,
         332, 11, 331, 333, 334, 0, 0, 0)
        """
        op2: OP2Geom = self.op2
        s = Struct(op2._endian + b'8i')
        ntotal = 32 * self.factor
        nelements = (len(data) - n) // ntotal  # 8*4
        assert (len(data) - n) % ntotal == 0
        elements = []
        for unused_i in range(nelements):
            edata = data[n:n + ntotal]
            out = s.unpack(edata)
            if op2.is_debug_file:
                op2.binary_debug.write('  CTRIAX=%s\n' % str(out))
            (eid, pid, n1, n2, n3, n4, n5, n6) = out
            nids = [n1, n2, n3, n4, n5, n6]
            assert n4 == 0, out
            assert n5 == 0, out
            assert n6 == 0, out
            elem = CTRIAX(eid, pid, nids, theta_mcid=0., comment='')
            elements.append(elem)
            n += ntotal
        return n, elements

    def _read_ctriax3fd_40(self, card_obj, data: bytes, n: int) -> int:
        """
        Word Name Type Description
        1 EID  I Element identification number
        2 PID  I Property identification number
        3 G(6) I Grid point identification numbers of connection points
        ?
        ?

        ints    = (16800, 168, 9978,
                   eid    pid    n1     n2     n3     4  5  6  ?   ?
                   16800, 16801, 16801, 16802, 16803, 0, 0, 0, 0, -1,
                   16801, 16801, 16801, 16803, 16804, 0, 0, 0, 0, -1)
        """
        op2: OP2Geom = self.op2
        ntotal = 40 * self.factor
        nelements = (len(data) - n) // ntotal  # 8*4
        assert (len(data) - n) % ntotal == 0

        elements = []
        s = Struct(op2._endian + b'10i')
        for unused_i in range(nelements):
            edata = data[n:n + ntotal]
            out = s.unpack(edata)
            if op2.is_debug_file:
                op2.binary_debug.write('  CTRIAX=%s\n' % str(out))
            (eid, pid, n1, n2, n3, n4, n5, n6, dunno, minus1) = out
            nids = [n1, n2, n3, n4, n5, n6]
            assert n4 == 0, out
            assert n5 == 0, out
            assert n6 == 0, out
            assert dunno == 0, dunno
            assert minus1 == -1, minus1
            elem = CTRIAX(eid, pid, nids, theta_mcid=0., comment='')
            elements.append(elem)
            n += ntotal
        return n, elements

    def read_ctria6fd(self, data: bytes, n: int) -> int:
        """
        Common method for reading CTRIA6s
        """
        card_name = 'CTRIA6'
        card_obj = CTRIA6
        methods = {
            32 : self._read_ctria6fd_32,
            40 : self._read_ctria6fd_40,
        }
        try:
            n = self._read_double_card(card_name, card_obj, self.add_op2_element,
                                       methods, data, n)
        except DoubleCardError:
            raise
        return n

    def _read_ctria6fd_40(self, card_obj, data: bytes, n: int) -> int:
        """
        Word Name Type Description
        1 EID  I Element identification number
        2 PID  I Property identification number
        3 G(6) I Grid point identification numbers of connection points
        ?
        ?

        ints    = (16700, 167, 9981,
                   eid    pid    n1     n2     n3     n4     n5     n6
                   16700, 16701, 16701, 16702, 16703, 16705, 16706, 16709, 0, -1,
                   16701, 16701, 16701, 16703, 16704, 16709, 16707, 16708, 0, -1)
        """
        op2: OP2Geom = self.op2
        ntotal = 40 * self.factor
        s = Struct(op2._endian + b'10i')
        nelements = (len(data) - n) // ntotal  # 8*4
        assert (len(data) - n) % ntotal == 0
        elements = []
        for unused_i in range(nelements):
            edata = data[n:n + ntotal]
            out = s.unpack(edata)
            if op2.is_debug_file:
                op2.binary_debug.write('  CTRIA6=%s\n' % str(out))
            (eid, pid, n1, n2, n3, n4, n5, n6, dunno, minus1) = out
            assert dunno == 0, dunno
            assert minus1 == -1, minus1
            nids = [n1, n2, n3, n4, n5, n6]
            #out = (eid, pid, n1, n2, n3, n4, n5, n6, theta, zoffs, t1, t2, t3, 0)
            elem = CTRIA6(eid, pid, nids,
                          theta_mcid=0., zoffset=0., tflag=0,
                          T1=None, T2=None, T3=None, comment='')
            elements.append(elem)
            n += ntotal
        return n, elements

    def _read_ctria6fd_32(self, card_obj, data: bytes, n: int) -> tuple[int, list[CTRIA6]]:
        """
        Word Name Type Description
        1 EID  I Element identification number
        2 PID  I Property identification number
        3 G(6) I Grid point identification numbers of connection points
        """
        op2: OP2Geom = self.op2
        ntotal = 32 * self.factor
        s = Struct(op2._endian + b'8i')
        nelements = (len(data) - n) // ntotal  # 8*4
        assert (len(data) - n) % ntotal == 0
        elements = []
        for unused_i in range(nelements):
            edata = data[n:n + ntotal]
            out = s.unpack(edata)
            if op2.is_debug_file:
                op2.binary_debug.write('  CTRIA6=%s\n' % str(out))
            (eid, pid, n1, n2, n3, n4, n5, n6) = out
            nids = [n1, n2, n3, n4, n5, n6]
            #out = (eid, pid, n1, n2, n3, n4, n5, n6, theta, zoffs, t1, t2, t3, 0)
            elem = CTRIA6(eid, pid, nids,
                          theta_mcid=0., zoffset=0., tflag=0,
                          T1=None, T2=None, T3=None, comment='')
            elements.append(elem)
            n += ntotal
        return n, elements

    def _read_ctria6_60(self, card_obj, data: bytes, n: int) -> tuple[int, list[CTRIA6]]:
        """
        CTRIA6(4801,48,327) - the marker for Record 96

        Record 90 -- CTRIA6(4801,48,327) # MSC 2005 - GEOM2
        Word Name Type Description
        1 EID     I Element identification number
        2 PID     I Property identification number
        3 G(6)    I Grid point identification numbers of connection points
        9 THETA  RS Material property orientation angle or coordinate system identification number
        10 ZOFFS RS Offset from the surface of grid points reference plane
        11 T(3)  RS Membrane thickness of element at grid points
        14 TFLAG  I Relative thickness flag
        -1
        """
        op2: OP2Geom = self.op2
        ntotal = 60 * self.factor # 15*4
        s = Struct(mapfmt(op2._endian + b'8i 5f i i', self.size))
        nelements = (len(data) - n) // ntotal
        assert (len(data) - n) % ntotal == 0
        elements = []
        for unused_i in range(nelements):
            edata = data[n:n + ntotal]
            out = s.unpack(edata)
            if op2.is_debug_file:
                op2.binary_debug.write('  CTRIA6=%s\n' % str(out))
            (eid, pid, n1, n2, n3, n4, n5, n6, theta, zoffs, t1, t2, t3, tflag, minus1) = out
            #print('eid=%s pid=%s nids[%s, %s %s] theta=%s zoffs=%s '
                  #'tflag=%s t1=%s t2=%s t3=%s' % (
                      #eid, pid, n1, n2, n3, theta, zoffs,
                      #tflag, t1, t2, t3))
            tflag, minus1 = out[-2:]
            assert minus1 == -1
            #op2.log.info('ctria6 tflag = %s' % tflag)
            #print(minus1)
            data_in = [eid, pid, n1, n2, n3, n4, n5, n6, theta, zoffs, t1, t2, t3, tflag]
            elem = CTRIA6.add_op2_data(data_in)
            self.add_op2_element(elem)
            assert tflag in [-1, 0, 1], tflag
            elements.append(elem)
            n += ntotal
        return n, elements

    def _read_ctria6_current_56(self, card_obj, data: bytes, n: int) -> tuple[int, list[CTRIA6]]:
        """
        CTRIA6(4801,48,327) - the marker for Record 96

        Record 90 -- CTRIA6(4801,48,327) # MSC 2005 - GEOM2
        Word Name Type Description
        1 EID     I Element identification number
        2 PID     I Property identification number
        3 G(6)    I Grid point identification numbers of connection points
        9 THETA  RS Material property orientation angle or coordinate system identification number
        10 ZOFFS RS Offset from the surface of grid points reference plane
        11 T(3)  RS Membrane thickness of element at grid points
        14 TFLAG  I Relative thickness flag
        """
        op2: OP2Geom = self.op2
        ntotal = 56 * self.factor # 14*4
        nelements = (len(data) - n) // ntotal
        assert (len(data) - n) % ntotal == 0
        elements = []
        s = Struct(mapfmt(op2._endian + b'8i 5f i', self.size))
        for unused_i in range(nelements):
            edata = data[n:n + ntotal]
            out = s.unpack(edata)
            if op2.is_debug_file:
                op2.binary_debug.write('  CTRIA6=%s\n' % str(out))
            #print('eid=%s pid=%s n1=%s n2=%s n3=%s theta=%s zoffs=%s '
                  #'blank1=%s blank2=%s tflag=%s t1=%s t2=%s t3=%s' % (
                      #eid, pid, n1, n2, n3, theta, zoffs,
                      #blank1, blank2, tflag, t1, t2, t3))
            #(eid, pid, n1, n2, n3, n4, n5, n6, theta, zoffs, t1, t2, t3, tflag) = out
            tflag = out[-1]
            #op2.log.info('ctria6 tflag = %s' % tflag)
            elem = CTRIA6.add_op2_data(out)
            self.add_op2_element(elem)
            assert tflag in [-1, 0, 1], tflag
            elements.append(elem)
            n += ntotal
        return n, elements

    def _read_ctria6_v2001_52(self, card_obj, data: bytes, n: int) -> int:
        """
        CTRIA6(4801,48,327) - the marker for Record 96

        1 EID     I Element identification number
        2 PID     I Property identification number
        3 G(6)    I Grid point identification numbers of connection points
        9 THETA  RS Material property orientation angle or coordinate system identification number
        10 ZOFFS RS Offset from the surface of grid points reference plane
        11 T(3)  RS Membrane thickness of element at grid points
        14 TFLAG  I Relative thickness flag
        """
        op2: OP2Geom = self.op2
        ntotal = 52 * self.factor
        s = Struct(op2._endian + b'8i 5f')
        nelements = (len(data) - n) // ntotal # 13*4
        assert (len(data) - n) % ntotal == 0
        elements = []
        for unused_i in range(nelements):
            edata = data[n:n + 52]
            out = s.unpack(edata)
            if op2.is_debug_file:
                op2.binary_debug.write('  CTRIA6=%s\n' % str(out))
            #print('eid=%s pid=%s n1=%s n2=%s n3=%s theta=%s zoffs=%s '
                  #'blank1=%s blank2=%s tflag=%s t1=%s t2=%s t3=%s' % (
                      #eid, pid, n1, n2, n3, theta, zoffs,
                      #blank1, blank2, tflag, t1, t2, t3))
            (eid, pid, n1, n2, n3, n4, n5, n6, theta, zoffs, t1, t2, t3) = out
            out = (eid, pid, n1, n2, n3, n4, n5, n6, theta, zoffs, t1, t2, t3, 0)
            elem = CTRIA6.add_op2_data(out)
            elements.append(elem)
            n += ntotal
        return n, elements

# CTRIA6FD
# CTRIAP

    def read_ctriar(self, data: bytes, n: int) -> int:
        """
        CTRIAR(9200,92,385) - the marker for Record 99
        """
        card_name = 'CTRIAR'
        card_obj = CTRIAR
        methods = {
            # nbytes
            52 : self._read_ctriar_13,
            56 : self._read_ctriar_14,
        }
        try:
            n = self._read_double_card(card_name, card_obj, self.add_op2_element,
                                       methods, data, n)
        except DoubleCardError:
            raise
            #self.op2.log.warning(f'try-except {card_name}')
            #n = self._read_split_card(data, n,
                                      #self._read_cquad8_current, self._read_cquad8_v2001,
                                      #card_name, op2.add_op2_element)
        #nelements = op2.card_count['CQUAD8']
        #op2.log.debug(f'nCQUAD8 = {nelements}')

        #n = self._read_dual_card(data, n, self._read_ctriax_8, self._read_ctriax_9,
                                 #'CTRIAX', op2.add_op2_element)
        return n

    def _read_ctriar_13(self, element: CTRIAR, data: bytes, n: int) -> tuple[int, list[CTRIAR]]:
        op2: OP2Geom = self.op2
        ntotal = 52 * self.factor  # 13*4
        s = Struct(mapfmt(op2._endian + b'5iff3i3f', self.size))
        nelements = (len(data) - n) // ntotal
        assert (len(data) - n) % ntotal == 0
        elements = []
        for unused_i in range(nelements):
            edata = data[n:n+ntotal]
            #self.show_data(data[n:n+ntotal+40])
            out = s.unpack(edata)
            (eid, pid, n1, n2, n3, theta, zoffs, unused_blank1,
             unused_blank2, tflag, t1, t2, t3) = out
            #print('eid=%s pid=%s nodes=(%s,%s,%s) theta=%s zoffs=%s '
                  #'blank1=%s blank2=%s tflag=%s t1-t3=(%s,%s,%s)' % (
                      #eid, pid, n1, n2, n3, theta, zoffs,
                      #unused_blank1, unused_blank2, tflag, t1, t2, t3))
            if op2.is_debug_file:
                op2.binary_debug.write('  CTRIAR=%s\n' % str(out))
            data_in = [eid, pid, n1, n2, n3, theta, zoffs, tflag, t1, t2, t3]
            elem = CTRIAR.add_op2_data(data_in)
            elements.append(elem)
            n += ntotal
        return n, elements

    def _read_ctriar_14(self, element: CTRIAR, data: bytes, n: int) -> tuple[int, list[CTRIAR]]:
        """same as ``read_ctriar_13`` but with a -1 to change the format"""
        op2: OP2Geom = self.op2
        ntotal = 56 * self.factor  # 14*4
        s = Struct(mapfmt(op2._endian + b'5iff3i3f i', self.size))
        nelements = (len(data) - n)// ntotal
        assert (len(data) - n) % ntotal == 0
        elements = []
        for unused_i in range(nelements):
            edata = data[n:n+ntotal]
            out = s.unpack(edata)
            (eid, pid, n1, n2, n3, theta, zoffs, unused_blank1,
             unused_blank2, tflag, t1, t2, t3, minus1) = out
            #print('eid=%s pid=%s nodes=(%s,%s,%s) theta=%s zoffs=%s '
                  #'blank1=%s blank2=%s tflag=%s t1-t3=(%s,%s,%s)' % (
                      #eid, pid, n1, n2, n3, theta, zoffs,
                      #unused_blank1, unused_blank2, tflag, t1, t2, t3))
            assert minus1 == -1, minus1
            if op2.is_debug_file:
                op2.binary_debug.write('  CTRIAR=%s\n' % str(out))
            data_in = [eid, pid, n1, n2, n3, theta, zoffs, tflag, t1, t2, t3]
            elem = CTRIAR.add_op2_data(data_in)
            elements.append(elem)
            n += ntotal
        return n, elements

    #def _read_ctriax_b(self, data: bytes, n: int) -> int:  # pragma: no cover
        #"""
        #CTRIAX  341     11      341     342     343     345     346     349     +TX1
        #+TX1    12
        #CTRIAX  342     111     341     343     344     349     347     348
        #(341, 11,  341, 342, 343, 345, 346, 349,
         #342, 111, 341, 343, 344, 349, 347, 348)
        #"""
        #ntotal = 32  # 8*4
        #nentries = (len(data) - n) // ntotal
        #struc = Struct(op2._endian + b'8i')
        #for unused_i in range(nentries):
            #edata = data[n:n + 32]
            #out = struc.unpack(edata)
            #if op2.is_debug_file:
                #op2.binary_debug.write('  CTRIAX=%s\n' % str(out))
            #eid, pid, n1, n2, n3, n4, n5, n6 = out
            #nids = [n1, n2, n3, n4, n5, n6]
            #elem = CTRIAX(eid, pid, nids, theta_mcid=0., comment='')
            #op2.add_op2_element(elem)
            #n += 32
        #op2.card_count['CTRIAX'] = nentries
        #return n

    def read_ctriax(self, data: bytes, n: int) -> int:
        """common method for reading CTRIAXs"""
        card_name = 'CTRIAX'
        card_obj = CTRIAX
        methods = {
            32 : self._read_ctriax_8,
            36 : self._read_ctriax_9,
            40 : self._read_ctriax_10,
        }
        try:
            n = self._read_double_card(card_name, card_obj, self.add_op2_element,
                                       methods, data, n)
        except DoubleCardError:
            raise
            #self.op2.log.warning(f'try-except {card_name}')
            #n = self._read_split_card(data, n,
                                      #self._read_cquad8_current, self._read_cquad8_v2001,
                                      #card_name, op2.add_op2_element)
        #nelements = op2.card_count['CQUAD8']
        #op2.log.debug(f'nCQUAD8 = {nelements}')

        #n = self._read_dual_card(data, n, self._read_ctriax_8, self._read_ctriax_9,
                                 #'CTRIAX', op2.add_op2_element)
        return n

    def _read_ctriax_8(self, card_obj, data: bytes, n: int) -> tuple[int, list[CTRIAX]]:
        """(10108, 101, 512)"""
        op2: OP2Geom = self.op2
        ntotal = 32 * self.factor  # 9*4
        struc = Struct(op2._endian + b'8i')

        nentries = (len(data) - n) // ntotal
        assert (len(data) - n) % ntotal == 0
        assert nentries > 0

        elems = []
        for unused_i in range(nentries):
            edata = data[n:n + ntotal]
            out = struc.unpack(edata)
            if op2.is_debug_file:
                op2.binary_debug.write('  CTRIAX=%s\n' % str(out))
            eid, pid, n1, n2, n3, n4, n5, n6 = out
            nids = [n1, n2, n3, n4, n5, n6]
            elem = CTRIAX(eid, pid, nids, theta_mcid=0., comment='no theta set')
            elems.append(elem)
            n += ntotal
        return n, elems

    def _read_ctriax_9(self, card_obj, data: bytes, n: int) -> tuple[int, list[CTRIAX]]:
        """(10108, 101, 512)"""
        op2: OP2Geom = self.op2
        ntotal = 36 * self.factor  # 9*4
        struc = Struct(op2._endian + b'9i')

        nentries = (len(data) - n) // ntotal
        assert (len(data) - n) % ntotal == 0
        assert nentries > 0

        elems = []
        for unused_i in range(nentries):
            edata = data[n:n + ntotal]
            out = struc.unpack(edata)
            if op2.is_debug_file:
                op2.binary_debug.write('  CTRIAX=%s\n' % str(out))
            eid, pid, n1, n2, n3, n4, n5, n6, unused_undef1 = out
            nids = [n1, n2, n3, n4, n5, n6]
            elem = CTRIAX(eid, pid, nids, theta_mcid=0., comment='no theta set')
            elems.append(elem)
            n += ntotal
        return n, elems

    def _read_ctriax_10(self, card_obj, data: bytes, n: int) -> tuple[int, list[CTRIAX]]:
        r"""(10108, 101, 512)

        C:\MSC.Software\msc_nastran_runs\el705ce.op2
        data = (16900, 169, 9977,
        eid    pid    n1     n2     n3     n4     n5     n6     ?   ?
        16900, 16901, 16901, 16902, 16903, 16905, 16906, 16909, 0, -1,
        16901, 16901, 16901, 16903, 16904, 16909, 16907, 16908, 0, -1
        """
        op2: OP2Geom = self.op2
        ntotal = 40 * self.factor  # 10*4
        struc = Struct(op2._endian + b'10i')

        nentries = (len(data) - n) // ntotal
        assert (len(data) - n) % ntotal == 0
        assert nentries > 0

        elems = []
        for unused_i in range(nentries):
            edata = data[n:n + ntotal]
            out = struc.unpack(edata)
            if op2.is_debug_file:
                op2.binary_debug.write('  CTRIAX=%s\n' % str(out))
            eid, pid, n1, n2, n3, n4, n5, n6, dunno, minus1 = out
            assert dunno == 0, dunno
            assert minus1 == -1, minus1
            nids = [n1, n2, n3, n4, n5, n6]
            elem = CTRIAX(eid, pid, nids, theta_mcid=0., comment='no theta set')
            elems.append(elem)
            n += ntotal
        return n, elems

    def read_ctriax6(self, data: bytes, n: int) -> int:  # 101
        """(6108, 61, 107)"""
        op2: OP2Geom = self.op2
        ntotal = 44 * self.factor  # 11*4
        nentries = (len(data) - n) // ntotal
        assert (len(data) - n) % ntotal == 0
        assert nentries > 0
        struc = Struct(mapfmt(op2._endian + b'8i f ii', self.size))
        for unused_i in range(nentries):
            edata = data[n:n + ntotal]
            out = struc.unpack(edata)
            if op2.is_debug_file:
                op2.binary_debug.write('  CTRIAX6=%s\n' % str(out))
            elem = CTRIAX6.add_op2_data(out)
            self.add_op2_element(elem)
            n += ntotal
        op2.card_count['CTRIAX6'] = nentries
        return n

# CTRIX3FD - 102
# CTRIX6FD - 103

    def read_gmbndc(self, data: bytes, n: int) -> int:
        """
        GMBNDC(3201,32,478)

        Word Name      Type  Description
        1    BID       I     Boundary identification number
        2    GRIDI     I     Initial grid identification number for boundary
        3    GRIDF     I     Final grid identification number for boundary
        4    ENTITY(2) CHAR4 Entity type for defining boundary
        6    EID       I     Entity identification numbers for boundary of subdomain
        Word 6 repeats until End of Record
        """
        raise UnsupportedCard('GMNBDC')
        #op2: OP2Geom = self.op2
        #op2.log.info('geom skipping GMBNDC in GEOM2')
        ##self.show_data(data)
        ##(1, 31, 32, GRID____, -1,
         ##2, 41, 42, GRID____, -1)

        ##ints= (3201, 32, 478,
        ## 2, 41, 42, 1145390406, 538985799, 41, -1,
        ## 990003, 101000045, 101000046, 1145655879, 538976288, 101000045, 101000046, -1)
        #ints = np.frombuffer(data[n:], op2.idtype) # .tolist()
        #isplit = np.where(ints == -1)[0]
        #nelements = len(isplit)

        #i0 = 0
        #for ispliti in isplit:
            #eid, gridi, gridf = ints[i0:i0+3]
            ##print(eid, gridi, gridf)
            #s0 = n + (i0 + 3) * 4
            #s1 = s0 + 8
            #entity = data[s0:s1].decode('latin1').rstrip()
            #eids = ints[i0+5:ispliti]
            #assert entity in ['FEEDGE', 'GRID', 'GMCURV', 'GMCURVE'], f'entity={entity!r}'
            ##print(eids)
            #i0 = ispliti + 1
        #op2.card_count['GMBNDC'] = nelements
        #return len(data)

    def read_ctube(self, data: bytes, n: int) -> int:
        """
        CTUBE(3701,37,49) - the marker for Record 104
        """
        op2: OP2Geom = self.op2
        ntotal = 16 * self.factor  # 4*4
        struct_4i = Struct(mapfmt(op2._endian + b'4i', self.size))
        nelements = (len(data) - n) // ntotal
        assert (len(data) - n) % ntotal == 0
        assert nelements > 0
        for unused_i in range(nelements):
            edata = data[n:n + ntotal]
            out = struct_4i.unpack(edata)
            if op2.is_debug_file:
                op2.binary_debug.write('  CTUBE=%s\n' % str(out))
            #(eid, pid, n1, n2) = out
            elem = CTUBE.add_op2_data(out)
            self.add_op2_element(elem)
            n += ntotal
        op2.card_count['CTUBE'] = nelements
        return n

    def read_cvisc(self, data: bytes, n: int) -> int:
        """CVISC(3901,39,50) - the marker for Record 105"""
        op2: OP2Geom = self.op2
        ntotal = 16 * self.factor  # 4*4
        struct_4i = Struct(mapfmt(op2._endian + b'4i', self.size))
        nelements = (len(data) - n) // ntotal
        assert (len(data) - n) % ntotal == 0
        assert nelements > 0
        for unused_i in range(nelements):
            edata = data[n:n + ntotal]
            out = struct_4i.unpack(edata)
            if op2.is_debug_file:
                op2.binary_debug.write('  CVISC=%s\n' % str(out))
            #(eid, pid, n1, n2) = out
            element = CVISC.add_op2_data(out)
            self.add_op2_element(element)
            n += ntotal
        op2.card_count['CVISC'] = nelements
        return n

    def read_cweld(self, data: bytes, n: int) -> int:
        """
        CWELD(11701,117,559) - Record 106
        same as CFAST
        """
        op2: OP2Geom = self.op2
        op2.log.info('geom skipping CWELD in GEOM2')
        if op2.is_debug_file:
            op2.binary_debug.write('geom skipping CWELD in GEOM2\n')
        return len(data)

    def read_cweldc(self, data: bytes, n: int) -> int:  # 107
        op2: OP2Geom = self.op2
        op2.log.info('geom skipping CWELDC in GEOM2')
        if op2.is_debug_file:
            op2.binary_debug.write('geom skipping CWELDC in GEOM2\n')
        return len(data)

    def read_cweldg(self, data: bytes, n: int) -> int:  # 108
        op2: OP2Geom = self.op2
        op2.log.info('geom skipping CWELDG in GEOM2')
        if op2.is_debug_file:
            op2.binary_debug.write('geom skipping CWELDG in GEOM2\n')
        return len(data)

# TDOO: above are checked by DMAP...
#-------------------------------
# CWSEAM
    def read_genel(self, data: bytes, n: int) -> int:
        r"""
        Word Name Type Description
        1 EID I Element identification number
        2 UI  I Independent grid point identification number
        3 CI  I Component number
        Words 2 through 3 repeat until End of Record

        4 M(C) I Number of rows and columns in K or Z and rows in S
        5 UD   I Dependent grid point identification number
        6 CD   I Component number
        Words 5 through 6 repeat until End of Record

        7 N(C)  I Number of columns in S
        8 F     I 1 means Z, 2 means K
        9 KZIJ RS Lower triangular terms of the K or Z
        matrix. See Notes.
        Word 9 repeats MM times
        10 NZERO(C) I
        NZERO =1 Actually " 0"
        11 SIJ RS Terms of the S matrix
        Word 11 repeats M times
        Word 11 repeats N times
        NZERO =0
        End NZERO
        12 UNDEF none
        Word 12 repeats until End of Record

        # C:\NASA\m4\formats\git\examples\move_tpl\ha145c.op2
        (432, # eid
        1, 3, 2, 3, 3, 3, 4, 3, 5, 3, 6, 3, 7, 3, 8, 3, 9, 3, 10, 3, -1, # (ui,ci)
        10, # M(c)
        11, 3, 11, 4, 11, 5, 11, 6, -1, # (ud,cd)
        4, 1, # N(c), f, KZij...floats...)

        (6.05360936588321e-43,
        1, 3,
        2.802596928649634e-45, 3,
        3, 3,
        4, 3,
        7.006492321624085e-45, 3,
        8.407790785948902e-45, 3,
        9.80908925027372e-45, 3,
        1.1210387714598537e-44, 3,
        1.2611686178923354e-44, 3,
        1.401298464324817e-44, 3,
        nan,
        1.401298464324817e-44, 1.5414283107572988e-44, 3, 1.5414283107572988e-44, 4, 1.5414283107572988e-44, 7.006492321624085e-45, 1.5414283107572988e-44, 8.407790785948902e-45,
        nan,
        4, 1,
        8.71720021677902e-06, 1.3361000128497835e-06, 1.2778000382240862e-05, 6.272000064200256e-06, 1.6251000488409773e-05, 1.0492000001249835e-05, 2.0478000806178898e-05, 1.562999932502862e-05, 2.428500010864809e-05, 2.0403000235091895e-05,
        3.086099968641065e-05, 6.272000064200256e-06, 3.229700087103993e-05, 1.0492000001249835e-05, 3.352899875608273e-05, 1.562999932502862e-05, 3.502099934848957e-05,
        2.025700086960569e-05, 3.578500036383048e-05, 2.7731999580282718e-05, 1.572600012877956e-05, 4.825499854632653e-05, 3.762800042750314e-05, 7.328399806283414e-05,
        6.433799717342481e-05, 9.580999903846532e-05, 8.837800123728812e-05, 6.374900112859905e-05, 3.762800042750314e-05, 8.013600017875433e-05, 6.433799717342481e-05,
        0.00010011999984271824, 8.837800123728812e-05, 0.00011811000149464235, 0.00012758000229950994, 0.00011344000085955486, 0.00019350000366102904, 0.0001816000003600493,
        0.0002528300101403147, 0.00024294000468216836, 0.0001699900021776557, 0.0001816000003600493, 0.000229199999012053, 0.00024294000468216836, 0.0002824899856932461,
        0.00036862000706605613, 0.00035051998565904796, 0.0005267499946057796, 0.0005117100081406534, 0.00042292001307941973, 0.0005117100081406534, 0.0005718700122088194,
        0.0008483999990858138, 0.0008233999833464622, 0.0009233999880962074, 4, 1.0, 90.0, -20.25, 45.0, 1.0, 90.0, 81.0, 45.0, 1.0, 186.0, -17.850000381469727,
        141.0, 1.0, 186.0, 71.4000015258789, 141.0, 1.0, 268.0, -15.800000190734863, 223.0, 1.0, 268.0, 63.20000076293945, 223.0, 1.0, 368.0, -13.300000190734863, 323.0, 1.0,
        368.0, 53.20000076293945, 323.0, 1.0, 458.0, -11.050000190734863, 413.0, 1.0, 458.0, 44.20000076293945, 413.0)
        """
        self.op2.log.info('geom skipping GENEL in GEOM2')
        #op2.log.info(f'geom skipping GENEL in GEOM2; len(data)={len(data)-12}')
        #print(n)
        ints = np.frombuffer(data[n:], dtype='int32').copy()
        #floats = np.frombuffer(data[n:], dtype='float32').copy()
        i = 0
        while i < len(ints):
            #1 EID I Element identification number
            eid = ints[i]

            iminus1 = np.where(ints[i+1:] == -1)[0]
            #print('iminus1', iminus1)
            idelta0 = iminus1[0]
            idelta1 = iminus1[1]
            # print('idelta0', idelta0)
            uc = ints[i+1:i+1+idelta0].reshape(idelta0//2, 2)
            print(uc)
            j = i + 1 + idelta0 + 1
            #2 UI  I Independent grid point identification number
            #3 CI  I Component number
            #Words 2 through 3 repeat until End of Record

            nrows = ints[j]
            print('nrows=', nrows)
            mucd = ints[j:i+1+idelta1]
            mc = mucd[0]
            nucd = len(mucd) - 1
            ucd = mucd[1:].reshape(nucd//2, 2)
            print(f'M(c) = {mc}')
            print(ucd)
            i = i + 1 + idelta1 + 1
            #4 M(C) I Number of rows and columns in K or Z and rows in S
            #5 UD   I Dependent grid point identification number
            #6 CD   I Component number
            #Words 5 through 6 repeat until End of Record

            # ---------------
            print('-----------------')

            #7 N(C)  I Number of columns in S (4)
            #8 F     I 1 means Z, 2 means K -> Z

            #9 KZIJ RS Lower triangular terms of the K or Z matrix. See Notes.
            #Word 9 repeats MM times
            #  10 NZERO(C) I
            #  NZERO =1 Actually " 0"
            #    11 SIJ RS Terms of the S matrix
            #    Word 11 repeats M times
            # Word 11 repeats N times
            # NZERO =0
            #End NZERO
            #12 UNDEF none
            #Word 12 repeats until End of Record

            nc = ints[i]
            f = ints[i+1]
            i += 2
            print(f'nc={nc} f={f}')
            #print(ints[i:].min())
            #print(ints[i+55])
            #print(floats[i:i+55].tolist())
            #print(ints[i+55:].tolist())
            #print(ints[i:].tolist())
            #print(floats[i:])
            #print(len(floats[i:]))
            break
        #self.show_data(data[12:])
        return len(data)
# GMDNDC
# GMBNDS
# GMINTC
# GMINTS
    def read_micpnt(self, data: bytes, n: int) -> int:
        """
        RECORD – MICPNT(2801,28,630)
        Word Name Type Description
        1 EID I Element identification number
        2 GID I Fluid grid identification number
        3 DESC(12) CHAR4 Description - 48 characters maximum
        """
        op2: OP2Geom = self.op2
        #size = self.size
        struc = Struct(op2._endian + b'2i 48s')
        ntotal = 8 + 48
        nelements = (len(data) - n) // ntotal
        assert (len(data) - n) % ntotal == 0
        assert nelements > 0
        for unused_i in range(nelements):
            edata = data[n:n + ntotal]  # 4*4
            out = struc.unpack(edata)
            if op2.is_debug_file:
                op2.binary_debug.write('  MICPNT=%s\n' % str(out))
            #(eid,n1,n2) = out
            eid, node_id, name_bytes = out
            name = name_bytes.decode('latin1').rstrip()
            op2.add_micpnt(eid, node_id, name)
            n += ntotal
        op2.card_count['MICPNT'] = nelements
        return n

    def read_plotel(self, data: bytes, n: int) -> int:  # 114
        """(5201, 52, 11)"""
        op2: OP2Geom = self.op2
        struct_3i = Struct(op2._endian + b'3i')
        ntotal = 12 * self.factor
        nelements = (len(data) - n) // ntotal
        assert (len(data) - n) % ntotal == 0
        assert nelements > 0
        for unused_i in range(nelements):
            edata = data[n:n + ntotal]  # 4*4
            out = struct_3i.unpack(edata)
            if op2.is_debug_file:
                op2.binary_debug.write('  PLOTEL=%s\n' % str(out))
            #(eid,n1,n2) = out
            elem = PLOTEL.add_op2_data(out)
            op2._add_methods._add_plotel_object(elem)
            n += ntotal
        op2.card_count['PLOTEL'] = nelements
        return n

    def read_radbc(self, data: bytes, n: int) -> int:
        """
        RADBC(12801,128,417)

        Word Name Type Description
        1 EID      I Element identification number
        2 FAMB    RS Radiation view factor between the face and the ambient point
        3 CNTRLND  I Control point for radiation boundary condition
        4 NODAMB   I
        """
        op2: OP2Geom = self.op2
        #C:\NASA\m4\formats\git\examples\move_tpl\ht15339.op2
        #(-99, 1.0, 0, 101)
        #radbc   101     1.0             -99
        #RADBC NODAMB   FAMB CNTRLND     EID1 EID2 EID3
        structi = Struct(op2._endian + b'ifii')
        ntotal = 16 * self.factor
        nelements = (len(data) - n) // ntotal
        assert (len(data) - n) % ntotal == 0
        assert nelements > 0
        for unused_i in range(nelements):
            edata = data[n:n + ntotal]  # 4*4
            out = structi.unpack(edata)
            if op2.is_debug_file:
                op2.binary_debug.write('  RADBC=%s\n' % str(out))
            eid, famb, cntrlnd, nodamb = out
            eids = eid
            boundary_condition = RADBC(nodamb, famb, cntrlnd, eids)
            op2._add_methods._add_thermal_bc_object(boundary_condition, boundary_condition.nodamb)
            n += ntotal
        op2.card_count['RADBC'] = nelements
        return n

# RADINT
    def read_sint(self, data: bytes, n: int) -> int:
        """
        Word Name Type Description
        1 EID    I Element identification number
        2 PID    I Property identification number
        3 PTELE  I Pointer to element identification number
        4 NSEG   I Number of segments
        5 STSC   I Stride for segment displacement data
        6 PTSC   I Pointer to segment displacements
        7 NBOUND I Number of boundaries
        8 BID    I Boundary identification number
        9 NFACE  I Number of faces
        10 STBC  I Stride for boundary displacement data
        11 NSEG  I Number of segments
        12 STLC1 I Stride for Boundary Lagrange Multiplier data
        13 PTBND I Pointer to boundary identification number
        14 PTBC  I Pointer to boundary displacements
        15 PTLC  I Pointer to boundary Lagrange
        Multipliers
        Words 8 through 15 repeat 5 times
        16 UNDEF(3 ) none
        """
        self.op2.log.info('geom skipping SINT in GEOM2')
        # C:\NASA\m4\formats\git\examples\move_tpl\ifscp88.op2
        # doesn't seem to be a card, more of a general info on the geometry...
        #ints = np.frombuffer(data[n:], dtype=op2.idtype).copy()
        #print(ints.tolist())
        return len(data)

    def read_spoint(self, data: bytes, n: int) -> int:
        """
        (5551,49,105)    - the marker for Record 118
        """
        op2: OP2Geom = self.op2
        ntotal = 4 * self.factor
        npoints = (len(data) - n) // ntotal
        nids = np.frombuffer(data[n:], op2.idtype8).tolist()
        if op2.is_debug_file:
            op2.binary_debug.write('SPOINT=%s\n' % nids)
        spoint = SPOINTs.add_op2_data(nids)
        op2._add_methods._add_spoint_object(spoint)
        op2.card_count['SPOINT'] = npoints
        return len(data)

    def read_vubeam(self, data: bytes, n: int) -> int:  # 119
        """(11601, 116, 9942)"""
        deltae = 100000000
        #deltan = 111000000 # 111001002
        def element(eid, pid, nids):
            x = [1., 0., 0.]
            g0 = None
            nids2 = [nid - deltae for nid in nids]
            elem = CBEAM(eid-deltae, pid, nids2, x, g0)
            return elem
        element.type = 'CBEAM' #  I love python
        self._run_2nodes(element, data, n)
        #assert len(self.elements) > 0, self.elements
        return n

# VUHEXA
# VUQUAD4
# VUPENTA
# VUTETRA
# VUTRIA
# VUBEAM
# VUHEXA
# VUQUAD4
# WELDP

    def read_ctrax3(self, data: bytes, n: int) -> int:
        """
        RECORD - CTRAX3(6111,61,996)

        Word Name Type Description
        1 EID    I Element identification number
        2 PID    I Property identification number
        3 G(3)   I Grid point identification numbers of connection points
        4 THETA RS Material property orientation angle
        """
        op2: OP2Geom = self.op2
        ntotal = 24 * self.factor  # 6*4
        s = Struct(mapfmt(op2._endian + b'5if', self.size))
        nelements = (len(data) - n)// ntotal
        for unused_i in range(nelements):
            edata = data[n:n+ntotal]
            out = s.unpack(edata)
            (eid, pid, n1, n2, n3, theta) = out
            if op2.is_debug_file:
                op2.binary_debug.write('  CTRAX3=%s\n' % str(out))
            #data_in = [eid, pid, n1, n2, n3, theta]
            elem = CTRAX3(eid, pid, [n1, n2, n3], theta)
            self.add_op2_element(elem)
            n += ntotal
        op2.card_count['CTRAX3'] = nelements
        return n

    #def _read_cquadx4(self, data: bytes, n: int) -> int:
        #op2.log.info('geom skipping CQUADX4 in GEOM2')
        #return len(data)

    def read_ctrax6(self, data: bytes, n: int) -> int:
        """
        RECORD - CTRAX6(6113, 61, 998)

        Word Name Type Description
        1 EID    I Element identification number
        2 PID    I Property identification number
        3 G(6)   I Grid point identification numbers of connection points
        4 THETA RS Material property orientation angle
        """
        op2: OP2Geom = self.op2
        ntotal = 36 * self.factor  # 9*4
        s = Struct(mapfmt(op2._endian + b'8if', self.size))
        nelements = (len(data) - n)// ntotal
        for unused_i in range(nelements):
            edata = data[n:n+ntotal]
            out = s.unpack(edata)
            (eid, pid, n1, n2, n3, n4, n5, n6, theta) = out
            if op2.is_debug_file:
                op2.binary_debug.write('  CTRAX6=%s\n' % str(out))
            #data_in = [eid, pid, n1, n2, n3, n4, n5, n6, theta]
            elem = CTRAX6(eid, pid, [n1, n2, n3, n4, n5, n6], theta)
            self.add_op2_element(elem)
            n += ntotal
        op2.card_count['CTRAX6'] = nelements
        return n

    def read_cquadx8(self, data: bytes, n: int) -> int:
        """
        RECORD - CQUADX8(6114, 61, 999)

        Word Name Type Description
        1 EID    I Element identification number
        2 PID    I Property identification number
        3 G(8)   I Grid point identification numbers of connection points
        4 THETA RS Material property orientation angle
        """
        op2: OP2Geom = self.op2
        ntotal = 44 * self.factor  # 11*4
        s = Struct(mapfmt(op2._endian + b'10if', self.size))
        nelements = (len(data) - n)// ntotal
        for unused_i in range(nelements):
            edata = data[n:n+ntotal]
            out = s.unpack(edata)
            (eid, pid, n1, n2, n3, n4, n5, n6, n7, n8, theta) = out
            if op2.is_debug_file:
                op2.binary_debug.write('  CQUADX8=%s\n' % str(out))
            #data_in = [eid, pid, n1, n2, n3, n4, n5, n6, n7, n8, theta]
            elem = CQUADX8(eid, pid, [n1, n2, n3, n4, n5, n6, n7, n8], theta)
            self.add_op2_element(elem)
            n += ntotal
        op2.card_count['CQUADX8'] = nelements
        return n

    def read_feface_pval(self, data: bytes, n: int) -> int:
        r"""
        C:\NASA\m4\formats\git\examples\move_tpl\pshp02.bdf
        ints = (
           100001001, 100, 111001001, 111001002, 111001004, 111001007, 111001008, 111001010,
           100001002, 100, 111001002, 111001003, 111001005, 111001008, 111001009, 111001011,
           100001003, 100, 111001005, 111001004, 111001002, 111001011, 111001010, 111001008,
           100001004, 100, 111001004, 111001005, 111001006, 111001010, 111001011, 111001012,
           100001005, 100, 111001007, 111001008, 111001010, 111001013, 111001014, 111001016,
           100001006, 100, 111001008, 111001009, 111001011, 111001014, 111001015, 111001017,
           100001007, 100, 111001011, 111001010, 111001008, 111001017, 111001016, 111001014,
           100001008, 100, 111001010, 111001011, 111001012, 111001016, 111001017, 111001018,
           100002001, 100, 111002001, 111002002, 111002004, 111002007, 111002008, 111002010,
           100002002, 100, 111002002, 111002003, 111002005, 111002008, 111002009, 111002011,
           100002003, 100, 111002005, 111002004, 111002002, 111002011, 111002010, 111002008,
           100002004, 100, 111002004, 111002005, 111002006, 111002010, 111002011, 111002012,
           100002005, 100, 111002007, 111002008, 111002010, 111002013, 111002014, 111002016,
           100002006, 100, 111002008, 111002009, 111002011, 111002014, 111002015, 111002017,
           100002007, 100, 111002011, 111002010, 111002008, 111002017, 111002016, 111002014,
           100002008, 100, 111002010, 111002011, 111002012, 111002016, 111002017, 111002018)

        C:\NASA\m4\formats\git\examples\move_tpl\pet1126.op2
        ints = (
           100001001, 1, 111001001, 111001002, 111001004, 111001007, 111001008, 111001010,
           100001002, 1, 111001002, 111001003, 111001005, 111001008, 111001009, 111001011,
           100001003, 1, 111001005, 111001004, 111001002, 111001011, 111001010, 111001008,
           100001004, 1, 111001004, 111001005, 111001006, 111001010, 111001011, 111001012,
           100001005, 1, 111001007, 111001008, 111001010, 111001013, 111001014, 111001016,
           100001006, 1, 111001008, 111001009, 111001011, 111001014, 111001015, 111001017,
           100001007, 1, 111001011, 111001010, 111001008, 111001017, 111001016, 111001014,
           100001008, 1, 111001010, 111001011, 111001012, 111001016, 111001017, 111001018)
        """
        #self.show_data(data[12:])
        return len(data)

    def read_feedge2(self, data: bytes, n: int) -> int:
        """
        (2901, 29, 9601)

        Word Name Type Description
        1 EDGEID     I Edge identification number
        2 GRID1      I Identification number of end GRID 1
        3 GRID2      I Identification number of end GRID 2
        4 CID        I Coordinate system identification number
        5 GEOMIN CHAR4 Type of referencing entry: "GMCURV" or "POINT"
        6 GEOMID1    I Identification number of a POINT or GMCURV entry
        7 GEOMID2    I Identification number of a POINT or GMCURV entry
        """
        op2: OP2Geom = self.op2
        # C:\NASA\m4\formats\git\examples\move_tpl\phsflux4.op2
        #(200000002, 3, 1002, 6, 12, 0, 0)
        # FEEDGE EDGEID GRID1 GRID2 CIDBC GEOMIN ID1 ID2
        #FEEDGE    1002    6     12
        # no FEFACE...

        # weird because the order is wrong and there are two extra FEFACE lines
        #C:\NASA\m4\formats\git\examples\move_tpl\phs19332.op2
        #(200000002, 3, 2, 74, 57, 24, 0)
        #(200000003, 3, 1, 7, 24, 57, 0)
        #FEFACE   1       7       24      57
        #FEFACE   2       74      57      24
        #FEFACE   3       1       51      18
        #FEFACE   4       68      18      51

        # C:\NASA\m4\formats\git\examples\move_tpl\ptsahd.op2
        #feedge,311, 311,331, ,gmcurv,31
        #feedge,321, 331,321, ,gmcurv,31
        #feface,21, 121,122,126,125,,12
        #feface,22, 125,126,123,  ,,12
        #feface,23, 124,123,126,  ,,12
        #(200000004, 5, 12, 121, 122, 126, 125)
        #(200000004, 5, 12, 123, 124, 126, 0)
        #(200000004, 5, 12, 123, 125, 126, 0)
        #(200000002, 4, 22, 221, 222, 0, 0)
        #(200000003, 1, 321, 321, 0, 0, 0)
        #(200000005, 3, 22, 125, 126, 123, 0)
        #(200000006, 3, 23, 124, 123, 126, 0)

        # C:\NASA\m4\formats\git\examples\move_tpl\phscvhg6.op2
        #C:\NASA\m4\formats\git\examples\move_tpl\phsconv4.op2


        #self.show_data(data[12:])
        ntotal = 28 * self.factor  # 7*4
        # s = Struct(op2._endian + b'4i 4s 2i') #expected
        s = Struct(op2._endian + b'7i')
        ndatai = len(data) - n
        nelements = ndatai // ntotal  # 7*4
        assert ndatai % ntotal == 0
        for unused_i in range(nelements):
            edata = data[n:n+ntotal]
            n += ntotal
            out = s.unpack(edata)
            #print(out)
            #edge_id, n1, n2, cid, geomin, geom1, geom2 = out # expected
            dunno, nfields, edge_id, n1, n2, n3, n4 = out
            assert nfields in [1, 2, 3, 4, 5], out
            #assert zero1 == 0, f'zero1={zero1} out={out}'
            #assert zero2 == 0, f'zero2={zero2} out={out}'
            if op2.is_debug_file:
                op2.binary_debug.write('  FEEDGE=%s\n' % str(out))

            geomin_str = 'POINT' # ???
            cid = 0
            geom1 = 0
            geom2 = 0
            if nfields == 2:
                edge = FEEDGE(edge_id, [n1, n2], cid, [geom1, geom2], geomin=geomin_str)
                if edge_id in op2.feedge:
                    edge_old = op2.feedge[edge_id]
                    if edge != edge_old:
                        msg = f'Duplicate FEEDGE\nold:\n{edge_old}\nnew:\n{edge}'
                        raise RuntimeError(msg)
                    continue
            feedge = op2.add_feedge(edge_id, [n1, n2], cid, [geom1, geom2], geomin=geomin_str)
            str(feedge)
            #elif nfields in [3, 4, 5]:
                #if nfields == 3:
                    #nids = [n1, n2]
                #elif nfields == 4:
                    #nids = [n1, n2, n3]
                #elif nfields == 5:
                    #nids = [n1, n2, n3, n4]
                ##elem = FEFACE(edge_id, nids)

            #data_in = [eid, pid, n1, n2, n3, n4, n5, n6, n7, n8, theta]
            # elem = CQUADX8(eid, pid, [n1, n2, n3, n4, n5, n6, n7, n8], theta)
            # self.add_op2_element(elem)
        op2.card_count['FEEDGE'] = nelements
        return n

    def read_gmbnds(self, data: bytes, n: int) -> int:
        """
        Word Name Type Description
        1 BID           I Boundary identification number
        2 GRIDC(4)      I Corner grid 1
        6 ENTITY(2) CHAR4 Entity type for defining boundary
        8 EID           I Entity identification numbers for boundary of subdomain
        Word 8 repeats until End of Record
        """
        raise UnsupportedCard('GMBNDS')
        #op2: OP2Geom = self.op2
        #op2.log.info('geom skipping GMBNDS in GEOM2')
        ##(1, 0, 0, 0, 0, 'FEFACE  ', 31, -1)
        #ints = np.frombuffer(data[n:], dtype=op2.idtype).copy()
        #iminus1 = np.where(ints == -1)[0]
        #i0 = 0
        #for iminus1i in iminus1:
            #bid, n1, n2, n3, n4 = ints[i0:i0+5]
            #s0 = n + (i0 + 5) * 4
            #s1 = s0 + 8
            #entity = data[s0:s1].decode('latin1').rstrip()
            #assert entity in ['FEFACE', 'GMSURF', 'GRID'], (bid, n1, n2, n3, n4, entity)
            #assert bid >= 0, (bid, n1, n2, n3, n4, entity)
            #eids = ints[i0+7:iminus1i]
            ##print(bid, n1, n2, n3, n4)
            ##print('entity = %r' % entity)
            ##print(eid)
            ##print('-----')
            #i0 = iminus1i + 1
        #return len(data)

    def read_cngret(self, data: bytes, n: int) -> int:
        # C:\NASA\m4\formats\git\examples\move_tpl\bpas101.op2
        # C:\NASA\m4\formats\git\examples\move_tpl\pass8.op2
        return len(data)

    def read_adapt(self, data: bytes, n: int) -> int:  # pragma: no cover
        raise UnsupportedCard('GMCORD')
        self.op2.log.info('geom skipping ADAPT card in GEOM2')
        return len(data)

    def read_cseam_maybe(self, data: bytes, n: int) -> int:
        """
        564 / 4 = 141
        141
        3, 47

          ints    = (
          77, 2011808,
          77, 2007308,
          8, 100001002, 4007101, 100001003, 4007101, 100001004, 4007101, 100001005, 4007101, 100001006, 4007101, 100001007, 4007101, 100001008, 4007101, 100001009, 4007101,
          78, 2011808,
          78, 2007308,
          4, 100001010, 4007101, 100001011, 4007101, 100001012, 4007101, 100001013, 4007101,
          79, 2011808,
          79, 2007308,
          4, 100001014, 4007101, 100001015, 4007101, 100001016, 4007101, 100001017, 4007101,
          87, 2011808,
          87, 2007308,
          8, 100001018, 4007101, 100001019, 4007101, 100001020, 4007101, 100001021, 4007101, 100001022, 4007101, 100001023, 4007101, 100001024, 4007101, 100001025, 4007101,
          88, 2011808,
          88, 2007308,
          4, 100001026, 4007101, 100001027, 4007101, 100001028, 4007101, 100001029, 4007101,
          89, 2011808,
          89, 2007308,
          4, 100001030, 4007101, 100001031, 4007101, 100001032, 4007101, 100001033, 4007101,
          97, 2011808,
          97, 2007308,
          8, 100001034, 4007101, 100001035, 4007101, 100001036, 4007101, 100001037, 4007101, 100001038, 4007101, 100001039, 4007101, 100001040, 4007101, 100001041, 4007101, # 2*8
          98, 2011808,
          98, 2007308,
          4, 100001042, 4007101, 100001043, 4007101, 100001044, 4007101, 100001045, 4007101,  # 2*4
          99, 2011808,
          99, 2007308,
          4, 100001046, 4007101, 100001047, 4007101, 100001048, 4007101, 100001049, 4007101)  # 2*4
        """
        self.op2.log.info('geom skipping CSEAM? card in GEOM2')
        #self.show_data(data[n:], types='ifqds')
        #self.show_data(data[n+4:], types='qd')
        #ntotal = 44 * self.factor  # 11*4
        ##s = Struct(mapfmt(op2._endian + b'10if', self.size))
        #nelements = (len(data) - n)// ntotal
        #for unused_i in range(nelements):
            #edata = data[n:n+ntotal]
            ##out = s.unpack(edata)
            #n += ntotal
        #aaa
        return len(data)

def convert_theta_to_mcid(theta):
    """odd function..."""
    # sort of guessed at this number...it seems reasonable-ish
    if theta > 511.:
        # per DMAP...you couldn't make a new record number?
        # theta = 512. * (cid + 1)
        # theta/512 = cid + 1
        # cid = theta/512. - 1
        #
        cid_float = theta / 512. - 1
        cid = int(cid_float)
        assert np.allclose(cid, cid_float), 'theta=%s cid=%s cid_float=%s' % (theta, cid, cid_float)
        theta = cid
    return theta

def get_minus_4_index(idata):
    """helper for ``get_minus_4_index``"""
    #print('idata =', idata)
    i = np.where((idata == -4) | (idata == -3))[0]
    if len(i) == 0:
        return len(idata)
    return i[0]
