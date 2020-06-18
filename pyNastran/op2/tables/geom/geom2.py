"""
defines readers for BDF objects in the OP2 GEOM2/GEOM2S table
"""
# pylint: disable=C0103
from struct import Struct
from functools import partial
from typing import Tuple, List, Union, Any
import numpy as np

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
from pyNastran.bdf.cards.elements.bars import CBAR, CBEND
from pyNastran.bdf.cards.elements.beam import CBEAM
from pyNastran.bdf.cards.elements.mass import (CONM1, CONM2, CMASS1, CMASS2,
                                               CMASS3, CMASS4)
from pyNastran.bdf.cards.elements.solid import (CTETRA4, CPYRAM5, CPENTA6, CHEXA8,
                                                CTETRA10, CPYRAM13, CPENTA15, CHEXA20,)
from pyNastran.bdf.cards.thermal.thermal import CHBDYG, CONV, CHBDYP, CHBDYE, CONVM
from pyNastran.bdf.cards.thermal.radiation import RADBC # , RADM, RADCAV, RADLST, RADMTX, VIEW, VIEW3D
from pyNastran.bdf.cards.nodes import SPOINTs
from pyNastran.bdf.cards.elements.bush import CBUSH
from pyNastran.bdf.cards.parametric.geometry import FEEDGE
from pyNastran.bdf.cards.elements.acoustic import CHACAB, CHACBR, CAABSF
from pyNastran.op2.errors import MixedVersionCard
from pyNastran.op2.tables.geom.geom_common import GeomCommon
from pyNastran.op2.op2_interface.op2_reader import mapfmt # , reshape_bytes_block


class GEOM2(GeomCommon):
    """defines methods for reading op2 elements"""

    def _read_geom2_4(self, data, ndata):
        return self._read_geom_4(self._geom2_map, data, ndata)

    def __init__(self):
        GeomCommon.__init__(self)
        self._geom2_map = {
            # per dmap-nx-10.pdf or nx12.pdf
            (15200, 152, 9912): ['ACFACE3', self._read_fake],
            (15500, 155, 9913): ['ACFACE4', self._read_fake],
            (15600, 156, 9914): ['ACFACE6', self._read_fake],
            (15700, 157, 9915): ['ACFACE8', self._read_fake],
            (2601, 26, 0): ['BEAMAERO', self._read_fake],
            (2708, 27, 59): ['CAABSF', self._read_caabsf],
            (2108, 21, 224): ['CAXIF2', self._read_fake],
            (2208, 22, 225): ['CAXIF3', self._read_fake],
            (2308, 23, 226): ['CAXIF4', self._read_fake],

            (2408, 24, 180): ['CBAR', self._read_cbar],         # record 8
            (4001, 40, 275): ['CBARAO', self._read_cbarao],     # record 9  - not done
            (5408, 54, 261): ['CBEAM', self._read_cbeam],       # record 10
            (11401, 114, 9016): ['CBEAMP', self._read_cbeamp],  # record 11 - not done
            (4601, 46, 298): ['CBEND', self._read_cbend],     # record 12 - not done
            (2608, 26, 60): ['CBUSH', self._read_cbush],      # record 13
            (5608, 56, 218): ['CBUSH1D', self._read_cbush1d], # record 14 - not done
            (5609, 60, 9899): ['CBUSH1DNL', self._read_fake],
            (14801, 148, 956): ['CCHOCK3', self._read_fake],
            (14901, 149, 957): ['CCHOCK4', self._read_fake],
            (15001, 150, 958): ['CCHOCK6', self._read_fake],
            (15101, 151, 959): ['CCHOCK8', self._read_fake],
            (2315, 23, 146): ['CCONE-10', self._read_ccone],  # nx10
            (2315, 23, 0): ['CCONE-12', self._read_fake],     # nx12
            (201, 2, 69): ['CDAMP1', self._read_cdamp1],
            (301, 3, 70): ['CDAMP2', self._read_cdamp2],
            (401, 4, 71): ['CDAMP3', self._read_cdamp3],
            (501, 5, 72): ['CDAMP4', self._read_cdamp4],
            (10608, 106, 404): ['CDAMPS', self._read_cdamp5],
            (6208, 62, 108): ['CDUM2', self._read_fake],
            (6308, 63, 109): ['CDUM3', self._read_fake],
            (6408, 64, 110): ['CDUM4', self._read_fake],
            (6508, 65, 111): ['CDUM5', self._read_fake],
            (6608, 66, 112): ['CDUM6', self._read_fake],
            (6708, 67, 113): ['CDUM7', self._read_fake],
            (6808, 68, 114): ['CDUM8', self._read_cdum8],
            (6908, 69, 115): ['CDUM9', self._read_cdum9],
            (601, 6, 73): ['CELAS1', self._read_celas1],
            (6010, 53, 9900): ['CELAS1NL', self._read_fake],
            (701, 7, 74): ['CELAS2', self._read_celas2],
            (7010, 20, 9898): ['CELAS2NL', self._read_fake],
            (801, 8, 75): ['CELAS3', self._read_celas3],
            (901, 9, 76): ['CELAS4', self._read_celas4],

            (9801, 98, 506): ['CFAST-nx10', self._read_cfast],  # nx10
            (13801, 138, 566): ['CFAST-nx12', self._read_cfast], # nx12

            (8515, 85, 209): ['CFLUID2-nx10', self._read_cfluid2],  # nx10
            (8615, 86, 210): ['CFLUID3-nx10', self._read_cfluid3],  # nx10
            (8715, 87, 211): ['CFLUID4-nx10', self._read_cfluid4],  # nx10
            (8515, 85, 0): ['CFLUID2-nx12', self._read_cfluid2], # nx12
            (8615, 86, 0): ['CFLUID3-nx12', self._read_cfluid3], # nx12
            (8715, 87, 0): ['CFLUID4-nx12', self._read_cfluid4], # nx12

            (1908, 19, 104): ['CGAP', self._read_cgap],
            #(13101, 131, 9901): ['CGPLSTN3', self._read_fake],
            #(13401, 134, 9904): ['CGPLSTN8', self._read_fake],
            #(13301, 133, 9903): ['CGPLSTN6', self._read_fake],
            #(13201, 132, 9902): ['CGPLSTN4', self._read_fake],
            (8308, 83, 405): ['CHBDYE', self._read_chbdye],
            (10808, 108, 406): ['CHBDYG', self._read_chbdyg],
            (10908, 109, 407): ['CHBDYP', self._read_chbdyp],
            #(7308, 73, 253): ['CHEXA', self._read_fake],
            #(16300, 163, 9989): ['CHEXA20F', self._read_fake],
            #(14100, 141, 9990): ['CHEXAF', self._read_fake],
            #(14000, 140, 9990): ['CHEXAFD', self._read_fake],
            #(7708, 77, 9944): ['CHEXAL', self._read_fake],
            #(11801, 118, 907): ['CHEXCZ', self._read_fake],
            #(12001, 120, 9011): ['CHEXP', self._read_fake],
            #(7409, 74, 9991): ['CHEXPR', self._read_fake],
            (1001, 10, 65): ['CMASS1', self._read_cmass1],    # record 52
            (1101, 11, 66): ['CMASS2', self._read_cmass2],    # record 53
            (1201, 12, 67): ['CMASS3', self._read_cmass3],    # record 54
            (1301, 13, 68): ['CMASS4', self._read_cmass4],    # record 55
            (2508, 25, 0): ['CMFREE', self._read_cmfree],     # record 56 - not done
            (1401, 14, 63): ['CONM1', self._read_conm1],      # record 57
            (1501, 15, 64): ['CONM2', self._read_conm2],      # record 58
            (1601, 16, 47): ['CONROD', self._read_conrod],    # record 59
            (12701, 127, 408): ['CONV', self._read_conv],
            (8908, 89, 422): ['CONVM', self._read_convm],
            (12101, 121, 9012): ['CPENP', self._read_fake],
            (4108, 41, 280): ['CPENTA', self._read_cpenta],
            #(14200, 142, 9906): ['CPENTAF', self._read_fake],
            #(7108, 71, 9943): ['CPENTAL', self._read_fake],
            (7509, 75, 9992): ['CPENPR', self._read_fake],
            #(16500, 165, 9987): ['CPENT15F', self._read_fake],
            #(16000, 160, 9988): ['CPENT6FD', self._read_fake],
            #(11901, 119, 908): ['CPENTCZ', self._read_fake],
            (1701, 17, 980): ['CPLSTN3', self._read_fake_nx],
            (5701, 57, 981): ['CPLSTN4', self._read_fake_nx],
            (5801, 58, 982): ['CPLSTN6', self._read_fake_nx],
            (7201, 72, 983): ['CPLSTN8', self._read_fake_nx],
            (8801, 88, 984): ['CPLSTS3', self._read_fake_nx],
            (8401, 84, 985): ['CPLSTS4', self._read_fake_nx],
            (1801, 18, 986): ['CPLSTS6', self._read_fake_nx],
            (3601, 36, 987): ['CPLSTS8', self._read_fake_nx],
            (17200, 172, 1000) : ['CPYRAM', self._read_cpyram], # nx-specific
            #(14400, 144, 9908): ['CPYRAMF', self._read_fake], # nx-specific
            (25700, 257, 9948) : ['CPYRA5FD', self._read_cpyram], # nx-specific
            (25800, 258, 9947) : ['CPYRA13F', self._read_cpyram], # nx-specific
            (7909, 79, 9946) : ['CPYRAMPR', self._read_cpyram], # nx-specific
            (17000, 170, 9980): ['CQDX4FD', self._read_cquad],
            (17100, 171, 9979): ['CQDX9FD', self._read_cquadx],
            #(25110, 170, 9951): ['CQDX4FDN', self._read_fake],
            #(25310, 171, 9949): ['CQDX8FDN', self._read_fake],
            (9108, 91, 507): ['CQUAD', self._read_cquad],
            (2958, 51, 177): ['CQUAD4', self._read_cquad4],
            #(14600, 146, 9910): ['CQUADF', self._read_fake],
            (13900, 139, 9984): ['CQUAD4FD', self._read_cquad],
            #(4701, 47, 326): ['CQUAD8', self._read_fake],
            #(3302, 33, 1694): ['CQUAD8L', self._read_fake],
            #(15901, 159, 9956): ['CQUAD8N', self._read_fake],
            (16400, 164, 9983): ['CQUAD9FD', self._read_cquad],
            #(11101, 111, 9014): ['CQUADP', self._read_fake],
            #(8009, 80, 367): ['CQUADR', self._read_fake],
            #(13002, 130, 1692): ['CQUADRL', self._read_fake],
            #(15401, 154, 9954): ['CQUADRN', self._read_fake],
            #(9008, 90, 508): ['CQUADX', self._read_fake],
            (6112, 61, 997): ['CQUADX4', self._read_fake],
            (6114, 61, 999): ['CQUADX8', self._read_cquadx8],
            (3001, 30, 48): ['CROD', self._read_crod],         # record 81
            #(14500, 145, 9909): ['CRODF', self._read_fake],
            (3501, 35, 1): ['CSBOLT', self._read_fake],
            (3101, 31, 61): ['CSHEAR', self._read_cshear],     # record 84
            (4408, 44, 227): ['CSLOT3', self._read_fake],
            (4508, 45, 228): ['CSLOT4', self._read_fake],
            #(12201, 122, 9013): ['CTETP', self._read_fake],
            #(5508, 55, 217): ['CTETRA', self._read_fake],
            #(14300, 143, 9907): ['CTETRAF', self._read_fake],
            #(7609, 76, 9993): ['CTETPR', self._read_fake],
            #(16600, 166, 9985): ['CTETR10F', self._read_fake],
            #(16100, 161, 9986): ['CTETR4FD', self._read_fake],
            #(25010, 168, 9952): ['CTRAX3FDN', self._read_fake],
            #(25210, 169, 9950): ['CTRAX6FDN', self._read_fake],
            #(5959, 59, 282): ['CTRIA3', self._read_fake],
            #(14700, 147, 9911): ['CTRIAF', self._read_fake],
            #(16200, 162, 9982): ['CTRIA3FD', self._read_fake],
            #(4801, 48, 327): ['CTRIA6', self._read_fake],
            #(16700, 167, 9981): ['CTRIA6FD', self._read_fake],
            #(3202, 32, 1693): ['CTRIA6L', self._read_fake],
            #(15801, 158, 9955): ['CTRIA6N', self._read_fake],
            #(11301, 113, 9015): ['CTRIAP', self._read_fake],
            #(9200, 92, 385): ['CTRIAR', self._read_fake],
            #(12902, 129, 1691): ['CTRIARL', self._read_fake],
            #(15301, 153, 9953): ['CTRIARN', self._read_fake],
            #(10108, 101, 512): ['CTRIAX', self._read_fake],
            #(6108, 61, 107): ['CTRIAX6', self._read_fake],
            (6111, 61, 996): ['CTRAX3', self._read_ctrax3],
            (6113, 61, 998): ['CTRAX6', self._read_ctrax6],
            #(16800, 168, 9978): ['CTRIX3FD', self._read_fake],
            #(16900, 169, 9977): ['CTRIX6FD', self._read_fake],
            (3701, 37, 49): ['CTUBE', self._read_ctube], # record 104
            (3901, 39, 50): ['CVISC', self._read_cvisc], # record 105
            #(11701, 117, 559): ['CWELD', self._read_fake],
            #(13501, 135, 564): ['CWELDC', self._read_fake],
            (13701, 137, 565): ['CWELDP', self._read_fake],
            #(13601, 136, 562): ['CWELDG', self._read_fake], # This record is no longer used
            (4301, 43, 28): ['GENEL', self._read_genel],
            #(3201, 32, 478): ['GMBNDC', self._read_fake],
            (12901, 129, 482): ['GMBNDS', self._read_gmbnds],
            (3301, 33, 479): ['GMINTC', self._read_fake],
            (13001, 130, 483): ['GMINTS', self._read_fake],
            #(2801, 28, 630): ['MICPNT', self._read_fake],
            (5201, 52, 11): ['PLOTEL', self._read_plotel],
            #(5202, 52, 669): ['PLOTEL3', self._read_fake],
            #(5203, 52, 670): ['PLOTEL4', self._read_fake],
            #(5204, 52, 671): ['PLOTEL6', self._read_fake],
            #(5205, 52, 672): ['PLOTEL8', self._read_fake],
            #(5206, 52, 673): ['PLOTHEX', self._read_fake],
            #(5208, 52, 675): ['PLOTPEN', self._read_fake],
            #(5209, 52, 676): ['PLOTPYR', self._read_fake],
            #(5207, 52, 674): ['PLOTTET', self._read_fake],
            #(3002, 46, 0): ['Q4AERO', self._read_fake],
            (12801, 128, 417): ['RADBC', self._read_radbc],
            (7801, 78, 8883): ['SINT', self._read_sint],
            (5551, 49, 105): ['SPOINT', self._read_spoint], # record 119
            #(2701, 27, 0): ['T3AERO', self._read_fake],
            #(11601, 116, 9942): ['VUBEAM', self._read_fake],
            #(12301, 123, 145): ['VUHEXA', self._read_fake],
            #(12401, 124, 146): ['VUPENTA', self._read_fake],
            #(11201, 112, 9940): ['VUQUAD4', self._read_fake],
            #(12501, 125, 147): ['VUTETRA', self._read_fake],
            #(11501, 115, 9941): ['VUTRIA3', self._read_fake],
            #(65535, 65535, 65535): ['EODB', self._read_fake],


            #------------------------------------------------------------------------------
            # MSC - DMAP 2016.1
            (2002, 20, 0): ['AEROQ4', self._read_fake], # 2
            (1801, 18, 0): ['AEROT3', self._read_fake], # 3
            (1701, 17, 0): ['BEAMAERO', self._read_fake], # 4
            (2708, 27, 59): ['CAABSF', self._read_caabsf], # 5
            (2108, 21, 224): ['CAXIF2', self._read_fake], # 6
            (2208, 22, 225): ['CAXIF3', self._read_fake], # 7
            (2308, 23, 226): ['CAXIF4', self._read_fake], # 8
            #(2408, 24, 180): ['CBAR', self._read_fake], # record 9

            #(4001, 40, 275): ['CBARAO', self._read_fake],    # record 10
            #(5408, 54, 261): ['CBEAM', self._read_fake],     # record 11
            #(11401, 114, 9016): ['CBEAMP', self._read_fake], # record 12
            #(4601, 46, 298): ['CBEND', self._read_fake],   # record 13
            #(2608, 26, 60): ['CBUSH', self._read_fake],    # record 14
            #(5608, 56, 218): ['CBUSH1D', self._read_fake], # record 15
            (2315, 23, 0): ['CCONE-msc', self._read_fake],  # record 16
            #(201, 2, 69): ['CDAMP1', self._read_cdamp1], # record 17
            #(301, 3, 70): ['CDAMP2', self._read_cdamp2], # record 18
            #(401, 4, 71): ['CDAMP3', self._read_cdamp3], # record 19
            #(501, 5, 72): ['CDAMP4', self._read_cdamp4],      # record 20
            #(10608, 106, 404): ['CDAMP5', self._read_cdamp5], # record 21
            #(6208, 62, 108): ['CDUM2', self._read_fake], # 22
            #(6308, 63, 109): ['CDUM3', self._read_fake], # 23
            #(6408, 64, 110): ['CDUM4', self._read_fake], # 24
            #(6508, 65, 111): ['CDUM5', self._read_fake], # 25
            #(6608, 66, 112): ['CDUM6', self._read_fake], # 26
            #(6708, 67, 113): ['CDUM7', self._read_fake], # 27
            #(6808, 68, 114): ['CDUM8', self._read_cdum8], # 28
            #(6908, 69, 115): ['CDUM9', self._read_cdum9], # 29

            #(601, 6, 73): ['CELAS1', self._read_fake], # record 30
            #(701, 7, 74): ['CELAS2', self._read_fake], # record 31
            #(801, 8, 75): ['CELAS3', self._read_fake], # record 32
            #(901, 9, 76): ['CELAS4', self._read_fake], # record 33
            #(9801, 98, 506): ['CFAST', self._read_fake], # record 34
            (9301, 93, 607): ['CFASTP', self._read_fake],    # 35
            #(8515, 85, 0): ['CFLUID2', self._read_cfluid2],  # record 36 - not done
            #(8615, 86, 0): ['CFLUID3', self._read_cfluid3],  # record 37 - not done
            #(8715, 87, 0): ['CFLUID4', self._read_cfluid4],  # record 38 - not done
            (7701, 77, 8881): ['CINT', self._read_cint],      # record 39 - not done

            #(1908, 19, 104): ['CGAP', self._read_cgap],      # record 40 - buggy
            (8100, 81, 381): ['CHACAB', self._read_chacab],    # 41
            (8200, 82, 383): ['CHACBR', self._read_chacbr],    # 42 - not done
            #(8308, 83, 405): ['CHBDYE', self._read_chbdye],  # record 43
            #(10808, 108, 406): ['CHBDYG', self._read_fake],  # 44
            #(10908, 109, 407): ['CHBDYP', self._read_fake],  # 45
            (7308, 73, 253): ['CHEXA', self._read_chexa],    # record 46
            (16300, 163, 9999): ['CHEXA20F', self._read_fake], # 47
            (14000, 140, 9990): ['CHEXAFD', self._read_chexa], # record 48
            (7908, 79, 369): ['CHEXAL', self._read_fake],     # 49

            (12001, 120, 9011): ['CHEXP', self._read_fake],  # record 50
            (7409, 74, 9991): ['CHEXPR', self._read_chexpr],  # record 51
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
            #(12101, 121, 9012): ['CPENP', self._read_fake],   # 62
            #(4108, 41, 280): ['CPENTA', self._read_cpenta],   # record 63
            #(7509, 75, 9992): ['CPENPR', self._read_fake],    # 64
            (16500, 165, 9999): ['CPENT15F', self._read_fake], # 65
            (16000, 160, 9999): ['CPENT6FD', self._read_fake], # 66
            (17000, 170, 9999): ['CQDX4FD', self._read_fake],  # 67
            (17100, 171, 9999) : ['CQDX9FD', self._read_cquadx],  # record 68
            #(9108, 91, 507): ['CQUAD', self._read_cquad],        # record 69 - not tested

            #(2958, 51, 177): ['CQUAD4', self._read_cquad4],     # record 70
            (13900, 139, 9989): ['CQUAD4FD', self._read_cquad4], # record 71
            (4701, 47, 326): ['CQUAD8', self._read_cquad8], # record 72
            (16400, 164, 9999): ['CQUAD9FD', self._read_fake], # 73
            (11101, 111, 9014): ['CQUADP', self._read_fake], # 74
            (8009, 80, 367): ['CQUADR', self._read_cquadr], # record 75
            (9008, 90, 508): ['CQUADX', self._read_cquadx], # record 76
            (14700, 147, 6662): ['CRBAR', self._read_crbar], # 77
            (17300, 173, 6664): ['CRBE1', self._read_crbe1], # 78
            (17200, 172, 6663): ['CRBE3', self._read_crbe3], # 79

            (11000, 110, 6667): ['CRJOINT', self._read_crjoint],  # 80
            #(3001, 30, 48): ['CROD', self._read_crod],         # record 81
            (12600, 126, 6661): ['CRROD', self._read_crrod],    # 82
            (13801, 138, 570): ['CSEAM', self._read_fake],     # 83
            #(3101, 31, 61): ['CSHEAR', self._read_cshear],     # record 84
            #(4408, 44, 227): ['CSLOT3', self._read_fake],      # 85
            #(4508, 45, 228): ['CSLOT4', self._read_fake],      # 86
            (12201, 122, 9013): ['CTETP', self._read_ctetrap], # record 87
            (5508, 55, 217): ['CTETRA', self._read_ctetra],    # record 88
            (7609, 76, 9993): ['CTETPR', self._read_ctetra],   # record 89

            (16600, 166, 9999): ['CTETR10F', self._read_ctetra], # record 90
            (16100, 161, 9999): ['CTETR4FD', self._read_ctetra], # record 91
            (14801, 148, 643): ['CTQUAD', self._read_fake], # 92
            (14901, 149, 644): ['CTTRIA', self._read_fake], # 93
            (5959, 59, 282): ['CTRIA3', self._read_ctria3], # record 94
            (16200, 162, 9999): ['CTRIA3FD', self._read_fake], # 95
            (4801, 48, 327): ['CTRIA6', self._read_ctria6], # record 96 - buggy
            (16700, 167, 9999): ['CTRIA6FD', self._read_ctria6], # 97
            (11301, 113, 9015): ['CTRIAP', self._read_fake], # 98
            (9200, 92, 385): ['CTRIAR', self._read_ctriar], # record 99

            (6108, 61, 107): ['CTRIAX6', self._read_ctriax6], # 101
            (16800, 168, 9978): ['CTRIX3FD', self._read_fake], # 102
            (16900, 169, 9977): ['CTRIX6FD', self._read_fake], # 103
            #(3701, 37, 49): ['CTUBE', self._read_ctube], # record 104
            #(3901, 39, 50): ['CVISC', self._read_cvisc], # record 105
            (11701, 117, 559): ['CWELD', self._read_fake], # 106; same as cfast
            (13501, 135, 564): ['CWELDC', self._read_fake], # 107
            (13601, 136, 562): ['CWELDG', self._read_fake], # 108
            (14600, 146, 630): ['CWSEAM', self._read_fake], # 109

            #(4301, 43, 28): ['GENEL', self._read_genel],     # 110
            (3201, 32, 478): ['GMBNDC', self._read_gmbndc],   # 111
            #(12901, 129, 482): ['GMBNDS', self.read_gmbnds], # 112
            #(3301, 33, 479): ['GMINTC', self._read_fake],   # 113
            #(13001, 130, 483): ['GMINTS', self._read_fake], # 114
            #(5201, 52, 11): ['PLOTEL', self._read_plotel],  # record 115
            #(12801, 128, 417): ['RADBC', self._read_radbc], # record 116
            (15501, 155, 634): ['RADINT', self._read_fake], # 117
            #(7801, 78, 8883): ['SINT', self._read_fake],    # 118
            #(5551, 49, 105): ['SPOINT', self._read_spoint], # record 119

            (11601, 116, 9942): ['VUBEAM', self._read_vubeam], # record 120
            (12301, 123, 145): ['VUHEXA', self._read_fake],    # 121
            (11201, 112, 9940): ['VUQUAD4', self._read_vuquad4],  # 122
            (12401, 124, 146): ['VUPENTA', self._read_fake],   # 123
            (12501, 125, 147): ['VUTETRA', self._read_fake],   # 124
            (11501, 115, 9941): ['VUTRIA3', self._read_vutria3],  # 125
            (13701, 137, 569): ['WELDP', self._read_fake],     # 126; same as CFASTP

            #----------------------------------------------

            # unorganized
            (6113, 61, 998): ['CTRAX6', self._read_ctrax6],
            (10108, 101, 512) : ['CTRIAX', self._read_ctriax],
            (2108, 21, 224): ['CAXIF2', self._read_fake],
            (5601, 56, 296): ['SESET', self._read_fake],
            (7509, 75, 9992): ['CPENPR', self._read_cpenta],
            (16000, 160, 9988): ['CPENTA6FD', self._read_cpenta],
            (16100, 161, 9986): ['CTETRAFD', self._read_ctetra],
            (16300, 163, 9989): ['CHEXA20F', self._read_chexa],
            (16700, 167, 9981): ['CTRI6FD', self._read_ctria6fd],
            (16800, 168, 9978): ['CTRIAX3FD', self._read_ctrix3fd],  # same as ctria6fd
            (16500, 165, 9987): ['CPENT15F', self._read_cpenta],
            (5008, 50, 258): ['CNGRET', self._read_cngret],
            (12301, 123, 9921): ['ADAPT card', self._read_adapt],
            (12401, 124, 9922): ['FEFACE/PVAL?', self._read_feface_pval],
            (7309, 73, 0): ['CaseControl SET?', self._read_fake],
            (12501, 125, 9923): ['ADAPT card 2', self._read_adapt],    # record
            (3401, 34, 9600): ['GMCONV?', self._read_fake],    # record
            (2901, 29, 9601): ['FEEDGE', self._read_feedge2],  # record
            (16600, 166, 9985) : ['CTETRA?', self._read_ctetra],  # record
            (16200, 162, 9982) : ['CTRIA3', self._read_ctria3fd],  # record
            (16900, 169, 9977) : ['CTRIAX', self._read_ctriax],  # record
            (23500, 235, 6662) : ['', self._read_fake],  # record
            (23800, 238, 6665) : ['', self._read_fake],  # record
            (23900, 239, 6666) : ['', self._read_fake],  # record

            (1976, 1, 1996) : ['', self._read_fake],  # record
            (6120, 1, 60434) : ['', self._read_fake],  # record
            (2024, 1001, 2024) : ['', self._read_fake],  # record
            (801, 1, 572) : ['', self._read_fake],  # record

            (1001, 100, 10000) : ['', self._read_fake],  # record
            (1118, 1, 1874) : ['', self._read_fake],  # record
        }

    def add_op2_element(self, elem):
        """checks that eids are positive and that -1 node ids become None"""
        if elem.eid <= 0:
            self.log.debug(str(elem))
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
        self._add_element_object(elem, allow_overwrites=False)
        #print(str(elem)[:-1])

# 1-AEROQ4 (???)
# AEROT3   (???)
# 1-BEAMAERO (1701,17,0)
    def _read_caabsf(self, data: bytes, n: int) -> int:
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

    def _read_cbar(self, data: bytes, n: int) -> int:
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
        ntotal = 64 * self.factor # 16*4
        fe1 = 28 * self.factor
        fe2 = 32 * self.factor
        nelements = (len(data) - n) // ntotal
        struct_i = self.struct_i if self.size == 4 else self.struct_q
        s1 = Struct(mapfmt(self._endian + b'4i3f3i6f', self.size))
        #s2 = Struct(self._endian + b'4i3f3i6f')
        s2 = s1
        s3 = Struct(mapfmt(self._endian + b'7ii2i6f', self.size))
        for unused_i in range(nelements):
            edata = data[n:n + ntotal]
            # we need this flag before we can figure out how to read f
            # per DMAP: F = FE bit-wise AND with 3
            fe, = struct_i.unpack(edata[fe1:fe2])
            f = fe & 3

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
            assert f == fe, 'f=%s type(f)=%s fe=%s\n%s' % (f, type(f), fe, elem)

            self.add_op2_element(elem)
            n += ntotal
        self.card_count['CBAR'] = nelements
        return n

    def _read_cbarao(self, data: bytes, n: int) -> int:
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
        nelements = (len(data) - n) // 36
        s = Struct(self._endian + b'2i7f')
        for unused_i in range(nelements):
            edata = data[n:n + 36]  # 9*4
            out = s.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  CBARAO=%s\n' % str(out))
            (eid, scale, x1, x2, x3, x4, x5, x6, unused_null) = out
            if scale == 2:
                scale = 'FR'
            else:
                NotImplementedError('CBARAO scale=%r; 2=FR' % scale)
            x = [x1, x2, x3, x4, x5, x6]
            self.add_cbarao(eid, scale, x, comment='')
            n += 36
        self.card_count['CBARAO'] = nelements
        return n

    def _read_cbeam(self, data: bytes, n: int) -> int:
        """
        CBEAM(5408,54,261) - the marker for Record 10
        """
        ntotal = 72 * self.factor  # 18*4
        fe1 = 40 * self.factor
        fe2 = 44 * self.factor
        nelements = (len(data) - n) // ntotal
        struct_i = self.struct_i if self.size == 4 else self.struct_q
        s1 = Struct(mapfmt(self._endian + b'6i3f3i6f', self.size))
        s3 = Struct(mapfmt(self._endian + b'12i6f', self.size))
        for unused_i in range(nelements):
            edata = data[n:n + ntotal]
            fe, = struct_i.unpack(edata[fe1:fe2])
            # per DMAP: F = FE bit-wise AND with 3
            f = fe & 3
            if f == 0:  # basic cid
                out = s1.unpack(edata)
                (eid, pid, ga, gb, sa, sb, x1, x2, x3, fe,
                 pa, pb, w1a, w2a, w3a, w1b, w2b, w3b) = out
                #self.log.info('CBEAM: eid=%s fe=%s f=%s; basic cid' % (eid, fe, f))

                data_in = [[eid, pid, ga, gb, sa, sb, pa, pb, w1a, w2a, w3a, w1b, w2b, w3b],
                           [f, x1, x2, x3]]
            elif f == 1:  # global cid
                # CBEAM    89616   5       384720  384521  0.      0.     -1.
                out = s1.unpack(edata)
                (eid, pid, ga, gb, sa, sb, x1, x2, x3, fe,
                 pa, pb, w1a, w2a, w3a, w1b, w2b, w3b) = out
                #self.log.info('CBEAM: eid=%s fe=%s f=%s; global cid' % (eid, fe, f))
                data_in = [[eid, pid, ga, gb, sa, sb, pa, pb, w1a, w2a, w3a, w1b, w2b, w3b],
                           [f, x1, x2, x3]]
            elif f == 2:  # grid option
                out = s3.unpack(edata)
                (eid, pid, ga, gb, sa, sb, g0, xxa, xxb, fe,
                 pa, pb, w1a, w2a, w3a, w1b, w2b, w3b) = out
                #self.log.info('CBEAM: eid=%s fe=%s f=%s; grid option '
                              #'(g0=%s xxa=%s xxb=%s)' % (eid, fe, f, g0, xxa, xxb))
                if g0 <= 0 or g0 >= 100000000 or xxa != 0 or xxb != 0:
                    # Nastran set this wrong...MasterModelTaxi
                    #CBEAM    621614  2672    900380  900379 .197266 -.978394.0600586
                    #    6


                    f = 1
                    out = s1.unpack(edata)
                    (eid, pid, ga, gb, sa, sb, x1, x2, x3, fe, pa,
                     pb, w1a, w2a, w3a, w1b, w2b, w3b) = out
                    #self.log.info('CBEAM: eid=%s fe=%s f=%s; global cid' % (eid, fe, f))
                    data_in = [[eid, pid, ga, gb, sa, sb, pa, pb, w1a, w2a, w3a, w1b, w2b, w3b],
                               [f, x1, x2, x3]]
                    #self.log.info('   (x1=%s x2=%s x3=%s)' % (x1, x2, x3))
                else:
                    data_in = [[eid, pid, ga, gb, sa, sb, pa, pb, w1a, w2a, w3a, w1b, w2b, w3b],
                               [f, g0]]
            else:
                raise RuntimeError('invalid f value...f=%r' % f)
            if self.is_debug_file:
                self.binary_debug.write('  CBEAM eid=%s f=%s fe=%s %s\n' % (
                    eid, f, fe, str(data_in)))

            elem = CBEAM.add_op2_data(data_in, f)
            self.add_op2_element(elem)
            n += ntotal
        self.card_count['CBEAM'] = nelements
        return n

    def _read_cbeamp(self, data: bytes, n: int) -> int:
        """
        CBEAMP(11401,114,9016) - the marker for Record 11
        """
        self.log.info('skipping CBEAMP in GEOM2')
        return len(data)

    def _read_cbend(self, data: bytes, n: int) -> int:
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
        F = 1 XYZ option - global cooridnate system
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
        ntotal = 52 # 4*13
        nentries = (len(data) - n) // ntotal
        fstruc = Struct(self._endian + b'4i 3f 6i')
        istruc = Struct(self._endian + b'4i 3i 6i')

        for unused_i in range(nentries):
            edata = data[n:n + 52]  # 13*4
            fe, = self.struct_i.unpack(edata[28:32])
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
        self.increase_card_count('CBEND', nentries)
        return n

    def _read_cbush(self, data: bytes, n: int) -> int:
        """
        CBUSH(2608,26,60) - the marker for Record 13
        """
        ntotal = 56 * self.factor
        nelements = (len(data) - n) // ntotal
        struct_obj1 = Struct(mapfmt(self._endian + b'4i iii i ifi3f', self.size))
        struct_obj2 = Struct(mapfmt(self._endian + b'4i fff i ifi3f', self.size))
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
        self.card_count['CBUSH'] = nelements
        return n

    def _read_cbush1d(self, data: bytes, n: int) -> int:
        """
        CBUSH1D(5608,56,218) - the marker for Record 14

        1 EID  I Element identification number
        2 PID  I Property identification number
        3 G(2) I Grid point identification numbers
        5 CID  I Coordinate system identification number
        6 UNDEF(3) none

        """
        ntotal = 32 *  self.factor # 4*8
        nelements = (len(data) - n) // ntotal
        struct_6i = Struct(mapfmt(self._endian + b'8i', self.size))
        for unused_i in range(nelements):
            edata = data[n:n + ntotal]
            out = struct_6i.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  CBUSH1D=%s\n' % str(out))
            (eid, pid, g1, g2, cid, unused_a, unused_b, unused_c) = out
            if cid == -1:
                cid = None
            self.add_cbush1d(eid, pid, [g1, g2], cid=cid)
            n += ntotal
        self.card_count['CBUSH1D'] = nelements
        return n

    def _read_ccone(self, data: bytes, n: int) -> int:
        """
        CCONE(2315,23,0) - the marker for Record 15
        """
        self.log.info('skipping CCONE in GEOM2')
        if self.is_debug_file:
            self.binary_debug.write('skipping CCONE in GEOM2\n')
        return len(data)

    def _read_cdamp1(self, data: bytes, n: int) -> int:
        """
        CDAMP1(201,2,69) - the marker for Record 16
        """
        ntotal = 24 * self.factor # 6*4
        nelements = (len(data) - n) // ntotal
        struct_6i = Struct(mapfmt(self._endian + b'6i', self.size))
        for unused_i in range(nelements):
            edata = data[n:n + ntotal]
            out = struct_6i.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  CDAMP1=%s\n' % str(out))
            #(eid, pid, g1, g2, c1, c2) = out
            elem = CDAMP1.add_op2_data(out)
            self.add_op2_element(elem)
            n += ntotal
        self.card_count['CDAMP1'] = nelements
        return n

    def _read_cdamp2(self, data: bytes, n: int) -> int:
        """
        CDAMP2(301,3,70) - the marker for Record 17
        """
        ntotal = 24 * self.factor # 6*4
        nelements = (len(data) - n) // ntotal
        s = Struct(mapfmt(self._endian + b'if4i', self.size))
        for unused_i in range(nelements):
            edata = data[n:n + ntotal]
            out = s.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  CDAMP2=%s\n' % str(out))
            #(eid, bdamp, g1, g2, c1, c2) = out
            elem = CDAMP2.add_op2_data(out)
            self.add_op2_element(elem)
            n += ntotal
        self.card_count['CDAMP2'] = nelements
        return n

    def _read_cdamp3(self, data: bytes, n: int) -> int:
        """
        CDAMP3(401,4,71) - the marker for Record 18
        """
        struct_4i = Struct(self._endian + b'4i')
        nelements = (len(data) - n) // 16
        for unused_i in range(nelements):
            edata = data[n:n + 16]  # 4*4
            out = struct_4i.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  CDAMP3=%s\n' % str(out))
            #(eid, pid, s1, s2) = out
            elem = CDAMP3.add_op2_data(out)
            self.add_op2_element(elem)
            n += 16
        self.card_count['CDAMP3'] = nelements
        return n

    def _read_cdamp4(self, data: bytes, n: int) -> int:
        """
        CDAMP4(501,5,72) - the marker for Record 19
        """
        s = Struct(self._endian + b'ifii')
        nelements = (len(data) - n) // 16
        for unused_i in range(nelements):
            edata = data[n:n + 16]  # 4*4
            out = s.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  CDAMP4=%s\n' % str(out))
            #(eid, bdamp, s1, s2) = out
            elem = CDAMP4.add_op2_data(out)
            self.add_op2_element(elem)
            n += 16
        self.card_count['CDAMP4'] = nelements
        return n

    def _read_cdamp5(self, data: bytes, n: int) -> int:
        """
        CDAMP5(10608,106,404) - the marker for Record 20
        """
        s = Struct(self._endian + b'4i')
        nelements = (len(data) - n) // 16
        for unused_i in range(nelements):
            edata = data[n:n + 16]  # 4*4
            out = s.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  CDAMP5=%s\n' % str(out))
            #(eid, pid, s1, s2) = out
            elem = CDAMP5.add_op2_data(out)
            self.add_op2_element(elem)
            n += 16
        self.card_count['CDAMP5'] = nelements
        return n

# CDUM2
# CDUM3
# CDUM4
# CDUM5
# CDUM6
# CDUM7

    def _read_cdum8(self, data: bytes, n: int) -> int:
        self.log.info('skipping CDUM9 in GEOM2')
        #ints = np.frombuffer(data[n:], dtype='int32').copy()
        #print('CDUM8', ints)
        return n

    def _read_cdum9(self, data: bytes, n: int) -> int:
        self.log.info('skipping CDUM9 in GEOM2')
        #ints = np.frombuffer(data[n:], dtype='int32').copy()
        #print('CDUM9', ints)
        return n

    def _read_celas1(self, data: bytes, n: int) -> int:
        """
        CELAS1(601,6,73) - the marker for Record 29
        """
        ntotal = 24 * self.factor  # 6*4
        struct_4i = Struct(mapfmt(self._endian + b'6i', self.size))
        nelements = (len(data) - n) // ntotal
        for unused_i in range(nelements):
            edata = data[n:n+ntotal]
            out = struct_4i.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  CELAS1=%s\n' % str(out))
            #(eid, pid, g1, g2, c1, c2) = out
            elem = CELAS1.add_op2_data(out)
            self.add_op2_element(elem)
            n += ntotal
        self.card_count['CELAS1'] = nelements
        return n

    def _read_celas2(self, data: bytes, n: int) -> int:
        """
        CELAS2(701,7,74) - the marker for Record 30
        """
        s1 = Struct(mapfmt(self._endian + b'if4iff', self.size))
        ntotal = 32 * self.factor
        nelements = (len(data) - n) // ntotal
        for unused_i in range(nelements):
            edata = data[n:n+ntotal]
            out = s1.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  CELAS2=%s\n' % str(out))
            #(eid, k, g1, g2, c1, c2, ge, s) = out
            elem = CELAS2.add_op2_data(out)
            self.add_op2_element(elem)
            n += ntotal
        self.card_count['CELAS2'] = nelements
        return n

    def _read_celas3(self, data: bytes, n: int) -> int:
        """
        CELAS3(801,8,75) - the marker for Record 31
        """
        ntotal = 16 * self.factor  # 4*4
        struct_4i = Struct(mapfmt(self._endian + b'4i', self.size))
        nelements = (len(data) - n) // ntotal
        for unused_i in range(nelements):
            edata = data[n:n+ntotal]
            out = struct_4i.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  CELAS3=%s\n' % str(out))
            #(eid, pid, s1, s2) = out
            elem = CELAS3.add_op2_data(out)
            self.add_op2_element(elem)
            n += ntotal
        self.card_count['CELAS3'] = nelements
        return n

    def _read_celas4(self, data: bytes, n: int) -> int:
        """
        CELAS4(901,9,76) - the marker for Record 32
        """
        ntotal = 16 * self.factor  # 4*4
        s = Struct(mapfmt(self._endian + b'ifii', self.size))
        nelements = (len(data) - n) // ntotal
        for unused_i in range(nelements):
            edata = data[n:n + ntotal]
            out = s.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  CELAS4=%s\n' % str(out))
            #(eid, k, s1, s2) = out
            elem = CELAS4.add_op2_data(out)
            self.add_op2_element(elem)
            n += ntotal
        self.card_count['CELAS4'] = nelements
        return n

    def _read_cfast(self, data: bytes, n: int) -> int:
        """
        CFAST(9801,98,506) - the marker for Record ???
        """
        self.log.info('skipping CFAST in GEOM2')
        return len(data)

# CFASTP

    def _read_cfluid2(self, data: bytes, n: int) -> int:
        """
        CFLUID2(8515,85,209) - the marker for Record 35

        1 EID       I Element identification number
        2 IDF1      I RINGFL point 1 identification number
        3 IDF2      I RINGFL point 2 identification number
        4 RHO      RS Mass density
        5 B        RS Bulk modulus
        6 HARMINDX  I Harmonic index
        """
        s = Struct(self._endian + b'3i2fi')
        nelements = (len(data) - n) // 24
        for unused_i in range(nelements):
            edata = data[n:n + 24]  # 6*4
            out = s.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  CFLUID2=%s\n' % str(out))
            eid, idf1, idf2, rho, bi, harmonic = out
            self.add_cfluid2(eid, [idf1, idf2], rho, bi, harmonic)
            n += 24
        self.card_count['CFLUID2'] = nelements
        return n

    def _read_cfluid3(self, data: bytes, n: int) -> int:
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
        s = Struct(self._endian + b'4i2fi')
        nelements = (len(data) - n) // 28
        for unused_i in range(nelements):
            edata = data[n:n + 28]  # 7*4
            out = s.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  CFLUID3=%s\n' % str(out))
            eid, idf1, idf2, idf3, rho, b, harmonic = out
            self.add_cfluid3(eid, [idf1, idf2, idf3], rho, b, harmonic)
            n += 28
        self.card_count['CFLUID3'] = nelements
        return n

    def _read_cfluid4(self, data: bytes, n: int) -> int:
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
        s = Struct(self._endian + b'5i2fi')
        nelements = (len(data) - n) // 32
        for unused_i in range(nelements):
            edata = data[n:n + 32]  # 8*4
            out = s.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  CFLUID4=%s\n' % str(out))
            eid, idf1, idf2, idf3, idf4, rho, bi, harmonic = out
            self.add_cfluid4(eid, [idf1, idf2, idf3, idf4], rho, bi, harmonic)
            n += 32
        self.card_count['CFLUID4'] = nelements
        return n

    def _read_cint(self, data: bytes, n: int) -> int:
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
        self.log.info('skipping CINT in GEOM2')
        # C:\NASA\m4\formats\git\examples\move_tpl\ifcq12p.op2
        # doesn't seem to be a card, more of a general info on the geometry...
        #ints = np.frombuffer(data[n:], dtype=self.idtype).copy()
        return len(data)

    def _read_cgap(self, data: bytes, n: int) -> int:
        """
        CGAP(1908,19,104) - the marker for Record 39
        """
        ntotal = 36 * self.factor  # 9*4
        s1 = Struct(mapfmt(self._endian + b'4i3fii', self.size))
        struct_i = self.struct_i if self.size == 4 else self.struct_q
        f2a = 28 * self.factor
        f2b = 32 * self.factor
        g0a = 16 * self.factor
        g0b = 20 * self.factor

        nelements = (len(data) - n) // ntotal
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
                if self.is_debug_file:
                    self.binary_debug.write('  CGAP eid=%s pid=%s gab0=[%s,%s,%s] x123=[%s,%s,%s] '
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
        self.card_count['CGAP'] = nelements
        return n

    def _read_chacab(self, data: bytes, n: int) -> int:
        """CHACAB(8100,81,381)"""
        return self._run_20nodes(CHACAB, data, n)

    def _read_chacbr(self, data: bytes, n: int) -> int:
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
        nelements = (len(data) - n) // 88
        s = Struct(self._endian + b'22i')
        #if self.is_debug_file:
            #self.binary_debug.write('ndata=%s\n' % (nelements * 44))

        #if self.is_debug_file:
            #self.binary_debug.write(f'  {element.type}=(eid, pid, [n1, n2]')

        for unused_i in range(nelements):
            edata = data[n:n + 88]  # 22*4
            out = s.unpack(edata)
            (eid, pid, *nodes) = out
            nodes = list(nodes)
            if self.is_debug_file:
                self.binary_debug.write(f'  {element.type}=({eid}, {pid}, {nodes}')

            elem = element(eid, pid, nodes)
            self.add_op2_element(elem)
            n += 88
        #if stop:
            #raise RuntimeError('theta is too large...make the quad wrong')
        self.card_count[elem.type] = nelements
        return n

# CHACBR

    def _read_chbdye(self, data: bytes, n: int) -> int:
        """
        CHBDYE(8308,83,405) - the marker for Record ???
        """
        ntotal = 28  # 7*4
        s = Struct(self._endian + b'7i')
        nelements = (len(data) - n) // ntotal
        for unused_i in range(nelements):
            edata = data[n:n+28]
            out = s.unpack(edata)
            (eid, eid2, side, iviewf, iviewb, radmidf, radmidb) = out
            if self.is_debug_file:
                self.binary_debug.write('  CHBDYE=%s\n' % str(out))
            #self.log.debug('  CHBDYE=%s' % str(out))
            data_in = [eid, eid2, side, iviewf, iviewb, radmidf, radmidb]
            elem = CHBDYE.add_op2_data(data_in)
            self._add_thermal_element_object(elem)
            n += ntotal
        self.card_count['CHBDYE'] = nelements
        return n

    def _read_chbdyg(self, data: bytes, n: int) -> int:
        """
        CHBDYG(10808,108,406) - the marker for Record 43
        """
        ntotal = 64  # 16*4
        s = Struct(self._endian + b'16i')
        nelements = (len(data) - n) // ntotal
        for unused_i in range(nelements):
            edata = data[n:n+64]
            out = s.unpack(edata)
            (eid, unused_blank, Type, iviewf, iviewb, radmidf, radmidb, unused_blank2,
             g1, g2, g3, g4, g5, g6, g7, g8) = out
            if self.is_debug_file:
                self.binary_debug.write('  CHBDYG=%s\n' % str(out))
            #self.log.debug('  CHBDYG=%s' % str(out))
            data_in = [eid, Type, iviewf, iviewb, radmidf, radmidb,
                       g1, g2, g3, g4, g5, g6, g7, g8]
            elem = CHBDYG.add_op2_data(data_in)
            self._add_thermal_element_object(elem)
            n += ntotal
        self.card_count['CHBDYG'] = nelements
        return n

    def _read_chbdyp(self, data: bytes, n: int) -> int:
        """
        CHBDYP(10908,109,407)
        """
        ntotal = 60  # 16*4
        s = Struct(self._endian + b'12i 3f')
        nelements = (len(data) - n) // ntotal
        for unused_i in range(nelements):
            edata = data[n:n+60]
            out = s.unpack(edata)
            (eid, pid, Type, iviewf, iviewb, g1, g2, g0, radmidf, radmidb,
             dislin, ce, e1, e2, e3) = out
            if self.is_debug_file:
                self.binary_debug.write('  CHBDYP=%s\n' % str(out))
            #self.log.debug('  CHBDYP=%s' % str(out))
            data_in = [eid, pid, Type, iviewf, iviewb, g1, g2, g0, radmidf, radmidb,
                       dislin, ce, e1, e2, e3]
            elem = CHBDYP.add_op2_data(data_in)
            self._add_thermal_element_object_safe(elem)
            n += ntotal
        self.card_count['CHBDYP'] = nelements
        return n

    def _add_thermal_element_object_safe(self, obj):
        if obj.eid in self.elements:
            self.reject_lines.append(obj.write_card(size=16))
        else:
            self._add_element_object(obj)
        #raise RuntimeError('this should be overwritten by the BDF class')

    def _read_chexa(self, data: bytes, n: int) -> int:
        """
        CHEXA(7308,73,253) - the marker for Record 45
        """
        s = Struct(mapfmt(self._endian + b'22i', self.size))
        ntotal = 88 * self.factor  # 22*4
        nelements = (len(data) - n) // ntotal
        for unused_i in range(nelements):
            edata = data[n:n+ntotal]
            out = s.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  CHEXA=%s\n' % str(out))
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
        self.card_count['CHEXA'] = nelements
        return n

# CHEXA20F
# CHEXAFD
# CHEXAL
# CHEXP
    #def _read_chexp(self, data: bytes, n: int) -> int:
        #"""
        #CHEXP(12001,120,9011) - the marker for Record 50
        #"""
        #self.log.info('skipping CHEXP in GEOM2')
        #if self.is_debug_file:
            #self.binary_debug.write('skipping CHEXP in GEOM2\n')
        #return len(data)

    def _read_chexpr(self, data: bytes, n: int) -> int:
        """
        CHEXPR(7409,74,9991) - the marker for Record 48
        """
        return self._read_chexa(data, n)

    def _read_cmass1(self, data: bytes, n: int) -> int:
        """
        CMASS1(1001,10,65) - the marker for Record 51
        """
        ntotal =  24 * self.factor  # 6*4
        struct_6i = Struct(mapfmt(self._endian + b'6i', self.size))
        nelements = (len(data) - n) // ntotal
        for unused_i in range(nelements):
            edata = data[n:n + ntotal]
            out = struct_6i.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  CMASS1=%s\n' % str(out))
            #(eid, pid, g1, g2, c1, c2) = out
            elem = CMASS1.add_op2_data(out)
            self._add_mass_object(elem)
            n += ntotal
        self.card_count['CMASS1'] = nelements
        return n

    def _read_cmass2(self, data: bytes, n: int) -> int:
        """
        CMASS2(1101,11,66) - the marker for Record 52
        """
        ntotal = 24 * self.factor  # 6*4
        s = Struct(mapfmt(self._endian + b'if4i', self.size))
        nelements = (len(data) - n) // ntotal
        for unused_i in range(nelements):
            edata = data[n:n + ntotal]
            out = s.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  CMASS2=%s\n' % str(out))
            #(eid, m, g1, g2, c1, c2) = out
            elem = CMASS2.add_op2_data(out)
            self._add_mass_object(elem)
            n += ntotal
        self.card_count['CMASS2'] = nelements
        return n

    def _read_cmass3(self, data: bytes, n: int) -> int:
        """
        CMASS3(1201,12,67) - the marker for Record 53
        """
        struct_4i = Struct(self._endian + b'4i')
        nelements = (len(data) - n) // 16
        for unused_i in range(nelements):
            edata = data[n:n + 16]  # 4*4
            out = struct_4i.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  CMASS3=%s\n' % str(out))
            #(eid, pid, s1, s2) = out
            elem = CMASS3.add_op2_data(out)
            self._add_mass_object(elem)
            n += 16
        self.card_count['CMASS3'] = nelements
        return n

    def _read_cmass4(self, data: bytes, n: int) -> int:
        """
        CMASS4(1301,13,68) - the marker for Record 54
        """
        nelements = (len(data) - n) // 16
        struct_if2i = Struct(self._endian + b'ifii')
        for unused_i in range(nelements):
            edata = data[n:n + 16]  # 4*4
            out = struct_if2i.unpack(edata)
            #(eid, m,s 1, s2) = out
            elem = CMASS4.add_op2_data(out)
            self._add_mass_object(elem)
            n += 16
        self.card_count['CMASS4'] = nelements
        return n

    def _read_cmfree(self, data: bytes, n: int) -> int:
        """
        CMFREE(2508,25,0) - the marker for Record 55

        1 EID  I Element identification number
        2   S  I
        3  S2  I
        4   Y RS
        5   N  I
        """
        assert n == 12, n
        nelements = (len(data) - n) // 20
        assert (len(data) - n) % 20 == 0
        struct_3ifi = Struct(self._endian + b'3ifi')
        for unused_i in range(nelements):
            edata = data[n:n + 20]  # 5*4
            out = struct_3ifi.unpack(edata)
            eid, s, s2, y, ncm = out
            self.add_cmfree(eid, s, s2, y, ncm)
            n += 20
        self.card_count['CMFREE'] = nelements
        return n

    def _read_conm1(self, data: bytes, n: int) -> int:
        """
        CONM1(1401,14,63) - the marker for Record 56
        """
        ntotal = 96 * self.factor   # 24*4
        s = Struct(mapfmt(self._endian + b'3i21f', self.size))
        nelements = (len(data) - n) // ntotal
        for unused_i in range(nelements):
            edata = data[n:n + ntotal]
            out = s.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  CONM1=%s\n' % str(out))
            #(eid, g, cid, m1, m2a, m2b, m3a, m3b, m3c, m4a, m4b, m4c, m4d,
             #m5a, m5b, m5c, m5d, m5e, m6a, m6b, m6c, m6d, m6e, m6f) = out
            elem = CONM1.add_op2_data(out)
            self._add_mass_object(elem)
            n += ntotal
        self.card_count['CONM1'] = nelements
        return n

    def _read_conm2(self, data: bytes, n: int) -> int:
        """
        CONM2(1501,15,64) - the marker for Record 57
        """
        ntotal = 52 * self.factor  # 13*4
        s = Struct(mapfmt(self._endian + b'3i10f', self.size))
        nelements = (len(data) - n) // ntotal
        for unused_i in range(nelements):
            edata = data[n:n+ntotal]
            out = s.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  CONM2=%s\n' % str(out))
            #(eid, g, cid, m, x1, x2, x3, i1, i2a, i2b, i3a, i3b, i3c) = out
            elem = CONM2.add_op2_data(out)
            self._add_mass_object(elem)
            n += ntotal
        self.card_count['CONM2'] = nelements
        return n

    def _read_conrod(self, data: bytes, n: int) -> int:
        """
        CONROD(1601,16,47) - the marker for Record 58
        """
        ntotal = 32 *self.factor  # 8*4
        s = Struct(mapfmt(self._endian + b'4i4f', self.size))
        nelements = (len(data) - n) // ntotal
        for unused_i in range(nelements):
            edata = data[n:n+ntotal]
            out = s.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  CONROD=%s\n' % str(out))
            #(eid, n1, n2, mid, a, j, c, nsm) = out
            elem = CONROD.add_op2_data(out)
            self.add_op2_element(elem)
            n += ntotal
        self.card_count['CONROD'] = nelements
        return n

    def _read_conv(self, data: bytes, n: int) -> int:
        """
        The CONV card is different between MSC and NX Nastran.
        The MSC version is 8 fields longer.
        """
        n0 = n
        if self.is_nx:
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
            self._add_thermal_bc_object(elem, elem.eid)
        self.card_count['CONV'] = nelements
        return n

    def _read_split_card(self, data, n, read1, read2, card_name, add_method):
        """
        generalization of multi read methods for different
        versions of MSC Nastran
        """
        n0 = n
        try:
            n, elements = read1(data, n)
        except AssertionError:
            self.log.info('AssertionError...try again reading %r' % card_name)
            n, elements = read2(data, n0)

        nelements = len(elements)
        for elem in elements:
            add_method(elem)
        self.card_count[card_name] = nelements
        return n

    def _read_dual_card(self, data, n, nx_read, msc_read, card_name, add_method):
        """
        generalization of multi read methods (MSC, NX)
        """
        n0 = n
        if self.is_nx:
            try:
                n, elements = nx_read(data, n)
            except (AssertionError, MixedVersionCard):
                #raise
                n, elements = msc_read(data, n0)
        else:
            try:
                n, elements = msc_read(data, n)
            except (AssertionError, MixedVersionCard):
                #raise
                n, elements = nx_read(data, n0)

        nelements = len(elements)
        assert n is not None
        for elem in elements:
            add_method(elem)
        self.card_count[card_name] = nelements
        return n

    def _read_conv_nx(self, data: bytes, n: int) -> int:
        """
        CONV(12701,127,408) - the marker for Record 59
        """
        ntotal = 48  # 12*4
        s = Struct(self._endian + b'4i 8i')
        nelements = (len(data) - n) // ntotal
        elements = []
        for unused_i in range(nelements):
            edata = data[n:n+ntotal]
            out = s.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  CONV=%s; len=%s\n' % (str(out), len(out)))
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
        ntotal = 80  # 20*4
        s = Struct(self._endian + b'12i 8f')
        nelements = (len(data) - n) // ntotal
        elements = []
        for unused_i in range(nelements):
            edata = data[n:n+80]
            out = s.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  CONV=%s; len=%s\n' % (str(out), len(out)))
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

    def _read_convm(self, data: bytes, n: int) -> int:
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
        #C:\Users\sdoyle\Dropbox\move_tpl\ht15330.op2
        ntotal6 = 24  # 6*4
        ntotal7 = 28  # 7*4
        ndata = len(data)
        nelements6 = (ndata - n) // ntotal6
        nelements7 = (ndata - n) // ntotal7

        is_six = (ndata - n) % ntotal6 == 0
        is_seven = (ndata - n) % ntotal7 == 0

        if is_six and is_seven:
            try:
                self.log.warning('CONVM: assuming 6 fields')
                elements, n = self.read_convm6(data, nelements6, n)
                self.log.debug('CONVM: 6 fields')
            except RuntimeError:  # eid < 0
                self.log.warning('CONVM: assuming 7 fields')
                elements, n = self.read_convm7(data, nelements7, n)
                self.log.debug('CONVM: 7 fields')
        elif is_six:
            elements, n = self.read_convm6(data, nelements6, n)
        elif is_seven:
            elements, n = self.read_convm7(data, nelements7, n)
        else:  # pragma: no cover
            raise RuntimeError('CONVM is_six=%s is_seven=%s' % (is_six, is_seven))
        nelements = len(elements)
        for element in elements:
            self._add_thermal_bc_object(element, element.eid)
        self.card_count['CONVM'] = nelements
        return n

    def read_convm6(self, data, nelements, n):
        structi = Struct(self._endian + b'6i')
        elements = []
        for unused_i in range(nelements):
            edata = data[n:n+24]
            out = structi.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  CONVM=%s\n' % str(out))
            (eid, pcon_id, flmnd, cntrlnd, ta1, ta2) = out
            #if eid <= 0:
            if eid <= 0 or pcon_id <= 0 or flmnd < 0 or cntrlnd <= 0 or ta1 <= 0 or ta2 <= 0:
                self.show_data(data, 'if')
                # TODO: I'm not sure that this really has 6 fields...
                raise RuntimeError(f'eid={eid} pconid={pcon_id} flmnd={flmnd} cntrlnd={cntrlnd} ta1={ta1} ta2={ta2} < 0')
            mdot = 0.
            data_in = [eid, pcon_id, flmnd, cntrlnd, ta1, ta2, mdot]
            elem = CONVM.add_op2_data(data_in)
            n += 24
            elements.append(elem)
        return elements, n

    def read_convm7(self, data, nelements, n):
        structi = Struct(self._endian + b'6if')
        elements = []
        for unused_i in range(nelements):
            edata = data[n:n+28]
            out = structi.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  CONVM=%s\n' % str(out))
            (eid, pcon_id, flmnd, cntrlnd, ta1, ta2, mdot) = out

            if eid <= 0 or pcon_id <= 0 or flmnd < 0 or cntrlnd <= 0 or ta1 <= 0 or ta2 <= 0:
                self.show_data(data, 'if')
                # TODO: I'm not sure that this really has 7 fields...
                raise RuntimeError(f'eid={eid} pconid={pcon_id} flmnd={flmnd} '
                                   f'cntrlnd={cntrlnd} ta1={ta1} ta2={ta2} < 0')
            data_in = [eid, pcon_id, flmnd, cntrlnd, ta1, ta2, mdot]
            elem = CONVM.add_op2_data(data_in)
            n += 28
            elements.append(elem)
        return elements, n

    def _read_cpyram(self, data: bytes, n: int) -> int:
        """
        CPYRAM(17200,172,1000) - the marker for Record ???

        Specific to NX Nastran
        """
        self.to_nx()
        ntotal = 64 * self.factor  # 16*4
        struct_16i = Struct(mapfmt(self._endian + b'16i', self.size))
        nelements = (len(data) - n) // ntotal
        for unused_i in range(nelements):
            edata = data[n:n + ntotal]
            out = struct_16i.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  CPENTA=%s\n' % str(out))
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
        self.card_count['CPYRAM'] = nelements
        self.to_nx()
        return n

# CPENP

    def _read_cpenta(self, data: bytes, n: int) -> int:
        """
        CPENTA(4108,41,280)      - the marker for Record 63
        CPENPR(7509,75,9992)     - the marker for Record 64
        CPENT15F(16500,165,9999) - the marker for Record 65
        CPENT6FD(16000,160,9999) - the marker for Record 66
        """
        ntotal = 68 * self.factor  # 17*4
        s = Struct(mapfmt(self._endian + b'17i', self.size))
        nelements = (len(data) - n) // ntotal
        for unused_i in range(nelements):
            edata = data[n:n + ntotal]
            out = s.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  CPENTA=%s\n' % str(out))
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
        self.card_count['CPENTA'] = nelements
        return n

# CQDX4FD
# CQDX9FD - same as CQDX4FD

    def _read_cquad(self, data: bytes, n: int) -> int:
        """CQUAD(9108,91,507)  - the marker for Record 69"""
        return self._run_cquad(CQUAD, data, n)

    def _read_cquad4(self, data: bytes, n: int) -> int:
        """
        CQUAD4(2958,51,177)    - the marker for Record 70
        CQUAD4(13900,139,9989) - the marker for Record 71
        """
        nx_method = partial(self._run_cquad4_nx, CQUAD4)
        msc_method = partial(self._run_cquad4_msc, CQUAD4)
        n = self._read_dual_card(
            data, n,
            nx_method, msc_method,
            'CQUAD4', self.add_op2_element)
        return n

    def _read_vutria3(self, data: bytes, n: int) -> int:
        return self._run_3nodes(CTRIA3, data, n)

    def _read_vuquad4(self, data: bytes, n: int) -> int:
        """VUQUAD4(11201,112,9940)"""
        return self._run_4nodes(CQUAD4, data, n)

    def _run_cquad(self, element: Union[CQUAD, CQUADX], data: bytes, n: int) -> int:
        """common method for CQUAD, CQUADX"""
        ntotal = 44 * self.factor  # 11*4
        s = Struct(mapfmt(self._endian + b'11i', self.size))
        nelements = (len(data) - n) // ntotal
        if self.is_debug_file:
            self.binary_debug.write('ndata=%s\n' % (nelements * ntotal))
        for unused_i in range(nelements):
            edata = data[n:n + ntotal]
            out = s.unpack(edata)
            (eid, pid, n1, n2, n3, n4, n5, n6, n7, n8, n9) = out
            if self.is_debug_file:
                self.binary_debug.write('  %s=%s\n' % (element.type, str(out)))
            #print('CQUAD eid=%s pid=%s n1=%s n2=%s n3=%s n4=%s n5=%s n6=%s n7=%s n8=%s' % (
                #eid, pid, n1, n2, n3, n4, n5, n6, n7, n8))
            datai = [eid, pid, n1, n2, n3, n4, n5, n6, n7, n8, n9]
            elem = element.add_op2_data(datai)
            self.add_op2_element(elem)
            n += ntotal
        self.card_count[element.type] = nelements
        return n

    def _run_2nodes(self, element, data: bytes, n: int) -> int:
        """common method for VUBEAM"""
        nelements = (len(data) - n) // 16
        s = Struct(self._endian + b'4i')
        #if self.is_debug_file:
            #self.binary_debug.write('ndata=%s\n' % (nelements * 44))

        if self.is_debug_file:
            self.binary_debug.write(f'  {element.type}=(eid, pid, [n1, n2]')

        for unused_i in range(nelements):
            edata = data[n:n + 16]  # 4*4
            out = s.unpack(edata)
            (eid, pid, n1, n2) = out
            if self.is_debug_file:
                self.binary_debug.write(
                    f'  {element.type}=({eid}, {pid}, [{n1}, {n2}]')

            nids = [n1, n2]
            elem = element(eid, pid, nids)
            self.add_op2_element(elem)
            n += 16
        #if stop:
            #raise RuntimeError('theta is too large...make the quad wrong')
        self.card_count[elem.type] = nelements
        return n

    def _run_3nodes(self, element: CTRIA3, data: bytes, n: int) -> int:
        """common method for CTRIA3, VUTRIA3"""
        nelements = (len(data) - n) // 20
        s = Struct(self._endian + b'5i')
        #if self.is_debug_file:
            #self.binary_debug.write('ndata=%s\n' % (nelements * 44))

        if self.is_debug_file:
            self.binary_debug.write(f'  {element.type}=(eid, pid, [n1, n2, n3]')

        for unused_i in range(nelements):
            edata = data[n:n + 20]  # 5*4
            out = s.unpack(edata)
            (eid, pid, n1, n2, n3) = out
            if self.is_debug_file:
                self.binary_debug.write(
                    f'  {element.type}=({eid}, {pid}, [{n1}, {n2}, {n3}]')

            nids = [n1, n2, n3]
            elem = element(eid, pid, nids)
            self.add_op2_element(elem)
            n += 20
        #if stop:
            #raise RuntimeError('theta is too large...make the quad wrong')
        self.card_count[element.type] = nelements
        return n

    def _run_4nodes(self, element: Union[CQUAD4, CAABSF], data: bytes, n: int) -> int:
        """common method for CQUAD4, CQUADR"""
        nelements = (len(data) - n) // 24
        s = Struct(self._endian + b'6i')
        #if self.is_debug_file:
            #self.binary_debug.write('ndata=%s\n' % (nelements * 44))

        if self.is_debug_file:
            self.binary_debug.write(f'  {element.type}=(eid, pid, [n1, n2, n3, n4]')

        for unused_i in range(nelements):
            edata = data[n:n + 24]  # 6*4
            out = s.unpack(edata)
            (eid, pid, n1, n2, n3, n4) = out
            if self.is_debug_file:
                self.binary_debug.write(
                    f'  {element.type}=({eid}, {pid}, [{n1}, {n2}, {n3}, {n4}]')

            nids = [n1, n2, n3, n4]
            elem = element(eid, pid, nids)
            self.add_op2_element(elem)
            n += 24
        #if stop:
            #raise RuntimeError('theta is too large...make the quad wrong')
        self.card_count[element.type] = nelements
        return n

    def _run_cquad4_msc(self, element: CQUAD4, data: bytes, n: int) -> Tuple[int, Any]:
        """
        buggy MSC 2018.2 version

        #TODO is this right?  What's with the intermediate zeros?
        data = (
            2958, 51, 177,
            1, 1, 1, 2, 73, 72, 0, 0, 0, 0, -1.0, -1.0, -1.0, -1.0, -1,
            2, 1, 2, 3, 74, 73, 0, 0, 0, 0, -1.0, -1.0, -1.0, -1.0, -1,
            3, 1, 3, 4, 75, 74, 0, 0, 0, 0, -1.0, -1.0, -1.0, -1.0, -1,
        )
        """
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
        s = Struct(mapfmt(self._endian + b'6i 4f 3fi i', self.size))
        #if self.is_debug_file:
            #self.binary_debug.write('ndata=%s\n' % (nelements * 44))

        if self.is_debug_file:
            self.binary_debug.write(f'  {element.type}=(eid, pid, [n1, n2, n3, n4], theta, zoffs, '
                                    'unused_blank, [tflag, t1, t2, t3, t4]); theta_mcid\n')

        for unused_i in range(nelements):
            edata = data[n:n + ntotal]
            out = s.unpack(edata)
            # theta, zoffs, blank, tflag, t1, t2, t3, t4
            (eid, pid, n1, n2, n3, n4,
             theta, zoffs, blank, tflag,
             t1, t2, t3, t4,
             minus1) = out
            minus1 = out[-1]

            msg = ('eid=%s pid=%s nodes=%s '
                   'theta=%s zoffs=%s '
                   'blank=%s tflag=%s '
                   't1-t4=%s minus1=%s' % (
                       eid, pid, (n1, n2, n3, n4),
                       theta, zoffs,
                       blank, tflag,
                       (t1, t2, t3, t4), minus1))

            assert theta == 0, msg
            assert zoffs == 0, msg
            assert blank == 0, msg
            assert tflag == 0, msg

            theta_mcid = convert_theta_to_mcid(theta)
            if self.is_debug_file:
                self.binary_debug.write(
                    f'  {element.type}=({eid}, {pid}, [{n1}, {n2}, {n3}, {n4}], '
                    f'{theta}, {zoffs}, {blank}, [{tflag}, {t1}, {t2}, {t3}, {t4})]; {theta_mcid}\n')
            assert minus1 == -1, minus1

            data_init = [
                eid, pid, n1, n2, n3, n4, theta_mcid, zoffs,
                tflag, t1, t2, t3, t4]
            elem = element.add_op2_data(data_init)
            elements.append(elem)
            n += ntotal
        #if stop:
            #raise RuntimeError('theta is too large...make the quad wrong')
        #self.card_count[element.type] = nelements
        return n, elements

    def _run_cquad4_nx(self, element: Union[CQUAD4, CQUADR],
                       data: bytes, n: int) -> Tuple[int, Any]:
        """
        common method for CQUAD4, CQUADR

        data = (
            2958, 51, 177,
            1, 1, 1, 2, 8, 7, 0, 0, 0, 0, -1.0, -1.0, -1.0, -1.0
            2, 1, 2, 3, 9, 8, 0, 0, 0, 0, -1.0, -1.0, -1.0, -1.0
            3, 1, 3, 4, 10, 9, 0, 0, 0, 0, -1.0, -1.0, -1.0, -1.0
        )
        """
        elements = []
        ntotal = 56 * self.factor  # 14*4
        nelements = (len(data) - n) // ntotal
        leftover = (len(data) - n) % ntotal
        assert leftover == 0, leftover
        s = Struct(mapfmt(self._endian + b'6i ff ii 4f', self.size))
        if self.is_debug_file:
            self.binary_debug.write('ndata=%s\n' % (nelements * 44))

        if self.is_debug_file:
            self.binary_debug.write(f'  {element.type}=(eid, pid, [n1, n2, n3, n4], theta, zoffs, '
                                    'unused_blank, [tflag, t1, t2, t3, t4]); theta_mcid\n')

        for unused_i in range(nelements):
            edata = data[n:n + ntotal]
            out = s.unpack(edata)
            (eid, pid, n1, n2, n3, n4, theta, zoffs, unused_blank, tflag,
             t1, t2, t3, t4) = out
            theta_mcid = convert_theta_to_mcid(theta)
            if self.is_debug_file:
                self.binary_debug.write(
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
            #self.add_op2_element(elem)
            elements.append(elem)
            n += ntotal
        #if stop:
            #raise RuntimeError('theta is too large...make the quad wrong')
        #self.card_count[element.type] = nelements
        return n, elements

# CQUAD4FD

    def _read_cquad8(self, data: bytes, n: int) -> int:
        """common method for reading CQUAD8s"""
        n = self._read_split_card(data, n,
                                  self._read_cquad8_current, self._read_cquad8_v2001,
                                  'CQUAD8', self.add_op2_element)
        return n

    def _read_cquad8_current(self, data: bytes, n: int) -> int:
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
        ntotal = 68 * self.factor # 17*4
        ndatai = len(data) - n
        nelements = ndatai // ntotal
        assert ndatai % ntotal == 0
        s = Struct(mapfmt(self._endian + b'10i 6f i', self.size))
        elements = []
        for unused_i in range(nelements):
            edata = data[n:n + ntotal]
            out = s.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  CQUAD8=%s\n' % str(out))
            #(eid, pid, n1, n2, n3, n4, n5, n6, n7, n8, t1, t2,
             #t3, t4, theta, zoffs, tflag) = out
            tflag = out[-1]
            #self.log.info('cquad8 tflag = %s' % tflag)
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

    def _read_cquad8_v2001(self, data: bytes, n: int) -> int:
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
        ntotal = 64 * self.factor  # 16*4
        ndatai = (len(data) - n)
        nelements = ndatai // ntotal
        assert ndatai % ntotal == 0
        s = Struct(mapfmt(self._endian + b'10i 6f', self.size))
        elements = []
        for unused_i in range(nelements):
            edata = data[n:n + ntotal]
            out = s.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  CQUAD8=%s\n' % str(out))
            (eid, pid, n1, n2, n3, n4, n5, n6, n7, n8, t1, t2,
             t3, t4, theta, zoffs) = out
            tflag = None
            out = (eid, pid, n1, n2, n3, n4, n5, n6, n7, n8, t1, t2,
                   t3, t4, theta, zoffs, tflag)
            #print('eid=%s pid=%s n1=%s n2=%s n3=%s n4=%s theta=%s zoffs=%s '
                  #'tflag=%s t1=%s t2=%s t3=%s t4=%s' % (
                      #eid, pid, n1, n2, n3, n4, theta, zoffs, tflag, t1, t2, t3, t4))
            #data_init = [eid,pid,n1,n2,n3,n4,theta,zoffs,tflag,t1,t2,t3,t4]
            elem = CQUAD8.add_op2_data(out)
            elements.append(elem)
            self.add_op2_element(elem)
            n += ntotal
        return n, elements

# CQUAD9FD
# CQUADP

    def _read_cquadr(self, data: bytes, n: int) -> int:
        """CQUADR(8009,80,367)  - the marker for Record 75"""
        nx_method = partial(self._run_cquad4_nx, CQUADR)
        msc_method = partial(self._run_cquad4_msc, CQUADR)
        n = self._read_dual_card(
            data, n,
            nx_method, msc_method,
            'CQUADR', self.add_op2_element)
        return n

    def _read_cquadx(self, data: bytes, n: int) -> int:
        """CQUADX(9008,90,508)  - the marker for Record 76"""
        return self._run_cquad(CQUADX, data, n)

    def _read_crbar(self, data: bytes, n: int) -> int:
        # C:\NASA\m4\formats\git\examples\move_tpl\nrgd20c.op2
        self.log.info('skipping RBAR in GEOM2')
        return len(data)

    def _read_crbe1(self, data: bytes, n: int) -> int:
        # C:\NASA\m4\formats\git\examples\move_tpl\nrgd406a.op2
        self.log.info('skipping RBE1 in GEOM2')
        return len(data)

    def _read_crbe3(self, data: bytes, n: int) -> int:
        # C:\NASA\m4\formats\git\examples\move_tpl\ngd720a.op2
        self.log.info('skipping RBE3 in GEOM2')
        return len(data)

    def _read_crjoint(self, data: bytes, n: int) -> int:
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
        self.log.info('skipping RJOINT in GEOM2')
        return len(data)

    def _read_crod(self, data: bytes, n: int) -> int:
        """
        CROD(3001,30,48)    - the marker for Record 81
        """
        ntotal = 16 * self.factor # 4*4
        struct_4i = Struct(mapfmt(self._endian + b'4i', self.size))
        nelements = (len(data) - n) // ntotal
        #is_long_ids = False
        for unused_i in range(nelements):
            edata = data[n:n + ntotal]
            out = struct_4i.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  CROD=%s\n' % str(out))
            #(eid, pid, n1, n2) = out
            #if n1 > 100000000 or n2 > 100000000:
                #is_long_ids = True
            elem = CROD.add_op2_data(out)
            self.add_op2_element(elem)
            n += ntotal
        #self._is_long_ids = is_long_ids
        self.card_count['CROD'] = nelements
        return n

    def _read_crrod(self, data: bytes, n: int) -> int:
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
        structi = Struct(self._endian + b'6if')
        nelements = (len(data) - n) // 28  # 7*4
        #is_long_ids = False
        for unused_i in range(nelements):
            edata = data[n:n + 28]
            out = structi.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  RROD=%s\n' % str(out))

            #print(out)
            (eid, n1, n2, unused_lagrange_multiplier_id, cma, cmb, alpha) = out
            assert cma == 0, (eid, n1, n2, unused_lagrange_multiplier_id, cma, cmb, alpha)
            assert cmb == 0, (eid, n1, n2, unused_lagrange_multiplier_id, cma, cmb, alpha)
            #if n1 > 100000000 or n2 > 100000000:
                #is_long_ids = True

            #eid, nids, cma='', cmb='', alpha=0.0
            elem = self.add_rrod(eid, [n1, n2], cma=cma, cmb=cmb, alpha=alpha)
            self.add_op2_element(elem)
            n += 28
        #self._is_long_ids = is_long_ids
        self.card_count['RROD'] = nelements
        return n

# CSEAM

    def _read_cshear(self, data: bytes, n: int) -> int:
        """
        CSHEAR(3101,31,61)    - the marker for Record 84
        """
        ntotal = 24 * self.factor  # 6*4
        struct_6i = Struct(mapfmt(self._endian + b'6i', self.size))
        nelements = (len(data) - n) // ntotal
        for unused_i in range(nelements):
            edata = data[n:n + ntotal]
            out = struct_6i.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  CSHEAR=%s\n' % str(out))
            #(eid, pid, n1, n2, n3, n4) = out
            elem = CSHEAR.add_op2_data(out)
            self.add_op2_element(elem)
            n += ntotal
        self.card_count['CSHEAR'] = nelements
        return n

# CSLOT3
# CSLOT4

    def _read_ctetrap(self, data: bytes, n: int) -> int:
        """
        CTETP(12201,122,9013)    - the marker for Record 87
        .. todo:: needs work
        """
        self.log.info('poor reading of CTETRAP in GEOM2')
        nelements = (len(data) - n) // 108  # 27*4
        struct_27i = Struct(self._endian + b'27i')
        for unused_i in range(nelements):
            edata = data[n:n+108]
            out = struct_27i.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  CTETP=%s\n' % str(out))

            print(out)
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

    def _read_ctetra(self, data: bytes, n: int) -> int:
        """
        CTETRA(5508,55,217)      - the marker for Record 88
        CTETPR(7609,76,9993)     - the marker for Record 89
        CTETR10F(16600,166,9999) - the marker for Record 90
        CTETR4FD(16100,161,9999) - the marker for Record 91
        """
        ntotal = 48 * self.factor  # 12*4
        s = Struct(mapfmt(self._endian + b'12i', self.size))
        nelements = (len(data) - n) // ntotal
        for unused_i in range(nelements):
            edata = data[n:n + ntotal]
            out = s.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  CTETRA=%s\n' % str(out))
            (eid, pid, n1, n2, n3, n4, n5, n6, n7, n8, n9, n10) = out
            #print("out = ",out)

            data_in = [eid, pid, n1, n2, n3, n4]
            big_nodes = [n5, n6, n7, n8, n9, n10]
            if sum(big_nodes) > 0:
                elem = CTETRA10.add_op2_data(data_in + big_nodes)
            else:
                elem = CTETRA4.add_op2_data(data_in)
            self.add_op2_element(elem)
            n += ntotal
        self.card_count['CTETRA'] = nelements
        return n

# CTQUAD - 92
# CTTRIA - 93

    def _read_ctria3(self, data: bytes, n: int) -> int:
        """
        CTRIA3(5959,59,282)    - the marker for Record 94
        """
        ntotal = 52 * self.factor  # 13*4
        s = Struct(mapfmt(self._endian + b'5iff3i3f', self.size))
        nelements = (len(data) - n)// ntotal
        for unused_i in range(nelements):
            edata = data[n:n+ntotal]
            out = s.unpack(edata)
            #print('eid=%s pid=%s n1=%s n2=%s n3=%s theta=%s zoffs=%s '
                  #'blank1=%s blank2=%s tflag=%s t1=%s t2=%s t3=%s' % (
                      #eid, pid, n1, n2, n3, theta, zoffs,
                      #blank1, blank2, tflag, t1, t2, t3))
            (eid, pid, n1, n2, n3, theta, zoffs, unused_blank1,
             unused_blank2, tflag, t1, t2, t3) = out
            if self.is_debug_file:
                self.binary_debug.write('  CTRIA3=%s\n' % str(out))

            theta_mcid = convert_theta_to_mcid(theta)
            data_in = [eid, pid, n1, n2, n3, theta_mcid, zoffs, tflag, t1, t2, t3]
            elem = CTRIA3.add_op2_data(data_in)
            self.add_op2_element(elem)
            n += ntotal
        self.card_count['CTRIA3'] = nelements
        return n


# CTRIAFD - 95

    def _read_ctria6(self, data: bytes, n: int) -> int:
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
        n = self._read_split_card(data, n,
                                  self._read_ctria6_current, self._read_ctria6_v2001,
                                  'CTRIA6', self.add_op2_element)
        return n

    def _read_ctria3fd(self, data: bytes, n: int) -> int:
        """
        Word Name Type Description
        1 EID  I Element identification number
        2 PID  I Property identification number
        3 G(6) I Grid point identification numbers of connection points
        """
        s = Struct(self._endian + b'8i')
        nelements = (len(data) - n) // 32  # 8*4
        assert (len(data) - n) % 32 == 0
        elements = []
        for unused_i in range(nelements):
            edata = data[n:n + 32]
            out = s.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  CTRIA3=%s\n' % str(out))
            (eid, pid, n1, n2, n3, n4, n5, n6) = out
            assert n4 == 0, out
            assert n5 == 0, out
            assert n6 == 0, out
            nids = [n1, n2, n3]
            elem = CTRIA3(eid, pid, nids,
                          theta_mcid=0., zoffset=0., tflag=0,
                          T1=None, T2=None, T3=None, comment='')
            elements.append(elem)
            n += 32

        nelements = len(elements)
        for elem in elements:
            self.add_op2_element(elem)
        self.card_count['CTRIA6'] = nelements
        return n

    def _read_ctrix3fd(self, data: bytes, n: int) -> int:
        """
        Word Name Type Description
        1 EID  I Element identification number
        2 PID  I Property identification number
        3 G(6) I Grid point identification numbers of connection points

        (331, 111, 331, 332, 333, 0, 0, 0,
         332, 11, 331, 333, 334, 0, 0, 0)
        """
        s = Struct(self._endian + b'8i')
        nelements = (len(data) - n) // 32  # 8*4
        assert (len(data) - n) % 32 == 0
        elements = []
        for unused_i in range(nelements):
            edata = data[n:n + 32]
            out = s.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  CTRIA6=%s\n' % str(out))
            (eid, pid, n1, n2, n3, n4, n5, n6) = out
            nids = [n1, n2, n3, n4, n5, n6]
            assert n4 == 0, out
            assert n5 == 0, out
            assert n6 == 0, out
            elem = CTRIAX(eid, pid, nids, theta_mcid=0., comment='')
            elements.append(elem)
            n += 32

        nelements = len(elements)
        for elem in elements:
            self.add_op2_element(elem)
        self.card_count['CTRIAX'] = nelements
        return n

    def _read_ctria6fd(self, data: bytes, n: int) -> int:
        """
        Word Name Type Description
        1 EID  I Element identification number
        2 PID  I Property identification number
        3 G(6) I Grid point identification numbers of connection points
        """
        s = Struct(self._endian + b'8i')
        nelements = (len(data) - n) // 32  # 8*4
        assert (len(data) - n) % 32 == 0
        elements = []
        for unused_i in range(nelements):
            edata = data[n:n + 32]
            out = s.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  CTRIA6=%s\n' % str(out))
            (eid, pid, n1, n2, n3, n4, n5, n6) = out
            nids = [n1, n2, n3, n4, n5, n6]
            #out = (eid, pid, n1, n2, n3, n4, n5, n6, theta, zoffs, t1, t2, t3, 0)
            elem = CTRIA6(eid, pid, nids,
                          theta_mcid=0., zoffset=0., tflag=0,
                          T1=None, T2=None, T3=None, comment='')
            elements.append(elem)
            n += 32

        nelements = len(elements)
        for elem in elements:
            self.add_op2_element(elem)
        self.card_count['CTRIA6'] = nelements
        return n

    def _read_ctria6_current(self, data: bytes, n: int) -> int:
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
        ntotal = 56 * self.factor # 14*4
        s = Struct(mapfmt(self._endian + b'8i 5f i', self.size))
        nelements = (len(data) - n) // ntotal
        assert (len(data) - n) % ntotal == 0
        elements = []
        for unused_i in range(nelements):
            edata = data[n:n + ntotal]
            out = s.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  CTRIA6=%s\n' % str(out))
            #print('eid=%s pid=%s n1=%s n2=%s n3=%s theta=%s zoffs=%s '
                  #'blank1=%s blank2=%s tflag=%s t1=%s t2=%s t3=%s' % (
                      #eid, pid, n1, n2, n3, theta, zoffs,
                      #blank1, blank2, tflag, t1, t2, t3))
            #(eid, pid, n1, n2, n3, n4, n5, n6, theta, zoffs, t1, t2, t3, tflag) = out
            tflag = out[-1]
            #self.log.info('ctria6 tflag = %s' % tflag)
            elem = CTRIA6.add_op2_data(out)
            self.add_op2_element(elem)
            assert tflag in [-1, 0, 1], tflag
            elements.append(elem)
            n += ntotal
        return n, elements

    def _read_ctria6_v2001(self, data: bytes, n: int) -> int:
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
        s = Struct(self._endian + b'8i 5f')
        nelements = (len(data) - n) // 52  # 13*4
        assert (len(data) - n) % 52 == 0
        elements = []
        for unused_i in range(nelements):
            edata = data[n:n + 52]
            out = s.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  CTRIA6=%s\n' % str(out))
            #print('eid=%s pid=%s n1=%s n2=%s n3=%s theta=%s zoffs=%s '
                  #'blank1=%s blank2=%s tflag=%s t1=%s t2=%s t3=%s' % (
                      #eid, pid, n1, n2, n3, theta, zoffs,
                      #blank1, blank2, tflag, t1, t2, t3))
            (eid, pid, n1, n2, n3, n4, n5, n6, theta, zoffs, t1, t2, t3) = out
            out = (eid, pid, n1, n2, n3, n4, n5, n6, theta, zoffs, t1, t2, t3, 0)
            elem = CTRIA6.add_op2_data(out)
            elements.append(elem)
            n += 52
        return n, elements

# CTRIA6FD
# CTRIAP

    def _read_ctriar(self, data: bytes, n: int) -> int:
        """
        CTRIAR(9200,92,385)    - the marker for Record 99
        """
        ntotal = 52 * self.factor  # 13*4
        s = Struct(mapfmt(self._endian + b'5iff3i3f', self.size))
        nelements = (len(data) - n)// ntotal
        for unused_i in range(nelements):
            edata = data[n:n+ntotal]
            out = s.unpack(edata)
            #print('eid=%s pid=%s n1=%s n2=%s n3=%s theta=%s zoffs=%s '
                  #'blank1=%s blank2=%s tflag=%s t1=%s t2=%s t3=%s' % (
                      #eid, pid, n1, n2, n3, theta, zoffs,
                      #blank1, blank2, tflag, t1, t2, t3))
            (eid, pid, n1, n2, n3, theta, zoffs, unused_blank1,
             unused_blank2, tflag, t1, t2, t3) = out
            if self.is_debug_file:
                self.binary_debug.write('  CTRIAR=%s\n' % str(out))
            data_in = [eid, pid, n1, n2, n3, theta, zoffs, tflag, t1, t2, t3]
            elem = CTRIAR.add_op2_data(data_in)
            self.add_op2_element(elem)
            n += ntotal
        self.card_count['CTRIAR'] = nelements
        return n

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
        #struc = Struct(self._endian + b'8i')
        #for unused_i in range(nentries):
            #edata = data[n:n + 32]
            #out = struc.unpack(edata)
            #if self.is_debug_file:
                #self.binary_debug.write('  CTRIAX=%s\n' % str(out))
            #eid, pid, n1, n2, n3, n4, n5, n6 = out
            #nids = [n1, n2, n3, n4, n5, n6]
            #elem = CTRIAX(eid, pid, nids, theta_mcid=0., comment='')
            #self.add_op2_element(elem)
            #n += 32
        #self.card_count['CTRIAX'] = nentries
        #return n

    def _read_ctriax(self, data: bytes, n: int) -> int:
        """common method for reading CTRIAXs"""
        n = self._read_dual_card(data, n, self._read_ctriax_8, self._read_ctriax_9,
                                 'CTRIAX', self.add_op2_element)
        return n

    def _read_ctriax_8(self, data: bytes, n: int) -> Tuple[int, List[CTRIAX]]:
        """(10108, 101, 512)"""
        ntotal = 32  # 9*4
        struc = Struct(self._endian + b'8i')

        nentries = (len(data) - n) // ntotal
        assert (len(data) - n) % ntotal == 0
        assert nentries > 0

        elems = []
        for unused_i in range(nentries):
            edata = data[n:n + ntotal]
            out = struc.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  CTRIAX=%s\n' % str(out))
            eid, pid, n1, n2, n3, n4, n5, n6 = out
            nids = [n1, n2, n3, n4, n5, n6]
            elem = CTRIAX(eid, pid, nids, theta_mcid=0., comment='no theta set')
            elems.append(elem)
            n += ntotal
        return n, elems

    def _read_ctriax_9(self, data: bytes, n: int) -> Tuple[int, List[CTRIAX]]:
        """(10108, 101, 512)"""
        ntotal = 36  # 9*4
        struc = Struct(self._endian + b'9i')

        nentries = (len(data) - n) // ntotal
        assert (len(data) - n) % ntotal == 0
        assert nentries > 0

        elems = []
        for unused_i in range(nentries):
            edata = data[n:n + ntotal]
            out = struc.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  CTRIAX=%s\n' % str(out))
            eid, pid, n1, n2, n3, n4, n5, n6, unused_undef1 = out
            nids = [n1, n2, n3, n4, n5, n6]
            elem = CTRIAX(eid, pid, nids, theta_mcid=0., comment='no theta set')
            elems.append(elem)
            n += ntotal
        return n, elems

    def _read_ctriax6(self, data: bytes, n: int) -> int:  # 101
        """(6108, 61, 107)"""
        ntotal = 44 * self.factor  # 11*4
        nentries = (len(data) - n) // ntotal
        assert (len(data) - n) % ntotal == 0
        assert nentries > 0
        struc = Struct(mapfmt(self._endian + b'8i f ii', self.size))
        for unused_i in range(nentries):
            edata = data[n:n + ntotal]
            out = struc.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  CTRIAX6=%s\n' % str(out))
            elem = CTRIAX6.add_op2_data(out)
            self.add_op2_element(elem)
            n += ntotal
        self.card_count['CTRIAX6'] = nentries
        return n

# CTRIX3FD - 102
# CTRIX6FD - 103

    def _read_gmbndc(self, data: bytes, n: int) -> int:
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
        self.log.info('skipping GMBNDC in GEOM2')
        #self.show_data(data)
        #(1, 31, 32, GRID____, -1,
         #2, 41, 42, GRID____, -1)

        #ints= (3201, 32, 478,
        # 2, 41, 42, 1145390406, 538985799, 41, -1,
        # 990003, 101000045, 101000046, 1145655879, 538976288, 101000045, 101000046, -1)
        ints = np.frombuffer(data[n:], self.idtype) # .tolist()
        isplit = np.where(ints == -1)[0]
        nelements = len(isplit)

        i0 = 0
        for ispliti in isplit:
            eid, gridi, gridf = ints[i0:i0+3]
            #print(eid, gridi, gridf)
            s0 = n + (i0 + 3) * 4
            s1 = s0 + 8
            entity = data[s0:s1].decode('latin1').rstrip()
            eids = ints[i0+5:ispliti]
            assert entity in ['FEEDGE', 'GRID', 'GMCURV'], f'entity={entity!r}'
            #print(eids)
            i0 = ispliti + 1

        self.card_count['GMBNDC'] = nelements
        return len(data)

    def _read_ctube(self, data: bytes, n: int) -> int:
        """
        CTUBE(3701,37,49) - the marker for Record 104
        """
        ntotal = 16 * self.factor  # 4*4
        struct_4i = Struct(mapfmt(self._endian + b'4i', self.size))
        nelements = (len(data) - n) // ntotal
        assert (len(data) - n) % ntotal == 0
        assert nelements > 0
        for unused_i in range(nelements):
            edata = data[n:n + ntotal]
            out = struct_4i.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  CTUBE=%s\n' % str(out))
            #(eid, pid, n1, n2) = out
            elem = CTUBE.add_op2_data(out)
            self.add_op2_element(elem)
            n += ntotal
        self.card_count['CTUBE'] = nelements
        return n

    def _read_cvisc(self, data: bytes, n: int) -> int:
        """CVISC(3901,39,50) - the marker for Record 105"""
        ntotal = 16 * self.factor  # 4*4
        struct_4i = Struct(mapfmt(self._endian + b'4i', self.size))
        nelements = (len(data) - n) // ntotal
        assert (len(data) - n) % ntotal == 0
        assert nelements > 0
        for unused_i in range(nelements):
            edata = data[n:n + ntotal]
            out = struct_4i.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  CVISC=%s\n' % str(out))
            #(eid, pid, n1, n2) = out
            element = CVISC.add_op2_data(out)
            self.add_op2_element(element)
            n += ntotal
        self.card_count['CVISC'] = nelements
        return n

    def _read_cweld(self, data: bytes, n: int) -> int:
        """
        CWELD(11701,117,559) - Record 106
        same as CFAST
        """
        self.log.info('skipping CWELD in GEOM2')
        if self.is_debug_file:
            self.binary_debug.write('skipping CWELD in GEOM2\n')
        return len(data)

    def _read_cweldc(self, data: bytes, n: int) -> int:  # 107
        self.log.info('skipping CWELDC in GEOM2')
        if self.is_debug_file:
            self.binary_debug.write('skipping CWELDC in GEOM2\n')
        return len(data)

    def _read_cweldg(self, data: bytes, n: int) -> int:  # 108
        self.log.info('skipping CWELDG in GEOM2')
        if self.is_debug_file:
            self.binary_debug.write('skipping CWELDG in GEOM2\n')
        return len(data)

# TDOO: above are checked by DMAP...
#-------------------------------
# CWSEAM
    def _read_genel(self, data: bytes, n: int) -> int:
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
        1.401298464324817e-45, 4.203895392974451e-45,
        2.802596928649634e-45, 4.203895392974451e-45,
        4.203895392974451e-45, 4.203895392974451e-45,
        5.605193857299268e-45, 4.203895392974451e-45,
        7.006492321624085e-45, 4.203895392974451e-45,
        8.407790785948902e-45, 4.203895392974451e-45,
        9.80908925027372e-45, 4.203895392974451e-45,
        1.1210387714598537e-44, 4.203895392974451e-45,
        1.2611686178923354e-44, 4.203895392974451e-45,
        1.401298464324817e-44, 4.203895392974451e-45,
        nan,
        1.401298464324817e-44, 1.5414283107572988e-44, 4.203895392974451e-45, 1.5414283107572988e-44, 5.605193857299268e-45, 1.5414283107572988e-44, 7.006492321624085e-45, 1.5414283107572988e-44, 8.407790785948902e-45,
        nan,
        5.605193857299268e-45, 1.401298464324817e-45,
        8.71720021677902e-06, 1.3361000128497835e-06, 1.2778000382240862e-05, 6.272000064200256e-06, 1.6251000488409773e-05, 1.0492000001249835e-05, 2.0478000806178898e-05, 1.562999932502862e-05, 2.428500010864809e-05, 2.0403000235091895e-05,
        3.086099968641065e-05, 6.272000064200256e-06, 3.229700087103993e-05, 1.0492000001249835e-05, 3.352899875608273e-05, 1.562999932502862e-05, 3.502099934848957e-05,
        2.025700086960569e-05, 3.578500036383048e-05, 2.7731999580282718e-05, 1.572600012877956e-05, 4.825499854632653e-05, 3.762800042750314e-05, 7.328399806283414e-05,
        6.433799717342481e-05, 9.580999903846532e-05, 8.837800123728812e-05, 6.374900112859905e-05, 3.762800042750314e-05, 8.013600017875433e-05, 6.433799717342481e-05,
        0.00010011999984271824, 8.837800123728812e-05, 0.00011811000149464235, 0.00012758000229950994, 0.00011344000085955486, 0.00019350000366102904, 0.0001816000003600493,
        0.0002528300101403147, 0.00024294000468216836, 0.0001699900021776557, 0.0001816000003600493, 0.000229199999012053, 0.00024294000468216836, 0.0002824899856932461,
        0.00036862000706605613, 0.00035051998565904796, 0.0005267499946057796, 0.0005117100081406534, 0.00042292001307941973, 0.0005117100081406534, 0.0005718700122088194,
        0.0008483999990858138, 0.0008233999833464622, 0.0009233999880962074, 5.605193857299268e-45, 1.0, 90.0, -20.25, 45.0, 1.0, 90.0, 81.0, 45.0, 1.0, 186.0, -17.850000381469727,
        141.0, 1.0, 186.0, 71.4000015258789, 141.0, 1.0, 268.0, -15.800000190734863, 223.0, 1.0, 268.0, 63.20000076293945, 223.0, 1.0, 368.0, -13.300000190734863, 323.0, 1.0,
        368.0, 53.20000076293945, 323.0, 1.0, 458.0, -11.050000190734863, 413.0, 1.0, 458.0, 44.20000076293945, 413.0)
        """
        self.log.info('skipping GENEL in GEOM2')
        #self.log.info(f'skipping GENEL in GEOM2; len(data)={len(data)-12}')
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

    def _read_plotel(self, data: bytes, n: int) -> int:  # 114
        """(5201, 52, 11)"""
        struct_3i = Struct(self._endian + b'3i')
        ntotal = 12
        nelements = (len(data) - n) // ntotal
        assert (len(data) - n) % ntotal == 0
        assert nelements > 0
        for unused_i in range(nelements):
            edata = data[n:n + ntotal]  # 4*4
            out = struct_3i.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  PLOTEL=%s\n' % str(out))
            #(eid,n1,n2) = out
            elem = PLOTEL.add_op2_data(out)
            self._add_plotel_object(elem)
            n += ntotal
        self.card_count['PLOTEL'] = nelements
        return n

    def _read_radbc(self, data: bytes, n: int) -> int:
        """
        RADBC(12801,128,417)

        Word Name Type Description
        1 EID      I Element identification number
        2 FAMB    RS Radiation view factor between the face and the ambient point
        3 CNTRLND  I Control point for radiation boundary condition
        4 NODAMB   I
        """
        #C:\NASA\m4\formats\git\examples\move_tpl\ht15339.op2
        #(-99, 1.0, 0, 101)
        #radbc   101     1.0             -99
        #RADBC NODAMB   FAMB CNTRLND     EID1 EID2 EID3
        structi = Struct(self._endian + b'ifii')
        ntotal = 16
        nelements = (len(data) - n) // ntotal
        assert (len(data) - n) % ntotal == 0
        assert nelements > 0
        for unused_i in range(nelements):
            edata = data[n:n + ntotal]  # 4*4
            out = structi.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  RADBC=%s\n' % str(out))
            eid, famb, cntrlnd, nodamb = out
            eids = eid
            boundary_condition = RADBC(nodamb, famb, cntrlnd, eids)
            self._add_thermal_bc_object(boundary_condition, boundary_condition.nodamb)
            n += ntotal
        self.card_count['RADBC'] = nelements
        return n

# RADINT
    def _read_sint(self, data: bytes, n: int) -> int:
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
        self.log.info('skipping SINT in GEOM2')
        # C:\NASA\m4\formats\git\examples\move_tpl\ifscp88.op2
        # doesn't seem to be a card, more of a general info on the geometry...
        #ints = np.frombuffer(data[n:], dtype=self.idtype).copy()
        #print(ints.tolist())
        return len(data)

    def _read_spoint(self, data: bytes, n: int) -> int:
        """
        (5551,49,105)    - the marker for Record 118
        """
        ntotal = 4 * self.factor
        npoints = (len(data) - n) // ntotal
        nids = np.frombuffer(data[n:], self.idtype8).tolist()
        if self.is_debug_file:
            self.binary_debug.write('SPOINT=%s\n' % nids)
        spoint = SPOINTs.add_op2_data(nids)
        self._add_spoint_object(spoint)
        self.card_count['SPOINT'] = npoints
        return len(data)

    def _read_vubeam(self, data: bytes, n: int) -> int:  # 119
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

    def _read_ctrax3(self, data: bytes, n: int) -> int:
        """
        RECORD - CTRAX3(6111,61,996)

        Word Name Type Description
        1 EID    I Element identification number
        2 PID    I Property identification number
        3 G(3)   I Grid point identification numbers of connection points
        4 THETA RS Material property orientation angle
        """
        ntotal = 24 * self.factor  # 6*4
        s = Struct(mapfmt(self._endian + b'5if', self.size))
        nelements = (len(data) - n)// ntotal
        for unused_i in range(nelements):
            edata = data[n:n+ntotal]
            out = s.unpack(edata)
            (eid, pid, n1, n2, n3, theta) = out
            if self.is_debug_file:
                self.binary_debug.write('  CTRAX3=%s\n' % str(out))
            #data_in = [eid, pid, n1, n2, n3, theta]
            elem = CTRAX3(eid, pid, [n1, n2, n3], theta)
            self.add_op2_element(elem)
            n += ntotal
        self.card_count['CTRAX3'] = nelements
        return n

    #def _read_cquadx4(self, data: bytes, n: int) -> int:
        #self.log.info('skipping CQUADX4 in GEOM2')
        #return len(data)

    def _read_ctrax6(self, data: bytes, n: int) -> int:
        """
        RECORD - CTRAX6(6113, 61, 998)

        Word Name Type Description
        1 EID    I Element identification number
        2 PID    I Property identification number
        3 G(6)   I Grid point identification numbers of connection points
        4 THETA RS Material property orientation angle
        """
        ntotal = 36 * self.factor  # 9*4
        s = Struct(mapfmt(self._endian + b'8if', self.size))
        nelements = (len(data) - n)// ntotal
        for unused_i in range(nelements):
            edata = data[n:n+ntotal]
            out = s.unpack(edata)
            (eid, pid, n1, n2, n3, n4, n5, n6, theta) = out
            if self.is_debug_file:
                self.binary_debug.write('  CTRAX6=%s\n' % str(out))
            #data_in = [eid, pid, n1, n2, n3, n4, n5, n6, theta]
            elem = CTRAX6(eid, pid, [n1, n2, n3, n4, n5, n6], theta)
            self.add_op2_element(elem)
            n += ntotal
        self.card_count['CTRAX6'] = nelements
        return n

    def _read_cquadx8(self, data: bytes, n: int) -> int:
        """
        RECORD - CQUADX8(6114, 61, 999)

        Word Name Type Description
        1 EID    I Element identification number
        2 PID    I Property identification number
        3 G(8)   I Grid point identification numbers of connection points
        4 THETA RS Material property orientation angle
        """
        ntotal = 44 * self.factor  # 11*4
        s = Struct(mapfmt(self._endian + b'10if', self.size))
        nelements = (len(data) - n)// ntotal
        for unused_i in range(nelements):
            edata = data[n:n+ntotal]
            out = s.unpack(edata)
            (eid, pid, n1, n2, n3, n4, n5, n6, n7, n8, theta) = out
            if self.is_debug_file:
                self.binary_debug.write('  CQUADX8=%s\n' % str(out))
            #data_in = [eid, pid, n1, n2, n3, n4, n5, n6, n7, n8, theta]
            elem = CQUADX8(eid, pid, [n1, n2, n3, n4, n5, n6, n7, n8], theta)
            self.add_op2_element(elem)
            n += ntotal
        self.card_count['CQUADX8'] = nelements
        return n

    def _read_feface_pval(self, data: bytes, n: int) -> int:
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

    def _read_feedge2(self, data: bytes, n: int) -> int:
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
        ntotal = 28  # 7*4
        # s = Struct(self._endian + b'4i 4s 2i') #expected
        s = Struct(self._endian + b'7i')
        nelements = (len(data) - n)// 28  # 7*4
        for unused_i in range(nelements):
            edata = data[n:n+28]
            out = s.unpack(edata)
            #print(out)
            #edge_id, n1, n2, cid, geomin, geom1, geom2 = out # expected
            dunno, nfields, edge_id, n1, n2, n3, n4 = out
            assert nfields in [1, 2, 3, 4, 5], out
            #assert zero1 == 0, f'zero1={zero1} out={out}'
            #assert zero2 == 0, f'zero2={zero2} out={out}'
            if self.is_debug_file:
                self.binary_debug.write('  FEEDGE=%s\n' % str(out))

            geomin_str = 'POINT' # ???
            cid = 0
            geom1 = 0
            geom2 = 0
            if nfields == 2:
                self.add_feedge(edge_id, [n1, n2], cid, [geom1, geom2], geomin=geomin_str)
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
            n += ntotal
        self.card_count['FEEDGE'] = nelements
        return n

    def _read_gmbnds(self, data: bytes, n: int) -> int:
        """
        Word Name Type Description
        1 BID           I Boundary identification number
        2 GRIDC(4)      I Corner grid 1
        6 ENTITY(2) CHAR4 Entity type for defining boundary
        8 EID           I Entity identification numbers for boundary of subdomain
        Word 8 repeats until End of Record
        """
        self.log.info('skipping GMBNDS in GEOM2')
        #(1, 0, 0, 0, 0, 'FEFACE  ', 31, -1)
        ints = np.frombuffer(data[n:], dtype=self.idtype).copy()
        iminus1 = np.where(ints == -1)[0]
        i0 = 0
        for iminus1i in iminus1:
            bid, n1, n2, n3, n4 = ints[i0:i0+5]
            s0 = n + (i0 + 5) * 4
            s1 = s0 + 8
            entity = data[s0:s1].decode('latin1').rstrip()
            assert entity in ['FEFACE', 'GMSURF', 'GRID'], (bid, n1, n2, n3, n4, entity)
            assert bid >= 0, (bid, n1, n2, n3, n4, entity)
            eids = ints[i0+7:iminus1i]
            #print(bid, n1, n2, n3, n4)
            #print('entity = %r' % entity)
            #print(eid)
            #print('-----')
            i0 = iminus1i + 1
        return len(data)

    def _read_cngret(self, data: bytes, n: int) -> int:
        # C:\NASA\m4\formats\git\examples\move_tpl\bpas101.op2
        # C:\NASA\m4\formats\git\examples\move_tpl\pass8.op2
        return len(data)

    def _read_adapt(self, data: bytes, n: int) -> int:
        self.log.info('skipping adapt card in GEOM2')
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
