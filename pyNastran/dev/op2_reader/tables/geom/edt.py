"""
defines readers for BDF objects in the OP2 EDT/EDTS table
"""
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)

#from pyNastran import is_release
from pyNastran.op2.tables.geom.geom_common import GeomCommon


class EDT(GeomCommon):
    """defines methods for reading aero and element deformations"""

    def _read_edt_4(self, data, ndata):
        """
        3.21 EDT
        Aero and element deformations.

        """
        return self._read_geom_4(self._edt_map, data, ndata)

    def __init__(self):
        GeomCommon.__init__(self)

        # F:\Program Files\Siemens\NXNastran\nxn10p1\nxn10p1\nast\tpl\fsw_eng.op2
        # F:\work\pyNastran\pyNastran\master2\pyNastran\bdf\test\nx_spike\out_boltld04i.op2
        # F:\work\pyNastran\pyNastran\master2\pyNastran\bdf\test\nx_spike\out_eliter17.op2
        # F:\work\pyNastran\pyNastran\master2\pyNastran\bdf\test\nx_spike\out_weld01i.op2
        # F:\work\pyNastran\examples\Dropbox\move_tpl\ac10901a_new.op2
        self._edt_map = {
            (5201, 52, 373) : ['ACMODL', self._read_fake],
            (6301, 63, 397) : ['ADAPT', self._read_fake],
            (7801, 78, 582) : ['AECOMP', self._read_fake],
            (7901, 79, 583) : ['AECOMPL', self._read_fake],
            (7301, 73, 574) : ['AEDW', self._read_fake],
            (4002, 40, 273) : ['AEFACT', self._read_fake],
            (7501, 75, 576) : ['AEFORCE', self._read_fake],
            (2602, 26, 386) : ['AELINK', self._read_fake],
            (2302, 23, 341) : ['AELIST', self._read_fake],
            (7001, 70, 571) : ['AEPARM', self._read_fake],
            (7401, 74, 575) : ['AEPRESS', self._read_fake],
            (3202, 32, 265) : ['AERO', self._read_fake],
            (2202, 22, 340) : ['AEROS', self._read_fake],
            (2102, 21, 339) : ['AESTAT', self._read_fake],
            (2002, 20, 338) : ['AESURF', self._read_fake],
            (7701, 77, 581) : ['AESURFS', self._read_fake],
            (3002, 30, 263) : ['CAERO1', self._read_fake],
            (4301, 43, 167) : ['CAERO2', self._read_fake],
            (4401, 44, 168) : ['CAERO3', self._read_fake],
            (4501, 45, 169) : ['CAERO4', self._read_fake],
            (5001, 50, 175) : ['CAERO5', self._read_fake],
            (6201, 62, 143) : ['CLOAD', self._read_fake],
            (6401, 64, 307) : ['CSSCHD', self._read_fake],
            (104, 1, 81) : ['DEFORM', self._read_fake],
            (2702, 27, 387) : ['DIVERG', self._read_fake],
            (4102, 41, 274) : ['FLFACT', self._read_fake],
            (3902, 39, 272) : ['FLUTTER', self._read_fake],
            (17400, 174, 616) : ['GROUP', self._read_fake],
            (3802, 38, 271) : ['MKAERO1', self._read_fake],
            (3702, 37, 270) : ['MKAERO2', self._read_fake],
            (7601, 76, 577) : ['MONPNT1', self._read_fake],
            (3102, 31, 264) : ['PAERO1', self._read_fake],
            (4601, 46, 170) : ['PAERO2', self._read_fake],
            (4701, 47, 171) : ['PAERO3', self._read_fake],
            (4801, 48, 172) : ['PAERO4', self._read_fake],
            (5101, 51, 176) : ['PAERO5', self._read_fake],
            (5301, 53, 378) : ['PANEL', self._read_fake],
            (3502, 35, 268) : ['SET1', self._read_fake],
            (3602, 36, 269) : ['SET2', self._read_fake],
            (4302, 43, 607) : ['SET3', self._read_fake],
            (3302, 33, 266) : ['SPLINE1', self._read_fake],
            (3402, 34, 267) : ['SPLINE2', self._read_fake],
            (4901, 49, 173) : ['SPLINE3', self._read_fake],
            (6501, 65, 308) : ['SPLINE4', self._read_fake],
            (6601, 66, 309) : ['SPLINE5', self._read_fake],
            (2402, 24, 342) : ['TRIM', self._read_fake],
            (7201, 72, 573) : ['UXVEC', self._read_fake],
            (7108, 822, 51) : ['BOLT', self._read_fake],
            (7108, 71, 251) : ['???', self._read_fake],
            (5808, 58, 220) : ['ITER', self._read_fake],
            (14000, 140, 568) : ['SWLDPRM', self._read_fake],
            (11001, 110, 581) : ['???', self._read_fake],
            #(10500, 105, 14) : ['???', self._read_fake],
            #(10500, 105, 14) : ['???', self._read_fake],
            #(10500, 105, 14) : ['???', self._read_fake],
            #(10500, 105, 14) : ['???', self._read_fake],
            #(10500, 105, 14) : ['???', self._read_fake],
            #(10500, 105, 14) : ['???', self._read_fake],
            #(10500, 105, 14) : ['???', self._read_fake],
            #(10500, 105, 14) : ['???', self._read_fake],
            #(10500, 105, 14) : ['???', self._read_fake],
            #(10500, 105, 14) : ['???', self._read_fake],
        }
