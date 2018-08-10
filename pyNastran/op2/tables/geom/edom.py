"""
defines readers for BDF objects in the OP2 EDOM/EDOMS table
"""
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)

#from pyNastran import is_release
from pyNastran.op2.tables.geom.geom_common import GeomCommon


class EDOM(GeomCommon):
    """defines methods for reading op2 properties"""

    def _read_edom4_4(self, data, ndata):
        """
        reads the EDOM table
        SOL 200 design optimization and sensitivity analysis bulk entries.

        """
        return self._read_geom_4(self._edom_map, data, ndata)

    def __init__(self):
        GeomCommon.__init__(self)

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
        self._edom_map = {
            # are these 3 really EDOM?
            #MAT1DOM(103,1,9944)
            #MAT10DOM(2801,28,9945)
            #MODTRAK(6006,60,477)
            (103, 1, 9944) : ['???', self._read_fake],
            (304, 3, 276) : ['???', self._read_fake],
            (404, 4, 277) : ['???', self._read_fake],
            (504, 5, 246) : ['???', self._read_fake],

            (4106, 41, 362) : ['DCONSTR', self._read_fake],
            #DDVAL(7000,70,563)
            #DRESP3(6700,67,433)

            (504, 5, 246) : ['???', self._read_fake],
            (504, 5, 246) : ['???', self._read_fake],
            (504, 5, 246) : ['???', self._read_fake],
            (504, 5, 246) : ['???', self._read_fake],
            (504, 5, 246) : ['???', self._read_fake],

            (3106, 31, 352) : ['DESVAR', self._read_fake],
            (3206, 32, 353) : ['DLINK', self._read_fake],
            (3306, 33, 354) : ['DVPREL1', self._read_fake],
            (3406, 34, 355) : ['DVPREL2', self._read_fake],
            #DOPTPRM(4306,43,364)
            (3706, 37, 358) : ['DTABLE', self._read_fake],
            (3806, 38, 359) : ['DRESP1', self._read_fake],
            (3906, 39, 360) : ['DRESP2', self._read_fake],
            (4106, 41, 362) : ['DCONSTR', self._read_fake],
            (4206, 42, 363) : ['DSCREEN', self._read_fake],
            (4306, 43, 364) : ['DOPTPRM', self._read_fake],
            (4406, 44, 372) : ['DVGRID', self._read_fake],
            #DVSHAP(5006,50,470)
            (5106, 51, 471) : ['DCONADD', self._read_fake],
            #DVBSHAP(5806,58,474)
            #DVGEOM(5906,59,356)
            (6006, 60, 477) : ['MODTRAK', self._read_fake],
            #DRESP3(6700,67,433)
            (6100, 61, 429) : ['DVCREL1', self._read_fake],
            (6200, 62, 430) : ['DVCREL2', self._read_fake],
            (6300, 63, 431) : ['DVMREL1', self._read_fake],
            (6400, 64, 432) : ['DVMREL2', self._read_fake],
            (6006, 60, 477) : ['???', self._read_fake],
            (7000, 70, 563) : ['DCONSTR/DDVAL?', self._read_fake],
        }
