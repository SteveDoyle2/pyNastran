"""
defines readers for BDF objects in the OP2 EDOM/EDOMS table
"""
from struct import Struct
from pyNastran.op2.tables.geom.geom_common import GeomCommon
from pyNastran.op2.op2_interface.op2_reader import mapfmt, reshape_bytes_block


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
            (103, 1, 9944) : ['MAT1DOM', self._read_fake],
            (304, 3, 276) : ['DSCONS', self._read_fake],
            (404, 4, 277) : ['DVAR', self._read_fake],
            (504, 5, 246) : ['DVSET', self._read_fake],

            (4106, 41, 362) : ['DCONSTR', self._read_fake],
            #DDVAL(7000,70,563)
            #DRESP3(6700,67,433)

            #(504, 5, 246) : ['???', self._read_fake],
            #(504, 5, 246) : ['???', self._read_fake],
            #(504, 5, 246) : ['???', self._read_fake],
            #(504, 5, 246) : ['???', self._read_fake],
            #(504, 5, 246) : ['???', self._read_fake],

            (3106, 31, 352) : ['DESVAR', self._read_desvar],
            (3206, 32, 353) : ['DLINK', self._read_fake],
            (3306, 33, 354) : ['DVPREL1', self._read_fake],
            (3406, 34, 355) : ['DVPREL2', self._read_fake],
            #DOPTPRM(4306,43,364)
            (3706, 37, 358) : ['DTABLE', self._read_fake],
            (3806, 38, 359) : ['DRESP1', self._read_fake],
            (3906, 39, 360) : ['DRESP2', self._read_fake],
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
        if self.size == 4:
            ntotal = 32  # 8*4
            structi = Struct(self._endian + b'i8s ffff i')
        else:
            ntotal = 64
            structi = Struct(self._endian + b'q16s dddd q')

        ncards = (len(data) - n) // ntotal
        for unused_i in range(ncards):
            edata = data[n:n + ntotal]
            desvar_id, blabel, xinit, xlb, xub, delx, ddval = structi.unpack(edata)
            label = blabel.decode('ascii')
            if delx == 0:
                delx = None
            if ddval == 0:
                ddval = None
            if desvar_id not in self.desvars:
                desvar = self.add_desvar(desvar_id, label, xinit, xlb=xlb, xub=xub,
                                         delx=delx, ddval=ddval, comment='')
            else:
                # duplicate DESVAR
                desvar_temp = self.add_desvar(1.0, label, xinit, xlb=xlb, xub=xub,
                                              delx=delx, ddval=ddval, comment='')
                del self.desvars[1.0]
                desvar_temp.desvar_id = desvar_id
                assert desvar_temp == self.desvars[desvar_id]
            n += ntotal
            #print(desvar)
        self.card_count['DESVAR'] = ncards
        return n
