"""
defines readers for BDF objects in the OP2 EDOM/EDOMS table
"""
from struct import Struct
from pyNastran.op2.tables.geom.geom_common import GeomCommon
#from pyNastran.op2.op2_interface.op2_reader import mapfmt, reshape_bytes_block


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
            (304, 3, 276) : ['DSCONS', self._read_dscons],
            (404, 4, 277) : ['DVAR', self._read_dvar],
            (504, 5, 246) : ['DVSET', self._read_dvset],

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
        """
        #self.show_data(data[12:50], types='ifs')
        structi = Struct(self._endian + b'iii ii ff ii')

        import numpy as np
        ints = np.frombuffer(data[n:], self.idtype8).copy()
        floats = np.frombuffer(data[n:], self.fdtype8).copy()
        iminus1 = np.where(ints == -1)[0]

        istart = [0] + list(iminus1[:-1] + 1)
        iend = iminus1
        for (i0, i1) in zip(istart, iend):
            assert ints[i1] == -1, ints[i1]
            #edata = data[n:n + ntotal]
            #out = structi.unpack(edata)
            #print(out)

            dvset_id, dvset_ptype1, dvset_ptype2, field, flag = ints[i0:i0+5]
            if flag == 1:
                mini, maxi = ints[i0+5:i0+7]
            elif flag == 2:
                mini, maxi = floats[i0+5:i0+7]
            elif flag == 3:
                #print(dvset_id, dvset_ptype1, dvset_ptype2, field, flag)
                print('  ? =', ints[i0+5:i0+7], floats[i0+5:i0+7])
                mini, maxi = '???', '???'
            else:
                print(dvset_id, dvset_ptype1, dvset_ptype2, field, flag)
                raise NotImplementedError(flag)
            pids = ints[i0+7:i1]
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
            print(dvset_id, (dvset_ptype1, dvset_ptype2), ptype, field, flag, (mini, maxi), pids)
        self.log.info(f'skipping {self.card_name} in {self.table_name}; ndata={len(data)-12}')
        return len(data)

    def _read_dvar(self, data: bytes, n: int) -> int:
        """
        DVAR    13013   SPARPNL .01     13013

          data = (404, 4, 277,
          11013, 'SPARPNL ', 0.01, 11013, -1)

        """
        ntotal = 24
        ndatai = len(data) - n
        ncards = ndatai // ntotal
        assert ndatai % ntotal == 0
        structi = Struct(self._endian + b'i 8s fii')
        for unused_i in range(ncards):
            edata = data[n:n + ntotal]
            #self.show_data(edata, types='ifs')
            #(11013, b'SPRCAPS ', 0.01, 11013, -1)
            #(11014, b'SPRCAPS ', 0.01, 11014, -1)
            #(11015, b'SPRCAPS ', 0.01, 11015, -1)
            #(11016, b'SPRCAPS ', 0.01, 11016, -1)
            out = structi.unpack(edata)
            #print(out)

            #idi, word, two_stress_three_force, idb, two_2, value, min_max = out
            #if two_stress_three_force == 2:
                #res_type = 'STRESS'
            #elif two_stress_three_force == 3:
                #res_type = 'FORCE'
            #else:
                #raise NotImplementedError(two_stress_three_force)
            #assert min_max in [0, 1], min_max
            #print(out)
            n += ntotal
        self.log.info(f'skipping {self.card_name} in {self.table_name}; ndata={len(data)-12}')
        return n

    def _read_dscons(self, data: bytes, n: int) -> int:
        """DSCONS

        DSCONS  110131  SPRCAPS STRESS  11013   2       25000.  MAX
        """
        ndatai = len(data) - n
        # !12
        ntotal = 32
        ncards = ndatai // ntotal
        assert ndatai % ntotal == 0
        structi = Struct(self._endian + b'i 8s i 2i fi')
        disp_stress_force_map = {
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
            idi, word, disp_stress_force, idb, two_2, value, min_max = out
            try:
                res_type = disp_stress_force_map[disp_stress_force]
            except KeyError:
                raise NotImplementedError(f'disp_stress_force={disp_stress_force} out={out}')
            assert min_max in [0, 1], min_max
            #print(out)
            n += ntotal
        self.log.info(f'skipping {self.card_name} in {self.table_name}; ndata={len(data)-12}')
        return n

    def _read_desvar(self, data: bytes, n: int) -> int:
        """
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
            s = Struct(self._endian + b'i8s ffff i')
        else:
            ntotal = 64
            s = Struct(self._endian + b'q16s dddd q')

        ncards = (len(data) - n) // ntotal
        for unused_i in range(ncards):
            edata = data[n:n + ntotal]
            desvar_id, blabel, xinit, xlb, xub, delx, ddval = s.unpack(edata)
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
