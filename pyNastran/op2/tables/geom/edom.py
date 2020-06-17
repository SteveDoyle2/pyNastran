"""
defines readers for BDF objects in the OP2 EDOM/EDOMS table
"""
from struct import Struct
from typing import Union
import numpy as np

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
            (103, 1, 9944) : ['MAT1DOM', self._read_mat1dom],
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
            (3306, 33, 354) : ['DVPREL1', self._read_dvprel1],
            (3406, 34, 355) : ['DVPREL2', self._read_fake],
            #DOPTPRM(4306,43,364)
            (3706, 37, 358) : ['DTABLE', self._read_fake],
            #(3806, 38, 359) : ['DRESP1', self._read_dresp1],
            (3806, 38, 359) : ['DRESP1', self._read_fake],
            (3906, 39, 360) : ['DRESP2', self._read_fake],
            (4106, 41, 362) : ['DCONSTR', self._read_fake],
            (4206, 42, 363) : ['DSCREEN', self._read_fake],
            (4306, 43, 364) : ['DOPTPRM', self._read_fake],
            (4406, 44, 372) : ['DVGRID', self._read_dvgrid],
            #DVSHAP(5006,50,470)
            (5106, 51, 471) : ['DCONADD', self._read_fake],
            #DVBSHAP(5806,58,474)
            #DVGEOM(5906,59,356)
            (6006, 60, 477) : ['MODTRAK', self._read_fake],
            #DRESP3(6700,67,433)
            (6100, 61, 429) : ['DVCREL1', self._read_fake],
            (6200, 62, 430) : ['DVCREL2', self._read_fake],
            (6300, 63, 431) : ['DVMREL1', self._read_dvmrel1],
            (6400, 64, 432) : ['DVMREL2', self._read_fake],
            (6006, 60, 477) : ['???', self._read_fake],
            (7000, 70, 563) : ['DCONSTR/DDVAL?', self._read_fake],
        }

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
        ntotal = 28 * self.factor # 7*4
        struct1 = Struct(mapfmt(self._endian + b'3i 4f', self.size))

        ncards = (len(data) - n) // ntotal
        for unused_i in range(ncards):
            edata = data[n:n + ntotal]
            dvgrid_id, nid, cid, coeff, *dxyz = struct1.unpack(edata)
            assert len(dxyz) == 3, dxyz

            dvgrid = self.add_dvgrid(dvgrid_id, nid, dxyz,
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
        ints = np.frombuffer(data[n:], self.idtype8).copy()
        floats = np.frombuffer(data[n:], self.fdtype8).copy()
        iminus1 = np.where(ints == -1)[0]

        #if self.size == 4:
            #struct1 = Struct(self._endian + b'i 8s i')
            #strs = np.frombuffer(data[n:], dtype='|S4')
        #else:
            #struct1 = Struct(self._endian + b'q 16s q')
            #strs = np.frombuffer(data[n:], dtype='|S8')
        #6i
        #ntotal1 = 16 * self.factor # 4*4

        istart = [0] + list(iminus1[:-1] + 1)
        iend = iminus1
        size = self.size
        for (i0, i1) in zip(istart, iend):
            #self.show_data(data[n+i0*size:n+i1*size], types='ifs')
            assert ints[i1] == -1, ints[i1]
            #print(i0, i1)
            dvprel_id = ints[i0]
            type_bytes = data[n+size:n+3*size]
            property_name_bytes = data[n+8*size:n+10*size]
            if size == 4:
                prop_type = type_bytes.decode('latin1').rstrip()
            else:
                prop_type = reshape_bytes_block(type_bytes).decode('latin1').rstrip()
            pid, fid = ints[i0+3:i0+5]
            pmin, pmax, c0 = floats[i0+5:i0+8]
            if fid == 0:
                fid = None
                if size == 4:
                    fid = property_name_bytes.decode('latin1').rstrip()
                else:
                    fid = reshape_bytes_block(property_name_bytes).decode('latin1').rstrip()

            # fid = fidi
            #print(dvprel_id, prop_type, pid, fid, (pmin, pmax, c0))
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
            dvprel = self.add_dvprel1(dvprel_id, prop_type, pid, fid,
                                      desvar_ids, coeffs,
                                      p_min=pmin, p_max=pmax, c0=c0,
                                      validate=True)
            dvprel.write_card_16()
            n += (i1 - i0 + 1) * size
        return n

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
        ints = np.frombuffer(data[n:], self.idtype8).copy()
        floats = np.frombuffer(data[n:], self.fdtype8).copy()
        iminus1 = np.where(ints == -1)[0]

        istart = [0] + list(iminus1[:-1] + 1)
        iend = iminus1
        size = self.size
        for (i0, i1) in zip(istart, iend):
            #self.show_data(data[n+i0*size:n+i1*size], types='ifs')
            assert ints[i1] == -1, ints[i1]
            dvmrel_id = ints[i0]
            mat_type_bytes = data[n+size:n+3*size]
            mid = ints[i0+3]
            mp_name_bytes = data[n+8*size:n+10*size]
            if size == 4:
                mat_type = mat_type_bytes.decode('latin1').rstrip()
                mp_name = mp_name_bytes.decode('latin1').rstrip()
            else:
                mat_type = reshape_bytes_block(mat_type_bytes).decode('latin1').rstrip()
                mp_name = reshape_bytes_block(mp_name_bytes).decode('latin1').rstrip()
            pid, fid = ints[i0+3:i0+5]
            mp_min, mp_max, c0 = floats[i0+5:i0+8]
            assert fid == 0, (dvmrel_id, mid, mat_type_bytes, mp_name, pid, fid)
            #print(dvmrel_id, mid, mat_type_bytes, mp_name, pid, fid)

            desvar_ids = ints[i0+10:i1:2]
            coeffs = floats[i0+11:i1:2]
            dvmrel = self.add_dvmrel1(dvmrel_id, mat_type, mid, mp_name,
                                      desvar_ids, coeffs,
                                      mp_min=mp_min, mp_max=mp_max, c0=c0,
                                      validate=True)
            dvmrel.write_card_16()
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
        FLAG = 31 PSDACCL
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

        """
        flag_to_resp = {
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
        }

        #self.show_data(data[n:], types='qds')
        ints = np.frombuffer(data[n:], self.idtype8).copy()
        floats = np.frombuffer(data[n:], self.fdtype8).copy()
        iminus1 = np.where(ints == -1)[0]

        #if self.size == 4:
            #struct1 = Struct(self._endian + b'i 8s i')
            #strs = np.frombuffer(data[n:], dtype='|S4')
        #else:
            #struct1 = Struct(self._endian + b'q 16s q')
            #strs = np.frombuffer(data[n:], dtype='|S8')
        #6i
        #ntotal1 = 16 * self.factor # 4*4

        istart = [0] + list(iminus1[:-1] + 1)
        iend = iminus1
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
                assert attb > -1.0e+8, attb
                #print(attbf)
                #ddd
            return attb

        for (i0, i1) in zip(istart, iend):
            assert ints[i1] == -1, ints[i1]
            #print(i0, i1)
            dresp_id = ints[i0]
            label_bytes = data[n+size:n+3*size]
            label = reshape_bytes_block(label_bytes).decode('latin1').rstrip()
            flag = ints[i0+3]
            response_type = flag_to_resp[flag]
            #print(dresp_id, flag, label)
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
                print(response_type, ints[i0+6:i1], floats[i0+6:i1])
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

            elif flag in [6, 7, 11]:
                # FLAG = 6 STRESS
                #  5 PTYPE(2) CHAR4 Element flag (ELEM) or property entry name
                #  7 REGION I Region identifier for constraint screening
                #  8 ATTA   I Response attribute
                #  9 ATTB   I Response attribute
                #  10 ATTi  I Element numbers (if Word 5 is ELEM) or property IDs
                #  Word 10 repeats until -1 occurs
                property_type_bytes = data[n+4*size:n+6*size]
                if self.size == 4:
                    property_type = property_type_bytes.decode('latin1').rstrip()
                else:
                    property_type = reshape_bytes_block(property_type_bytes).decode('latin1').rstrip()
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

            elif flag == 12:
                # FLAG = 12 FREQ
                #   5 UNDEF(2) None
                #   7 REGION I Region identifier for constraint screening
                #   8 ATTA   I Response attribute
                #   9 ATTB   I Response attribute
                #   10 MONE  I Entry is -1
                property_type = None
                region, atta, attb = ints[i0+6:i0+9]
                atti = None
            elif flag == 15:
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
            elif flag == 19:
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

            elif flag == 20:
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
            elif flag == 22:
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

            elif flag in [24, 25]:
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
                property_type_bytes = data[n+4*size:n+6*size]
                if self.size == 4:
                    property_type = property_type_bytes.decode('latin1').rstrip()
                else:
                    property_type = reshape_bytes_block(property_type_bytes).decode('latin1').rstrip()

                region, atta, attbi = ints[i0+6:i0+9]
                attbf = floats[i0+8]
                attb = _pick_attbi_attbf(attbi, attbf)
                atti = ints[i0+9:i1].tolist()
                #print(property_type, region, atta, attb, atti)
            #elif flag == 25:
                # FRDISP
                # 5 PTYPE(2) CHAR4 Element flag (ELEM) or property entry name
                # 7 REGION I Region identifier for constraint screening
                # 8 ATTA   I Response attribute
                # 9 ATTB  RS Frequency value; -1 (integer) spawn for all
                # frequencies in set; -1.10000E+08 for SUM;
                # -1.20000E+08 for AVG; -1.30000E+08 for SSQ;
                # -1.40000E+08 for RSS; -1.50000E+08 for MAX;
                # -1.60000E+08 for MIN
                # 10 ATTi I Element numbers (if Word 5 is ELEM) or property IDs
                # Word 10 repeats until -1 occurs
            elif flag == 29:
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
                print(ints[i0+4:i1+5])
                print(floats[i0+4:i1+5])
                attbf = floats[i0+8]
                attb = _pick_attbi_attbf(attbi, attbf)
                atti = ints[i0+9:i1].tolist()
            else:
                raise NotImplementedError(flag)

            print(response_type)
            if atta == 0:
                atta = None
            if attb == 0:
                attb = None
            if atta is not None:
                atta = int(atta)

            print(dresp_id, label,
                  response_type, property_type, region,
                  atta, attb, atti)
            dresp1 = self.add_dresp1(dresp_id, label,
                                     response_type, property_type, region,
                                     atta, attb, atti, validate=True)
            dresp1.write_card_16()
            n += (i1 - i0 + 1) * self.size
            del dresp_id, label, response_type, property_type, region, atta, attb, atti

        #for i in range(10):
            #edata = data[n:n+ntotal1]
            #dresp1_id, label_bytes, flag = struct1.unpack(edata)
            ##, undef1, undef2, atta, attb, mone
            #label = reshape_bytes_block(label_bytes).decode('latin1').rstrip()
            #print(dresp1_id, label, flag)
            #sss
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
