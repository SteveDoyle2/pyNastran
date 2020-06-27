"""
defines readers for BDF objects in the OP2 EPT/EPTS table
"""
#pylint: disable=C0103,R0914
from struct import unpack, Struct
from typing import Tuple, List

import numpy as np

#from pyNastran import is_release
from pyNastran.bdf.cards.properties.mass import PMASS, NSM, NSML
from pyNastran.bdf.cards.properties.bars import PBAR, PBARL, PBEND
from pyNastran.bdf.cards.properties.beam import PBEAM, PBEAML, PBCOMP
from pyNastran.bdf.cards.properties.bush import PBUSH, PBUSHT
from pyNastran.bdf.cards.properties.damper import PDAMP, PVISC
from pyNastran.bdf.cards.properties.properties import PFAST, PGAP
from pyNastran.bdf.cards.properties.rods import PROD, PTUBE
from pyNastran.bdf.cards.properties.shell import PSHEAR, PSHELL, PCOMP
from pyNastran.bdf.cards.properties.solid import PSOLID
from pyNastran.bdf.cards.properties.springs import PELAS, PELAST

from pyNastran.bdf.cards.thermal.thermal import PCONV, PHBDY, PCONVM
# PCOMPG, PBUSH1D, PBEAML, PBEAM3
from pyNastran.op2.tables.geom.geom_common import GeomCommon
from pyNastran.op2.op2_interface.op2_reader import (
    mapfmt, reshape_bytes_block, reshape_bytes_block_size)


class EPT(GeomCommon):
    """defines methods for reading op2 properties"""

    def _read_ept_4(self, data, ndata):
        return self._read_geom_4(self._ept_map, data, ndata)

    def __init__(self):
        GeomCommon.__init__(self)
        self.big_properties = {}
        self._ept_map = {
            (3201, 32, 55): ['NSM', self._read_nsm],          # record 2
            (52, 20, 181): ['PBAR', self._read_pbar],         # record 11 - buggy
            (9102, 91, 52): ['PBARL', self._read_pbarl],      # record 12 - almost there...
            (2706, 27, 287): ['PCOMP', self._read_pcomp],     # record 22 - buggy
            (302, 3, 46): ['PELAS', self._read_pelas],        # record 39
            (2102, 21, 121): ['PGAP', self._read_pgap],       # record 42
            (902, 9, 29): ['PROD', self._read_prod],          # record 49
            (1002, 10, 42): ['PSHEAR', self._read_pshear],    # record 50
            (2402, 24, 281): ['PSOLID', self._read_psolid],   # record 51
            (2302, 23, 283): ['PSHELL', self._read_pshell],   # record 52
            (1602, 16, 30): ['PTUBE', self._read_ptube],      # record 56

            (5402, 54, 262): ['PBEAM', self._read_pbeam],      # record 14 - not done
            (9202, 92, 53): ['PBEAML', self._read_pbeaml],     # record 15
            (2502, 25, 248): ['PBEND', self._read_pbend],      # record 16 - not done
            (1402, 14, 37): ['PBUSH', self._read_pbush],       # record 19 - not done
            (3101, 31, 219): ['PBUSH1D', self._read_pbush1d],  # record 20 - not done
            (152, 19, 147): ['PCONEAX', self._read_pconeax],   # record 24 - not done
            (11001, 110, 411): ['PCONV', self._read_pconv],    # record 25 - not done
            # record 26
            (202, 2, 45): ['PDAMP', self._read_pdamp],      # record 27 - not done
            (2802, 28, 236): ['PHBDY', self._read_phbdy],   # record 43 - not done
            (402, 4, 44): ['PMASS', self._read_pmass],      # record 48
            (1802, 18, 31): ['PVISC', self._read_pvisc],    # record 59
            (10201, 102, 400): ['PVAL', self._read_pval],   # record 58 - not done
            (2606, 26, 289): ['VIEW', self._read_view],     # record 62 - not done
            (3201, 32, 991) : ['NSM', self._read_fake],  # record
            (3301, 33, 992) : ['NSM1', self._read_fake],  # record
            (3701, 37, 995) : ['NSML1', self._read_fake],    # record
            (15006, 150, 604): ['PCOMPG', self._read_pcompg],  # record

            (702, 7, 38): ['PBUSHT', self._read_pbusht],  # record 1
            (3301, 33, 56): ['NSM1', self._read_fake],  # record 3
            (3401, 34, 57) : ['NSMADD', self._read_fake],    # record 5
            (3501, 35, 58): ['NSML', self._read_fake],  # record 6
            (3501, 35, 994) : ['NSML', self._read_fake],
            (3601, 36, 62): ['NSML1', self._read_fake],  # record 7
            (1502, 15, 36): ['PAABSF', self._read_fake],  # record 8
            (8300, 83, 382): ['PACABS', self._read_fake],  # record 9
            (8500, 85, 384): ['PACBAR', self._read_fake],  # record 10
            (5403, 55, 349): ['PBCOMP', self._read_pbcomp],  # record 13
            (13301, 133, 509): ['PBMSECT', self._read_fake],  # record 17
            (2902, 29, 420): ['PCONVM', self._read_pconvm],  # record 26
            (1202, 12, 33): ['PDAMPT', self._read_pdampt],  # record 28
            (8702, 87, 412): ['PDAMP5', self._read_pdamp5],  # record 29
            (6802, 68, 164): ['PDUM8', self._read_fake],  # record 37
            (6902, 69, 165): ['PDUM9', self._read_fake],  # record 38
            (1302, 13, 34): ['PELAST', self._read_pelast],  # record 41
            (12001, 120, 480): ['PINTC', self._read_fake],  # record 44
            (12101, 121, 484): ['PINTS', self._read_fake],  # record 45
            (4606, 46, 375): ['PLPLANE', self._read_plplane],  # record 46
            (4706, 47, 376): ['PLSOLID', self._read_plsolid],  # record 47
            (10301, 103, 399): ['PSET', self._read_pset],  # record 57
            (3002, 30, 415): ['VIEW3D', self._read_fake],  # record 63

            (13501, 135, 510) : ['PFAST', self._read_pfast_msc],  # MSC-specific

            # NX-specific
            (3601, 36, 55) : ['PFAST', self._read_pfast_nx],  # NX-specific
            (3801, 38, 979) : ['PPLANE', self._read_fake],
            (11801, 118, 560) : ['PWELD', self._read_fake],
            (3401, 34, 993) : ['NSMADD', self._read_fake],
        }

    def _add_op2_property(self, prop):
        """helper method for op2"""
        #if prop.pid > 100000000:
            #raise RuntimeError('bad parsing; pid > 100000000...%s' % str(prop))
        #print(str(prop)[:-1])
        ntables = self.table_names.count(b'EPT') + self.table_names.count(b'EPTS')
        pid = prop.pid
        allow_overwrites = (
            ntables > 1 and
            pid in self.properties and
            self.properties[pid].type == prop.type)
        self._add_property_object(prop, allow_overwrites=allow_overwrites)

    def _add_op2_property_mass(self, prop):
        """helper method for op2"""
        #if prop.pid > 100000000:
            #raise RuntimeError('bad parsing; pid > 100000000...%s' % str(prop))
        #print(str(prop)[:-1])
        ntables = self.table_names.count(b'EPT') + self.table_names.count(b'EPTS')
        pid = prop.pid
        allow_overwrites = (
            ntables > 1 and
            pid in self.properties_mass and
            self.properties_mass[pid].type == prop.type)
        self._add_property_mass_object(prop, allow_overwrites=allow_overwrites)

    def _add_pconv(self, prop):
        if prop.pconid > 100000000:
            raise RuntimeError('bad parsing pconid > 100000000...%s' % str(prop))
        self._add_convection_property_object(prop)

# HGSUPPR

    def _read_nsm(self, data: bytes, n: int) -> int:
        """NSM"""
        n = self._read_dual_card(data, n, self._read_nsm_nx, self._read_nsm_msc,
                                 'NSM', self._add_nsm_object)
        return n

    def _read_nsm_msc(self, data: bytes, n: int) -> int:
        """
        NSM(3201,32,55) - the marker for Record 2

        MSC
        1 SID       I Set identification number
        2 PROP  CHAR4 Set of property or elements
        3 ID        I Property or element identification number
        4 VALUE    RS Nonstructural mass value
        ORIGIN=0 NSM Bulk Data entry
          5 ID I Property or element ID
          6 VALUE RS Nonstructural mass value
          Words 5 through 6 repeat until End of Record
        ORIGIN=2 NSML Bulk Data entry
          5 ID I Property or element ID
          6 VALUE RS Nonstructural mass value
          Words 5 through 6 repeat until End of Record
        Words 3 through 4 repeat until End of Record
        """
        properties = []
        struct1 = Struct(self._endian + b'i 4s if')
        ndelta = 16

        i = 0
        ints = np.frombuffer(data[n:], self.idtype).copy()
        floats = np.frombuffer(data[n:], self.fdtype).copy()

        while n < len(data):
            edata = data[n:n+ndelta]
            out = struct1.unpack(edata)
            (sid, prop_set, pid, value) = out
            #            538976312
            assert pid < 100000000
            i += 4
            n += ndelta

            prop_set = prop_set.decode('utf8').rstrip(' ') # \x00
            values = [value]
            #print('ints[i:]=', ints[i:])
            while ints[i] != -1:
                value2 = floats[i]
                values.append(value2)
                n += 4
                i += 1
            self.log.info("MSC: NSM-sid=%s prop_set=%s pid=%s values=%s" % (
                sid, prop_set, pid, values))
            prop = NSM.add_op2_data([sid, prop_set, pid, value])
            #self._add_nsm_object(prop)
            properties.append(prop)

            # handle the trailing -1
            i += 1
            n += 4
        return n, properties

    def _read_nsm_nx(self, data: bytes, n: int) -> int:
        """
        NSM(3201,32,55) - the marker for Record 2

        1 SID         I Set identification number
        2 PROP(2) CHAR4 Set of properties or elements
        4 ORIGIN      I  Entry origin
        5 ID          I  Property or element identification number
        6 VALUE      RS Nonstructural mass value
        Words 5 through 6 repeat until End of Record
        """
        properties = []

        #NX: C:\Users\sdoyle\Dropbox\move_tpl\nsmlcr2s.op2
        struct1 = Struct(self._endian + b'i 8s ii f')
        ndelta = 24
        #self.show_data(data[12:], 'ifs')

        i = 0
        ints = np.frombuffer(data[n:], self.idtype).copy()
        floats = np.frombuffer(data[n:], self.fdtype).copy()

        unused_packs = break_by_minus1(ints)
        #for pack in packs:
            #print(pack)

        #ipack = 0
        while n < len(data):
            #print('ints[i:]=', ints[i:].tolist())
            #i1, i2 = packs[ipack]
            #print('idata=%s' % idata[i1:i2])
            #print('fdata=%s' % fdata[i1:i2])
            #print(idata[i1:i2])
            edata = data[n:n+ndelta]
            out = struct1.unpack(edata)
            (sid, prop_set, origin, pid, value) = out
            #            538976312
            assert pid < 100000000
            i += 6
            n += ndelta

            prop_set = prop_set.decode('utf8').rstrip(' ') # \x00
            pids = [pid]
            values = [value]
            #print('ints[i:]=', ints[i:].tolist())
            while ints[i] != -1:
                pid = ints[i]
                value2 = floats[i+1]
                assert pid != -1
                pids.append(pid)
                values.append(value2)
                n += 8
                i += 2

            for pid, value in zip(pids, values):
                if origin == 0:
                    #self.log.info("NX: NSM-sid=%s prop_set=%s pid=%s values=%s" % (
                        #sid, prop_set, pid, values))
                    prop = NSM.add_op2_data([sid, prop_set, pid, value])
                elif origin == 2:
                    #self.log.info("NX: NSML-sid=%s prop_set=%s pid=%s values=%s" % (
                        #sid, prop_set, pid, values))
                    prop = NSML.add_op2_data([sid, prop_set, pid, value])

                #print(prop.rstrip(), pid, value)
                #self._add_nsm_object(prop)
                properties.append(prop)
            #print('----')

            # handle the trailing -1
            i += 1
            n += 4
            #ipack += 1
        return n, properties

# NSM1
# NSML1
# NSMADD
# NSML
# NSML1
# PAABSF
# PACABS
# PACBAR

    def _read_pbar(self, data: bytes, n: int) -> int:
        """
        PBAR(52,20,181) - the marker for Record 11
        .. warning:: this makes a funny property...

        MSC 2016/NX10

        Word Name Type Description
        1  PID  I Property identification number
        2  MID  I Material identification number
        3  A   RS Area
        4  I1  RS Area moment of inertia in plane 1
        5  I2  RS Area moment of inertia in plane 2
        6  J   RS Torsional constant
        7  NSM RS Nonstructural mass per unit length
        8  FE  RS
        9  C1  RS Stress recovery location at point C in element y-axis
        10 C2  RS Stress recovery location at point C in element z-axis
        11 D1  RS Stress recovery location at point D in element y-axis
        12 D2  RS Stress recovery location at point D in element z-axis
        13 E1  RS Stress recovery location at point E in element y-axis
        14 E2  RS Stress recovery location at point E in element z-axis
        15 F1  RS Stress recovery location at point F in element y-axis
        16 F2  RS Stress recovery location at point F in element z-axis
        17 K1  RS Area factor for shear in plane 1
        18 K2  RS Area factor for shear in plane 2
        19 I12 RS Area product of inertia for plane 1 and 2
        """
        ntotal = 76 * self.factor  # 19*4
        struct1 = Struct(mapfmt(self._endian + b'2i17f', self.size))
        nentries = (len(data) - n) // ntotal
        for unused_i in range(nentries):
            edata = data[n:n+ntotal]
            out = struct1.unpack(edata)
            #(pid, mid, a, I1, I2, J, nsm, fe, c1, c2, d1, d2,
             #e1, e2, f1, f2, k1, k2, I12) = out
            prop = PBAR.add_op2_data(out)
            self._add_op2_property(prop)
            n += ntotal
        self.card_count['PBAR'] = nentries
        return n

    def _read_pbarl(self, data: bytes, n: int) -> int:
        """
        PBARL(9102,91,52) - the marker for Record 12
        TODO: buggy
        It's possible to have a PBARL and a PBAR at the same time.
        NSM is at the end of the element.
        """
        valid_types = {
            "ROD": 1,
            "TUBE": 2,
            "I": 6,
            "CHAN": 4,
            "T": 4,
            "BOX": 4,
            "BAR": 2,
            "CROSS": 4,
            "H": 4,
            "T1": 4,
            "I1": 4,
            "CHAN1": 4,
            "Z": 4,
            "CHAN2": 4,
            "T2": 4,
            "BOX1": 6,
            "HEXA": 3,
            "HAT": 4,
            "HAT1": 5,
            "DBOX": 10,  # was 12
            #'MLO TUBE' : 2,
        }  # for GROUP="MSCBML0"

        size = self.size
        ntotal = 28 * self.factor  # 7*4 - ROD - shortest entry...could be buggy... # TODO fix this
        if size == 4:
            struct1 = Struct(self._endian + b'2i 8s 8s f')
        else:
            struct1 = Struct(self._endian + b'2q 16s 16s d')

        #nentries = (len(data) - n) // ntotal
        #print(self.show_ndata(80))
        ndata = len(data)

        while ndata - n > ntotal:
            edata = data[n:n+ntotal]
            n += ntotal

            out = struct1.unpack(edata)
            (pid, mid, group, beam_type, value) = out
            if pid > 100000000 or pid < 1:
                self.log.debug("  pid=%s mid=%s group=%r beam_type=%r value=%s" % (
                    pid, mid, group, beam_type, value))
                raise RuntimeError('bad parsing...')

            beam_type = reshape_bytes_block_size(beam_type, size=size)
            group = reshape_bytes_block_size(group, size=size)
            data_in = [pid, mid, group, beam_type, value]

            expected_length = valid_types[beam_type]
            iformat = self._endian + b'%if' % expected_length

            ndelta = expected_length * 4
            dims_nsm = list(unpack(iformat, data[n:n+ndelta]))
            data_in += dims_nsm
            #print("  pid=%s mid=%s group=%r beam_type=%r value=%s dims_nsm=%s" % (
                #pid, mid, group, beam_type, value, dims_nsm))

            # TODO why do i need the +4???
            #      is that for the nsm?
            #min_len =  expected_length * 4 + 4
            #if len(data)
            #data = data[n + expected_length * 4 + 4:]
            n += ndelta

            #prin( "len(out) = ",len(out)))
            #print("PBARL = %s" % data_in)
            prop = PBARL.add_op2_data(data_in)  # last value is nsm
            pid = prop.pid
            if pid in self.properties:
                #self.log.debug("removing:\n%s" % self.properties[pid])
                self._type_to_id_map['PBAR'].remove(pid)
                del self.properties[pid]
            self._add_op2_property(prop)
            #self.properties[pid] = prop
            #print(prop.get_stats())
            #print(self.show_data(data[n-8:-100]))

            # the PBARL ends with a -1 flag
            #value, = unpack(self._endian + b'i', data[n:n+4])
            n += 4 * self.factor
        if len(self._type_to_id_map['PBAR']) == 0 and 'PBAR' in self.card_count:
            del self._type_to_id_map['PBAR']
            del self.card_count['PBAR']
        self.increase_card_count('PBARL')
        #assert len(data) == n
        return n

    def _read_pbcomp(self, data: bytes, n: int) -> int:
        """
        PBCOMP(5403, 55, 349)

                    pid      mid  A      I1      I2         I12 J           NSM
        PBCOMP         3       2 2.00E-4 6.67E-9 1.67E-9    0.0 4.58E-9     0.0 +
                 pid mid
        floats = (3, 2, 0.0002, 6.67e-09, 1.67e-09, 0.0, 4.58e-09, 0.0, 1.0, 1.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
        ints   = (3, 2, 0.0002, 6.67E-9, 1.67E-9, 0, 4.58E-9, 0, 1.0, 1.0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)

        """
        struct1 = Struct(mapfmt(self._endian + b'2i 12f i', self.size))
        struct2 = Struct(mapfmt(self._endian + b'3f 2i', self.size))
        nproperties = 0
        ntotal1 = 60 * self.factor  # 4*15
        ntotal2 = 20 * self.factor

        ndata = len(data)
        #print(ntotal1, ntotal2)
        if self.factor == 2:
            self.show_data(data[12*self.factor:], types='qd')
        #print(len(data[12*self.factor:]))
        while n < ndata:
            self.log.debug(f"n={n} ndata={ndata}")
            edata = data[n:n+ntotal1]
            #if len(edata) == ntotal1:
            data1 = struct1.unpack(edata)
            #else:
                #self.show_data(edata, types='qdi')
                #n += ntotal2
                #continue
            nsections = data1[-1]
            if self.is_debug_file:
                (pid, mid, a, i1, i2, i12, j, nsm, k1, k2,
                 m1, m2, n1, n2, unused_nsections) = data1
                self.log.info(f'PBCOMP pid={pid} mid={mid} nsections={nsections} '
                              f'k1={k1} k2={k2} m=({m1},{m2}) n=({n1},{n2})\n')
            #if pid > 0 and nsections == 0:
                #print('n1')
                #n += ntotal1
                #continue
            #if pid == 0 and nsections == 0:
                #print('n2')
                #n += ntotal2
                #continue

            data2 = []
            n += ntotal1
            if nsections in [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20]:
                # 16 Y   RS    Lumped area location along element's y-axis
                # 17 Z   RS    Lumped area location along element's z-axis
                # 18 C   RS    Fraction of the total area for the lumped area
                # 19 MID I     Material identification number
                # 20     UNDEF None
                # Words 16 through 20 repeat NSECT times
                for unused_i in range(nsections):
                    datai = data[n:n+ntotal2]
                    xi, yi, ci, mid, unused_null = struct2.unpack(datai)
                    data2.append((xi, yi, ci, mid))
                    n += ntotal2
            else:
                raise NotImplementedError('PBCOMP nsections=%r' % nsections)

            if self.is_debug_file:
                self.binary_debug.write('     PBCOMP: %s\n' % str([data1, data2]))
                msg = (
                    '    i=%-2s so=%s xxb=%.1f a=%g i1=%g i2=%g i12=%g j=%g nsm=%g '
                    'c=[%s,%s] d=[%s,%s] e=[%s,%s] f=[%s,%s]' % (
                        nsections, None, -9999., a, i1, i2, i12, j, nsm,
                        None, None, None, None, None, None, None, None,)
                )
                self.log.debug(msg)
            self.log.debug(data1)
            self.log.debug(data2)

            data_in = [data1, data2]
            prop = PBCOMP.add_op2_data(data_in)
            pid = data1[0]
            if pid in self.properties:
                self._type_to_id_map['PBEAM'].remove(pid)
                del self.properties[pid]

            self._add_op2_property(prop)
            nproperties += 1
        #print(f"n={n} ndata={ndata}")
        assert nproperties > 0, 'PBCOMP nproperties=%s' % (nproperties)
        if len(self._type_to_id_map['PBEAM']) == 0 and 'PBEAM' in self.card_count:
            del self._type_to_id_map['PBEAM']
            del self.card_count['PBEAM']
        self.card_count['PBCOMP'] = nproperties
        return n

    def _read_pbeam(self, data: bytes, n: int) -> int:
        """
        PBEAM(5402,54,262) - the marker for Record 14
        .. todo:: add object
        """
        cross_section_type_map = {
            0 : 'variable',
            1 : 'constant',
            2 : '???',
        }

        struct1 = Struct(mapfmt(self._endian + b'4if', self.size))
        struct2 = Struct(mapfmt(self._endian + b'16f', self.size))
        struct3 = Struct(mapfmt(self._endian + b'16f', self.size))
        unused_ntotal = 768 # 4*(5+16*12)
        #nproperties = (len(data) - n) // ntotal
        #assert nproperties > 0, 'ndata-n=%s n=%s datai\n%s' % (len(data)-n, n, self.show_data(data[n:100+n]))
        ndata = len(data)
        #self.show_data(data[12:], 'if')
        #assert ndata % ntotal == 0, 'ndata-n=%s n=%s ndata%%ntotal=%s' % (len(data)-n, n, ndata % ntotal)
        nproperties = 0

        ntotal1 = 20 * self.factor
        ntotal2 = 64 * self.factor
        while n < ndata:
        #while 1: #for i in range(nproperties):
            edata = data[n:n+ntotal1]
            n += ntotal1
            data_in = list(struct1.unpack(edata))
            #if self.is_debug_file:
                #self.log.info('PBEAM pid=%s mid=%s nsegments=%s ccf=%s x=%s\n' % tuple(data_in))
            (pid, unused_mid, unused_nsegments, ccf, unused_x) = data_in
            #self.log.info('PBEAM pid=%s mid=%s nsegments=%s ccf=%s x=%s' % tuple(data_in))

            # Constant cross-section flag: 1=yes and 0=no
            # what is 2?
            if ccf not in [0, 1, 2]:
                msg = ('  PBEAM pid=%s mid=%s nsegments=%s ccf=%s x=%s; '
                       'ccf must be in [0, 1, 2]\n' % tuple(data_in))
                raise ValueError(msg)

            cross_section_type = cross_section_type_map[ccf]
            #print('cross_section_type = %s' % cross_section_type)

            is_pbcomp = False
            for i in range(11):
                edata = data[n:n+ntotal2]
                if len(edata) != ntotal2:
                    endpack = []
                    raise RuntimeError('PBEAM unexpected length i=%s...' % i)
                n += ntotal2
                pack = struct2.unpack(edata)
                (soi, xxb, a, i1, i2, i12, j, nsm, c1, c2,
                 d1, d2, e1, e2, f1, f2) = pack

                if soi == 0.0:
                    so_str = 'NO'
                elif soi == 1.0:
                    so_str = 'YES'
                else:
                    if soi < 0.:
                        msg = 'PBEAM pid=%s i=%s x/xb=%s soi=%s; soi not in 0.0 or 1.0; assuming PBCOMP & dropping' % (
                            pid, i, xxb, soi)
                        self.log.error(msg)
                        is_pbcomp = True

                    so_str = str(soi)
                    #msg = 'PBEAM pid=%s i=%s x/xb=%s soi=%s; soi not in 0.0 or 1.0' % (
                        #pid, i, xxb, soi)
                    #raise NotImplementedError(msg)

                #if xxb != 0.0:
                    #msg = 'PBEAM pid=%s i=%s x/xb=%s soi=%s; xxb not in 0.0 or 1.0' % (
                        #pid, i, xxb, soi)
                    #raise NotImplementedError(msg)

                pack2 = (so_str, xxb, a, i1, i2, i12, j, nsm, c1, c2,
                         d1, d2, e1, e2, f1, f2)
                data_in.append(pack2)
                if self.is_debug_file:
                    self.binary_debug.write('     %s\n' % str(pack))
                    msg = (
                        '    i=%-2s' % i + ' so=%s xxb=%.1f a=%g i1=%g i2=%g i12=%g j=%g nsm=%g '
                        'c=[%s,%s] d=[%s,%s] e=[%s,%s] f=[%s,%s]' % (tuple(pack2))
                    )
                    self.binary_debug.write(msg)
                #msg = (
                    #'    i=%-2s' % i + ' so=%s xxb=%.1f a=%g i1=%g i2=%g i12=%g j=%g nsm=%g '
                    #'c=[%s,%s] d=[%s,%s] e=[%s,%s] f=[%s,%s]' % (tuple(pack2))
                #)
                #print(msg)

            edata = data[n:n+ntotal2]
            if len(edata) != ntotal2:
                endpack = []
                raise RuntimeError('PBEAM unexpected length 2...')
            endpack = struct3.unpack(edata)
            n += ntotal2

            assert len(endpack) == 16, endpack
            #(k1, k2, s1, s2, nsia, nsib, cwa, cwb, # 8
             #m1a, m2a, m1b, m2b, n1a, n2a, n1b, n2b) = endpack # 8 -> 16
            if self.is_debug_file:
                self.binary_debug.write('    k=[%s,%s] s=[%s,%s] nsi=[%s,%s] cw=[%s,%s] '
                                        'ma=[%s,%s] mb=[%s,%s] na=[%s,%s] nb=[%s,%s]' % (
                                            tuple(endpack)))
            data_in.append(endpack)

            if is_pbcomp:
                continue
            if pid in self.properties:
                if self.properties[pid].type == 'PBCOMP':
                    continue
            prop = PBEAM.add_op2_data(data_in)
            nproperties += 1
            self._add_op2_property(prop)
        if nproperties:
            self.card_count['PBEAM'] = nproperties
        return n

    def _read_pbeaml(self, data: bytes, n: int) -> int:
        """
        PBEAML(9202,92,53)

        Word Name Type Description
        1 PID        I   Property identification number
        2 MID        I   Material identification number
        3 GROUP(2) CHAR4 Cross-section group name
        5 TYPE(2)  CHAR4 Cross section type
        7 VALUE      RS  Cross section values for XXB, SO, NSM, and dimensions
        Word 7 repeats until (-1) occurs
        """
        #strs = numpy.core.defchararray.reshapesplit(data, sep=",")
        #ints = np.frombuffer(data[n:], self._uendian + 'i').copy()
        #floats = np.frombuffer(data[n:], self._uendian + 'f').copy()
        ints = np.frombuffer(data[n:], self.idtype8).copy()
        floats = np.frombuffer(data[n:], self.fdtype8).copy()
        iminus1 = np.where(ints == -1)[0]

        istart = [0] + list(iminus1[:-1] + 1)
        iend = iminus1

        size = self.size
        nproperties = len(istart)
        if size == 4:
            struct1 = Struct(self._endian + b'2i 8s 8s')
        else:
            struct1 = Struct(self._endian + b'2q 16s 16s')

        for unused_i, (istarti, iendi) in enumerate(zip(istart, iend)):
            idata = data[n+istarti*size : n+(istarti+6)*size]
            pid, mid, group, beam_type = struct1.unpack(idata)
            group = group.decode('latin1').strip()
            beam_type = beam_type.decode('latin1').strip()
            fvalues = floats[istarti+6: iendi]
            if self.is_debug_file:
                self.binary_debug.write('     %s\n' % str(fvalues))
                self.log.debug(f'pid={pid:d} mid={mid:d} group={group} beam_type={beam_type}')
                self.log.debug(fvalues)
            #self.log.debug(f'pid={pid:d} mid={mid:d} group={group} beam_type={beam_type}')
            data_in = [pid, mid, group, beam_type, fvalues]
            prop = PBEAML.add_op2_data(data_in)
            if pid in self.properties:
                # this is a fake PSHELL
                propi = self.properties[pid]
                assert propi.type in ['PBEAM'], propi.get_stats()
                nproperties -= 1
                continue
            self._add_op2_property(prop)
        if nproperties:
            self.card_count['PBEAML'] = nproperties
        return len(data)

    def _read_pbend(self, data: bytes, n: int) -> int:
        """PBEND"""
        n = self._read_dual_card(data, n, self._read_pbend_nx, self._read_pbend_msc,
                                 'PBEND', self._add_property_object)
        return n

    def _read_pbend_msc(self, data: bytes, n: int) -> int:
        """
        PBEND

        1 PID     I  Property identification number
        2 MID     I  Material identification number
        3 A       RS Area
        4 I1      RS Area moment of inertia in plane 1
        5 I2      RS Area moment of inertia in plane 2
        6 J       RS Torsional constant
        7 FSI     I  flexibility and stress intensification factors
        8 RM      RS Mean cross-sectional radius of the curved pipe
        9 T       RS Wall thickness of the curved pipe
        10 P      RS Internal pressure
        11 RB     RS Bend radius of the line of centroids
        12 THETAB RS Arc angle of element
        13 C1     RS Stress recovery location at point C in element y-axis
        14 C2     RS Stress recovery location at point C in element z-axis
        15 D1     RS Stress recovery location at point D in element y-axis
        16 D2     RS Stress recovery location at point D in element z-axis
        17 E1     RS Stress recovery location at point E in element y-axis
        18 E2     RS Stress recovery location at point E in element z-axis
        19 F1     RS Stress recovery location at point F in element y-axis
        20 F2     RS Stress recovery location at point F in element z-axis
        21 K1     RS Area factor for shear in plane 1
        22 K2     RS Area factor for shear in plane 2
        23 NSM    RS Nonstructural mass per unit length
        24 RC     RS Radial offset of the geometric centroid
        25 ZC     RS Offset of the geometric centroid
        26 DELTAN  I Radial offset of the neutral axis from the geometric
                     centroid
        """
        ntotal = 104  # 26*4
        struct1 = Struct(self._endian + b'2i 4f i 18f f')  # delta_n is a float, not an integer
        nproperties = (len(data) - n) // ntotal
        assert (len(data) - n) % ntotal == 0
        assert nproperties > 0, 'table=%r len=%s' % (self.table_name, len(data) - n)
        properties = []
        for unused_i in range(nproperties):
            edata = data[n:n+104]
            out = struct1.unpack(edata)
            (pid, mid, area, i1, i2, j, fsi, rm, t, p, rb, theta_b,
             c1, c2, d1, d2, e1, e2, f1, f2, k1, k2, nsm, rc, zc,
             delta_n) = out
            beam_type = fsi

            if (area, rm, t, p) == (0., 0., 0., 0.):
                area = None
                rm = None
                t = None
                p = None
                delta_n = None
                beam_type = 2
            if delta_n == 0:
                #: Radial offset of the neutral axis from the geometric
                #: centroid, positive is toward the center of curvature
                delta_n = None
            pbend = PBEND(pid, mid, beam_type, area, i1, i2, j,
                          c1, c2, d1, d2, e1, e2, f1, f2, k1, k2,
                          nsm, rc, zc, delta_n, fsi, rm, t, p, rb, theta_b)
            #print(pbend)
            pbend.validate()

            properties.append(pbend)
            n += ntotal
        return n, properties

    def _read_pbend_nx(self, data: bytes, n: int) -> int:
        """
        PBEND

        1 PID     I  Property identification number
        2 MID     I  Material identification number
        3 A       RS Area
        4 I1      RS Area moment of inertia in plane 1
        5 I2      RS Area moment of inertia in plane 2
        6 J       RS Torsional constant
        7 FSI     I  Flexibility and stress intensification factors
        8 RM      RS Mean cross-sectional radius of the curved pipe
        9 T       RS Wall thickness of the curved pipe
        10 P      RS Internal pressure
        11 RB     RS Bend radius of the line of centroids
        12 THETAB RS Arc angle of element
        13 C1     RS Stress recovery location at point C in element y-axis
        14 C2     RS Stress recovery location at point C in element z-axis
        15 D1     RS Stress recovery location at point D in element y-axis
        16 D2     RS Stress recovery location at point D in element z-axis
        17 E1     RS Stress recovery location at point E in element y-axis
        18 E2     RS Stress recovery location at point E in element z-axis
        19 F1     RS Stress recovery location at point F in element y-axis
        20 F2     RS Stress recovery location at point F in element z-axis
        21 K1     RS Area factor for shear in plane 1
        22 K2     RS Area factor for shear in plane 2
        23 NSM    RS Nonstructural mass per unit length
        24 RC     RS Radial offset of the geometric centroid
        25 ZC     RS Offset of the geometric centroid
        26 DELTAN RS Radial offset of the neutral axis from the geometric
                     centroid
        27 SACL   RS Miter spacing at center line.
        28 ALPHA  RS One-half angle between the adjacent miter axis
                     (Degrees).
        29 FLANGE I  For FSI=5, defines the number of flanges attached.
        30 KX     RS For FSI=6, the user defined flexibility factor for the
                  torsional moment.
        31 KY     RS For FSI=6, the user defined flexibility factor for the
                  out-of-plane bending moment.
        32 KZ     RS For FSI=6, the user defined flexbility factor for the
                  in-plane bending moment.
        33 Not used
        """
        #self.log.info('skipping PBEND in EPT')
        #return len(data)
        ntotal = 132  # 33*4
        struct1 = Struct(self._endian + b'2i 4f i 21f i 4f')
        nproperties = (len(data) - n) // ntotal
        assert (len(data) - n) % ntotal == 0
        assert nproperties > 0, 'table=%r len=%s' % (self.table_name, len(data) - n)
        properties = []
        for unused_i in range(nproperties):
            edata = data[n:n+132]
            out = struct1.unpack(edata)
            (pid, mid, area, i1, i2, j, fsi, rm, t, p, rb, theta_b,
             c1, c2, d1, d2, e1, e2, f1, f2, k1, k2, nsm, rc, zc,
             delta_n, unused_sacl, unused_alpha, unused_flange,
             unused_kx, unused_ky, unused_kz, unused_junk,) = out
            beam_type = fsi

            pbend = PBEND(pid, mid, beam_type, area, i1, i2, j,
                          c1, c2, d1, d2, e1, e2, f1, f2, k1, k2,
                          nsm, rc, zc, delta_n, fsi, rm, t, p, rb, theta_b)
            pbend.validate()
            properties.append(pbend)
            n += ntotal
        return n, properties

# PBMSECT
# PBRSECT

    def _read_pbush(self, data: bytes, n: int) -> int:
        """
        The PBUSH card is different between MSC and NX Nastran.

        DMAP NX 11
        ----------
        NX has 23 fields in NX 11 (supported)
        NX has 18 fields in the pre-2001 format (supported)

        DMAP MSC 2005
        -------------
        MSC has 23 fields in 2005 (supported)
        MSC has 18 fields in the pre-2001 format (supported)

        DMAP MSC 2016
        -------------
        TODO: MSC has 24 fields in 2016.1 (not supported)
        MSC has 18 fields in the pre-2001 format (supported)
        """
        # we're listing nx twice because NX/MSC used to be consistent
        # the new form for MSC is not supported
        n = self._read_dual_card(data, n, self._read_pbush_nx, self._read_pbush_msc,
                                 'PBUSH', self._add_op2_property)
        return n

    def _read_pbush_nx(self, data: bytes, n: int) -> int:
        """PBUSH(1402,14,37) - 18 fields"""
        ntotal = 72 * self.factor
        struct1 = Struct(mapfmt(self._endian + b'i17f', self.size))
        ndata = len(data) - n
        nentries = ndata // ntotal
        assert nentries > 0, 'table={self.table_name} len={ndata}'
        assert ndata % ntotal == 0, f'table={self.table_name} leftover = {ndata} % {ntotal} = {ndata % ntotal}'
        props = []
        for unused_i in range(nentries):
            edata = data[n:n+ntotal]
            out = struct1.unpack(edata)
            (pid, k1, k2, k3, k4, k5, k6, b1, b2, b3, b4, b5, b6,
             g1, sa, st, ea, et) = out
            #self.log.debug(out)
            g2 = g3 = g4 = g5 = g6 = g1
            data_in = (pid, k1, k2, k3, k4, k5, k6, b1, b2, b3, b4, b5, b6,
                       g1, g2, g3, g4, g5, g6, sa, st, ea, et)
            prop = PBUSH.add_op2_data(data_in)
            props.append(prop)
            n += ntotal
        return n, props

    def _read_pbush_msc(self, data: bytes, n: int) -> int:
        """PBUSH(1402,14,37) - 23 fields"""
        ntotal = 92 * self.factor # 23*4
        struct1 = Struct(mapfmt(self._endian + b'i22f', self.size))

        ndata = len(data) - n
        nentries = ndata // ntotal
        assert nentries > 0, 'table={self.table_name} len={ndata}'
        assert ndata % ntotal == 0, f'table={self.table_name} leftover = {ndata} % {ntotal} = {ndata % ntotal}'

        props = []
        for unused_i in range(nentries):
            edata = data[n:n+ntotal]
            out = struct1.unpack(edata)
            #(pid, k1, k2, k3, k4, k5, k6, b1, b2, b3, b4, b5, b6,
             #g1, g2, g3, g4, g5, g6, sa, st, ea, et) = out
            prop = PBUSH.add_op2_data(out)
            props.append(prop)
            n += ntotal
        return n, props

    def _read_pbush1d(self, data: bytes, n: int) -> int:
        """
        Record 18 -- PBUSH1D(3101,31,219)

        1  PID    I  Property identification number
        2  K      RS Stiffness
        3  C      RS Viscous Damping
        4  M      RS Mass
        5  ALPHA  RS Temperature coefficient
        6  SA     RS Stress recovery coefficient
        7  EA/SE  RS Strain recovery coefficient

        8  TYPEA  I  Shock data type:0=Null, 1=Table, 2=Equation
        9  CVT    RS Coefficient of translation velocity tension
        10 CVC    RS Coefficient of translation velocity compression
        11 EXPVT  RS Exponent of velocity tension
        12 EXPVC  RS Exponent of velocity compression
        13 IDTSU  I  TABLEDi or DEQATN entry identification number for scale factor vs displacement
        14 IDTCU  I  DEQATN entry identification number for scale factor vs displacement
        15 IDTSUD I  DEQATN entry identification number for derivative tension
        16 IDCSUD I  DEQATN entry identification number for derivative compression

        17 TYPES  I  Spring data type: 0=Null, 1=Table, 2=Equation
        18 IDTS   I  TABLEDi or DEQATN entry identification number for tension compression
        19 IDCS   I  DEQATN entry identification number for compression
        20 IDTDU  I  DEQATN entry identification number for scale factor vs displacement
        21 IDCDU  I  DEQATN entry identification number for force vs displacement

        22 TYPED  I  Damper data type: 0=Null, 1=Table, 2=Equation
        23 IDTD   I  TABLEDi or DEQATN entry identification number for tension compression
        24 IDCD   I  DEQATN entry identification number for compression
        25 IDTDV  I  DEQATN entry identification number for scale factor versus velocity
        26 IDCDV  I  DEQATN entry identification number for force versus velocity

        27 TYPEG  I  General data type: 0=Null, 1=Table, 2=Equation
        28 IDTG   I  TABLEDi or DEQATN entry identification number for tension compression
        29 IDCG   I  DEQATN entry identification number for compression
        30 IDTDU  I  DEQATN entry identification number for scale factor versus displacement
        31 IDCDU  I  DEQATN entry identification number for force versus displacement
        32 IDTDV  I  DEQATN entry identification number for scale factor versus velocity
        33 IDCDV  I  DEQATN entry identification number for force vs velocity

        34 TYPEF  I  Fuse data type: 0=Null, 1=Table
        35 IDTF   I  TABLEDi entry identification number for tension
        36 IDCF   I  TABLEDi entry identification number for compression

        37 UT     RS Ultimate tension
        38 UC     RS Ultimate compression
        """
        type_map = {
            0 : None,  # NULL
            1 : 'EQUAT',
            2 : 'TABLE',
        }
        ntotal = 152 * self.factor  # 38*4
        struct1 = Struct(mapfmt(self._endian + b'i 6f i 4f 24i 2f', self.size))
        nentries = (len(data) - n) // ntotal
        for unused_i in range(nentries):
            edata = data[n:n+ntotal]
            out = struct1.unpack(edata)
            (pid, k, c, m, unused_alpha, sa, se,
             typea, cvt, cvc, expvt, expvc, idtsu, idtcu, idtsud, idcsud,
             types, idts, idcs, idtdus, idcdus,
             typed, idtd, idcd, idtdvd, idcdvd,
             typeg, idtg, idcg, idtdug, idcdug, idtdvg, idcdvg,
             typef, idtf, idcf,
             unused_ut, unused_uc) = out
            #  test_op2_other_05
            #pbush1d, 204, 1.e+5, 1000., , , , , , +pb1
            #+pb1, spring, table, 205, , , , , , +pb2
            #+pb2, damper, table, 206
            #pid=204 k=100000.0 c=1000.0 m=0.0 sa=nan se=nan


            msg = f'PBUSH1D pid={pid} k={k} c={c} m={m} sa={sa} se={se}'
            optional_vars = {}
            typea_str = type_map[typea]
            types_str = type_map[types]
            typed_str = type_map[typed]
            unused_typeg_str = type_map[typeg]
            unused_typef_str = type_map[typef]

            if min([typea, types, typed, typeg, typef]) < 0:
                raise RuntimeError(f'typea={typea} types={types} typed={typed} typeg={typeg} typef={typef}')
            if typea in [1, 2]:  # SHOCKA?
                #pbush1d, 204, 1.e+5, 1000., , , , , , +pb4
                #+pb4, shocka, table, 1000., , 1., , 214, , +pb41
                #+pb41, spring, table, 205

                idts = idtsu if typea_str == 'TABLE' else 0
                idets = idtsu if typea_str == 'EQUAT' else 0
                optional_vars['SHOCKA'] = [typea_str, cvt, cvc, expvt, expvc,
                                           idts, idets, idtcu, idtsud, idcsud]
                #(shock_type, shock_cvt, shock_cvc, shock_exp_vt, shock_exp_vc,
                 #shock_idts, shock_idets, shock_idecs, shock_idetsd, shock_idecsd
                #)
                #print('shock_idts, shock_idets', typea_str, idtsu, idtsu)
                msg += (
                    f'  SHOCKA type={typea} cvt={cvt} cvc={cvc} expvt={expvt} expvc={expvc}\n'
                    f'    idtsu={idtsu} (idts={idts} idets={idets}) idtcu={idtcu} idtsud={idtsud} idcsud={idcsud}')
            if types in [1, 2]: # SPRING: Spring data type: 0=Null, 1=Table, 2=Equation
                #(spring_type, spring_idt, spring_idc, spring_idtdu, spring_idcdu) = values
                # SPRING, TYPE IDT IDC IDTDU IDCDU
                optional_vars['SPRING'] = [types_str, idts, idcs, idtdus, idcdus]
                msg += f'  SPRING type={types} idt={idts} idc={idcs} idtdu={idtdus} idcdu={idcdus}'
            if typed in [1, 2]: # Damper data type: 0=Null, 1=Table, 2=Equation
                optional_vars['DAMPER'] = [typed_str, idtd, idcd, idtdvd, idcdvd]
                msg += f'  DAMPER type={typed} idt={idtd} idc={idtd} idtdv={idtdvd} idcdv={idcdvd}'
            if typeg in [1, 2]: # general, GENER?: 0=Null, 1=Table 2=Equation
                # C:\NASA\m4\formats\git\examples\move_tpl\ar29scbt.bdf
                #pbush1d, 206, 1.e+3, 10., , , , , , +pb6
                #+pb6, gener, equat, 315, , 3015, , 3016
                msg += f'  GENER  type={typeg} idt={idtg} idc={idcg} idtdu={idtdug} idcdu={idcdug} idtdv={idtdvg} idcdv={idcdvg}'
                optional_vars['GENER'] = [idtg, idcg, idtdug, idcdug, idtdvg, idcdvg]
            if typef in [1, 2]: # Fuse data type: 0=Null, 1=Table
                raise NotImplementedError(f'typef={typef} idtf={idtf} idcf={idcf}')

            if self.is_debug_file:
                self.binary_debug.write(msg)

            self.add_pbush1d(pid, k=k, c=c, m=m, sa=sa, se=se,
                             optional_vars=optional_vars,)
            n += ntotal
        self.card_count['PBUSH1D'] = nentries
        return n

    def _read_pbusht(self, data: bytes, n: int) -> int:
        """reads the PBUSHT"""
        n, props = self._read_pbusht_nx(data, n)
        for prop in props:
            #print(prop)
            self._add_pbusht_object(prop)
        return n

    def _read_pbusht_nx(self, data: bytes, n: int) -> int:
        """
        NX 12 / MSC 2005
        Word Name Type Description
        1 PID       I Property identification number
        2 TKID(6)   I TABLEDi entry identification numbers for stiffness
        8 TBID(6)   I TABLEDi entry identification numbers for viscous damping
        14 TGEID(6) I TABLEDi entry identification number for structural damping
        20 TKNID(6) I TABLEDi entry identification numbers for force versus deflection

        old style
        Word Name Type Description
        1 PID       I Property identification number
        2 TKID(6)   I TABLEDi entry identification numbers for stiffness
        8 TBID(6)   I TABLEDi entry identification numbers for viscous damping
        14 TGEID    I TABLEDi entry identification number for structural damping
        15 TKNID(6) I TABLEDi entry IDs for force versus deflection
        """
        #self.show_data(data[12:])
        ndata = (len(data) - n) // self.factor

        if ndata % 100 == 0 and ndata % 80 == 0:
            self.log.warning(f"skipping PBUSHT in EPT because nfields={ndata//4}, which is "
                             'nproperties*25 or nproperties*20')
            return len(data), []
        if ndata % 100 == 0:
            n, props = self._read_pbusht_100(data, n)
        elif ndata % 80 == 0:
            n, props = self._read_pbusht_80(data, n)
        return n, props

    def _read_pbusht_80(self, data: bytes, n: int) -> int:
        ntotal = 80 * self.factor
        struct1 = Struct(self._endian + b'20i')
        nentries = (len(data) - n) // ntotal
        assert nentries > 0, 'table=%r len=%s' % (self.table_name, len(data) - n)

        props = []
        for unused_i in range(nentries):
            edata = data[n:n+ntotal]
            out = struct1.unpack(edata)
            #(pid,
             #k1, k2, k3, k4, k5, k6,
             #b1, b2, b3, b4, b5, b6,
             #g1, sa, st, ea, et) = out
            (pid,
             k1, k2, k3, k4, k5, k6,
             b1, b2, b3, b4, b5, b6,
             g1,
             n1, n2, n3, n4, n5, n6) = out
            g2 = g3 = g4 = g5 = g6 = g1
            k_tables = [k1, k2, k3, k4, k5, k6]
            b_tables = [b1, b2, b3, b4, b5, b6]
            ge_tables = [g1, g2, g3, g4, g5, g6]
            kn_tables = [n1, n2, n3, n4, n5, n6]
            prop = PBUSHT(pid, k_tables, b_tables, ge_tables, kn_tables)
            props.append(prop)
            n += ntotal
        return n, props

    def _read_pbusht_100(self, data: bytes, n: int) -> int:
        props = []
        ntotal = 100 * self.factor
        struct1 = Struct(mapfmt(self._endian + b'25i', self.size))
        nentries = (len(data) - n) // ntotal
        assert nentries > 0, 'table=%r len=%s' % (self.table_name, len(data) - n)
        for unused_i in range(nentries):
            edata = data[n:n+ntotal]
            out = struct1.unpack(edata)
            (pid,
             k1, k2, k3, k4, k5, k6,
             b1, b2, b3, b4, b5, b6,
             g1, g2, g3, g4, g5, g6,
             n1, n2, n3, n4, n5, n6) = out
            k_tables = [k1, k2, k3, k4, k5, k6]
            b_tables = [b1, b2, b3, b4, b5, b6]
            ge_tables = [g1, g2, g3, g4, g5, g6]
            kn_tables = [n1, n2, n3, n4, n5, n6]
            prop = PBUSHT(pid, k_tables, b_tables, ge_tables, kn_tables)
            props.append(prop)
            n += ntotal
        return n, props

    def _read_pcomp(self, data: bytes, n: int) -> int:
        """PCOMP(2706,27,287) - the marker for Record 22

        standard:
          EPTS; 64-bit: C:\MSC.Software\simcenter_nastran_2019.2\tpl_post1\cqrdbxdra3lg.op2
        """
        if self.size == 4:
            n2, props = self._read_pcomp_32_bit(data, n)
            nproperties = len(props)
            for prop in props:
                self._add_op2_property(prop)
            self.card_count['PCOMP'] = nproperties
        else:
            n2 = self._read_dual_card(data, n, self._read_pcomp_32_bit,
                                      self._read_pcomp_64_bit,
                                      'PCOMP', self._add_op2_property)
        return n2

    def _read_pcomp_64_bit(self, data: bytes, n: int) -> Tuple[int, List[PCOMP]]:
        """
        PCOMP(2706,27,287) - the marker for Record 22

        1  PID   I  Property identification number
        2  N(C)  I  Number of plies
        3  Z0    RS Distance from the reference plane to the bottom surface
        4  NSM   RS Nonstructural mass per unit area
        5  SB    RS Allowable shear stress of the bonding material
        6  FT    I  Failure theory
        7  TREF  RS Reference temperature
        8  GE    RS Damping coefficient

        9  MID   I  Material identification number
        10 T     RS Thicknesses of the ply
        11 THETA RS Orientation angle of the longitudinal direction of the ply
        12 SOUT  I Stress or strain output request of the ply
        Words 9 through 12 repeat N times

        TODO:
           64-bit bug: why is the number of plies 0???

          doubles (float64) = (
          1, 0.0, 1.7368e-18, 0.0, 1.0, 1.5e-323, 0.0, 0.0,
            1, 0.11, 0, 1,
            1, 0.11, 0, 1,
            1, 0.11, 0, 1,
          -1, -1, -1, -1,
          21, 0.0, 1.7368e-18, 0.0, 1.0, 1.5e-323, 0.0, 0.0,
            1, 0.11, 0, 1,
            1, 0.11, 0, 1,
            1, 0.11, 0, 1,
            1, 0.11, 0, 1,
          -1, -1, -1, -1)
          long long (int64) = (
          1, 0,   1.7368e-18, 0,   1.0, 3, 0, 0, 1, 4592590756007337001, 0, 1,
            1, 0.11, 0, 1,
            1, 0.11, 0, 1,
            1, 0.11, 0, 1,
          -1, -1, -1, -1,
          21, 0, 4341475431749739292, 0, 4607182418800017408, 3, 0, 0,
            1, 0.11, 0, 1,
            1, 0.11, 0, 1,
            1, 0.11, 0, 1,
            1, 0.11, 0, 1,
          -1, -1, -1, -1)

          doubles (float64) = (5e-324, 0.0, -0.005, 0.0, 0.0, 0.0, 0.0, 0.0,
                               4e-323, 0.005, 0.0, 5e-324,
                               4e-323, 0.005, 0.0, 5e-324,
                               nan, nan, nan, nan)
          long long (int64) = (1, 0, -4650957407178058629, 0, 0, 0, 0, 0,
                                  8, 4572414629676717179, 0, 1,
                                  8, 4572414629676717179, 0, 1,
                               -1, -1, -1, -1)

        """
        self.to_nx()
        nproperties = 0
        s1 = Struct(mapfmt(self._endian + b'2i3fi2f', self.size))
        ntotal1 = 32 * self.factor
        s2 = Struct(mapfmt(self._endian + b'i2fi', self.size))

        four_minus1 = Struct(mapfmt(self._endian + b'4i', self.size))
        ndata = len(data)
        ntotal2 = 16 * self.factor
        props = []
        while n < (ndata - ntotal1):
            out = s1.unpack(data[n:n+ntotal1])
            (pid, nlayers, z0, nsm, sb, ft, tref, ge) = out
            assert pid > 0
            if self.binary_debug:
                self.binary_debug.write(f'PCOMP pid={pid} nlayers={nlayers} z0={z0} nsm={nsm} '
                                        f'sb={sb} ft={ft} Tref={tref} ge={ge}')
            assert isinstance(nlayers, int), out
            #print(f'PCOMP pid={pid} nlayers={nlayers} z0={z0} nsm={nsm} '
                  #f'sb={sb} ft={ft} Tref={tref} ge={ge}')
            n += ntotal1

            # None, 'SYM', 'MEM', 'BEND', 'SMEAR', 'SMCORE', 'NO'
            is_symmetrical = 'NO'
            #if nlayers < 0:
                #is_symmetrical = 'SYM'
                #nlayers = abs(nlayers)

            mids = []
            T = []
            thetas = []
            souts = []
            edata2 = data[n:n+ntotal2]
            idata = four_minus1.unpack(edata2)
            while idata != (-1, -1, -1, -1):
                (mid, t, theta, sout) = s2.unpack(edata2)
                mids.append(mid)
                T.append(t)
                thetas.append(theta)
                souts.append(sout)
                if self.is_debug_file:
                    self.binary_debug.write(f'      mid={mid} t={t} theta={theta} sout={sout}\n')
                n += ntotal2
                #print(f'      mid={mid} t={t} theta={theta} sout={sout}')
                edata2 = data[n:n+ntotal2]
                idata = four_minus1.unpack(edata2)

            if self.size == 4:
                assert 0 < nlayers < 100, 'pid=%s nlayers=%s z0=%s nms=%s sb=%s ft=%s Tref=%s ge=%s' % (
                    pid, nlayers, z0, nsm, sb, ft, tref, ge)
            else:
                assert nlayers == 0, nlayers
                nlayers = len(mids)

            data_in = [
                pid, z0, nsm, sb, ft, tref, ge,
                is_symmetrical, mids, T, thetas, souts]
            prop = PCOMP.add_op2_data(data_in)
            nproperties += 1
            n += ntotal2
            props.append(prop)
        return n, props

    def _read_pcomp_32_bit(self, data: bytes, n: int) -> Tuple[int, List[PCOMP]]:  # pragma: no cover
        """PCOMP(2706,27,287) - the marker for Record 22"""
        nproperties = 0
        s1 = Struct(mapfmt(self._endian + b'2i3fi2f', self.size))
        ntotal1 = 32 * self.factor
        s2 = Struct(mapfmt(self._endian + b'i2fi', self.size))

        ndata = len(data)
        ntotal2 = 16 * self.factor
        props = []
        while n < (ndata - ntotal1):
            out = s1.unpack(data[n:n+ntotal1])
            (pid, nlayers, z0, nsm, sb, ft, tref, ge) = out
            assert pid > 0

            if self.binary_debug:
                self.binary_debug.write(f'PCOMP pid={pid} nlayers={nlayers} z0={z0} nsm={nsm} '
                                        f'sb={sb} ft={ft} Tref={tref} ge={ge}')
            assert isinstance(nlayers, int), out
            #print(f'PCOMP pid={pid} nlayers={nlayers} z0={z0} nsm={nsm} '
                  #f'sb={sb} ft={ft} Tref={tref} ge={ge}')
            n += ntotal1

            mids = []
            T = []
            thetas = []
            souts = []

            # None, 'SYM', 'MEM', 'BEND', 'SMEAR', 'SMCORE', 'NO'
            is_symmetrical = 'NO'
            if nlayers < 0:
                is_symmetrical = 'SYM'
                nlayers = abs(nlayers)
            assert nlayers > 0, out

            assert 0 < nlayers < 100, 'pid=%s nlayers=%s z0=%s nms=%s sb=%s ft=%s Tref=%s ge=%s' % (
                pid, nlayers, z0, nsm, sb, ft, tref, ge)

            if self.is_debug_file:
                self.binary_debug.write('    pid=%s nlayers=%s z0=%s nms=%s sb=%s ft=%s Tref=%s ge=%s\n' % (
                    pid, nlayers, z0, nsm, sb, ft, tref, ge))
            for unused_ilayer in range(nlayers):
                (mid, t, theta, sout) = s2.unpack(data[n:n+ntotal2])
                mids.append(mid)
                assert mid > 0

                T.append(t)
                thetas.append(theta)
                souts.append(sout)
                if self.is_debug_file:
                    self.binary_debug.write(f'      mid={mid} t={t} theta={theta} sout={sout}\n')
                n += ntotal2
                #print(f'      mid={mid} t={t} theta={theta} sout={sout}\n')

            data_in = [
                pid, z0, nsm, sb, ft, tref, ge,
                is_symmetrical, mids, T, thetas, souts]
            prop = PCOMP.add_op2_data(data_in)
            props.append(prop)
            nproperties += 1
        return n, props

    def _read_pcompg(self, data: bytes, n: int) -> int:
        """
        PCOMP(2706,27,287)

        1 PID      I  Property identification number
        2 LAMOPT   I  Laminate option
        3 Z0       RS Distance from the reference plane to the bottom surface
        4 NSM      RS Nonstructural mass per unit area
        5 SB       RS Allowable shear stress of the bonding material
        6 FT       I  Failure theory
        7 TREF     RS Reference temperature
        8 GE       RS Damping coefficient

        9  GPLYIDi I  Global ply IDs.
        10 MID     I  Material identification number
        11 T       RS Thicknesses of the ply
        12 THETA   RS Orientation angle of the longitudinal direction of the ply
        13 SOUT    I  Stress or strain output request of the ply
        Words 9 through 13 repeat N times (until -1, -1, -1, -1, -1 as Nplies doesn't exist...)

        float = (15006, 150, 604,
                 5, 0.0, 1.7368e-18, 0.0, 0.0, 0.0, 20.0, 0.0,
                     5e-324, 5e-324, 2.0, 0.0, 0.0,
                     1e-323, 1e-323, 3.0, 0.0, 0.0,
                     1.5e-323, 1e-323, 3.0, 0.0, 0.0,
                     2e-323, 5e-324, 2.0, 0.0, 0.0,
                     nan, nan, nan, nan, nan)
        int   = (15006, 150, 604,
                 5, 0,   1.7368e-18, 0,   0,   0,   20.0, 0,
                     1, 1, 4611686018427387904, 0, 0,
                     2, 2, 4613937818241073152, 0, 0,
                     3, 2, 4613937818241073152, 0, 0,
                     4, 1, 4611686018427387904, 0, 0,
                     -1, -1, -1, -1, -1)

        """
        nproperties = 0
        s1 = Struct(mapfmt(self._endian + b'2i 3f i 2f', self.size))
        s2 = Struct(mapfmt(self._endian + b'2i 2f i', self.size))
        struct_i5 = Struct(mapfmt(self._endian + b'5i', self.size))

        # lam - SYM, MEM, BEND, SMEAR, SMCORE, None
        lam_map = {
            0 : None,
            # MEM
            # BEND
            # SMEAR
            # SMCORE
        }

        # ft - HILL, HOFF, TSAI, STRN, None
        ft_map = {
            0 : None,
            # HILL
            # HOFF
            3 : 'TSAI',
            # STRN
        }
        # sout - YES, NO
        sout_map = {
            0 : 'NO',
            1 : 'YES',
        }
        ndata = len(data)
        #self.show_data(data, types='qd')
        ntotal1 = 32 * self.factor
        ntotal2 = 20 * self.factor
        while n < (ndata - ntotal1):
            out = s1.unpack(data[n:n+ntotal1])
            (pid, lam_int, z0, nsm, sb, ft_int, tref, ge) = out
            if self.binary_debug:
                self.binary_debug.write(f'PCOMPG pid={pid} lam_int={lam_int} z0={z0} nsm={nsm} '
                                        f'sb={sb} ft_int={ft_int} tref={tref} ge={ge}')
            #print(f'PCOMPG pid={pid} lam_int={lam_int} z0={z0} nsm={nsm} sb={sb} '
                  #f'ft_int={ft_int} tref={tref} ge={ge}')
            assert isinstance(lam_int, int), out
            assert pid > -1, out
            n += ntotal1

            mids = []
            thicknesses = []
            thetas = []
            souts = []
            global_ply_ids = []

            # None, 'SYM', 'MEM', 'BEND', 'SMEAR', 'SMCORE', 'NO'
            #is_symmetrical = 'NO'
            #if nlayers < 0:
                #is_symmetrical = 'SYM'
                #nlayers = abs(nlayers)
            #assert nlayers > 0, out

            #assert 0 < nlayers < 100, 'pid=%s nlayers=%s z0=%s nms=%s sb=%s ft=%s tref=%s ge=%s' % (
                #pid, nlayers, z0, nsm, sb, ft, tref, ge)

            #if self.is_debug_file:
                #self.binary_debug.write('    pid=%s nlayers=%s z0=%s nms=%s sb=%s ft=%s tref=%s ge=%s\n' % (
                    #pid, nlayers, z0, nsm, sb, ft, tref, ge))
            ilayer = 0
            while ilayer < 1000:
                ints5 = struct_i5.unpack(data[n:n+ntotal2])
                if ints5 == (-1, -1, -1, -1, -1):
                    if self.is_debug_file:
                        self.binary_debug.write('      global_ply=%-1 mid=%-1 t=%-1 theta=%-1 sout=-1\n')
                    break
                (global_ply, mid, t, theta, sout_int) = s2.unpack(data[n:n+ntotal2])
                #print('  ', (global_ply, mid, t, theta, sout_int))
                try:
                    sout = sout_map[sout_int]
                except KeyError:
                    self.log.error('cant parse global_ply=%s sout=%s; assuming 0=NO' % (
                        global_ply, sout_int))
                    sout = 'NO'

                global_ply_ids.append(global_ply)
                mids.append(mid)
                thicknesses.append(t)
                thetas.append(theta)
                souts.append(sout)
                if self.is_debug_file:
                    self.binary_debug.write('      global_ply=%s mid=%s t=%s theta=%s sout_int=%s sout=%r\n' % (
                        global_ply, mid, t, theta, sout_int, sout))
                n += ntotal2
                ilayer += 1
            n += ntotal2

            try:
                ft = ft_map[ft_int]
            except KeyError:
                self.log.error('pid=%s cant parse ft=%s; should be HILL, HOFF, TSAI, STRN'
                               '...skipping' % (pid, ft_int))
                continue

            try:
                lam = lam_map[lam_int]
            except KeyError:
                self.log.error('pid=%s cant parse lam=%s; should be HILL, HOFF, TSAI, STRN'
                               '...skipping' % (pid, lam_int))
                continue

            # apparently Nastran makes duplicate property ids...
            if pid in self.properties and self.properties[pid].type == 'PCOMP':
                del self.properties[pid]

            self.add_pcompg(pid, global_ply_ids, mids, thicknesses, thetas=thetas, souts=souts,
                            nsm=nsm, sb=sb, ft=ft, tref=tref, ge=ge, lam=lam, z0=z0, comment='')
            nproperties += 1
        self.card_count['PCOMPG'] = nproperties
        return n

# PCOMPA

    def _read_pconeax(self, data: bytes, n: int) -> int:
        """
        (152,19,147) - Record 24
        """
        self.log.info('skipping PCONEAX in EPT')
        return len(data)

    def _read_pconv(self, data: bytes, n: int) -> int:
        """common method for reading PCONVs"""
        n = self._read_dual_card(data, n, self._read_pconv_nx, self._read_pconv_msc,
                                 'PCONV', self._add_pconv)
        return n

    def _read_pconv_nx(self, data: bytes, n: int) -> int:
        """
        (11001,110,411)- NX version
        """
        ntotal = 16  # 4*4
        struct_3if = Struct(self._endian + b'3if')
        nentries = (len(data) - n) // ntotal
        assert (len(data) - n) % ntotal == 0
        props = []
        for unused_i in range(nentries):
            out = struct_3if.unpack(data[n:n+ntotal])
            (pconid, mid, form, expf) = out
            ftype = tid = chlen = gidin = ce = e1 = e2 = e3 = None
            data_in = (pconid, mid, form, expf, ftype, tid, chlen,
                       gidin, ce, e1, e2, e3)

            prop = PCONV.add_op2_data(data_in)
            props.append(prop)
            n += ntotal
        return n, props

    def _read_pconv_msc(self, data: bytes, n: int) -> int:
        """
        (11001,110,411)- MSC version - Record 25
        """
        ntotal = 56  # 14*4
        s = Struct(self._endian + b'3if 4i fii 3f')
        nentries = (len(data) - n) // ntotal
        assert (len(data) - n) % ntotal == 0
        props = []
        for unused_i in range(nentries):
            out = s.unpack(data[n:n+ntotal])
            (pconid, mid, form, expf, ftype, tid, unused_undef1, unused_undef2, chlen,
             gidin, ce, e1, e2, e3) = out
            data_in = (pconid, mid, form, expf, ftype, tid, chlen,
                       gidin, ce, e1, e2, e3)

            prop = PCONV.add_op2_data(data_in)
            props.append(prop)
            n += ntotal
        return n, props

    def _read_pconvm(self, data: bytes, n: int) -> int:
        """Record 24 -- PCONVM(2902,29,420)

        1 PID    I Property identification number
        2 MID    I Material identification number
        3 FORM   I Type of formula used for free convection
        4 FLAG   I Flag for mass flow convection
        5 COEF  RS Constant coefficient used for forced convection
        6 EXPR  RS Reynolds number convection exponent
        7 EXPPI RS Prandtl number convection exponent into the working fluid
        8 EXPPO RS Prandtl number convection exponent out of the working fluid
        """
        ntotal = 32  # 8*4
        structi = Struct(self._endian + b'4i 4f')
        nentries = (len(data) - n) // ntotal
        for unused_i in range(nentries):
            out = structi.unpack(data[n:n+ntotal])
            if out != (0, 0, 0, 0, 0., 0., 0., 0.):
                (pconid, mid, form, flag, coeff, expr, expri, exppo) = out
                #print(out)
                prop = PCONVM(pconid, mid, coeff, form=form, flag=flag,
                              expr=expr, exppi=expri, exppo=exppo, comment='')
                self._add_convection_property_object(prop)
            n += ntotal
        self.card_count['PCONVM'] = nentries
        return n

    def _read_pdamp(self, data: bytes, n: int) -> int:
        """
        PDAMP(202,2,45) - the marker for Record ???
        """
        ntotal = 8 * self.factor # 2*4
        struct_if = Struct(mapfmt(self._endian + b'if', self.size))
        nentries = (len(data) - n) // ntotal
        for unused_i in range(nentries):
            out = struct_if.unpack(data[n:n+ntotal])
            #(pid, b) = out
            prop = PDAMP.add_op2_data(out)
            self._add_op2_property(prop)
            n += ntotal
        self.card_count['PDAMP'] = nentries
        return n

    def _read_pdampt(self, data: bytes, n: int) -> int:  # 26
        self.log.info('skipping PDAMPT in EPT')
        return len(data)

    def _read_pdamp5(self, data: bytes, n: int) -> int:  # 26
        self.log.info('skipping PDAMP5 in EPT')
        return len(data)

# PDUM1
# PDUM2
# PDUM3
# PDUM4
# PDUM5
# PDUM6
# PDUM7
# PDUM8
# PDUM9

    def _read_pelas(self, data: bytes, n: int) -> int:
        """PELAS(302,3,46) - the marker for Record 39"""
        struct_i3f = Struct(mapfmt(self._endian + b'i3f', self.size))
        ntotal = 16 * self.factor # 4*4
        nproperties = (len(data) - n) // ntotal
        for unused_i in range(nproperties):
            edata = data[n:n+ntotal]
            out = struct_i3f.unpack(edata)
            #(pid, k, ge, s) = out
            if self.is_debug_file:
                self.binary_debug.write('  PELAS=%s\n' % str(out))
            prop = PELAS.add_op2_data(out)
            self._add_op2_property(prop)
            n += ntotal
        self.card_count['PELAS'] = nproperties
        return n

    def _read_pfast_msc(self, data: bytes, n: int) -> int:
        ntotal = 92 # 23*4
        struct1 = Struct(self._endian + b'ifii 4f')
        nproperties = (len(data) - n) // ntotal
        for unused_i in range(nproperties):
            edata = data[n:n+ntotal]
            out = struct1.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  PFAST=%s\n' % str(out))
            (pid, d, mcid, unused_connbeh, unused_conntype, unused_extcon,
             unused_condtype, unused_weldtype, unused_minlen, unused_maxlen,
             unused_gmcheck, unused_spcgs, mass, ge,
             unused_aa, unused_bb, unused_cc, mcid, mflag,
             kt1, kt2, kt3, kr1, kr2, kr3) = out

            data_in = (pid, d, mcid, mflag, kt1, kt2, kt3,
                       kr1, kr2, kr3, mass, ge)
            prop = PFAST.add_op2_data(data_in)
            self._add_op2_property(prop)
            n += ntotal
        self.card_count['PFAST'] = nproperties
        return n

    def _read_pfast_nx(self, data: bytes, n: int) -> int:
        """
        PFAST(3601,36,55)
        NX only
        """
        self.to_nx()
        ntotal = 48
        struct1 = Struct(self._endian + b'ifii 8f')
        nproperties = (len(data) - n) // ntotal
        delta = (len(data) - n) % ntotal
        assert delta == 0, 'len(data)-n=%s n=%s' % (len(data) - n, (len(data) - n) / 48.)
        for unused_i in range(nproperties):
            edata = data[n:n+ntotal]
            out = struct1.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  PFAST=%s\n' % str(out))
            (pid, d, mcid, mflag, kt1, kt2, kt3, kr1, kr2, kr3, mass, ge) = out

            data_in = (pid, d, mcid, mflag, kt1, kt2, kt3,
                       kr1, kr2, kr3, mass, ge)
            prop = PFAST.add_op2_data(data_in)
            self._add_op2_property(prop)
            n += ntotal
        self.card_count['PFAST'] = nproperties
        return n

    def _read_pelast(self, data: bytes, n: int) -> int:
        """
        Record 41 -- PELAST(1302,13,34)

        1 PID   I Property identification number
        2 TKID  I TABLEDi entry identification number for stiffness
        3 TGEID I TABLEDi entry identification number for structural
                  damping
        4 TKNID I TABLEDi entry
        """
        ntotal = 16 * self.factor
        struct_4i = Struct(mapfmt(self._endian + b'4i', self.size))
        nproperties = (len(data) - n) // ntotal
        for unused_i in range(nproperties):
            edata = data[n:n+ntotal]
            out = struct_4i.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  PELAST=%s\n' % str(out))
            #(pid, tkid, tgeid, tknid) = out
            prop = PELAST.add_op2_data(out)
            self._add_pelast_object(prop)
            n += ntotal
        self.card_count['PELAST'] = nproperties
        return n

    def _read_pgap(self, data: bytes, n: int) -> int:
        """
        PGAP(2102,21,121) - the marker for Record 42
        """
        ntotal = 44 * self.factor
        struct_i10f = Struct(mapfmt(self._endian + b'i10f', self.size))
        nproperties = (len(data) - n) // ntotal
        for unused_i in range(nproperties):
            edata = data[n:n+ntotal]
            out = struct_i10f.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  PGAP=%s\n' % str(out))
            #(pid,u0,f0,ka,kb,kt,mu1,mu2,tmax,mar,trmin) = out
            prop = PGAP.add_op2_data(out)
            self._add_op2_property(prop)
            n += ntotal
        self.card_count['PGAP'] = nproperties
        return n

    def _read_phbdy(self, data: bytes, n: int) -> int:
        """
        PHBDY(2802,28,236) - the marker for Record 43
        """
        struct_i3f = Struct(self._endian + b'ifff')
        nproperties = (len(data) - n) // 16
        for unused_i in range(nproperties):
            edata = data[n:n+16]
            out = struct_i3f.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  PHBDY=%s\n' % str(out))
            #(pid, af, d1, d2) = out
            prop = PHBDY.add_op2_data(out)
            self._add_phbdy_object(prop)
            n += 16
        self.card_count['PHBDY'] = nproperties
        return n

    def _read_pintc(self, data: bytes, n: int) -> int:
        self.log.info('skipping PINTC in EPT')
        return len(data)

    def _read_pints(self, data: bytes, n: int) -> int:
        self.log.info('skipping PINTS in EPT')
        return len(data)

    def _read_plplane(self, data: bytes, n: int) -> int:
        """
        PLPLANE(4606,46,375)

        NX 10
        1 PID     I Property identification number
        2 MID     I Material identification number
        3 CID     I Coordinate system identification number
        4 STR CHAR4 Location of stress and strain output
        5 T      RS Default membrane thickness for Ti on the connection entry
        6 CSOPT  I  Reserved for coordinate system definition of plane
        7 UNDEF(5) None

        MSC 2016
        PID       I Property identification number
        2 MID     I Material identification number
        3 CID     I Coordinate system identification number
        4 STR CHAR4 Location of stress and strain output
        5 UNDEF(7 ) none Not used

        .. warning:: CSOPT ad T are not supported
        """
        ntotal = 44 * self.factor  # 4*11
        if self.size == 4:
            s = Struct(self._endian + b'3i 4s f 6i')
        else:
            s = Struct(self._endian + b'3q 8s d 6q')
        nentries = (len(data) - n) // ntotal
        for unused_i in range(nentries):
            out = s.unpack(data[n:n+ntotal])
            pid, mid, cid, location, unused_t, unused_csopt = out[:6]
            location = location.decode('latin1')
            #self.show_data(data[n:n+ntotal], 'ifs')
            self.add_plplane(pid, mid, cid=cid, stress_strain_output_location=location)
            n += ntotal
        self.card_count['PLPLANE'] = nentries
        return n

    def _read_plsolid(self, data: bytes, n: int) -> int:
        """
        MSC 2016
        1 PID I Property identification number
        2 MID I Material identification number
        3 STR CHAR4 Location of stress and strain output
        4 UNDEF(4 ) none Not used

        NX 10
        1 PID I Property identification number
        2 MID I Material identification number
        3 STR CHAR4 Location of stress and strain output
        4 CSOPT I Reserved for coordinate system definition of plane
        5 UNDEF(3) None

        .. warning:: CSOPT is not supported
        """
        ntotal = 28 * self.factor  # 4*7
        if self.size == 4:
            struct1 = Struct(self._endian + b'2i 4s 4i')
        else:
            struct1 = Struct(self._endian + b'2q 8s 4q')
        nentries = (len(data) - n) // ntotal
        for unused_i in range(nentries):
            out = struct1.unpack(data[n:n+ntotal])
            pid, mid, location, unused_csopt, unused_null_a, unused_null_b, unused_null_c = out
            location = location.decode('latin1')
            #self.show_data(data[n:n+ntotal], 'ifs')
            self.add_plsolid(pid, mid, stress_strain=location, ge=0.)
            n += ntotal
        self.card_count['PLSOLID'] = nentries
        return n

    def _read_pmass(self, data: bytes, n: int) -> int:
        """
        PMASS(402,4,44) - the marker for Record 48
        """
        ntotal = 8 * self.factor # 2*4
        struct_if = Struct(mapfmt(self._endian + b'if', self.size))
        nentries = (len(data) - n) // ntotal
        for unused_i in range(nentries):
            edata = data[n:n + ntotal]
            out = struct_if.unpack(edata)
            #out = (pid, mass)
            if self.is_debug_file:
                self.binary_debug.write('  PMASS=%s\n' % str(out))
            prop = PMASS.add_op2_data(out)
            self._add_op2_property_mass(prop)
            n += ntotal
        return n

    def _read_prod(self, data: bytes, n: int) -> int:
        """
        PROD(902,9,29) - the marker for Record 49
        """
        ntotal = 24 * self.factor  # 6*4
        struct_2i4f = Struct(mapfmt(self._endian + b'2i4f', self.size))
        nproperties = (len(data) - n) // ntotal
        for unused_i in range(nproperties):
            edata = data[n:n+ntotal]
            out = struct_2i4f.unpack(edata)
            #(pid, mid, a, j, c, nsm) = out
            prop = PROD.add_op2_data(out)
            if self.is_debug_file:
                self.binary_debug.write('  PROD=%s\n' % str(out))
            self._add_op2_property(prop)
            n += ntotal
        self.card_count['PROD'] = nproperties
        return n

    def _read_pshear(self, data: bytes, n: int) -> int:
        """
        PSHEAR(1002,10,42) - the marker for Record 50
        """
        ntotal = 24 * self.factor
        struct_2i4f = Struct(mapfmt(self._endian + b'2i4f', self.size))
        nproperties = (len(data) - n) // ntotal
        for unused_i in range(nproperties):
            edata = data[n:n+ntotal]
            out = struct_2i4f.unpack(edata)
            #(pid, mid, t, nsm, f1, f2) = out
            if self.is_debug_file:
                self.binary_debug.write('  PSHEAR=%s\n' % str(out))
            prop = PSHEAR.add_op2_data(out)
            self._add_op2_property(prop)
            n += ntotal
        self.card_count['PSHEAR'] = nproperties
        return n

    def _read_pshell(self, data: bytes, n: int) -> int:
        """
        PSHELL(2302,23,283) - the marker for Record 51
        """
        ntotal = 44 * self.factor  # 11*4
        s = Struct(mapfmt(self._endian + b'iififi4fi', self.size))
        nproperties = (len(data) - n) // ntotal
        for unused_i in range(nproperties):
            edata = data[n:n+ntotal]
            out = s.unpack(edata)
            (pid, mid1, unused_t, mid2, unused_bk, mid3, unused_ts,
             unused_nsm, unused_z1, unused_z2, mid4) = out
            if self.is_debug_file:
                self.binary_debug.write('  PSHELL=%s\n' % str(out))
            prop = PSHELL.add_op2_data(out)
            n += ntotal

            if pid in self.properties:
                # this is a fake PSHELL
                propi = self.properties[pid]
                if prop == propi:
                    self.log.warning('Fake PSHELL (skipping):\n%s' % propi)
                    nproperties -= 1
                    continue
                assert propi.type in ['PCOMP', 'PCOMPG'], propi.get_stats()
                self.log.warning(f'PSHELL is also PCOMP (skipping PSHELL):\n{propi}{prop}')
                nproperties -= 1
                continue

            if max(pid, mid1, mid2, mid3, mid4) > 1e8:
                self.big_properties[pid] = prop
            else:
                self._add_op2_property(prop)
        if nproperties:
            self.card_count['PSHELL'] = nproperties
        return n

    def _read_psolid(self, data: bytes, n: int) -> int:
        """
        PSOLID(2402,24,281) - the marker for Record 52
        """
        #print("reading PSOLID")
        if self.size == 4:
            ntotal = 28  # 7*4
            struct_6i4s = Struct(self._endian + b'6i4s')
        else:
            ntotal = 28 * 2
            struct_6i4s = Struct(self._endian + b'6q8s')

        nproperties = (len(data) - n) // ntotal
        for unused_i in range(nproperties):
            edata = data[n:n+ntotal]
            out = struct_6i4s.unpack(edata)
            #(pid, mid, cid, inp, stress, isop, fctn) = out
            #data_in = [pid, mid, cid, inp, stress, isop, fctn]
            if self.is_debug_file:
                self.binary_debug.write('  PSOLID=%s\n' % str(out))
            prop = PSOLID.add_op2_data(out)
            self._add_op2_property(prop)
            n += ntotal
        self.card_count['PSOLID'] = nproperties
        return n

# PSOLIDL
# PTRIA6
# PTSHELL

    def _read_ptube(self, data: bytes, n: int) -> int:
        """
        PTUBE(1602,16,30) - the marker for Record 56

        .. todo:: OD2 only exists for heat transfer...
                  how do i know if there's heat transfer at this point?
                  I could store all the tubes and add them later,
                  but what about themal/non-thermal subcases?

        .. warning:: assuming OD2 is not written (only done for thermal)
        """
        struct_2i3f = Struct(self._endian + b'2i3f')
        nproperties = (len(data) - n) // 20
        for unused_i in range(nproperties):
            edata = data[n:n+20]  # or 24???
            out = struct_2i3f.unpack(edata)
            (pid, mid, OD, t, nsm) = out
            data_in = [pid, mid, OD, t, nsm]
            if self.is_debug_file:
                self.binary_debug.write('  PTUBE=%s\n' % str(out))
            prop = PTUBE.add_op2_data(data_in)
            self._add_op2_property(prop)
            n += 20
        self.card_count['PTUBE'] = nproperties
        return n

    def _read_pset(self, data: bytes, n: int) -> int:
        struct_5i4si = Struct(self._endian + b'5i4si')
        nentries = 0
        while  n < len(data):
            edata = data[n:n+28]
            out = struct_5i4si.unpack(edata)
            #print(out)
            idi, poly1, poly2, poly3, cid, typei, typeid = out
            typei = typei.rstrip().decode('latin1')
            assert typei in ['SET', 'ELID'], (idi, poly1, poly2, poly3, cid, typei, typeid)
            if self.is_debug_file:
                self.binary_debug.write('  PVAL=%s\n' % str(out))
            #print(idi, poly1, poly2, poly3, cid, typei, typeid)
            typeids = []
            n += 28
            while typeid != -1:
                typeids.append(typeid)
                typeid, = self.struct_i.unpack(data[n:n+4])
                n += 4
                #print(val)
            #print(typeids)
            # PSET ID POLY1 POLY2 POLY3 CID SETTYP ID
            self.add_pset(idi, poly1, poly2, poly3, cid, typei, typeids)
        self.card_count['PSET'] = nentries
        return n

    def _read_pval(self, data: bytes, n: int) -> int:
        """
        PVAL(10201,102,400)

        Word Name Type Description
        1 ID       I p-value set identification number
        2 POLY1    I Polynomial order in 1 direction of the CID system
        3 POLY2    I Polynomial order in 2 direction of the CID system
        4 POLY3    I Polynomial order in 2 direction of the CID system
        5 CID      I Coordinate system identification number
        6 TYPE CHAR4 Type of set provided: "SET" or "ELID"
        7 TYPEID   I SET identification number or element identification
                     number with this p-value specification.
        Words 1 through 7 repeat until End of Record
        """
        #self.show_data(data[n:])
        if self.size == 4:
            struct_5i4si = Struct(self._endian + b'5i 4s i')
            struct_i = self.struct_i
        else:
            struct_5i4si = Struct(self._endian + b'5q 8s q')
            struct_i = self.struct_q

        nentries = 0
        ntotal = 28 * self.factor
        size = self.size
        while  n < len(data):
            edata = data[n:n+ntotal]
            out = struct_5i4si.unpack(edata)
            #print(out)
            idi, poly1, poly2, poly3, cid, typei, typeid = out
            typei = typei.rstrip().decode('latin1')
            assert typei in ['SET', 'ELID'], f'idi={idi} poly1={poly1} poly2={poly2} poly3={poly3} cid={cid} typei={typei} typeid={typeid}'
            if self.is_debug_file:
                self.binary_debug.write('  PVAL=%s\n' % str(out))
            #print(idi, poly1, poly2, poly3, cid, typei, typeid)
            typeids = []
            n += ntotal
            while typeid != -1:
                typeids.append(typeid)
                typeid, = struct_i.unpack(data[n:n+size])
                n += size
                #print(val)
            #print(typeids)
            # PVAL ID POLY1 POLY2 POLY3 CID SETTYP ID
            self.add_pval(idi, poly1, poly2, poly3, cid, typei, typeids)
        self.card_count['PVAL'] = nentries
        return n

    def _read_pvisc(self, data: bytes, n: int) -> int:
        """PVISC(1802,18,31) - the marker for Record 39"""
        struct_i2f = Struct(self._endian + b'i2f')
        nproperties = (len(data) - n) // 12
        for unused_i in range(nproperties):
            edata = data[n:n+12]
            out = struct_i2f.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  PVISC=%s\n' % str(out))
            #(pid, ce, cr) = out
            prop = PVISC.add_op2_data(out)
            self._add_op2_property(prop)
            n += 12
        self.card_count['PVISC'] = nproperties
        return n

# PWELD
# PWSEAM
    def _read_view(self, data: bytes, n: int) -> int:
        self.log.info('skipping VIEW in EPT')
        return len(data)

    def _read_view3d(self, data: bytes, n: int) -> int:
        self.log.info('skipping VIEW3D in EPT')
        return len(data)

def break_by_minus1(idata):
    """helper for ``read_nsm_nx``"""
    i1 = 0
    i = 0
    i2 = None
    packs = []
    for idatai in idata:
        #print('data[i:] = ', data[i:])
        if idatai == -1:
            i2 = i
            packs.append((i1, i2))
            i1 = i2 + 1
            i += 1
            continue
        i += 1
    #print(packs)
    return packs
