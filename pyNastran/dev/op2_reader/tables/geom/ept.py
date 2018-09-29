"""
defines readers for BDF objects in the OP2 EPT/EPTS table
"""
#pylint: disable=C0301,W0612,C0111,R0201,C0103,W0613,R0914
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
from struct import unpack, Struct
from six import b

import numpy as np

from pyNastran import is_release
from pyNastran.bdf.cards.properties.mass import PMASS, NSM, NSML
from pyNastran.bdf.cards.properties.bars import PBAR, PBARL, PBEND
from pyNastran.bdf.cards.properties.beam import PBEAM, PBEAML, PBCOMP
from pyNastran.bdf.cards.properties.bush import PBUSH
from pyNastran.bdf.cards.properties.damper import PDAMP, PVISC
from pyNastran.bdf.cards.properties.properties import PFAST, PGAP
from pyNastran.bdf.cards.properties.rods import PROD, PTUBE
from pyNastran.bdf.cards.properties.shell import PSHEAR, PSHELL, PCOMP
from pyNastran.bdf.cards.properties.solid import PSOLID
from pyNastran.bdf.cards.properties.springs import PELAS, PELAST

from pyNastran.bdf.cards.thermal.thermal import PCONV, PHBDY
# PCOMPG, PBUSH1D, PBEAML, PBEAM3
from pyNastran.op2.tables.geom.geom_common import GeomCommon


class EPT(GeomCommon):
    """defines methods for reading op2 properties"""

    def _read_ept_4(self, data, ndata):
        return self._read_geom_4(self._ept_map, data, ndata)

    def __init__(self):
        GeomCommon.__init__(self)
        self.big_properties = {}
        self._ept_map = {
            (3201, 32, 55): ['NSM', self._read_nsm],          # record 2  - needs an object holder (e.g. self.elements/self.properties)
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
            (3501, 35, 994) : ['NSM?', self._read_fake],
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
            (10301, 103, 399): ['PSET', self._read_fake],  # record 57
            (3002, 30, 415): ['VIEW3D', self._read_fake],  # record 63

            (13501, 135, 510) : ['PFAST', self._read_pfast_msc],  # MSC-specific

            # NX-specific
            (3601, 36, 55) : ['PFAST', self._read_pfast_nx],  # NX-specific
            (3801, 38, 979) : ['PPLANE', self._read_fake],
            (11801, 118, 560) : ['PWELD', self._read_fake],
            (3401, 34, 993) : ['NSMADD', self._read_fake],
        }

    def _add_op2_property(self, prop):
        #if prop.pid > 100000000:
            #raise RuntimeError('bad parsing; pid > 100000000...%s' % str(prop))
        self._add_property_object(prop, allow_overwrites=True)
        #print(str(prop)[:-1])

    def _add_pconv(self, prop):
        if prop.pconid > 100000000:
            raise RuntimeError('bad parsing pconid > 100000000...%s' % str(prop))
        self._add_convection_property_object(prop)

# HGSUPPR

    def _read_nsm(self, data, n):
        """NSM"""
        n = self._read_dual_card(data, n, self._read_nsm_nx, self._read_nsm_msc,
                                 'NSM', self._add_nsm_object)
        return n

    def _read_nsm_msc(self, data, n):
        """
        NSM(3201,32,55) - the marker for Record 2

        MSC
        1 SID       I Set identification number
        2 PROP  CHAR4 Set of property or elements
        3 ID        I Property or element identification number
        4 VALUE    RS Nonstructural mass value
        ORIGIN =0 NSM Bulk Data entry
        5 ID I Property or element ID
        6 VALUE RS Nonstructural mass value
        Words 5 through 6 repeat until End of Record
        ORIGIN =2 NSML Bulk Data entry
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

    def _read_nsm_nx(self, data, n):
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
        idata = np.frombuffer(data[n:], self.idtype).copy()
        fdata = np.frombuffer(data[n:], self.fdtype).copy()
        packs = break_by_minus1(idata)
        #for pack in packs:
            #print(pack)

        ipack = 0
        while n < len(data):
            #print('ints[i:]=', ints[i:].tolist())
            i1, i2 = packs[ipack]
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
            ipack += 1
        return n, properties

# NSM1
# NSML1
# NSMADD
# NSML
# NSML1
# PAABSF
# PACABS
# PACBAR

    def _read_pbar(self, data, n):
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
        ntotal = 76  # 19*4
        struct1 = Struct(self._endian + b'2i17f')
        nentries = (len(data) - n) // ntotal
        for i in range(nentries):
            edata = data[n:n+76]
            out = struct1.unpack(edata)
            (pid, mid, a, I1, I2, J, nsm, fe, c1, c2, d1, d2,
             e1, e2, f1, f2, k1, k2, I12) = out
            prop = PBAR.add_op2_data(out)
            self._add_op2_property(prop)
            n += ntotal
        self.card_count['PBAR'] = nentries
        return n

    def _read_pbarl(self, data, n):
        """
        PBARL(9102,91,52) - the marker for Record 12
        TODO: buggy
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

        ntotal = 28  # 7*4 - ROD - shortest entry...could be buggy... # TODO fix this
        struct1 = Struct(self._endian + b'2i8s8sf')
        #nentries = (len(data) - n) // ntotal
        #print(self.show_ndata(80))
        ndata = len(data)
        while ndata - n > ntotal:
            edata = data[n:n+28]
            n += 28

            out = struct1.unpack(edata)
            (pid, mid, group, beam_type, value) = out
            beam_type = beam_type.strip().decode('latin1')
            group = group.strip().decode('latin1')
            data_in = [pid, mid, group, beam_type, value]
            #self.log.debug("  pid=%s mid=%s group=%r beam_type=%r value=%s" % (
                #pid, mid, group, beam_type, value))
            if pid > 100000000:
                raise RuntimeError('bad parsing...')
            expected_length = valid_types[beam_type]
            iformat = b(self._uendian + '%if' % expected_length)

            ndelta = expected_length * 4
            data_in += list(unpack(iformat, data[n:n+ndelta]))
            # TODO why do i need the +4???
            #min_len =  expected_length * 4 + 4
            #if len(data)
            #data = data[n + expected_length * 4 + 4:]
            n += ndelta

            #prin( "len(out) = ",len(out)))
            #print("PBARL = %s" % data_in)
            prop = PBARL.add_op2_data(data_in)  # last value is nsm
            self._add_op2_property(prop)
            #print(self.show_data(data[n-8:-100]))
            break
        self.increase_card_count('PBARL')
        #assert len(data) == n
        return n

    def _read_pbcomp(self, data, n):
        struct1 = Struct(self._endian + b'2i 12f i')
        struct2 = Struct(self._endian + b'3f 2i')
        nproperties = 0
        ndata = len(data)
        while n < ndata:
            edata = data[n:n+60]  # 4*15
            n += 60
            data1 = struct1.unpack(edata)
            nsections = data1[-1]
            if self.is_debug_file:
                (pid, mid, a, i1, i2, i12, j, nsm, k1, k2, m1, m2, n1, n2, nsections) = data1
                self.log.info('PBCOMP pid=%s mid=%s nsections=%s\n' % (pid, mid, nsections))

            data2 = []
            if nsections in [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20]:
                # 16 Y   RS    Lumped area location along element's y-axis
                # 17 Z   RS    Lumped area location along element's z-axis
                # 18 C   RS    Fraction of the total area for the lumped area
                # 19 MID I     Material identification number
                # 20     UNDEF None
                # Words 16 through 20 repeat NSECT times
                for i in range(nsections):
                    datai = data[n:n+20]
                    xi, yi, ci, mid, null = struct2.unpack(datai)
                    data2.append((xi, yi, ci, mid))
                    n += 20
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
            self._add_op2_property(prop)
            nproperties += 1
        assert nproperties > 0, 'PBCOMP nproperties=%s' % (nproperties)
        self.card_count['PBCOMP'] = nproperties
        return n

    def _read_pbeam(self, data, n):
        """
        PBEAM(5402,54,262) - the marker for Record 14
        .. todo:: add object
        """
        struct1 = Struct(self._endian + b'4if')
        struct2 = Struct(self._endian + b'16f')
        struct3 = Struct(self._endian + b'16f')
        ntotal = 1072  # 44+12*84+20
        nproperties = (len(data) - n) // ntotal
        #assert nproperties > 0, 'ndata-n=%s n=%s datai\n%s' % (len(data)-n, n, self.show_data(data[n:100+n]))
        ndata = len(data)
        #self.show_data(data[12:], 'if')
        #assert ndata % ntotal == 0, 'ndata-n=%s n=%s ndata%%ntotal=%s' % (len(data)-n, n, ndata % ntotal)
        while n < ndata:
        #while 1: #for i in range(nproperties):
            edata = data[n:n+20]
            n += 20
            data_in = list(struct1.unpack(edata))
            #if self.is_debug_file:
                #self.log.info('PBEAM pid=%s mid=%s nsegments=%s ccf=%s x=%s\n' % tuple(data_in))
            (pid, mid, nsegments, ccf, x) = data_in
            #self.log.info('PBEAM pid=%s mid=%s nsegments=%s ccf=%s x=%s' % tuple(data_in))

            # Constant cross-section flag: 1=yes and 0=no
            # what is 2?
            if ccf not in [0, 1, 2]:
                msg = ('  PBEAM pid=%s mid=%s nsegments=%s ccf=%s x=%s; '
                       'ccf must be in [0, 1, 2]\n' % tuple(data_in))
                raise ValueError(msg)

            for i in range(11):
                edata = data[n:n+64]
                if len(edata) != 64:
                    endpack = []
                    raise RuntimeError('PBEAM unexpected length i=%s...' % i)
                n += 64
                pack = struct2.unpack(edata)
                (soi, xxb, a, i1, i2, i12, j, nsm, c1, c2,
                 d1, d2, e1, e2, f1, f2) = pack

                if soi == 0.0:
                    so_str = 'NO'
                elif soi == 1.0:
                    so_str = 'YES'
                else:
                    so_str = str(soi)
                    #msg = 'PBEAM pid=%s i=%s x/xb=%s soi=%s; soi not in 0.0 or 1.0' % (pid, i, xxb, soi)
                    #raise NotImplementedError(msg)

                #if xxb != 0.0:
                    #msg = 'PBEAM pid=%s i=%s x/xb=%s soi=%s; xxb not in 0.0 or 1.0' % (pid, i, xxb, soi)
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
                    self.log.debug(msg)
            edata = data[n:n+64]
            if len(edata) != 64:
                endpack = []
                raise RuntimeError('PBEAM unexpected length 2...')
                #break
            else:
                endpack = struct3.unpack(edata)
                n += 64

            assert len(endpack) == 16, endpack
            (k1, k2, s1, s2, nsia, nsib, cwa, cwb, # 8
             m1a, m2a, m1b, m2b, n1a, n2a, n1b, n2b) = endpack # 8 -> 16
            if self.is_debug_file:
                self.log.debug('    k=[%s,%s] s=[%s,%s] nsi=[%s,%s] cw=[%s,%s] '
                               'ma=[%s,%s] mb=[%s,%s] na=[%s,%s] nb=[%s,%s]' % (tuple(endpack)))
            data_in.append(endpack)

            if pid in self.properties:
                if self.properties[pid].type == 'PBCOMP':
                    continue
            prop = PBEAM.add_op2_data(data_in)
            self._add_op2_property(prop)
        self.card_count['PBEAM'] = nproperties
        return n

    def _read_pbeaml(self, data, n):
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
        ints = np.frombuffer(data[n:], self._uendian + 'i').copy()
        floats = np.frombuffer(data[n:], self._uendian + 'f').copy()
        iminus1 = np.where(ints == -1)[0]

        istart = [0] + list(iminus1[:-1] + 1)
        iend = iminus1

        struct1 = Struct(self._endian + b'2i8s8s')
        for i, (istarti, iendi) in enumerate(zip(istart, iend)):
            idata = data[n+istarti*4 : n+(istarti+6)*4]
            pid, mid, group, beam_type = struct1.unpack(idata)
            group = group.decode('latin1').strip()
            beam_type = beam_type.decode('latin1').strip()
            fvalues = floats[istarti+6: iendi]
            if self.is_debug_file:
                self.binary_debug.write('     %s\n' % str(fvalues))
                self.log.debug('pid=%i mid=%i group=%r beam_type=%r' % (pid, mid, group, beam_type))
                self.log.debug(fvalues)
            #self.log.debug('pid=%i mid=%i group=%s beam_type=%s' % (pid, mid, group, beam_type))
            data_in = [pid, mid, group, beam_type, fvalues]
            prop = PBEAML.add_op2_data(data_in)
            self._add_op2_property(prop)
        nproperties = len(istart)
        self.card_count['PBEAML'] = nproperties
        return len(data)

    def _read_pbend(self, data, n):
        """PBEND"""
        n = self._read_dual_card(data, n, self._read_pbend_nx, self._read_pbend_msc,
                                 'PBEND', self._add_property_object)
        return n

    def _read_pbend_msc(self, data, n):
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
        assert nproperties > 0, 'table=%r len=%s' % (self.table_name, len(data) - n)
        properties = []
        for i in range(nproperties):
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

    def _read_pbend_nx(self, data, n):
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
        assert nproperties > 0, 'table=%r len=%s' % (self.table_name, len(data) - n)
        properties = []
        for i in range(nproperties):
            edata = data[n:n+132]
            out = struct1.unpack(edata)
            (pid, mid, area, i1, i2, j, fsi, rm, t, p, rb, theta_b,
             c1, c2, d1, d2, e1, e2, f1, f2, k1, k2, nsm, rc, zc,
             delta_n, sacl, alpha, flange, kx, ky, kz, junk,) = out
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

    def _read_pbush(self, data, n):
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
        n = self._read_dual_card(data, n, self._read_pbush_nx, self._read_pbush_nx,
                                 'PBUSH', self._add_op2_property)
        return n

    def _read_pbush_nx(self, data, n):
        """PBUSH(1402,14,37) - 18 fields"""
        #if self.table_name == ['EPTS', 'EPT']:
        ntotal = 72
        struct1 = Struct(self._endian + b'i17f')
        nentries = (len(data) - n) // ntotal
        assert nentries > 0, 'table=%r len=%s' % (self.table_name, len(data) - n)
        props = []
        for i in range(nentries):
            edata = data[n:n+72]
            out = struct1.unpack(edata)
            (pid, k1, k2, k3, k4, k5, k6, b1, b2, b3, b4, b5, b6,
             g1, sa, st, ea, et) = out
            g2 = g3 = g4 = g5 = g6 = g1
            data_in = (pid, k1, k2, k3, k4, k5, k6, b1, b2, b3, b4, b5, b6,
                       g1, g2, g3, g4, g5, g6, sa, st, ea, et)
            prop = PBUSH.add_op2_data(data_in)
            props.append(prop)
            n += ntotal
        return n, props

    def _read_pbush_msc(self, data, n):
        """PBUSH(1402,14,37)"""
        ntotal = 92  # 23*4
        struct1 = Struct(self._endian + b'i22f')
        nentries = (len(data) - n) // ntotal
        assert nentries > 0, 'table=%r len=%s' % (self.table_name, len(data) - n)
        props = []
        for i in range(nentries):
            edata = data[n:n+92]
            out = struct1.unpack(edata)
            (pid, k1, k2, k3, k4, k5, k6, b1, b2, b3, b4, b5, b6,
             g1, g2, g3, g4, g5, g6, sa, st, ea, et) = out
            prop = PBUSH.add_op2_data(out)
            props.append(prop)
            n += ntotal
        return n, props

    def _read_pbush1d(self, data, n):
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
        24 IDTD   I  DEQATN entry identification number for compression
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
            0 : None,
            1 : 'EQUAT',
            2 : 'TABLE',
        }
        ntotal = 152  # 38*4
        struct1 = Struct(self._endian + b'i 6f i 4f 24i 2f')
        nentries = (len(data) - n) // ntotal
        for i in range(nentries):
            edata = data[n:n+152]
            out = struct1.unpack(edata)
            (pid, k, c, m, unused_alpha, sa, se,
             typea, cvt, cvc, expvt, expvc, idtsu, idtcu, idtsud, idcsud,
             types, idts, idcs, idtdu, idcdu, typed, idtd, idtd, idtdv, idcdv,
             typeg, idtg, idcg, idtdu, idcdu, idtdv, idcdv,
             typef, idtf, idcf,
             ut, uc) = out
            if not is_release:
                if typea in [1, 2]:
                    raise NotImplementedError(str((typea, cvt, cvc, expvt, expvc, idtsu, idtcu, idtsud, idcsud)))
                if types in [1, 2]:
                    raise NotImplementedError(str((types, idts, idcs, idtdu, idcdu, typed, idtd, idtd, idtdv, idcdv)))
                if typeg in [1, 2]:
                    raise NotImplementedError(str((typeg, idtg, idcg, idtdu, idcdu, idtdv, idcdv)))
                if typef in [1, 2]:
                    raise NotImplementedError(str((idtf, idcf)))
                typea_str = type_map[typea]
                types_str = type_map[types]
                typeg_str = type_map[typeg]
                typef_str = type_map[typef]
            self.add_pbush1d(pid, k=k, c=c, m=m, sa=sa, se=se,
                             optional_vars=None,)
            n += ntotal
        self.card_count['PBUSH1D'] = nentries
        return n

    def _read_pbusht(self, data, n):
        self.log.info('skipping PBUSHT in EPT')
        return len(data)

    def _read_pcomp(self, data, n):
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
        """
        nproperties = 0
        s1 = Struct(self._endian + b'2i3fi2f')
        s2 = Struct(self._endian + b'i2fi')

        ndata = len(data)
        while n < (ndata - 32):
            out = s1.unpack(data[n:n+32])
            (pid, nlayers, z0, nsm, sb, ft, Tref, ge) = out
            if self.binary_debug:
                self.log.debug('PCOMP pid=%s nlayers=%s z0=%s nsm=%s sb=%s ft=%s Tref=%s ge=%s' % tuple(out))
            assert isinstance(nlayers, int), out
            n += 32

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
                pid, nlayers, z0, nsm, sb, ft, Tref, ge)

            if self.is_debug_file:
                self.binary_debug.write('    pid=%s nlayers=%s z0=%s nms=%s sb=%s ft=%s Tref=%s ge=%s\n' % (
                    pid, nlayers, z0, nsm, sb, ft, Tref, ge))
            for ilayer in range(nlayers):
                (mid, t, theta, sout) = s2.unpack(data[n:n+16])
                mids.append(mid)
                T.append(t)
                thetas.append(theta)
                souts.append(sout)
                if self.is_debug_file:
                    self.binary_debug.write('      mid=%s t=%s theta=%s sout=%s\n' % (
                        mid, t, theta, sout))
                n += 16

            data_in = [
                pid, z0, nsm, sb, ft, Tref, ge,
                is_symmetrical, mids, T, thetas, souts]
            prop = PCOMP.add_op2_data(data_in)
            self._add_op2_property(prop)
            nproperties += 1
        self.card_count['PCOMP'] = nproperties
        return n

    def _read_pcompg(self, data, n):
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
        """
        nproperties = 0
        s1 = Struct(self._endian + b'2i 3f i 2f')
        s2 = Struct(self._endian + b'2i2fi')
        struct_i5 = Struct(self._endian + b'5i')

        # lam - SYM, MEM, BEND, SMEAR, SMCORE, None
        lam_map = {
            0 : None,
        }

        # ft - HILL, HOFF, TSAI, STRN, None
        ft_map = {
            0 : None,
        }
        # sout - YES, NO
        sout_map = {
            0 : 'NO',
        }
        ndata = len(data)
        while n < (ndata - 32):
            out = s1.unpack(data[n:n+32])
            (pid, lam_int, z0, nsm, sb, ft_int, tref, ge) = out
            if self.binary_debug:
                self.binary_debug.write('PCOMPG pid=%s lam_int=%s z0=%s nsm=%s sb=%s ft_int=%s tref=%s ge=%s' % tuple(out))
            assert isinstance(lam_int, int), out
            n += 32

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
                ints5 = struct_i5.unpack(data[n:n+20])
                if ints5 == (-1, -1, -1, -1, -1):
                    if self.is_debug_file:
                        self.binary_debug.write('      global_ply=%-1 mid=%-1 t=%-1 theta=%-1 sout=-1\n')
                    break
                (global_ply, mid, t, theta, sout_int) = s2.unpack(data[n:n+20])
                try:
                    sout = sout_map[sout_int]
                except KeyError:
                    self.log.error('cant parse iply=%s sout=%s; assuming 0=NO' % (iply, sout_int))
                    sout = 'NO'

                global_ply_ids.append(global_ply)
                mids.append(mid)
                thicknesses.append(t)
                thetas.append(theta)
                souts.append(sout)
                if self.is_debug_file:
                    self.binary_debug.write('      global_ply=%s mid=%s t=%s theta=%s sout_int=%s sout=%r\n' % (
                        global_ply, mid, t, theta, sout_int, sout))
                n += 20
                ilayer += 1

            try:
                ft = ft_map[ft_int]
            except KeyError:
                self.log.error('pid=%s cant parse ft=%s; should be HILL, HOFF, TSAI, STRN...skipping' % (pid, ft_int))
                continue

            try:
                lam = lam_map[lam_int]
            except KeyError:
                self.log.error('pid=%s cant parse lam=%s; should be HILL, HOFF, TSAI, STRN...skipping' % (pid, lam_int))
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

    def _read_pconeax(self, data, n):
        """
        (152,19,147) - Record 24
        """
        self.log.info('skipping PCONEAX in EPT')
        return len(data)

    def _read_pconv(self, data, n):
        """common method for reading PCONVs"""
        n = self._read_dual_card(data, n, self._read_pconv_nx, self._read_pconv_msc,
                                 'PCONV', self._add_pconv)
        return n

    def _read_pconv_nx(self, data, n):
        """
        (11001,110,411)- NX version
        """
        ntotal = 16  # 4*4
        struct_3if = Struct(self._endian + b'3if')
        nentries = (len(data) - n) // ntotal
        assert (len(data) - n) % ntotal == 0
        props = []
        for i in range(nentries):
            out = struct_3if.unpack(data[n:n+ntotal])
            (pconid, mid, form, expf) = out
            ftype = tid = chlen = gidin = ce = e1 = e2 = e3 = None
            data_in = (pconid, mid, form, expf, ftype, tid, chlen,
                       gidin, ce, e1, e2, e3)

            prop = PCONV.add_op2_data(data_in)
            props.append(prop)
            n += ntotal
        return n, props

    def _read_pconv_msc(self, data, n):
        """
        (11001,110,411)- MSC version - Record 25
        """
        ntotal = 56  # 14*4
        s = Struct(self._endian + b'3if 4i fii 3f')
        nentries = (len(data) - n) // ntotal
        assert (len(data) - n) % ntotal == 0
        props = []
        for i in range(nentries):
            out = s.unpack(data[n:n+ntotal])
            (pconid, mid, form, expf, ftype, tid, undef1, undef2, chlen,
             gidin, ce, e1, e2, e3) = out
            data_in = (pconid, mid, form, expf, ftype, tid, chlen,
                       gidin, ce, e1, e2, e3)

            prop = PCONV.add_op2_data(data_in)
            props.append(prop)
            n += ntotal
        return n, props

    def _read_pconvm(self, data, n):  # 26
        self.log.info('skipping PCONVM in EPT')
        return len(data)

    def _read_pdamp(self, data, n):
        """
        PDAMP(202,2,45) - the marker for Record ???
        """
        ntotal = 8  # 2*4
        struct_if = Struct(self._endian + b'if')
        nentries = (len(data) - n) // ntotal
        for i in range(nentries):
            out = struct_if.unpack(data[n:n+ntotal])
            #(pid, b) = out
            prop = PDAMP.add_op2_data(out)
            self._add_op2_property(prop)
            n += ntotal
        self.card_count['PDAMP'] = nentries
        return n

    def _read_pdampt(self, data, n):  # 26
        self.log.info('skipping PDAMPT in EPT')
        return len(data)

    def _read_pdamp5(self, data, n):  # 26
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

    def _read_pelas(self, data, n):
        """PELAS(302,3,46) - the marker for Record 39"""
        struct_i3f = Struct(self._endian + b'i3f')
        ntotal = 16  # 4*4
        nproperties = (len(data) - n) // ntotal
        for i in range(nproperties):
            edata = data[n:n+16]
            out = struct_i3f.unpack(edata)
            #(pid, k, ge, s) = out
            if self.is_debug_file:
                self.binary_debug.write('  PELAS=%s\n' % str(out))
            prop = PELAS.add_op2_data(out)
            self._add_op2_property(prop)
            n += ntotal
        self.card_count['PELAS'] = nproperties
        return n

    def _read_pfast_msc(self, data, n):
        ntotal = 92 # 23*4
        struct1 = Struct(self._endian + b'ifii 4f')
        nproperties = (len(data) - n) // ntotal
        for i in range(nproperties):
            edata = data[n:n+ntotal]
            out = struct1.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  PFAST=%s\n' % str(out))
            (pid, d, mcid, connbeh, conntype, extcon, condtype, weldtype,
             minlen, maxlen, gmcheck, spcgs, mass, ge, aa, bb, cc, mcid, mflag,
             kt1, kt2, kt3, kr1, kr2, kr3) = out

            data_in = (pid, d, mcid, mflag, kt1, kt2, kt3,
                       kr1, kr2, kr3, mass, ge)
            prop = PFAST.add_op2_data(data_in)
            self._add_op2_property(prop)
            n += ntotal
        self.card_count['PFAST'] = nproperties
        return n

    def _read_pfast_nx(self, data, n):
        """
        PFAST(3601,36,55)
        NX only
        """
        ntotal = 48
        struct1 = Struct(self._endian + b'ifii 8f')
        nproperties = (len(data) - n) // ntotal
        delta = (len(data) - n) % ntotal
        assert delta == 0, 'len(data)-n=%s n=%s' % (len(data) - n, (len(data) - n) / 48.)
        for i in range(nproperties):
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

    def _read_pelast(self, data, n):
        """
        Record 41 -- PELAST(1302,13,34)

        1 PID   I Property identification number
        2 TKID  I TABLEDi entry identification number for stiffness
        3 TGEID I TABLEDi entry identification number for structural
                  damping
        4 TKNID I TABLEDi entry
        """
        ntotal = 16
        struct_4i = Struct(self._endian + b'4i')
        nproperties = (len(data) - n) // ntotal
        for i in range(nproperties):
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

    def _read_pgap(self, data, n):
        """
        PGAP(3201,32,55) - the marker for Record 42
        """
        ntotal = 44
        struct_i10f = Struct(self._endian + b'i10f')
        nproperties = (len(data) - n) // ntotal
        for i in range(nproperties):
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

    def _read_phbdy(self, data, n):
        """
        PHBDY(2802,28,236) - the marker for Record 43
        """
        struct_i3f = Struct(self._endian + b'ifff')
        nproperties = (len(data) - n) // 16
        for i in range(nproperties):
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

    def _read_pintc(self, data, n):
        self.log.info('skipping PINTC in EPT')
        return len(data)

    def _read_pints(self, data, n):
        self.log.info('skipping PINTS in EPT')
        return len(data)

    def _read_plplane(self, data, n):
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
        ntotal = 44  # 4*11
        s = Struct(self._endian + b'3i 4s f 6i')
        nentries = (len(data) - n) // ntotal
        for i in range(nentries):
            out = s.unpack(data[n:n+ntotal])
            pid, mid, cid, location, t, csopt = out[:6]
            location = location.decode('latin1')
            #self.show_data(data[n:n+ntotal], 'ifs')
            plplane = self.add_plplane(pid, mid, cid=cid, stress_strain_output_location=location)
            #print(plplane)
            n += ntotal
        self.card_count['PLPLANE'] = nentries
        return n

    def _read_plsolid(self, data, n):
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
        ntotal = 28  # 4*7
        struct1 = Struct(self._endian + b'2i 4s 4i')
        nentries = (len(data) - n) // ntotal
        for i in range(nentries):
            out = struct1.unpack(data[n:n+ntotal])
            pid, mid, location, csopt, null_a, null_b, null_c = out
            location = location.decode('latin1')
            #self.show_data(data[n:n+ntotal], 'ifs')
            plsolid = self.add_plsolid(pid, mid, stress_strain=location, ge=0.)
            n += ntotal
        self.card_count['PLSOLID'] = nentries
        return n

    def _read_pmass(self, data, n):
        """
        PMASS(402,4,44) - the marker for Record 48
        """
        struct_if = Struct(self._endian + b'if')
        nentries = (len(data) - n) // 8  # 2*4
        for i in range(nentries):
            edata = data[n:n + 8]
            out = struct_if.unpack(edata)
            #out = (pid, mass)
            if self.is_debug_file:
                self.binary_debug.write('  PMASS=%s\n' % str(out))
            prop = PMASS.add_op2_data(out)
            self._add_op2_property(prop)
            n += 8
        return n

    def _read_prod(self, data, n):
        """
        PROD(902,9,29) - the marker for Record 49
        """
        ntotal = 24  # 6*4
        struct_2i4f = Struct(self._endian + b'2i4f')
        nproperties = (len(data) - n) // ntotal
        for i in range(nproperties):
            edata = data[n:n+24]
            out = struct_2i4f.unpack(edata)
            #(pid, mid, a, j, c, nsm) = out
            prop = PROD.add_op2_data(out)
            if self.is_debug_file:
                self.binary_debug.write('  PROD=%s\n' % str(out))
            self._add_op2_property(prop)
            n += ntotal
        self.card_count['PROD'] = nproperties
        return n

    def _read_pshear(self, data, n):
        """
        PSHEAR(1002,10,42) - the marker for Record 50
        """
        struct_2i4f = Struct(self._endian + b'2i4f')
        nproperties = (len(data) - n) // 24
        for i in range(nproperties):
            edata = data[n:n+24]
            out = struct_2i4f.unpack(edata)
            #(pid, mid, t, nsm, f1, f2) = out
            if self.is_debug_file:
                self.binary_debug.write('  PSHEAR=%s\n' % str(out))
            prop = PSHEAR.add_op2_data(out)
            self._add_op2_property(prop)
            n += 24
        self.card_count['PSHEAR'] = nproperties
        return n

    def _read_pshell(self, data, n):
        """
        PSHELL(2302,23,283) - the marker for Record 51
        """
        ntotal = 44  # 11*4
        s = Struct(self._endian + b'iififi4fi')
        nproperties = (len(data) - n) // ntotal
        for i in range(nproperties):
            edata = data[n:n+44]
            out = s.unpack(edata)
            (pid, mid1, t, mid2, bk, mid3, ts, nsm, z1, z2, mid4) = out
            if self.is_debug_file:
                self.binary_debug.write('  PSHELL=%s\n' % str(out))
            prop = PSHELL.add_op2_data(out)

            if max(pid, mid1, mid2, mid3, mid4) > 1e8:
                #print("PSHELL = ",out)
                self.big_properties[pid] = prop
            else:
                self._add_op2_property(prop)
            n += ntotal
        self.card_count['PSHELL'] = nproperties
        return n

    def _read_psolid(self, data, n):
        """
        PSOLID(2402,24,281) - the marker for Record 52
        """
        #print("reading PSOLID")
        ntotal = 28  # 7*4
        struct_6i4s = Struct(self._endian + b'6i4s')
        nproperties = (len(data) - n) // ntotal
        for i in range(nproperties):
            edata = data[n:n+28]
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

    def _read_ptube(self, data, n):
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
        for i in range(nproperties):
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

    def _read_pset(self, data, n):
        self.log.info('skipping PSET in EPT')
        return len(data)

    def _read_pval(self, data, n):
        self.log.info('skipping PVAL in EPT')
        return len(data)

    def _read_pvisc(self, data, n):
        """PVISC(1802,18,31) - the marker for Record 39"""
        struct_i2f = Struct(self._endian + b'i2f')
        nproperties = (len(data) - n) // 12
        for i in range(nproperties):
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
    def _read_view(self, data, n):
        self.log.info('skipping VIEW in EPT')
        return len(data)

    def _read_view3d(self, data, n):
        self.log.info('skipping VIEW3D in EPT')
        return len(data)
