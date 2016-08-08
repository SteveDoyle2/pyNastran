#pylint: disable=C0301,W0612,C0111,R0201,C0103,W0613,R0914
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
from struct import unpack, Struct
from six import b
from six.moves import range

from pyNastran.bdf.bdf import NSM, PMASS

from pyNastran.bdf.cards.properties.bars import PBAR, PBARL
from pyNastran.bdf.cards.properties.beam import PBEAM
from pyNastran.bdf.cards.properties.bush import PBUSH
from pyNastran.bdf.cards.properties.damper import PDAMP, PVISC
from pyNastran.bdf.cards.properties.properties import PFAST, PGAP
from pyNastran.bdf.cards.properties.rods import PROD, PTUBE
from pyNastran.bdf.cards.properties.shell import PSHEAR, PSHELL, PCOMP
from pyNastran.bdf.cards.properties.solid import PSOLID
from pyNastran.bdf.cards.properties.springs import PELAS

from pyNastran.bdf.cards.thermal.thermal import PCONV, PHBDY
# PCOMPG, PBUSH1D, PBEAML, PBEAM3
from pyNastran.op2.tables.geom.geom_common import GeomCommon


class EPT(GeomCommon):
    """defines methods for reading op2 properties"""

    def _read_ept_4(self, data, ndata):
        return self._read_geom_4(self._ept_map, data, ndata)

    def __init__(self):
        GeomCommon.__init__(self)
        self.bigProperties = {}
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
            (9202, 92, 53): ['PBEAML', self._read_pbeaml],     # record 15 - not done
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
            (15006, 150, 604): ['PCOMPG', self._read_fake],  # record

            (702, 7, 38): ['PBUSHT', self._read_fake],  # record 1
            (3301, 33, 56): ['NSM1', self._read_fake],  # record 3
            (3401, 34, 57) : ['NSMADD', self._read_fake],    # record 5
            (3501, 35, 58): ['NSML', self._read_fake],  # record 6
            (3601, 36, 62): ['NSML1', self._read_fake],  # record 7
            (1502, 15, 36): ['PAABSF', self._read_fake],  # record 8
            (8300, 83, 382): ['PACABS', self._read_fake],  # record 9
            (8500, 85, 384): ['PACBAR', self._read_fake],  # record 10
            (5403, 55, 349): ['PBCOMP', self._read_fake],  # record 13
            (13301, 133, 509): ['PBMSECT', self._read_fake],  # record 17
            (2902, 29, 420): ['PCONVM', self._read_fake],  # record 26
            (1202, 12, 33): ['PDAMPT', self._read_fake],  # record 28
            (8702, 87, 412): ['PDAMP5', self._read_fake],  # record 29
            (6802, 68, 164): ['PDUM8', self._read_fake],  # record 37
            (6902, 69, 165): ['PDUM9', self._read_fake],  # record 38
            (1302, 13, 34): ['PELAST', self._read_fake],  # record 41
            (12001, 120, 480): ['PINTC', self._read_fake],  # record 44
            (12101, 121, 484): ['PINTS', self._read_fake],  # record 45
            (4606, 46, 375): ['PLPLANE', self._read_fake],  # record 46
            (4706, 47, 376): ['PLSOLID', self._read_fake],  # record 47
            (10301, 103, 399): ['PSET', self._read_fake],  # record 57
            (3002, 30, 415): ['VIEW3D', self._read_fake],  # record 63

            (13501, 135, 510) : ['PFAST', self._read_pfast_msc],  # MSC-specific

            # NX-specific
            (3601, 36, 55) : ['PFAST', self._read_pfast_nx],  # NX-specific
        }

    def _add_op2_property(self, prop):
        if prop.pid > 100000000:
            raise RuntimeError('bad parsing...%s' % str(prop))
        self.add_property(prop, allow_overwrites=True)
        #print(str(prop)[:-1])

    def _add_pconv(self, prop):
        if prop.pconid > 100000000:
            raise RuntimeError('bad parsing...%s' % str(prop))
        self.add_convection_property(prop)

# HGSUPPR

    def _read_nsm(self, data, n):
        """
        NSM(3201,32,55) - the marker for Record 2
        .. todo:: this isnt a property...
        """
        self.log.debug('skipping NSM in EPT\n')
        return len(data)
        s = Struct(b(self._endian + 'i4sif'))
        while len(data) >= 16:  # 4*4
            edata = data[:16]
            data = data[16:]
            out = s.unpack(edata)
            (sid, prop_set, ID, value) = out
            #print("sid=%s propSet=%s ID=%s value=%s" %(sid,propSet,ID,value))
            prop = NSM.add_op2_data([sid, prop_set, ID, value])
            #self._add_op2_property(prop)
        return n

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
        """
        ntotal = 76  # 19*4
        s = Struct(b(self._endian + '2i17f'))
        nentries = (len(data) - n) // ntotal
        for i in range(nentries):
            edata = data[n:n+76]
            out = s.unpack(edata)
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
        s = Struct(b(self._endian + '2i8s8sf'))
        #nentries = (len(data) - n) // ntotal
        #print(self.show_ndata(80))
        ndata = len(data)
        while ndata - n > ntotal:
            edata = data[n:n+28]
            n += 28

            out = s.unpack(edata)
            (pid, mid, group, Type, value) = out
            Type = Type.strip().decode('latin1')
            group = group.strip().decode('latin1')
            data_in = [pid, mid, group, Type, value]
            #print("pid=%s mid=%s group=%r Type=%r value=%s" % (
                #pid, mid, group, Type, value))
            if pid > 100000000:
                raise RuntimeError('bad parsing...')
            expected_length = valid_types[Type]
            iformat = b'%if' % expected_length

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
        self._increase_card_count('PBARL')
        #assert len(data) == n
        return n

# PBCOMP

    def _read_pbeam(self, data, n):
        """
        PBEAM(5402,54,262) - the marker for Record 14
        .. todo:: add object
        """
        struct1 = Struct(b(self._endian + '4if'))
        struct2 = Struct(b(self._endian + '16f'))
        struct3 = Struct(b(self._endian + '16f'))
        ntotal = 1072  # 44+12*84+20
        nproperties = (len(data) - n) // ntotal
        #assert nproperties > 0, 'ndata-n=%s n=%s datai\n%s' % (len(data)-n, n, self.show_data(data[n:100+n]))
        ndata = len(data)
        while n < ndata:
        #while 1: #for i in range(nproperties):
            edata = data[n:n+20]
            n += 20
            data_in = list(struct1.unpack(edata))
            if self.is_debug_file:
                self.log.info('PBEAM pid=%s mid=%s nsegments=%s ccf=%s x=%s\n' % tuple(data_in))
            (pid, mid, nsegments, ccf, x) = data_in
            #self.log.info('PBEAM pid=%s mid=%s nsegments=%s ccf=%s x=%s' % tuple(data_in))

            # Constant cross-section flag: 1=yes and 0=no
            # what is 2?
            assert ccf in [0, 1, 2], '  PBEAM pid=%s mid=%s nsegments=%s ccf=%s x=%s\n' % tuple(data_in)
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
                    raise NotImplementedError('PBEAM pid=%s i=%s x/xb=%s soi=%s' % (pid, i, xxb, soi))

                pack2 = (so_str, xxb, a, i1, i2, i12, j, nsm, c1, c2,
                         d1, d2, e1, e2, f1, f2)
                data_in.append(pack2)
                if self.is_debug_file:
                    self.binary_debug.write('     %s\n' % str(pack))
                    self.log.info('    i=%-2s' % i + ' so=%s xxb=%.1f a=%g i1=%g i2=%g i12=%g j=%g nsm=%g '
                                  'c=[%s,%s] d=[%s,%s] e=[%s,%s] f=[%s,%s]' % (tuple(pack2)))
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
             m1a, m2a, m1b, m2b, n1a, n2a, n1b, n2b ) = endpack # 8 -> 16
            self.log.info('    k=[%s,%s] s=[%s,%s] nsi=[%s,%s] cw=[%s,%s] '
                          'ma=[%s,%s] mb=[%s,%s] na=[%s,%s] nb=[%s,%s]' % (tuple(endpack)))
            data_in.append(endpack)

            prop = PBEAM.add_op2_data(data_in)
            self._add_op2_property(prop)
        self.card_count['PBEAM'] = nproperties
        return n

    def _read_pbeaml(self, data, n):
        self.log.debug('skipping PBEAML in EPT\n')
        if self.is_debug_file:
            self.binary_debug.write('skipping PBEAML in EPT\n')
        return len(data)

    def _read_pbend(self, data, n):
        self.log.debug('skipping PBEND in EPT\n')
        if self.is_debug_file:
            self.binary_debug.write('skipping PBEND in EPT\n')
        return len(data)

# PBMSECT
# PBRSECT

    def _read_pbush(self, data, n):
        """
        The PBUSH card is different between MSC and NX Nastran.
        The NX version has 23 fields.
        The MSC version has 24 fields.
        There's also an 18 field version from pre-MSC/NX 2001.
        """
        if self.is_nx:
            return self._read_pbush_nx(data, n)
        return self._read_pbush_msc(data, n)

    def _read_pbush_nx(self, data, n):
        """PBUSH(1402,14,37)"""
        #if self.table_name == ['EPTS', 'EPT']:
        ntotal = 72
        s = Struct(b(self._endian + 'i17f'))
        nentries = (len(data) - n) // ntotal
        assert nentries > 0, 'table=%r len=%s' % (self.table_name, len(data) - n)
        for i in range(nentries):
            edata = data[n:n+72]
            out = s.unpack(edata)
            (pid, k1, k2, k3, k4, k5, k6, b1, b2, b3, b4, b5, b6,
             g1, sa, st, ea, et) = out
            g2 = g3 = g4 = g5 = g6 = g1
            data_in = (pid, k1, k2, k3, k4, k5, k6, b1, b2, b3, b4, b5, b6,
                       g1, g2, g3, g4, g5, g6, sa, st, ea, et)
            prop = PBUSH.add_op2_data(data_in)
            self._add_op2_property(prop)
            n += ntotal
        #else:
            #ntotal = 92  # 23*4
            #s = Struct(b(self._endian + 'i22f'))
            #nentries = (len(data) - n) // ntotal
            #assert nentries > 0, 'table=%r len=%s' % (self.table_name, len(data) - n)
            #for i in range(nentries):
                #edata = data[n:n+92]
                #out = s.unpack(edata)
                #(pid, k1, k2, k3, k4, k5, k6, b1, b2, b3, b4, b5, b6,
                 #g1, g2, g3, g4, g5, g6, sa, st, ea, et) = out
                #prop = PBUSH.add_op2_data(out)
                #self._add_op2_property(prop)
                #n += ntotal
        self.card_count['PBUSH'] = nentries
        return n

    def _read_pbush_msc(self, data, n):
        """PBUSH(1402,14,37)"""
        ntotal = 92  # 23*4
        s = Struct(b(self._endian + 'i22f'))
        nentries = (len(data) - n) // ntotal
        assert nentries > 0, 'table=%r len=%s' % (self.table_name, len(data) - n)
        for i in range(nentries):
            edata = data[n:n+92]
            out = s.unpack(edata)
            (pid, k1, k2, k3, k4, k5, k6, b1, b2, b3, b4, b5, b6,
             g1, g2, g3, g4, g5, g6, sa, st, ea, et) = out
            prop = PBUSH.add_op2_data(out)
            self._add_op2_property(prop)
            n += ntotal
        self.card_count['PBUSH'] = nentries
        return n

    def _read_pbush1d(self, data, n):
        self.log.debug('skipping PBUSH1D in EPT\n')
        if self.is_debug_file:
            self.binary_debug.write('skipping PBUSH1D in EPT\n')
        return len(data)

    def _read_pbusht(self, data, n):
        self.log.debug('skipping PBUSHT in EPT\n')
        if self.is_debug_file:
            self.binary_debug.write('skipping PBUSHT in EPT\n')
        return len(data)

    def _read_pcomp(self, data, n):
        """
        PCOMP(2706,27,287) - the marker for Record 22
        """
        nproperties = 0
        s1 = Struct(b(self._endian + '2i3fi2f'))
        s2 = Struct(b(self._endian + 'i2fi'))

        ndata = len(data)
        while n < (ndata - 32):
            out = s1.unpack(data[n:n+32])
            (pid, nlayers, z0, nsm, sb, ft, Tref, ge) = out
            if self.binary_debug:
                self.log.debug('PCOMP pid=%s nlayers=%s z0=%s nsm=%s sb=%s ft=%s Tref=%s ge=%s' % tuple(out))
            assert isinstance(nlayers, int), out
            n += 32

            Mid = []
            T = []
            Theta = []
            Sout = []

            # None, 'SYM', 'MEM', 'BEND', 'SMEAR', 'SMCORE', 'NO'
            is_symmetrical = 'NO'
            if nlayers < 0:
                is_symmetrical = 'SYM'
                nlayers = abs(nlayers)
            assert nlayers > 0, out

            assert 0 < nlayers < 100, 'pid=%s nlayers=%s z0=%s nms=%s sb=%s ft=%s Tref=%s ge=%s' % (
                pid, nlayers, z0, nsm, sb, ft, Tref, ge)
            for ilayer in range(nlayers):
                (mid, t, theta, sout) = s2.unpack(data[n:n+16])
                Mid.append(mid)
                T.append(t)
                Theta.append(theta)
                Sout.append(sout)
                if self.is_debug_file:
                    self.binary_debug.write('  mid=%s t=%s theta=%s sout=%s' % (mid, t, theta, sout))
                n += 16

            data_in = [
                pid, z0, nsm, sb, ft, Tref, ge,
                is_symmetrical, Mid, T, Theta, Sout]
            prop = PCOMP.add_op2_data(data_in)
            self._add_op2_property(prop)
            nproperties += 1
        self.card_count['PCOMP'] = nproperties
        return n

# PCOMPA
    def _read_pconeax(self, data, n):
        """
        (152,19,147) - Record 24
        """
        self.log.debug('skipping PCONEAX in EPT\n')
        if self.is_debug_file:
            self.binary_debug.write('skipping PCONEAX\n')
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
        s = Struct(b(self._endian + '3if'))
        nentries = (len(data) - n) // ntotal
        assert (len(data) - n) % ntotal == 0
        props = []
        for i in range(nentries):
            out = s.unpack(data[n:n+ntotal])
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
        s = Struct(b(self._endian + '3if 4i fii 3f'))
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
        self.log.debug('skipping PCONVM in EPT\n')
        if self.is_debug_file:
            self.binary_debug.write('skipping PCONVM\n')
        return len(data)

    def _read_pdamp(self, data, n):
        """
        PDAMP(202,2,45) - the marker for Record ???
        """
        ntotal = 8  # 2*4
        s = Struct(b(self._endian + 'if'))
        nentries = (len(data) - n) // ntotal
        for i in range(nentries):
            out = s.unpack(data[n:n+ntotal])
            #(pid, b) = out
            prop = PDAMP.add_op2_data(out)
            self._add_op2_property(prop)
            n += ntotal
        self.card_count['PDAMP'] = nentries
        return n

# PDAMPT
# PDAMP5
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
        s = Struct(b(self._endian + 'i3f'))
        ntotal = 16  # 4*4
        nproperties = (len(data) - n) // ntotal
        for i in range(nproperties):
            edata = data[n:n+16]
            out = s.unpack(edata)
            #(pid,k,ge,s) = out
            if self.is_debug_file:
                self.binary_debug.write('  PELAS=%s\n' % str(out))
            prop = PELAS.add_op2_data(out)
            self._add_op2_property(prop)
            n += ntotal
        self.card_count['PELAS'] = nproperties
        return n

    def _read_pfast_msc(self, data, n):
        ntotal = 92 # 23*4
        s = Struct(b(self._endian + 'ifii 4f'))
        nproperties = (len(data) - n) // ntotal
        for i in range(nproperties):
            edata = data[n:n+ntotal]
            out = s.unpack(edata)
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
        s = Struct(b(self._endian + 'ifii 4f'))
        nproperties = (len(data) - n) // ntotal
        for i in range(nproperties):
            edata = data[n:n+ntotal]
            out = s.unpack(edata)
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
# PELAST

    def _read_pgap(self, data, n):
        """
        PGAP(3201,32,55) - the marker for Record 42
        """
        ntotal = 44
        s = Struct(b(self._endian + 'i10f'))
        nproperties = (len(data) - n) // ntotal
        for i in range(nproperties):
            edata = data[n:n+ntotal]
            out = s.unpack(edata)
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
        s = Struct(b(self._endian + 'ifff'))
        nproperties = (len(data) - n) // 16
        for i in range(nproperties):
            edata = data[n:n+16]
            out = s.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  PHBDY=%s\n' % str(out))
            #(pid, af, d1, d2) = out
            prop = PHBDY.add_op2_data(out)
            self.add_PHBDY(prop)
            n += 16
        self.card_count['PHBDY'] = nproperties
        return n

    def _read_pintc(self, data, n):
        self.log.debug('skipping PINTC in EPT\n')
        return len(data)

    def _read_pints(self, data, n):
        self.log.debug('skipping PINTS in EPT\n')
        return len(data)

    def _read_plplane(self, data, n):
        self.log.debug('skipping PLPLANE in EPT\n')
        return len(data)

    def _read_plsolid(self, data, n):
        self.log.debug('skipping PLSOLID in EPT\n')
        return len(data)

    def _read_pmass(self, data, n):
        """
        PMASS(402,4,44) - the marker for Record 48
        """
        n = 0
        s = self.struct_2i
        nentries = (len(data) - n) // 8  # 2*4
        for i in range(nentries):
            edata = data[n:n + 8]
            out = s.unpack(edata)
            #out = (pid,mass)
            if self.is_debug_file:
                self.binary_debug.write('  PMASS=%s\n' % str(out))
            prop = PMASS(data=out)
            self._add_op2_property(prop)
            n += 8
        return n

    def _read_prod(self, data, n):
        """
        PROD(902,9,29) - the marker for Record 49
        """
        ntotal = 24  # 6*4
        s = Struct(b(self._endian + '2i4f'))
        nproperties = (len(data) - n) // ntotal
        for i in range(nproperties):
            edata = data[n:n+24]
            out = s.unpack(edata)
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
        s = Struct(b(self._endian + '2i4f'))
        nproperties = (len(data) - n) // 24
        for i in range(nproperties):
            edata = data[n:n+24]
            out = s.unpack(edata)
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
        s = Struct(b(self._endian + 'iififi4fi'))
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
                self.bigProperties[pid] = prop
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
        s = Struct(b(self._endian + '6i4s'))
        nproperties = (len(data) - n) // ntotal
        for i in range(nproperties):
            edata = data[n:n+28]
            out = s.unpack(edata)
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
        .. todo:: OD2 only exists for heat transfer...how do i know if there's heat transfer at this point...
        .. todo:: I could store all the tubes and add them later, but what about themal/non-thermal subcases
        .. warning:: assuming OD2 is not written (only done for thermal)
        """
        s = Struct(b(self._endian + '2i3f'))
        nproperties = (len(data) - n) // 20
        for i in range(nproperties):
            edata = data[n:n+20]  # or 24???
            out = s.unpack(edata)
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
        self.log.debug('skipping PSET in EPT\n')
        return len(data)

    def _read_pval(self, data, n):
        self.log.debug('skipping PVAL in EPT\n')
        return len(data)

    def _read_pvisc(self, data, n):
        """PVISC(1802,18,31) - the marker for Record 39"""
        s = Struct(b(self._endian + 'i2f'))
        nproperties = (len(data) - n) // 12
        for i in range(nproperties):
            edata = data[n:n+12]
            out = s.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  PVISC=%s\n' % str(out))
            #(pid,ce,cr) = out
            prop = PVISC.add_op2_data(out)
            self._add_op2_property(prop)
            n += 12
        self.card_count['PVISC'] = nproperties
        return n

# PWELD
# PWSEAM
    def _read_view(self, data, n):
        self.log.debug('skipping VIEW in EPT\n')
        return len(data)

    def _read_view3d(self, data, n):
        self.log.debug('skipping VIEW3D in EPT\n')
        return len(data)
