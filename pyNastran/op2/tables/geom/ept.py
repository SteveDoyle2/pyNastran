#pylint: disable=C0301,W0612,C0111,R0201,C0103,W0613,R0914
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
from struct import unpack, Struct
from six import b
from six.moves import range

from pyNastran.bdf.bdf import (NSM, PBAR, PBARL, #PBEAM,
                               PROD, PSHELL, PSHEAR,
                               PCOMP, PSOLID,
                               PVISC, PELAS, PMASS,
                               PTUBE, PGAP)
# PCOMPG, PBUSH1D, PBEAML, PBEAM3, PBUSH,
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
            (702, 7, 38): ['', self._read_fake],            # record
            (10301, 103, 399): ['', self._read_fake],
            (5403, 55, 349): ['', self._read_fake],
            (6902, 69, 165): ['', self._read_fake],
            (3002, 30, 415): ['', self._read_fake],
            (13301, 133, 509): ['', self._read_fake],
            (6802, 68, 164): ['', self._read_fake],
            (4606, 46, 375): ['', self._read_fake],
            (1302, 13, 34): ['', self._read_fake],
            (4706, 47, 376): ['', self._read_fake],
            (8702, 87, 412): ['', self._read_fake],
            (2902, 29, 420): ['', self._read_fake],
            (1502, 15, 36): ['', self._read_fake],
            (3201, 32, 991) : ['', self._read_fake],  # record
            (3301, 33, 992) : ['', self._read_fake],  # record
            (3301, 33, 56): ['', self._read_fake],
            (3401, 34, 57) : ['', self._read_fake],    # record
            (3701, 37, 995) : ['', self._read_fake],    # record
            (1202, 12, 33): ['', self._read_fake],  # record
            (12001, 120, 480): ['', self._read_fake],  # record
            (12101, 121, 484): ['', self._read_fake],  # record
            (3501, 35, 58): ['', self._read_fake],  # record
            (3601, 36, 62): ['', self._read_fake],  # record
            (8300, 83, 382): ['', self._read_fake],  # record
            (8500, 85, 384): ['', self._read_fake],  # record
            (15006, 150, 604): ['', self._read_fake],  # record
        }

    def _add_op2_property(self, prop):
        self.add_property(prop, allow_overwrites=True)
        #print(str(prop)[:-1])

# HGSUPPR

    def _read_nsm(self, data, n):
        """
        NSM(3201,32,55) - the marker for Record 2
        .. todo:: this isnt a property...
        """
        return n
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
            Type = Type.strip()
            data_in = [pid, mid, group, Type, value]
            print("pid=%s mid=%s group=%r Type=%r value=%s" %(pid, mid, group, Type, value))
            if pid > 100000000:
                raise RuntimeError('bad parsing...')
            expectedLength = valid_types[Type]
            iFormat = b'%if' % expectedLength

            ndelta = expectedLength * 4
            data_in += list(unpack(iFormat, data[n:n+ndelta]))
            # TODO why do i need the +4???
            #min_len =  expectedLength * 4 + 4
            #if len(data)
            #data = data[n + expectedLength * 4 + 4:]
            n += ndelta

            #prin( "len(out) = ",len(out)))
            #print("PBARL = %s" % data_in)
            prop = PBARL.add_op2_data(data_in)  # last value is nsm
            self._add_op2_property(prop)
            #print(self.show_data(data[n-8:-100]))
            break
        self._increase_card_count('PBARL')
        #assert len(data) == n
        return len(data)

# PBCOMP

    def _read_pbeam(self, data, n):
        """
        PBEAM(5402,54,262) - the marker for Record 14
        .. todo:: add object
        """
        s1 = Struct(b(self._endian + '4if'))
        s2 = Struct(b(self._endian + '16f'))
        s3 = Struct(b(self._endian + '11f'))
        ntotal = 1072  # 44+12*84+20
        nproperties = (len(data) - n) // ntotal
        for i in range(nproperties):
            edata = data[n:n+20]
            n += 20
            data_in = list(s1.unpack(edata))
            if self.is_debug_file:
                self.binary_debug.write('  PBEAM=%s\n' % str(data_in))
            (pid, mid, nsegs, ccf, x) = data_in

            for i in range(12):
                edata = data[n:n+64]
                n += 64
                pack = s2.unpack(edata)
                (so, xxb, a, i1, i2, i12, j, nsm, c1, c2,
                 d1, d2, e1, e2, f1, f2) = pack
                data_in.append(pack)
                if self.is_debug_file:
                    self.binary_debug.write('     %s\n' % str(pack))
            edata = data[n:n+44]

            data_in = list(s3.unpack(edata))
            #(k1,k2,s1,s2,nsia,nsib,cwa,cwb,m1a,m2a,m1b,m2b,n1a,n2a,n1b,n2b) = pack

            # prop = PBEAM(None, data_in)
            # self._add_op2_property(prop)
            #sys.exit('ept-PBEAM')
        self.card_count['PBEAM'] = nproperties
        return n

    def _read_pbeaml(self, data, n):
        if self.is_debug_file:
            self.binary_debug.write('skipping PBEAML in EPT\n')
        return len(data)

    def _read_pbend(self, data, n):
        if self.is_debug_file:
            self.binary_debug.write('skipping PBEND in EPT\n')
        return len(data)

# PBMSECT
# PBRSECT

    def _read_pbush(self, data, n):
        if self.is_debug_file:
            self.binary_debug.write('skipping PBUSH in EPT\n')
        return len(data)

    def _read_pbush1d(self, data, n):
        if self.is_debug_file:
            self.binary_debug.write('skipping PBUSH1D in EPT\n')
        return len(data)

    def _read_pbusht(self, data, n):
        if self.is_debug_file:
            self.binary_debug.write('skipping PBUSHT in EPT\n')
        return len(data)

    def _read_pcomp(self, data, n):
        """
        PCOMP(2706,27,287) - the marker for Record 22
        """
        nproperties = 0
        n2 = n
        s1 = Struct(b(self._endian + '2i3fi2f'))
        s2 = Struct(b(self._endian + 'i2fi'))
        while n2 < n:  #len(data) >= 32:  # 8*4 - dynamic
            #print("len(data) = %s" % len(data))
            #print(self.print_block(data[0:200]))
            isSymmetrical = 'NO'
            edata = data[n:n+32]
            out = s1.unpack(edata)
            self.binary_debug.write('  PCOMP=%s\n' % str(out))
            (pid, nlayers, z0, nsm, sb, ft, Tref, ge,) = out

            edata = data[n:n+16 * (nlayers)]
            Mid = []
            T = []
            Theta = []
            Sout = []
            if nlayers < 0:
                is_symmetrical = 'YES'
                nlayers = abs(nlayers)
            #print("nlayers = ",nlayers)
            assert 0 < nlayers < 100, 'pid=%s nlayers=%s z0=%s nms=%s sb=%s ft=%s Tref=%s ge=%s' % (
                pid, nlayers, z0, nsm, sb, ft, Tref, ge)

            idata = 0
            for ilayer in range(nlayers):
                (mid, t, theta, sout) = s2.unpack(edata[idata:idata+16])
                Mid.append(mid)
                T.append(t)
                Theta.append(theta)
                Sout.append(sout)
                idata += 16

            data_in = [
                pid, z0, nsm, sb, ft, Tref, ge,
                is_symmetrical, Mid, T, Theta, Sout]
            #print("PCOMP = %s" % (data_in))
            prop = PCOMP.add_op2_data(data_in)
            self._add_op2_property(prop)
            nproperties += 1
        self.card_count['PCOMP'] = nproperties
        return n

# PCOMPA
    def _read_pconeax(self, data, n):  # 24
        if self.is_debug_file:
            self.binary_debug.write('skipping PCONEAX\n')
        return len(data)
    def _read_pconv(self, data, n):  # 25
        if self.is_debug_file:
            self.binary_debug.write('skipping PCONV\n')
        return len(data)
    def _read_pconvm(self, data, n):  # 26
        if self.is_debug_file:
            self.binary_debug.write('skipping PCONVM\n')
        return len(data)
    def _read_pdamp(self, data, n):
        if self.is_debug_file:
            self.binary_debug.write('skipping PDAMP\n')
        return len(data)

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

# PFAST
# PELAST

    def _read_pgap(self, data, n):
        """
        PGAP(3201,32,55) - the marker for Record 42
        """
        s = Struct(b(self._endian + 'i10f'))
        nproperties = (len(data) - n) // 44
        for i in range(nproperties):
            edata = data[n:n+44]
            out = s.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  PGAP=%s\n' % str(out))
            #(pid,u0,f0,ka,kb,kt,mu1,mu2,tmax,mar,trmin) = out
            prop = PGAP.add_op2_data(out)
            self._add_op2_property(prop)
        return n

    def _read_phbdy(self, data, n):
        return len(data)
    def _read_pintc(self, data, n):
        return len(data)
    def _read_pints(self, data, n):
        return len(data)
    def _read_plplane(self, data, n):
        return len(data)
    def _read_plsolid(self, data, n):
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
            (pid, mid, a, j, c, nsm) = out
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
            (pid, mid, t, nsm, f1, f2) = out
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
        return len(data)
    def _read_pval(self, data, n):
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
        return len(data)
    def _read_view3d(self, data, n):
        return len(data)
