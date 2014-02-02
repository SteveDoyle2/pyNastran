from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
from struct import unpack, Struct

from pyNastran.bdf.bdf import (NSM, PBAR, PBARL, PBEAM,
                               PROD, PSHELL, PSHEAR,
                               PCOMP, PSOLID,
                               PVISC, PELAS, PMASS,
                               PTUBE, PGAP)
# PCOMPG,PBUSH1D,PBEAML,PBEAM3,PBUSH,


class EPT(object):

    def __init__(self):
        self.bigProperties = {}
        self._ept_map = {
            #(3201,32,55):    self.readNSM,     # record 2  - needs an object holder (e.g. self.elements/self.properties)
            (52, 20, 181): self.readPBAR,    # record 11 - buggy
            (9102, 91, 52): self.readPBARL,   # record 12 - almost there...
            #(2706,27,287):   self.readPCOMP,   # record 22 - buggy
            (302, 3, 46): self.readPELAS,   # record 39
            (2102, 21, 121): self.readPGAP,    # record 42
            (902, 9, 29): self.readPROD,    # record 49
            (1002, 10, 42): self.readPSHEAR,  # record 50
            (2402, 24, 281): self.readPSOLID,  # record 51
            (2302, 23, 283): self.readPSHELL,  # record 52
            (1602, 16, 30): self.readPTUBE,   # record 56

            #(5402,54,262):   self.readPBEAM,   # record 14 - not done
            #(9202, 92,  53): self.readPBEAML,  # record 15 - not done
            (2502, 25, 248): self.readPBEND,   # record 16 - not done
            (3101, 31, 219): self.readPBUSH1D,  # record 20 - not done
            #(152,  19, 147): self.readPCONEAX, # record 24 - not done
            #(11001,110,411): self.readPCONV,   # record 25 - not done
            (202, 2, 45): self.readPDAMP,   # record 27 - not done
            #(2802, 28, 236): self.readPHBDY,   # record 43 - not done
            (402, 4, 44): self.readPMASS,   # record 48
            (1802, 18, 31): self.readPVISC,   # record 59
            #(10201,102,400): self.readPVAL,    # record 58 - not done
            #(2606, 26, 289): self.readVIEW,    # record 62 - not done
            #(1402, 14, 37):   self.readFake,    # record
        }

    def addOp2Property(self, prop):
        self.add_property(prop, allowOverwrites=True)
        print(str(prop)[:-1])

# HGSUPPR

    def readNSM(self, data, n):
        """
        NSM(3201,32,55) - the marker for Record 2
        .. todo:: this isnt a property...
        """
        #print "reading NSM"
        return
        while len(data) >= 16:  # 4*4
            eData = data[:16]
            data = data[16:]
            out = unpack(b'i4sif', eData)
            (sid, propSet, ID, value) = out
            #print "sid=%s propSet=%s ID=%s value=%s" %(sid,propSet,ID,value)
            prop = NSM(None, None, [sid, propSet, ID, value])
            #self.addOp2Property(prop)
        return n

# NSM1
# NSML1
# NSMADD
# NSML
# NSML1
# PAABSF
# PACABS
# PACBAR
    def readPBAR(self, data, n):
        """
        PBAR(52,20,181) - the marker for Record 11
        .. warning:: this makes a funny property...
        """
        #print "reading PBAR"
        ntotal = 76  # 19*4
        nentries = (len(data) - n) // ntotal
        for i in xrange(nentries):
            eData = data[n:n+76]
            out = unpack(b'2i17f', eData)
            #print out
            (pid, mid, a, I1, I2, J, nsm, fe, c1, c2, d1, d2,
                e1, e2, f1, f2, k1, k2, I12) = out
            prop = PBAR(None, out)
            self.addOp2Property(prop)
            n += ntotal
        self.card_count['PBAR'] = nentries
        return n

    def readPBARL(self, data, n):
        """
        PBARL(9102,91,52) - the marker for Record 12
        """
        validTypes = {
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
        }  # for GROUP="MSCBML0"

        #print "reading PBARL"
        ntotal = 28# 7*4 - ROD - shortest entry...could be buggy... # TODO fix this
        nEntries = (len(data) - n) // ntotal
        for i in xrange(nEntries):
            eData = data[n:n+28]
            out = unpack(b'2i8s8sf', eData)
            (pid, mid, group, Type, value) = out
            Type = Type.strip()
            dataIn = [pid, mid, group, Type, value]
            #print "pid=%s mid=%s group=|%s| Type=|%s| value=%s" %(pid,mid,group,Type,value)
            expectedLength = validTypes[Type]
            iFormat = 'f' * expectedLength
            iFormat = bytes(iFormat)
            dataIn += list(unpack(iFormat, data[:expectedLength * 4]))

            data = data[expectedLength * 4 + 4:]
                # TODO why do i need the +4???

            #print "len(out) = ",len(out)
            #print "PBARL = ",dataIn
            prop = PBARL(None, dataIn)
            self.addOp2Property(prop)
            #print self.print_section(20)
        self.card_count['PBARL'] = nEntries
        return n

# PBCOMP

    def readPBEAM(self, data, n):
        """
        PBEAM(5402,54,262) - the marker for Record 14
        .. todo:: add object
        """
        #print "reading PBEAM"
        while len(data) >= 1072:  # 44+12*84+20
            eData = data[:20]
            data = data[20:]
            dataIn = list(unpack(b'4if', eData))
            #print "len(out) = ",len(out)
            #print out
            (pid, mid, nsegs, ccf, x) = dataIn

            for i in xrange(12):
                eData = data[64:]
                data = data[:64]
                pack = unpack(b'16f', eData)
                (so, xxb, a, i1, i2, i12, j, nsm, c1, c2,
                    d1, d2, e1, e2, f1, f2) = pack
                dataIn.append(pack)

            eData = data[:44]
            data = data[44:]

            dataIn = list(unpack(b'11f', eData))
            #(k1,k2,s1,s2,nsia,nsib,cwa,cwb,m1a,m2a,m1b,m2b,n1a,n2a,n1b,n2b) = pack

            prop = PBEAM(None, dataIn)
            self.addOp2Property(prop)
            #sys.exit('ept-PBEAM')

# PBEAML

    def readPBEND(self, data, n):
        self.skippedCardsFile.write('skipping PBEND in EPT\n')
        return n

# PBMSECT
# PBRSECT

    def readPBUSH(self, data, n):
        self.skippedCardsFile.write('skipping PBUSH in EPT\n')
        return n

    def readPBUSH1D(self, data, n):
        self.skippedCardsFile.write('skipping PBUSH1D in EPT\n')
        return n

    def readPBUSHT(self, data, n):
        self.skippedCardsFile.write('skipping PBUSHT in EPT\n')
        return n

    def readPCOMP(self, data, n):
        """
        PCOMP(2706,27,287) - the marker for Record 22
        """
        #print "reading PCOMP"
        while len(data) >= 32:  # 8*4 - dynamic
            #print "len(data) = ",len(data)
            #print self.print_block(data[0:200])
            isSymmetrical = 'NO'
            eData = data[:32]
            data = data[32:]
            out = unpack(b'2i3fi2f', eData)
            (pid, nLayers, z0, nsm, sb, ft, Tref, ge,) = out

            eData = data[:16 * (nLayers)]
            data = data[16 * (nLayers):]

            Mid = []
            T = []
            Theta = []
            Sout = []
            if nLayers < 0:
                isSymmetrical = 'YES'
                nLayers = abs(nLayers)
            #print "nLayers = ",nLayers
            assert 0 < nLayers < 100, 'pid=%s nLayers=%s z0=%s nms=%s sb=%s ft=%s Tref=%s ge=%s' % (pid, nLayers, z0, nsm, sb, ft, Tref, ge)

            for n in xrange(nLayers):
                #print "len(eData) = ",len(eData)
                (mid, t, theta, sout) = unpack(b'i2fi', eData[0:16])
                Mid.append(mid)
                T.append(t)
                Theta.append(theta)
                Sout.append(sout)
                eData = eData[16:]

            dataIn = [pid, z0, nsm, sb, ft, Tref, ge,
                      isSymmetrical, Mid, T, Theta, Sout]
            #print "PCOMP = %s" %(dataIn)
            prop = PCOMP(None, dataIn)
            self.addOp2Property(prop)

# PCOMPA
# PCONEAX
# PCONV
# PCONVM

    def readPDAMP(self, data, n):
        self.skippedCardsFile.write('skipping PDAMP\n')
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

    def readPELAS(self, data, n):
        """PELAS(302,3,46) - the marker for Record 39"""
        #print "reading PELAS"
        while len(data) >= 16:  # 4*4
            eData = data[:16]
            data = data[16:]
            out = unpack(b'i3f', eData)
            #(pid,k,ge,s) = out
            prop = PELAS(data=out)
            self.addOp2Property(prop)
        return n

# PFAST
# PELAST

    def readPGAP(self, data, n):
        """
        PGAP(3201,32,55) - the marker for Record 42
        """
        #print "reading PGAP"
        while len(data) >= 44:  # 11*4
            eData = data[:44]
            data = data[44:]
            out = unpack(b'i10f', eData)
            #(pid,u0,f0,ka,kb,kt,mu1,mu2,tmax,mar,trmin) = out
            prop = PGAP(None, out)
            self.addOp2Property(prop)
        return n

# PHBDY
# PINTC
# PINTS
# PLPLANE
# PLSOLID

    def readPMASS(self, data, n):
        """
        PMASS(402,4,44) - the marker for Record 48
        """
        n = 0
        nEntries = len(data) // 8  # 2*4
        for i in xrange(nEntries):
            eData = data[n:n + 8]
            out = unpack(b'ii', eData)
            #out = (pid,mass)
            prop = PMASS(data=out)
            self.addOp2Property(prop)
            n += 8
        return n

    def readPROD(self, data, n):
        """
        PROD(902,9,29) - the marker for Record 49
        """
        #print "reading PROD"
        ntotal = 24  # 6*4
        nproperties = (len(data) - n) // ntotal
        for i in xrange(nproperties):
            eData = data[n:n+24]
            out = unpack(b'2i4f', eData)
            (pid, mid, a, j, c, nsm) = out
            prop = PROD(None, out)
            self.addOp2Property(prop)
            n += ntotal
        self.card_count['PROD'] = nproperties

    def readPSHEAR(self, data, n):
        """
        PSHEAR(1002,10,42) - the marker for Record 50
        """
        #print "reading PSHEAR"
        while len(data) >= 24:  # 6*4
            eData = data[:24]
            data = data[24:]
            out = unpack(b'2i4f', eData)
            (pid, mid, t, nsm, f1, f2) = out
            prop = PSHEAR(data=out)
            self.addOp2Property(prop)

    def readPSHELL(self, data, n):
        """
        PSHELL(2302,23,283) - the marker for Record 51
        """
        #print "reading PSHELL"
        ntotal = 44  # 11*4
        nproperties = (len(data) - n) // ntotal
        for i in xrange(nproperties):
            eData = data[n:n+44]
            out = unpack(b'iifififfffi', eData)
            (pid, mid1, t, mid2, bk, mid3, ts, nsm, z1, z2, mid4) = out
            prop = PSHELL(None, out)

            if max(pid, mid1, mid2, mid3, mid4) > 1e8:
                #print "PSHELL = ",out
                self.bigProperties[pid] = prop
            else:
                self.addOp2Property(prop)
            n += ntotal
        self.card_count['PSHELL'] = nproperties


    def readPSOLID(self, data, n):
        """
        PSOLID(2402,24,281) - the marker for Record 52
        """
        #print "reading PSOLID"
        ntotal = 28  # 7*4
        nproperties = (len(data) - n) // ntotal
        for i in xrange(nproperties):
            eData = data[n:n+28]
            (pid, mid, cid, inp, stress, isop, fctn) = unpack(b'6i4s', eData)
            dataIn = [pid, mid, cid, inp, stress, isop, fctn]
            prop = PSOLID(None, dataIn)
            self.addOp2Property(prop)
            n += ntotal
        self.card_count['PSOLID'] = nproperties

# PSOLIDL
# PTRIA6
# PTSHELL

    def readPTUBE(self, data, n):
        """
        PTUBE(1602,16,30) - the marker for Record 56
        .. todo:: OD2 only exists for heat transfer...how do i know if there's heat transfer at this point...
        .. todo:: I could store all the tubes and add them later, but what about themal/non-thermal subcases
        .. warning:: assuming OD2 is not written (only done for thermal)
        """
        #print "reading PTUBE"
        while len(data) >= 20:  # 5*4
            eData = data[:20]
            data = data[20:]  # or 24
            (pid, mid, OD, t, nsm) = unpack(b'2i3f', eData)
            dataIn = [pid, mid, OD, t, nsm]
            prop = PTUBE(None, dataIn)
            self.addOp2Property(prop)

# PSET
# PVAL

    def readPVISC(self, data, n):
        """PVISC(1802,18,31) - the marker for Record 39"""
        #print "reading PVISC"
        while len(data) >= 12:  # 3*4
            eData = data[:12]
            data = data[12:]
            out = unpack(b'i2f', eData)
            #(pid,ce,cr) = out
            prop = PVISC(data=out)
            self.addOp2Property(prop)

# PWELD
# PWSEAM
# PVIEW
# PVIEW3D
