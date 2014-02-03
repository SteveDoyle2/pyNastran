#pylint: disable=C0301,W0612,C0111,R0201,C0103,W0613,R0914,C0326
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
from struct import unpack, Struct
import StringIO

from pyNastran.bdf.bdf import (NSM, PBAR, PBARL, PBEAM,
                               PROD, PSHELL, PSHEAR,
                               PCOMP, PSOLID,
                               PVISC, PELAS, PMASS,
                               PTUBE, PGAP)
# PCOMPG, PBUSH1D, PBEAML, PBEAM3, PBUSH,


class EPT(object):
    def add_property(self, card, allowOverwrites=True):
        raise RuntimeError('this should be overwritten')
    def readFake(self, data, n):
        raise RuntimeError('this should be overwritten')

    def __init__(self):
        self.skippedCardsFile = StringIO.StringIO()
        self.card_count = {}
        self.bigProperties = {}
        self._ept_map = {
            (3201,32,55):    ['NSM',    self.readNSM],     # record 2  - needs an object holder (e.g. self.elements/self.properties)
            (52,   20, 181): ['PBAR',   self.readPBAR],    # record 11 - buggy
            (9102, 91,  52): ['PBARL',  self.readPBARL],   # record 12 - almost there...
            (2706, 27, 287): ['PCOMP',  self.readPCOMP],   # record 22 - buggy
            (302,   3,  46): ['PELAS',  self.readPELAS],   # record 39
            (2102, 21, 121): ['PGAP',   self.readPGAP],    # record 42
            (902,   9,  29): ['PROD',   self.readPROD],    # record 49
            (1002, 10,  42): ['PSHEAR', self.readPSHEAR],  # record 50
            (2402, 24, 281): ['PSOLID', self.readPSOLID],  # record 51
            (2302, 23, 283): ['PSHELL', self.readPSHELL],  # record 52
            (1602, 16,  30): ['PTUBE',  self.readPTUBE],   # record 56

            (5402,  54, 262): ['PBEAM',   self.readPBEAM],   # record 14 - not done
            (9202,  92,  53): ['PBEAML',  self.readPBEAML],  # record 15 - not done
            (2502,  25, 248): ['PBEND',   self.readPBEND],   # record 16 - not done
            (3101,  31, 219): ['PBUSH1D', self.readPBUSH1D], # record 20 - not done
            (152,   19, 147): ['PCONEAX', self.readPCONEAX], # record 24 - not done
            (11001,110, 411): ['PCONV',   self.readPCONV],   # record 25 - not done
            # record 26
            (202,    2,  45): ['PDAMP',   self.readPDAMP],   # record 27 - not done
            (2802,  28, 236): ['PHBDY',   self.readPHBDY],   # record 43 - not done
            (402,    4,  44): ['PMASS',   self.readPMASS],   # record 48
            (1802,  18,  31): ['PVISC',   self.readPVISC],   # record 59
            (10201,102, 400): ['PVAL',   self.readPVAL],     # record 58 - not done
            (2606,  26, 289): ['VIEW',   self.readVIEW],     # record 62 - not done
            (1402,   14,  37): ['', self.readFake],    # record
            (2706,   27, 287): ['', self.readFake],    # record
            (702,     7,  38): ['', self.readFake],    # record
            (10301, 103, 399): ['', self.readFake],
            (5403, 55, 349): ['', self.readFake],
            (6902, 69, 165): ['', self.readFake],
            (3002, 30, 415): ['', self.readFake],
            (13301, 133, 509): ['', self.readFake],
            (6802, 68, 164): ['', self.readFake],
            (4606, 46, 375): ['', self.readFake],
            (1302, 13, 34): ['', self.readFake],
            (4706, 47, 376): ['', self.readFake],
            (8702, 87, 412): ['', self.readFake],
            (2902, 29, 420): ['', self.readFake],
            (1502, 15, 36): ['', self.readFake],
            (3301, 33, 56): ['', self.readFake],
            (3401, 34, 57): ['', self.readFake],    # record
            (1202, 12, 33): ['', self.readFake],  # record
            (12001, 120, 480): ['', self.readFake],  # record
            (12101, 121, 484): ['', self.readFake],  # record
            (3501, 35, 58): ['', self.readFake],  # record
            (3601, 36, 62): ['', self.readFake],  # record
            (8300, 83, 382): ['', self.readFake],  # record
            (8500, 85, 384): ['', self.readFake],  # record
        }

    def addOp2Property(self, prop):
        self.add_property(prop, allowOverwrites=True)
        #print(str(prop)[:-1])

# HGSUPPR

    def readNSM(self, data, n):
        """
        NSM(3201,32,55) - the marker for Record 2
        .. todo:: this isnt a property...
        """
        #print "reading NSM"
        return n
        s = Struct(b'i4sif')
        while len(data) >= 16:  # 4*4
            eData = data[:16]
            data = data[16:]
            out = s.unpack(eData)
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
        s = Struct(b'2i17f')
        nentries = (len(data) - n) // ntotal
        for i in xrange(nentries):
            eData = data[n:n+76]
            out = s.unpack(eData)
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
        s = Struct(b'2i8s8sf')
        nentries = (len(data) - n) // ntotal
        for i in xrange(nentries):
            eData = data[n:n+28]
            out = s.unpack(eData)
            (pid, mid, group, Type, value) = out
            Type = Type.strip()
            dataIn = [pid, mid, group, Type, value]
            #print "pid=%s mid=%s group=|%s| Type=|%s| value=%s" %(pid,mid,group,Type,value)
            expectedLength = validTypes[Type]
            iFormat = b'%if' % expectedLength
            dataIn += list(unpack(iFormat, data[:expectedLength * 4]))

            # TODO why do i need the +4???
            data = data[expectedLength * 4 + 4:]

            #print "len(out) = ",len(out)
            #print "PBARL = ",dataIn
            prop = PBARL(None, dataIn)
            self.addOp2Property(prop)
            #print self.print_section(20)
        self.card_count['PBARL'] = nentries
        return n

# PBCOMP

    def readPBEAM(self, data, n):
        """
        PBEAM(5402,54,262) - the marker for Record 14
        .. todo:: add object
        """
        #print "reading PBEAM"
        s1 = Struct(b'4if')
        s2 = Struct(b'16f')
        s3 = Struct(b'11f')
        ntotal = 1072  # 44+12*84+20
        nproperties = (len(data) - n) // ntotal
        for i in xrange(nproperties):
            eData = data[n:n+20]
            n += 20
            dataIn = list(s1.unpack(eData))
            #print "len(out) = ",len(out)
            #print out
            (pid, mid, nsegs, ccf, x) = dataIn

            for i in xrange(12):
                eData = data[n:n+64]
                n += 64
                pack = s2.unpack(eData)
                (so, xxb, a, i1, i2, i12, j, nsm, c1, c2,
                    d1, d2, e1, e2, f1, f2) = pack
                dataIn.append(pack)

            eData = data[n:n+44]

            dataIn = list(s3.unpack(eData))
            #(k1,k2,s1,s2,nsia,nsib,cwa,cwb,m1a,m2a,m1b,m2b,n1a,n2a,n1b,n2b) = pack

            prop = PBEAM(None, dataIn)
            self.addOp2Property(prop)
            #sys.exit('ept-PBEAM')
        self.card_count['PBEAM'] = nproperties
        return n

    def readPBEAML(self, data, n):
        self.skippedCardsFile.write('skipping PBEAML in EPT\n')
        return n

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
        nproperties = 0
        n2 = n
        s1 = Struct(b'2i3fi2f')
        s2 = Struct(b'i2fi')
        while n2 < n:  #len(data) >= 32:  # 8*4 - dynamic
            #print "len(data) = ",len(data)
            #print self.print_block(data[0:200])
            isSymmetrical = 'NO'
            eData = data[n:n+32]
            out = s1.unpack(eData)
            (pid, nLayers, z0, nsm, sb, ft, Tref, ge,) = out

            eData = data[n:n+16 * (nLayers)]
            Mid = []
            T = []
            Theta = []
            Sout = []
            if nLayers < 0:
                isSymmetrical = 'YES'
                nLayers = abs(nLayers)
            #print "nLayers = ",nLayers
            assert 0 < nLayers < 100, 'pid=%s nLayers=%s z0=%s nms=%s sb=%s ft=%s Tref=%s ge=%s' % (pid, nLayers, z0, nsm, sb, ft, Tref, ge)

            idata = 0
            for ilayer in xrange(nLayers):
                (mid, t, theta, sout) = s2.unpack(eData[idata:idata+16])
                Mid.append(mid)
                T.append(t)
                Theta.append(theta)
                Sout.append(sout)
                idata += 16

            dataIn = [pid, z0, nsm, sb, ft, Tref, ge,
                      isSymmetrical, Mid, T, Theta, Sout]
            #print "PCOMP = %s" % (dataIn)
            prop = PCOMP(None, dataIn)
            self.addOp2Property(prop)
            nproperties += 1
        self.card_count['PCOMP'] = nproperties
        return n

# PCOMPA
    def readPCONEAX(self, data, n):  # 24
        self.skippedCardsFile.write('skipping PCONEAX\n')
    def readPCONV(self, data, n):  # 25
        self.skippedCardsFile.write('skipping PCONV\n')
    def readPCONVM(self, data, n):  # 26
        self.skippedCardsFile.write('skipping PCONVM\n')

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
        s = Struct(b'i3f')
        ntotal = 16  # 4*4
        nproperties = (len(data) - n) // ntotal
        for i in xrange(nproperties):
            eData = data[n:n+16]
            out = s.unpack(eData)
            #(pid,k,ge,s) = out
            prop = PELAS(data=out)
            self.addOp2Property(prop)
            n += ntotal
        self.card_count['PELAS'] = nproperties
        return n

# PFAST
# PELAST

    def readPGAP(self, data, n):
        """
        PGAP(3201,32,55) - the marker for Record 42
        """
        #print "reading PGAP"
        s = Struct(b'i10f')
        nproperties = (len(data) - n) // 44
        for i in xrange(nproperties):
            eData = data[n:n+44]
            out = s.unpack(eData)
            #(pid,u0,f0,ka,kb,kt,mu1,mu2,tmax,mar,trmin) = out
            prop = PGAP(None, out)
            self.addOp2Property(prop)
        return n

    def readPHBDY(self, data, n):
        return n
    def readPINTC(self, data, n):
        return n
    def readPINTS(self, data, n):
        return n
    def readPLPLANE(self, data, n):
        return n
    def readPLSOLID(self, data, n):
        return n

    def readPMASS(self, data, n):
        """
        PMASS(402,4,44) - the marker for Record 48
        """
        n = 0
        s = Struct(b'ii')
        nEntries = (len(data) - n) // 8  # 2*4
        for i in xrange(nEntries):
            eData = data[n:n + 8]
            out = s.unpack(eData)
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
        s = Struct(b'2i4f')
        nproperties = (len(data) - n) // ntotal
        for i in xrange(nproperties):
            eData = data[n:n+24]
            out = s.unpack(eData)
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
        s = Struct(b'2i4f')
        nproperties = (len(data) - n) // 24
        for i in xrange(nproperties):
            eData = data[n:n+24]
            out = s.unpack(eData)
            (pid, mid, t, nsm, f1, f2) = out
            prop = PSHEAR(data=out)
            self.addOp2Property(prop)
            n += 24
        self.card_count['PSHEAR'] = nproperties
        return n

    def readPSHELL(self, data, n):
        """
        PSHELL(2302,23,283) - the marker for Record 51
        """
        #print "reading PSHELL"
        ntotal = 44  # 11*4
        s = Struct(b'iififi4fi')
        nproperties = (len(data) - n) // ntotal
        for i in xrange(nproperties):
            eData = data[n:n+44]
            out = s.unpack(eData)
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
        s = Struct(b'6i4s')
        nproperties = (len(data) - n) // ntotal
        for i in xrange(nproperties):
            eData = data[n:n+28]
            (pid, mid, cid, inp, stress, isop, fctn) = s.unpack(eData)
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
        s = Struct(b'2i3f')
        nproperties = (len(data) - n) // 20
        for i in xrange(nproperties):
            eData = data[n:n+20]  # or 24???
            (pid, mid, OD, t, nsm) = s.unpack(eData)
            dataIn = [pid, mid, OD, t, nsm]
            prop = PTUBE(None, dataIn)
            self.addOp2Property(prop)
            n += 20
        self.card_count['PTUBE'] = nproperties

    def readPSET(self, data, n):
        return n
    def readPVAL(self, data, n):
        return n

    def readPVISC(self, data, n):
        """PVISC(1802,18,31) - the marker for Record 39"""
        #print "reading PVISC"
        s = Struct(b'i2f')
        nproperties = (len(data) - n) // 12
        for i in xrange(nproperties):
            eData = data[n:n+12]
            out = s.unpack(eData)
            #(pid,ce,cr) = out
            prop = PVISC(data=out)
            self.addOp2Property(prop)
            n += 12
        self.card_count['PVISC'] = nproperties
        return n

# PWELD
# PWSEAM
    def readVIEW(self, data, n):
        return n
    def readVIEW3D(self, data, n):
        return n
