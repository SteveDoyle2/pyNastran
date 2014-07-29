from numpy import array, zeros, arange, concatenate, searchsorted, where, unique

from pyNastran.bdf.fieldWriter import print_card_8
from pyNastran.bdf.bdfInterface.assign_type import (integer, #integer_or_blank,
    double, double_string_or_blank, integer_double_string_or_blank,
    string,
    double_or_blank)#, integer_double_or_blank, blank)


class PBEAM(object):
    type = 'PBEAM'
    def __init__(self, model):
        """
        Defines the PCOMP object.

        :param self: the PCOMP object
        :param model: the BDF object
        :param cards: the list of PCOMP cards
        """
        self.model = model
        self.n = 0
        self._cards = []
        self._comments = []

    def add(self, card, comment):
        self._cards.append(card)
        self._comments.append(comment)

    def build(self):
        cards = self._cards
        ncards = len(cards)
        self.n = ncards

        if ncards:
            #: Property ID
            self.property_id = zeros(ncards, 'int32')
            self.material_id = zeros(ncards, 'int32')
            self.area = zeros(ncards, 'float64')
            self.I1 = zeros(ncards, 'float64')
            self.I2 = zeros(ncards, 'float64')
            self.J = zeros(ncards, 'float64')
            self.nsm = zeros(ncards, 'float64')

            #: Shear stiffness factor K in K*A*G for plane 1.
            self.k1 = zeros(ncards, 'float32')
            #: Shear stiffness factor K in K*A*G for plane 2.
            self.k2 = zeros(ncards, 'float32')

            #: Shear relief coefficient due to taper for plane 1.
            self.s1 = zeros(ncards, 'float32')
            #: Shear relief coefficient due to taper for plane 2.
            self.s2 = zeros(ncards, 'float32')

            #: non structural mass moment of inertia per unit length
            #: about nsm center of gravity at Point A.
            self.nsia = zeros(ncards, 'float32')
            #: non structural mass moment of inertia per unit length
            #: about nsm center of gravity at Point B.
            self.nsib = zeros(ncards, 'float32')

            #: warping coefficient for end A.
            self.cwa = zeros(ncards, 'float32')
            #: warping coefficient for end B.
            self.cwb = zeros(ncards, 'float32')

            #: y coordinate of center of gravity of
            #: nonstructural mass for end A.
            self.m1a = zeros(ncards, 'float32')
            #: z coordinate of center of gravity of
            #: nonstructural mass for end A.
            self.m2a = zeros(ncards, 'float32')

            #: y coordinate of center of gravity of
            #: nonstructural mass for end B.
            self.m1b = zeros(ncards, 'float32')
            #: z coordinate of center of gravity of
            #: nonstructural mass for end B.
            self.m2b = zeros(ncards, 'float32')

            #: y coordinate of neutral axis for end A.
            self.n1a = zeros(ncards, 'float32')
            #: z coordinate of neutral axis for end A.
            self.n2a = zeros(ncards, 'float32')

            #: y coordinate of neutral axis for end B.
            self.n1b = zeros(ncards, 'float32')
            #: z coordinate of neutral axis for end B.
            self.n2b = zeros(ncards, 'float32')

            for i, card in enumerate(cards):
                self.add_card(card, i)

            i = self.property_id.argsort()
            self.property_id = self.property_id[i]
            self.material_id = self.material_id[i]

            self.area = self.area[i]
            self.I1 = self.I1[i]
            self.I2 = self.I2[i]
            self.J = self.J[i]
            self.nsm = self.nsm[i]

            unique_pids = unique(self.property_id)

            if len(unique_pids) != len(self.property_id):
                raise RuntimeError('There are duplicate PBEAM IDs...')
            self._cards = []
            self._comments = []
        #import sys
        #sys.exit()

    #=========================================================================
    def add_card(self, card, i):
        #: property ID
        self.property_id[i] = integer(card, 1, 'property_id')

        #: material ID
        self.material_id[i] = integer(card, 2, 'material_id')
        SO = []
        XXB = []
        A = []
        I1 = []
        I2 = []
        I12 = []
        J = []
        NSM = []

        # at least one cross section are required
        # so[0] and xxb[0] aren't used
        #: Output flag
        SO = ['YES']
        #: Section position
        XXB = [0.]
        A.append(double(card, 3, 'A'))
        I1.append(double_or_blank(card, 4, 'I1', 0.0))
        I2.append(double_or_blank(card, 5, 'I2', 0.0))
        I12.append(double_or_blank(card, 6, 'I12', 0.0))
        J.append(double_or_blank(card, 7, 'J', 0.0))
        NSM.append(double_or_blank(card, 8, 'nsm', 0.0))
        print(card)
        assert I1[0] > 0, I1[0]
        assert I2[0] > 0, I2[0]

        isCDEF = False
        field9 = double_string_or_blank(card, 9, 'field9')
        field17 = double_string_or_blank(card, 17, 'field17')
        try:
            isinstance(field17, float)
            isCDEF = True
            isFooter = True
        except SyntaxError:
            pass
        #print("f9=%s f17=%s" % (field9, field17))

        #nlines = nfields // 8

        if field9 in ['YES', 'YESA', 'NO']:
            isCDEF = False
            isContinue = True
        elif field17 in ['YES', 'YESA', 'NO']:
            isCDEF = True
            isContinue = True
        else:
            isContinue = False
            isCDEF = True

        C1 = [double_or_blank(card, 9, 'c1', 0.0)]
        C2 = [double_or_blank(card, 10, 'c2', 0.0)]
        D1 = [double_or_blank(card, 11, 'd1', 0.0)]
        D2 = [double_or_blank(card, 12, 'd2', 0.0)]
        E1 = [double_or_blank(card, 13, 'e1', 0.0)]
        E2 = [double_or_blank(card, 14, 'e2', 0.0)]
        F1 = [double_or_blank(card, 15, 'f1', 0.0)]
        F2 = [double_or_blank(card, 16, 'f2', 0.0)]

        # ----------------------------------------------------------------
        # figure out how many YES/YESA/NO fields there are
        # and if there is a footer

        # counting continuation cards
        nfields = len(card) - 1  # -1 for PBEAM field
        if isCDEF:
            nmajor = nfields // 16
            nleftover = nfields % 16
            if nleftover == 0:
                nmajor -= 1

            if nmajor == 0:
                nmajor = 1

            # jump to the last block of 16
            x = nmajor * 16 + 1

            # If it's an SO field, we don't read the footer
            # remark 6:
            # The fourth and fifth continuation entries, which
            # contain fields K1 through N2(B), are optional
            # and may be omitted if the default values are appropriate.
            val = integer_double_string_or_blank(card, x, 'YES/YESA/NO')
            if val in ['YES', 'YESA', 'NO']:  # there is no footer
                nmajor += 1
                x += 16
            else:
                # read the footer
                pass
        else:
            nmajor = nfields // 8
            nleftover = nfields % 8
            if nleftover == 0:
                nmajor -= 1
            if nmajor == 0:
                nmajor = 1
            x = nmajor * 8 + 1

            val = integer_double_string_or_blank(card, x, 'YES/YESA/NO')
            if val in ['YES', 'YESA', 'NO']:  # there is no footer
                nmajor += 1
                x += 8
            else:
                # read the footer
                pass

        # ----------------------------------------------------------------
        for nRepeated in xrange(1, nmajor): # start at 1 to drop the header
            # field 17 is the first possible so
            if isCDEF:
                nStart = nRepeated * 16 + 1
            else:
                nStart = nRepeated * 8 + 1

            so  = string(card,  nStart, 'SO%i' % nRepeated)
            xxb = double(card,  nStart + 1, 'x/xb%i' % nRepeated)
            a   = double_or_blank(card,  nStart + 2, 'Area%i' % nRepeated, 0.0)
            i1  = double_or_blank(card, nStart + 3, 'I1 %i' % nRepeated, 0.0)
            i2  = double_or_blank(card, nStart + 4, 'I2 %i' % nRepeated, 0.0)
            i12 = double_or_blank(card, nStart + 5, 'I12 %i' % nRepeated, 0.0)
            j   = double_or_blank(card, nStart + 6, 'J%i' % nRepeated, 0.0)
            nsm = double_or_blank(card, nStart + 7, 'nsm%i' % nRepeated, 0.0)

            SO.append(so)
            XXB.append(xxb)
            A.append(a)
            I1.append(i1)
            I2.append(i2)
            I12.append(i12)
            J.append(j)
            NSM.append(nsm)

            if isCDEF:
                c1 = double_or_blank(card, nStart + 8, 'c1 %i' % nRepeated, 0.0)
                c2 = double_or_blank(card, nStart + 9, 'c2 %i' % nRepeated, 0.0)
                d1 = double_or_blank(card, nStart + 10, 'd1 %i' % nRepeated, 0.0)
                d2 = double_or_blank(card, nStart + 11, 'd2 %i' % nRepeated, 0.0)
                e1 = double_or_blank(card, nStart + 12, 'e1 %i' % nRepeated, 0.0)
                e2 = double_or_blank(card, nStart + 13, 'e2 %i' % nRepeated, 0.0)
                f1 = double_or_blank(card, nStart + 14, 'f1 %i' % nRepeated, 0.0)
                f2 = double_or_blank(card, nStart + 15, 'f2 %i' % nRepeated, 0.0)
                C1.append(c1)
                C2.append(c2)
                D1.append(d1)
                D2.append(d2)
                E1.append(e1)
                E2.append(e2)
                F1.append(f1)
                F2.append(f2)
            else:
                # YESA or NO, values MUST be omitted; remark 5
                C1.append(None)
                C2.append(None)
                D1.append(None)
                D2.append(None)
                E1.append(None)
                E2.append(None)
                F1.append(None)
                F2.append(None)
        if len(XXB) > 1:
            assert min(XXB) == 0.0, 'min=%s, but should be 0.0\nxxb=%s' % (min(XXB), XXB)
            assert max(XXB) == 1.0, 'max=%s, but should be 1.0\nxxb=%s' % (max(XXB), XXB)


        # footer fields
        #: Shear stiffness factor K in K*A*G for plane 1.
        self.k1[i] = double_or_blank(card, x, 'k1', 1.0)
        #: Shear stiffness factor K in K*A*G for plane 2.
        self.k2[i] = double_or_blank(card, x + 1, 'k2', 1.0)

        #: Shear relief coefficient due to taper for plane 1.
        self.s1[i] = double_or_blank(card, x + 2, 's1', 0.0)
        #: Shear relief coefficient due to taper for plane 2.
        self.s2[i] = double_or_blank(card, x + 3, 's2', 0.0)

        #: non structural mass moment of inertia per unit length
        #: about nsm center of gravity at Point A.
        self.nsia[i] = double_or_blank(card, x + 4, 'nsia', 0.0)
        #: non structural mass moment of inertia per unit length
        #: about nsm center of gravity at Point B.
        self.nsib[i] = double_or_blank(card, x + 5, 'nsib', self.nsia)

        #: warping coefficient for end A.
        self.cwa[i] = double_or_blank(card, x + 6, 'cwa', 0.0)
        #: warping coefficient for end B.
        self.cwb[i] = double_or_blank(card, x + 7, 'cwb', self.cwa)

        #: y coordinate of center of gravity of
        #: nonstructural mass for end A.
        self.m1a[i] = double_or_blank(card, x + 8, 'm1a', 0.0)
        #: z coordinate of center of gravity of
        #: nonstructural mass for end A.
        self.m2a[i] = double_or_blank(card, x + 9, 'm2a', self.m1a)

        #: y coordinate of center of gravity of
        #: nonstructural mass for end B.
        self.m1b[i] = double_or_blank(card, x + 10, 'm1b', 0.0)
        #: z coordinate of center of gravity of
        #: nonstructural mass for end B.
        self.m2b[i] = double_or_blank(card, x + 11, 'm2b', self.m1b)

        #: y coordinate of neutral axis for end A.
        self.n1a[i] = double_or_blank(card, x + 12, 'n1a', 0.0)
        #: z coordinate of neutral axis for end A.
        self.n2a[i] = double_or_blank(card, x + 13, 'n2a', self.n1a)

        #: y coordinate of neutral axis for end B.
        self.n1b[i] = double_or_blank(card, x + 14, 'n1a', 0.0)
        #: z coordinate of neutral axis for end B.
        self.n2b[i] = double_or_blank(card, x + 15, 'n2b', self.n1b)
    #=========================================================================
    def get_index(self, property_ids):
        if isinstance(property_ids, int):
            property_ids = array([property_ids])
        if property_ids is None:
            return arange(self.n)

        indexs = searchsorted(self.property_id, property_ids)
        assert len(indexs) == len(property_ids), 'indexs=%s pids=%s' % (indexs, property_ids)
        return indexs

    #=========================================================================
    def write_bdf(self, f, size=8, property_ids=None):
        if self.n:
            if property_ids is None:
                i = arange(self.n)
            else:
                i = searchsorted(self.property_id, property_ids)

            for (pid, mid, area, I1, I2, J) in zip(self.property_id[i], self.material_id[i],
                    self.area[i], self.I1[i], self.I2[i], self.J[i]):
                card = ['PBEAM', pid, mid, area, I1, I2, J]
                f.write(print_card_8(card))
