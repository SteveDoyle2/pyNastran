
from numpy import arange, sqrt, linspace, where, zeros, array
from struct import pack, unpack
from pyNastran.op2.fortranFile import FortranFile
from pyNastran.converters.usm3d.usm3d_reader import Usm3dReader
from math import floor

def factors(nraw):
    result = []
    n = nraw
    for i in range(2,n+1): # test all integers between 2 and n
        s = 0;
        while n % i == 0: # is n/i an integer?
            n = n / float(i)
            s += 1
        if s > 0:
            for k in range(s):
                result.append(i) # i is a pf s times
            if n == 1:
                return result
    return [nraw]

class IFace(FortranFile):
    def __init__(self, log=None, debug=None):
        FortranFile.__init__(self)


    def read_poin1(self, poin1_filename):
        f = open(poin1_filename, 'r')
        ipoin1, ilines = f.readline().strip().split()
        ipoin1 = int(ipoin1)
        ilines = int(ilines)

        poin1 = zeros(ilines, 'int32')
        for i, line in enumerate(f.readlines()):
            ii, poin1i = line.split()
            ii = int(ii)
            assert i + 1 == ii, 'i=%s ii=%s' % (i+1, ii)
            poin1[i] = poin1i
        assert poin1.max() == ipoin1
        f.close()
        return poin1

    def read_m2(self, m2_filename):
        m2 = open(m2_filename)
        self.op2 = m2

        self.print_section2(5000, '>')

    def read_iface(self, iface_filename):
        """
        BC File...
        nFaces    nBouc,      nRegions, ???

        nFaces   - number of faces on the surface
        nBouc    - ???
        nRegions - number of independent surfaces that are set in the mapbc file
        ????     -

        #--------------------------------
        New2:

        BC File...
           56864     944       7          2
        nFaces     nBouc       nRegions ???

        Cogsg File...
            {'dummy': 6266912, 'nBoundPts': 28434,
            'nPoints': 79734,  'nElements': 391680,
            'nViscPts': 26304, 'nViscElem': 130560,
            'tc': 0.0, 'inew': -1, }

        IFace File...
        A=16235848 B=811792 C=56864
        A = [2, 2, 2, 1051, 1931]
        B = [2, 2, 2, 2, 113, 449]
        C = nfaces = [2, 2, 2, 2, 2, 1777]

        Flo File...
        nPoints = nlines = 79734

        --------------------------------------------------------------
        Box:

        BC File...
        6810     284       6          2
        nFaces   nBouc     nRegions ???

        Cogsg File...
            {'dummy': 1006496, 'nBoundPts': 3407,
            'nPoints': 12283, 'nElements': 62904,
            'nViscPts': 10907, 'nViscElem': 51564,
            'tc': 0.0, 'inew': -1,}

        Front File...
        npoit      nelet    npoif    nface    nboun    nfacs    nbouc    nfaci      npoiv      nelev
        12283      62904        0        0     3407     6810      284     6236      10907      51564

             6

        Poin1 File...
        3407       10907
       1       1
       2       2
       3       3
       4       4
       10907 lines

        IFace File...

        A=2584268 B=129213 C=6810
        # A = [2, 2, 646067]
        # B = [3, 3, 7, 7, 293]
        # C = [2, 3, 5, 227]
        A = ???
        B =
        C = nfaces
        """
        f = open(iface_filename, 'rb')

        # makes it so FortranFile works...
        self.op2 = f
        n = 0

        data = f.read(4 * 3)  # A, B, C
        A, B, C = unpack('>3i', data)

        print "A=%s B=%s C=%s" % (A, B, C)
        nints = C
        Format = '>%ii' % nints
        data = f.read(4 * nints)  # read nints ints
        ints = unpack(Format, data)
        self.n += 4 * nints

        assert max(ints) == ints[-1], 'max(ints)=%i ints[-1]=%i'  % (max(ints), ints[-1])

        # A=16235848 B=811792 C=56864
        # A = 2^3 * 3^4
        # B = 2^4 * 50737
        # C = 2^5 *

        # A=2584268 B=129213 C=6810
        # A = [2, 2, 646067]
        # B = [3, 3, 7, 7]
        # C = [2, 3, 5, 227]

        nints = B
        Format = '>%ii' % nints
        data = f.read(4 * nints)  # read nint ints
        ints = unpack(Format, data)
        assert max(ints) < nints, 'max(ints)=%i nints=%i'  % (max(ints), nints)
        self.n += 4 * nints
        #print ints

        #print factors(A)
        #print factors(B)
        #print factors(C)

        print self.print_section2(n, '>')
        pass

    def read_flo(self, flo_filename, n=None, nvars=5):
        """
        ipltqn is a format code where:
         - ipltqn = 0  (no printout)
         - ipltqn = 1  (unformatted)
         - ipltqn = 2  (formatted) - default

        nvars = 5
          - (nodeID,rho,rhoU,rhoV,rhoW) = sline
            (e) = line

        nvars = 6
          - (nodeID,rho,rhoU,rhoV,rhoW,e) = line
        """
        formatCode = 2

        f = open(flo_filename, 'r')
        mach = float(f.readline().strip())

        node_id = zeros(n, 'int32')
        rho = zeros(n, 'float32')
        rhoU = zeros(n, 'float32')
        rhoV = zeros(n, 'float32')
        rhoW = zeros(n, 'float32')
        e = zeros(n, 'float32')

        if n:
            if nvars == 6:
                for i in xrange(n):
                    sline1 = f.readline().strip().split()
                    #print sline1
                    rhoi, rhoui, rhovi, rhowi, ei = Float(sline1[1:], 5)
                    node_id[i] = sline1[0]
                    rho[i] = rhoi
                    rhoU[i] = rhoui
                    rhoV[i] = rhovi
                    rhoW[i] = rhowi
                    e[i] = ei
                    assert len(sline1) == 6, 'len(sline1)=%s' % len(sline1)
            else:
                for i in xrange(n):
                    sline1 = f.readline().strip().split()
                    node_id[i] = sline1[0]
                    rhoi, rhoui, rhovi, rhowi = Float(sline1[1:], 4)

                    assert len(line) == 5, 'len(sline1)=%s' % len(sline1)

                    sline2 = f.readline().strip().split()
                    ei = Float(sline2, 1)

                    rho[i] = rhoi
                    rhoU[i] = rhoui
                    rhoV[i] = rhovi
                    rhoW[i] = rhowi
                    e[i] = ei
                    assert len(sline2) == 1, 'len(sline2)=%s' % len(sline2)
        else:
            asdf

        f.close()
        gamma = 1.4
        gm1 = gamma - 1
        rhoVV = (rhoU**2+rhoV**2+rhoW**2) / rho
        if 'p' in result_names or 'Mach' in result_names:
            pND = gm1*(e - rhoVV/2. )
            if 'p' in result_names:
                loads['p'] = pND
        if 'Mach' in result_names:
            mach = (rhoVV/(gamma*pND))**0.5
            loads['Mach'] = Mach

        T = gamma*pND/rho # =a^2 as well
        if 'T' in result_names:
            loads['T'] = T
        return loads

def Float(sline, n):
    vals = []
    for val in sline:
        try:
            vals.append(float(val))
        except:
            vals.append(0.0)
    return vals

if __name__ == '__main__':
    cogsg_obj = Usm3dReader()
    model = 'new2'
    cogsg_obj.read_cogsg(model + '.cogsg')
    iface_obj = IFace()
    if model in ['box']:
        iface_obj.read_poin1(model + '.poin1')

    #iface_obj.read_m2(model + '.m2')
    #iface_obj.read_iface(model + '.iface')
    if model == 'new2':
        iface_obj.read_flo(model + '.flo', n=79734, nvars=6)
    elif model == 'box':
        iface_obj.read_flo(model + '.flo', n=12440, nvars=5)