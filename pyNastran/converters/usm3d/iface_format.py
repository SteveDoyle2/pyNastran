from struct import unpack
from numpy import zeros
#from pyNastran.converters.usm3d.usm3d_reader import Usm3d

def factors(nraw):
    """
    function for getting the primes factors of a number.

    This is supposed to help with figuring out the file format.
    I'm sure there's a better method, but it doesn't matter too much.
    """
    result = []
    n = nraw
    for i in range(2, n + 1): # test all integers between 2 and n
        s = 0
        while n % i == 0: # is n/i an integer?
            n = n / float(i)
            s += 1
        if s > 0:
            for k in range(s):
                result.append(i) # i is a pf s times
            if n == 1:
                return result
    return [nraw]


class IFace:
    def __init__(self, log=None, debug=None):
        self.n = 0
        self.log = log
        self.debug = debug

    def read_poin1(self, poin1_filename):
        with open(poin1_filename, 'r') as poin1_file:
            ipoin1, ilines = poin1_file.readline().strip().split()
            ipoin1 = int(ipoin1)
            ilines = int(ilines)

            poin1 = zeros(ilines, 'int32')
            for i, line in enumerate(poin1_file.readlines()):
                ii, poin1i = line.split()
                ii = int(ii)
                assert i + 1 == ii, 'i=%s ii=%s' % (i+1, ii)
                poin1[i] = poin1i
            assert poin1.max() == ipoin1
        return poin1

    #def read_m2(self, m2_filename):
        #m2 = open(m2_filename)
        #self.print_section2(5000, '>')

    def read_iface(self, iface_filename):
        """
        BC File...
        nFaces    nBouc,      nRegions, ???

        nFaces   - number of faces on the surface
        nBouc    - ???
        nRegions - number of independent surfaces that are set in the mapbc file
        ????     -

        **Cogsg File...**
        header = {
            'dummy': 6266912, 'nBoundPts': 28434,
            'nPoints': 79734,  'nElements': 391680,
            'nViscPts': 26304, 'nViscElem': 130560,
            'tc': 0.0, 'inew': -1,
        }
        """
        with open(iface_filename, 'rb') as iface_file:

            data = iface_file.read(4 * 3)  # A, B, C
            A, B, C = unpack('>3i', data)

            self.log.debug("A=%s B=%s C=%s" % (A, B, C))
            nints = C
            sformat = '>%ii' % nints
            data = iface_file.read(4 * nints)  # read nints ints
            ints = unpack(sformat, data)
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
            data = iface_file.read(4 * nints)  # read nint ints
            ints = unpack(Format, data)
            assert max(ints) < nints, 'max(ints)=%i nints=%i'  % (max(ints), nints)
            self.n += 4 * nints
            #print(ints

            #print(factors(A))
            #print(factors(B))
            #print(factors(C))

            #print(self.print_section2(n, '>'))


#if __name__ == '__main__':  # pragma: no cover
    #cogsg_obj = Usm3d()
    #model = 'new2'
    #cogsg_obj.read_cogsg(model + '.cogsg')
    #iface_obj = IFace()
    #if model in ['box']:
        #iface_obj.read_poin1(model + '.poin1')

    ##iface_obj.read_m2(model + '.m2')
    ##iface_obj.read_iface(model + '.iface')
    #if model == 'new2':
        #iface_obj.read_flo(model + '.flo', n=79734)
    #elif model == 'box':
        #iface_obj.read_flo(model + '.flo', n=12440)
