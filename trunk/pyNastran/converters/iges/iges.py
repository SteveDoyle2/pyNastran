import pyNastran
from collections import defaultdict
from docopt import docopt

class IGES(object):
    supported  = [106, 110, 112, 114, 124, 126, 128, 142, 144,]
    maybe_supported = [100]
    def __init__(self):
        """
        http://www.transcendata.com/support/cadfix/faq/IGESimport.pdf
        http://www.okino.com/conv/imp_iges.htm

        Supported in VGRID
        ------------------
        100 - Circular Arc (???)
        102 - Composite Curve (???)
        106 - Copious Data
        110 - Line
        112 - Parametric Spline
        114 - Parametric Surface Spline
        116 - Point (???)
        124 - Matrix
        126 - NURBS curve (Rational B-Spline curve)
        128 - NURBS surface
        142 - Curve on a Parametric surface
        144 - Trimmed Surface

        Not Supported in VGRID (???)
        104 - Conic Arc
        108 - Bounded plane
        118 - Ruled Surface
        120 - Surface of Revolution
        122 - Tabulated Surface
        128 - Rational B-Spline surface
        143 - Bounded Surface
        186 - Manifold Solid

        124 - Tansformation Matrix
        308 - Subfigure definition
        314 - Color definition
        402 - Group associativity
        408 - Subfigure instance
        416 - External file reference

        130 - Offset Curve
        140 - Offset Surface
        141/143 - Bounded Surface

        150-169 - Constructive Solid Geometry entities
        180-186 - Manifold solid B-Rep entities
        """
        pass

    def read_iges(self, igs_name):
        f = open(igs_name)
        lines = f.readlines()
        f.close()

        'SolidWorks IGES file using analytic representation for surfaces         S      1'
        Start1 = []
        Global1 = []
        Directory1 = []
        Parameter1 = []
        Terminate1 = []
        for line in lines:
            SGDPT = line[72:73]
            if SGDPT in 'S':
                Start1.append(line)
            elif SGDPT in 'G':
                #print '%r' % line[:72]
                Global1.append(line[:72])
            elif SGDPT in 'D':
                Directory1.append(line)
            elif SGDPT in 'P':
                Parameter1.append(line)
            elif SGDPT in 'T':
                Terminate1.append(line)
            else:
                raise NotImplementedError(SGDPT)
            #print '%r' % line[72:73]
        del lines
        #print S1
        #for
        G_global = ''.join(Global1)
        print "G =", G_global.split(',')
        print "----P----"
        P_parameter_data = self.combine(Parameter1, True)
        #for key, line in sorted(P_parameter_data.iteritems()):
            #print key, line

        print "----D----"
        for line in Directory1:
            dline = line[0:72].split()
            print '%r' % dline
            dType = int(dline[0])
            assert dType in self.supported+self.maybe_supported, '%i is not supported' % dType

        #D_directory_entry = self.combine(Directory1)
        #for key, line in sorted(D_directory_entry.iteritems()):
            #print '***', key, line

    def combine(self, P1, debug=False):
        P = defaultdict(str)
        for pline in P1:
            i = pline[63:72].lstrip()
            #if debug:
                #print "i =", i
            #print '%r' % pline[63:72].lstrip()
            P[i] += pline[:70]
        return P
    def write_iges(self, fname):
        pass

def run_arg_parse():
    msg  = "Usage:\n"
    msg += "  iges (IGS_FILENAME)"
    ver = str(pyNastran.__version__)
    data = docopt(msg, version=ver)

    igs_filename = data['IGS_FILENAME']
    return igs_filename

def main():
    igs_name = run_arg_parse()
    igs_name2 = 'out.IGS'
    igs = IGES()
    igs.read_iges(igs_name)
    igs.write_iges(igs_name2)

if __name__ == '__main__':
    main()