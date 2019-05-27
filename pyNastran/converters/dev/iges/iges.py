from collections import defaultdict
from docopt import docopt
import pyNastran

class IGES:
    supported = [106, 110, 112, 114, 124, 126, 128, 142, 144,]
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
        with open(igs_name, 'r') as iges_file:
            lines = iges_file.readlines()

        # 'SolidWorks IGES file using analytic representation for surfaces         S      1'
        start1 = []
        global1 = []
        directory1 = []
        parameter1 = []
        terminate1 = []
        for line in lines:
            sgdpt = line[72:73]
            if sgdpt in 'S':
                start1.append(line)
            elif sgdpt in 'G':
                #print('%r' % line[:72])
                global1.append(line[:72])
            elif sgdpt in 'D':
                directory1.append(line)
            elif sgdpt in 'P':
                parameter1.append(line)
            elif sgdpt in 'T':
                terminate1.append(line)
            else:
                raise NotImplementedError(sgdpt)
            #print('%r' % line[72:73])
        del lines
        #print(S1)
        #for
        g_global = ''.join(global1)
        print("G =", g_global.split(','))
        print("----P----")
        p_parameter_data = self.combine(parameter1, True)
        #for key, line in sorted(p_parameter_data.items()):
            #print(key, line)

        print("----D----")
        for line in directory1:
            dline = line[0:72].split()
            print('%r' % dline)
            dtype = int(dline[0])
            assert dtype in self.supported+self.maybe_supported, '%i is not supported' % dtype

        #D_directory_entry = self.combine(Directory1)
        #for key, line in sorted(D_directory_entry.items()):
            #print('***', key, line)

    def combine(self, input_dict, debug=False):
        combined_dict = defaultdict(str)
        for pline in input_dict:
            i = pline[63:72].lstrip()
            #if debug:
                #print("i = %s" % i)
            #print('%r' % pline[63:72].lstrip())
            combined_dict[i] += pline[:70]
        return combined_dict

    def write_iges(self, fname):
        pass

def run_arg_parse():
    msg = "Usage:\n"
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

if __name__ == '__main__':  # pragma: no cover
    main()
