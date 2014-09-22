from numpy import array, zeros, arange, ravel, ones, cross
from numpy.linalg import norm

from pyNastran.converters.shabp.shabp_results import ShabpOut

class SHABP(ShabpOut):
    def __init__(self, log=None, debug=False):
        ShabpOut.__init__(self, log, debug)
        self.X = {}
        self.Y = {}
        self.Z = {}
        self.trailer = None

        self.component_name_to_patch = {}
        self.patch_to_component_num = {}
        self.component_to_params = {}
        self.component_num_to_name = {}
        self.component_name_to_num = {}

    def write_shabp(self, out_filename):
        pass

    def get_area_xlength_by_component(self, components=None):
        if components is None:
            components = self.component_name_to_patch.keys()
        ncomponents = len(components)

        # we're not using a vector because component name is
        # probably more useful
        areas = {}
        lengths = {}
        xmax = None
        xmin = None
        for name in components:
            patches = self.component_name_to_patch[name]

            A = 0.
            for i, ipatch in enumerate(patches):
                X = self.X[ipatch-1]
                nrows, ncols = X.shape

                npoints = nrows * ncols
                XYZ = zeros((npoints, 3), dtype='float32')
                XYZ[:, 0] = self.X[ipatch-1].ravel()
                XYZ[:, 1] = self.Y[ipatch-1].ravel()
                XYZ[:, 2] = self.Z[ipatch-1].ravel()

                if xmin is None:
                    xmin = XYZ[:, 0].min()
                    xmax = XYZ[:, 0].max()
                else:
                    xmin = min(xmin, XYZ[:, 0].min())
                    xmax = max(xmax, XYZ[:, 0].max())

                # TODO: can we vectorize this efficiently?
                for irow in xrange(nrows-1):
                    for jcol in xrange(ncols-1):
                        i1 = irow*ncols +jcol,
                        i2 = (irow+1)*ncols +(jcol  ),
                        i3 = (irow+1)*ncols +(jcol+1),
                        i4 = irow*ncols +(jcol+1),
                        a = XYZ[i3, :] - XYZ[i1, :]
                        b = XYZ[i4, :] - XYZ[i2, :]
                        A += 0.5 * norm(cross(a, b))

            areas[name] = A
            lengths[name] = xmax - xmin
        print "areas =", areas
        print "lengths =", lengths
        return areas

    def get_area_by_component(self, components=None):
        if components is None:
            components = self.component_name_to_patch.keys()
        ncomponents = len(components)

        # we're not using a vector because component name is
        # probably more useful
        areas = {}
        for name in components:
            patches = self.component_name_to_patch[name]
            A = self.get_area_by_patch(patches)
            areas[name] = A.sum()
        print "areas =", areas
        return areas

    def get_area_by_patch(self, ipatches=None):
        if ipatches is None:
            ipatches = arange(ipatches)

        areas = zeros(len(ipatches), dtype='float64')
        for i, ipatch in enumerate(ipatches):
            X = self.X[ipatch-1]
            nrows, ncols = X.shape

            npoints = nrows * ncols
            XYZ = zeros((npoints, 3), dtype='float32')
            XYZ[:, 0] = self.X[ipatch-1].ravel()
            XYZ[:, 1] = self.Y[ipatch-1].ravel()
            XYZ[:, 2] = self.Z[ipatch-1].ravel()

            A = 0.
            # TODO: can we vectorize this efficiently?
            for irow in xrange(nrows-1):
                for jcol in xrange(ncols-1):
                    i1 = irow*ncols +jcol,
                    i2 = (irow+1)*ncols +(jcol  ),
                    i3 = (irow+1)*ncols +(jcol+1),
                    i4 = irow*ncols +(jcol+1),
                    a = XYZ[i3, :] - XYZ[i1, :]
                    b = XYZ[i4, :] - XYZ[i2, :]
                    A += 0.5 * norm(cross(a, b))
            areas[i] = A
        return areas

    def read_shabp(self, infilename):
        lines = open(infilename).readlines()
        if(infilename.lower().endswith(".geo")):
            i = 0
        else:  # this supports standard .inp and .mk5 files
            i = 3
        self.header = lines[0:i]

        patches = []
        while i < len(lines)-1:
            patch = []
            # C  1001         10        3  patch 1
            header = lines[i]
            #print "%r" % header.strip()
            name = header[:7].strip()
            #num = int(header[4:7])
            nrows = int(header[7:20])
            ncols = int(header[20:29])
            name2 = header[29:].strip()
            #print lines[i].strip()
            flag = name[4]
            if flag == '1':
                is_last_patch = True
            elif flag == '0':
                is_last_patch = False
            else:
                raise RuntimeError('last patch flag = %r; must be 0 or 1' % flag)
            #print "name=%r nrows=%i ncols=%i name2=%r" % (name, nrows, ncols, name2)

            i += 2
            row = []
            while 1:
                #-1071.0480   77.2500  -66.94202 -987.7440   77.2500  -66.94200         3
                line = lines[i]
                i += 1

                x1,y1,z1,t1 = line[:10], line[10:20],line[20:30],line[30:31]
                x1 = float(x1)
                y1 = float(y1)
                z1 = float(z1)
                t1 = int(t1)
                if t1 == 1:
                    patch.append(row)
                    row = []
                    row.append([x1, y1, z1])
                elif t1 in [0, 2]:
                    row.append([x1, y1, z1])
                    pass
                elif t1 == 3:
                    row.append([x1, y1, z1])
                    patch.append(row)
                    #print "breaking1"
                    break
                else:
                    aaaa

                x2,y2,z2,t2 = line[31:41], line[41:51],line[51:61],line[61:62]
                x2 = float(x2)
                y2 = float(y2)
                z2 = float(z2)
                t2 = int(t2)
                if t2 == 1:
                    patch.append(row)
                    row = []
                    row.append([x2, y2, z2])
                elif t2 in [0, 2]:
                    row.append([x2, y2, z2])
                    pass
                elif t2 == 3:
                    row.append([x2, y2, z2])
                    patch.append(row)
                    #print "breaking2"
                    break
                else:
                    bbb

            patches.append(array(patch, dtype='float32'))
            if is_last_patch:
                #print "***last patch - lines[%i] = %r" % (i, lines[i])
                break
            #print "lines[%i] = %r" % (i, lines[i])
        #print "********len(patches) =", len(patches)
        self.trailer = lines[i:]
        self.build_patches(patches)
        self.parse_trailer()

    def build_patches(self, patches):
        X = []
        Y = []
        Z = []
        for ipatch, patch in enumerate(patches):
            nrows = len(patch)
            ncols = len(patch[0])
            x = zeros((nrows, ncols), dtype='float32')
            y = zeros((nrows, ncols), dtype='float32')
            z = zeros((nrows, ncols), dtype='float32')
            for irow, row in enumerate(patch):
                for icol, col in enumerate(row):
                    x[irow, icol] = col[0]
                    y[irow, icol] = col[1]
                    z[irow, icol] = col[2]
            X.append(x)
            Y.append(y)
            Z.append(z)
        self.X = X
        self.Y = Y
        self.Z = Z
        print "npatches = ", len(self.X)


    def getPointsElementsRegions(self):
        npatches = len(self.X)
        npoints = 0
        nelements = 0
        for ipatch in xrange(npatches):
            X = self.X[ipatch]
            nrows, ncols = X.shape
            npoints += nrows * ncols
            nelements += (nrows-1) * (ncols-1)

        ipoint = 0
        ielement = 0
        j = 0
        XYZ = zeros((npoints, 3), dtype='float32')
        elements2 = zeros((nelements, 4), dtype='float32')
        components = ones(nelements, dtype='int32')
        patches = ones(nelements, dtype='int32')
        impact = ones(nelements, dtype='int32')
        shadow = ones(nelements, dtype='int32')

        for ipatch in xrange(npatches):
            #print "---ipatch=%s---" % ipatch
            if ipatch not in self.patch_to_component_num:
                comp_num = 0
                #raise RuntimeError('ipatch=%s keys=%s' % (ipatch, sorted(self.patch_to_component_num)))
                impact_val = 0
                shadow_val = 0
                #continue
            else:
                comp_num = self.patch_to_component_num[ipatch]
                data = self.component_to_params[comp_num-1]
                impact_val = data[0]
                shadow_val = data[1]
            nrows, ncols = self.X[ipatch].shape
            npointsi = nrows * ncols
            nelementsi = (nrows-1) * (ncols-1)
            XYZ[ipoint:ipoint+npointsi, 0] = self.X[ipatch].ravel()
            XYZ[ipoint:ipoint+npointsi, 1] = self.Y[ipatch].ravel()
            XYZ[ipoint:ipoint+npointsi, 2] = self.Z[ipatch].ravel()

            #if comp_num == 0:
                #continue
            elements = []

            for irow in xrange(nrows-1):
                for jcol in xrange(ncols-1):
                    element = [
                         irow*ncols + jcol,
                         (irow+1)*ncols + (jcol  ),
                         (irow+1)*ncols + (jcol+1),
                         irow*ncols + (jcol+1),
                         ]
                    elements.append(element)
            elements = array(elements, dtype='int32')

            patches[ielement:ielement+nelementsi] *= (ipatch+1)
            components[ielement:ielement+nelementsi] *= comp_num
            impact[ielement:ielement+nelementsi] *= impact_val
            shadow[ielement:ielement+nelementsi] *= shadow_val
            elements2[ielement:ielement+nelementsi,:] = elements[:,:] + ipoint
            #print "  ipatch=%i Cp[%i:%i]" % (ipatch+1, ielement, ielement+nelementsi)
            ipoint += npointsi
            ielement += nelementsi
        return XYZ, elements2, patches, components, impact, shadow

    def parse_trailer(self):
        #for line in self.trailer:
            #print line.rstrip()
        self.title = self.trailer[0].strip()
        line2 = self.trailer[1]
        npatches = line2[:2]

        mach_line = self.trailer[2].rstrip()
        mach = float(mach_line[0 :10].strip())
        alt = float(mach_line[10:20].strip())
        pstag = float(mach_line[20:30].strip())
        tstag = float(mach_line[30:40].strip())
        igas = int(mach_line[40:41].strip())
        nalpha_beta = int(mach_line[41:43].strip())

        ref_line = self.trailer[3].rstrip()
        Sref = float(ref_line[0:10].strip())
        Cref = float(ref_line[10:20].strip())
        Bref = float(ref_line[20:30].strip())
        Xcg = float(ref_line[30:40].strip())
        Ycg = float(ref_line[40:50].strip())
        Zcg = float(ref_line[50:60].strip())

        i = 4
        for n in xrange(nalpha_beta):
            alpha_line = self.trailer[i+n].rstrip()
            #print "alpha_line =", alpha_line
            i += 1

        #self.getMethods()
        ncomponents_line = self.trailer[i].rstrip()
        #print "comp line =", ncomponents_line
        ncomponents = int(ncomponents_line[:2])
        zero = ncomponents_line[2:3]
        assert zero == '0', 'zero=%r; %s' %(zero, ncomponents_line)

        i += 1
        for icomponent in xrange(ncomponents):
            #print "lines[%i] = %s" % (i, self.trailer[i].rstrip())
            #print "lines[%i] = %s\n" % (i+1, self.trailer[i+1].rstrip())

            line2 = self.trailer[i+1].rstrip()
            impact = int(line2[0:2].strip())
            shadow = int(line2[2:4].strip())
            iprint = int(line2[4:5].strip())
            ipin   = int(line2[5:6].strip())
            isave  = int(line2[6:7].strip())
            pdata1 = float(line2[10:20].strip())
            pdata2 = float(line2[20:30].strip())
            pdata3 = float(line2[30:40].strip())
            pdata4 = float(line2[40:50].strip())
            pdata5 = float(line2[50:60].strip())
            pdata6 = float(line2[60:70].strip())
            self.component_to_params[icomponent] = [impact, shadow, iprint,
                ipin, isave, pdata1, pdata2, pdata3, pdata4, pdata5, pdata6]
            i += 2

        print "ncomps =", len(self.component_to_params)
        print "keys =", sorted(self.component_to_params.keys())
        print "**lines[%i] = %s\n" % (i+1, self.trailer[i].rstrip())
        #i += 1

        methods = []
        print "methods ", self.trailer[i].strip()
        for v in self.trailer[i].strip():
            v2 = int(v)
            if v2 > 0:
                methods.append(v2)
        print "methods =", methods

        val1 = self.trailer[i][0:2]
        print "val1 = ", val1
        i += 1
        if val1 == '10':
            print "****10****"
            pass
            i += 1
        elif val1 == '13':
            print "****13****"
            val2 = int(self.trailer[i][2:4])
            assert val2 == ncomponents, 'val2=%r ncomponents=%r' % (val2, ncomponents)
            i += 1
            for v in xrange(val2):
                print "  lines[%i] = %s" % (i+1, self.trailer[i].rstrip())
                i += 1
        else:
            print "*lines[%i] = %s\n" % (i+1, self.trailer[i].rstrip())
            aaaa

        # component names   7
        comp_names_line = self.trailer[i].rstrip()
        print "comp_names_line = %r" % comp_names_line
        ncomps = int(comp_names_line.strip().split()[2])
        i += 1

        for icomp in xrange(ncomps):
            line = self.trailer[i]
            print "line =", line.strip()
            npatches = int(line[:2])
            name = line[2:].strip()

            line = self.trailer[i+1] + ' '
            patches = []
            for n in xrange(0, len(line), 2):
                #print "n =", n
                ipatch = line[n:n+2].strip()
                if len(ipatch) == 0:
                    break
                int_ipatch = int(ipatch)
                #print "int_ipatch =", int_ipatch
                patches.append(int_ipatch)
                self.patch_to_component_num[int_ipatch-1] = icomp+1
            print "patches =", patches, '\n'
            self.component_name_to_patch[name] = patches
            self.component_num_to_name[icomp] = name
            self.component_name_to_num[name] = icomp
            i += 2

        # 2noseconeright
        #3142
        print 'done with trailer'

if __name__ == '__main__':
    import sys
    s = SHABP()
    s.read_shabp(sys.argv[1])
    s.read_shabp_out('SHABP.OUT')
    #s.get_area_by_component()
    #s.get_area_xlength_by_component()

