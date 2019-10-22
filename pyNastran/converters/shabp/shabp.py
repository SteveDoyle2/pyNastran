"""
dfeines:
 - read_shabp(shabp_filename, log=None, debug=False)
 - SHABP(log=None, debug=False)

"""
from numpy import array, zeros, arange, ones, cross
from numpy.linalg import norm  # type: ignore
from cpylog import get_logger2

from pyNastran.converters.shabp.shabp_results import ShabpOut
#from pyNastran.converters.shabp.parse_trailer import parse_trailer

def read_shabp(shabp_filename, read_special_routines=False, log=None, debug=False):
    """reads an S/HABP file"""
    model = SHABP(log=log, debug=debug)
    model.read_shabp(shabp_filename, read_special_routines=read_special_routines)
    return model

class SHABP(ShabpOut):
    """defines the SHABP class"""
    def __init__(self, log=None, debug=False):
        """
        Initializes the SHABP object

        Parameters
        ----------
        debug : bool/None; default=True
            used to set the logger if no logger is passed in
                True:  logs debug/info/error messages
                False: logs info/error messages
                None:  logs error messages
        log : logging module object / None
            if log is set, debug is ignored and uses the
            settings the logging object has

        """
        #self.xyz = {}
        self.X = {}
        self.Y = {}
        self.Z = {}
        self.trailer = None

        self.component_name_to_patch = {}
        self.patch_to_component_num = {}
        self.component_to_params = {}
        self.component_num_to_name = {}
        self.component_name_to_num = {}
        self.log = get_logger2(log, debug=debug)

        self.title = ''
        self.header = ''
        self.shabp_cases = {}

    #def write_shabp(self, out_filename):
        #"""writes an S/HABP file"""
        #pass

    def get_area_xlength_by_component(self, components=None):
        """gets the area and length of a set of components"""
        if components is None:
            components = self.component_name_to_patch.keys()
        #ncomponents = len(components)

        # we're not using a vector because component name is
        # probably more useful
        areas = {}
        lengths = {}
        xmax = None
        xmin = None
        for name in components:
            patches = self.component_name_to_patch[name]

            area = 0.
            for unused_i, ipatch in enumerate(patches):
                X = self.X[ipatch-1]
                nrows, ncols = X.shape

                npoints = nrows * ncols
                xyz = zeros((npoints, 3), dtype='float32')
                xyz[:, 0] = self.X[ipatch-1].ravel()
                xyz[:, 1] = self.Y[ipatch-1].ravel()
                xyz[:, 2] = self.Z[ipatch-1].ravel()

                if xmin is None:
                    xmin = xyz[:, 0].min()
                    xmax = xyz[:, 0].max()
                else:
                    xmin = min(xmin, xyz[:, 0].min())
                    xmax = max(xmax, xyz[:, 0].max())

                # TODO: can we vectorize this efficiently?
                for irow in range(nrows-1):
                    for jcol in range(ncols-1):
                        i1 = irow*ncols +jcol,
                        i2 = (irow+1)*ncols + jcol,
                        i3 = (irow+1)*ncols + (jcol+1),
                        i4 = irow*ncols +(jcol+1),
                        a = xyz[i3, :] - xyz[i1, :]
                        b = xyz[i4, :] - xyz[i2, :]
                        area += 0.5 * norm(cross(a, b))

            areas[name] = area
            lengths[name] = xmax - xmin
        return areas, lengths

    def get_area_by_component(self, components=None):
        """gets the area of a set of components"""
        if components is None:
            components = self.component_name_to_patch.keys()
        #ncomponents = len(components)

        # we're not using a vector because component name is
        # probably more useful
        areas = {}
        for name in components:
            patches = self.component_name_to_patch[name]
            area = self.get_area_by_patch(patches)
            areas[name] = area.sum()
        return areas

    def get_area_by_patch(self, ipatches=None):
        """gets the area of a set of patches"""
        if ipatches is None:
            npatches = len(self.X)
            ipatches = arange(npatches)

        areas = zeros(len(ipatches), dtype='float32')
        for i, ipatch in enumerate(ipatches):
            X = self.X[ipatch-1]
            nrows, ncols = X.shape

            npoints = nrows * ncols
            xyz = zeros((npoints, 3), dtype='float32')
            xyz[:, 0] = self.X[ipatch-1].ravel()
            xyz[:, 1] = self.Y[ipatch-1].ravel()
            xyz[:, 2] = self.Z[ipatch-1].ravel()

            area = 0.
            # TODO: can we vectorize this efficiently?
            for irow in range(nrows-1):
                for jcol in range(ncols-1):
                    i1 = irow*ncols +jcol,
                    i2 = (irow+1)*ncols +(jcol),
                    i3 = (irow+1)*ncols +(jcol+1),
                    i4 = irow*ncols + (jcol+1),
                    a = xyz[i3, :] - xyz[i1, :]
                    b = xyz[i4, :] - xyz[i2, :]
                    area += 0.5 * norm(cross(a, b))
            areas[i] = area
        return areas

    def read_shabp(self, shabp_filename, read_special_routines=False):
        """reads an SHABP.INP / SHABP.mk5 file"""
        with open(shabp_filename) as shabp_file:
            lines = shabp_file.readlines()

        if shabp_filename.lower().endswith(".geo"):
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

            if 1:
                name = header[:7].strip()
                #num = int(header[4:7])
                #nrows = int(header[7:20])
                #ncols = int(header[20:29])
                #name2 = header[29:].strip()
                #print lines[i].strip()
                flag = name[4]
                self.log.debug(f'name={name} flag={flag}')
            else:
                stream_cross = header[5]       # heat/tk A-10
                symmetry = header[16]          # heat/tk A-10
                scale_factor = header[17]      # heat/tk A-10
                #nadj1 = header[18:20]  # not used
                #nadj2 = header[20:22]  # not used
                #nadj3 = header[22:24]  # not used
                #nadj4 = header[24:26]  # not used
                self.log.debug(f'stream_cross={stream_cross} symmetry={symmetry} scale_factor={scale_factor}')

                vis_type = int(header[26])
                #if vis_type == 0:
                    #print("inviscid_A")
                #if vis_type == 1:
                    #print("neither")
                #if vis_type == 2:
                    #print("viscous_A")
                #if vis_type == 3:
                    #print("inviscid_B")

                name_old = header[:7].strip()
                self.log.info(f'name={name!r} stream_cross={stream_cross!r} symmetry={symmetry!r} '
                              f'scale_factor={scale_factor!r} vis_type={vis_type!r}')

                header2 = lines[i+1]
                xsc = header2[0:10].strip()
                ysc = header2[10:20].strip()
                zsc = header2[20:30].strip()

                delx = header2[30:40].strip()
                dely = header2[40:50].strip()
                delz = header2[50:60].strip()
                self.log.info('sc=%s,%s,%s del=%r,%r,%r' % (xsc, ysc, zsc, delx, dely, delz))
                #print(lines[i].strip())
                flag = name_old[4]

            if flag == '1':
                is_last_patch = True
            elif flag == '0':
                is_last_patch = False
            else:
                raise RuntimeError('last patch flag = %r; must be 0 or 1' % flag)
            #print(f'name={name} nrows={nrows} ncols={ncols:d} name2={name2!r}')

            i += 2
            row = []
            while 1:
                #-1071.0480   77.2500  -66.94202 -987.7440   77.2500  -66.94200         3
                line = lines[i]
                #print(i, line, end='')
                i += 1

                #t1 = int(line[30])    # STAT - status flag - p. 27
                x1, y1, z1, t1 = line[:10], line[10:20], line[20:30], line[30:31]
                t1 = line[30]
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
                elif t1 == 3:
                    row.append([x1, y1, z1])
                    patch.append(row)
                    #assert len(row) == nrows, 'len(row)=%s nrows=%s ncols=%s' % (len(row), nrows, ncols)
                    break
                else:
                    raise RuntimeError()

                if t1 == 3:
                    break
                x2, y2, z2, t2 = line[31:41], line[41:51], line[51:61], line[61:62]
                t2 = line[61]
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
                elif t2 == 3:
                    row.append([x2, y2, z2])
                    patch.append(row)
                    break
                else:
                    raise RuntimeError()

            try:
                patchi = array(patch, dtype='float32')
            except:
                print('patch =', patch)
                for i, patchi in enumerate(patch):
                    print('i=%s n=%s patch=%s' % (i, len(patchi), patchi))
                raise

            patches.append(patchi)
            if is_last_patch:
                #print "***last patch - lines[%i] = %r" % (i, lines[i])
                break
            #print "lines[%i] = %r" % (i, lines[i])
        self.trailer = lines[i:]
        self.build_patches(patches)
        try:
            self.parse_trailer(read_special_routines)
        except (RuntimeError, ValueError):
            self.log.error('failed parsing trailer')
            #raise

    def build_patches(self, patches):
        X = []
        Y = []
        Z = []
        #XYZ = []
        for patch in patches:
            nrows = len(patch)
            ncols = len(patch[0])
            #xyz = zeros((nrows, ncols, 3), dtype='float32')
            x = zeros((nrows, ncols), dtype='float32')
            y = zeros((nrows, ncols), dtype='float32')
            z = zeros((nrows, ncols), dtype='float32')
            for irow, row in enumerate(patch):
                for icol, col in enumerate(row):
                    x[irow, icol] = col[0]
                    y[irow, icol] = col[1]
                    z[irow, icol] = col[2]
                    #xyz[irow, icol, :] = [col[0], col[1], col[2]]
            X.append(x)
            Y.append(y)
            Z.append(z)
            #XYZ.append(xyz)
        self.X = X
        self.Y = Y
        self.Z = Z
        #self.XYZ = XYZ

    def get_points_elements_regions(self):
        npatches = len(self.X)
        npoints = 0
        nelements = 0
        for ipatch in range(npatches):
            X = self.X[ipatch]
            nrows, ncols = X.shape
            npoints += nrows * ncols
            nelements += (nrows-1) * (ncols-1)

        ipoint = 0
        ielement = 0
        xyz = zeros((npoints, 3), dtype='float32')
        elements2 = zeros((nelements, 4), dtype='int32')
        components = ones(nelements, dtype='int32')
        patches = ones(nelements, dtype='int32')
        impact = ones(nelements, dtype='int32')
        shadow = ones(nelements, dtype='int32')

        for ipatch in range(npatches):
            comp_num, impact_val, shadow_val = self.get_impact_shadow(ipatch)

            nrows, ncols = self.X[ipatch].shape
            npointsi = nrows * ncols
            nelementsi = (nrows-1) * (ncols-1)
            xyz[ipoint:ipoint+npointsi, 0] = self.X[ipatch].ravel()
            xyz[ipoint:ipoint+npointsi, 1] = self.Y[ipatch].ravel()
            xyz[ipoint:ipoint+npointsi, 2] = self.Z[ipatch].ravel()

            #if comp_num == 0:
                #continue
            elements = []

            for irow in range(nrows-1):
                for jcol in range(ncols-1):
                    element = [
                        irow*ncols + jcol,
                        (irow+1)*ncols + (jcol),
                        (irow+1)*ncols + (jcol+1),
                        irow*ncols + (jcol+1),
                    ]
                    elements.append(element)
            elements = array(elements, dtype='int32')

            patches[ielement:ielement+nelementsi] *= (ipatch+1)
            components[ielement:ielement+nelementsi] *= comp_num
            impact[ielement:ielement+nelementsi] *= impact_val
            shadow[ielement:ielement+nelementsi] *= shadow_val
            elements2[ielement:ielement+nelementsi, :] = elements[:, :] + ipoint
            #print("  ipatch=%i Cp[%i:%i]" % (ipatch+1, ielement, ielement+nelementsi))
            ipoint += npointsi
            ielement += nelementsi
        return xyz, elements2, patches, components, impact, shadow

    def get_impact_shadow(self, ipatch):
        """We may not have inviscid pressure"""
        comp_num = 0
        impact_val = 0
        shadow_val = 0
        if ipatch not in self.patch_to_component_num:
            self.log.debug(f'skipping ipatch={ipatch} comp_num={comp_num}')
            return comp_num, impact_val, shadow_val
            #raise RuntimeError('ipatch=%s keys=%s' % (
                #ipatch, sorted(self.patch_to_component_num)))

        comp_num = self.patch_to_component_num[ipatch]
        if comp_num not in self.component_num_to_name:
            self.log.debug(f'skipping ipatch={ipatch} comp_num={comp_num}')
            return comp_num, impact_val, shadow_val

        name = self.component_num_to_name[comp_num]
        patch = self.component_name_to_patch[name]
        if comp_num-1 not in self.component_to_params:
            self.log.debug(f'skipping ipatch={ipatch} comp_num={comp_num} name={name!r} patch={patch}')
            return comp_num, impact_val, shadow_val

        data = self.component_to_params[comp_num-1]
        #except KeyError:
            #name = self.component_num_to_name[comp_num]
            #self.log.error(f'ipatch={ipatch} comp_num={comp_num} name={name}')
            #print(self.component_to_params)
            #raise
        impact_val = data[0]
        shadow_val = data[1]
        self.log.debug(f'loading ipatch={ipatch} comp_num={comp_num} name={name!r} patch={patch}')
        return comp_num, impact_val, shadow_val

    def parse_trailer(self, read_special_routines=False):
        """parses the case information (e.g., number of angles of attack, control surface info)"""
        #out = parse_trailer(self.trailer)
        #order, component_names, cases, components = out
        #print('order = %s' % order)
        #print('component_names = %s' % component_names)
        #print('cases = %s' % cases)
        #print('components = %s' % components)
        #return out

        #for line in self.trailer:
            #print line.rstrip()
        self.title = self.trailer[0].strip()
        #print('title = %r' % self.title)
        line2 = self.trailer[1]

        #print('line2 = %r' % line2[:20])
        methods, summation_flag = _get_methods(line2)
        #print('methods =', methods)
        #print('summation_flag =', summation_flag)

        unused_npatches = line2[:2]
        #print('npatches =', npatches)

        mach_line = self.trailer[2].rstrip()
        mach, alt, pstag, tstag, igas, nalpha_beta, ideriv, dangle = _parse_flight_condition(mach_line)
        #print(f'mach={mach} alt={alt} pstag={pstag} tstag={tstag} igas={igas} '
              #f'nalpha_beta={nalpha_beta} ideriv={ideriv} dangle={dangle}')

        ref_line = self.trailer[3].rstrip()
        unused_reference_quantities = _parse_reference_conditions(ref_line)
        #sref, cref, bref, xcg, ycg, zcg = reference_quantities
        #print(f'sref={sref} cref={cref} bref={bref} cg=[{xcg},{ycg},{zcg}]')

        i = 4
        for n in range(nalpha_beta):
            alpha_line = self.trailer[i].rstrip()
            alpha, beta, roll, cdelta, qi, ri, pi = _read_alpha_beta_line(alpha_line)
            #print(f'alpha={alpha} beta={beta} roll={roll} cdelta={cdelta} pqr=[{qi},{ri},{pi}]')

            self.shabp_cases[n] = [mach, alpha, beta]
            i += 1

        #self.getMethods()
        self.log.debug(f'methods = {methods}')
        for method in methods:
            line = self.trailer[i].rstrip()
            if method == 3:
                self.log.debug(f'inviscid_pressure: {line}')
                i, ncomponents, ifsave, params = parse_inviscid_pressure(self.trailer, line, i)
                self.component_to_params = params
            elif method == 4:
                i = parse_viscous(self.trailer, line, i)
            elif method == 5:
                i = parse_special_routines(self.trailer, line, i,
                                           read_special_routines=read_special_routines)
            else:
                raise NotImplementedError(method)

        #print "ncomps =", len(self.component_to_params)
        #print "keys =", sorted(self.component_to_params.keys())
        #print "**lines[%i] = %s\n" % (i+1, self.trailer[i].rstrip())
        #i += 1

        # component names   7
        comp_names_line = self.trailer[i].rstrip()
        #print("comp_names_line = %r" % comp_names_line)
        ncomps = int(comp_names_line.strip().split()[2])
        i += 1

        for icomp in range(ncomps):
            line = self.trailer[i]
            #print("line =", line.strip())
            unused_npatches = int(line[:2])
            name = line[2:].strip()

            line = self.trailer[i+1] + ' '
            patches = []
            for n in range(0, len(line), 2):
                #print("n =", n)
                ipatch = line[n:n+2].strip()
                if len(ipatch) == 0:
                    break
                int_ipatch = int(ipatch)
                #print("int_ipatch =", int_ipatch)
                patches.append(int_ipatch)
                self.patch_to_component_num[int_ipatch-1] = icomp+1
            #print("patches =", patches, '\n')
            self.component_name_to_patch[name] = patches
            self.component_num_to_name[icomp] = name
            self.component_name_to_num[name] = icomp
            i += 2

        # 2noseconeright
        #3142
        self.log.debug('done with trailer')


def _parse_flight_condition(mach_line):
    """reads a flight condition line"""
    mach = float(mach_line[0 :10].strip())
    alt = float(mach_line[10:20].strip())
    pstag = float(mach_line[20:30].strip())
    tstag = float(mach_line[30:40].strip())
    igas = int(mach_line[40:41].strip())
    nalpha_beta = int(mach_line[41:43].strip())

    # 0-no derivatives
    # 1: pitch static stability derivatives
    # 2: pitch control derivatives
    # 3: pitch dynamic derivatives
    # 4: lat/dir static stability derivatives
    # 5: lat/dir control derivatives
    # 6: lat/dir dynamic derivatives
    # 7: ideriv=1,2,3 (all pitch)
    # 8: ideriv=4,5,6 (all lat/dir)
    # 9: ideriv=1-6 (all)
    try:
        ideriv = int(mach_line[43].strip()) # 0-9
        dangle = float(mach_line[45:55].strip())  # used for ideriv=1,2,4,5
    except:  # old SHABP
        raise # old SHABP
        #iDeriv = 0
        #dAngle = 0.0
    assert igas in [0, 1], igas # 0=air, 1=helium
    return mach, alt, pstag, tstag, igas, nalpha_beta, ideriv, dangle

def _get_methods(line):
    """
    per vecc A-53
    1 - Flow Field Analysis
    2 - Shielding Analysis
    3 - Inviscid Pressures
    4 - Viscous Methods
    5 - Special Routines
    """
    fields = line[:20]
    assert len(fields) == 20, fields
    methods = []
    for inti in line[:20]:
        if inti == '0':
            continue
        methodi = int(inti)
        methods.append(methodi)
        summation_flag = int(line[20])
    return methods, summation_flag

def _read_alpha_beta_line(alpha_line):
    """reads the alpha line"""
    #print("alpha_line =", alpha_line)
    alpha = float(alpha_line[0:10])
    beta = float(alpha_line[10:20])
    roll = float(alpha_line[20:30])
    cdelta = float(alpha_line[30:40]) # unused
    qi = float(alpha_line[40:50])
    ri = float(alpha_line[50:60])
    #pi = float(alpha_line[60:70])
    pi = None
    return alpha, beta, roll, cdelta, qi, ri, pi

def _parse_reference_conditions(ref_line):
    """reads the reference quantities line"""
    sref = float(ref_line[0:10].strip())
    cref = float(ref_line[10:20].strip())
    bref = float(ref_line[20:30].strip())
    xcg = float(ref_line[30:40].strip())
    ycg = float(ref_line[40:50].strip())
    zcg = float(ref_line[50:60].strip())
    return sref, cref, bref, xcg, ycg, zcg

def parse_inviscid_pressure(lines, line, i):
    """reads inviscid pressure (method=3)"""
    ncomponents = int(line[0:2])
    # 0 : save component force data; use new Unit 9
    # 1 : save component force data; use old Unit 9
    # 2 : dont save force data
    ifsave = int(line[2:3])
    unused_title = line[6:66].strip()
    assert ifsave == 0, 'ifsave=%r; %s' % (ifsave, line)
    #print(line)
    #print(f'ncomponents={ncomponents} ifsave={ifsave} title={title!r}')
    i += 1
    #return i, ncomponents, ifsave

    params = {}
    for icomponent in range(ncomponents):
        hinge_line = lines[i].rstrip()
        #print('hinge_line = %r' % hinge_line)
        #print "icount=",icount
        # get deflection angle, check for a hinge
        #n+=1
        #print "n = ",n
        #print "trailer[",n+1,"] = ",trailer[n+1].strip()
        hinge_method = int(hinge_line[7:8])
        unused_dangle = hinge_line[12:18]
        unused_hinge_patch = int(hinge_line[1:3])
        #print(f'hinge_method={hinge_method} dangle={dangle} hinge_patch={hinge_patch}')
        i += 1

        line2 = lines[i].rstrip()
        #print(line2)
        impact = int(line2[0:2].strip())
        shadow = int(line2[2:4].strip())
        iprint = int(line2[4:5].strip())
        ipin = int(line2[5:6].strip())
        isave = int(line2[6:7].strip())
        pdata1 = float(line2[10:20].strip())
        pdata2 = float(line2[20:30].strip())
        pdata3 = float(line2[30:40].strip())
        pdata4 = float(line2[40:50].strip())
        pdata5 = float(line2[50:60].strip())
        pdata6 = float(line2[60:70].strip())
        #print(f'icomponent={icomponent} impact={impact} shadow={shadow} iprint={iprint}')
        params[icomponent] = [
            impact, shadow, iprint,
            ipin, isave, pdata1, pdata2, pdata3, pdata4, pdata5, pdata6
        ]
        i += 1

        if hinge_method == 0:
            pass
        elif hinge_method == 3:
            hinge_line = lines[i].rstrip()
            #print('hinge_line2 = %r' % hinge_line)
            x1 = float(hinge_line[0:10])
            y1 = float(hinge_line[10:20])
            z1 = float(hinge_line[20:30])
            x2 = float(hinge_line[30:40])
            y2 = float(hinge_line[40:50])
            z2 = float(hinge_line[50:60])
            unused_hinge = (x1, y1, z1, x2, y2, z2)
            i += 1
            #print('hinge =', hinge)
        else:
            print(hinge_line)
            raise NotImplementedError(hinge_method)
    return i, ncomponents, ifsave, params

def parse_viscous(lines, line, i):
    """
    Parses the viscous table (method=4)

    140   viscous_run (level 2)
      11  1
    000
     5 5000       0.0000    1.0000    0.0000    1.0000    3.0000    3.0000
     2 0 4 1 1 1 1 1 0 0 0 1 1 0 1
       1.500000E2   3.000000E2   1.600000E1   3.000000E1   9.000000E1
       5.000000E2  8.000000E-1   1.200000E1   0.000000E0   0.000000E0
      21  1
      ...
    """
    line = lines[i].rstrip()
    # Skin Friction Base & Title card; A-73
    ncomp = int(line[0:2])
    unused_ifsave = line[2:3]
    unused_title = line[6:66].strip()
    #print(f'ncomp={ncomp} ifsave={ifsave} title={title}')
    i += 1
    #print('***', line)

    for unused_icompi in range(ncomp):
        #geometry data source card; A-73
        line = lines[i].rstrip()
        #print(line)
        unused_a, unused_b = line.split()
        unused_icomp = int(line[:3])
        # 0=skin friction method cards will be read for each alpha/beta
        # 1=read only one skin friction method card and assume it applies to all alpha/betas for this component
        # 2=use the same skin friction data as for the prepvious component (skip reading data)
        unused_isk = line[3]

        unused_nskin_friction_elements = line[4:7]
        #print(f'icomp={icomp} isk={isk} nskin_friction_elements={nskin_friction_elements}')
        i += 6
        # skin friction method cards; A-74
        # 0=level 2 viscous
        # 1=level 1 viscous
        #isfmeth =

        # 0=dont print
        # 1=print detailed skin friction intermediate results
        #iprint =

        # 0=don't save
        # 1=write skin friction coefficient for each element
        #isave =
    line = lines[i].rstrip()
    #print(line)
    return i

def parse_special_routines(lines, line, i, read_special_routines=False):
    """
    Parses the special routines (method=5)

    10000000000000000000
    01100   1 3 0 0 0 0
    01100   2 4 0 0 0 0
    11100   3 1 2 5 6 7

    """
    #if not read_special_routines:
        #raise RuntimeError('special routines')
    line = lines[i].rstrip()
    #print('***', line)
    special_lines = []
    while line[0:1] in ['0', '1']:
        special_lines.append(line)
        i += 1
        line = lines[i].rstrip()
    #print('***', line)
    #print('**2', lines[i])

    return i
