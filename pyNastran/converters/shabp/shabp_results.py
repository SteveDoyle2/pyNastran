from numpy import zeros

class ShabpOut:
    def __init__(self, model, log=None, debug=False):
        self.model = model
        self.X = model.X
        self.npatches = len(model.X)

    #@property
    #def component_name_to_num(self):
        #return self.model.component_name_to_num

    def readline(self, f, i):
        i += 1
        return f.readline(), i

    def readline_n(self, f, i, n):
        i += n
        for unused_j in range(n-1):
            f.readline()
        return f.readline(), i

    def read_shabp_out(self, out_filename):
        #npatches = len(self.X)
        istart = zeros(self.npatches, dtype='int32')
        nelements = 0
        for ipatch in range(self.npatches):
            X = self.X[ipatch]
            nrows, ncols = X.shape
            nelementsi = (nrows-1) * (ncols-1)
            istart[ipatch] = nelements
            nelements += nelementsi
        #print('istart =', istart)

        Cp_out = {}
        delta_out = {}

        # all of case 1, then all of case 2
        Cp_dict, delta_dict, ncases = self._parse_results(out_filename)
        #unused_ncomps = len(self.component_num_to_name)
        #print('Cp_dict =', Cp_dict)
        components = self.model.component_name_to_patch.keys()
        for icase in range(ncases):
            #print(f'ncases={ncases} components={components} nelements={nelements}')
            Cp = zeros(nelements, dtype='float32')
            delta = zeros(nelements, dtype='float32')
            for name in sorted(components):
                icomp = self.model.component_name_to_num[name]
                patches = self.model.component_name_to_patch[name]
                #print(f'  name={name} icomp={icomp} patches={patches}')

                Cp_array = Cp_dict[icomp]
                #print(f'  Cp_array={Cp_array}')
                delta_array = delta_dict[icomp]
                ndata = len(Cp_array) // ncases
                jelement_start = ndata * icase
                for unused_i, ipatch in enumerate(patches):  # ipatch starts at 1
                    X = self.X[ipatch-1]
                    iistart = istart[ipatch-1]
                    #print(f'    X.shape={X.shape} iistart={iistart}')

                    nrows, ncols = X.shape
                    nelementsi = (nrows-1) * (ncols-1)
                    Cp[iistart:iistart+nelementsi] = Cp_array[jelement_start:jelement_start+nelementsi]
                    delta[iistart:iistart+nelementsi] = delta_array[jelement_start:jelement_start+nelementsi]
                    jelement_start += nelementsi
            #print('Cp =', Cp)
            #print('delta =', delta)
            Cp_out[icase] = Cp
            delta_out[icase] = delta
            #break
        return Cp_out, delta_out

    def _parse_results(self, out_filename):
        with open(out_filename, 'r') as resfile:
            i = 0
            line, i = self.readline(resfile, i)
            while '******** MAIN PROGRAM NOW HAS CONTROL OF SYSTEM ********' not in line:
                line, i = self.readline(resfile, i)

            while '*** PRESSURE CALCULATION PROGRAM' not in line:
                line, i = self.readline(resfile, i)

            out = self._read_inviscid_pressure(resfile, i)
        line, i, Cp_dict_components, delta_dict_components, ncases = out
        #self._read_viscous2(line, i)
        return Cp_dict_components, delta_dict_components, ncases

    def read_viscous2(self, f, i):
        line = ''
        while '*** SKIN FRICTION FORCE PROGRAM' not in line:
            line, i = self.readline(f, i)
        line, i = self.readline_n(f, i, 2)
        assert 'Alpha' in line, line

        nstreamlines = 19
        streamlines = []
        for unused_istreamline in range(nstreamlines):
            # streamline 1
            streamline = []
            line, i = self.readline(f, i)

            #                Skin          Heat Adiabatic
            #    Wetted  Friction      Transfer      Wall      Wall                  Wall  Flow
            #  Distance    Coeff.        Coeff.  Enthalpy     Temp.   Heat Flux  Pressure  Regime           X         Y         Z  RE_TH/ML     MACHL
            #      (in)   (---)     (lbm/ft2/s)  (Btu/lb)       (F) (Btu/ft2-s) (lbf/ft2)  (---)           (in)      (in)      (in)
            line, i = self.readline_n(f, i, 4)
            while 'Streamline:' not in line and line.strip() != '':
                line, i = self.readline(f, i)
                streamline.append(line.strip())
            streamlines.append(streamline)
        return streamlines

    def _read_inviscid_pressure(self, f, i):
        #npatches = 44
        #6
        #0ELEMENT DATA   MACH=  6.000  ALT =  50000.  S REF =196272.0  SPAN =  669.6  IMPACT =  1  IMPACI =  3
        #                XCG = -713.4  YCG =    0.0    ZCG  =    0.0    MAC =  240.0   ISHAD =  1  ISHADI =  3
        #      ANGLE OF ATTACK = 45.00   YAW ANGLE =  0.00   K = 1.83000   ETAC =  1.0000   DELTA E =  0.00
        Cp_dict_components = {}
        delta_dict_components = {}
        _line, i = self.readline_n(f, i, 5)
        for icomponent in range(self.npatches):
            mach_line, i = self.readline(f, i)
            if 'Summation Number' in mach_line:
                break
            if 'MACH' not in mach_line:
                raise RuntimeError('couldnt find MACH in line[%i]=%r' % (i, mach_line.strip()))

            unused_xcg_line, i = self.readline(f, i)
            unused_alpha_line, i = self.readline(f, i)

            unused_xcent_line, i = self.readline_n(f, i, 3)
            #0  L      DEL CA        DEL CY        DEL CN       DEL CLL       DEL CLM       DEL CLN       CP             AREA
            #            CA            CY            CN           CLL           CLM            CLN        DELTA
            #           XCENT         YCENT         ZCENT      NX             NY             NZ

            line, i = self.readline(f, i)
            Cp = []
            Delta = []
            while 1:
                while line[0] == '0':
                    (del_ca, del_cy, del_cn, del_cll, del_clm, del_cln, cp, area) = (
                        line[8 :20], line[21:33],
                        line[35:47], line[49:63],
                        line[63:75], line[77:89],
                        line[91:103], line[105:117],
                    )
                    del del_ca, del_cy, del_cn, del_cll, del_clm, del_cln, area

                    line, i = self.readline(f, i)
                    unused_ca, unused_cy, unused_cn, unused_cll, unused_clm, unused_cln, delta = (
                        line[8 :20], line[21:33],
                        line[35:47], line[49:63],
                        line[63:75], line[77:89],
                        line[91:103],
                    )

                    line, i = self.readline(f, i)
                    unused_xc, unused_yc, unused_zc, unused_nx, unused_ny, unused_nz = (
                        line[8 :20], line[21:33],
                        line[35:47], line[49:63],
                        line[63:75], line[77:89],
                    )

                    Cp.append(cp)
                    Delta.append(delta)
                    line, i = self.readline(f, i)

                #4 check on COMPONENT
                line, i = self.readline_n(f, i, 4)
                if '*AIR' in line:
                    # keep going
                    line, i = self.readline_n(f, i, 7)
                elif 'COMPONENT' in line:
                    #npanels = None
                    while 1:
                        #PANELS:  patch32             patch33             patch41
                        line, i = self.readline(f, i)
                        panels = line.strip().split()
                        if len(panels) == 0:
                            break

                    while 'Mach   Velocity   Reynolds #   Altitude   -Freestream Conditions-   Gas' not in line:
                        line, i = self.readline(f, i)

                    while 'ALPHA  BETA     C N        C A        C M        C L        C D        L/D        C Y        C LN       C LL' not in line:
                        line, i = self.readline(f, i)

                    # this is a weird way to read this...
                    ncases = 0
                    while 'S/HABP' not in line:
                        #ncases += 1
                        line, i = self.readline(f, i)
                        #print(line.rstrip())
                    ncases += 1
                    #ncases -= 2  # correct for reading too many lines

                    line, i = self.readline_n(f, i, 3)
                    break
            Cp_dict_components[icomponent] = Cp
            delta_dict_components[icomponent] = Delta
        return line, i, Cp_dict_components, delta_dict_components, ncases
