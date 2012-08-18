class PanairWrite(object):
    def __init__(self):
        pass

    def printAbutments(self):
        msg = ''
        msg += '               SUMMARY OF FACING SURFACES (+:upper, -:lower)\n'
        msg += ' abutment   nw-ident  ntd  knet.edge    nw-ident  ntd  knet.edge\n'
        msg += '        1   winga      12     1.3+      wingwk     18     7.1+\n'
        msg += '            wingwk     18     7.1-      winga      12     1.1+\n'
        msg += '            winga      12     1.1-      winga      12     1.3-\n'
        msg += '        9   1st p-o-s   0    -1.0       bodyl      12     3.4+\n'
        msg += '            bodyl      12     3.4-      1st p-o-s   0    -1.0\n'

        for patchID, patch in self.patches.items():
            (p1, x1, y1, z1) = patch.getEdges()
            self.log.debug("p[%s] = %s" % (patchID, p1))
            #print "x = ",x1
            #print "y = ",y1
            #print "z = ",z1

        return msg

    def printOptions(self):
        msg = ''
        msg += '0               print options\n'
        msg += '            %i = singularity grid print flag\n' % (
            self.isings)
        msg += '            %i = panel geometry print flag\n' % (self.igeomp)
        msg += '            %i = spline data flag  ( 0 ==> off, nonzero ==> on )\n' % (self.isingp)
        msg += '            %i = control point information print flag\n' % (
            self.icontp)
        msg += '            %i = boundary condition data print flag \n' % (
            self.ibconp)
        msg += '            %i = edge matching information print flag\n' % (
            self.iedgep)
        msg += '            %i = index of control point for which aic-s are printed\n' % (self.ipraic)
        msg += '            %i = edge control point flow properties print flag\n' % (self.nexdgn)
        msg += '            %i = output control flag (-1 ==> no surface flow properties, 0 ==> standard output, 1 ==> short form output )\n' % (self.ioutpr)
        msg += '            %s = force/moment control flag (-1 ==> no force and moment data, 0 ==> standard output, 1 ==> nw totals only )\n' % (self.ifmcpr)
        msg += '            %i = print flag for detailed cost information during execution of job\n' % (self.icostp)
        msg += '            1 = print flag for singularity parameter maps\n'
        msg += '0               abutment processing options\n'
        msg += '   %10s = global edge abutment tolerance specified by user.  if this value is zero, a default value will be calculated\n' % (fortranValue(self.epsgeo))
        msg += '                later.   this default value is taken as:  .001  * (minimum panel diameter)\n'
        msg += '            %i = print flag controlling geometry printout  b e f o r e  the abutment processing.  ( nonzero ==> do print )\n' % (self.igeoin)
        msg += '            %i = print flag controlling geometry printout  a f t e r    the abutment processing.  ( nonzero ==> do print )\n' % (self.igeout)
        msg += '            %i = network/abutment/abutment-intersection print flag.  ( nonzero ==> generate the cross referenced abutment listing\n' % (self.nwxref)
        msg += '            %i = control index for panel intersection checking.  ( nonzero ==> do perform the check. )\n' % (self.triint)
        msg += '            %i = abutment/abutment-intersection (short listing) print flag ( 0 ==> suppress, nonzero ==> generate usual print )\n' % (2)
        msg += ' \n'
        msg += '                force and moment reference parameters\n'
        msg += '   %10s = reference area for force and moment calculations.    (sref)\n' % (fortranValue(self.sref))
        msg += '   %10s = rolling moment reference length  (bref)\n' % (
            fortranValue(self.bref))
        msg += '   %10s = pitching moment reference length (cref)\n' % (
            fortranValue(self.cref))
        msg += '   %10s = yawing moment reference length   (dref)\n' % (
            fortranValue(self.dref))
        msg += '   %10s = x - coordinate for the point about which moments will be calculated  (xref)\n' % (fortranValue(self.xref))
        msg += '   %10s = y - coordinate for the point about which moments will be calculated  (yref)\n' % (fortranValue(self.yref))
        msg += '   %10s = z - coordinate for the point about which moments will be calculated  (zref)\n' % (fortranValue(self.zref))
        msg += '            3 = pressure coefficient index (nprcof) (1=linear, 2=slenderbody, 3=2nd, 4=isentropic)\n'
        msg += '1\n'
        return msg

    def printOutHeader(self):
        msg = """\n
         ****************************************************************************************************

          dynamic memory management initialization

          max no. levels         15   max no. arrays        200   maximum scratch storage     900000   total storage provided  900000
             addr(maplev)          0     addr(maplws)          0     addr(scratch storage)          1

        ****************************************************************************************************
          wopen call on unit    1  blocks:    10   status:     0
          wopen call on unit    2  blocks:    10   status:     0
          wopen call on unit    3  blocks:    10   status:     0
        1


        *****************************************************************************
         *                                                                           *
         *                     a502 - pan-air technology program                     *
         *                                                                           *
         *               potential flow about arbitrary configurations               *
         *               version id = ht2 (12 feb 92) boeing ver i00                 *
         *                                                                           *
         *                                   02/12/92                                *
         *                                                                           *
         *                                                                           *
         * simple wing-body with composite panel. (run with a502i)
         * saaris  865-6209  m/s 7c-36
         *                                                                           *
         *                                                                           *
         *                                                                           *
         *****************************************************************************
        1
        0*b*input-da






                               - list of a502 input data cards -\n"""
        return msg

    def printGridSummary(self):
        msg = ''
        msg2 = ''
        msg3 = ''

        self.nOffBodyPoints = 0
        self.nStreamlines = 0
        msg += '        1               ***  quick summary of a502 input  ***\n'

        for i, titleLine in enumerate(self.titleLines):
            msg += "title%s:%s\n" % (i + 1, titleLine)

        msg += '0               processing options\n'
        msg += '            %i = datacheck.   (0=regular run,1=full datacheck,2=short datacheck)\n' % (self.dataCheck)
        msg += '            0 = s.p. flag.   (0 ==> no s.p. file (ft09) provided, 1 ==> local file ft09 with singularity values is provided)\n'
        msg += '            0 = aic flag.    (0 ==>  no aic file (ft04) provided, 1 ==> local file ft04 with aic-s is provided by the user)\n'
        msg += '            0 = b.l. flag    (0 ==> no boundary layer file requested, 1 ==> boundary layer data will be written to file ft17)\n'
        msg += '            0 = velocity correction index.  (0 ==> no correction, 1 ==> mclean correction, 2 ==>  boctor correction)\n'
        msg += '            0 = flow visualization flag.  (nonzero ==> off-body and streamline processing will be performed)\n'
        msg += '            0 = off-body calculation type. (0 ==> mass flux, nonzero ==> velocity)\n'
        msg += '            0 = streamline calculation type. (0 ==> mass flux, nonzero ==> velocity)\n'
        msg += '           %2i = number of off-body points.\n' % (
            self.nOffBodyPoints)
        msg += '           %2i = number of streamlines to be traced.\n' % (
            self.nStreamlines)
        msg += '0               case summary\n'
        msg += '           %2i = number of cases\n' % (self.ncases)
        msg += '     %f = mach number\n' % (self.mach)
        msg += '     %f = compressibility axis angle of attack (alpc)\n' % (
            self.alphaC)
        msg += '     %f = compressibliity axis angle of sideslip (betc)\n' % (
            self.betaC)

        msg3 += '0network id&index   #rows   #cols  kt  src  dblt  nlopt1  nropt1  nlopt2  nropt2    ipot   # pts  # pans  cpnorm  cum pt  cum pn\n'
        msg3 += '---------- -----   -----   -----  --  ---  ----  ------  ------  ------  ------    ----    ----    ----  ------  ------  ------\n'

        totalPoints = 0
        totalPanels = 0
        for patchID in xrange(self.nPatches()):
            patch = self.patch(patchID)
            msg3 += patch.quickSummary(totalPoints, totalPanels)
            totalPanels += patch.nPanels()
            totalPoints += patch.nPoints()

        msg2 += '0  case       alpha          beta      mag(f-s-v)\n'
        msg2 += ' ------    ----------    ----------   -----------\n'
        for (iCase) in xrange(self.ncases):
            alpha = self.alphas[iCase]
            beta = self.betas[iCase]
            msg2 += '     %2s      %f      %f      1.000000\n' % (
                iCase + 1, alpha, beta)

        msg2 += '0               symmetry options\n'
        msg2 += '            %s = number of planes of symmetry\n' % (
            self.nSymmetryPlanes)
        msg2 += '            %s = x-z plane of symmetry flag (0 ==> no symmetry, 1==> flow symmetry, -1 ==> flow antisymmetry)\n' % (self.XZsymmetry)
        msg2 += '            %s = x-y plane of symmetry flag (0 ==> no symmetry, 1==> flow symmetry, -1 ==> flow antisymmetry)\n' % (self.XYsymmetry)
        msg2 += '0               configuration summary\n'
        msg2 += '          %3s = total number of networks read in\n' % (
            self.nNetworks)
        msg2 += '         %4s = total number of mesh points\n' % (totalPoints)
        msg2 += '         %4s = total number of panels\n' % (totalPanels)

        return msg + msg2 + msg3

    def writeDataCheck(self):
        msg = '$datacheck\n'
        msg += '%s.\n' % (self.dataCheck)
        return msg

    def writePrintout(self):
        msg = '$printout options\n'
        msg += "%-10s%-10s%-10s%-10s%-10s%-10s\n" % (self.isings, self.igeomp,
                                                     self.isingp, self.icontp, self.ibconp, self.iedgep)
        msg += "%-10s%-10s%-10s%-10s%-10s\n" % (self.ipraic,
                                                self.nexdgn, self.ioutpr, self.ifmcpr, self.icostp)
        return msg
