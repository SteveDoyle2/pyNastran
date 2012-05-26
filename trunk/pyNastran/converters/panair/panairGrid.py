import os
import sys
import copy
from numpy import array,zeros,ones,transpose,cross,dot
from numpy.linalg import norm

from panairGridPatch import PanairPatch,PanairWakePatch,PanairGridHelper,sInt
from pyNastran.general.general import ListPrint

#CL = -Fx*sin(alpha)*cos(beta) + Fy*sin(alpha)*sin(beta) +Fz*cos(alpha)
#CD =  Fx*cos(alpha)*cos(beta) - Fy*cos(alpha)*sin(beta) +Fz*sin(alpha)
#CY =  Fx*sin(beta) +Fy*cos(beta)


def fortranValue(value):
    return "%8.4E" %(value)

class PanairGrid(PanairGridHelper):
    def __init__(self,infileName):
        self.infileName = infileName
        self.nNetworks = 0
        self.patches = {}
        
        self.alphas = [0.]
        self.ncases = 0.
        self.betas = [0.]
        self.alphaC = 0.
        self.betaC  = 0.


        self.sref=1.
        self.bref=1.
        self.cref=1.
        self.dref=1.

        self.xref=0.
        self.yref=0.
        self.zref=0.

        self.xyzSection = ''
        self.streamlineSection = ''
        self.flowSection = ''
        self.sectionalPropSection = ''
        self.gridSection = ''

        self.msg = ''
    
    def printAbutments(self):
        msg = ''
        msg += '               SUMMARY OF FACING SURFACES (+:upper, -:lower)\n'
        msg += ' abutment   nw-ident  ntd  knet.edge    nw-ident  ntd  knet.edge\n'
        msg += '        1   winga      12     1.3+      wingwk     18     7.1+\n'
        msg += '            wingwk     18     7.1-      winga      12     1.1+\n'
        msg += '            winga      12     1.1-      winga      12     1.3-\n'
        msg += '        9   1st p-o-s   0    -1.0       bodyl      12     3.4+\n'
        msg += '            bodyl      12     3.4-      1st p-o-s   0    -1.0\n'

        for patchID,patch in self.patches.items():
            (p1,x1,y1,z1) = patch.getEdges()
            print "p[%s] = %s" %(patchID,p1)
            #print "x = ",x1
            #print "y = ",y1
            #print "z = ",z1

        return msg

    def printOptions(self):
        msg=''
        msg += '0               print options\n'
        msg += '            %i = singularity grid print flag\n'  %(self.isings)
        msg += '            %i = panel geometry print flag\n' %(self.igeomp)
        msg += '            %i = spline data flag  ( 0 ==> off, nonzero ==> on )\n' %(self.isingp)
        msg += '            %i = control point information print flag\n' %(self.icontp)
        msg += '            %i = boundary condition data print flag \n' %(self.ibconp)
        msg += '            %i = edge matching information print flag\n' %(self.iedgep)
        msg += '            %i = index of control point for which aic-s are printed\n' %(self.ipraic)
        msg += '            %i = edge control point flow properties print flag\n' %(self.nexdgn)
        msg += '            %i = output control flag (-1 ==> no surface flow properties, 0 ==> standard output, 1 ==> short form output )\n' %(self.ioutpr)
        msg += '            %s = force/moment control flag (-1 ==> no force and moment data, 0 ==> standard output, 1 ==> nw totals only )\n' %(self.ifmcpr)
        msg += '            %i = print flag for detailed cost information during execution of job\n' %(self.icostp)
        msg += '            1 = print flag for singularity parameter maps\n'
        msg += '0               abutment processing options\n'
        msg += '   %10s = global edge abutment tolerance specified by user.  if this value is zero, a default value will be calculated\n' %(fortranValue(self.epsgeo))
        msg += '                later.   this default value is taken as:  .001  * (minimum panel diameter)\n'
        msg += '            %i = print flag controlling geometry printout  b e f o r e  the abutment processing.  ( nonzero ==> do print )\n' %(self.igeoin)
        msg += '            %i = print flag controlling geometry printout  a f t e r    the abutment processing.  ( nonzero ==> do print )\n' %(self.igeout)
        msg += '            %i = network/abutment/abutment-intersection print flag.  ( nonzero ==> generate the cross referenced abutment listing\n' %(self.nwxref)
        msg += '            %i = control index for panel intersection checking.  ( nonzero ==> do perform the check. )\n' %(self.triint)
        msg += '            %i = abutment/abutment-intersection (short listing) print flag ( 0 ==> suppress, nonzero ==> generate usual print )\n' %(2)
        msg += ' \n'
        msg += '                force and moment reference parameters\n'
        msg += '   %10s = reference area for force and moment calculations.    (sref)\n' %(fortranValue(self.sref))
        msg += '   %10s = rolling moment reference length  (bref)\n' %(fortranValue(self.bref))
        msg += '   %10s = pitching moment reference length (cref)\n' %(fortranValue(self.cref))
        msg += '   %10s = yawing moment reference length   (dref)\n' %(fortranValue(self.dref))
        msg += '   %10s = x - coordinate for the point about which moments will be calculated  (xref)\n' %(fortranValue(self.xref))
        msg += '   %10s = y - coordinate for the point about which moments will be calculated  (yref)\n' %(fortranValue(self.yref))
        msg += '   %10s = z - coordinate for the point about which moments will be calculated  (zref)\n' %(fortranValue(self.zref))
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
    def printFile(self):
        msg = ''
        for i,line in enumerate(self.lines):
            msg += "%6s %s" %(i+1,line )
        msg += '                                      record of input processing\n\n\n'
        return msg

    def printGridSummary(self):
        msg='';msg2='';msg3=''

        self.nOffBodyPoints = 0
        self.nStreamlines = 0
        msg += '        1               ***  quick summary of a502 input  ***\n'
        
        for i,titleLine in enumerate(self.titleLines):
            msg += "title%s:%s\n" %(i+1,titleLine)
        
        msg += '0               processing options\n'
        msg += '            %i = datacheck.   (0=regular run,1=full datacheck,2=short datacheck)\n' %(self.dataCheck)
        msg += '            0 = s.p. flag.   (0 ==> no s.p. file (ft09) provided, 1 ==> local file ft09 with singularity values is provided)\n'
        msg += '            0 = aic flag.    (0 ==>  no aic file (ft04) provided, 1 ==> local file ft04 with aic-s is provided by the user)\n'
        msg += '            0 = b.l. flag    (0 ==> no boundary layer file requested, 1 ==> boundary layer data will be written to file ft17)\n'
        msg += '            0 = velocity correction index.  (0 ==> no correction, 1 ==> mclean correction, 2 ==>  boctor correction)\n'
        msg += '            0 = flow visualization flag.  (nonzero ==> off-body and streamline processing will be performed)\n'
        msg += '            0 = off-body calculation type. (0 ==> mass flux, nonzero ==> velocity)\n'
        msg += '            0 = streamline calculation type. (0 ==> mass flux, nonzero ==> velocity)\n'
        msg += '           %2i = number of off-body points.\n' %(self.nOffBodyPoints)
        msg += '           %2i = number of streamlines to be traced.\n' %(self.nStreamlines)
        msg += '0               case summary\n'
        msg += '           %2i = number of cases\n' %(self.ncases)
        msg += '     %f = mach number\n' %(self.mach)
        msg += '     %f = compressibility axis angle of attack (alpc)\n' %(self.alphaC)
        msg += '     %f = compressibliity axis angle of sideslip (betc)\n' %(self.betaC)


        msg3 += '0network id&index   #rows   #cols  kt  src  dblt  nlopt1  nropt1  nlopt2  nropt2    ipot   # pts  # pans  cpnorm  cum pt  cum pn\n'
        msg3 += '---------- -----   -----   -----  --  ---  ----  ------  ------  ------  ------    ----    ----    ----  ------  ------  ------\n'

        totalPoints = 0
        totalPanels = 0
        for patchID in range(self.nPatches()):
            patch = self.patch(patchID)
            msg3 += patch.quickSummary(totalPoints,totalPanels)
            totalPanels += patch.nPanels()
            totalPoints += patch.nPoints()
        ###

        msg2 += '0  case       alpha          beta      mag(f-s-v)\n'
        msg2 += ' ------    ----------    ----------   -----------\n'
        for (iCase) in range(self.ncases):
            alpha=self.alphas[iCase]; beta=self.betas[iCase]
            msg2 += '     %2s      %f      %f      1.000000\n' %(iCase+1,alpha,beta)

        msg2 += '0               symmetry options\n'
        msg2 += '            %s = number of planes of symmetry\n' %(self.nSymmetryPlanes)
        msg2 += '            %s = x-z plane of symmetry flag (0 ==> no symmetry, 1==> flow symmetry, -1 ==> flow antisymmetry)\n' %(self.XZsymmetry)
        msg2 += '            %s = x-y plane of symmetry flag (0 ==> no symmetry, 1==> flow symmetry, -1 ==> flow antisymmetry)\n' %(self.XYsymmetry)
        msg2 += '0               configuration summary\n'
        msg2 += '          %3s = total number of networks read in\n' %(self.nNetworks)
        msg2 += '         %4s = total number of mesh points\n' %(totalPoints)
        msg2 += '         %4s = total number of panels\n' %(totalPanels)

        return msg+msg2+msg3

    def getPanelPoints(self,iPanel):
        r = iPanel%(self.nRows-1)
        c = iPanel/(self.nRows-1)

        #print "r=%s c=%s" %(r,c)
        p1 = self.getPoint(r,  c  )
        p2 = self.getPoint(r,  c+1)
        p3 = self.getPoint(r+1,c+1)
        p4 = self.getPoint(r+1,c  )
        return (p1,p2,p3,p4)

    def nPanels(self):
        totalNPanels = 0
        for patchID in range(self.nPatches()):
            patch = self.patch(patchID)
            if patch.kt==1:
                totalNPanels += patch.nPanels()
                print "nPanels = ",totalNPanels
        return totalNPanels

    def nPatches(self):
        return len(self.patches)

    def patch(self,ID):
        return self.patches[ID]

    def updateCases(self):
        """
        reduces confusion by only printing cases that will run
        """
        if len(self.alphas)>self.ncases:
           self.alphas = self.alphas[:self.ncases] 
        if len(self.betas)>self.ncases:
           self.betas = self.alphas[:self.ncases] 

    def writeGrid(self,outfileName):
        self.updateCases()
        outfile = open(outfileName,'wb')
        outfile.write(self.titleSection)
        outfile.write(self.writeDataCheck())
        outfile.write(self.symmetrySection)
        
        outfile.write(self.writeMach())
        outfile.write(self.writeCases())
        outfile.write(self.writeAlphas())
        outfile.write(self.writeReferenceQuantities())


        #outfile.write(self.alphaSection)
        #outfile.write(self.caseSection)
        
        for patchName,patch in sorted(self.patches.items()):
            outfile.write(str(patch))
        
        outfile.write(self.xyzSection)
        outfile.write(self.streamlineSection)
        outfile.write(self.flowSection)
        outfile.write(self.sectionalPropSection)
        outfile.write(self.gridSection)
        outfile.write(self.writePrintout())

        outfile.write(self.peaSection)
        outfile.write(self.writeLiberalizedAbutments())
        outfile.write(self.writeEnd())

        outfile.close()

    def removeComments(self,lines):
        lines2 = []
        for line in lines:
            line = line.rstrip().lower()
            if '=' in line:
                #print "line -> |%s|" %(line)
                if '=' is not line[0]:
                    print "line[0] -> ",line[0]
                    line = line.split('=')[0]
                    print "******"
                    lines2.append(line)
                # else: skip
            else:
                lines2.append(line)
        #sys.exit('stopping')
        
        return lines2

    def getTitle(self,section):
        #print "hi"
        #self.title = section[1:]
        self.titleSection = '\n'.join(section)+'\n'
        self.titleLines = section[1:]
        return True

    def writeDataCheck(self):
        msg = '$datacheck\n'
        msg += '%s.\n' %(self.dataCheck)
        return msg

    def getDataCheck(self,section):
        self.dataCheck = int(float(section[1][0:10]))
        self.dataCheck = 2
        return True

    def getSymmetry(self,section):
        """
        $symmetry - xz plane of symmetry
        =misymm   mjsymm
        1.        0.
        
        @warning
            doesnt consider antisymmetryic
        """
        self.XZsymmetry = int(float(section[1][0 :10]))  # doesnt consider antisymmetric
        self.XYsymmetry = int(float(section[1][10:20]))  # doesnt consider antisymmetric
        self.nSymmetryPlanes = self.XZsymmetry+self.XYsymmetry
        self.symmetrySection = '\n'.join(section)+'\n'
        print "XZsymmetry=%s XYsymmetry=%s" %(self.XZsymmetry,self.XYsymmetry)
        return True

    def splitPoints(self,lines,nActual,nRemainder):
        """
        reads the points
        """
        points = []
        for n in range(nActual):
            #(x1,y1,z1) = lines[n][0 :10].strip(),lines[n][10:20].strip(),lines[n][20:30].strip()
            #print "x1=%s y1=%s z1=%s" %(x1,y1,z1)
            (x1,y1,z1) = float(lines[n][0 :10]),float(lines[n][10:20]),float(lines[n][20:30])
            (x2,y2,z2) = float(lines[n][30:40]),float(lines[n][40:50]),float(lines[n][50:60])
            point1 = array([x1,y1,z1])
            point2 = array([x2,y2,z2])
            #point1 = [x1,y1,z1]
            #point2 = [x2,y2,z2]
            #print ListPrint(point1)
            #print ListPrint(point2)
            points.append(point1); points.append(point2)
        ###

        if nRemainder:
            n+=1
            #print "***"
            (x1,y1,z1) = float(lines[n][0 :10]),float(lines[n][10:20]),float(lines[n][20:30])
            point1 = array([x1,y1,z1])
            #print ListPrint(point1)
            points.append(point1)
        #print "points = ",ListPrint(points)
        return points

    def addWakePatch(self,netName,options,x,y,z):
        #print "self.nNetworks = ",self.nNetworks
        patch = PanairWakePatch(self.nNetworks,netName,options,x,y,z)
        self.msg += patch.process()
        #print "patch = ",patch
        self.patches[patch.iNetwork] = patch # deepcopy?
        self.nNetworks+=1

    def addPatch(self,netName,kt,cpNorm,x,y,z):
        #print "self.nNetworks = ",self.nNetworks
        patch = PanairPatch(self.nNetworks,netName,kt,cpNorm,x,y,z)
        self.msg += patch.process()
        #print "patch = ",patch
        self.patches[patch.iNetwork] = patch # deepcopy?
        self.nNetworks+=1
    
    def findPatchByName(self,netName):
        names = []
        for patchID,patch in self.patches.items():
            print "patchID=%s" %(patchID)
            print "*patch = ",patch
            print "patch.netName=%s" %(patch.netName)
            if patch.netName==netName:
                return patch
            names.append(patch.netName)
        raise Exception('couldnt findPatchbyName name=|%s| names=%s' %(netName,names))

    def getPoints(self,section):
        """
        $points - wing-body  with composite panels
        =kn                                               cpnorm
        4.                                                2.
        =kt
        1.
        =nm       nn                                                          netname
        11.       3.                                                          winga
        =x(1,1)   y(1,1)    z(1,1)    x(*,*)    y(*,*)    z(*,*)
           69.4737    9.2105    0.0000   63.7818    9.5807    0.7831
        """
        nNetworks = int(float(section[1][0 :10]))
        cpNorm = section[1][50:60].strip()
        #    cpNorm    = float(section[1][50:60])
        
        kt = int(float(section[2][0 :10]))
        if cpNorm:
            print "nNetworks=%s cpNorm=%s" %(nNetworks,cpNorm)
        else:
            print "nNetworks=%s" %(nNetworks)
        
        n=4
        self.msg += '      kn,kt            %i          %i\n' %(nNetworks,kt)

        for iNetwork in range(nNetworks):
            print "lines[* %s] = %s" %(n-1,section[n-1])
            nm = int(float(section[n-1][0 :10]))
            nn = int(float(section[n-1][10:20]))
            netName = section[n-1][70:80]
            print "kt=%s nm=%s nn=%s netname=%s" %(kt,nm,nn,netName)
            #sys.exit('asdf')
            x = zeros([nm,nn])
            y = zeros([nm,nn])
            z = zeros([nm,nn])
            nFullLines = nm/2
            nPartialLines = nm%2
            nLines = nFullLines+nPartialLines
            #print "nFullLines=%s nPartialLines=%s nLines=%s" %(nFullLines,nPartialLines,nLines)
            for j in range(nn):
                lines = section[n:n+nLines]
                n+=nLines
                #print '\n'.join(lines)
                points = self.splitPoints(lines,nFullLines,nPartialLines)
                
                for i,point in enumerate(points):
                    x[i][j] = point[0]
                    y[i][j] = point[1]
                    z[i][j] = point[2]

            #print "--X--"
            #print x
            self.addPatch(netName,kt,cpNorm,x,y,z)
            n+=1
        return True

    def getGrid(self,section):
        self.gridSection = '\n'.join(section)+'\n'
        return True

    def getXYZ(self,section):
        self.xyzSection = '\n'.join(section)+'\n'
        return True

    def getStreamlines(self,section):
        self.streamlineSection = '\n'.join(section)+'\n'
        return True

    def writePrintout(self):
        msg = '$printout options\n'
        msg += "%-10s%-10s%-10s%-10s%-10s%-10s\n" %(self.isings, self.igeomp, self.isingp, self.icontp, self.ibconp, self.iedgep)
        msg += "%-10s%-10s%-10s%-10s%-10s\n"      %(self.ipraic, self.nexdgn, self.ioutpr, self.ifmcpr, self.icostp)
        return msg

    def getPrintout(self,section):
        """
        isings  igeomp  isingp  icontp  ibconp  iedgep
        ipraic  nexdgn  ioutpr  ifmcpr  icostp
        """
        #self.printoutSection = '\n'.join(section)+'\n'

        ########  ROW 1 #######
        self.isings = float(section[1][ 0:10])
        #self.isings = 5.
        #2.0 output wake singularity grid-singularity strength and gradients at 9 pts/panel
        #3.0 add table of singularity parameter values
        #4.0 add network array of singularity values
        #5.0 add singularity grid for all network
        
        self.igeomp = float(section[1][10:20])
        self.igeomp = 1.
        # 1.0 outputs panel data - geometry diagnostic data
        
        self.isingp = float(section[1][20:30])
        #self.isignp = 1.  # ???
        # 1.0 outputs singularity spline data
        
        self.icontp = float(section[1][30:40])
        self.icontp = 1.
        # 1.0 outputs control pt location, upper surface unit normals, control point 
        #     diagnositc data, control pt maps and singularity parameter maps
        # 2.0 add additional control point map for each network
        
        self.ibconp = float(section[1][40:50])
        #self.ibconp = 0.0
        # 1.0 outputs boundary condition coeffs, diagnositc data, 
        #     boundary condition maps, control point maps, 
        #     and singularity paramter maps

        self.iedgep = float(section[1][50:60])
        self.iedgep = .0
        # 1.0 outputs edge-matching diagnositic data
        #print "igeomp=|%s| isingp=|%s| icontp=|%s| ibconp=|%s| iedgep=|%s|" %(igeomp,isingp,icontp,ibconp,iedgep)
        
        ########  ROW 2 #######

        self.ipraic = float(section[2][ 0:10])
        #self.ipraic = 3.
        # xxx inputs control point sequence number for
        #     which AIC values are to be pritned (rarely used)

        self.nexdgn = float(section[2][10:20])
        self.nexdgn = 0.
        # 1.0 outputs edge and extra control point data

        self.ioutpr = float(section[2][20:30])
        self.ioutpr = -1.0
        # -1.0 omits flow parameter output
        #  0.0 outputs 49 flow parameters
        #  1.0 outputs 13 flow parameters
        
        self.ifmcpr = section[2][30:40].strip()
        self.ifmcpr = 0.
        # 0.0 omits network force and moment output
        # 1.0 outputs network forces and moments summary per column, 
        #     per network, and accumulation over all previous networks
        
        self.icostp = section[2][40:50].strip()
        self.icostp = 0.
        #self.icostp = 1.0
        #print "ipraic=|%s| nexdgn=|%s| ioutpr=|%s| ifmcpr=|%s| ifmcpr=|%s|" %(ipraic,nexdgn,ioutpr,ifmcpr,iedgep)
        # 1.0 outputs detailed job cost for different segments of the
        #     program (rarely used)
        

        #$printout options
        #=isings   igeomp    isingp    icontp    ibconp    iedgep
        #4.        0.        0.        1.        1.        0.
        #=ipraic   nexdgn    ioutpr    ifmcpr
        #.0        .0        1.        0.                  3.

        return True

    def getTrailingWakes(self,section):
        """
        $trailing wakes from body
        =kn                                               cpnorm
        2.
        =kt       matcw
        18.       1.
        =inat     insd      xwake     twake                                   netname
        bodyl     3.        100.      .0                                      bodylwk 
        bodyu     3.        100.      .0                                      bodyuwk 
        """
        nNetworks = int(float(section[1][0 :10]))
        cpNorm = section[1][50:60].strip()
        if cpNorm:
            cpNorm = float(cpNorm)
        else:
            cpNorm = 0
        
        #print "section[2] = ",section[2]
        kt    = int(float(section[2][0 :10]))
        
        matchw = section[2][10:20].strip()
        if matchw:
            matchw = float(matchw)
        #if cpNorm:
        #    print "nNetworks=%s cpNorm=%s" %(nNetworks,cpNorm)
        #else:
        #    print "nNetworks=%s" %(nNetworks)
        #print ""
        #matcw = int(float(section[2][10:20]))
        
        n=3
        self.msg += '      kn,kt            %i          %i\n' %(nNetworks,kt)
        print 'kt=%s cpNorm=%s matchw=%s' %(kt,cpNorm,matchw)
        for iNetwork in range(nNetworks):
            print "lines[* %s] = %s" %(n,section[n])
            
            trailedPanel = section[n][ 0:10].strip()
            edgeNumber   = int(float(section[n][10:20]))
            xWake = float(section[n][20:30]) # x distance
            tWake = float(section[n][30:40]) # 0-wake parallel to x axis;  1-wake in direction of compressibility
            netName = section[n][70:80].strip()
            print 'trailedPanel=%s edgeNumber=%s xWake=%s tWake=%s netName=%s' %(trailedPanel,edgeNumber,xWake,tWake,netName)
            try:
                patch = self.findPatchByName(trailedPanel)
            except KeyError:
                print 'trailedPanel isnt defined...trailedPanel=|%s|' %(trailedPanel)
                raise

            #xPoints = patch.x
            #yPoints = patch.y
            #print "xPoints = ",xPoints
            #print "yPoints = ",yPoints
            (p1,x1,y1,z1) = patch.getEdge(edgeNumber)

            nPoints = len(x1)
            X = zeros([2,nPoints])
            Y = zeros([2,nPoints])
            Z = zeros([2,nPoints])
            
            X[0][:] = x1[0]
            Y[0][:] = y1[0]
            Z[0][:] = z1[0]

            if tWake==0.:
                X[1][:] = ones([nPoints])*xWake
                Y[1][:] = y1[0]
                Z[1][:] = z1[0]
            else:
                #alphaC, betaC
                raise Exception('tWake isnt supported')
            print "--X---"
            print X
            print "--Y---"
            print Y
            print "--Z---"
            print Z
            #nm = int(float(section[n-1][0 :10]))
            #nn = int(float(section[n-1][10:20]))
            options = [kt,cpNorm,matchw,trailedPanel,edgeNumber,xWake,tWake]
            patch = self.addWakePatch(netName,options,X,Y,Z)
            n+=1
        ###
        return True

    def getPartialEdgeAbutments(self,section):
        self.peaSection = '\n'.join(section)+'\n'
        return True

    def writeLiberalizedAbutments(self):
        msg = '$eat - liberalized abutments\n'
        msg += '%10s'*6 %(self.epsgeo,self.igeoin,self.igeout,self.nwxref,self.triint,self.iabsum) + '\n'
        return msg

    def getLiberalizedAbutments(self,section):
        print "section[1] = ",section[1]
        self.epsgeo = section[1][ 0:10].strip()
        if self.epsgeo:
            self.epsgeo = float(self.epsgeo)
        else:
            self.epsgeo = 0.0
        # absolute value of tolerance used to establith equivalent network
        # edge points; minimum panel diagonal times 0.001 <= espgeo <= minimum
        # panel diagonal times 0.03; to not be restricted use negative numbers

        self.igeoin = float(section[1][10:20])
        self.igeoin = -1.
        # prints network geometry before $EAT
        #  0.0 prints geometry
        # -1.0 omits printout

        self.igeout = float(section[1][20:30])
        self.igeout = 0.
        # prints network geometry after $EAT
        # 0.0 omits printout
        # 1.0 prints geometry

        self.nwxref = float(section[1][30:40])
        self.nwxref = 0.
        #self.nwxref = 1.0
        # prints abutment cross-reference; for each network, all abutments and 
        # abutment intersections are described in a format similar to that
        # used in the abutment summary and abutment intersection summary output sections
        # 0.0 omits printout
        # 1.0 prints geometry
        
        self.triint = float(section[1][40:50])
        self.triint = 0.0
        # optionally checks for intersection of network edges through other panels
        # for final geometry; used only if needed on complex configurations; can
        # increase the cost of a computer run by 1 or 2 times that of a datachek (2.) run
        # 0.0 omits intersection checking
        # 1.0 checks intersection
        
        self.iabsum = float(section[1][50:60])
        self.iabsum = -0.0
        #  0.0 abutment summary
        #  1.0 plus intersection summary
        #  2.0 plus complete list of each network and associated abutments
        #  3.0 plus diagnostic data for program developer
        # -1.0 no abutment printout
        return True

    def getSectionalProperties(self,section):
        self.sectionalPropSection = '\n'.join(section)+'\n'
        return True

    def getFlowfieldProperties(self,section):
        self.flowSection = '\n'.join(section)+'\n'
        return True

    #def get(self,section):
    #    pass

    #def get(self,section):
    #    pass

    #def get(self,section):
    #    pass

    def groupSections(self,sections,sectionNames):
        #self.Points = []
        #self.Streamlines = []
        #self.Trailing = []

        #self.Title = None
        #self.dataCheck = None
        #self.symmetry = None
        #self.mach = None
        #self.cases = None
        #self.alphas = None

        #self.pea = None
        #self.eat = None

        #self.printout = None
        #self.end = None

        #self.Grid = []
        
        sectionMap = {
            'tit':self.getTitle,
            'ref':self.getReferenceQuantities,
            'dat':self.getDataCheck,
            'sym':self.getSymmetry,
            'mac':self.getMach,
            'cas':self.getCases,

            'ang':self.getAlphas,
            'yaw':self.getBetas,
            'pri':self.getPrintout,
            'poi':self.getPoints,
            'tra':self.getTrailingWakes,
            'pea':self.getPartialEdgeAbutments,
            'eat':self.getLiberalizedAbutments,
            'sec':self.getSectionalProperties,
            'flo':self.getFlowfieldProperties,
            'xyz':self.getXYZ,
            'gri':self.getGrid,
            'str':self.getStreamlines,
            'end':self.getEnd,
        }
        #validMaps = ['tit','ref','dat','sym','mac','cas','ang','poi','tra','end']
        #validMaps += ['xyz','str','flo','sec','pri','pea','eat','gri']
        #validMaps = ['yaw','poi','tra']
        validMaps = sectionMap.keys()

        for section,sectionName in zip(sections,sectionNames): # 1st line
            self.msg += '  $%s\n' %(sectionName)
            #print "section = ",len(section)
            print "sectionName=%s" %(sectionName)
            if sectionName in validMaps:
                print "section[0] = ",section[0]
                functionMap = sectionMap[sectionName]
                ran = functionMap(section)
                assert ran==True,'%s didnt run' %(sectionName)
                print ""
            if sectionName=='end':
                break
        return sections

    def splitIntoSections(self,lines):
        sections = []
        sectionNames = []
        section = None
        sectionName = None
        for line in lines:
            if '$' in line:
                sections.append(section)
                sectionNames.append(sectionName)
                section = []
                sectionName = line[1:4]
                if sectionName=='end':
                    break
                #print "new section...name=%s line - %s" %(sectionName,line)
            section.append(line)
        ###

        # end section
        section.append(line)
        sections.append(section)
        sectionNames.append(sectionName)
        
        # get rid of first dummy section
        sections = sections[1:]
        sectionNames = sectionNames[1:]
        
        assert len(sections)==len(sectionNames),"%s %s" %(len(sections),len(sectionNames))

        #for section in sections: # 1st line
        #    #print "section = ",len(section)
        #    print "section[0] = ",section[0]
        return (sections,sectionNames)

    def readGrid(self):
        infile = open(self.infileName,'r')
        self.lines = infile.readlines()
        infile.close()
        
        lines = self.removeComments(self.lines)
        (sections,sectionNames) = self.splitIntoSections(lines)
        groups = self.groupSections(sections,sectionNames)
        print "nPatches = ",self.nPatches()
        # split into headings
        #for panel in panels:
        #    points = readPoints()
        print "self.msg = ",self.msg

if __name__=='__main__':
    infileName  = 'SWB.INP'
    outfileName = 'SWB_new.INP'
    #infileName = 'ELLIP.INP'
    grid = PanairGrid(infileName)
    grid.readGrid()
    grid.writeGrid(outfileName)
    print "\ninfileName=%s" %(infileName)

