import os
#import sys

from pyNastran.converters.panair.panairGridPatch import PanairGridHelper
from pyNastran.converters.cart3d.cart3d_reader import genericCart3DReader

title = 'simple wing-body with composite panel. (run with a502i)'
alphas = 4.
alphaCompressibility = 4.

beta = 0.
betaCompressibility = 0.
xySym = True
yzSym = False
mach = 0.6
Sref = 2400.
Bref = 60.
Cref = 40.
Dref = 90.
xref = 46.
yref = 0.
zref = 0.
bcMap = {
    1: [1., None],  # kt,cpnorm
    #2: [1.,2.], # kt,cpnorm
    #3: [1.,2.], # kt,cpnorm
}

#$title
#simple wing-body with composite panel. (run with a502i)
#saaris  865-6209  m/s 7c-36
#$datacheck
# 0.
#$symmetry - xz plane of symmetry
#=misymm   mjsymm
#1.        0.
#$mach number
#=amach
#.6
#$cases - no. of solutions
#=nacase
#1.
#$angles-of-attack
#=alpc
#4.
#=alpha(1) alpha(2)  alpha(3)
#4.        10.       0.
#$printout options
#=isings   igeomp    isingp    icontp    ibconp    iedgep
#4.        0.        0.        1.        1.        0.
#=ipraic   nexdgn    ioutpr    ifmcpr
#.0        .0        1.        0.                  3.
#$references for accumulated forces and moments
#=xref     yref      zref      nref
#46.       0.        0.
#=sref     bref      cref      dref
#2400.     60.       40.       90.


def sInt(value):
    """
    int represented as a short float
    """
    value = "%f" % (value)
    return value.rstrip('0')


class Cart3dToPanair(PanairGridHelper):
    def __init__(self, cart3dGeom, oname, BCMap):
        self.printout = ("$printout options\n"
                         "=isings   igeomp    isingp    icontp    ibconp    iedgep\n"
                         "4.        0.        0.        1.        1.        0.\n"
                         "=ipraic   nexdgn    ioutpr    ifmcpr\n"
                         ".0        .0        1.        0.                  3.\n")
        self.run(cart3dGeom, oname, BCMap)

    def writePoints(self, point1, point2):
        point1 = self.fixPoint(point1)
        point2 = self.fixPoint(point2)
        #print point1
        #print point2
        out = "%-10s" * 6 % (point1[0], point1[1], point1[2],
                             point2[0], point2[1], point2[2])
        return out + '\n'

    def writePoint(self, point1):
        point1 = self.fixPoint(point1)
        out = "%-10s" * 3 % (point1[0], point1[1], point1[2])
        return out + '\n'

    def fixPoint(self, pointIn):
        pointOut = []
        for value in pointIn:
            sValue = '%s' % (value)
            if len(sValue) > 10:
                sValue = sValue[0:9]
            pointOut.append(sValue.rstrip('0'))
            #print "sValue=%s len=%s" %(sValue,len(sValue))
        #print "pointOut = ",pointOut
        return pointOut

    def run(self, cart3dGeom, oname, BCMap):
            f = open(oname, 'wb')
            print "oname", oname
            self.mach = mach
            self.ncases = 1
            self.alphaC = alphaCompressibility
            self.alphas = [alphas]

            self.betaC = betaCompressibility
            self.betas = [beta]
            self.xref = xref
            self.yref = yref
            self.zref = zref
            self.sref = Sref
            self.bref = Bref
            self.cref = Cref
            self.dref = Dref
            self.isEnd = True
            msg = ''
            #msg += self.writeTitle()
            msg += self.writeMach()
            msg += self.writeCases()
            msg += self.writeAlphas()
            msg += self.writeBetas()
            msg += self.writeReferenceQuantities()
            msg += self.printout
            f.write(msg)

            cart = genericCart3DReader(cart3dGeom)
            (points, elements, regions, loads) = cart.readCart3d(cart3dGeom)

            for pid, point in sorted(points.iteritems()):
                #if pid<85:
                #    print pid,point
                pass
            for eid, element in sorted(elements.iteritems()):
                region = regions[eid]
                if region not in BCMap:
                    continue
                (kt, cpNorm) = BCMap[region]
                if cpNorm is None:
                    cpNorm = ''

                #print "****"
                #print "element =",element
                #print "region  =",region
                #if eid==2:
                #    print "points = ",points
                nid1, nid2, nid3 = element
                n1, n2, n3 = points[nid1], points[nid2], points[nid3]
                #print "n1=%s" %(n1)
                #print "n2=%s" %(n2)
                #print "n3=%s" %(n3)
                #p1 =
                #sys.exit()

                netName = 'e%s' % (eid)

                header = '$points - surface panels\n'

                header += '%-10s%-10s\n' % ('1.', cpNorm)  # nNetworks is 1
                header += '%-10s\n' % (sInt(kt))
                header += '%-10s%-10s%50s%-10s\n' % (
                    sInt(2), sInt(2), '', netName)
                pointsOut = self.writePoints(n1, n2)
                pointsOut += self.writePoints(n3, n3)
                f.write(header + pointsOut)
                #break
            #print points
            #print outfilename

            f.write('$end of panair inputs\n')
            #sys.exit()


if __name__ == '__main__':
    cart3dGeom = os.path.join('models', 'threePlugs.tri')
    #cart3dGeom  = os.path.join('models','spike.a.tri')
    outfilename = os.path.join('models', 'panair.inp')
    Cart3dToPanair(cart3dGeom, outfilename, bcMap)
    print "done..."
