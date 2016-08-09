from __future__ import print_function
from six import iteritems
import os
import sys

#from pyNastran.converters.panair.panairGridPatch import PanairGridHelper
from pyNastran.converters.cart3d.cart3d import Cart3d

def load_panair_file(fname='panair.in'):
    if not os.path.exists(fname):
        raise IOError('%s does not exist' % fname)
    execfile(fname)
    varnames = {
                   'title': 'default title',
                   'alpha': 0.,
                   'alpha_compressibility': 0.,
                   'beta': 0.,
                   'beta_compressibility': 0.,
                   'xy_sym': True,
                   'yz_sym': False,
                   'mach':  0.8,
                   'Sref':  1.,
                   'Bref':  1.,
                   'Cref':  1.,
                   'Dref':  1.,
                   'Xref':  0.,
                   'Yref':  0.,
                   'Zref':  0.,
                   'bcMap': {},
                }
    varmap = {}
    localvars = locals()
    for varname, default in sorted(iteritems(varnames)):
        if varname in localvars:
            print("type(%s)=%s type(localvars[varname])=%s" % (varname, type(default), type(localvars[varname])))
            if type(default) != type(localvars[varname]):
                #msg = 'type(%s) != type(%s)\n' %(default, localvars[varname])
                msg = 'type(%s)=%s type(%s)=%s' %(default, type(default),
                                                   localvars[varname], type(localvars[varname]))
                raise RuntimeError(msg)
            assert type(default)==type(localvars[varname])
            varmap[varname] = localvars[varname]
        else:
            if default is not {}:
                varmap[varname] = default
            else:
                raise RuntimeError('variable %s is not defined' % varname)

    if 'bcMap' not in varmap:
        raise RuntimeError('variable bcMap is not defined')
    return varmap

if 0:
    title = 'simple wing-body with composite panel. (run with a502i)'
    alphas = 4.
    alpha_compressibility = 4.

    beta = 0.
    beta_compressibility = 0.
    xy_sym = True
    yz_sym = False
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
    def __init__(self, cart3d_geom_filename, oname, varmap):
        self.printout = ("$printout options\n"
                         "=isings   igeomp    isingp    icontp    ibconp    iedgep\n"
                         "4.        0.        0.        1.        1.        0.\n"
                         "=ipraic   nexdgn    ioutpr    ifmcpr\n"
                         ".0        .0        1.        0.                  3.\n")
        self.run(cart3d_geom_filename, oname, varmap)

    def write_points(self, point1, point2):
        point1 = self.fix_point(point1)
        point2 = self.fix_point(point2)
        out = "%-10s" * 6 % (point1[0], point1[1], point1[2],
                             point2[0], point2[1], point2[2])
        return out + '\n'

    def write_point(self, point1):
        point1 = self.fix_point(point1)
        out = "%-10s" * 3 % (point1[0], point1[1], point1[2])
        return out + '\n'

    def fix_point(self, pointIn):
        pointOut = []
        for value in pointIn:
            sValue = '%s' % (value)
            if len(sValue) > 10:
                sValue = sValue[0:9]
            pointOut.append(sValue.rstrip('0'))
            #print("sValue=%s len=%s" %(sValue,len(sValue)))
        #print("pointOut = ",pointOut)
        return pointOut

    def run(self, cart3d_geom_filename, oname, varmap):
            f = open(oname, 'wb')
            print("oname", oname)
            self.title = varmap['title']
            self.mach = varmap['mach']
            self.ncases = 1
            self.alphaC = varmap['alpha_compressibility']
            self.alphas = [varmap['alpha']]

            self.betaC = varmap['beta_compressibility']
            self.betas = [varmap['beta']]
            self.xref = varmap['Xref']
            self.yref = varmap['Yref']
            self.zref = varmap['Zref']
            self.sref = varmap['Sref']
            self.bref = varmap['Bref']
            self.cref = varmap['Cref']
            self.dref = varmap['Dref']
            self.isEnd = True
            msg = ''
            msg += self.write_title()
            msg += self.write_mach()
            msg += self.write_cases()
            msg += self.write_alphas()
            msg += self.write_betas()
            msg += self.write_reference_quantities()
            msg += self.printout
            f.write(msg)

            BCMap = varmap['bcMap']

            cart = Cart3d()
            cart.read_cart3d(cart3d_geom_filename)
            points = cart.points
            elements = cart.elements
            regions = cart.regions

            #for pid, point in sorted(iteritems(points)):
                #if pid<85:
                #    print(pid,point)
                #pass
            region_old = 9
            for eid, element in sorted(iteritems(elements)):
                header = ''
                region = regions[eid]
                if region not in BCMap:
                    #continue
                    msg = 'regionID=%s is not defined in the BCMap' % region
                    raise RuntimeError(msg)

                (kt, cpNorm) = BCMap[region]
                if cpNorm is None:
                    cpNorm = ''
                if region != region_old:
                    header += '=region %s\n' % region
                region_old = region

                #print("****")
                #print('element = %s' % element)
                #print('region  = %s' % region)
                #if eid==2:
                    #print('points = %s' % points)
                nid1, nid2, nid3 = element
                n1, n2, n3 = points[nid1], points[nid2], points[nid3]
                #print("n1=%s" % (n1))
                #print("n2=%s" % (n2))
                #print("n3=%s" % (n3))
                #p1 =
                #sys.exit()

                net_name = 'e%s' % eid

                header += '$points - surface panels\n'

                header += '%-10s%-10s\n' % ('1.', cpNorm)  # nNetworks is 1
                header += '%-10s\n' % sInt(kt)
                header += '%-10s%-10s%50s%-10s\n' % (
                    sInt(2), sInt(2), '', net_name)
                points_out = self.write_points(n1, n2)
                points_out += self.write_points(n3, n3)
                f.write(header + points_out)
                #break
            #print(points)
            #print(outfilename)

            f.write('$end of panair inputs\n')
            #sys.exit()


def main():
    panair_in = sys.argv[1]
    cart3d_geom_filename = sys.argv[2]
    panair_inp_filename = sys.argv[3]
    varmap = load_panair_file(panair_in)
    #cart3dGeom = os.path.join('models', 'threePlugs.tri')
    #cart3dGeom  = os.path.join('models','spike.a.tri')
    #outfilename = os.path.join('models', 'panair.inp')
    Cart3dToPanair(cart3d_geom_filename, panair_inp_filename, varmap)
    print("done...")

if __name__ == '__main__':  # pragma: no cover
    main()
