## GNU Lesser General Public License
## 
## Program pyNastran - a python interface to NASTRAN files
## Copyright (C) 2011-2012  Steven Doyle, Al Danial
## 
## Authors and copyright holders of pyNastran
## Steven Doyle <mesheb82@gmail.com>
## Al Danial    <al.danial@gmail.com>
## 
## This file is part of pyNastran.
## 
## pyNastran is free software: you can redistribute it and/or modify
## it under the terms of the GNU Lesser General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
## 
## pyNastran is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU Lesser General Public License
## along with pyNastran.  If not, see <http://www.gnu.org/licenses/>.
## 
# pylint: disable=R0904,R0902

from __future__ import division, print_function
from numpy import array, cross

from pyNastran.bdf.cards.loads.staticLoads import Moment, Force
from pyNastran.bdf.cards.elements.shell import ShellElement


class BDFMethods(object):
    def __init__(self):
        pass
    #--------------------
    # METHODS

    def MassProperties(self):
        """
        Caclulates mass properties in the global system about <0,0,0>
        I   = mass*centroid*centroid
        Ixx =  mass*(y^2+z^2)
        Iyz = -mass*y*z
        http://en.wikipedia.org/wiki/Moment_of_inertia#Moment_of_inertia_tensor
        """
                 #Ixx Iyy Izz, Ixy, Ixz Iyz
        I  = array([0., 0., 0.,  0.,  0., 0.,])
        cg = array([0., 0., 0.])
        mass = 0.
        for element in self.elements.itervalues():
            try:
                p = element.Centroid()  # not really coded across the board
                m = element.Mass()
                (x, y, z) = p
                x2 = x*x
                y2 = y*y
                z2 = z*z
                I[0] += m*(y2+z2)  # Ixx
                I[1] += m*(x2+z2)  # Iyy
                I[2] += m*(x2+y2)  # Izz
                I[3] -= m*x*y      # Ixy
                I[4] -= m*x*z      # Ixz
                I[5] -= m*y*z      # Iyz
                mass += m
                cg += m*p
            except:
                self.log().warning("could not get inertia for element"
                                   "...\n%s" %(element))
            ###
        ###
        cg = cg/mass
        return (mass, cg, I)

    def Mass(self):
        """Caclulates mass in the global coordinate system"""
        mass = 0.
        for element in self.elements.itervalues():
            m = element.Mass()
            mass += m
        return (mass)

    def resolveGrids(self, cid=0):
        """
        Puts all nodes in a common coordinate system (mainly for cid testing)
        @param self the object pointer
        @param cid the cid to resolve the nodes to
        @note loses association with previous coordinate systems so to go back
        requires another fem
        """
        assert cid in self.coords, ('cannot resolve nodes to '
                                    'cid=|%s| b/c it doesnt exist' %(cid))
        for nid,node in sorted(self.nodes.iteritems()):
            p = node.PositionWRT(self, cid)
            #p = node.Position(self)
            #print "p = ",p
            node.UpdatePosition(self, p, cid)
        ###

    def unresolveGrids(self, femOld):
        """
        Puts all nodes back to original coordinate system.
        @param self
          the object pointer
        @param femOld
          the old model that hasnt lost it's connection to the node cids
        @warning
          hasnt been tested well...
        """
        debug = False
        for (nid, nodeOld) in femOld.nodes.iteritems():
            coord = femOld.node.cp
            (p, matrix) = coord.transformToGlobal(self.xyz, debug=debug)
            p2 = coord.transformToLocal(p,matrix, debug=debug)
            self.nodes[nid].UpdatePosition(self, p2, coord.cid)
        ###

    def sumForces(self):
        """
        Sums applied forces for all load cases.
        Considers FORCE, FORCE1, FORCE2.

        @retval Forces the forces as a numpy array
        @warning not validated
        """
        for (key, loadCase) in self.loads.iteritems():
            F = array([0., 0., 0.])
            #print "loadCase = ",loadCase
            for load in loadCase:
                #print "load = ",load
                if isinstance(load, Force):
                    f = load.mag*load.xyz
                    #print "f = ",f
                    F += f
                ###
            self.log.info("case=%s F=%s\n\n" %(key, F))
        ###
        return F

    def sumMoments(self, p0):
        """
        Sums applied forces & moments about a reference point p0 
        for all load cases.
        Considers FORCE, FORCE1, FORCE2, MOMENT, MOMENT1, MOMENT2.

        @param p0 the reference point
        @retval Moments the moments as a numpy array
        @retval Forces the forces as a numpy array
        @warning not validated
        """
        p = array(p0)
        for (key, loadCase) in self.loads.iteritems():
            M = array([0., 0., 0.])
            F = array([0., 0., 0.])
            #print "loadCase = ",loadCase
            for load in loadCase:
                #print "load = ",load
                if isinstance(load, Force):
                    f = load.mag*load.xyz
                    node = self.Node(load.node)
                    #print "node = ",node
                    r = node.Position() - p
                    m = cross(r,f)
                    #print "m    = ",m
                    M += m
                    F += f
                elif isinstance(load, Moment):
                    m = load.mag*load.xyz
                    M += m
                ###
            print("case=%s F=%s M=%s\n\n" %(key, F, M))
        ###
        return (M, F)

