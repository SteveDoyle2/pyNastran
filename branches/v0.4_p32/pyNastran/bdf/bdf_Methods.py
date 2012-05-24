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
from numpy import array
from pyNastran.bdf.cards.loads import *
from pyNastran.bdf.errors import *

# 3rd party
#import numpy
#from numpy import any,cross


class bdfMethods(object):
    def __init__(self):
        pass
    #--------------------
    # METHODS

    def MassProperties(self):
        """
        Caclulates mass properties in the global system about 0,0,0
        I   = mass*centroid*centroid
        Ixx =  mass*(y^2+z^2)
        Iyz = -mass*y*z
        http://en.wikipedia.org/wiki/Moment_of_inertia#Moment_of_inertia_tensor
        """
                 #Ixx Iyy Izz, Ixy, Ixz Iyz
        I  = array([0., 0., 0.,  0.,  0., 0.,])
        cg = array([0., 0., 0.])
        mass = 0.
        for eid,element in self.elements.items():
            try:
                p = e.Centroid()  # not really coded across the board
                m = e.Mass()
                (x,y,z) = p
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
                self.log().warning("could not get inertia for element...\n%s" %(element))
            ###
        ###
        cg = cg/mass
        return (mass,cg,I)

    def Mass(self):
        """Caclulates mass in the global coordinate system"""
        mass = 0.
        for element in self.elements:
            m = e.Mass()
            mass += m
        return (mass)

    def resolveGrids(self,cid=0):
        """
        puts all nodes in a common coordinate system (mainly for cid testing)
        @param self the object pointer
        @param cid the cid to resolve the node to
        @note loses association with previous coordinate systems so to go back
        requires another fem
        """
        assert cid in self.coords,'cannot resolve nodes to cid=|%s| b/c it doesnt exist' %(cid)
        for nid,node in sorted(self.nodes.items()):
            p = node.PositionWRT(self,cid)
            #p = node.Position(self)
            #print "p = ",p
            node.UpdatePosition(self,p,cid)
        ###

    def unresolveGrids(self,femOld):
        """
        puts all nodes back to original coordinate system
        @param self the object pointer
        @param femOld the old model that hasnt lost it's connection to the node cids
        @warning hasnt been tested well...
        """
        for nid,nodeOld in femOld.nodes.items():
            coord = femOld.node.cp
            p,matrix  = coord.transformToGlobal(self.xyz,debug=debug)
            p2 = coord.transformToLocal(p,matrix,debug=debug)
            node.UpdatePosition(self,p2,cid)
        ###

    def sumForces(self):
        for key,loadCase in self.loads.items():
            F = array([0.,0.,0.])
            #print "loadCase = ",loadCase
            for load in loadCase:
                #print "load = ",load
                if isinstance(load,Force):
                    f = load.mag*load.xyz
                    #print "f = ",f
                    F += f
                ###
            self.log.info("case=%s F=%s\n\n" %(key,F))
        ###

    def sumMoments(self,p0):
        p = array(p0)
        for key,loadCase in self.loads.items():
            M = array([0.,0.,0.])
            F = array([0.,0.,0.])
            #print "loadCase = ",loadCase
            for load in loadCase:
                #print "load = ",load
                if isinstance(load,Force):
                    f = load.mag*load.xyz
                    node = self.Node(load.node)
                    #print "node = ",node
                    r = node.Position() - p
                    m = cross(r,f)
                    #print "m    = ",m
                    M += m
                    F += f
                elif isinstance(load,Moment):
                    m = load.mag*load.xyz
                    M += m
                ###
            print("case=%s F=%s M=%s\n\n" %(key,F,M))
        ###
