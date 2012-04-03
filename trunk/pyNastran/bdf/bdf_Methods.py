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
        I = mass*centroid*centroid
        Ixx = mass*x*x
        Iyz = mass*y*z
        @warning centroid isnt coded across the board
        """
                 #Ixx Iyy Izz, Ixy, Ixz Iyz
        I = array([0., 0., 0.,  0.,  0., 0.,])
        cg = array([0., 0., 0.])
        mass = 0.
        for element in self.elements:
            p = e.Centroid()  # not really coded across the board
            m = e.Mass()
            (x,y,z) = p
            I[0] = m*x*x  # Ixx
            I[1] = m*y*y  # Iyy
            I[2] = m*z*z  # Izz
            I[3] = m*x*y  # Ixy
            I[4] = m*x*z  # Ixz
            I[5] = m*y*z  # Iyz
            cg += m*p
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
        @warning hasnt been tested...
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
                    print "f = ",f
                    F += f
                ###
            self.log.info("case=%s F=%s\n\n" %(key,F))
        ###

    def sumMoments(self):
        p = array([0.,0.5,0.])
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
            print "case=%s F=%s M=%s\n\n" %(key,F,M)
        ###
