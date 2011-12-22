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
                 #Ixx Iyy Izz, Ixy, Ixz Iyz
        I = array(0., 0., 0.,  0.,  0., 0.,)
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
        ###
        return I
            
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
