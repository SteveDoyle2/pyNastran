from numpy import array
from pyNastran.bdf.cards.loads import *
from pyNastran.bdf.errors import *

from pyNastran.bdf.cards.plates.elementsShell import ShellElement
from pyNastran.general.generalMath import buildSpline
from scipy.interpolate import splev

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
        for eid,element in self.elements.iteritems():
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

    def flipNormals(self,starterEid,eids=None,flipStarter=False):
        """
        Takes the normals of SHELL elements and flips it to a common direction
        This method follows the contour of the body, so assuming
        no internal elements, all the normals on the outside will point
        outwards (or inwards).

        @param starterEid the element to copy the normal of
        @param eids       the element IDs to flip to the common direction (default=None -> all)
        @param flipStarter should the staring element be flipped (default=False)
        
        @todo finish method...think i need to build a edge list...
              that'd be a lot easier to loop through stuff...
        """
        normals   = {}
        validNids = set([])
        isCorrectNormal = set([])

        allEids = eids
        if allEids==None:
            allEids==self.elements.keys()
        ###
        setAllEids = set(allEids)


        if flipStarter:
            elem = self.Element(starterEid)
            elem.flipNormal()
        normals[starterEid] = elem.Normal()
        
        for eid in allEids:
            if isinstance(element,ShellElement):
                elem = self.Element(starterEid)
                normals[starterEid] = elem.Normal()
                validNids = validNids.union(set(elem.nodeIDs()))
            ###
        ###
        
        ## clean up the elements that will be considered
        elemsToCheck = set([])
        nidToEidMap = self.getNodeIDToElementIDsMap()
        for nid,eidsMap in sorted(nidToEidMap.iteritems()):
            if nid not in validNids:  # clean up extra nodes
                del nidToEidMap[nid]
            else:
                eids = list(set(eids))  # do i need this?
                for eid in eids:  # clean up ROD/SOLID elements
                    eids2 = []
                    if eid in setAllEids:
                        eids2.append(eid)
                    ###
                    elemsToCheck = elemsToCheck.union(set(eids2))
                nidToEidMap[nid] = eids2
            ###
        ###
        
        ## starts with the starter element, loops thru adjacent elements
        ## and checks to see if the normal is 'close' to the elements
        ## normal from before
        goEid = starterEid
        
        # no recursion to avoid recursion limit
        while 1:
            elem = self.Element(goEid)
            nids = elem.getNodeIDs()
            normals = self.getAdjacentNormals(nids,nidToEidMap)
            normal = normals[goEid]
        ###

    def getAdjacentElements(self,nids,nidToEidMap):
        """
        @todo doesnt work...
        """
        normals = {}
        #for nid in 

    def resolveGrids(self,cid=0):
        """
        puts all nodes in a common coordinate system (mainly for cid testing)
        @param self the object pointer
        @param cid the cid to resolve the node to
        @note loses association with previous coordinate systems so to go back
        requires another fem
        """
        assert cid in self.coords,'cannot resolve nodes to cid=|%s| b/c it doesnt exist' %(cid)
        for nid,node in sorted(self.nodes.iteritems()):
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
        for nid,nodeOld in femOld.nodes.iteritems():
            coord = femOld.node.cp
            p,matrix  = coord.transformToGlobal(self.xyz,debug=debug)
            p2 = coord.transformToLocal(p,matrix,debug=debug)
            node.UpdatePosition(self,p2,cid)
        ###

    def sumForces(self):
        for key,loadCase in self.loads.iteritems():
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
        for key,loadCase in self.loads.iteritems():
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

    def tankFill(self,mFuel,percentStart=.50,rhoFuel=51.19,
                      tankElements=[],gravity=None,
                      mTol=0.05,nIterMax=10,addElements=True):
        """
        Fills a single fuel tank in consistent units
        @param percentStart the percentage to start at. percentStart = z/(zMax-zMin); default=0.50
        @param mFuel mass (or weight) of fuel to add to the tank
        @param rhoFuel density of fuel (default = 51.19 lb/ft^3 = 6.6 lb/gal)
        @param tankElements list of elements defining the boundary of the tank
        @param gravity vector defining the direction of gravity in cid=0; <0,0,-32.2>=default
        @param massTol tolerance on the mass (massError = massTol*mFuel; default=0.05)
        @param nIterMax the maximum number of iterations (default=10)
        @param addElements create CONM2 elements if True

        @retval nodalMasses the masses on each node
        @retval percentStarts the percentages used in the interpolation
        @retval massToli = (mFuel-totalMass)/mFuel

        @note massTol should be valid for any tank size
              (adjust the percent error from 5% if necessary)

        Method:
           1.  Create a new CORD2R coordinate system (cid=-2)
           2.  Rotate the geometry into the gravity coordinate system
           3.  find the z0 (zero fill line) by taking z0=percentStart*(zMax-zMin)
           4.  find the nodes below the waterline (z0)
           5.  apply mass to the nodes based on their depth; mass = (z-z0)*g*rho
           6.  compare this to the required mass (or weight) for convergence test
           7.  interpolate using a spline to get the next point to check
        
        Requirements:
           1.  No cid=-2 is already being used
           2.  Tank doesnt have to be closed (but it's probably a good thing)
           3.  No requirement on normals
           4.  percentStart < 1.
           5.  mTol < 0.20
        """
        assert percentStart<1.
        assert mTol<0.20
        percentStarts.append(0.)  # x; empty tank
        massErrors.append(-mFuel) # y; the mass is too low by mFuel

        # find the vector with the maximum difference with the gravity vector
        if gravity==None:
            gravity = array([0.,0.,-32.2])
        magGravity = abs(array) # magnitude
        A=array([1.,0.,0.]) # global x
        B=array([0.,1.,0.]) # global y
        C=array([0.,0.,1.]) # global z
        
        AdotGravity = dot(A,gravity)
        BdotGravity = dot(B,gravity)
        CdotGravity = dot(C,gravity)
        ABC = [A,B,C]
        ABCdot = [AdotGravity,CdotGravity,CdotGravity]
        ABCmax = max(ABCdot)
        i = ABCdot.index(ABCmax)

        (bx,by,bz) = ABC[i]
        gx,gy,gz = gravity

        # Create a new CORD2R coordinate system (cid=-2)
        # origin is at [0.,0.,0.] and doesnt matter
        # z axis is the gravity direction (gravity)
        # point on x-z plane is the max unit vector (ABC)
        coordCard = ['CORD2R',-2,0,  0.,0.,0., gx,gy,gz, bx,by,bz]
        cardName = 'CORD2R'
        self.addCard(self,coordCard,cardName)
        #cardObj = BDF_Card(coordCard)
        #coord = CORD2R(cardObj)
        
        # convert all points into the gravity frame
        cid = -2
        elementNodeIDs = {} ## CQUAD4, CTRIA3
        nodeLocations = {}
        
        for eid in tankElements:
            elem = self.elements[eid]
            if (elem.type=='CQUAD4') or (elem.type=='CTRIA3'):
                nodes = elem.nodes

                elementNodeIDs[eid] = []
                for node in nodes:
                    nid = node.nid
                    elementNodeIDs[eid].append(nid)
                    if nid not in nodeLocations:
                        p = node.PositionWRT(self,cid)
                        nodeLocations[nid] = p
                    ###
                ###
            ###
        ###
        
        zMax = nodeLocations[nid][2]
        zMin = zMax
        for nid,node in sorted(nodeLocations.iteritems()):
            zMax = max(zMax,node[2])
            zMin = min(zMin,node[2])
        
        # max sure to go into the while looop
        massToli = 2. # 2*mFuel
        percentFill = percentStart

        nIter = 0
        while massToli>massTol and nIter<nIterMax:
            # find the z0 (zero fill line) by taking z0=percentStart*(zMax-zMin)
            z0 = percentFill*(zMax-zMin)

            aboveNodes=set(); belowNodes=set()
            for nid,node in sorted(nodeLocations.iteritems()):
                if node[2]>=z0:
                    aboveNodes.add(nid)
                else:
                    belowNodes.add(nid)
                ###
            ###

            if 0:
                belowElements=set(); partialElements=set()
                for eid,nodeIDs in sorted(elementNodeIDs.items()):
                    elem = self.elements[eid]

                    isAboveBelow = set() # True=Above False=Below
                    for nid in nodeIDs:
                        if nid in aboveNodes:
                            isAboveBelow.add(True)
                        else:
                            isAboveBelow.add(False)
                        ###
                    ###
                    if   True in isAboveBelow and False not in isAboveBelow:  # all nodes are above
                        pass
                    elif True not in isAboveBelow and False in isAboveBelow:  # all nodes are below
                        belowElements.add(eid)
                    elif True in isAboveBelow and False in isAboveBelow:  # some nodes are above, some below
                        partialElements.add(eid)
                    else;
                        raise RuntimeError('not above, not below, not partial...')
                    ###
                ###

            if 0:
                for eid in belowElements:
                    elem = self.elements[eid]
                    nodeIDs = elementNodeIDs[eid]
                    if elem.type=='CQUAD4':
                        pass
                    elif elem.type=='CTRIA3':
                        pass
                    else:
                        raise NotImplementedError()

            # compute the total and elemental masses
            nodalMasses = {}
            totalMass = 0.
            for nid in belowNodes:
                # mass = g*rho*(z0-z)
                # it's (z0-z) b/c magGravity is always positive and z0 is higher than z
                mass = magGravity*(z0-nodeLocations[nid][2])
                nodalMasses[nid] = mass
                totalMass += mass
            ###
            massError = mFuel-totalMass

            percentStarts.append(percentStart) # x
            massErrors.append(massError)       # y
            massToli = massError/mFuel
            
            #x=[]; y=[]
            for xi,yi in massFound:
                x.append(xi)
                y.append(yi)
            i = argsort(x)
            X = array(percentStarts)[i]  # sorted x
            Y = array(massToli)[i]       # sorted y
            
            spline = buildSpline(Y,X) # reverse interpolation
            yi = 0. # find 0. mass
            xi = splev(yi,spline) # the percentFill for 0. mass
            percentFill = xi
            
            nIter+=1
        ###

        if addElements:
            maxEid = max(self.elements)+1 # get the next available eid
            for nid,mass in sorted(nodalMasses.iteritems()):
                card = ['CONM2',maxEid,nid,0,mass]
                self.addCard(self,card,'CONM2')
                maxEid+=1

        del self.coords[cid]
        return masses,X,Y
