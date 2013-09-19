#imort sys
#from math import sqrt

class Node(object):
    pass

class KdTree(object):
    def __init__(self,treeType,pointList,nClose=0):
        self.nClose = nClose
        self.treeType = treeType
        if treeType!='node' and treeType!='element':  # verifies you're calling the right
            msg = 'Error!  Invalid treeType\n'
            msg += "treeType=|%s| valid='node','element'" %(treeType)
            raise Exception(msg)
        
        nodes = []
        for nid,nodeLoc in sorted(pointList.items()):
            #print "nid     = ",nid
            #print "nodeLoc = ",nodeLoc
            #print "pointList[15] = ",pointList[15]
            n = list(nodeLoc)+[nid]
            #print "n = ",n
            nodes.append(n)
        self.tree = self.buildTree(nodes,nClose)

    def getCloseElementIDs(self,point):
        #print "point = ",point
        closeNodesDists = self.nNearestPoints(point,self.nClose)
        #print "closeNodesDists = ",closeNodesDists
        
        closeIDs = []
        dists    = []
        for nodeDist in closeNodesDists:
            ID = nodeDist[0][3]
            dist = nodeDist[1]
            #print "ID = ",ID
            closeIDs.append(ID)  # nid
            dists.append(dist)
        return closeIDs,dists

    def buildTree(self,pointList, depth=0):
        if not pointList:
            return None

        # Select axis based on depth so that axis cycles through all valid values
        k = len(pointList[0])-1 # Assumes all points have the same dimension, and that
        # the last element in a point is an identifier
        axis = depth % k

        # Sort point list to select median
        pointList.sort(key=lambda x: x[axis])
        median = len(pointList) // 2 # Choose median

        # Create node and construct subtrees
        node = Node()
        node.location   = pointList[median]
        node.leftChild  = self.buildTree(pointList[:median], depth+1)
        node.rightChild = self.buildTree(pointList[median+1:], depth+1)
        return node

    def InCircle(self,p,q,radius):
        d = 0.
        for i in range(len(p)):
            d = d+(p[i]-q[i])**2
        if d<=radius**2:
            return True
        return False

    def Distance(self,p,q):
        d = 0.
        for i in range(len(p)):
            d += (p[i]-q[i])**2
        return d**0.5

    def PointsInSphere(self,tree,p,radius,depth=0,ptlist=[]):
        if depth==0:
            #print "PointsInSphere"
            #print tree
            #print p
            #print radius
            #print depth
            #print ptlist
            ptlist = []
        
        k = len(p)
        axis = depth%k
        distance = self.Distance(p,tree.location)
        if distance<=radius:
            ptlist.append((tree.location,distance))
        if tree.leftChild:
            if tree.location[axis]>=p[axis]-radius:
                self.PointsInSphere(tree.leftChild,p,radius,depth+1,ptlist)
        if tree.rightChild:
            if tree.location[axis]<=p[axis]+radius:
                self.PointsInSphere(tree.rightChild,p,radius,depth+1,ptlist)
        return ptlist


    def nNearestPoints(self,p,n):
        tree = self.tree
        ptlist = None
        self.radius_guess = 0.
        npts = 0
        while npts<n:
            if self.radius_guess>0.:
                self.radius_guess = self.radius_guess*2.
            else:
                self.radius_guess = 0.001
            #print "radius_guess = ",self.radius_guess
            ptlist = self.PointsInSphere(tree,p,self.radius_guess)
            npts = len(ptlist)
            #print "npts = ",npts
        ptlist.sort(key=lambda x: x[1])
        self.radius_guess = ptlist[n-1][1]
        return ptlist[:n]

