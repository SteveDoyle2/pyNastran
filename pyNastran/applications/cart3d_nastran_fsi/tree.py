from numpy import array,argsort
from numpy.linalg import norm
from mathFunctions import distance


def loadList(listA):
    """"Turns a dictionary into a list"""
    dictA = {}
    nItems = len(listA)
    for i in range(1,nItems+1):
        dictA[i] = listA.pop(0)
    return dictA


#nodes is a dictionary
class badTree(object):
    """
    treeType = 'node','element'
    links the node number in the fromNodes to:
      - list of NClose toNodes that are closest to each fromNode
      - list of corresponding distances
    """
    def __init__(self,treeType,nClose=2):
        self.tree = {}
        self.nClose = nClose
        self.treeType = treeType
        if treeType!='node' and treeType!='element':  # verifies you're calling the right
            msg = 'Error!  Invalid treeType\n'
            msg += "treeType=|%s| valid='node','element'" %(treeType)
            raise Exception(msg)

    def reduceToClose(self,dFT,sortList):
        #print "len(dFT)=%s len(sortList)=%s" %(len(dFT),len(sortList))
        nClose = min(len(sortList),self.nClose)
        #print "dFT      = ",dFT
        dFTshort = [dFT[sortList[n]] for n in range(nClose)]
        #print "dFTshort = ",dFTshort
        return dFTshort

    def getCloseElementIDs(self,eid):
        assert self.treeType=='element'
        #print "self.tree = ",self.tree
        closeElements = self.tree[eid]
        assert len(closeElements)>0
        return closeElements

    def getCloseNodeIDs(self,nid):
        assert self.treeType=='node'
        closeNodes = self.tree[nid]
        assert len(closeNodes)>0
        return closeNodes

    def buildTree(self,fromNodes,toNodes):
        #print "fromNodes = ",fromNodes
        #print "toNodes = ",toNodes
        fromKeys = fromNodes.keys()
        nFromNodes = len(fromKeys)

        toKeys = toNodes.keys()
        nToNodes = len(toKeys)
        
        nMax = len(fromKeys)
        n = 0
        
        for fromKey in fromKeys:
            if n%1000==0:
                print "n/%s = %s; %.2f%%" %(n,nMax,n*100./nMax)
            #print "fromKey = ",fromKey
            fromNode = fromNodes[fromKey] #.xyz
            #print "fromNode = ",fromNode
            dFT = []  # distances from-to
            for toKey in toKeys:
                toNode = toNodes[toKey] #.xyz
                dist = distance(fromNode,toNode)
                dFT.append(dist)
                #print "dFT = ",dFT
            dFT = array(dFT)

            #print "dFT = ",dFT
            sortList = argsort(dFT)
            #print "sortList = ",sortList

            dFTshort = self.reduceToClose(dFT,sortList)
            nIDshort = self.reduceToClose(toKeys,sortList)
            
            #print "dFTshort = ",dFTshort
            #print "nIDshort = ",nIDshort
            #print "node[%s]=%s" %(nIDshort[0],toNodes[nIDshort[0]])
            #print "node[%s]=%s" %(nIDshort[1],toNodes[nIDshort[1]])
            self.tree[fromKey] = [nIDshort,dFTshort]
            n+=1
        print "finished constructing distance tree"
        return self.tree


if __name__=='__main__':
    n1 = array([0.,0.,0.])
    n2 = array([1.,1.,1.])
    n3 = array([1.,0.,0.])
    n4 = array([5.,3.,0.])
    n5 = array([2.,0.,4.])
    
    nodes = {
        1:n1,
        2:n2,
        3:n3,
        4:n4,
        5:n5,}
    treeObj = Tree(treeType='node',nClose=3)
    tree = treeObj.buildTree(nodes,nodes)
    
    for nkey,dist in tree.items():
        print nkey,dist[0],dist[1]
    

