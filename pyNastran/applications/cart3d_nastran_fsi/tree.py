from __future__ import print_function
from six.moves import range
from numpy import array, argsort
#from numpy.linalg import norm
from pyNastran.applications.cart3d_nastran_fsi.math_functions import distance


#def load_list(listA):
    #""""Turns a dictionary into a list"""
    #dictA = {}
    #nItems = len(listA)
    #for i in range(1, nItems+1):
        #dictA[i] = listA.pop(0)
    #return dictA


#nodes is a dictionary
class BadTree(object):
    """
    tree_type = 'node','element'
    links the node number in the fromNodes to:
      - list of NClose toNodes that are closest to each fromNode
      - list of corresponding distances
    """
    def __init__(self, tree_type, nClose=2):
        self.tree = {}
        self.nClose = nClose
        self.tree_type = tree_type
        if tree_type not in ['node', 'element']:
            # verifies you're calling the right tree
            msg = 'Error!  Invalid tree_type\n'
            msg += "tree_type=%r valid='node','element'" % tree_type
            raise Exception(msg)

    def reduce_to_close(self, dFT, sort_list):
        #print "len(dFT)=%s len(sort_list)=%s" % (len(dFT), len(sort_list))
        nClose = min(len(sort_list), self.nClose)
        #print "dFT      = ",dFT
        dFTshort = [dFT[sort_list[n]] for n in range(nClose)]
        #print "dFTshort = ", dFTshort
        return dFTshort

    def getCloseElementIDs(self, eid):
        assert self.tree_type == 'element'
        #print "self.tree = ", self.tree
        close_elements = self.tree[eid]
        assert len(close_elements) > 0
        return close_elements

    def getCloseNodeIDs(self, nid):
        assert self.tree_type == 'node'
        close_nodes = self.tree[nid]
        assert len(close_nodes) > 0
        return close_nodes

    def build_tree(self, from_nodes, to_nodes):
        #print "from_nodes = ", from_nodes
        #print "to_nodes = ", to_nodes
        from_keys = from_nodes.keys()
        #n_from_nodes = len(from_keys)

        to_keys = to_nodes.keys()
        #n_to_nodes = len(to_keys)

        nMax = len(from_keys)
        n = 0

        for from_key in from_keys:
            if n % 1000 == 0:
                print("n/%s = %s; %.2f%%" % (n, nMax, n * 100. / nMax))
            #print "from_key = ",from_key
            from_node = from_nodes[from_key] #.xyz
            #print "from_node = ",from_node
            dFT = []  # distances from-to
            for to_key in to_keys:
                to_node = to_nodes[to_key] #.xyz
                dist = distance(from_node, to_node)
                dFT.append(dist)
                #print "dFT = ", dFT
            dFT = array(dFT)

            #print "dFT = ",dFT
            sort_list = argsort(dFT)
            #print "sort_list = ",sort_list

            dFTshort = self.reduce_to_close(dFT, sort_list)
            nIDshort = self.reduce_to_close(to_keys, sort_list)

            #print "dFTshort = ", dFTshort
            #print "nIDshort = ", nIDshort
            #print "node[%s]=%s" %(nIDshort[0] ,toNodes[nIDshort[0]])
            #print "node[%s]=%s" %(nIDshort[1], toNodes[nIDshort[1]])
            self.tree[from_key] = [nIDshort, dFTshort]
            n += 1
        print("finished constructing distance tree")
        return self.tree


def main():
    n1 = array([0., 0., 0.])
    n2 = array([1., 1., 1.])
    n3 = array([1., 0., 0.])
    n4 = array([5., 3., 0.])
    n5 = array([2., 0., 4.])

    nodes = {
        1 : n1,
        2 : n2,
        3 : n3,
        4 : n4,
        5 : n5,
    }
    tree_obj = BadTree(tree_type='node', nClose=3)
    tree = tree_obj.build_tree(nodes, nodes)

    for nkey, dist in tree.items():
        print(nkey, dist[0], dist[1])


if __name__ == '__main__':
    main()
