from __future__ import print_function
from six import iteritems
from six.moves import range
from numpy import array
from scipy.spatial import KDTree as scipyKDTree

class Node(object):
    pass

class KdTree2(object):
    def __init__(self, tree_type, point_list, nclose=0):
        self.tree_type = tree_type
        self.nclose = nclose
        if tree_type not in ['node', 'element']:
            # verifies you're calling the right
            msg = 'Error!  Invalid tree_type\n'
            msg += "tree_type=%r valid='node','element'" % tree_type
            raise RuntimeError(msg)
        npoints = len(point_list)
        nodes = array((npoints, 3), dtype='float64')
        self.node_ids = array(npoints, dtype='int32')

        i = 0
        for node_id, xyz in sorted(iteritems(point_list)):
            nodes[i, :] = xyz
            self.node_ids[i] = node_id
        self.tree = scipyKDTree(point_list)

    def get_close_element_ids(self, point):
        dists, i = self.tree.query(point, k=self.nclose)
        close_ids = self.node_ids[i]
        return close_ids, dists


class KdTree(object):
    """this si the class we use"""
    def __init__(self, tree_type, point_list, nclose=0):
        self.nclose = nclose
        self.tree_type = tree_type
        if tree_type not in ['node', 'element']:
            # verifies you're calling the right
            msg = 'Error!  Invalid tree_type\n'
            msg += "tree_type=%r valid='node','element'" % tree_type
            raise RuntimeError(msg)

        nodes = []
        for node_id, xyz in sorted(iteritems(point_list)):
            n = list(xyz) + [node_id]
            nodes.append(n)
        self.tree = self.build_tree(nodes, nclose)

    def get_close_element_ids(self, point):
        close_nodes_dists = self.nnearest_points(point, self.nclose)
        close_ids = []
        dists = []
        for node_dist in close_nodes_dists:
            ID = node_dist[0][3]
            dist = node_dist[1]
            close_ids.append(ID)  # nid
            dists.append(dist)
        return close_ids, dists

    def build_tree(self, point_list, depth=0):
        if not point_list:
            return None

        # Select axis based on depth so that axis cycles through all valid values
        k = len(point_list[0]) - 1 # Assumes all points have the same dimension, and that
        # the last element in a point is an identifier
        axis = depth % k

        # Sort point list to select median
        point_list.sort(key=lambda x: x[axis])
        median = len(point_list) // 2 # Choose median

        # Create node and construct subtrees
        node = Node()
        node.location = point_list[median]
        node.left_child = self.build_tree(point_list[:median], depth+1)
        node.right_child = self.build_tree(point_list[median+1:], depth+1)
        return node

    def in_circle(self, p, q, radius):
        d = 0.
        for i in range(len(p)):
            d = d + (p[i] - q[i])**2
        if d <= radius**2:
            return True
        return False

    def distance(self, p, q):
        d = 0.
        for i in range(len(p)):
            d += (p[i] - q[i])**2
        return d**0.5

    def points_in_sphere(self, tree, p, radius, depth=0, ptlist=None):
        if ptlist is None:
            ptlist = []
        if depth == 0:
            ptlist = []

        k = len(p)
        axis = depth % k
        distance = self.distance(p, tree.location)
        if distance <= radius:
            ptlist.append((tree.location, distance))
        if tree.left_child:
            if tree.location[axis] >= p[axis] - radius:
                self.points_in_sphere(tree.left_child, p, radius, depth + 1, ptlist)
        if tree.right_child:
            if tree.location[axis] <= p[axis] + radius:
                self.points_in_sphere(tree.right_child, p, radius, depth + 1, ptlist)
        return ptlist

    def nnearest_points(self, p, n):
        tree = self.tree
        ptlist = None
        self.radius_guess = 0.
        npts = 0
        while npts < n:
            if self.radius_guess > 0.:
                self.radius_guess = self.radius_guess * 2.
            else:
                self.radius_guess = 0.001
            ptlist = self.points_in_sphere(tree, p, self.radius_guess)
            npts = len(ptlist)
        ptlist.sort(key=lambda x: x[1])
        self.radius_guess = ptlist[n-1][1]
        return ptlist[:n]
