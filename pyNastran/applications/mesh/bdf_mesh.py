from pyNastran.bdf.bdf import BDF
from collections import defaultdict

class MeshTools(BDF):
    def __init__(self, log=None, debug=False):
        BDF.__init__(self, debug=debug, log=log)

    def get_elemental_nodeids(self):
        nid_to_eids = defaultdict(list)
        eid_to_nids = {}
        for eid, element in self.elements.iteritems():
            nids = element.nodeIDs()
            eid_to_nids[eid] = nids
            for nid in nids:
                nid_to_eids[nid].append(eid)

        for eid, element in self.rigidElements.iteritems():
            nids = element.nodeIDs()
            #print nids
            eid_to_nids[eid] = nids
            for nid in nids:
                nid_to_eids[nid].append(eid)

        #print('mpc', self.mpcObject)
        for (eid, mpcs) in sorted(self.mpcs.iteritems()):
            for mpc in mpcs:
                #print str(mpc)
                nids = mpc.nodeIDs()
                eid_to_nids[eid] = nids
                for nid in nids:
                    nid_to_eids[nid].append(eid)
        return nid_to_eids, eid_to_nids

    def get_free_nodes(self, eids):
        """
        ===================================
        Example 1 - rods
        ===================================
        1        2         3
        *---11---*----22---*
        nids_to_eids = {
            1 : [11],
            2 : [11, 22],
            3 : [22],
        }
        eid_to_nids = {
            11 : [1, 2],
            22 : [2, 3],
        }
        eids = [11]
        all_used_nids = {1, 2}
        all_used_touching_eids = {11, 22}
        free_nodes = {2}

        ===================================
        Example 2 - quads
        ===================================

        1     3      5     7
        *-----*------*-----*
        |     |      |     |
        | 11  |  22  |  33 |
        |     |      |     |
        *-----*------*-----*
        2     4      6     8

        nids_to_eids = {
            1 : [11],
            2 : [11],
            3 : [11, 22],
            4 : [11, 22],
            5 : [22, 33],
            6 : [22, 33],
            7 : [33],
            8 : [33],
        }
        eid_to_nids = {
            11 : [1, 3, 4, 2],
            22 : [3, 5, 6, 4],
            33 : [5, 7, 8, 6],
        }

        eids = [11]
        used_nids = {1, 2, 3, 4}
        used_touching_eids = {11, 22}
        free_nids = {3, 4}

        eids = [22]
        used_nids = {3, 5, 6, 4}
        used_touching_eids = {11, 22, 33}
        free_nids = {1, 2, 7, 8}
        """
        nid_to_eids, eid_to_nids = self.get_elemental_nodeids()
        used_nids = set([])
        for eid in eids:
            #elem = self.elements[eid]
            e_nids = eid_to_nids[eid]
            #used_nids.add(e_nids)
            used_nids.update(e_nids)

        used_free_eids = set([])
        for nid in used_nids:
            n_eids = nid_to_eids[nid]
            used_free_eids.update(n_eids)

        used_free_nids = set([])
        for eid in used_free_eids:
            e_nids = eid_to_nids[eid]
            used_free_nids.update(e_nids)

        free_nids = used_free_nids.difference(used_nids)
        free_eids = used_free_eids.difference(eids)
        #all_eids = set(self.elements.keys())
        return used_free_eids, used_free_nids, free_eids, free_nids

