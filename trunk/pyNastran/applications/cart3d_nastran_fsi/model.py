class Model(object):
    def __init__(self):
        pass

    def getNodeIDLocations(self,nIDs):
        nodes = []
        for nid in nIDs:
            node = self.Node(nid)
            nodes.append(node)
        return nodes

    def get_element_properties(self, eid):
        raise NotImplementedError('overwrite this method...')

    def NodeIDs(self):
        raise NotImplementedError('overwrite this method...')

    def Node(self,nid):
        raise NotImplementedError('overwrite this method...')
