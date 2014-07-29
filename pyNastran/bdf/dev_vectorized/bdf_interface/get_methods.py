class GetMethods(object):
    def __init__(self):
        pass

    def Coord(self, cid, msg=''):
        try:
            return self.coords[cid]
        except KeyError:
            raise KeyError('cid=%s not found%s.  Allowed Cids=%s'
                           % (cid, msg, self.coordIDs()))
