from pyNastran.bdf.cards.base_card import BaseCard

class FEEDGE(BaseCard):
    type = 'FEEDGE'
    def __init__(self, edge_id, nids, cid, geomin, geom_ids, comment=''):
        if comment:
            self.comment = comment
        self.edge_id = edge_id
        self.nids = nids
        self.cid = cid
        self.geomin = geomin
        self.geom_ids = geom_ids
        assert geomin in ['POINT', 'GMCURV'], geomin

    def raw_fields(self):
        return ['FEEDGE', self.edge_id] + self.nids + [self.cid, self.geomin] + self.geom_ids

class FEFACE(BaseCard):
    type = 'FEFACE'
    def __init__(self, face_id, nids, cid, surf_ids, comment=''):
        if comment:
            self.comment = comment
        self.face_id = face_id
        self.nids = nids
        self.cid = cid
        self.surf_ids = surf_ids

    def raw_fields(self):
        return ['FEFACE', self.face_id] + self.nids + [self.cid] + self.surf_ids
