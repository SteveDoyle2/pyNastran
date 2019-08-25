from pyNastran.bdf.cards.base_card import BaseCard
from pyNastran.bdf.field_writer_8 import print_card_8

from pyNastran.bdf.bdf_interface.assign_type import (
    integer, integer_or_blank, # integer_or_string,
    #parse_components, components_or_blank as fcomponents_or_blank,
    string, string_or_blank, # integer_string_or_blank,
)

class PSET(BaseCard):
    type = 'PSET'
    """
    PSET SID POLY1 POLY2 POLY3 CID SETTYP ID
    PSET 127 1     2     1                12
    """
    def __init__(self, sid, poly1, poly2, poly3, cid, typei, idi, comment=''):
        if comment:
            self.comment = comment
        self.sid = sid
        self.poly1 = poly1
        self.poly2 = poly2
        self.poly3 = poly3
        self.cid = cid
        self.Type = typei
        #self.typeids = typeids
        self.idi = idi

    @classmethod
    def add_card(cls, card, comment=''):
        sid = integer(card, 1, 'sid')
        poly1 = integer(card, 2, 'poly1')
        poly2 = integer(card, 3, 'poly2')
        poly3 = integer(card, 4, 'poly3')
        cid = integer_or_blank(card, 5, 'cid')
        typei = string_or_blank(card, 6, 'typei', 'SET')
        idi = integer_or_blank(card, 7, 'id', 999999)
        return PSET(sid, poly1, poly2, poly3, cid, typei, idi, comment=comment)

    def raw_fields(self):
        return ['PVAL', self.sid, self.poly1, self.poly2, self.poly3, self.cid, self.Type, self.idi]

    def write_card(self, size=8, is_double=False):
        fields = self.raw_fields()
        return print_card_8(fields)


class PVAL(BaseCard):
    type = 'PVAL'
    """
    PVAL ID  POLY1 POLY2 POLY3 CID SETTYP ID
    PVAL 127 1     2     1
    """
    def __init__(self, idi, poly1, poly2, poly3, cid, typei, typeids, comment=''):
        if comment:
            self.comment = comment
        if isinstance(typeids, int):
            typeids = [typeids]
        self.idi = idi
        self.poly1 = poly1
        self.poly2 = poly2
        self.poly3 = poly3
        self.cid = cid
        self.Type = typei
        self.typeids = typeids

    @classmethod
    def add_card(cls, card, comment=''):
        idi = integer(card, 1, 'idi')
        poly1 = integer(card, 2, 'nid1')
        poly2 = integer(card, 3, 'poly2')
        poly3 = integer(card, 4, 'poly3')
        cid = integer_or_blank(card, 5, 'cid')
        typei = string_or_blank(card, 6, 'typei', 'SET')
        typeids = integer(card, 7, 'typeids') # 999999
        return PVAL(idi, poly1, poly2, poly3, cid, typei, typeids, comment=comment)

    def raw_fields(self):
        return ['PVAL', self.idi, self.poly1, self.poly2, self.poly3, self.cid, self.Type] + self.typeids

    def write_card(self, size=8, is_double=False):
        fields = self.raw_fields()
        return print_card_8(fields)


class FEEDGE(BaseCard):
    type = 'FEEDGE'
    def __init__(self, edge_id, nids, cid, geom_ids, geomin='POINT', comment=''):
        if comment:
            self.comment = comment
        self.edge_id = edge_id
        self.nids = nids
        self.cid = cid
        self.geomin = geomin
        self.geom_ids = geom_ids
        assert geomin in ['POINT', 'GMCURV'], geomin
        if geom_ids[0] == 0:
            assert geom_ids[1] is None, f'geom1 must be 0 when geom_id1=0; geom_ids={geom_ids}'

    @classmethod
    def add_card(cls, card, comment=''):
        edge_id = integer(card, 1, 'edge_id')
        nid1 = integer(card, 2, 'nid1')
        nid2 = integer(card, 3, 'nid2')
        cid = integer_or_blank(card, 4, 'cid')
        geomin = string_or_blank(card, 5, 'geomin', 'POINT')
        geom_id1 = integer_or_blank(card, 6, 'geom_id1', 0)
        geom_id2 = integer_or_blank(card, 7, 'geom_id2', None)

        nids = [nid1, nid2]
        geom_ids = [geom_id1, geom_id2]
        return FEEDGE(edge_id, nids, cid, geom_ids, geomin=geomin, comment=comment)

    def raw_fields(self):
        return ['FEEDGE', self.edge_id] + self.nids + [self.cid, self.geomin] + self.geom_ids

    def write_card(self, size=8, is_double=False):
        fields = self.raw_fields()
        return print_card_8(fields)


class FEFACE(BaseCard):
    type = 'FEFACE'
    """
    FEFACE FACEID GRID1 GRID2 GRID3 GRID4 CIDBC SURFID
    FEFACE 101    123   547   243   295   12
    """
    def __init__(self, face_id, nids, cid, surf_ids, comment=''):
        if comment:
            self.comment = comment
        if isinstance(surf_ids, int):
            surf_ids = [surf_ids]
        self.face_id = face_id
        self.nids = nids
        self.cid = cid
        self.surf_ids = surf_ids

    @classmethod
    def add_card(cls, card, comment=''):
        face_id = integer(card, 1, 'face_id')
        nid1 = integer(card, 2, 'nid1')
        nid2 = integer(card, 3, 'nid2')
        nid3 = integer(card, 4, 'nid3')
        nid4 = integer_or_blank(card, 5, 'nid4')
        cid_bc = integer_or_blank(card, 6, 'cid_bc')
        surf_id = integer_or_blank(card, 7, 'surf_id', 0)

        nids = [nid1, nid2, nid3, nid4]
        return FEFACE(face_id, nids, cid_bc, surf_id, comment=comment)

    def raw_fields(self):
        return ['FEFACE', self.face_id] + self.nids + [self.cid] + self.surf_ids

    def write_card(self, size=8, is_double=False):
        fields = self.raw_fields()
        return print_card_8(fields)


class GMCURV(BaseCard):
    type = 'GMCURV'
    def __init__(self, curve_id, group, cid_in, cid_bc, data, comment=''):
        if comment:
            self.comment = comment
        self.curve_id = curve_id
        self.group = group
        self.cid_in = cid_in
        self.cid_bc = cid_bc
        self.data = data
        assert group in ['MSCGRP0'], group

    @classmethod
    def add_card(cls, card_lines, comment=''):
        #print(card_lines)
        from pyNastran.bdf.bdf_interface.utils import to_fields
        from pyNastran.bdf.bdf_interface.bdf_card import BDFCard

        card = BDFCard(to_fields([card_lines[0]], 'GMCURV'))
        curve_id = integer(card, 1, 'curve_id')
        group = string(card, 2, 'group')
        cid_in = integer(card, 3, 'cid_in')
        cid_bc = integer(card, 4, 'cid_bc')
        #print(card.fields())
        data = card_lines[1:]
        return GMCURV(curve_id, group, cid_in, cid_bc, data, comment=comment)

    def raw_fields(self):
        return ['GMCURV', self.curve_id, self.group, self.cid_in, self.cid_bc, self.data]

    def write_card(self, size=8, is_double=False):
        data = self.data.strip()
        #print(repr(data))
        data_split = ['        %s\n' % data[i:i+64].strip() for i in range(0, len(data), 64)]
        #print(data_split)
        msg = 'GMCURV  %8i%8s%8i%8i\n' % (
            self.curve_id, self.group, self.cid_in, self.cid_bc)
        msg += ''.join(data_split)
        print(msg)
        return msg
