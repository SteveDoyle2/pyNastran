from pyNastran.bdf.cards.base_card import BaseCard
from pyNastran.bdf.bdf_interface.bdf_card import BDFCard
from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.bdf_interface.utils import to_fields


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
    def _init_from_empty(cls):
        sid = 1
        poly1 = 1
        poly2 = 1
        poly3 = 1
        cid = None
        typei = 'SET'
        idi = 42
        return PSET(sid, poly1, poly2, poly3, cid, typei, idi, comment='')

    @classmethod
    def add_card(cls, card, comment=''):
        sid = integer(card, 1, 'sid')
        poly1 = integer(card, 2, 'poly1')
        poly2 = integer_or_blank(card, 3, 'poly2', poly1)
        poly3 = integer_or_blank(card, 4, 'poly3', poly1)
        cid = integer_or_blank(card, 5, 'cid')
        typei = string_or_blank(card, 6, 'typei', 'SET')
        idi = integer_or_blank(card, 7, 'id', 999999)
        return PSET(sid, poly1, poly2, poly3, cid, typei, idi, comment=comment)

    def raw_fields(self):
        return ['PSET', self.sid, self.poly1, self.poly2, self.poly3, self.cid, self.Type, self.idi]

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
    def _init_from_empty(cls):
        idi = 1
        poly1 = 1
        poly2 = 1
        poly3 = 1
        cid = None
        typei = 'SET'
        typeids = 42
        return PVAL(idi, poly1, poly2, poly3, cid, typei, typeids)

    @classmethod
    def add_card(cls, card, comment=''):
        idi = integer(card, 1, 'idi')
        poly1 = integer(card, 2, 'nid1')
        poly2 = integer_or_blank(card, 3, 'poly2', poly1)
        poly3 = integer_or_blank(card, 4, 'poly3', poly1)
        cid = integer_or_blank(card, 5, 'cid')
        typei = string_or_blank(card, 6, 'typei', 'SET')
        typeids = integer_or_blank(card, 7, 'typeids', 999999)
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
        geom_ids = [idi if idi is not None else 0
                    for  idi in geom_ids]
        #if geom_ids[0] is None:
            #geom_ids[0] = 0
        #if geom_ids[1] is None:
            #geom_ids[1] = 0
        if geom_ids[0] == 0:
            # ID1      ID2      Shape of the FEEDGE
            # =======  =======  ===================
            # Blank/0  Blank/0  Linear
            # >0       Blank/0  Quadratic
            # >0      >0        Cubic
            # Blank/0 >0        Not allowed

            # 0        0        Linear
            # >0       0        Quadratic
            # >0      >0        Cubic
            # 0       >0        Not allowed
            assert geom_ids[1] == 0, f'geom1 must be 0 when geom_id1=0; geom_ids={geom_ids}'

    @classmethod
    def _init_from_empty(cls):
        edge_id = 1
        nids = [1, 2]
        cid = None
        geom_ids = [0, 0]
        return FEEDGE(edge_id, nids, cid, geom_ids, geomin='POINT', comment='')

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
    def _init_from_empty(cls):
        face_id = 1
        nids = [1, 2]
        cid = None
        surf_id = 0
        return FEFACE(face_id, nids, cid, surf_id, comment='')

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
    def __init__(self, curve_id, group, data, cid_in=0, cid_bc=0, comment=''):
        """
        | GMCURV | CURVID | GROUP | CIDIN | CIDBC |
        |Evaluator Specific Data and Format |
        """
        if comment:
            self.comment = comment
        self.curve_id = curve_id
        self.group = group
        self.cid_in = cid_in
        self.cid_bc = cid_bc
        self.data = data
        assert isinstance(data, list), type(data)
        assert group in ['MSCGRP0', 'MSCGRP1', 'MSCGRP2'], group

    @classmethod
    def _init_from_empty(cls):
        curve_id = 1
        group = 'MSCGRP0'
        data = ['CAT']
        return GMCURV(curve_id, group, data, cid_in=0, cid_bc=0, comment='')

    @classmethod
    def add_card(cls, card_lines, comment=''):
        card = BDFCard(to_fields([card_lines[0]], 'GMCURV'))
        curve_id = integer(card, 1, 'curve_id')
        group = string(card, 2, 'group')
        cid_in = integer_or_blank(card, 3, 'cid_in', 0)
        cid_bc = integer_or_blank(card, 4, 'cid_bc', 0)
        data = card_lines[1:]
        return GMCURV(curve_id, group, data, cid_in=cid_in, cid_bc=cid_bc, comment=comment)

    def raw_fields(self):
        return ['GMCURV', self.curve_id, self.group, self.cid_in, self.cid_bc, self.data]

    def write_card(self, size=8, is_double=False):
        #data = self.data.strip()
        #print(repr(data))
        #data_split = ['        %s\n' % data[i:i+64].strip() for i in range(0, len(data), 64)]

        data_split = self.data
        #print(data_split)
        msg = 'GMCURV  %8i%8s%8i%8i\n' % (
            self.curve_id, self.group, self.cid_in, self.cid_bc)
        #msg += ''.join(data_split)
        msg += '\n'.join(data_split) + '\n'
        return msg

    def __repr__(self):
        return self.write_card(size=8, is_double=False)


class GMSURF(BaseCard):
    type = 'GMSURF'
    def __init__(self, surf_id, group, data, cid_in=0, cid_bc=0, comment=''):
        """
        | GMSURF | SURFID | GROUP | CIDIN | CIDBC |
        |Evaluator Specific Data and Format |

        """
        if comment:
            self.comment = comment
        self.surf_id = surf_id
        self.group = group
        self.cid_in = cid_in
        self.cid_bc = cid_bc
        self.data = data
        assert isinstance(data, list), type(data)
        assert group in ['MSCGRP0', 'MSCGRP1', 'MSCGRP2'], group # 'MSCGRP0', 'MSCGRP1', 'MSCGRP2'

    @classmethod
    def _init_from_empty(cls):
        surf_id = 1
        group = 'MSCGRP0'
        data = ['CAT']
        return GMSURF(surf_id, group, data, cid_in=0, cid_bc=0, comment='')

    @classmethod
    def add_card(cls, card_lines, comment=''):
        card = BDFCard(to_fields([card_lines[0]], 'GMSURF'))
        surf_id = integer(card, 1, 'surf_id')
        group = string(card, 2, 'group')
        cid_in = integer_or_blank(card, 3, 'cid_in', 0)
        cid_bc = integer_or_blank(card, 4, 'cid_bc', 0)
        data = card_lines[1:]
        return GMSURF(surf_id, group, data, cid_in=cid_in, cid_bc=cid_bc, comment=comment)

    def raw_fields(self):
        return ['GMSURF', self.surf_id, self.group, self.cid_in, self.cid_bc, self.data]

    def write_card(self, size=8, is_double=False):
        #data = self.data.strip()
        #print(repr(data))
        #data_split = ['        %s\n' % data[i:i+64].strip() for i in range(0, len(data), 64)]

        data_split = self.data
        #print(data_split)
        msg = 'GMSURF  %8i%8s%8i%8i\n' % (
            self.surf_id, self.group, self.cid_in, self.cid_bc)
        #msg += ''.join(data_split)
        msg += '\n'.join(data_split) + '\n'
        return msg

    def __repr__(self):
        return self.write_card(size=8, is_double=False)
