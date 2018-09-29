"""
 - BCRPARA
 - BCTADD
 - BCTSET
 - BSURF
 - BSURFS
"""
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)

from pyNastran.bdf.cards.base_card import BaseCard, expand_thru_by
from pyNastran.bdf.bdf_interface.assign_type import (
    integer, integer_or_blank, integer_string_or_blank, double_or_blank,
    integer_double_or_blank, string, string_or_blank)
from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.field_writer_16 import print_card_16

class BLSEG(BaseCard):
    """
    3D Contact Region Definition by Shell Elements (SOLs 101, 601 and 701)

    Defines a 3D contact region by shell element IDs.

    BLSEG ID G1 G2 G3 G4 G5 G6 G7
    BLSEG ID G1 THRU G2 BY INC
    """
    type = 'BLSEG'
    def __init__(self, line_id, nodes, comment=''):
        if comment:
            self.comment = comment

        self.line_id = line_id
        self.nodes = expand_thru_by(nodes)

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a BLSEG card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        line_id = integer(card, 1, 'line_id')
        #: Number (float)
        nfields = card.nfields
        i = 2
        nodes = []
        while i < nfields:
            d = integer_string_or_blank(card, i, 'field_%s' % i)
            if d is not None:
                nodes.append(d)
            i += 1
        return BLSEG(line_id, nodes, comment=comment)

    #def raw_fields(self):
        #fields = ['BSURF', self.sid]
        #return fields + list(self.eids)

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        return self.comment + print_card_8(card)

class BCONP(BaseCard):
    """
    3D Contact Region Definition by Shell Elements (SOLs 101, 601 and 701)

    Defines a 3D contact region by shell element IDs.

    +-------+----+-------+--------+-----+------+--------+-------+-----+
    |   1   |  2 |   3   |   4    |  5  |   6  |   7    |   8   |  9  |
    +=======+====+=======+========+=====+======+========+=======+=====+
    | BCONP | ID | SLAVE | MASTER |     | SFAC | FRICID | PTYPE | CID |
    +-------+----+-------+--------+-----+------+--------+-------+-----+
    | BCONP | 95 |   10  |   15   |     |  1.0 |   33   |   1   |     |
    +-------+----+-------+--------+-----+------+--------+-------+-----+

    """
    type = 'BCONP'
    def __init__(self, contact_id, slave, master, sfac, fric_id, ptype, cid, comment=''):
        if comment:
            self.comment = comment

        self.contact_id = contact_id
        self.slave = slave
        self.master = master
        self.sfac = sfac
        self.fric_id = fric_id
        self.ptype = ptype
        self.cid = cid

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a BCONP card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        contact_id = integer(card, 1, 'contact_id')
        slave = integer(card, 2, 'slave')
        master = integer(card, 3, 'master')
        sfac = double_or_blank(card, 5, 'sfac', 1.0)
        fric_id = integer_or_blank(card, 6, 'fric_id')
        ptype = integer_or_blank(card, 7, 'ptype', 1)
        cid = integer_or_blank(card, 8, 'cid', 0)
        return BCONP(contact_id, slave, master, sfac, fric_id, ptype, cid, comment=comment)

    #def raw_fields(self):
        #fields = ['BSURF', self.sid]
        #return fields + list(self.eids)

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        return self.comment + print_card_8(card)

class BSURF(BaseCard):
    """
    3D Contact Region Definition by Shell Elements (SOLs 101, 601 and 701)

    Defines a 3D contact region by shell element IDs.

    +-------+------+------+-------+-------+--------+------+------+------+
    |   1   |   2  |   3  |   4   |   5   |    6   |  7   |   8  |   9  |
    +=======+======+======+=======+=======+========+======+======+======+
    | BSURF |  ID  | EID1 | EID2  | EID3  |  EID4  | EID5 | EID6 | EID7 |
    +-------+------+------+-------+-------+--------+------+------+------+
    |       | EID8 | EID9 | EID10 | -etc- |        |      |      |      |
    +-------+------+------+-------+-------+--------+------+------+------+

    +-------+------+------+-------+-------+--------+------+------+------+
    | BSURF |  ID  | EID1 |  THRU | EID2  |   BY   | INC  |      |      |
    +-------+------+------+-------+-------+--------+------+------+------+
    |       | EID8 | EID9 | EID10 | EID11 | -etc.- |      |      |      |
    +-------+------+------+-------+-------+--------+------+------+------+
    |       | EID8 | THRU | EID9  |  BY   |  INC   |      |      |      |
    +-------+------+------+-------+-------+--------+------+------+------+

    +-------+------+------+-------+-------+--------+------+------+------+
    | BSURF |  15  |  5   | THRU  |  21   |   BY   |  4   |      |      |
    +-------+------+------+-------+-------+--------+------+------+------+
    |       |  27  |  30  |  32   |  33   |        |      |      |      |
    +-------+------+------+-------+-------+--------+------+------+------+
    |       |  35  | THRU |  44   |       |        |      |      |      |
    +-------+------+------+-------+-------+--------+------+------+------+
    |       |  67  |  68  |  70   |  85   |   92   |      |      |      |
    +-------+------+------+-------+-------+--------+------+------+------+
    """
    type = 'BSURF'
    def __init__(self, sid, eids, comment=''):
        if comment:
            self.comment = comment
        #: Set identification number. (Unique Integer > 0)
        self.sid = sid
        #: Element identification numbers of shell elements. (Integer > 0)
        self.eids = eids

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a BSURF card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        sid = integer(card, 1, 'sid')
        #: Number (float)
        nfields = card.nfields
        i = 2
        eid_data = []
        while i < nfields:
            d = integer_string_or_blank(card, i, 'field_%s' % i)
            if d is not None:
                eid_data.append(d)
            i += 1
        eids = expand_thru_by(eid_data)
        return BSURF(sid, eids, comment=comment)

    def raw_fields(self):
        fields = ['BSURF', self.sid]
        return fields + list(self.eids)

        #fields = ['BSURF', self.sid, None, None, None, None, None, None, None]
        ## is this right???
        #packs = collapse_thru_by(self.eids, get_packs=True)

        #pack = packs[0]
        #if len(pack) == 3:
            #minv, maxv, dv = pack
            #if dv == 1:
                #fields[2:5] = [minv, 'THRU', maxv]
            #else:
                #fields[2:7] = [minv, 'THRU', maxv, 'BY', dv]
        #else:
            #fields[3:3+len(pack)] = pack

        #for pack in packs[1:]:
            ##fields += pack + [None, None, None]
            #if len(pack) == 3:
                #minv, maxv, dv = pack
                #if dv == 1:
                    #fields += [minv, 'THRU', maxv, None, None, None, None]
                #else:
                    #fields += [minv, 'THRU', maxv, 'BY', dv, None, None]
            #else:
                #fields += pack + [None] * (8 - len(pack))
        ##for sid, tid, fric, mind, maxd in zip(self.sids, self.tids, self.frictions,
        ##                                      self.min_distances, self.max_distances):
        ##    fields += [sid, tid, fric, mind, maxd, None, None]
        #return fields

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class BSURFS(BaseCard):
    """
    Defines a 3D contact region by the faces of the CHEXA, CPENTA or CTETRA
    elements.

    Notes
    -----
    1. The continuation field is optional.
    2. BSURFS is a collection of one or more element faces on solid elements.
       BSURFS defines a contact region which may act as a contact source
       (contactor) or target.
    3. The ID must be unique with respect to all other BSURFS, BSURF, and
       BCPROP entries.
    """
    type = 'BSURFS'
    def __init__(self, id, eids, g1s, g2s, g3s, comment=''):
        if comment:
            self.comment = comment
        #: Identification number of a contact region. See Remarks 2 and 3.
        #: (Integer > 0)
        self.id = id

        #: Element identification numbers of solid elements. (Integer > 0)
        self.eids = eids

        #: Identification numbers of 3 corner grid points on the face (triangular
        #: or quadrilateral) of the solid element. (Integer > 0)
        self.g1s = g1s
        self.g2s = g2s
        self.g3s = g3s

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a BSURFS card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        id = integer(card, 1, 'id')
        eids = []
        g1s = []
        g2s = []
        g3s = []

        n = card.nfields - 5
        i = 0
        j = 1
        while i < n:
            eid = integer(card, 5 + i, 'eid%s' % j)
            g1 = integer(card, 5 + i + 1, 'g3_%s' % j)
            g2 = integer(card, 5 + i + 2, 'g2_%s' % j)
            g3 = integer(card, 5 + i + 3, 'g1_%s' % j)
            j += 1
            i += 4
            eids.append(eid)
            g1s.append(g1)
            g2s.append(g2)
            g3s.append(g3)
        return BSURFS(id, eids, g1s, g2s, g3s, comment=comment)

    def raw_fields(self):
        fields = ['BSURFS', self.id, None, None, None]
        for eid, g1, g2, g3 in zip(self.eids, self.g1s, self.g2s, self.g3s):
            fields += [eid, g1, g2, g3]
        return fields

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class BCTSET(BaseCard):
    """
    3D Contact Set Definition (SOLs 101, 601 and 701 only)
    Defines contact pairs of a 3D contact set.

    +--------+-------+------+-------+-------+-------+-------+
    |   1    |   2   | 3    |  4    |   5   |   6   |   7   |
    +========+=======+======+=======+=======+=======+=======+
    | BCTSET | CSID  | SID1 | TID1  | FRIC1 | MIND1 | MAXD1 |
    +--------+-------+------+-------+-------+-------+-------+
    |        | SID2  | TID2 | FRIC2 | MIND2 | MAXD2 |       |
    +--------+-------+------+-------+-------+-------+-------+
    |        | -etc- |      |       |       |       |       |
    +--------+-------+------+-------+-------+-------+-------+
    """
    type = 'BCTSET'
    def __init__(self, csid, sids, tids, frictions, min_distances, max_distances,
                 comment='', sol=101):
        if comment:
            self.comment = comment
        #: CSID Contact set identification number. (Integer > 0)
        self.csid = csid
        #: SIDi Source region (contactor) identification number for contact pair i.
        #: (Integer > 0)
        self.sids = sids

        #: TIDi Target region identification number for contact pair i. (Integer > 0)
        self.tids = tids

        #: FRICi Static coefficient of friction for contact pair i. (Real; Default=0.0)
        self.frictions = frictions

        #: MINDi Minimum search distance for contact. (Real) (Sol 101 only)
        self.min_distances = min_distances

        #: MAXDi Maximum search distance for contact. (Real) (Sol 101 only)
        self.max_distances = max_distances

    @classmethod
    def add_card(cls, card, comment='', sol=101):
        csid = integer(card, 1, 'csid')
        sids = []
        tids = []
        frictions = []
        min_distances = []
        max_distances = []

        nfields = card.nfields
        i = 2
        j = 1
        while i < nfields:
            sids.append(integer(card, i, 'sid%s' % j))
            tids.append(integer(card, i + 1, 'tid%s' % j))
            frictions.append(double_or_blank(card, i + 2, 'fric%s' % j, 0.0))
            if sol == 101:
                min_distances.append(double_or_blank(card, i + 3, 'mind%s' % j, 0.0))
                max_distances.append(double_or_blank(card, i + 4, 'maxd%s' % j, 0.0))
            else:
                min_distances.append(None)
                max_distances.append(None)
            i += 8
            j += 1
        return BCTSET(csid, sids, tids, frictions, min_distances,
                      max_distances, comment=comment,
                      sol=sol)

    def raw_fields(self):
        fields = ['BCTSET', self.csid]
        for sid, tid, fric, mind, maxd in zip(self.sids, self.tids, self.frictions,
                                              self.min_distances, self.max_distances):
            fields += [sid, tid, fric, mind, maxd, None, None, None]
        return fields

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)


class BCRPARA(BaseCard):
    """
    +---------+------+------+--------+------+-----+---+---+---+----+
    |    1    |   2  |   3  |   4    |   5  |  6  | 7 | 8 | 9 | 10 |
    +=========+======+======+========+======+=====+===+===+===+====+
    | BCRPARA | CRID | SURF | OFFSET | TYPE | GP  |   |   |   |    |
    +---------+------+------+--------+------+-----+---+---+---+----+
    """
    type = 'BCRPARA'
    def __init__(self, crid, offset=None, surf='TOP', Type='FLEX', grid_point=0,
                 comment=''):
        """
        Creates a BCRPARA card

        Parameters
        ----------
        crid : int
            CRID Contact region ID.
        offset : float; default=None
            Offset distance for the contact region (Real > 0.0).
            None : OFFSET value in BCTPARA entry
        surf : str; default='TOP'
            SURF Indicates the contact side. See Remark 1.  {'TOP', 'BOT'; )
        Type : str; default='FLEX'
            Indicates whether a contact region is a rigid surface if it
            is used as a target region. {'RIGID', 'FLEX'}.
            This is not supported for SOL 101.
        grid_point : int; default=0
            Control grid point for a target contact region with TYPE=RIGID
            or when the rigid-target algorithm is used.  The grid point
            may be used to control the motion of a rigid surface.
            (Integer > 0).  This is not supported for SOL 101.
        comment : str; default=''
            a comment for the card

        """
        if comment:
            self.comment = comment

        #: CRID Contact region ID. (Integer > 0)
        self.crid = crid

        #: SURF Indicates the contact side. See Remark 1. (Character = "TOP" or
        #: "BOT"; Default = "TOP")
        self.surf = surf

        #: Offset distance for the contact region. See Remark 2. (Real > 0.0,
        #: Default =OFFSET value in BCTPARA entry)
        self.offset = offset

        #: Indicates whether a contact region is a rigid surface if it is used as a
        #: target region. See Remarks 3 and 4. (Character = "RIGID" or "FLEX",
        #: Default = "FLEX"). This is not supported for SOL 101.
        self.Type = Type

        #: Control grid point for a target contact region with TYPE=RIGID or
        #: when the rigid-target algorithm is used. The grid point may be
        #: used to control the motion of a rigid surface. (Integer > 0)
        #: This is not supported for SOL 101.
        self.grid_point = grid_point

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a BCRPARA card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        crid = integer(card, 1, 'crid')
        surf = string_or_blank(card, 2, 'surf', 'TOP')
        offset = double_or_blank(card, 3, 'offset', None)
        Type = string_or_blank(card, 4, 'type', 'FLEX')
        grid_point = integer_or_blank(card, 5, 'grid_point', 0)
        return BCRPARA(crid, surf=surf, offset=offset, Type=Type,
                       grid_point=grid_point, comment=comment)

    def raw_fields(self):
        fields = ['BCRPARA', self.crid, self.surf, self.offset, self.Type, self.grid_point]
        return fields

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)


class BCTPARA(BaseCard):
    """
    Defines parameters for a surface-to-surface contact region.

    +---------+--------+--------+--------+--------+--------+---------+--------+
    |    1    |   2    |    3   |   4    |   5    |   6    |    7    |    8   |
    +=========+========+========+========+========+========+=========+========+
    | BCTPARA | CSID   | Param1 | Value1 | Param2 | Value2 | Param3  | Value3 |
    +---------+--------+--------+--------+--------+--------+---------+--------+
    |         | Param4 | Value4 | Param5 | Value5 | -etc-  |         |        |
    +---------+--------+--------+--------+--------+--------+---------+--------+
    | BCTPARA |   1    | TYPE   |   0    | NSIDE  |   2    | SEGNORM |  -1    |
    +---------+--------+--------+--------+--------+--------+---------+--------+
    |         | CSTIFF |   1    | OFFSET | 0.015  |        |         |        |
    +---------+--------+--------+--------+--------+--------+---------+--------+
    """
    type = 'BCTPARA'
    def __init__(self, csid, params, comment=''):
        """
        Creates a BCTPARA card

        Parameters
        ----------
        csid : int
            Contact set ID. Parameters defined in this command apply to
            contact set CSID defined by a BCTSET entry. (Integer > 0)
        params : dict[key] : value
            the optional parameters
        comment : str; default=''
            a comment for the card

        """
        if comment:
            self.comment = comment

        #: Contact set ID. Parameters defined in this command apply to
        #: contact set CSID defined by a BCTSET entry. (Integer > 0)
        self.csid = csid
        self.params = params

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a BCTPARA card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        csid = integer(card, 1, 'csid')
        i = 2
        j = 1
        params = {}
        while i < card.nfields:
            param = string(card, i, 'param%s' % j)
            i += 1
            if param == 'TYPE':
                value = integer_or_blank(card, i, 'value%s' % j, 0)
                assert value in [0, 1, 2], 'TYPE must be [0, 1, 2]; TYPE=%r' % value
            #elif param == 'TYPE': # NX
                #value = string_or_blank(card, i, 'value%s' % j, 'FLEX').upper()
                #assert value in ['FLEX', 'RIGID', 'COATING'], 'TYPE must be [FLEX, RIGID, COATING.]; CSTIFF=%r' % value

            elif param == 'NSIDE':
                value = integer_or_blank(card, i, 'value%s' % j, 1)
                assert value in [1, 2], 'NSIDE must be [1, 2]; NSIDE=%r' % value
            elif param == 'TBIRTH':
                value = double_or_blank(card, i, 'value%s' % j, 0.0)
            elif param == 'TDEATH':
                value = double_or_blank(card, i, 'value%s' % j, 0.0)
            elif param == 'INIPENE':
                value = integer_or_blank(card, i, 'value%s' % j, 0)
                assert value in [0, 1, 2, 3], 'INIPENE must be [0, 1, 2]; INIPENE=%r' % value
            elif param == 'PDEPTH':
                value = double_or_blank(card, i, 'value%s' % j, 0.0)
            elif param == 'SEGNORM':
                value = integer_or_blank(card, i, 'value%s' % j, 0)
                assert value in [-1, 0, 1], 'SEGNORM must be [-1, 0, 1]; SEGNORM=%r' % value
            elif param == 'OFFTYPE':
                value = integer_or_blank(card, i, 'value%s' % j, 0)
                assert value in [0, 1, 2], 'OFFTYPE must be [0, 1, 2]; OFFTYPE=%r' % value
            elif param == 'OFFSET':
                value = double_or_blank(card, i, 'value%s' % j, 0.0)
            elif param == 'TZPENE':
                value = double_or_blank(card, i, 'value%s' % j, 0.0)

            elif param == 'CSTIFF':
                value = integer_or_blank(card, i, 'value%s' % j, 0)
                assert value in [0, 1], 'CSTIFF must be [0, 1]; CSTIFF=%r' % value
            elif param == 'TIED':
                value = integer_or_blank(card, i, 'value%s' % j, 0)
                assert value in [0, 1], 'TIED must be [0, 1]; TIED=%r' % value
            elif param == 'TIEDTOL':
                value = double_or_blank(card, i, 'value%s' % j, 0.0)
            elif param == 'EXTFAC':
                value = double_or_blank(card, i, 'value%s' % j, 0.001)
                assert 1.0E-6 <= value <= 0.1, 'EXTFAC must be 1.0E-6 < EXTFAC < 0.1; EXTFAC=%r' % value
            else:
                # FRICMOD, FPARA1/2/3/4/5, EPSN, EPST, CFACTOR1, PENETOL
                # NCMOD, TCMOD, RFORCE, LFORCE, RTPCHECK, RTPMAX, XTYPE
                # ...
                value = integer_double_or_blank(card, i, 'value%s' % j)
                assert value is not None, '%s%i must not be None' % (param, j)

            params[param] = value
            i += 1
            j += 1
            if j == 4:
                i += 1
        return BCTPARA(csid, params, comment=comment)

    def raw_fields(self):
        fields = ['BCTPARA', self.csid]
        i = 0
        for key, value in sorted(self.params.items()):
            if i == 3:
                fields.append(None)
                i = 0
            fields.append(key)
            fields.append(value)
            i += 1
        return fields

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)


class BCTADD(BaseCard):
    """
    +--------+------+----+-------+----+----+----+----+----+
    |   1    |  2   | 3  |   4   |  5 | 6  |  7 | 8  |  9 |
    +========+======+====+=======+====+====+====+====+====+
    | BCTADD | CSID | SI |  S2   | S3 | S4 | S5 | S6 | S7 |
    +--------+------+----+-------+----+----+----+----+----+
    |        |  S8  | S9 |  etc. |    |    |    |    |    |
    +--------+------+----+-------+----+----+----+----+----+

    Remarks:
    1. To include several contact sets defined via BCTSET entries in a model,
       BCTADD must be used to combine the contact sets. CSID in BCTADD is
       then selected with the Case Control command BCSET.
    2. Si must be unique and may not be the identification of this or any other
       BCTADD entry.
    """
    type = 'BCTADD'
    def __init__(self, csid, contact_sets, comment=''):
        if comment:
            self.comment = comment
        #: Contact set identification number. (Integer > 0)
        self.csid = csid

        #: Identification numbers of contact sets defined via BCTSET entries.
        #: (Integer > 0)
        self.contact_sets = contact_sets

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a BCTADD card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        csid = integer(card, 1, 'csid')
        contact_sets = []

        i = 1
        j = 1
        while i < card.nfields:
            contact_set = integer(card, i, 'S%i' % j)
            contact_sets.append(contact_set)
            i += 1
            j += 1
        return BCTADD(csid, contact_sets, comment=comment)

    def raw_fields(self):
        fields = ['BCTADD'] + self.contact_sets
        return fields

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)
