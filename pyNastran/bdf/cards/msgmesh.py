# EQUIV
# DELETE
# ACTIVE
# RSPLNG

# CGEN
# SPCG
# SETG
# TEMPG
# ETEMP
# LIST
# EDGER
# INSECT
# GRIDMOD
# DISTORT

# TRICON
# CBARG
# PGEN

# EIDL
# EIDH

# GRIDMOD
# GRIDD
# GRIDG
# GRIDU
# EGRID

# CURVE
# EDGER

# MPOINT
# MDIM
# MLINE
# MAREA
# MTAB
# MFUN

# MESHOFF
# MESHON

# PLOTG
# PLOTG
# PLOTOPT

# PLOADG

class CGEN(object):
    """
    +------+------+-----------+-----+----------+-----+-------------+------+------+
    |   1  |   2  |     3     | 4   |     5    |  6  |      7      |  8   |  9   |
    +======+======+===========+=====+==========+=====+=============+======+======+
    | CGEN | TYPE | FIELD_EID | PID | FIELD_ID | DIR | TH_GEOM_OPT | EIDL | EIDH |
    +------+------+-----------+-----+----------+-----+-------------+------+------+
    """
    type = 'CGEN'

    def __init__(self, Type, field_eid, pid, field_id, th_geom_opt, eidl, eidh, cdir='L', comment=''):
        """
        Creates the GRID card in a functional way

        Parameters
        ----------
        Type : str
            LINE, TRIA, QUAD, HEX
        field_eid : int
           starting element id
        pid : int
           property id
        field_id : int
            GRIDG id
        cdir : str
            L, M, N
        th_geom_opt :
            TH : ???
                ???
            GEOM : ???
                only when Type = BEND
            OPT : ???
                ???
        eidl : int
            ???
        eidh : int
            ???
        """
        if comment:
            self._comment = comment
        self.Type = Type
        self.field_eid = field_eid
        self.pid = pid
        self.field_id = field_id
        self.cdir = cdir

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Creates the GRID card

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        data : List[int/float]; default=None
            a list with the GRID fields defined in OP2 format
        comment : str; default=''
          a comment for the card
        """
        nfields = len(card)
        #: Node ID
        nid = integer(card, 1, 'nid')

        #: Grid point coordinate system
        cp = integer_or_blank(card, 2, 'cp', 0)

        #: node location in local frame
        xyz = np.array([
            double_or_blank(card, 3, 'x1', 0.),
            double_or_blank(card, 4, 'x2', 0.),
            double_or_blank(card, 5, 'x3', 0.)], dtype='float64')

        if nfields > 6:
            #: Analysis coordinate system
            cd = integer_or_blank(card, 6, 'cd', 0)

            #: SPC constraint
            ps = u(integer_or_blank(card, 7, 'ps', ''))

            #: Superelement ID
            seid = integer_or_blank(card, 8, 'seid', 0)
            assert len(card) <= 9, 'len(GRID card) = %i\ncard=%s' % (len(card), card)
        else:
            cd = 0
            ps = ''
            seid = 0
        return CGEN(nid, cp, xyz, cd, ps, seid, comment=comment)

    def _validate_input(self):
        assert self.nid > 0, 'nid=%s' % (self.nid)
        assert self.cp >= 0, 'cp=%s' % (self.cp)
        assert self.cd >= -1, 'cd=%s' % (self.cd)
        assert self.seid >= 0, 'seid=%s' % (self.seid)
        assert len(self.xyz) == 3



    def _verify(self, xref):
        """
        Verifies all methods for this object work

        Parameters
        ----------
        xref : bool
            has this model been cross referenced
        """
        pass

    def cross_reference(self, model, grdset=None):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        grdset : GRDSET / None; default=None
            a GRDSET if available (default=None)

        .. note::  The gridset object will only update the fields that
                   have not been set
        """
        if grdset:
            # update using a gridset object
            if not self.cp:
                self.cp = grdset.cp
                self.cp_ref = self.cp
            if not self.cd:
                self.cd = grdset.cd
                self.cd_ref = self.cd
            if not self.ps:
                self.ps = grdset.ps
                self.ps_ref = self.ps
            if not self.seid:
                self.seid = grdset.seid
                self.seid_ref = self.seid
        msg = ' which is required by %s nid=%s' % (self.type, self.nid)
        self.cp = model.Coord(self.cp, msg=msg)
        self.cp_ref = self.cp
        if self.cd != -1:
            self.cd = model.Coord(self.cd, msg=msg)
            self.cd_ref = self.cd

    def uncross_reference(self):
        self.cp = self.Cp()
        self.cd = self.Cd()
        del self.cp_ref, self.cd_ref
        if hasattr(self, 'elements'):
            del self.elements

    def raw_fields(self):
        """
        Gets the fields in their unmodified form

        Returns
        -------
        fields : List[int/float/str]
            the fields that define the card
        """
        list_fields = ['GRID', self.nid, self.Cp()] + list(self.xyz) + \
                      [self.Cd(), self.ps, self.SEid()]
        return list_fields

    def repr_fields(self):
        """
        Gets the fields in their simplified form

        Returns
        -------
        fields : List[int/float/str]
            the fields that define the card
        """
        cp = set_blank_if_default(self.Cp(), 0)
        cd = set_blank_if_default(self.Cd(), 0)
        seid = set_blank_if_default(self.SEid(), 0)
        list_fields = ['CGEN', self.nid, cp] + list(self.xyz) + [cd, self.ps,
                                                                 seid]
        return list_fields

    #def write_card(self, size=8, is_double=False):
        #"""
        #The writer method used by BDF.write_card

        #Parameters
        #----------
        #size : int; default=8
            #the size of the card (8/16)
        #is_double : bool; default=False
            #should this card be written with double precision

        #Returns
        #-------
        #msg : str
            #the card as a string
        #"""
        #if size == 8:
            #return self.write_card_8()
        #else:
            #return self.write_card_16(is_double)
